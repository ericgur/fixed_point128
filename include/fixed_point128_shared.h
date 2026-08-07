/***********************************************************************************
    MIT License

    Copyright (c) 2025 Eric Gur (ericgur@iname.com)

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.
************************************************************************************/
#ifndef FP128_FIXED_POINT128_SHARED_H
#define FP128_FIXED_POINT128_SHARED_H

/**
 * @file fixed_point128_shared.h
 * @brief Shared definitions, intrinsics, and helper functions for 128-bit fixed-point arithmetic.
 *
 * Provides compiler-specific intrinsic wrappers (MSVC and GCC/Clang), bit
 * manipulation utilities (128-bit shifts, leading zero counts, population counts),
 * multi-word division algorithms (32-bit and 64-bit word-based, derived from
 * Hacker's Delight), IEEE 754 bit-layout views, and build configuration macros.
 *
 * This header is consumed by fixed_point128.h and should not be included directly.
 */

#include <array>
#include <bit>  // countl_zero, popcount - the constant evaluated stand-ins for the lzcnt/popcnt intrinsics
#include <cstdint>
#include <cassert>
#include <cctype>   // tolower, isspace
#include <cstdio>   // snprintf
#include <cstring>  // strlen, strncmp
#include <stdexcept>
#include <memory>
#include <string>
#include <type_traits>

/***********************************************************************************
 *                                  Build Options
 ************************************************************************************/
// Note that under VS 2022/2026, both __clang__ and _MSC_VER are pre-defined when using the Clang toolset
#if defined(__GNUC__) || defined(__clang__)
#define FP128_CLANG
#elif defined(_MSC_VER)
#define FP128_MSVC
#endif

// Detect AArch64 (ARMv8 64 bit, e.g. Apple Silicon). Used to select hand written
// assembly implementations of the x86 intrinsics emulated in the Clang/GCC section below.
#if defined(FP128_CLANG) && (defined(__aarch64__) || defined(_M_ARM64))
#define FP128_ARM64
#endif

// Set to TRUE to disable function inlining - useful for profiling a specific function. Default 0
#ifndef FP128_DISABLE_INLINE
#define FP128_DISABLE_INLINE 0
#endif

#ifndef FP128_NO_INLINE
#if defined(FP128_MSVC)
#define FP128_NO_INLINE __declspec(noinline)
#else
#define FP128_NO_INLINE __attribute__((noinline))
#endif
#endif  // FP128_NO_INLINE

/**
 * @def FP128_FORCE_INLINE
 * @brief Marks a function the compiler must inline rather than merely may.
 *
 * Which of the two markers a member of the 128 bit types carries follows a single rule:
 * <UL>
 * <LI>FP128_FORCE_INLINE for the trivial members - accessors, comparisons, conversions between the
 *     builtin types - and for every member that only converts its operand and delegates to another
 *     member, such as the binary operators that forward to their compound assignment counterparts
 *     and the mixed type overloads that wrap their argument in the class type.</LI>
 * <LI>FP128_INLINE for everything with a body of its own: division and modulo, the string
 *     conversions and the transcendental functions.</LI>
 * </UL>
 *
 * The distinction matters because all four types are 16 bytes, which the x64 ABI passes and returns
 * in memory. Left to its own judgement MSVC declines to inline these forwarders under /GL (whole
 * program optimization), and the round trip through memory then costs more than the operation being
 * forwarded to. Forcing the wrapper open does not force the callee open with it, so the code that
 * does the actual work still gets outlined when the optimizer thinks that is better.
 *
 * The by value shift operators of float128 and fixed_point128 are the exception, and say why at
 * their definitions: what they forward to is large enough that expanding the wrapper crowds out the
 * arithmetic around it. Treat any addition to the forced set the same way - measure it, and record
 * the number if the answer is surprising.
 */
#if FP128_DISABLE_INLINE != 0
#define FP128_INLINE       FP128_NO_INLINE
#define FP128_FORCE_INLINE FP128_NO_INLINE
#else

#if defined(FP128_MSVC)
#define FP128_INLINE inline
#define FP128_FORCE_INLINE __forceinline
#else
#define FP128_INLINE inline
#define FP128_FORCE_INLINE inline __attribute__((always_inline))
#endif
#endif

static constexpr bool FP128_CPP_STYLE_MODULO = true;  ///< Use C++ modulo semantics (false = Python-style).

/**
 * @def FP128_USE_RECIPROCAL_FOR_DIVISION
 * @brief Selects how fixed_point128 divides by a value that is neither a power of two nor an integer.
 *
 * Non zero (the default) computes <tt>a / b</tt> as <tt>a * reciprocal(b)</tt>, where reciprocal()
 * refines a double precision estimate with Newton iterations. Zero selects the hand written long
 * division instead. Only the general case is affected either way: a power of two divisor is still
 * turned into a shift, and an integral divisor that fits in 64 bit still goes through div_64bit().
 *
 * Override it on the command line (<tt>/DFP128_USE_RECIPROCAL_FOR_DIVISION=0</tt> or
 * <tt>-DFP128_USE_RECIPROCAL_FOR_DIVISION=0</tt>) or by defining it before including any header of
 * this library, exactly as with FP128_DISABLE_INLINE.
 *
 * The default is non zero because the reciprocal is the faster of the two. Measured over a table of
 * 256 random divisors, it runs at 1.4x to 1.7x the rate of the long division for fixed_point128<10>
 * and 1.8x for fixed_point128<32>, on both MSVC and clang-cl. What it buys with that is accuracy:
 * the two algorithms disagree on roughly 40% of those divisors, by up to 1.7 ulp, and comparing the
 * residual <tt>|a - q * b|</tt> of each puts the long division closer to the exact quotient every
 * single time it differs. Divide with this off when the last two bits have to be right.
 *
 * The flag is deliberately specific to fixed_point128. The other three types are not built the same
 * way, and measuring them says to leave them alone:
 * <UL>
 * <LI>float128 loses on both counts - multiplying by a reciprocal runs at 0.47x (MSVC) to 0.55x
 *     (clang-cl) of its long division, and is the less accurate of the two by the same residual test.
 *     reciprocal() is the slower half: it normalizes its operand and then runs two or three float128
 *     multiplications, each of which renormalizes and rounds, where the fixed_point128 equivalent
 *     multiplies raw 128 bit words. So float128::operator/=() always divides.</LI>
 * <LI>uint128_t and int128_t have no reciprocal to multiply by - 1 / b is zero for every |b| > 1 in
 *     an integer representation. Their equivalent is a scaled reciprocal, floor(2^k / b) followed by
 *     a multiply and a correction, and computing that scaled reciprocal is itself a division of the
 *     kind it is meant to replace. It only pays when one divisor is reused across many dividends,
 *     which operator/=() has no way to know, so those two divide directly.</LI>
 * </UL>
 *
 * @note Setting this to zero used to be enough to fail a unit test, which is worth recording because
 *       the cause was not the division. asin() and acos() refine a double precision estimate with
 *       Newton's method, and at the ends of their domain the derivative they divide by is ~6.1e-17
 *       while the numerator is nothing but the error of sin(). The exact quotient of those two is
 *       a correction of ~7e-5 applied to an estimate that was already right, and the iteration never
 *       recovered. The reciprocal path had been hiding it rather than avoiding it: 1 / 6.1e-17
 *       overflows the type, so reciprocal() saturated and the bad correction came out clamped to
 *       something harmless. Both functions now reject a step that does not improve their residual,
 *       and produce identical results under either setting of this flag.
 */
#ifndef FP128_USE_RECIPROCAL_FOR_DIVISION
#define FP128_USE_RECIPROCAL_FOR_DIVISION 1
#endif

/***********************************************************************************
 *                                  Macros
 ***********************************************************************************/
/** @name Utility Macros
 *  @brief Bit manipulation and exception macros for fixed-point arithmetic.
 *  @{
 */

#define FP128_ONE_SHIFT(x) (1ull << (x))
#define FP128_MAX_VALUE_64(x) (UINT64_MAX >> (64 - (x)))
#define FP128_GET_BIT(x, n) (((x) >> (n)) & 1)
#define FP128_GET_BITS(x, b, count) (((x) >> (b)) & FP128_MAX_VALUE_64(count))
#define FP128_INT_DIVIDE_BY_ZERO_EXCEPTION throw std::logic_error("Integer divide by zero!")
#define FP128_FLOAT_DIVIDE_BY_ZERO_EXCEPTION throw std::logic_error("Floating point divide by zero!")
#define FP128_NOT_IMPLEMENTED_EXCEPTION throw std::runtime_error("Not implemented!")
#if defined _DEBUG || defined DEBUG
#define FP128_ASSERT assert
#define FP128_THROW_ONLY_IN_DEBUG
#else
#define FP128_ASSERT(x)
#define FP128_THROW_ONLY_IN_DEBUG noexcept
#endif  // _DEBUG

/// @name Platform-Specific Intrinsics
/// @brief Compiler-specific intrinsic wrappers for bit manipulation and extended arithmetic.
/// @{

//
// MSVC specific section - keep intrinsics as-is
//
#if defined(FP128_MSVC)
#pragma warning(disable : 4996)

#include <intrin.h>
#include <immintrin.h>

#define udiv64 _udiv64
#define udiv128 _udiv128
#define alloca _alloca

// The 128 bit funnel shifts, which map to a single SHRD/SHLD instruction here.
// See the note on the naming in the Clang section below.
#define FP128_SHIFTRIGHT128 fp128_shiftright128
#define FP128_SHIFTLEFT128 fp128_shiftleft128

//
// The bit counting and extended arithmetic intrinsics are wrapped in constexpr functions rather
// than aliased with a macro, so that the operations built on them (fixed_point128's addition,
// multiplication and get_exponent among them) can be evaluated at compile time.
//
// An intrinsic is never a constant expression, so each wrapper serves a constant evaluated call
// from a portable implementation of the same operation and leaves every other call to the
// intrinsic. std::is_constant_evaluated() is folded away by the optimizer, so the runtime path is
// the bare intrinsic it was before, and it is the only path a runtime call can take: the portable
// code never reaches code generation.
//

/**
 * @brief Right funnel shift of a 128 bit value held in two QWORDs.
 *
 * Wraps __shiftright128, a single SHRD instruction. The shift count is taken modulo 64 and a
 * count of zero returns the low QWORD unchanged, which is what the instruction does and what the
 * constant evaluated path below reproduces.
 *
 * @param l Low QWORD.
 * @param h High QWORD.
 * @param shift Bits to shift.
 * @return Lower 64 bits of the result.
 */
FP128_FORCE_INLINE constexpr uint64_t fp128_shiftright128(uint64_t l, uint64_t h, unsigned char shift) noexcept
{
    if (std::is_constant_evaluated()) {
        const unsigned char s = shift & 63;
        // shifting a 64 bit value by 64 is undefined in C++, which is what the zero case would reach
        return (s == 0) ? l : (l >> s) | (h << (64 - s));
    }
    return __shiftright128(l, h, shift);
}

/**
 * @brief Left funnel shift of a 128 bit value held in two QWORDs.
 *
 * Wraps __shiftleft128, a single SHLD instruction. The shift count is taken modulo 64 and a count
 * of zero returns the high QWORD unchanged.
 *
 * @param l Low QWORD.
 * @param h High QWORD.
 * @param shift Bits to shift.
 * @return Upper 64 bits of the result.
 */
FP128_FORCE_INLINE constexpr uint64_t fp128_shiftleft128(uint64_t l, uint64_t h, unsigned char shift) noexcept
{
    if (std::is_constant_evaluated()) {
        const unsigned char s = shift & 63;
        return (s == 0) ? h : (h << s) | (l >> (64 - s));
    }
    return __shiftleft128(l, h, shift);
}

/**
 * @brief Count leading zeros in a 32-bit value.
 * @param x Value to inspect.
 * @return Number of leading zero bits, 32 for a zero operand.
 */
FP128_FORCE_INLINE constexpr uint32_t lzcnt32(uint32_t x) noexcept
{
    if (std::is_constant_evaluated()) {
        return static_cast<uint32_t>(std::countl_zero(x));
    }
    return __lzcnt(x);
}

/**
 * @brief Count leading zeros in a 64-bit value.
 * @param x Value to inspect.
 * @return Number of leading zero bits, 64 for a zero operand.
 */
FP128_FORCE_INLINE constexpr uint64_t lzcnt64(uint64_t x) noexcept
{
    if (std::is_constant_evaluated()) {
        return static_cast<uint64_t>(std::countl_zero(x));
    }
    return __lzcnt64(x);
}

/**
 * @brief Count set bits in a 32-bit value.
 * @param x Value to inspect.
 * @return Number of 1 bits in x.
 */
FP128_FORCE_INLINE constexpr uint32_t popcnt32(uint32_t x) noexcept
{
    if (std::is_constant_evaluated()) {
        return static_cast<uint32_t>(std::popcount(x));
    }
    return __popcnt(x);
}

/**
 * @brief Count set bits in a 64-bit value.
 * @param x Value to inspect.
 * @return Number of 1 bits in x.
 */
FP128_FORCE_INLINE constexpr uint64_t popcnt64(uint64_t x) noexcept
{
    if (std::is_constant_evaluated()) {
        return static_cast<uint64_t>(std::popcount(x));
    }
    return __popcnt64(x);
}

/**
 * @brief 64x64 -> 128-bit unsigned multiply, wrapping the MULX instruction.
 *
 * MSVC has no 128 bit integer type to fall back on, so the constant evaluated path assembles the
 * product from four 32x32 -> 64 bit partial products, the schoolbook algorithm on 32 bit limbs.
 *
 * @param a First operand.
 * @param b Second operand.
 * @param hi Pointer to receive the upper 64 bits of the product.
 * @return Lower 64 bits of the product.
 */
FP128_FORCE_INLINE constexpr uint64_t mulx_u64(uint64_t a, uint64_t b, uint64_t* hi) noexcept
{
    FP128_ASSERT(hi != nullptr);  // Caller must provide a valid pointer for the high part.
    if (std::is_constant_evaluated()) {
        const uint64_t a_lo = a & UINT32_MAX, a_hi = a >> 32;
        const uint64_t b_lo = b & UINT32_MAX, b_hi = b >> 32;
        const uint64_t p_ll = a_lo * b_lo;
        const uint64_t p_lh = a_lo * b_hi;
        const uint64_t p_hl = a_hi * b_lo;
        const uint64_t p_hh = a_hi * b_hi;
        // The middle column holds the carry out of the low 32 bits plus the low halves of both
        // cross products. Three values below 2^32 cannot overflow a 64 bit accumulator.
        const uint64_t mid = (p_ll >> 32) + (p_lh & UINT32_MAX) + (p_hl & UINT32_MAX);
        *hi = p_hh + (p_lh >> 32) + (p_hl >> 32) + (mid >> 32);
        return (mid << 32) | (p_ll & UINT32_MAX);
    }
    return _mulx_u64(a, b, hi);
}

/**
 * @brief 64-bit add with carry, wrapping the ADCX instruction.
 * @param c Input carry (0 or 1).
 * @param a First operand.
 * @param b Second operand.
 * @param out Pointer to receive the 64-bit sum.
 * @return Output carry (0 or 1).
 */
FP128_FORCE_INLINE constexpr unsigned char addcarryx_u64(unsigned char c, uint64_t a, uint64_t b, uint64_t* out) noexcept
{
    FP128_ASSERT(out != nullptr);  // Caller must provide a valid pointer for the result.
    if (std::is_constant_evaluated()) {
        // Unsigned addition wraps, so a sum that came out smaller than an operand is exactly the
        // carry out. The total is at most 2^65-1, so at most one of the two additions can carry.
        const uint64_t sum = a + b;
        const uint64_t res = sum + c;
        *out = res;
        return static_cast<unsigned char>((sum < a) | (res < sum));
    }
    return _addcarryx_u64(c, a, b, out);
}

/**
 * @brief 64-bit subtract with borrow, wrapping the SBB instruction.
 * @param b Input borrow (0 or 1).
 * @param a Minuend.
 * @param c Subtrahend.
 * @param out Pointer to receive the 64-bit difference.
 * @return Output borrow (0 or 1).
 */
FP128_FORCE_INLINE constexpr unsigned char subborrow_u64(unsigned char b, uint64_t a, uint64_t c, uint64_t* out) noexcept
{
    FP128_ASSERT(out != nullptr);  // Caller must provide a valid pointer for the result.
    if (std::is_constant_evaluated()) {
        // Unsigned subtraction wraps, so a difference that came out larger than the minuend is
        // exactly the borrow out. The result is at least -(2^65-1), so at most one of the two
        // subtractions can borrow.
        const uint64_t diff = a - c;
        const uint64_t res = diff - b;
        *out = res;
        return static_cast<unsigned char>((diff > a) | (res > diff));
    }
    return _subborrow_u64(b, a, c, out);
}

//
// GCC/Clang portable fallback implementations
//
#elif defined(FP128_CLANG)
// The funnel shifts fall back to this library's own implementations.
//
// These deliberately do NOT reuse the __shiftright128 / __shiftleft128 spelling. Those names
// belong to the implementation, and Microsoft's STL calls them for real inside __msvc_int128.hpp,
// which <algorithm> and friends pull in. Defining them as macros here rewrites those calls to an
// unqualified name that is not visible at that point, so any translation unit that reached an
// fp128 header before <algorithm> failed to compile. Prefixing the macro keeps the two apart.
#define FP128_SHIFTRIGHT128 shift_right128
#define FP128_SHIFTLEFT128 shift_left128

#if defined(FP128_ARM64)
#include <arm_neon.h>  // for the uint8x8_t type used as an operand of the NEON population count assembly
#endif

/**
 * @brief Portable 32-bit unsigned division with remainder (GCC/Clang fallback).
 *
 * On AArch64 this is a single UDIV followed by an MSUB (multiply-subtract) which computes
 * the remainder without a second division. Both operands are widened to 64 bit, matching the
 * x86 _udiv64 semantics of a 64 bit numerator; unlike x86 there is no trap when the quotient
 * exceeds 32 bits, the value is simply truncated.
 *
 * @param dividend Numerator.
 * @param divisor Denominator.
 * @param remainder Pointer to receive the remainder (may be nullptr).
 * @return Quotient.
 */
FP128_FORCE_INLINE constexpr uint32_t udiv64(uint64_t dividend, uint32_t divisor, uint32_t* remainder)
{
#if defined(FP128_ARM64)
    // Inline assembly is not allowed during constant evaluation, use the portable path instead.
    if (!std::is_constant_evaluated()) {
        const uint64_t divisor64 = divisor;
        uint64_t quot64 = 0, rem64 = 0;
        // quot64 = dividend / divisor64
        // rem64  = dividend - quot64 * divisor64
        __asm__("udiv %[q], %[n], %[d]\n\t"
                "msub %[r], %[q], %[d], %[n]"
                : [q] "=&r"(quot64), [r] "=r"(rem64)  // q is early-clobber: MSUB reads n and d after q was written
                : [n] "r"(dividend), [d] "r"(divisor64));
        if (remainder) {
            *remainder = static_cast<uint32_t>(rem64);
        }
        return static_cast<uint32_t>(quot64);
    }
#endif
    uint32_t quot = dividend / divisor;
    if (remainder) {
        *remainder = dividend % divisor;
    }
    return quot;
}

/**
 * @brief Portable 128-bit by 64-bit unsigned division (GCC/Clang fallback).
 *
 * @note AArch64 has no 128/64 bit divide instruction (UDIV is at most 64/64), so there is no
 *       assembly variant of this function. The __uint128_t expression below lowers to a call to
 *       the compiler runtime helper (__udivti3 / __umodti3) which is faster than any short
 *       hand written long division. Callers on ARM64 should prefer the reciprocal based
 *       division path (see FP128_USE_RECIPROCAL_FOR_DIVISION) built on mulx_u64.
 *
 * @param hi_dividend Upper 64 bits of the 128-bit dividend.
 * @param lo_dividend Lower 64 bits of the 128-bit dividend.
 * @param divisor 64-bit divisor.
 * @param remainder Pointer to receive the 64-bit remainder (may be nullptr).
 * @return 64-bit quotient.
 */
FP128_FORCE_INLINE constexpr uint64_t udiv128(uint64_t hi_dividend, uint64_t lo_dividend, uint64_t divisor, uint64_t* remainder)
{
    __uint128_t dividend = (static_cast<__uint128_t>(hi_dividend) << 64) | lo_dividend;
    uint64_t quot = static_cast<uint64_t>(dividend / divisor);
    if (remainder) {
        *remainder = static_cast<uint64_t>(dividend % divisor);
    }
    return quot;
}

/**
 * @brief 64x64 -> 128-bit unsigned multiply (GCC/Clang fallback for _mulx_u64).
 *
 * The AArch64 variant uses the MUL/UMULH instruction pair, the direct equivalent of the x86
 * MULX instruction. The two instructions are independent and dual issue on Apple Silicon.
 *
 * @param a First operand.
 * @param b Second operand.
 * @param hi Pointer to receive the upper 64 bits of the product.
 * @return Lower 64 bits of the product.
 */
FP128_FORCE_INLINE static constexpr uint64_t mulx_u64(uint64_t a, uint64_t b, uint64_t* hi) noexcept
{
    FP128_ASSERT(hi != nullptr);  // Caller must provide a valid pointer for the high part. Compatibility with MSVC intrinsic.
#if defined(FP128_ARM64)
    // Inline assembly is not allowed during constant evaluation, use the portable path instead.
    if (!std::is_constant_evaluated()) {
        uint64_t lo_res = 0, hi_res = 0;
        __asm__("umulh %[hi], %[a], %[b]\n\t"
                "mul   %[lo], %[a], %[b]"
                : [hi] "=&r"(hi_res), [lo] "=r"(lo_res)  // hi is early-clobber: MUL reads a and b after hi was written
                : [a] "r"(a), [b] "r"(b));
        *hi = hi_res;
        return lo_res;
    }
#endif
    __uint128_t r = (__uint128_t)a * b;
    *hi = (uint64_t)(r >> 64);
    return (uint64_t)r;
}

/**
 * @brief 64-bit add with carry (GCC/Clang fallback for _addcarryx_u64).
 *
 * The AArch64 variant mirrors the x86 ADCX sequence: CMP materializes the incoming carry into
 * the C flag (the unsigned compare c - 1 sets C when c is non zero), ADCS performs the addition
 * with that carry and updates the flags, and CSET extracts the outgoing carry.
 *
 * @param c Input carry (0 or 1).
 * @param a First operand.
 * @param b Second operand.
 * @param out Pointer to receive the 64-bit sum.
 * @return Output carry (0 or 1).
 */
FP128_FORCE_INLINE static constexpr unsigned char addcarryx_u64(unsigned char c, uint64_t a, uint64_t b, uint64_t* out) noexcept
{
    FP128_ASSERT(out != nullptr);  // Caller must provide a valid pointer for the result. Compatibility with MSVC intrinsic.
#if defined(FP128_ARM64)
    // Inline assembly is not allowed during constant evaluation, use the portable path instead.
    if (!std::is_constant_evaluated()) {
        uint64_t sum = 0;
        uint32_t carry_out = 0;
        __asm__("cmp  %w[cin], #1\n\t"      // C = (cin != 0)
                "adcs %[sum], %[a], %[b]\n\t"
                "cset %w[cout], cs"
                : [sum] "=&r"(sum), [cout] "=&r"(carry_out)  // both outputs are early-clobber, they must not alias the inputs
                : [cin] "r"((uint32_t)c), [a] "r"(a), [b] "r"(b)
                : "cc");
        *out = sum;
        return (unsigned char)carry_out;
    }
#endif
    __uint128_t r = (__uint128_t)a + b + c;
    *out = (uint64_t)r;
    return (unsigned char)(r >> 64);
}

/**
 * @brief 64-bit subtract with borrow (GCC/Clang fallback for _subborrow_u64).
 *
 * Unlike its addcarryx_u64 neighbour this has no hand written AArch64 variant. The subtraction of a
 * borrow is expressed naturally by the __uint128_t expression below, which Clang already lowers to
 * the SUBS/SBCS pair an assembly version would have spelled out, so there is nothing left to win.
 *
 * @param b Input borrow (0 or 1).
 * @param a Minuend.
 * @param c Subtrahend.
 * @param out Pointer to receive the 64-bit difference.
 * @return Output borrow (0 or 1).
 */
FP128_FORCE_INLINE static constexpr unsigned char subborrow_u64(unsigned char b, uint64_t a, uint64_t c, uint64_t* out) noexcept
{
    FP128_ASSERT(out != nullptr);  // Caller must provide a valid pointer for the result. Compatibility with MSVC intrinsic.
    __uint128_t r = (__uint128_t)a - c - b;
    *out = (uint64_t)r;
    // The difference wraps when it borrows, so bit 64 of the 128-bit result is the borrow out.
    return (unsigned char)((r >> 64) & 1);
}

/**
 * @brief Count leading zeros in a 32-bit value (GCC/Clang fallback).
 *
 * The AArch64 CLZ instruction is defined for a zero operand (it returns the operand width),
 * so unlike the x86 BSR based lowering no zero test is needed.
 */
FP128_FORCE_INLINE static constexpr uint32_t lzcnt32(uint32_t x) noexcept
{
#if defined(FP128_ARM64)
    // Inline assembly is not allowed during constant evaluation, use the portable path instead.
    if (!std::is_constant_evaluated()) {
        uint32_t res = 0;
        __asm__("clz %w[res], %w[val]" : [res] "=r"(res) : [val] "r"(x));
        return res;
    }
#endif
    return (x == 0) ? 32u : (uint32_t)__builtin_clz(x);
}

/**
 * @brief Count leading zeros in a 64-bit value (GCC/Clang fallback).
 *
 * The AArch64 CLZ instruction returns 64 for a zero operand, no zero test is needed.
 */
FP128_FORCE_INLINE static constexpr uint64_t lzcnt64(uint64_t x) noexcept
{
#if defined(FP128_ARM64)
    // Inline assembly is not allowed during constant evaluation, use the portable path instead.
    if (!std::is_constant_evaluated()) {
        uint64_t res = 0;
        __asm__("clz %[res], %[val]" : [res] "=r"(res) : [val] "r"(x));
        return res;
    }
#endif
    return (x == 0) ? 64u : (uint64_t)__builtin_clzll(x);
}

/**
 * @brief Count set bits in a 64-bit value (GCC/Clang fallback).
 *
 * AArch64 has no general purpose register population count before ARMv8.9 (FEAT_CSSC), which
 * Apple Silicon does not implement. The value is moved to a NEON register, CNT counts the bits
 * of each of the 8 bytes in parallel and ADDV sums the 8 byte lanes into a single byte.
 */
FP128_FORCE_INLINE static constexpr uint64_t popcnt64(uint64_t x) noexcept
{
#if defined(FP128_ARM64)
    // Inline assembly is not allowed during constant evaluation, use the portable path instead.
    if (!std::is_constant_evaluated()) {
        uint32_t res = 0;
        uint8x8_t tmp;
        __asm__("fmov %d[tmp], %[val]\n\t"
                "cnt  %[tmp].8b, %[tmp].8b\n\t"
                "addv %b[tmp], %[tmp].8b\n\t"
                "fmov %w[res], %s[tmp]"
                : [res] "=r"(res), [tmp] "=&w"(tmp)
                : [val] "r"(x));
        return res;
    }
#endif
    return (uint64_t)__builtin_popcountll(x);
}

/**
 * @brief Count set bits in a 32-bit value (GCC/Clang fallback).
 *
 * Same NEON sequence as popcnt64. FMOV of a W register zeroes the upper lanes of the NEON
 * register, so the 4 unused bytes contribute nothing to the ADDV sum.
 */
FP128_FORCE_INLINE static constexpr uint32_t popcnt32(uint32_t x) noexcept
{
#if defined(FP128_ARM64)
    // Inline assembly is not allowed during constant evaluation, use the portable path instead.
    if (!std::is_constant_evaluated()) {
        uint32_t res = 0;
        uint8x8_t tmp;
        __asm__("fmov %s[tmp], %w[val]\n\t"
                "cnt  %[tmp].8b, %[tmp].8b\n\t"
                "addv %b[tmp], %[tmp].8b\n\t"
                "fmov %w[res], %s[tmp]"
                : [res] "=r"(res), [tmp] "=&w"(tmp)
                : [val] "r"(x));
        return res;
    }
#endif
    return (uint32_t)__builtin_popcount(x);
}

#endif  // #if defined (FP128_CLANG)

// portable alignment macro
#if defined(FP128_MSVC)
#define FP128_ALIGN(_a) __declspec(align(_a))
//
// Clang/GCC specific macros
//
#elif defined(FP128_CLANG)
#define FP128_ALIGN(_a) __attribute__((aligned(_a)))
#endif

#define FP128_ALIGN16 FP128_ALIGN(16)
/** @} */  // end Utility Macros

namespace fp128
{

/***********************************************************************************
 *                                  Constants
 ************************************************************************************/

/// @name IEEE 754 Layout Constants
/// @{
static constexpr int32_t flt_frac_bits = 23;  ///< Mantissa bit count of an IEEE 754 float.
static constexpr int32_t flt_exp_bits = 8;    ///< Exponent bit count of an IEEE 754 float.
static constexpr int32_t dbl_frac_bits = 52;  ///< Mantissa bit count of an IEEE 754 double.
static constexpr int32_t dbl_exp_bits = 11;   ///< Exponent bit count of an IEEE 754 double.
/// @}

/***********************************************************************************
 *                                  Containers
 ************************************************************************************/
/**
 * @struct Double
 * @brief Bit level view of an IEEE 754 double-precision value.
 *
 * Gives access to the mantissa, exponent and sign of a double without manual bit shifting, and
 * assembles a double back out of those three fields.
 *
 * The value is held as a raw bit pattern and converted with std::bit_cast rather than being
 * overlaid with a union. The two produce identical code - a single move between register classes -
 * but reading the inactive member of a union is not allowed during constant evaluation, which
 * would keep every conversion between these types and a double out of a constant expression.
 */
struct Double {
    uint64_t bits;  ///< The raw bit pattern of the double.

    /// @brief Constructs from a double, or from positive zero when no value is given.
    constexpr Double(double v = 0) noexcept : bits(std::bit_cast<uint64_t>(v)) {}

    /// @brief Mantissa (fraction) bits.
    [[nodiscard]] constexpr uint64_t f() const noexcept { return bits & FRAC_MASK; }
    /// @brief Biased exponent bits.
    [[nodiscard]] constexpr uint64_t e() const noexcept { return (bits >> dbl_frac_bits) & EXP_MASK; }
    /// @brief Sign bit (0 = positive, 1 = negative).
    [[nodiscard]] constexpr uint64_t s() const noexcept { return bits >> SIGN_SHIFT; }
    /// @brief The double these bits encode.
    [[nodiscard]] constexpr double val() const noexcept { return std::bit_cast<double>(bits); }

    /// @brief Sets the mantissa. Bits above the field width are dropped, as a bit-field assignment would.
    constexpr void set_f(uint64_t v) noexcept { bits = (bits & ~FRAC_MASK) | (v & FRAC_MASK); }
    /// @brief Sets the biased exponent. Bits above the field width are dropped.
    constexpr void set_e(uint64_t v) noexcept { bits = (bits & ~(EXP_MASK << dbl_frac_bits)) | ((v & EXP_MASK) << dbl_frac_bits); }
    /// @brief Sets the sign bit.
    constexpr void set_s(uint64_t v) noexcept { bits = (bits & ~(1ull << SIGN_SHIFT)) | ((v & 1) << SIGN_SHIFT); }

    /**
     * @brief Assembles a double out of its three fields.
     *
     * Preferred over three set_ calls wherever the value is built from scratch: those each mask
     * the old field out before merging the new one, which is wasted work when there is nothing
     * there yet.
     *
     * @param s Sign bit.
     * @param e Biased exponent. Bits above the field width are dropped.
     * @param f Mantissa. Bits above the field width are dropped.
     * @return The assembled double.
     */
    [[nodiscard]] static constexpr double make(uint64_t s, uint64_t e, uint64_t f) noexcept
    {
        return std::bit_cast<double>(((s & 1) << SIGN_SHIFT) | ((e & EXP_MASK) << dbl_frac_bits) | (f & FRAC_MASK));
    }

private:
    static constexpr uint64_t FRAC_MASK = (1ull << dbl_frac_bits) - 1;   ///< Mask of the mantissa field.
    static constexpr uint64_t EXP_MASK = (1ull << dbl_exp_bits) - 1;     ///< Mask of the exponent field, once shifted down.
    static constexpr int32_t SIGN_SHIFT = dbl_frac_bits + dbl_exp_bits;  ///< Bit position of the sign.
};
static_assert(sizeof(Double) == sizeof(double), "The Double view should have the same size as a double variable!");

/**
 * @struct Float
 * @brief Bit level view of an IEEE 754 single-precision value.
 *
 * The single precision counterpart of @ref Double, see its documentation.
 */
struct Float {
    uint32_t bits;  ///< The raw bit pattern of the float.

    /// @brief Constructs from a float, or from positive zero when no value is given.
    constexpr Float(float v = 0) noexcept : bits(std::bit_cast<uint32_t>(v)) {}

    /// @brief Mantissa (fraction) bits.
    [[nodiscard]] constexpr uint32_t f() const noexcept { return bits & FRAC_MASK; }
    /// @brief Biased exponent bits.
    [[nodiscard]] constexpr uint32_t e() const noexcept { return (bits >> flt_frac_bits) & EXP_MASK; }
    /// @brief Sign bit (0 = positive, 1 = negative).
    [[nodiscard]] constexpr uint32_t s() const noexcept { return bits >> SIGN_SHIFT; }
    /// @brief The float these bits encode.
    [[nodiscard]] constexpr float val() const noexcept { return std::bit_cast<float>(bits); }

    /// @brief Sets the mantissa. Bits above the field width are dropped, as a bit-field assignment would.
    constexpr void set_f(uint32_t v) noexcept { bits = (bits & ~FRAC_MASK) | (v & FRAC_MASK); }
    /// @brief Sets the biased exponent. Bits above the field width are dropped.
    constexpr void set_e(uint32_t v) noexcept { bits = (bits & ~(EXP_MASK << flt_frac_bits)) | ((v & EXP_MASK) << flt_frac_bits); }
    /// @brief Sets the sign bit.
    constexpr void set_s(uint32_t v) noexcept { bits = (bits & ~(1u << SIGN_SHIFT)) | ((v & 1) << SIGN_SHIFT); }

    /// @brief Assembles a float out of its three fields, see Double::make().
    /// @param s Sign bit. @param e Biased exponent. @param f Mantissa. @return The assembled float.
    [[nodiscard]] static constexpr float make(uint32_t s, uint32_t e, uint32_t f) noexcept
    {
        return std::bit_cast<float>(((s & 1) << SIGN_SHIFT) | ((e & EXP_MASK) << flt_frac_bits) | (f & FRAC_MASK));
    }

private:
    static constexpr uint32_t FRAC_MASK = (1u << flt_frac_bits) - 1;     ///< Mask of the mantissa field.
    static constexpr uint32_t EXP_MASK = (1u << flt_exp_bits) - 1;       ///< Mask of the exponent field, once shifted down.
    static constexpr int32_t SIGN_SHIFT = flt_frac_bits + flt_exp_bits;  ///< Bit position of the sign.
};
static_assert(sizeof(Float) == sizeof(float), "The Float view should have the same size as a float variable!");

/***********************************************************************************
 *                                  Functions
 ************************************************************************************/

/**
 * @brief Copy a string to lowercase, up to n characters.
 * @param dest Destination buffer.
 * @param src Source string.
 * @param n Maximum number of characters to copy.
 * @return 0 on success, EINVAL if dest or src is nullptr.
 */
FP128_INLINE static errno_t strnlwr(char* dest, const char* src, size_t n) noexcept
{
    if (src == nullptr || dest == nullptr)
        return EINVAL;
    size_t i = 0;
    for (; i < n && src[i] != '\0'; ++i) {
        dest[i] = static_cast<char>(::tolower(static_cast<unsigned char>(src[i])));
    }
    if (i < n)
        dest[i] = '\0';
    return 0;
}

/**
 * @brief Calculates the element count of a C style array at build time
 * Example:
 * int a[5];
 * constexpr int len = array_length(a); // returns 5 at build time
 * @tparam T array type
 * @param a array instance
 * @return Element count in array
 */
template <typename T>[[nodiscard]] constexpr uint32_t array_length(const T& a)
{
    static_assert(sizeof(a[0]) != 0, "Requires an array of non-zero sized elements!");
    return sizeof(a) / sizeof(a[0]);
}
/**
 * @brief Right-shift a 64-bit value with rounding.
 *
 * Undefined behavior when shift is outside the range [1, 63].
 *
 * @param x Value to shift.
 * @param shift Number of bits to shift.
 * @return The rounded result of x >> shift.
 */
[[nodiscard]] FP128_INLINE constexpr uint64_t shift_right64_round(uint64_t x, int shift) noexcept
{
    FP128_ASSERT(shift > 0 && shift < 64);
    x += 1ull << (shift - 1);
    return x >> shift;
}

/**
 * @brief Right shift a 128 bit unsigned integer (inplace).
 * Limited range, inplace and no parameter checks (checked via assert in debug builds)
 * @param l Low QWORD
 * @param h High QWORD
 * @param shift Bits to shift, between 1-63
 * @return void
 */
FP128_INLINE constexpr void shift_right128_inplace(uint64_t& l, uint64_t& h, int shift) noexcept
{
    FP128_ASSERT(shift > 0 && shift < 64);
    l = (l >> shift) | (h << (64 - shift));
    h >>= shift;
}
/**
 * @brief Left shift a 128 bit integer (inplace).
 * Limited range, inplace and no parameter checks (checked via assert in debug builds)
 * @param l Low QWORD
 * @param h High QWORD
 * @param shift Bits to shift, between 1-63
 * @return void
 */
FP128_INLINE constexpr void shift_left128_inplace(uint64_t& l, uint64_t& h, int shift) noexcept
{
    FP128_ASSERT(shift > 0 && shift < 64);
    h = (h << shift) | (l >> (64 - shift));
    l <<= shift;
}
/**
 * @brief Right shift a 128 bit integer (inplace) with rounding.
 * Handles any positive shift value.
 * @param l Low QWORD
 * @param h High QWORD
 * @param shift Bits to shift, between 1-inf
 * @return void
 */
FP128_FORCE_INLINE constexpr void shift_right128_inplace_safe(uint64_t& l, uint64_t& h, int shift) noexcept
{
    FP128_ASSERT(shift >= 0);
    if (shift == 0)
        return;
    uint64_t lsb = 0;
    switch (shift >> 6) {
    case 0:  // 1-63 bit
        lsb = (shift == 1) ? (l & 3) << 1 : (l >> (shift - 2)) & 7;
        l = (l >> shift) | (h << (64 - shift));
        h >>= shift;
        break;
    case 1:  // 64-127 bit
        shift -= 64;
        switch (shift) {
        case 0:
            lsb = ((h & 1) << 2) | ((l >> 63) << 1) |
                  ((l & 0x7FFFFFFFFFFFFFFFull) != 0 ? 1 : 0);  // the last clause checks if any of the bits that got shifted away are 1, if so we need to round
                                                               // up in case of a tie (when h's lsb is 1 and the rest of the bits are zero)
            break;
        case 1:
            lsb = ((h & 3) << 1) | (l != 0 ? 1 : 0);
            break;
        default:
            lsb = (h >> (shift - 2)) & 7;
        }

        l = h >> shift;
        h = 0;
        break;
    default:  // >127 bit or negative
        h = l = 0;
    }

    // Use rounding half to even
    // Middle bit is the bit that got shifted away.
    // It get rounded up in 2 cases:
    //   1) The 2 rightmost bits are b11 (lsb == 3 or 7), this equal to 0.75
    //   2) The value's msb is 1 (odd number) and the right bits are b10 (0.5) so the result will be an even number
    if (lsb >= 6 || lsb == 3) {
        ++l;  // low will wrap around to zero if overflowed
        h += l == 0;
    }
}
/**
 * @brief Left shift a 128 bit integer (inplace).
 * Handles any positive shift value.
 * @param l Low QWORD
 * @param h High QWORD
 * @param shift Bits to shift, between 1-inf
 * @return void
 */
FP128_INLINE constexpr void shift_left128_inplace_safe(uint64_t& l, uint64_t& h, int shift) noexcept
{
    FP128_ASSERT(shift >= 0);
    if (shift == 0)
        return;

    switch (shift >> 6) {
    case 0:  // 1-63 bit
        h = (h << shift) | (l >> (64 - shift));
        l <<= shift;
        break;
    case 1:  // 64-127 bit
        h = l << (shift - 64);
        l = 0;
        break;
    default:  // >127 bit or negative
        h = l = 0;
    }
}
/**
 * @brief Right shift a 128 bit integer. When shift is a compile time constant, this function generates optimal code for all shift values.
 * @param l Low QWORD
 * @param h High QWORD
 * @
 * @return Lower 64 bit of the result
 */
template <int shift> [[nodiscard]] FP128_INLINE constexpr uint64_t shift_right128(uint64_t l, uint64_t h) noexcept
{
    FP128_ASSERT(shift >= 0 && shift < 128);
    if constexpr (shift == 0) {
        return l;
    } else if constexpr (shift < 64) {
        return (l >> shift) | (h << (64 - shift));
    } else if constexpr (shift < 128) {
        return h >> (shift - 64);
    } else {
        return 0;
    }
 }
/**
 * @brief Left shift a 128 bit integer. When shift is a compile time constant, this function generates optimal code for all shift values.
 * @param l Low QWORD
 * @param h High QWORD
 * @
 * @return Upper 64 bit of the result
 */
 template <int shift> [[nodiscard]] FP128_FORCE_INLINE constexpr uint64_t shift_left128(uint64_t l, uint64_t h) noexcept
 {
    FP128_ASSERT(shift >= 0 && shift < 128);
    if constexpr (shift == 0) {
        return h;
    } else if constexpr (shift < 64) {
        return (h << shift) | (l >> (64 - shift));
    } else if constexpr (shift < 128) {
        return l << (shift - 64);
    } else {
        return 0;
    }
}
/**
 * @brief Right shift a 128 bit integer.
 * @param l Low QWORD
 * @param h High QWORD
 * @param shift Bits to shift, between 0-127
 * @return Lower 64 bit of the result
 */
[[nodiscard]] FP128_FORCE_INLINE constexpr uint64_t shift_right128(uint64_t l, uint64_t h, int shift) noexcept
{
    FP128_ASSERT(shift >= 0 && shift < 128);
    switch (shift >> 6) {
    case 0:  // 0-63 bit
        return (l >> shift) | (h << (64 - shift));
    case 1:  // 64-127 bit
        return h >> (shift ^ 64);
    default:
        return 0;
    }
}
/**
 * @brief Right shift a 128 bit integer with rounding.
 * @param l Low QWORD
 * @param h High QWORD
 * @param shift Bits to shift, between 0-127
 * @return Lower 64 bit of the result
 */
[[nodiscard]] FP128_FORCE_INLINE constexpr uint64_t shift_right128_round(uint64_t l, uint64_t h, int shift) noexcept
{
    shift_right128_inplace_safe(l, h, shift);
    return l;
}
/**
 * @brief Left shift a 128 bit integer.
 * @param l Low QWORD
 * @param h High QWORD
 * @param shift Bits to shift, between 0-127
 * @return Upper 64 bit of the result
 */
[[nodiscard]] FP128_FORCE_INLINE constexpr uint64_t shift_left128(uint64_t l, uint64_t h, int shift) noexcept
{
    FP128_ASSERT(shift >= 0 && shift < 128);
    switch (shift >> 6) {
        case 0:  // 1-63 bit
            return (h << shift) | (l >> (64 - shift));
        case 1:
            return l << (shift - 64);
        default:
            return 0;
    }
}
/**
 * @brief converts a 128 integer to negative via 2's complement.
 * @param l Low QWORD (ref)
 * @param h High QWORD (ref)
 * @return void
 */
FP128_INLINE constexpr void twos_complement128(uint64_t& l, uint64_t& h) noexcept
{
    l = ~l + 1ull;
    h = ~h + (l == 0);
}
/**
 * @brief 32 bit words unsigned divide function. Variation of the code from the book Hacker's Delight.
 * @param q (output) Pointer to receive the quotient
 * @param r (output, optional) Pointer to receive the remainder. Can be nullptr
 * @param u Pointer Numerator, an array of uint32_t
 * @param v denominator (uint32_t)
 * @param m Count of elements in u
 * @return 0 for success
 */
FP128_INLINE static int32_t div_32bit(uint32_t* q, uint32_t* r, const uint32_t* u, uint32_t v, int64_t m) noexcept
{
    if (u == nullptr || q == nullptr || v == 0)
        return 1;

    while (m > 0 && u[m - 1] == 0)
        --m;

    uint32_t k = 0;
    for (auto j = m - 1; j >= 0; --j) {
        q[j] = udiv64((((uint64_t)k) << 32) + u[j], v, &k);
    }

    if (r != nullptr)
        *r = k;
    return 0;
}
/**
 * @brief 32 bit words unsigned divide function. Variation of the code from the book Hacker's Delight.
 * @param q (output) Pointer to receive the quotient
 * @param r (output, optional) Pointer to receive the remainder. Can be nullptr
 * @param u Pointer numerator, an array of uint32_t
 * @param v Pointer denominator, an array of uint32_t
 * @param m Count of elements in u
 * @param n Count of elements in v
 * @return 0 for success
 */
inline static int div_32bit(uint32_t* q, uint32_t* r, const uint32_t* u, const uint32_t* v, int m, int n) noexcept
{
    if (q == nullptr || u == nullptr || v == nullptr)
        return 1;

    constexpr uint64_t WORD_WIDTH = 32ull;         // bit width of a word
    constexpr uint64_t BASE = 1ull << WORD_WIDTH;  // Number base (32 bits).
    constexpr uint64_t MASK = BASE - 1;            // 32 bit mask
    uint32_t *un, *vn;                             // Normalized form of u, v.
    uint64_t qhat;                                 // Estimated quotient digit.
    uint64_t rhat;                                 // A remainder.
    uint64_t p;                                    // Product of two digits.
    int64_t t, k;                                  // Temporary variables
    int32_t i, j;                                  // Indexes
    // disable various warnings, some are bogus in VS2022.
    // the below code relies on the implied truncation (to 32 bit) of several expressions.
#if defined(FP128_MSVC)
#pragma warning(push)
#pragma warning(disable : 6255)
#pragma warning(disable : 4244)
#pragma warning(disable : 6297)
#pragma warning(disable : 6385)
#pragma warning(disable : 6386)
#pragma warning(disable : 26451)
#pragma warning(disable : 26493)
#pragma warning(disable : 26438)
#endif

    // shrink the arrays to avoid extra work on small numbers
    while (m > 0 && u[m - 1] == 0)
        --m;
    while (n > 0 && v[n - 1] == 0)
        --n;

    if (m < n || n <= 0 || v[n - 1] == 0)
        return 1;  // Return if invalid param.

    // Take care of the case of a single-digit divisor here.
    if (n == 1)
        return div_32bit(q, r, u, v[0], m);

    /* Normalize by shifting v left just enough so that its high-order
    bit is on, and shift u left the same amount. We may have to append a
    high-order digit on the dividend; we do that unconditionally. */

    const int32_t s = lzcnt32(v[n - 1]);  // 0 <= s <= WORD_WIDTH-1.
    const int32_t s_comp = WORD_WIDTH - s;
    vn = (uint32_t*)alloca(sizeof(uint32_t) * n);
    for (i = n - 1; i > 0; --i) {
        vn[i] = (v[i] << s) | ((uint64_t)v[i - 1] >> s_comp);
    }
    vn[0] = v[0] << s;

    un = (uint32_t*)alloca(sizeof(uint32_t) * (m + 1));
    un[m] = (uint64_t)u[m - 1] >> s_comp;
    for (i = m - 1; i > 0; --i)
        un[i] = (u[i] << s) | ((uint64_t)u[i - 1] >> s_comp);
    un[0] = u[0] << s;

    for (j = m - n; j >= 0; --j) {  // Main loop.
        // Compute estimate qhat of q[j].
        qhat = udiv128(0, ((uint64_t)un[j + n] << WORD_WIDTH) | un[j + n - 1], vn[n - 1], &rhat);
        // qhat = (un[j + n] * BASE + un[j + n - 1]) / vn[n - 1];
        // rhat = (un[j + n] * BASE + un[j + n - 1]) - qhat * vn[n - 1];
again:
        if (qhat >= BASE || qhat * vn[n - 2] > ((rhat << WORD_WIDTH) | un[j + n - 2])) {
            --qhat;
            rhat += vn[n - 1];
            if (rhat < BASE)
                goto again;
        }

        // Multiply and subtract.
        k = 0;
        for (i = 0; i < n; ++i) {
            p = qhat * vn[i];
            t = un[i + j] - k - (p & MASK);
            un[i + j] = t;
            k = (p >> WORD_WIDTH) - (t >> WORD_WIDTH);
        }
        t = un[j + n] - k;
        un[j + n] = t;

        q[j] = qhat;          // Store quotient digit.
        if (t < 0) {          // If we subtracted too
            q[j] = q[j] - 1;  // much, add back.
            k = 0;
            for (i = 0; i < n; ++i) {
                t = (uint64_t)un[i + j] + vn[i] + k;
                un[i + j] = t;
                k = t >> WORD_WIDTH;
            }
            un[j + n] = un[j + n] + k;
        }
    }  // End j.
    // If the caller wants the remainder, unnormalize
    // it and pass it back.
    if (r != nullptr) {
        for (i = 0; i < n - 1; ++i)
            r[i] = (un[i] >> s) | ((uint64_t)un[i + 1] << s_comp);

        r[n - 1] = un[n - 1] >> s;
    }
    return 0;
#if defined(FP128_MSVC)
#pragma warning(pop)
#endif
}
/**
 * @brief 64 bit words unsigned divide function. Variation of the code from the book Hacker's Delight.
 * @param q (output) Pointer to receive the quotient. Expected to be initialized to zero
 * @param r (output, optional) Pointer to receive the remainder. Can be nullptr
 * @param u Pointer to Numerator, an array of uint64_t
 * @param v denominator (uint64_t)
 * @param m Count of elements in u
 * @return 0 for success
 */
FP128_INLINE static int32_t div_64bit(uint64_t* q, uint64_t* r, const uint64_t* u, uint64_t v, int64_t m) noexcept
{
    if (u == nullptr || q == nullptr)
        return 1;
    uint64_t dummy_reminder {};
    if (r == nullptr) {
        r = &dummy_reminder;
    }

    if (v == 0)  // error case
        return 1;

    while (m > 0 && u[m - 1] == 0)
        --m;

    // Trivial cases
    if (m < 2) {
        if (m == 0 || u[0] < v) {
            *r = (m == 0) ? 0 : u[0];
            return 0;
        }
        if (u[0] == v) {
            *q = 1;
            *r = 0;
            return 0;
        }
    }

    uint64_t k[2] = {};
    for (auto j = m - 1; j >= 0; --j) {
        k[0] = u[j];
        q[j] = udiv128(k[1], k[0], v, &k[1]);
    }

    // Remainder
    *r = k[1];
    return 0;
}
/**
 * @brief Counts the number of 1 bits (population count) in a 128-bit unsigned integer.
 * @param x input value.
 * @return Number of 1 bits in x.
 */
[[nodiscard]] FP128_INLINE constexpr uint64_t popcnt128(uint64_t l, uint64_t h) noexcept
{
    return popcnt64(l) + popcnt64(h);
}
/**
 * @brief Left zero count 128 bit
 * @param l Low QWORD
 * @param h High QWORD
 * @return Left zero count
 */
[[nodiscard]] FP128_INLINE constexpr uint64_t lzcnt128(uint64_t l, uint64_t h) noexcept
{
    return (h != 0) ? lzcnt64(h) : 64 + lzcnt64(l);
}
/**
 * @brief Calculates the Log base 2 of x: log2(x)
 * Rounding is always towards zero so the maximum error is close to 1.
 * @param l Lower QWORD of the value
 * @param h High QWORD of the value
 * @return log2(x). Returns zero when x is zero.
 */
[[nodiscard]] FP128_INLINE constexpr uint64_t log2(uint64_t l, uint64_t h) noexcept
{
    return (h != 0 || l != 0) ? 127 - lzcnt128(l, h) : 0;
}
/**
 * @brief Calculates the Log base 2 of x: log2(x)
 * Rounding is always towards zero so the maximum error is close to 1.
 * @param x The number to perform log2 on.
 * @return log2(x). Returns zero when x is zero.
 */
[[nodiscard]] FP128_INLINE constexpr uint64_t log2(uint64_t x) noexcept
{
    return (x) ? 63ull - lzcnt64(x) : 0;
}
/**
 * @brief Calculates the Log base 2 of x: log2(x)
 * Rounding is always towards zero so the maximum error is close to 1.
 * @param x The number to perform log2 on.
 * @return log2(x). Returns zero when x is zero.
 */
[[nodiscard]] FP128_INLINE constexpr uint32_t log2(uint32_t x) noexcept
{
    return (x) ? 31ull - lzcnt32(x) : 0;
}

}  // namespace fp128

#endif  // FP128_FIXED_POINT128_SHARED_H
