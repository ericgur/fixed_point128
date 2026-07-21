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
 * Hacker's Delight), IEEE 754 bit-layout unions, and build configuration macros.
 *
 * This header is consumed by fixed_point128.h and should not be included directly.
 */

#include <cstdint>
#include <cassert>
#include <stdexcept>
#include <memory>
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

static constexpr bool FP128_CPP_STYLE_MODULO = true;             ///< Use C++ modulo semantics (false = Python-style).
static constexpr bool FP128_USE_RECIPROCAL_FOR_DIVISION = true;  ///< Use reciprocal approximation for division.

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

#define lzcnt32 __lzcnt
#define lzcnt64 __lzcnt64
#define udiv64 _udiv64
#define udiv128 _udiv128
#define mulx_u64 _mulx_u64
#define addcarryx_u64 _addcarryx_u64
#define popcnt32 __popcnt
#define popcnt64 __popcnt64
#define alloca _alloca

//
// GCC/Clang portable fallback implementations
//
#elif defined(FP128_CLANG)
#define __shiftright128 shift_right128
#define __shiftleft128 shift_left128

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
FP128_FORCE_INLINE static uint64_t mulx_u64(uint64_t a, uint64_t b, uint64_t* hi) noexcept
{
    FP128_ASSERT(hi != nullptr);  // Caller must provide a valid pointer for the high part. Compatibility with MSVC intrinsic.
#if defined(FP128_ARM64)
    uint64_t lo_res = 0, hi_res = 0;
    __asm__("umulh %[hi], %[a], %[b]\n\t"
            "mul   %[lo], %[a], %[b]"
            : [hi] "=&r"(hi_res), [lo] "=r"(lo_res)  // hi is early-clobber: MUL reads a and b after hi was written
            : [a] "r"(a), [b] "r"(b));
    *hi = hi_res;
    return lo_res;
#else
    __uint128_t r = (__uint128_t)a * b;
    *hi = (uint64_t)(r >> 64);
    return (uint64_t)r;
#endif
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
FP128_FORCE_INLINE static unsigned char addcarryx_u64(unsigned char c, uint64_t a, uint64_t b, uint64_t* out) noexcept
{
    FP128_ASSERT(out != nullptr);  // Caller must provide a valid pointer for the result. Compatibility with MSVC intrinsic.
#if defined(FP128_ARM64)
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
#else
    __uint128_t r = (__uint128_t)a + b + c;
    *out = (uint64_t)r;
    return (unsigned char)(r >> 64);
#endif
}

/**
 * @brief Count leading zeros in a 32-bit value (GCC/Clang fallback).
 *
 * The AArch64 CLZ instruction is defined for a zero operand (it returns the operand width),
 * so unlike the x86 BSR based lowering no zero test is needed.
 */
FP128_FORCE_INLINE static uint32_t lzcnt32(uint32_t x) noexcept
{
#if defined(FP128_ARM64)
    uint32_t res = 0;
    __asm__("clz %w[res], %w[val]" : [res] "=r"(res) : [val] "r"(x));
    return res;
#else
    return (x == 0) ? 32u : (uint32_t)__builtin_clz(x);
#endif
}

/**
 * @brief Count leading zeros in a 64-bit value (GCC/Clang fallback).
 *
 * The AArch64 CLZ instruction returns 64 for a zero operand, no zero test is needed.
 */
FP128_FORCE_INLINE static uint64_t lzcnt64(uint64_t x) noexcept
{
#if defined(FP128_ARM64)
    uint64_t res = 0;
    __asm__("clz %[res], %[val]" : [res] "=r"(res) : [val] "r"(x));
    return res;
#else
    return (x == 0) ? 64u : (uint64_t)__builtin_clzll(x);
#endif
}

/**
 * @brief Count set bits in a 64-bit value (GCC/Clang fallback).
 *
 * AArch64 has no general purpose register population count before ARMv8.9 (FEAT_CSSC), which
 * Apple Silicon does not implement. The value is moved to a NEON register, CNT counts the bits
 * of each of the 8 bytes in parallel and ADDV sums the 8 byte lanes into a single byte.
 */
FP128_FORCE_INLINE static uint64_t popcnt64(uint64_t x) noexcept
{
#if defined(FP128_ARM64)
    uint32_t res = 0;
    uint8x8_t tmp;
    __asm__("fmov %d[tmp], %[val]\n\t"
            "cnt  %[tmp].8b, %[tmp].8b\n\t"
            "addv %b[tmp], %[tmp].8b\n\t"
            "fmov %w[res], %s[tmp]"
            : [res] "=r"(res), [tmp] "=&w"(tmp)
            : [val] "r"(x));
    return res;
#else
    return (uint64_t)__builtin_popcountll(x);
#endif
}

/**
 * @brief Count set bits in a 32-bit value (GCC/Clang fallback).
 *
 * Same NEON sequence as popcnt64. FMOV of a W register zeroes the upper lanes of the NEON
 * register, so the 4 unused bytes contribute nothing to the ADDV sum.
 */
FP128_FORCE_INLINE static uint32_t popcnt32(uint32_t x) noexcept
{
#if defined(FP128_ARM64)
    uint32_t res = 0;
    uint8x8_t tmp;
    __asm__("fmov %s[tmp], %w[val]\n\t"
            "cnt  %[tmp].8b, %[tmp].8b\n\t"
            "addv %b[tmp], %[tmp].8b\n\t"
            "fmov %w[res], %s[tmp]"
            : [res] "=r"(res), [tmp] "=&w"(tmp)
            : [val] "r"(x));
    return res;
#else
    return (uint32_t)__builtin_popcount(x);
#endif
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
#if defined(FP128_MSVC)
#pragma warning(push)
#pragma warning(disable : 4201)  // nameless union/structs
#endif

/**
 * @struct Double
 * @brief Union for accessing IEEE 754 double-precision bit fields.
 *
 * Allows direct access to the mantissa, exponent, and sign of a double
 * without manual bit shifting.
 */
struct Double {
    Double(double v = 0) noexcept : val(v) {}
    union {
        struct {
            uint64_t f : dbl_frac_bits;  ///< Mantissa (fraction) bits.
            uint64_t e : dbl_exp_bits;   ///< Biased exponent bits.
            uint64_t s : 1;              ///< Sign bit (0 = positive, 1 = negative).
        };
        double val;  ///< The raw double value.
    };
};
static_assert(sizeof(Double) == sizeof(double), "The Double union should have the same size as a double variable!");

/**
 * @struct Float
 * @brief Union for accessing IEEE 754 single-precision bit fields.
 *
 * Allows direct access to the mantissa, exponent, and sign of a float
 * without manual bit shifting.
 */
struct Float {
    Float(float v = 0) noexcept : val(v) {}
    union {
        struct {
            uint32_t f : flt_frac_bits;  ///< Mantissa (fraction) bits.
            uint32_t e : flt_exp_bits;   ///< Biased exponent bits.
            uint32_t s : 1;              ///< Sign bit (0 = positive, 1 = negative).
        };
        float val;  ///< The raw float value.
    };
};
static_assert(sizeof(Float) == sizeof(float), "The Float union should have the same size as a float variable!");

#if defined(FP128_MSVC)
#pragma warning(pop)
#endif

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
[[nodiscard]] FP128_INLINE uint64_t shift_right64_round(uint64_t x, int shift) noexcept
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
FP128_INLINE void shift_right128_inplace(uint64_t& l, uint64_t& h, int shift) noexcept
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
FP128_INLINE void shift_left128_inplace(uint64_t& l, uint64_t& h, int shift) noexcept
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
FP128_INLINE void shift_right128_inplace_safe(uint64_t& l, uint64_t& h, int shift) noexcept
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
FP128_INLINE void shift_left128_inplace_safe(uint64_t& l, uint64_t& h, int shift) noexcept
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
template <int shift> [[nodiscard]] FP128_INLINE uint64_t shift_right128(uint64_t l, uint64_t h) noexcept
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
template <int shift> [[nodiscard]] FP128_INLINE uint64_t shift_left128(uint64_t l, uint64_t h) noexcept
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
[[nodiscard]] FP128_INLINE uint64_t shift_right128(uint64_t l, uint64_t h, int shift) noexcept
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
[[nodiscard]] FP128_FORCE_INLINE uint64_t shift_right128_round(uint64_t l, uint64_t h, int shift) noexcept
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
[[nodiscard]] FP128_FORCE_INLINE uint64_t shift_left128(uint64_t l, uint64_t h, int shift) noexcept
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
FP128_INLINE void twos_complement128(uint64_t& l, uint64_t& h) noexcept
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
[[nodiscard]] FP128_INLINE uint64_t popcnt128(uint64_t l, uint64_t h) noexcept
{
    return popcnt64(l) + popcnt64(h);
}
/**
 * @brief Left zero count 128 bit
 * @param l Low QWORD
 * @param h High QWORD
 * @return Left zero count
 */
[[nodiscard]] FP128_INLINE uint64_t lzcnt128(uint64_t l, uint64_t h) noexcept
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
[[nodiscard]] FP128_INLINE uint64_t log2(uint64_t l, uint64_t h) noexcept
{
    return (h != 0 || l != 0) ? 127 - lzcnt128(l, h) : 0;
}
/**
 * @brief Calculates the Log base 2 of x: log2(x)
 * Rounding is always towards zero so the maximum error is close to 1.
 * @param x The number to perform log2 on.
 * @return log2(x). Returns zero when x is zero.
 */
[[nodiscard]] FP128_INLINE uint64_t log2(uint64_t x) noexcept
{
    return (x) ? 63ull - lzcnt64(x) : 0;
}
/**
 * @brief Calculates the Log base 2 of x: log2(x)
 * Rounding is always towards zero so the maximum error is close to 1.
 * @param x The number to perform log2 on.
 * @return log2(x). Returns zero when x is zero.
 */
[[nodiscard]] FP128_INLINE uint32_t log2(uint32_t x) noexcept
{
    return (x) ? 31ull - lzcnt32(x) : 0;
}

}  // namespace fp128

#endif  // FP128_FIXED_POINT128_SHARED_H
