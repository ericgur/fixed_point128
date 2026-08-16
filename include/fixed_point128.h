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

/***********************************************************************************
                                Acknologements
    The function div_32bit is derived from the book "Hacker's Delight" 2nd Edition
    by Henry S. Warren Jr. It was converted to 32 bit operations + a bugfix.

    The functions log, log2, log10 are derived from Dan Moulding's code:
    https://github.com/dmoulding/log2fix

    The function sqrt is based on the book "Math toolkit for real time programming"
    by Jack W. Crenshaw. The sin/cos/atan functions use some ideas from the book.

************************************************************************************/

/**
 * @file fixed_point128.h
 * @brief Header-only 128-bit fixed-point arithmetic library.
 *
 * Provides the `fixed_point128<I>` template class where I is the number of
 * integer bits (range [1, 63]), with one bit used for the sign and the
 * remaining 127-I for the fractional part. Supports full arithmetic (+, -, *, /, %, shifts),
 * comparison operators, type conversions (int, float, double, string), and
 * a comprehensive set of math functions (sqrt, sin, cos, tan, atan, atan2,
 * sinh, cosh, tanh, exp, log, pow, etc.).
 *
 * Platform-specific intrinsics are used for MSVC and GCC/Clang to maximize
 * performance. All methods are inline. The library is 64-bit only.
 *
 * @tparam I Number of integer bits in [1, 63]. The fraction uses 127-I bits and the sign takes one.
 *
 * @see fp128_shared.h for supporting intrinsics and utilities.
 */

#pragma once

// override some static analysis checks
#if defined(_MSC_VER)
#pragma warning(push)
#pragma warning(disable : 26472)  // Don't use a static_cast for arithmetic conversions. Use brace initialization
#pragma warning(disable : 26485)  // No array to pointer decay
#pragma warning(disable : 26481)  // Don't use pointer arithmetic. Use span instead
#pragma warning(disable : 26446)  // Prefer to use gsl::at() instead of unchecked subscript operator
#pragma warning(disable : 26482)  // Only index into arrays using constant expressions
#pragma warning(disable : 26408)  // Avoid malloc() and free(), prefer the nothrow version of new with delete
#endif

#ifndef FP128_FIXED_POINT128_T_H
#define FP128_FIXED_POINT128_T_H

#include <format>      // std::formatter
#include <functional>  // std::hash
#include <limits>
#include <istream>
#include <ostream>
#include "fp128_shared.h"

namespace fp128
{

/***********************************************************************************
 *                                  Forward declarations
 ************************************************************************************/
class fp128_gtest;                          ///< Google test class (friend).
template <int32_t I> class fixed_point128;  ///< Forward declaration of the main template class.

// Note: Release builds will fail without these forward declarations. Hints towards compiler a bug (VS2022 v17.4)
// The compiler and IntelliSense don't match these functions in some cases and try to use the CRT versions which
// causes a compilation error.

/// @name CRT-Style Math Functions (Forward Declarations)
/// @brief Free functions providing standard math library equivalents for fixed_point128.
/// @{
template <int32_t I> [[nodiscard]] fixed_point128<I> fabs(const fixed_point128<I>& x) noexcept;
template <int32_t I> [[nodiscard]] fixed_point128<I> floor(const fixed_point128<I>& x) noexcept;
template <int32_t I> [[nodiscard]] fixed_point128<I> ceil(const fixed_point128<I>& x) noexcept;
template <int32_t I> [[nodiscard]] fixed_point128<I> trunc(const fixed_point128<I>& x) noexcept;
template <int32_t I> [[nodiscard]] fixed_point128<I> round(const fixed_point128<I>& x) noexcept;
template <int32_t I> [[nodiscard]] int32_t ilogb(const fixed_point128<I>& x) noexcept;
template <int32_t I> [[nodiscard]] fixed_point128<I> copysign(const fixed_point128<I>& x, const fixed_point128<I>& y) noexcept;
template <int32_t I> [[nodiscard]] fixed_point128<I> fmod(const fixed_point128<I>& x, const fixed_point128<I>& y);
template <int32_t I> [[nodiscard]] fixed_point128<I> modf(const fixed_point128<I>& x, fixed_point128<I>* iptr) noexcept;
template <int32_t I> [[nodiscard]] fixed_point128<I> fdim(const fixed_point128<I>& x, const fixed_point128<I>& y) noexcept;
template <int32_t I> [[nodiscard]] fixed_point128<I> fmin(const fixed_point128<I>& x, const fixed_point128<I>& y) noexcept;
template <int32_t I> [[nodiscard]] fixed_point128<I> fmax(const fixed_point128<I>& x, const fixed_point128<I>& y) noexcept;
template <int32_t I> [[nodiscard]] fixed_point128<I> hypot(const fixed_point128<I>& x, const fixed_point128<I>& y) noexcept;
template <int32_t I> [[nodiscard]] fixed_point128<I> sqrt(const fixed_point128<I>& x, uint32_t iterations = 3) noexcept;
template <int32_t I> [[nodiscard]] fixed_point128<I> sin(fixed_point128<I> x) noexcept;
template <int32_t I> [[nodiscard]] fixed_point128<I> asin(fixed_point128<I> x) noexcept;
template <int32_t I> [[nodiscard]] fixed_point128<I> cos(fixed_point128<I> x) noexcept;
template <int32_t I> [[nodiscard]] fixed_point128<I> acos(fixed_point128<I> x) noexcept;
template <int32_t I> [[nodiscard]] fixed_point128<I> tan(fixed_point128<I> x) noexcept;
template <int32_t I> [[nodiscard]] fixed_point128<I> atan(fixed_point128<I> x) noexcept;
template <int32_t I> [[nodiscard]] fixed_point128<I> atan2(fixed_point128<I> y, fixed_point128<I> x) noexcept;
template <int32_t I> [[nodiscard]] fixed_point128<I> sinh(const fixed_point128<I>& x) noexcept;
template <int32_t I> [[nodiscard]] fixed_point128<I> asinh(const fixed_point128<I>& x) noexcept;
template <int32_t I> [[nodiscard]] fixed_point128<I> cosh(const fixed_point128<I>& x) noexcept;
template <int32_t I> [[nodiscard]] fixed_point128<I> acosh(const fixed_point128<I>& x) noexcept;
template <int32_t I> [[nodiscard]] fixed_point128<I> tanh(const fixed_point128<I>& x) noexcept;
template <int32_t I> [[nodiscard]] fixed_point128<I> atanh(const fixed_point128<I>& x) noexcept;
template <int32_t I> [[nodiscard]] fixed_point128<I> exp(const fixed_point128<I>& x) noexcept;
template <int32_t I> [[nodiscard]] fixed_point128<I> exp2(const fixed_point128<I>& x) noexcept;
template <int32_t I> [[nodiscard]] fixed_point128<I> expm1(const fixed_point128<I>& x) noexcept;
template <int32_t I> [[nodiscard]] fixed_point128<I> pow(const fixed_point128<I>& x, const fixed_point128<I>& y) noexcept;
template <int32_t I> [[nodiscard]] fixed_point128<I> log(fixed_point128<I> x) noexcept;
template <int32_t I> [[nodiscard]] fixed_point128<I> log2(fixed_point128<I> x) noexcept;
template <int32_t I> [[nodiscard]] fixed_point128<I> log10(fixed_point128<I> x) noexcept;
template <int32_t I> [[nodiscard]] fixed_point128<I> logb(fixed_point128<I> x) noexcept;
template <int32_t I> [[nodiscard]] fixed_point128<I> log1p(fixed_point128<I> x) noexcept;
/// @}

/// @name Non-CRT Utility Functions (Forward Declarations)
/// @{
template <int32_t I> [[nodiscard]] uint64_t lzcnt128(const fixed_point128<I>& x) noexcept;
template <int32_t I> [[nodiscard]] fixed_point128<I> reciprocal(const fixed_point128<I>& x) noexcept;
template <int32_t I> void fact_reciprocal(int x, fixed_point128<I>& res) noexcept;
/// @}

/***********************************************************************************
 *                                  Main Code
 ************************************************************************************/

/**
 * @brief 128 bit fixed point class template.
 *
 * This template provides a floating point-like type that performs various math operations quickly compared to traditional high precision libraries.
 * The only template parameter <B>I</B> is the bit count of the integer part, not counting the sign.
 * The sign takes one bit and the fraction takes the remaining 127 - <B>I</B>.<BR>
 * <B>I</B> is limited to the range [1,63] in order to simplify the implementation and increase performance.<BR>
 * This restriction is enforced at compile time.
 * All of fixed_point128's methods are inline for maximum performance.
 *
 * <B>Representation:</B><BR>
 * The object is exactly 128 bits wide and holds nothing but the value. Those bits are a two's
 * complement integer, held as the QWORD pair (high:low) with the sign in the MSB of <B>high</B>,
 * and the value they stand for is that integer divided by 2^<B>F</B>:
 * @code
 *     value = (int128)(high:low) / 2^F,  where F = 127 - I
 * @endcode
 * The representable range is therefore [-2^I, 2^I - 2^-F], which is the asymmetric range of every
 * two's complement type: the most negative value has no positive counterpart, so negating it (and
 * fabs of it) wraps back to itself.
 *
 * <B>Minimum integer bits:</B><BR>
 * Some of the constants and functions cannot work in the narrowest instantiations, and each one
 * says so with a static_assert in its body rather than by overflowing quietly or failing to link.
 * The bound is the smallest <B>I</B> that works, not a safe margin:
 * <UL>
 * <LI><B>I >= 2</B>: pi(), e(), and the exponential family exp(), exp2(), expm1(), pow() and
 *     tanh(), which reach e() or - for tanh - form e^x + e^-x, never below 2.</LI>
 * <LI><B>I >= 3</B>: pi2(), which is 6.28.</LI>
 * <LI><B>I >= 4</B>: the trigonometric group sin1(), cos1(), sin(), cos(), tan(), asin(), acos(),
 *     atan(), atan2(), and the hyperbolic sinh() and cosh().</LI>
 * </UL>
 * Everything else works for every <B>I</B>, including half_pi(), sqrt_2() and golden_ratio(), which
 * are all below two. A function whose <I>result</I> leaves the range - exp() of a large argument,
 * log() of a tiny one - is a matter of the argument rather than of the type, and overflows silently
 * like any other operation.
 *
 * <B>Implementation notes:</B>
 * <UL>
 * <LI>Overflow is handled silently, similar to builtin integer operations.</LI>
 * <LI>The sign is part of the value rather than a separate field, so addition, subtraction, the
 *     comparisons and the bitwise operators are the plain 128 bit integer operations and there is
 *     no negative zero to normalize away.</LI>
 * <LI>Multiplication and division work on magnitudes internally and restore the sign afterwards.</LI>
 * <LI>A fixed_point128 object is not thread safe. Accessing a const object from multiple threads is safe.</LI>
 * <LI>fixed_point128 is <B>conditionally safe</B>, 2 different non const objects can be accessed concurrently.</LI>
 * <LI>Only 64 bit builds are supported.</LI>
 * </UL>
 *
 * <B>Compile time evaluation:</B><BR>
 * Everything that stays within integer arithmetic is constexpr:
 * <UL>
 * <LI>Construction from and conversion to the integer types and the floating point types, copy,
 *     move and the cross-template conversions between two fixed_point128 instantiations.</LI>
 * <LI>Addition, subtraction, multiplication, square(), the increment and decrement operators,
 *     the shifts, the bitwise operators, the unary operators and the comparisons.</LI>
 * <LI>The queries: is_int(), is_positive(), is_negative(), is_zero(), get_bit() and
 *     get_exponent(), and the exact constants one(), half() and epsilon().</LI>
 * <LI>The math functions fabs, floor, ceil, trunc, round, copysign, modf, fmin, fmax, fdim, sqr,
 *     ilogb, lzcnt128, log2 and logb.</LI>
 * </UL>
 *
 * The bit counting and extended arithmetic intrinsics these rest on are not constant expressions,
 * so fp128_shared.h wraps each one in a constexpr function that serves a constant
 * evaluated call from a portable implementation of the same operation. A runtime call still
 * reaches the bare intrinsic and generates the same code it did before.
 *
 * The rest cannot be constexpr, for one of two reasons:
 * <UL>
 * <LI>Division and modulo: div_32bit needs alloca and a goto, neither of which C++20 permits in
 *     a constexpr function, and div_64bit rests on the _udiv128 intrinsic.</LI>
 * <LI>Function local statics, which a constexpr function may not declare. The transcendental
 *     constants (pi, e, sqrt_2, ...) and the factorial reciprocal table are parsed from strings
 *     into them, and reciprocal holds its bounds the same way. sqrt is doubly out, it calls the
 *     CRT's sqrt() for its initial estimate.</LI>
 * </UL>
 */
template <int32_t I> class fixed_point128
{
    // build time validation of template parameters
    static_assert(1 <= I && I <= 63, "Template parameter <I> must be in the [1,63] range!");
    static_assert(sizeof(void*) == 8, "fixed_point128 is supported in 64 bit builds only!");

    // friends
    template <int32_t I2> friend class fixed_point128;  // this class is a friend of all its template instances. Avoids awkward getter/setter functions.
    friend class fp128_gtest;

private:
    /// @name Data Members
    /// @{
    uint64_t low;   ///< Lower 64-bit QWORD of the 128-bit two's complement value.
    uint64_t high;  ///< Upper 64-bit QWORD of the 128-bit two's complement value. Bit 63 holds the sign.
    /// @}

public:
    static constexpr int32_t F = 127 - I;  ///< Fraction bit count (127 - I). The remaining bit is the sign.

private:
    /// @name Internal Constants
    /// @{
    static constexpr int32_t upper_frac_bits = F - 64;                   ///< Fraction bits residing in the upper QWORD (63 - I).
    static constexpr uint64_t unity = 1ull << upper_frac_bits;           ///< Upper QWORD value representing 1.0.
    static inline const double upper_unity = ::pow(2, 64 - F);           ///< Scale factor for upper QWORD to double conversion.
    static inline const double lower_unity_l = ::pow(2, -F);             ///< Scale factor for lower QWORD (low 32 bits) to double.
    static inline const double lower_unity_h = ::pow(2, 32 - F);         ///< Scale factor for lower QWORD (high 32 bits) to double.
    static constexpr uint64_t int_mask = UINT64_MAX << upper_frac_bits;  ///< Bitmask for the sign and integer bits in the upper QWORD.
    static constexpr int32_t max_frac_digits = (int)(F / 3.3);           ///< Maximum meaningful base-10 fractional digits.
    static constexpr uint64_t sign_mask = 1ull << 63;                    ///< Bitmask of the sign bit within the upper QWORD.
    /// @}

    /**
     * @brief Reads the magnitude of this object into a QWORD pair.
     *
     * The multiply, divide and floating point conversions all need |x| as an unsigned 128 bit
     * value; this is the one place that knows how to produce it. The most negative value is its own
     * two's complement, and 2^127 is exactly what its magnitude is as an unsigned pair, so that
     * case needs no special handling.
     *
     * The sign is returned as a mask rather than a bool so that two of them combine with xor. The
     * negation is left as a branch on purpose: two branchless spellings were measured, the two's
     * complement as (x ^ -1) + 1 against (x ^ 0) + 0, and computing both forms for a conditional
     * move to pick between. Both are slower - by 30% and 20% respectively on the multiply and on
     * everything built on it - MSVC putting the pair through memory rather than keeping it in
     * registers.
     *
     * @param l Receives the low QWORD of the magnitude
     * @param h Receives the high QWORD of the magnitude
     * @return The sign as a mask, all ones when the object was negative and all zeros otherwise.
     *         Two of them can be combined with xor to give the sign of a quotient, and NegateIf()
     *         takes one.
     */
    [[nodiscard]] FP128_FORCE_INLINE constexpr uint64_t GetMagnitude(uint64_t& l, uint64_t& h) const noexcept
    {
        l = low;
        h = high;
        const uint64_t sign_bits = static_cast<uint64_t>(static_cast<int64_t>(high) >> 63);
        if (sign_bits != 0)
            twos_complement128(l, h);

        return sign_bits;
    }

    /**
     * @brief Negates this object in place, in the two's complement sense.
     * The most negative value has no positive counterpart and is left unchanged, which is the same
     * silent wrap the builtin integer types perform.
     */
    FP128_FORCE_INLINE constexpr void Negate() noexcept { twos_complement128(low, high); }

    /**
     * @brief Negates this object when the mask says to, without branching.
     * The counterpart of GetMagnitude(), which produces the masks this takes.
     * @param sign_bits All ones to negate, all zeros to leave the object alone.
     */
    FP128_FORCE_INLINE constexpr void NegateIf(uint64_t sign_bits) noexcept
    {
        if (sign_bits != 0)
            twos_complement128(low, high);
    }

    /**
     * @brief Replaces this object with its magnitude, in place.
     * What the division does to its numerator: the work below it is unsigned and the sign goes back
     * on at the end.
     * @return The sign as a mask, in the form NegateIf() takes.
     */
    [[nodiscard]] FP128_FORCE_INLINE constexpr uint64_t TakeMagnitude() noexcept
    {
        const uint64_t sign_bits = static_cast<uint64_t>(static_cast<int64_t>(high) >> 63);
        NegateIf(sign_bits);
        return sign_bits;
    }

    /**
     * @brief Extracts the mantissa of a floating point representation of this object, rounded.
     *
     * Keeps the frac_bits+1 most significant bits of the magnitude and rounds the dropped
     * remainder to nearest, ties going to even. This is what IEEE 754 mandates and what the
     * builtin conversions do.<BR>
     * Rounding needs three pieces of information: the bits that survive, the highest dropped bit
     * (which chooses the direction) and whether anything below it is set (the sticky bit, which
     * distinguishes an exact tie from a remainder larger than half). The caller passes the
     * magnitude, so the sign plays no part here.
     *
     * @param low Low QWORD of the magnitude.
     * @param high High QWORD of the magnitude.
     * @param msb Bit position of the most significant set bit. Must be larger than frac_bits.
     * @param frac_bits Fraction bit count of the target type, 23 for float or 52 for double.
     * @return The rounded mantissa. Holds frac_bits+1 bits, or frac_bits+2 when rounding carried
     *         into the next power of 2.
     */
    [[nodiscard]] FP128_INLINE static constexpr uint64_t RoundedMantissa(uint64_t low, uint64_t high, int32_t msb, int32_t frac_bits) noexcept
    {
        const int32_t shift = msb - frac_bits;  // count of dropped bits, at least 1
        uint64_t mant = shift_right128(low, high, shift);

        // the highest dropped bit sits just below the mantissa
        const uint64_t round_bit = (shift <= 64) ? FP128_GET_BIT(low, shift - 1) : FP128_GET_BIT(high, shift - 65);
        if (round_bit == 0)
            return mant;  // remainder is below half, round down

        // is anything set below the rounding bit? Each branch below shifts the bits of interest
        // to the top of the QWORD, which avoids building a mask of a variable width.
        bool sticky;
        if (shift < 2) {
            sticky = false;  // the rounding bit is the only dropped bit
        } else if (shift <= 64) {
            sticky = (low << (65 - shift)) != 0;  // bits [shift-2:0] of low
        } else if (shift == 65) {
            sticky = low != 0;  // all of low, nothing of high
        } else {
            sticky = low != 0 || (high << (129 - shift)) != 0;  // all of low plus bits [shift-66:0] of high
        }

        // above half always rounds up, an exact tie rounds towards the even mantissa
        if (sticky || (mant & 1))
            ++mant;

        return mant;
    }

public:
    /**
     * @brief Reads the raw representation: the two QWORDs of the 128 bit two's complement value.
     * @param l Receives the low QWORD
     * @param h Receives the high QWORD, sign included
     */
    FP128_FORCE_INLINE constexpr void get_components(uint64_t& l, uint64_t& h) const noexcept
    {
        l = low;
        h = high;
    }

    static constexpr uint64_t max_int_value = FP128_MAX_VALUE_64(I);  ///< Largest representable positive integer, 2^I - 1.
    typedef fixed_point128<I> type;                                   ///< Self type alias.
    typedef fixed_point128<I>* ptr_type;                              ///< Pointer type alias.
    typedef fixed_point128<I>& ref_type;                              ///< Reference type alias.

    /// @name Constructors
    /// @{

    /**
     * @brief Default constructor, creates an instance with a value of zero.
     */
    FP128_FORCE_INLINE constexpr fixed_point128() noexcept : low(0), high(0) {}
    /**
     * @brief Copy constructor
     * @param rhs Object to copy from
     */
    FP128_FORCE_INLINE constexpr fixed_point128(const fixed_point128& rhs) noexcept : low(rhs.low), high(rhs.high) {}
    /**
     * @brief cross-template Copy constructor, can be used between two different fixed_point128 templates
     * A value too large for the destination wraps silently, as everywhere else in this class.
     * @param rhs fixed_point128 instance with from a different template instance.
     */
    template <int32_t I2> FP128_FORCE_INLINE constexpr fixed_point128(const fixed_point128<I2>& rhs) noexcept
    {
        if constexpr (I == I2) {
            high = rhs.high;
            low = rhs.low;
        }
        // rhs has less integer bits and more fraction bits
        else if constexpr (I < I2) {
            // shift left by I2 - I bits
            constexpr int shift = I2 - I;
            low = rhs.low << shift;
            high = shift_left128(rhs.low, rhs.high, shift);
        }
        // rhs has more integer bits and less fraction bits
        else {  // I > I2
            // shift right by I - I2 bits, replicating the sign bit so a negative value stays negative
            constexpr int shift = I - I2;
            low = rhs.low;
            high = rhs.high;
            shift_right128_inplace_safe_signed(low, high, shift);
        }
    }

    /**
     * @brief Move constructor
     * Doesn't modify the right hand side object. Acts like a copy constructor.
     * @param rhs Object to copy from
     */
    FP128_FORCE_INLINE constexpr fixed_point128(fixed_point128&& rhs) noexcept : low(rhs.low), high(rhs.high) {}
    /**
     * @brief Constructor from the double type
     * Underflow goes to zero. Overflow, NaN and +-INF saturate to the extreme value carrying the
     * sign of the input: max() for a positive one and lowest() for a negative one.
     * Not constexpr: the exponent and mantissa are read through the Double union, and reading the
     * inactive member of a union is not allowed during constant evaluation.
     * @param x Input value
     */
    FP128_INLINE constexpr fixed_point128(double x) noexcept
    {
        // very common case
        if (x == 0) {
            low = high = 0;
            return;
        }

        // hack the double bit fields
        const Double d(x);

        const uint32_t negative = static_cast<uint32_t>(d.s());
        const int32_t e = static_cast<int32_t>(d.e()) - 1023;
        uint64_t f = d.f();

        // overflow which also catches NaN and Inf
        if (e >= I) {
            low = (negative) ? 0 : UINT64_MAX;
            high = (negative) ? sign_mask : ~sign_mask;
            return;
        }

        // normal number, produces non zero value
        if (e >= -F) {
            // bit 52 in f is the unity value of the float. it needs to move to the unity position in fixed point
            f |= FP128_ONE_SHIFT(dbl_frac_bits);
            // the leading one of f sits at bit 52 and belongs at bit F + e, which is this many bits
            // into the upper QWORD. The magnitude is built first and negated below if needed.
            int32_t bits_to_shift = F - 64 - dbl_frac_bits + e;

            // f fits in high QWORD
            if (bits_to_shift >= 0) {
                high = f << bits_to_shift;
                low = 0;
            }
            // shift right
            else {
                bits_to_shift = -bits_to_shift;
                // f has some bits in high QWORD
                if (bits_to_shift <= 53) {
                    high = f >> bits_to_shift;
                    low = FP128_GET_BITS(f, 0, bits_to_shift);
                    low <<= 64ll - bits_to_shift;
                }
                // shift f into low QWORD
                else {
                    high = 0;
                    // f has bit 52 high, shift it left to moves to bit 63
                    f <<= 63 - dbl_frac_bits;
                    bits_to_shift -= 64 - (63 - dbl_frac_bits);
                    low = f >> bits_to_shift;
                }
            }

            if (negative)
                Negate();
        }
        // too small to be represented, no need to bother.
        else {
            high = low = 0;
        }
    }
    /**
     * @brief Constructor from uint64_t type
     * A value of 2^I or above overflows into the sign and wraps silently.
     * @param x Input value
     */
    FP128_FORCE_INLINE constexpr fixed_point128(uint64_t x) noexcept
    {
        low = 0;
        high = x << upper_frac_bits;
    }
    /**
     * @brief Constructor from int64_t type
     * The shift is done in the unsigned domain, where it wraps rather than overflowing, and a
     * negative value carries its sign bits along with it: scaling a two's complement value by a
     * power of two is the same left shift its magnitude needs.
     * @param x Input value
     */
    FP128_FORCE_INLINE constexpr fixed_point128(int64_t x) noexcept
    {
        low = 0;
        high = static_cast<uint64_t>(x) << upper_frac_bits;
    }
    /**
     * @brief Constructor from uint32_t type
     * @param x Input value
     */
    FP128_FORCE_INLINE constexpr fixed_point128(uint32_t x) noexcept
    {
        low = 0;
        high = static_cast<uint64_t>(x) << upper_frac_bits;
    }
    /**
     * @brief Constructor from int32_t type
     * @param x Input value
     */
    FP128_FORCE_INLINE constexpr fixed_point128(int32_t x) noexcept
    {
        low = 0;
        // sign extended to 64 bits first, see the int64_t constructor
        high = static_cast<uint64_t>(static_cast<int64_t>(x)) << upper_frac_bits;
    }
    /**
     * @brief Constructor from const char* (C string).
     * Accurate up to 37 digits after the decimal point.
     * Allows creating very high precision values. Much slower than the rhs constructors.
     * A value too large for the type wraps silently, which unlike the saturation the double
     * constructor performs can land on either sign.
     * @param x Input string
     */
    fixed_point128(const char* x) noexcept
    {
        // 10^-n in fixed_point128<1> form, which holds 126 fraction bits whatever I is
        static fixed_point128<1> base10_table[] = {
            //  {low QWORD,             high QWORD}
            {0x0000000000000000ull, 0x4000000000000000ull},  // 10^0, not used, makes the code simpler
            {0x6666666666666666ull, 0x0666666666666666ull},  // 10^-1
            {0x0A3D70A3D70A3D71ull, 0x00A3D70A3D70A3D7ull},  // 10^-2
            {0xE76C8B4395810625ull, 0x0010624DD2F1A9FBull},  // 10^-3
            {0xCA57A786C226809Dull, 0x0001A36E2EB1C432ull},  // 10^-4
            {0xE108C3F3E0370CDDull, 0x000029F16B11C6D1ull},  // 10^-5
            {0x634DAD31FCD24E16ull, 0x00000431BDE82D7Bull},  // 10^-6
            {0xBD215E1E99483B02ull, 0x0000006B5FCA6AF2ull},  // 10^-7
            {0x461CEFCFDC20D2B3ull, 0x0000000ABCC77118ull},  // 10^-8
            {0x6D694B2E62D01512ull, 0x0000000112E0BE82ull},  // 10^-9
            {0xD7BDBAB7D6AE6882ull, 0x000000001B7CDFD9ull},  // 10^-10
            {0x2F2C92ABFBDE3DA7ull, 0x0000000002BFAFFCull},  // 10^-11
            {0x04B7A84465FC9FC4ull, 0x0000000000465E66ull},  // 10^-12
            {0x9A125DA07099432Dull, 0x0000000000070970ull},  // 10^-13
            {0xDC35095CD80F5385ull, 0x000000000000B424ull},  // 10^-14
            {0xAF9EE756159B21F4ull, 0x0000000000001203ull},  // 10^-15
            {0x2B297D889BC2B698ull, 0x00000000000001CDull},  // 10^-16
            {0x1DEA8C8DA92D1242ull, 0x000000000000002Eull},  // 10^-17
            {0x9C97747490EAE83Aull, 0x0000000000000004ull},  // 10^-18
            {0x760F253EDB4AB0D3ull, 0x0000000000000000ull},  // 10^-19
            {0x0BCE5086492111AFull, 0x0000000000000000ull},  // 10^-20
            {0x012E3B40A0E9B4F8ull, 0x0000000000000000ull},  // 10^-21
            {0x001E392010175EE6ull, 0x0000000000000000ull},  // 10^-22
            {0x000305B66802564Aull, 0x0000000000000000ull},  // 10^-23
            {0x00004D5F0A66A23Bull, 0x0000000000000000ull},  // 10^-24
            {0x000007BCB43D769Full, 0x0000000000000000ull},  // 10^-25
            {0x000000C612062576ull, 0x0000000000000000ull},  // 10^-26
            {0x00000013CE9A36F2ull, 0x0000000000000000ull},  // 10^-27
            {0x00000001FB0F6BE5ull, 0x0000000000000000ull},  // 10^-28
            {0x0000000032B4BDFDull, 0x0000000000000000ull},  // 10^-29
            {0x0000000005121300ull, 0x0000000000000000ull},  // 10^-30
            {0x000000000081CEB3ull, 0x0000000000000000ull},  // 10^-31
            {0x00000000000CFB12ull, 0x0000000000000000ull},  // 10^-32
            {0x0000000000014C4Full, 0x0000000000000000ull},  // 10^-33
            {0x000000000000213Bull, 0x0000000000000000ull},  // 10^-34
            {0x0000000000000353ull, 0x0000000000000000ull},  // 10^-35
            {0x0000000000000055ull, 0x0000000000000000ull},  // 10^-36
            {0x0000000000000009ull, 0x0000000000000000ull},  // 10^-37
            {0x0000000000000001ull, 0x0000000000000000ull}   // 10^-38
        };
        constexpr uint64_t base10_table_size = array_length(base10_table);
        static_assert(max_frac_digits < base10_table_size);

        low = high = 0;
        if (x == nullptr)
            return;

        const auto x_len = 1 + strlen(x);
        // make_unique_for_overwrite throws bad_alloc, which this noexcept constructor must not
        // propagate. The nothrow form reports a failed allocation as a null pointer instead.
        auto str_ptr = std::unique_ptr<char[]>(new (std::nothrow) char[x_len]);
        char* p = str_ptr.get();
        if (p == nullptr)
            return;
        strnlwr(p, x, x_len);

        // skip leading white space. The cast keeps a negative char from reaching isspace, which
        // only accepts values representable as unsigned char (or EOF).
        while (*p && isspace(static_cast<unsigned char>(*p)))
            ++p;

        if (*p == '\0') {
            *this = fixed_point128();
            return;
        }

        // note a negative sign if present. The magnitude is built first and negated at the end.
        bool negative = false;
        if (*p == '-') {
            negative = true;
            ++p;
        } else if (*p == '+')
            ++p;

        char* dec = strchr(p, '.');
        // number is an integer
        if (dec == nullptr) {
            high = std::strtoull(p, nullptr, 10) << upper_frac_bits;
            if (negative)
                Negate();
            return;
        }

        // number is a float, get the integer part using strtoull()
        *dec = '\0';
        const uint64_t int_val = std::strtoull(p, nullptr, 10) << upper_frac_bits;

        p = dec + 1;
        uint32_t digits = 0;
        fixed_point128<1> frac;
        // multiply each digits by 10^-n
        while (++digits < base10_table_size && isdigit(static_cast<unsigned char>(*p))) {
            uint32_t d = static_cast<uint64_t>(p[0] - '0');
            frac += base10_table[digits] * d;
            ++p;
        }

        // the fraction was accumulated with 126 bits of it, this instantiation keeps F of them
        frac >>= (I - 1);
        low = frac.low;
        high = frac.high + int_val;
        if (negative)
            Negate();
    }
    /**
     * @brief Constructor from std::string.
     * Accurate to 37 digits after the decimal point.
     * Allows creating very high precision values. Much slower than the rhs constructors.
     * @param x Input string
     */
    FP128_INLINE fixed_point128(const std::string& x) noexcept
    {
        fixed_point128 temp = x.c_str();
        *this = temp;
    }
    /**
     * @brief Constructor from the two QWORDs of the raw value, useful for creating very small constants.
     * The pair is read as a 128 bit two's complement integer scaled by 2^-F, so a negative value is
     * spelled with the MSB of h set rather than with a separate sign.
     * @param l Low QWORD
     * @param h High QWORD, sign included
     */
    FP128_FORCE_INLINE constexpr fixed_point128(uint64_t l, uint64_t h) noexcept : low(l), high(h) {}

    /**
     * @brief Destructor
     */
    constexpr ~fixed_point128() noexcept = default;
    /**
     * @brief Assignment operator
     * @param rhs Object to copy from
     * @return This object.
     */
    FP128_FORCE_INLINE constexpr fixed_point128& operator=(const fixed_point128& rhs) noexcept
    {
        high = rhs.high;
        low = rhs.low;
        return *this;
    }
    /**
     * @brief Move assignment operator
     * @param rhs Object to copy from
     * @return This object.
     */
    FP128_FORCE_INLINE constexpr fixed_point128& operator=(fixed_point128&& rhs) noexcept
    {
        high = rhs.high;
        low = rhs.low;
        return *this;
    }
    /**
     * @brief cross-template assignment operator, can be used between two different fixed_point128 templates
     * @param rhs fixed_point128 instance with from a different template instance.
     * @return This object.
     */
    template <int32_t I2> FP128_FORCE_INLINE constexpr fixed_point128<I>& operator=(const fixed_point128<I2>& rhs) noexcept
    {
        if constexpr (I == I2) {
            high = rhs.high;
            low = rhs.low;
        }
        // rhs has less integer bits and more fraction bits
        else if constexpr (I < I2) {
            // shift left by I2 - I bits
            constexpr int shift = I2 - I;
            low = rhs.low << shift;
            high = shift_left128(rhs.low, rhs.high, shift);
        }
        // rhs has more integer bits and less fraction bits
        else {  // I > I2
            // shift right by I - I2 bits, replicating the sign bit
            constexpr int shift = I - I2;
            low = rhs.low;
            high = rhs.high;
            shift_right128_inplace_safe_signed(low, high, shift);
        }

        return *this;
    }

    //
    /// @}

    /// @name Conversion Operators
    /// @{

    /**
     * @brief operator int64_t - converts to a int64_t
     * The fraction is truncated towards zero, as the builtin floating point to integer conversions
     * do, rather than shifted away - which would round a negative value towards minus infinity.
     * @return Object value.
     */
    [[nodiscard]] FP128_FORCE_INLINE constexpr operator int64_t() const noexcept
    {
        uint64_t l = 0, h = 0;
        const uint64_t sign_bits = GetMagnitude(l, h);
        const uint64_t magnitude = h >> upper_frac_bits;
        // negated in the unsigned domain, which wraps instead of overflowing for the most negative
        // value of a 64 bit result. (x ^ -1) - -1 is that negation, and (x ^ 0) - 0 leaves it.
        return static_cast<int64_t>((magnitude ^ sign_bits) - sign_bits);
    }
    /**
     * @brief operator uint64_t - converts to a uint64_t
     * A negative value wraps modulo 2^64, the way a builtin signed to unsigned conversion does.
     * @return Object value.
     */
    [[nodiscard]] FP128_FORCE_INLINE constexpr operator uint64_t() const noexcept { return static_cast<uint64_t>(operator int64_t()); }
    /**
     * @brief operator int32_t - converts to a int32_t
     * @return Object value, truncated towards zero and wrapped to 32 bits.
     */
    [[nodiscard]] FP128_FORCE_INLINE constexpr operator int32_t() const noexcept { return static_cast<int32_t>(operator int64_t()); }
    /**
     * @brief operator uint32_t - converts to a uint32_t
     * @return Object value, truncated towards zero and wrapped to 32 bits.
     */
    [[nodiscard]] FP128_FORCE_INLINE constexpr operator uint32_t() const noexcept { return static_cast<uint32_t>(operator int64_t()); }
    /**
     * @brief operator float - converts to a float
     * Not constexpr: the result is assembled through the Float union, and writing one member of a
     * union and reading another is not allowed during constant evaluation.
     * @return Object value.
     */
    [[nodiscard]] FP128_INLINE constexpr operator float() const noexcept
    {
        if (!*this)
            return 0.0f;

        // the conversion is done on the magnitude, the sign is put back by Float::make()
        uint64_t l = 0, h = 0;
        const uint32_t s = static_cast<uint32_t>(GetMagnitude(l, h) & 1);
        const int32_t msb = 127 - static_cast<int32_t>(lzcnt128(l, h));
        const int32_t expo = msb - F;  // base 2 exponent of the value

        // A float only reaches down to 2^-126 before going denormal. The smallest value of any
        // instantiation is 2^-F, and F peaks at 126, so no value here is denormal as a float.

        // the value fits in the fraction, no bits are lost so no rounding is needed.
        // float doesn't hold the msb, it's implicit, so bit [22:0] hold the rest of the value
        if (msb <= flt_frac_bits) {
            return Float::make(s, static_cast<uint32_t>(expo + 127), static_cast<uint32_t>(l << (flt_frac_bits - msb)));
        }

        // more bits than the fraction can hold, drop the extra ones with rounding.
        // rounding up can carry into the next power of 2, costing a fraction bit and
        // incrementing the exponent
        const uint64_t mant = RoundedMantissa(l, h, msb, flt_frac_bits);
        const uint64_t carry = mant >> (flt_frac_bits + 1);
        return Float::make(s, static_cast<uint32_t>(expo + static_cast<int32_t>(carry) + 127), static_cast<uint32_t>(mant >> carry));
    }
    /**
     * @brief operator double - converts to a double
     * The whole value range is within reach of a double's exponent, so nothing overflows or
     * goes denormal.
     * Not constexpr: the result is assembled through the Double union, and writing one member of a
     * union and reading another is not allowed during constant evaluation.
     * @return Object value, rounded to nearest with ties going to even.
     */
    [[nodiscard]] FP128_INLINE constexpr operator double() const noexcept
    {
        if (!*this)
            return 0.0;

        // the conversion is done on the magnitude, the sign is put back by Double::make()
        uint64_t l = 0, h = 0;
        const uint64_t s = GetMagnitude(l, h) & 1;
        const int32_t msb = 127 - static_cast<int32_t>(lzcnt128(l, h));
        const int32_t expo = msb - F;  // base 2 exponent of the value

        // the value fits in the fraction, no bits are lost so no rounding is needed.
        // double doesn't hold the msb, it's implicit, so bit [51:0] hold the rest of the value
        if (msb <= dbl_frac_bits) {
            return Double::make(s, static_cast<uint64_t>(expo + 1023), l << (dbl_frac_bits - msb));
        }

        // more bits than the fraction can hold, drop the extra ones with rounding.
        // rounding up can carry into the next power of 2, costing a fraction bit and
        // incrementing the exponent
        const uint64_t mant = RoundedMantissa(l, h, msb, dbl_frac_bits);
        const uint64_t carry = mant >> (dbl_frac_bits + 1);
        return Double::make(s, static_cast<uint64_t>(expo + static_cast<int32_t>(carry) + 1023), mant >> carry);
    }
    /**
     * @brief operator long double - converts to a long double
     * @return Object value.
     */
    [[nodiscard]] FP128_FORCE_INLINE constexpr operator long double() const noexcept { return operator double(); }
    /**
     * @brief Converts to a std::string (slow) string holds all meaningful fraction bits.
     * @return object string representation
     */
    [[nodiscard]] FP128_FORCE_INLINE operator std::string() const noexcept { return operator char*(); }
    /**
     * @brief Converts to a C string (slow) string holds all meaningful fraction bits.
     * @return object string representation
     */
    [[nodiscard]] explicit FP128_INLINE operator char*() const noexcept
    {
        static thread_local char str[128];  // need roughly a (meaningful) decimal digits per 3.2 bits

        char* p = &str[0];
        // the digits are produced from the magnitude, the sign is written out separately
        uint64_t l = 0, h = 0;
        if (GetMagnitude(l, h))
            *p++ = '-';

        uint64_t integer = FP128_GET_BITS(h, upper_frac_bits, I + 1);
        // the size argument is what is left of the buffer, not its total size plus the offset
        p += snprintf(p, sizeof(str) - static_cast<size_t>(p - str), "%llu", integer);
        // what is left is the fraction, which is below one and so always a positive value of this type
        fixed_point128 temp(l, h & ~int_mask);
        // check if temp has additional digits (not zero)
        if (temp) {
            *p++ = '.';
        }
        // the faster way, requires temp *= 10 not overflowing
        int digits = 0;
        uint64_t res[2] {};
        while (digits++ < max_frac_digits && temp) {
            if constexpr (I < 4) {
                res[0] = mulx_u64(temp.high, 10ull, &res[1]);  // multiply by 10
                // extract the integer part
                integer = shift_right128_round(res[0], res[1], upper_frac_bits);
                temp *= 10;  // move another digit to the integer area
            } else {
                temp *= 10;  // move another digit to the integer area
                integer = FP128_GET_BITS(temp.high, upper_frac_bits, I + 1);
            }
            *p++ = '0' + (char)integer;
            temp.high &= ~int_mask;
        }
        // round: if we stopped due to max_frac_digits and there's still a remaining
        // fraction, peek at the next digit and round up if >= 5
        if (temp) {
            uint64_t next_digit;
            if constexpr (I < 4) {
                res[0] = mulx_u64(temp.high, 10ull, &res[1]);
                next_digit = shift_right128_round(res[0], res[1], upper_frac_bits);
            } else {
                fixed_point128 t = temp;
                t *= 10;
                next_digit = FP128_GET_BITS(t.high, upper_frac_bits, I + 1);
            }
            if (next_digit >= 5) {
                // carry-propagate backwards through the digit string
                char* q = p - 1;
                bool carry = true;
                while (carry && q >= &str[0]) {
                    if (*q >= '0' && *q <= '9') {
                        if (*q == '9') {
                            *q-- = '0';
                        } else {
                            ++(*q);
                            carry = false;
                        }
                    } else {
                        --q;  // skip '.' or '-'
                    }
                }
                if (carry) {
                    // all digits were 9 (e.g. "9.999" → "10.000")
                    // insert '1' before the integer digits, after any '-'
                    char* start = str + (str[0] == '-' ? 1 : 0);
                    memmove(start + 1, start, p - start + 1);  // +1 for '\0'
                    *start = '1';
                    ++p;
                }
            }
        }
        *p = '\0';
        // remove trailing zeros after the decimal point
        if (char* dot = strchr(str + (str[0] == '-' ? 1 : 0), '.')) {
            char* last = p - 1;
            while (last > dot && *last == '0')
                --last;
            if (last == dot)
                --last;  // remove the dot too if no fraction remains
            *(last + 1) = '\0';
        }
        return str;
    }
    /// @}

    /// @name Arithmetic Operators
    /// @{

    /**
     * @brief Performs right shift operation.
     *
     * Hinted rather than forced, unlike the other one line forwarders (see FP128_FORCE_INLINE).
     * operator>>= is itself forced, so forcing this one too would expand the whole shift, its
     * rounding included, at every use of the by value form. The equivalent pair in float128
     * measurably hurt a loop that shifts alongside other arithmetic; this one is left hinted for
     * the same reason. The callee being forced means the shift is inlined either way once the
     * compiler decides to inline this wrapper.
     *
     * @param shift bits to shift
     * @return Temporary object with the result of the operation
     */
    template <typename T> [[nodiscard]] FP128_INLINE constexpr fixed_point128 operator>>(T shift) const noexcept
    {
        fixed_point128 temp(*this);
        return temp >>= static_cast<int32_t>(shift);
    }
    /**
     * @brief Performs left shift operation.
     * Hinted rather than forced, see operator>> above.
     * @param shift bits to shift
     * @return Temporary object with the result of the operation
     */
    template <typename T> [[nodiscard]] FP128_INLINE constexpr fixed_point128 operator<<(T shift) const noexcept
    {
        fixed_point128 temp(*this);
        return temp <<= static_cast<int32_t>(shift);
    }
    /**
     * @brief Add a value to this object
     * @param rhs Right hand side operand
     * @return This object.
     */
    FP128_FORCE_INLINE constexpr fixed_point128& operator+=(const fixed_point128& rhs) noexcept
    {
        // a two's complement addition needs no attention to the signs at all: the sum of the two
        // 128 bit patterns is the pattern of the sum, and a result outside the range wraps silently
        const uint8_t carry = addcarryx_u64(0, low, rhs.low, &low);
        addcarryx_u64(carry, high, rhs.high, &high);
        return *this;
    }
    /**
     * @brief Add a value to this object
     * @param rhs Right hand side operand
     * @return This object.
     */
    template <typename T> FP128_FORCE_INLINE constexpr fixed_point128& operator+=(const T& rhs) noexcept { return operator+=(fixed_point128(rhs)); }
    /**
     * @brief Subtract a value to this object
     * @param rhs Right hand side operand
     * @return This object.
     */
    FP128_FORCE_INLINE constexpr fixed_point128& operator-=(const fixed_point128& rhs) noexcept
    {
        // subtracted directly rather than by adding the two's complement of rhs, which would cost
        // a negation that the borrow chain gets for free. Both wrap on overflow.
        const uint8_t borrow = subborrow_u64(0, low, rhs.low, &low);
        subborrow_u64(borrow, high, rhs.high, &high);
        return *this;
    }
    /**
     * @brief Subtract a value to this object
     * @param rhs Right hand side operand
     * @return This object.
     */
    template <typename T> FP128_FORCE_INLINE constexpr fixed_point128& operator-=(const T& rhs) noexcept { return operator-=(fixed_point128(rhs)); }
    /**
     * @brief Multiplies a value to this object
     * @param rhs Right hand side operand
     * @return This object.
     */
    FP128_FORCE_INLINE constexpr fixed_point128& operator*=(const fixed_point128& rhs) noexcept
    {
        // The 256 bit intermediate product lives in four named locals rather than an array. Both
        // spellings hand the address of a QWORD to the extended arithmetic intrinsics, but MSVC
        // promotes the scalars into registers more readily than array elements and stops parking
        // a couple of the partial products on the stack.
        // They are left uninitialized to get 10-50% extra performance. Zero initialization is a
        // 10% penalty and using a thread_local static variable lowers performance by >50%.
        uint64_t r0, r1, r2, r3;  // 256 bit of result, low QWORD first

        // The magnitudes are multiplied and the sign is put back at the end. Multiplying the two's
        // complement patterns directly is possible - the signed product is the unsigned one less
        // 2^128 times each negative operand, and the window kept below reaches that far - but it
        // rounds the tie of a negative product towards positive infinity, and correcting that back
        // to away from zero costs more than it saves. Three spellings of the correction were
        // measured against this one and every one of them was slower, by 20% here and by up to 20%
        // again in the functions built on it.
        uint64_t lhs_low = 0, lhs_high = 0, rhs_low = 0, rhs_high = 0;
        const uint64_t result_sign = GetMagnitude(lhs_low, lhs_high) ^ rhs.GetMagnitude(rhs_low, rhs_high);

        // multiply low QWORDs
        r0 = mulx_u64(lhs_low, rhs_low, &r1);

        // multiply high QWORDs (overflow can happen)
        r2 = mulx_u64(lhs_high, rhs_high, &r3);

        // Each cross product is folded in as one unbroken carry chain across r1, r2 and r3.
        // Spelling the last link as a third addcarryx_u64 of zero rather than as `r3 += carry`
        // is what keeps it a chain: MSVC cannot fold an added carry return value back into an
        // ADC and emits SETB, MOVZX and ADD in its place, three instructions and three cycles of
        // latency where the chain wants one. Clang folds either spelling into ADC and is
        // unaffected. Together with the named locals above this took the Mandelbrot escape-time
        // loop from 212 to 195 instructions (25 to 19 stack accesses, 9 to 15 ADCX, 6 to 0 MOVZX)
        // and about 6% off the fp128 render under MSVC, on both P and E cores. The results are
        // bit identical either way.

        // multiply low this and high rhs
        uint64_t cross_high = 0;
        uint64_t cross_low = mulx_u64(lhs_low, rhs_high, &cross_high);
        uint8_t carry = addcarryx_u64(0, r1, cross_low, &r1);
        carry = addcarryx_u64(carry, r2, cross_high, &r2);
        addcarryx_u64(carry, r3, 0, &r3);

        // multiply high this and low rhs
        cross_low = mulx_u64(lhs_high, rhs_low, &cross_high);
        carry = addcarryx_u64(0, r1, cross_low, &r1);
        carry = addcarryx_u64(carry, r2, cross_high, &r2);
        addcarryx_u64(carry, r3, 0, &r3);

        // extract the bits from the product keeping the precision the same as this object
        // shift result by F
        constexpr int32_t index = (F == 64) ? 0 : F / 64;
        constexpr int32_t lsb = (F == 64) ? 64 : (F & FP128_MAX_VALUE_64(6)); // bit within the 64bit QWORD at `index`
        constexpr uint64_t half = 1ull << (lsb - 1);                          // used for rounding

        // The three consecutive QWORDs the result is read from, starting at `index`. I is
        // constrained to [1,63], so F is in [64,126] and `index` is 0 or 1: what was an array
        // subscript is a compile time choice between two triples.
        const uint64_t w0 = (index == 0) ? r0 : r1;
        const uint64_t w1 = (index == 0) ? r1 : r2;
        const uint64_t w2 = (index == 0) ? r2 : r3;
        const bool need_rounding = (w0 & half) != 0;

        // copy block #1 (lowest)
        low  = shift_right128<lsb>(w0, w1);
        high = shift_right128<lsb>(w1, w2);

        FP128_ADD_ROUND_BIT(low, high, need_rounding);
        NegateIf(result_sign);
        return *this;
    }
    /**
     * @brief Squares this object in place.
     *
     * Cheaper than operator*=(*this): a square is symmetric, so the two cross products
     * low*high and high*low are the same value and only one 64x64->128 bit multiply is
     * needed instead of two. The sign is also known in advance - a magnitude times itself is
     * the answer as it stands - which removes the sign combination and the negation.
     *
     * The 256 bit intermediate product is accumulated exactly as in operator*=, so the
     * result is bit identical to (*this) * (*this), including rounding.
     *
     * @return This object.
     */
    FP128_FORCE_INLINE constexpr fixed_point128& square() noexcept
    {
        // Named locals rather than arrays, and one unbroken carry chain per cross product, for
        // the reasons spelled out in operator*=. They are left uninitialized to get extra
        // performance, same as there.
        uint64_t r0, r1, r2, r3;  // 256 bit of result, low QWORD first

        // squared on the magnitude, as operator*= multiplies on it, which is what makes the two
        // agree bit for bit for a negative value as well as a positive one
        uint64_t l = 0, h = 0;
        (void)GetMagnitude(l, h);

        // multiply the low QWORD by itself
        r0 = mulx_u64(l, l, &r1);

        // multiply the high QWORD by itself (overflow can happen)
        r2 = mulx_u64(h, h, &r3);

        // the low * high cross product, which appears twice in the sum
        uint64_t cross_high = 0;
        const uint64_t cross_low = mulx_u64(l, h, &cross_high);

        uint8_t carry = addcarryx_u64(0, r1, cross_low, &r1);
        carry = addcarryx_u64(carry, r2, cross_high, &r2);
        addcarryx_u64(carry, r3, 0, &r3);

        carry = addcarryx_u64(0, r1, cross_low, &r1);
        carry = addcarryx_u64(carry, r2, cross_high, &r2);
        addcarryx_u64(carry, r3, 0, &r3);

        // extract the bits from the product keeping the precision the same as this object
        // shift result by F
        constexpr int32_t index = (F == 64) ? 0 : F / 64;
        constexpr int32_t lsb = (F == 64) ? 64 : (F & FP128_MAX_VALUE_64(6)); // bit within the 64bit QWORD at `index`
        constexpr uint64_t half = 1ull << (lsb - 1);                          // used for rounding

        // the three consecutive QWORDs the result is read from, see the note in operator*=
        const uint64_t w0 = (index == 0) ? r0 : r1;
        const uint64_t w1 = (index == 0) ? r1 : r2;
        const uint64_t w2 = (index == 0) ? r2 : r3;
        const bool need_rounding = (w0 & half) != 0;

        // copy block #1 (lowest)
        low  = shift_right128<lsb>(w0, w1);
        high = shift_right128<lsb>(w1, w2);

        FP128_ADD_ROUND_BIT(low, high, need_rounding);
        return *this;
    }
    /**
     * @brief Multiplies a value to this object
     * @param x Right hand side operand
     * @return This object.
     */
    template <typename T> FP128_FORCE_INLINE constexpr fixed_point128& operator*=(T x) noexcept
    {
        // floating point
        if constexpr (std::is_floating_point_v<T>) {
            return operator*=(fixed_point128(x));
        } else {
            // integers: convert to uint64 for a simpler operation. The magnitude is taken in the
            // unsigned domain, where the negation wraps; negating the signed value is undefined
            // for the most negative one.
            uint64_t magnitude = static_cast<uint64_t>(x);
            bool negate_result = false;
            if constexpr (std::is_signed_v<T>) {
                // always do positive multiplication
                if (x < 0) {
                    magnitude = 0ull - magnitude;
                    negate_result = true;
                }
            }

            operator*=(magnitude);
            if (negate_result)
                Negate();

            return *this;
        }
    }
    /**
     * @brief Multiplies a 64 bit value to this object
     * The product is taken modulo 2^128, which is the signed product of a two's complement value by
     * an unsigned one: multiplying is repeated addition and addition already wraps correctly.
     * @param x Right hand side operand
     * @return This object.
     */
    template <> FP128_FORCE_INLINE constexpr fixed_point128& operator*= <uint64_t>(uint64_t x) noexcept
    {
        uint64_t temp;

        // multiply low QWORDs
        low = mulx_u64(low, x, &temp);
        high = high * x + temp;
        return *this;
    }
    /**
     * @brief Divide this object by x.
     * @param rhs Right hand side operator (denominator)
     * @return this object.
     * @throws std::logic_error when rhs is zero, leaving this object in an unspecified state - the
     *         magnitude of what it held, the sign not yet put back.
     */
    FP128_INLINE fixed_point128& operator/=(const fixed_point128& rhs)
    {
        bool need_rounding = false;
        // trivial case, this object is zero
        if (!*this)
            return *this;

        // The sign of the quotient is settled here and everything below works on magnitudes: the
        // shift, div_64bit and the long division are all unsigned. Taking the magnitude of the
        // numerator in place leaves this object holding it until the sign is put back at the end.
        uint64_t denom_low = 0, denom_high = 0;
        const uint64_t result_sign = TakeMagnitude() ^ rhs.GetMagnitude(denom_low, denom_high);

        // exponent of 2, convert to a much faster shift operation
        if (1 == popcnt128(denom_low, denom_high)) {
            // the leading one of the magnitude sits at bit 127 - lzcnt and stands for 2^(that - F)
            const int32_t expo = I - static_cast<int32_t>(lzcnt128(denom_low, denom_high));
            if (expo > 0)
                shift_right128_inplace_safe(low, high, expo);
            else if (expo < 0)
                shift_left128_inplace_safe(low, high, -expo);
        }
        // optimization for when dividing by an integer. The integer part of a magnitude is at most
        // 2^63, so it always fits in the QWORD div_64bit takes.
        else if (rhs.is_int()) {
            uint64_t q[2] {};
            const uint64_t nom[2] = {low, high};
            const uint64_t denom = denom_high >> upper_frac_bits;
            uint64_t r;
            if (0 == div_64bit((uint64_t*)q, &r, (uint64_t*)nom, denom, 2)) {
                need_rounding = r > (denom >> 1);
                high = q[1];
                low = q[0];
            } else {  // error
                FP128_FLOAT_DIVIDE_BY_ZERO_EXCEPTION;
            }
        } else {
            bool use_long_division = true;
            if constexpr (FP128_USE_RECIPROCAL_FOR_DIVISION != 0) {
                // A numerator magnitude of 2^127 is what lowest() leaves behind. It reads as a
                // negative value of this type, so the multiply cannot be used for it and the long
                // division, which works on the raw QWORDs, takes over. Every other magnitude, and
                // the divisor's - which is not a power of two here, so it is not 2^127 either -
                // is an ordinary positive value.
                use_long_division = high >= sign_mask;
                if (!use_long_division)
                    *this *= reciprocal(fixed_point128(denom_low, denom_high));
            }

            if (use_long_division)
                need_rounding = DivideMagnitude(denom_low, denom_high);
        }

        // Left as a branch, unlike the carry propagating form operator*= and square() use. Both
        // spellings were measured here: this one is 9% faster on the 128 bit divisor benchmark,
        // because a division reaches it once rather than once per result bit and the branch is
        // predictable enough to be free, while the unconditional add sits on the dependency chain.
        if (need_rounding) {
            ++low;
            high += low == 0;
        }
        NegateIf(result_sign);
        return *this;
    }
    /**
     * @brief Divide this object by x.
     * @param x Denominator.
     * @return This object.
     */
    template <typename T> FP128_FORCE_INLINE fixed_point128& operator/=(T x)
    {
        if constexpr (std::is_floating_point_v<T>) {
            return operator/=(static_cast<double>(x));
        } else {
            // integers: convert to uint64 for a simpler operation. The magnitude is taken in the
            // unsigned domain, where the negation wraps; negating the signed value is undefined
            // for the most negative one.
            uint64_t magnitude = static_cast<uint64_t>(x);
            bool negate_result = false;
            if constexpr (std::is_signed_v<T>) {
                // always do positive division
                if (x < 0) {
                    magnitude = 0ull - magnitude;
                    negate_result = true;
                }
            }

            operator/=(magnitude);
            if (negate_result)
                Negate();

            return *this;
        }
    }
    /**
     * @brief Divide this object by x.
     * @param x Denominator.
     * @return This object.
     */
    template <> FP128_FORCE_INLINE fixed_point128& operator/= <double>(double x)
    {
        if (0 == x)
            FP128_FLOAT_DIVIDE_BY_ZERO_EXCEPTION;

        // Simple and common case, the value is an exponent of 2
        // Convert to a much faster shift operation
        const Double d(x);
        if (0 == d.f()) {
            if (d.s())
                Negate();
            const int32_t e = static_cast<int32_t>(d.e()) - 1023;
            return (e >= 0) ? *this >>= e : *this <<= -e;
        }

        // normal division
        return *this /= fixed_point128(x);
    }
    /**
     * @brief Divide this object by x.
     * @param x Denominator.
     * @return This object.
     */
    template <> FP128_INLINE fixed_point128& operator/= <uint64_t>(uint64_t x)
    {
        if (0 == x)
            FP128_FLOAT_DIVIDE_BY_ZERO_EXCEPTION;

        // div_64bit is unsigned, so the magnitude is divided and the sign put back afterwards
        const uint64_t sign_bits = TakeMagnitude();

        const uint64_t nom[2] = {low, high};
        // div_64bit shortcuts a numerator that is smaller than the divisor and returns without
        // writing the quotient at all, so the destination has to start at zero. Otherwise a value
        // whose raw 128 bit form is below the divisor would be left completely unchanged.
        low = high = 0;
        // the results is stored in low and high, the function returns non zero if error (divide by zero or overflow)
        if (0 != div_64bit(&low, nullptr, (uint64_t*)nom, x, 2)) {
            low = high = 0;
        } else {
            NegateIf(sign_bits);
        }

        return *this;
    }
    /**
     * @brief %= operator
     * @param rhs Modulo operand.
     * @return This object.
     */
    FP128_INLINE fixed_point128& operator%=(const fixed_point128& rhs)
    {
        // trivial cases, this object is zero or the rhs is zero
        if (!*this)
            return *this;
        if (!rhs)
            FP128_FLOAT_DIVIDE_BY_ZERO_EXCEPTION;

        // simple case, both are integers (fractions is zero)
        if (is_int() && rhs.is_int()) {
            return operator=(static_cast<int64_t>(*this) % static_cast<int64_t>(rhs));
        }
        // num or denom are fractions
        // x mod y =  x - y * trunc(x/y)
        fixed_point128 x_div_y = *this / rhs;
        // Integer result - remainder is zero.
        // Avoid the extra computation and precision loss with the standard equation.
        if (x_div_y.is_int()) {
            *this = 0;
        }
        // Fraction result - remainder is non zero.
        else {
            *this -= rhs * trunc(x_div_y);
        }

        // the result keeps the sign of the numerator, as fmod() does
        return *this;
    }
    /**
     * @brief Performs modulo with a generic type by converting to fixed_point128.
     * @tparam T Right hand side type
     * @param rhs Divisor
     * @return This object.
     */
    template <typename T> FP128_FORCE_INLINE fixed_point128& operator%=(T rhs) { return operator%=(fixed_point128(rhs)); }
    /**
     * @brief Shift right this object.
     * The shift is arithmetic: the sign bit is replicated, so the result is the value divided by
     * 2^shift and rounded to nearest, with a tie going to the even neighbour.
     * @param shift Bits to shift. Negative or very high values cause undefined behavior.
     * @return This object.
     */
    FP128_FORCE_INLINE constexpr fixed_point128& operator>>=(int32_t shift) noexcept
    {
        shift_right128_inplace_safe_signed(low, high, shift);
        return *this;
    }
    /**
     * @brief Shift left this object.
     * A value shifted past the range of the type wraps silently, sign included.
     * @param shift Bits to shift. Negative or very high values cause undefined behavior.
     * @return This object.
     */
    FP128_FORCE_INLINE constexpr fixed_point128& operator<<=(int32_t shift) noexcept
    {
        shift_left128_inplace_safe(low, high, shift);
        return *this;
    }
    /**
     * @brief Bitwise AND=
     * The operators below work on all 128 bits, the sign bit included.
     * @param rhs AND mask.
     * @return This object.
     */
    FP128_FORCE_INLINE constexpr fixed_point128& operator&=(const fixed_point128& rhs) noexcept
    {
        low &= rhs.low;
        high &= rhs.high;
        return *this;
    }
    /**
     * @brief Bitwise AND=
     * @param rhs Right hand side operand
     * @return This object.
     */
    template <typename T> FP128_FORCE_INLINE constexpr fixed_point128& operator&=(const T& rhs) { return operator&=(fixed_point128(rhs)); }
    /**
     * @brief Bitwise OR=
     * @param rhs OR mask.
     * @return This object.
     */
    FP128_FORCE_INLINE constexpr fixed_point128& operator|=(const fixed_point128& rhs) noexcept
    {
        low |= rhs.low;
        high |= rhs.high;
        return *this;
    }
    /**
     * @brief Bitwise OR=
     * @param rhs Right hand side operand
     * @return This object.
     */
    template <typename T> FP128_FORCE_INLINE constexpr fixed_point128& operator|=(const T& rhs) { return operator|=(fixed_point128(rhs)); }
    /**
     * @brief Bitwise XOR=
     * @param rhs XOR mask.
     * @return This object.
     */
    FP128_FORCE_INLINE constexpr fixed_point128& operator^=(const fixed_point128& rhs) noexcept
    {
        low ^= rhs.low;
        high ^= rhs.high;
        return *this;
    }
    /**
     * @brief Bitwise XOR=
     * @param rhs Right hand side operand
     * @return This object.
     */
    template <typename T> FP128_FORCE_INLINE constexpr fixed_point128& operator^=(const T& rhs) { return operator^=(fixed_point128(rhs)); }
    /**
     * @brief Prefix ++ operation (++a)
     * @return This object.
     */
    FP128_FORCE_INLINE constexpr fixed_point128& operator++() noexcept
    {
        *this += one();
        return *this;
    }
    /**
     * @brief Postfix ++ operation (a++)
     * @return This object.
     */
    FP128_FORCE_INLINE constexpr fixed_point128 operator++(int32_t) noexcept
    {
        fixed_point128 temp(*this);
        ++*this;  // call the prefix implementation
        return temp;
    }
    /**
     * @brief Prefix -- operation (--a)
     * @return This object.
     */
    FP128_FORCE_INLINE constexpr fixed_point128& operator--() noexcept
    {
        *this -= one();
        return *this;
    }
    /**
     * @brief Postfix -- operation (a--)
     * @return This object.
     */
    FP128_FORCE_INLINE constexpr fixed_point128 operator--(int32_t) noexcept
    {
        fixed_point128 temp(*this);
        --*this;  // call the prefix implementation
        return temp;
    }

    //
    // unary operations
    //
    /**
     * @brief Convert to bool
     */
    [[nodiscard]] FP128_FORCE_INLINE constexpr operator bool() const noexcept { return high != 0 || low != 0; }
    /**
     * @brief Logical not (!). Opposite of operator bool.
     */
    [[nodiscard]] FP128_FORCE_INLINE constexpr bool operator!() const noexcept { return high == 0 && low == 0; }
    /**
     * @brief Bitwise not (~). Flips all 128 bits, the sign included, which makes it -x - epsilon().
     */
    [[nodiscard]] FP128_FORCE_INLINE constexpr fixed_point128 operator~() const noexcept
    {
        fixed_point128 temp(*this);
        temp.high = ~high;
        temp.low = ~low;
        return temp;
    }
    /**
     * @brief Unary +. Returns a copy of the object.
     */
    [[nodiscard]] FP128_FORCE_INLINE constexpr fixed_point128 operator+() const noexcept
    {
        fixed_point128 temp(*this);
        return temp;
    }
    /**
     * @brief Unary -. Returns the two's complement of the object.
     * The most negative value has no positive counterpart and is returned unchanged, the same
     * silent wrap the builtin integer types perform.
     */
    [[nodiscard]] FP128_FORCE_INLINE constexpr fixed_point128 operator-() const noexcept
    {
        fixed_point128 temp(*this);
        temp.Negate();
        return temp;
    }

    /// @}

    /// @name Query and Utility Functions
    /// @{

    /**
     * @brief Returns true if the value is an int (fraction is zero)
     * @return True when the fraction is zero.
     */
    [[nodiscard]] FP128_FORCE_INLINE constexpr bool is_int() const noexcept {
        if constexpr (I == 63) {
            // the fraction is exactly the low QWORD, shifting high by 64 would be undefined
            return 0 == low;
        } else {
            // the fraction bits of the upper QWORD are the ones below the integer part and the sign
            return 0 == low && 0 == (high << (I + 1));
        }
    }
    /**
     * @brief Returns true if the value is positive (including zero)
     * @return True when the value is positive
     */
    [[nodiscard]] FP128_FORCE_INLINE constexpr bool is_positive() const noexcept { return 0 == (high & sign_mask); }
    /**
     * @brief Returns true if the value negative (smaller than zero)
     * @return True when the value is negative
     */
    [[nodiscard]] FP128_FORCE_INLINE constexpr bool is_negative() const noexcept { return 0 != (high & sign_mask); }
    /**
     * @brief Returns true if the value is zero
     * @return Returns true if the value is zero
     */
    [[nodiscard]] FP128_FORCE_INLINE constexpr bool is_zero() const noexcept { return 0 == low && 0 == high; }
    /**
     * @brief get a specific bit within the 128 fixed point data
     * @param bit bit to get [0,127]
     * @return 0 or 1. Undefined when bit > 127
     */
    [[nodiscard]] FP128_FORCE_INLINE constexpr int32_t get_bit(uint32_t bit) const noexcept
    {
        if (bit < 64) {
            return FP128_GET_BIT(low, bit);
        }
        return FP128_GET_BIT(high, bit - 64);
    }
    /**
     * @brief Returns the exponent of the object - like the base 2 exponent of a floating point
     * A value of 2.1 would return 1, values in the range [0.5,1.0) would return -1.
     * Reads the magnitude, so the sign of the value makes no difference.
     * @return Exponent of the number
     */
    [[nodiscard]] FP128_FORCE_INLINE constexpr int32_t get_exponent() const noexcept
    {
        const int32_t s = static_cast<int32_t>(lzcnt128(*this));
        return I - s;
    }
    /**
     * @brief Returns an instance of fixed_point128 with the value of pi
     * Needs <B>I</B> to be 2 or more.
     * @return pi
     */
    [[nodiscard]] FP128_INLINE static const fixed_point128& pi() noexcept
    {
        static_assert(I >= 2, "pi() needs template parameter <I> to be at least 2: 3.14159... does not fit in one integer bit!");
        static const fixed_point128 pi = "3.14159265358979323846264338327950288419716939937510";  // 50 first digits of pi
        return pi;
    }
    /**
     * @brief Returns an instance of fixed_point128 with the value of pi * 2
     * Needs <B>I</B> to be 3 or more.
     * @return pi * 2
     */
    [[nodiscard]] FP128_INLINE static const fixed_point128& pi2() noexcept
    {
        static_assert(I >= 3, "pi2() needs template parameter <I> to be at least 3: 6.28318... does not fit in fewer integer bits!");
        static const fixed_point128 pi2 = "6.28318530717958647692528676655900576839433879875021";  // 50 first digits of pi * 2
        return pi2;
    }
    /**
     * @brief Returns an instance of fixed_point128 with the value of pi / 2
     * @return pi / 2
     */
    [[nodiscard]] FP128_INLINE static const fixed_point128& half_pi() noexcept
    {
        static const fixed_point128 half_pi = "1.57079632679489661923132169163975144209858469968755";  // 50 first digits of pi / 2
        return half_pi;
    }
    /**
     * @brief Returns an instance of fixed_point128 with the value of the golden ratio
     * @return golden ratio constant
     */
    [[nodiscard]] FP128_INLINE static const fixed_point128& golden_ratio() noexcept
    {
        static const fixed_point128 golden_ratio = "1.6180339887498948482045868343656381177203";  // 40 first digits of the golden ratio
        return golden_ratio;
    }
    /**
     * @brief Returns an instance of fixed_point128 with the value of e
     * Needs <B>I</B> to be 2 or more.
     * @return e
     */
    [[nodiscard]] FP128_INLINE static const fixed_point128& e() noexcept
    {
        static_assert(I >= 2, "e() needs template parameter <I> to be at least 2: 2.71828... does not fit in one integer bit!");
        static const fixed_point128 e = "2.71828182845904523536028747135266249775724709369";  // 50 first digits of e
        return e;
    }
    /**
     * @brief Returns an instance of fixed_point128 with the value of sqrt(2)
     * @return e
     */
    [[nodiscard]] FP128_INLINE static const fixed_point128& sqrt_2() noexcept
    {
        static const fixed_point128 sqrt_2 = "1.41421356237309504880168872420969807856967187537";  // 50 first digits of sqrt(2)
        return sqrt_2;
    }
    /**
     * @brief Return an instance of fixed_point128 with the value of 1
     * These three exact constants are returned by value rather than by reference like the
     * transcendental ones above: they are representable in every instantiation without a string
     * conversion, so a constexpr function can build them on the spot and callers get a compile
     * time constant instead of a lazily initialized static.
     * @return 1
     */
    [[nodiscard]] FP128_FORCE_INLINE static constexpr fixed_point128 one() noexcept { return fixed_point128(1); }
    /**
     * @brief Return an instance of fixed_point128 with the value of 0.5
     * Shifting one() right is exact for any <B>I</B>. Note that 0.5 lands in the low QWORD when
     * <B>I</B> is 63, which rules out spelling it as a constant of the upper QWORD alone.
     * @return 0.5
     */
    [[nodiscard]] FP128_FORCE_INLINE static constexpr fixed_point128 half() noexcept { return one() >> 1; }
    /**
     * @brief Return an instance of fixed_point128 with the smallest positive value possible
     * @return 1
     */
    [[nodiscard]] FP128_FORCE_INLINE static constexpr fixed_point128 epsilon() noexcept { return fixed_point128(1, 0); }

private:
    /**
     * @brief Number of Mercator series terms log2_fraction() needs for this instantiation.
     *
     * The reduction leaves |z| <= 2^-6, so term n is bounded by 2^(-6n), and the series is cut off
     * once that is below the last bit of the result with eight bits to spare.
     */
    static constexpr int32_t LOG2_TERMS = (F + 8 + log2_reduction_bits - 1) / log2_reduction_bits;

    /**
     * @brief Reads entry n of log2_inv_n_table as a fixed_point128<1>.
     *
     * The table holds 1/(n*ln2) with 127 fraction bits, which is the form float128's log2 wants;
     * fixed_point128<1> carries 126, so each entry is shifted down by a bit as it is read.
     *
     * The shift cannot be hoisted out of the series loop by rescaling the whole recurrence. Doing
     * that is sound - the accumulator is linear in the table entries, so leaving them doubled just
     * doubles the answer - but the accumulator peaks around 1.47 and a doubled one would not fit
     * the range of fixed_point128<1>, which stops just under 2.
     *
     * @param n Index into the table.
     * @return The entry, scaled for the type the series runs in.
     */
    [[nodiscard]] FP128_FORCE_INLINE static fixed_point128<1> Log2InvN(size_t n) noexcept
    {
        const uint64_t high = log2_inv_n_table[n][0], low = log2_inv_n_table[n][1];
        return fixed_point128<1>(shift_right128<1>(low, high), high >> 1);
    }

    /**
     * @brief Fraction part of log2(x) for x in [1,2), by argument reduction and a short series.
     *
     * Three steps:
     * -# Reduction. The leading six fraction bits of x select a tabulated reciprocal, and
     *    multiplying by it brings the value to within 2^-6 of one. Because log2_value_table holds
     *    the logarithm of the reciprocal that is actually stored rather than of the round number it
     *    approximates, the identity log2(x) = log2(x * recip) - log2(recip) is exact and the
     *    reduction introduces no error of its own. The multiply is done on the raw QWORDs as a pure
     *    128 bit fraction, so it neither knows nor cares what I is.
     * -# Series. log2(1+z) = (z - z^2/2 + z^3/3 - ...) / ln(2), evaluated by Horner. With
     *    |z| <= 2^-6 each term buys six more bits, so LOG2_TERMS of them suffice.
     * -# Scaling. The result is converted to this instantiation's scaling exactly once, at the end.
     *
     * The series runs in fixed_point128<1>, which has 126 fraction bits whatever the caller's I is.
     * That is what keeps the accuracy: z is exact to 2^-127 coming out of the reduction, and
     * rounding it onto a grid of F bits before the series would cost 0.72 ulp of the result all on
     * its own. Holding it at 126 bits instead leaves the table entry and the final conversion as
     * the only meaningful error terms. The accumulator stays below 1.02 and 1/ln(2) is 1.443, so
     * both fit the range of fixed_point128<1>, which is just under 2.
     *
     * The caller passes the raw significand rather than a value already brought to [1,2). Dividing
     * down to that range is a right shift, and for an argument whose leading one sits above the
     * binary point it pushes up to I-1 of the significant bits off the bottom - worth about 1.4 ulp
     * of the answer, which was the largest single error term in the previous implementation.
     * Normalizing upwards here instead is a left shift, so the whole significand survives.
     *
     * @param low Low QWORD of the argument's raw significand.
     * @param high High QWORD of the argument's raw significand.
     * @param leading_zeros Count of leading zero bits in that significand.
     * @return The fraction part of log2, in [0,1).
     */
    [[nodiscard]] static FP128_INLINE fixed_point128 log2_fraction(uint64_t low, uint64_t high, int32_t leading_zeros) noexcept
    {
        using work_t = fixed_point128<1>;
        constexpr int32_t K = log2_reduction_bits;

        // Bring the leading one to bit 127, turning the significand into a fraction in [0.5,1).
        uint64_t m_low = low, m_high = high;
        shift_left128_inplace_safe(m_low, m_high, leading_zeros);

        // The leading one sits at bit 127; the six bits below it choose the reciprocal.
        const size_t j = static_cast<size_t>((m_high >> (63 - K)) & ((1ull << K) - 1));

        uint64_t p_low = 0, p_high = 0;
        mul128_high(m_low, m_high, log2_recip_table[j][1], log2_recip_table[j][0], p_low, p_high);

        // z = 2 * (x/2 * recip) - 1, which is the product with its leading bit removed. Subtracting
        // one half from the raw fraction wraps in the unsigned domain into exactly the two's
        // complement of a negative difference, which is what a reciprocal that rounded down
        // produces. The difference stands for z scaled by 2^127, one bit more than the 126 the raw
        // form of a fixed_point128<1> carries, so it is shifted down by one - arithmetically, z
        // having a sign now.
        constexpr uint64_t one_half = 1ull << 63;
        uint64_t z_low = p_low, z_high = p_high - one_half;
        shift_right128_inplace_safe_signed(z_low, z_high, 1);
        const work_t z(z_low, z_high);

        // Horner over 1/(n*ln2), from the last term down, so the result is already base two.
        work_t acc = Log2InvN(LOG2_TERMS - 1);
        for (int32_t n = LOG2_TERMS - 1; n >= 1; --n) {
            acc = Log2InvN(n - 1) - z * acc;
        }
        const work_t series = z * acc;

        // -log2(recip), the part of the answer the reduction removed, converted from a raw 128 bit
        // fraction to this scaling. The shift is by I + 1, which reaches 64 for the widest integer
        // part, so it goes through the variant that handles a shift of a whole QWORD.
        uint64_t table_low = log2_value_table[j][1], table_high = log2_value_table[j][0];
        shift_right128_inplace_safe(table_low, table_high, I + 1);

        return fixed_point128(table_low, table_high) + fixed_point128(series);
    }

    /**
     * @brief Fraction part of log2(x) for x in [1,2), one bit at a time.
     *
     * Squaring the argument repeatedly and recording whether each square left the [1,2) range
     * yields the bits of the answer from the top down: F squarings for F bits. log2_fraction() is
     * about five times faster and more accurate, so this survives only for constant evaluation,
     * where neither the tables nor the raw QWORD reduction of the fast path can be reached.
     *
     * The iteration runs on x/2 in [0.5,1) rather than on x itself, which keeps every square
     * inside [0.25,1). Squaring x directly reaches 4, and a magnitude of 4 needs three integer
     * bits, so the earlier form of this loop returned a completely wrong answer for
     * fixed_point128<1> - it compared against a constant 2 that instantiation cannot hold.
     * Halving costs the lowest bit of the argument, which is why this is the slow path's problem
     * to have and not the fast one's.
     *
     * @param x Value in [1,2).
     * @return log2(x), in [0,1).
     */
    [[nodiscard]] static constexpr fixed_point128 log2_fraction_bitwise(fixed_point128 x) noexcept
    {
        const fixed_point128 half = fixed_point128::one() >> 1;  // 0.5
        fixed_point128 u = x >> 1;                               // in [0.5,1)
        fixed_point128 b = half;
        fixed_point128 fy;  // fraction part of the result
        for (auto i = 0; i < fixed_point128::F; ++i) {
            u.square();  // in [0.25,1)
            // the square of x left [1,2) exactly when the square of x/2 reached 0.5
            if (u >= half) {
                fy |= b;  // ORing is identical (in this case) but faster than addition.
            } else {
                shift_left128_inplace(u.low, u.high, 1);
            }
            // divide base by 2 using inplace shifts
            shift_right128_inplace(b.low, b.high, 1);
        }

        return fy;
    }

    /**
     * @brief Divides the magnitude this object holds by another magnitude, in place.
     *
     * The long division path of operator/=, which two of that function's branches reach. Both
     * operands are unsigned 128 bit magnitudes rather than values of this type, and that is what
     * lets this serve the one numerator the reciprocal path cannot: 2^127, which is the magnitude
     * of lowest() and reads as a negative value when held in a fixed_point128.
     *
     * The numerator is extended with 128 zero bits below it, so the quotient arrives with 128
     * fraction bits and is brought down to the F this type keeps.
     *
     * @param denom_low Low QWORD of the divisor's magnitude
     * @param denom_high High QWORD of the divisor's magnitude
     * @return True when the dropped part of the quotient reached half a bit, which the caller
     *         rounds up for.
     */
    FP128_INLINE bool DivideMagnitude(uint64_t denom_low, uint64_t denom_high)
    {
        uint64_t q[4] {};
        const uint64_t nom[4] = {0, 0, low, high};
        const uint64_t denom[2] = {denom_low, denom_high};

        if (0 != div_32bit((uint32_t*)q, nullptr, (uint32_t*)nom, (uint32_t*)denom, 2ll * array_length(nom), 2ll * array_length(denom))) {
            FP128_FLOAT_DIVIDE_BY_ZERO_EXCEPTION;
        }

        // the quotient carries 128 fraction bits and this type keeps F of them, so it is shifted
        // right by 128 - F, which is I + 1
        constexpr uint64_t half = 1ull << I;  // the highest bit dropped, used for rounding
        const bool need_rounding = (q[0] & half) != 0;
        if constexpr (I == 63) {
            // a shift of exactly one QWORD, which the intrinsic behind FP128_SHIFTRIGHT128 takes
            // modulo 64 and would leave unshifted
            low = q[1];
            high = q[2];
        } else {
            high = FP128_SHIFTRIGHT128(q[1], q[2], I + 1);
            low = FP128_SHIFTRIGHT128(q[0], q[1], I + 1);
        }

        return need_rounding;
    }

    /// @}

    /// @name Binary Math Operators (Friend)
    /// @{

    // Each of these applies the compound assignment to the by value left hand side and then returns
    // it on a line of its own, rather than the shorter `return lhs OP= rhs;`. The two are
    // equivalent, but the short form makes the return statement copy construct from an lvalue, which
    // MSVC compiles into a store forwarding stall. int128_base's binary operators carry the full
    // explanation and the measurement.

    /**
     * @brief Adds 2 values and returns the result.
     * @param lhs left hand side operand
     * @param rhs Right hand side operand
     * @return Result of the operation
     */
    template <typename T>
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr fixed_point128 operator+(fixed_point128 lhs, const T& rhs) noexcept
    {
        lhs += rhs;
        return lhs;
    }
    /**
     * @brief subtracts the right hand side operand to this object to and returns the result.
     * @param lhs left hand side operand
     * @param rhs Right hand side operand
     * @return The fixed_point128 result
     */
    template <typename T>
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr fixed_point128 operator-(fixed_point128 lhs, const T& rhs) noexcept
    {
        lhs -= rhs;
        return lhs;
    }
    /**
     * @brief Multiplies the right hand side operand with this object to and returns the result.
     * @param lhs left hand side operand
     * @param rhs Right hand side operand
     * @return The fixed_point128 result
     */
    template <typename T>
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr fixed_point128 operator*(fixed_point128 lhs, const T& rhs) noexcept
    {
        lhs *= rhs;
        return lhs;
    }
    /**
     * @brief Divides this object by the right hand side operand and returns the result.
     * @param lhs left hand side operand
     * @param rhs Right hand side operand
     * @return The fixed_point128 result
     */
    template <typename T> [[nodiscard]] friend FP128_FORCE_INLINE fixed_point128 operator/(fixed_point128 lhs, const T& rhs)
    {
        lhs /= rhs;
        return lhs;
    }
    /**
     * @brief Calculates modulo.
     * @param lhs left hand side operand
     * @param rhs Right hand side operand
     * @return The fixed_point128 result
     */
    template <typename T> [[nodiscard]] friend FP128_FORCE_INLINE fixed_point128 operator%(fixed_point128 lhs, const T& rhs)
    {
        lhs %= rhs;
        return lhs;
    }

    /// @}

    /// @name Binary Bitwise Operators (Friend)
    /// @{

    /**
     * @brief Performs bitwise AND (&)
     * @param lhs left hand side operand
     * @param rhs Right hand side operand
     * @return The fixed_point128 result
     */
    template <typename T> [[nodiscard]] friend FP128_FORCE_INLINE constexpr fixed_point128 operator&(fixed_point128 lhs, const T& rhs)
    {
        lhs &= rhs;
        return lhs;
    }
    /**
     * @brief Performs bitwise OR (|)
     * @param lhs left hand side operand
     * @param rhs Right hand side operand
     * @return The fixed_point128 result
     */
    template <typename T> [[nodiscard]] friend FP128_FORCE_INLINE constexpr fixed_point128 operator|(fixed_point128 lhs, const T& rhs)
    {
        lhs |= rhs;
        return lhs;
    }
    /**
     * @brief Performs bitwise XOR (^)
     * @param lhs left hand side operand
     * @param rhs Right hand side operand
     * @return The fixed_point128 result
     */
    template <typename T> [[nodiscard]] friend FP128_FORCE_INLINE constexpr fixed_point128 operator^(fixed_point128 lhs, const T& rhs)
    {
        lhs ^= rhs;
        return lhs;
    }

    /// @}

    /// @name Binary Operators With The Scalar On The Left Hand Side (Friend)
    /// @brief Overloads that accept a builtin arithmetic type as the left operand.
    ///
    /// Without these, an expression like (1 + x) is ambiguous: converting the literal to
    /// fixed_point128 and converting x to a builtin type are both one user defined conversion, so
    /// neither overload wins. Restricting the left operand to the arithmetic types keeps these from
    /// competing with the fixed_point128 on the left versions above, which would otherwise be an
    /// equally good match. The comparison operators already carry the same pair of overloads.
    ///
    /// The scalar is widened to fixed_point128 rather than the object being narrowed to a builtin
    /// type, so (1 - x) keeps the fraction of x instead of truncating it away.
    /// @{

    /// @brief Adds a scalar and a fixed_point128, in that order. @param lhs Left operand @param rhs Right operand @return The fixed_point128 result
    template <typename T>
        requires std::is_arithmetic_v<T>
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr fixed_point128 operator+(const T& lhs, const fixed_point128& rhs) noexcept
    {
        return fixed_point128(lhs) += rhs;
    }
    /// @brief Subtracts a fixed_point128 from a scalar. @param lhs Left operand @param rhs Right operand @return The fixed_point128 result
    template <typename T>
        requires std::is_arithmetic_v<T>
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr fixed_point128 operator-(const T& lhs, const fixed_point128& rhs) noexcept
    {
        return fixed_point128(lhs) -= rhs;
    }
    /// @brief Multiplies a scalar and a fixed_point128, in that order. @param lhs Left operand @param rhs Right operand @return The fixed_point128 result
    template <typename T>
        requires std::is_arithmetic_v<T>
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr fixed_point128 operator*(const T& lhs, const fixed_point128& rhs) noexcept
    {
        return fixed_point128(lhs) *= rhs;
    }
    /// @brief Divides a scalar by a fixed_point128. @param lhs Left operand @param rhs Right operand @return The fixed_point128 result
    template <typename T>
        requires std::is_arithmetic_v<T>
    [[nodiscard]] friend FP128_FORCE_INLINE fixed_point128 operator/(const T& lhs, const fixed_point128& rhs)
    {
        return fixed_point128(lhs) /= rhs;
    }
    /// @brief Performs modulo of a scalar by a fixed_point128. @param lhs Left operand @param rhs Right operand @return The fixed_point128 result
    template <typename T>
        requires std::is_arithmetic_v<T>
    [[nodiscard]] friend FP128_FORCE_INLINE fixed_point128 operator%(const T& lhs, const fixed_point128& rhs)
    {
        return fixed_point128(lhs) %= rhs;
    }
    /// @brief Performs bitwise AND of a scalar and a fixed_point128. @param lhs Left operand @param rhs Right operand @return The fixed_point128 result
    template <typename T>
        requires std::is_arithmetic_v<T>
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr fixed_point128 operator&(const T& lhs, const fixed_point128& rhs)
    {
        return fixed_point128(lhs) &= rhs;
    }
    /// @brief Performs bitwise OR of a scalar and a fixed_point128. @param lhs Left operand @param rhs Right operand @return The fixed_point128 result
    template <typename T>
        requires std::is_arithmetic_v<T>
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr fixed_point128 operator|(const T& lhs, const fixed_point128& rhs)
    {
        return fixed_point128(lhs) |= rhs;
    }
    /// @brief Performs bitwise XOR of a scalar and a fixed_point128. @param lhs Left operand @param rhs Right operand @return The fixed_point128 result
    template <typename T>
        requires std::is_arithmetic_v<T>
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr fixed_point128 operator^(const T& lhs, const fixed_point128& rhs)
    {
        return fixed_point128(lhs) ^= rhs;
    }
    /**
     * @brief Shifts a scalar left by a fixed_point128 shift count.
     *
     * The left operand is widened to fixed_point128 first, so the result is 128 bit rather than the
     * builtin type of lhs.
     *
     * The count is rhs converted to int32_t, the same conversion the fixed_point128 on the left
     * overload applies. That conversion is a shift, so the fraction of rhs is truncated toward zero
     * rather than rounded.
     *
     * @param lhs Left operand, the value being shifted
     * @param rhs Right operand, the shift count
     * @return The fixed_point128 result
     */
    template <typename T>
        requires std::is_arithmetic_v<T>
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr fixed_point128 operator<<(const T& lhs, const fixed_point128& rhs) noexcept
    {
        return fixed_point128(lhs) <<= static_cast<int32_t>(rhs);
    }
    /**
     * @brief Shifts a scalar right by a fixed_point128 shift count.
     * The left operand is widened to fixed_point128 first, see operator<< above.
     * @param lhs Left operand, the value being shifted
     * @param rhs Right operand, the shift count
     * @return The fixed_point128 result
     */
    template <typename T>
        requires std::is_arithmetic_v<T>
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr fixed_point128 operator>>(const T& lhs, const fixed_point128& rhs) noexcept
    {
        return fixed_point128(lhs) >>= static_cast<int32_t>(rhs);
    }

    /// @}

    /// @name Comparison Operators
    /// @{

    /**
     * @brief Compare logical/bitwise equal.
     * @param lhs left hand side operand
     * @param rhs Right hand side operand
     * @return True if this and rhs are equal.
     */
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr bool operator==(const fixed_point128& lhs, const fixed_point128& rhs) noexcept
    {
        return lhs.high == rhs.high && lhs.low == rhs.low;
    }
    /// @overload
    template <typename T> [[nodiscard]] friend FP128_FORCE_INLINE constexpr bool operator==(const fixed_point128& lhs, const T& rhs) noexcept
    {
        return lhs == fixed_point128(rhs);
    }
    /// @overload
    template <typename T> [[nodiscard]] friend FP128_FORCE_INLINE constexpr bool operator==(const T& lhs, const fixed_point128& rhs) noexcept
    {
        return rhs == fixed_point128(lhs);
    }
    /**
     * @brief Return true when objects are not equal. Can be used as logical XOR.
     * @param lhs left hand side operand
     * @param rhs Right hand side operand
     * @return True if not equal.
     */
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr bool operator!=(const fixed_point128& lhs, const fixed_point128& rhs) noexcept
    {
        return lhs.high != rhs.high || lhs.low != rhs.low;
    }
    /// @overload
    template <typename T> [[nodiscard]] friend FP128_FORCE_INLINE constexpr bool operator!=(const fixed_point128& lhs, const T& rhs) noexcept
    {
        return lhs != fixed_point128(rhs);
    }
    /// @overload
    template <typename T> [[nodiscard]] friend FP128_FORCE_INLINE constexpr bool operator!=(const T& lhs, const fixed_point128& rhs) noexcept
    {
        return rhs != fixed_point128(lhs);
    }
    /**
     * @brief Return true if this object is small than the rhs
     * @param lhs left hand side operand
     * @param rhs Right hand side operand
     * @return True when this object is smaller.
     */
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr bool operator<(const fixed_point128& lhs, const fixed_point128& rhs) noexcept
    {
        // the upper QWORDs are compared as signed values, which settles the sign and the magnitude
        // above bit 64 in one step; the lower ones are pure magnitude and compare unsigned
        if (lhs.high != rhs.high)
            return static_cast<int64_t>(lhs.high) < static_cast<int64_t>(rhs.high);

        return lhs.low < rhs.low;
    }
    /// @overload
    template <typename T> [[nodiscard]] friend FP128_FORCE_INLINE constexpr bool operator<(const fixed_point128& lhs, const T& rhs) noexcept
    {
        return lhs < fixed_point128(rhs);
    }
    /// @overload
    template <typename T> [[nodiscard]] friend FP128_FORCE_INLINE constexpr bool operator<(const T& lhs, const fixed_point128& rhs) noexcept
    {
        return fixed_point128(lhs) < rhs;
    }
    /**
     * @brief Return true this object is small or equal than the rhs
     * @param lhs left hand side operand
     * @param rhs Right hand side operand
     * @return True when this object is smaller or equal.
     */
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr bool operator<=(const fixed_point128& lhs, const fixed_point128& rhs) noexcept { return !(lhs > rhs); }
    /// @overload
    template <typename T> [[nodiscard]] friend FP128_FORCE_INLINE constexpr bool operator<=(const fixed_point128& lhs, const T& rhs) noexcept
    {
        return !(lhs > fixed_point128(rhs));
    }
    /// @overload
    template <typename T> [[nodiscard]] friend FP128_FORCE_INLINE constexpr bool operator<=(const T& lhs, const fixed_point128& rhs) noexcept
    {
        return !(fixed_point128(lhs) > rhs);
    }
    /**
     * @brief Return true this object is larger than the rhs
     * @param lhs left hand side operand
     * @param rhs Right hand side operand
     * @return True when this object is larger.
     */
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr bool operator>(const fixed_point128& lhs, const fixed_point128& rhs) noexcept
    {
        // see operator< above
        if (lhs.high != rhs.high)
            return static_cast<int64_t>(lhs.high) > static_cast<int64_t>(rhs.high);

        return lhs.low > rhs.low;
    }
    /// @overload
    template <typename T> [[nodiscard]] friend FP128_FORCE_INLINE constexpr bool operator>(const fixed_point128& lhs, const T& rhs) noexcept
    {
        return lhs > fixed_point128(rhs);
    }
    /// @overload
    template <typename T> [[nodiscard]] friend FP128_FORCE_INLINE constexpr bool operator>(const T& lhs, const fixed_point128& rhs) noexcept
    {
        return fixed_point128(lhs) > rhs;
    }
    /**
     * @brief Return true this object is larger or equal than the rhs
     * @param lhs left hand side operand
     * @param rhs Right hand side operand
     * @return True when this objext is larger or equal.
     */
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr bool operator>=(const fixed_point128& lhs, const fixed_point128& rhs) noexcept { return !(lhs < rhs); }
    /// @overload
    template <typename T> [[nodiscard]] friend FP128_FORCE_INLINE constexpr bool operator>=(const fixed_point128& lhs, const T& rhs) noexcept
    {
        return !(lhs < fixed_point128(rhs));
    }
    /// @overload
    template <typename T> [[nodiscard]] friend FP128_FORCE_INLINE constexpr bool operator>=(const T& lhs, const fixed_point128& rhs) noexcept
    {
        return !(fixed_point128(lhs) < rhs);
    }

    /// @}

    /// @name Friend Math Functions
    /// @brief CRT-style math functions implemented as friend functions.
    /// @{

    /**
     * @brief Returns the absolute value
     * @param x Input value
     * @return A copy of x with the sign removed. lowest() has no positive counterpart and comes
     *         back unchanged, the same silent wrap the builtin integer types perform.
     */
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr fixed_point128 fabs(const fixed_point128& x) noexcept
    {
        fixed_point128 temp = x;
        if (temp.is_negative())
            temp.Negate();

        return temp;
    }
    /**
     * @brief Performs the floor() function, similar to libc's floor(), rounds down towards -infinity.
     * @param x Input value
     * @return A fixed_point128 holding the integer value. Overflow is not reported.
     */
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr fixed_point128 floor(const fixed_point128& x) noexcept
    {
        // dropping the fraction bits of a two's complement value moves it towards minus infinity,
        // which is what floor does, so the mask is the whole of it
        return fixed_point128(0, x.high & int_mask);
    }

    /**
     * @brief Performs the ceil() function, similar to libc's ceil(), rounds up towards infinity.
     * @param x Input value
     * @return A fixed_point128 holding the integer value. Overflow is not reported.
     */
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr fixed_point128 ceil(const fixed_point128& x) noexcept
    {
        if (x.is_int())
            return x;

        // the integer above the floor, which for a value with a fraction is the ceiling
        return fixed_point128(0, (x.high & int_mask) + unity);
    }
    /**
     * @brief Rounds towards zero
     * @param x Value to truncate
     * @return Integer value, rounded towards zero.
     */
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr fixed_point128 trunc(const fixed_point128& x) noexcept
    {
        // towards zero is up for a negative value and down for a positive one
        return (x.is_negative()) ? ceil(x) : floor(x);
    }
    /**
     * @brief Rounds towards the nearest integer.
     * The halfway value (0.5) is rounded away from zero.
     * @param x Value to round
     * @return Integer value, rounded towards the nearest integer.
     */
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr fixed_point128 round(const fixed_point128& x) noexcept
    {
        // the magnitude is rounded and the sign put back, which is what rounding away from zero on
        // a tie means
        const fixed_point128 res = floor(fabs(x) + fixed_point128::half());
        return (x.is_negative()) ? -res : res;
    }
    /**
     * @brief Retrieves an integer that represents the base-2 exponent of the specified value.
     * @param x The specified value.
     * @return Integer value, rounded towards the nearest integer.
     */
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr int32_t ilogb(const fixed_point128& x) noexcept { return x.get_exponent(); }
    /**
     * @brief returns the value of x with the sign of y.
     * @param x The value that's returned as the magnitude of the result.
     * @param y The sign of the result.
     * @return The copysign functions return a floating-point value that combines the magnitude of x and the sign of y.
     */
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr fixed_point128 copysign(const fixed_point128& x, const fixed_point128& y) noexcept
    {
        const fixed_point128 magnitude = fabs(x);
        return (y.is_negative()) ? -magnitude : magnitude;
    }
    /**
     * @brief Performs the fmod() function, similar to libc's fmod(), returns the remainder of a division x/root.
     * @param x Numerator
     * @param y Denominator
     * @return A fixed_point128 holding the modulo value.
     */
    [[nodiscard]] friend FP128_FORCE_INLINE fixed_point128 fmod(const fixed_point128& x, const fixed_point128& y) { return x % y; }
    /**
     * @brief Split into integer and fraction parts.
     * Both results carry the sign of the input variable.
     * @param x Input value
     * @param iptr Pointer to fixed_point128 holding the integer part of x.
     * @return The fraction part of x. Undefined when iptr is nullptr.
     */
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr fixed_point128 modf(const fixed_point128& x, fixed_point128* iptr) noexcept
    {
        if (iptr == nullptr)
            return 0;

        *iptr = trunc(x);  // the integer part, carrying the sign of x
        // What is left is the fraction, which carries that sign as well. The subtraction is exact:
        // the two parts have no bits in common.
        return x - *iptr;
    }
    /**
     * @brief Determines the positive difference between the first and second values.
     * @param x First value
     * @param y Second value
     * @return If x > y returns x - y. Otherwise zero.
     */
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr fixed_point128 fdim(const fixed_point128& x, const fixed_point128& y) noexcept
    {
        return (x > y) ? x - y : fixed_point128(0);
    }
    /**
     * @brief Returns the minimum between x and y.
     * @param x First value
     * @param y Second value
     * @return If x < y returns x. Otherwise y.
     */
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr fixed_point128 fmin(const fixed_point128& x, const fixed_point128& y) noexcept { return (x < y) ? x : y; }
    /**
     * @brief Returns the maximum between x and y.
     * @param x First value
     * @param y Second value
     * @return If x > y returns x. Otherwise y.
     */
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr fixed_point128 fmax(const fixed_point128& x, const fixed_point128& y) noexcept { return (x > y) ? x : y; }
    /**
     * @brief Calculates the hypotenuse. i.e. sqrt(x^2 + y^2)
     * @param x First value
     * @param y Second value
     * @return sqrt(x^2 + y^2).
     */
    [[nodiscard]] friend FP128_FORCE_INLINE fixed_point128 hypot(const fixed_point128& x, const fixed_point128& y) noexcept { return sqrt(sqr(x) + sqr(y)); }
    /**
     * @brief Calculates the square of a value. i.e. x^2
     *
     * Faster than x * x and bit identical to it. Prefer this function wherever a value is
     * multiplied by itself, see fixed_point128::square().
     *
     * @param x Value to square
     * @return x^2, which is never negative.
     */
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr fixed_point128 sqr(fixed_point128 x) noexcept { return x.square(); }
    /**
     * @brief Calculates the left zero count of value x, ignoring the sign.
     * Counted on the magnitude, so a negative value gives the same answer as its opposite rather
     * than the zero its two's complement bits would produce.
     * @param x input value.
     * @return lzc (uint32_t) of th result.
     */
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr uint64_t lzcnt128(const fixed_point128& x) noexcept
    {
        uint64_t l = 0, h = 0;
        (void)x.GetMagnitude(l, h);
        return (h != 0) ? lzcnt64(h) : 64 + lzcnt64(l);
    }
    /**
     * @brief Calculates the square root using Newton's method.
     * Based on the book "Math toolkit for real time programming" by Jack W. Crenshaw
     * @param x Value to calculate the root of
     * @param iterations how many iterations to perform (more is more accurate). Sensible values are 0-5.
     * @return Square root of (x), zero when x <= 0.
     */
    [[nodiscard]] friend FP128_INLINE fixed_point128 sqrt(const fixed_point128& x, uint32_t iterations = 3) noexcept
    {
        static const fixed_point128 factor = "0.70710678118654752440084436210484903928483593768847403658833981";  // sqrt(2) / 2
        if (x.is_negative() || !x)
            return 0;

        // normalize the input to the range [0.5, 1)
        int32_t expo = x.get_exponent() + 1;
        fixed_point128 norm_x = (expo > 0) ? x >> expo : x << -expo;

        // use existing HW to provide an excellent first estimate.
        fixed_point128 root = ::sqrt(static_cast<double>(norm_x));

        // iterate several times via Newton's method
        //                  X
        //   Xn+1 = 0.5 * (---- + Xn )
        //                  Xn
        for (auto i = iterations; i != 0; --i) {
            root = (norm_x * reciprocal(root) + root) >> 1;
        }

        if (expo & 1) {
            root *= factor;
            ++expo;
        }

        // Denormalize the result
        if (expo > 0)
            root <<= (expo + 1) / 2;
        else if (expo < 0)
            root >>= -expo / 2;

        return root;
    }
    /**
     * @brief Factorial reciprocal (inverse). Calculates 1 / x!
     * Maximum value of x that may produce non zero values is 34.
     * This value depends on the amount of fraction bits.
     * @param x Input value
     * @param res Result of the function
     */
    friend FP128_INLINE void fact_reciprocal(int x, fixed_point128& res) noexcept
    {
        static const fixed_point128 c[] = {
            "1",                                          // 1 /  0!
            "1",                                          // 1 /  1!
            "0.5",                                        // 1 /  2!
            "0.166666666666666666666666666666666666667",  // 1 /  3!
            "0.041666666666666666666666666666666666667",  // 1 /  4!
            "0.008333333333333333333333333333333333333",  // 1 /  5!
            "0.001388888888888888888888888888888888889",  // 1 /  6!
            "0.000198412698412698412698412698412698413",  // 1 /  7!
            "0.000024801587301587301587301587301587302",  // 1 /  8!
            "0.000002755731922398589065255731922398589",  // 1 /  9!
            "0.000000275573192239858906525573192239859",  // 1 / 10!
            "0.000000025052108385441718775052108385442",  // 1 / 11!
            "0.000000002087675698786809897921009032120",  // 1 / 12!
            "0.000000000160590438368216145993923771702",  // 1 / 13!
            "0.000000000011470745597729724713851697979",  // 1 / 14!
            "0.000000000000764716373181981647590113199",  // 1 / 15!
            "0.000000000000047794773323873852974382075",  // 1 / 16!
            "0.000000000000002811457254345520763198946",  // 1 / 17!
            "0.000000000000000156192069685862264622164",  // 1 / 18!
            "0.000000000000000008220635246624329716956",  // 1 / 19!
            "0.000000000000000000411031762331216485848",  // 1 / 20!
            "0.000000000000000000019572941063391261231",  // 1 / 21!
            "0.000000000000000000000889679139245057329",  // 1 / 22!
            "0.000000000000000000000038681701706306840",  // 1 / 23!
            "0.000000000000000000000001611737571096118",  // 1 / 24!
            "0.000000000000000000000000064469502843845",  // 1 / 25!
            "0.000000000000000000000000002479596263225",  // 1 / 26!
            "0.000000000000000000000000000091836898638",  // 1 / 27!
            "0.000000000000000000000000000003279889237",  // 1 / 28!
            "0.000000000000000000000000000000113099628",  // 1 / 29!
            "0.000000000000000000000000000000003769988",  // 1 / 30!
            "0.000000000000000000000000000000000121613",  // 1 / 31!
            "0.000000000000000000000000000000000003800",  // 1 / 32!
            "0.000000000000000000000000000000000000115",  // 1 / 33!
            "0.000000000000000000000000000000000000003"   // 1 / 34!
        };
        constexpr int series_len = array_length(c);

        if (x >= 0 && x < series_len) {
            res = c[x];
        } else {
            res = 0;
        }
    }
    /**
     * @brief Calculates the reciprocal of a value. y = 1 / x
     * Using newton iterations: Yn+1 = Yn(2 - x * Yn)
     * @param x Input value
     * @return 1 / x. Returns zero on overflow or division by zero
     */
    [[nodiscard]] friend FP128_INLINE fixed_point128 reciprocal(const fixed_point128& x) noexcept
    {
        static const fixed_point128 one = 1, two = 2;
        static const fixed_point128 xy_max = one + (fixed_point128::epsilon() << 2);
        static const fixed_point128 xy_min = one - (fixed_point128::epsilon() << 2);
        constexpr int max_iterations = 3;
        constexpr int debug = false;
        // 1/0 is infinity, which the double constructor saturates to the largest representable
        // value rather than zero. Catch it here so the documented result actually holds.
        if (!x)
            return fixed_point128(0);

        fixed_point128 y = 1.0 / static_cast<double>(x);

        if (!y)
            return y;

        fixed_point128 xy, y_prev;
        // Newton iterations:
        int i = 0;
        for (; i < max_iterations && (y_prev != y) && (xy < xy_min || xy > xy_max); ++i) {
            y_prev = y;
            xy = x * y;
            // y = y * (two - x * y);
            y *= two - xy;
        }

        if constexpr (debug) {
            static int debug_max_iter = 0;
            if (i > debug_max_iter || i == max_iterations) {
                debug_max_iter = i;
                printf("reciprocal took %i iterations for %.10lf\n", i, static_cast<double>(x));
            }
        }

        return y;
    }
    /**
     * @brief Calculate Sine and Cosine using CORDIC using a limited range of [-pi/2, pi/2]
     * @param x
     * @param sin_x Result of the Sine of x.
     * @param cos_x Result of the Cosine of x.
     * @param apply_scale_factor Apply a scale factor when 60 or more iterations are performed. This is required to get accurate results when the number of
     * iterations is high, but it can be disabled for performance reasons when the number of iterations is low.
     */
    friend FP128_INLINE void _sincos_cordic(fixed_point128 x, fixed_point128& sin_x, fixed_point128& cos_x, bool apply_scale_factor) noexcept
    {
        static const fixed_point128 angles[] = {
            "0.7853981633974483096156608458198757210492",  // arctan(2^-0)
            "0.4636476090008061162142562314612144020285",  // arctan(2^-1)
            "0.2449786631268641541720824812112758109141",  // arctan(2^-2)
            "0.1243549945467614350313548491638710255731",  // arctan(2^-3)
            "0.0624188099959573484739791129855051136062",  // arctan(2^-4)
            "0.0312398334302682762537117448924909770324",  // arctan(2^-5)
            "0.0156237286204768308028015212565703189111",  // arctan(2^-6)
            "0.0078123410601011112964633918421992816212",  // arctan(2^-7)
            "0.0039062301319669718276286653114243871403",  // arctan(2^-8)
            "0.0019531225164788186851214826250767139316",  // arctan(2^-9)
            "0.0009765621895593194304034301997172908516",  // arctan(2^-10)
            "0.0004882812111948982754692396256448486661",  // arctan(2^-11)
            "0.0002441406201493617640167229432596599862",  // arctan(2^-12)
            "0.0001220703118936702042390586461179563009",  // arctan(2^-13)
            "0.0000610351561742087750216625691738291537",  // arctan(2^-14)
            "0.0000305175781155260968618259534385360197",  // arctan(2^-15)
            "0.0000152587890613157621072319358126978851",  // arctan(2^-16)
            "0.0000076293945311019702633884823401050905",  // arctan(2^-17)
            "0.0000038146972656064962829230756163729937",  // arctan(2^-18)
            "0.0000019073486328101870353653693059172441",  // arctan(2^-19)
            "0.0000009536743164059608794206706899231123",  // arctan(2^-20)
            "0.0000004768371582030888599275838214492470",  // arctan(2^-21)
            "0.0000002384185791015579824909479772189326",  // arctan(2^-22)
            "0.0000001192092895507806853113684971379221",  // arctan(2^-23)
            "0.0000000596046447753905544139210621417888",  // arctan(2^-24)
            "0.0000000298023223876953036767401327677095",  // arctan(2^-25)
            "0.0000000149011611938476551470925165959632",  // arctan(2^-26)
            "0.0000000074505805969238279871365645744953",  // arctan(2^-27)
            "0.0000000037252902984619140452670705718119",  // arctan(2^-28)
            "0.0000000018626451492309570290958838214764",  // arctan(2^-29)
            "0.0000000009313225746154785153557354776845",  // arctan(2^-30)
            "0.0000000004656612873077392577788419347105",  // arctan(2^-31)
            "0.0000000002328306436538696289020427418388",  // arctan(2^-32)
            "0.0000000001164153218269348144525990927298",  // arctan(2^-33)
            "0.0000000000582076609134674072264967615912",  // arctan(2^-34)
            "0.0000000000291038304567337036132730326989",  // arctan(2^-35)
            "0.0000000000145519152283668518066395978373",  // arctan(2^-36)
            "0.0000000000072759576141834259033201841046",  // arctan(2^-37)
            "0.0000000000036379788070917129516601402005",  // arctan(2^-38)
            "0.0000000000018189894035458564758300761188",  // arctan(2^-39)
            "0.0000000000009094947017729282379150388117",  // arctan(2^-40)
            "0.0000000000004547473508864641189575194999",  // arctan(2^-41)
            "0.0000000000002273736754432320594787597617",  // arctan(2^-42)
            "0.0000000000001136868377216160297393798823",  // arctan(2^-43)
            "0.0000000000000568434188608080148696899413",  // arctan(2^-44)
            "0.0000000000000284217094304040074348449706",  // arctan(2^-45)
            "0.0000000000000142108547152020037174224853",  // arctan(2^-46)
            "0.0000000000000071054273576010018587112426",  // arctan(2^-47)
            "0.0000000000000035527136788005009293556213",  // arctan(2^-48)
            "0.0000000000000017763568394002504646778106",  // arctan(2^-49)
            "0.0000000000000008881784197001252323389053",  // arctan(2^-50)
            "0.0000000000000004440892098500626161694526",  // arctan(2^-51)
            "0.0000000000000002220446049250313080847263",  // arctan(2^-52)
            "0.0000000000000001110223024625156540423631",  // arctan(2^-53)
            "0.0000000000000000555111512312578270211815",  // arctan(2^-54)
            "0.0000000000000000277555756156289135105907",  // arctan(2^-55)
            "0.0000000000000000138777878078144567552953",  // arctan(2^-56)
            "0.0000000000000000069388939039072283776476",  // arctan(2^-57)
            "0.0000000000000000034694469519536141888238",  // arctan(2^-58)
            "0.0000000000000000017347234759768070944119",  // arctan(2^-59)
            "0.0000000000000000008673617379884035472059",  // arctan(2^-60)
            "0.0000000000000000004336808689942017736029",  // arctan(2^-61)
            "0.0000000000000000001084202172485504434007",  // arctan(2^-63)
            //"0.0000000000000000000542101086242752217003", // arctan(2^-64)
            //"0.0000000000000000000271050543121376108501", // arctan(2^-65)
            //"0.0000000000000000000135525271560688054250", // arctan(2^-66)
            //"0.0000000000000000000067762635780344027125", // arctan(2^-67)
            //"0.0000000000000000000033881317890172013562", // arctan(2^-68)
            //"0.0000000000000000000016940658945086006781", // arctan(2^-69)
            //"0.0000000000000000000008470329472543003390", // arctan(2^-70)
            //"0.0000000000000000000004235164736271501695", // arctan(2^-71)
            //"0.0000000000000000000002117582368135750847", // arctan(2^-72)
            //"0.0000000000000000000001058791184067875423", // arctan(2^-73)
            //"0.0000000000000000000000529395592033937711", // arctan(2^-74)
            //"0.0000000000000000000000264697796016968855", // arctan(2^-75)
            //"0.0000000000000000000000132348898008484427", // arctan(2^-76)
            //"0.0000000000000000000000066174449004242213", // arctan(2^-77)
            //"0.0000000000000000000000033087224502121106", // arctan(2^-78)
            //"0.0000000000000000000000016543612251060553", // arctan(2^-79)
            //"0.0000000000000000000000008271806125530276", // arctan(2^-80)
            //"0.0000000000000000000000004135903062765138", // arctan(2^-81)
            //"0.0000000000000000000000002067951531382569", // arctan(2^-82)
            //"0.0000000000000000000000001033975765691284", // arctan(2^-83)
        };
        constexpr uint32_t angles_count = array_length(angles);
        static const fixed_point128 K = "0.607252935008881256169446752504928263";  // scaling factor for 60 or more iterations

        sin_x = 0.0;  // sin x result
        cos_x = 1.0;  // cos x result
        int32_t base = 0;
        auto angle = angles[0];

        for (auto j = 1; j <= angles_count; ++j) {
            if (x.is_positive()) {  // positive sigma
                fixed_point128 s2 = sin_x + (cos_x >> base);
                fixed_point128 c2 = cos_x - (sin_x >> base);
                cos_x = c2;
                sin_x = s2;
                //
                //  Update the remaining angle.
                //
                x -= angle;
            } else {  // negative sigma
                fixed_point128 s2 = sin_x - (cos_x >> base);
                fixed_point128 c2 = cos_x + (sin_x >> base);
                cos_x = c2;
                sin_x = s2;
                //
                //  Update the remaining angle.
                //
                x += angle;
            }

            ++base;
            //
            //  Update the angle from table, or eventually by just dividing by two.
            //
            if (angles_count < j + 1) {
                angle >>= 1;
            } else {
                angle = angles[j];
            }
        }

        // apply scaling factor for 60 iterations.
        if (apply_scale_factor) {
            sin_x *= K;
            cos_x *= K;
        }
    }
    /**
     * @brief Calculate the sine function over a limited range [-0.5pi, 0.5pi]
     * Using the Maclaurin series expansion, the formula is:
     *              x^3   x^5   x^7
     * sin(x) = x - --- + --- - --- + ...
     *               3!    5!    7!
     * @param x value in Radians in the range [-0.5pi, 0.5pi]
     * @return Sine of x
     */
    [[nodiscard]] friend FP128_INLINE fixed_point128 sin1(fixed_point128 x) noexcept
    {
        static_assert(I >= 4, "sin1() needs template parameter <I> to be at least 4: the series it evaluates carries powers of the argument!");
        assert(fabs(x) <= fixed_point128::half_pi());

        // first part of the series is just 'x'
        const fixed_point128 xx = sqr(x);
        fixed_point128 elem_denom, elem_nom = x;

        // compute the rest of the series, starting with: -(x^3 / 3!)
        for (int i = 3, sign = 1;; i += 2, sign = 1 - sign) {
            elem_nom *= xx;
            fact_reciprocal(i, elem_denom);
            fixed_point128 elem = elem_nom * elem_denom;  // next element in the series
            // precision limit has been hit
            if (!elem)
                break;
            x += (sign) ? -elem : elem;
        }

        return x;
    }
    /**
     * @brief Calculate the cosine function over a limited range [-0.5pi, 0.5pi]
     * Since the sin1 function converges faster, call it with the modified angle.
     * @param x value in Radians in the range [-0.5pi, 0.5pi]
     * @return Cosine of x
     */
    [[nodiscard]] friend FP128_INLINE fixed_point128 cos1(const fixed_point128& x) noexcept
    {
        static_assert(I >= 4, "cos1() needs template parameter <I> to be at least 4: it evaluates sin1()!");
        static const fixed_point128& half_pi = fixed_point128::half_pi();
        assert(fabs(x) <= half_pi);
        return (x.is_positive()) ? sin1(half_pi - x) : -sin1(-half_pi - x);
    }
    /**
     * @brief Calculate the Sine function
     * Ultimately uses sin() with a reduced range of [-pi/4, pi/4]
     * @param x value in Radians
     * @return Sine of x
     */
    [[nodiscard]] friend fixed_point128 sin(fixed_point128 x) noexcept
    {
        static_assert(I >= 4, "sin() needs template parameter <I> to be at least 4: it evaluates sin1() and cos1()!");
        static const fixed_point128& half_pi = fixed_point128::half_pi();  // pi / 2
        double round = (x.is_positive()) ? 0.5 : -0.5;

        int64_t n = static_cast<int64_t>((x / half_pi) + round);
        x -= half_pi * n;
        n = n & 3ull;  // n mod 4
        switch (n) {
        case 0:  // [-45-45) degrees
            return sin1(x);
        case 1:  // [45-135) degrees
            return cos1(x);
        case 2:  // [135-225) degrees
            return -sin1(x);
        case 3:  // [225-315) degrees
        default:
            return -cos1(x);
        }
    }
    /**
     * @brief Calculate the inverse sine function
     * Uses Newton's method to converge quickly.
     * @param x value in radians in the range [-1,1]
     * @return Inverse sine of x
     */
    [[nodiscard]] friend fixed_point128 asin(fixed_point128 x) noexcept
    {
        static_assert(I >= 4, "asin() needs template parameter <I> to be at least 4: it evaluates sin() and cos()!");
        static const fixed_point128 eps = fixed_point128::epsilon() << 1;
        constexpr int max_iterations = 6;
        if (x < -1 || x > 1)
            return 0;

        //              sin(Xn) - a
        // Xn+1 = Xn - -------------
        //                cos(Xn)
        // where 'a' is the argument, each iteration will converge on the result if the initial
        //  estimate is close enough.
        //
        // A step is only taken when it moves sin(Xn) closer to the argument. Newton's method needs
        // that guard here because the derivative it divides by vanishes at the ends of the domain:
        // asin(1) starts at pi/2, where cos(Xn) is ~6.1e-17 while the numerator is nothing but the
        // error of sin(). Dividing the one by the other produces a correction of ~7e-5 to an
        // estimate that was already good to the last bits, and the iteration never recovers.
        // The guard is close to free: sin(candidate) is the value the next iteration needs anyway,
        // and the converged case below leaves before evaluating it at all, so the whole function
        // measures ~1% slower than the unguarded loop it replaced.
        const bool negative = x.is_negative();
        x = fabs(x);
        fixed_point128 res = ::asin(static_cast<double>(x));
        fixed_point128 sin_res = sin(res);
        fixed_point128 residual = fabs(sin_res - x);
        for (int i = 0; i < max_iterations; ++i) {
            const fixed_point128 cos_res = cos(res);
            // the derivative vanished. Dividing by it would throw out of this noexcept function,
            // and there is nothing left to refine anyway.
            if (!cos_res)
                break;

            const fixed_point128 e = (sin_res - x) / cos_res;
            const fixed_point128 candidate = res - e;
            // A step of a couple of ulp has converged and cannot be the runaway case, so it is
            // taken without checking it. Leaving now also skips the sin() below, which is what
            // keeps the guard from costing an extra evaluation per call.
            if (fabs(e) <= eps) {
                res = candidate;
                break;
            }

            const fixed_point128 sin_candidate = sin(candidate);
            const fixed_point128 candidate_residual = fabs(sin_candidate - x);
            // the step made the estimate worse, or no better: it is noise, so keep what we have
            if (candidate_residual >= residual)
                break;

            res = candidate;
            sin_res = sin_candidate;
            residual = candidate_residual;
        }

        return (negative) ? -res : res;
    }
    /**
     * @brief Calculate the cosine function
     * Ultimately uses sin1() with a reduced range of [-pi/4, pi/4]
     * Sine's Maclaurin series converges faster than Cosine's.
     * @param x value in Radians
     * @return Cosine of x
     */
    [[nodiscard]] friend fixed_point128 cos(fixed_point128 x) noexcept
    {
        static_assert(I >= 4, "cos() needs template parameter <I> to be at least 4: it evaluates sin1() and cos1()!");
        static const fixed_point128& half_pi = fixed_point128::half_pi();  // pi / 2
        double round = (x.is_positive()) ? 0.5 : -0.5;

        int64_t n = static_cast<int64_t>((x / half_pi) + round);
        x -= half_pi * n;
        n = n & 3ull;  // n mod 4
        switch (n) {
        case 0:  // [-45-45) degrees
            return cos1(x);
        case 1:  // [45-135) degrees
            return -sin1(x);
        case 2:  // [135-225) degrees
            return -cos1(x);
        case 3:  // [225-315) degrees
        default:
            return sin1(x);
        }
    }
    /**
     * @brief Calculate the inverse cosine function
     * Uses Newton's method to converge quickly.
     * @param x value in radians in the range [-1,1]
     * @return Inverse cosine of x
     */
    [[nodiscard]] friend fixed_point128 acos(fixed_point128 x) noexcept
    {
        static_assert(I >= 4, "acos() needs template parameter <I> to be at least 4: it evaluates sin() and cos()!");
        static const fixed_point128 eps = fixed_point128::epsilon() << 1;
        constexpr int max_iterations = 6;
        if (x < -1 || x > 1)
            return 0;
        //              cos(Xn) - a           a - cos(Xn)
        // Xn+1 = Xn - ------------- = Xn -  ------------
        //                -sin(Xn)              sin(Xn)
        // where 'a' is the argument, each iteration will converge on the result if the initial
        //  estimate is close enough.
        // A step is only taken when it moves cos(Xn) closer to the argument, for the reason spelled
        // out at asin(): at the ends of the domain the derivative divided by is all but zero, and
        // an unguarded step there destroys an estimate that was already correct.
        fixed_point128 res = ::acos(static_cast<double>(x));
        fixed_point128 cos_res = cos(res);
        fixed_point128 residual = fabs(cos_res - x);
        for (int i = 0; i < max_iterations; ++i) {
            const fixed_point128 sin_res = sin(res);
            // At the ends of the domain the derivative vanishes: acos(1) starts at res == 0 where
            // sin is exactly zero. Dividing by it would throw out of this noexcept function, and
            // there is nothing left to refine anyway.
            if (!sin_res)
                break;

            const fixed_point128 e = (x - cos_res) / sin_res;
            const fixed_point128 candidate = res - e;
            // as in asin(): a converged step is taken unchecked, which also skips the cos() below
            if (fabs(e) <= eps) {
                res = candidate;
                break;
            }

            const fixed_point128 cos_candidate = cos(candidate);
            const fixed_point128 candidate_residual = fabs(cos_candidate - x);
            // the step made the estimate worse, or no better: it is noise, so keep what we have
            if (candidate_residual >= residual)
                break;

            res = candidate;
            cos_res = cos_candidate;
            residual = candidate_residual;
        }

        return res;
    }
    /**
     * @brief Calculate the tangent function
     * tan(x) = sin(x)/cos(x)
     * @param x value
     * @return Tangent of x. Zero at the poles (odd multiples of pi/2), where the true value is
     *         unbounded and the cosine comes out exactly zero.
     */
    [[nodiscard]] friend FP128_INLINE fixed_point128 tan(fixed_point128 x) noexcept
    {
        static_assert(I >= 4, "tan() needs template parameter <I> to be at least 4: it evaluates sin() and cos()!");
        constexpr bool use_cordic = false;  // CORDIC is currently slower and less accurate
        if constexpr (use_cordic) {
            fixed_point128 sin_x, cos_x;
            _sincos_cordic(x, sin_x, cos_x, false);

            if (cos_x)
                return sin_x / cos_x;
            return 0;
        } else {
            const fixed_point128 cos_x = cos(x);
            // the range reduction lands exactly on zero at the poles. Dividing by it would throw
            // out of this noexcept function, so report zero like the CORDIC branch above does.
            if (!cos_x)
                return 0;
            return sin(x) / cos_x;
        }
    }
    /**
     * @brief Calculate the inverse tangent function
     * @param x value
     * @return Arctangent of x
     */
    [[nodiscard]] friend fixed_point128 atan(fixed_point128 x) noexcept
    {
        static_assert(I >= 4, "atan() needs template parameter <I> to be at least 4: it evaluates sin() and cos()!");
        // constants for segmentation
        static const fixed_point128& half_pi = fixed_point128::half_pi();  // pi / 2
        static const fixed_point128 eps = fixed_point128::epsilon() << 1;
        bool comp = false;
        constexpr int max_iterations = 6;

        // make argument positive, save the sign
        const bool negative = x.is_negative();
        x = fabs(x);

        // limit argument to 0..1
        if (x > 1) {
            comp = true;
            x = reciprocal(x);
        }

        // initial step uses the CRT function.
        fixed_point128 res = ::atan(static_cast<double>(x));
        //
        // Xn+1 =  Xn - cos(Xn) * ( sin(Xn) - a * cos(Xn))
        //
        // where 'a' is the argument, each iteration will converge on the result if the initial
        //  estimate is close enough.
        for (int i = 0; i < max_iterations; ++i) {
            fixed_point128 cos_xn = cos(res);
            fixed_point128 sin_xn = sin(res);
            fixed_point128 e = cos_xn * (sin_xn - x * cos_xn);  // this is the iteration estimated error
            res -= e;
            if (fabs(e) <= eps)
                break;
        }

        // restore complement if needed
        if (comp)
            res = half_pi - res;

        // restore sign if needed
        return (negative) ? -res : res;
    }
    /**
     * @brief Calculate the inverse tangent function of the ratio y / x
     * @param y value
     * @param x value
     * @return Arctangent of y / x in the range [-pi, pi]
     */
    [[nodiscard]] friend fixed_point128 atan2(fixed_point128 y, fixed_point128 x) noexcept
    {
        static_assert(I >= 4, "atan2() needs template parameter <I> to be at least 4: it evaluates atan()!");
        // constants for segmentation
        static const fixed_point128& pi = fixed_point128::pi();
        static const fixed_point128& half_pi = fixed_point128::half_pi();          // pi / 2

        // x == 0
        if (!x) {
            if (!y)
                return fixed_point128(0);

            return (y.is_negative()) ? -half_pi : half_pi;
        }
        // y == 0
        if (!y)
            return (x.is_negative()) ? pi : fixed_point128(0);

        fixed_point128 res;
        // save the signs of x, y
        bool comp = fabs(y) > fabs(x);
        fixed_point128 ratio;

        // calculate the ratio keeping it below 1.0
        ratio = (comp) ? x / y : y / x;
        res = atan(ratio);

        if (comp)
            res = (res.is_negative()) ? -half_pi - res : half_pi - res;

        if (x > 0)
            return res;

        // x < 0
        return (y < 0) ? res - pi : res + pi;
    }
    /**
     * @brief Calculate the hyperbolic sine function
     * Use the exponent function which produces more accurate results than the power series.
     *           e^x - e^(-x)
     * sinh(x) = ------------
     *                2
     * @param x value
     * @return Sine of x
     */
    [[nodiscard]] friend FP128_FORCE_INLINE fixed_point128 sinh(const fixed_point128& x) noexcept
    {
        static_assert(I >= 4, "sinh() needs template parameter <I> to be at least 4!");
        return (exp(x) - exp(-x)) >> 1;
        // the below code while faster, produces lower precision results
        //    if (fabs(x) > 1) {
        //        return (exp(x) - exp(-x)) >> 1;
        //    }
        //    else {
        //        // Using the Maclaurin series expansion, the formula is :
        //        //                x^3   x^5   x^7
        //        //  sinh(x) = x + --- + --- + ---  +...
        //        //                 2!    4!    6!
        //
        //        // first part of the series is just 'x'
        //        const fixed_point128 xx = x * x;
        //        fixed_point128 elem_denom, elem_nom = x;
        //        fixed_point128 res = x;

        //        // compute the rest of the series, starting with: -(x^3 / 2!)
        //        for (int i = 3; ; i += 2) {
        //            elem_nom *= xx;
        //            fact_reciprocal(i - 1, elem_denom);
        //            fixed_point128 elem = elem_nom * elem_denom; // next element in the series
        //            // precision limit has been hit
        //            if (!elem)
        //                break;
        //            res += elem;
        //        }
        //        return res;
        //    }
    }
    /**
     * @brief Calculates the inverse hyperbolic sine
     * For positive x:
     * asinh(x) = log(x + sqrt(x^2 + 1))
     * For negative x, the function returns the result with the sign inverted
     * @param x value
     * @return Inverse hyperbolic sine of x
     */
    [[nodiscard]] friend FP128_INLINE fixed_point128 asinh(const fixed_point128& x) noexcept
    {
        fixed_point128 absx = fabs(x);
        fixed_point128 res = log(absx + sqrt(sqr(absx) + fixed_point128::one()));

        return (x.is_positive()) ? res : -res;
    }
    /**
     * @brief Calculate the hyperbolic cosine function over a limited range [-0.5pi, 0.5pi]
     *           e^x + e^(-x)
     * cosh(x) = ------------
     *                2
     * @param x value in Radians in the range [-0.5pi, 0.5pi]
     * @return Cosine of x
     */
    [[nodiscard]] friend FP128_FORCE_INLINE fixed_point128 cosh(const fixed_point128& x) noexcept
    {
        static_assert(I >= 4, "cosh() needs template parameter <I> to be at least 4!");
        return (exp(x) + exp(-x)) >> 1;

        // Using the Maclaurin series expansion, the formula is:
        //               x^2   x^4   x^6
        // sinh(x) = 1 + --- + --- + --- + ...
        //                2!    4!    6!

        // const fixed_point128 xx = x * x;
        // fixed_point128 elem_denom, elem_nom = fixed_point128::one();
        //// first part of the series is 1
        // fixed_point128 res = fixed_point128::one();

        //// compute the rest of the series, starting with: -(x^2 / 2!)
        // for (int i = 2; ; i += 2) {
        //     elem_nom *= xx;
        //     fact_reciprocal(i, elem_denom);
        //     fixed_point128 elem = elem_nom * elem_denom; // next element in the series
        //     // precision limit has been hit
        //     if (!elem)
        //         break;
        //     res += elem;
        // }

        // return res;
    }
    /**
     * @brief Calculates the inverse hyperbolic cosine
     * For x >= 1:
     * acosh(x) = log(x + sqrt(x^2 - 1))
     * For x < 1, the function return zero
     * @param x value in the range [1, inf]
     * @return Inverse hyperbolic cosine of x
     */
    [[nodiscard]] friend FP128_INLINE fixed_point128 acosh(const fixed_point128& x) noexcept
    {
        if (x < 1)
            return 0;

        fixed_point128 res = log(x + sqrt(sqr(x) - fixed_point128::one()));
        return res;
    }
    /**
     * @brief Calculates the hyperbolic tangent
     *           e^x - e^(-x)
     * tanh(x) = ------------
     *           e^x + e^(-x)
     * @param x value
     * @return hyperbolic tangent of x
     */
    [[nodiscard]] friend FP128_INLINE fixed_point128 tanh(const fixed_point128& x) noexcept
    {
        // e^x + e^-x never goes below 2, so one integer bit cannot hold the denominator for any argument
        static_assert(I >= 2, "tanh() needs template parameter <I> to be at least 2: e^x + e^-x is never below 2!");
        fixed_point128 ex = exp(x);     // e^x
        fixed_point128 exm1 = exp(-x);  // e^(-x)
        //
        //           e^x - e^(-x)
        // tanh(x) = ------------
        //           e^x + e^(-x)
        //
        return (ex - exm1) / (ex + exm1);
    }
    /**
     * @brief Calculates the inverse hyperbolic tangent
     *                       1 + x
     * atanh(x) = 0.5 * log( -----)
     *                       1 - x
     * @param x value in the range (-1, 1)
     * @return Inverse hyperbolic tangent of x
     */
    [[nodiscard]] friend FP128_INLINE fixed_point128 atanh(const fixed_point128& x) noexcept
    {
        auto one = fixed_point128::one();
        if (fabs(x) >= 1)
            return 0;

        return log((one + x) / (one - x)) >> 1;
    }
    /**
     *                                       x
     * @brief Calculates the exponent of x: e
     * Using the Maclaurin series expansion, the formula is:
     *                  1       2       3
     *                 x       x       x
     * exp(x) = 1  +  ---  +  ---  +  --- + ...
     *                 1!      2!      3!
     *
     * The Maclaurin series will quickly overflow as x's power increases rapidly.
     *
     *                     x   ix   fx
     * Using the equality e = e  * e
     * Where ix is the integer part of x and fx is the fraction part.
     * ix is computed via multiplications which won't overflow if the result value can be held.
     * fx is computed via Maclaurin series expansion, but since fx < 1, it won't overflow.
     * @param x A number specifying a power.
     * @return Exponent of x
     */
    [[nodiscard]] friend FP128_INLINE fixed_point128 exp(const fixed_point128& x) noexcept
    {
        static_assert(I >= 2, "exp() needs template parameter <I> to be at least 2: it multiplies by e(), which needs 2 integer bits!");
        static const fixed_point128 e = fixed_point128::e();
        fixed_point128 _ix, exp_ix;  // integer part of x
        fixed_point128 fx = modf(fabs(x), &_ix);
        uint64_t ix = static_cast<uint64_t>(_ix);  // 64 bit is an overkill to hold the exponent
        fixed_point128 res;

        // compute e^ix (integer part of x)
        if (ix > 0) {
            exp_ix = 1;            // result
            fixed_point128 b = e;  // value of e^1
            while (ix > 0) {
                if (ix & 1)
                    exp_ix *= b;
                ix >>= 1;
                b.square();
            }
        } else {
            exp_ix = 1;
        }

        // compute e^fx (fraction part of x)
        if (fx) {
            fixed_point128 exp_fx = fixed_point128::one() + fx;  // first and second elements of the series
            fixed_point128 elem_denom, elem_nom = fx;

            for (int i = 2;; ++i) {
                elem_nom *= fx;
                fact_reciprocal(i, elem_denom);
                fixed_point128 elem = elem_nom * elem_denom;
                if (!elem)
                    break;
                exp_fx += elem;  // next element in the series
            }

            res = exp_ix * exp_fx;
        } else {
            res = exp_ix;
        }

        return (x.is_positive()) ? res : reciprocal(res);
    }
    /**
     * @brief Calculates the exponent of x and reduces 1 from the result: (e^x) - 1
     * @param x A number specifying a power.
     * @return Exponent of x
     */
    [[nodiscard]] friend FP128_FORCE_INLINE fixed_point128 expm1(const fixed_point128& x) noexcept
    {
        static_assert(I >= 2, "expm1() needs template parameter <I> to be at least 2: it evaluates exp()!");
        return exp(x) - fixed_point128::one();
    }
    /**
     * @brief Computes 2 to the power of x
     * @param x Exponent value
     * @return 2^x
     */
    [[nodiscard]] friend FP128_INLINE fixed_point128 exp2(const fixed_point128& x) noexcept
    {
        static_assert(I >= 2, "exp2() needs template parameter <I> to be at least 2: it evaluates exp()!");
        //
        // Based on exponent law: (x^n)^m = x^(m*n)
        // Convert the exponent x (function parameter) to produce an exponent that will work with exp()
        // y = log(2)
        // 2^x = e^(y*x) = exp(y*x)
        //
        static const fixed_point128 lan2 = "0.693147180559945309417232121458176575";
        return exp(x * lan2);
    }
    /**
     * @brief Computes x to the power of y
     * @param x Base value, must be positive
     * @param y Exponent value
     * @return x^y
     */
    [[nodiscard]] friend FP128_INLINE fixed_point128 pow(const fixed_point128& x, const fixed_point128& y) noexcept
    {
        static_assert(I >= 2, "pow() needs template parameter <I> to be at least 2: it evaluates exp()!");
        //
        // Based on exponent law: (x^n)^m = x^(m * n)
        // Convert the exponent y (function parameter) to produce an exponent that will work with exp()
        // z = log(x)
        // pow(x, y) = x^y = e^(y * z) = exp(y * z)
        //
        if (x.is_negative())
            return 0;

        // log() throws on a zero input and this function is noexcept, so zero has to be handled
        // here. Follows the CRT: pow(0, 0) is 1 and pow(0, y) is zero for any other y.
        if (x.is_zero())
            return (y.is_zero()) ? fixed_point128::one() : fixed_point128(0);

        fixed_point128 lan_x = log(x);
        if (!lan_x)
            return lan_x;

        return exp(y * lan_x);
    }
    /**
     * @brief Calculates the Log base 2 of x: y = log2(x)
     * @param x The number to perform log2 on.
     * @return log2(x)
     */
    [[nodiscard]] friend FP128_INLINE constexpr fixed_point128 log2(fixed_point128 x)
    {
        if (x.is_zero() || x.is_negative()) {
            throw std::domain_error("Math domain error! Function accepts positive, non-zero values only.");
        }

        // Calculate the log in 2 steps:
        // - The integer part (iy) is the position of the leading one, which lzcnt gives directly.
        // - The fraction part is where the work is; see log2_fraction() below.
        // The result is the sum of the two. Based on the identity:
        // log(x + y) = log(x) + log(y)
        const int32_t leading_zeros = static_cast<int32_t>(lzcnt128(x.low, x.high));
        const int32_t expo = 127 - leading_zeros - F;
        const fixed_point128 iy = expo;  // integer part of the result

        // x is an exponent of 2: there is nothing below its leading one, so the fraction is zero.
        const uint64_t below_low = x.low - 1;
        const uint64_t below_high = x.high - ((x.low == 0) ? 1ull : 0ull);
        if (((x.low & below_low) | (x.high & below_high)) == 0) {
            return iy;
        }

        if (std::is_constant_evaluated()) {
            // The bitwise path wants the value brought to [1,2) first.
            if (expo < 0) {
                x <<= -expo;
            } else if (expo > 0) {
                x >>= expo;
            }
            return iy + log2_fraction_bitwise(x);
        }

        return iy + log2_fraction(x.low, x.high, leading_zeros);
    }
    /**
     * @brief Calculates the natural Log (base e) of x: log(x)
     * @param x The number to perform log on.
     * @return log(x)
     */
    [[nodiscard]] friend FP128_INLINE fixed_point128 log(fixed_point128 x)
    {
        static const fixed_point128 inv_log2_e = "0.693147180559945309417232121458176575";
        if (x.is_zero() || x.is_negative()) {
            throw std::domain_error("Math domain error! Function accepts positive, non-zero values only.");
        }
        if (x == 1)
            return 0;

        fixed_point128 y = log2(x);
        return y * inv_log2_e;
    }
    /**
     * @brief Calculates the natural Log (base e) of 1 + x: log(1 + x)
     * Not noexcept: it forwards to log(), which rejects a non positive argument. Marking it
     * noexcept would turn that domain_error into a call to std::terminate.
     * @param x The number to perform log on.
     * @return log1p(x)
     * @throws std::domain_error when x is -1 or below.
     */
    [[nodiscard]] friend FP128_FORCE_INLINE fixed_point128 log1p(fixed_point128 x) { return log(fixed_point128::one() + x); }
    /**
     * @brief Calculates Log base 10 of x: log10(x)
     * @param x The number to perform log on.
     * @return log10(x)
     */
    [[nodiscard]] friend FP128_INLINE fixed_point128 log10(fixed_point128 x)
    {
        static const fixed_point128 inv_log2_10 = "0.301029995663981195213738894724493068";
        if (x.is_zero() || x.is_negative()) {
            throw std::domain_error("Math domain error! Function accepts positive, non-zero values only.");
        }

        if (x == 1)
            return 0;

        fixed_point128 y = log2(x);
        return y * inv_log2_10;
    }
    /**
     * @brief Calculates Log base 2 of x as an integer ignoring the sign of x.
     * Similar to: floor(log2(fabs(x)))
     * @param x The number to perform log on.
     * @return logb(x)
     */
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr fixed_point128 logb(fixed_point128 x)
    {
        if (x.is_zero() || x.is_negative()) {
            throw std::domain_error("Math domain error! Function accepts positive, non-zero values only.");
        }

        return x.get_exponent();
    }

    /// @}
};  // class fixed_point128


/**
 * @brief Writes a fixed_point128 to a stream, honouring its width, fill and adjustment.
 * @tparam I Integer bit count of the type
 */
template <int32_t I> inline std::ostream& operator<<(std::ostream& os, const fixed_point128<I>& value)
{
    std::string text = static_cast<std::string>(value);
    if ((os.flags() & std::ios_base::showpos) && !text.empty() && text[0] != '-')
        text.insert(text.begin(), '+');

    const std::streamsize width = os.width();
    os.width(0);
    if (static_cast<std::streamsize>(text.size()) < width) {
        const size_t padding = static_cast<size_t>(width) - text.size();
        const std::ios_base::fmtflags adjust = os.flags() & std::ios_base::adjustfield;
        if (adjust == std::ios_base::left)
            text.append(padding, os.fill());
        else if (adjust == std::ios_base::internal && !text.empty() && (text[0] == '-' || text[0] == '+'))
            text.insert(1, padding, os.fill());
        else
            text.insert(0, padding, os.fill());
    }
    return os << text;
}

/**
 * @brief Reads a fixed_point128 from a stream, accepting what the constructor from a string does.
 * @tparam I Integer bit count of the type
 */
template <int32_t I> inline std::istream& operator>>(std::istream& is, fixed_point128<I>& value)
{
    std::string token;
    if (!(is >> token))
        return is;

    value = fixed_point128<I>(token.c_str());
    return is;
}

}  // namespace fp128


namespace std
{
/**
 * @brief Numeric properties of `fixed_point128<I>`.
 *
 * A fixed point type is exact on its own grid rather than approximate like a floating point one,
 * so is_exact is true and epsilon is the spacing of that grid - the same everywhere, unlike a
 * float's, which is the spacing near one.
 */
template <int32_t I> class numeric_limits<fp128::fixed_point128<I>>
{
    using value_type = fp128::fixed_point128<I>;

public:
    static constexpr bool is_specialized = true;
    static constexpr bool is_signed = true;
    static constexpr bool is_integer = false;
    /// @brief Every representable value is exactly the number it stands for.
    static constexpr bool is_exact = true;
    static constexpr bool has_infinity = false;
    static constexpr bool has_quiet_NaN = false;
    static constexpr bool has_signaling_NaN = false;
    static constexpr bool has_denorm_loss = false;
    static constexpr float_denorm_style has_denorm = denorm_absent;
    static constexpr float_round_style round_style = round_toward_zero;
    static constexpr bool is_iec559 = false;
    static constexpr bool is_bounded = true;
    static constexpr bool is_modulo = false;
    static constexpr bool traps = false;
    static constexpr bool tinyness_before = false;

    /// @brief Value bits. One of the 128 is the sign, so 127 of them carry magnitude.
    static constexpr int digits = 127;
    /// @brief Decimal digits that survive a round trip: floor(127 * log10(2)).
    static constexpr int digits10 = 38;
    static constexpr int max_digits10 = 39;
    static constexpr int radix = 2;
    /// @brief The last place is 2^(I-127), and the largest magnitude is just under 2^I.
    static constexpr int min_exponent = I - 127;
    static constexpr int max_exponent = I;
    static constexpr int min_exponent10 = static_cast<int>((I - 127) * 0.30102999566398119521);
    static constexpr int max_exponent10 = static_cast<int>(I * 0.30102999566398119521);

    /// @brief Smallest positive value, which is one unit in the last place.
    [[nodiscard]] static value_type min() noexcept { return value_type::epsilon(); }
    /// @brief Largest value, every bit below the sign set: 2^I - epsilon().
    [[nodiscard]] static value_type max() noexcept { return value_type(UINT64_MAX, 0x7FFFFFFFFFFFFFFFull); }
    /// @brief Most negative value, exactly -2^I. Two's complement makes the range asymmetric, so
    ///        this one has no positive counterpart and negating it wraps back to itself.
    [[nodiscard]] static value_type lowest() noexcept { return value_type(0, 1ull << 63); }
    /// @brief Spacing of the grid, the same at every magnitude.
    [[nodiscard]] static value_type epsilon() noexcept { return value_type::epsilon(); }
    [[nodiscard]] static value_type round_error() noexcept { return value_type::half(); }
    [[nodiscard]] static value_type infinity() noexcept { return value_type(); }
    [[nodiscard]] static value_type quiet_NaN() noexcept { return value_type(); }
    [[nodiscard]] static value_type signaling_NaN() noexcept { return value_type(); }
    [[nodiscard]] static value_type denorm_min() noexcept { return value_type(); }
};

/// @brief const, volatile and cv qualified fixed_point128 have the same numeric properties.
template <int32_t I> class numeric_limits<const fp128::fixed_point128<I>> : public numeric_limits<fp128::fixed_point128<I>>
{
};
template <int32_t I> class numeric_limits<volatile fp128::fixed_point128<I>> : public numeric_limits<fp128::fixed_point128<I>>
{
};
template <int32_t I> class numeric_limits<const volatile fp128::fixed_point128<I>> : public numeric_limits<fp128::fixed_point128<I>>
{
};

/// @brief Hash support, so a fixed_point128 can be a key in an unordered container.
template <int32_t I> struct hash<fp128::fixed_point128<I>>
{
    [[nodiscard]] size_t operator()(const fp128::fixed_point128<I>& value) const noexcept
    {
        // every value has exactly one representation, so the two QWORDs are all there is to mix
        uint64_t low = 0, high = 0;
        value.get_components(low, high);

        uint64_t state = low + 0x9E3779B97F4A7C15ull;
        state = (state ^ (state >> 30)) * 0xBF58476D1CE4E5B9ull;
        state = (state ^ (state >> 27)) * 0x94D049BB133111EBull;
        state ^= high + 0x9E3779B97F4A7C15ull + (state << 6) + (state >> 2);
        state = (state ^ (state >> 30)) * 0xBF58476D1CE4E5B9ull;
        return static_cast<size_t>(state ^ (state >> 31));
    }
};

/**
 * @brief std::format support for fixed_point128.
 *
 * Accepts fill and alignment, a sign, zero padding and a width. The type produces its own decimal
 * text, which carries every digit the grid distinguishes, so no precision or presentation type is
 * offered.
 */
template <int32_t I> struct formatter<fp128::fixed_point128<I>, char>
{
    constexpr auto parse(basic_format_parse_context<char>& context)
    {
        auto it = context.begin();
        const auto end = context.end();
        if (it == end || *it == '}')
            return it;

        const auto is_align = [](char c) { return c == '<' || c == '>' || c == '^'; };
        if (it + 1 != end && is_align(*(it + 1))) {
            fill = *it;
            align = *(it + 1);
            it += 2;
        } else if (is_align(*it)) {
            align = *it++;
        }

        if (it != end && (*it == '+' || *it == '-' || *it == ' '))
            sign = *it++;

        if (it != end && *it == '0') {
            zero_pad = true;
            ++it;
        }

        while (it != end && *it >= '0' && *it <= '9')
            width = width * 10 + (*it++ - '0');

        if (it != end && *it != '}')
            throw format_error("invalid format specification for fixed_point128");
        return it;
    }

    template <typename FormatContext> auto format(const fp128::fixed_point128<I>& value, FormatContext& context) const
    {
        string text = static_cast<string>(value);
        if (sign != '-' && !text.empty() && text[0] != '-')
            text.insert(text.begin(), sign);

        if (static_cast<int>(text.size()) < width) {
            const size_t padding = static_cast<size_t>(width) - text.size();
            if (zero_pad && align == 0) {
                const size_t offset = (!text.empty() && (text[0] == '-' || text[0] == '+' || text[0] == ' ')) ? 1u : 0u;
                text.insert(offset, padding, '0');
            } else if (align == '<') {
                text.append(padding, fill);
            } else if (align == '^') {
                text.insert(0, padding / 2, fill);
                text.append(padding - padding / 2, fill);
            } else {
                text.insert(0, padding, fill);
            }
        }
        return std::copy(text.begin(), text.end(), context.out());
    }

    char fill = ' ';        ///< Character the padding is made of.
    char align = 0;         ///< One of < > ^, or zero when none was given.
    char sign = '-';        ///< One of + - space.
    bool zero_pad = false;  ///< The 0 flag: pad with zeros after the sign.
    int width = 0;          ///< Minimum field width.
};
}  // namespace std

#endif  // FP128_FIXED_POINT128_T_H
