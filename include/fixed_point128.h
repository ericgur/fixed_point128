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
    The function div_128bit is derived from the book "Hacker's Delight" 2nd Edition
    by Henry S. Warren Jr. It implements Knuth's "Algorithm D", specialized to a 128 bit
    divisor and written in 64 bit limbs.

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

#include <cmath>       // std::sqrt, std::ldexp, std::frexp, FP_ILOGB0
#include <format>      // std::formatter
#include <functional>  // std::hash
#include <limits>
#include <istream>
#include <ostream>
#include "fp128_shared.h"
#include "fp128_decimal.h"

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

/**
 * @brief Divides a 256 bit numerator by a 128 bit divisor: the one call of div_128bit() in this header.
 *
 * The long division, the remainder and the integer square root of fixed_point128 all divide by a
 * divisor that needs both QWORDs, and all of them come through here. The word count is fixed at four
 * by the template argument of div_128bit<4>(), so its loops are unrolled for this shape whether the
 * compiler inlines it or not. That used to depend on inlining: with the count an ordinary argument
 * and callers in int128_shared.h and float128.h sharing one out of line copy, clang++ compiled a
 * generic division and the long division here lost 15%.
 *
 * @param q (output) Quotient, three words, q[0] lowest.
 * @param r (output, optional) Remainder, two words. Can be nullptr.
 * @param u Numerator, four words, u[0] lowest.
 * @param v Divisor, two words; v[1] must not be zero.
 * @return true on success, false when v[1] is zero or a pointer is null.
 */
[[nodiscard]] inline bool Divide256By128(uint64_t* q, uint64_t* r, const uint64_t* u, const uint64_t* v) noexcept
{
    return div_128bit<4>(q, r, u, v);
}

/**
 * @brief 1 / (n(n+1)) for n = 1..40, with 125 fraction bits rounded to nearest, {high, low}.
 *
 * The ratio of consecutive terms of the sine and cosine series: the sine term after x^n / n! is the
 * previous one times x^2 / ((n+1)(n+2)), entries 2, 4, 6, ... of this table, and the cosine terms use
 * entries 1, 3, 5, ... Generated with Python fractions (exact), see fixed_point128::TrigSeries().
 */
inline constexpr uint64_t fp128_trig_series_table[][2] = {
    {0x1000000000000000ull, 0x0000000000000000ull},  // 1 / (1 * 2)
    {0x0555555555555555ull, 0x5555555555555555ull},  // 1 / (2 * 3)
    {0x02AAAAAAAAAAAAAAull, 0xAAAAAAAAAAAAAAABull},  // 1 / (3 * 4)
    {0x0199999999999999ull, 0x999999999999999Aull},  // 1 / (4 * 5)
    {0x0111111111111111ull, 0x1111111111111111ull},  // 1 / (5 * 6)
    {0x00C30C30C30C30C3ull, 0x0C30C30C30C30C31ull},  // 1 / (6 * 7)
    {0x0092492492492492ull, 0x4924924924924925ull},  // 1 / (7 * 8)
    {0x0071C71C71C71C71ull, 0xC71C71C71C71C71Cull},  // 1 / (8 * 9)
    {0x005B05B05B05B05Bull, 0x05B05B05B05B05B0ull},  // 1 / (9 * 10)
    {0x004A7904A7904A79ull, 0x04A7904A7904A790ull},  // 1 / (10 * 11)
    {0x003E0F83E0F83E0Full, 0x83E0F83E0F83E0F8ull},  // 1 / (11 * 12)
    {0x0034834834834834ull, 0x8348348348348348ull},  // 1 / (12 * 13)
    {0x002D02D02D02D02Dull, 0x02D02D02D02D02D0ull},  // 1 / (13 * 14)
    {0x0027027027027027ull, 0x0270270270270270ull},  // 1 / (14 * 15)
    {0x0022222222222222ull, 0x2222222222222222ull},  // 1 / (15 * 16)
    {0x001E1E1E1E1E1E1Eull, 0x1E1E1E1E1E1E1E1Eull},  // 1 / (16 * 17)
    {0x001AC5701AC5701Aull, 0xC5701AC5701AC570ull},  // 1 / (17 * 18)
    {0x0017F405FD017F40ull, 0x5FD017F405FD017Full},  // 1 / (18 * 19)
    {0x00158ED2308158EDull, 0x2308158ED2308159ull},  // 1 / (19 * 20)
    {0x0013813813813813ull, 0x8138138138138138ull},  // 1 / (20 * 21)
    {0x0011BB4A4046ED29ull, 0x011BB4A4046ED290ull},  // 1 / (21 * 22)
    {0x00103091B51F5E1Aull, 0x4EECC652F8EAC041ull},  // 1 / (22 * 23)
    {0x000ED7303B5CC0EDull, 0x7303B5CC0ED7303Bull},  // 1 / (23 * 24)
    {0x000DA740DA740DA7ull, 0x40DA740DA740DA74ull},  // 1 / (24 * 25)
    {0x000C9A633FCD9673ull, 0x00C9A633FCD96730ull},  // 1 / (25 * 26)
    {0x000BAB656100BAB6ull, 0x56100BAB656100BBull},  // 1 / (26 * 27)
    {0x000AD602B580AD60ull, 0x2B580AD602B580ADull},  // 1 / (27 * 28)
    {0x000A16B312EA8FC3ull, 0x77CD8E80A16B312Full},  // 1 / (28 * 29)
    {0x00096A850096A850ull, 0x096A850096A85009ull},  // 1 / (29 * 30)
    {0x0008CF008CF008CFull, 0x008CF008CF008CF0ull},  // 1 / (30 * 31)
    {0x0008421084210842ull, 0x1084210842108421ull},  // 1 / (31 * 32)
    {0x0007C1F07C1F07C1ull, 0xF07C1F07C1F07C1Full},  // 1 / (32 * 33)
    {0x00074D1E92F0074Dull, 0x1E92F0074D1E92F0ull},  // 1 / (33 * 34)
    {0x0006E25006E25006ull, 0xE25006E25006E250ull},  // 1 / (34 * 35)
    {0x0006806806806806ull, 0x8068068068068068ull},  // 1 / (35 * 36)
    {0x000626703D806267ull, 0x03D80626703D8062ull},  // 1 / (36 * 37)
    {0x0005D392C1005D39ull, 0x2C1005D392C1005Dull},  // 1 / (37 * 38)
    {0x0005871513005871ull, 0x5130058715130058ull},  // 1 / (38 * 39)
    {0x0005405405405405ull, 0x4054054054054054ull},  // 1 / (39 * 40)
    {0x0004FEC04FEC04FEull, 0xC04FEC04FEC04FECull},  // 1 / (40 * 41)
};

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
 *     tanh(), which all reach e().</LI>
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
 * <LI><B>Rounding:</B> every result that does not fit the grid is rounded to the nearest grid value,
 *     and an exact tie goes away from zero (the magnitude is rounded half up). This holds for the
 *     multiplication, the division by a fixed_point128 or by an integer, the right shifts, and the
 *     conversion to an instantiation with more integer bits. The conversion from a floating point
 *     value is the exception: it truncates towards zero, as the conversion to an integer does.</LI>
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
 * <LI>Division and modulo: div_64bit rests on the _udiv128 intrinsic. div_128bit does not, and
 *     nothing in it is barred from a constant expression, but the function that calls it is still
 *     bound by the rest of this list.</LI>
 * <LI>Function local statics, which a constexpr function may not declare. The transcendental
 *     constants (pi, e, sqrt_2, ...) and the factorial reciprocal table are parsed from strings
 *     into them. reciprocal is a division. sqrt and hypot seed their integer square root with the
 *     CRT's sqrt().</LI>
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
    // the stream extractor reads with ParseDecimal(), which reports a value out of range
    template <int32_t J> friend std::istream& operator>>(std::istream& is, fixed_point128<J>& value);

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
    /// Decimal places the text form shows: ceil(F * log10(2)), the fewest that tell every pair of
    /// neighbouring values apart, 10^-max_frac_digits being the first power of ten below the grid
    /// spacing 2^-F. F * log10(2) is never an integer, so the floor plus one is that ceiling. It used
    /// to be F / 3.3, one short for most I: 1 - epsilon printed as 1 and max() as 2^I.
    static constexpr int32_t max_frac_digits = static_cast<int32_t>(F * 0.30102999566398119521) + 1;
    /// Fraction digits the string constructor reads. A point halfway between two grid values has F + 1
    /// decimal places, so 130 settle every rounding; the digits after them cannot reach a halfway
    /// point and, ties going away from zero, cannot change the outcome at one either.
    static constexpr int32_t parsed_frac_digits = 130;
    static constexpr uint64_t sign_mask = 1ull << 63;                    ///< Bitmask of the sign bit within the upper QWORD.
    /// @}

    /**
     * @brief Reads a decimal number to the nearest grid value, a tie going away from zero.
     *
     * Grammar: an optional sign, decimal digits, and an optional point followed by more digits; at
     * least one digit in all. The integer part is kept modulo 2^64, which is all that reaches the
     * 128 bit result: shifted into place it is the integer part modulo 2^(I+1). The fraction is
     * kept exactly, as the integer N of its first parsed_frac_digits digits over 10^D, and rounded
     * as q = floor(N * 2^(F+1) / 10^D), the grid value being (q + 1) / 2 rounded down, which
     * rounds half away from zero. Both steps are exact big integer arithmetic.
     *
     * @param first Start of the text, at the sign or the first digit
     * @param last One past the end of the text
     * @param l Receives the low QWORD of the raw value, wrapped modulo 2^128
     * @param h Receives the high QWORD
     * @param end Receives a pointer past the last character read, first when nothing was read
     * @return 0 when there were no digits (the value is zero), 1 when the value is in range, 2 when
     *         it is outside [lowest(), max()] and wrapped.
     */
    [[nodiscard]] static int32_t ParseDecimal(const char* first, const char* last, uint64_t& l, uint64_t& h, const char*& end) noexcept
    {
        l = h = 0;
        end = first;
        const char* p = first;
        bool negative = false;
        if (p < last && (*p == '-' || *p == '+')) {
            negative = (*p == '-');
            ++p;
        }

        // The integer part, modulo 2^64, and whether it reached 2^64 at all.
        uint64_t integer = 0;
        bool integer_overflow = false;
        bool any_digit = false;
        while (p < last && *p >= '0' && *p <= '9') {
            const uint64_t digit = static_cast<uint64_t>(*p - '0');
            if (integer > (UINT64_MAX - digit) / 10) {
                integer_overflow = true;
            }
            integer = integer * 10 + digit;
            any_digit = true;
            ++p;
        }

        // The fraction as an exact integer over 10^digits.
        detail::big_uint fraction;
        int32_t digits = 0;
        if (p < last && *p == '.') {
            ++p;
            while (p < last && *p >= '0' && *p <= '9') {
                if (digits < parsed_frac_digits) {
                    fraction.mul_small(10u);
                    fraction.add_small(static_cast<uint32_t>(*p - '0'));
                    ++digits;
                }
                any_digit = true;
                ++p;
            }
        }
        if (!any_digit)
            return 0;
        end = p;

        // q = floor(fraction * 2^(F+1) / 10^digits), below 2^(F+1), so it fits the 128 bits read out.
        bool remainder = false;
        fraction.shift_left(F + 1);
        fraction.div_pow10(digits, remainder);
        uint64_t q_low = static_cast<uint64_t>(fraction.at(0)) | (static_cast<uint64_t>(fraction.at(1)) << 32);
        uint64_t q_high = static_cast<uint64_t>(fraction.at(2)) | (static_cast<uint64_t>(fraction.at(3)) << 32);
        // (q + 1) / 2: half away from zero, applied to the magnitude
        const uint8_t carry = addcarryx_u64(0, q_low, 1, &q_low);
        addcarryx_u64(carry, q_high, 0, &q_high);
        shift_right128_inplace(q_low, q_high, 1);

        // magnitude = integer * 2^F + fraction, which F >= 64 puts in the upper QWORD for the integer
        uint64_t m_low = q_low, m_high = 0;
        const uint8_t carry2 = addcarryx_u64(0, q_high, integer << upper_frac_bits, &m_high);
        // Out of range: an integer part above 2^I can only wrap; at 2^I and below the 128 bit
        // magnitude is exact and compares against 2^127 - 1, or 2^127 for a negative value.
        const bool integer_too_large = integer_overflow || integer > (1ull << I);
        const uint64_t limit_high = negative ? (1ull << 63) : (UINT64_MAX >> 1);
        const uint64_t limit_low = negative ? 0 : UINT64_MAX;
        const bool too_large = integer_too_large || carry2 != 0 || m_high > limit_high || (m_high == limit_high && m_low > limit_low);
        if (negative) {
            twos_complement128(m_low, m_high);
        }
        l = m_low;
        h = m_high;
        return too_large ? 2 : 1;
    }

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
     * @brief The decimal digit that ten times a fraction carries into the integer part.
     *
     * Ten times a value below one needs four integer bits, which fixed_point128<1> to <3> do not
     * have: their multiply by ten wraps the integer part away. The digit is therefore read from the
     * full 192 bit product of the raw fraction and ten, and truncated - a digit is the floor of the
     * product, never rounded, or 0.25 would print as 0.35 and 0.96 would produce the character after
     * '9'.
     *
     * @param low Low QWORD of the fraction, a non negative value below one.
     * @param high High QWORD of the fraction.
     * @return floor(10 * fraction), in [0, 9].
     */
    [[nodiscard]] FP128_INLINE static constexpr uint64_t DigitTimesTen(uint64_t low, uint64_t high) noexcept
    {
        uint64_t carry_low = 0, carry_high = 0;
        (void)mulx_u64(low, 10ull, &carry_low);
        uint64_t middle = mulx_u64(high, 10ull, &carry_high);
        const uint64_t top = carry_high + addcarryx_u64(0, middle, carry_low, &middle);
        // the product carries F fraction bits; what is above them is the digit
        return shift_right128(middle, top, upper_frac_bits);
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
     * A value too large for the destination wraps silently, as everywhere else in this class. Fraction
     * bits the destination cannot hold are rounded to nearest, a tie going away from zero.
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
            // shift right by I - I2 bits, replicating the sign bit so a negative value stays negative,
            // and round the dropped fraction bits to nearest
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
     *
     * Reads optional leading white space, an optional sign, decimal digits and an optional point
     * followed by more digits, and stops at the first character that is not part of that. The value
     * is the grid value nearest the decimal one, a tie going away from zero - correctly rounded
     * for any number of digits and any I, which is what lets every value survive the trip through
     * its own text form.
     *
     * The digits after the point are read as one exact integer and scaled by 2^F / 10^digits in big
     * integer arithmetic, see ParseDecimal(). They used to be summed as a series of powers of ten held
     * to 126 bits each, which was 22.8 ulp off at I = 1 and 3 ulp off at I = 4.
     *
     * Much slower than the other constructors. A value too large for the type wraps silently, modulo
     * 2^(I+1) in its integer part, which unlike the saturation the double constructor performs can
     * land on either sign; the stream extractor reports it instead.
     * @param x Input string
     */
    fixed_point128(const char* x) noexcept
    {
        low = high = 0;
        if (x == nullptr)
            return;

        // skip leading white space. The cast keeps a negative char from reaching isspace, which
        // only accepts values representable as unsigned char (or EOF).
        while (*x && isspace(static_cast<unsigned char>(*x)))
            ++x;

        const char* end = nullptr;
        (void)ParseDecimal(x, x + strlen(x), low, high, end);
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
     * Behaves as the cross-template constructor: a value too large for this instantiation wraps, and
     * fraction bits it cannot hold are rounded to nearest, a tie going away from zero.
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
            // shift right by I - I2 bits, replicating the sign bit and rounding to nearest
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
    [[nodiscard]] FP128_INLINE operator std::string() const noexcept { return operator char*(); }
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
        while (digits++ < max_frac_digits && temp) {
            if constexpr (I < 4) {
                integer = DigitTimesTen(temp.low, temp.high);
                temp *= 10;  // move another digit to the integer area, the integer part wraps
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
                next_digit = DigitTimesTen(temp.low, temp.high);
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
     * The quotient is rounded to nearest, a tie going away from zero, whichever of the three paths
     * below computes it.
     * @param rhs Right hand side operator (denominator)
     * @return this object.
     * @throws std::logic_error when rhs is zero - zero divided by zero included - leaving this object
     *         in an unspecified state: the magnitude of what it held, the sign not yet put back.
     */
    FP128_INLINE fixed_point128& operator/=(const fixed_point128& rhs)
    {
        bool need_rounding = false;
        // A zero numerator needs no shortcut: every path below produces zero for it. Taking one here
        // used to return zero for 0 / 0 too, where the division by zero has to throw.

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
            if (div_64bit((uint64_t*)q, &r, (uint64_t*)nom, denom, 2)) {
                // round half up on the magnitude, the rule the long division and the multiply use:
                // 2r >= denom, spelled so that it cannot overflow
                need_rounding = r >= denom - r;
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
     * The quotient is rounded to nearest, a tie going away from zero: the same result as dividing by
     * fixed_point128(x), which the integer divisor only makes cheaper.
     * @param x Denominator.
     * @return This object.
     * @throws std::logic_error when x is zero.
     */
    template <> FP128_INLINE fixed_point128& operator/= <uint64_t>(uint64_t x)
    {
        if (0 == x)
            FP128_FLOAT_DIVIDE_BY_ZERO_EXCEPTION;

        // div_64bit is unsigned, so the magnitude is divided and the sign put back afterwards
        const uint64_t sign_bits = TakeMagnitude();

        const uint64_t nom[2] = {low, high};
        // div_64bit shortcuts a numerator that is smaller than the divisor and returns without
        // writing the quotient at all, so the quotient has to start at zero. Otherwise a value whose
        // raw 128 bit form is below the divisor would be left completely unchanged.
        uint64_t q[2] {};
        uint64_t r = 0;
        // x is not zero and the pointers are valid, which are the only ways div_64bit can fail
        (void)div_64bit(q, &r, nom, x, 2);

        // round half up on the magnitude: 2r >= x, spelled so that it cannot overflow
        low = q[0];
        high = q[1];
        FP128_ADD_ROUND_BIT(low, high, r >= x - r);
        NegateIf(sign_bits);
        return *this;
    }
    /**
     * @brief %= operator, the remainder of a truncating division: x - y * trunc(x / y).
     *
     * The result is exact and keeps the sign of the numerator, as fmod() does. Two values on the same
     * grid of 2^-F are integers scaled by the same power of two, so their remainder is the remainder
     * of the raw 128 bit integers, scaled the same way: no quotient is formed, so it can neither
     * overflow nor be rounded onto the wrong integer.
     *
     * @param rhs Modulo operand.
     * @return This object.
     * @throws std::logic_error when rhs is zero, zero modulo zero included.
     */
    FP128_INLINE fixed_point128& operator%=(const fixed_point128& rhs)
    {
        if (!rhs)
            FP128_FLOAT_DIVIDE_BY_ZERO_EXCEPTION;

        // magnitudes, the remainder takes the sign of this object. lowest() has the magnitude 2^127,
        // which an unsigned QWORD pair holds as it is.
        uint64_t denom_low = 0, denom_high = 0;
        (void)rhs.GetMagnitude(denom_low, denom_high);
        const uint64_t sign_bits = TakeMagnitude();

        uint64_t rem[2] {};
        if (denom_high != 0) {
            // a divisor that needs both QWORDs; the numerator is widened to the four words
            // Divide256By128() takes, which costs two passes over leading zero words
            const uint64_t nom[4] = {low, high, 0, 0};
            const uint64_t denom[2] = {denom_low, denom_high};
            uint64_t q[3] {};
            (void)Divide256By128(q, rem, nom, denom);
        } else {
            const uint64_t nom[2] = {low, high};
            uint64_t q[2] {};
            (void)div_64bit(q, rem, nom, denom_low, 2);
        }

        low = rem[0];
        high = rem[1];
        NegateIf(sign_bits);
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
     * 2^shift and rounded to nearest, a tie going away from zero - the rounding rule of the multiply
     * and the divisions.
     *
     * That rule makes +-epsilon() a fixed point of a shift by one: epsilon() >> 1 is exactly half an
     * epsilon, a tie, and rounds back up to epsilon(). A loop that halves a value until it reaches
     * zero has to bound its iterations or test for epsilon() instead.
     * @param shift Bits to shift. A negative count causes undefined behavior; 128 or more leaves
     *              zero, apart from lowest() >> 128, which is the tie -epsilon() / 2 and so -epsilon().
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

    /// @brief ln(2) as a 128 bit fraction, {high, low}: round(ln(2) * 2^128). The decimal string it
    ///        replaces had wrong 35th and 36th digits.
    static constexpr uint64_t ln2_fraction[2] = {0xB17217F7D1CF79ABull, 0xC9E3B39803F2F6AFull};
    /// @brief log10(2) as a 128 bit fraction, {high, low}: round(log10(2) * 2^128). The decimal
    ///        string it replaces had lost a digit, an error of 4e-35.
    static constexpr uint64_t log10_2_fraction[2] = {0x4D104D427DE7FBCCull, 0x47C4ACD605BE48BCull};

    /**
     * @brief Multiplies this object by a constant in [0, 1) held as a 128 bit fraction.
     *
     * The full 256 bit product of the magnitude and the constant is formed and its upper half kept,
     * rounded to nearest with a tie going away from zero. The constant therefore contributes 2^-128
     * of relative error instead of the 2^-F a fixed_point128 of this instantiation would carry, which
     * matters once the other factor is large: log() multiplies a log2 of up to 63 by ln(2).
     *
     * @param c_low Low QWORD of the constant, scaled by 2^128.
     * @param c_high High QWORD of the constant.
     * @return This object times the constant.
     */
    [[nodiscard]] FP128_INLINE constexpr fixed_point128 MulByFraction128(uint64_t c_low, uint64_t c_high) const noexcept
    {
        uint64_t l = 0, h = 0;
        const uint64_t sign_bits = GetMagnitude(l, h);

        uint64_t r1 = 0, r2 = 0, r3 = 0;
        (void)mulx_u64(l, c_low, &r1);  // only its carry into r1 matters
        r2 = mulx_u64(h, c_high, &r3);
        uint64_t cross_high = 0;
        uint64_t cross_low = mulx_u64(l, c_high, &cross_high);
        uint8_t carry = addcarryx_u64(0, r1, cross_low, &r1);
        carry = addcarryx_u64(carry, r2, cross_high, &r2);
        addcarryx_u64(carry, r3, 0, &r3);
        cross_low = mulx_u64(h, c_low, &cross_high);
        carry = addcarryx_u64(0, r1, cross_low, &r1);
        carry = addcarryx_u64(carry, r2, cross_high, &r2);
        addcarryx_u64(carry, r3, 0, &r3);

        // the product is scaled by 2^128 too much: keep r3:r2, rounded by the highest bit of r1
        fixed_point128 res(r2, r3);
        FP128_ADD_ROUND_BIT(res.low, res.high, (r1 >> 63) != 0);
        res.NegateIf(sign_bits);
        return res;
    }

    /// @brief 1/sqrt(2) as a 128 bit fraction, {high, low}: round(2^-0.5 * 2^128).
    static constexpr uint64_t inv_sqrt2_fraction[2] = {0xB504F333F9DE6484ull, 0x597D89B3754ABE9Full};

    /**
     * @brief sqrt((1 - x) / 2) for x in [1/2, 1], the argument asin() and acos() reduce to.
     * Formed as sqrt(1 - x) / sqrt(2): 1 - x is exact and the root correctly rounded, where halving
     * 1 - x first would drop its last bit and the root would magnify that near x = 1.
     * @param x Value in [1/2, 1].
     * @return sqrt((1 - x) / 2), in [0, 1/2].
     */
    [[nodiscard]] FP128_INLINE static fixed_point128 RootOfHalfComplement(const fixed_point128& x) noexcept
    {
        return sqrt(fixed_point128::one() - x).MulByFraction128(inv_sqrt2_fraction[1], inv_sqrt2_fraction[0]);
    }

    /**
     * @brief asin(x) for x in [0, 1/2], by Newton's method on sin().
     *
     *              sin(Xn) - a
     * Xn+1 = Xn - -------------
     *                cos(Xn)
     * where 'a' is the argument, starting from the double precision asin(). The cosine stays above
     * 0.86 on this range, so the iteration converges quadratically. A step is only taken when it
     * moves sin(Xn) closer to the argument: a correction below the precision of sin() is noise, and
     * the guard keeps the iteration from wandering on it. The converged case leaves before evaluating
     * sin(candidate), so the guard costs close to nothing.
     *
     * @param x Value in [0, 1/2].
     * @return asin(x).
     */
    [[nodiscard]] static fixed_point128 AsinNewton(const fixed_point128& x) noexcept
    {
        constexpr fixed_point128 eps = fixed_point128::epsilon() << 1;
        constexpr int max_iterations = 6;
        fixed_point128 res = ::asin(static_cast<double>(x));
        fixed_point128 sin_res = sin(res);
        fixed_point128 residual = fabs(sin_res - x);
        for (int i = 0; i < max_iterations; ++i) {
            const fixed_point128 e = (sin_res - x) / cos(res);
            const fixed_point128 candidate = res - e;
            // a step of a couple of ulp has converged; taking it unchecked skips the sin() below
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
        return res;
    }

    /// @brief e^-1 as a 128 bit fraction, {high, low}: round(e^-1 * 2^128).
    static constexpr uint64_t inv_e_fraction[2] = {0x5E2D58D8B3BCDF1Aull, 0xBADEC7829054F90Eull};

    /**
     * @brief Converts a constant in [0, 1) held as a 128 bit fraction to this instantiation.
     * @param high High QWORD of the fraction, scaled by 2^128.
     * @param low Low QWORD of the fraction.
     * @return The constant, rounded to nearest.
     */
    [[nodiscard]] FP128_INLINE static constexpr fixed_point128 FromFraction128(uint64_t high, uint64_t low) noexcept
    {
        // 128 fraction bits down to F is a right shift by I + 1, rounded
        shift_right128_inplace_safe(low, high, I + 1);
        return fixed_point128(low, high);
    }

    /**
     * @brief sqrt(x^2 + 1), or sqrt(x^2 - 1) for |x| >= 1, rounded to nearest.
     *
     * The square is formed as an exact 256 bit integer and the root taken with RoundedSqrt256(), so
     * neither function forms x^2 in this type. That matters twice over: x^2 overflows once |x|
     * reaches 2^(I/2), and x^2 - 1 for an x just above one is a tiny value whose rounding error the
     * root magnifies - acosh() was out by 5e8 ulp at fixed_point128<1> and by 1e4 ulp at I = 63.
     *
     * @param x The value to square.
     * @param minus True for x^2 - 1, which needs |x| >= 1; false for x^2 + 1.
     * @return The root. A result of 2^I or more wraps like any overflow.
     */
    [[nodiscard]] FP128_INLINE static fixed_point128 RootOfSquareAndOne(const fixed_point128& x, bool minus) noexcept
    {
        uint64_t l = 0, h = 0;
        (void)x.GetMagnitude(l, h);
        uint64_t s[4] {};
        Square256(l, h, s);

        // one is 2^F in raw units, so its square is the single bit 2F of the 256 bit square
        uint64_t one_square[4] {};
        one_square[(2 * F) / 64] = 1ull << ((2 * F) % 64);
        uint8_t carry = 0;
        for (int32_t i = 0; i < 4; ++i) {
            carry = (minus) ? subborrow_u64(carry, s[i], one_square[i], &s[i]) : addcarryx_u64(carry, s[i], one_square[i], &s[i]);
        }

        fixed_point128 res;
        RoundedSqrt256(s, res.low, res.high);
        return res;
    }

    /**
     * @brief log(x + r) for the magnitude x and the root r that asinh() and acosh() add to it.
     *
     * The sum reaches 2^I once x reaches 2^(I-1). It used to wrap there, into a negative value that
     * made log() throw out of the two noexcept functions and terminate the process. A sum beyond the
     * range is now formed in 129 bits, divided by four, and log(4) added back to its logarithm.
     *
     * @param x_low Low QWORD of the magnitude of x, read as unsigned: 2^127 for lowest().
     * @param x_high High QWORD of the magnitude of x.
     * @param r The root, read as an unsigned raw value too: for lowest() it is a little above 2^127.
     * @return log(x + r).
     */
    [[nodiscard]] FP128_INLINE static fixed_point128 LogOfSum(uint64_t x_low, uint64_t x_high, const fixed_point128& r)
    {
        uint64_t s_low = 0, s_high = 0;
        const uint8_t carry = addcarryx_u64(addcarryx_u64(0, x_low, r.low, &s_low), x_high, r.high, &s_high);
        if (carry == 0 && (s_high >> 63) == 0) {
            return log(fixed_point128(s_low, s_high));
        }

        // a quarter of the 129 bit sum, rounded half up, is below 2^126 and fits every instantiation
        const uint64_t round_bit = (s_low >> 1) & 1;
        s_low = (s_low >> 2) | (s_high << 62);
        s_high = (s_high >> 2) | (static_cast<uint64_t>(carry) << 62);
        s_high += addcarryx_u64(0, s_low, round_bit, &s_low);
        const fixed_point128 ln2 = FromFraction128(ln2_fraction[0], ln2_fraction[1]);
        return log(fixed_point128(s_low, s_high)) + ln2 + ln2;
    }

    /**
     * @brief Maclaurin series of the sine or the cosine of an argument in [0, pi/2].
     *
     * Each term is the previous one times x^2 / (n(n+1)), the ratio read from fp128_trig_series_table,
     * so no power of x and no factorial is ever formed on its own. The whole series runs in
     * fixed_point128<2>, which has 125 fraction bits whatever I is: x^2 reaches 2.47, beyond the range
     * of fixed_point128<1>, while every term and partial sum stays at or below one in magnitude for the
     * cosine and below pi/2 for the sine. The caller rounds the sum to its own grid once, at the end.
     *
     * The earlier evaluation raised x to the n-th power and multiplied by 1/n! in the caller's own
     * instantiation: x^25 is 8e4 for x = pi/2 and wrapped in fixed_point128<16>, and 1/n! held to F bits
     * lost most of its digits in a wide instantiation. cos(0) came out 7026 ulp short of one at I = 63
     * and every cosine was 8e-21 off at I = 16.
     *
     * The series stops once a term drops below an eighth of the caller's last place - the rest of it
     * adds up to less than that - so a narrow instantiation evaluates more terms than a wide one.
     *
     * @param x Argument in [0, pi/2].
     * @param cosine True for the cosine series, false for the sine series.
     * @return The sum of the series, in fixed_point128<2>.
     */
    [[nodiscard]] static FP128_INLINE fixed_point128<2> TrigSeries(const fixed_point128<2>& x, bool cosine) noexcept
    {
        using work_t = fixed_point128<2>;
        // this instantiation's last place is 2^(I-2) units of the 2^-125 grid; an eighth of it
        // (one unit at I = 4, where an eighth would be below the grid)
        constexpr uint64_t cutoff = 1ull << ((I >= 5) ? (I - 5) : 0);

        const work_t xx = sqr(x);
        work_t term = (cosine) ? work_t::one() : x;
        work_t sum = term;
        bool subtract = true;
        // the cosine multiplies by 1/(1*2), 1/(3*4), ... and the sine by 1/(2*3), 1/(4*5), ...
        for (size_t n = (cosine) ? 0 : 1; n < array_length(fp128_trig_series_table); n += 2) {
            term *= xx;
            term *= work_t(fp128_trig_series_table[n][1], fp128_trig_series_table[n][0]);
            uint64_t term_low = 0, term_high = 0;
            term.get_components(term_low, term_high);
            if (term_high == 0 && term_low < cutoff) {
                break;
            }
            if (subtract) {
                sum -= term;
            } else {
                sum += term;
            }
            subtract = !subtract;
        }
        return sum;
    }

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

        // A divisor small enough to fit in one QWORD is div_64bit's job; div_128bit rejects it
        // rather than handle a case that has a cheaper route. Both fill the same first words of q.
        const bool ok = (denom_high != 0) ? Divide256By128(q, nullptr, nom, denom)
                                          : div_64bit(q, nullptr, nom, denom_low, array_length(nom));
        if (!ok) {
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

    /**
     * @brief Squares an unsigned 128 bit integer into its full 256 bit product.
     * @param l Low QWORD of the value
     * @param h High QWORD of the value
     * @param r Receives the product, r[0] lowest.
     */
    FP128_INLINE static constexpr void Square256(uint64_t l, uint64_t h, uint64_t (&r)[4]) noexcept
    {
        r[0] = mulx_u64(l, l, &r[1]);
        r[2] = mulx_u64(h, h, &r[3]);

        // the cross product l * h appears twice: it is added doubled, a 129 bit value, at bit 64
        uint64_t cross_high = 0;
        const uint64_t cross_low = mulx_u64(l, h, &cross_high);
        const uint64_t doubled_low = cross_low << 1;
        const uint64_t doubled_high = (cross_high << 1) | (cross_low >> 63);
        const uint64_t doubled_top = cross_high >> 63;
        uint8_t carry = addcarryx_u64(0, r[1], doubled_low, &r[1]);
        carry = addcarryx_u64(carry, r[2], doubled_high, &r[2]);
        addcarryx_u64(carry, r[3], doubled_top, &r[3]);
    }

    /**
     * @brief Compares two unsigned 256 bit integers.
     * @param a Left hand side, a[0] lowest
     * @param b Right hand side, b[0] lowest
     * @return Negative, zero or positive as a is below, equal to or above b.
     */
    [[nodiscard]] FP128_INLINE static constexpr int32_t Compare256(const uint64_t (&a)[4], const uint64_t (&b)[4]) noexcept
    {
        for (int32_t i = 3; i >= 0; --i) {
            if (a[i] != b[i]) {
                return (a[i] < b[i]) ? -1 : 1;
            }
        }
        return 0;
    }

    /**
     * @brief Square root of an unsigned 256 bit integer, rounded to the nearest integer.
     *
     * A double precision estimate of the root, nudged above it, seeds the integer Newton iteration
     * r' = floor((r + floor(s / r)) / 2), which never drops below floor(sqrt(s)) once it starts above
     * it. The estimate is good to about 44 bits, so two steps leave r at floor(sqrt(s)) or one above
     * it. Comparing s with r^2 - r and r^2 + r then settles the rounded root exactly: it is r when
     * r^2 - r < s <= r^2 + r, the integer form of (r - 1/2)^2 <= s < (r + 1/2)^2. A root exactly
     * halfway between two integers is impossible, so there is no tie to break.
     *
     * Two 256 by 128 bit divisions do the work, the cost of the three Newton steps of the fixed point
     * iteration this replaced. The result no longer depends on how far the argument had to be
     * normalized: the old form shifted bits of a large argument away before taking the root and
     * shifted the root left afterwards, which lost up to half of the integer bits.
     *
     * @param s Radicand, s[0] lowest. At most 2^255, so the root fits in 128 bits.
     * @param res_low Receives the low QWORD of the rounded root.
     * @param res_high Receives the high QWORD of the rounded root.
     */
    FP128_INLINE static void RoundedSqrt256(const uint64_t (&s)[4], uint64_t& res_low, uint64_t& res_high) noexcept
    {
        res_low = res_high = 0;
        int32_t top = 3;
        while (top >= 0 && s[top] == 0) {
            --top;
        }
        if (top < 0) {
            return;
        }

        // s is about lead * 2^k, lead holding its 64 leading bits
        const int32_t lz = static_cast<int32_t>(lzcnt64(s[top]));
        uint64_t lead = s[top] << lz;
        if (lz != 0 && top != 0) {
            lead |= s[top - 1] >> (64 - lz);
        }
        const int32_t k = 64 * top - lz;
        // k >> 1 rounds towards minus infinity and k & 1 is the matching remainder, so the split is
        // right for a negative k as well. The factor puts the estimate safely above the root.
        const double estimate = std::ldexp(std::sqrt(std::ldexp(static_cast<double>(lead), k & 1)), k >> 1) * (1.0 + 0x1p-45);

        uint64_t r_low = 0, r_high = 0;
        if (estimate < 0x1p64) {
            r_low = static_cast<uint64_t>(estimate) + 2;  // makes up for the truncation, and then some
        } else {
            int exponent = 0;
            const double mantissa = std::frexp(estimate, &exponent);  // estimate = mantissa * 2^exponent
            r_low = static_cast<uint64_t>(std::ldexp(mantissa, 64));  // exact, it has 53 significant bits
            shift_left128_inplace_safe(r_low, r_high, exponent - 64);
            r_high += addcarryx_u64(0, r_low, 1, &r_low);
        }

        // Two Newton steps from above. The quotient fits in 128 bits: r stays at or above
        // floor(sqrt(s)), so s / r is at most two above it.
        for (int32_t i = 0; i < 2; ++i) {
            uint64_t q[4] {};
            if (r_high != 0) {
                const uint64_t divisor[2] = {r_low, r_high};
                (void)Divide256By128(q, nullptr, s, divisor);
            } else {
                (void)div_64bit(q, nullptr, s, r_low, 4);
            }
            // (r + q) / 2, the carry out of the sum shifted back in at the top
            const uint8_t sum_carry = addcarryx_u64(addcarryx_u64(0, r_low, q[0], &r_low), r_high, q[1], &r_high);
            r_low = (r_low >> 1) | (r_high << 63);
            r_high = (r_high >> 1) | (static_cast<uint64_t>(sum_carry) << 63);
        }

        // r is floor(sqrt(s)) or one more, so the rounded root is r - 1, r or r + 1
        uint64_t square[4] {};
        Square256(r_low, r_high, square);
        uint64_t below[4] {}, above[4] {};  // r^2 - r and r^2 + r
        uint8_t borrow = subborrow_u64(0, square[0], r_low, &below[0]);
        borrow = subborrow_u64(borrow, square[1], r_high, &below[1]);
        borrow = subborrow_u64(borrow, square[2], 0, &below[2]);
        (void)subborrow_u64(borrow, square[3], 0, &below[3]);
        uint8_t carry = addcarryx_u64(0, square[0], r_low, &above[0]);
        carry = addcarryx_u64(carry, square[1], r_high, &above[1]);
        carry = addcarryx_u64(carry, square[2], 0, &above[2]);
        (void)addcarryx_u64(carry, square[3], 0, &above[3]);

        if (Compare256(s, below) <= 0) {
            r_high -= subborrow_u64(0, r_low, 1, &r_low);
        } else if (Compare256(s, above) > 0) {
            r_high += addcarryx_u64(0, r_low, 1, &r_low);
        }

        res_low = r_low;
        res_high = r_high;
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
     * floor(log2(|x|)), as the C function: the sign of x makes no difference.
     * @param x The specified value.
     * @return The exponent. FP_ILOGB0 for zero, which is what the C function returns for it.
     */
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr int32_t ilogb(const fixed_point128& x) noexcept
    {
        return (x.is_zero()) ? FP_ILOGB0 : x.get_exponent();
    }
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
     * @brief Performs the fmod() function, similar to libc's fmod(): x - y * trunc(x / y).
     * Exact, and the result keeps the sign of x. See operator%=().
     * @param x Numerator
     * @param y Denominator
     * @return A fixed_point128 holding the modulo value.
     * @throws std::logic_error when y is zero. This type has no NaN to return instead.
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
     *
     * Correctly rounded, and free of intermediate overflow: hypot(x, y) * 2^F is the square root of
     * the sum of the squared raw values, which is formed as an exact 256 bit integer. Spelled as
     * sqrt(sqr(x) + sqr(y)) instead, it overflowed as soon as |x| or |y| reached 2^(I/2), long before
     * the result itself leaves the range.
     *
     * @param x First value
     * @param y Second value
     * @return sqrt(x^2 + y^2), rounded to nearest. A result of 2^I or more wraps, like any overflow.
     */
    [[nodiscard]] friend FP128_INLINE fixed_point128 hypot(const fixed_point128& x, const fixed_point128& y) noexcept
    {
        uint64_t x_low = 0, x_high = 0, y_low = 0, y_high = 0;
        (void)x.GetMagnitude(x_low, x_high);
        (void)y.GetMagnitude(y_low, y_high);

        // each magnitude is at most 2^127, so the sum of the squares is at most 2^255
        uint64_t sum[4] {}, y_square[4] {};
        Square256(x_low, x_high, sum);
        Square256(y_low, y_high, y_square);
        uint8_t carry = 0;
        for (int32_t i = 0; i < 4; ++i) {
            carry = addcarryx_u64(carry, sum[i], y_square[i], &sum[i]);
        }

        fixed_point128 res;
        RoundedSqrt256(sum, res.low, res.high);
        return res;
    }
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
     * @brief Calculates the square root, correctly rounded.
     *
     * sqrt(x) * 2^F is the square root of the raw value scaled up by another F bits, which makes it
     * the integer square root of a 256 bit number: RoundedSqrt256() computes it rounded to nearest,
     * so the result is the representable value closest to the exact root at every I and every
     * magnitude. (A tie between two representable values cannot occur.)
     *
     * @param x Value to calculate the root of
     * @param iterations Ignored. It used to set the step count of an iterative refinement, which lost
     *        up to half of the integer bits for a large argument; it is kept so that existing calls
     *        still compile.
     * @return Square root of (x), zero when x <= 0.
     */
    [[nodiscard]] friend FP128_INLINE fixed_point128 sqrt(const fixed_point128& x, [[maybe_unused]] uint32_t iterations = 3) noexcept
    {
        if (x.is_negative() || !x)
            return 0;

        // raw * 2^F: F is 64 + upper_frac_bits, so the raw value moves up one QWORD and that many bits
        uint64_t s[4] {};
        if constexpr (upper_frac_bits == 0) {
            s[1] = x.low;
            s[2] = x.high;
        } else {
            s[1] = x.low << upper_frac_bits;
            s[2] = (x.high << upper_frac_bits) | (x.low >> (64 - upper_frac_bits));
            s[3] = x.high >> (64 - upper_frac_bits);
        }

        fixed_point128 root;
        RoundedSqrt256(s, root.low, root.high);
        return root;
    }
    /**
     * @brief Factorial reciprocal (inverse). Calculates 1 / x!
     * Maximum value of x that may produce non zero values is 34.
     * This value depends on the amount of fraction bits.
     * The table is parsed from decimal strings, so it carries the accuracy of the string constructor:
     * within 0.5 + 23 / 2^(I-1) ulp of the exact 1 / x!, measured 17.5 ulp at I = 1 and 0.5 from
     * I = 16 up.
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
     *
     * Correctly rounded: the quotient comes from the same long division operator/=() uses, rounded to
     * nearest with a tie going away from zero. The division is the faster of the two since
     * div_128bit, and the Newton refinement of a double estimate this replaced lost accuracy as 1 / x
     * grew - tens of thousands of ulp at fixed_point128<16>. It is called directly rather than
     * through operator/=(), which divides by multiplying with this function when
     * FP128_USE_RECIPROCAL_FOR_DIVISION is set.
     *
     * @param x Input value
     * @return 1 / x. Zero when x is zero or when 1 / x does not fit the range [-2^I, 2^I), that is
     *         for |x| < 2^-I and for x == 2^-I; x == -2^-I gives lowest().
     */
    [[nodiscard]] friend FP128_INLINE fixed_point128 reciprocal(const fixed_point128& x) noexcept
    {
        uint64_t l = 0, h = 0;
        const uint64_t sign_bits = x.GetMagnitude(l, h);
        if ((l | h) == 0) {
            return fixed_point128(0);
        }

        // |x| = 2^(msb - F) and up, so 1 / |x| reaches 2^I exactly when msb is F - I or below.
        const int32_t msb = 127 - static_cast<int32_t>(lzcnt128(l, h));
        if (msb < F - I) {
            return fixed_point128(0);
        }
        if (msb == F - I && popcnt128(l, h) == 1) {
            // exactly 2^-I: the reciprocal is 2^I, which only the negative side of the range holds
            return (sign_bits != 0) ? fixed_point128(0, 1ull << 63) : fixed_point128(0);
        }

        fixed_point128 res = one();
        const bool need_rounding = res.DivideMagnitude(l, h);
        FP128_ADD_ROUND_BIT(res.low, res.high, need_rounding);
        res.NegateIf(sign_bits);
        return res;
    }
    /**
     * @brief Calculate Sine and Cosine using CORDIC using a limited range of [-pi/2, pi/2]
     * 64 iterations, one per entry of the arctangent table, which give an absolute accuracy of about
     * 2^-63: measured within 1.1e-19 of mpmath for |x| < 1. The table used to skip arctan(2^-62),
     * which held the error near 1.8e-19. Not used by tan(), which divides sin() by cos().
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
            "0.0000000000000000002168404344971008868015",  // arctan(2^-62)
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

        for (uint32_t j = 1; j <= angles_count; ++j) {
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
     * The series runs in a wider fraction than this instantiation has, see TrigSeries().
     * @param x value in Radians in the range [-0.5pi, 0.5pi]
     * @return Sine of x. See sin() for the accuracy.
     */
    [[nodiscard]] friend FP128_INLINE fixed_point128 sin1(fixed_point128 x) noexcept
    {
        static_assert(I >= 4, "sin1() needs template parameter <I> to be at least 4: its callers reduce the argument with pi, which needs more than one integer bit!");
        assert(fabs(x) <= fixed_point128::half_pi());

        // sine is odd: the series runs on the magnitude, which keeps every term non negative
        const bool negative = x.is_negative();
        const fixed_point128 res(TrigSeries(fixed_point128<2>(fabs(x)), false));
        return (negative) ? -res : res;
    }
    /**
     * @brief Calculate the cosine function over a limited range [-0.5pi, 0.5pi]
     * Using the Maclaurin series expansion, the formula is:
     *              x^2   x^4   x^6
     * cos(x) = 1 - --- + --- - --- + ...
     *               2!    4!    6!
     * Evaluated directly rather than as sin1(pi/2 - x), which carried the rounding error of pi/2 into
     * every result; cos1(0) is now exactly one. See TrigSeries().
     * @param x value in Radians in the range [-0.5pi, 0.5pi]
     * @return Cosine of x. See cos() for the accuracy.
     */
    [[nodiscard]] friend FP128_INLINE fixed_point128 cos1(const fixed_point128& x) noexcept
    {
        static_assert(I >= 4, "cos1() needs template parameter <I> to be at least 4: its callers reduce the argument with pi, which needs more than one integer bit!");
        assert(fabs(x) <= fixed_point128::half_pi());
        // cosine is even
        return fixed_point128(TrigSeries(fixed_point128<2>(fabs(x)), true));
    }
    /**
     * @brief Calculate the Sine function
     * Ultimately uses sin1() or cos1() with a reduced range of [-pi/4, pi/4].
     *
     * Accuracy, measured against mpmath at I = 1, 2, 4, 16, 32, 62 and 63 where the function exists: within 2.5 ulp
     * for |x| < 4 (1.4 ulp from I = 16 up). Beyond that the reduction x - n * pi/2 adds n times the
     * rounding error of pi/2 in this instantiation, at most |x| / 3 ulp more: 320 ulp for |x| near
     * 1000 at I = 32.
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
     * Newton's method on sin() for |x| up to one half, see AsinNewton(). Beyond that the identity
     * asin(x) = pi/2 - 2 asin(sqrt((1 - x) / 2)) keeps the iteration away from pi/2, where the
     * derivative it divides by vanishes and the refinement stalled at double precision: asin(1) was
     * 4e16 ulp off at fixed_point128<16>. The square root is correctly rounded, so the identity costs
     * a few ulp at most.
     *
     * Accuracy, measured against mpmath at I = 4, 16, 32, 62 and 63: within 3.5 ulp.
     * @param x value in radians in the range [-1,1]
     * @return Inverse sine of x, zero outside [-1,1]
     */
    [[nodiscard]] friend fixed_point128 asin(fixed_point128 x) noexcept
    {
        static_assert(I >= 4, "asin() needs template parameter <I> to be at least 4: it evaluates sin() and cos()!");
        if (x < -1 || x > 1)
            return 0;

        const bool negative = x.is_negative();
        x = fabs(x);
        fixed_point128 res;
        if (x <= fixed_point128::half()) {
            res = AsinNewton(x);
        } else {
            res = fixed_point128::half_pi() - (AsinNewton(RootOfHalfComplement(x)) << 1);
        }
        return (negative) ? -res : res;
    }
    /**
     * @brief Calculate the cosine function
     * Ultimately uses sin1() or cos1() with a reduced range of [-pi/4, pi/4]. cos(0) is exactly one.
     *
     * Accuracy, measured against mpmath at I = 1, 2, 4, 16, 32, 62 and 63 where the function exists: within 1.8 ulp
     * for |x| < 4 (1.3 ulp from I = 16 up), plus at most |x| / 3 ulp of argument reduction beyond, see
     * sin().
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
     * pi/2 - asin(x) for |x| up to one half. Beyond that acos(x) = 2 asin(sqrt((1 - x) / 2)) for a
     * positive x and pi - 2 asin(sqrt((1 + x) / 2)) for a negative one, which keep Newton's method
     * away from the ends of the domain: acos(-1) was 8e16 ulp off at fixed_point128<16>.
     *
     * Accuracy, measured against mpmath at I = 4, 16, 32, 62 and 63: within 4 ulp.
     * @param x value in radians in the range [-1,1]
     * @return Inverse cosine of x, zero outside [-1,1]
     */
    [[nodiscard]] friend fixed_point128 acos(fixed_point128 x) noexcept
    {
        static_assert(I >= 4, "acos() needs template parameter <I> to be at least 4: it evaluates sin() and cos()!");
        if (x < -1 || x > 1)
            return 0;

        const fixed_point128 magnitude = fabs(x);
        if (magnitude <= fixed_point128::half()) {
            const fixed_point128 res = AsinNewton(magnitude);
            return fixed_point128::half_pi() + ((x.is_negative()) ? res : -res);
        }

        const fixed_point128 twice = AsinNewton(RootOfHalfComplement(magnitude)) << 1;
        return (x.is_negative()) ? fixed_point128::pi() - twice : twice;
    }
    /**
     * @brief Calculate the tangent function
     * tan(x) = sin(x)/cos(x)
     *
     * The error of the cosine is divided by the cosine itself, so it grows towards the poles. Measured
     * against mpmath for |x| < 4 at I = 4, 16, 32, 62 and 63: within 45 * max(1, |tan(x)|) ulp.
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
     * Accuracy, measured against mpmath at I = 1, 2, 4, 16, 32, 62 and 63 where the function exists: within 2.2 ulp.
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

        // limit argument to 0..1. x is greater than one here, so the division cannot be by zero.
        if (x > 1) {
            comp = true;
            x = one() / x;
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
     * Accuracy, measured against mpmath at I = 1, 2, 4, 16, 32, 62 and 63 where the function exists: within 3.3 ulp.
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

        // The ratio is kept at or below one in magnitude. Steeper than the diagonal, the angle is
        // +-pi/2 - atan(x / y) with the sign of y, in every quadrant. That used to take its sign
        // from atan(x / y) instead, which lost it whenever x / y rounded to zero: atan2(1000, -1e-30)
        // came out as 3pi/2 at fixed_point128<16>.
        if (fabs(y) > fabs(x)) {
            return ((y.is_negative()) ? -half_pi : half_pi) - atan(x / y);
        }

        // flatter than the diagonal: atan(y / x), moved into the left half plane when x is negative
        const fixed_point128 res = atan(y / x);
        if (x.is_positive()) {
            return res;
        }

        return (y.is_negative()) ? res - pi : res + pi;
    }
    /**
     * @brief Calculate the hyperbolic sine function
     * Use the exponent function which produces more accurate results than the power series.
     *           e^x - e^(-x)
     * sinh(x) = ------------
     *                2
     * Accuracy, measured against mpmath at I = 1, 2, 4, 16, 32, 62 and 63 where the function exists: within
     * 10 * max(1, |sinh(x)|) ulp at I = 4 and 3 * max(1, |sinh(x)|) ulp from I = 16 up. e^x itself has
     * to fit the range, so a result between 2^(I-1) and 2^I overflows.
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
     * For negative x, the function returns the result with the sign inverted.
     * The root is exact before it is rounded (RootOfSquareAndOne()) and the sum is kept from overflowing
     * (LogOfSum()), so every argument works, lowest() included. Accuracy, measured against mpmath at
     * I = 1, 2, 4, 16, 32, 62 and 63 over the whole range: within 2.8 ulp.
     * @param x value
     * @return Inverse hyperbolic sine of x
     */
    [[nodiscard]] friend FP128_INLINE fixed_point128 asinh(const fixed_point128& x) noexcept
    {
        uint64_t l = 0, h = 0;
        (void)x.GetMagnitude(l, h);
        const fixed_point128 res = LogOfSum(l, h, RootOfSquareAndOne(x, false));
        return (x.is_negative()) ? -res : res;
    }
    /**
     * @brief Calculate the hyperbolic cosine function
     *           e^x + e^(-x)
     * cosh(x) = ------------
     *                2
     * Accuracy, measured against mpmath at I = 1, 2, 4, 16, 32, 62 and 63 where the function exists: within
     * 10 * max(1, cosh(x)) ulp at I = 4 and 3 * max(1, cosh(x)) ulp from I = 16 up. e^|x| itself has to
     * fit the range, see sinh().
     * @param x value
     * @return Hyperbolic cosine of x
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
     * Accuracy, measured against mpmath at I = 1, 2, 4, 16, 32, 62 and 63 over the whole range: within
     * 2.9 ulp. It was out by up to 5e8 ulp just above one, where x^2 - 1 cancels.
     * @param x value in the range [1, inf]
     * @return Inverse hyperbolic cosine of x
     */
    [[nodiscard]] friend FP128_INLINE fixed_point128 acosh(const fixed_point128& x) noexcept
    {
        if (x < 1)
            return 0;

        // The root of x^2 - 1 is exact before rounding, see RootOfSquareAndOne(): formed in this type
        // the difference cancels for an x near one and the root magnifies its rounding error. The sum
        // is kept from overflowing by LogOfSum().
        return LogOfSum(x.low, x.high, RootOfSquareAndOne(x, true));
    }
    /**
     * @brief Calculates the hyperbolic tangent
     *           e^x - e^(-x)     1 - e^(-2|x|)
     * tanh(x) = ------------ = ---------------, with the sign of x
     *           e^x + e^(-x)     1 + e^(-2|x|)
     * The second form is the one evaluated, for every argument. Accuracy, measured against mpmath at I = 1, 2, 4, 16, 32, 62 and 63 where the function exists:
     * within 4.2 ulp.
     * @param x value
     * @return hyperbolic tangent of x
     */
    [[nodiscard]] friend FP128_INLINE fixed_point128 tanh(const fixed_point128& x) noexcept
    {
        static_assert(I >= 2, "tanh() needs template parameter <I> to be at least 2: it evaluates exp()!");
        // (1 - e^(-2|x|)) / (1 + e^(-2|x|)) with the sign of x: only e^(-2|x|), which is at most one,
        // is ever formed. The quotient of e^x - e^-x and e^x + e^-x this replaced overflowed with e^x,
        // which happens once |x| reaches I ln(2) - tanh(12) was garbage at fixed_point128<16>.
        const fixed_point128 magnitude = fabs(x);
        fixed_point128 e2;
        if (magnitude.is_positive() && magnitude.get_exponent() < I - 1) {
            e2 = exp(-(magnitude << 1));
        } else {
            // 2|x| does not fit the range (lowest() included, whose fabs() wraps): square e^-|x| instead
            e2 = sqr(exp(-magnitude));
        }
        const fixed_point128 res = (fixed_point128::one() - e2) / (fixed_point128::one() + e2);
        return (x.is_negative()) ? -res : res;
    }
    /**
     * @brief Calculates the inverse hyperbolic tangent
     *                       1 + x
     * atanh(x) = 0.5 * log( -----) = 0.5 * (log(1 + x) - log(1 - x))
     *                       1 - x
     * The second form is the one evaluated. Both of its logarithms take an argument in (0, 2), inside
     * the range of every instantiation, while the quotient of the first form leaves the range as x
     * approaches one - beyond 0.88 at fixed_point128<4>, within 2^-15 of one at fixed_point128<16> -
     * and the wrapped value made log() throw out of this noexcept function, terminating the process.
     * Accuracy, measured against mpmath at I = 4, 16, 32, 62 and 63: within 1.4 ulp. Below I = 4 the
     * same holds where log(1 - x) fits the range, which log() needs for its log2 (see log()): for
     * |x| up to 0.75 at I = 1 and 0.9375 at I = 2.
     * @param x value in the range (-1, 1)
     * @return Inverse hyperbolic tangent of x, zero outside (-1, 1). A result beyond the range of the
     *         instantiation, which only a narrow one can meet, wraps like any overflow.
     */
    [[nodiscard]] friend FP128_INLINE fixed_point128 atanh(const fixed_point128& x) noexcept
    {
        const fixed_point128 one = fixed_point128::one();
        // compared without fabs(), which leaves lowest() negative
        if (x <= -one || x >= one)
            return 0;

        return (log(one + x) - log(one - x)) >> 1;
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
     * A negative x works with 1/e and the alternating series of e^-fx, so it never forms e^|x|.
     *
     * Accuracy, measured against mpmath at I = 1, 2, 4, 16, 32, 62 and 63 where the function exists: within
     * 6.2 * max(1, e^x) ulp. The argument itself is only known to half an ulp, which alone moves a
     * large result by e^x / 2 ulp.
     * @param x A number specifying a power.
     * @return Exponent of x
     */
    [[nodiscard]] friend FP128_INLINE fixed_point128 exp(const fixed_point128& x) noexcept
    {
        static_assert(I >= 2, "exp() needs template parameter <I> to be at least 2: it multiplies by e(), which needs 2 integer bits!");
        static const fixed_point128 e = fixed_point128::e();
        // A negative exponent works with 1/e and an alternating series rather than inverting e^|x|:
        // e^|x| overflows long before e^-|x| leaves the grid - exp(-12) needs e^12 = 162755, beyond
        // the range of fixed_point128<16>, and came out as garbage.
        const bool negative = x.is_negative();
        fixed_point128 int_part;
        const fixed_point128 fx = modf(fabs(x), &int_part);  // fraction part of |x|
        uint64_t ix = static_cast<uint64_t>(int_part);        // 64 bit is an overkill to hold the exponent

        // e^ix, or e^-ix, by binary powering
        fixed_point128 exp_ix = fixed_point128::one();
        fixed_point128 b = (negative) ? FromFraction128(inv_e_fraction[0], inv_e_fraction[1]) : e;
        while (ix > 0) {
            if (ix & 1) {
                exp_ix *= b;
            }
            ix >>= 1;
            if (ix > 0) {
                b.square();
            }
        }

        // e^fx, or e^-fx, by the Maclaurin series: fx < 1, so the terms shrink and none overflows
        if (!fx)
            return exp_ix;

        fixed_point128 exp_fx = (negative) ? fixed_point128::one() - fx : fixed_point128::one() + fx;
        fixed_point128 elem_denom, elem_nom = fx;
        for (int i = 2;; ++i) {
            elem_nom *= fx;
            fact_reciprocal(i, elem_denom);
            const fixed_point128 elem = elem_nom * elem_denom;
            if (!elem)
                break;
            if (negative && (i & 1)) {
                exp_fx -= elem;
            } else {
                exp_fx += elem;
            }
        }

        return exp_ix * exp_fx;
    }
    /**
     * @brief Calculates the exponent of x and reduces 1 from the result: (e^x) - 1
     * Evaluated as exp(x) - 1, which in a fixed point type loses nothing to cancellation: the
     * absolute error is that of exp(), within 9.3 * max(1, |e^x - 1|) ulp measured at I = 4 and
     * 6.2 * max(1, e^x) ulp from I = 16 up.
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
     * exp(x * ln(2)), the constant applied with 128 fraction bits. Accuracy, measured against mpmath at I = 1, 2, 4, 16, 32, 62 and 63 where the function exists:
     * within 5.5 * max(1, 2^x) ulp.
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
        return exp(x.MulByFraction128(ln2_fraction[1], ln2_fraction[0]));
    }
    /**
     * @brief Computes x to the power of y
     * exp(y * log(x)). Accuracy, measured against mpmath at I = 1, 2, 4, 16, 32, 62 and 63 where the function exists: within
     * 7.3 * max(1, x^y) ulp at I = 4 and 4.4 * max(1, x^y) ulp from I = 16 up. Exact powers are not
     * exact: pow(2, 3) is 8 - 15 ulp at I = 16. y * log(x) has to fit the range.
     * @param x Base value, must be positive
     * @param y Exponent value
     * @return x^y, one when x is one, zero when x is negative, and for x = 0 one when y = 0 and zero
     *         otherwise
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

        // log(1) is zero, and so is the log of a value within an ulp or two of one: x^y is one there
        // for every y. Returning the zero itself, as this used to, made pow(1, y) zero.
        const fixed_point128 lan_x = log(x);
        if (!lan_x)
            return fixed_point128::one();

        return exp(y * lan_x);
    }
    /**
     * @brief Calculates the Log base 2 of x: y = log2(x)
     * Accuracy, measured against mpmath at I = 1, 2, 4, 16, 32, 62 and 63 where the function exists: within 2 ulp
     * at I = 1, 1.7 ulp at I = 2 and 1.1 ulp from I = 4 up. The compile time evaluation takes a
     * different, bitwise path that truncates instead: it can be one ulp below the runtime result.
     * A result outside the range wraps, which a narrow instantiation meets for small arguments:
     * log2(x) needs x >= 2^(-2^I).
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
     * log2(x) * ln(2), the constant applied with 128 fraction bits by MulByFraction128() - a constant
     * held to this instantiation's F bits, multiplied by a log2 of up to 63, would cost tens of ulp.
     * Accuracy, measured against mpmath at I = 1, 2, 4, 16, 32, 62 and 63 where the function exists: within 1.9 ulp.
     * log2(x) has to fit the range as well, so a narrow instantiation cannot take the log of an
     * argument below 2^(-2^I) even when its log would fit: at I = 1 that is x below 0.25.
     * @param x The number to perform log on.
     * @return log(x)
     * @throws std::domain_error when x is zero or negative.
     */
    [[nodiscard]] friend FP128_INLINE fixed_point128 log(fixed_point128 x)
    {
        if (x.is_zero() || x.is_negative()) {
            throw std::domain_error("Math domain error! Function accepts positive, non-zero values only.");
        }
        if (x == 1)
            return 0;

        return log2(x).MulByFraction128(ln2_fraction[1], ln2_fraction[0]);
    }
    /**
     * @brief Calculates the natural Log (base e) of 1 + x: log(1 + x)
     * Not noexcept: it forwards to log(), which rejects a non positive argument. Marking it
     * noexcept would turn that domain_error into a call to std::terminate. 1 + x is exact in a fixed
     * point type, so the error is that of log(): within 1.2 ulp measured, with log()'s limit on small
     * arguments at I < 7.
     * @param x The number to perform log on.
     * @return log1p(x)
     * @throws std::domain_error when x is -1 or below.
     */
    [[nodiscard]] friend FP128_FORCE_INLINE fixed_point128 log1p(fixed_point128 x) { return log(fixed_point128::one() + x); }
    /**
     * @brief Calculates Log base 10 of x: log10(x)
     * log2(x) * log10(2), the constant applied with 128 fraction bits, see log(). Accuracy, measured against mpmath at I = 1, 2, 4, 16, 32, 62 and 63 where the function exists:
     * within 1 ulp, with log()'s limit on small arguments at I < 7.
     * @param x The number to perform log on.
     * @return log10(x)
     */
    [[nodiscard]] friend FP128_INLINE fixed_point128 log10(fixed_point128 x)
    {
        if (x.is_zero() || x.is_negative()) {
            throw std::domain_error("Math domain error! Function accepts positive, non-zero values only.");
        }

        if (x == 1)
            return 0;

        // log2(x) * log10(2), the constant applied with 128 fraction bits, see log()
        return log2(x).MulByFraction128(log10_2_fraction[1], log10_2_fraction[0]);
    }
    /**
     * @brief Calculates Log base 2 of x as an integer ignoring the sign of x.
     * floor(log2(fabs(x))), as the C function: logb(-8) is 3. A narrow instantiation cannot hold the
     * exponent of its smallest values - logb(epsilon()) is -126 for fixed_point128<1> - and the result
     * then wraps like any other overflow.
     * @param x The number to perform log on.
     * @return logb(x)
     * @throws std::domain_error when x is zero. The C function returns minus infinity for it, which
     *         this type cannot represent.
     */
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr fixed_point128 logb(fixed_point128 x)
    {
        if (x.is_zero()) {
            throw std::domain_error("Math domain error! logb() of zero is minus infinity.");
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
 * @brief Reads a fixed_point128 from a stream, the way the builtin arithmetic extractors do.
 *
 * Reads what the constructor from a string reads - an optional sign, digits, an optional point and
 * more digits - correctly rounded to the grid, and stops at the first character that is not part
 * of it, leaving it in the stream: "2.5,3" gives 2.5 and leaves ",3".
 *
 * Failure sets failbit: no digits at all stores zero, and a value outside [lowest(), max()] stores
 * the nearer of the two. Leading white space is skipped unless skipws is off. This used to read a
 * whole white space delimited token through the string constructor, which never failed: "xyz" read
 * as zero and an overflow wrapped silently.
 * @tparam I Integer bit count of the type
 */
template <int32_t I> inline std::istream& operator>>(std::istream& is, fixed_point128<I>& value)
{
    const std::istream::sentry sentry(is);
    if (!sentry)
        return is;

    using traits = std::char_traits<char>;
    std::streambuf* const buffer = is.rdbuf();
    std::string text;
    int c = buffer->sgetc();
    const auto take = [&]() {
        text += traits::to_char_type(c);
        buffer->sbumpc();
        c = buffer->sgetc();
    };
    if (c == '+' || c == '-') {
        take();
    }
    while (c != traits::eof() && c >= '0' && c <= '9') {
        take();
    }
    if (c == '.') {
        take();
        while (c != traits::eof() && c >= '0' && c <= '9') {
            take();
        }
    }
    if (c == traits::eof()) {
        is.setstate(std::ios_base::eofbit);
    }

    uint64_t low = 0, high = 0;
    const char* end = nullptr;
    const int32_t status = fixed_point128<I>::ParseDecimal(text.data(), text.data() + text.size(), low, high, end);
    if (status == 0) {
        value = fixed_point128<I>();
        is.setstate(std::ios_base::failbit);
    } else if (status == 2) {
        value = (!text.empty() && text[0] == '-') ? std::numeric_limits<fixed_point128<I>>::lowest() : std::numeric_limits<fixed_point128<I>>::max();
        is.setstate(std::ios_base::failbit);
    } else {
        value = fixed_point128<I>(low, high);
    }
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
 *
 * The exponent members read the grid the way the floating point definitions read a significand:
 * the smallest positive value epsilon() = 2^(I-127) is 2^(min_exponent-1), and max() is just under
 * 2^max_exponent = 2^I. min() follows the floating point reading too and is the smallest positive
 * value; lowest() is the most negative one.
 *
 * Every member is usable in a constant expression, as the standard requires of a specialization.
 * The const, volatile and const volatile forms need no specialization of their own: \<limits\>
 * already defines numeric_limits<cv T> to have the members of numeric_limits<T>.
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
    // has_denorm and has_denorm_loss are deprecated since C++23 but remain part of the interface a
    // generic caller may read, so they are provided.
    FP128_SUPPRESS_DEPRECATED_BEGIN
    static constexpr bool has_denorm_loss = false;
    static constexpr float_denorm_style has_denorm = denorm_absent;
    FP128_SUPPRESS_DEPRECATED_END
    /// @brief Multiplication, division (by a fixed_point128 or an integer), the right shift and the
    ///        narrowing conversions round to the nearest grid value, a tie going away from zero.
    static constexpr float_round_style round_style = round_to_nearest;
    static constexpr bool is_iec559 = false;
    static constexpr bool is_bounded = true;
    /**
     * @brief +, - and * wrap around: they keep the low 128 bits of the two's complement result,
     *        which differs from the exact one by a multiple of 2^(I+1) = max() - lowest() + epsilon().
     *
     * As for int128_t, whose builtin counterparts report false only because their overflow is
     * undefined. Here it is defined: max() + epsilon() is lowest(). This used to say false.
     */
    static constexpr bool is_modulo = true;
    /// @brief Division by zero throws std::logic_error rather than producing a value.
    static constexpr bool traps = true;
    static constexpr bool tinyness_before = false;

    /// @brief Value bits. One of the 128 is the sign, so 127 of them carry magnitude.
    static constexpr int digits = 127;
    /// @brief Decimal digits that survive a round trip: floor(127 * log10(2)).
    static constexpr int digits10 = 38;
    static constexpr int radix = 2;
    /// @brief 2^(min_exponent-1) = 2^(I-127) is the smallest positive value, epsilon().
    static constexpr int min_exponent = I - 126;
    /// @brief 2^(max_exponent-1) = 2^(I-1) is the largest power of two in range; max() is 2^I - epsilon().
    static constexpr int max_exponent = I;
    static constexpr int min_exponent10 = static_cast<int>((I - 127) * 0.30102999566398119521);
    static constexpr int max_exponent10 = static_cast<int>(I * 0.30102999566398119521);
    /**
     * @brief Significant decimal digits that tell every pair of neighbouring values apart.
     *
     * The top decade, from 10^max_exponent10 to max(), needs the most. Its values have
     * max_exponent10 + 1 integer digits, and the fraction needs floor((127 - I) * log10(2)) + 1
     * decimal places before their spacing falls below the grid spacing 2^(I-127). The total is 40
     * for the I where the two floors add up to 38 (4, 7, 10, 14, ..., 57, 60) and 39 for the rest.
     */
    static constexpr int max_digits10 = max_exponent10 + static_cast<int>((127 - I) * 0.30102999566398119521) + 2;

    /// @brief Smallest positive value, which is one unit in the last place.
    [[nodiscard]] static constexpr value_type min() noexcept { return value_type::epsilon(); }
    /// @brief Largest value, every bit below the sign set: 2^I - epsilon().
    [[nodiscard]] static constexpr value_type max() noexcept { return value_type(UINT64_MAX, 0x7FFFFFFFFFFFFFFFull); }
    /// @brief Most negative value, exactly -2^I. Two's complement makes the range asymmetric, so
    ///        this one has no positive counterpart and negating it wraps back to itself.
    [[nodiscard]] static constexpr value_type lowest() noexcept { return value_type(0, 1ull << 63); }
    /// @brief Spacing of the grid, the same at every magnitude.
    [[nodiscard]] static constexpr value_type epsilon() noexcept { return value_type::epsilon(); }
    /// @brief Largest rounding error of an operation in units in the last place, one half.
    [[nodiscard]] static constexpr value_type round_error() noexcept { return value_type::half(); }
    [[nodiscard]] static constexpr value_type infinity() noexcept { return value_type(); }
    [[nodiscard]] static constexpr value_type quiet_NaN() noexcept { return value_type(); }
    [[nodiscard]] static constexpr value_type signaling_NaN() noexcept { return value_type(); }
    /// @brief No value is subnormal, so this is the smallest positive one, the same as min().
    [[nodiscard]] static constexpr value_type denorm_min() noexcept { return value_type::epsilon(); }
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
