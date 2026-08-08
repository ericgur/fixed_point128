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
 * @file float128.h
 * @brief 128-bit IEEE 754-2008 quadruple-precision floating-point type.
 *
 * Provides the @ref fp128::float128 class with the same bit layout as
 * IEEE 754 binary128: 112-bit fraction, 15-bit exponent, and 1-bit sign.
 * Supports full arithmetic, comparison, and conversion operators, as well as
 * a comprehensive set of math functions (sqrt, cbrt, sin, cos, tan, asin,
 * acos, atan, atan2, sinh, cosh, tanh, asinh, acosh, atanh, exp, exp2, log,
 * log2, log10, pow, erf, erfc, etc.).
 * All methods are inline for maximum performance.
 *
 * @see fp128_shared.h for supporting intrinsics and utilities.
 */

#ifndef FP128_FLOAT128_H
#define FP128_FLOAT128_H

#include <algorithm>
#include <array>
#include <cmath>    // FP_NAN and friends, and the double seeds of sqrt/cbrt
#include <climits>  // INT_MAX, the ilogb answer for an infinity
#include <cstdlib>  // strtoull, for the nan() payload
#include <charconv>   // to_chars_result, from_chars_result, chars_format
#include <format>     // std::formatter
#include <functional> // std::hash
#include <istream>
#include <ostream>
#include "fp128_shared.h"
#include "fp128_decimal.h"
#include "uint128_t.h"

namespace fp128
{
/***********************************************************************************
 *                                  Forward declarations
 ************************************************************************************/

class fp128_gtest;  // Google test class
class float128;

/// @brief Reads a decimal string into a float128. Defined below, declared here because the
///        constructor from a string delegates to it.
std::from_chars_result from_chars(const char* first, const char* last, float128& value);

namespace detail
{
/// @brief Renders a float128 as text. Defined below, declared here because the conversions to a
///        string delegate to it.
[[nodiscard]] std::string render(const float128& value, char style, int32_t precision, bool uppercase, bool alternate, char sign_char);
}  // namespace detail

/// @brief User-defined literal for constructing float128 from a string.
float128 operator""_f128(const char*);

/// @name CRT-Style Math Functions (Forward Declarations)
/// @brief Free functions providing standard math library equivalents for float128.
/// @{
constexpr float128 fabs(const float128& x) noexcept;
constexpr float128 floor(const float128& x) noexcept;
constexpr float128 ceil(const float128& x) noexcept;
constexpr float128 trunc(const float128& x) noexcept;
constexpr float128 round(const float128& x) noexcept;
constexpr int64_t llrint(const float128& x) noexcept;
constexpr int64_t llround(const float128& x) noexcept;
constexpr int32_t lrint(const float128& x) noexcept;
constexpr int32_t lround(const float128& x) noexcept;
constexpr int32_t ilogb(const float128& x) noexcept;
constexpr float128 copysign(const float128& x, const float128& y) noexcept;
float128 fmod(const float128& x, const float128& y) noexcept;
constexpr float128 modf(const float128& x, float128* iptr) noexcept;
constexpr float128 fdim(const float128& x, const float128& y) noexcept;
constexpr float128 fmin(const float128& x, const float128& y) noexcept;
constexpr float128 fmax(const float128& x, const float128& y) noexcept;
constexpr float128 fma(float128 x, float128 y, float128 z) noexcept;
float128 hypot(const float128& x, const float128& y) noexcept;
float128 cbrt(const float128 x, uint32_t iterations = 2) noexcept;
float128 sqrt(const float128& x, uint32_t iterations = 3) noexcept;
float128 erf(float128 x) noexcept;
float128 erfc(float128 x) noexcept;
float128 sin(float128 x) noexcept;
float128 asin(float128 x) noexcept;
float128 cos(float128 x) noexcept;
float128 acos(float128 x) noexcept;
float128 tan(float128 x) noexcept;
float128 atan(float128 x) noexcept;
float128 atan2(float128 y, float128 x) noexcept;
float128 sinh(float128 x) noexcept;
float128 asinh(float128 x) noexcept;
float128 cosh(float128 x) noexcept;
float128 acosh(float128 x) noexcept;
float128 tanh(float128 x) noexcept;
float128 atanh(float128 x) noexcept;
float128 exp(const float128& x) noexcept;
float128 exp2(const float128& x) noexcept;
float128 expm1(const float128& x) noexcept;
float128 pow(const float128& x, const float128& y) noexcept;
float128 pow(const float128& x, int32_t y) noexcept;
float128 log(float128 x) noexcept;
float128 log2(float128 x) noexcept;
float128 log10(float128 x) noexcept;
float128 logb(float128 x) noexcept;
float128 log1p(float128 x) noexcept;
constexpr float128 frexp(float128 x, int* expptr) noexcept;
constexpr float128 ldexp(float128 x, int exp) noexcept;
constexpr bool isfinite(const float128& x) noexcept;
constexpr bool signbit(const float128& x) noexcept;
constexpr bool isnormal(const float128& x) noexcept;
constexpr int fpclassify(const float128& x) noexcept;
constexpr bool isunordered(const float128& x, const float128& y) noexcept;
constexpr bool isgreater(const float128& x, const float128& y) noexcept;
constexpr bool isgreaterequal(const float128& x, const float128& y) noexcept;
constexpr bool isless(const float128& x, const float128& y) noexcept;
constexpr bool islessequal(const float128& x, const float128& y) noexcept;
constexpr bool islessgreater(const float128& x, const float128& y) noexcept;
constexpr float128 rint(const float128& x) noexcept;
constexpr float128 nearbyint(const float128& x) noexcept;
constexpr float128 scalbn(const float128& x, int n) noexcept;
constexpr float128 scalbln(const float128& x, long n) noexcept;
constexpr float128 nextafter(const float128& x, const float128& y) noexcept;
constexpr float128 nexttoward(const float128& x, const float128& y) noexcept;
float128 remainder(const float128& x, const float128& y) noexcept;
float128 remquo(const float128& x, const float128& y, int* quo) noexcept;
float128 hypot(const float128& x, const float128& y, const float128& z) noexcept;
float128 lgamma(float128 x) noexcept;
float128 tgamma(float128 x) noexcept;
constexpr float128 abs(const float128& x) noexcept;
/// The only one of these that a call cannot reach through argument dependent lookup: its argument
/// is a plain character pointer, so the name has to be visible at namespace scope - and a call has
/// to be qualified as fp128::nan(...) wherever <cmath> puts its own nan(const char*) in scope too.
float128 nan(const char* payload) noexcept;
/// @}

/// @name Non-CRT Utility Functions (Forward Declarations)
/// @{
float128 reciprocal(const float128& x) noexcept;
void fact_reciprocal(int x, float128& res) noexcept;
float128 double_factorial(int x) noexcept;
/// @}

/***********************************************************************************
 *                                  Main Code
 ************************************************************************************/

/**
 * @brief 128 bit floating point class.
 *
 * This class implements the standard operators a floating point data type.<BR>
 * All of float128's methods are inline for maximum performance.
 *
 * <B>Implementation notes:</B>
 * <UL> Same bit layout as binary128: 112 bit fraction, 15 bit exponent and 1 bit for the sign.
 * <LI>A float128 object is not thread safe. Accessing a const object from multiple threads is safe.</LI>
 * <LI>Only 64 bit builds are supported.</LI>
 * </UL>
 *
 * <B>Compile time evaluation:</B><BR>
 * Everything that stays within the 128 bit encoding is constexpr:
 * <UL>
 * <LI>Construction from any builtin arithmetic type and from the raw QWORDs, copy, move,
 *     assignment and the conversions to the integer and floating point types.</LI>
 * <LI>Addition, subtraction, multiplication, square(), the shifts, the unary operators and the
 *     comparisons.</LI>
 * <LI>The queries and the component accessors: is_zero(), is_normal(), is_nan(), is_inf(),
 *     is_int(), get_exponent(), get_class(), get_components(), set_components() and the rest,
 *     along with every named constant (one(), pi(), inf(), nan(), ...).</LI>
 * <LI>The math functions fabs, floor, ceil, trunc, round, llround, lround, llrint, lrint, ilogb,
 *     copysign, modf, fdim, fmin, fmax, fma, frexp, ldexp, sqr, isnan, isinf, isfinite, and the
 *     nextUp/nextDown/exp10 helpers.</LI>
 * </UL>
 *
 * Two things make that possible. The bit counting and extended arithmetic intrinsics are not
 * constant expressions, so fp128_shared.h wraps each one in a constexpr function that
 * serves a constant evaluated call from a portable implementation of the same operation. And the
 * fields of the high QWORD are read through the shift and mask accessors rather than the
 * _float128_bits view, because reading the inactive member of a union is not allowed during
 * constant evaluation. A runtime call is unaffected by either.
 *
 * The rest cannot be constexpr, for one of two reasons:
 * <UL>
 * <LI>Division and modulo, and everything built on them: div_32bit needs alloca and a goto,
 *     neither of which C++20 permits in a constexpr function.</LI>
 * <LI>The string conversions allocate, and the transcendental functions parse their constants
 *     from strings held in function local statics, which a constexpr function may not declare.</LI>
 * </UL>
 */

class FP128_ALIGN16 float128
{
    // build time validation of template parameters
    static_assert(sizeof(void*) == 8, "float128 is supported in 64 bit builds only!");
    friend class fp128_gtest;

    static constexpr int32_t EXP_BITS = 15;            ///< Number of exponent bits in binary128.
    static constexpr int32_t EXP_BIAS = 0x3FFF;         ///< Exponent bias (16383) for binary128.
    static constexpr int32_t ZERO_EXP_BIASED = -EXP_BIAS;        ///< Biased exponent value representing zero.
    static constexpr int32_t ZERO_EXP_UNBIASED = 0;              ///< Unbiased exponent value for zero encoding.
    static constexpr int32_t SUBNORM_EXP_BIASED = 0;             ///< Biased exponent value for subnormal numbers.
    static constexpr int32_t SUBNORM_EXP_UNBIASED = -EXP_BIAS;   ///< Unbiased exponent value for subnormal numbers.
    static constexpr int32_t INF_EXP_BIASED = 0x7FFF;            ///< Biased exponent value for infinity/NaN.
    static constexpr int32_t INF_EXP_UNBIASED = INF_EXP_BIASED - EXP_BIAS; ///< Unbiased exponent value for infinity/NaN.
    static constexpr uint64_t EXP_MASK = INF_EXP_BIASED;         ///< Bitmask for the exponent field.
    static constexpr int32_t FRAC_BITS = 112;                    ///< Number of fraction (mantissa) bits.
    static constexpr int32_t EXP_SHIFT = FRAC_BITS - 64;         ///< Bit position of the exponent field within the high QWORD.
    static constexpr uint64_t UPPER_FRAC_MASK = FP128_MAX_VALUE_64(FRAC_BITS - 64); ///< Bitmask for upper fraction bits within the high QWORD.
    static constexpr uint64_t FRAC_UNITY = FP128_ONE_SHIFT(FRAC_BITS - 64);         ///< The implicit unity bit position in the high QWORD.
    static constexpr uint64_t SIGN_MASK = 1ull << 63;            ///< Bitmask for the sign bit.
    static constexpr uint64_t QUIET_NAN_BIT = 1ull << (EXP_SHIFT - 1);  ///< Leading fraction bit, set on a quiet NaN.

    /// @brief Bit-field view of the upper 64 bits of a float128.
    ///
    /// Retained for the debugger visualizer (fixed_point128.natvis), which reads the three fields
    /// by name. The code itself reaches the same fields through the shift and mask accessors
    /// below: reading the inactive member of a union is not allowed during constant evaluation,
    /// and every operation on this type funnels through those accessors.
    struct _float128_bits {
        uint64_t f : 48; ///< Upper 48 bits of the 112-bit fraction.
        uint64_t e : 15; ///< 15-bit biased exponent.
        uint64_t s : 1;  ///< Sign bit (0 = positive, 1 = negative).
    };

    /// @brief IEEE 754 classification categories for float128 values.
    enum float128_class_t {
        signalingNaN,       ///< Signaling NaN.
        quietNaN,           ///< Quiet NaN.
        negativeInfinity,   ///< Negative infinity.
        negativeNormal,     ///< Negative normal number.
        negativeSubnormal,  ///< Negative subnormal (denormalized) number.
        negativeZero,       ///< Negative zero.
        positiveZero,       ///< Positive zero.
        positiveSubnormal,  ///< Positive subnormal (denormalized) number.
        positiveNormal,     ///< Positive normal number.
        positiveInfinity    ///< Positive infinity.
    };

    /// @brief Internal storage: 128-bit value split into low and high QWORDs.
    struct {
        uint64_t low;  ///< Lower 64 bits of the float128 encoding.
        union {
            uint64_t high;              ///< Upper 64 bits (raw).
            _float128_bits high_bits;   ///< Upper 64 bits (bit-field view).
        };
    };

public:
    /**
     * @brief Default constructor, creates an instance with a value of zero.
     */
    FP128_FORCE_INLINE constexpr float128() noexcept : low(0), high(0) {}
    /**
     * @brief Copy constructor
     * @param other Object to copy from
     */
    FP128_FORCE_INLINE constexpr float128(const float128& other) noexcept : low(other.low), high(other.high) {}
    /**
     * @brief Move constructor
     * Doesn't modify the right hand side object. Acts like a copy constructor.
     * @param other Object to copy from
     */
    FP128_FORCE_INLINE constexpr float128(float128&& other) noexcept : low(other.low), high(other.high) {}
    /**
     * @brief Low level constructor
     * @param l Low QWORD
     * @param h High QWORD
     */
    FP128_FORCE_INLINE constexpr float128(uint64_t l, uint64_t h) noexcept : low(l), high(h) {}
    /**
     * @brief Low level constructor
     * @param lf Low fraction part (bits 63:0)
     * @param hf High fraction part (bits 111:64)
     * @param e Exponent
     * @param s sign
     */
    FP128_FORCE_INLINE constexpr float128(uint64_t lf, uint64_t hf, uint32_t e, uint32_t s) noexcept :
        low(lf), high((hf & FP128_MAX_VALUE_64(48)) | (((0x7FFFull & e) << 48)) | ((1ull & s) << 63))
    {
    }
    /**
     * @brief Constructor from the double type
     * @param x Input value
     */
    FP128_INLINE constexpr float128(double x) noexcept
    {
        low = high = 0;
        // very common case
        if (x == 0)
            return;

        // hack the double bit fields
        const Double d(x);

        // subnormal numbers
        if (d.e() == 0) {
            // the exponent is -1022 (1-1023)
            auto msb = 64 - static_cast<int32_t>(lzcnt64(d.f()));
            // exponent
            int32_t x_expo = static_cast<int32_t>(d.e()) - 1023;
            int32_t expo = x_expo + msb - dbl_frac_bits;
            // fraction
            low = d.f() & ~(1ull << (msb - 1));  // clear the msb
            auto shift = static_cast<int32_t>(FRAC_BITS - msb + 1);
            shift_left128_inplace_safe(low, high, shift);
            set_exponent(expo);
        }
        // NaN & INF
        else if (d.e() == 0x7FF) {
            set_exponent_bits(INF_EXP_BIASED);
            // zero for +- INF, non-zero for NaN
            set_fraction_bits((d.f()) ? 1 : 0);
        }
        // normal numbers
        else {
            low = d.f() << 60;
            high = d.f() >> 4;
            set_exponent(static_cast<int32_t>(d.e()) - 1023);
        }

        // copy the sign
        set_sign(d.s());
    }
    /**
     * @brief Generic constructor for integral and floating-point types.
     * Floating-point types are converted via double. Integral types are converted directly.
     * Character pointer types delegate to the const char* constructor.
     * @tparam T Source type (must be arithmetic or a character pointer type)
     * @param x Input value
     */
    template <typename T> constexpr float128(T x) noexcept : low(0), high(0)
    {
        if constexpr (std::is_floating_point_v<T>) {
            *this = float128(static_cast<double>(x));
            return;
        } else if constexpr (std::is_same_v<char*, T> || std::is_same_v<unsigned char*, T> || std::is_same_v<const unsigned char*, T>) {
            *this = float128(static_cast<const char*>(x));
            return;
        } else if constexpr (std::is_integral_v<T>) {
            uint64_t sign = 0;
            // The magnitude is taken in the unsigned domain, where the conversion sign extends and
            // the negation wraps. Negating the signed value instead is undefined for the most
            // negative one, which has no positive counterpart, and produces the same bit pattern
            // for every other value.
            low = static_cast<uint64_t>(x);
            if constexpr (std::is_signed_v<T>) {
                // always do positive multiplication
                if (x < 0) {
                    low = 0ull - low;
                    sign = 1;
                }
            }

            if (low == 0)
                return;

            auto expo = log2(low);  // this is the index of the msb as well
            auto shift = static_cast<int32_t>(FRAC_BITS - expo);
            shift_left128_inplace_safe(low, high, shift);
            set_sign(sign);
            set_exponent(static_cast<int32_t>(expo));
            return;
        }
    }
    /**
     * @brief Construct from a string
     * Allows creating very high precision values, approximately 34 decimal digits.
     * Much slower than the other constructors.
     * @param x Input string
     */
    float128(const char* x) noexcept
    {
        low = high = 0;
        if (x == nullptr)
            return;

        // The decimal form is read by from_chars(), which builds the exact value and rounds it
        // once to the nearest representable one. The code below reads a hexadecimal literal, where
        // every digit lands on a bit and no rounding arises.
        //
        // The decimal path used to live here too: it accumulated the digits as a float128 and
        // divided by a power of ten held in the type. Negative powers of ten are not
        // representable in binary, so that division rounded, and the value came back about an ulp
        // away from the one the string named - close enough to look right and far enough that
        // writing it back out took 35 digits instead of one.
        const char* probe = x;
        while (*probe != 0 && isspace(static_cast<unsigned char>(*probe)))
            ++probe;
        const char* body = probe;
        if (*body == '-' || *body == '+')
            ++body;
        if (!(body[0] == '0' && (body[1] == 'x' || body[1] == 'X'))) {
            from_chars(probe, probe + strlen(probe), *this);
            return;
        }

        constexpr uint64_t base16_max_digits = (112 + 4) / 4;  // 29 hex digits. 28 for the fraction (112 bit) and another for the unity
        constexpr uint64_t base10_max_digits = 35;             // maximum for 112 bit of mantissa/fraction is 34, read one extra to get maximum precision
        uint32_t sign = 0;
        uint32_t base = 10;
        int32_t expo2 = 0;   // base2 exponent
        int32_t expo10 = 0;  // base10 exponent

        // convert the input string to lowercase for simpler processing.
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
            *this = float128();
            return;
        }
        // set negative sign if needed
        if (*p == '-') {
            sign = 1;
            ++p;
        } else if (*p == '+')
            ++p;

        // check for infinity
        if (0 == strncmp(p, "inf", 3)) {
            *this = inf();
            return;
        }
        // check for nan
        else if (0 == strncmp(p, "nan", 3)) {
            *this = nan();
            return;
        }

        int32_t max_digits = base10_max_digits;
        // check for a hex string
        if (0 == strncmp(p, "0x", 2)) {
            base = 16;
            p += 2;
            max_digits = base16_max_digits;
        }

        // skip leading zeros
        while (*p == '0')
            ++p;

        int32_t int_digits = 0, frac_digits = 0;
        char* int_start = p;
        char* frac_start = nullptr;

        // count the integer digits
        while (isdigit(static_cast<unsigned char>(*p)) || (base == 16 && *p >= 'a' && *p <= 'f')) {
            ++int_digits;
            ++p;
        }

        // got a hex unsigned int literal
        // every digit is 4 bits, need to keep at most 112 bits after the msb.
        if (base == 16) {
            // zero value
            if (int_digits == 0) {
                set_sign(sign);
                return;
            }

            uint64_t* frac_bits = &low;  // fill the internal data structure directly
            // start at the leftmost digit and iterate right

            int32_t digits_consumed = std::min(int_digits, max_digits);
            char* cur_digit = int_start;
            char* const end_digit = int_start + digits_consumed;

            // fill the internal structure starting at the top bits of high
            while (cur_digit < end_digit) {
                uint64_t d = *cur_digit;
                if (d >= '0' && d <= '9')
                    d -= '0';
                else
                    d = 10ull + d - 'a';

                ++cur_digit;
                auto index = 1 - (expo2 >> 6);  // fill high part first
                assert(index >= 0 && index <= 1);
                auto shift = 60 - (expo2 & 63);
                frac_bits[index] |= d << shift;
                expo2 += 4;
            }
            // shift the result into position and fix the exponent
            uint64_t left_digit = high >> 60;
            assert(left_digit != 0);
            static constexpr int32_t digit_msb_lut[16] = {0, 0, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3};
            const auto digit_msb = digit_msb_lut[left_digit];
            expo2 -= 4 - digit_msb;
            expo2 += 4 * (int_digits - digits_consumed);  // account for digits which do not fit in the fraction.

            // overflow
            if (expo2 > INF_EXP_UNBIASED) {
                *this = inf();
            } else {
                shift_right128_inplace_safe(low, high, 124 + digit_msb - FRAC_BITS);  // move digit's 2nd  msb to bit 111
                set_exponent(expo2);
                set_sign(sign);
            }

            // a hex input value has no exponent or fraction. see note below about exponent support for hex.
            return;
        }

        // Note: fraction and exponent are valid only with base 10 until 'p' style strings are supported (base2 hex exponents). e.g. "1.EDp5F"
        assert(base == 10);

        // check for the optional decimal point
        if (*p == '.') {
            *p = '\0';
            ++p;
            frac_start = (isdigit(static_cast<unsigned char>(*p))) ? p : nullptr;

            // count the fraction digits if they exist
            if (frac_start) {
                while (isdigit(static_cast<unsigned char>(*p))) {
                    ++frac_digits;
                    ++p;
                }

                // back track and erase the trailing zeros
                char* pp = p - 1;
                while (*pp == '0') {
                    --frac_digits;
                    *pp-- = '\0';
                }
            }

            // TODO: optimize small numbers
            // integer part is zero - skip the leading zeros in the fraction and adjust the exponent
            // example 0.01 == 0.1E-1
            // if (int_digits == 0 && frac_start != nullptr) {
            //    while (*frac_start == '0') {
            //        ++frac_start;
            //        --frac_digits;
            //        --expo10_adjust;
            //    }
            //}
        }

        // check for the optional exponent
        if (*p == 'e') {
            *p = '\0';  // terminate the fraction string
            ++p;
            // convert the exponent
            expo10 = strtol(p, nullptr, 10);

            // underflow
            if (expo10 < -4965) {
                return;
            }

            // overflow
            if (expo10 > 4932) {
                *this = inf();
            }
        }

        // compute the integer part
        if (base == 10) {
            uint128_t int_part;
            int32_t digits_consumed = std::min(int_digits, max_digits);
            char* const end_digit = int_start + digits_consumed;
            *end_digit = '\0';
            int_part = int_start;
            int32_t shift_bits = 0;

            // mark the extra exponent that may exist if we had enough bits to represent the entire value
            auto extra_digits = int_digits - digits_consumed;

            // multiply by 10 for each digit
            while (extra_digits > 0) {
                --extra_digits;
                uint64_t msb = int_part ? log2(int_part) : 0;
                if (msb >= 123) {
                    int_part >>= 4;
                    shift_bits += 4;
                    assert(msb - log2(int_part) == 4);
                }
                int_part *= 10;
            }

            int32_t msb = int_part ? static_cast<int32_t>(log2(int_part)) : 0;
            expo2 = msb + shift_bits;

            // shift the integer value into position: msb at bit 112
            int32_t shift = 0;
            shift = FRAC_BITS - msb;
            if (shift > 0)
                int_part <<= shift;
            else
                int_part >>= -shift;

            // set the integer part if non-zero
            if (int_part) {
                int_part.get_components(low, high);  // unity bit is erased by set_exponent()
                set_exponent(expo2);
                set_sign(sign);
            }
        }

        // if both integer & fraction are zero, the result is zero regardless of the exponent e.g. 0E9999
        if (int_digits == 0 && frac_digits == 0) {
            set_sign(sign);
            return;
        }

        // The fraction part is relevant if:
        // 1) Adding more digits actually changes the value
        // 2) Fraction digits exist and not all zero
        float128 frac_part;
        if (frac_digits > 0) {
            // A fraction below 0.1 begins with zeros that carry no information. Skipping them and
            // compensating through the divisor's exponent means the digit budget is spent on
            // significant digits only. Charging the leading zeros against it, as this used to,
            // left a value such as 1.1e-7 with barely 29 significant digits, too few to read back.
            int32_t leading_zeros = 0;
            if (int_digits == 0) {
                while (leading_zeros < frac_digits && frac_start[leading_zeros] == '0')
                    ++leading_zeros;
            }

            // take the minimum of the actual digits in the string versus what is the maximum possible to hold in 112 bit.
            const int32_t digits = std::min(frac_digits - leading_zeros, max_digits - int_digits);

            // Read the digits as a plain integer and scale down with a single division.
            //
            // The accumulation is exact: a float128 holds every integer below 2^113, which is a
            // little above 10^34, so all the digits kept here survive without rounding. Only the
            // final divide rounds.
            //
            // The previous implementation instead multiplied by a stored 10^-9 constant once per
            // group of 9 digits and summed the pieces. Negative powers of ten are not
            // representable in binary, so that constant was already inexact, its running product
            // compounded the error, and every partial sum rounded again. Even "0.5", which is
            // exactly representable, came out one unit in the last place low.
            float128 numerator;
            for (int32_t i = 0; i < digits; ++i) {
                numerator = numerator * 10 + (frac_start[leading_zeros + i] - '0');
            }

            // the skipped zeros are put back through the exponent of the divisor
            frac_part = numerator / exp10(digits + leading_zeros);
            frac_part.set_sign(sign);
        }

        // integer only, no fraction or exponent
        if (frac_digits == 0 && expo10 == 0) {
            return;
        }

        //
        // assemble the integer and fraction to a single value
        //
        *this += frac_part;

        // adjust the result based on the exponent
        if (expo10 != 0) {
            float128 e = exp10(expo10);
            // let the below handle overflow/underflow
            *this *= e;
        }
    }
    /**
     * @brief Constructor from std::string.
     * Allows creating very high precision values. Much slower than the other constructors.
     * @param x Input string
     */
    FP128_INLINE float128(const std::string& x) noexcept
    {
        // delegate to the char* c'tor
        *this = float128(x.c_str());
    }
    /**
     * @brief Destructor
     */
    constexpr ~float128() = default;

    /**
     * @brief Assignment operator
     * @param other Object to copy from
     * @return This object.
     */
    constexpr FP128_FORCE_INLINE float128& operator=(const float128& other) noexcept
    {
        high = other.high;
        low = other.low;
        return *this;
    }
    /**
     * @brief Move assignment operator
     * @param other Object to copy from
     * @return This object.
     */
    constexpr FP128_FORCE_INLINE float128& operator=(float128&& other) noexcept
    {
        high = other.high;
        low = other.low;
        return *this;
    }

    //
    // conversion operators
    //

    /**
     * @brief Operator double
     */
    [[nodiscard]] FP128_INLINE constexpr operator double() const noexcept
    {
        Double d {};

        // nan and inf
        if (get_exponent_bits() == INF_EXP_BIASED) {
            // zero fraction means inf, otherwise nan
            if (low == 0 && get_fraction_bits() == 0) {
                return (get_sign()) ? -HUGE_VAL : HUGE_VAL;
            }
            return NAN;
        }

        int32_t expo = get_exponent();
        // subnormal and underflow
        if (expo < -1022) {
            int32_t shift = (int32_t)(FRAC_BITS - dbl_frac_bits) - 1022 - expo;

            // underflow
            if (shift >= FRAC_BITS) {
                return 0;
            }

            // add the msb back
            uint64_t h = get_fraction_bits() | (1ull << (FRAC_BITS - 64));
            return Double::make(get_sign(), 0, shift_right128_round(low, h, shift));
        }
        // too big for double
        else if (expo > 1023) {
            return HUGE_VAL;
        }
        // normal numbers
        else {
            d.set_e(1023ull + expo);
            d.set_f(shift_right128_round(low, get_fraction_bits(), FRAC_BITS - dbl_frac_bits));
            d.set_s(get_sign());

            // fraction caused a round up
            if (d.f() == 0 && get_fraction_bits() != 0)
                d.set_e(d.e() + 1);

            return d.val();
        }
    }
    /**
     * @brief operator float converts to a float
     */
    [[nodiscard]] FP128_FORCE_INLINE constexpr operator float() const noexcept
    {
        // TODO: write proper function
        double v = static_cast<double>(*this);
        return static_cast<float>(v);
    }
    /**
     * @brief operator uint64_t converts to a uint64_t
     */
    [[nodiscard]] FP128_INLINE constexpr operator uint64_t() const noexcept
    {
        uint64_t l, h;
        int32_t e;
        uint32_t s;
        get_components(l, h, e, s);

        if (e > 63)
            return (s) ? static_cast<uint64_t>(INT64_MIN) : UINT64_MAX;
        if (e < 0)
            return 0;

        auto shift = static_cast<int>(FRAC_BITS) - e;
        shift_right128_inplace_safe(l, h, shift);
        return l;
    }
    /**
     * @brief operator int64_t converts to a int64_t
     */
    [[nodiscard]] FP128_INLINE constexpr operator int64_t() const noexcept
    {
        uint64_t l, h;
        int32_t e;
        uint32_t s;
        get_components(l, h, e, s);

        if (e > 62)
            return (s) ? INT64_MIN : INT64_MAX;
        if (e < 0)
            return 0;

        auto shift = static_cast<int>(FRAC_BITS) - e;
        shift_right128_inplace_safe(l, h, shift);
        int64_t res = static_cast<int64_t>(l);
        return (s) ? -res : res;
    }
    /**
     * @brief operator uint32_t converts to a uint32_t
     */
    [[nodiscard]] FP128_INLINE constexpr operator uint32_t() const noexcept
    {
        uint64_t l, h;
        int32_t e;
        uint32_t s;
        get_components(l, h, e, s);

        if (e > 31)
            return (s) ? static_cast<uint32_t>(INT32_MIN) : UINT32_MAX;
        if (e < 0)
            return 0;

        auto shift = static_cast<int>(FRAC_BITS) - e;
        shift_right128_inplace_safe(l, h, shift);
        return static_cast<uint32_t>(l);
    }
    /**
     * @brief operator int32_t converts to a int32_t
     */
    [[nodiscard]] FP128_INLINE constexpr operator int32_t() const noexcept
    {
        uint64_t l, h;
        int32_t e;
        uint32_t s;
        get_components(l, h, e, s);

        if (e > 30)
            return (s) ? INT32_MIN : INT32_MAX;
        if (e < 0)
            return 0;

        auto shift = static_cast<int>(FRAC_BITS) - e;
        shift_right128_inplace_safe(l, h, shift);
        int32_t res = static_cast<int32_t>(l);
        return (s) ? -res : res;
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
        // The buffer is thread local so that the pointer stays valid after the call returns, which
        // is what an implicit conversion to char* has to provide. It is large enough for the
        // longest shortest-form string: a sign, 36 digits, a point and a five digit exponent.
        static thread_local char str[64];
        const std::string text = detail::render(*this, '\0', -1, false, false, '-');
        const size_t length = (text.size() < sizeof(str) - 1) ? text.size() : sizeof(str) - 1;
        for (size_t i = 0; i < length; ++i)
            str[i] = text[i];
        str[length] = '\0';
        return str;
    }
    /**
     * @brief Writes the value in scientific notation into a caller supplied buffer.
     *
     * Kept for the callers that predate std::format support. Everything it does is also reachable
     * as fp128::to_chars(..., std::chars_format::scientific) or as std::format("{:e}", value),
     * both of which say what they produce more clearly.
     *
     * @param str Output buffer
     * @param buff_size Output buffer size in bytes
     */
    void to_e_format(char* str, int32_t buff_size) const
    {
        if (str == nullptr || buff_size < 1)
            return;

        const std::string text = detail::render(*this, 'e', -1, false, false, '-');
        const size_t length = (text.size() < static_cast<size_t>(buff_size) - 1) ? text.size() : static_cast<size_t>(buff_size) - 1;
        for (size_t i = 0; i < length; ++i)
            str[i] = text[i];
        str[length] = '\0';
    }

    //
    // math operators
    //
    /**
     * @brief Shift right this object.
     * @param shift Bits to shift. Values less than 1 do nothing, high values can cause the value to reach zero.
     * @return This object.
     */
    FP128_INLINE constexpr float128& operator>>=(int32_t shift) noexcept
    {
        if (shift < 1 || is_special())
            return *this;

        uint64_t l, h;
        int32_t e;
        uint32_t s;
        get_components(l, h, e, s);
        e -= shift;
        set_components(l, h, e, s);
        return *this;
    }
    /**
     * @brief Shift left this object.
     * @param shift Bits to shift. Values less than 1 do nothing, high values can cause the value to reach infinity.
     * @return This object.
     */
    FP128_INLINE constexpr float128& operator<<=(int32_t shift) noexcept
    {
        if (shift < 1 || is_special())
            return *this;

        uint64_t l, h;
        int32_t e;
        uint32_t s;
        get_components(l, h, e, s);
        e += shift;
        set_components(l, h, e, s);
        return *this;
    }
    /**
     * @brief Performs right shift operation.
     *
     * Hinted rather than forced, unlike the other one line forwarders (see FP128_FORCE_INLINE).
     * What it forwards to is not cheap: operator>>= splits the value into components and
     * reassembles it, which brings the subnormal handling of get_components() and
     * set_components() along with it. Forcing this pair open cost ~9% of the float128 Mandelbrot
     * benchmark, which shifts inside a loop that already inlines an add, a multiply and two
     * squares. Left as a hint, the shift is still inlined wherever the caller has room for it.
     *
     * @param shift bits to shift
     * @return Temporary object with the result of the operation
     */
    template <typename T> FP128_INLINE constexpr float128 operator>>(T shift) const noexcept
    {
        float128 temp(*this);
        return temp >>= static_cast<int32_t>(shift);
    }
    /**
     * @brief Performs left shift operation.
     * Hinted rather than forced, see operator>> above.
     * @param shift bits to shift
     * @return Temporary object with the result of the operation
     */
    template <typename T> FP128_INLINE constexpr float128 operator<<(T shift) const noexcept
    {
        float128 temp(*this);
        return temp <<= static_cast<int32_t>(shift);
    }

    /**
     * @brief Add a value to this object
     * @param other Right hand side operand
     * @return This object.
     */
    FP128_INLINE constexpr float128& operator+=(const float128& other) noexcept
    {
        // check trivial cases
        if (is_special() || other.is_special()) {
            // A NaN operand propagates. Adding infinities of opposite signs is the invalid
            // operation and produces a NaN as well.
            if (is_nan() || other.is_nan() || (is_inf() && other.is_inf() && get_sign() != other.get_sign())) {
                *this = nan();
                return *this;
            }

            // Either a single operand is infinite, or both are infinite with the same sign.
            // In both cases the result is that infinity, keeping its own sign.
            if (other.is_inf())
                *this = other;
            return *this;
        }

        uint32_t sign, other_sign;
        int32_t expo, other_expo;
        uint64_t l1, h1, l2, h2;
        get_components(l1, h1, expo, sign);
        other.get_components(l2, h2, other_expo, other_sign);
        constexpr int32_t shift_left_bits = 127 - 2 - FRAC_BITS;  // move bit 112 left, keep 1 bit for additional exponent and one for the sign

        if (expo > other_expo) {
            int32_t shift = expo - other_expo - shift_left_bits;  // how many bits to shift right
            // exponents are too far apart, result will stay the same
            if (shift > FRAC_BITS)
                return *this;

            shift_left128_inplace_safe(l1, h1, shift_left_bits);
            if (shift >= 0)
                shift_right128_inplace_safe(l2, h2, shift);
            else
                shift_left128_inplace_safe(l2, h2, -shift);

            // fix the exponent
            expo -= shift_left_bits;
        } else if (expo < other_expo) {
            int32_t shift = other_expo - expo - shift_left_bits;  // how many bits to shift right
            // exponents are too far apart, use the other value
            if (shift > FRAC_BITS) {
                *this = other;
                return *this;
            }
            shift_left128_inplace_safe(l2, h2, shift_left_bits);

            if (shift >= 0)
                shift_right128_inplace_safe(l1, h1, shift);
            else
                shift_left128_inplace_safe(l1, h1, -shift);

            // result base exponent comes from the other value
            expo = other_expo - shift_left_bits;
        }

        // same sign: the simple case
        if (other.get_sign() == get_sign()) {
            // add the other value
            const uint8_t carry = addcarryx_u64(0, l1, l2, &l1);
            addcarryx_u64(carry, h1, h2, &h1);
        }
        // different sign: invert the sign for other and subtract
        else {
            // this value is negative
            if (is_negative())
                twos_complement128(l1, h1);
            // other value is negative
            else
                twos_complement128(l2, h2);

            // add the other value, results stored in l1
            const uint8_t carry = addcarryx_u64(0, l1, l2, &l1);
            addcarryx_u64(carry, h1, h2, &h1);

            // bit 63 is high - got a negative result
            // flip the bits and invert the sign
            sign = FP128_GET_BIT(h1, 63);
            if (sign) {
                twos_complement128(l1, h1);
            }
        }

        norm_fraction(l1, h1, expo);
        set_components(l1, h1, expo, sign);
        return *this;
    }
    /**
     * @brief Add a value to this object
     * @param other Right hand side operand
     * @return This object.
     */
    template <typename T> FP128_FORCE_INLINE constexpr float128& operator+=(const T& other) { return operator+=(float128(other)); }
    /**
     * @brief Subtract a value from this object
     * @param other Right hand side operand
     * @return This object.
     */
    FP128_FORCE_INLINE constexpr float128& operator-=(const float128& other) noexcept { return *this += (-other); }
    /**
     * @brief Subtract a value from this object
     * @param other Right hand side operand
     * @return This object.
     */
    template <typename T> FP128_FORCE_INLINE constexpr float128& operator-=(const T& other) { return operator+=(-float128(other)); }
    /**
     * @brief Multiply a value to this object
     * @param other Right hand side operand
     * @return This object.
     */
    FP128_INLINE constexpr float128& operator*=(const float128& other) noexcept
    {
        // check trivial cases
        if (is_special() || other.is_special()) {
            // A NaN operand propagates, and so does the invalid operation inf * zero.
            // Note the zero tests are only reachable for the operand that is not special.
            if (is_nan() || other.is_nan() || (is_inf() && other.is_zero()) || (is_zero() && other.is_inf())) {
                *this = nan();
                return *this;
            }

            // at least one operand is infinite, the result is infinite with the combined sign
            const uint32_t res_sign = get_sign() ^ other.get_sign();
            *this = inf();
            set_sign(res_sign);
            return *this;
        } else if (is_zero() || other.is_zero()) {
            // the sign of a zero product is the combination of both operand signs
            const uint32_t res_sign = get_sign() ^ other.get_sign();
            *this = 0;
            set_sign(res_sign);
            return *this;
        }
        // extract fractions and exponents
        uint32_t sign, other_sign;
        int32_t expo, other_expo;
        uint64_t l1, h1, l2, h2;
        get_components(l1, h1, expo, sign);
        other.get_components(l2, h2, other_expo, other_sign);
        // Deliberately asked of the objects rather than derived from the mantissas above, even
        // though get_components() has just produced the bits this needs. Testing the raw value
        // leaves the question independent of the extraction, so it can be answered alongside it;
        // phrasing it as (h1 == FRAC_UNITY && l1 == 0) chains it behind instead and costs Clang
        // 15% of the multiply benchmark, against no gain on MSVC.
        bool is_exp2 = is_exponent_of_2();
        bool other_is_exp2 = other.is_exponent_of_2();

        // add the exponents
        expo += other_expo;

        // optimize for exponents of 2
        if (is_exp2 || other_is_exp2) {
            // copy the fraction as needed
            if (is_exp2) {
                l1 = l2;
                h1 = h2;
            }

            set_components(l1, h1, expo, sign ^ other_sign);
            return *this;
        }

        // multiply the fractions
        // the fractions are in u16.112 precision
        // the result will be u32.224 precision and will be shifted-right by 112 bit
        uint64_t res[4];  // 256 bit of result
        uint64_t temp1[2], temp2[2];

        // multiply low QWORDs
        res[0] = mulx_u64(l1, l2, &res[1]);

        // multiply high QWORDs (overflow can happen)
        res[2] = mulx_u64(h1, h2, &res[3]);

        // multiply low this and high other
        temp1[0] = mulx_u64(l1, h2, &temp1[1]);
        uint8_t carry = addcarryx_u64(0, res[1], temp1[0], &res[1]);
        res[3] += addcarryx_u64(carry, res[2], temp1[1], &res[2]);

        // multiply high this and low other
        temp2[0] = mulx_u64(h1, l2, &temp2[1]);
        carry = addcarryx_u64(0, res[1], temp2[0], &res[1]);
        res[3] += addcarryx_u64(carry, res[2], temp2[1], &res[2]);

        // extract the bits from res[] keeping the precision the same as this object
        // shift result by F
        constexpr int32_t index = 1;
        // constexpr int32_t lsb = FRAC_BITS & 63;            // bit within the 64bit data pointed by res[index]
        constexpr int32_t lsb = (FRAC_BITS & 63) - 1;  // bit within the 64bit data pointed by res[index] minus 1 to improve rounding
        // constexpr uint64_t half = 1ull << (lsb - 1);       // used for rounding
        // const bool need_rounding = (res[index] & half) != 0;

        const uint64_t sticky_bits = res[0] | (res[1] & FP128_MAX_VALUE_64(lsb));
        l1 = shift_right128(res[index], res[index + 1], lsb);  // custom function is 20% faster in Mandelbrot than the intrinsic
        h1 = shift_right128(res[index + 1], res[index + 2], lsb);
        --expo;

        // if (need_rounding) {
        //     ++l1; // low will wrap around to zero if overflowed
        //     h1 += l1 == 0;
        // }

        norm_product(l1, h1, expo, sticky_bits != 0);
        set_components(l1, h1, expo, sign ^ other_sign);
        return *this;
    }
    /**
     * @brief Squares this object in place.
     *
     * Cheaper than operator*=(*this) on several counts:
     * - a square is symmetric, so the two cross products of the fraction multiply are the
     *   same value and only three 64x64->128 bit multiplies are needed instead of four
     * - the result sign is always positive, so no sign combining is needed
     * - only one operand has to be tested for the special cases and for being an exponent
     *   of 2, and in the exponent of 2 case the fraction does not have to be copied
     * - the exponents are added to each other, which is a doubling
     *
     * The result is identical to (*this) * (*this) for every value, specials included.
     *
     * @return This object.
     */
    FP128_INLINE constexpr float128& square() noexcept
    {
        // check trivial cases
        if (is_special()) {
            // a NaN stays a NaN, +/-inf squared is +inf
            *this = is_nan() ? nan() : inf();
            return *this;
        } else if (is_zero()) {
            *this = 0;
            return *this;
        }
        // extract the fraction and exponent, the sign of a square is always positive
        uint32_t sign;
        int32_t expo;
        uint64_t l, h;
        get_components(l, h, expo, sign);

        // add the exponent to itself
        expo += expo;

        // optimize for exponents of 2, the fraction is unchanged by the multiply
        if (is_exponent_of_2()) {
            set_components(l, h, expo, 0);
            return *this;
        }

        // square the fraction
        // the fraction is in u16.112 precision
        // the result will be u32.224 precision and will be shifted-right by 112 bit
        uint64_t res[4];  // 256 bit of result
        uint64_t cross[2];

        // multiply the low QWORD by itself
        res[0] = mulx_u64(l, l, &res[1]);

        // multiply the high QWORD by itself (overflow can happen)
        res[2] = mulx_u64(h, h, &res[3]);

        // the low * high cross product, which appears twice in the sum
        cross[0] = mulx_u64(l, h, &cross[1]);

        uint8_t carry = addcarryx_u64(0, res[1], cross[0], &res[1]);
        res[3] += addcarryx_u64(carry, res[2], cross[1], &res[2]);

        carry = addcarryx_u64(0, res[1], cross[0], &res[1]);
        res[3] += addcarryx_u64(carry, res[2], cross[1], &res[2]);

        // extract the bits from res[] keeping the precision the same as this object
        constexpr int32_t index = 1;
        constexpr int32_t lsb = (FRAC_BITS & 63) - 1;  // bit within the 64bit data pointed by res[index] minus 1 to improve rounding

        const uint64_t sticky_bits = res[0] | (res[1] & FP128_MAX_VALUE_64(lsb));
        l = shift_right128(res[index], res[index + 1], lsb);
        h = shift_right128(res[index + 1], res[index + 2], lsb);
        --expo;

        norm_product(l, h, expo, sticky_bits != 0);
        set_components(l, h, expo, 0);
        return *this;
    }
    /**
     * @brief Multiply a value to this object
     * @param other Right hand side operand
     * @return This object.
     */
    template <typename T> FP128_FORCE_INLINE constexpr float128& operator*=(const T& other) { return operator*=(float128(other)); }
    /**
     * @brief Divide this object by a value
     * @param other Right hand side operand
     * @return This object.
     */
    FP128_INLINE float128& operator/=(const float128& other)
    {
        // check trivial cases
        // A NaN operand propagates, and so do the invalid operations zero / zero and inf / inf.
        if (is_nan() || other.is_nan() || (is_zero() && other.is_zero()) || (is_inf() && other.is_inf())) {
            *this = nan();
            return *this;
        }

        // the sign of any of the results below is the combination of both operand signs
        const uint32_t res_sign = get_sign() ^ other.get_sign();

        // An infinite dividend or a zero divisor produce an infinity, a zero dividend or an
        // infinite divisor produce a zero. The combinations where both apply, which are the
        // two invalid operations, were already handled above.
        if (is_inf() || other.is_zero()) {
            *this = inf();
            set_sign(res_sign);
            return *this;
        } else if (is_zero() || other.is_inf()) {
            *this = 0;
            set_sign(res_sign);
            return *this;
        }

        // extract fractions and exponents
        uint32_t sign, other_sign;
        int32_t expo, other_expo;
        uint64_t l1, h1, l2, h2;
        get_components(l1, h1, expo, sign);
        other.get_components(l2, h2, other_expo, other_sign);

        // subtract the exponents
        expo -= other_expo;

        // optimize for other value is an exponent of 2
        if (other.is_exponent_of_2()) {
            set_components(l1, h1, expo, sign ^ other_sign);
            return *this;
        }

        // divide the fractions
        uint64_t q[4] {};
        const uint64_t nom[4] = {0, 0, l1, h1};
        const uint64_t denom[2] = {l2, h2};

        uint64_t rem[2] {};
        if (0 == div_32bit((uint32_t*)q, (uint32_t*)rem, (uint32_t*)nom, (uint32_t*)denom, 2ll * array_length(nom), 2ll * array_length(denom))) {
            // 128 bit were added to the dividend, 112 were lost:
            // need to shift right 16 bit (128 - 112) but we don't go all the way so norm_fraction()
            //  can produce accurate rounding
            constexpr int32_t dshift = 127 - FRAC_BITS;
            // A non zero remainder means the quotient was cut short, and the bits shifted out of
            // q[0] are lost as well. Both have to reach the rounding step as a sticky bit.
            const uint64_t sticky_bits = rem[0] | rem[1] | (q[0] & FP128_MAX_VALUE_64(dshift));
            l1 = shift_right128(q[0], q[1], dshift);
            h1 = shift_right128(q[1], q[2], dshift);
            --expo;
            norm_fraction_sticky(l1, h1, expo, sticky_bits != 0);
        } else {  // error
            *this = inf();
            return *this;
        }

        set_components(l1, h1, expo, sign ^ other_sign);
        return *this;
    }
    /**
     * @brief Divide this object by a value
     * @param other Right hand side operand
     * @return This object.
     */
    template <typename T> FP128_FORCE_INLINE float128& operator/=(const T& other) { return operator/=(float128(other)); }

    //
    // unary operations
    //
    /**
     * @brief Convert to bool
     */
    [[nodiscard]] FP128_FORCE_INLINE constexpr operator bool() const noexcept { return !is_zero(); }
    /**
     * @brief Logical not (!). Opposite of operator bool.
     * Uses is_zero() rather than testing the raw words: negative zero has its sign bit set and
     * would otherwise be reported as a non zero value.
     */
    [[nodiscard]] FP128_FORCE_INLINE constexpr bool operator!() const noexcept { return is_zero(); }
    /**
     * @brief Unary +. Returns a copy of the object.
     */
    [[nodiscard]] FP128_FORCE_INLINE constexpr float128 operator+() const noexcept { return *this; }
    /**
     * @brief Unary -. Returns a copy of the object with sign inverted.
     */
    [[nodiscard]] FP128_FORCE_INLINE constexpr float128 operator-() const noexcept { return float128(low, high ^ SIGN_MASK); }

    //
    // useful public functions
    //
    /**
     * @brief Returns true if the value is positive (including zero and NaN)
     * @return True when the sign is 0
     */
    [[nodiscard]] FP128_FORCE_INLINE constexpr bool is_positive() const noexcept { return 0 == (high & SIGN_MASK); }
    /**
     * @brief Returns true if the value is negative (including zero and NaN).
     * @return True when the sign is 1
     */
    [[nodiscard]] FP128_FORCE_INLINE constexpr bool is_negative() const noexcept { return 0 != (high & SIGN_MASK); }
    /**
     * @brief Returns true if and only if the value is ±0.
     * @return Returns true if the value is zero
     */
    [[nodiscard]] FP128_FORCE_INLINE constexpr bool is_zero() const noexcept { return 0 == low && 0 == (high & ~SIGN_MASK); }
    /**
     * @brief Returns true if and only if x is zero, subnormal or normal (not infinite or NaN).
     * @return True if and only if x is zero, subnormal or normal (not infinite or NaN).
     */
    [[nodiscard]] FP128_FORCE_INLINE constexpr bool is_finite() const noexcept { return !is_special(); }
    /**
     * @brief Tests if the value is subnormal
     * @return True when the value is subnormal
     */
    [[nodiscard]] FP128_FORCE_INLINE constexpr bool is_subnormal() const noexcept { return get_exponent_bits() == 0; }
    /**
     * @brief Tests if the value is normal (not zero, subnormal, infinite, or NaN)
     * @return True if and only if the value is normal
     */
    [[nodiscard]] FP128_FORCE_INLINE constexpr bool is_normal() const noexcept { return get_exponent_bits() != 0 && get_exponent_bits() != INF_EXP_BIASED; }
    /**
     * @brief Tests if this value is a NaN
     * @return True when the value is a NaN
     */
    [[nodiscard]] FP128_FORCE_INLINE constexpr bool is_nan() const
    {
        // fraction is zero for +- INF, non-zero for NaN
        return get_exponent_bits() == INF_EXP_BIASED && (get_fraction_bits() != 0 || low != 0);
    }
    /**
     * @brief Tests if this value is a signaling NaN
     * @return True if this value is a signaling NaN
     */
    [[nodiscard]] FP128_FORCE_INLINE constexpr bool is_signaling() const
    {
        // A NaN is quiet when the leading fraction bit is set and signaling when it is clear,
        // which is what every binary interchange format in IEEE 754-2008 uses.
        return is_nan() && (high & QUIET_NAN_BIT) == 0;
    }
    /**
     * @brief Tests if this value is an Infinite (negative or positive)
     * @return True when the value is an Infinite
     */
    [[nodiscard]] FP128_FORCE_INLINE constexpr bool is_inf() const
    {
        // fraction is zero for +- INF, non-zero for NaN
        return get_exponent_bits() == INF_EXP_BIASED && get_fraction_bits() == 0;
    }
    /**
     * @brief Tests if the value is an exponent of 2 (fraction part is zero)
     * @return True when the value is an exponent of 2
     */
    [[nodiscard]] FP128_FORCE_INLINE constexpr bool is_exponent_of_2() const
    {
        // fraction is zero for +- INF, non-zero for NaN
        return get_fraction_bits() == 0 && low == 0;
    }
    /**
     * @brief return true when the value is either an inf or nan
     * @return true for inf and nan
     */
    [[nodiscard]] FP128_FORCE_INLINE constexpr bool is_special() const
    {
        // fraction is zero for +- INF, non-zero for NaN
        return get_exponent_bits() == INF_EXP_BIASED;
    }
    /**
     * @brief Returns if the value is an integer (fraction is zero).
     * @return True when the value is an integer.
     */
    [[nodiscard]] FP128_FORCE_INLINE constexpr bool is_int() const
    {
        int32_t expo = get_exponent();
        if (expo < 0)
            return false;
        if (expo >= FRAC_BITS)
            return true;
        return get_fraction().is_zero();
    }
    /**
     * @brief get a specific bit within the float128 data
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
     * @brief Return the fraction part as a float128
     */
    [[nodiscard]] FP128_INLINE constexpr float128 get_fraction() const
    {
        auto expo = get_exponent();
        int32_t frac_bits = static_cast<int32_t>(FRAC_BITS) - expo;
        // all the bits are fraction
        if (frac_bits > FRAC_BITS)
            return *this;
        // the exponent is too large to hold a fraction
        if (frac_bits <= 0)
            return 0;

        uint64_t l = low, h = get_fraction_bits();
        if (frac_bits <= 64) {
            h = 0;
            l &= FP128_MAX_VALUE_64(frac_bits);
        } else {
            h &= FP128_MAX_VALUE_64(frac_bits - 64);
        }

        // no fraction bits are high
        if (l == 0 && h == 0) {
            return 0;
        }

        // find the msb and shift to bit 112
        int32_t msb = static_cast<int32_t>(log2(l, h));
        int32_t shift = FRAC_BITS - msb;  // how many bits to shift left
        shift_left128_inplace_safe(l, h, shift);
        expo -= shift;
        float128 res(l, h, expo + EXP_BIAS, get_sign());
        return res;
    }
    /**
     * @brief Reads the biased exponent field out of the high QWORD.
     *
     * This accessor and the three below it are the only places that know where the fields sit
     * within the high QWORD, which is the same layout the (l, h, e, s) constructor assembles.
     * They deliberately shift and mask rather than going through the _float128_bits view: that
     * would be a read of the inactive member of a union, which constant evaluation rejects.
     *
     * @return The biased exponent, in [0, 0x7FFF].
     */
    [[nodiscard]] FP128_FORCE_INLINE constexpr uint64_t get_exponent_bits() const noexcept { return (high >> EXP_SHIFT) & EXP_MASK; }
    /**
     * @brief Sets the biased exponent field, leaving the fraction and the sign alone.
     * @param e Biased exponent. Bits above the 15 the field holds are dropped.
     */
    FP128_FORCE_INLINE constexpr void set_exponent_bits(uint64_t e) noexcept { high = (high & ~(EXP_MASK << EXP_SHIFT)) | ((e & EXP_MASK) << EXP_SHIFT); }
    /**
     * @brief Reads the upper 48 bits of the fraction out of the high QWORD.
     * @return Bits [111:64] of the fraction.
     */
    [[nodiscard]] FP128_FORCE_INLINE constexpr uint64_t get_fraction_bits() const noexcept { return high & UPPER_FRAC_MASK; }
    /**
     * @brief Sets the upper 48 bits of the fraction, leaving the exponent and the sign alone.
     * @param f Fraction bits. Bits above the 48 the field holds are dropped.
     */
    FP128_FORCE_INLINE constexpr void set_fraction_bits(uint64_t f) noexcept { high = (high & ~UPPER_FRAC_MASK) | (f & UPPER_FRAC_MASK); }
    /**
     * @brief Inverts the sign
     */
    FP128_FORCE_INLINE constexpr void invert_sign() noexcept { high ^= SIGN_MASK; }
    /**
     * @brief Sets the sign
     */
    FP128_FORCE_INLINE constexpr void set_sign(uint64_t s) noexcept { high = (high & ~SIGN_MASK) | ((s & 1) << 63); }
    /**
     * @brief Gets the sign
     */
    [[nodiscard]] FP128_FORCE_INLINE constexpr uint32_t get_sign() const noexcept { return static_cast<uint32_t>(high >> 63); }
    /**
     * @brief Reads the raw binary128 encoding.
     *
     * The counterpart of the float128(low, high) constructor: the two together are what a
     * std::bit_cast to and from the storage would be, without needing the layout to be spelled
     * out at the call site. Nothing is interpreted, so a NaN keeps its payload and a subnormal is
     * not normalized - use get_components() for the value the encoding stands for.
     *
     * @param l Receives the low QWORD (fraction bits 63:0)
     * @param h Receives the high QWORD (sign, exponent, fraction bits 111:64)
     */
    FP128_FORCE_INLINE constexpr void get_bits(uint64_t& l, uint64_t& h) const noexcept
    {
        l = low;
        h = high;
    }
    /**
     * @brief Classifies the float128 value according to IEEE 754.
     * @return The classification category (normal, subnormal, zero, inf, or NaN).
     */
    [[nodiscard]] FP128_INLINE constexpr float128_class_t get_class() const noexcept
    {
        // inf and Nan
        if (get_exponent_bits() == INF_EXP_BIASED) {
            if (is_nan())
                // TODO: support signalling NaN
                return quietNaN;
            return is_positive() ? positiveInfinity : negativeInfinity;
        }
        if (is_zero()) {
            return is_positive() ? positiveZero : negativeZero;
        }
        if (is_subnormal()) {
            return is_positive() ? positiveSubnormal : negativeSubnormal;
        }

        return is_positive() ? positiveNormal : negativeNormal;
    }
    /**
     * @brief Returns the exponent of the object - like the base 2 exponent of a floating point
     * A value of 2.1 would return 1, values in the range [0.5,1.0) would return -1.
     * @return Exponent of the number
     */
    [[nodiscard]] FP128_FORCE_INLINE constexpr int32_t get_exponent() const noexcept { return static_cast<int32_t>(get_exponent_bits()) - EXP_BIAS; }
    /**
     * @brief Set the exponent
     * @param e Exponent value
     */
    FP128_FORCE_INLINE constexpr void set_exponent(int32_t e) noexcept
    {
        e += EXP_BIAS;
        assert(e >= 0);
        assert(e <= INF_EXP_BIASED);
        set_exponent_bits(static_cast<uint64_t>(e));
    }
    /**
     * @brief break the float into its components.
     * Normalizes subnormal values
     * @param l Reference to receive the low fraction
     * @param h Reference to receive the high fraction
     * @param e Reference to receive the unbiased exponent
     * @param s Reference to receive the sign
     */
    FP128_INLINE constexpr void get_components(uint64_t& l, uint64_t& h, int32_t& e, uint32_t& s) const noexcept
    {
        l = low;
        h = get_fraction_bits();
        e = get_exponent();
        s = get_sign();

        if (is_subnormal()) {
            // shift the bits to the left so the msb is on bit 112
            auto shift = FRAC_BITS - static_cast<int32_t>(log2(l, h));
            shift_left128_inplace_safe(l, h, shift);

            // A subnormal has no implicit leading one, so its stored exponent field says nothing
            // about its magnitude - every subnormal holds the same one. What the value is worth is
            // decided entirely by where its leading bit sits: the fraction is scaled by 2^-16494,
            // and normalizing it left by `shift` puts the exponent that many places below the
            // smallest normal one. Adding the shift to the stored exponent instead, as this used
            // to, put the smallest subnormal 223 binades above where it belongs and made every
            // operation on a subnormal operand return an unrelated value. set_components() has
            // always used this convention, which is why the two did not round trip.
            e = SUBNORM_EXP_UNBIASED + 1 - shift;
        }
        // normal numbers
        else if (is_normal()) {
            // add the unity value
            h |= FRAC_UNITY;
        }
    }
    /**
     * @brief Sets the internal components handling al special cases
     * The fraction is expected to have 113 bits (includes the unity bit)
     * It will store the float128 value taking into account infinity ands subnormals.
     * @param l Low part of the fraction
     * @param h High part of the fraction
     * @param e Unbiased exponent, can be any value.
     * @param s Sign (1 is negative)
     */
    /**
     * @brief Number of Mercator series terms log2() needs.
     *
     * The reduction leaves |z| <= 2^-6, so term n is bounded by 2^(-6n), and the series is cut off
     * once that is below the last bit of a 113 bit mantissa with eight bits to spare.
     */
    static constexpr int32_t LOG2_TERMS = (FRAC_BITS + 1 + 8 + log2_reduction_bits - 1) / log2_reduction_bits;

    /**
     * @brief Builds a float128 from a 128 bit fraction, a value in [0,1).
     *
     * The counterpart of reading a mantissa out with get_components(): log2() does its argument
     * reduction and series on plain 128 bit fractions, and this is how the two results it ends up
     * with re-enter floating point.
     *
     * @param low Low QWORD of the fraction.
     * @param high High QWORD of the fraction.
     * @return The value the fraction stands for, or zero if it is zero.
     */
    [[nodiscard]] static FP128_INLINE float128 from_fraction128(uint64_t low, uint64_t high) noexcept
    {
        if ((low | high) == 0) {
            return float128();
        }

        // Bring the leading one to bit 127, then keep the top 113 bits as the mantissa.
        const int32_t leading_zeros = static_cast<int32_t>(lzcnt128(low, high));
        shift_left128_inplace_safe(low, high, leading_zeros);
        shift_right128_inplace_safe(low, high, 127 - FRAC_BITS);
        int32_t expo = -(leading_zeros + 1);

        // That last shift rounds, and rounding up can carry into a 114th bit.
        if (high > (FRAC_UNITY | UPPER_FRAC_MASK)) {
            shift_right128_inplace_safe(low, high, 1);
            ++expo;
        }

        float128 res;
        res.set_components(low, high, expo, 0);

        return res;
    }

    FP128_INLINE constexpr void set_components(uint64_t l, uint64_t h, int32_t e, uint32_t s) noexcept
    {
        // overflow
        if (e >= INF_EXP_UNBIASED) {
            *this = inf();
            set_sign(s);
            return;
        }
        // sub normals
        if (e <= SUBNORM_EXP_UNBIASED) {
            // fix the fraction, remove the bits 112+ and shift to the right
            int32_t shift = SUBNORM_EXP_UNBIASED + 1 - e;
            // The shift is applied unconditionally. Zeroing the fraction once it reached the
            // mantissa's width, as this used to, threw away the smallest subnormal itself: its
            // leading bit sits at bit 112 and a shift of exactly 112 lands it on bit 0 rather than
            // off the end. The shift helper handles every width, and rounds, which is what decides
            // whether a value just under half of the smallest subnormal survives as one.
            shift_right128_inplace_safe(l, h, shift);
            // override the exponent to mark the value as subnormal
            e = SUBNORM_EXP_UNBIASED;
        }

        low = l;
        set_fraction_bits(h);
        e += EXP_BIAS;
        set_exponent_bits(static_cast<uint64_t>(e));
        set_sign(s != 0);
    }
    /**
     * @brief Normalize the fraction so the msb (unity bit) is on bit 112.
     * The fraction value must contain the unity value
     * The exponent component is adjusted based on the bit shift direction required.
     * @param l Low part of the fraction
     * @param h High part of the fraction
     * @param e Unbiased exponent, can be any value.
     */
    /// @brief Word count of the fixed point expansion used to print the fraction digits.
    static constexpr int32_t FIXED_WORDS = 4;

    /**
     * @brief Expands the fraction part of this value into an exact fixed point number.
     *
     * The result is the fraction scaled by 2^256 and held in four words, frac[0] being the least
     * significant. Nothing is lost: a fraction bit of a binary128 has weight 2^(e-112) with e no
     * smaller than -111 on the path that calls this, so scaling by 2^256 always lands on an
     * integer, and the fraction is below one so the product stays under 2^256.
     *
     * Printing the digits from this expansion is exact. Deriving them by repeatedly multiplying a
     * float128 by 100000 instead, as this used to, rounds once per group and the error compounds
     * over the eight or so groups that a full precision value needs.
     *
     * @param frac Output, four words holding the fraction scaled by 2^256
     * @return True when a fraction exists, false when the value is an integer
     */
    [[nodiscard]] FP128_INLINE constexpr bool fraction_to_fixed(uint64_t frac[FIXED_WORDS]) const noexcept
    {
        for (int32_t i = 0; i < FIXED_WORDS; ++i)
            frac[i] = 0;

        uint64_t l, h;
        int32_t e;
        uint32_t s;
        get_components(l, h, e, s);

        // no bits sit below the binary point
        if (e >= FRAC_BITS)
            return false;

        // keep only the bits whose weight is below one. A value under one contributes all of them.
        if (e >= 0) {
            const int32_t keep = FRAC_BITS - e;  // 1 - 112
            if (keep <= 64) {
                l &= FP128_MAX_VALUE_64(keep);
                h = 0;
            } else {
                h &= FP128_MAX_VALUE_64(keep - 64);
            }
        }

        if (l == 0 && h == 0)
            return false;

        // bit i of the surviving fraction has weight 2^(e-112+i), which scaled by 2^256 lands on
        // bit e+144+i. The shift is non negative for every exponent that reaches this function.
        const int32_t shift = e + 2 * 64 + (FIXED_WORDS * 64 - 2 * 64) - FRAC_BITS;  // e + 144
        FP128_ASSERT(shift >= 0 && shift < FIXED_WORDS * 64);
        const int32_t word = shift >> 6;
        const int32_t bit = shift & 63;

        if (bit == 0) {
            frac[word] = l;
            if (word + 1 < FIXED_WORDS)
                frac[word + 1] = h;
        } else {
            frac[word] = l << bit;
            if (word + 1 < FIXED_WORDS)
                frac[word + 1] = (l >> (64 - bit)) | (h << bit);
            if (word + 2 < FIXED_WORDS)
                frac[word + 2] = h >> (64 - bit);
        }
        return true;
    }
    /**
     * @brief Multiplies the fixed point fraction by a small value and returns the part that
     *        crossed the binary point.
     * @param frac Fraction scaled by 2^256, updated in place to hold the new fraction
     * @param mul Multiplier, small enough that the result stays below 2^64
     * @return The integer part produced by the multiply, which is the next group of digits
     */
    [[nodiscard]] FP128_INLINE static constexpr uint64_t fixed_mul_extract(uint64_t frac[FIXED_WORDS], uint64_t mul) noexcept
    {
        uint64_t carry = 0;
        for (int32_t i = 0; i < FIXED_WORDS; ++i) {
            uint64_t hi;
            const uint64_t lo = mulx_u64(frac[i], mul, &hi);
            const uint8_t c = addcarryx_u64(0, lo, carry, &frac[i]);
            carry = hi + c;
        }
        // the fraction was below one, so the overflow is below mul and fits in a single word
        return carry;
    }
    /// @brief True when every word of the fixed point fraction is zero.
    [[nodiscard]] FP128_INLINE static constexpr bool fixed_is_zero(const uint64_t frac[FIXED_WORDS]) noexcept
    {
        uint64_t acc = 0;
        for (int32_t i = 0; i < FIXED_WORDS; ++i)
            acc |= frac[i];
        return acc == 0;
    }

    /**
     * @name Wide accumulator
     *
     * A fixed 384 bit unsigned integer, used by fma() to hold the exact sum of a 226 bit product
     * and a 113 bit addend before a single rounding is applied to it. The width is what the two
     * together need in the worst case: whichever of them is larger is placed with its leading bit
     * at the top of the accumulator, which leaves 383-225 = 158 bits below the product and
     * 383-112 = 271 bits below the addend. An operand that does not fit under the other one is
     * more than 158 bits smaller than it, so it cannot influence anything but the sticky bit, and
     * is folded into that instead of being placed.
     *
     * These are implementation details of fma() rather than part of the interface.
     * @{
     */
    static constexpr int32_t WIDE_WORDS = 6;               ///< Words in the accumulator.
    static constexpr int32_t WIDE_BITS = WIDE_WORDS * 64;  ///< Bits in the accumulator.
    /// @brief Bit the leading one of the larger operand is placed on.
    ///
    /// One below the top, so that adding two operands of the same magnitude cannot carry out of
    /// the accumulator.
    static constexpr int32_t WIDE_TOP = WIDE_BITS - 2;

    /**
     * @brief Full 256 bit product of two 128 bit values.
     * @param al Low QWORD of the first operand
     * @param ah High QWORD of the first operand
     * @param bl Low QWORD of the second operand
     * @param bh High QWORD of the second operand
     * @param res Output, four words with res[0] the least significant
     */
    FP128_INLINE static constexpr void mul128to256(uint64_t al, uint64_t ah, uint64_t bl, uint64_t bh, uint64_t res[4]) noexcept
    {
        uint64_t hi = 0;
        res[0] = mulx_u64(al, bl, &hi);
        res[1] = hi;
        res[2] = 0;
        res[3] = 0;

        uint64_t lo = mulx_u64(al, bh, &hi);
        unsigned char c = addcarryx_u64(0, res[1], lo, &res[1]);
        c = addcarryx_u64(c, res[2], hi, &res[2]);
        res[3] += c;

        lo = mulx_u64(ah, bl, &hi);
        c = addcarryx_u64(0, res[1], lo, &res[1]);
        c = addcarryx_u64(c, res[2], hi, &res[2]);
        res[3] += c;

        lo = mulx_u64(ah, bh, &hi);
        c = addcarryx_u64(0, res[2], lo, &res[2]);
        res[3] += hi + c;
    }

    /**
     * @brief Index of the highest set bit of a multi word value, or -1 when it is zero.
     * @param a Words, a[0] the least significant
     * @param words Word count
     */
    [[nodiscard]] FP128_INLINE static constexpr int32_t wide_msb(const uint64_t* a, int32_t words) noexcept
    {
        for (int32_t i = words - 1; i >= 0; --i) {
            if (a[i] != 0)
                return i * 64 + 63 - static_cast<int32_t>(lzcnt64(a[i]));
        }
        return -1;
    }

    /**
     * @brief Writes a value into a zeroed accumulator, shifted to a given bit position.
     *
     * A negative offset puts part or all of the value below the accumulator's least significant
     * bit. Those bits are not representable there but still decide the rounding, so whether any of
     * them was set is reported back rather than dropped.
     *
     * @param acc Accumulator, must be zero on entry
     * @param src Value to place, src[0] the least significant word
     * @param words Word count of src
     * @param offset Bit position the least significant bit of src lands on, may be negative
     * @return True when a set bit fell below the accumulator.
     */
    [[nodiscard]] FP128_INLINE static constexpr bool wide_place(uint64_t* acc, const uint64_t* src, int32_t words, int32_t offset) noexcept
    {
        if (offset >= 0) {
            const int32_t word = offset >> 6;
            const int32_t bit = offset & 63;
            for (int32_t i = 0; i < words; ++i) {
                const int32_t dst = i + word;
                if (dst >= WIDE_WORDS)
                    break;
                acc[dst] |= src[i] << bit;
                if (bit != 0 && dst + 1 < WIDE_WORDS)
                    acc[dst + 1] |= src[i] >> (64 - bit);
            }
            return false;
        }

        const int32_t drop = -offset;
        if (drop >= words * 64) {
            uint64_t any = 0;
            for (int32_t i = 0; i < words; ++i)
                any |= src[i];
            return any != 0;
        }

        const int32_t word = drop >> 6;
        const int32_t bit = drop & 63;
        uint64_t lost = 0;
        for (int32_t i = 0; i < word; ++i)
            lost |= src[i];
        if (bit != 0)
            lost |= src[word] & FP128_MAX_VALUE_64(bit);

        for (int32_t i = 0; word + i < words; ++i) {
            uint64_t value = src[word + i] >> bit;
            if (bit != 0 && word + i + 1 < words)
                value |= src[word + i + 1] << (64 - bit);
            if (i < WIDE_WORDS)
                acc[i] |= value;
        }
        return lost != 0;
    }

    /// @brief Adds b into a. The accumulator is wide enough that the sum cannot carry out.
    FP128_INLINE static constexpr void wide_add(uint64_t* a, const uint64_t* b) noexcept
    {
        unsigned char c = 0;
        for (int32_t i = 0; i < WIDE_WORDS; ++i)
            c = addcarryx_u64(c, a[i], b[i], &a[i]);
    }

    /// @brief Subtracts b from a, which must not be smaller.
    FP128_INLINE static constexpr void wide_sub(uint64_t* a, const uint64_t* b) noexcept
    {
        unsigned char borrow = 0;
        for (int32_t i = 0; i < WIDE_WORDS; ++i) {
            const uint64_t rhs = b[i];
            const uint64_t diff = a[i] - rhs - borrow;
            borrow = (a[i] < rhs || (borrow && a[i] == rhs)) ? 1 : 0;
            a[i] = diff;
        }
    }

    /// @brief Subtracts one from the accumulator, which must not be zero.
    FP128_INLINE static constexpr void wide_dec(uint64_t* a) noexcept
    {
        for (int32_t i = 0; i < WIDE_WORDS; ++i) {
            if (a[i]-- != 0)
                break;
        }
    }

    /// @brief Three way comparison of two accumulators.
    [[nodiscard]] FP128_INLINE static constexpr int32_t wide_cmp(const uint64_t* a, const uint64_t* b) noexcept
    {
        for (int32_t i = WIDE_WORDS - 1; i >= 0; --i) {
            if (a[i] != b[i])
                return (a[i] > b[i]) ? 1 : -1;
        }
        return 0;
    }
    /// @}
    /**
     * @brief Normalize the fraction to bit 112, rounding half to even.
     *
     * The multiply and the divide both produce more bits than the fraction can hold and then throw
     * the surplus away. Rounding correctly needs two pieces of information about what was
     * discarded: the highest dropped bit, which chooses the direction, and whether anything below
     * it was set, which separates an exact tie from a remainder above half. The caller passes the
     * latter for the words it dropped before calling; this function collects the rest.
     *
     * norm_fraction() below rounds from a three bit window instead, which cannot see the discarded
     * low words at all, so its ties resolve arbitrarily.
     *
     * @param l Low part of the fraction
     * @param h High part of the fraction
     * @param e Unbiased exponent, adjusted to match the normalized fraction
     * @param sticky True when the caller already dropped one or more set bits below l
     */
    /**
     * @brief Normalizes the product of two mantissas to bit 112, rounding half to even.
     *
     * The multiply and the square both arrive here with (h:l) equal to their 256 bit product shifted
     * right by 111. Both operands were normalized to [2^112, 2^113) by get_components(), so the
     * product is in [2^224, 2^226) and (h:l) is in [2^113, 2^115): its leading one is at bit 113 or
     * bit 114 and nowhere else, which makes the remaining shift 1 or 2.
     *
     * norm_fraction_sticky() would reach the same answer, but it has to find the leading one with a
     * count of leading zeros first and then handle a shift that could be anything from negative to
     * past the end of a QWORD. Here the choice is a single bit test and the two shift widths are
     * small enough that none of its range handling applies. The rounding is deliberately identical
     * to it, down to the carry out case, so the two produce the same bits for every input.
     *
     * @param l Low part of the shifted product, replaced by the normalized fraction
     * @param h High part of the shifted product, replaced by the normalized fraction
     * @param e Unbiased exponent, adjusted to match the normalized fraction
     * @param sticky True when the caller already dropped one or more set bits below l
     */
    FP128_INLINE static constexpr void norm_product(uint64_t& l, uint64_t& h, int32_t& e, bool sticky) noexcept
    {
        // bit 114 of the product decides between the two, and it is bit 50 of the high QWORD
        const int32_t shift = 1 + static_cast<int32_t>(h >> 50);
        e += shift;

        // the highest dropped bit decides the direction, everything under it is sticky
        const uint64_t guard = (l >> (shift - 1)) & 1;
        const bool below = sticky || (shift == 2 && (l & 1) != 0);

        l = shift_right128(l, h, shift);
        h >>= shift;

        if (guard && (below || (l & 1))) {
            if (++l == 0)
                ++h;
            // a carry out of the fraction moves to the next power of two
            if ((h >> 48) != 1)
                ++e;
        }
    }
    FP128_INLINE constexpr void norm_fraction_sticky(uint64_t& l, uint64_t& h, int32_t& e, bool sticky) const noexcept
    {
        if (l == 0 && h == 0) {
            e = ZERO_EXP_BIASED;
            return;
        }
        const int32_t msb = static_cast<int32_t>(log2(l, h));
        const int32_t shift = msb - FRAC_BITS;
        e += shift;
        if (shift <= 0) {
            shift_left128_inplace_safe(l, h, -shift);
            return;
        }

        // the highest dropped bit decides the direction, everything under it is sticky
        const uint64_t guard = (shift <= 64) ? FP128_GET_BIT(l, shift - 1) : FP128_GET_BIT(h, shift - 65);
        bool below = sticky;
        if (!below && shift >= 2) {
            below = (shift <= 64) ? ((l << (65 - shift)) != 0)
                                  : (l != 0 || (shift >= 66 && (h << (129 - shift)) != 0));
        }

        const uint64_t new_l = shift_right128(l, h, shift);
        h = (shift < 64) ? (h >> shift) : 0;
        l = new_l;

        if (guard && (below || (l & 1))) {
            if (++l == 0)
                ++h;
            // a carry out of the fraction moves to the next power of two
            if ((h >> 48) != 1)
                ++e;
        }
    }
    FP128_INLINE constexpr void norm_fraction(uint64_t& l, uint64_t& h, int32_t& e) const noexcept
    {
        // l and h are both zero
        if (l == 0 && h == 0) {
            e = ZERO_EXP_BIASED;
            return;
        }

        // fix the exponent
        auto msb = static_cast<int32_t>(log2(l, h));

        // if the msb is exactly msb == FRAC_BITS the exponent stays the same
        auto shift = msb - FRAC_BITS;

        e += shift;
        if (shift > 0) {
            shift_right128_inplace_safe(l, h, shift);
            // rounding up may have happened, expect the upper 16 bit to be exactly 1
            if ((h >> 48) != 1) {
                ++e;
            }
        } else {
            // assert(shift == 0);
            shift_left128_inplace_safe(l, h, -shift);
        }
    }
    /**
     * @brief Produces the closest value larger than x
     * nextUp(x) is the least floating-point number in the format of x that compares greater than x.
     * If x is the negative number of least magnitude in x’s format, nextUp(x) is −0.
     * nextUp(±0) is the positive number of least magnitude in x’s format.
     * nextUp(+∞) is +∞, and nextUp(−∞) is the finite negative number largest in magnitude.
     * When x is NaN, then the result is according to 6.2. nextUp(x) is quiet except for sNaNs.
     * @param x Source value
     * @return Higher value closest to x
     */
    [[nodiscard]] FP128_INLINE static constexpr float128 nextUp(float128 x)
    {
        switch (x.get_class()) {
        case positiveInfinity:
        case quietNaN:
        case signalingNaN:
            break;
        case negativeInfinity:
            return float128(UINT64_MAX, UINT64_MAX, EXP_MASK - 1, 1);
        case negativeZero:
        case positiveZero:
            x.low = 1;
            x.set_sign(0);
            break;
        case positiveSubnormal:
        case positiveNormal:
            ++x.low;
            if (x.low == 0) {
                ++x.high;  // takes care of raising the exponent if needed as well
            }
            break;
        case negativeSubnormal:
        case negativeNormal:
            --x.low;
            if (x.low == UINT64_MAX) {
                --x.high;  // takes care of decreasing the exponent if needed as well
            }
            break;
        }
        return x;
    }
    /**
     * @brief Produces the closest value smaller than x
     * @param x Source value
     * @return Lower value closest to x
     */
    [[nodiscard]] FP128_INLINE static constexpr float128 nextDown(float128 x)
    {
        switch (x.get_class()) {
        case negativeInfinity:
        case quietNaN:
        case signalingNaN:
            break;
        case positiveInfinity:
            return float128(UINT64_MAX, UINT64_MAX, EXP_MASK - 1, 0);
        case negativeZero:
        case positiveZero:
            x.low = 1;
            x.set_sign(1);
            break;
        case negativeSubnormal:
        case negativeNormal:
            ++x.low;
            if (x.low == 0) {
                ++x.high;  // takes care of raising the exponent if needed as well
            }
            break;
        case positiveSubnormal:
        case positiveNormal:
            --x.low;
            if (x.low == UINT64_MAX) {
                --x.high;  // takes care of decreasing the exponent if needed as well
            }
            break;
        }
        return x;
    }

    /**
     * @brief Return the infinite constant
     * @return INF
     */
    [[nodiscard]] FP128_FORCE_INLINE static constexpr float128 inf() { return float128(0, 0, INF_EXP_BIASED, 0); }
    /**
     * @brief Return the quiet (non-signaling) NaN constant
     * @return NaN
     */
    [[nodiscard]] FP128_FORCE_INLINE static constexpr float128 nan() { return float128(0, QUIET_NAN_BIT, INF_EXP_BIASED, 0); }
    /**
     * @brief Return a signaling NaN
     *
     * Nothing in the library raises an exception on encountering one - there is no floating point
     * status to raise it in - but the encoding is distinguishable, which is what
     * numeric_limits::signaling_NaN() and is_signaling() need.
     * @return A signaling NaN
     */
    [[nodiscard]] FP128_FORCE_INLINE static constexpr float128 signaling_nan() { return float128(1, 0, INF_EXP_BIASED, 0); }
    /**
     * @brief Return the value of pi
     * @return pi
     */
    [[nodiscard]] FP128_FORCE_INLINE static constexpr float128 pi() { return float128(0x8469898CC51701B8, 0x921FB54442D1, 0x4000, 0); }
    /**
     * @brief Return the value of pi / 2
     * @return pi / 2
     */
    [[nodiscard]] FP128_FORCE_INLINE static constexpr float128 half_pi() { return float128(0x8469898CC51701B8, 0x921FB54442D1, 0x3FFF, 0); }
    /**
     * @brief Return the value of pi / 4
     * @return pi / 4
     */
    [[nodiscard]] FP128_FORCE_INLINE static constexpr float128 quarter_pi() { return float128(0x8469898CC51701B8, 0x921FB54442D1, 0x3FFE, 0); }
    /**
     * @brief Return the value of e
     * @return e
     */
    [[nodiscard]] FP128_FORCE_INLINE static constexpr float128 e()
    {
        // static const float128 e = "2.71828182845904523536028747135266249775724709369"; // 50 first digits of e
        // return e;
        return float128(0x95355FB8AC404E7A, 0x5BF0A8B14576, 0x4000, 0);
    }
    /**
     * @brief Returns a value of sqrt(2)
     * @return
     */
    [[nodiscard]] FP128_FORCE_INLINE static constexpr float128 sqrt_2() noexcept
    {
        // static const float128 sqrt_2 = "1.41421356237309504880168872420969807856967187537"; // 50 first digits of sqrt(2)
        // return sqrt_2;
        return float128(0xC908B2FB1366EA95, 0x00006A09E667F3BC, 0x3FFF, 0);
    }
    /**
     * @brief  Returns a value of 1
     * @return 1
     */
    [[nodiscard]] FP128_FORCE_INLINE static constexpr float128 one() noexcept { return float128(0, 0, 0x3FFF, 0); }
    /**
     * @brief  Returns a value of 0.5
     * @return 0.5
     */
    [[nodiscard]] FP128_FORCE_INLINE static constexpr float128 half() noexcept { return float128(0, 0, 0x3FFE, 0); }
    /**
     * @brief Return 0.1 using maximum precision
     * @return
     */
    [[nodiscard]] FP128_FORCE_INLINE static constexpr float128 tenth() noexcept
    {
        // 0.1 using maximum precision
        return float128(0x999999999999999A, 0x999999999999, EXP_BIAS - 4, 0u);
    }

    /**
     * @name Mathematical constants
     *
     * The binary128 nearest each constant, given as the encoding rather than parsed from a decimal
     * string: the string constructor is accurate to about an ulp, which is one bit too few to
     * define a constant the rest of the library computes from. The set mirrors <numbers>, so a
     * generic function template written against std::numbers has the same names available here.
     * @{
     */
    /// @brief 2*pi
    [[nodiscard]] FP128_FORCE_INLINE static constexpr float128 two_pi() noexcept { return float128(0x8469898CC51701B8, 0x921FB54442D1, 0x4001, 0); }
    /// @brief 1/pi
    [[nodiscard]] FP128_FORCE_INLINE static constexpr float128 inv_pi() noexcept { return float128(0x2A53F84EAFA3EA6A, 0x45F306DC9C88, 0x3FFD, 0); }
    /// @brief 1/sqrt(pi)
    [[nodiscard]] FP128_FORCE_INLINE static constexpr float128 inv_sqrt_pi() noexcept { return float128(0xD11AE3A914FED7FE, 0x20DD750429B6, 0x3FFE, 0); }
    /// @brief Natural logarithm of 2
    [[nodiscard]] FP128_FORCE_INLINE static constexpr float128 ln2() noexcept { return float128(0xF35793C7673007E6, 0x62E42FEFA39E, 0x3FFE, 0); }
    /// @brief Natural logarithm of 10
    [[nodiscard]] FP128_FORCE_INLINE static constexpr float128 ln10() noexcept { return float128(0x582DD4ADAC5705A6, 0x26BB1BBB5551, 0x4000, 0); }
    /// @brief Base 2 logarithm of e
    [[nodiscard]] FP128_FORCE_INLINE static constexpr float128 log2_e() noexcept { return float128(0xE1777D0FFDA0D23A, 0x71547652B82F, 0x3FFF, 0); }
    /// @brief Base 10 logarithm of e
    [[nodiscard]] FP128_FORCE_INLINE static constexpr float128 log10_e() noexcept { return float128(0xE32A6AB7555F5A68, 0xBCB7B1526E50, 0x3FFD, 0); }
    /// @brief Base 10 logarithm of 2
    [[nodiscard]] FP128_FORCE_INLINE static constexpr float128 log10_2() noexcept { return float128(0xEF311F12B35816F9, 0x34413509F79F, 0x3FFD, 0); }
    /// @brief Square root of 3
    [[nodiscard]] FP128_FORCE_INLINE static constexpr float128 sqrt_3() noexcept { return float128(0xA73B25742D7078B8, 0xBB67AE8584CA, 0x3FFF, 0); }
    /// @brief 1/sqrt(3)
    [[nodiscard]] FP128_FORCE_INLINE static constexpr float128 inv_sqrt_3() noexcept { return float128(0xC4D218F81E4AFB25, 0x279A74590331, 0x3FFE, 0); }
    /// @brief The Euler-Mascheroni constant
    [[nodiscard]] FP128_FORCE_INLINE static constexpr float128 egamma() noexcept { return float128(0x8F49A37C7F0202A6, 0x2788CFC6FB61, 0x3FFE, 0); }
    /// @brief The golden ratio
    [[nodiscard]] FP128_FORCE_INLINE static constexpr float128 phi() noexcept { return float128(0x7C15F39CC0605CEE, 0x9E3779B97F4A, 0x3FFF, 0); }
    /// @}
    /**
     * @brief calculates 10^e
     * @param e integer exponent, in the range
     * @return 10^e
     */
    [[nodiscard]] FP128_INLINE static constexpr float128 exp10(int32_t e) noexcept
    {
        // check the limits first
        if (e < -4965) {
            return 0;
        } else if (e > 4932) {
            return inf();
        }

        // calculate the exponent optimally
        float128 res = 1;
        float128 b;
        if (e > 0)
            b = 10;  // 10^1
        else if (e < 0) {
            b = tenth();
            e = -e;
        }

        while (e > 0) {
            if (e & 1)
                res *= b;
            e >>= 1;
            b.square();
        }
        return res;
    }

    /// @brief Returns false; this type does not conform to IEEE 754-1985.
    static constexpr bool is754version1985(void) { return false; }
    /// @brief Returns true; this type conforms to IEEE 754-2008 (binary128).
    static constexpr bool is754version2008(void) { return true; }
    //
    // End of class method implementation
    //

    //
    // Binary math operators
    //
    // Each of these applies the compound assignment to the by value left hand side and then returns
    // it on a line of its own, rather than the shorter `return lhs OP= rhs;`. The two are
    // equivalent, but the short form makes the return statement copy construct from an lvalue, which
    // MSVC compiles into a store forwarding stall. int128_base's binary operators carry the full
    // explanation and the measurement.
    //
    /**
     * @brief Adds 2 values and returns the result.
     * @param lhs left hand side operand
     * @param rhs Right hand side operand
     * @return Result of the operation
     */
    template <typename T> [[nodiscard]] friend FP128_FORCE_INLINE constexpr float128 operator+(float128 lhs, const T& rhs) noexcept
    {
        lhs += rhs;
        return lhs;
    }
    /**
     * @brief subtracts the right hand side operand to this object to and returns the result.
     * @param lhs left hand side operand
     * @param rhs Right hand side operand
     * @return The float128 result
     */
    template <typename T> [[nodiscard]] friend FP128_FORCE_INLINE constexpr float128 operator-(float128 lhs, const T& rhs) noexcept
    {
        lhs -= rhs;
        return lhs;
    }
    /**
     * @brief Multiplies the right hand side operand with this object to and returns the result.
     * @param lhs left hand side operand
     * @param rhs Right hand side operand
     * @return The float128 result
     */
    template <typename T> [[nodiscard]] friend FP128_FORCE_INLINE constexpr float128 operator*(float128 lhs, const T& rhs) noexcept
    {
        lhs *= rhs;
        return lhs;
    }
    /**
     * @brief Divides this object by the right hand side operand and returns the result.
     * @param lhs left hand side operand
     * @param rhs Right hand side operand
     * @return The float128 result
     */
    template <typename T> [[nodiscard]] friend FP128_FORCE_INLINE float128 operator/(float128 lhs, const T& rhs)
    {
        lhs /= rhs;
        return lhs;
    }

    //
    // Binary math operators with the scalar on the left hand side
    //
    // Without these, an expression like (1 + x) is ambiguous: converting the literal to float128
    // and converting x to a builtin type are both one user defined conversion, so neither
    // overload wins. Restricting the left operand to the arithmetic types keeps these from
    // competing with the float128 on the left versions above, which would otherwise be an equally
    // good match. The comparison operators already carry the same pair of overloads.
    //

    /// @brief Adds a scalar and a float128, in that order. @param lhs Left operand @param rhs Right operand @return The float128 result
    template <typename T>
        requires std::is_arithmetic_v<T>
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr float128 operator+(const T& lhs, const float128& rhs) noexcept
    {
        return float128(lhs) += rhs;
    }
    /// @brief Subtracts a float128 from a scalar. @param lhs Left operand @param rhs Right operand @return The float128 result
    template <typename T>
        requires std::is_arithmetic_v<T>
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr float128 operator-(const T& lhs, const float128& rhs) noexcept
    {
        return float128(lhs) -= rhs;
    }
    /// @brief Multiplies a scalar and a float128, in that order. @param lhs Left operand @param rhs Right operand @return The float128 result
    template <typename T>
        requires std::is_arithmetic_v<T>
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr float128 operator*(const T& lhs, const float128& rhs) noexcept
    {
        return float128(lhs) *= rhs;
    }
    /// @brief Divides a scalar by a float128. @param lhs Left operand @param rhs Right operand @return The float128 result
    template <typename T>
        requires std::is_arithmetic_v<T>
    [[nodiscard]] friend FP128_FORCE_INLINE float128 operator/(const T& lhs, const float128& rhs)
    {
        return float128(lhs) /= rhs;
    }
    /**
     * @brief Shifts a scalar left by a float128 shift count, scaling it by a power of two.
     *
     * The left operand is widened to float128 first, so the result is 128 bit rather than the
     * builtin type of lhs.
     *
     * The count is rhs converted to int32_t, the same conversion the float128 on the left overload
     * applies. Note that conversion rounds to nearest rather than truncating, so a count of 3.75
     * shifts by 4. fixed_point128 truncates in the same place, its conversion being a plain shift.
     *
     * @param lhs Left operand, the value being shifted
     * @param rhs Right operand, the shift count
     * @return The float128 result
     */
    template <typename T>
        requires std::is_arithmetic_v<T>
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr float128 operator<<(const T& lhs, const float128& rhs) noexcept
    {
        return float128(lhs) <<= static_cast<int32_t>(rhs);
    }
    /**
     * @brief Shifts a scalar right by a float128 shift count, scaling it by a power of two.
     * The left operand is widened to float128 first, see operator<< above.
     * @param lhs Left operand, the value being shifted
     * @param rhs Right operand, the shift count
     * @return The float128 result
     */
    template <typename T>
        requires std::is_arithmetic_v<T>
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr float128 operator>>(const T& lhs, const float128& rhs) noexcept
    {
        return float128(lhs) >>= static_cast<int32_t>(rhs);
    }

    //
    // Comparison operators
    //

    /**
     * @brief Compare logical/bitwise equal.
     * @param lhs left hand side operand
     * @param rhs Right hand side operand
     * @return True if this and other are equal.
     */
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr bool operator==(const float128& lhs, const float128& rhs) noexcept
    {
        // A NaN compares equal to nothing, not even to another NaN with the same bits.
        if (lhs.is_nan() || rhs.is_nan())
            return false;
        // Positive and negative zero are numerically equal even though their bits differ.
        if (lhs.is_zero() && rhs.is_zero())
            return true;
        return lhs.high == rhs.high && lhs.low == rhs.low;
    }
    /// @overload
    // Constrained to a builtin type. Unconstrained, T deduces to float128 for a comparison of two
    // of them, which makes this template an exact match alongside the non-template above - and
    // under the C++20 rewriting of == and != the reversed forms join in, which MSVC reports as an
    // ambiguity from inside any standard container that compares keys.
    template <typename T>
        requires std::is_arithmetic_v<T>
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr bool operator==(const float128& lhs, const T& rhs) noexcept { return lhs == float128(rhs); }
    /// @overload
    // Constrained to a builtin type, see operator== above.
    template <typename T>
        requires std::is_arithmetic_v<T>
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr bool operator==(const T& lhs, const float128& rhs) noexcept
    {
        return rhs == float128(lhs);
    }
    /**
     * @brief Return true when objects are not equal. Can be used as logical XOR.
     * @param lhs left hand side operand
     * @param rhs Right hand side operand
     * @return True if not equal.
     */
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr bool operator!=(const float128& lhs, const float128& rhs) noexcept { return !(lhs == rhs); }
    /// @overload
    // Constrained to a builtin type, see operator== above.
    template <typename T>
        requires std::is_arithmetic_v<T>
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr bool operator!=(const float128& lhs, const T& rhs) noexcept { return lhs != float128(rhs); }
    /// @overload
    /// @copydoc operator!=
    template <typename T>
        requires std::is_arithmetic_v<T>
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr bool operator!=(const T& lhs, const float128& rhs) noexcept { return rhs != float128(lhs); }
    /**
     * @brief Return true if this object is small than the other
     * @param lhs left hand side operand
     * @param rhs Right hand side operand
     * @return True when this object is smaller.
     */
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr bool operator<(const float128& lhs, const float128& rhs) noexcept
    {
        // A NaN is unordered with everything, so every relational test involving one is false.
        if (lhs.is_nan() || rhs.is_nan())
            return false;
        // The two zeros are numerically equal, so neither is smaller than the other. Without this
        // the differing sign bits would make -0 compare smaller than +0.
        if (lhs.is_zero() && rhs.is_zero())
            return false;

        auto rhs_sign = rhs.get_sign();
        auto lhs_sign = lhs.get_sign();

        // signs are different
        if (lhs_sign != rhs_sign)
            return lhs_sign > rhs_sign;  // true when lhs_sign is 1 and rhs.sign is 0

        // MSB is the same, check the LSB, implies the exponent is identical
        if (lhs.high == rhs.high)
            return (lhs_sign) ? lhs.low > rhs.low : lhs.low < rhs.low;

        return (lhs_sign) ? lhs.high > rhs.high : lhs.high < rhs.high;
    }
    /// @overload
    template <typename T>
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr bool operator<(const float128& lhs, const T& rhs) noexcept { return lhs < float128(rhs); }
    /// @overload
    template <typename T>
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr bool operator<(const T& lhs, const float128& rhs) noexcept { return float128(lhs) < rhs; }
    /**
     * @brief Return true this object is small or equal than the other
     * @param lhs left hand side operand
     * @param rhs Right hand side operand
     * @return True when this object is smaller or equal.
     */
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr bool operator<=(const float128& lhs, const float128& rhs) noexcept
    {
        // Not simply !(lhs > rhs): a NaN makes every relational test false, so negating the
        // opposite test would wrongly report that a NaN is less than or equal to everything.
        if (lhs.is_nan() || rhs.is_nan())
            return false;
        return !(lhs > rhs);
    }
    /// @overload
    template <typename T>
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr bool operator<=(const float128& lhs, const T& rhs) noexcept { return lhs <= float128(rhs); }
    /// @overload
    template <typename T>
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr bool operator<=(const T& lhs, const float128& rhs) noexcept { return float128(lhs) <= rhs; }
    /**
     * @brief Return true this object is larger than the other
     * @param lhs left hand side operand
     * @param rhs Right hand side operand
     * @return True when this object is larger.
     */
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr bool operator>(const float128& lhs, const float128& rhs) noexcept
    {
        // A NaN is unordered with everything, so every relational test involving one is false.
        if (lhs.is_nan() || rhs.is_nan())
            return false;
        // the two zeros are numerically equal, so neither is larger than the other
        if (lhs.is_zero() && rhs.is_zero())
            return false;

        auto rhs_sign = rhs.get_sign();
        auto lhs_sign = lhs.get_sign();

        // signs are different
        if (lhs_sign != rhs_sign)
            return lhs_sign < rhs_sign;  // true when lhs_sign is 1 and rhs.sign is 0

        // MSB is the same, check the LSB, implies the exponent is identical
        if (lhs.high == rhs.high)
            return (lhs_sign) ? lhs.low < rhs.low : lhs.low > rhs.low;

        return (lhs_sign) ? lhs.high < rhs.high : lhs.high > rhs.high;
    }
    /// @overload
    template <typename T>
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr bool operator>(const float128& lhs, const T& rhs) noexcept { return lhs > float128(rhs); }
    /// @overload
    template <typename T>
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr bool operator>(const T& lhs, const float128& rhs) noexcept { return float128(lhs) > rhs; }
    /**
     * @brief Return true this object is larger or equal than the other
     * @param lhs left hand side operand
     * @param rhs Right hand side operand
     * @return True when this objext is larger or equal.
     */
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr bool operator>=(const float128& lhs, const float128& rhs) noexcept
    {
        // see the note on operator<= about why this is not simply !(lhs < rhs)
        if (lhs.is_nan() || rhs.is_nan())
            return false;
        return !(lhs < rhs);
    }
    /// @overload
    template <typename T>
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr bool operator>=(const float128& lhs, const T& rhs) noexcept { return lhs >= float128(rhs); }
    /// @overload
    template <typename T>
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr bool operator>=(const T& lhs, const float128& rhs) noexcept { return float128(lhs) >= rhs; }

    /**
     * @brief Return the NaN constant
     * @param
     * @return
     */
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr float128 nan(const float128&) { return float128::nan(); }
    /**
     * @brief Tests if the value is a NaN
     * @param x Value to test
     * @return True when the value is a NaN
     */
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr bool isnan(const float128& x)
    {
        // zero for +- INF, non-zero for NaN
        return x.is_nan();
    }
    /**
     * @brief Tests if the value is an Infinite (negative or positive)
     * @param x Value to test
     * @return True when the value is an Infinite
     */
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr bool isinf(const float128& x) { return x.is_inf(); }

    /**
     * @brief Returns the absolute value of x.
     * @param x Input value
     * @return |x|
     */
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr float128 fabs(const float128& x) noexcept
    {
        float128 temp = x;
        temp.set_sign(0);
        return temp;
    }
    /**
     * @brief Performs the floor() function, similar to libc's floor(), rounds down towards -infinity.
     * @param x Input value
     * @return A float128 holding the integer value. Overflow is not reported.
     */
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr float128 floor(const float128& x) noexcept
    {
        float128 fraction = x.get_fraction();
        if (fraction.is_zero())
            return x;

        float128 res = x - fraction;
        if (fraction.is_negative())
            return res - 1;
        return res;
    }
    /**
     * @brief Performs the ceil() function, similar to libc's ceil(), rounds up towards infinity.
     * @param x Input value
     * @return A float128 holding the integer value. Overflow is not reported.
     */
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr float128 ceil(const float128& x) noexcept
    {
        float128 fraction = x.get_fraction();
        if (fraction.is_zero())
            return x;

        float128 res = x - fraction;
        if (fraction.is_positive())
            return res + 1;
        return res;
    }
    /**
     * @brief Rounds towards zero
     * @param x Value to truncate
     * @return Integer value, rounded towards zero.
     */
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr float128 trunc(const float128& x) noexcept
    {
        float128 fraction = x.get_fraction();
        if (fraction.is_zero())
            return x;

        return x - fraction;
    }
    /**
     * @brief Rounds towards the nearest integer.
     * The halfway value (0.5) is rounded away from zero.
     * @param x Value to round
     * @return Integer value, rounded towards the nearest integer.
     */
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr float128 round(const float128& x) noexcept
    {
        float128 h = (x.is_positive()) ? half() : -half();
        return trunc(x + h);
    }
    /**
     * @brief Rounds x to the nearest integer and returns the result as int64_t.
     * @param x Input value
     * @return Nearest integer as int64_t. Returns 0 on overflow.
     */
    // rint rounds ties to even, round() rounds them away from zero, so these cannot forward to
    // llround/lround the way they used to: llrint(2.5) is 2 and llround(2.5) is 3.
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr int64_t llrint(const float128& x) noexcept { return static_cast<int64_t>(rint(x)); }
    /**
     * @brief Rounds towards the nearest integer.
     * The halfway value (0.5) is rounded away from zero.
     * @param x Value to round
     * @return Integer value, rounded towards the nearest integer.
     */
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr int64_t llround(const float128& x) noexcept
    {
        float128 res = round(x);
        if (res.is_special() || res > INT64_MAX || res < INT64_MIN)
            return 0;
        return static_cast<int64_t>(res);
    }
    /**
     * @brief Rounds x to the nearest integer and returns the result as int32_t.
     * @param x Input value
     * @return Nearest integer as int32_t. Returns 0 on overflow.
     */
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr int32_t lrint(const float128& x) noexcept { return static_cast<int32_t>(rint(x)); }
    /**
     * @brief Rounds towards the nearest integer.
     * The halfway value (0.5) is rounded away from zero.
     * @param x Value to round
     * @return Integer value, rounded towards the nearest integer.
     */
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr int32_t lround(const float128& x) noexcept
    {
        float128 res = round(x);
        if (res.is_special() || res > INT32_MAX || res < INT32_MIN)
            return 0;
        return static_cast<int32_t>(res);
    }

    /**
     * @brief Retrieves an integer that represents the base-2 exponent of the specified value.
     * @param x The specified value.
     * @return Integer value, rounded towards the nearest integer.
     */
    [[nodiscard]] friend FP128_INLINE constexpr int32_t ilogb(const float128& x) noexcept
    {
        // The three answers <cmath> reserves for the arguments that have no exponent. A subnormal
        // does have one, but not the one its exponent field holds, so it goes through
        // get_components() with the rest.
        if (x.is_zero())
            return FP_ILOGB0;
        if (x.is_nan())
            return FP_ILOGBNAN;
        if (x.is_inf())
            return INT_MAX;

        uint64_t l = 0, h = 0;
        int32_t expo = 0;
        uint32_t sign = 0;
        x.get_components(l, h, expo, sign);
        return expo;
    }
    /**
     * @brief returns the value of x with the sign of y.
     * @param x The value that's returned as the magnitude of the result.
     * @param y The sign of the result.
     * @return The copysign functions return a floating-point value that combines the magnitude of x and the sign of y.
     */
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr float128 copysign(const float128& x, const float128& y) noexcept
    {
        float128 temp = x;
        temp.set_sign(y.get_sign());
        return temp;
    }
    /**
     * @brief Performs the fmod() function, similar to libc's fmod(), returns the remainder of a division x/root.
     * @param x Numerator
     * @param y Denominator
     * @return The modulo value.
     */
    [[nodiscard]] friend float128 fmod(const float128& x, const float128& y) noexcept
    {
        // a NaN operand propagates
        if (x.is_nan() || y.is_nan())
            return nan();

        // trivial case, x is zero
        if (x.is_zero())
            return x;

        // fmod(x, 0) is the invalid operation and produces a NaN, matching the CRT.
        // An infinite dividend is invalid for the same reason.
        if (y.is_zero() || x.is_inf())
            return nan();

        // an infinite divisor leaves the dividend untouched
        if (y.is_inf())
            return x;

        // fmod is exact: the result is x - n*y for an integer n, and it is representable. Reaching
        // it through trunc(x/y), as this used to, is only exact while the quotient is: once x/y
        // needs more than 113 bits the truncated quotient is a rounded value and the subtraction
        // returns something unrelated to the remainder. The mantissas are integers, so the whole
        // operation is a long division on them, done here one 64 bit block of the exponent
        // difference at a time.
        uint64_t lx = 0, hx = 0, ly = 0, hy = 0;
        int32_t ex = 0, ey = 0;
        uint32_t sx = 0, sy = 0;
        x.get_components(lx, hx, ex, sx);
        y.get_components(ly, hy, ey, sy);

        // |x| < |y| leaves the dividend as the remainder
        if (ex < ey || (ex == ey && uint128_t(lx, hx) < uint128_t(ly, hy)))
            return x;

        // Both mantissas hold 113 bits with the leading one at bit 112, so the values compared
        // here are x = rem * 2^(ex-112) and y = mod * 2^(ey-112).
        uint128_t rem(lx, hx);
        const uint128_t mod(ly, hy);
        rem %= mod;
        int32_t shift = ex - ey;

        // A remainder is below the divisor and so holds at most 113 bits; shifting 15 more in
        // keeps the dividend inside the 128 the modulo can take.
        constexpr int32_t block = 15;
        while (shift > 0 && !rem.is_zero()) {
            const int32_t step = (shift < block) ? shift : block;
            rem <<= step;
            shift -= step;
            rem %= mod;
        }

        if (rem.is_zero()) {
            float128 zero;
            zero.set_sign(sx);
            return zero;
        }

        // Rebuild the value: the remainder is an integer scaled by the divisor's last place.
        uint64_t rl = 0, rh = 0;
        rem.get_components(rl, rh);
        const int32_t msb = static_cast<int32_t>(log2(rl, rh));
        shift_left128_inplace_safe(rl, rh, FRAC_BITS - msb);

        float128 res;
        res.set_components(rl, rh, ey - FRAC_BITS + msb, sx);
        return res;
    }
    /**
     * @brief Split into integer and fraction parts.
     * Both results carry the sign of the input variable.
     * @param x Input value
     * @param iptr Pointer to float128 holding the integer part of x.
     * @return The fraction part of x. Undefined when iptr is nullptr.
     */
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr float128 modf(const float128& x, float128* iptr) noexcept
    {
        if (iptr == nullptr)
            return 0;

        // fraction
        float128 res = x.get_fraction();
        // integer
        *iptr = x - res;
        return res;
    }
    /**
     * @brief Determines the positive difference between the first and second values.
     * @param x First value
     * @param y Second value
     * @return If x > y returns x - y. Otherwise zero.
     */
    [[nodiscard]] friend FP128_INLINE constexpr float128 fdim(const float128& x, const float128& y) noexcept
    {
        if (x.is_nan() || y.is_nan())
            return nan();
        return (x > y) ? x - y : float128();
    }
    /**
     * @brief Returns the minimum between x and y.
     * @param x First value
     * @param y Second value
     * @return If x < y returns x. Otherwise y.
     */
    [[nodiscard]] friend FP128_INLINE constexpr float128 fmin(const float128& x, const float128& y) noexcept
    {
        // A NaN operand is treated as missing rather than as a value, so the other one wins. Every
        // comparison against a NaN is false, which made the plain conditional return whichever
        // operand happened to sit on the false branch.
        if (x.is_nan())
            return y;
        if (y.is_nan())
            return x;
        // The comparison operators treat the two zeros as equal, and the standard asks for the
        // negative one here.
        if (x.is_zero() && y.is_zero())
            return x.is_negative() ? x : y;
        return (x < y) ? x : y;
    }
    /**
     * @brief Returns the maximum between x and y.
     * @param x First value
     * @param y Second value
     * @return If x > y returns x. Otherwise y.
     */
    [[nodiscard]] friend FP128_INLINE constexpr float128 fmax(const float128& x, const float128& y) noexcept
    {
        if (x.is_nan())
            return y;
        if (y.is_nan())
            return x;
        if (x.is_zero() && y.is_zero())
            return x.is_positive() ? x : y;
        return (x > y) ? x : y;
    }
    /**
     * @brief Calculates the hypotenuse. i.e. sqrt(x^2 + y^2)
     * @param x First value
     * @param y Second value
     * @return sqrt(x^2 + y^2).
     */
    [[nodiscard]] friend FP128_FORCE_INLINE float128 hypot(const float128& x, const float128& y) noexcept { return sqrt(sqr(x) + sqr(y)); }
    /**
     * @brief Calculates the square of a value. i.e. x^2
     *
     * Faster than x * x and identical to it, see float128::square().
     *
     * @param x Value to square
     * @return x^2, which is never negative.
     */
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr float128 sqr(float128 x) noexcept { return x.square(); }
    /**
     * @brief Calculates the square root using Newton's method.
     * Based on the book "Math toolkit for real time programming" by Jack W. Crenshaw
     * @param x Value to calculate the root of
     * @param iterations how many iterations to perform (more is more accurate). Sensible values are 0-5.
     * @return Square root of (x), zero when x <= 0.
     */
    [[nodiscard]] friend float128 sqrt(const float128& x, uint32_t iterations) noexcept
    {
        if (x.is_nan())
            return x;
        if (x.is_zero())
            return x;  // sqrt(-0) is -0
        if (x.is_negative())
            return float128::nan();
        if (x.is_inf())
            return x;

        // Split off an even power of two, leaving a mantissa in [1, 4).
        //
        // Halving an even exponent is exact, so the result needs no correction factor afterwards.
        // Normalizing to [0.5, 1) instead, as this used to, leaves an odd exponent half the time
        // and has to multiply the root by sqrt(2)/2 to make up for it - an irrational factor, and
        // therefore a rounding, applied to a value that was otherwise correct to the last bit.
        //
        // Reading the mantissa through get_components() also normalizes a subnormal argument,
        // which the previous exponent arithmetic did not account for at all.
        uint64_t l = 0, h = 0;
        int32_t expo = 0;
        uint32_t sign = 0;
        x.get_components(l, h, expo, sign);
        const bool odd = (expo & 1) != 0;
        const float128 norm_x(l, h, static_cast<uint32_t>(EXP_BIAS + (odd ? 1 : 0)), 0);
        const int32_t half_expo = (expo - (odd ? 1 : 0)) >> 1;

        // The hardware double gives 53 correct bits to start from, and every Newton step doubles
        // them, so three passes cover the 113 the mantissa holds several times over.
        float128 root = ::sqrt(static_cast<double>(norm_x));

        // iterate several times via Newton's method
        //                  X
        //   Xn+1 = 0.5 * (---- + Xn )
        //                  Xn
        for (auto i = iterations; i != 0; --i) {
            root = (norm_x / root + root) >> 1;
        }

        return ldexp(root, half_expo);
    }
    /**
     * @brief Calculates the cube root.
     * Uses the Halley's method.
     * @param x Floating point value
     * @param iterations how many Halley to perform, usually 1 is enough
     * @return cube root of x
     */
    [[nodiscard]] friend FP128_INLINE float128 cbrt(const float128 x, uint32_t iterations) noexcept
    {
        if (x.is_nan())
            return x;
        // The real cube root is defined for a negative argument and is an odd function, which is
        // what the C library's cbrt computes. Returning a NaN, as this used to, made cbrt(-8)
        // undefined where the standard has it equal to -2.
        if (x.is_zero() || x.is_inf())
            return x;

        // Split off a power of two whose exponent is a multiple of three, leaving a mantissa in
        // [1, 8). Dividing that exponent by three is exact, so no correction factor is needed.
        //
        // The factors this used to multiply by, cbrt(2)/2 and cbrt(4)/2, were written out to the
        // 17 digits a double round trips through and no further, so two thirds of all arguments
        // came back with 53 correct bits instead of 113.
        uint64_t l = 0, h = 0;
        int32_t expo = 0;
        uint32_t sign = 0;
        x.get_components(l, h, expo, sign);
        int32_t rem = expo % 3;
        if (rem < 0)
            rem += 3;
        const float128 norm_x(l, h, static_cast<uint32_t>(EXP_BIAS + rem), 0);
        const int32_t third_expo = (expo - rem) / 3;

        // Halley's method triples the correct digits per pass, so two passes take the 53 bits the
        // hardware seed provides past the 113 the mantissa holds.
        float128 root = ::cbrt(static_cast<double>(norm_x));

        //                3
        //              Xn  + 2X
        //   Xn+1 = Xn ----------
        //                3
        //              2Xn  + X
        const auto x2 = norm_x << 1;
        for (auto i = iterations; i != 0; --i) {
            float128 r_cube = root * root * root;
            root = root * (r_cube + x2) / ((r_cube << 1) + norm_x);
        }

        root = ldexp(root, third_expo);
        root.set_sign(sign);
        return root;
    }
    /**
     * @brief Calculates the reciprocal of a value. y = 1 / x
     * Using newton iterations: Yn+1 = Yn(2 - x * Yn)
     * @param x Input value
     * @return 1 / x. Returns zero on overflow or division by zero
     */
    [[nodiscard]] friend FP128_INLINE float128 reciprocal(const float128& x) noexcept
    {
        static const float128 two = 2;
        constexpr int max_iterations = 3;
        constexpr int debug = false;
        const auto x_sign = x.get_sign();
        if (x.is_nan())
            return x;
        if (x.is_inf()) {
            float128 zero;
            zero.set_sign(x_sign);
            return zero;
        }
        if (x.is_zero())
            return (x_sign) ? -inf() : inf();

        // get_components() normalizes a subnormal, so the mantissa below is always in [1, 2) and
        // the Newton iteration sees the same problem whatever the magnitude of x was. Rebuilding
        // the value from the mantissa alone is also what keeps the scaling exact.
        uint64_t l = 0, h = 0;
        int32_t expo = 0;
        uint32_t sign = 0;
        x.get_components(l, h, expo, sign);
        const float128 norm_x(l, h, EXP_BIAS, 0);

        float128 y = 1.0 / static_cast<double>(norm_x);

        if (!y)
            return y;

        float128 xy, y_prev;
        // Newton iterations:
        int i = 0;
        for (; i < max_iterations && (y_prev != y); ++i) {
            y_prev = y;
            xy = norm_x * y;
            // y = y * (two - xy);
            y *= two - xy;
        }

        if constexpr (debug) {
            static int debug_max_iter = 0;
            if (i > debug_max_iter || i == max_iterations) {
                debug_max_iter = i;
                printf("reciprocal took %i iterations for %.10lf\n", i, static_cast<double>(x));
            }
        }

        // The mantissa was divided out above, so undoing the exponent is a scaling rather than a
        // replacement. Overwriting it, as this used to, was right only while y kept an exponent of
        // -1, which is every value except the one case y comes out exactly 1: reciprocal() of any
        // power of two was off by a factor of two, and exp() of a small negative argument with it.
        y = ldexp(y, -expo);
        y.set_sign(x_sign);
        return y;
    }
    /**
     * @brief Factorial reciprocal (inverse). Calculates 1 / x!
     * Maximum supported value of x is 50.
     * @param x Input value
     * @param res Result of the function
     */
    friend FP128_INLINE void fact_reciprocal(int x, float128& res) noexcept
    {
        // Given as encodings rather than decimal strings. Several of the strings this replaced
        // were the 17 digit expansion of a double - 1/6! among them - which capped the accuracy of
        // every series built on the table at 53 bits. The exponential's own reduction is good to
        // better than 2^-110, and a single term formed from a double sized constant was throwing
        // 40 of those bits away.
        static constexpr float128 c[] = {
            float128(0x0000000000000000, 0x000000000000, 0x3FFF, 0),  // 1 / 0!
            float128(0x0000000000000000, 0x000000000000, 0x3FFF, 0),  // 1 / 1!
            float128(0x0000000000000000, 0x000000000000, 0x3FFE, 0),  // 1 / 2!
            float128(0x5555555555555555, 0x555555555555, 0x3FFC, 0),  // 1 / 3!
            float128(0x5555555555555555, 0x555555555555, 0x3FFA, 0),  // 1 / 4!
            float128(0x1111111111111111, 0x111111111111, 0x3FF8, 0),  // 1 / 5!
            float128(0x6C16C16C16C16C17, 0x6C16C16C16C1, 0x3FF5, 0),  // 1 / 6!
            float128(0xA01A01A01A01A01A, 0xA01A01A01A01, 0x3FF2, 0),  // 1 / 7!
            float128(0xA01A01A01A01A01A, 0xA01A01A01A01, 0x3FEF, 0),  // 1 / 8!
            float128(0x38FAAC1C88E50017, 0x71DE3A556C73, 0x3FEC, 0),  // 1 / 9!
            float128(0xC72EF016D3EA6679, 0x27E4FB7789F5, 0x3FE9, 0),  // 1 / 10!
            float128(0x38FE747E4B837DC7, 0xAE64567F544E, 0x3FE5, 0),  // 1 / 11!
            float128(0x7B544DA987ACFE85, 0x1EED8EFF8D89, 0x3FE2, 0),  // 1 / 12!
            float128(0x97CA38331D23AF68, 0x6124613A86D0, 0x3FDE, 0),  // 1 / 13!
            float128(0xD20BADF145DFA3E5, 0x93974A8C07C9, 0x3FDA, 0),  // 1 / 14!
            float128(0xF11D8656B0EE8CB0, 0xAE7F3E733B81, 0x3FD6, 0),  // 1 / 15!
            float128(0xF11D8656B0EE8CB0, 0xAE7F3E733B81, 0x3FD2, 0),  // 1 / 16!
            float128(0xA6B2605197771B00, 0x952C77030AD4, 0x3FCE, 0),  // 1 / 17!
            float128(0x77BB004886A2C2AB, 0x6827863B97D9, 0x3FCA, 0),  // 1 / 18!
            float128(0x724CA1EC3B7B9675, 0x2F49B4681415, 0x3FC6, 0),  // 1 / 19!
            float128(0x507A9CAD2BF8F0BB, 0xE542BA402022, 0x3FC1, 0),  // 1 / 20!
            float128(0x18BEF146FCEE6E45, 0x71B8EF6DCF57, 0x3FBD, 0),  // 1 / 21!
            float128(0x29450C90B7F338EC, 0x0CE396DB7F85, 0x3FB9, 0),  // 1 / 22!
            float128(0x9D97B8704DD7F628, 0x761B41316381, 0x3FB4, 0),  // 1 / 23!
            float128(0x7CCA4B4067CA9D8A, 0xF2CF01972F57, 0x3FAF, 0),  // 1 / 24!
            float128(0x8D4E44A419776F11, 0x3F3CCDD165FA, 0x3FAB, 0),  // 1 / 25!
            float128(0x9A38F2050BA6B015, 0x88E85FC6A4E5, 0x3FA6, 0),  // 1 / 26!
            float128(0x320A9A18F15D4277, 0xD1AB1C2DCCEA, 0x3FA1, 0),  // 1 / 27!
            float128(0xD373C5C51C354A8D, 0x0A18A2635085, 0x3F9D, 0),  // 1 / 28!
            float128(0xD7ABE30E7766F129, 0x259F98B4358A, 0x3F98, 0),  // 1 / 29!
            float128(0xE60CADED4C2989C5, 0x3932C5047D60, 0x3F93, 0),  // 1 / 30!
            float128(0xC42E1EE46FA6BFC4, 0x434D2E783F5B, 0x3F8E, 0),  // 1 / 31!
            float128(0xC42E1EE46FA6BFC4, 0x434D2E783F5B, 0x3F89, 0),  // 1 / 32!
            float128(0x1B5382CDFFA97422, 0x3981254DD0D5, 0x3F84, 0),  // 1 / 33!
            float128(0xA13F8A2B4AF9D6B7, 0x2710231C0FD7, 0x3F7F, 0),  // 1 / 34!
            float128(0xF2833C7F5A7E0624, 0x0DC59C716D91, 0x3F7A, 0),  // 1 / 35!
            float128(0x92B06B8D12A7275C, 0xDF983290C2CA, 0x3F74, 0),  // 1 / 36!
            float128(0xAF4C78B15C3D89D3, 0x9EC8D1C94E85, 0x3F6F, 0),  // 1 / 37!
            float128(0xAE913D370A4EC4E8, 0x5D4ACB9C0C3A, 0x3F6A, 0),  // 1 / 38!
            float128(0xDE0104476AEB4C3B, 0x1E99449A4BAC, 0x3F65, 0),  // 1 / 39!
            float128(0x3001A07244ABAD2B, 0xCA8ED42A12AE, 0x3F5F, 0),  // 1 / 40!
            float128(0x0C7E25CFD1B1B2DD, 0x65E61C39D024, 0x3F5A, 0),  // 1 / 41!
            float128(0x836C4D9225DCB90A, 0x10AF527530DE, 0x3F55, 0),  // 1 / 42!
            float128(0x22DCBAE56DEF3720, 0x95DB45257E51, 0x3F4F, 0),  // 1 / 43!
            float128(0xA4FD9F327E7F6DE9, 0x272B1B03FEC6, 0x3F4A, 0),  // 1 / 44!
            float128(0x78E02C5E91C64CAC, 0xA3CB87222064, 0x3F44, 0),  // 1 / 45!
            float128(0x062CA46E4F25C608, 0x240804F65951, 0x3F3F, 0),  // 1 / 46!
            float128(0x9B78B454D8B628E5, 0x8DA8E0A127EB, 0x3F39, 0),  // 1 / 47!
            float128(0x67A5CD8DE5CEC5EE, 0x091B406B6FF2, 0x3F34, 0),  // 1 / 48!
            float128(0x5D94A3FD410E1232, 0x5A42F0DFEB08, 0x3F2E, 0),  // 1 / 49!
            float128(0x8205F0A053453602, 0xBB36F6E12CD7, 0x3F28, 0),  // 1 / 50!
        };
        constexpr int series_len = array_length(c);
        static_assert(series_len == 51);

        res = (x >= 0 && x < series_len) ? c[x] : float128();
    }
    /**
     * @brief returns the double factorial of a number
     * @param x
     * @param res
     */
    [[nodiscard]] friend FP128_INLINE float128 double_factorial(int x) noexcept
    {
        constexpr int32_t arr_size = 100;
        static float128 c[arr_size];
        if (c[0].is_zero()) {
            c[0] = 1;
            c[1] = 1;

            // compute the odd and even double factorials
            for (int i = 2; i < arr_size - 1; i += 2) {
                c[i] = c[i - 2] * i;
                c[i + 1] = c[i - 1] * (i + 1);
            }
        }
        if (x < arr_size)
            return c[x];

        // TODO: compute the following members
        return float128::inf();
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
     *                     x   ix   fx
     * Using the equality e = e  * e
     * Where ix is the integer part of x and fx is the fraction part.
     * ix is computed via multiplications which won't overflow if the result value can be held.
     * fx is computed via Maclaurin series expansion, but since fx < 1, it won't overflow.
     * @param x A number specifying a power.
     * @return Exponent of x
     */
    /// @brief Terms of the Maclaurin series exp_reduced() runs. |r| <= ln2/2, so r^28/28! is past the last bit.
    static constexpr int32_t EXP_TERMS = 28;

    /**
     * @brief e to the power of a reduced argument, |r| <= ln2/2.
     * @param r Reduced argument
     * @return exp(r)
     */
    [[nodiscard]] FP128_INLINE static float128 exp_reduced(const float128& r) noexcept
    {
        float128 acc, factorial;
        fact_reciprocal(EXP_TERMS, acc);
        for (int32_t n = EXP_TERMS - 1; n >= 0; --n) {
            fact_reciprocal(n, factorial);
            acc = factorial + r * acc;
        }
        return acc;
    }

    /**
     * @brief Computes e to the power of x.
     *
     * x is written as k*ln2 + r with k an integer and |r| <= ln2/2, which turns the answer into
     * 2^k * exp(r): the power of two is an exact scaling and only the reduced argument goes
     * through a series.
     *
     * The reduction is where the accuracy of the whole function is decided. ln2 is irrational, so
     * k*ln2 cannot be represented and subtracting a rounded one leaves an error of k ulp - at the
     * top of the range k reaches 16000 and fourteen bits of the answer are gone before the series
     * even starts. Splitting ln2 into a high part whose low 14 bits are zero, so that k times it
     * is exact, and a low part carrying the rest, keeps r accurate to far below its last bit.
     *
     * The previous implementation raised e to the integer part by repeated squaring, which doubles
     * its own error at every step, and reached the negative half of the range through reciprocal().
     *
     * @param x A number specifying a power.
     * @return e^x
     */
    [[nodiscard]] friend FP128_INLINE float128 exp(const float128& x) noexcept
    {
        if (x.is_nan())
            return x;
        if (x.is_inf())
            return x.is_negative() ? float128() : x;
        if (x.is_zero())
            return float128::one();

        // log of the largest finite value is 11356.52, log of the smallest subnormal is -11433.46
        if (x > float128(11357))
            return inf();
        if (x < float128(-11434))
            return float128();

        // ln2 to 321 bits, in three pieces. The first two have their low 16 mantissa bits zeroed,
        // so multiplying either by a k of up to 2^16 is exact and the products can be removed from
        // the argument without a rounding of their own.
        //
        // Two pieces would only carry ln2 to the 113 bits a float128 holds, which is not enough:
        // the reduction multiplies whatever error the constant has by k, and at the top of the
        // range k reaches 16000, so a constant good to 2^-114 leaves the reduced argument good to
        // only 2^-100 and the answer 8000 ulp wide.
        constexpr float128 ln2_hi(0xF35793C767300000, 0x62E42FEFA39E, 0x3FFE, 0);
        constexpr float128 ln2_mid(0x93394C5B16C50000, 0xF97B57A079A1, 0x3F98, 0);
        constexpr float128 ln2_lo(0x7CF70EC40DBD7593, 0xA2EB71755F45, 0x3F32, 0);

        const int32_t k = static_cast<int32_t>(llround(x * float128::log2_e()));
        const float128 kf(k);
        // The first subtraction is exact, cancelling two values within a factor of two of each
        // other; the others remove quantities far below the last place of the result.
        const float128 r = ((x - kf * ln2_hi) - kf * ln2_mid) - kf * ln2_lo;

        return ldexp(exp_reduced(r), k);
    }
    /**
     * @brief Computes 2 to the power of x
     * @param x Exponent value
     * @return 2^x
     */
    [[nodiscard]] friend FP128_INLINE float128 exp2(const float128& x) noexcept
    {
        //
        // Based on exponent law: (x^n)^m = x^(m*n)
        // Convert the exponent x (function parameter) to produce an exponent that will work with exp()
        // y = log(2)
        // 2^x = e^(y*x) = exp(y*x)
        //
        if (x.is_nan() || x.is_inf())
            return (x.is_inf() && x.is_negative()) ? float128() : x;
        if (x.is_zero())
            return float128::one();
        if (x > float128(16384))
            return inf();
        if (x < float128(-16495))
            return float128();

        // The integer part scales exactly, so only the fraction reaches exp() and the argument it
        // is handed stays small enough that the rounding of the multiply below is the only error.
        // That rounding is then recovered with an exact fma and folded back in: exp(p + e) is
        // exp(p) * (1 + e) to well past the last bit for an e this small.
        const int32_t k = static_cast<int32_t>(llround(x));
        const float128 f = x - float128(k);
        const float128 p = f * float128::ln2();
        const float128 correction = fma(f, float128::ln2(), -p);

        return ldexp(exp(p) * (float128::one() + correction), k);
    }
    /**
     * @brief Calculates the exponent of x and reduces 1 from the result: (e^x) - 1
     * @param x A number specifying a power.
     * @return Exponent of x
     */
    [[nodiscard]] friend FP128_INLINE float128 expm1(const float128& x) noexcept
    {
        if (x.is_nan() || x.is_zero())
            return x;
        if (x.is_inf())
            return x.is_negative() ? -float128::one() : x;

        // Half is the point where exp(x) - 1 stops throwing away significant bits: below it the
        // subtraction cancels most of what exp() produced, and at 2^-60 there is nothing left of
        // the answer at all. The Maclaurin series has no cancellation - every term is formed from
        // the argument itself - and needs 30 terms to reach the last bit at this magnitude.
        if (fabs(x) <= float128::half()) {
            constexpr int32_t terms = 30;
            float128 acc, factorial;
            fact_reciprocal(terms + 1, acc);
            for (int32_t n = terms - 1; n >= 0; --n) {
                fact_reciprocal(n + 1, factorial);
                acc = factorial + x * acc;
            }
            return x * acc;
        }

        return exp(x) - float128::one();
    }
    /**
     * @brief Computes x to the power of y
     * @param x Base value
     * @param y Exponent value (integer)
     * @return x^y
     */
    [[nodiscard]] friend FP128_INLINE float128 pow(const float128& x, int32_t y) noexcept
    {
        static const float128 max_exponent = 11355;  // log(16382) / log2
        float128 res = 1;
        // check the trivial cases
        if (y == 1) {
            return x;
        } else if (y == 0) {
            return 1;
        } else if (x == 1) {
            return x;
        }
        // The magnitude is taken in the unsigned domain: negating the most negative int32_t is
        // undefined, and the name abs now resolves to this namespace's float128 overload rather
        // than to the one from <cstdlib>.
        uint32_t expo = (y < 0) ? (0u - static_cast<uint32_t>(y)) : static_cast<uint32_t>(y);
        // check if the value isn't too large
        if (y > max_exponent) {
            res = inf();
        }
        // check if the value isn't too small
        else if (y < -max_exponent)
            res = 0;
        // compute x^y
        else if (expo > 0) {
            float128 b = x;  // value of e^1
            while (expo > 0) {
                if (expo & 1)
                    res *= b;
                expo >>= 1;
                b.square();
            }
        }

        return (y >= 0) ? res : reciprocal(res);
    }
    /**
     * @brief Computes x to the power of y
     * @param x Base value
     * @param y Exponent value
     * @return x^y
     */
    [[nodiscard]] friend FP128_INLINE float128 pow(const float128& x, const float128& y) noexcept
    {
        //
        // Based on exponent law: (x^n)^m = x^(m * n)
        // Convert the exponent y (function parameter) to produce an exponent that will work with exp()
        // z = log(x)
        // pow(x, y) = x^y = e^(y * z) = exp(y * z)
        //
        if (y.is_int()) {
            return pow(x, static_cast<int32_t>(y));
        } else if (x.is_negative()) {
            return -nan();
        }

        float128 lan_x = log(x);
        if (!lan_x)
            return lan_x;

        return exp(y * lan_x);
    }
    /**
     * @brief Calculates the natural Log (base e) of x: log(x)
     * @param x The number to perform log on.
     * @return log(x)
     */
    /// @brief Largest |t| the log1p_small() series is used for. Chosen so twelve terms suffice.
    [[nodiscard]] FP128_FORCE_INLINE static constexpr float128 log1p_small_limit() noexcept { return float128(0, 0, EXP_BIAS - 4, 0); }
    /// @brief Terms of the log1p_small() series. |s| stays below 2^-4.9, so s^24 is past the last bit.
    static constexpr int32_t LOG1P_TERMS = 12;

    /**
     * @brief Natural logarithm of 1+t for |t| <= 2^-4, accurate relative to the result.
     *
     * log(1+t) = 2 * atanh(t / (2+t)). The substitution replaces the Mercator series in t, which
     * would need a hundred terms at this precision, with one in s^2 where |s| < 2^-4.9, so twelve
     * terms reach the last bit of the mantissa.
     *
     * What matters more than the term count is that the result stays accurate relative to its own
     * size. Near one the logarithm is proportional to t, and t reaches this function exactly - the
     * subtraction that produced it cancelled two values within a factor of two of each other, which
     * is exact. Computing the same thing as a difference of two numbers near one, which is what
     * log2()'s argument reduction does, loses a bit for every power of two the argument sits closer
     * to one: at t = 2^-60 barely half the mantissa survives.
     *
     * @param t Offset from one, must satisfy |t| <= 2^-4
     * @return log(1+t)
     */
    [[nodiscard]] FP128_INLINE static float128 log1p_small(const float128& t) noexcept
    {
        if (t.is_zero())
            return t;

        const float128 s = t / (float128(2) + t);
        const float128 w = sqr(s);

        // 1/(2k+1) for the Horner pass below. Built once: the divisions are far more expensive
        // than the series itself, and the values do not depend on the argument.
        static const auto odd_reciprocal = [] {
            std::array<float128, LOG1P_TERMS + 1> table {};
            for (int32_t k = 0; k <= LOG1P_TERMS; ++k)
                table[static_cast<size_t>(k)] = float128::one() / float128(2 * k + 1);
            return table;
        }();

        float128 acc = odd_reciprocal[LOG1P_TERMS];
        for (int32_t k = LOG1P_TERMS - 1; k >= 0; --k)
            acc = odd_reciprocal[static_cast<size_t>(k)] + w * acc;

        return (s * acc) << 1;
    }
    [[nodiscard]] friend FP128_INLINE float128 log(float128 x) noexcept
    {
        // Sterbenz guarantees the subtraction is exact over the range the branch covers, so the
        // series below sees the offset from one with every bit it has.
        const float128 t = x - float128::one();
        if (fabs(t) <= log1p_small_limit())
            return log1p_small(t);

        return log2(x) * float128::ln2();
    }
    /**
     * @brief Calculates the Log base 2 of x: y = log2(x)
     *
     * Accurate relative to its own result, which matters most for an argument close to one: the
     * answer is then proportional to the mantissa's fraction, and that fraction is carried through
     * exactly rather than being rounded onto a grid it would barely register on. The earlier
     * implementation built the answer as a fixed point fraction of 112 bits before returning it as
     * a float, which left log2(1 + 2^-112) with no correct significant bits at all.
     *
     * @param x The number to perform log2 on.
     * @return log2(x)
     */
    [[nodiscard]] friend FP128_INLINE float128 log2(float128 x) noexcept
    {
        if (x.is_negative() || x.is_zero()) {
            return -inf();
        }

        // Near one the reduction below forms the answer as a difference of two values close to
        // one, which costs it a bit for every power of two the argument sits closer to one. The
        // series has no such cancellation, see log1p_small().
        const float128 offset = x - float128::one();
        if (fabs(offset) <= log1p_small_limit())
            return log1p_small(offset) * float128::log2_e();

        // Calculate the log in 2 steps:
        // - The integer part is simple and fast via the get_exponent() function.
        // - The fraction part is log2 of the mantissa, by argument reduction and a short series.
        // The result is the sum of the two. Based on the identity:
        // log(x + y) = log(x) + log(y)
        const int32_t expo = x.get_exponent();

        // x is an exponent of 2, so the mantissa contributes nothing. This also keeps the earlier
        // handling of infinity, whose fraction bits are zero.
        if (x.is_exponent_of_2()) {
            return float128(expo);
        }

        uint64_t frac_low = 0, frac_high = 0;
        int32_t mantissa_expo = 0;
        uint32_t mantissa_sign = 0;
        x.get_components(frac_low, frac_high, mantissa_expo, mantissa_sign);
        frac_high &= UPPER_FRAC_MASK;  // drop the unity bit, leaving f in [0,1) scaled by 2^112

        // The leading log2_reduction_bits fraction bits choose the reciprocal to reduce with.
        const size_t j = static_cast<size_t>(frac_high >> (FRAC_BITS - 64 - log2_reduction_bits));

        // z, the reduced argument, as a 128 bit fraction. The series below needs |z| <= 2^-6.
        uint64_t z_low = 0, z_high = 0;
        uint32_t z_sign = 0;
        if (j == 0) {
            // Already inside the series' range, so no reduction is applied - and none is wanted.
            // f is exact here, and multiplying it by a rounded reciprocal would cost it every
            // significant bit it has when it is small. That is what makes log2 near 1 accurate:
            // the answer is then proportional to f, and f survives intact.
            z_low = frac_low << 16;
            z_high = (frac_high << 16) | (frac_low >> 48);
        } else {
            // mantissa/2 as a 128 bit fraction, which is in [0.5,1)
            const uint64_t m_low = frac_low << 15;
            const uint64_t m_high = (1ull << 63) | (frac_high << 15) | (frac_low >> 49);
            uint64_t p_low = 0, p_high = 0;
            mul128_high(m_low, m_high, log2_recip_table[j][1], log2_recip_table[j][0], p_low, p_high);

            // z = 2 * (mantissa/2 * recip) - 1. Normally non negative; a reciprocal that rounded
            // down can take it just below zero.
            constexpr uint64_t one_half = 1ull << 63;
            uint64_t d_low = 0, d_high = 0;
            if (p_high >= one_half) {
                d_low = p_low;
                d_high = p_high - one_half;
            } else {
                d_low = 0ull - p_low;
                d_high = one_half - p_high - ((p_low != 0) ? 1ull : 0ull);
                z_sign = 1;
            }
            z_low = d_low << 1;
            z_high = (d_high << 1) | (d_low >> 63);
        }

        // Horner over 1/(n*ln2), from the last term down, giving A = log2(1+z)/z.
        //
        // Every value here is halved relative to the mathematics - log2_inv_n_table already holds
        // 1/(2n*ln2) - so that all of them stay below one and can be held as plain 128 bit
        // fractions. A is about 1.44 and would not fit otherwise. The accumulator stays inside
        // [0.019, 0.734] throughout, so neither the subtraction nor the addition can leave range.
        uint64_t acc_low = log2_inv_n_table[LOG2_TERMS - 1][1];
        uint64_t acc_high = log2_inv_n_table[LOG2_TERMS - 1][0];
        for (int32_t n = LOG2_TERMS - 1; n >= 1; --n) {
            uint64_t term_low = 0, term_high = 0;
            mul128_high(z_low, z_high, acc_low, acc_high, term_low, term_high);
            const uint64_t q_low = log2_inv_n_table[n - 1][1], q_high = log2_inv_n_table[n - 1][0];
            if (z_sign) {
                const uint8_t carry = addcarryx_u64(0, q_low, term_low, &acc_low);
                addcarryx_u64(carry, q_high, term_high, &acc_high);
            } else {
                const uint8_t borrow = subborrow_u64(0, q_low, term_low, &acc_low);
                subborrow_u64(borrow, q_high, term_high, &acc_high);
            }
        }

        // The product is formed in float128 rather than in the fraction arithmetic above. A fraction
        // grid is absolute, and log2(1+z) is proportional to z, so holding the product on that grid
        // would leave a small result with only as many significant bits as it has room above 2^-128.
        // Multiplying an exact z by an A that is accurate in its own right keeps the answer accurate
        // relative to its own size, which is what a floating point type is expected to do.
        float128 series = from_fraction128(z_low, z_high) * from_fraction128(acc_low, acc_high);
        series <<= 1;  // undo the halving the table carries
        series.set_sign(z_sign);

        // -log2(recip), the part of the answer the reduction removed. Zero when nothing was reduced.
        const float128 table_value = (j != 0) ? from_fraction128(log2_value_table[j][1], log2_value_table[j][0]) : float128();

        // The two fraction parts are summed before the exponent is added in. Both are below one, so
        // that first addition rounds against a small value; adding the exponent first would round
        // twice against a number as large as 16000 and cost the answer a bit for nothing.
        return float128(expo) + (table_value + series);
    }
    /**
     * @brief Calculates Log base 10 of x: log10(x)
     * @param x The number to perform log on.
     * @return log10(x)
     */
    [[nodiscard]] friend FP128_INLINE float128 log10(float128 x) noexcept
    {
        const float128 t = x - float128::one();
        if (fabs(t) <= log1p_small_limit())
            return log1p_small(t) * float128::log10_e();

        return log2(x) * float128::log10_2();
    }
    /**
     * @brief Calculates Log base 2 of x as an integer ignoring the sign of x.
     * Similar to: floor(log2(fabs(x)))
     * @param x The number to perform log on.
     * @return logb(x)
     */
    [[nodiscard]] friend FP128_INLINE float128 logb(float128 x) noexcept
    {
        if (x.is_nan() || x.is_inf())
            return fabs(x);
        if (x.is_zero())
            return -inf();

        // get_components() normalizes a subnormal, whose stored exponent field is the same for
        // every one of them and therefore says nothing about the value.
        uint64_t l = 0, h = 0;
        int32_t expo = 0;
        uint32_t sign = 0;
        x.get_components(l, h, expo, sign);
        return float128(expo);
    }
    /**
     * @brief Calculates the natural Log (base e) of 1 + x: log(1 + x)
     * @param x The number to perform log on.
     * @return log1p(x)
     */
    [[nodiscard]] friend FP128_INLINE float128 log1p(float128 x) noexcept
    {
        if (x.is_nan() || x.is_zero() || x.is_inf())
            return (x.is_inf() && x.is_negative()) ? nan() : x;

        const float128 one_value = float128::one();
        if (x < -one_value)
            return nan();
        if (x == -one_value)
            return -inf();

        if (fabs(x) <= log1p_small_limit())
            return log1p_small(x);

        // Above the series' range but still below one, 1+x rounds and the logarithm of a value
        // that close to one magnifies the error. Both subtractions below are exact - the first by
        // Sterbenz, the second because the operands agree to within a rounding - so the correction
        // puts back exactly what the rounding of 1+x took away.
        if (fabs(x) < one_value) {
            const float128 sum = one_value + x;
            const float128 rounding = (sum - one_value) - x;
            return log(sum) - rounding / sum;
        }

        return log(one_value + x);
    }

    //
    // Trigonometric functions
    //

    /**
     * @brief Calculate the sine function over a limited range [-0.5pi, 0.5pi]
     * Using the Maclaurin series expansion, the formula is:
     *              x^3   x^5   x^7
     * sin(x) = x - --- + --- - --- + ...
     *               3!    5!    7!
     * @param x value in Radians in the range [-0.5pi, 0.5pi]
     * @return Sine of x
     */
    [[nodiscard]] friend FP128_INLINE float128 sin1(float128 x) noexcept
    {
        assert(fabs(x) <= float128::half_pi());

        // first part of the series is just 'x'
        const float128 xx = x * x;
        float128 elem_denom, elem_nom = x;

        // compute the rest of the series, starting with: -(x^3 / 3!)
        for (int i = 3, sign = 1;; i += 2, sign = 1 - sign) {
            elem_nom *= xx;
            fact_reciprocal(i, elem_denom);
            float128 elem = elem_nom * elem_denom;  // next element in the series
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
    [[nodiscard]] friend FP128_INLINE float128 cos1(const float128& x) noexcept
    {
        constexpr float128 half_pi = float128::half_pi();
        assert(fabs(x) <= half_pi);
        return (x.is_positive()) ? sin1(half_pi - x) : -sin1(-half_pi - x);
    }
    /**
     * @brief Calculate the Sine function
     * Ultimately uses sin() with a reduced range of [-pi/4, pi/4]
     * @param x value in Radians
     * @return Sine of x
     */
    /**
     * @name Unevaluated sums
     *
     * A value held as a pair of float128 whose sum is the quantity meant, with the second one
     * below the last place of the first. That is 226 bits of significand out of two 113 bit
     * values, which is what the trigonometric argument reduction needs: the residual it produces
     * has to be known to the absolute precision of the *result*, and near a zero of the sine the
     * result is many orders of magnitude below the argument it came from.
     * @{
     */

    /**
     * @brief Exact sum of two values, as a rounded sum and the error it made.
     *
     * Dekker's algorithm. Nothing here is an approximation: a + b is exactly sum + err for any two
     * finite operands.
     *
     * @param a First operand
     * @param b Second operand
     * @param sum Receives the rounded sum
     * @param err Receives a + b - sum, exactly
     */
    FP128_INLINE static void two_sum(const float128& a, const float128& b, float128& sum, float128& err) noexcept
    {
        sum = a + b;
        const float128 b_part = sum - a;
        err = (a - (sum - b_part)) + (b - b_part);
    }

    /**
     * @brief Exact sum of two values whose magnitudes are ordered.
     * @param a First operand, must not be smaller in magnitude than b
     * @param b Second operand
     * @param sum Receives the rounded sum
     * @param err Receives a + b - sum, exactly
     */
    FP128_INLINE static void quick_two_sum(const float128& a, const float128& b, float128& sum, float128& err) noexcept
    {
        sum = a + b;
        err = b - (sum - a);
    }
    /// @}

    /**
     * @brief Nearest integer to x / (pi/2), the quadrant the argument falls in.
     *
     * Beyond 2^62 the quadrant no longer fits in the integer the reduction counts with, and an
     * argument that large has fewer bits below the binary point than a quadrant is wide - its
     * sine is not determined by the value in any useful sense. Everything up to there reduces
     * exactly, see reduce_half_pi().
     *
     * @param x Argument
     * @return Quadrant index, zero for an argument too large to reduce.
     */
    [[nodiscard]] FP128_INLINE static int64_t quadrant_of(const float128& x) noexcept
    {
        constexpr float128 two_over_pi(0x2A53F84EAFA3EA6A, 0x45F306DC9C88, 0x3FFE, 0);
        const float128 scaled = x * two_over_pi;
        if (fabs(scaled) >= ldexp(float128::one(), 62))
            return 0;
        return llround(scaled);
    }

    /**
     * @brief Subtracts n*(pi/2) from x, leaving the residual as an unevaluated sum.
     *
     * pi/2 is irrational, so no single float128 multiple of it can be subtracted without leaving
     * an error of n ulp behind - and the sine of the residual is only as accurate as the residual
     * itself. Where the answer is small, which is exactly where an argument lands near a multiple
     * of pi/2, that error is the whole answer: sin() used to return a value with 55 correct bits
     * for an argument a few ulp away from pi.
     *
     * Two things fix it. pi/2 is carried in three pieces covering 342 bits rather than 113, and
     * each product n*piece is split into its rounded value and its exact error with fma(), so no
     * multiplication is lost either. The residual accumulates as a pair, giving it 226 bits.
     *
     * @param x Argument
     * @param n Multiple of pi/2 to remove
     * @param hi Receives the residual
     * @param lo Receives what did not fit in hi
     */
    FP128_INLINE static void reduce_half_pi(const float128& x, int64_t n, float128& hi, float128& lo) noexcept
    {
        constexpr float128 half_pi_parts[3] = {
            float128(0x8469898CC51701B8, 0x921FB54442D1, 0x3FFF, 0),
            float128(0xA67CC74020BBEA64, 0xCD129024E088, 0x3F8C, 0),
            float128(0xE19C72FEC8841ABA, 0x3B19376BAD7D, 0x3F1A, 1),
        };

        hi = x;
        lo = float128();
        if (n == 0)
            return;

        const float128 nf(n);
        for (const float128& part : half_pi_parts) {
            // n * part is exactly product + product_err
            const float128 product = nf * part;
            const float128 product_err = fma(nf, part, -product);

            float128 sum, err;
            two_sum(hi, -product, sum, err);
            // Everything that did not fit is collected and folded back in, which keeps the pair
            // normalized for the next piece.
            quick_two_sum(sum, (err + lo) - product_err, hi, lo);
        }
    }

    [[nodiscard]] friend float128 sin(float128 x) noexcept
    {
        if (x.is_nan() || x.is_zero())
            return x;
        if (x.is_inf())
            return nan();

        const int64_t n = quadrant_of(x);
        float128 hi, lo;
        reduce_half_pi(x, n, hi, lo);

        switch (n & 3) {
        case 0:  // [-45-45) degrees
            return sin1(hi) + cos1(hi) * lo;
        case 1:  // [45-135) degrees
            return cos1(hi) - sin1(hi) * lo;
        case 2:  // [135-225) degrees
            return -(sin1(hi) + cos1(hi) * lo);
        case 3:  // [225-315) degrees
        default:
            return -(cos1(hi) - sin1(hi) * lo);
        }
    }
    /**
     * @brief Calculate the inverse sine function
     * Uses Newton's method to converge quickly.
     * @param x value in radians in the range [-1,1]
     * @return Inverse sine of x
     */
    [[nodiscard]] friend float128 asin(float128 x) noexcept
    {
        constexpr int max_iterations = 6;
        if (x < -1 || x > 1)
            return 0;

        //              sin(Xn) - a
        // Xn+1 = Xn - -------------
        //                cos(Xn)
        // where 'a' is the argument, each iteration will converge on the result if the initial
        //  estimate is close enough.
        auto sign = x.get_sign();
        x.set_sign(0);

        // initial estimate using the standard library
        float128 res = ::asin(static_cast<double>(x));
        const float128 eps = fabs(res >> 110);
        for (int i = 0; i < max_iterations; ++i) {
            float128 e = (sin(res) - x) / cos(res);
            res -= e;
            if (fabs(e) <= eps)
                break;
        }

        res.set_sign(sign);
        return res;
    }
    /**
     * @brief Calculate the cosine function
     * Ultimately uses sin1() with a reduced range of [-pi/4, pi/4]
     * Sine's Maclaurin series converges faster than Cosine's.
     * @param x value in Radians
     * @return Cosine of x
     */
    [[nodiscard]] friend float128 cos(float128 x) noexcept
    {
        if (x.is_nan())
            return x;
        if (x.is_inf())
            return nan();

        const int64_t n = quadrant_of(x);
        float128 hi, lo;
        reduce_half_pi(x, n, hi, lo);

        switch (n & 3) {
        case 0:  // [-45-45) degrees
            return cos1(hi) - sin1(hi) * lo;
        case 1:  // [45-135) degrees
            return -(sin1(hi) + cos1(hi) * lo);
        case 2:  // [135-225) degrees
            return -(cos1(hi) - sin1(hi) * lo);
        case 3:  // [225-315) degrees
        default:
            return sin1(hi) + cos1(hi) * lo;
        }
    }
    /**
     * @brief Calculate the inverse cosine function
     * Uses Newton's method to converge quickly.
     * @param x value in radians in the range [-1,1]
     * @return Inverse cosine of x
     */
    [[nodiscard]] friend float128 acos(float128 x) noexcept
    {
        constexpr int max_iterations = 6;
        if (x < -1 || x > 1)
            return 0;
        //              cos(Xn) - a           a - cos(Xn)
        // Xn+1 = Xn - ------------- = Xn -  ------------
        //                -sin(Xn)              sin(Xn)
        // where 'a' is the argument, each iteration will converge on the result if the initial
        //  estimate is close enough.
        float128 res = ::acos(static_cast<double>(x));
        const float128 eps = fabs(res >> 110);

        for (int i = 0; i < max_iterations; ++i) {
            float128 cos_xn = cos(res);
            float128 sin_xn = sin(res);
            float128 e = (x - cos_xn) / sin_xn;
            res -= e;
            if (fabs(e) <= eps)
                break;
        }

        return res;
    }
    /**
     * @brief Calculate the tangent function
     * tan(x) = sin(x)/cos(x)
     * @param x value
     * @return Tangent of x
     */
    [[nodiscard]] friend FP128_FORCE_INLINE float128 tan(float128 x) noexcept { return sin(x) / cos(x); }
    /**
     * @brief Calculate the inverse tangent function
     * @param x value
     * @return Arctangent of x
     */
    [[nodiscard]] friend float128 atan(float128 x) noexcept
    {
        // constants for segmentation
        constexpr float128 half_pi = float128::half_pi();  // pi / 2
        bool comp = false;
        constexpr int max_iterations = 6;

        // make argument positive, save the sign
        auto sign = x.get_sign();
        x.set_sign(0);

        // limit argument to 0..1
        if (x > 1) {
            comp = true;
            x = reciprocal(x);
        }

        // initial step uses the CRT function.
        float128 res = ::atan(static_cast<double>(x));
        const float128 eps = fabs(res >> 110);

        //
        // Xn+1 =  Xn - cos(Xn) * ( sin(Xn) - a * cos(Xn))
        //
        // where 'a' is the argument, each iteration will converge on the result if the initial
        //  estimate is close enough.
        for (int i = 0; i < max_iterations; ++i) {
            float128 cos_xn = cos(res);
            float128 sin_xn = sin(res);
            float128 e = cos_xn * (sin_xn - x * cos_xn);  // this is the iteration estimated error
            res -= e;
            if (fabs(e) <= eps)
                break;
        }

        // restore complement if needed
        if (comp)
            res = half_pi - res;
        // restore sign if needed
        res.set_sign(sign);
        return res;
    }
    /**
     * @brief Calculate the inverse tangent function of the ratio y / x
     * @param y value
     * @param x value
     * @return Arctangent of y / x in the range [-pi, pi]
     */
    [[nodiscard]] friend float128 atan2(float128 y, float128 x) noexcept
    {
        // constants for segmentation
        constexpr float128 pi = float128::pi();
        constexpr float128 half_pi = float128::half_pi();  // pi / 2

        // x == 0
        if (!x) {
            if (!y)
                return 0;

            return (y.is_negative()) ? -half_pi : half_pi;
        }
        // y == 0
        if (!y)
            return (x.is_negative()) ? -pi : pi;

        float128 res;
        // save the signs of x, y
        bool comp = fabs(y) > fabs(x);
        float128 ratio;

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
    [[nodiscard]] friend FP128_INLINE float128 sinh(const float128 x) noexcept
    {
        if (x.is_nan() || x.is_inf() || x.is_zero())
            return x;

        const float128 a = fabs(x);
        float128 res;
        // Below one, exp(a) and exp(-a) agree to within a factor of e and their difference cancels
        // away most of the mantissa - at a = 2^-60 all of it. expm1() holds the same quantity with
        // no cancellation at all, and t/(t+1) is exp(-a)-1 expressed through it.
        if (a < float128::one()) {
            const float128 t = expm1(a);
            res = (t + t / (t + float128::one())) >> 1;
        } else {
            res = (exp(a) - exp(-a)) >> 1;
        }

        return (x.is_negative()) ? -res : res;
    }
    /**
     * @brief Calculates the inverse hyperbolic sine
     * For positive x:
     * asinh(x) = log(x + sqrt(x^2 + 1))
     * For negative x, the function returns the result with the sign inverted
     * @param x value
     * @return Inverse hyperbolic sine of x
     */
    [[nodiscard]] friend FP128_INLINE float128 asinh(const float128 x) noexcept
    {
        if (x.is_nan() || x.is_inf() || x.is_zero())
            return x;

        const float128 one_value = float128::one();
        const float128 a = fabs(x);
        float128 res;

        // Three regimes, because the closed form is only usable in the middle one.
        //
        // Below 2^-57 the answer is the argument to within half an ulp, and the closed form would
        // square it to zero and hand back log(1) = 0 instead.
        //
        // Above 2^57 the square overflows long before the argument does, and asinh(a) is log(2a)
        // to well past the last bit.
        //
        // In between, adding a to sqrt(a*a+1) cancels for a negative argument, which is why the
        // sign is taken out first, and the a*a/(1+sqrt(1+a*a)) form is used rather than
        // sqrt(a*a+1)-a: near zero the latter subtracts two values that agree to every bit the
        // result needs.
        if (a < ldexp(one_value, -57)) {
            res = a;
        } else if (a > ldexp(one_value, 57)) {
            res = log(a) + float128::ln2();
        } else {
            const float128 a2 = sqr(a);
            res = log1p(a + a2 / (one_value + sqrt(one_value + a2)));
        }

        return (x.is_positive()) ? res : -res;
    }
    /**
     * @brief Calculate the hyperbolic cosine function over a limited range [-0.5pi, 0.5pi]
     *           e^x + e^(-x)
     * cosh(x) = ------------
     *                2
     * @param x value in Radians in the range [-0.5pi, 0.5pi]
     * @return Sine of x
     */
    [[nodiscard]] friend FP128_FORCE_INLINE float128 cosh(const float128 x) noexcept { return (exp(x) + exp(-x)) >> 1; }
    /**
     * @brief Calculates the inverse hyperbolic cosine
     * For x >= 1:
     * acosh(x) = log(x + sqrt(x^2 - 1))
     * For x < 1, the function return zero
     * @param x value in the range [1, inf]
     * @return Inverse hyperbolic cosine of x
     */
    [[nodiscard]] friend FP128_INLINE float128 acosh(const float128 x) noexcept
    {
        if (x.is_nan())
            return x;
        // Outside the domain the result is undefined rather than zero, which is what the previous
        // version returned for every argument below one.
        const float128 one_value = float128::one();
        if (x < one_value)
            return nan();
        if (x.is_inf())
            return x;

        // Near one, x*x-1 cancels: at x = 1+2^-100 the square rounds back to one and the root
        // comes out zero. The argument of log1p below is built from t = x-1, which is exact over
        // that range, so the answer keeps its significant bits.
        if (x < float128(2)) {
            const float128 t = x - one_value;
            return log1p(t + sqrt((t << 1) + sqr(t)));
        }
        // Beyond 2^57 the square overflows before the argument does and acosh(x) is log(2x).
        if (x > ldexp(one_value, 57))
            return log(x) + float128::ln2();

        return log(x + sqrt(sqr(x) - one_value));
    }
    /**
     * @brief Calculates the hyperbolic tangent
     *           e^x - e^(-x)
     * tanh(x) = ------------
     *           e^x + e^(-x)
     * @param x value
     * @return hyperbolic tangent of x
     */
    [[nodiscard]] friend FP128_INLINE float128 tanh(const float128 x) noexcept
    {
        if (x.is_nan() || x.is_zero())
            return x;

        const float128 a = fabs(x);
        float128 res;
        // Past 40 the two exponentials differ by more than the mantissa can express and the
        // quotient is one to the last bit.
        if (x.is_inf() || a > float128(40)) {
            res = float128::one();
        } else {
            //           e^x - e^(-x)         e^(2x) - 1
            // tanh(x) = ------------  =  ----------------
            //           e^x + e^(-x)      (e^(2x) - 1) + 2
            //
            // The right hand form is written in terms of expm1, which holds e^(2x)-1 without the
            // cancellation that subtracting the two exponentials suffers for a small argument.
            const float128 t = expm1(a << 1);
            res = t / (t + float128(2));
        }

        return (x.is_negative()) ? -res : res;
    }
    /**
     * @brief Calculates the inverse hyperbolic tangent
     *                       1 + x
     * atanh(x) = 0.5 * log( -----)
     *                       1 - x
     * @param x value in the range (-1, 1)
     * @return Inverse hyperbolic tangent of x
     */
    [[nodiscard]] friend FP128_INLINE float128 atanh(const float128 x) noexcept
    {
        if (x.is_nan() || x.is_zero())
            return x;

        constexpr auto one_value = float128::one();
        const float128 a = fabs(x);
        // The endpoints are poles and outside them the function is undefined. Returning zero for
        // all three, as this used to, gave atanh(1) a finite value and atanh(2) a defined one.
        if (a > one_value)
            return nan();
        if (a == one_value)
            return x.is_negative() ? -inf() : inf();

        // 2*atanh(a) = log1p(2a/(1-a)). Splitting the argument as below keeps the numerator from
        // rounding away for a small a, where the answer is proportional to it: the previous form,
        // log((1+x)/(1-x)), formed a quotient that rounds to exactly one and then took its
        // logarithm, which is zero.
        float128 res;
        if (a < float128::half())
            res = log1p((a << 1) + (sqr(a) << 1) / (one_value - a)) >> 1;
        else
            res = log1p((a << 1) / (one_value - a)) >> 1;

        return (x.is_negative()) ? -res : res;
    }
    /**
     * @brief Calculates the Maclaurin series constants for the erf function.
     * The array will hold  1 / (n! * (2n + 1))
     * @param a pointer to array that receives the results. The array must be preallocated.
     * @param count Element count in the array
     */
    friend void erf_constants(float128* a, int32_t count)
    {
        if (a == nullptr)
            return;

        a[0] = 1;
        float128 f = 1;  // value of 0!
        for (int32_t i = 1; i < count; ++i) {
            f *= i;
            a[i] = float128::one() / (f * (2 * i + 1));
        }
    }
    /**
     * @brief Computes the error function of a value.
     * The erf function return a value in the range -1.0 to 1.0.
     * There's no error return.
     * @param x A floating-point value.
     * @return The erf functions return the Gauss error function of x.
     */
    [[nodiscard]] friend FP128_INLINE float128 erf(float128 x) noexcept
    {
        if (x.is_nan() || x.is_zero())
            return x;

        const uint32_t sign = x.get_sign();
        const float128 a = fabs(x);
        float128 res;

        // erfc(9) is 4.1e-37, which is below half the last place of one, so erf is one from here
        // on and the series need never be run outside the range it converges quickly over.
        if (x.is_inf() || a >= float128(9)) {
            res = float128::one();
        } else {
            //                        2                   n
            //           2x        -x     inf         (2x^2)
            // erf(x) = ------ * e     *  sum   ------------------
            //          sqrt(pi)          n=0    1*3*5*...*(2n+1)
            //
            // Every term of this series is positive, unlike the Maclaurin series in x^(2n+1) that
            // this replaced: that one alternates, and its terms peak around n = x^2 at a value
            // 2^24 times the sum by the time x reaches 4, so the answer was formed by cancelling
            // away 24 bits of it. Past x = 2 the old code did not use a series at all - it called
            // the CRT's double erf and returned 53 correct bits.
            const float128 xx = sqr(a);
            // The square rounds, and exp() magnifies that by x^2, which reaches 81 here. Keeping
            // the exact error of the multiply and applying it as a factor puts those bits back.
            const float128 xx_err = fma(a, a, -xx);
            const float128 two_xx = xx << 1;

            // The sum runs to eighty terms near the top of the range and every one of them rounds
            // against a running total, which alone would cost the answer seven bits. two_sum()
            // hands back what each addition lost, so it can be carried and added in at the end.
            float128 term = float128::one();
            float128 sum = term;
            float128 lost;
            for (int32_t n = 1; n < 1024; ++n) {
                term = term * two_xx / float128(2 * n + 1);
                float128 rounded, err;
                two_sum(sum, term, rounded, err);
                sum = rounded;
                lost += err;
                if (term.get_exponent() + FRAC_BITS + 8 < sum.get_exponent())
                    break;
            }
            sum += lost;

            const float128 gaussian = exp(-xx) * (float128::one() - xx_err);
            res = ((a * sum) << 1) * float128::inv_sqrt_pi() * gaussian;
        }

        res.set_sign(sign);
        return res;
    }
    /**
     * @brief Computes the complementary error function, 1 - erf(x).
     *
     * Above one the answer is a small difference of two values close to one, so it cannot be
     * reached through erf() - at x = 6 the subtraction would leave 55 of the 113 bits, and past
     * 13 nothing at all. It is computed directly instead, from Laplace's continued fraction
     *
     *                       2
     *                     -x
     *                    e                 1
     *      erfc(x) = ---------- * ------------------------
     *                 sqrt(pi)     x + (1/2)/(x + 1/(x + ...))
     *
     * evaluated with the modified Lentz algorithm. The fraction converges for every x above one,
     * and fastest where the series form is weakest.
     *
     * @param x Input value
     * @return 1 - erf(x)
     */
    [[nodiscard]] friend FP128_INLINE float128 erfc(float128 x) noexcept
    {
        if (x.is_nan())
            return x;
        if (x.is_inf())
            return x.is_negative() ? float128(2) : float128();
        if (x.is_zero())
            return float128::one();

        const float128 one_value = float128::one();
        // erfc is symmetric about one: the negative half is where the answer approaches two and
        // keeps every bit, so it is the reflection that is well conditioned there.
        if (x.is_negative())
            return float128(2) - erfc(-x);
        if (x < one_value)
            return one_value - erf(x);
        // exp(-x*x) reaches the smallest subnormal at 106.9
        if (x > float128(107))
            return float128();

        constexpr int32_t max_iterations = 4096;
        const float128 tiny = ldexp(one_value, -16000);
        const float128 epsilon = ldexp(one_value, -FRAC_BITS - 8);

        float128 f = tiny;
        float128 c = f;
        float128 d;
        for (int32_t n = 1; n < max_iterations; ++n) {
            // a_1 is one, every a after it is (n-1)/2
            const float128 a = (n == 1) ? one_value : (float128(n - 1) >> 1);

            d = x + a * d;
            if (d.is_zero())
                d = tiny;
            c = x + a / c;
            if (c.is_zero())
                c = tiny;
            d = one_value / d;

            const float128 delta = c * d;
            f *= delta;
            if (fabs(delta - one_value) < epsilon)
                break;
        }

        const float128 xx = sqr(x);
        const float128 xx_err = fma(x, x, -xx);
        const float128 gaussian = exp(-xx) * (one_value - xx_err);
        return gaussian * float128::inv_sqrt_pi() * f;
    }
    /**
     * @brief Tests whether x is finite (not infinite and not NaN).
     * @param x Input value
     * @return Non-zero if finite, zero otherwise.
     */
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr bool isfinite(const float128& x) noexcept { return x.is_finite(); }
    /**
     * @brief Computes (x * y) + z with a single rounding at the end.
     *
     * The whole point of a fused multiply add is that the product is not rounded before the
     * addend reaches it, which is what makes the result correct to the last bit and what lets a
     * caller recover the error of a multiplication as fma(x, y, -x*y). Writing it as x * y + z
     * rounds twice and loses everything the addend cancels against.
     *
     * The exact product of the two 113 bit mantissas is 226 bits, and the addend is another 113.
     * Both are placed in a wide accumulator with their leading bits aligned, added or subtracted
     * there, and the single rounding is applied when the result is read back out. An operand too
     * far below the other one to reach the accumulator cannot change any bit of the result except
     * through the tie break, so it is folded into the sticky bit.
     *
     * @param x The first value to multiply.
     * @param y The second value to multiply.
     * @param z The value to add to the product.
     * @return (x * y) + z, correctly rounded.
     */
    [[nodiscard]] friend FP128_INLINE constexpr float128 fma(float128 x, float128 y, float128 z) noexcept
    {
        const uint32_t product_sign = x.get_sign() ^ y.get_sign();

        // A NaN operand propagates, and so does the invalid product of a zero and an infinity.
        if (x.is_nan() || y.is_nan() || z.is_nan())
            return nan();
        if ((x.is_inf() && y.is_zero()) || (y.is_inf() && x.is_zero()))
            return nan();

        // An infinite product decides the result on its own, unless the addend is the opposite
        // infinity and the sum is undefined.
        if (x.is_inf() || y.is_inf()) {
            if (z.is_inf() && z.get_sign() != product_sign)
                return nan();
            float128 res = inf();
            res.set_sign(product_sign);
            return res;
        }
        if (z.is_inf())
            return z;

        // A zero product leaves the addend. Two zeros give a zero that is negative only when both
        // of them are, which is what round to nearest requires of a sum of opposite signs.
        if (x.is_zero() || y.is_zero()) {
            if (!z.is_zero())
                return z;
            float128 res;
            res.set_sign((z.get_sign() == product_sign) ? product_sign : 0);
            return res;
        }
        // A zero addend leaves the product, which the multiply already rounds correctly.
        if (z.is_zero())
            return x * y;

        uint64_t lx = 0, hx = 0, ly = 0, hy = 0, lz = 0, hz = 0;
        int32_t ex = 0, ey = 0, ez = 0;
        uint32_t sx = 0, sy = 0, sz = 0;
        x.get_components(lx, hx, ex, sx);
        y.get_components(ly, hy, ey, sy);
        z.get_components(lz, hz, ez, sz);

        // The exact 226 bit product of the two 113 bit mantissas. Both were normalized to
        // [2^112, 2^113), so the product is in [2^224, 2^226) and its value is
        // product * 2^(ex + ey - 224).
        uint64_t product[4] {};
        mul128to256(lx, hx, ly, hy, product);
        const int32_t product_bits = wide_msb(product, 4);

        // Weight of the leading bit of each operand, which is what they are aligned by.
        const int32_t product_msb = ex + ey - 2 * FRAC_BITS + product_bits;
        const int32_t top = (product_msb > ez) ? product_msb : ez;

        // Bit WIDE_TOP of an accumulator carries the weight 2^top.
        uint64_t product_acc[WIDE_WORDS] {};
        uint64_t addend_acc[WIDE_WORDS] {};
        const uint64_t addend[2] = {lz, hz};
        bool sticky = wide_place(product_acc, product, 4, WIDE_TOP - (top - product_msb) - product_bits);
        sticky = wide_place(addend_acc, addend, 2, WIDE_TOP - (top - ez) - FRAC_BITS) || sticky;

        // Only the operand with the smaller leading bit can have been placed below the
        // accumulator, so a sticky bit always belongs to the smaller of the two.
        uint64_t* acc = product_acc;
        uint32_t sign = product_sign;
        if (product_sign == sz) {
            wide_add(product_acc, addend_acc);
        } else {
            const int32_t cmp = wide_cmp(product_acc, addend_acc);
            // Equal accumulators mean equal values: a dropped tail puts its operand more than 157
            // bits below the other one, which no accumulator of the larger one can match.
            if (cmp == 0)
                return float128();

            if (cmp > 0) {
                wide_sub(product_acc, addend_acc);
            } else {
                wide_sub(addend_acc, product_acc);
                acc = addend_acc;
                sign = sz;
            }
            // What the subtrahend lost below the accumulator makes the difference smaller by less
            // than one of its last places. Borrowing that one place and marking the remainder
            // sticky says exactly that.
            if (sticky)
                wide_dec(acc);
        }

        const int32_t msb = wide_msb(acc, WIDE_WORDS);
        int32_t expo = top - WIDE_TOP + msb;

        // Read the 113 bit mantissa out of the accumulator, with the bit below it and whether
        // anything below that was set.
        uint64_t l = 0, h = 0;
        uint64_t guard = 0;
        const int32_t lsb = msb - FRAC_BITS;
        if (lsb <= 0) {
            // Fewer bits survived than the mantissa holds, so the result is exact. Cancellation
            // that deep cannot coexist with a dropped tail, which needs the operands far apart.
            l = acc[0];
            h = acc[1];
            shift_left128_inplace_safe(l, h, -lsb);
        } else {
            const int32_t word = lsb >> 6;
            const int32_t bit = lsb & 63;
            l = acc[word] >> bit;
            h = (word + 1 < WIDE_WORDS) ? (acc[word + 1] >> bit) : 0;
            if (bit != 0) {
                if (word + 1 < WIDE_WORDS)
                    l |= acc[word + 1] << (64 - bit);
                if (word + 2 < WIDE_WORDS)
                    h |= acc[word + 2] << (64 - bit);
            }
            h &= FRAC_UNITY | UPPER_FRAC_MASK;

            guard = FP128_GET_BIT(acc[(lsb - 1) >> 6], (lsb - 1) & 63);
            for (int32_t i = 0; i < ((lsb - 1) >> 6); ++i)
                sticky = sticky || acc[i] != 0;
            const int32_t below = (lsb - 1) & 63;
            if (below != 0)
                sticky = sticky || (acc[(lsb - 1) >> 6] & FP128_MAX_VALUE_64(below)) != 0;
        }

        // round half to even
        if (guard != 0 && (sticky || (l & 1) != 0)) {
            if (++l == 0)
                ++h;
            // Rounding up out of the mantissa lands on the next power of two.
            if ((h >> EXP_SHIFT) != 1) {
                h >>= 1;
                ++expo;
            }
        }

        float128 res;
        res.set_components(l, h, expo, sign);
        return res;
    }
    /**
     * @brief Gets the mantissa and exponent of a floating-point number.
     * @param x Floating-point value.
     * @param expptr Floating-point value.
     * @return Mantissa in the [0.5,1) range.
     */
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr float128 frexp(float128 x, int* expptr) noexcept
    {
        if (x.is_special() || x.is_zero()) {
            *expptr = 0;
            return x;
        }

        uint64_t l, h;
        int32_t e;
        uint32_t s;
        x.get_components(l, h, e, s);
        *expptr = e + 1;
        float128 res(l, h, EXP_BIAS - 1, s);
        return res;
    }
    /**
     * @brief Multiplies a floating-point number by an integral power of two.
     * @param x Floating-point value.
     * @param exp Integer exponent.
     * @return The ldexp functions return the value of x * 2^exp if successful. On overflow, and depending on the sign of x, ldexp returns +/- inf
     */
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr float128 ldexp(float128 x, int exp) noexcept
    {
        if (x.is_zero())
            return x;
        if (x.is_special())
            return x;

        if (exp > 0)
            return x << exp;
        return x >> -exp;
    }

    /**
     * @name Classification and comparison
     *
     * The <cmath> predicates. They are functions rather than the macros the C header defines, in
     * the same way the standard library's own C++ overloads are, so they take part in overload
     * resolution and can be found by argument dependent lookup.
     * @{
     */

    /// @brief True when the sign bit is set, including for -0 and for a negative NaN.
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr bool signbit(const float128& x) noexcept { return x.get_sign() != 0; }
    /// @brief True when x is neither zero, subnormal, infinite nor NaN.
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr bool isnormal(const float128& x) noexcept { return x.is_normal() && !x.is_zero(); }
    /// @brief One of FP_NAN, FP_INFINITE, FP_ZERO, FP_SUBNORMAL or FP_NORMAL.
    [[nodiscard]] friend FP128_INLINE constexpr int fpclassify(const float128& x) noexcept
    {
        if (x.is_nan())
            return FP_NAN;
        if (x.is_inf())
            return FP_INFINITE;
        if (x.is_zero())
            return FP_ZERO;
        return x.is_subnormal() ? FP_SUBNORMAL : FP_NORMAL;
    }
    /// @brief True when either operand is a NaN, so that the two do not compare in any order.
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr bool isunordered(const float128& x, const float128& y) noexcept
    {
        return x.is_nan() || y.is_nan();
    }
    /// @brief x > y, false rather than unordered when either is a NaN.
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr bool isgreater(const float128& x, const float128& y) noexcept
    {
        return !isunordered(x, y) && x > y;
    }
    /// @brief x >= y, false when either is a NaN.
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr bool isgreaterequal(const float128& x, const float128& y) noexcept
    {
        return !isunordered(x, y) && x >= y;
    }
    /// @brief x < y, false when either is a NaN.
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr bool isless(const float128& x, const float128& y) noexcept
    {
        return !isunordered(x, y) && x < y;
    }
    /// @brief x <= y, false when either is a NaN.
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr bool islessequal(const float128& x, const float128& y) noexcept
    {
        return !isunordered(x, y) && x <= y;
    }
    /// @brief x < y or x > y, false when either is a NaN.
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr bool islessgreater(const float128& x, const float128& y) noexcept
    {
        return !isunordered(x, y) && (x < y || x > y);
    }
    /// @}

    /**
     * @brief A quiet NaN carrying the payload spelled out in the argument.
     *
     * Mirrors the C library's nan(): the string is read as an unsigned integer, decimal by
     * default, hexadecimal after an 0x prefix, and its low bits become the NaN's fraction. An
     * empty or unparsable string gives the default NaN.
     *
     * @param payload Digits of the payload, may be empty
     * @return A quiet NaN.
     */
    [[nodiscard]] friend FP128_INLINE float128 nan(const char* payload) noexcept
    {
        uint64_t bits = 1;
        if (payload != nullptr && *payload != '\0') {
            const uint64_t parsed = strtoull(payload, nullptr, 0);
            // A zero fraction is an infinity rather than a NaN, so the default payload stands in.
            if (parsed != 0)
                bits = parsed;
        }
        // The leading fraction bit marks the NaN quiet; the payload fills the bits below it.
        float128 res = float128(bits, 0, INF_EXP_BIASED, 0);
        res.high |= QUIET_NAN_BIT;
        return res;
    }

    /// @brief True when x is an integer whose last place is one, which is the parity a tie needs.
    [[nodiscard]] FP128_INLINE static constexpr bool is_odd_int(const float128& x) noexcept
    {
        uint64_t l = 0, h = 0;
        int32_t expo = 0;
        uint32_t sign = 0;
        x.get_components(l, h, expo, sign);
        // Below one there is no units bit; above 2^112 the last place is two or more, so every
        // value that far out is even.
        if (expo < 0 || expo > FRAC_BITS)
            return false;
        const int32_t bit = FRAC_BITS - expo;
        return (bit < 64) ? ((l >> bit) & 1) != 0 : ((h >> (bit - 64)) & 1) != 0;
    }

    /**
     * @name Rounding to the current mode
     *
     * The library rounds to nearest with ties to even and offers no way to change that, so
     * nearbyint and rint are the same function and neither can raise the inexact exception rint is
     * otherwise allowed to. round() differs from both: it breaks a tie away from zero.
     * @{
     */
    [[nodiscard]] friend FP128_INLINE constexpr float128 rint(const float128& x) noexcept
    {
        if (x.is_special() || x.is_zero() || x.is_int())
            return x;

        const float128 truncated = trunc(x);
        const float128 fraction = fabs(x - truncated);  // exact
        const float128 half_value = float128::half();
        if (fraction < half_value)
            return truncated;

        float128 step = float128::one();
        step.set_sign(x.get_sign());
        if (fraction > half_value)
            return truncated + step;

        // Exactly halfway: the tie goes to the even neighbour.
        return is_odd_int(truncated) ? truncated + step : truncated;
    }
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr float128 nearbyint(const float128& x) noexcept { return rint(x); }
    /// @}

    /**
     * @name Exponent manipulation
     * @{
     */
    /// @brief x * 2^n, the same operation ldexp performs.
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr float128 scalbn(const float128& x, int n) noexcept { return ldexp(x, n); }
    /// @brief x * 2^n for a long exponent, clamped to what the format can reach.
    [[nodiscard]] friend FP128_INLINE constexpr float128 scalbln(const float128& x, long n) noexcept
    {
        // Anything past the format's range saturates the same way the shift itself would, and
        // clamping keeps the conversion to int well defined.
        constexpr long limit = 2 * (EXP_BIAS + FRAC_BITS);
        const long clamped = (n > limit) ? limit : ((n < -limit) ? -limit : n);
        return ldexp(x, static_cast<int>(clamped));
    }
    /// @}

    /**
     * @brief The representable value next to x in the direction of y.
     * @param x Starting value
     * @param y Value giving the direction
     * @return The neighbour of x towards y, or y itself when the two are equal.
     */
    [[nodiscard]] friend FP128_INLINE constexpr float128 nextafter(const float128& x, const float128& y) noexcept
    {
        if (x.is_nan() || y.is_nan())
            return nan();
        if (x == y)
            return y;  // the sign of y is what the standard hands back here
        return (y > x) ? nextUp(x) : nextDown(x);
    }
    /// @brief The neighbour of x towards y. The same as nextafter: there is no wider type to convert from.
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr float128 nexttoward(const float128& x, const float128& y) noexcept
    {
        return nextafter(x, y);
    }

    /**
     * @brief |x| mod |y|, along with the low bits of the quotient it implies.
     *
     * The shared engine behind fmod, remainder and remquo. All three are exact operations on the
     * mantissas, and all three need the same long division; only what they do with the quotient
     * differs.
     *
     * @param x Dividend, must be finite and non zero
     * @param y Divisor, must be finite and non zero
     * @param quotient Receives the low bits of the integer quotient
     * @return The remainder, with a positive sign.
     */
    [[nodiscard]] FP128_INLINE static float128 fmod_magnitude(const float128& x, const float128& y, uint64_t& quotient) noexcept
    {
        quotient = 0;

        uint64_t lx = 0, hx = 0, ly = 0, hy = 0;
        int32_t ex = 0, ey = 0;
        uint32_t sx = 0, sy = 0;
        x.get_components(lx, hx, ex, sx);
        y.get_components(ly, hy, ey, sy);

        if (ex < ey || (ex == ey && uint128_t(lx, hx) < uint128_t(ly, hy)))
            return fabs(x);

        uint128_t rem(lx, hx);
        const uint128_t mod(ly, hy);
        uint128_t block_quotient = rem / mod;
        rem -= block_quotient * mod;
        uint64_t low = 0, high = 0;
        block_quotient.get_components(low, high);
        quotient = low;

        int32_t shift = ex - ey;
        constexpr int32_t block = 15;
        while (shift > 0) {
            const int32_t step = (shift < block) ? shift : block;
            rem <<= step;
            shift -= step;
            block_quotient = rem / mod;
            rem -= block_quotient * mod;
            block_quotient.get_components(low, high);
            // Only the low bits of the quotient are ever wanted, so the overflow is dropped.
            quotient = (quotient << step) + low;
        }

        if (rem.is_zero())
            return float128();

        uint64_t rl = 0, rh = 0;
        rem.get_components(rl, rh);
        const int32_t msb = static_cast<int32_t>(log2(rl, rh));
        shift_left128_inplace_safe(rl, rh, FRAC_BITS - msb);

        float128 res;
        res.set_components(rl, rh, ey - FRAC_BITS + msb, 0);
        return res;
    }

    /**
     * @brief IEEE 754 remainder: x - n*y with n the quotient rounded to the nearest even integer.
     *
     * Unlike fmod, whose result keeps the sign of x, this one lands in [-|y|/2, |y|/2].
     *
     * @param x Dividend
     * @param y Divisor
     * @return The remainder.
     */
    [[nodiscard]] friend FP128_INLINE float128 remainder(const float128& x, const float128& y) noexcept
    {
        int quotient = 0;
        return remquo(x, y, &quotient);
    }

    /**
     * @brief The IEEE remainder together with the low bits of the quotient.
     * @param x Dividend
     * @param y Divisor
     * @param quo Receives at least the low three bits of the quotient, with the sign of x/y
     * @return The remainder.
     */
    [[nodiscard]] friend FP128_INLINE float128 remquo(const float128& x, const float128& y, int* quo) noexcept
    {
        if (quo != nullptr)
            *quo = 0;
        if (x.is_nan() || y.is_nan() || x.is_inf() || y.is_zero())
            return nan();
        if (y.is_inf() || x.is_zero())
            return x;

        uint64_t bits = 0;
        const float128 magnitude = fabs(y);
        float128 rem = fmod_magnitude(x, y, bits);

        // fmod leaves the remainder in [0, |y|); the IEEE one is the nearer of that and
        // |y| - it, with a tie going to whichever leaves an even quotient.
        const float128 twice = rem << 1;
        const bool round_up = (twice > magnitude) || (twice == magnitude && (bits & 1) != 0);
        if (round_up) {
            rem -= magnitude;
            ++bits;
        }

        if (quo != nullptr) {
            // The standard asks for the low bits of the magnitude of the quotient, signed by the
            // signs of the operands. Seven bits is the customary amount and fits any int.
            const int low_bits = static_cast<int>(bits & 0x7F);
            *quo = (x.get_sign() != y.get_sign()) ? -low_bits : low_bits;
        }

        rem.set_sign(rem.is_zero() ? x.get_sign() : (x.get_sign() ^ rem.get_sign()));
        return rem;
    }

    /**
     * @brief Length of the diagonal of a box, computed without an intermediate overflow.
     * @param x First side
     * @param y Second side
     * @param z Third side
     * @return sqrt(x*x + y*y + z*z)
     */
    [[nodiscard]] friend FP128_INLINE float128 hypot(const float128& x, const float128& y, const float128& z) noexcept
    {
        if (x.is_inf() || y.is_inf() || z.is_inf())
            return inf();
        if (x.is_nan() || y.is_nan() || z.is_nan())
            return nan();

        // Scaling by the largest term keeps every square inside the format's range, which squaring
        // the values as they came would not: a side above 2^8192 overflows on its own.
        const float128 largest = fmax(fabs(x), fmax(fabs(y), fabs(z)));
        if (largest.is_zero())
            return largest;

        const int32_t expo = ilogb(largest);
        const float128 a = ldexp(x, -expo);
        const float128 b = ldexp(y, -expo);
        const float128 c = ldexp(z, -expo);
        return ldexp(sqrt(sqr(a) + sqr(b) + sqr(c)), expo);
    }

    /// @brief Terms of the Stirling series lgamma() runs. Sixteen reach 2^-119 for an argument of 40.
    static constexpr int32_t STIRLING_TERMS = 16;

    /**
     * @brief Natural logarithm of the absolute value of the gamma function.
     *
     * Stirling's asymptotic series converges usefully only for a large argument, so a small one is
     * walked up to 40 with the recurrence gamma(x+1) = x*gamma(x) and the logarithms of the
     * factors are subtracted off afterwards. A negative argument goes through the reflection
     * formula, which is why the result is the logarithm of the absolute value: the gamma function
     * alternates sign between the poles at the non positive integers.
     *
     * @param x Input value
     * @return log(|gamma(x)|)
     */
    [[nodiscard]] friend FP128_INLINE float128 lgamma(float128 x) noexcept
    {
        if (x.is_nan())
            return x;
        if (x.is_inf())
            return inf();
        // The non positive integers are the poles of the gamma function
        if (x.is_zero() || (x.is_negative() && x.is_int()))
            return inf();
        // gamma(1) and gamma(2) are both one. The recurrence and the series below would reach the
        // logarithm of 31! and subtract it from itself, which cancels to a small non zero value
        // rather than to the exact answer.
        if (x == float128::one() || x == float128(2))
            return float128();

        constexpr float128 coefficients[STIRLING_TERMS] = {
            float128(0x5555555555555555, 0x555555555555, 0x3FFB, 0),  // B2 / (2*1)
            float128(0x6C16C16C16C16C17, 0x6C16C16C16C1, 0x3FF6, 1),  // B4 / (4*3)
            float128(0xA01A01A01A01A01A, 0xA01A01A01A01, 0x3FF4, 0),  // B6 / (6*5)
            float128(0x3813813813813814, 0x381381381381, 0x3FF4, 1),  // B8 / (8*7)
            float128(0x3570EA73806E5479, 0xB951E2B18FF2, 0x3FF4, 0),  // B10 / (10*9)
            float128(0xC81F6AB0D9993C7D, 0xF6AB0D9993C7, 0x3FF5, 1),  // B12 / (12*11)
            float128(0xA41A41A41A41A41A, 0xA41A41A41A41, 0x3FF7, 0),  // B14 / (14*13)
            float128(0x7DC2064A8ED3175C, 0xE4286CB0F539, 0x3FF9, 1),  // B16 / (16*15)
            float128(0xFFA1876FE96381E0, 0x6FE96381E067, 0x3FFC, 0),  // B18 / (18*17)
            float128(0x9EDBDB9CE625987D, 0x6476701181F3, 0x3FFF, 1),  // B20 / (20*19)
            float128(0x5A74F53910C8B380, 0xACE44322CE00, 0x4002, 0),  // B22 / (22*21)
            float128(0xAAB67EE25D73C0F9, 0x39B2525CCCC1, 0x4006, 1),  // B24 / (24*23)
            float128(0x1B4E81B4E81B4E82, 0x12234E81B4E8, 0x400A, 0),  // B26 / (26*25)
            float128(0x7EB3FEDDD8496920, 0x1A198AE1C4AB, 0x400E, 1),  // B28 / (28*27)
            float128(0xA38433DC9FB888D4, 0x51A2089A6E11, 0x4012, 0),  // B30 / (30*29)
            float128(0x77880C2D3577880C, 0xD1089B142D35, 0x4016, 1),  // B32 / (32*31)
        };
        constexpr float128 half_log_two_pi(0x4A69297920028832, 0xD67F1C864BEB, 0x3FFE, 0);
        constexpr float128 stirling_limit(0, 0, EXP_BIAS + 5, 0);  // 32, where sixteen terms suffice

        // gamma(x) * gamma(1-x) = pi / sin(pi*x) reflects the negative half onto the positive one,
        // where the series lives.
        if (x.is_negative()) {
            const float128 reflected = sin(float128::pi() * x);
            return log(float128::pi() / fabs(reflected)) - lgamma(float128::one() - x);
        }

        // Walk up to where the asymptotic series is accurate, remembering what to divide out.
        float128 scale_log;
        while (x < stirling_limit) {
            scale_log += log(x);
            x += float128::one();
        }

        const float128 inv_xx = float128::one() / sqr(x);
        float128 series;
        for (int32_t i = STIRLING_TERMS - 1; i >= 1; --i)
            series = (series + coefficients[i]) * inv_xx;
        series = (series + coefficients[0]) / x;

        return (x - float128::half()) * log(x) - x + half_log_two_pi + series - scale_log;
    }

    /**
     * @brief The gamma function.
     *
     * Exact for the small integer arguments, where the answer is a factorial the format holds
     * without rounding, and exp(lgamma) elsewhere. The exponential is what limits the accuracy: a
     * logarithm as large as 7000 carries its own last bit into the seventh place of the result.
     *
     * @param x Input value
     * @return gamma(x)
     */
    [[nodiscard]] friend FP128_INLINE float128 tgamma(float128 x) noexcept
    {
        if (x.is_nan())
            return x;
        if (x.is_inf())
            return x.is_negative() ? nan() : x;
        // The poles, and the two zeros of 1/gamma that a signed zero argument picks out
        if (x.is_zero())
            return x.is_negative() ? -inf() : inf();
        if (x.is_negative() && x.is_int())
            return nan();
        // gamma(1756) overflows
        if (x > float128(1756))
            return inf();

        // 34! is the largest factorial a binary128 holds exactly, so every integer argument up to
        // 35 comes back with no rounding at all.
        if (x.is_int() && x <= float128(35)) {
            float128 result = float128::one();
            for (float128 i = float128(2); i < x; i += float128::one())
                result *= i;
            return result;
        }

        const float128 magnitude = exp(lgamma(x));
        if (x.is_positive())
            return magnitude;

        // Between two poles the gamma function keeps one sign, alternating with every step: it is
        // negative on (-1, 0), positive on (-2, -1), and so on, which is the parity of floor(x).
        const float128 floor_x = floor(x);
        return is_odd_int(floor_x) ? -magnitude : magnitude;
    }

    /// @brief Absolute value, the name <cmath> gives the floating point overload alongside fabs.
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr float128 abs(const float128& x) noexcept { return fabs(x); }

    /// @brief User-defined literal implementation for constructing float128 from a string.
    friend float128 operator""_f128(const char* literal) { return float128(literal); }
};

static_assert(sizeof(float128) == sizeof(uint64_t) * 2);

/***********************************************************************************
 *                        Text conversion and standard library integration
 ************************************************************************************/

namespace detail
{
/// @brief Default precision the e, f and g styles use when none is given, as in printf.
inline constexpr int32_t DEFAULT_PRECISION = 6;

/// @brief Appends the decimal exponent of a scientific form, with a sign and at least two digits.
inline void append_exponent(std::string& out, int32_t exponent, char marker)
{
    out += marker;
    out += (exponent < 0) ? '-' : '+';
    // The magnitude is taken in the unsigned domain, where the negation wraps. Negating the signed
    // value is undefined for the most negative one and produces the same bit pattern for the rest.
    uint32_t magnitude = (exponent < 0) ? (0u - static_cast<uint32_t>(exponent)) : static_cast<uint32_t>(exponent);

    // Ten digits is the widest a uint32_t can be. No exponent this function is handed comes close -
    // the largest is 16494, from the hexadecimal form of the smallest subnormal - but the buffer is
    // sized for the type rather than for the callers, so the loops below cannot run past it.
    //
    // The digits are produced least significant first and so are written from the end of the
    // buffer backwards, which leaves them already in order and lets the whole run be appended at
    // once. Filling forwards and reversing afterwards needs three loops over one shared counter,
    // and MSVC's code analysis cannot see that the counter stays in range across all of them.
    constexpr int32_t max_digits = 10;
    char digits[max_digits];
    int32_t index = max_digits;
    while (index > 0) {
        digits[--index] = static_cast<char>('0' + (magnitude % 10));
        magnitude /= 10;
        if (magnitude == 0)
            break;
    }

    // The exponent always shows at least two digits, as it does for a double.
    while (index > max_digits - 2)
        digits[--index] = '0';

    out.append(digits + index, static_cast<size_t>(max_digits - index));
}

/**
 * @brief Fewest significant digits that read back as the same value.
 *
 * The search is a bisection rather than a scan because the property is monotone: a decimal with
 * more digits is at least as close to the value as one with fewer, so if some length reads back
 * correctly then every longer one does too.
 *
 * @param low Low QWORD of the mantissa
 * @param high High QWORD of the mantissa
 * @param exponent Unbiased exponent of the value
 * @return Digit count in [1, 36].
 */
[[nodiscard]] inline int32_t shortest_digit_count(uint64_t low, uint64_t high, int32_t exponent)
{
    char digits[40];
    int32_t lower = 1;
    int32_t upper = 36;  // max_digits10, which always reads back

    while (lower < upper) {
        const int32_t middle = (lower + upper) / 2;
        const int32_t exponent10 = to_decimal_digits(low, high, exponent, middle, digits);

        uint64_t back_low = 0, back_high = 0;
        int32_t back_exponent = 0;
        const bool parsed = from_decimal_digits(digits, middle, exponent10 - middle, back_low, back_high, back_exponent);
        if (parsed && back_low == low && back_high == high && back_exponent == exponent)
            upper = middle;
        else
            lower = middle + 1;
    }
    return lower;
}

/**
 * @brief Renders a float128 as text.
 *
 * The one place the character form of a value is produced. std::format, the stream inserter,
 * to_chars and the conversion to std::string all come through here so that they cannot disagree
 * with each other.
 *
 * @param value Value to render
 * @param style One of 'a', 'e', 'f' or 'g', or '\0' for the shortest form that reads back exactly
 * @param precision Digits after the point, or significant digits for 'g'. Negative selects the
 *        default for the style
 * @param uppercase Emit the digits, the exponent marker and the special values in upper case
 * @param alternate Keep the decimal point and, for 'g', the trailing zeros
 * @param sign_char Character to emit for a non negative value: '-' means emit nothing
 * @return The rendered text.
 */
[[nodiscard]] inline std::string render(const float128& value, char style, int32_t precision, bool uppercase, bool alternate, char sign_char)
{
    std::string out;
    if (value.is_negative())
        out += '-';
    else if (sign_char != '-')
        out += sign_char;

    if (value.is_nan()) {
        out += uppercase ? "NAN" : "nan";
        return out;
    }
    if (value.is_inf()) {
        out += uppercase ? "INF" : "inf";
        return out;
    }

    const char* const hex_digits = uppercase ? "0123456789ABCDEF" : "0123456789abcdef";

    // The hexadecimal form is a direct reading of the encoding, so it needs none of the decimal
    // machinery and is exact by construction at any precision.
    if (style == 'a') {
        out += uppercase ? "0X" : "0x";
        if (value.is_zero()) {
            out += '0';
            if (alternate || precision > 0) {
                out += '.';
                for (int32_t i = 0; i < precision; ++i)
                    out += '0';
            }
            append_exponent(out, 0, uppercase ? 'P' : 'p');
            return out;
        }

        uint64_t low = 0, high = 0;
        int32_t exponent = 0;
        uint32_t sign = 0;
        value.get_components(low, high, exponent, sign);

        // 28 hex digits hold the whole fraction, and the leading one sits above the point.
        char fraction[28];
        for (int32_t i = 0; i < 28; ++i) {
            const int32_t shift = 108 - 4 * i;
            const uint64_t nibble = (shift >= 64) ? ((high >> (shift - 64)) & 0xF) : ((low >> shift) & 0xF);
            fraction[i] = hex_digits[nibble];
        }

        const auto hex_value = [](char c) {
            if (c >= '0' && c <= '9')
                return c - '0';
            return ((c >= 'a') ? (c - 'a') : (c - 'A')) + 10;
        };

        int32_t kept = (precision >= 0) ? precision : 28;
        if (kept > 28)
            kept = 28;
        // Rounding a shortened fraction: the first dropped digit decides, ties go up. A carry can
        // run through the fraction and into the leading one, which is the next power of two.
        bool leading_two = false;
        if (precision >= 0 && precision < 28) {
            if (hex_value(fraction[precision]) >= 8) {
                int32_t index = precision - 1;
                while (index >= 0) {
                    const int32_t digit = hex_value(fraction[index]);
                    if (digit != 15) {
                        fraction[index] = hex_digits[digit + 1];
                        break;
                    }
                    fraction[index] = '0';
                    --index;
                }
                if (index < 0)
                    leading_two = true;
            }
        }
        if (precision < 0) {
            // Trailing zeros carry no information in the shortest form.
            while (kept > 0 && fraction[kept - 1] == '0')
                --kept;
        }

        out += leading_two ? '2' : '1';
        if (kept > 0 || alternate) {
            out += '.';
            out.append(fraction, static_cast<size_t>(kept));
        }
        append_exponent(out, exponent, uppercase ? 'P' : 'p');
        return out;
    }

    if (value.is_zero()) {
        out += '0';
        // The fixed and scientific layouts always show the digits asked for; the general one shows
        // none unless the alternate form wanted them.
        int32_t zeros = 0;
        if (style == 'f' || style == 'e')
            zeros = (precision >= 0) ? precision : DEFAULT_PRECISION;
        else if (alternate)
            zeros = (precision > 0) ? precision - 1 : DEFAULT_PRECISION - 1;

        if (zeros > 0 || alternate) {
            out += '.';
            out.append(static_cast<size_t>(zeros), '0');
        }
        if (style == 'e')
            append_exponent(out, 0, uppercase ? 'E' : 'e');
        return out;
    }

    uint64_t low = 0, high = 0;
    int32_t exponent = 0;
    uint32_t sign = 0;
    value.get_components(low, high, exponent, sign);

    char digits[MAX_OUTPUT_DIGITS];
    int32_t significant = 0;
    int32_t exponent10 = 0;

    if (style == '\0') {
        significant = shortest_digit_count(low, high, exponent);
        exponent10 = to_decimal_digits(low, high, exponent, significant, digits);
    } else if (style == 'e') {
        significant = ((precision >= 0) ? precision : DEFAULT_PRECISION) + 1;
        if (significant > MAX_OUTPUT_DIGITS)
            significant = MAX_OUTPUT_DIGITS;
        exponent10 = to_decimal_digits(low, high, exponent, significant, digits);
    } else if (style == 'g') {
        significant = (precision > 0) ? precision : ((precision == 0) ? 1 : DEFAULT_PRECISION);
        if (significant > MAX_OUTPUT_DIGITS)
            significant = MAX_OUTPUT_DIGITS;
        exponent10 = to_decimal_digits(low, high, exponent, significant, digits);
    } else {
        // Fixed notation asks for a count of digits after the point rather than significant ones,
        // and how many that is depends on where the value sits. One digit is generated first to
        // find that out, then the real request is made.
        const int32_t after_point = (precision >= 0) ? precision : DEFAULT_PRECISION;
        bool exact = false;
        int32_t probe_exponent10 = to_decimal_digits(low, high, exponent, 1, digits, &exact);
        significant = probe_exponent10 + after_point;

        // A single digit probe rounds, and a value just below a power of ten rounds up into the
        // next one: 99 comes back as 1e2 and claims one more digit than it has. Generating the
        // digits reveals the true exponent, and one retry with it is always enough because the
        // rounding can only move the exponent by one.
        if (significant > 0 && significant <= MAX_OUTPUT_DIGITS) {
            const int32_t actual_exponent10 = to_decimal_digits(low, high, exponent, significant, digits);
            if (actual_exponent10 != probe_exponent10) {
                probe_exponent10 = actual_exponent10;
                significant = actual_exponent10 + after_point;
            }
        }

        if (significant <= 0) {
            // The value is below half of the last place asked for, or exactly on it. A tie goes to
            // the even digit, which is the zero already there.
            const bool round_up = (significant == 0) && (digits[0] > '5' || (digits[0] == '5' && !exact));
            out += '0';
            if (after_point > 0 || alternate) {
                out += '.';
                for (int32_t i = 0; i < after_point; ++i)
                    out += (round_up && i == after_point - 1) ? '1' : '0';
            }
            return out;
        }

        if (significant > MAX_OUTPUT_DIGITS) {
            // Past MAX_OUTPUT_DIGITS the exact expansion is no longer produced and the rest of
            // the request is filled with zeros. The expansion of a subnormal runs to 16494 places
            // after the point; the cap is set at what writing the widest finite value out in full
            // needs, which is 4933.
            exponent10 = to_decimal_digits(low, high, exponent, MAX_OUTPUT_DIGITS, digits);
            const int32_t padding = significant - MAX_OUTPUT_DIGITS;
            significant = MAX_OUTPUT_DIGITS;
            std::string body(digits, static_cast<size_t>(significant));
            body.append(static_cast<size_t>(padding), '0');
            // Lay the point out from the exponent the generator returned.
            if (exponent10 <= 0) {
                out += "0.";
                out.append(static_cast<size_t>(-exponent10), '0');
                out += body;
            } else {
                out.append(body, 0, static_cast<size_t>(exponent10));
                out += '.';
                out.append(body, static_cast<size_t>(exponent10), std::string::npos);
            }
            return out;
        }

        exponent10 = to_decimal_digits(low, high, exponent, significant, digits);
    }

    const int32_t scientific_exponent = exponent10 - 1;

    int32_t kept = significant;
    if ((style == 'g' || style == '\0') && !alternate) {
        while (kept > 1 && digits[kept - 1] == '0')
            --kept;
    }

    bool scientific = (style == 'e');
    if (style == 'g') {
        // printf's rule for %g, which std::format keeps for the g type.
        scientific = (scientific_exponent < -4) || (scientific_exponent >= significant);
    } else if (style == '\0') {
        // The default is whichever layout is shorter, with a tie going to the fixed one, which is
        // what std::to_chars produces and therefore what std::format gives a double. The rule for
        // %g would print 100 as 1e+02: it compares the exponent against the digit count, and the
        // shortest form of a round number has very few digits.
        int32_t fixed_length = 0;
        if (exponent10 <= 0)
            fixed_length = 2 - exponent10 + kept;
        else if (exponent10 >= kept)
            fixed_length = exponent10;
        else
            fixed_length = kept + 1;

        int32_t exponent_digits = 2;
        for (int32_t magnitude = (scientific_exponent < 0) ? -scientific_exponent : scientific_exponent; magnitude >= 100; magnitude /= 10)
            ++exponent_digits;
        const int32_t scientific_length = kept + ((kept > 1) ? 1 : 0) + 2 + exponent_digits;

        scientific = scientific_length < fixed_length;
    }

    if (scientific) {
        out += digits[0];
        if (kept > 1 || alternate) {
            out += '.';
            out.append(digits + 1, static_cast<size_t>(kept - 1));
        }
        append_exponent(out, scientific_exponent, uppercase ? 'E' : 'e');
        return out;
    }

    if (exponent10 <= 0) {
        out += "0.";
        out.append(static_cast<size_t>(-exponent10), '0');
        out.append(digits, static_cast<size_t>(kept));
    } else if (exponent10 >= kept) {
        out.append(digits, static_cast<size_t>(kept));
        out.append(static_cast<size_t>(exponent10 - kept), '0');
        if (alternate)
            out += '.';
    } else {
        out.append(digits, static_cast<size_t>(exponent10));
        out += '.';
        out.append(digits + exponent10, static_cast<size_t>(kept - exponent10));
    }

    // Fixed notation pads out to the requested number of digits after the point.
    if (style == 'f') {
        const int32_t after_point = (precision >= 0) ? precision : DEFAULT_PRECISION;
        const size_t point = out.find('.');
        const size_t written = (point == std::string::npos) ? 0 : (out.size() - point - 1);
        if (after_point > 0) {
            if (point == std::string::npos)
                out += '.';
            for (size_t i = written; i < static_cast<size_t>(after_point); ++i)
                out += '0';
        }
    }
    return out;
}
}  // namespace detail

/**
 * @brief The shortest decimal string that reads back as this exact value.
 * @param value Value to convert
 * @return Decimal text, in scientific notation when the exponent is far from zero.
 */
[[nodiscard]] inline std::string to_string(const float128& value)
{
    return detail::render(value, '\0', -1, false, false, '-');
}

/**
 * @brief Converts text to a float128, in the shape of std::from_chars.
 *
 * Reads the longest prefix of [first, last) that forms a number: an optional sign, decimal digits
 * with an optional point and an optional exponent, or one of inf, infinity and nan. The result is
 * the representable value nearest the one the text names, correctly rounded.
 *
 * @param first Start of the text
 * @param last One past the end of the text
 * @param value Receives the parsed value, untouched when nothing was parsed
 * @return ptr points past what was consumed; ec is invalid_argument when no number was found and
 *         result_out_of_range when the value is beyond the format's range.
 */
inline std::from_chars_result from_chars(const char* first, const char* last, float128& value)
{
    std::from_chars_result result {first, std::errc {}};
    const char* cursor = first;
    if (cursor == last) {
        result.ec = std::errc::invalid_argument;
        return result;
    }

    bool negative = false;
    if (*cursor == '-' || *cursor == '+') {
        negative = (*cursor == '-');
        ++cursor;
    }

    const auto lower = [](char c) { return static_cast<char>((c >= 'A' && c <= 'Z') ? (c - 'A' + 'a') : c); };
    const auto matches = [&](const char* word, size_t length) {
        if (static_cast<size_t>(last - cursor) < length)
            return false;
        for (size_t i = 0; i < length; ++i) {
            if (lower(cursor[i]) != word[i])
                return false;
        }
        return true;
    };

    if (matches("inf", 3)) {
        cursor += 3;
        if (matches("inity", 5))
            cursor += 5;
        value = negative ? -float128::inf() : float128::inf();
        result.ptr = cursor;
        return result;
    }
    if (matches("nan", 3)) {
        cursor += 3;
        value = float128::nan();
        value.set_sign(negative ? 1 : 0);
        result.ptr = cursor;
        return result;
    }

    // The digits are collected without the point, which only decides the exponent.
    char digits[detail::MAX_SIGNIFICANT_DIGITS];
    int32_t count = 0;
    int32_t exponent10 = 0;
    bool any_digit = false;
    bool seen_significant = false;

    // What is being built is `digits` read as an integer, times 10^exponent10. A digit before the
    // point that does not fit raises the exponent, since its place is still there. One after the
    // point lowers it when it is kept, and a leading zero after the point lowers it as well even
    // though nothing is stored for it.
    const auto take_integer = [&](char digit) {
        any_digit = true;
        seen_significant = seen_significant || (digit != '0');
        if (!seen_significant)
            return;
        if (count < detail::MAX_SIGNIFICANT_DIGITS)
            digits[count++] = digit;
        else
            ++exponent10;
    };
    const auto take_fraction = [&](char digit) {
        any_digit = true;
        seen_significant = seen_significant || (digit != '0');
        if (!seen_significant) {
            --exponent10;
            return;
        }
        if (count < detail::MAX_SIGNIFICANT_DIGITS) {
            digits[count++] = digit;
            --exponent10;
        }
    };

    while (cursor < last && *cursor >= '0' && *cursor <= '9')
        take_integer(*cursor++);

    if (cursor < last && *cursor == '.') {
        ++cursor;
        while (cursor < last && *cursor >= '0' && *cursor <= '9')
            take_fraction(*cursor++);
    }

    if (!any_digit) {
        result.ec = std::errc::invalid_argument;
        return result;
    }
    result.ptr = cursor;

    // An exponent is only consumed when it is well formed, so that "1e" reads as 1.
    if (cursor < last && (*cursor == 'e' || *cursor == 'E')) {
        const char* probe = cursor + 1;
        bool exponent_negative = false;
        if (probe < last && (*probe == '-' || *probe == '+')) {
            exponent_negative = (*probe == '-');
            ++probe;
        }
        if (probe < last && *probe >= '0' && *probe <= '9') {
            int64_t magnitude = 0;
            while (probe < last && *probe >= '0' && *probe <= '9') {
                if (magnitude < 1000000)
                    magnitude = magnitude * 10 + (*probe - '0');
                ++probe;
            }
            exponent10 += static_cast<int32_t>(exponent_negative ? -magnitude : magnitude);
            result.ptr = probe;
        }
    }

    uint64_t low = 0, high = 0;
    int32_t exponent = 0;
    if (!detail::from_decimal_digits(digits, count, exponent10, low, high, exponent)) {
        if (exponent > 100000) {
            value = negative ? -float128::inf() : float128::inf();
            result.ec = std::errc::result_out_of_range;
        } else {
            value = float128();
            value.set_sign(negative ? 1 : 0);
            if (exponent < -100000)
                result.ec = std::errc::result_out_of_range;
        }
        return result;
    }

    value.set_components(low, high, exponent, negative ? 1u : 0u);
    return result;
}

/**
 * @brief Converts a float128 to text, in the shape of std::to_chars.
 * @param first Start of the output range
 * @param last One past the end of the output range
 * @param value Value to convert
 * @param fmt Layout to use
 * @param precision Digits after the point, or significant digits for `general`
 * @return ptr points past the last character written; ec is value_too_large when the range was too
 *         small, in which case nothing was written.
 */
inline std::to_chars_result to_chars(char* first, char* last, const float128& value, std::chars_format fmt, int precision)
{
    char style = '\0';
    switch (fmt) {
    case std::chars_format::scientific: style = 'e'; break;
    case std::chars_format::fixed:      style = 'f'; break;
    case std::chars_format::hex:        style = 'a'; break;
    default:                            style = 'g'; break;
    }

    const std::string text = detail::render(value, style, precision, false, false, '-');
    if (text.size() > static_cast<size_t>(last - first))
        return {last, std::errc::value_too_large};

    for (size_t i = 0; i < text.size(); ++i)
        first[i] = text[i];
    return {first + text.size(), std::errc {}};
}

/// @brief Converts a float128 to text with the default precision for the layout.
inline std::to_chars_result to_chars(char* first, char* last, const float128& value, std::chars_format fmt)
{
    return to_chars(first, last, value, fmt, -1);
}

/// @brief Converts a float128 to the shortest text that reads back as the same value.
inline std::to_chars_result to_chars(char* first, char* last, const float128& value)
{
    const std::string text = detail::render(value, '\0', -1, false, false, '-');
    if (text.size() > static_cast<size_t>(last - first))
        return {last, std::errc::value_too_large};

    for (size_t i = 0; i < text.size(); ++i)
        first[i] = text[i];
    return {first + text.size(), std::errc {}};
}

/**
 * @brief Writes a float128 to a stream, honouring its formatting state.
 *
 * The stream's precision, its fixed, scientific and hexfloat flags, showpos, showpoint, uppercase,
 * width, fill and adjustfield are all applied, so a float128 behaves in a stream the way a double
 * does.
 */
inline std::ostream& operator<<(std::ostream& os, const float128& value)
{
    const std::ios_base::fmtflags flags = os.flags();
    const std::ios_base::fmtflags style_flags = flags & std::ios_base::floatfield;

    char style = '\0';
    int32_t precision = static_cast<int32_t>(os.precision());
    if (style_flags == std::ios_base::fixed) {
        style = 'f';
    } else if (style_flags == std::ios_base::scientific) {
        style = 'e';
    } else if (style_flags == (std::ios_base::fixed | std::ios_base::scientific)) {
        style = 'a';
        precision = -1;  // hexfloat ignores the precision
    } else {
        // The default floatfield is the general layout, where the precision counts significant
        // digits and zero means one.
        style = 'g';
    }

    const char sign_char = (flags & std::ios_base::showpos) ? '+' : '-';
    const bool uppercase = (flags & std::ios_base::uppercase) != 0;
    const bool alternate = (flags & std::ios_base::showpoint) != 0;
    std::string text = detail::render(value, style, precision, uppercase, alternate, sign_char);

    const std::streamsize width = os.width();
    os.width(0);
    if (static_cast<std::streamsize>(text.size()) < width) {
        const size_t padding = static_cast<size_t>(width) - text.size();
        const std::ios_base::fmtflags adjust = flags & std::ios_base::adjustfield;
        if (adjust == std::ios_base::left) {
            text.append(padding, os.fill());
        } else if (adjust == std::ios_base::internal && !text.empty() && (text[0] == '-' || text[0] == '+')) {
            // internal puts the fill between the sign and the digits
            text.insert(1, padding, os.fill());
        } else {
            text.insert(0, padding, os.fill());
        }
    }
    return os << text;
}

/**
 * @brief Reads a float128 from a stream.
 *
 * Accepts what the constructor from a string does. Failure leaves the value untouched and sets
 * failbit, the way the builtin extractors do.
 */
inline std::istream& operator>>(std::istream& is, float128& value)
{
    std::string token;
    if (!(is >> token))
        return is;

    float128 parsed;
    const std::from_chars_result result = from_chars(token.data(), token.data() + token.size(), parsed);
    if (result.ec == std::errc::invalid_argument || result.ptr != token.data() + token.size()) {
        is.setstate(std::ios_base::failbit);
        return is;
    }

    value = parsed;
    return is;
}

}  // namespace fp128

namespace std
{
/**
 * @brief Numeric properties of fp128::float128, the binary128 interchange format.
 *
 * Specialized so that a function template written against a builtin floating point type - anything
 * reaching for numeric_limits<T>::epsilon() to size a tolerance, or for max() to seed a minimum -
 * compiles and behaves correctly when instantiated with float128.
 *
 * The values are given as encodings rather than computed, which keeps every one of them usable in
 * a constant expression and independent of the string parser.
 */
template <> class numeric_limits<fp128::float128>
{
public:
    static constexpr bool is_specialized = true;
    static constexpr bool is_signed = true;
    static constexpr bool is_integer = false;
    static constexpr bool is_exact = false;
    static constexpr bool has_infinity = true;
    static constexpr bool has_quiet_NaN = true;
    static constexpr bool has_signaling_NaN = true;
    static constexpr bool is_bounded = true;
    static constexpr bool is_modulo = false;
    /// @brief The format is binary128 exactly as IEEE 754-2008 defines it.
    static constexpr bool is_iec559 = true;
    /// @brief No floating point status word exists, so nothing can trap or be flagged.
    static constexpr bool traps = false;
    static constexpr bool tinyness_before = false;
    static constexpr float_round_style round_style = round_to_nearest;

    /// @brief Mantissa bits including the implicit leading one.
    static constexpr int digits = 113;
    /// @brief Decimal digits that survive a round trip through the type.
    static constexpr int digits10 = 33;
    /// @brief Decimal digits needed to distinguish every value of the type.
    static constexpr int max_digits10 = 36;
    static constexpr int radix = 2;
    static constexpr int min_exponent = -16381;
    static constexpr int max_exponent = 16384;
    static constexpr int min_exponent10 = -4931;
    static constexpr int max_exponent10 = 4932;

    // has_denorm and has_denorm_loss are deprecated in C++23 but remain part of the interface a
    // generic caller may read, so they are provided.
    static constexpr float_denorm_style has_denorm = denorm_present;
    static constexpr bool has_denorm_loss = false;

    /// @brief Smallest positive normal value, 2^-16382.
    [[nodiscard]] static constexpr fp128::float128 min() noexcept { return fp128::float128(0, 0x0001000000000000ull); }
    /// @brief Largest finite value, (2 - 2^-112) * 2^16383.
    [[nodiscard]] static constexpr fp128::float128 max() noexcept { return fp128::float128(UINT64_MAX, 0x7FFEFFFFFFFFFFFFull); }
    /// @brief Most negative finite value.
    [[nodiscard]] static constexpr fp128::float128 lowest() noexcept { return fp128::float128(UINT64_MAX, 0xFFFEFFFFFFFFFFFFull); }
    /// @brief Difference between one and the next larger value, 2^-112.
    [[nodiscard]] static constexpr fp128::float128 epsilon() noexcept { return fp128::float128(0, 0x3F8F000000000000ull); }
    /// @brief Largest rounding error in units in the last place, one half.
    [[nodiscard]] static constexpr fp128::float128 round_error() noexcept { return fp128::float128(0, 0x3FFE000000000000ull); }
    /// @brief Smallest positive subnormal, 2^-16494.
    [[nodiscard]] static constexpr fp128::float128 denorm_min() noexcept { return fp128::float128(1, 0); }
    [[nodiscard]] static constexpr fp128::float128 infinity() noexcept { return fp128::float128::inf(); }
    [[nodiscard]] static constexpr fp128::float128 quiet_NaN() noexcept { return fp128::float128::nan(); }
    [[nodiscard]] static constexpr fp128::float128 signaling_NaN() noexcept { return fp128::float128::signaling_nan(); }
};

/// @brief const, volatile and cv qualified float128 have the same numeric properties.
template <> class numeric_limits<const fp128::float128> : public numeric_limits<fp128::float128>
{
};
template <> class numeric_limits<volatile fp128::float128> : public numeric_limits<fp128::float128>
{
};
template <> class numeric_limits<const volatile fp128::float128> : public numeric_limits<fp128::float128>
{
};

/**
 * @brief std::format support for fp128::float128.
 *
 * Accepts the whole floating point format specification: fill and alignment, a sign, the alternate
 * form, zero padding, a width, a precision, and any of the a, A, e, E, f, F, g and G types. An
 * empty specification gives the shortest text that reads back as the same value, which is what
 * std::format produces for a double.
 *
 * A width or precision given as a nested replacement field is not supported; both have to be
 * written out as digits.
 */
template <> struct formatter<fp128::float128, char>
{
    constexpr auto parse(basic_format_parse_context<char>& context)
    {
        auto it = context.begin();
        const auto end = context.end();
        if (it == end || *it == '}')
            return it;

        // [[fill]align]
        const auto is_align = [](char c) { return c == '<' || c == '>' || c == '^'; };
        if (it + 1 != end && is_align(*(it + 1))) {
            fill = *it;
            align = *(it + 1);
            it += 2;
        } else if (is_align(*it)) {
            align = *it++;
        }

        // [sign]
        if (it != end && (*it == '+' || *it == '-' || *it == ' '))
            sign = *it++;

        // [#]
        if (it != end && *it == '#') {
            alternate = true;
            ++it;
        }

        // [0], which is an alignment rather than a fill in its own right
        if (it != end && *it == '0') {
            zero_pad = true;
            ++it;
        }

        // [width]
        while (it != end && *it >= '0' && *it <= '9')
            width = width * 10 + (*it++ - '0');

        // [.precision]
        if (it != end && *it == '.') {
            ++it;
            precision = 0;
            while (it != end && *it >= '0' && *it <= '9')
                precision = precision * 10 + (*it++ - '0');
        }

        // [type]
        if (it != end && *it != '}') {
            switch (*it) {
            case 'a':
            case 'A':
            case 'e':
            case 'E':
            case 'f':
            case 'F':
            case 'g':
            case 'G':
                type = *it++;
                break;
            default:
                throw format_error("invalid type in the format specification for float128");
            }
        }

        if (it != end && *it != '}')
            throw format_error("unmatched brace in the format specification for float128");
        return it;
    }

    template <typename FormatContext> auto format(const fp128::float128& value, FormatContext& context) const
    {
        const string text = render(value);
        return std::copy(text.begin(), text.end(), context.out());
    }

    /// @brief Builds the text, including the padding the specification asks for.
    [[nodiscard]] string render(const fp128::float128& value) const
    {
        const bool uppercase = (type >= 'A' && type <= 'Z');
        const char style = static_cast<char>(uppercase ? (type - 'A' + 'a') : type);
        string text = fp128::detail::render(value, style, precision, uppercase, alternate, sign);

        if (static_cast<int>(text.size()) >= width)
            return text;

        const size_t padding = static_cast<size_t>(width) - text.size();
        // Zero padding goes between the sign and the digits, and only when no alignment was given.
        // A special value is never zero padded, the way it is not for a double either.
        if (zero_pad && align == 0 && !value.is_nan() && !value.is_inf()) {
            const size_t offset = (!text.empty() && (text[0] == '-' || text[0] == '+' || text[0] == ' ')) ? 1u : 0u;
            text.insert(offset, padding, '0');
            return text;
        }

        switch (align) {
        case '<':
            text.append(padding, fill);
            break;
        case '^':
            text.insert(0, padding / 2, fill);
            text.append(padding - padding / 2, fill);
            break;
        case '>':
        default:
            // A number aligns right by default, as every arithmetic type does.
            text.insert(0, padding, fill);
            break;
        }
        return text;
    }

    char fill = ' ';       ///< Character the padding is made of.
    char align = 0;        ///< One of < > ^, or zero when none was given.
    char sign = '-';       ///< One of + - space.
    char type = 0;         ///< One of a A e E f F g G, or zero for the shortest form.
    bool alternate = false;///< The # flag: keep the point and the trailing zeros.
    bool zero_pad = false; ///< The 0 flag: pad with zeros after the sign.
    int width = 0;         ///< Minimum field width.
    int precision = -1;    ///< Digits after the point, or -1 when none was given.
};

/**
 * @brief Hash support, so a float128 can be a key in an unordered container.
 *
 * The two zeros compare equal and therefore have to hash equal, which the raw encoding would not
 * do: they differ in the sign bit.
 */
template <> struct hash<fp128::float128>
{
    [[nodiscard]] size_t operator()(const fp128::float128& value) const noexcept
    {
        uint64_t low = 0, high = 0;
        value.get_bits(low, high);
        if (value.is_zero())
            low = high = 0;

        // splitmix64's finalizer, applied to the two halves in turn
        uint64_t state = low + 0x9E3779B97F4A7C15ull;
        state = (state ^ (state >> 30)) * 0xBF58476D1CE4E5B9ull;
        state = (state ^ (state >> 27)) * 0x94D049BB133111EBull;
        state ^= high + 0x9E3779B97F4A7C15ull + (state << 6) + (state >> 2);
        state = (state ^ (state >> 30)) * 0xBF58476D1CE4E5B9ull;
        return static_cast<size_t>(state ^ (state >> 31));
    }
};

/// @name Common type
///
/// float128 is wider than every builtin arithmetic type, so a mixed expression produces a
/// float128. Without these the default common_type would try to form the ternary operator over the
/// two, which is ambiguous: float128 converts to double and double converts to float128.
/// @{
template <typename T>
    requires is_arithmetic_v<T>
struct common_type<fp128::float128, T>
{
    using type = fp128::float128;
};
template <typename T>
    requires is_arithmetic_v<T>
struct common_type<T, fp128::float128>
{
    using type = fp128::float128;
};
template <> struct common_type<fp128::float128, fp128::float128>
{
    using type = fp128::float128;
};
/// @}
}  // namespace std

#endif  // FP128_FLOAT128_H
