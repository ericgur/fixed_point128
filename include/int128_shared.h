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

************************************************************************************/

/**
 * @file int128_shared.h
 * @brief Implementation shared by the 128 bit signed and unsigned integer types.
 *
 * Provides @ref fp128::int128_base, a single class template parameterized on the
 * signedness of the type. @ref fp128::uint128_t and @ref fp128::int128_t are
 * aliases of its two instantiations, declared by uint128_t.h and int128_t.h
 * respectively. Include one of those headers rather than this one.
 *
 * @see uint128_t.h for the unsigned alias and its user-defined literal.
 * @see int128_t.h for the signed alias, its user-defined literal and abs().
 * @see fp128_shared.h for supporting intrinsics and utilities.
 */

#ifndef FP128_INT128_SHARED_H
#define FP128_INT128_SHARED_H

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
class fp128_gtest;  // Google test class

/***********************************************************************************
 *                                  Main Code
 ************************************************************************************/

/**
 * @brief 128 bit integer class, signed or unsigned depending on the template parameter.
 *
 * This class implements the standard operators a builtin integer type has.<BR>
 * All of int128_base's methods are inline for maximum performance.<BR>
 * Both instantiations store the same two QWORDs and only differ where the interpretation of
 * the msb matters, which is why a single template covers them. The signed instantiation uses
 * a two's complement representation, so addition, subtraction, multiplication, the bitwise
 * operators, the left shift and equality are bit identical between the two and are written
 * once, without a test on the signedness.
 *
 * <B>Implementation notes:</B>
 * <UL>
 * <LI>The signedness dependent code is selected with <TT>if constexpr (is_signed)</TT>, or by
 *     calling is_negative()/is_positive(), which are compile time constants for the unsigned
 *     instantiation and let the surrounding tests fold away.</LI>
 * <LI>Overflow is handled silently, similar to builtin integer operations.</LI>
 * <LI>An int128_base object is not thread safe. Accessing a const object from multiple threads is safe.</LI>
 * <LI>int128_base is <B>conditionally safe</B>, 2 different non const objects can be accessed concurrently.</LI>
 * <LI>Only 64 bit builds are supported.</LI>
 * </UL>
 *
 * <B>Compile time evaluation:</B><BR>
 * Everything that stays within integer arithmetic is constexpr:
 * <UL>
 * <LI>Construction from any builtin arithmetic type and from the two QWORDs, copy, move,
 *     assignment and the conversions to the integer and floating point types.</LI>
 * <LI>Addition, subtraction, multiplication, square(), the increment and decrement operators,
 *     the shifts, the bitwise operators, the unary operators and the comparisons.</LI>
 * <LI>The queries is_positive(), is_negative(), is_zero(), get_bit() and get_components(),
 *     and the constant one().</LI>
 * <LI>The math functions abs, sqr, lzcnt128, log2 and pow.</LI>
 * </UL>
 *
 * The bit counting and extended arithmetic intrinsics these rest on are not constant expressions,
 * so fp128_shared.h wraps each one in a constexpr function that serves a constant
 * evaluated call from a portable implementation of the same operation. A runtime call still
 * reaches the bare intrinsic and generates the same code it did before.
 *
 * The rest cannot be constexpr, for one of two reasons:
 * <UL>
 * <LI>Division and modulo, and sqrt which is built on them. A 128 bit divisor goes through
 *     div_32bit, which needs alloca and a goto, neither of which C++20 permits in a constexpr
 *     function. The paths for a divisor that fits in 64 bit are no better off: they rest on the
 *     _udiv128 intrinsic, and one of them writes both QWORDs through a pointer to the first,
 *     which is out of bounds as far as constant evaluation is concerned.</LI>
 * <LI>The string conversions allocate, and log/log10 look their result up in a function local
 *     static table, which a constexpr function may not declare.</LI>
 * </UL>
 *
 * @tparam IsSigned True for a two's complement signed type, false for an unsigned one.
 */
template <bool IsSigned> class FP128_ALIGN16 int128_base
{
    // build time validation of template parameters
    static_assert(sizeof(void*) == 8, "int128_base is supported in 64 bit builds only!");
    friend class fp128_gtest;

private:
    //
    // members
    //
    uint64_t low;   ///< Lower 64 bits (QWORD).
    uint64_t high;  ///< Upper 64 bits (QWORD). For the signed type the sign is in the MSB (two's complement).

    /**
     * @brief Produces the upper QWORD when sign extending a builtin integral value.
     * The test is compiled out for unsigned types, which also keeps the compiler from warning
     * about an unsigned or bool value being compared against zero.
     * @tparam T Any builtin integral type
     * @param x Input value
     * @return All ones when x is negative, zero otherwise.
     */
    template <typename T> [[nodiscard]] static constexpr uint64_t HighFromIntegral(T x) noexcept
    {
        if constexpr (std::is_signed_v<T>) {
            return (x < T(0)) ? UINT64_MAX : 0;
        } else {
            return 0;
        }
    }
    /**
     * @brief Compares the magnitude of this object with another, treating both as unsigned.
     *
     * The division code converts its operands to absolute values and from that point on they are
     * magnitudes, not signed numbers. The signed comparison operators would misread the magnitude
     * of the most negative value, whose absolute value keeps the sign bit set.
     *
     * @param rhs Right hand side operand
     * @return True when this object's bit pattern is larger when read as an unsigned integer.
     */
    [[nodiscard]] FP128_FORCE_INLINE constexpr bool MagnitudeGreaterThan(const int128_base& rhs) const noexcept
    {
        return high > rhs.high || (high == rhs.high && low > rhs.low);
    }
    /**
     * @brief Extracts the mantissa of a floating point representation of this object, rounded.
     *
     * Keeps the frac_bits+1 most significant bits and rounds the dropped remainder to nearest,
     * ties going to even. This is what IEEE 754 mandates and what the builtin integer to floating
     * point conversions do, so the result matches converting via a wider builtin type.<BR>
     * Rounding needs three pieces of information: the bits that survive, the highest dropped bit
     * (which chooses the direction) and whether anything below it is set (the sticky bit, which
     * distinguishes an exact tie from a remainder larger than half).<BR>
     * Expects this object to hold a magnitude (the conversion operators take the absolute value
     * before calling it), so the sign bit plays no part.
     *
     * @param expo Bit position of the msb of this object. Must be larger than frac_bits.
     * @param frac_bits Fraction bit count of the target type, 23 for float or 52 for double.
     * @return The rounded mantissa. Holds frac_bits+1 bits, or frac_bits+2 when rounding carried
     *         into the next power of 2.
     */
    [[nodiscard]] FP128_INLINE constexpr uint64_t RoundedMantissa(uint64_t expo, int32_t frac_bits) const noexcept
    {
        const int32_t shift = static_cast<int32_t>(expo) - frac_bits;  // count of dropped bits, at least 1
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
    static constexpr bool is_signed = IsSigned;  ///< True when the msb of the high QWORD is a sign bit.

    typedef int128_base type;      ///< Self type alias.
    typedef int128_base* ptr_type; ///< Pointer type alias.
    typedef int128_base& ref_type; ///< Reference type alias.

    //
    // Constructors
    //

    /**
     * @brief Default constructor, creates an instance with a value of zero.
     * Explicitly zeroing both QWORDs rather than defaulting the constructor keeps the type usable
     * in a constant expression, which a default initialized object with indeterminate members is not.
     */
    FP128_FORCE_INLINE constexpr int128_base() noexcept : low(0), high(0) {}
    /**
     * @brief Copy constructor
     * Defaulted, which is what makes the type trivially copyable.
     */
    constexpr int128_base(const int128_base&) noexcept = default;
    /**
     * @brief Move constructor
     * Doesn't modify the right hand side object. Acts like a copy constructor.
     */
    constexpr int128_base(int128_base&&) noexcept = default;
    /**
     * @brief Constructor from the double type
     * Underflow goes to zero. Overflow, NaN and +-INF go to the largest supported magnitude of
     * the matching sign.<BR>
     * For the unsigned type a negative value wraps around via 2's complement, matching both the
     * integer constructors and the way builtin unsigned types behave. For example -5.0 yields
     * 2^128 - 5.
     * @param x Input value
     */
    constexpr int128_base(double x) noexcept
    {
        const Double d(x);
        // very common case
        if (x == 0) {
            low = high = 0;
            return;
        }

        const int32_t e = static_cast<int32_t>(d.e()) - 1023;
        uint64_t f = d.f();

        // overflow which catches NaN and Inf
        if constexpr (is_signed) {
            // -(2**127) is represented correctly (doesn't overflow)
            if (e > 126) {
                if (d.s() == 0) {
                    high = 0x8000000000000000ull - 1;
                    low = UINT64_MAX;
                } else {
                    high = 0x8000000000000000ull;
                    low = 0;
                }
                return;
            }
        } else {
            if (e > 127) {
                high = low = UINT64_MAX;
                return;
            }
        }

        // normal number, produces non zero value
        if (e >= 0) {
            // bit 52 in f is the unity value of the float. it needs to move to the unity position in fixed point
            f |= FP128_ONE_SHIFT(dbl_frac_bits);
            const int32_t bits_to_shift = e - dbl_frac_bits;
            low = f;
            high = 0;

            // f fits in high QWORD
            if (bits_to_shift > 0) {
                *this <<= bits_to_shift;
            }
            // shift right
            else {
                *this >>= -bits_to_shift;
            }

            // negative number, the unsigned type wraps around just like the integer constructors do
            if (d.s()) {
                twos_complement128(low, high);
            }
        }
        // too small to be represented, no need to bother.
        else {
            high = low = 0;
        }
    }
    /**
     * @brief Constructor from any builtin integral type.
     *
     * A single constrained template is used instead of one overload per fixed width type.
     * Overloads would only cover the types the fixed width aliases happen to name on a given
     * platform, leaving the rest ambiguous: on MSVC, long and unsigned long are distinct from
     * both int and long long, so int128_base(1ul) could not pick a best match. Which of long and
     * long long int64_t names also differs between LP64 and LLP64, so no fixed set of overloads
     * covers every platform.<BR>
     * Signed values are sign extended, unsigned values are zero extended. For the unsigned type a
     * negative value therefore wraps around via 2's complement, exactly like the builtin unsigned
     * types do.
     * @tparam T Any builtin integral type (including bool and the character types)
     * @param x Input value
     */
    template <typename T>
        requires std::is_integral_v<T>
    FP128_FORCE_INLINE constexpr int128_base(T x) noexcept : low(static_cast<uint64_t>(x)), high(HighFromIntegral(x))
    {
    }
    /**
     * @brief Constructor from const char* (C string).
     * Allows creating 128 bit values from a string. Much slower than the other constructors.<BR>
     * Input string can be decimal or hex. Hex initialization is faster. Commas and apostrophes are
     * ignored, so both conventional digit grouping and the C++ digit separator are accepted, the
     * latter being what reaches the type's user defined literals. A nullptr, an empty string or a
     * leading illegal character all
     * produce zero; parsing otherwise stops at the first character that isn't a digit of the
     * detected base. Overflow wraps around silently, like the builtin integer types.<BR>
     * An optional leading '+' or '-' is accepted. For the unsigned type a negative value wraps
     * around via 2's complement, so "-5" and -5.0 both produce 2^128 - 5, matching strtoull.<BR>
     * Never throws: an allocation failure produces zero.
     * @param x Input string
     */
    int128_base(const char* x) noexcept
    {
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

        // trim leading white space. The cast keeps a negative char from reaching isspace, which
        // only accepts values representable as unsigned char (or EOF).
        while (*p && isspace(static_cast<unsigned char>(*p)))
            ++p;

        // An optional sign, which has to be consumed before the base prefix so that "-0x10" is
        // recognized as hex. The magnitude is parsed first and negated at the end.
        bool negative = false;
        if (*p == '-' || *p == '+') {
            negative = (*p == '-');
            ++p;
        }

        // base 10 or 16 are supported
        const uint32_t base = (0 == strncmp("0x", p, 2)) ? 16u : 10u;
        if (base == 16)
            p += 2;

        // trim leading zeros
        while (*p == '0')
            ++p;

        // convert one digit at a time
        while (*p && !isspace(static_cast<unsigned char>(*p)) && *p != '.') {
            uint64_t d = static_cast<unsigned char>(*p);
            if (d >= '0' && d <= '9')
                d -= '0';
            else if (base == 16 && (d >= 'a' && d <= 'f'))
                d = 10ull + d - 'a';
            // Digit grouping separators are skipped. The apostrophe matters for the user defined
            // literals: a raw literal operator receives the source characters of its token and the
            // compiler does not strip the C++ digit separator out of them, so 1'000'000_uint128
            // arrives here as "1'000'000" and would otherwise stop parsing at the first separator,
            // silently yielding 1.
            else if (d == ',' || d == '\'') {
                ++p;
                continue;
            } else
                break;  // treat an unknown char as end of string

            // 4 bits per digit
            if (base == 16) {
                *this <<= 4;
                low |= d;
            } else {
                *this *= base;
                *this += d;
            }

            ++p;
        }

        // negate the magnitude. For the unsigned type this wraps around, the same way the numeric
        // constructors and strtoull do.
        if (negative) {
            twos_complement128(low, high);
        }
    }
    /**
     * @brief Constructor from std::string.
     * Allows creating very high precision values. Much slower than the other constructors.
     * @param x Input string
     */
    int128_base(const std::string& x) noexcept { *this = x.c_str(); }
    /**
     * @brief Constructor from the 2 int128_base elements, useful for creating special constants.
     * @param l Low QWORD
     * @param h High QWORD
     */
    FP128_FORCE_INLINE constexpr int128_base(uint64_t l, uint64_t h) noexcept : low(l), high(h) {}
    /**
     * @brief Destructor
     */
    constexpr ~int128_base() noexcept = default;
    /**
     * @brief Assignment operator
     * Defaulted, which is what makes the type trivially copyable.
     */
    FP128_INLINE constexpr int128_base& operator=(const int128_base&) noexcept = default;
    /**
     * @brief Move assignment operator
     * Defaulted, which is what makes the type trivially copyable.
     */
    FP128_INLINE constexpr int128_base& operator=(int128_base&&) noexcept = default;
    /**
     * @brief Assignment operator
     * @param rhs Value to copy from
     * @return This object.
     */
    template <typename T> FP128_FORCE_INLINE constexpr int128_base& operator=(const T& rhs) noexcept
    {
        *this = int128_base(rhs);
        return *this;
    }

    //
    // conversion operators
    //
    /**
     * @brief operator uint64_t - converts to a uint64_t
     * @return Object value.
     */
    [[nodiscard]] FP128_FORCE_INLINE constexpr operator uint64_t() const noexcept { return low; }
    /**
     * @brief operator int64_t - converts to a int64_t
     * @return Object value.
     */
    [[nodiscard]] FP128_FORCE_INLINE constexpr operator int64_t() const noexcept { return static_cast<int64_t>(low); }
    /**
     * @brief operator uint32_t - converts to a uint32_t
     * @return Object value.
     */
    [[nodiscard]] FP128_FORCE_INLINE constexpr operator uint32_t() const noexcept { return static_cast<uint32_t>(low); }
    /**
     * @brief operator int32_t - converts to a int32_t
     * @return Object value.
     */
    [[nodiscard]] FP128_FORCE_INLINE constexpr operator int32_t() const noexcept { return static_cast<int32_t>(low); }
    /**
     * @brief operator float - converts to a float
     * The signed range fits within a float's exponent, so no signed value overflows. Unsigned
     * values above FLT_MAX become infinity, as they do with the builtin unsigned conversions:
     * the carry out of the rounding step below pushes the exponent to its all ones encoding.
     * @return Object value, rounded to nearest with ties going to even.
     */
    [[nodiscard]] FP128_INLINE constexpr operator float() const noexcept
    {
        if (is_zero())
            return 0;

        // the magnitude is what gets converted, the sign is carried in the result's sign bit.
        // Both tests are a compile time false for the unsigned instantiation.
        int128_base temp = *this;
        const bool sign = is_negative();
        if (sign) {
            twos_complement128(temp.low, temp.high);
        }

        // log2 of the magnitude. The free function form is used because the absolute value of the
        // most negative value keeps its sign bit set, which the int128_base overload rejects.
        const uint64_t expo = log2(temp.low, temp.high);  // returns the bit location of the msb
        // the value fits in the fraction, no bits are lost so no rounding is needed.
        // float doesn't hold the msb, it's implicit, so bit [22:0] hold the rest of the value
        if (expo <= flt_frac_bits) {
            return Float::make(sign, static_cast<uint32_t>(expo) + 127, static_cast<uint32_t>(temp.low << (flt_frac_bits - expo)));
        }

        // more bits than the fraction can hold, drop the extra ones with rounding.
        // rounding up can carry into the next power of 2, costing a fraction bit and
        // incrementing the exponent. An exponent of 255 is the encoding of infinity, which is
        // exactly what an unsigned value above FLT_MAX must produce.
        const uint64_t mant = temp.RoundedMantissa(expo, flt_frac_bits);
        const uint64_t carry = mant >> (flt_frac_bits + 1);
        return Float::make(sign, static_cast<uint32_t>(expo + carry) + 127, static_cast<uint32_t>(mant >> carry));
    }
    /**
     * @brief operator double - converts to a double
     * The whole 128 bit range is within reach of a double's exponent, so no value overflows.
     * @return Object value, rounded to nearest with ties going to even.
     */
    [[nodiscard]] FP128_INLINE constexpr operator double() const noexcept
    {
        if (is_zero())
            return 0;

        // the magnitude is what gets converted, the sign is carried in the result's sign bit.
        // Both tests are a compile time false for the unsigned instantiation.
        int128_base temp = *this;
        const bool sign = is_negative();
        if (sign) {
            twos_complement128(temp.low, temp.high);
        }

        // log2 of the magnitude. The free function form is used because the absolute value of the
        // most negative value keeps its sign bit set, which the int128_base overload rejects.
        // It is also noexcept, which the overload that throws on a zero input is not.
        const uint64_t expo = log2(temp.low, temp.high);  // returns the bit location of the msb
        // the value fits in the fraction, no bits are lost so no rounding is needed.
        // double doesn't hold the msb, it's implicit, so bit [51:0] hold the rest of the value
        if (expo <= dbl_frac_bits) {
            return Double::make(sign, expo + 1023, temp.low << (dbl_frac_bits - expo));
        }

        // more bits than the fraction can hold, drop the extra ones with rounding.
        // rounding up can carry into the next power of 2, costing a fraction bit and
        // incrementing the exponent
        const uint64_t mant = temp.RoundedMantissa(expo, dbl_frac_bits);
        const uint64_t carry = mant >> (dbl_frac_bits + 1);
        return Double::make(sign, expo + carry + 1023, mant >> carry);
    }
    /**
     * @brief operator long double - converts to a long double
     * @return Object value.
     */
    [[nodiscard]] FP128_FORCE_INLINE constexpr operator long double() const noexcept { return operator double(); }
    /**
     * @brief Converts to a char* (slow) string holds all meaningful fraction bits.
     * @return object string representation
     */
    [[nodiscard]] FP128_FORCE_INLINE char* to_string() const noexcept { return operator char*(); }
    /**
     * @brief Converts to a std::string (slow) string holds all meaningful fraction bits.
     * @return object string representation
     */
    [[nodiscard]] FP128_FORCE_INLINE operator std::string() const { return operator char*(); }
    /**
     * @brief Converts to a C string (slow) string holds all meaningful fraction bits.
     * The returned string is a statically, thread-allocated buffer.
     * Additional calls to this function from the same thread, overwrite the previous result.
     * @return C string with describing the value of the object.
     * @return object string representation
     */
    [[nodiscard]] explicit FP128_INLINE operator char*() const noexcept
    {
        constexpr int buff_size = 45;
        static thread_local char str[buff_size];

        // Small numbers, use the fast snprintf method.
        if constexpr (is_signed) {
            // The value round trips through an int64_t only when the high QWORD is nothing but
            // copies of the sign bit of the low QWORD. Testing high == 0 alone would hand %lld a
            // negative number for every value in [2^63, 2^64), printing 2^64-1 as "-1". It also
            // lets small negatives take this path, whose high QWORD is all ones.
            if (high == static_cast<uint64_t>(static_cast<int64_t>(low) >> 63)) {
                snprintf(str, buff_size, "%lld", static_cast<int64_t>(low));
                return str;
            }
        } else {
            if (high == 0) {
                snprintf(str, buff_size, "%llu", low);
                return str;
            }
        }

        // the digits are produced from the magnitude, both tests below are a compile time false
        // for the unsigned instantiation
        int128_base temp = *this;
        const bool sign = is_negative();
        if (sign) {
            twos_complement128(temp.low, temp.high);
        }

        // The digits are written backwards starting at the terminator, which sits at the very end
        // of the buffer. The longest result is 39 digits plus a sign, well within the 44 available.
        str[buff_size - 1] = 0;
        char* p = str + buff_size - 1;  // writing the string in reverse
        uint64_t q[2] {};
        uint64_t digit {};
        // as long as the intermediate value is >64 bit, use the more expensive long division.
        while (temp.high) {
            const uint64_t nom[2] = {temp.low, temp.high};
            div_64bit((uint64_t*)q, &digit, (uint64_t*)nom, 10ull, 2);
            temp.low = q[0];
            temp.high = q[1];
            *--p = static_cast<char>(digit + '0');
        }
        // the intermediate value fits in 64 bit, use the much faster native-64 bit division
        while (temp.low) {
            uint64_t r = temp.low % 10ull;
            temp.low /= 10ull;
            *--p = static_cast<char>(r + '0');
        }

        if (sign) {
            *--p = '-';
        }

        return p;
    }

    //
    // math operators
    //

    /**
     * @brief Performs right shift operation.
     * @param shift bits to shift
     * @return Temporary object with the result of the operation
     */
    [[nodiscard]] FP128_FORCE_INLINE constexpr int128_base operator>>(int32_t shift) const noexcept
    {
        int128_base temp(*this);
        return temp >>= shift;
    }
    /**
     * @brief Performs left shift operation.
     * @param shift bits to shift
     * @return Temporary object with the result of the operation
     */
    [[nodiscard]] FP128_FORCE_INLINE constexpr int128_base operator<<(int32_t shift) const noexcept
    {
        int128_base temp(*this);
        return temp <<= shift;
    }
    /**
     * @brief Add a value to this object
     * like other builtin integer types, overflow wraps around silently
     * @param rhs Right hand side operand
     * @return This object.
     */
    FP128_FORCE_INLINE constexpr int128_base& operator+=(const int128_base& rhs) noexcept
    {
        const uint8_t carry = addcarryx_u64(0, low, rhs.low, &low);
        addcarryx_u64(carry, high, rhs.high, &high);
        return *this;
    }
    /**
     * @brief Add a value to this object
     * like other builtin integer types, overflow wraps around silently
     * @param rhs Right hand side operand
     * @return This object.
     */
    template <typename T> FP128_FORCE_INLINE constexpr int128_base& operator+=(const T& rhs) noexcept { return operator+=(int128_base(rhs)); }
    /**
     * @brief Subtract a value from this object
     * like other builtin integer types, underflow wraps around silently
     * @param rhs Right hand side operand
     * @return This object.
     */
    FP128_FORCE_INLINE constexpr int128_base& operator-=(const int128_base& rhs) noexcept
    {
        // Subtracting directly rather than adding the two's complement of rhs. Both are the same
        // difference modulo 2^128, but neither compiler folds the negate-and-add form back into a
        // borrow chain: MSVC emits the whole neg/not/sete sequence and lands on twice the
        // instruction count of operator+=. Writing the difference back in place is safe when rhs
        // aliases this object, each QWORD of rhs being read before its counterpart is overwritten.
        const uint8_t borrow = subborrow_u64(0, low, rhs.low, &low);
        subborrow_u64(borrow, high, rhs.high, &high);
        return *this;
    }
    /**
     * @brief Subtract a value from this object
     * like other builtin integer types, underflow wraps around silently
     * @param rhs Right hand side operand
     * @return This object.
     */
    template <typename T> FP128_FORCE_INLINE constexpr int128_base& operator-=(const T& rhs) noexcept { return operator-=(int128_base(rhs)); }
    /**
     * @brief Multiplies a value to this object
     *
     * The same code serves both instantiations. A truncated 128 bit product is a multiplication
     * modulo 2^128, and in that ring the two's complement bit pattern of a negative value is its
     * value, so multiplying the raw bits gives the signed result without converting the operands
     * to magnitudes and reapplying the sign. int128_base::square() relies on the same identity.
     *
     * @param rhs Right hand side operand
     * @return This object.
     */
    FP128_FORCE_INLINE constexpr int128_base& operator*=(const int128_base& rhs) noexcept
    {
        // Snapshot both operands before writing anything. The multiply below overwrites both
        // QWORDs of this object, and rhs is allowed to alias it (e.g. a *= a), so reading
        // rhs after that point would pick up the partial result instead of the operand.
        const uint64_t l = low, h = high;
        const uint64_t rhs_l = rhs.low, rhs_h = rhs.high;

        // multiply low QWORDs
        low = mulx_u64(l, rhs_l, &high);

        // multiply low this and high rhs; multiply high this and low rhs
        high += l * rhs_h + h * rhs_l;
        return *this;
    }
    /**
     * @brief Squares this object in place.
     *
     * Cheaper than `operator*=(*this)`: a square is symmetric, so the `low*high` and `high*low`
     * cross products are the same value. Computing it once and doubling it replaces one of
     * the three multiplies with an addition. The temporary copy that `operator*=` needs to
     * guard against aliasing is not needed either.
     *
     * The result is truncated to 128 bit and is bit identical to `(*this) * (*this)`, which for
     * the signed type means it is negative in the same overflow cases.
     *
     * @return This object.
     */
    FP128_FORCE_INLINE constexpr int128_base& square() noexcept
    {
        const uint64_t l = low, h = high;

        // multiply the low QWORD by itself
        low = mulx_u64(l, l, &high);

        // the low * high cross product appears twice in the sum, so double it.
        // Both operands only contribute their lower 64 bit to the result.
        high += 2 * (l * h);
        return *this;
    }
    /**
     * @brief Multiplies a value to this object
     *
     * A non negative operand that fits in 64 bit skips one of the three partial products. It is
     * valid for both instantiations for the same reason operator*=(const int128_base&) is: the
     * product is taken modulo 2^128, where the two's complement bits of this object are its value.
     *
     * @param x Right hand side operand
     * @return This object.
     */
    template <typename T> FP128_FORCE_INLINE constexpr int128_base& operator*=(T x) noexcept
    {
        // floating point values are converted first. Casting them straight to uint64_t below
        // would be undefined behavior for anything at or above 2^64 and would silently drop the
        // upper bits of operands the int128_base constructor handles exactly.
        if constexpr (std::is_floating_point_v<T>) {
            return operator*=(int128_base(x));
        }
        // check if the type is signed or not
        // for negative values, convert to int128_base and multiply.
        else if constexpr (std::is_signed_v<T>) {
            if (x < 0)
                return operator*=(int128_base(x));
        }

        uint64_t temp;
        const uint64_t uval = static_cast<uint64_t>(x);
        // multiply low QWORDs
        low = mulx_u64(low, uval, &temp);
        high = high * uval + temp;

        return *this;
    }
    /**
     * @brief Divide this object by rhs.
     *
     * Rounds towards zero, matching the builtin integer types: -7 / 2 is -3, not -4. The two
     * operands are reduced to magnitudes, divided as unsigned values (which truncates), and the
     * sign is applied afterwards, which is exactly truncation towards zero. Every sign related
     * step below folds away for the unsigned instantiation.
     *
     * @param rhs_in Right hand side operator (denominator)
     * @return this object.
     */
    FP128_INLINE int128_base& operator/=(const int128_base& rhs_in)
    {
        // The signed path negates a negative divisor before dividing, so it needs its own copy.
        // The copy also makes a divisor that aliases this object (a /= a) safe.
        int128_base rhs = rhs_in;

        // check some trivial cases
        if (rhs.is_zero()) {
            FP128_INT_DIVIDE_BY_ZERO_EXCEPTION;
        }
        if (is_zero()) {
            low = high = 0;
            return *this;
        }
        if (rhs == *this) {
            low = 1;
            high = 0;
            return *this;
        }

        // convert both to absolute values. From here on the two are magnitudes, so every
        // comparison and shift below has to treat them as unsigned.
        const bool sign = is_negative();
        const bool rhs_sign = rhs.is_negative();
        const bool res_is_negative = sign != rhs_sign;
        if (sign) {
            twos_complement128(low, high);
        }
        if (rhs_sign) {
            twos_complement128(rhs.low, rhs.high);
        }

        // a smaller magnitude divided by a larger one truncates to zero, whatever the signs
        if (rhs.MagnitudeGreaterThan(*this)) {
            low = high = 0;
            return *this;
        }

        // exponent of 2, convert to a much faster shift operation
        if (1 == popcnt128(rhs.low, rhs.high)) {
            // the free function form of log2 is used because the magnitude of the most negative
            // value keeps its sign bit set, which the int128_base overload rejects
            const auto bits = static_cast<int32_t>(log2(rhs.low, rhs.high));
            // operator>>= would replicate the sign bit, this is a magnitude and needs a logical
            // shift. A shift of zero (dividing by 1) would also be undefined below.
            if (bits > 0) {
                low = shift_right128(low, high, bits);
                high = (bits < 64) ? (high >> bits) : 0;
            }
        } else {
            uint64_t q[2] {};
            const uint64_t nom[2] = {low, high};

            // optimization for when dividing by a small (<= 64 bit) integer
            if (rhs.high == 0) {
                if (div_64bit((uint64_t*)q, nullptr, (uint64_t*)nom, rhs.low, 2)) {
                    FP128_INT_DIVIDE_BY_ZERO_EXCEPTION;
                }
            }
            // divide by a 128 bit divisor
            else {
                const uint64_t denom[2] = {rhs.low, rhs.high};
                if (div_32bit((uint32_t*)q, nullptr, (uint32_t*)nom, (uint32_t*)denom, 2ll * array_length(nom), 2ll * array_length(denom))) {
                    FP128_INT_DIVIDE_BY_ZERO_EXCEPTION;
                }
            }
            low = q[0];
            high = q[1];
        }

        // apply sign if needed
        if (res_is_negative) {
            twos_complement128(low, high);
        }

        return *this;
    }
    /**
     * @brief Divide this object by x.
     *
     * The unsigned type keeps a dedicated path for a divisor that fits in 64 bit, which writes the
     * quotient straight into this object. The signed type delegates to the int128_base overload
     * instead: the shortcut for a small numerator relies on an unsigned comparison, which would
     * zero every negative value, and the overload it defers to already special cases a divisor
     * that fits in 64 bit.
     *
     * @param x Right hand side operator (denominator)
     * @return this object.
     */
    template <typename T> FP128_FORCE_INLINE int128_base& operator/=(T x)
    {
        if constexpr (is_signed) {
            return operator/=(int128_base(x));
        } else {
            if (x == 0) {
                FP128_INT_DIVIDE_BY_ZERO_EXCEPTION;
            }

            // floating point values are converted first. Casting them straight to uint64_t below
            // would be undefined behavior for anything at or above 2^64, and a divisor between zero
            // and one would truncate to zero and raise a bogus divide by zero.
            if constexpr (std::is_floating_point_v<T>) {
                return operator/=(int128_base(x));
            }
            // check if the type is signed or not
            // for negative values only, convert to int128_base and divide.
            else if constexpr (std::is_signed_v<T>) {
                if (x < 0)
                    return operator/=(int128_base(x));
            }

            uint64_t uval = static_cast<uint64_t>(x);

            // check some trivial cases
            if (is_zero() || *this < uval) {
                low = high = 0;
                return *this;
            }

            if (*this == uval) {
                low = 1;
                high = 0;
                return *this;
            }

            // exponent of 2, convert to a much faster shift operation
            if (1 == popcnt64(uval)) {
                return *this >>= (int32_t)log2(uval);
            }

            uint64_t nom[2] = {low, high};
            uint64_t* words = &low;
            if (div_64bit(words, nullptr, (uint64_t*)nom, uval, 2)) {
                FP128_INT_DIVIDE_BY_ZERO_EXCEPTION;
            }
            return *this;
        }
    }
    /**
     * @brief %= operator
     *
     * For the signed type the remainder takes the sign of the dividend, matching the builtin
     * integer types: -7 % 2 is -1. It is derived from the quotient so the two can never disagree.
     * The unsigned type gets the remainder straight out of the division instead, which is cheaper.
     *
     * @param rhs Modulo operand.
     * @return This object.
     */
    FP128_INLINE int128_base& operator%=(const int128_base& rhs)
    {
        if (rhs.is_zero()) {
            FP128_INT_DIVIDE_BY_ZERO_EXCEPTION;
        }

        if constexpr (is_signed) {
            // x mod y == x - y * trunc(x / y), which keeps the identity (x/y)*y + x%y == x
            const int128_base quotient = *this / rhs;
            *this -= rhs * quotient;
            return *this;
        } else {
            // check some trivial cases
            if (*this < rhs)
                return *this;

            if (*this == rhs) {
                low = 0;
                high = 0;
                return *this;
            }

            uint64_t q[2] {};
            const uint64_t nom[2] = {low, high};

            // optimization for when dividing by a small integer
            if (rhs.high == 0) {
                if (div_64bit((uint64_t*)q, &low, (uint64_t*)nom, rhs.low, 2)) {
                    FP128_INT_DIVIDE_BY_ZERO_EXCEPTION;
                }
                high = 0;
            } else {
                const uint64_t denom[2] = {rhs.low, rhs.high};
                // div_32bit shrinks the denominator past its leading zero words and fills only that
                // many words of the remainder. Collect the result in a zeroed buffer so the words it
                // leaves untouched read as zero instead of retaining the numerator's bits.
                uint64_t r[2] {};
                if (div_32bit((uint32_t*)q, (uint32_t*)r, (uint32_t*)nom, (uint32_t*)denom, 2ll * array_length(nom), 2ll * array_length(denom))) {
                    FP128_INT_DIVIDE_BY_ZERO_EXCEPTION;
                }
                low = r[0];
                high = r[1];
            }
            return *this;
        }
    }
    /**
     * @brief %= operator
     * @param x Modulo operand.
     * @return This object.
     */
    template <typename T> FP128_FORCE_INLINE int128_base& operator%=(T x) { return operator%=(int128_base(x)); }
    /**
     * @brief Shift right this object.
     * The signed type shifts arithmetically, replicating the sign bit; the unsigned type shifts
     * logically. Both are done without rounding, to match int64_t and uint64_t.
     * @param shift Bits to shift. Zero and negative values leave the object unchanged. Values of
     *              128 and above produce zero, or -1 for a negative value of the signed type.
     * @return This object.
     */
    FP128_INLINE constexpr int128_base& operator>>=(int32_t shift) noexcept
    {
        if (shift < 1)
            return *this;

        if constexpr (is_signed) {
            const int64_t temp = static_cast<int64_t>(high);
            const uint64_t sign_bits = static_cast<uint64_t>(temp >> 63);  // all zeros or all ones

            // 1-63 bit shift - most common
            if (shift < 64) {
                low = FP128_SHIFTRIGHT128(low, high, static_cast<uint8_t>(shift));
                high = static_cast<uint64_t>(temp >> shift);
            }
            // the whole value shifts out, leaving nothing but copies of the sign bit. Without this
            // check the shift below would be taken modulo 64 by the hardware.
            else if (shift > 127) {
                low = high = sign_bits;
            } else {
                low = static_cast<uint64_t>(temp >> (shift - 64));
                high = sign_bits;
            }
        } else {
            // 1-63 bit shift - most common
            if (shift < 64) {
                low = FP128_SHIFTRIGHT128(low, high, static_cast<uint8_t>(shift));
                high >>= shift;
            }
            // the whole value shifts out. Without this check the shift below would be taken modulo 64
            // by the hardware and produce a non zero result.
            else if (shift > 127) {
                low = high = 0;
            } else {
                low = high >> (shift - 64);
                high = 0;
            }
        }
        return *this;
    }
    /**
     * @brief Shift left this object.
     * @param shift Bits to shift. Zero and negative values leave the object unchanged,
     *              values of 128 and above produce zero.
     * @return This object.
     */
    FP128_INLINE constexpr int128_base& operator<<=(int32_t shift) noexcept
    {
        if (shift < 1)
            return *this;
        if (shift < 64) {
            high = FP128_SHIFTLEFT128(low, high, static_cast<uint8_t>(shift));
            low <<= shift;
        }
        // the whole value shifts out. Without this check the shift below would be taken modulo 64
        // by the hardware and produce a non zero result.
        else if (shift > 127) {
            low = high = 0;
        } else {
            high = low << (shift - 64);
            low = 0;
        }
        return *this;
    }
    /**
     * @brief Bitwise AND=
     * @param rhs AND mask.
     * @return This object.
     */
    FP128_FORCE_INLINE constexpr int128_base& operator&=(const int128_base& rhs) noexcept
    {
        low &= rhs.low;
        high &= rhs.high;
        return *this;
    }
    /**
     * @brief Bitwise AND=
     * @param rhs AND mask.
     * @return This object.
     */
    template <typename T> FP128_FORCE_INLINE constexpr int128_base& operator&=(const T& rhs) noexcept { return operator&=(int128_base(rhs)); }
    /**
     * @brief Bitwise OR=
     * @param rhs OR mask.
     * @return This object.
     */
    FP128_FORCE_INLINE constexpr int128_base& operator|=(const int128_base& rhs) noexcept
    {
        low |= rhs.low;
        high |= rhs.high;
        return *this;
    }
    /**
     * @brief Bitwise OR=
     * @param rhs OR mask.
     * @return This object.
     */
    template <typename T> FP128_FORCE_INLINE constexpr int128_base& operator|=(const T& rhs) noexcept { return operator|=(int128_base(rhs)); }
    /**
     * @brief Bitwise XOR=
     * @param rhs XOR mask.
     * @return This object.
     */
    FP128_FORCE_INLINE constexpr int128_base& operator^=(const int128_base& rhs) noexcept
    {
        low ^= rhs.low;
        high ^= rhs.high;
        return *this;
    }
    /**
     * @brief Bitwise XOR=
     * @param rhs XOR mask.
     * @return This object.
     */
    template <typename T> FP128_FORCE_INLINE constexpr int128_base& operator^=(const T& rhs) noexcept { return operator^=(int128_base(rhs)); }
    /**
     * @brief Prefix ++ operation (++a)
     * @return This object.
     */
    FP128_FORCE_INLINE constexpr int128_base& operator++() noexcept
    {
        *this += 1;
        return *this;
    }
    /**
     * @brief Postfix ++ operation (a++)
     * @return This object.
     */
    FP128_FORCE_INLINE constexpr int128_base operator++(int32_t) noexcept
    {
        int128_base temp(*this);
        ++*this;  // call the prefix implementation
        return temp;
    }
    /**
     * @brief Prefix -- operation (--a)
     * @return This object.
     */
    FP128_FORCE_INLINE constexpr int128_base& operator--() noexcept
    {
        *this -= 1;
        return *this;
    }
    /**
     * @brief Postfix -- operation (a--)
     * @return This object.
     */
    FP128_FORCE_INLINE constexpr int128_base operator--(int32_t) noexcept
    {
        int128_base temp(*this);
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
     * @brief Bitwise not (~).
     */
    [[nodiscard]] FP128_FORCE_INLINE constexpr int128_base operator~() const noexcept
    {
        int128_base temp(*this);
        temp.high = ~high;
        temp.low = ~low;
        return temp;
    }
    /**
     * @brief Unary +. Returns a copy of the object.
     */
    [[nodiscard]] FP128_FORCE_INLINE constexpr int128_base operator+() const noexcept
    {
        int128_base temp(*this);
        return temp;
    }
    /**
     * @brief Unary -. Returns a copy of the object with sign inverted.
     * Performs a 2's complement operation just like the native integer types
     */
    [[nodiscard]] FP128_FORCE_INLINE constexpr int128_base operator-() const noexcept
    {
        int128_base temp = *this;
        twos_complement128(temp.low, temp.high);
        return temp;
    }

    //
    // Comparison operators
    //
    /**
     * @brief Compare logical/bitwise equal.
     * @param rhs Righthand operand
     * @return True if this and rhs are equal.
     */
    [[nodiscard]] FP128_FORCE_INLINE constexpr bool operator==(const int128_base& rhs) const noexcept { return high == rhs.high && low == rhs.low; }
    /**
     * @brief Compare logical/bitwise equal.
     * @param rhs Righthand operand
     * @return True if this and rhs are equal.
     */
    template <typename T> [[nodiscard]] FP128_FORCE_INLINE constexpr bool operator==(const T& rhs) const noexcept { return *this == int128_base(rhs); }
    /**
     * @brief Return true when objects are not equal. Can be used as logical XOR.
     * @param rhs Righthand operand.
     * @return True of not equal.
     */
    [[nodiscard]] FP128_FORCE_INLINE constexpr bool operator!=(const int128_base& rhs) const noexcept { return low != rhs.low || high != rhs.high; }
    /**
     * @brief Return true when objects are not equal. Can be used as logical XOR.
     * @param rhs Righthand operand.
     * @return True of not equal.
     */
    template <typename T> [[nodiscard]] FP128_FORCE_INLINE constexpr bool operator!=(const T& rhs) const noexcept { return *this != int128_base(rhs); }
    /**
     * @brief Return true if this object is small than the rhs
     * @param rhs Righthand operand.
     * @return True when this object is smaller.
     */
    [[nodiscard]] FP128_FORCE_INLINE constexpr bool operator<(const int128_base& rhs) const noexcept
    {
        if constexpr (is_signed) {
            const bool sign = is_negative();
            const bool rhs_sign = rhs.is_negative();
            // objects have a different sign
            if (sign != rhs_sign) {
                return sign;
            }
        }

        // objects have the same sign, or the type is unsigned
        return high < rhs.high || (high == rhs.high && low < rhs.low);
    }
    /**
     * @brief Return true if this object is small than the rhs
     * @param rhs Righthand operand.
     * @return True when this object is smaller.
     */
    template <typename T> [[nodiscard]] FP128_FORCE_INLINE constexpr bool operator<(const T& rhs) const noexcept { return *this < int128_base(rhs); }
    /**
     * @brief Return true this object is small or equal than the rhs
     * @param rhs Righthand operand.
     * @return True when this object is smaller or equal.
     */
    [[nodiscard]] FP128_FORCE_INLINE constexpr bool operator<=(const int128_base& rhs) const noexcept { return !(*this > rhs); }
    /**
     * @brief Return true this object is small or equal than the rhs
     * @param rhs Righthand operand.
     * @return True when this object is smaller or equal.
     */
    template <typename T> [[nodiscard]] FP128_FORCE_INLINE constexpr bool operator<=(const T& rhs) const noexcept { return !(*this > int128_base(rhs)); }
    /**
     * @brief Return true this object is larger than the rhs
     * @param rhs Righthand operand.
     * @return True when this objext is larger.
     */
    [[nodiscard]] FP128_FORCE_INLINE constexpr bool operator>(const int128_base& rhs) const noexcept
    {
        if constexpr (is_signed) {
            const bool sign = is_negative();
            const bool rhs_sign = rhs.is_negative();
            // objects have a different sign
            if (sign != rhs_sign) {
                return rhs_sign;
            }
        }

        // objects have the same sign, or the type is unsigned
        return high > rhs.high || (high == rhs.high && low > rhs.low);
    }
    /**
     * @brief Return true this object is larger than the rhs
     * @param rhs Righthand operand.
     * @return True when this objext is larger.
     */
    template <typename T> [[nodiscard]] FP128_FORCE_INLINE constexpr bool operator>(const T& rhs) const noexcept { return *this > int128_base(rhs); }
    /**
     * @brief Return true this object is larger or equal than the rhs
     * @param rhs Righthand operand.
     * @return True when this objext is larger or equal.
     */
    [[nodiscard]] FP128_FORCE_INLINE constexpr bool operator>=(const int128_base& rhs) const noexcept { return !(*this < rhs); }
    /**
     * @brief Return true this object is larger or equal than the rhs
     * @param rhs Righthand operand.
     * @return True when this objext is larger or equal.
     */
    template <typename T> [[nodiscard]] FP128_FORCE_INLINE constexpr bool operator>=(const T& rhs) const noexcept { return !(*this < int128_base(rhs)); }

    //
    // useful public functions
    //
    /**
     * @brief Returns true if the value is an int (fraction is zero)
     * @return True when the fraction is zero.
     */
    [[nodiscard]] FP128_FORCE_INLINE constexpr bool is_int() const noexcept { return true; }
    /**
     * @brief Returns true if the value is positive (including zero)
     * Always true for the unsigned type, which lets the sign tests in the shared code fold away.
     * @return True when the value is positive
     */
    [[nodiscard]] FP128_FORCE_INLINE constexpr bool is_positive() const noexcept
    {
        if constexpr (is_signed) {
            return 0ull == high >> 63;
        } else {
            return true;
        }
    }
    /**
     * @brief Returns true if the value negative (smaller than zero)
     * Always false for the unsigned type, which lets the sign tests in the shared code fold away.
     * @return True when the value is negative
     */
    [[nodiscard]] FP128_FORCE_INLINE constexpr bool is_negative() const noexcept
    {
        if constexpr (is_signed) {
            return 1ull == high >> 63;
        } else {
            return false;
        }
    }
    /**
     * @brief Returns true if the value is zero
     * @return Returns true if the value is zero
     */
    [[nodiscard]] FP128_FORCE_INLINE constexpr bool is_zero() const noexcept { return 0 == low && 0 == high; }
    /**
     * @brief get a specific bit within the 128 data
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
     * @brief Extracts the low and high QWORD
     * @param l Reference to the low QWORD
     * @param h Reference to the high QWORD
     */
    FP128_FORCE_INLINE constexpr void get_components(uint64_t& l, uint64_t& h) const noexcept
    {
        l = low;
        h = high;
    }
    /**
     * @brief Return an instance of int128_base with the value of 1
     * @return 1
     */
    [[nodiscard]] FP128_FORCE_INLINE static constexpr int128_base one() noexcept { return int128_base(1u); }
    /**
     * @brief Converts this object to a hex C string.
     * The returned string is a statically thread-allocated buffer.
     * Additional calls to this function from the same thread overwrite the previous result.
     * @return C string containing the HEX value of the object.
     */
    [[nodiscard]] FP128_INLINE char* hex() const noexcept
    {
        constexpr int buff_size = 35;
        static thread_local char str[buff_size];
        snprintf(str, buff_size, "0x%016llX%016llX", high, low);
        return str;
    }

    //
    // Binary math operators
    //
    // Each of these applies the compound assignment to the by value left hand side and then returns
    // it on a line of its own. The shorter `return lhs OP= rhs;` is equivalent and was what these
    // used to say, but it costs MSVC a factor of six on the multiplication: the compound assignment
    // returns int128_base&, so the return statement copy constructs the result from an lvalue, and
    // MSVC emits that copy as a 16 byte XMM load of the object whose two halves it has just written
    // with 8 byte GPR stores. A load that overlaps two narrower stores cannot forward, so every
    // evaluation ends in a store forwarding stall of about fifteen cycles - which is five times the
    // multiplication itself. Returning the named parameter instead keeps it in registers.
    // Clang generates identical code either way. Measured: uint128_t * uint128_t goes from 259M/s
    // to 1.73G/s on MSVC 19.51, matching Clang.
    //

    /**
     * @brief Adds 2 values and returns the result.
     * @param lhs left hand side operand
     * @param rhs Right hand side operand
     * @return Result of the operation
     */
    template <typename T> [[nodiscard]] friend FP128_FORCE_INLINE constexpr int128_base operator+(int128_base lhs, const T& rhs) noexcept
    {
        lhs += rhs;
        return lhs;
    }
    /**
     * @brief Subtracts 2 values and returns the result.
     * @param lhs left hand side operand
     * @param rhs Right hand side operand
     * @return Result of the operation
     */
    template <typename T> [[nodiscard]] friend FP128_FORCE_INLINE constexpr int128_base operator-(int128_base lhs, const T& rhs) noexcept
    {
        lhs -= rhs;
        return lhs;
    }
    /**
     * @brief Multiplies 2 values and returns the result.
     * @param lhs left hand side operand
     * @param rhs Right hand side operand
     * @return Result of the operation
     */
    template <typename T> [[nodiscard]] friend FP128_FORCE_INLINE constexpr int128_base operator*(int128_base lhs, const T& rhs) noexcept
    {
        lhs *= rhs;
        return lhs;
    }
    /**
     * @brief Divides 2 values and returns the result.
     *
     * Forced inline like every other forwarder here, see FP128_FORCE_INLINE. This one is where the
     * cost was first measured: leaving the decision to MSVC under /GL (whole program optimization)
     * lands on a round trip through memory worth ~10% of a 128 bit by 64 bit division.
     *
     * @param lhs left hand side operand
     * @param rhs Right hand side operand
     * @return Result of the operation
     */
    template <typename T> [[nodiscard]] friend FP128_FORCE_INLINE int128_base operator/(int128_base lhs, const T& rhs)
    {
        lhs /= rhs;
        return lhs;
    }
    /**
     * @brief Performs modulo and returns the result.
     *
     * @param lhs left hand side operand
     * @param rhs Right hand side operand
     * @return Result of the operation
     */
    template <typename T> [[nodiscard]] friend FP128_FORCE_INLINE int128_base operator%(int128_base lhs, const T& rhs)
    {
        lhs %= rhs;
        return lhs;
    }

    //
    // Binary math operators
    //
    /**
     * @brief Performs bitwise AND (&).
     * @param lhs left hand side operand
     * @param rhs Right hand side operand
     * @return Result of the operation
     */
    template <typename T> [[nodiscard]] friend FP128_FORCE_INLINE constexpr int128_base operator&(int128_base lhs, const T& rhs)
    {
        lhs &= rhs;
        return lhs;
    }
    /**
     * @brief Performs bitwise OR (|).
     * @param lhs left hand side operand
     * @param rhs Right hand side operand
     * @return Result of the operation
     */
    template <typename T> [[nodiscard]] friend FP128_FORCE_INLINE constexpr int128_base operator|(int128_base lhs, const T& rhs)
    {
        lhs |= rhs;
        return lhs;
    }
    /**
     * @brief Performs bitwise XOR (^).
     * @param lhs left hand side operand
     * @param rhs Right hand side operand
     * @return Result of the operation
     */
    template <typename T> [[nodiscard]] friend FP128_FORCE_INLINE constexpr int128_base operator^(int128_base lhs, const T& rhs)
    {
        lhs ^= rhs;
        return lhs;
    }

    //
    // Binary math operators with the scalar on the left hand side
    //
    // Without these, an expression like (1 + x) is ambiguous: converting the literal to
    // int128_base and converting x to a builtin type are both one user defined conversion, so
    // neither overload wins. Restricting the left operand to the arithmetic types keeps these from
    // competing with the int128_base on the left versions above, which would otherwise be an
    // equally good match.
    //

    /// @brief Adds a scalar and an int128_base, in that order. @param lhs Left operand @param rhs Right operand @return Result of the operation
    template <typename T>
        requires std::is_arithmetic_v<T>
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr int128_base operator+(const T& lhs, const int128_base& rhs) noexcept
    {
        return int128_base(lhs) += rhs;
    }
    /// @brief Subtracts an int128_base from a scalar. @param lhs Left operand @param rhs Right operand @return Result of the operation
    template <typename T>
        requires std::is_arithmetic_v<T>
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr int128_base operator-(const T& lhs, const int128_base& rhs) noexcept
    {
        return int128_base(lhs) -= rhs;
    }
    /// @brief Multiplies a scalar and an int128_base, in that order. @param lhs Left operand @param rhs Right operand @return Result of the operation
    template <typename T>
        requires std::is_arithmetic_v<T>
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr int128_base operator*(const T& lhs, const int128_base& rhs) noexcept
    {
        return int128_base(lhs) *= rhs;
    }
    /// @brief Divides a scalar by an int128_base. @param lhs Left operand @param rhs Right operand @return Result of the operation
    template <typename T>
        requires std::is_arithmetic_v<T>
    [[nodiscard]] friend FP128_FORCE_INLINE int128_base operator/(const T& lhs, const int128_base& rhs)
    {
        return int128_base(lhs) /= rhs;
    }
    /// @brief Performs modulo of a scalar by an int128_base. @param lhs Left operand @param rhs Right operand @return Result of the operation
    template <typename T>
        requires std::is_arithmetic_v<T>
    [[nodiscard]] friend FP128_FORCE_INLINE int128_base operator%(const T& lhs, const int128_base& rhs)
    {
        return int128_base(lhs) %= rhs;
    }
    /// @brief Performs bitwise AND of a scalar and an int128_base. @param lhs Left operand @param rhs Right operand @return Result of the operation
    template <typename T>
        requires std::is_arithmetic_v<T>
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr int128_base operator&(const T& lhs, const int128_base& rhs)
    {
        return int128_base(lhs) &= rhs;
    }
    /// @brief Performs bitwise OR of a scalar and an int128_base. @param lhs Left operand @param rhs Right operand @return Result of the operation
    template <typename T>
        requires std::is_arithmetic_v<T>
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr int128_base operator|(const T& lhs, const int128_base& rhs)
    {
        return int128_base(lhs) |= rhs;
    }
    /// @brief Performs bitwise XOR of a scalar and an int128_base. @param lhs Left operand @param rhs Right operand @return Result of the operation
    template <typename T>
        requires std::is_arithmetic_v<T>
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr int128_base operator^(const T& lhs, const int128_base& rhs)
    {
        return int128_base(lhs) ^= rhs;
    }
    /**
     * @brief Shifts a scalar left by an int128_base shift count.
     *
     * The left operand is widened to int128_base first, so the result is 128 bit rather than the
     * builtin type of lhs. This is what makes (1 << n) with a large n produce the expected power of
     * two instead of the undefined behavior a builtin shift past the operand width would have.
     *
     * @param lhs Left operand, the value being shifted
     * @param rhs Right operand, the shift count. Only its low 32 bits are used.
     * @return Result of the operation
     */
    template <typename T>
        requires std::is_arithmetic_v<T>
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr int128_base operator<<(const T& lhs, const int128_base& rhs) noexcept
    {
        return int128_base(lhs) <<= static_cast<int32_t>(rhs);
    }
    /**
     * @brief Shifts a scalar right by an int128_base shift count.
     * The left operand is widened to int128_base first, see operator<< above.
     * @param lhs Left operand, the value being shifted
     * @param rhs Right operand, the shift count. Only its low 32 bits are used.
     * @return Result of the operation
     */
    template <typename T>
        requires std::is_arithmetic_v<T>
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr int128_base operator>>(const T& lhs, const int128_base& rhs) noexcept
    {
        return int128_base(lhs) >>= static_cast<int32_t>(rhs);
    }

    //
    // Comparison operators with the scalar on the left hand side
    //
    // The member comparisons above only cover an int128_base on the left. Without these, an
    // expression like (1 == x) has nothing but the builtin comparisons to choose from, and the
    // seven conversion operators of int128_base make every one of them equally good, so the call is
    // ambiguous. Widening the scalar keeps the comparison exact at 128 bits instead of narrowing
    // the object down to whichever builtin type the compiler happened to pick.
    //

    /// @brief Compares a scalar and an int128_base for equality. @param lhs Left operand @param rhs Right operand @return True when the two are equal
    template <typename T>
        requires std::is_arithmetic_v<T>
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr bool operator==(const T& lhs, const int128_base& rhs) noexcept
    {
        return int128_base(lhs) == rhs;
    }
    /// @brief Compares a scalar and an int128_base for inequality. @param lhs Left operand @param rhs Right operand @return True when the two differ
    template <typename T>
        requires std::is_arithmetic_v<T>
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr bool operator!=(const T& lhs, const int128_base& rhs) noexcept
    {
        return int128_base(lhs) != rhs;
    }
    /// @brief Returns true when a scalar is smaller than an int128_base. @param lhs Left operand @param rhs Right operand @return True when lhs is smaller
    template <typename T>
        requires std::is_arithmetic_v<T>
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr bool operator<(const T& lhs, const int128_base& rhs) noexcept
    {
        return int128_base(lhs) < rhs;
    }
    /// @brief Returns true when a scalar is smaller or equal to an int128_base. @param lhs Left operand @param rhs Right operand @return True when lhs is smaller or equal
    template <typename T>
        requires std::is_arithmetic_v<T>
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr bool operator<=(const T& lhs, const int128_base& rhs) noexcept
    {
        return int128_base(lhs) <= rhs;
    }
    /// @brief Returns true when a scalar is larger than an int128_base. @param lhs Left operand @param rhs Right operand @return True when lhs is larger
    template <typename T>
        requires std::is_arithmetic_v<T>
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr bool operator>(const T& lhs, const int128_base& rhs) noexcept
    {
        return int128_base(lhs) > rhs;
    }
    /// @brief Returns true when a scalar is larger or equal to an int128_base. @param lhs Left operand @param rhs Right operand @return True when lhs is larger or equal
    template <typename T>
        requires std::is_arithmetic_v<T>
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr bool operator>=(const T& lhs, const int128_base& rhs) noexcept
    {
        return int128_base(lhs) >= rhs;
    }

    //
    // Various math functions, implemented as friend functions
    //
    // These are hidden friends: they are found by argument dependent lookup only, never by
    // ordinary unqualified lookup. That is what keeps abs() below from hiding ::abs() for the
    // builtin types inside namespace fp128.
    //

    /**
     * @brief Calculates the absolute value of x.
     * Only declared for the signed instantiation. An unqualified abs() on an unsigned value keeps
     * resolving the way it does today.
     * @param x input value.
     * @return |x|
     */
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr int128_base abs(const int128_base& x) noexcept
        requires IsSigned
    {
        return x.is_positive() ? x : -x;
    }
    /**
     * @brief Calculates the square of a value. i.e. x^2
     *
     * Faster than x * x and bit identical to it, see int128_base::square().
     *
     * @param x Value to square
     * @return x^2, truncated to 128 bit
     */
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr int128_base sqr(int128_base x) noexcept { return x.square(); }
    /**
     * @brief Calculates the left zero count of value x.
     * @param x input value.
     * @return lzc (uint32_t) of the result.
     */
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr uint64_t lzcnt128(const int128_base& x) noexcept
    {
        return (x.high != 0) ? lzcnt64(x.high) : 64 + lzcnt64(x.low);
    }
    /**
     * @brief Calculates the square root using Newton's method.
     * Based on the book "Math toolkit for real time programming" by Jack W. Crenshaw
     * @param x Value to calculate the root of
     * @return Square root of (x) rounded down to the nearest integer, zero when x <= 0.
     */
    [[nodiscard]] friend FP128_INLINE uint64_t sqrt(const int128_base& x) noexcept
    {
        // zero has to be rejected here as well as the negatives. log2 throws on a zero input and
        // this function is noexcept, so letting it through would terminate the process.
        // The negative test is a compile time false for the unsigned instantiation.
        if (x.is_negative() || x.is_zero())
            return 0;

        // an msb at bit zero means the value is 1, which is its own root and would make the
        // initial guess below zero
        const auto expo = static_cast<uint32_t>(log2(x));
        if (expo == 0) {
            return 1;
        }

        int128_base root = int128_base::one();
        int128_base e, temp;
        root <<= ((expo + 1) >> 1);

        // Newton iterations to reduce the error
        do {
            temp = x / root;
            // working with unsigned numbers, must keep the below positive at all times
            e = ((root > temp) ? (root - temp) : (temp - root)) >> 1;
            root = (root + temp) >> 1;
        } while (e);

        return root;
    }
    /**
     * @brief Calculates the Log base 2 of x: log2(x)
     * Rounding is always towards zero so the maximum error is close to 1.
     * @param x The number to perform log2 on.
     * @return log2(x), which is also the bit position of the msb.
     * @throws std::domain_error when x is zero, or negative for the signed type.
     */
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr uint64_t log2(const int128_base& x)
    {
        if (x.is_negative() || x.is_zero()) {
            throw std::domain_error("Math domain error! Function accepts positive, non-zero values only.");
        }

        return 127ull - lzcnt128(x);
    }
    /**
     * @brief Calculates the natural Log (base e) of x: log(x), rounded down to the nearest integer.
     * @param x The number to perform log on.
     * @return log(x)
     * @throws std::domain_error when x is zero, or negative for the signed type.
     */
    [[nodiscard]] friend FP128_FORCE_INLINE uint64_t log(const int128_base& x)
    {
        if (x.is_negative() || x.is_zero()) {
            throw std::domain_error("Math domain error! Function accepts positive, non-zero values only.");
        }

        // The table below holds the value of various powers of e (~2.718), index i holds
        // ceil(pow(e, i)). That is the smallest integer whose log is i: e^i is never itself an
        // integer, so an integer x satisfies x >= e^i exactly when x >= ceil(e^i). Every entry
        // therefore has to be exact to the last digit, a value one too large misclassifies the
        // boundary itself and one too small misclassifies everything up to the true boundary.
        // The entries are produced with exact arithmetic rather than pow(), which loses the low
        // digits well before e^88.<BR>
        // Being a function local static with an initializer, the language guarantees its
        // construction happens once even when several threads race to get here.
        static const int128_base lan_table[] = {
            "1",                                       // e^0
            "3",                                       // e^1
            "8",                                       // e^2
            "21",                                      // e^3
            "55",                                      // e^4
            "149",                                     // e^5
            "404",                                     // e^6
            "1097",                                    // e^7
            "2981",                                    // e^8
            "8104",                                    // e^9
            "22027",                                   // e^10
            "59875",                                   // e^11
            "162755",                                  // e^12
            "442414",                                  // e^13
            "1202605",                                 // e^14
            "3269018",                                 // e^15
            "8886111",                                 // e^16
            "24154953",                                // e^17
            "65659970",                                // e^18
            "178482301",                               // e^19
            "485165196",                               // e^20
            "1318815735",                              // e^21
            "3584912847",                              // e^22
            "9744803447",                              // e^23
            "26489122130",                             // e^24
            "72004899338",                             // e^25
            "195729609429",                            // e^26
            "532048240602",                            // e^27
            "1446257064292",                           // e^28
            "3931334297145",                           // e^29
            "10686474581525",                          // e^30
            "29048849665248",                          // e^31
            "78962960182681",                          // e^32
            "214643579785917",                         // e^33
            "583461742527455",                         // e^34
            "1586013452313431",                        // e^35
            "4311231547115196",                        // e^36
            "11719142372802612",                       // e^37
            "31855931757113757",                       // e^38
            "86593400423993747",                       // e^39
            "235385266837019986",                      // e^40
            "639843493530054950",                      // e^41
            "1739274941520501048",                     // e^42
            "4727839468229346562",                     // e^43
            "12851600114359308276",                    // e^44
            "34934271057485095349",                    // e^45
            "94961194206024488746",                    // e^46
            "258131288619006739624",                   // e^47
            "701673591209763173866",                   // e^48
            "1907346572495099690526",                  // e^49
            "5184705528587072464088",                  // e^50
            "14093490824269387964493",                 // e^51
            "38310080007165768493036",                 // e^52
            "104137594330290877971835",                // e^53
            "283075330327469390044207",                // e^54
            "769478526514201713818275",                // e^55
            "2091659496012996153907072",               // e^56
            "5685719999335932222640349",               // e^57
            "15455389355901039303530767",              // e^58
            "42012104037905142549565935",              // e^59
            "114200738981568428366295719",             // e^60
            "310429793570191990870734215",             // e^61
            "843835666874145448907332949",             // e^62
            "2293783159469609879099352841",            // e^63
            "6235149080811616882909238709",            // e^64
            "16948892444103337141417836115",           // e^65
            "46071866343312915426773184429",           // e^66
            "125236317084221378051352196075",          // e^67
            "340427604993174052137690718701",          // e^68
            "925378172558778760024239791669",          // e^69
            "2515438670919167006265781174253",         // e^70
            "6837671229762743866755892826678",         // e^71
            "18586717452841279803403701812546",        // e^72
            "50523936302761041945570383321858",        // e^73
            "137338297954017618778418852980854",       // e^74
            "373324199679900164025490831726471",       // e^75
            "1014800388113888727832461784131717",      // e^76
            "2758513454523170206286469819902662",      // e^77
            "7498416996990120434675630591224061",      // e^78
            "20382810665126687668323137537172633",     // e^79
            "55406223843935100525711733958316613",     // e^80
            "150609731458503054835259413016767499",    // e^81
            "409399696212745469666091422932782905",    // e^82
            "1112863754791759412087071478183940806",   // e^83
            "3025077322201142338266566396443428743",   // e^84
            "8223012714622913510304328016407774696",   // e^85
            "22352466037347150474430657323327147399",  // e^86
            "60760302250568721495223289381302760753",  // e^87
            "165163625499400185552832979626485876707"  // e^88
        };
        static constexpr uint64_t lan_table_len = array_length(lan_table);
        if (x < 3)
            return 0;
        if (x >= lan_table[lan_table_len - 1])
            return lan_table_len - 1;

        // binary search the result
        uint64_t l = 0, h = lan_table_len - 1, res = (h + l) >> 1;

        while (l < h) {
            if (x < lan_table[res]) {
                h = res;
            } else if (l == res)
                break;
            else {
                l = res;
            }
            res = (h + l) >> 1;
        }
        return res;
    }
    /**
     * @brief Calculates Log base 10 of x: log10(x), rounded down to the nearest integer.
     * @param x The number to perform log on.
     * @return log10(x)
     * @throws std::domain_error when x is zero, or negative for the signed type.
     */
    [[nodiscard]] friend FP128_FORCE_INLINE uint64_t log10(const int128_base& x)
    {
        if (x.is_negative() || x.is_zero()) {
            throw std::domain_error("Math domain error! Function accepts positive, non-zero values only.");
        }

        static constexpr uint64_t log10_max = 38;  // 10^38 is the largest power of ten below 2^127
        // Holds every power of ten that fits in 128 bit. Building it inside the initializer of a
        // function local static makes the first call thread safe: the language guarantees that
        // concurrent callers block until the initialization finishes. Filling a plain static array
        // behind an 'is it still zero' check does not, both threads can pass the check at once.
        static const auto log10_table = []() {
            std::array<int128_base, log10_max + 1> table;
            table[0] = 1ull;
            for (uint32_t i = 1u; i <= log10_max; ++i) {
                table[i] = table[i - 1] * 10ull;
            }
            return table;
        }();

        if (x >= log10_table[log10_max])
            return log10_max;

        // binary search the result
        uint64_t l = 0, h = log10_max, res = (h + l) >> 1;

        while (l < h) {
            if (x < log10_table[res]) {
                h = res;
            } else if (l == res)
                break;
            else {
                l = res;
            }
            res = (h + l) >> 1;
        }
        return res;
    }
    /**
     * @brief Calculates x to the power of y.
     * @param x Base value
     * @param y Exponent
     * @return x to the power of y
     */
    [[nodiscard]] friend FP128_INLINE constexpr int128_base pow(const int128_base& x, uint32_t y) noexcept
    {
        // zero power always yields 1
        // Even if x is zero! Same behavior as pow(double, double) but different from Python which returns zero.
        if (y == 0)
            return 1;

        // special case where base is zero
        if (!x)
            return x;

        int128_base res = 1;
        int128_base b = x;
        while (y > 0) {
            if (y & 1)
                res *= b;
            y >>= 1;
            b.square();
        }
        return res;
    }
};  // class int128_base


/**
 * @brief Writes a 128 bit integer to a stream, honouring its width, fill and adjustment.
 *
 * The value is always written in decimal: that is the conversion the type provides, and a stream's
 * hex and oct flags have no counterpart in it.
 */
template <bool IsSigned> inline std::ostream& operator<<(std::ostream& os, const int128_base<IsSigned>& value)
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
 * @brief Reads a 128 bit integer from a stream.
 *
 * Accepts what the constructor from a string does: an optional sign followed by decimal digits, or
 * a hexadecimal literal after an 0x prefix.
 */
template <bool IsSigned> inline std::istream& operator>>(std::istream& is, int128_base<IsSigned>& value)
{
    std::string token;
    if (!(is >> token))
        return is;

    value = int128_base<IsSigned>(token.c_str());
    return is;
}

}  // namespace fp128



namespace std
{
/**
 * @brief Numeric properties of the 128 bit integer types.
 *
 * Specialized so that generic code written against a builtin integer - anything reaching for
 * numeric_limits<T>::max() to seed a minimum, or for digits10 to size a buffer - compiles and
 * behaves correctly when instantiated with int128_t or uint128_t.
 */
template <bool IsSigned> class numeric_limits<fp128::int128_base<IsSigned>>
{
    using value_type = fp128::int128_base<IsSigned>;

public:
    static constexpr bool is_specialized = true;
    static constexpr bool is_signed = IsSigned;
    static constexpr bool is_integer = true;
    static constexpr bool is_exact = true;
    static constexpr bool has_infinity = false;
    static constexpr bool has_quiet_NaN = false;
    static constexpr bool has_signaling_NaN = false;
    static constexpr bool has_denorm_loss = false;
    static constexpr float_denorm_style has_denorm = denorm_absent;
    static constexpr float_round_style round_style = round_toward_zero;
    static constexpr bool is_iec559 = false;
    static constexpr bool is_bounded = true;
    /// @brief Arithmetic wraps around, which is what the truncated 128 bit operators do.
    static constexpr bool is_modulo = !IsSigned;
    static constexpr bool traps = true;
    static constexpr bool tinyness_before = false;

    /// @brief Value bits, which excludes the sign bit for the signed type.
    static constexpr int digits = IsSigned ? 127 : 128;
    /// @brief Decimal digits that can be represented without change: floor(digits * log10(2)).
    static constexpr int digits10 = IsSigned ? 38 : 38;
    static constexpr int max_digits10 = 0;
    static constexpr int radix = 2;
    static constexpr int min_exponent = 0;
    static constexpr int max_exponent = 0;
    static constexpr int min_exponent10 = 0;
    static constexpr int max_exponent10 = 0;

    [[nodiscard]] static constexpr value_type min() noexcept
    {
        // The signed minimum is the sign bit on its own; the unsigned one is zero.
        return IsSigned ? value_type(0ull, 1ull << 63) : value_type(0ull, 0ull);
    }
    [[nodiscard]] static constexpr value_type max() noexcept
    {
        return IsSigned ? value_type(UINT64_MAX, UINT64_MAX >> 1) : value_type(UINT64_MAX, UINT64_MAX);
    }
    [[nodiscard]] static constexpr value_type lowest() noexcept { return min(); }
    [[nodiscard]] static constexpr value_type epsilon() noexcept { return value_type(); }
    [[nodiscard]] static constexpr value_type round_error() noexcept { return value_type(); }
    [[nodiscard]] static constexpr value_type infinity() noexcept { return value_type(); }
    [[nodiscard]] static constexpr value_type quiet_NaN() noexcept { return value_type(); }
    [[nodiscard]] static constexpr value_type signaling_NaN() noexcept { return value_type(); }
    [[nodiscard]] static constexpr value_type denorm_min() noexcept { return value_type(); }
};

/// @brief const, volatile and cv qualified 128 bit integers have the same numeric properties.
template <bool IsSigned> class numeric_limits<const fp128::int128_base<IsSigned>> : public numeric_limits<fp128::int128_base<IsSigned>>
{
};
template <bool IsSigned> class numeric_limits<volatile fp128::int128_base<IsSigned>> : public numeric_limits<fp128::int128_base<IsSigned>>
{
};
template <bool IsSigned>
class numeric_limits<const volatile fp128::int128_base<IsSigned>> : public numeric_limits<fp128::int128_base<IsSigned>>
{
};

/// @brief Hash support, so a 128 bit integer can be a key in an unordered container.
template <bool IsSigned> struct hash<fp128::int128_base<IsSigned>>
{
    [[nodiscard]] size_t operator()(const fp128::int128_base<IsSigned>& value) const noexcept
    {
        uint64_t low = 0, high = 0;
        value.get_components(low, high);

        // splitmix64's finalizer over the two halves
        uint64_t state = low + 0x9E3779B97F4A7C15ull;
        state = (state ^ (state >> 30)) * 0xBF58476D1CE4E5B9ull;
        state = (state ^ (state >> 27)) * 0x94D049BB133111EBull;
        state ^= high + 0x9E3779B97F4A7C15ull + (state << 6) + (state >> 2);
        state = (state ^ (state >> 30)) * 0xBF58476D1CE4E5B9ull;
        return static_cast<size_t>(state ^ (state >> 31));
    }
};

/**
 * @brief std::format support for the 128 bit integer types.
 *
 * Accepts fill and alignment, a sign, a width and the d type. The other integer presentations that
 * `<format>` defines for a builtin - b, o, x and c - are not offered, because the conversion the
 * type provides is decimal.
 */
template <bool IsSigned> struct formatter<fp128::int128_base<IsSigned>, char>
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

        if (it != end && *it == 'd')
            ++it;

        if (it != end && *it != '}')
            throw format_error("invalid format specification for a 128 bit integer");
        return it;
    }

    template <typename FormatContext> auto format(const fp128::int128_base<IsSigned>& value, FormatContext& context) const
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
#endif  // FP128_INT128_SHARED_H
