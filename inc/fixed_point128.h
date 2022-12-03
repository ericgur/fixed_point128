/***********************************************************************************
    MIT License

    Copyright (c) 2022 Eric Gur (ericgur@iname.com)

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
    by Jack W. Crenshaw

************************************************************************************/

#ifndef FIXED_POINT128_H
#define FIXED_POINT128_H

#include "fixed_point128_shared.h"

namespace fp128 {

/***********************************************************************************
*                                  Main Code
************************************************************************************/

/**
 * @brief 128 bit fixed point class template.
 * 
 * This template provides a floating point-like type that performs various math operations quickly compared to traditional high precision libraries.
 * The only template parameter <B>I</B> is the bit count of the integer part. The fraction part complements <B>I</B> to 128 bit.<BR>
 * <B>I</B> is limited to the range [1,64] in order to simplify the implementation and increase preformance.<BR>
 * This restriction is enforced at compile time.
 * All of fixed_point128's methods are inline for maximum performance.
 *
 * <B>Implementation notes:</B>
 * <UL>
 * <LI>Overflow is handled silently, similar to builtin integer operations.</LI>
 * <LI>Unlike integer operations, the sign is held in a separate data member which simplifies the implementation.</LI>
 * <LI>A fixed_point128 object is not thread safe. Accessing a const object from multiple threads is safe.</LI>
 * <LI>fixed_point128 is <B>conditionally safe</B>, 2 different non const objects can be accessed concurrently.</LI>
 * <LI>Only 64 bit builds are supported.</LI>
 * </UL>
*/
template<int32_t I>
class fixed_point128
{
    // build time validation of template parameters
    static_assert(1 <= I && I <= 64, "Template parameter <I> must be in the [1,64] range!");
    static_assert(sizeof(void*) == 8, "fixed_point128 is supported in 64 bit builds only!");

    // friends
    friend class fixed_point128; // this class is a friend of all its template instances. Avoids awkward getter/setter functions.
private:
    //
    // members
    //
    uint64_t low;  // lower QWORD
    uint64_t high; // upper QWORD
    unsigned sign; // 0 = positive, 1 negative

    // useful const calculations
    static constexpr int32_t F = 128 - I;                               // fraction bit count
    static constexpr int32_t upper_frac_bits = F - 64;                  // how many bits of the fraction exist in the upper QWORD
    static constexpr uint64_t unity = 1ull << (64 - I);                 // upper QWORD value equal to '1'
    static inline const double upper_unity = pow(2, 64 - F);     // convert upper QWORD to floating point
    static inline const double lower_unity = pow(2, -F);         // convert lower QWORD to floating point
    static constexpr uint64_t int_mask = UINT64_MAX << upper_frac_bits; // mask of the integer bits in the upper QWORD
    static constexpr int32_t max_frac_digits = (int)(1 + F / 3.1);      // meaningful base 10 digits for the fraction
public:
    typedef fixed_point128<I> type;
    typedef fixed_point128<I>* ptr_type;
    typedef fixed_point128<I>& ref_type;

    //
    // ctors
    //

    /**
     * @brief Default constructor, creates an instance with a value of zero.
    */
    fixed_point128() noexcept :
        low(0), high(0), sign(0) {}
    /**
     * @brief Copy constructor
     * @param other Object to copy from
    */
    fixed_point128(const fixed_point128& other) noexcept :
        low(other.low), high(other.high), sign(other.sign) {}
    /**
     * @brief Move constructor
     * Doesn't modify the right hand side object. Acts like a copy constructor.
     * @param other Object to copy from
    */
    fixed_point128(const fixed_point128&& other) noexcept :
        low(other.low), high(other.high), sign(other.sign) {}
    /**
     * @brief Constructor from the double type
     * Underflow goes to zero. Overflow, NaN and +-INF go to max supported positive value.
     * @param x Input value
    */
    fixed_point128(double x) noexcept {
        // brute convert to uint64_t for easy bit manipulation
        const uint64_t i = *reinterpret_cast<uint64_t*>(&x);
        // very common case
        if (i == 0) {
            low = high = 0;
            sign = 0;
            return;
        }

        sign = FP128_GET_BIT(i, 63);
        const int32_t e = FP128_GET_BITS(i, dbl_frac_bits, dbl_exp_bits) - 1023;
        uint64_t f = (i & FP128_MAX_VALUE_64(dbl_frac_bits));
        
        // overflow which catches NaN and Inf
        if (e >= I) {
            high = low = UINT64_MAX;
            sign = 0;
            return;
        }

        // normal number, produces non zero value
        if (e >= -F) {
            // bit 52 in f is the unity value of the float. it needs to move to the unity position in fixed point
            f |= FP128_ONE_SHIFT(dbl_frac_bits);
            int32_t bits_to_shift = 64 - I - dbl_frac_bits + e;

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
                    bits_to_shift -= dbl_frac_bits;
                    f <<= 63 - dbl_frac_bits;
                    low = f >> bits_to_shift;
                }
            }
        }
        // too small to be represented, no need to bother.
        else {
            high = low = 0;
            sign = 0;
        }
    }
    /**
     * @brief Constructor from uint64_t type
     * @param x Input value
    */
    fixed_point128(uint64_t x) noexcept {
        low = 0;
        sign = 0;
        high = x << upper_frac_bits;
    }
    /**
     * @brief Constructor from int64_t type
     * @param x Input value
    */
    fixed_point128(int64_t x) noexcept {
        low = 0;
        sign = FP128_GET_BIT(x, 63);
        high = ((sign != 0) ? -x : x) << upper_frac_bits;
    }
    /**
     * @brief Constructor from uint32_t type
     * @param x Input value
   */
    fixed_point128(uint32_t x) noexcept {
        low = 0;
        sign = 0;
        high = (uint64_t)x << upper_frac_bits;
    }
    /**
     * @brief Constructor from int32_t type
     * @param x Input value
   */
    fixed_point128(int32_t x) noexcept {
        low = 0;
        sign = FP128_GET_BIT(x, 31);
        high = (uint64_t)((sign != 0) ? -x : x) << upper_frac_bits;
    }
    /**
     * @brief Constructor from const char* (C string).
     * Accurate up to 37 digits after the decimal point.
     * Allows creating very high precision values. Much slower than the other constructors.
     * @param x Input string
    */
    fixed_point128(const char* x) noexcept {
        low = high = 0;
        sign = 0;
        if (x == nullptr) return;
        
        char* str = _strdup(x);
        char* p = str;
        if (p == nullptr) return;

        // skip white space
        while (*p == ' ') ++p;

        if (*p == '\0') {
            *this = fixed_point128();
            free(str);
            return;
        }

        // set negative sign if needed
        if (p[0] == '-') {
            sign = 1;
            ++p;
        }
        char* dec = strchr(p, '.');
        // number is an integer
        if (dec == nullptr) {
            high = std::strtoull(p, nullptr, 10) << upper_frac_bits;
            free(str);
            return;
        }

        // number is a float, get the integer part using strtoull()
        *dec = '\0';
        uint64_t int_val = std::strtoull(p, nullptr, 10) << upper_frac_bits;

        p = dec + 1;
        int32_t digits = 0;
        fixed_point128<1> base(0xCCCCCCCCCCCCCCCD, 0x0CCCCCCCCCCCCCCC, 0); // maximum precision to represent 0.1
        fixed_point128<1> step = base;
        fixed_point128<1> frac;
        // multiply each digits by 0.1**n
        while (digits++ < max_frac_digits && *p != '\0' && base) {
            fixed_point128<1> temp = base * (p[0] - '0');
            unsigned char carry = _addcarry_u64(0, frac.low, temp.low, &frac.low);
            frac.high += temp.high + carry;
            base *= step;
            ++p;
        }        
        frac >>= (I - 1);
        low = frac.low;
        high = frac.high + int_val;
        free(str);
    }
    /**
     * @brief Constructor from std::string.
     * Accurate to 37 digits after the decimal point.
     * Allows creating very high precision values. Much slower than the other constructors.
     * @param x Input string
    */
    fixed_point128(const std::string& x) noexcept {
        fixed_point128 temp = x.c_str();
        *this = temp;
    }
    /**
     * @brief Constructor from the 3 fixed_point128 base elements, useful for creating very small constants.
     * @param l Low QWORD
     * @param h High QWORD
     * @param s Sign - zero for positive, 1 for negative.
    */
    fixed_point128(uint64_t l, uint64_t h, uint32_t s) noexcept:
        low(l), high(h) ,sign(s) {
        sign = (sign != 0);
    }
    
    /**
     * @brief Destructor
    */
    ~fixed_point128() noexcept {}
    /**
     * @brief Assignment operator
     * @param other Object to copy from 
     * @return This object.
    */
    FP128_INLINE fixed_point128& operator=(const fixed_point128& other) noexcept {
        high = other.high;
        low = other.low;
        sign = other.sign;
        return *this;
    }
    /**
     * @brief Move assignment operator
     * @param other Object to copy from
     * @return This object.
    */
    FP128_INLINE fixed_point128& operator=(const fixed_point128&& other) noexcept {
        high = other.high;
        low = other.low;
        sign = other.sign;
        return *this;
    }
    /**
     * @brief cross-template assignment operator, can be used between two different fixed_point128 templates
     * @param other fixed_point128 instance with from a different template instance.
     * @return This object.
    */
    template<int32_t I2>
    FP128_INLINE fixed_point128<I>& operator=(const fixed_point128<I2>& other)
    {
        sign = other.sign;
        if constexpr (I == I2) {
            high = other.high;
            low = other.low;
        }
        // other has less integer bits and more fraction bits
        else if constexpr (I < I2) {
            // shift left by I2 - I bits
            const int shift = I2 - I;
            low = other.low << shift;
            high = shift_left128(other.low, other.high, (uint8_t)(64 - shift));
        }
        // other has more integer bits and less fraction bits
        else { // I > I2
            // shift right by I - I2 bits
            const int shift = I - I2;
            low = shift_right128_round(other.low, other.high, (uint8_t)shift);
            high = other.high >> shift;
        }

        return *this;
    }

    //
    // conversion operators
    //
    /**
     * @brief operator uint64_t - converts to a uint64_t
     * @return Object value.
    */
    FP128_INLINE operator uint64_t() const noexcept {
        return (high >> upper_frac_bits) & UINT64_MAX;
    }
    /**
     * @brief operator int64_t - converts to a int64_t
     * @return Object value.
    */
    FP128_INLINE operator int64_t() const noexcept {
        int64_t res = (sign) ? -1ll : 1ll;
        return res * ((high >> upper_frac_bits) & UINT64_MAX);
    }
    /**
     * @brief operator uint32_t - converts to a uint32_t
     * @return Object value.
    */
    FP128_INLINE operator uint32_t() const noexcept {
        return (high >> upper_frac_bits) & UINT32_MAX;
    }
    /**
     * @brief operator int32_t - converts to a int32_t
     * @return Object value.
    */
    FP128_INLINE operator int32_t() const noexcept {
        int32_t res = (sign) ? -1 : 1;
        return res * ((int32_t)((int64_t)high >> upper_frac_bits) & (UINT32_MAX));
    }
    /**
     * @brief operator float - converts to a float
     * @return Object value.
    */
    FP128_INLINE operator float() const noexcept {
        double res = (double)(high * upper_unity); // bits [64:127]
        res += (double)(low * lower_unity);        // bits [0:63]
        return (float)((sign) ? -res : res);
    }
    /**
     * @brief operator double - converts to a double
     * @return Object value.
    */
    FP128_INLINE operator double() const noexcept {
        double res = (double)(high * upper_unity); // bits [64:127]
        res += (double)(low * lower_unity);        // bits [0:63]
        return (sign) ? -res : res;
    }
    /**
     * @brief operator long double - converts to a long double
     * @return Object value.
    */
    FP128_INLINE operator long double() const noexcept {
        long double res = (long double)(high * upper_unity); // bits [64:127]
        res += (long double)(low * lower_unity);             // bits [0:63]
        return (sign) ? -res : res;
    }
    /**
     * @brief Converts to a std::string (slow) string holds all meaningful fraction bits.
     * @return object string representation
    */
    FP128_INLINE operator std::string() const {
        return fp2s();
    }
    /**
     * @brief Converts to a C string (slow) string holds all meaningful fraction bits.
     * @return object string representation
    */
    explicit FP128_INLINE operator char*() const {
        return fp2s();
    }
    //
    // math operators
    //
    /**
     * @brief Adds the right hand side operand to this object to and returns the result.
     * @param other Right hand side operand
     * @return Temporary object with the result of the operation
    */
    FP128_INLINE fixed_point128 operator+(const fixed_point128& other) const {
        fixed_point128 temp(*this);
        return temp += other;
    }
    /**
     * @brief subtracts the right hand side operand to this object to and returns the result.
     * @param other Right hand side operand
     * @return Temporary object with the result of the operation
    */
    FP128_INLINE fixed_point128 operator-(const fixed_point128& other) const {
        fixed_point128 temp(*this);
        return temp -= other;
    }
    /**
     * @brief Multiplies the right hand side operand with this object to and returns the result.
     * @param other Right hand side operand
     * @return Temporary object with the result of the operation
    */
    FP128_INLINE fixed_point128 operator*(const fixed_point128& other) const {
        fixed_point128 temp(*this);
        return temp *= other;
    }
    /**
     * @brief Multiplies the right hand side operand with this object to and returns the result.
     * @param x Right hand side operand
     * @return Temporary object with the result of the operation
    */
    FP128_INLINE fixed_point128 operator*(double x) const {
        fixed_point128 temp(*this);
        return temp *= fixed_point128(x);
    }
    /**
     * @brief Multiplies the right hand side operand with this object to and returns the result.
     * @param x Right hand side operand
     * @return Temporary object with the result of the operation
    */
    FP128_INLINE fixed_point128 operator*(int64_t x) const {
        fixed_point128 temp(*this);
        return temp *= x;
    }
    /**
     * @brief Multiplies the right hand side operand with this object to and returns the result.
     * @param x Right hand side operand
     * @return Temporary object with the result of the operation
    */
    FP128_INLINE fixed_point128 operator*(int32_t x) const {
        fixed_point128 temp(*this);
        return temp *= x;
    }
    /**
     * @brief Divides this object by the right hand side operand and returns the result.
     * @param other Right hand side operand (denominator)
     * @return Temporary object with the result of the operation
    */
    FP128_INLINE fixed_point128 operator/(const fixed_point128& other) const {
        fixed_point128 temp(*this);
        return temp /= other;
    }
    /**
     * @brief Divides this object by the right hand side operand and returns the result.
     * @param x Right hand side operand (denominator)
     * @return Temporary object with the result of the operation
    */
    FP128_INLINE fixed_point128 operator/(double x) const {
        if (x == 0)
            FP128_FLOAT_DIVIDE_BY_ZERO_EXCEPTION;

        fixed_point128 temp(*this);
        temp /= x;
        return temp;
    }
    /**
     * @brief Calculates modulo.
     * @param other Right hand side operand (denominator)
     * @return Temporary object with the result of the operation
    */
    FP128_INLINE fixed_point128 operator%(const fixed_point128& other) const {
        fixed_point128 temp(*this);
        return temp %= other;
    }
    /**
     * @brief Performs right shift operation.
     * @param shift bits to shift
     * @return Temporary object with the result of the operation
    */
    FP128_INLINE fixed_point128 operator>>(int32_t shift) const {
        fixed_point128 temp(*this);
        return temp >>= shift;
    }
    /**
     * @brief Performs left shift operation.
     * @param shift bits to shift
     * @return Temporary object with the result of the operation
    */
    FP128_INLINE fixed_point128 operator<<(int32_t shift) const {
        fixed_point128 temp(*this);
        return temp <<= shift;
    }
    /**
     * @brief Performs bitwise AND (&)
     * @param other Right hand side operand
     * @return Temporary object with the result of the operation
    */
    FP128_INLINE fixed_point128 operator&(const fixed_point128& other) const {
        fixed_point128 temp(*this);
        return temp &= other;
    }
    /**
     * @brief Performs bitwise OR (|)
     * @param other Right hand side operand
     * @return Temporary object with the result of the operation
    */
    FP128_INLINE fixed_point128 operator|(const fixed_point128& other) const {
        fixed_point128 temp(*this);
        return temp |= other;
    }
    /**
     * @brief Performs bitwise XOR (^)
     * @param other Right hand side operand
     * @return Temporary object with the result of the operation
    */
    FP128_INLINE fixed_point128 operator^(const fixed_point128& other) const {
        fixed_point128 temp(*this);
        return temp ^= other;
    }
    /**
     * @brief Add a value to this object
     * @param other Right hand side operand
     * @return This object.
    */
    FP128_INLINE fixed_point128& operator+=(const fixed_point128& other) {
        if (!other) {
            return *this;
        }

        // different sign: convert other to negative and use operator -=
        if (other.sign != sign) {
            fixed_point128 temp = -other;
            return subtract(temp);
        }

        return add(other);
    }
    /**
     * @brief Subtract a value to this object
     * @param other Right hand side operand
     * @return This object.
    */
    FP128_INLINE fixed_point128& operator-=(const fixed_point128& other) {
        if (!other) {
            return *this;
        }

        // different sign: convert other to negative and use operator +=
        if (other.sign != sign) {
            fixed_point128 temp = -other;
            return add(temp);
        }

        return subtract(other);
    }
    /**
     * @brief Multiplies a value to this object
     * @param other Right hand side operand
     * @return This object.
    */
    FP128_INLINE fixed_point128& operator*=(const fixed_point128& other) noexcept{
        // Temporary arrays to store the result. They are uninitialzied to get 10-50% extra performance.
        // Zero initialization is a 10% penalty and using a thread_local static varible lowers
        //  performance by >50%.

        uint64_t res[4]; // 256 bit of result
        uint64_t temp1[2], temp2[2];
        unsigned char carry;

        // multiply low QWORDs
        res[0] = _umul128(low, other.low, &res[1]);

        // multiply high QWORDs (overflow can happen)
        res[2] = _umul128(high, other.high, &res[3]);

        // multiply low this and high other
        temp1[0] = _umul128(low, other.high, &temp1[1]);
        carry = _addcarry_u64(0, res[1], temp1[0], &res[1]);
        res[3] += _addcarry_u64(carry, res[2], temp1[1], &res[2]);

        // multiply high this and low other
        temp2[0] = _umul128(high, other.low, &temp2[1]);
        carry = _addcarry_u64(0, res[1], temp2[0], &res[1]);
        res[3] += _addcarry_u64(carry, res[2], temp2[1], &res[2]);

        // extract the bits from root[] keeping the precision the same as this object
        // shift result by F
        static constexpr int32_t index = F >> 6; // / 64;
        static constexpr int32_t lsb = F & FP128_MAX_VALUE_64(6); // bit within the 64bit data pointed by root[index]
        static constexpr uint64_t half = 1ull << (lsb - 1);       // used for rounding
        const bool need_rounding = (res[index] & half) != 0;

        // copy block #1 (lowest)
        low = shift_right128(res[index], res[index + 1], lsb);
        high = shift_right128(res[index+1], res[index + 2], lsb);

        if (need_rounding) {
            ++low; // low will wrap around to zero if overflowed
            high += low == 0;
        }
        // set the sign
        sign ^= other.sign;
        // set sign to 0 when both low and high are zero (avoid having negative zero value)
        sign &= (0 != low || 0 != high);
        return *this;
    }
    /**
     * @brief Multiplies a value to this object
     * @param x Right hand side operand
     * @return This object.
    */
    FP128_INLINE fixed_point128& operator*=(int32_t x) {
        // alway do positive multiplication
        if (x < 0) {
            x = -x;
            sign ^= 1;
        }
        uint64_t uval = (uint64_t)x;
        uint64_t temp;

        // multiply low QWORDs
        low = _umul128(low, uval, &temp);
        high = high * uval + temp;
        // set sign to 0 when both low and high are zero (avoid having negative zero value)
        sign &= (0 != low || 0 != high);
        return *this;
    }
    /**
     * @brief Multiplies a value to this object
     * @param x Right hand side operand
     * @return This object.
    */
    FP128_INLINE fixed_point128& operator*=(uint32_t x) {
        uint64_t temp;

        // multiply low QWORDs
        low = _umul128(low, x, &temp);
        high = high * x + temp;
        // set sign to 0 when both low and high are zero (avoid having negative zero value)
        sign &= (0 != low || 0 != high);
        return *this;
    }
    /**
     * @brief Multiplies a value to this object
     * @param x Right hand side operand
     * @return This object.
    */
    FP128_INLINE fixed_point128& operator*=(int64_t x) {
        // alway do positive multiplication
        if (x < 0) {
            x = -x;
            sign ^= 1;
        }
        uint64_t uval = (uint64_t)x;
        uint64_t temp;

        // multiply low QWORDs
        low = _umul128(low, uval, &temp);
        high = high * uval + temp;
        // set sign to 0 when both low and high are zero (avoid having negative zero value)
        sign &= (0 != low || 0 != high);
        return *this;
    }
    /**
     * @brief Multiplies a value to this object
     * @param x Right hand side operand
     * @return This object.
    */
    FP128_INLINE fixed_point128& operator*=(uint64_t x) {
        uint64_t temp;

        // multiply low QWORDs
        low = _umul128(low, x, &temp);
        high = high * x + temp;
        // set sign to 0 when both low and high are zero (avoid having negative zero value)
        sign &= (0 != low || 0 != high);
        return *this;
    }
    /**
     * @brief Divide this object by x.
     * @param other Right hand side operator (denominator)
     * @return this object.
    */
    inline fixed_point128& operator/=(const fixed_point128& other) {
        uint64_t q[4]{}; 
        bool need_rounding = false;

        // optimization for when dividing by an integer
        if (other.is_int() && (uint64_t)other <= UINT64_MAX) {
            uint64_t nom[2] = { low, high };
            uint64_t denom = (uint64_t)other;
            uint64_t r;
            if (0 == div_64bit((uint64_t*)q, &r, (uint64_t*)nom, denom, 2)) {
                need_rounding = r > (denom >> 1);
                high = q[1];
                low = q[0];
            }
            else { // error
                FP128_INT_DIVIDE_BY_ZERO_EXCEPTION;
            }
        }
        else {
            uint64_t nom[4] = {0, 0, low, high};
            uint64_t denom[2] = {other.low, other.high};
            if (0 == div_32bit((uint32_t*)q, nullptr, (uint32_t*)nom, (uint32_t*)denom, 2ll * array_length(nom), 2ll * array_length(denom))) {
                static constexpr uint64_t half = 1ull << (I - 1);  // used for rounding
                need_rounding = (q[0] & half) != 0;
                // result in q needs to shifted left by F
                // shifting right by 128-F is simpler.
                high = shift_right128(q[1], q[2], I);
                low = shift_right128(q[0], q[1], I);
            }
            else { // error
                FP128_INT_DIVIDE_BY_ZERO_EXCEPTION;
            }
        }

        if (need_rounding) {
            ++low;
            high += low == 0;
        }
        sign ^= other.sign;
        // set sign to 0 when both low and high are zero (avoid having negative zero value)
        sign &= (0 != low || 0 != high);
        return *this;
    }
    /**
     * @brief Divide this object by x.
     * @param x Denominator.
     * @return This object.
    */
    FP128_INLINE fixed_point128& operator/=(double x) {
        const uint64_t i = *((uint64_t*)(&x));
        // infinity
        if (0 == i) FP128_FLOAT_DIVIDE_BY_ZERO_EXCEPTION;

        uint64_t f = (i & FP128_MAX_VALUE_64(dbl_frac_bits));
        // simple and common case, the value is an exponent of 2
        if (0 == f) {
            sign ^= int32_t(i >> 63);
            int32_t e = FP128_GET_BITS(i, dbl_frac_bits, dbl_exp_bits) - 1023;
            return (e >= 0) ? *this >>= e  : *this <<= e;
        }
        
        // normal division
        *this /= fixed_point128(x);
        return *this;
    }
    /**
     * @brief %= operator
     * @param other Modulo operand.
     * @return This object.
    */
    inline fixed_point128& operator%=(const fixed_point128& other) {
        uint64_t nom[4] = {0, 0, low, high};
        uint64_t denom[2] = {other.low, other.high};
        uint64_t q[4] = {0}, r[4] = {0};
        
        //do the division in with positive numbers
        if (0 == div_32bit((uint32_t*)q, (uint32_t*)r, (uint32_t*)nom, (uint32_t*)denom, 2ll * array_length(nom), 2ll * array_length(denom))) {
            // simple case, both are integers (fractions is zero)
            if (is_int() && other.is_int()) {
                // result is in r (remainder) needs to shift left by F
                low = r[0];
                high = r[1];
            }
            // nom or denom are fractions
            // x mod root =  x - root * floor(x/root)
            else { 
                fixed_point128 x_div_y; // x / root. 
                x_div_y.high = (q[2] << upper_frac_bits) | (q[1] >> I);
                x_div_y.low = (q[1] << upper_frac_bits) | (q[0] >> I);
                x_div_y.sign = sign ^ other.sign;
                *this -= other * floor(x_div_y);
            }
            
            // common case (fractions and integers) where one of the values is negative
            if (sign != other.sign) {
                // the remainder + denominator
                *this += other;
            }

            // Note if signs are the same, for nom/denom, the result keeps the sign.
            // set sign to 0 when both low and high are zero (avoid having negative zero value)
            sign &= (0 != low || 0 != high);
        }
        else { // error
            FP128_INT_DIVIDE_BY_ZERO_EXCEPTION;
        }
        return *this;
    }
    /**
     * @brief Shift right this object.
     * @param shift Bits to shift. Negative or very high values cause undefined behavior.
     * @return This object.
    */
    FP128_INLINE fixed_point128& operator>>=(int32_t shift) noexcept {
        if (shift < 1)
            return *this;
        // 0-64 bit shift - most common
        if (shift <= 64) {
            low = shift_right128_round(low, high, (uint8_t)shift);
            high >>= shift;
        }
        else if (shift >= 128) {
            low = high = 0;
        }
        else if (shift >= 64) {
            low = shift_right64_round(high, shift - 64);
            high = 0;
        }
        // set sign to 0 when both low and high are zero (avoid having negative zero value)
        sign &= (0 != low || 0 != high);
        return *this;
    }
    /**
     * @brief Shift left this object.
     * @param shift Bits to shift. Negative or very high values cause undefined behavior. 
     * @return This object.
    */
    FP128_INLINE fixed_point128& operator<<=(int32_t shift) noexcept {
        if (shift < 1)
            return *this;
        if (shift <= 64) {
            high = shift_left128(low, high, (unsigned char)shift);
            low <<= shift;
        }
        else if (shift < 128) {
            high = low << (shift - 64);
            low = 0;
        }
        else {
            low = high = 0;
        }
        // set sign to 0 when both low and high are zero (avoid having negative zero value)
        sign &= (0 != low || 0 != high);
        return *this;
    }
    /**
     * @brief Bitwise AND=
     * @param other AND mask.
     * @return This object.
    */
    FP128_INLINE fixed_point128& operator&=(const fixed_point128& other) noexcept {
        low &= other.low;
        high &= other.high;
        sign &= other.sign;
        // set sign to 0 when both low and high are zero (avoid having negative zero value)
        sign &= (0 != low || 0 != high);
        return *this;
    }
    /**
     * @brief Bitwise OR=
     * @param other OR mask.
     * @return This object.
    */
    FP128_INLINE fixed_point128& operator|=(const fixed_point128& other) noexcept {
        low |= other.low;
        high |= other.high;
        sign |= other.sign;
        return *this;
    }
    /**
     * @brief Bitwise XOR=
     * @param other XOR mask.
     * @return This object.
    */
    FP128_INLINE fixed_point128& operator^=(const fixed_point128& other) noexcept {
        low ^= other.low;
        high ^= other.high;
        sign ^= other.sign;
        return *this;
    }
    /**
     * @brief Prefix ++ operation (++a)
     * @return This object.
    */
    FP128_INLINE fixed_point128& operator++() noexcept {
        *this += one();
        return *this;
    }
    /**
     * @brief Postfix ++ operation (a++)
     * @return This object.
    */
    FP128_INLINE fixed_point128 operator++(int32_t) noexcept {
        fixed_point128 temp(*this);
        ++*this; // call the prefix implementation
        return temp;
    }
    /**
     * @brief Prefix -- operation (--a)
     * @return This object.
    */
    FP128_INLINE fixed_point128& operator--() noexcept {
        *this -= one();
        return *this;
    }
    /**
     * @brief Postfix -- operation (a--)
     * @return This object.
    */
    FP128_INLINE fixed_point128 operator--(int32_t) noexcept {
        fixed_point128 temp(*this);
        --*this; // call the prefix implementation
        return temp;
    }
    
    //
    // unary operations
    //     
    /**
     * @brief Convert to bool
    */
    FP128_INLINE operator bool() const noexcept {
        return high != 0 || low != 0;
    }
    /**
     * @brief Logical not (!). Opposite of operator bool.
    */
    FP128_INLINE bool operator!() const noexcept {
        return high == 0 && low == 0;
    }
    /**
     * @brief Bitwise not (~).
    */
    FP128_INLINE fixed_point128 operator~() const noexcept {
        fixed_point128 temp(*this);
        temp.high = ~high;
        temp.low = ~low;
        // set sign to 0 when both low and high are zero (avoid having negative zero value)
        temp.sign &= (0 != low || 0 != high);
        return temp;
    }
    /**
     * @brief Unary +. Returns a copy of the object.
    */
    FP128_INLINE fixed_point128 operator+() const noexcept {
        fixed_point128 temp(*this);
        return temp;
    }
    /**
     * @brief Unary -. Returns a copy of the object with sign inverted.
    */
    FP128_INLINE fixed_point128 operator-() const noexcept {
        fixed_point128 temp(*this);

        // set sign to 0 when both low and high are zero (avoid having negative zero value)
        temp.sign ^= 1;
        temp.sign &= (0 != low || 0 != high);
        return temp;
    }

    //
    // Comparison operators
    //
    /**
     * @brief Compare logical/bitwise equal.
     * @param other Righthand operand
     * @return True if this and other are equal.
    */
    FP128_INLINE bool operator==(const fixed_point128& other) const noexcept {
        return sign == other.sign && high == other.high && low == other.low;
    }
    /**
     * @brief Return true when objects are not equal. Can be used as logical XOR.
     * @param other Righthand operand.
     * @return True of not equal.
    */
    FP128_INLINE bool operator!=(const fixed_point128& other) const noexcept {
        return sign != other.sign || high != other.high || low != other.low;
    }
    /**
     * @brief Return true if this object is small than the other
     * @param other Righthand operand.
     * @return True when this object is smaller.
    */
    FP128_INLINE bool operator<(const fixed_point128& other) const noexcept {
        // signs are different
        if (sign != other.sign)
            return sign > other.sign; // true when sign is 1 and other.sign is 0

        // MSB is the same, check the LSB
        if (high == other.high)
            return (sign) ? low > other.low : low < other.low;

        return (sign) ? high > other.high : high < other.high;
    }
    /**
     * @brief Return true this object is small or equal than the other
     * @param other Righthand operand.
     * @return True when this object is smaller or equal.
    */
    FP128_INLINE bool operator<=(const fixed_point128& other) const noexcept {
        return !(*this > other);
    }
    /**
     * @brief Return true this object is larger than the other
     * @param other Righthand operand.
     * @return True when this objext is larger.
    */
    FP128_INLINE bool operator>(const fixed_point128& other) const noexcept {
        // signs are different
        if (sign != other.sign)
            return sign < other.sign; // true when sign is 0 and other.sign is 1

        // MSB is the same, check the LSB
        if (high == other.high)
            return (sign) ? low < other.low : low > other.low;

        return (sign) ? high < other.high : high > other.high;
    }
    /**
     * @brief Return true this object is larger or equal than the other
     * @param other Righthand operand.
     * @return True when this objext is larger or equal.
    */
    FP128_INLINE bool operator>=(const fixed_point128& other) const noexcept {
        return !(*this < other);
    }
    //
    // useful public functions
    //
    /**
     * @brief Returns true if the value is an int (fraction is zero)
     * @return True when the fraction is zero.
    */
    FP128_INLINE bool is_int() const noexcept
    {
        return 0 == low && 0 == (high << I);
    }
    /**
     * @brief Returns true if the value positive (incuding zero)
     * @return True when the the value positive
    */
    FP128_INLINE bool is_positive() const noexcept
    {
        return 0 == sign;
    }
    /**
     * @brief Returns true if the value negative (smaller than zero)
     * @return True when the the value negative
    */
    FP128_INLINE bool is_negative() const noexcept
    {
        return 1 == sign;
    }
    /**
     * @brief Returns true if the value is zero
     * @return Returns true if the value is zero 
    */
    FP128_INLINE bool is_zero() const noexcept
    {
        return 0 == low && 0 == high;
    }
    /**
     * @brief get a specific bit within the 128 fixed point data
     * @param bit bit to get [0,127]
     * @return 0 or 1. Undefined when bit > 127
    */
    FP128_INLINE int32_t get_bit(unsigned bit) const noexcept
    {
        if (bit < 64) {
            return FP128_GET_BIT(low, bit);
        }
        return FP128_GET_BIT(high, bit-64);
    }
    /**
     * @brief Returns the exponent of the object - like the base 2 exponent of a floating point
     * A value of 2.1 would return 1, [0.5,1.0) would return -1.
     * @return Exponent of the number
    */
    FP128_INLINE int32_t get_exponent() const
    {
        int32_t s = (int32_t)lzcnt128(*this);
        s = I - s - 1;
        return s;
    }
    /**
     * @brief Returns an instance of fixed_point128 with the value of pi
     * @return pi
    */
    FP128_INLINE static const fixed_point128& pi() noexcept {
        static const fixed_point128 pi = "3.14159265358979323846264338327950288419716939937510"; // 50 first digits of pi
        return pi;
    }
    /**
     * @brief Returns an instance of fixed_point128 with the value of the golden ratio
     * @return golden ratio constant
    */
    FP128_INLINE static const fixed_point128& golden_ratio() noexcept {
        static const fixed_point128 golden_ratio = "1.6180339887498948482045868343656381177203"; // 40 first digits of the golden ratio
        return golden_ratio;
    }
    /**
     * @brief Returns an instance of fixed_point128 with the value of e
     * @return e
    */
    FP128_INLINE static const fixed_point128& e() noexcept {
        static const fixed_point128 e = "2.71828182845904523536028747135266249775724709369"; // 50 first digits of e
        return e;
    }
    /**
     * @brief Return an instance of fixed_point128 with the value of 1
     * @return 1
    */
    FP128_INLINE static const fixed_point128& one() noexcept {
        static const fixed_point128 one = 1;
        return one;
    }
    /**
     * @brief Return an instance of fixed_point128 with the smallest positive value possible
     * @return 1
    */
    FP128_INLINE static const fixed_point128& epsilon() noexcept {
        static const fixed_point128 epsilon(1, 0, 0);
        return epsilon;
    }

private:
    /**
     * @brief Converts this object to a C string.
     * The returned string is a statically thread-allocated buffer.
     * Additional calls to this function from the same thread, overwrite the previous result.
     * @return C string with describing the value of the object.
    */
    FP128_INLINE char* fp2s() const {
        static thread_local char str[128]; // need roughly a (meaningful) decimal digit per 3.2 bits

        char* p = str;
        fixed_point128 temp = *this;

        //number is negative
        if (temp.sign)
            *p++ = '-';

        uint64_t integer = FP128_GET_BITS(temp.high, upper_frac_bits, 63);
        p += snprintf(p, sizeof(str) + p - str, "%llu", integer);
        temp.high &= ~int_mask; // remove the integer part
        // check if temp has additional digits (not zero)
        if (temp) {
            *p++ = '.';
        }
        // the faster way, requires temp *= 10 not overflowing
        int digits = 0;
        uint64_t res[2]{};
        while (digits++ < max_frac_digits && temp) {
            if constexpr (I < 4) {
                
                res[0] = _umul128(high, 10ull, &res[1]); // multiply by 10
                // extract the integer part
                integer = shift_right128_round(res[0], res[1], upper_frac_bits);
                temp *= 10; // move another digit to the integer area
            }
            else {
                temp *= 10; // move another digit to the integer area
                integer = FP128_GET_BITS(temp.high, upper_frac_bits, 63);
            }
            *p++ = '0' + (char)integer;
            temp.high &= ~int_mask;
        }
        *p = '\0';
        return str;
    }
    /**
     * @brief Adds 2 fixed_point128 objects of the same sign. Throws exception otherwise. this = this + other.
     * @param other The right side of the addition operation
     * @return This object.
    */
    FP128_INLINE fixed_point128& add(const fixed_point128& other) {
        unsigned char carry;
        FP128_ASSERT(other.sign == sign); // bug if asserted, calling method should take care of this
        // equal sign: simple case
        carry = _addcarry_u64(0, low, other.low, &low);
        high += other.high + carry;

        // set sign to 0 when both low and high are zero (avoid having negative zero value)
        sign &= (0 != low || 0 != high);
        return *this;
    }
    /**
     * @brief Subtracts 2 fixed_point128 objects of the same sign. Throws exception otherwise. this = this + other.
     * @param other The right side of the subtraction operation.
     * @return This object.
    */
    FP128_INLINE fixed_point128& subtract(const fixed_point128& other) {
        unsigned char carry;
        FP128_ASSERT(other.sign == sign); // bug if asserted, calling method should take care of this

        // convert other high/low to 2's complement (flip bits, add +1)
        uint64_t other_low = static_cast<uint64_t>(-(int64_t)other.low);
        uint64_t other_high = ~other.high + (other_low == 0);

        //add the other value
        carry = _addcarry_u64(0, low, other_low, &low);
        high += other_high + carry;

        // if result is is negative, invert it along with the sign.
        if (high & FP128_ONE_SHIFT(63)) {
            sign ^= 1;
            low = static_cast<uint64_t>(-(int64_t)low);
            high = ~high + (low == 0);
        }

        // set sign to 0 when both low and high are zero (avoid having negative zero value)
        sign &= (0 != low || 0 != high);
        return *this;
    }
public:
    //
    // Floating point style functions
    //
    /**
     * @brief Returns the absolute value (sets sign to 0)
     * @param x Fixed_point128 object
     * @return A copy of x with sign removed
    */
    friend FP128_INLINE fixed_point128 fabs(const fixed_point128& x) noexcept
    {
        fixed_point128 temp = x;
        temp.sign = 0;
        return temp;
    }
    /**
     * @brief Performs the floor() function, similar to libc's floor(), rounds down towards -infinity.
     * @param x Input value
     * @return A fixed_point128 holding the integer value. Overflow is not reported.
    */
    friend FP128_INLINE fixed_point128 floor(const fixed_point128& x) noexcept
    {
        if (x.is_int()) return x;

        int64_t val = (int64_t)x;
        return (val < 0) ? --val : val;
    }

    /**
     * @brief Performs the ceil() function, similar to libc's ceil(), rounds up towards infinity.
     * @param x Input value
     * @return A fixed_point128 holding the integer value. Overflow is not reported.
    */
    friend FP128_INLINE fixed_point128 ceil(const fixed_point128& x) noexcept
    {
        if (x.is_int()) return x;

        int64_t val = (int64_t)x;
        return (val >= 0) ? ++val : val;
    }

    /**
     * @brief Performs the fmod() function, similar to libc's fmod(), returns the remainder of a division x/root.
     * @param x Numerator
     * @param y Denominator
     * @return A fixed_point128 holding the modulo value.
    */
    friend FP128_INLINE fixed_point128 fmod(const fixed_point128& x, const fixed_point128& y) noexcept
    {
        return x % y;
    }

    /**
     * @brief Split into integer and fraction parts.
     * Both results carry the sign of the input variable.
     * @param x Input value
     * @param iptr Pointer to fixed_point128 holding the integer part of x.
     * @return The fraction part of x. Undefined when iptr is nullptr.
    */
    friend FP128_INLINE fixed_point128 modf(const fixed_point128& x, fixed_point128* iptr) noexcept
    {
        if (iptr == nullptr) {
            return 0;
        }
        iptr->high = x.high & x.int_mask; // lose the fraction
        iptr->low = 0;
        iptr->sign = x.sign;

        fixed_point128 res = x;
        res.high &= ~x.int_mask; // lose the integer part
        return res;
    }
    /**
     * @brief Calculates the left zero count of value x, ignoring the sign.
     * @param x input value.
     * @return lzc (uint32_t) of th result.
    */
    friend FP128_INLINE uint64_t lzcnt128(const fixed_point128& x) noexcept
    {
        return (x.high != 0) ? __lzcnt64(x.high) : 64 + __lzcnt64(x.low);
    }
    /**
     * @brief Calculates the square root using Newton's method.
     * Based on the book "Math toolkit for real time programming" by Jack W. Crenshaw 
     * @param x Value to calculate the root of
     * @param iterations how many iterations to perform (more is more accurate). Sensible values are 0-5.
     * @return Square root of (x), zero when x <= 0.
    */
    friend FP128_INLINE fixed_point128 sqrt(const fixed_point128& x, uint32_t iterations = 3) noexcept
    {
        if (x.sign || !x)
            return 0;
        fixed_point128 root;
        static const fixed_point128 A = "0.41730759963886499890887997563983148154050544649511960252715318";
        static const fixed_point128 B = "0.59016206709064458299663118037100097432703047929336615665205464";
        static const fixed_point128 factor = "0.70710678118654752440084436210484903928483593768847403658833981"; // sqrt(2) / 2
        int32_t expo;
        /*
        * The A and B constants represent a polynom that approximates value of sqrt(x) when 0.5 <= x < 1.0
        * sqrt(x) ~= A + Bx
        * This allows a very good first guess for Newton's method to converge quickly.
        */

        // normalize the input to the range [0.5, 1)
        expo = x.get_exponent() + 1;
        fixed_point128 norm_x = (expo > 0) ? x >> expo : x << -expo;
        
        // get square root of the normalized number
        // generate a first guess
        root = A + B * norm_x;

        // iterate several times via Newton's method
        //
        //                X
        // Xn+1 = 0.5 * (---- + Xn )
        //                Xn
        for (auto i = iterations; i != 0; --i) {
            root = (norm_x / root + root) >> 1;
        }

        if (expo & 1) {
            root *= factor;
            ++expo;
        }

        // Denormalize the result
        if (expo > 0) {
            root <<= (expo + 1) / 2;
        }
        else if (expo < 0)
        {
            root >>= -expo / 2;
        }

        return root;
    }
    /**
     * @brief Calculates the square root using binary search.
     * @param x Value to calculate the root of
     * @return Square root of (x), zero when X <= 0.
    */
    friend FP128_INLINE fixed_point128 sqrt_slow(const fixed_point128& x) noexcept
    {
        if (x.sign || !x)
            return 0;

        fixed_point128 ul = 0, ll = 0, t = 0, e(1, 0, 0);
        int32_t s = 0;
        ul.low = 1ull;
        s = (x.high != 0) ? 128 - (int32_t)__lzcnt64(x.high) : 64 - (int32_t)__lzcnt64(x.low);

        // x >= 1
        if (s >= x.F) {
            ul = x;     // upper limit
            ll.low = 1; // lower limit
            ll <<= x.F + ((s - x.F - 1) >> 1);
        }
        // x < 1
        else {
            ul.low = 1; // upper limit
            ul <<= x.F;
            ll = x;     // lower limit
        }

        // yeh old binary search - need an int256 type to use Newton-Raphson or similar methods
        t = (ul + ll) >> 1;
        while (ul > ll + e) {
            // check if the guess (t) is too big
            if (t * t > x) {
                ul = t; // decrease upper limit
            }
            else {
                ll = t; // increase lower limit
            }
            t = (ul + ll) >> 1;
        }

        return t;
    }
    /**
     * @brief Factorial reciprocal (inverse). Calculates 1 / x!
     * Maximum value of x that may produce non zero values is 34. 
     * This value depends on the amount of fraction bits.
     * @param x Input value
     * @param root Result of the function
    */
    friend FP128_INLINE void fact_reciprocal(int x, fixed_point128& res) noexcept
    {
        static const fixed_point128 c[] = {
            "1",                                         // 1 /  0!
            "1",                                         // 1 /  1!
            "0.5",                                       // 1 /  2!
            "0.166666666666666666666666666666666666667", // 1 /  3!
            "0.041666666666666666666666666666666666667", // 1 /  4!
            "0.008333333333333333333333333333333333333", // 1 /  5!
            "0.001388888888888888888888888888888888889", // 1 /  6!
            "0.000198412698412698412698412698412698413", // 1 /  7!
            "0.000024801587301587301587301587301587302", // 1 /  8!
            "0.000002755731922398589065255731922398589", // 1 /  9!
            "0.000000275573192239858906525573192239859", // 1 / 10!
            "0.000000025052108385441718775052108385442", // 1 / 11!
            "0.000000002087675698786809897921009032120", // 1 / 12!
            "0.000000000160590438368216145993923771702", // 1 / 13!
            "0.000000000011470745597729724713851697979", // 1 / 14!
            "0.000000000000764716373181981647590113199", // 1 / 15!
            "0.000000000000047794773323873852974382075", // 1 / 16!
            "0.000000000000002811457254345520763198946", // 1 / 17!
            "0.000000000000000156192069685862264622164", // 1 / 18!
            "0.000000000000000008220635246624329716956", // 1 / 19!
            "0.000000000000000000411031762331216485848", // 1 / 20!
            "0.000000000000000000019572941063391261231", // 1 / 21!
            "0.000000000000000000000889679139245057329", // 1 / 22!
            "0.000000000000000000000038681701706306840", // 1 / 23!
            "0.000000000000000000000001611737571096118", // 1 / 24!
            "0.000000000000000000000000064469502843845", // 1 / 25!
            "0.000000000000000000000000002479596263225", // 1 / 26!
            "0.000000000000000000000000000091836898638", // 1 / 27!
            "0.000000000000000000000000000003279889237", // 1 / 28!
            "0.000000000000000000000000000000113099628", // 1 / 29!
            "0.000000000000000000000000000000003769988", // 1 / 30!
            "0.000000000000000000000000000000000121613", // 1 / 31!
            "0.000000000000000000000000000000000003800", // 1 / 32!
            "0.000000000000000000000000000000000000115", // 1 / 33!
            "0.000000000000000000000000000000000000003"  // 1 / 34!
        };
        constexpr int series_len = array_length(c);

        if (x >= 0 && x < series_len) {
            res = c[x];
        }
        else {
            res = 0;
        }
    }
    /**
     * @brief Calculate the sine function
     * Using the Maclaurin series expansion, the formula is:
     * sin(x) = x - (x^3 / 3!) + (x^5 / 5!) - (x^7 / 7!) + ...
     * 
     * @param x value in Radians
     * @return Sine of x
    */
    friend FP128_INLINE fixed_point128 sin(const fixed_point128& x) noexcept
    {
        static_assert(I >= 4, "fixed_point128 must have at least 4 integer bits to use sin()!");
        static const fixed_point128 pi = fixed_point128::pi();
        static const fixed_point128 pi2 = pi << 1; // 2 * pi
        static const fixed_point128 half_pi = pi >> 1; // pi / 2

        // first part of the series
        fixed_point128 res = x;
        // move to the range  [-pi, pi]
        if (res > pi || res < -pi) {
            res = fmod(x, pi2);
            if (res > pi)
                res -= pi2;
        }

        // bring closest to zero as possible to minimize the error
        // move to the range [-1/2pi, 1/2pi]
        if (res > half_pi)
            res = pi - res;
        else if (res < -half_pi)
            res = -(pi + res);
        
        const fixed_point128 xx = res * res;
        fixed_point128 elem_denom, elem_nom = res;
        
        // compute the rest of the series, starting with: -(x^3 / 3!)
        for (int i = 3, sign = 1; ; i += 2, sign = 1 - sign) {
            elem_nom *= xx;
            fact_reciprocal(i, elem_denom);

            // precision limit has been hit
            if (!elem_denom)
                break;
            fixed_point128 elem = elem_nom * elem_denom; // next element in the series
            res += (sign) ? -elem : elem;
        }

        return res;
    }
    /**
     * @brief Calculate the cosine function
     * Using the Maclaurin series expansion, the formula is:
     * cos(x) = 1 - (x^2 / 2!) + (x^4 / 4!) - (x^6 / 6!) + ...
     *
     * @param x value in Radians
     * @return Cosine of x
    */
    friend FP128_INLINE fixed_point128 cos(const fixed_point128& x) noexcept
    {
        static_assert(I >= 4, "fixed_point128 must have at least 4 integer bits to use cos()!");
        static const fixed_point128 pi = fixed_point128::pi();
        static const fixed_point128 pi2 = pi << 1; // 2 * pi
        static const fixed_point128 half_pi = pi >> 1; // pi / 2

        // first part of the series
        fixed_point128 res = x;
        // move to the range  [-pi, pi]
        if (res > pi || res < -pi) {
            res = fmod(x, pi2);
            if (res > pi)
                res -= pi2;
        }

        // bring closest to zero as possible to minimize the error
        // move to the range [-1/2pi, 1/2pi]
        if (res > half_pi)
            res = pi - res;
        else if (res < -half_pi)
            res = -(pi + res);

        const fixed_point128 xx = res * res;
        res = fixed_point128::one(); // first element in the series
        fixed_point128 elem_denom, elem_nom = res;

        // compute the rest of the series starting with: -(x^2 / 2!)
        for (int i = 2, sign = 1; ; i += 2, sign = 1 - sign) {
            elem_nom *= xx;
            fact_reciprocal(i, elem_denom);
            
            // precision limit has been hit
            if (!elem_denom)
                break;
            fixed_point128 elem = elem_nom * elem_denom; // next element in the series
            res += (sign) ? -elem : elem;
        }

        return res;
    }
    /**
     * @brief Calculates the exponent of x: e^x
     * Using the Maclaurin series expansion, the formula is:
     * exp(x) = x^0 + (x^1 / 1!) + (x^2 / 2!) + (x^3 / 3!) + ...
     *
     * @param x A number specifying a power. 
     * @return Exponent of x
    */
    friend FP128_INLINE fixed_point128 exp(const fixed_point128& x) noexcept
    {
        static_assert(I >= 4, "fixed_point128 must have at least 4 integer bits to use cos()!");
        static const fixed_point128 e = fixed_point128::e();
        static const fixed_point128 two(2);
        
        // first and second elements of the series
        fixed_point128 res = fixed_point128::one() + x;
        fixed_point128 elem_denom, elem_nom = x;

        for (int i = 2; ; ++i) {
            elem_nom *= x;
            fact_reciprocal(i, elem_denom);
            if (!elem_denom)
                break;

            res += elem_nom * elem_denom; // next element in the series
        }

        return res;
    }
    /**
     * @brief Calculates the Log base 2 of x: log2(x)
     * @param x The number to perform log2 on.
     * @return log2(x)
    */
    friend FP128_INLINE fixed_point128 log2(fixed_point128 x) noexcept
    {
        static const fixed_point128 two(2); 
        fixed_point128 b = fixed_point128::one() >> 1;
        fixed_point128 y = 0;

        if (!x.is_positive() || x.is_zero()) {
            return fixed_point128(UINT64_MAX, UINT64_MAX, 1); // represents negative infinity
        }

        // bring x to the range [1,2)
        while (x < fixed_point128::one()) {
            x <<= 1;
            --y;
        }

        while (x >= two) {
            x >>= 1;
            ++y;
        }

        // x is an exponent of 2.
        if (x == fixed_point128::one())
            return y;

        fixed_point128 z = x;
        for (size_t i = 0; i < fixed_point128::F; ++i) {
            z *= z;
            if (z >= two) {
                z >>= 1;
                y += b;
            }
            b >>= 1;
        }

        return y;
    }
    /**
     * @brief Calculates the natural Log (base e) of x: log(x)
     * @param x The number to perform log on.
     * @return log(x)
    */
    friend FP128_INLINE fixed_point128 log(fixed_point128 x) noexcept
    {
        static const fixed_point128 inv_log2_e = fixed_point128("0.693147180559945309417232121458176575");
        fixed_point128 y = log2(x);

        return y * inv_log2_e;
    }
    /**
     * @brief Calculates Log base 10 of x: log10(x)
     * @param x The number to perform log on.
     * @return log10(x)
    */
    friend FP128_INLINE fixed_point128 log10(fixed_point128 x) noexcept
    {
        static const fixed_point128 inv_log2_10 = fixed_point128("0.301029995663981195213738894724493068");
        fixed_point128 y = log2(x);

        return y * inv_log2_10;
    }
}; //class fixed_point128


} //namespace fp128

#endif // #ifndef FIXED_POINT128_H
