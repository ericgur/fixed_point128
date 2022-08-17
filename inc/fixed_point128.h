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
    The function div_32bit is derived from the book "Hacker's Delight" 2nd Edition by 
    Henry S. Warren Jr. It was converted to 32 bit operations + a bugfix.

************************************************************************************/

// TODO: implement faster _umul128

#ifndef FIXED_POINT128_H
#define FIXED_POINT128_H

#include <intrin.h>
#include <string>
#include <exception>
#include <stdexcept>

typedef __int64 int64;
typedef unsigned __int64 uint64;
typedef int int32;
typedef unsigned int uint32;
typedef short int16;
typedef unsigned short uint16;
typedef char int8;
typedef unsigned char uint8;

// useful macros
#define FP128_ONE_SHIFT(x)          (1ull << (x))
#define FP128_MAX_VALUE_64(x)       (((uint64)-1ll) >> (64 - x))
#define FP128_GET_BIT(x, n)         (((x) >> (n)) & 1)
#define FP128_GET_BITS(x, b, count) (((x) >> (b)) & FP128_MAX_VALUE_64(count))
#define FP128_INT_DIVIDE_BY_ZERO_EXCEPTION   throw std::logic_error("Integer divide by zero!")
#define FP128_FLOAT_DIVIDE_BY_ZERO_EXCEPTION throw std::logic_error("Floating point divide by zero!")
#define FP128_NOT_IMPLEMENTED_EXCEPTION throw std::exception("Not implemented!")
#if defined _DEBUG || defined DEBUG
    #define FP128_ASSERT(x) { if (!(x)) throw std::exception("FP128_ASSERT failed: "#x); } 
#else
    #define FP128_ASSERT(x)
#endif // _DEBUG

// uncomment this to force the functions to not become inline (useful for profiling a specific function)
//#define FP128_DISABLE_INLINE
#ifdef FP128_DISABLE_INLINE
#define FP128_INLINE __declspec(noinline)
#else
#define FP128_INLINE __forceinline
#endif

namespace fp128 {

// Forward declarations
FP128_INLINE int32 div_32bit(uint32* q, uint32* r, const uint32* u, const uint32* v, int64 m, int64 n);

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
template<int32 I>
class fixed_point128
{
    static_assert(1 <= I && I <= 64, "Template parameter <I> must be in the [1,64] range!");
    static_assert(sizeof(void*) == 8, "fixed_point128 is 64 bit only!");

    // friends
    friend class fixed_point128; // this class is a friend of all its template instances. Avoids awkward getter/setter functions.
    //
    // members
    //
    uint64 low; 
    uint64 high;
    unsigned sign; // 0 = positive, 1 negative

    // useful const calculations
    static constexpr int32 F = 128 - I;
    static constexpr int32 upper_frac_bits = F - 64;
    static constexpr uint64 unity = 1ull << (64 - I);
    static inline const double upper_unity = pow(2, 64 - F);
    static inline const double lower_unity = pow(2, -F);
    static constexpr uint32 max_dword_value = (uint32)(-1);
    static constexpr uint64 max_qword_value = (uint64)(-1);
    static constexpr uint64 int_mask = max_qword_value << upper_frac_bits;
    static constexpr int32 dbl_exp_bits = 11;
    static constexpr int32 dbl_frac_bits = 52;
    static constexpr int32 max_frac_digits = (int)(1 + F / 3.1);
public:
    typedef fixed_point128<I> type;
    typedef fixed_point128<I>* ptr_type;
    typedef fixed_point128<I>& ref_type;

    //
    // ctors
    //

    /**
     * @brief Default constructor
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
     * @brief Constructor from double type
     * @param x Input value
    */
    fixed_point128(double x) noexcept {
        // brute convert to uint64 for easy bit manipulation
        uint64 i = *((uint64*)(&x));
        // very common case
        if (i == 0) {
            low = high = 0;
            sign = 0;
            return;
        }

        sign = FP128_GET_BIT(i, 63);
        int32 e = FP128_GET_BITS(i, dbl_frac_bits, dbl_exp_bits) - 1023;
        uint64 f = (i & FP128_MAX_VALUE_64(dbl_frac_bits));
        
        // overflow which catches NaN and Inf
        if (e >= I) {
            high = low = max_qword_value;
            sign = 0;
            return;
        }

        // normal number, produces non zero value
        if (e >= -F) {
            // bit 52 in f is the unity value of the float. it needs to move to the unity position in fixed point
            f |= FP128_ONE_SHIFT(dbl_frac_bits);
            int32 bits_to_shift = 64 - I - dbl_frac_bits + e;

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
     * @brief Constructor from uint64 type
     * @param x Input value
    */
    fixed_point128(uint64 x) noexcept {
        low = 0;
        sign = 0;
        high = x << upper_frac_bits;
    }
    /**
     * @brief Constructor from int64 type
     * @param x Input value
    */
    fixed_point128(int64 x) noexcept {
        low = 0;
        sign = FP128_GET_BIT(x, 63);
        high = ((sign != 0) ? -x : x) << upper_frac_bits;
    }
    /**
     * @brief Constructor from uint32 type
     * @param x Input value
   */
    fixed_point128(uint32 x) noexcept {
        low = 0;
        sign = 0;
        high = (uint64)x << upper_frac_bits;
    }
    /**
     * @brief Constructor from int32 type
     * @param x Input value
   */
    fixed_point128(int32 x) noexcept {
        low = 0;
        sign = FP128_GET_BIT(x, 31);
        high = (uint64)((sign != 0) ? -x : x) << upper_frac_bits;
    }
    /**
     * @brief Constructor from const char* (C string).
     * Allows creating very high precision values. Much slower than the other constructors.
     * @param x Input string
    */
    fixed_point128(const char* x) noexcept {
        low = high = 0;
        sign = 0;
        if (x == nullptr) return;
        
        char* str = _strdup(x);

        char* p = str;
        // set negative sign if needed
        if (p[0] == '-') {
            sign = 1;
            ++p;
        }
        char* dec = strchr(p, '.');
        // number is an integer
        if (dec == nullptr) {
            uint64 int_val = std::strtoull(p, nullptr, 10);
            high = int_val << upper_frac_bits;
            free(str);
            return;
        }
        // number is a float, get the integer part using strtoull()
        *dec = '\0';
        uint64 int_val = std::strtoull(p, nullptr, 10);
        high = int_val << upper_frac_bits;
        p = dec + 1;
        fixed_point128<1> base = 0.1, step = 0.1;
        int32 digits = 0;
        while (digits++  < max_frac_digits && *p != '\0' && base) {
            fixed_point128<1> temp = base * (p[0] - '0');
            // make them the same precision
            temp >>= (I - 1);
            unsigned char carry = _addcarry_u64(0, low, temp.low, &low);
            high += temp.high + carry;
            base *= step;
            ++p;
        }        

        free(str);
    }
    /**
     * @brief Constructor from std::string.
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
    fixed_point128(uint64 l, uint64 h, uint32 s) :
        low(l), high(h) ,sign(s) {
        FP128_ASSERT(sign < 2);
    }

    /**
     * @brief Assignment operator
     * @param other Object to copy from 
     * @return This object.
    */
    FP128_INLINE fixed_point128& operator=(const fixed_point128& other) {
        high = other.high;
        low = other.low;
        sign = other.sign;
        return *this;
    }
    /**
     * @brief template assignment operator, can be used between two different fixed_point128 templates
     * @param other fixed_point128 instance with from a different template instance.
     * @return This object.
    */
    template<int32 I2>
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
            int shift = I2 - I;
            low = other.low << shift;
            high = shift_right128(other.low, other.high, (uint8)(64 - shift));
        }
        // other has more integer bits and less fraction bits
        else { // I > I2
            // shift right by I - I2 bits
            int shift = I - I2;
            low = shift_right128(other.low, other.high, (uint8)shift);
            high = other.high >> shift;
        }

        return *this;
    }

    //
    // conversion operators
    //
    /**
     * @brief operator uint64 - converts to a uint64
     * @return Object value.
    */
    FP128_INLINE operator uint64() const noexcept {
        return (high >> upper_frac_bits) & max_qword_value;
    }
    /**
     * @brief operator int64 - converts to a int64
     * @return Object value.
    */
    FP128_INLINE operator int64() const noexcept {
        int64 res = (sign) ? -1ll : 1ll;
        return res * ((high >> upper_frac_bits) & max_qword_value);
    }
    /**
     * @brief operator uint32 - converts to a uint32
     * @return Object value.
    */
    FP128_INLINE operator uint32() const noexcept {
        return (high >> upper_frac_bits) & max_dword_value;
    }
    /**
     * @brief operator int32 - converts to a int32
     * @return Object value.
    */
    FP128_INLINE operator int32() const noexcept {
        int32 res = (sign) ? -1 : 1;
        return res * ((int32)((int64)high >> upper_frac_bits) & (max_dword_value));
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
     * @brief operator std::string - converts to a std::string (slow)
     *                 string holds all fraction bits.
     * @return object string representation
    */
    FP128_INLINE operator std::string() const {
        char str[128]; // need roughly a (meaningful) digit per 3.2 bits
        
        char* p = str;
        fixed_point128 temp = *this;
        
        //number is negative
        if (temp.sign)
            *p++ = '-';

        uint64 integer = FP128_GET_BITS(temp.high, upper_frac_bits, 63);
        p += snprintf(p, sizeof(str) + p - str, "%lld", integer);
        temp.high &= ~int_mask; // remove the integer part
        // check if temp has additional digits (not zero)
        if (temp) {
            *p++ = '.';
        }
        // the faster way, requires temp *= 10 not overflowing
        int digits = 0;
        while (digits++ < max_frac_digits && temp) {
            if constexpr (I < 4) {
                uint64 res[2];
                res[0] = _umul128(high, 10ull, &res[1]); // multiply by 10
                // extract the integer part
                integer = shift_right128(res[0], res[1], (uint8)upper_frac_bits);
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
    FP128_INLINE fixed_point128 operator*(int64 x) const {
        fixed_point128 temp(*this);
        return temp *= x;
    }
    /**
     * @brief Multiplies the right hand side operand with this object to and returns the result.
     * @param x Right hand side operand
     * @return Temporary object with the result of the operation
    */
    FP128_INLINE fixed_point128 operator*(int32 x) const {
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
     * @param Shift bits to shift
     * @return Temporary object with the result of the operation
    */
    FP128_INLINE fixed_point128 operator>>(int32 shift) const {
        fixed_point128 temp(*this);
        return temp >>= shift;
    }
    /**
     * @brief Performs left shift operation.
     * @param Shift bits to shift
     * @return Temporary object with the result of the operation
    */
    FP128_INLINE fixed_point128 operator<<(int32 shift) const {
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
    FP128_INLINE fixed_point128& operator*=(const fixed_point128& other) {
        uint64 res[4]; // 256 bit of result
        uint64 temp1[2], temp2[2];
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

        // extract the bits from res[] keeping the precision the same as this object
        // shift result by F
        static constexpr int32 index = F >> 6; // / 64;
        static constexpr int32 lsb = F & FP128_MAX_VALUE_64(6);

        // copy block #1 (lowest)
        low = shift_right128(res[index], res[index + 1], lsb);
        high = shift_right128(res[index+1], res[index + 2], lsb);

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
    FP128_INLINE fixed_point128& operator*=(int32 x) {
        // alway do positive multiplication
        if (x < 0) {
            x = -x;
            sign ^= 1;
        }
        uint64 uval = (uint64)x;
        uint64 temp;

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
    FP128_INLINE fixed_point128& operator*=(uint32 x) {
        uint64 temp;

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
    FP128_INLINE fixed_point128& operator*=(int64 x) {
        // alway do positive multiplication
        if (x < 0) {
            x = -x;
            sign ^= 1;
        }
        uint64 uval = (uint64)x;
        uint64 temp;

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
    FP128_INLINE fixed_point128& operator*=(uint64 x) {
        uint64 temp;

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
        uint64 nom[4] = {0, 0, low, high};
        uint64 denom[2] = {other.low, other.high};
        uint64 q[4] = {0}, *r = nullptr; // don't need the reminder
        if (0 == div_32bit((uint32*)q, (uint32*)r, (uint32*)nom, (uint32*)denom, sizeof(nom) / sizeof(uint32), sizeof(denom) / sizeof(uint32))) {
            // result in q needs to shift left by F
            high = shift_right128(q[1], q[2], I);
            low = shift_right128(q[0], q[1], I);
            sign ^= other.sign;
            // set sign to 0 when both low and high are zero (avoid having negative zero value)
            sign &= (0 != low || 0 != high);
        }
        else { // error
            FP128_INT_DIVIDE_BY_ZERO_EXCEPTION;
        }
        return *this;
    }
    /**
     * @brief Divide this object by x.
     * @param x Denominator.
     * @return This object.
    */
    FP128_INLINE fixed_point128& operator/=(double x) {
        uint64 i = *((uint64*)(&x));
        // infinity
        if (0 == i) FP128_FLOAT_DIVIDE_BY_ZERO_EXCEPTION;

        uint64 f = (i & FP128_MAX_VALUE_64(dbl_frac_bits));
        // simple and common case, the value is an exponent of 2
        if (0 == f) {
            sign ^= int32(i >> 63);
            int32 e = FP128_GET_BITS(i, dbl_frac_bits, dbl_exp_bits) - 1023;
            return (e >= 0) ? *this >>= e  : *this <<= e;
        }

        *this /= fixed_point128(x);
        return *this;
    }
    /**
     * @brief %= operator
     * @param other Modulo operand.
     * @return This object.
    */
    inline fixed_point128& operator%=(const fixed_point128& other) {
        uint64 nom[4] = {0, 0, low, high};
        uint64 denom[2] = {other.low, other.high};
        uint64 q[4] = {0}, r[4] = {0};
        
        //do the division in with positive numbers
        if (0 == div_32bit((uint32*)q, (uint32*)r, (uint32*)nom, (uint32*)denom, sizeof(nom) / sizeof(uint32), sizeof(denom) / sizeof(uint32))) {
            // simple case, both are integers (fractions is zero)
            if (is_int() && other.is_int()) {
                // result is in r (remainder) needs to shift left by F
                low = r[0];
                high = r[1];
            }
            // nom or denom are fractions
            // x mod res =  x - res * floor(x/res)
            else { 
                fixed_point128 x_div_y; // x / res. 
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
     * @param shift Bits to shift. negative or very high values cause undefined behavior.
     * @return This object.

    */
    FP128_INLINE fixed_point128& operator>>=(int32 shift) {
        // 0-64 bit shift - most common
        if (shift <= 64) {
            low = shift_right128(low, high, (uint8)shift);
            high >>= shift;
        }
        else if (shift >= 128) {
            low = high = 0;
        }
        else if (shift >= 64) {
            low = high >> (shift - 64);
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
    FP128_INLINE fixed_point128& operator<<=(int32 shift) {
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
    FP128_INLINE fixed_point128& operator&=(const fixed_point128& other) {
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
    FP128_INLINE fixed_point128& operator|=(const fixed_point128& other) {
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
    FP128_INLINE fixed_point128& operator^=(const fixed_point128& other) {
        low ^= other.low;
        high ^= other.high;
        sign ^= other.sign;
    }
    /**
     * @brief Prefix ++ operation (++a)
     * @return This object.
    */
    FP128_INLINE fixed_point128& operator++() {
        high += unity;
        // set sign to 0 when both low and high are zero (avoid having negative zero value)
        sign &= (0 != low || 0 != high);
        return *this;
    }
    /**
     * @brief Postfix ++ operation (a++)
     * @return This object.
    */
    FP128_INLINE fixed_point128 operator++(int32) {
        fixed_point128 temp(*this);
        ++*this; // call the prefix implementation
        return temp;
    }
    /**
     * @brief Prefix -- operation (--a)
     * @return This object.
    */
    FP128_INLINE fixed_point128& operator--() {
        // unity is in the upper QWORD
        high -= unity;
        if (high == max_qword_value) {
            high = 0;
            --low;
        }
        // set sign to 0 when both low and high are zero (avoid having negative zero value)
        sign &= (0 != low || 0 != high);
        return *this;
    }
    /**
     * @brief Postfix -- operation (a--)
     * @return This object.
    */
    FP128_INLINE fixed_point128 operator--(int32) {
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
    FP128_INLINE operator bool() const {
        return high != 0 || low != 0;
    }
    /**
     * @brief Logical not (!). Opposite of operator bool.
    */
    FP128_INLINE bool operator!() const {
        return high == 0 && low == 0;
    }
    /**
     * @brief Bitwise not (~).
    */
    FP128_INLINE fixed_point128 operator~() const {
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
    FP128_INLINE fixed_point128 operator+() const {
        fixed_point128 temp(*this);
        return temp;
    }
    /**
     * @brief Unary -. Returns a copy of the object with sign inverted.
    */
    FP128_INLINE fixed_point128 operator-() const {
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
    FP128_INLINE bool operator==(const fixed_point128& other) const {
        return sign == other.sign && high == other.high && low == other.low;
    }
    /**
     * @brief Return true when objects are not equal. Can be used as logical XOR.
     * @param other Righthand operand.
     * @return True of not equal.
    */
    FP128_INLINE bool operator!=(const fixed_point128& other) const {
        return sign != other.sign || high != other.high || low != other.low;
    }
    /**
     * @brief Return true if this object is small than the other
     * @param other Righthand operand.
     * @return True when this object is smaller.
    */
    FP128_INLINE bool operator<(const fixed_point128& other) const {
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
    FP128_INLINE bool operator<=(const fixed_point128& other) const {
        return !(*this > other);
    }
    /**
     * @brief Return true this object is larger than the other
     * @param other Righthand operand.
     * @return True when this objext is larger.
    */
    FP128_INLINE bool operator>(const fixed_point128& other) const {
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
    FP128_INLINE bool operator>=(const fixed_point128& other) const {
        return !(*this < other);
    }
    //
    // useful public functions
    //
    /**
     * @brief Returns true if the value is an int (fraction is zero)
     * @return True when the fraction is zero.
    */
    FP128_INLINE bool is_int() const
    {
        return 0 == low && 0 == (high << I);
    }
    /**
     * @brief Returns true if the value positive (incuding zero)
     * @return True when the the value positive
    */
    FP128_INLINE bool is_positive() const
    {
        return 0 == sign;
    }
    /**
     * @brief Returns true if the value is zero
     * @return Returns true if the value is zero 
    */
    FP128_INLINE bool is_zero() const
    {
        return 0 == low && 0 == high;
    }
    /**
     * @brief get a specific bit wihtin the 128 fixed point data
     * @param bit bit to get [0,127]
     * @return 0 or 1. Undefined when bit > 127
    */
    FP128_INLINE int32 get_bit(unsigned bit) const
    {
        if (bit <= 64) {
            return FP128_GET_BIT(low, bit);
        }
        return FP128_GET_BIT(high, bit-64);
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
        uint64 other_low = static_cast<uint64>(-(int64)other.low);
        uint64 other_high = ~other.high + (other_low == 0);

        //add the other value
        carry = _addcarry_u64(0, low, other_low, &low);
        high += other_high + carry;

        // if result is is negative, invert it along with the sign.
        if (high & FP128_ONE_SHIFT(63)) {
            sign ^= 1;
            low = static_cast<uint64>(-(int64)low);
            high = ~high + (low == 0);
        }

        // set sign to 0 when both low and high are zero (avoid having negative zero value)
        sign &= (0 != low || 0 != high);
        return *this;
    }

    /**
     * @brief Right shift a 128 bit ineteger.
     * @param l Low QWORD
     * @param h High QWORD
     * @param shift Bits to shift
     * @return Lower 64 bit of the result
    */
    static inline uint64 shift_right128(uint64 l, uint64 h, int shift)
    {
        return (l >> shift) | (h << (64 - shift));
    }
    /**
     * @brief Left shift a 128 bit ineteger.
     * @param l Low QWORD
     * @param h High QWORD
     * @param shift Bits to shift
     * @return Upper 64 bit of the result
    */
    static inline uint64 shift_left128(uint64 l, uint64 h, int shift)
    {
        return (h << shift) | (l >> (64 - shift));
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
        auto temp = x;
        temp.low = 0;
        uint64 frac = temp.high & ~temp.int_mask;
        temp.high &= temp.int_mask;
        // floor always rounds towards -infinity
        if (0 != temp.sign and 0 != frac) {
            ++temp;
        }
        return temp;
    }

    /**
     * @brief Performs the ciel() function, similar to libc's ciel(), rounds up towards infinity.
     * @param x Input value
     * @return A fixed_point128 holding the integer value. Overflow is not reported.
    */
    friend FP128_INLINE fixed_point128 ciel(const fixed_point128& x) noexcept
    {
        auto temp = x;
        temp.low = 0;
        uint64 frac = temp.high & ~temp.int_mask;
        temp.high &= temp.int_mask;
        // ciel always rounds towards infinity
        if (0 != temp.sign and 0 != frac) {
            ++temp;
        }
        return temp;
    }

    /**
     * @brief Performs the fmod() function, similar to libc's fmod(), returns the remainder of a division x/res.
     * @param x Nominator
     * @param res Denominator
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
     * @brief Calculates the square root.
     * @param x Value to calculate the root of
     * @return Square root of (x), zero when X <= 0.
    */
    friend FP128_INLINE fixed_point128 sqrt(const fixed_point128& x) noexcept
    {
        if (x.sign || !x)
            return 0;

        fixed_point128 ul = 0, ll = 0, t = 0, e(1, 0, 0);
        int32 s = 0;
        ul.low = 1ull;
        s = (x.high != 0) ? 128 - (int32)__lzcnt64(x.high) : 64 - (int32)__lzcnt64(x.low);

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
            // printf("g0: %0.15lf\n", (double)ul);
            // printf("g1: %0.15lf\n", (double)ll);
            // printf("t: %0.15lf\n", (double)t);
            
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
     * @brief Calculate the sine function
     * Using the Maclaurin series, the formula is:
     * sin(x) = x - (x^3 / 3!) + (x^5 / 5!) - (x^7 / 7!)...
     * 
     * @param x value in Radians
     * @param precision maximum error bits, default 0 means masimum precision
     * @return Sine of x
    */
    friend FP128_INLINE fixed_point128 sin(const fixed_point128& x, int precision = 0) noexcept
    {
    //#define FAST_SIN
        static_assert(I >= 4, "fixed_point128 must have at least 4 integer bits to use sin()!");
        static const fixed_point128 pi = fixed_point128::pi();
        static const fixed_point128 pi2 = pi << 1; // 2 * pi
        static const fixed_point128 half_pi = pi >> 1; // pi / 2
    #ifdef FAST_SIN
        UNREFERENCED_PARAMETER(precision);

        static const fixed_point128 c[] = {
            fixed_point128(1.0 /   (2 * 3)), // 1 / 3!
            fixed_point128(1.0 /   (4 * 5)), // 1 / 5!
            fixed_point128(1.0 /   (6 * 7)),  // 1 / 7!
            fixed_point128(1.0 /   (8 * 9)),  // 1 / 9!
            fixed_point128(1.0 / (10 * 11)),  // 1 / 11!
            fixed_point128(1.0 / (12 * 13)),  // 1 / 13!
            fixed_point128(1.0 / (14 * 15)),  // 1 / 15!
            fixed_point128(1.0 / (16 * 17)),  // 1 / 17!
            fixed_point128(1.0 / (18 * 19)),  // 1 / 19!
            fixed_point128(1.0 / (20 * 21)),  // 1 / 21!
        };
        constexpr int series_len = sizeof(c) / sizeof(fixed_point128);
    #endif        

        // first part of the series
        fixed_point128 res = fmod(x, pi2);
        if (res > pi)
            res -= pi2;
        // bring closest to zero as possible to minimize the error
        if (res > half_pi)
            res = pi - res;
        else if (res < -half_pi)
            res = -(pi + res);
        
        const fixed_point128 xx = res * res;
        fixed_point128 temp = res;
    #ifdef FAST_SIN
        for (int i = 0; i < series_len; ++i) {
            temp *= c[i] * xx;
            res += (i & 1) ? temp : -temp;
        }
    #else 
        fixed_point128 max_error = (precision == 0) ? fixed_point128::epsilon() : fixed_point128::epsilon() << (precision);
        int64 k = 0;
        for (int i = 0; temp >= max_error; ++i) {
            k += 2;
            double fact = 1.0 / (k * (k + 1));
            temp *= xx * fixed_point128(fact);
            res += (i & 1) ? temp : -temp;
        }

    #endif
        return res;
    }
};
/**
 * @brief 32 bit words unsigned divide function. Variation of the code from the book Hacker's Delight.
 * @param q (output) Pointer to receive the quote
 * @param r (output, optional) Pointer to receive the remainder. Can be nullptr
 * @param u Pointer nominator, an array of uint32
 * @param v Pointer denominator, an array of uint32
 * @param m Count of elements in u
 * @param n Count of elements in v
 * @return 0 for success
*/
int32 div_32bit(uint32* q, uint32* r, const uint32* u, const uint32* v, int64 m, int64 n)
{
    const uint64 b = 1ull << 32; // Number base (32 bits).
    const uint64 mask = b - 1;
    uint32* un, * vn;            // Normalized form of u, v.
    uint64 qhat;                 // Estimated quotient digit.
    uint64 rhat;                 // A remainder.
    uint64 p;                    // Product of two digits.
    int64 s, s_comp, i, j, t, k;
    // shrink the arrays to avoid extra work on small numbers
    while (m >= 0 && u[m - 1] == 0) --m;
    while (n >= 0 && v[n - 1] == 0) --n;

    if (m < n || n <= 0 || v[n - 1] == 0)
        return 1; // Return if invalid param.

    // Take care of the case of a single-digit divisor here.
    if (n == 1) {
        k = 0;
        for (j = m - 1; j >= 0; j--) {
            q[j] = (uint32)((k * b + u[j]) / v[0]);
            k = (k * b + u[j]) - (uint64)q[j] * v[0];
        }

        if (r != nullptr)
            r[0] = (uint32)k;
        return 0;
    }
    // Normalize by shifting v left just enough so that its high-order bit is on, and shift u left the same amount.
    // We may have to append a high-order digit on the dividend; we do that unconditionally.
    s = (uint64)__lzcnt(v[n - 1]); // 0 <= s <= 32. 
    s_comp = 32 - s;
    
    vn = (uint32*)_malloca(sizeof(uint32) * n);
    if (nullptr == vn) return 1;

    for (i = n - 1; i > 0; i--)
        vn[i] = (uint32)(v[i] << s) | (v[i - 1] >> s_comp);
    vn[0] = v[0] << s;
    un = (uint32*)_malloca(sizeof(uint32) * (m + 1));
    if (nullptr == un) return 1;

    un[m] = u[m - 1] >> s_comp;
    for (i = m - 1; i > 0; i--)
        un[i] = (u[i] << s) | (u[i - 1] >> s_comp);

    un[0] = u[0] << s;
    for (j = m - n; j >= 0; j--) { // Main loop. 
                                   // Compute estimate qhat of q[j]. 
        qhat = ((uint64)un[j + n] * b + (uint64)un[j + n - 1]) / (uint64)vn[n - 1];
        rhat = ((uint64)un[j + n] * b + (uint64)un[j + n - 1]) - qhat * (uint64)vn[n - 1];
        while (qhat >= b || qhat * (uint64)vn[n - 2] > b * rhat + (uint64)un[j + n - 2]) {
            qhat = qhat - 1;
            rhat = rhat + (uint64)vn[n - 1];
            if (rhat >= b)
                break;
        }
        // Multiply and subtract. 
        k = 0;
        for (i = 0; i < n; i++) {
            p = qhat * vn[i];
            t = (uint64)un[i + j] - k - (p & mask);
            un[i + j] = (uint32)t;
            k = (p >> 32) - (t >> 32);
        }

        t = (uint64)un[j + n] - k; un[j + n] = (uint32)t;
        q[j] = (uint32)qhat;    // Store quotient digit. 
        if (t < 0) {            // If we subtracted too
            q[j] = q[j] - 1;    // much, add back. 
            k = 0;
            for (i = 0; i < n; i++) {
                t = (uint64)un[i + j] + vn[i] + k;
                un[i + j] = (uint32)t;
                k = t >> 32;
            }
            un[j + n] = (uint32)((uint64)un[j + n] + k);
        }
    } // End j.
    // If the caller wants the remainder, unnormalize it and pass it back. 
    if (r != nullptr) {
        for (i = 0; i < n; i++) {
            r[i] = (un[i] >> s) | (un[i + 1] << s_comp);
        }
    }
    return 0;
}

} //namespace fp128

#endif // #ifndef FIXED_POINT128_H
