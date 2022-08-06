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

    The sin, cos, exp  functions were adapted from the Fixed Point Math 
    Library. Copyright (c) 2007-2009: Peter Schregle

************************************************************************************/

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
#ifndef FP128_ONE_SHIFT
#define FP128_ONE_SHIFT(x)  (1ull << (x))
#endif

#ifndef FP128_MAX_VALUE_64
#define FP128_MAX_VALUE_64(x)   (((uint64)-1ll) >> (64 - x))
#endif

#ifndef FP128_GET_BIT
#define FP128_GET_BIT(x, n)  (((x) >> n) & 1)
#endif 

#ifndef FP128_GET_BITS
#define FP128_GET_BITS(x, b, count)   (((x) >> (b)) & FP128_MAX_VALUE_64(count))
#endif

#ifndef FP128_INT_DIVIDE_BY_ZERO_EXCEPTION
#define FP128_INT_DIVIDE_BY_ZERO_EXCEPTION   throw std::logic_error("Integer divide by zero!")
#endif

#ifndef FP128_FLOAT_DIVIDE_BY_ZERO_EXCEPTION 
#define FP128_FLOAT_DIVIDE_BY_ZERO_EXCEPTION throw std::logic_error("Floating point divide by zero!")
#endif

#ifndef FP128_NOT_IMPLEMENTED_EXCEPTION 
#define FP128_NOT_IMPLEMENTED_EXCEPTION throw std::exception("Not implemented!")
#endif

#ifndef FP128_ASSERT
#ifdef _DEBUG 
#define FP128_ASSERT(x) { if (!(x)) std::exception("FP128_ASSERT failed: "#x); } 
#else
#define FP128_ASSERT(x)
#endif // _DEBUG
#endif // FP128_ASSERT

namespace fp128 {

// Forward declarations
inline int32 div_32bit(uint32* q, uint32* r, const uint32* u, const uint32* v, int64 m, int64 n);

/**
 * @brief 128 bit fixed point type template.
 * The only template parameter I is the bit count of the integer part. The fraction part complements I to 128 bit.
 * I is limited to the range [1,64] in order to simplify the implementation and increase preformance. This restriction is enforced at compile time.
 * All of fixed_point128's methods are inline for maximum performance.
 * Overflow is handled silently, similar to builtin integer operations.
 * Unlike integer operations, the sign is held in a separate data member which simplifies the implementation.
*/
template<int32 I>
class fixed_point128
{
    static_assert(1 <= I && I <= 64, "Template parameter <I> must be in the [1,64] range!");
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
public:
    typedef fixed_point128<I> type;

    //
    // ctors
    //

    /**
     * @brief Default constructor
    */
    fixed_point128() noexcept { 
        low = high = 0ull; sign = 0; 
    }
    /**
     * @brief Copy constructor
     * @param other Object to copy from
    */
    fixed_point128(const fixed_point128& other) noexcept {
        low = other.low;
        high = other.high;
        sign = other.sign;
    }
    /**
     * @brief Constructor from double type
     * @param x input value
    */
    fixed_point128(double x) noexcept {
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

        // normal number, smaller exponents underflow silently to zero
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
     * @param x input value
    */
    fixed_point128(uint64 x) noexcept {
        low = 0;
        sign = 0;
        high = x << upper_frac_bits;
    }
    /**
     * @brief Constructor from int64 type
     * @param x input value
    */
    fixed_point128(int64 x) noexcept {
        low = 0;
        sign = FP128_GET_BIT(x, 63);
        high = ((sign != 0) ? -x : x) << upper_frac_bits;
    }
    /**
     * @brief Constructor from uint32 type
     * @param x input value
   */
    fixed_point128(uint32 x) noexcept {
        low = 0;
        sign = 0;
        high = (uint64)x << upper_frac_bits;
    }
    /**
     * @brief Constructor from int32 type
     * @param x input value
   */
    fixed_point128(int32 x) noexcept {
        low = 0;
        sign = FP128_GET_BIT(x, 31);
        high = (uint64)((sign != 0) ? -x : x) << upper_frac_bits;
    }
    /**
     * @brief Constructor from const char* (C string).
     * Allows creating very high precision values. Much slower than the other constructors.
     * @param x input value
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
            uint64 int_val = atoll(p);
            high = int_val << upper_frac_bits;
            free(str);
            return;
        }
        // number is a float, get the integer part using atoll()
        *dec = '\0';
        uint64 int_val = atoll(p);
        high = int_val << upper_frac_bits;
        p = dec + 1;
        fixed_point128<1> base = 0.1, step = 0.1;
        while (*p != '\0' && base) {
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
     * @brief Constructor from the 3 fixed_point128 base elements, useful for creating very small constants.
     * @param l Low QWORD
     * @param h High QWORD
     * @param s Sign - zero fo rpositive, 1 for negative.
    */
    fixed_point128(uint64 l, uint64 h, uint32 s) {
        FP128_ASSERT(sign < 2);
        low = l;
        high = h;
        sign = s;
    }

    /**
     * @brief assignment operator
     * @param other: object to copy from 
     * @return this
    */
    inline fixed_point128& operator=(const fixed_point128& other) {
        high = other.high;
        low = other.low;
        sign = other.sign;
        return *this;
    }
    //
    // conversion operators
    //
    /**
     * @brief operator uint64 - converts to a uint64
     * @return object value
    */
    inline operator uint64() const noexcept {
        return (high >> upper_frac_bits) & max_qword_value;
    }
    /**
     * @brief operator int64 - converts to a int64
     * @return object value
    */
    inline operator int64() const noexcept {
        int64 res = (sign) ? -1ll : 1ll;
        return res * ((high >> upper_frac_bits) & max_qword_value);
    }
    /**
     * @brief operator uint32 - converts to a uint32
     * @return object value
    */
    inline operator uint32() const noexcept {
        return (high >> upper_frac_bits) & max_dword_value;
    }
    /**
     * @brief operator int32 - converts to a int32
     * @return object value
    */
    inline operator int32() const noexcept {
        int32 res = (sign) ? -1 : 1;
        return res * ((int32)((int64)high >> upper_frac_bits) & (max_dword_value));
    }
    /**
     * @brief operator float - converts to a float
     * @return object value
    */
    inline operator float() const noexcept {
        double res = (double)(high * upper_unity); // bits [64:127]
        res += (double)(low * lower_unity);        // bits [0:63]
        return (float)((sign) ? -res : res);
    }
    /**
     * @brief operator double - converts to a double
     * @return object value
    */
    inline operator double() const noexcept {
        double res = (double)(high * upper_unity); // bits [64:127]
        res += (double)(low * lower_unity);        // bits [0:63]
        return (sign) ? -res : res;
    }
    /**
     * @brief operator long double - converts to a long double
     * @return object value
    */
    inline operator long double() const noexcept {
        double res = (long double)(high * upper_unity); // bits [64:127]
        res += (long double)(low * lower_unity);        // bits [0:63]
        return (sign) ? -res : res;
    }
    /**
     * @brief operator std::string - converts to a std::string (slow)
     *                 string holds all fraction bits.
     * @return object string representation
    */
    inline operator std::string() const {
        char str[128]; // need roughly a digit per 3.5 bits
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
        while (temp) {
            if constexpr (I < 4) {
                uint64 res[2];
                res[0] = _umul128(high, 10ull, &res[1]); // multiply by 10
                // extract the integer part
                integer = __shiftright128(res[0], res[1], (unsigned char)upper_frac_bits);
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
     * @brief adds the right hand side operand to this object to and returns the result.
     * @param other: right hand side operand
     * @return temporary object with the result of the operation
    */
    inline fixed_point128 operator+(const fixed_point128& other) const {
        fixed_point128 temp(*this);
        return temp += other;
    }
    /**
     * @brief subtracts the right hand side operand to this object to and returns the result.
     * @param other: right hand side operand
     * @return temporary object with the result of the operation
    */
    inline fixed_point128 operator-(const fixed_point128& other) const {
        fixed_point128 temp(*this);
        return temp -= other;
    }
    /**
     * @brief multiplies the right hand side operand with this object to and returns the result.
     * @param other: right hand side operand
     * @return temporary object with the result of the operation
    */
    inline fixed_point128 operator*(const fixed_point128& other) const {
        fixed_point128 temp(*this);
        return temp *= other;
    }
    /**
     * @brief multiplies the right hand side operand with this object to and returns the result.
     * @param x: right hand side operand
     * @return temporary object with the result of the operation
    */
    inline fixed_point128 operator*(double x) const {
        fixed_point128 temp(*this);
        return temp *= fixed_point128(x);
    }
    /**
     * @brief multiplies the right hand side operand with this object to and returns the result.
     * @param x: right hand side operand
     * @return temporary object with the result of the operation
    */
    inline fixed_point128 operator*(int64 x) const {
        fixed_point128 temp(*this);
        return temp *= x;
    }
    /**
     * @brief multiplies the right hand side operand with this object to and returns the result.
     * @param x: right hand side operand
     * @return temporary object with the result of the operation
    */
    inline fixed_point128 operator*(int32 x) const {
        fixed_point128 temp(*this);
        return temp *= (int64)x;
    }
    /**
     * @brief divides this object by the right hand side operand and returns the result.
     * @param other: right hand side operand (denominator)
     * @return temporary object with the result of the operation
    */
    inline fixed_point128 operator/(const fixed_point128& other) const {
        fixed_point128 temp(*this);
        return temp /= other;
    }
    /**
     * @brief divides this object by the right hand side operand and returns the result.
     * @param x: right hand side operand (denominator)
     * @return temporary object with the result of the operation
    */
    inline fixed_point128 operator/(double x) const {
        if (x == 0)
            FP128_FLOAT_DIVIDE_BY_ZERO_EXCEPTION;

        fixed_point128 temp(*this);
        temp /= x;
        return temp;
    }
    /**
     * @brief calculates modulo.
     * @param other: right hand side operand (denominator)
     * @return temporary object with the result of the operation
    */
    inline fixed_point128 operator%(const fixed_point128& other) const {
        fixed_point128 temp(*this);
        return temp %= other;
    }
    /**
     * @brief performs right shift operation.
     * @param shift: bits to shift
     * @return temporary object with the result of the operation
    */
    inline fixed_point128 operator>>(int32 shift) const {
        fixed_point128 temp(*this);
        return temp >>= shift;
    }
    /**
     * @brief performs left shift operation.
     * @param shift: bits to shift
     * @return temporary object with the result of the operation
    */
    inline fixed_point128 operator<<(int32 shift) const {
        fixed_point128 temp(*this);
        return temp <<= shift;
    }
    /**
     * @brief performs bitwise AND (&)
     * @param other: right hand side operand
     * @return temporary object with the result of the operation
    */
    inline fixed_point128 operator&(const fixed_point128& other) const {
        fixed_point128 temp(*this);
        return temp &= other;
    }
    /**
     * @brief performs bitwise OR (|)
     * @param other: right hand side operand
     * @return temporary object with the result of the operation
    */
    inline fixed_point128 operator|(const fixed_point128& other) const {
        fixed_point128 temp(*this);
        return temp |= other;
    }
    /**
     * @brief performs bitwise XOR (^)
     * @param other: right hand side operand
     * @return temporary object with the result of the operation
    */
    inline fixed_point128 operator^(const fixed_point128& other) const {
        fixed_point128 temp(*this);
        return temp ^= other;
    }
    /**
     * @brief add a value to this object
     * @param other: right hand side operand
     * @return this
    */
    inline fixed_point128& operator+=(const fixed_point128& other) {
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
     * @brief subtract a value to this object
     * @param other: right hand side operand
     * @return this
    */
    inline fixed_point128& operator-=(const fixed_point128& other) {
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
     * @brief multiplies a value to this object
     * @param other: right hand side operand
     * @return this
    */
    inline fixed_point128& operator*=(const fixed_point128& other) {
    //__declspec(noinline) fixed_point128& operator*=(const fixed_point128& other) { // use this line instead when profiling this function
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
        low = __shiftright128(res[index], res[index + 1], lsb);
        high = __shiftright128(res[index+1], res[index + 2], lsb);

        // set the sign
        sign ^= other.sign;
        // set sign to 0 when both low and high are zero (avoid having negative zero value)
        sign &= (0 != low || 0 != high);
        return *this;
    }
    /**
     * @brief multiplies a value to this object
     * @param x: right hand side operand
     * @return this
    */
    inline fixed_point128& operator*=(int32 x) {
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
     * @brief multiplies a value to this object
     * @param x: right hand side operand
     * @return this
    */
    inline fixed_point128& operator*=(uint32 x) {
        uint64 temp;

        // multiply low QWORDs
        low = _umul128(low, x, &temp);
        high = high * x + temp;
        // set sign to 0 when both low and high are zero (avoid having negative zero value)
        sign &= (0 != low || 0 != high);
        return *this;
    }
    /**
     * @brief multiplies a value to this object
     * @param x: right hand side operand
     * @return this
    */
    inline fixed_point128& operator*=(int64 x) {
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
     * @brief multiplies a value to this object
     * @param x: right hand side operand
     * @return this
    */
    inline fixed_point128& operator*=(uint64 x) {
        uint64 temp;

        // multiply low QWORDs
        low = _umul128(low, x, &temp);
        high = high * x + temp;
        // set sign to 0 when both low and high are zero (avoid having negative zero value)
        sign &= (0 != low || 0 != high);
        return *this;
    }
    /**
     * @brief divide this object by x.
     * @param other: right hand side operator (denominator)
     * @return this object.
    */
    inline fixed_point128& operator/=(const fixed_point128& other) {
        uint64 nom[4] = {0, 0, low, high};
        uint64 denom[2] = {other.low, other.high};
        uint64 q[4] = {0}, *r = nullptr; // don't need the reminder
        if (0 == div_32bit((uint32*)q, (uint32*)r, (uint32*)nom, (uint32*)denom, sizeof(nom) / sizeof(uint32), sizeof(denom) / sizeof(uint32))) {
            // result in q needs to shift left by F
            high = __shiftright128(q[1], q[2], I);
            low = __shiftright128(q[0], q[1], I);
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
     * @brief divide this object by x.
     * @param x: denominator.
     * @return this object.
    */
    inline fixed_point128& operator/=(double x) {
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
     * @param other: modulo operand.
     * @return this
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
            // x mod y =  x - y * floor(x/y)
            else { 
                fixed_point128 x_div_y; // x / y. 
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
     * @brief left right for this object.
     * @param shift: how many bits to shift. neagive or very high values cause undefined behavior.
     * @return this
    */
    inline fixed_point128& operator>>=(int32 shift) {
        // 0-64 bit shift - most common
        if (shift <= 64) {
            low = __shiftright128(low, high, (unsigned char)shift);
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
     * @brief left shift for this object.
     * @param shift: how many bits to shift. neagive or very high values cause undefined behavior. 
     * @return this
    */
    inline fixed_point128& operator<<=(int32 shift) {
        if (shift <= 64) {
            high = __shiftleft128(low, high, (unsigned char)shift);
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
     * @brief bitwise AND=
     * @param other mask
     * @return this
    */
    inline fixed_point128& operator&=(const fixed_point128& other) {
        low &= other.low;
        high &= other.high;
        sign &= other.sign;
        // set sign to 0 when both low and high are zero (avoid having negative zero value)
        sign &= (0 != low || 0 != high);
        return *this;
    }
    /**
     * @brief bitwise OR=
     * @param other mask
     * @return this
    */
    inline fixed_point128& operator|=(const fixed_point128& other) {
        low |= other.low;
        high |= other.high;
        sign |= other.sign;
        return *this;
    }
    /**
     * @brief bitwise XOR=
     * @param other mask
     * @return this
    */
    inline fixed_point128& operator^=(const fixed_point128& other) {
        low ^= other.low;
        high ^= other.high;
        sign ^= other.sign;
    }
    /**
     * @brief prefix ++ operation (++a)
     * @return this
    */
    inline fixed_point128& operator++() {
        high += unity;
        // set sign to 0 when both low and high are zero (avoid having negative zero value)
        sign &= (0 != low || 0 != high);
        return *this;
    }
    /**
     * @brief postfix ++ operation (a++)
     * @return this
    */
    inline fixed_point128 operator++(int32) {
        fixed_point128 temp(*this);
        ++*this; // call the prefix implementation
        return temp;
    }
    /**
     * @brief prefix -- operation (--a)
     * @return this
    */
    inline fixed_point128& operator--() {
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
     * @brief postfix -- operation (a--)
     * @return this
    */
    inline fixed_point128 operator--(int32) {
        fixed_point128 temp(*this);
        --*this; // call the prefix implementation
        return temp;
    }
    //
    // unary operations
    // 
    
    /**
     * @brief convert to bool
    */
    inline operator bool() const {
        return high != 0 || low != 0;
    }
    /**
     * @brief logical not (!). Opposite of operator bool.
    */
    inline bool operator!() const {
        return high == 0 && low == 0;
    }
    /**
     * @brief bitwise not (~).
    */
    inline fixed_point128 operator~() const {
        fixed_point128 temp(*this);
        temp.high = ~high;
        temp.low = ~low;
        // set sign to 0 when both low and high are zero (avoid having negative zero value)
        temp.sign &= (0 != low || 0 != high);
        return temp;
    }
    /**
     * @brief unary +. returns a copy of the object.
    */
    inline fixed_point128 operator+() const {
        fixed_point128 temp(*this);
        return temp;
    }
    /**
     * @brief unary -. returns a copy of the object with sign inverted.
    */
    inline fixed_point128 operator-() const {
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
     * @brief compare logical/bitwise equal.
     * @param other: righthand operand
     * @return true of equal
    */
    inline bool operator==(const fixed_point128& other) const {
        return sign == other.sign && high == other.high && low == other.low;
    }
    /**
     * @brief return true when objects are not equal. Can be used as logical XOR.
     * @param other: righthand operand
     * @return true of not equal
    */
    inline bool operator!=(const fixed_point128& other) const {
        return sign != other.sign || high != other.high || low != other.low;
    }
    /**
     * @brief return true this object is small than the other
     * @param other: righthand operand
     * @return true when this is smaller
    */
    inline bool operator<(const fixed_point128& other) const {
        // signs are different
        if (sign != other.sign)
            return sign > other.sign; // true when sign is 1 and other.sign is 0

        // MSB is the same, check the LSB
        if (high == other.high)
            return (sign) ? low > other.low : low < other.low;

        return (sign) ? high > other.high : high < other.high;
    }
    /**
     * @brief return true this object is small or equal than the other
     * @param other: righthand operand
     * @return true when this is smaller or equal
    */
    inline bool operator<=(const fixed_point128& other) const {
        return !(*this > other);
    }
    /**
     * @brief return true this object is larger than the other
     * @param other: righthand operand
     * @return true when this is larger
    */
    inline bool operator>(const fixed_point128& other) const {
        // signs are different
        if (sign != other.sign)
            return sign < other.sign; // true when sign is 0 and other.sign is 1

        // MSB is the same, check the LSB
        if (high == other.high)
            return (sign) ? low < other.low : low > other.low;

        return (sign) ? high < other.high : high > other.high;
    }
    /**
     * @brief return true this object is larger or equal than the other
     * @param other: righthand operand
     * @return true when this is larger or equal
    */
    inline bool operator>=(const fixed_point128& other) const {
        return !(*this < other);
    }
    //
    // useful public functions
    //

    /**
     * @brief returns true if the value is an int (fraction is zero)
     * @return true when the fraction is zero
    */
    inline bool is_int() const
    {
        return 0 == low && 0 == (high << I);
    }
    /**
     * @brief returns true if the value positive (incuding zero)
     * @return true when the the value positive
    */
    inline bool is_positive() const
    {
        return 0 == sign;
    }
    /**
     * @brief returns true if the value is zero
     * @return returns true if the value is zero 
    */
    inline bool is_zero() const
    {
        return 0 == low && 0 == high;
    }
    /**
     * @brief get a specific bit wihtin the 128 fixed point data
     * @param bit bit to get [0,127]
     * @return 0 or 1. undefined when bit is > 127
    */
    inline int32 get_bit(unsigned bit) const
    {
        if (bit <= 64) {
            return FP128_GET_BIT(low, bit);
        }
        return FP128_GET_BIT(high, bit-64);
    }
    /**
     * @brief Return an instance of fixed_point128 with the value of pi
     * @return pi
    */
    inline static const fixed_point128& pi() noexcept {
        static const fixed_point128 pi = "3.14159265358979323846264338327950288419716939937510"; // 50 first digits of pi
        return pi;
    }
    /**
     * @brief Return an instance of fixed_point128 with the value of e
     * @return e
    */
    inline static const fixed_point128& e() noexcept {
        static const fixed_point128 e = "2.71828182845904523536028747135266249775724709369"; // 50 first digits of e
        return e;
    }
    /**
     * @brief Return an instance of fixed_point128 with the value of 1
     * @return 1
    */
    inline static const fixed_point128& one() noexcept {
        static const fixed_point128 one = 1;
        return one;
    }

private:
    /**
     * @brief Adds 2 fixed_point128 objects of the same sign. Throws exception otherwise. this = this + other.
     * @param other The right side of the addition operation
     * @return This object
    */
    inline fixed_point128& add(const fixed_point128& other) {
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
     * @return This object
    */
    inline fixed_point128& subtract(const fixed_point128& other) {
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

public:
    /**
     * @brief returns the absolute value (sets sign to 0)
     * @param x - fixed_point128 element
     * @return - a copy of x with sign removed
    */
    friend inline fixed_point128 fabs(const fixed_point128& x) noexcept
    {
        fixed_point128 temp = x;
        temp.sign = 0;
        return temp;
    }
    /**
     * @brief peforms the floor() function, similar to libc's floor(), rounds down towards -infinity.
     * @param x - input value
     * @return a fixed_point128 holding the integer value. Overflow is not reported.
    */
    friend inline fixed_point128 floor(const fixed_point128& x) noexcept
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
     * @brief peforms the ciel() function, similar to libc's ciel(), rounds up towards infinity.
     * @param x - input value
     * @return a fixed_point128 holding the integer value. Overflow is not reported.
    */
    friend inline fixed_point128 ciel(const fixed_point128& x) noexcept
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
     * @brief peforms the fmod() function, similar to libc's fmod(), returns the remainder of a division x/y.
     * @param x - nominator
     * @param y - denominator
     * @return a fixed_point128 holding the modulo value.
    */
    friend inline fixed_point128 fmod(const fixed_point128& x, const fixed_point128& y) noexcept
    {
        return x % y;
    }

    /**
     * @brief Split into integer and fraction parts.
     * @param x - input value
     * @param iptr - pointer to fixed_point128 holding the integer part of x
     * @return the fraction part of x. Undefined when iptr is nullptr.
    */
    friend inline fixed_point128 modf(const fixed_point128& x, fixed_point128* iptr) noexcept
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
     * @param x - value to calculate the root of
     * @return square root of (x), zero for negative or zero values of x
    */
    friend inline fixed_point128 sqrt(const fixed_point128& x) noexcept
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
};
/**
 * @brief 32 bit words unsigned divide function. Variation of the code from the book Hacker's Delight.
 * @param q - (output) Pointer to receive the quote
 * @param r - (output, optional) Pointer to receive the remainder. Can be nullptr
 * @param u - Pointer nominator, an array of uint32
 * @param v - Pointer denominator, an array of uint32
 * @param m - count of elements in u
 * @param n - count of elements in v
 * @return 0 for success
*/
inline int32 div_32bit(uint32* q, uint32* r, const uint32* u, const uint32* v, int64 m, int64 n)
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
