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

#ifndef ONE_SHIFT
#define ONE_SHIFT(x)  (1ull << (x))
#endif

#ifndef MAX_BITS_VALUE_64
#define MAX_BITS_VALUE_64(x)   (((uint64)-1ll) >> (64 - x))
#endif

#ifndef GET_BIT
#define GET_BIT(x, n)  (((x) >> n) & 1)
#endif

#ifndef GET_BITS
#define GET_BITS(x, b, count)   (((x) >> (b)) & MAX_BITS_VALUE_64(count))
#endif

#ifndef INT_DIVIDE_BY_ZERO_EXCEPTION
#define INT_DIVIDE_BY_ZERO_EXCEPTION   throw std::logic_error("Integer divide by zero!")
#endif

#ifndef FLOAT_DIVIDE_BY_ZERO_EXCEPTION 
#define FLOAT_DIVIDE_BY_ZERO_EXCEPTION throw std::logic_error("Floating point divide by zero!")
#endif

#ifndef NOT_IMPLEMENTED_EXCEPTION 
#define NOT_IMPLEMENTED_EXCEPTION throw std::exception("Not implemented!")
#endif
#ifndef ASSERT
#ifdef _DEBUG 
#define ASSERT(x) { if (!(x)) std::exception("ASSERT failed: "#x); } 
#else
#define ASSERT(x)
#endif // _DEBUG
#endif // ASSERT

namespace fp128 {

// Forward declarations
inline int div_32bit(uint32* q, uint32* r, const uint32* u, const uint32* v, int64 m, int64 n);
template<int int_bits> class fixed_point128;
template<int int_bits> inline fixed_point128<int_bits> abs(const fixed_point128<int_bits>& val) noexcept;
template<int int_bits> inline fixed_point128<int_bits> floor(const fixed_point128<int_bits>& val) noexcept;
template<int int_bits> inline fixed_point128<int_bits> ciel(const fixed_point128<int_bits>& val) noexcept;
// Main fixed point type template
template<int int_bits = 16>
class fixed_point128
{
    static_assert(1 <= int_bits && int_bits <= 64, "Template parameter <int_bits> must be in the [1,64] range!");

private:
    // members
    uint64 low;
    uint64 high;
    unsigned sign; // 0 = positive, 1 negative

    // useful const calculations
    static constexpr int frac_bits = 128 - int_bits;
    static constexpr int upper_frac_bits = frac_bits - 64;
    static constexpr uint64 unity = 1ull << (64 - int_bits);
    static inline const double upper_unity = pow(2, -upper_frac_bits);
    static inline const double lower_unity = pow(2, -frac_bits);
    static constexpr unsigned max_dword_value = (unsigned)(-1);
    static constexpr uint64 max_qword_value = (uint64)(-1);
    static constexpr uint64 int_mask = max_qword_value << upper_frac_bits;
    static constexpr int dbl_exp_bits = 11;
    static constexpr int dbl_frac_bits = 52;
public:
    // ctors
    fixed_point128() noexcept { 
        low = high = 0ull; sign = 0; 
    }

    fixed_point128(const fixed_point128& other) noexcept {
        low = other.low;
        high = other.high;
        sign = other.sign;
    }

    fixed_point128(double val) noexcept {
        uint64 i = *((uint64*)(&val));
        // very common case
        if (i == 0) {
            low = high = 0;
            sign = 0;
            return;
        }

        sign = GET_BIT(i, 63);
        int e = GET_BITS(i, dbl_frac_bits, dbl_exp_bits) - 1023;
        uint64 f = (i & MAX_BITS_VALUE_64(dbl_frac_bits));
        
        // overflow which catches NaN and Inf
        if (e >= int_bits) {
            high = low = max_qword_value;
            sign = 0;
            return;
        }

        // normal number, smaller exponents underflow silently to zero
        if (e >= -frac_bits) {
            // bit 52 in f is the unity value of the float. it needs to move to the unity position in fixed point
            f |= ONE_SHIFT(dbl_frac_bits);
            int bits_to_shift = 64 - int_bits - dbl_frac_bits + e;

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
                    low = GET_BITS(f, 0, bits_to_shift);
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

    fixed_point128(uint64 val) noexcept {
        low = 0;
        sign = 0;
        high = val << upper_frac_bits;
    }

    fixed_point128(int64 val) noexcept {
        low = 0;
        sign = GET_BIT(val, 63);
        high = ((sign != 0) ? -val : val) << upper_frac_bits;
    }

    fixed_point128(unsigned val) noexcept {
        low = 0;
        sign = 0;
        high = (uint64)val << upper_frac_bits;
    }

    fixed_point128(int val) noexcept {
        low = 0;
        sign = GET_BIT(val, 31);
        high = (uint64)((sign != 0) ? -val : val) << upper_frac_bits;
    }

    // assignment operators
    inline fixed_point128& operator=(const fixed_point128& other) {
        high = other.high;
        low = other.low;
        sign = other.sign;
        return *this;
    }

    // conversion operators
    inline operator uint64() const noexcept {
        return (high >> upper_frac_bits) & max_qword_value;
    }

    inline operator int64() const noexcept {
        int64 res = (sign) ? -1ll : 1ll;
        return res * ((high >> upper_frac_bits) & max_qword_value);
    }

    inline operator uint32() const noexcept {
        return (high >> upper_frac_bits) & max_dword_value;
    }

    inline operator int32() const noexcept {
        int res = (sign) ? -1 : 1;
        return res * ((int)((int64)high >> upper_frac_bits) & (max_dword_value));
    }

    inline operator double() const noexcept {
        double res = (double)(high * upper_unity); // bits [64:127]
        res += (double)(low * lower_unity);        // bits [0:63]
        return (sign) ? -res : res;
    }

    inline operator std::string() const {
        char str[64]; // need roughly a digit per 3.5 bits
        char* p = str;
        fixed_point128 temp = *this;
        
        //number is negative
        if (temp.sign)
            *p++ = '-';

        uint64 integer = GET_BITS(temp.high, upper_frac_bits, 63);
        p += snprintf(p, sizeof(str) + p - str, "%lld", integer);
        temp.high &= ~int_mask; // remove the integer part
        // check if temp has additional digits (not zero)
        if (temp) {
            *p++ = '.';
        }
        while (temp) {
            temp *= 10; // move another digit to the integer area
            integer = GET_BITS(temp.high, upper_frac_bits, 63);
            *p++ = '0' + (char)integer;
            temp.high &= ~int_mask;
        }
        *p = '\0';
        return str;
    }

    // math operators
    inline fixed_point128 operator+(const fixed_point128& other) const {
        fixed_point128 temp(*this);
        return temp += other;
    }

    inline fixed_point128 operator-(const fixed_point128& other) const {
        fixed_point128 temp(*this);
        return temp -= other;
    }

    inline fixed_point128 operator*(const fixed_point128& other) const {
        fixed_point128 temp(*this);
        return temp *= other;
    }

    inline fixed_point128 operator*(double val) const {
        fixed_point128 temp(*this);
        return temp *= fixed_point128(val);
    }

    inline fixed_point128 operator*(int64 val) const {
        fixed_point128 temp(*this);
        return temp *= val;
    }

    inline fixed_point128 operator*(int val) const {
        fixed_point128 temp(*this);
        return temp *= (int64)val;
    }

    inline fixed_point128 operator/(const fixed_point128& other) const {
        fixed_point128 temp(*this);
        return temp /= other;
    }

    inline fixed_point128 operator/(double val) const {
        if (val == 0)
            FLOAT_DIVIDE_BY_ZERO_EXCEPTION;

        fixed_point128 temp(*this);
        temp /= val;
        return temp;
    }

    inline fixed_point128 operator%(const fixed_point128& other) const {
        fixed_point128 temp(*this);
        return temp %= other;
    }

    inline fixed_point128 operator>>(int shift) const {
        fixed_point128 temp(*this);
        return temp >>= shift;
    }

    inline fixed_point128 operator<<(int shift) const {
        fixed_point128 temp(*this);
        return temp <<= shift;
    }

    inline fixed_point128 operator&(const fixed_point128& other) const {
        fixed_point128 temp(*this);
        return temp &= other;
    }

    inline fixed_point128 operator|(const fixed_point128& other) const {
        fixed_point128 temp(*this);
        return temp |= other;
    }

    inline fixed_point128 operator^(const fixed_point128& other) const {
        fixed_point128 temp(*this);
        return temp ^= other;
    }

    inline fixed_point128& operator+=(const fixed_point128& other) {
        unsigned char carry;
        // different sign: convert other to negative and use operator -=
        if (other.sign != sign) {
            fixed_point128 temp(other);
            temp.sign ^= 1;
            *this -= temp;
            return *this;
        }

        // equal sign: simple case
        carry = _addcarry_u64(0, low, other.low, &low);
        _addcarry_u64(carry, high, other.high, &high);
        
        // set sign to 0 when both low and high are zero (avoid having negative zero value)
        sign &= (0 != low || 0 != high);
        return *this;
    }

    inline fixed_point128& operator-=(const fixed_point128& other) {
        unsigned char carry;
        // different sign: convert other to negative and use operator +=
        if (other.sign != sign) {
            fixed_point128 temp(other);
            temp.sign ^= 1;
            *this += temp;
            return *this;
        }
        
        // equal sign: simple case
        // convert other high/low to 2's complement (flip bits, add +1)
        uint64 other_low, other_high = ~other.high;
        other_high += _addcarry_u64(0, ~other.low, 1, &other_low);
        
        //add the other value
        carry = _addcarry_u64(0, low, other_low, &low); // convert low to 2's complement
        _addcarry_u64(carry, high, other_high, &high); // convert low to 2's complement

        // if result is is negative, invert it along with the sign.
        if (0 != (high >> 63)) {
            sign ^= 1;
            carry = _addcarry_u64(0, ~low, 1, &low);
            high = ~high + carry;
        }
        
        // set sign to 0 when both low and high are zero (avoid having negative zero value)
        sign &= (0 != low || 0 != high);
        return *this;
    }

    inline fixed_point128& operator*=(const fixed_point128& other) {
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
        // shift result by frac_bits
        static constexpr int index = frac_bits >> 6; // / 64;
        static constexpr int lsb = frac_bits & MAX_BITS_VALUE_64(6);
        static constexpr int lsb_comp = 64 - lsb;
        static constexpr uint64 lsb_comp_mask = MAX_BITS_VALUE_64(lsb_comp);
        static_assert(lsb <= 64);

        // copy block #1 (lowest)
        low = res[index + 1] << lsb_comp;
        low |= (res[index] >> lsb) & lsb_comp_mask;
        high = res[index + 2] << lsb_comp;
        high |= (res[index + 1] >> lsb) & lsb_comp_mask;

        // set the sign
        sign ^= other.sign;
        // set sign to 0 when both low and high are zero (avoid having negative zero value)
        sign &= (0 != low || 0 != high);
        return *this;
    }

    inline fixed_point128& operator*=(int32 val) {
        // alway do positive multiplication
        if (val < 0) {
            val = -val;
            sign ^= 1;
        }
        uint64 uval = (uint64)val;
        uint64 temp;

        // multiply low QWORDs
        low = _umul128(low, uval, &temp);
        high = high * uval + temp;
        // set sign to 0 when both low and high are zero (avoid having negative zero value)
        sign &= (0 != low || 0 != high);
        return *this;
    }

    inline fixed_point128& operator*=(uint32 val) {
        uint64 temp;

        // multiply low QWORDs
        low = _umul128(low, val, &temp);
        high = high * val + temp;
        // set sign to 0 when both low and high are zero (avoid having negative zero value)
        sign &= (0 != low || 0 != high);
        return *this;
    }

    inline fixed_point128& operator*=(int64 val) {
        // alway do positive multiplication
        if (val < 0) {
            val = -val;
            sign ^= 1;
        }
        uint64 uval = (uint64)val;
        uint64 temp;

        // multiply low QWORDs
        low = _umul128(low, uval, &temp);
        high = high * uval + temp;
        // set sign to 0 when both low and high are zero (avoid having negative zero value)
        sign &= (0 != low || 0 != high);
        return *this;
    }

    inline fixed_point128& operator*=(uint64 val) {
        uint64 temp;

        // multiply low QWORDs
        low = _umul128(low, val, &temp);
        high = high * val + temp;
        // set sign to 0 when both low and high are zero (avoid having negative zero value)
        sign &= (0 != low || 0 != high);
        return *this;
    }

    inline fixed_point128& operator/=(const fixed_point128& other) {
        uint64 nom[4] = {0, 0, low, high};
        uint64 denom[2] = {other.low, other.high};
        uint64 q[4] = {0}, *r = nullptr; // don't need the reminder
        if (0 == div_32bit((uint32*)q, (uint32*)r, (uint32*)nom, (uint32*)denom, sizeof(nom) / sizeof(uint32), sizeof(denom) / sizeof(uint32))) {
            // result in q needs to shift left by frac_bits
            high = (q[2] << upper_frac_bits) | (q[1] >> int_bits);
            low  = (q[1] << upper_frac_bits) | (q[0] >> int_bits);
            sign ^= other.sign;
            // set sign to 0 when both low and high are zero (avoid having negative zero value)
            sign &= (0 != low || 0 != high);
        }
        else { // error
            INT_DIVIDE_BY_ZERO_EXCEPTION;
        }
        return *this;
    }

    inline fixed_point128& operator/=(double val) {
        uint64 i = *((uint64*)(&val));
        // infinity
        if (0 == i) FLOAT_DIVIDE_BY_ZERO_EXCEPTION;

        uint64 f = (i & MAX_BITS_VALUE_64(dbl_frac_bits));
        // simple and common case, the value is an exponent of 2
        if (0 == f) {
            sign ^= int(i >> 63);
            int e = GET_BITS(i, dbl_frac_bits, dbl_exp_bits) - 1023;
            return (e >= 0) ? *this >>= e  : *this <<= e;
        }

        *this /= fixed_point128(val);
        return *this;
    }

    inline fixed_point128& operator%=(const fixed_point128& other) {
        uint64 nom[4] = {0, 0, low, high};
        uint64 denom[2] = {other.low, other.high};
        uint64 q[4] = {0}, r[4] = {0};
        
        //do the division in with positive numbers
        if (0 == div_32bit((uint32*)q, (uint32*)r, (uint32*)nom, (uint32*)denom, sizeof(nom) / sizeof(uint32), sizeof(denom) / sizeof(uint32))) {
            // simple case, both are integers (fractions is zero)
            if (is_int() && other.is_int()) {
                // result is in r (remainder) needs to shift left by frac_bits
                low = r[0];
                high = r[1];
            }
            // nom or denom are fractions
            // x mod y =  x - y * floor(x/y)
            else { 
                fixed_point128 x_div_y; // x / y. 
                x_div_y.high = (q[2] << upper_frac_bits) | (q[1] >> int_bits);
                x_div_y.low = (q[1] << upper_frac_bits) | (q[0] >> int_bits);
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
            INT_DIVIDE_BY_ZERO_EXCEPTION;
        }
        return *this;
    }

    inline fixed_point128& operator>>=(int shift) {
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

    inline fixed_point128& operator<<=(int shift) {
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

    inline fixed_point128& operator&=(const fixed_point128& other) {
        low &= other.low;
        high &= other.high;
        sign &= other.sign;
        // set sign to 0 when both low and high are zero (avoid having negative zero value)
        sign &= (0 != low || 0 != high);
        return *this;
    }

    inline fixed_point128& operator|=(const fixed_point128& other) {
        low |= other.low;
        high |= other.high;
        sign |= other.sign;
        return *this;
    }

    inline fixed_point128& operator^=(const fixed_point128& other) {
        low ^= other.low;
        high ^= other.high;
        sign ^= other.sign;
    }

    // prefix ++ operation
    inline fixed_point128& operator++() {
        high += unity;
        // set sign to 0 when both low and high are zero (avoid having negative zero value)
        sign &= (0 != low || 0 != high);
        return *this;
    }

    // postfix ++ operation
    inline fixed_point128 operator++(int) {
        fixed_point128 temp(*this);
        ++*this; // call the prefix implementation
        return temp;
    }

    // prefix -- operation
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

    // postfix -- operation
    inline fixed_point128 operator--(int) {
        fixed_point128 temp(*this);
        --*this; // call the prefix implementation
        return temp;
    }
    //
    // unary operations
    // 
    
    // conversion to bool
    inline operator bool() const {
        return high != 0 || low != 0;
    }
    // logical not
    inline bool operator!() const {
        return high == 0 && low == 0;
    }
    // bit wise not
    inline fixed_point128 operator~() const {
        fixed_point128 temp(*this);
        temp.high = ~high;
        temp.low = ~low;
        // set sign to 0 when both low and high are zero (avoid having negative zero value)
        temp.sign &= (0 != low || 0 != high);
        return temp;
    }
    // unary +
    inline fixed_point128 operator+() const {
        fixed_point128 temp(*this);
        return temp;
    }
    // unary -
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
    inline bool operator==(const fixed_point128& other) const {
        return sign == sign && high == high && low == low;
    }

    inline bool operator!=(const fixed_point128& other) const {
        return sign != sign || high != high || low != low;
    }

    inline bool operator<(const fixed_point128& other) const {
        // signs are different
        if (sign != other.sign)
            return sign > other.sign; // true when sign is 1 and other.sign is 0

        // MSB is the same, check the LSB
        if (high == other.high)
            return (sign) ? low > other.low : low < other.low;

        return (sign) ? high > other.high : high < other.high;
    }
    inline bool operator<=(const fixed_point128& other) const {
        // signs are different
        if (sign != other.sign)
            return sign > other.sign; // true when sign is 1 and other.sign is 0

        // MSB is the same, check the LSB
        if (high == other.high)
            return (sign) ? low >= other.low : low <= other.low;

        return (sign) ? high >= other.high : high <= other.high;
    }
    inline bool operator>(const fixed_point128& other) const {
        // signs are different
        if (sign != other.sign)
            return sign < other.sign; // true when sign is 0 and other.sign is 1

        // MSB is the same, check the LSB
        if (high == other.high)
            return (sign) ? low < other.low : low > other.low;

        return (sign) ? high < other.high : high > other.high;
    }
    inline bool operator>=(const fixed_point128& other) const {
        // signs are different
        if (sign != other.sign)
            return sign < other.sign; // true when sign is 0 and other.sign is 1

        // MSB is the same, check the LSB
        if (high == other.high)
            return (sign) ? low <= other.low : low >= other.low;

        return (sign) ? high <= other.high : high >= other.high;
    }

    // useful public functions
    inline bool is_int() const
    {
        return 0 == low && 0 == (high << int_bits);
    }
    inline bool is_positive() const
    {
        return 0 == sign;
    }

    // friends
    friend fixed_point128<int_bits> fp128::abs(const fixed_point128<int_bits>&) noexcept;
    friend fixed_point128<int_bits> fp128::floor(const fixed_point128<int_bits>&) noexcept;
    friend fixed_point128<int_bits> fp128::ciel(const fixed_point128<int_bits>&) noexcept;
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
inline int div_32bit(uint32* q, uint32* r, const uint32* u, const uint32* v, int64 m, int64 n)
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
    // Normalize by shifting v left just enough so that 
    // its high-order bit is on, and shift u left the 
    // same amount. We may have to append a high-order 
    // digit on the dividend; we do that unconditionally.
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
      // If the caller wants the remainder, unnormalize 
      // it and pass it back. 
    if (r != nullptr) {
        for (i = 0; i < n; i++) {
            r[i] = (un[i] >> s) | (un[i + 1] << s_comp);
        }
    }
    return 0;
}

/**
 * @brief returns the absolute value (sets sign to 0)
 * @param val - fixed_point128 element
 * @return - a copy of val with sign removed
*/
template<int int_bits>
inline fixed_point128<int_bits> abs(const fixed_point128<int_bits>& val) noexcept
{
    fixed_point128 temp = val;
    temp.sign = 0;
    return temp;
}

/**
 * @brief peforms the floor() function, similar to libc's floor(), rounds down towards -infinity.
 * @param val - input value
 * @return a fixed_point128 holding the integer value. Overflow is not reported.
*/
template<int int_bits>
inline fixed_point128<int_bits> floor(const fixed_point128<int_bits>& val) noexcept
{
    auto temp = val;
    temp.low = 0;
    temp.high &= temp.int_mask;
    // floor always rounds towards -infinity
    if (0 != temp.sign) {
        ++temp;
    }
    return temp;
}

/**
 * @brief peforms the ciel() function, similar to libc's ciel(), rounds up towards infinity.
 * @param val - input value
 * @return a fixed_point128 holding the integer value. Overflow is not reported.
*/
template<int int_bits>
inline fixed_point128<int_bits> ciel(const fixed_point128<int_bits>& val) noexcept
{
    auto temp = val;
    temp.low = 0;
    temp.high &= temp.int_mask;
    // ciel always rounds towards infinity
    if (0 == temp.sign) {
        ++temp;
    }
    return temp;
}


} //namespace fp128
#endif // #ifndef FIXED_POINT128_H
