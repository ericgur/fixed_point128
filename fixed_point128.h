#ifndef FIXED_POINT128_H
#define FIXED_POINT128_H

#include <intrin.h>
#include <exception>
#include <stdexcept>

typedef __int64 int64;
typedef unsigned __int64 uint64;
typedef int int32;
typedef unsigned int uint32;
typedef short int16;
typedef unsigned short uint16;

int div_short(unsigned short* q, unsigned short* r, const unsigned short* u, const unsigned short* v, int m, int n);
int div_32bit(uint32* q, uint32* r, const uint32* u, const uint32* v, int64 m, int64 n);

#pragma once
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

#define USE32BIT_DIVISION

template<int int_bits>
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
    static constexpr int upper_frac_bits = (frac_bits <= 64) ? 0 : frac_bits - 64;
    static constexpr uint64 unity = 1ull << (64 - int_bits);
    static inline double upper_unity = pow(2, -upper_frac_bits);
    static inline double lower_unity = pow(2, -frac_bits);
    static constexpr unsigned max_dword_value = (unsigned)(-1);
    static constexpr uint64 max_qword_value = (uint64)(-1);
public:
    fixed_point128() { 
        low = high = 0ull; sign = 0; 
    }
    fixed_point128(const fixed_point128& other) {
        low = other.low;
        high = other.high;
        sign = other.sign;
    }
    fixed_point128(double val) {
        uint64 i = *((uint64*)(&val));
        // very common case
        if (i == 0) {
            low = high = 0;
            sign = 0;
            return;
        }

        sign = GET_BIT(i, 63);
        int e = GET_BITS(i, 52, 11) - 1023;
        uint64 f = (i & MAX_BITS_VALUE_64(52));
        
        // overflow which catches NaN and Inf
        if (e >= int_bits) {
            high = low = max_qword_value;
            sign = 0;
            return;
        }

        // normal number
        if (e > -1023) {
            // bit 52 in f is the unity value of the float. it needs to move to the unity position in fixed point
            f |= ONE_SHIFT(52);
            int bits_to_shift = 64 - int_bits - 52 + e;

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
                    low <<= (int64)(64 - bits_to_shift);
                }
                // shift f into low QWORD
                else {
                    high = 0;
                    bits_to_shift -= 52;
                    f <<= (63 - 52);
                    low = f >> bits_to_shift;
                }
            }
        }
        // too small, no need to bother.
        else {
            high = low = 0;
            sign = 0;
        }
    }

    fixed_point128(uint64 val) {
        low = 0;
        sign = 0;
        high = val << upper_frac_bits;
    }

    fixed_point128(int64 val) {
        low = 0;
        sign = GET_BIT(val, 63);
        high = ((sign != 0) ? (-val) : val) << upper_frac_bits;
    }

    fixed_point128(unsigned val) {
        low = 0;
        sign = 0;
        high = (uint64)val << upper_frac_bits;
    }

    fixed_point128(int val) {
        low = 0;
        sign = GET_BIT(val, 31);
        high = (uint64)((sign != 0) ? (-val) : val) << upper_frac_bits;
    }

    // conversion operators
    inline operator unsigned int() const {
        return (high >> upper_frac_bits) & max_dword_value;
    }

    inline operator int() const {
        return (int)((int64)high >> upper_frac_bits) & (max_dword_value);
    }

    inline operator double() const {
        double res = (double)(high * upper_unity); // bits [64:127]
        res += (double)(low * lower_unity);        // bits [0:63]
        return (sign) ? -res : res;
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
        fixed_point128 temp(*this);
        temp /= val;
        return temp;
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
        uint64 other_low, other_high = ~other.high;
        carry = _addcarry_u64(0, ~other.low, 1, &other_low); // convert low to 2's complement
        other_high += carry;
        
        //add the other value
        carry = _addcarry_u64(0, low, other_low, &low); // convert low to 2's complement
        _addcarry_u64(carry, high, other_high, &high); // convert low to 2's complement

        // if result is is negative, invert it along with the sign.
        if (1 == GET_BIT(high, 63)) {
            sign ^= 1;
            carry = _addcarry_u64(0, ~low, 1, &low);
            high = ~high + carry;
        }
        
        return *this;
    }

    inline fixed_point128& operator*=(const fixed_point128& other) {
        uint64 res[4] = {0}; // 256 bit of result
        uint64 temp[2] = {0};
        unsigned char carry;

        // multiply low QWORDs
        res[0] = _umul128(low, other.low, &res[1]);

        // multiply high QWORDs (overflow can happen)
        res[2] = _umul128(high, other.high, &res[3]);

        // multiply low this and high other
        temp[0] = _umul128(low, other.high, &temp[1]);
        carry = _addcarry_u64(0, res[1], temp[0], &res[1]);
        carry = _addcarry_u64(carry, res[2], temp[1], &res[2]);
        res[3] += carry;

        // multiply high this and low other
        temp[0] = _umul128(high, other.low, &temp[1]);
        carry = _addcarry_u64(0, res[1], temp[0], &res[1]);
        carry = _addcarry_u64(carry, res[2], temp[1], &res[2]);
        res[3] += carry;

        // extract the bits from res[] keeping the precision the same as this object
        // shift result by frac_bits
        static constexpr int index = frac_bits >> 6; // / 64;
        static constexpr int lsb = frac_bits & MAX_BITS_VALUE_64(6);
        static constexpr int lsb_comp = 64 - lsb;
        //printf("index %i, lsb %i, lsb_comp %i", index, lsb, lsb_comp);
        // copy block #1 (lowest)
        low = res[index + 1] << lsb_comp;
        low |= (res[index] >> lsb) & MAX_BITS_VALUE_64(lsb_comp);

        high = res[index + 2] << lsb_comp;
        high |= (res[index + 1] >> lsb) & MAX_BITS_VALUE_64(lsb_comp);

        // set the sign
        sign ^= other.sign;
        return *this;
    }

    inline fixed_point128& operator*=(int64 val) {
        // alway do positive multiplication
        if (val < 0) {
            val = -val;
            sign ^= 1;
        }
        uint64 uval = uint64(val);
        uint64 temp;

        // multiply low QWORDs
        low = _umul128(low, uval, &temp);
        high = high * uval + temp;
        return *this;
    }

    inline fixed_point128& operator/=(const fixed_point128& other) {
        uint64 nom[4] = {0, 0, low, high};
        uint64 denom[2] = {other.low, other.high};
        uint64 q[4] = {0}, r = NULL; // don't need the reminder
    #ifdef USE32BIT_DIVISION
        if (0 == div_32bit((uint32*)q, (uint32*)r, (uint32*)nom, (uint32*)denom, sizeof(nom) / sizeof(uint32), sizeof(denom) / sizeof(uint32))) {
    #else
        if (0 == div_short((unsigned short*)q, (unsigned short*)r, (unsigned short*)nom, (unsigned short*)denom, sizeof(nom) / sizeof(short), sizeof(denom) / sizeof(short))) {
    #endif
            // result in q needs to shift left by frac_bits
            high = (q[2] << upper_frac_bits) | (q[1] >> int_bits);
            low  = (q[1] << upper_frac_bits) | (q[0] >> int_bits);
            sign ^= other.sign;
        }
        else { // error
            throw std::logic_error("Divide by zero!");
        }
        return *this;
    }

    inline fixed_point128& operator/=(double val) {
        uint64 i = *((uint64*)(&val));
        // infinity
        if (i == 0) {
            low = high = max_qword_value;
            sign = 0;
            return *this;
        }

        uint64 f = (i & MAX_BITS_VALUE_64(52));
        // simple and common case, the value is an exponent of 2
        if (0 == f) {
            sign ^= int(i >> 63);
            int e = GET_BITS(i, 52, 11) - 1023;
            return (e >= 0) ? *this >>= e  : *this <<= e;
        }
        *this /= fixed_point128(val);
        return *this;
    }

    inline fixed_point128& operator>>=(int shift) {
        if (shift >= 128) {
            low = high = 0;
            sign = 0;
        }
        else if (shift >= 64) {
            low = high >> (shift - 64);
            high = 0;
        }
        // 0-64 bit shift
        else {
            low = __shiftright128(low, high, (unsigned char)shift);
            high >>= shift;
        }
        return *this;
    }

    inline fixed_point128& operator<<=(int shift) {
        // 0-64 bit shift - most common
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
            sign = 0;
        }

        return *this;
    }

    inline fixed_point128& operator&=(const fixed_point128& other) {
        low &= other.low;
        high &= other.high;
        sign &= other.sign;
        return *this;
    }
    inline fixed_point128& operator|=(const fixed_point128& other) {
        low |= other.low;
        high |= other.high;
        sign |= other.sign;
        return *this;
    }

    // prefix ++ operation
    inline fixed_point128& operator++() {
        high += unity;
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
        temp.sign ^= 1;
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

};


inline int div_32bit(uint32* q, uint32* r, const uint32* u, const uint32* v, int64 m, int64 n)
{
    const uint64 b = 1ull << 32; // Number base (32 bits).
    const uint64 mask = b - 1;
    uint32* un, * vn;            // Normalized form of u, v.
    uint64 qhat;               // Estimated quotient digit.
    uint64 rhat;               // A remainder.
    uint64 p;                  // Product of two digits.
    int64 s, s_comp, i, j, t, k;
    if (m < n || n <= 0 || v[n - 1] == 0)
        return 1; // Return if invalid param.

    if (n == 1) {                           // Take care of 
        k = 0;                              // the case of a 
        for (j = m - 1; j >= 0; j--) {      // single-digit 
            q[j] = (uint32)((k * b + u[j]) / v[0]);   // divisor here. 
            k = (k * b + u[j]) - (uint64)q[j] * v[0];
        }

        if (r != NULL)
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
    if (NULL == vn) return 1;

    for (i = n - 1; i > 0; i--)
        vn[i] = (uint32)(v[i] << s) | (v[i - 1] >> s_comp);
    vn[0] = v[0] << s;
    un = (uint32*)_malloca(sizeof(uint32) * (m + 1));
    if (NULL == un) return 1;
    un[m] = u[m - 1] >> s_comp;
    for (i = m - 1; i > 0; i--)
        un[i] = (u[i] << s) | (u[i - 1] >> s_comp);
    un[0] = u[0] << s;
    for (j = m - n; j >= 0; j--) { // Main loop. 
                                   // Compute estimate qhat of q[j]. 
        qhat = ((uint64)un[j + n] * b + (uint64)un[j + n - 1]) / (uint64)vn[n - 1];
        rhat = ((uint64)un[j + n] * b + (uint64)un[j + n - 1]) - qhat * (uint64)vn[n - 1];
    again:
        if (qhat >= b || qhat * (uint64)vn[n - 2] > b * rhat + (uint64)un[j + n - 2]) {
            qhat = qhat - 1;
            rhat = rhat + (uint64)vn[n - 1];
            if (rhat < b)
                goto again;
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
        if (t < 0) {                    // If we subtracted too

            q[j] = q[j] - 1; // much, add back. 
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
    if (r != NULL) {
        for (i = 0; i < n; i++)
            r[i] = (un[i] >> s) | (un[i + 1] << s_comp);

    }
    return 0;
}

inline int div_short(unsigned short* q, unsigned short* r, const unsigned short* u, const unsigned short* v, int m, int n)
{
    const unsigned b = 65536; // Number base (16 bits).
    const unsigned mask = b - 1;
    unsigned short* un, * vn;  // Normalized form of u, v.
    unsigned qhat;            // Estimated quotient digit.
    unsigned rhat;            // A remainder.
    unsigned p;               // Product of two digits.
    int s, i, j, t, k;
    if (m < n || n <= 0 || v[n - 1] == 0)
        return 1; // Return if invalid param.

    if (n == 1) {                           // Take care of 
        k = 0;                              // the case of a 
        for (j = m - 1; j >= 0; j--) {      // single-digit 
            q[j] = (unsigned short)((k * b + u[j]) / v[0]);   // divisor here. 
            k = (k * b + u[j]) - q[j] * v[0];
        }

        if (r != NULL)
            r[0] = (unsigned short)k;
        return 0;
    }
    // Normalize by shifting v left just enough so that 
    // its high-order bit is on, and shift u left the 
    // same amount. We may have to append a high-order 
    // digit on the dividend; we do that unconditionally.
    s = __lzcnt16(v[n - 1]); // -16; // 0 <= s <= 16. 
    vn = (unsigned short*)_malloca(2 * n);
    if (NULL == vn)
        return 1;

    for (i = n - 1; i > 0; i--)
        vn[i] = (unsigned short)(v[i] << s) | (v[i - 1] >> (16 - s));
    vn[0] = v[0] << s;
    un = (unsigned short*)_malloca(2 * (m + 1));
    if (NULL == un)
        return 1;

    un[m] = u[m - 1] >> (16 - s);
    for (i = m - 1; i > 0; i--)
        un[i] = (u[i] << s) | (u[i - 1] >> (16 - s));
    un[0] = u[0] << s;
    for (j = m - n; j >= 0; j--) { // Main loop. 
                                   // Compute estimate qhat of q[j]. 
        qhat = (un[j + n] * b + un[j + n - 1]) / vn[n - 1];
        rhat = (un[j + n] * b + un[j + n - 1]) - qhat * vn[n - 1];
    again:
        if (qhat >= b || qhat * vn[n - 2] > b * rhat + un[j + n - 2]) {
            qhat = qhat - 1;
            rhat = rhat + vn[n - 1];
            if (rhat < b)
                goto again;
        }
        // Multiply and subtract. 
        k = 0;
        for (i = 0; i < n; i++) {
            p = qhat * vn[i];
            t = un[i + j] - k - (p & mask);
            un[i + j] = (unsigned short)t;
            k = (p >> 16) - (t >> 16);
        }

        t = un[j + n] - k; un[j + n] = (unsigned short)t;
        q[j] = (unsigned short)qhat;    // Store quotient digit. 
        if (t < 0) {                    // If we subtracted too

            q[j] = q[j] - 1; // much, add back. 
            k = 0;
            for (i = 0; i < n; i++) {
                t = un[i + j] + vn[i] + k;
                un[i + j] = (unsigned short)t;
                k = t >> 16;
            }
            un[j + n] = (unsigned short)(un[j + n] + k);
        }
    } // End j.
      // If the caller wants the remainder, unnormalize 
      // it and pass it back. 
    if (r != NULL) {
        for (i = 0; i < n; i++)
            r[i] = (un[i] >> s) | (un[i + 1] << (16 - s));

    }
    return 0;
}

#endif // #ifndef FIXED_POINT128_H