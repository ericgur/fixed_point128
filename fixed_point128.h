#include <intrin.h>
#include <exception>


typedef __int64 int64;
typedef unsigned __int64 uint64;

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


template<int int_bits>
class fixed_point128
{
    static_assert(1 <= int_bits && int_bits <= 64, "Template parameter <int_bits> must be in the [1,64] range!");

private:
    // members
    union
    {
        char bytes[16];
        struct
        {
            uint64 low;
            uint64 high;
            unsigned sign; // 0 = positive, 1 negative
        };
    };
    static constexpr int frac_bits = 128 - int_bits;
    static constexpr int upper_frac_bits = (frac_bits <= 64) ? 0 : frac_bits - 64;
    static constexpr uint64 unity = 1ull << (64 - int_bits);
    static constexpr double upper_unity = (0 != upper_frac_bits) ? (1.0 / (double)(1ull << upper_frac_bits)) : 0;
    static constexpr double lower_unity = upper_unity / (double)(1ull << (32)) / (double)(1ull << (32));
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
        // too small. TODO: handle double's subnormal numbers
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
        high = (sign != 0) ? (-val) << upper_frac_bits : val << upper_frac_bits;
    }

    // conversion operators
    inline operator unsigned int() const {
        return (high >> upper_frac_bits) & max_dword_value;
    }

    inline operator int() const {
        return (int)((int64)high >> upper_frac_bits) & (max_dword_value);
    }

    inline operator double() const {
        double res = (double)(high * upper_unity); // high 64 bit part - bits 64-127
        res += (double)(low * lower_unity);
        if (sign) res = -res;
        return res;
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

    inline fixed_point128 operator>>(int shift) const {
        fixed_point128 temp(*this);
        return temp >> shift;
    }

    inline fixed_point128 operator<<(int shift) const {
        fixed_point128 temp(*this);
        return temp << shift;
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

        // multiply high QWORDs (overflow can happen)
        temp[0] = _umul128(high, other.high, &temp[1]);
        carry = _addcarry_u64(0, res[2], temp[0], &res[2]);
        res[2] += carry;
        res[3] += temp[1];

        // extract the bits from res[] keeping the precision the same as this object
        // shift result by frac_bits
        static constexpr int index = frac_bits >> 6; // / 64;
        static constexpr int lsb = frac_bits & MAX_BITS_VALUE_64(6);
        static constexpr int lsb_comp = 64 - lsb;
        //printf("index %i, lsb %i, lsb_comp %i", index, lsb, lsb_comp);
        // copy block #1 (lowest)
        low = res[index] >> lsb;
        low |= GET_BITS(res[index + 1], 0, lsb) << lsb_comp;
        high = res[index + 1] >> lsb;
        high |= GET_BITS(res[index + 2], 0, lsb) << lsb_comp;

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
        throw std::exception("Not implemented!");
        // use div128
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
            low = __shiftright128(low, high, shift);
            high >>= shift;
        }
        return *this;
    }

    inline fixed_point128& operator<<=(int shift) {
        if (shift >= 128) {
            low = high = 0;
            sign = 0;
        }
        else if (shift >= 64) {
            high = low << (shift - 64);
            low = 0;
        }
        // 0-64 bit shift
        else {
            high = __shiftleft128(low, high, shift);
            low <<= shift;
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
