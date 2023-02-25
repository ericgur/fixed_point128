/***********************************************************************************
    MIT License

    Copyright (c) 2023 Eric Gur (ericgur@iname.com)

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

#ifndef FLOAT128_H
#define FLOAT128_H

// override some static analysis checks
#pragma warning(push)
#pragma warning(disable: 26472) // Don't use a static_cast for arithmetic conversions. Use brace initialization
#pragma warning(disable: 26485) // No array to pointer decay
#pragma warning(disable: 26481) // Don't use pointer arithmetic. Use span instead
#pragma warning(disable: 26446) // Prefer to use gsl::at() instead of unchecked subscript operator
#pragma warning(disable: 26482) // Only index into arrays using constant expressions
#pragma warning(disable: 26408) // Avoid malloc() and free(), prefer the nothrow version of new with delete
#pragma warning(disable: 6255)  // _alloca indicates failure by raising a stack overflow exception.  Consider using _malloca instead
#pragma warning(disable: 4996)  // This function or variable may be unsafe.Consider using strncpy_s instead.To disable deprecation, use _CRT_SECURE_NO_WARNINGS.See online help for details.


#include "fixed_point128_shared.h"

namespace fp128 {
/***********************************************************************************
*                                  Forward declarations
************************************************************************************/
class fp128_gtest; // Google test class
class float128;

// CRT style functions
float128 fabs(const float128& x) noexcept;
float128 floor(const float128& x) noexcept;
float128 ceil(const float128& x) noexcept;
float128 trunc(const float128& x) noexcept;
float128 round(const float128& x) noexcept;
int32_t ilogb(const float128& x) noexcept;
float128 copysign(const float128& x, const float128& y) noexcept;
float128 fmod(const float128& x, const float128& y) noexcept;
float128 modf(const float128& x, float128* iptr) noexcept;
float128 fdim(const float128& x, const float128& y) noexcept;
float128 fmin(const float128& x, const float128& y) noexcept;
float128 fmax(const float128& x, const float128& y) noexcept;
float128 hypot(const float128& x, const float128& y) noexcept;
float128 cbrt(const float128 x, uint32_t iterations = 1) noexcept;
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
float128 pow(const float128& x, const float128& y, int32_t f = 112) noexcept;
float128 pow(const float128& x, int32_t y) noexcept;
float128 log(float128 x, int32_t f = 112) noexcept;
float128 log2(float128 x, int32_t f = 112) noexcept;
float128 log10(float128 x, int32_t f = 112) noexcept;
float128 logb(float128 x, int32_t f = 112) noexcept;
float128 log1p(float128 x, int32_t f = 112) noexcept;
// non CRT function
float128 reciprocal(const float128& x) noexcept;
void fact_reciprocal(int x, float128& res) noexcept;
float128 double_factorial(int x) noexcept;

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
*/

struct _float128_bits {
    uint64_t f : 48;
    uint64_t e : 15;
    uint64_t s : 1;
};

enum float128_class_t {
    signalingNaN,
    quietNaN,
    negativeInfinity,
    negativeNormal,
    negativeSubnormal,
    negativeZero,
    positiveZero,
    positiveSubnormal,
    positiveNormal,
    positiveInfinity
};

class __declspec(align(16)) float128
{
    // build time validation of template parameters
    static_assert(sizeof(void*) == 8, "float128 is supported in 64 bit builds only!");
    friend class fp128_gtest;
    
    static constexpr int32_t EXP_BITS = 15;
    static constexpr int32_t EXP_BIAS = 0x3FFF;
    static constexpr int32_t ZERO_EXP_BIASED = -EXP_BIAS;
    static constexpr int32_t ZERO_EXP_UNBIASED = 0;
    static constexpr int32_t SUBNORM_EXP_BIASED = 0;
    static constexpr int32_t SUBNORM_EXP_UNBIASED = -EXP_BIAS;
    static constexpr int32_t INF_EXP_BIASED = 0x7FFF;
    static constexpr int32_t INF_EXP_UNBIASED = INF_EXP_BIASED - EXP_BIAS;
    static constexpr uint64_t EXP_MASK = INF_EXP_BIASED;
    static constexpr int32_t FRAC_BITS = 112;
    static constexpr uint64_t UPPER_FRAC_MASK = FP128_MAX_VALUE_64(FRAC_BITS - 64);
    static constexpr uint64_t FRAC_UNITY = FP128_ONE_SHIFT(FRAC_BITS - 64);
    static constexpr uint64_t SIGN_MASK = 1ull << 63;
#pragma warning(push)
#pragma warning(disable: 4201) // nameless union/structs
    struct {
        uint64_t low;
        union {
            uint64_t high;
            _float128_bits high_bits; 
        };
    };
#pragma warning(pop)
public:
    /**
     * @brief Default constructor, creates an instance with a value of zero.
    */
    __forceinline constexpr float128() noexcept : 
        low(0), high(0) {}
    /**
     * @brief Copy constructor
     * @param other Object to copy from
    */
    __forceinline float128(const float128& other) noexcept :
        low(other.low), high(other.high) {}
    /**
     * @brief Move constructor
     * Doesn't modify the right hand side object. Acts like a copy constructor.
     * @param other Object to copy from
    */
    __forceinline float128(const float128&& other) noexcept :
        low(other.low), high(other.high) {}
    /**
     * @brief Low level constructor
     * @param l Low QWORD
     * @param h High QWORD
    */
    __forceinline constexpr float128(uint64_t l, uint64_t h) noexcept :
        low(l), high(h) {}
    /**
     * @brief Low level constructor
     * @param lf Low fraction part (bits 63:0)
     * @param hf High fraction part (bits 111:64)
     * @param e Exponent
     * @param s sign
    */
    __forceinline float128(uint64_t lf, uint64_t hf, uint32_t e, uint32_t s) noexcept :
        low(lf) {
        high_bits.f = hf;
        high_bits.e = e;
        high_bits.s = s;
    }
    /**
     * @brief Constructor from the double type
     * @param x Input value
    */
    FP128_INLINE float128(double x) noexcept {
        low = high = 0;
        // very common case
        if (x == 0) return;

        // hack the double bit fields
        const Double d(x);
        
        // subnormal numbers
        if (d.e == 0) {
            // the exponent is -1022 (1-1023)
            auto msb = 64 - static_cast<int32_t>(_lzcnt_u64(d.f));
            // exponent
            int32_t x_expo = static_cast<int32_t>(d.e) - 1023;
            int32_t expo = x_expo + msb - dbl_frac_bits;
            //fraction
            low = d.f & ~(1ull << (msb - 1)); // clear the msb
            auto shift = static_cast<int32_t>(FRAC_BITS - msb + 1);
            shift_left128_inplace_safe(low, high, shift);
            set_exponent(expo);
        }
        // NaN & INF
        else if (d.e == 0x7FF)
        {
            high_bits.e = INF_EXP_BIASED;
            // zero for +- INF, non-zero for NaN 
            high_bits.f = (d.f) ? 1 : 0;
        }
        // normal numbers
        else {
            low = d.f << 60;
            high = d.f >> 4;
            set_exponent(static_cast<int32_t>(d.e) - 1023);
        }

        // copy the sign
        set_sign(d.s);
    }
    template<typename T>
    __forceinline float128(T x) noexcept {
        if constexpr (std::is_floating_point_v<T>) {
            new(this) float128(static_cast<double>(x));
            return;
        }
        else if constexpr (std::is_same_v<char*, T> ||
                           std::is_same_v<unsigned char*, T> ||
                           std::is_same_v<const unsigned char*, T>) {
            new(this) float128(static_cast<const char*>(x));
            return;
        }
        else if constexpr (std::is_integral_v<T>) {
            uint64_t sign = 0;
            if constexpr (std::is_signed_v<T>) {
                // alway do positive multiplication
                if (x < 0) {
                    x = -x;
                    sign = 1;
                }
            }

            // integers: convert to uint64 for a simpler operation.
            high = 0;
            low = static_cast<uint64_t>(x);
            if (low == 0)
                return;

            auto expo = log2(low); // this is the index of the msb as well
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
    FP128_INLINE float128(const char* x) noexcept {
        low = high = 0;
        if (x == nullptr) return;

        constexpr uint64_t base16_max_digits = (112 + 4) / 4;  // 29 hex digits. 28 for the fraction (112 bit) and another for the unity
        constexpr uint64_t base10_max_digits = 35;             // maximum for 112 bit of manstissa/fraction is 34, read one extra to get maximum precision
        uint32_t sign = 0;
        uint32_t base = 10;
        int32_t expo2 = 0;         // base2 exponent
        int32_t expo10 = 0;        // base10 exponent

        // convert the input string to lowercase for simpler processing.
        const auto x_len = 1 + strlen(x);
        char* str = static_cast<char*>(_alloca(x_len));
        strncpy(str, x, x_len);
        _strlwr_s(str, x_len);

        char* p = str;
        if (p == nullptr) return;

        // skip white space
        while (*p == ' ') ++p;

        if (*p == '\0') {
            *this = float128();
            return;
        }
        // set negative sign if needed
        if (*p == '-') {
            sign = 1;
            ++p;
        }
        else if (*p == '+')
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
        while (*p == '0') ++p;

        int32_t int_digits = 0, frac_digits = 0;
        char* int_start = p;
        char* frac_start = nullptr;

        // count the integer digits
        while (isdigit(*p) || (base == 16 && *p >= 'a' && *p <= 'f')) {
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

            uint64_t* frac_bits = &low; // fill the internal data structure directly
            // start at the leftmost digit and iterate right

            int32_t digits_consumed = min(int_digits, max_digits);
            char* cur_digit = int_start;
            char* const end_digit = int_start + digits_consumed;

            // fill the internal structure starting the the top bits of high
            while (cur_digit < end_digit) {
                uint64_t d = *cur_digit;
                if (d >= '0' && d <= '9')
                    d -= '0';
                else
                    d = 10ull + d - 'a';

                ++cur_digit;
                auto index = 1 - (expo2 >> 6); // fill high part first
                assert(index >= 0 && index <= 1);
                auto shift = 60 - (expo2 & 63);
                frac_bits[index] |= d << shift;
                expo2 += 4;
            }
            // shift the result into position and fix the exponent
            uint64_t left_digit = high >> 60;
            assert(left_digit != 0);
            static constexpr int32_t digit_msb_lut[16] = { 0, 0, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3 };
            const auto digit_msb = digit_msb_lut[left_digit];
            expo2 -= 4 - digit_msb;
            expo2 += 4 * (int_digits - digits_consumed); // account for digits which do not fit in the fraction.

            // overflow
            if (expo2 > INF_EXP_UNBIASED) {
                *this = inf();
            }
            else {
                shift_right128_inplace_safe(low, high, 124 + digit_msb - FRAC_BITS); // move digit's 2nd  msb to bit 111
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
            frac_start = (isdigit(*p)) ? p : nullptr;

            // count the fraction digits if they exist
            if (frac_start) {
                while (isdigit(*p)) {
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
            // integer part is zero - skip the leading zeros in the fraction and ajust the exponent
            // example 0.01 == 0.1E-1
            //if (int_digits == 0 && frac_start != nullptr) {
            //    while (*frac_start == '0') {
            //        ++frac_start;
            //        --frac_digits;
            //        --expo10_adjust;
            //    }
            //}
        }

        // check for the optional exponent
        if (*p == 'e') {
            *p = '\0'; // terminate the fraction string
            ++p;
            // convert the exponent
            expo10 = strtol(p, nullptr, 10);

            // underflow
            if (expo10 < -4965) {
                return;
            }

            // overflow
            if (expo10 > 4932) {
                *this == inf();
            }
        }

        // compute the integer part
        if (base == 10) {
            uint128_t int_part;
            int32_t digits_consumed = min(int_digits, max_digits);
            char* const end_digit = int_start + digits_consumed;
            *end_digit = '\0';
            int_part = int_start;
            int32_t shift_bits = 0;

            // mark the extra exponent that may exist if we had enough bits to represent the entire value
            auto extra_digits = int_digits - digits_consumed;

            // multiply by 10 for each digit
            while (extra_digits > 0) {
                --extra_digits;
                uint64_t msb = log2(int_part);
                if (msb >= 123) {
                    int_part >>= 4;
                    shift_bits += 4;
                    assert(msb - log2(int_part) == 4);
                }
                int_part *= 10;
            }

            int32_t msb = static_cast<int32_t>(log2(int_part));
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
                int_part.get_components(low, high); // unity bit is erased by set_exponent()
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
            // take the minimum of the actual digits in the string versus what is the maximum possible to hold in 112 bit.
            int32_t digits = min(frac_digits, max_digits - int_digits);
            constexpr int32_t digit_group = 9;
            static_assert(digit_group <= 9); // must fit in 32 bit
            int32_t i = 0;
            //const float128 group_base = exp10(-digit_group);
            const float128 group_base = float128(0x4b2e62d01511f12a, 0x12e0be826d69, EXP_BIAS - 30, 0); // 10^-9
            float128 current_base = 1;

            for (; i < digits; i += digit_group) {
                // convert to an int
                uint32_t val = 0;
                for (auto j = 0; j < digit_group; ++j) {
                    val *= 10;
                    if (i + j < digits)
                        val += frac_start[i + j] - '0';
                }
                
                current_base *= group_base;
                frac_part += current_base * val;
            }

            //// handle the remaining digits one at a time
            //for (; i < digits; ++i) {
            //    // convert to an int
            //    uint32_t val = frac_start[i] - '0';
            //    current_base *= tenth();
            //    if (val == 0)
            //        continue;
            //    frac_part += current_base * val;
            //}

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
    __forceinline float128(const std::string& x) noexcept {
        float128 temp = x.c_str();
        *this = temp;
    }
    /**
     * @brief Destructor
    */
    __forceinline ~float128() {}
    /**
     * @brief Assignment operator
     * @param other Object to copy from
     * @return This object.
    */
    __forceinline float128& operator=(const float128& other) noexcept {
        high = other.high;
        low = other.low;
        return *this;
    }
    /**
     * @brief Move assignment operator
     * @param other Object to copy from
     * @return This object.
    */
    __forceinline float128& operator=(const float128&& other) noexcept {
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
    FP128_INLINE operator double() const noexcept {
        Double d{};

        // nan and inf
        if (high_bits.e == INF_EXP_BIASED) {
            // zero fraction means inf, otherwise nan
            if (low == 0 && high_bits.f == 0) {
                return (get_sign()) ? -HUGE_VAL : HUGE_VAL;
            }
            return NAN;
        }

        int32_t expo = get_exponent();
        // subnormal and underflow
        if (expo < -1022) {
            int32_t shift = (int32_t)(FRAC_BITS - dbl_frac_bits) -1022 - expo;
            
            // underflow
            if (shift >= FRAC_BITS) {
                return 0;
            }

            d.e = 0;
            // add the msb back
            uint64_t h = high_bits.f | (1ull << (FRAC_BITS - 64));
            d.f = shift_right128_round(low, h, shift);
            d.s = high_bits.s;
            return d.val;
        }
        // too big for double
        else if (expo > 1023) {
            return HUGE_VAL;
        }
        // normal numbers
        else {
            d.e = 1023ull + expo;
            d.f = shift_right128_round(low, high_bits.f, FRAC_BITS - dbl_frac_bits);
            d.s = high_bits.s;
            
            // fraction caused a round up
            if (d.f == 0 && high_bits.f != 0)
                d.e += 1;

            return d.val;
        }
    }
    /**
     * @brief operator float converts to a float
    */
    FP128_INLINE operator float() const noexcept {
        // TODO: write proper function
        double v = static_cast<double>(*this);
        return static_cast<float>(v);
    }
    /**
     * @brief operator uint64_t converts to a uint64_t
    */
    FP128_INLINE operator uint64_t() const noexcept {
        uint64_t l, h;
        int32_t e;
        uint32_t s;
        get_components(l, h, e, s);

        if (e > 63) return (s) ? static_cast<uint64_t>(INT64_MIN) : UINT64_MAX;
        if (e < 0) return 0;
        
        auto shift = static_cast<int>(FRAC_BITS) - e;
        shift_right128_inplace_safe(l, h, shift);
        return l;
    }
    /**
     * @brief operator int64_t converts to a int64_t
    */
    FP128_INLINE operator int64_t() const noexcept {
        uint64_t l, h;
        int32_t e;
        uint32_t s;
        get_components(l, h, e, s);

        if (e > 62) return (s) ? INT64_MIN : INT64_MAX;
        if (e < 0) return 0;

        auto shift = static_cast<int>(FRAC_BITS) - e;
        shift_right128_inplace_safe(l, h, shift);
        int64_t res = static_cast<int64_t>(l);
        return  (s) ? -res : res;
    }
    /**
     * @brief operator uint32_t converts to a uint32_t
    */
    FP128_INLINE operator uint32_t() const noexcept {
        uint64_t l, h;
        int32_t e;
        uint32_t s;
        get_components(l, h, e, s);

        if (e > 31) return (s) ? static_cast<uint32_t>(INT32_MIN) : UINT32_MAX;
        if (e < 0) return 0;

        auto shift = static_cast<int>(FRAC_BITS) - e;
        shift_right128_inplace_safe(l, h, shift);
        return static_cast<uint32_t>(l);
    }
    /**
     * @brief operator int32_t converts to a int32_t
    */
    FP128_INLINE operator int32_t() const noexcept {
        uint64_t l, h;
        int32_t e;
        uint32_t s;
        get_components(l, h, e, s);

        if (e > 30) return (s) ? INT32_MIN : INT32_MAX;
        if (e < 0) return 0;

        auto shift = static_cast<int>(FRAC_BITS) - e;
        shift_right128_inplace_safe(l, h, shift);
        int32_t res = static_cast<int32_t>(l);
        return  (s) ? -res : res;
    }
    /**
      * @brief operator long double - converts to a long double
      * @return Object value.
     */
    FP128_INLINE operator long double() const noexcept {
        return operator double();
    }
    /**
     * @brief Converts to a std::string (slow) string holds all meaningful fraction bits.
     * @return object string representation
    */
    FP128_INLINE operator std::string() const noexcept {
        return operator char* ();
    }
    /**
     * @brief Converts to a C string (slow) string holds all meaningful fraction bits.
     * @return object string representation
    */
    explicit FP128_INLINE operator char* () const noexcept {
        constexpr int32_t buff_size = 128;
        static thread_local char str[buff_size]; // need roughly a (meaningful) decimals digit per 3.2 bits
        char tmp_buf[buff_size]; // used to write the integer part
        char* s = str; // needed for debugging

        if (is_special()) {
            if (is_nan()) {
                strcpy(str, "nan");
            }
            else if (is_inf()) {
                sprintf(str, "%sinf", (is_negative()) ? "-" : "");
            }
            return str;
        }
        if (is_zero()) {
            strcpy(str, "0");
            return str;
        }

        // TODO: find base 10 exponent via log10
        auto expo = get_exponent();
        constexpr int32_t max_exponent = 127;
        constexpr int32_t min_exponent = -FRAC_BITS;
        if (expo >= max_exponent || expo <= min_exponent) {
            to_e_format(str, buff_size);
            return str;
        }
        
        char* p = str; 
        if (get_sign()) {
            *p++ = '-';
        }
        // convert the integer part
        float128 int_part, frac_part, ten = 10;
        frac_part = modf(fabs(*this), &int_part);
        int32_t int_digits = 0;
        char* pp = tmp_buf;
        do {
            auto digit = static_cast<uint32_t>(fmod(int_part, ten));
            *pp++ = static_cast<char>(digit + '0');
            int_part /= 10;
            ++int_digits;
        } while (int_digits < buff_size && int_part >= 1);
        *pp = '\0';
        for (auto i = int_digits; i != 0; --i) {
            *p++ = *--pp;
        }

        // not enough digits.
        if (buff_size == int_digits) {
            str[buff_size - 4] = str[buff_size - 3] = str[buff_size - 2] = '.';
            str[buff_size - 1] = '\0';
            return str;
        }

        // convert the fraction part
        if (!frac_part.is_zero()) {
            *p++ = '.';
            ++int_digits;
            constexpr int32_t digit_group = 5, max_frac_digits = 35;
            static_assert(digit_group <= 9); // must fit in 32 bit
            //float128 base = exp10(digit_group);
            float128 base = 100000;
            char fmt[10];
            sprintf(fmt, "%%0%ii", digit_group);
            int32_t groups = min(max_frac_digits, (buff_size - int_digits - 1)) / digit_group;
            while (groups-- > 0) {
                frac_part *= base;
                frac_part = modf(frac_part, &int_part);
                uint32_t val = static_cast<uint32_t>(int_part);
                p += sprintf(p, fmt, val);
            }

            // back track and remove trailing zero
            while (p[-1] == '0')
                --p;

        }

        // convert the fraction part
        //if (!frac_part.is_zero()) {
        //    *p++ = '.';
        //    ++int_digits;
        //    for (auto i = 0; i < 35; ++i) {
        //        frac_part *= 10;
        //        frac_part = modf(frac_part, &int_part);
        //        uint32_t digit = static_cast<uint32_t>(int_part);
        //        *p++ = static_cast<char>(digit + '0');
        //    }

        //    // back track and remove trailing zero
        //    while (p[-1] == '0') 
        //        --p;

        //}

        //// there are fraction bits left to print
        //if (frac_bits > 0) {
        //    uint32_t digits = (uint32_t)((double)frac_bits / 3.29);
        //    *p++ = '.';
        //    uint128_t mask(UINT64_MAX, 0xFFFFFFFFFFFF); // keeps lower 112 bits
        //    frac_part <<= (128 - frac_bits); // erase the remaining integer bits
        //    frac_part >>= 16;
        //    while (frac_part != 0 && digits-- > 0) {
        //        frac_part *= 10;
        //        int32_t digit = frac_part >> FRAC_BITS;
        //        frac_part &= mask;
        //        if (digit > 0 || frac_part != 0)
        //            *p++ = static_cast<char>(digit + '0');
        //    }
        //}

        //auto msb = 

        *p = '\0';
        return s;
    }
    /**
     * @brief converts the stored value to a string with scientific notation
     * @param str Output buffer
     * @param buff_size Output buffer size in bytes
    */
    void to_e_format(char* str, int32_t buff_size) const {
        UNREFERENCED_PARAMETER(buff_size);
        strcpy(str, "e format not supported yet");
    }

    //
    // math operators
    //
    /**
     * @brief Shift right this object.
     * @param shift Bits to shift. Values less than 1 do nothing, high values can cause the value to reach zero.
     * @return This object.
    */
    FP128_INLINE float128& operator>>=(int32_t shift) noexcept {
        if (shift < 1)
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
    FP128_INLINE float128& operator<<=(int32_t shift) noexcept {
        if (shift < 1)
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
     * @param shift bits to shift
     * @return Temporary object with the result of the operation
    */
    template<typename T>
    __forceinline float128 operator>>(T shift) const noexcept {
        float128 temp(*this);
        return temp >>= static_cast<int32_t>(shift);
    }
    /**
     * @brief Performs left shift operation.
     * @param shift bits to shift
     * @return Temporary object with the result of the operation
    */
    template<typename T>
    __forceinline float128 operator<<(T shift) const noexcept {
        float128 temp(*this);
        return temp <<= static_cast<int32_t>(shift);
    }

    /**
     * @brief Add a value to this object
     * @param other Right hand side operand
     * @return This object.
    */
    FP128_INLINE float128& operator+=(const float128& other) noexcept {
        // check trivial cases
        if (high_bits.e == INF_EXP_BIASED || other.high_bits.e == INF_EXP_BIASED) {
            if (is_nan() || other.is_nan())
                *this = nan();
            // return inf with the right sign
            if (other.is_inf())
                *this = other;
            return *this;
        }

        uint32_t sign, other_sign;
        int32_t expo, other_expo;
        uint64_t l1, h1, l2, h2;
        get_components(l1, h1, expo, sign);
        other.get_components(l2, h2, other_expo, other_sign);
        constexpr int32_t shift_left_bits = 127 - 2 - FRAC_BITS; // move bit 112 left, keep 1 bit for additional exponent and one for the sign
        
        if (expo > other_expo) {
            int32_t shift = expo - other_expo - shift_left_bits; // how many bits to shift right
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
        }
        else if (expo < other_expo) {
            int32_t shift = other_expo - expo - shift_left_bits; // how many bits to shift right
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
            //add the other value
            const uint8_t carry = _addcarryx_u64(0, l1, l2, &l1);
            _addcarryx_u64(carry, h1, h2, &h1);
        }
        // different sign: invert the sign for other and subtract
        else {
            // this value is negative
            if (high_bits.s)
                twos_complement128(l1, h1);
            // other value is negative
            else
                twos_complement128(l2, h2);

            //add the other value, results stored in l1
            const uint8_t carry = _addcarryx_u64(0, l1, l2, &l1);
            _addcarryx_u64(carry, h1, h2, &h1);

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
    template<typename T>
    FP128_INLINE float128& operator+=(const T& other) {
        return operator+=(float128(other));
    }
    /**
     * @brief Subtract a value from this object
     * @param other Right hand side operand
     * @return This object.
    */
    FP128_INLINE float128& operator-=(const float128& other) noexcept {
        return *this+=(-other);
    }
    /**
     * @brief Subtract a value from this object
     * @param other Right hand side operand
     * @return This object.
    */
    template<typename T>
    FP128_INLINE float128& operator-=(const T& other) {
        return operator+=(-float128(other));
    }
    /**
     * @brief Multiply a value to this object
     * @param other Right hand side operand
     * @return This object.
    */
    FP128_INLINE float128& operator*=(const float128& other) noexcept {
        // check trivial cases
        if (is_special() || other.is_special()) {
            if (is_nan() || other.is_nan())
                *this = nan();
            // return inf with the right sign
            if (other.is_inf())
                *this = other;
            return *this;
        }
        else if (is_zero() || other.is_zero()) {
            *this = 0;
            return *this;
        }
        // extract fractions and exponents
        uint32_t sign, other_sign;
        int32_t expo, other_expo;
        uint64_t l1, h1, l2, h2;
        get_components(l1, h1, expo, sign);
        other.get_components(l2, h2, other_expo, other_sign);
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
        uint64_t res[4]; // 256 bit of result
        uint64_t temp1[2], temp2[2];

        // multiply low QWORDs
        res[0] = _mulx_u64(l1, l2, &res[1]);

        // multiply high QWORDs (overflow can happen)
        res[2] = _mulx_u64(h1, h2, &res[3]);

        // multiply low this and high other
        temp1[0] = _mulx_u64(l1, h2, &temp1[1]);
        uint8_t carry = _addcarryx_u64(0, res[1], temp1[0], &res[1]);
        res[3] += _addcarryx_u64(carry, res[2], temp1[1], &res[2]);

        // multiply high this and low other
        temp2[0] = _mulx_u64(h1, l2, &temp2[1]);
        carry = _addcarryx_u64(0, res[1], temp2[0], &res[1]);
        res[3] += _addcarryx_u64(carry, res[2], temp2[1], &res[2]);

        // extract the bits from res[] keeping the precision the same as this object
        // shift result by F
        constexpr int32_t index = 1;
        //constexpr int32_t lsb = FRAC_BITS & 63;            // bit within the 64bit data pointed by res[index]
        constexpr int32_t lsb = (FRAC_BITS & 63) - 1;            // bit within the 64bit data pointed by res[index] minus 1 to improve rounding
        //constexpr uint64_t half = 1ull << (lsb - 1);       // used for rounding
        //const bool need_rounding = (res[index] & half) != 0;

        l1 = shift_right128(res[index], res[index + 1], lsb); // custom function is 20% faster in Mandelbrot than the intrinsic
        h1 = shift_right128(res[index + 1], res[index + 2], lsb);
        --expo;

        //if (need_rounding) {
        //    ++l1; // low will wrap around to zero if overflowed
        //    h1 += l1 == 0;
        //}

        norm_fraction(l1, h1, expo);
        set_components(l1, h1, expo, sign ^ other_sign);
        return *this;
    }
    /**
     * @brief Multiply a value to this object
     * @param other Right hand side operand
     * @return This object.
    */
    template<typename T>
    FP128_INLINE float128& operator*=(const T& other) {
        return operator*=(float128(other));
    }
    /**
     * @brief Divide this object by a value
     * @param other Right hand side operand
     * @return This object.
    */
    FP128_INLINE float128& operator/=(const float128& other) {
        // check trivial cases
        if (other.is_zero()) {
            *this = inf();
            return *this;
        }
        else if (is_zero()) {
            return *this;
        }
        if (is_special() || other.is_special()) {
            if (is_nan() || other.is_nan())
                *this = nan();
            // return inf with the right sign
            if (other.is_inf()) {
                *this = other;
                set_sign(get_sign() ^ other.get_sign());
            }
                
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
        uint64_t q[4]{};
        const uint64_t nom[4] = { 0, 0, l1, h1 };
        const uint64_t denom[2] = { l2, h2 };

        if (0 == div_32bit((uint32_t*)q, nullptr, (uint32_t*)nom, (uint32_t*)denom, 2ll * array_length(nom), 2ll * array_length(denom))) {
            // 128 bit were added to the dividend, 112 were lost:
            // need to shift right 16 bit (128 - 112) but we don't go all the way so norm_fraction() 
            //  can produce accurate rounding
            l1 = shift_right128(q[0], q[1], 127 - FRAC_BITS);
            h1 = shift_right128(q[1], q[2], 127 - FRAC_BITS);
            --expo;
        }
        else { // error
            *this = inf();
            return *this;
        }

        norm_fraction(l1, h1, expo);
        set_components(l1, h1, expo, sign ^ other_sign);
        return *this;
    }
    /**
     * @brief Divide this object by a value
     * @param other Right hand side operand
     * @return This object.
    */
    template<typename T>
    FP128_INLINE float128& operator/=(const T& other) {
        return operator/=(float128(other));
    }

    //
    // unary operations
    //     
    /**
     * @brief Convert to bool
    */
    __forceinline operator bool() const noexcept {
        return high != 0 || low != 0;
    }
    /**
     * @brief Logical not (!). Opposite of operator bool.
    */
    __forceinline bool operator!() const noexcept {
        return high == 0 && low == 0;
    }
    /**
     * @brief Unary +. Returns a copy of the object.
    */
    __forceinline float128 operator+() const noexcept {
        return *this;
    }
    /**
     * @brief Unary -. Returns a copy of the object with sign inverted.
    */
    __forceinline float128 operator-() const noexcept {
        return float128(low, high ^ SIGN_MASK);
    }

    //
    // useful public functions
    //
    /**
     * @brief Returns true if the value is positive (including zero and NaN)
     * @return True when the sign is 0
    */
    __forceinline bool is_positive() const noexcept {
        return high_bits.s == 0;
    }
    /**
     * @brief Returns true if the value is negative (including zero and NaN).
     * @return True when the sign is 1
    */
    __forceinline bool is_negative() const noexcept {
        return high_bits.s == 1;
    }
    /**
     * @brief Returns true if and only if the value is ±0.
     * @return Returns true if the value is zero
    */
    __forceinline bool is_zero() const noexcept {
        return 0 == low && 0 == (high & ~SIGN_MASK);
    }
    /**
     * @brief Returns true if and only if x is zero, subnormal or normal (not infinite or NaN).
     * @return True if and only if x is zero, subnormal or normal (not infinite or NaN).
    */
    __forceinline bool is_finite() const noexcept {
        return !is_special();
    }
    /**
     * @brief Tests if the value is subnormal
     * @return True when the value is subnormal
    */
    __forceinline bool is_subnormal() const noexcept {
        return high_bits.e == 0;
    }
    /**
     * @brief Tests if the value is normal (not zero, subnormal, infinite, or NaN)
     * @return True if and only if the value is normal
    */
    __forceinline bool is_normal() const noexcept {
        return high_bits.e != 0 && high_bits.e != INF_EXP_BIASED;
    }
    /**
     * @brief Tests if this value is a NaN
     * @return True when the value is a NaN
    */
    __forceinline bool is_nan() const {
        // fraction is zero for +- INF, non-zero for NaN 
        return high_bits.e == INF_EXP_BIASED && (high_bits.f != 0 || low != 0);
    }
    /**
     * @brief Tests if this value is a signaling NaN
     * @return True if this value is a signaling NaN
    */
    __forceinline bool is_signaling() const {
        // TODO: supprot sNaN
        return false;
    }
    /**
     * @brief Tests if this value is an Infinite (negative or positive)
     * @return True when the value is an Infinite
    */
    __forceinline bool is_inf() const {
        // fraction is zero for +- INF, non-zero for NaN 
        return high_bits.e == INF_EXP_BIASED && high_bits.f == 0;
    }
    /**
     * @brief Tests if the value is an exponent of 2 (fraction part is zero)
     * @return True when the value is an exponent of 2
    */
    __forceinline bool is_exponent_of_2() const {
        // fraction is zero for +- INF, non-zero for NaN 
        return high_bits.f == 0 && low == 0;
    }
    /**
     * @brief return true when the value is either an inf or nan
     * @return true for inf and nan
    */
    __forceinline bool is_special() const {
        // fraction is zero for +- INF, non-zero for NaN 
        return high_bits.e == INF_EXP_BIASED;
    }
    /**
     * @brief Returns if the value is an integer (fraction is zero).
     * @return True when the value is an integer.
    */
    __forceinline bool is_int() const {
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
    __forceinline int32_t get_bit(uint32_t bit) const noexcept
    {
        if (bit < 64) {
            return FP128_GET_BIT(low, bit);
        }
        return FP128_GET_BIT(high, bit - 64);
    }
    /**
     * @brief Return the fraction part as a float128
    */
    FP128_INLINE float128 get_fraction() const {
        auto expo = get_exponent();
        int32_t frac_bits = static_cast<int32_t>(FRAC_BITS) - expo;
        // all the bits are fraction
        if (frac_bits > FRAC_BITS) 
            return *this;
        // the exponent is too large to hold a fraction
        if (frac_bits <= 0)
            return 0;

        uint64_t l = low, h = high_bits.f;
        if (frac_bits <= 64) {
            h = 0;
            l &= FP128_MAX_VALUE_64(frac_bits);
        }
        else {
            h &= FP128_MAX_VALUE_64(frac_bits - 64);
        }
        
        // no fraction bits are high
        if (l == 0 && h == 0) {
            return 0;
        }

        // find the msb and shift to bit 112
        int32_t msb = static_cast<int32_t>(log2(l, h));
        int32_t shift = FRAC_BITS - msb; // how many bits to shift left
        shift_left128_inplace_safe(l, h, shift);
        expo -= shift;
        float128 res(l, h, expo + EXP_BIAS, get_sign());
        return res;
    }
    /**
     * @brief Inverts the sign
    */
    __forceinline void invert_sign() noexcept
    {
        high_bits.s ^= 1;
    }
    /**
     * @brief Sets the sign
    */
    __forceinline void set_sign(uint64_t s) noexcept
    {
        high_bits.s = s;
    }
    /**
     * @brief Gets the sign
    */
    __forceinline uint32_t get_sign() const noexcept
    {
        return high_bits.s;
    }
    __forceinline float128_class_t get_class() const noexcept {
        // inf and Nan
        if (high_bits.e == INF_EXP_BIASED) {
            if (is_nan())
                // TODO: support signalling NaN
                return quietNaN;
            return (high_bits.s == 0) ? positiveInfinity : negativeInfinity;
        }
        if (is_zero()) {
            return (high_bits.s == 0) ? positiveZero : negativeZero;
        }
        if (is_subnormal()) {
            return (high_bits.s == 0) ? positiveSubnormal : negativeSubnormal;
        }

        return (high_bits.s == 0) ? positiveNormal : negativeNormal;
    }
    /**
     * @brief Returns the exponent of the object - like the base 2 exponent of a floating point
     * A value of 2.1 would return 1, values in the range [0.5,1.0) would return -1.
     * @return Exponent of the number
    */
    __forceinline int32_t get_exponent() const noexcept
    {
        return static_cast<int32_t>(high_bits.e) - EXP_BIAS;
    }
    /**
     * @brief Set the exponent
     * @param e Exponent value
    */
    __forceinline void set_exponent(int32_t e) noexcept
    {
        e += EXP_BIAS;
        assert(e >= 0);
        assert(e <= INF_EXP_BIASED);
        high_bits.e =  static_cast<uint64_t>(e);
    }
    /**
     * @brief break the float into its components.
     * Normalizes subnormal values
     * @param l Reference to receive the low fraction
     * @param h Reference to receive the high fraction
     * @param e Reference to receive the unbiased exponent
     * @param s Reference to receive the sign
    */
    __forceinline void get_components(uint64_t& l, uint64_t& h, int32_t& e, uint32_t& s) const noexcept
    {
        l = low;
        h = high_bits.f;
        e = get_exponent();
        s = get_sign();
        
        if (is_subnormal()) {
            // shift the bits to the lft so the msb is on bit 112
            auto shift = FRAC_BITS - static_cast<int32_t>(log2(l, h));
            shift_left128_inplace_safe(l, h, shift);

            // update the exponent
            e += shift;
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
    __forceinline void set_components(uint64_t l, uint64_t h, int32_t e, uint32_t s) noexcept {
        // overflow 
        if (e >= INF_EXP_UNBIASED) {
            e = INF_EXP_UNBIASED;
        }
        // sub normals
        if (e <= SUBNORM_EXP_UNBIASED) {
            // fix the fraction, remove the bits 112+ and shift to the right
            int32_t shift = SUBNORM_EXP_UNBIASED + 1 - e;
            if (shift >= FRAC_BITS) {
                l = h = 0;
            }
            // some bits stay in the low or high QWORD
            else 
                shift_right128_inplace_safe(l, h, shift);
            // override the exponent to mark the value as subnormal
            e = SUBNORM_EXP_UNBIASED;
        }

        low = l;
        high_bits.f = h;
        e += EXP_BIAS;
        high_bits.e = static_cast<uint64_t>(e);
        high_bits.s = s != 0;
    }
    /**
     * @brief Normalize the fraction so the msb (unity bit) is on bit 112.
     * The fraction value must contain the unity value
     * The exponent component is adjusted based on the bit shift direction required.
     * @param l Low part of the fraction
     * @param h High part of the fraction
     * @param e Unbiased exponent, can be any value.
    */
    __forceinline void norm_fraction(uint64_t& l, uint64_t& h, int32_t& e) const noexcept {
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
            // rounding up may have happended, expect the the upper 16 bit to be exactly 1
            if ((h >> 48) != 1) {
                ++e;
            }

        }
        else {
            //assert(shift == 0);
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
    __forceinline static float128 nextUp(float128 x) {
        FP128_NOT_IMPLEMENTED_EXCEPTION;
    }
    /**
     * @brief Produces the closest value smaller than x
     * @param x Source value
     * @return Lower value closest to x
    */
    __forceinline static float128 nextDown(float128 x) {
        FP128_NOT_IMPLEMENTED_EXCEPTION;
    }

    /**
     * @brief Return the infinite constant
     * @return INF
    */
    __forceinline static float128 inf() {
        return float128(0, 0, INF_EXP_BIASED, 0);
    }
    /**
     * @brief Return the quiet (non-signaling) NaN constant
     * @return NaN
    */
    __forceinline static float128 nan() {
        return float128(1, 0, INF_EXP_BIASED, 0);
    }
    /**
     * @brief Return the value of pi
     * @return pi
    */
    __forceinline static float128 pi() {
        return float128(0x8469898CC51701B8, 0x921FB54442D1, 0x4000, 0);
    }
    /**
     * @brief Return the value of pi / 2
     * @return pi / 2
    */
    __forceinline static float128 half_pi() {
        return pi() >> 1;
    }
    /**
     * @brief Return the value of pi / 4
     * @return pi / 4
    */
    __forceinline static float128 quarter_pi() {
        return pi() >> 2;
    }
    /**
     * @brief Return the value of e
     * @return e
    */
    __forceinline static float128 e() {
        static const float128 e = "2.71828182845904523536028747135266249775724709369"; // 50 first digits of e
        return e;
    }
    /**
     * @brief Returns a value of sqrt(2)
     * @return
    */
    __forceinline static float128 sqrt_2() noexcept {
        static const float128 sqrt_2 = "1.41421356237309504880168872420969807856967187537"; // 50 first digits of sqrt(2)
        return sqrt_2;
    }
    /**
     * @brief  Returns a value of 1
     * @return 1
    */
    __forceinline static float128 one() noexcept {
        static const float128 one = 1;
        return one;
    }
    /**
     * @brief  Returns a value of 0.5
     * @return 0.5
    */
    __forceinline static float128 half() noexcept
    {
        static const float128 half = 0.5;
        return half;
    }
    /**
     * @brief Return 0.1 using maximum precision
     * @return 
    */
    __forceinline static float128 tenth() noexcept
    {
        // 0.1 using maximum precision
        return float128(0x999999999999999A, 0x999999999999, EXP_BIAS - 4, 0u);
    }
    /**
     * @brief calculates 10^e
     * @param e integer exponent, in the range
     * @return 10^e
    */
    FP128_INLINE static float128 exp10(int32_t e) noexcept
    {
        // check the limits first
        if (e < -4965) {
            return 0;
        } 
        else if (e > 4932) {
            return inf();
        }

        // calculate the exponent optimally
        float128 res = 1;
        float128 b;
        if (e > 0)
            b = 10; // 10^1
        else if (e < 0) {
            b = tenth(); 
            e = -e;
        }

        while (e > 0) {
            if (e & 1)
                res *= b;
            e >>= 1;
            b *= b;
        }
        return res;
    }

    static constexpr bool is754version1985(void) { return false; }
    static constexpr bool is754version2008(void) { return true; }
    //
    // End of class method implementation
    //

    //
    // Binary math operators
    //
    /**
     * @brief Adds 2 values and returns the result.
     * @param lhs left hand side operand
     * @param rhs Right hand side operand
     * @return Result of the operation
    */
    template<typename T>
    friend __forceinline float128 operator+(float128 lhs, const T& rhs) noexcept {
        return lhs += rhs;
    }
    /**
     * @brief subtracts the right hand side operand to this object to and returns the result.
     * @param lhs left hand side operand
     * @param rhs Right hand side operand
     * @return The float128 result
    */
    template<typename T>
    friend __forceinline float128 operator-(float128 lhs, const T& rhs) noexcept {
        return lhs -= rhs;
    }
    /**
     * @brief Multiplies the right hand side operand with this object to and returns the result.
     * @param lhs left hand side operand
     * @param rhs Right hand side operand
     * @return The float128 result
    */
    template<typename T>
    friend __forceinline float128 operator*(float128 lhs, const T& rhs) noexcept {
        return lhs *= rhs;
    }
    /**
     * @brief Divides this object by the right hand side operand and returns the result.
     * @param lhs left hand side operand
     * @param rhs Right hand side operand
     * @return The float128 result
    */
    template<typename T>
    friend __forceinline float128 operator/(float128 lhs, const T& rhs) {
        return lhs /= rhs;
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
    friend __forceinline bool operator==(const float128& lhs, const float128& rhs) noexcept {
        return lhs.high == rhs.high && lhs.low == rhs.low;
    }
    template<typename T>
    friend __forceinline bool operator==(const float128& lhs, const T& rhs) noexcept {
        return lhs == float128(rhs);
    }
    template<typename T>
    friend __forceinline bool operator==(const T& lhs, const float128& rhs) noexcept {
        return rhs == float128(lhs);
    }
    /**
     * @brief Return true when objects are not equal. Can be used as logical XOR.
     * @param lhs left hand side operand
     * @param rhs Right hand side operand
     * @return True if not equal.
    */
    friend __forceinline bool operator!=(const float128& lhs, const float128& rhs) noexcept {
        return lhs.high != rhs.high || lhs.low != rhs.low;
    }
    template<typename T>
    friend __forceinline bool operator!=(const float128& lhs, const T& rhs) noexcept {
        return lhs != float128(rhs);
    }
    template<typename T>
    friend __forceinline bool operator!=(const T& lhs, const float128& rhs) noexcept {
        return rhs != float128(lhs);
    }
    /**
     * @brief Return true if this object is small than the other
     * @param lhs left hand side operand
     * @param rhs Right hand side operand
     * @return True when this object is smaller.
    */
    friend __forceinline bool operator<(const float128& lhs, const float128& rhs) noexcept {
        auto rhs_sign = rhs.get_sign();
        auto lhs_sign = lhs.get_sign();

        // signs are different
        if (lhs_sign != rhs_sign)
            return lhs_sign > rhs_sign; // true when lhs_sign is 1 and rhs.sign is 0

        // MSB is the same, check the LSB, implies the exponent is identical
        if (lhs.high == rhs.high)
            return (lhs_sign) ? lhs.low > rhs.low : lhs.low < rhs.low;

        return (lhs_sign) ? lhs.high > rhs.high : lhs.high < rhs.high;
    }
    template<typename T>
    friend __forceinline bool operator<(const float128& lhs, const T& rhs) noexcept {
        return lhs < float128(rhs);
    }
    template<typename T>
    friend __forceinline bool operator<(const T& lhs, const float128& rhs) noexcept {
        return float128(lhs) < rhs;
    }
    /**
     * @brief Return true this object is small or equal than the other
     * @param lhs left hand side operand
     * @param rhs Right hand side operand
     * @return True when this object is smaller or equal.
    */
    friend __forceinline bool operator<=(const float128& lhs, const float128& rhs) noexcept {
        return !(lhs > rhs);
    }
    template<typename T>
    friend __forceinline bool operator<=(const float128& lhs, const T& rhs) noexcept {
        return !(lhs > float128(rhs));
    }
    template<typename T>
    friend __forceinline bool operator<=(const T& lhs, const float128& rhs) noexcept {
        return !(float128(lhs) > rhs);
    }
    /**
     * @brief Return true this object is larger than the other
     * @param lhs left hand side operand
     * @param rhs Right hand side operand
     * @return True when this object is larger.
    */
    friend __forceinline bool operator>(const float128& lhs, const float128& rhs) noexcept {
        auto rhs_sign = rhs.get_sign();
        auto lhs_sign = lhs.get_sign();

        // signs are different
        if (lhs_sign != rhs_sign)
            return lhs_sign < rhs_sign; // true when lhs_sign is 1 and rhs.sign is 0

        // MSB is the same, check the LSB, implies the exponent is identical
        if (lhs.high == rhs.high)
            return (lhs_sign) ? lhs.low < rhs.low : lhs.low > rhs.low;

        return (lhs_sign) ? lhs.high < rhs.high : lhs.high > rhs.high;
    }
    template<typename T>
    friend __forceinline bool operator>(const float128& lhs, const T& rhs) noexcept {
        return lhs > float128(rhs);
    }
    template<typename T>
    friend __forceinline bool operator>(const T& lhs, const float128& rhs) noexcept {
        return float128(lhs) > rhs;
    }
    /**
     * @brief Return true this object is larger or equal than the other
     * @param lhs left hand side operand
     * @param rhs Right hand side operand
     * @return True when this objext is larger or equal.
    */
    friend __forceinline bool operator>=(const float128& lhs, const float128& rhs) noexcept {
        return !(lhs < rhs);
    }
    template<typename T>
    friend __forceinline bool operator>=(const float128& lhs, const T& rhs) noexcept {
        return !(lhs < float128(rhs));
    }
    template<typename T>
    friend __forceinline bool operator>=(const T& lhs, const float128& rhs) noexcept {
        return !(float128(lhs) < rhs);
    }

    /**
     * @brief Return the NaN constant
     * @param  
     * @return 
    */
    friend float128 nan(const float128&) {
        return float128::nan();
    }
    /**
     * @brief Tests if the value is a NaN 
     * @param x Value to test
     * @return True when the value is a NaN
    */
    friend bool isnan(const float128& x) {
        // zero for +- INF, non-zero for NaN 
        return x.is_nan();
    }
    /**
     * @brief Tests if the value is an Infinite (negative or positive)
     * @param x Value to test
     * @return True when the value is an Infinite
    */
    friend bool isinf(const float128& x) {
        return x.is_inf();
    }

    friend __forceinline float128 fabs(const float128& x) noexcept {
        float128 temp = x;
        temp.set_sign(0);
        return temp;
    }
    /**
     * @brief Performs the floor() function, similar to libc's floor(), rounds down towards -infinity.
     * @param x Input value
     * @return A float128 holding the integer value. Overflow is not reported.
    */
    friend __forceinline float128 floor(const float128& x) noexcept {
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
    friend __forceinline float128 ceil(const float128& x) noexcept {
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
    friend __forceinline float128 trunc(const float128& x) noexcept {
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
    friend __forceinline float128 round(const float128& x) noexcept {
        float128 h = (x.is_positive()) ? half() : -half();
        return trunc(x + h);
    }
    /**
     * @brief Retrieves an integer that represents the base-2 exponent of the specified value.
     * @param x The specified value.
     * @return Integer value, rounded towards the nearest integer.
    */
    friend __forceinline int32_t ilogb(const float128& x) noexcept {
        return x.get_exponent();
    }
    /**
     * @brief returns the value of x with the sign of y.
     * @param x The value that's returned as the magnitude of the result.
     * @param y The sign of the result.
     * @return The copysign functions return a floating-point value that combines the magnitude of x and the sign of y.
    */
    friend __forceinline float128 copysign(const float128& x, const float128& y) noexcept {
        float128 temp = x;
        temp.high_bits.s = y.high_bits.s;
        return temp;
    }
    /**
     * @brief Performs the fmod() function, similar to libc's fmod(), returns the remainder of a division x/root.
     * @param x Numerator
     * @param y Denominator
     * @return The modulo value.
    */
    friend float128 fmod(const float128& x, const float128& y) noexcept {
        // trivial case, x is zero
        if (x.is_zero())
            return x;

        if (y.is_zero())
            return copysign(inf(), y);

        // do the division in with positive numbers
        float128 x_div_y = x / y;

        // Integer result - remainder is zero. 
        // Avoid the extra computation and precision loss with the standard equation.
        if (x_div_y.is_int()) {
            return 0;
        }
        // Fraction result - remainder is non zero.
        float128 res = x - y * trunc(x_div_y);
        return res;
    }
    /**
     * @brief Split into integer and fraction parts.
     * Both results carry the sign of the input variable.
     * @param x Input value
     * @param iptr Pointer to float128 holding the integer part of x.
     * @return The fraction part of x. Undefined when iptr is nullptr.
    */
    friend float128 modf(const float128& x, float128* iptr) noexcept {
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
    friend __forceinline float128 fdim(const float128& x, const float128& y) noexcept {
        return (x > y) ? x - y : float128();
    }
    /**
     * @brief Returns the mimimun between x and y.
     * @param x First value
     * @param y Second value
     * @return If x < y returns x. Otherwise y.
    */
    friend __forceinline float128 fmin(const float128& x, const float128& y) noexcept {
        return (x < y) ? x : y;
    }
    /**
     * @brief Returns the maximum between x and y.
     * @param x First value
     * @param y Second value
     * @return If x > y returns x. Otherwise y.
    */
    friend __forceinline float128 fmax(const float128& x, const float128& y) noexcept {
        return (x > y) ? x : y;
    }
    /**
     * @brief Calculates the hypotenuse. i.e. sqrt(x^2 + y^2)
     * @param x First value
     * @param y Second value
     * @return sqrt(x^2 + y^2).
    */
    friend FP128_INLINE float128 hypot(const float128& x, const float128& y) noexcept
    {
        return sqrt(x * x + y * y);
    }
    /**
     * @brief Calculates the square root using Newton's method.
     * Based on the book "Math toolkit for real time programming" by Jack W. Crenshaw
     * @param x Value to calculate the root of
     * @param iterations how many iterations to perform (more is more accurate). Sensible values are 0-5.
     * @return Square root of (x), zero when x <= 0.
    */
    friend float128 sqrt(const float128& x, uint32_t iterations) noexcept {
        static const float128 factor = "0.70710678118654752440084436210484903928483593768847403658833981"; // sqrt(2) / 2
        if (x.is_negative())
            return float128::nan();
        if (x.is_zero())
            return 0;
        if (x.is_special())
            return x;

        // normalize the input to the range [0.5, 1)
        int32_t expo = x.get_exponent() + 1;
        float128 norm_x = x;
        norm_x.set_exponent(-1);

        // use existing HW to provide an excellent first estimate.
        // regardless of what x was, norm_x is within double's precision range
        double temp = static_cast<double>(norm_x);
        float128 root = ::sqrt(temp);

        // iterate several times via Newton's method
        //                  X
        //   Xn+1 = 0.5 * (---- + Xn )
        //                  Xn
        for (auto i = iterations; i != 0; --i) {
            root = (norm_x / root + root) >> 1;
        }

        if (expo & 1) {
            root *= factor;
            ++expo;
        }

        // Denormalize the result
        int32_t root_expo = root.get_exponent();
        root_expo += (expo > 0) ? (expo + 1) / 2 : expo / 2;
        root.set_exponent(root_expo);
        return root;
    }
    /**
     * @brief Calculates the cube root.
     * Uses the Halley's method.
     * @param x Floating point value
     * @param iterations how many Halley to perform, usually 1 is enough
     * @return cube root of x
    */
    friend FP128_INLINE float128 cbrt(const float128 x, uint32_t iterations) noexcept {
        static const float128 factor1 = "0.62996052494743659533327218014164827764034271240234"; // cbrt(2) / 2
        static const float128 factor2 = "0.79370052598409979172089379062526859343051910400391"; // cbrt(4) / 2

        if (x.is_negative())
            return float128::nan();
        if (x.is_zero())
            return 0;
        if (x.is_special())
            return x;

        // normalize the input to the range [0.5, 1)
        int32_t expo = x.get_exponent() + 1;
        float128 norm_x = x;
        norm_x.set_exponent(-1);

        // use existing HW to provide an excellent first estimate.
        // regardless of what x was, norm_x is within double's precision range
        double temp = static_cast<double>(norm_x);
        float128 root = ::cbrt(temp);

        // iterate several times via Halley's method
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
        
        // correct the result if the exponent was not a multiple of 3.
        // the offset is to have the expo always positive.
        switch ((expo + 300000) % 3) {
        case 1:
            root *= factor1;
            break;
        case 2:
            root *= factor2;
            break;
        default:
            break;
        }

        // Denormalize the result
        int32_t root_expo = root.get_exponent();
        root_expo += (expo > 0) ? (expo + 2) / 3 : expo / 3;
        root.set_exponent(root_expo);
        return root;
    }
    /**
     * @brief Calculates the reciprocal of a value. y = 1 / x
     * Using newton iterations: Yn+1 = Yn(2 - x * Yn)
     * @param x Input value
     * @return 1 / x. Returns zero on overflow or division by zero
    */
    friend FP128_INLINE float128 reciprocal(const float128& x) noexcept {
        static const float128 one = 1, two = 2;
        constexpr int max_iterations = 3;
        constexpr int debug = false;
        auto x_sign = x.get_sign();
        if (x.is_special()) return x;
        if (x.is_subnormal() || x.is_zero()) return (x_sign) ? -inf() : inf();
        
        float128 norm_x = x;
        const int32_t expo = 1 + x.get_exponent();
        norm_x.set_exponent(0);
        norm_x.set_sign(0);

        float128 y = 1.0 / static_cast<double>(norm_x);

        if (!y)
            return y;

        float128 xy, y_prev;
        // Newton iterations:
        int i = 0;
        for (; i < max_iterations && (y_prev != y); ++i) {
            y_prev = y;
            xy = norm_x * y;
            //y = y * (two - xy);
            y *= two - xy;
        }

        if constexpr (debug) {
            static int debug_max_iter = 0;
            if (i > debug_max_iter || i == max_iterations) {
                debug_max_iter = i;
                printf("reciprocal took %i iterations for %.10lf\n", i, static_cast<double>(x));
            }
        }

        y.set_exponent(-expo);
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
        static const float128 c[] = {
            "1.0", // 1 / 0!
            "1.0", // 1 / 1!
            "0.5", // 1 / 2!
            "1.6666666666666666666666666666666667e-1", // 1 / 3!
            "4.1666666666666666666666666666666667e-2", // 1 / 4!
            "8.3333333333333333333333333333333333e-3", // 1 / 5!
            "1.3888888888888889418943284326246612e-3", // 1 / 6!
            "1.9841269841269841269841269841269841e-4", // 1 / 7!
            "2.4801587301587301587301587301587302e-5", // 1 / 8!
            "2.7557319223985890652557319223985891e-6", // 1 / 9!
            "2.7557319223985890652557319223985891e-7", // 1 / 10!
            "2.5052108385441718775052108385441719e-8", // 1 / 11!
            "2.0876756987868098979210090321201432e-9", // 1 / 12!
            "1.6059043836821614599392377170154948e-10", // 1 / 13!
            "1.1470745597729724713851697978682106e-11", // 1 / 14!
            "7.6471637318198164759011319857880704e-13", // 1 / 15!
            "4.7794773323873852974382074911175440e-14", // 1 / 16!
            "2.8114572543455207631989455830103200e-15", // 1 / 17!
            "1.5619206968586226462216364350057333e-16", // 1 / 18!
            "8.2206352466243297169559812368722807e-18", // 1 / 19!
            "4.1103176233121648584779906184361404e-19", // 1 / 20!
            "1.9572941063391261230847574373505430e-20", // 1 / 21!
            "8.8967913924505732867488974425024683e-22", // 1 / 22!
            "3.8681701706306840377169119315228123e-23", // 1 / 23!
            "1.6117375710961183490487133048011718e-24", // 1 / 24!
            "6.4469502843844733961948532192046872e-26", // 1 / 25!
            "2.4795962632247974600749435458479566e-27", // 1 / 26!
            "9.1836898637955461484257168364739134e-29", // 1 / 27!
            "3.2798892370698379101520417273121119e-30", // 1 / 28!
            "1.1309962886447716931558764576938317e-31", // 1 / 29!
            "3.7699876288159056438529215256461057e-33", // 1 / 30!
            "1.2161250415535179496299746856922921e-34", // 1 / 31!
            "3.8003907548547435925936708927884130e-36", // 1 / 32!
            "1.1516335620771950280586881493298221e-37", // 1 / 33!
            "3.3871575355211618472314357333230062e-39", // 1 / 34!
            "9.6775929586318909920898163809228749e-41", // 1 / 35!
            "2.6882202662866363866916156613674652e-42", // 1 / 36!
            "7.2654601791530713153827450307228790e-44", // 1 / 37!
            "1.9119632050402819251007223765060208e-45", // 1 / 38!
            "4.9024697565135433976941599397590277e-47", // 1 / 39!
            "1.2256174391283858494235399849397569e-48", // 1 / 40!
            "2.9893108271424045107891219144872120e-50", // 1 / 41!
            "7.1174067312914393114026712249695524e-52", // 1 / 42!
            "1.6552108677421951886982956337138494e-53", // 1 / 43!
            "3.7618428812322617924961264402587486e-55", // 1 / 44!
            "8.3596508471828039833247254227972192e-57", // 1 / 45!
            "1.8173154015614791268097229179993955e-58", // 1 / 46!
            "3.8666285139605938868291976978710542e-60", // 1 / 47!
            "8.0554760707512372642274952038980296e-62", // 1 / 48!
            "1.6439747083165790335158153477342917e-63", // 1 / 49!
            "3.2879494166331580670316306954685835e-65", // 1 / 50!
        };
        constexpr int series_len = array_length(c);
        static_assert(series_len == 51);

        if (x >= 0 && x < series_len) {
            res = c[x];
        }
        else {
            res = 0;
        }
    }
    /**
     * @brief returns the double factorial of a number
     * @param x 
     * @param res 
    */
    friend FP128_INLINE float128 double_factorial(int x) noexcept {
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
    friend FP128_INLINE float128 exp(const float128& x) noexcept {
        static const float128 e = float128::e();
        static const float128 max_exponent = 11355; // log(16382) / log2

        // check if the value isn't too large
        if (x > max_exponent)
            return inf();
        // check if the value isn't too small
        if (x <  -max_exponent)
            return 0;

        float128 _ix, exp_ix; // integer part of x
        float128 fx = modf(fabs(x), &_ix);
        uint64_t ix = static_cast<uint64_t>(_ix); // 64 bit is an overkill to hold the exponent
        float128 res;

        // compute e^ix (integer part of x)
        if (ix > 0) {
            exp_ix = 1;     // result
            float128 b = e; // value of e^1
            while (ix > 0) {
                if (ix & 1)
                    exp_ix *= b;
                ix >>= 1;
                b *= b;
            }
        }
        else {
            exp_ix = 1;
        }

        float128 exp_fx;
        if (!fx.is_zero()) {
            // compute e^fx (fraction part of x)
            // first and second elements of the series
            if (!fx.is_zero()) {
                exp_fx = float128::one() + fx;
                float128 elem, elem_denom, elem_nom = fx;

                for (int i = 2; ; ++i) {
                    elem_nom *= fx;
                    fact_reciprocal(i, elem_denom);
                    elem = elem_nom * elem_denom;
                    // value is too small to add any bits as the result's exponent is either 0 or 1. exp_fx = [1, e)
                    if (elem.get_exponent() <= -FRAC_BITS)
                        break;
                    exp_fx += elem; // next element in the series
                }
            }

            res = exp_ix * exp_fx;
        }
        else {
            res = exp_ix;
        }

        return (x.is_positive()) ? res : reciprocal(res);
    }
    /**
     * @brief Computes 2 to the power of x
     * @param x Exponent value
     * @return 2^x
    */
    friend FP128_INLINE float128 exp2(const float128& x) noexcept {
        //
        // Based on exponent law: (x^n)^m = x^(m*n)
        // Convert the exponent x (function parameter) to produce an exponent that will work with exp()
        // y = log(2) 
        // 2^x = e^(y*x) = exp(y*x)
        //
        static const float128 lan2 = "0.693147180559945309417232121458176575";
        return exp(x * lan2);
    }
    /**
     * @brief Calculates the exponent of x and reduces 1 from the result: (e^x) - 1
     * @param x A number specifying a power.
     * @return Exponent of x
    */
    friend FP128_INLINE float128 expm1(const float128& x) noexcept {
        return exp(x) - float128::one();
    }
    /**
     * @brief Computes x to the power of y
     * @param x Base value
     * @param y Exponent value (integer)
     * @return x^y
    */
    friend FP128_INLINE float128 pow(const float128& x, int32_t y) noexcept {
        static const float128 max_exponent = 11355; // log(16382) / log2
        float128 res = 1;
        // check the trivial cases
        if (y == 1) {
            return x;
        }
        else if (y == 0) {
            return 1;
        }
        else if (x == 1) {
            return x;
        }
        auto expo = abs(y);
        // check if the value isn't too large
        if (y > max_exponent) {
            res = inf();
        }
        // check if the value isn't too small
        else if (y < -max_exponent)
            res = 0;
        // compute x^y 
        else if (expo > 0) {
            float128 b = x; // value of e^1
            while (expo > 0) {
                if (expo & 1)
                    res *= b;
                expo >>= 1;
                b *= b;
            }
        }

        return (y >= 0) ? res : reciprocal(res);
    }
    /**
     * @brief Computes x to the power of y
     * @param x Base value
     * @param y Exponent value
     * @param f Optional: how many fraction bits in the result. Default to all.
     * @return x^y
    */
    friend FP128_INLINE float128 pow(const float128& x, const float128& y, int32_t f) noexcept
    {
        //
        // Based on exponent law: (x^n)^m = x^(m * n)
        // Convert the exponent y (function parameter) to produce an exponent that will work with exp()
        // z = log(x) 
        // pow(x, y) = x^y = e^(y * z) = exp(y * z)
        //
        if (y.is_int()) {
            return pow(x, static_cast<int32_t>(y));
        }
        else if (x.is_negative()) {
            return -nan();
        }

        float128 lan_x = log(x, f);
        if (!lan_x)
            return lan_x;

        return exp(y * lan_x);
    }
    /**
     * @brief Calculates the natural Log (base e) of x: log(x)
     * @param x The number to perform log on.
     * @param f Optional: how many fraction bits in the result. Default to all.
     * @return log(x)
    */
    friend FP128_INLINE float128 log(float128 x, int32_t f) noexcept {
        static const float128 lan2 = "0.693147180559945309417232121458176575";
        float128 y = log2(x, f);
        return y * lan2;
    }
    /**
     * @brief Calculates the Log base 2 of x: y = log2(x)
     * @param x The number to perform log2 on.
     * @param f Optional: how many fraction bits in the result. Default to all.
     * @return log2(x)
    */
    friend FP128_INLINE float128 log2(float128 x, int32_t f) noexcept {
        if (x.is_negative() || x.is_zero()) {
            return -inf();
        }

        // Calculate the log in 2 steps:
        // - The integer part (iy) is simple and fast via the get_exponent() function.
        // - The fraction part (fy) is trickier. Uses Binary Logarithm
        // The result is the sum of the two. Based on the identity:
        // log(x + y) = log(x) + log(y)

        // bring x to the range [1,2)
        auto expo = x.get_exponent();
        float128 iy = expo; // integer part of the result
        // x is an exponent of 2.
        if (x.is_exponent_of_2())
            return iy;

        x.set_exponent(0);

        static const float128 two(2);
        float128 b = float128::half(); // 0.5
        float128 fy; // fraction part of the result
        for (size_t i = 0; i < f; ++i) {
            // x = x * x
            x *= x;
            // if x is greater than 2, we have another bit in the result
            if (x  >= two) {
                // divide x by 2 using shifts
                x >>= 1;
                fy += b;
            }
            // divide 2 using shifts
            b >>= 1;
        }

        return iy + fy;
    }
    /**
     * @brief Calculates Log base 10 of x: log10(x)
     * @param x The number to perform log on.
     * @param f Optional: how many fraction bits in the result. Default to all.
     * @return log10(x)
    */
    friend FP128_INLINE float128 log10(float128 x, int32_t f) noexcept {
        static const float128 log10_2 = "0.301029995663981195213738894724493068";  // log10(2)
        float128 y = log2(x, f);
        return y * log10_2;
    }
    /**
     * @brief Calculates Log base 2 of x as an integer ignoring the sign of x.
     * Similar to: floor(log2(fabs(x)))
     * @param x The number to perform log on.
     * @return logb(x)
    */
    friend FP128_INLINE float128 logb(float128 x, int32_t) noexcept {
        return x.get_exponent();
    }
    /**
     * @brief Calculates the natural Log (base e) of 1 + x: log(1 + x)
     * @param x The number to perform log on.
     * @param f Optional: how many fraction bits in the result. Default to all.
     * @return log1p(x)
    */
    friend FP128_INLINE float128 log1p(float128 x, int32_t f) noexcept {
        return log(float128::one() + x, f);
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
    friend FP128_INLINE float128 sin1(float128 x) noexcept
    {
        assert(fabs(x) <= float128::half_pi());

        // first part of the series is just 'x'
        const float128 xx = x * x;
        float128 elem_denom, elem_nom = x;

        // compute the rest of the series, starting with: -(x^3 / 3!)
        for (int i = 3, sign = 1; ; i += 2, sign = 1 - sign) {
            elem_nom *= xx;
            fact_reciprocal(i, elem_denom);
            float128 elem = elem_nom * elem_denom; // next element in the series
            // precision limit has been hit
            if (!elem)
                break;
            x += (sign) ? -elem : elem;
        }

        return x;
    }
    /**
         * @brief Calculate the cosine function over a limited range [-0.5pi, 0.5pi]
         * Since the sin1 function converges faster, call it with the modifed angle.
         * @param x value in Radians in the range [-0.5pi, 0.5pi]
         * @return Cosine of x
        */
    friend __forceinline float128 cos1(const float128& x) noexcept
    {
        const float128 half_pi = float128::half_pi();
        assert(fabs(x) <= half_pi);
        return (x.is_positive()) ?
            sin1(half_pi - x) :
            -sin1(-half_pi - x);
    }
    /**
     * @brief Calculate the Sine function
     * Ultimately uses sin() with a reduced range of [-pi/4, pi/4]
     * @param x value in Radians
     * @return Sine of x
    */
    friend float128 sin(float128 x) noexcept
    {
        const float128 half_pi = float128::half_pi();
        double round = (x.is_positive()) ? 0.5 : -0.5;

        int64_t n = static_cast<int64_t>((x / half_pi) + round);
        x -= half_pi * n;
        n = n & 3ull; // n mod 4
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
    friend float128 asin(float128 x) noexcept
    {
        constexpr int max_iterations = 6;
        if (x < -1 || x > 1) return 0;

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
    friend float128 cos(float128 x) noexcept
    {
        const float128 half_pi = float128::half_pi(); // pi / 2
        double round = (x.is_positive()) ? 0.5 : -0.5;

        int64_t n = static_cast<int64_t>((x / half_pi) + round);
        x -= half_pi * n;
        n = n & 3ull; // n mod 4
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
    friend float128 acos(float128 x) noexcept
    {
        constexpr int max_iterations = 6;
        if (x < -1 || x > 1) return 0;
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
    friend FP128_INLINE float128 tan(float128 x) noexcept
    {
        return sin(x) / cos(x);
    }
    /**
     * @brief Calculate the inverse tangent function
     * @param x value
     * @return Arctangent of x
    */
    friend float128 atan(float128 x) noexcept
    {
        // constants for segmentation
        const float128 half_pi = float128::half_pi(); // pi / 2
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
            float128 e = cos_xn * (sin_xn - x * cos_xn); // this is the iteration estimated error
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
    friend float128 atan2(float128 y, float128 x) noexcept
    {
        // constants for segmentation
        const float128 pi = float128::pi();
        const float128 half_pi = float128::half_pi(); // pi / 2
        const float128 quarter_pi = float128::quarter_pi(); // pi / 4

        // x == 0
        if (!x) {
            if (!y) return 0;

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
        const float128 eps = fabs(res >> 110);

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
    friend FP128_INLINE float128 sinh(const float128 x) noexcept
    {
        return (exp(x) - exp(-x)) >> 1;
    }
    /**
     * @brief Calculates the inverse hyperbolic sine
     * For positive x:
     * asinh(x) = log(x + sqrt(x^2 + 1))
     * For negative x, the function returns the result with the sign inverted
     * @param x value
     * @return Inverse hyperbolic sine of x
    */
    friend FP128_INLINE float128 asinh(const float128 x) noexcept
    {
        float128 absx = fabs(x);
        float128 res = log(absx + sqrt(absx * absx + float128::one()));

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
    friend FP128_INLINE float128 cosh(const float128 x) noexcept
    {
        return (exp(x) + exp(-x)) >> 1;
    }
    /**
     * @brief Calculates the inverse hyperbolic cosine
     * For x >= 1:
     * acosh(x) = log(x + sqrt(x^2 - 1))
     * For x < 1, the function return zero
     * @param x value in the range [1, inf]
     * @return Inverse hyperbolic cosine of x
    */
    friend FP128_INLINE float128 acosh(const float128 x) noexcept
    {
        if (x < 1) return 0;

        float128 res = log(x + sqrt(x * x - 1));
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
    friend FP128_INLINE float128 tanh(const float128 x) noexcept
    {
        float128 ex = exp(x); // e^x
        float128 exm1 = exp(-x); // e^(-x)
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
    friend FP128_INLINE float128 atanh(const float128 x) noexcept
    {
        auto one = float128::one();
        if (fabs(x) >= 1)
            return 0;

        return log((one + x) / (one - x)) >> 1;
    }
    /**
     * @brief Calculates the Maclaurin series constants for the erf function.
     * The array will hold  1 / (n! * (2n + 1))
     * @param a pointer to array that receives the results. The array must be preallocated.
     * @param count Element count in the array
    */
    friend void erf_constants(float128* a, int32_t count) {
        if (a == nullptr) return;
        
        a[0] = 1;
        float128 f = 1; // value of 0!
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
    friend FP128_INLINE float128 erf(float128 x) noexcept {
        static const float128 sqrt_pi = sqrt(float128::pi());
        static const float128 two_by_sqrt_pi = float128(2) / sqrt_pi;
        constexpr int32_t C_len = 30;
        static float128 C[C_len];
        if (C[0].is_zero()) {
            erf_constants(C, C_len);
        }
        if (x.is_zero())
            return 0;
        if (x.is_inf())
            return (x.is_positive()) ? 1 : -1;
        if (x.is_nan())
            return nan;

        auto x_sign = x.get_sign();
        x.set_sign(0);

        int32_t expo = x.get_exponent();
        float128 xx = x * x;

        if (x < 1) {
            // Maclaurin series:
            //                             3      5      7      9
            //             2              x      x      x      x
            // erf(x) = -------- * ( x - ---  + ---- - ---- + ----- - ... )
            //          sqrt(pi)          3      10     42     216
            //
            //
            // each element in the series is:
            //      n    2n+1
            //  (-1)  * x
            //  -------------
            //  n! * (2n + 1)
            //
    
            // the first series element is x
            // compute the rest of the series: 
            float128 elem_nom = x;
            for (int i = 1, sign = 1; i < C_len; ++i, sign = 1 - sign) {
                elem_nom *= xx;

                float128 elem = elem_nom * C[i]; // next element in the series
                // precision limit has been hit
                if (elem.get_exponent() + float128::FRAC_BITS < expo)
                    break;
                x += (sign) ? -elem : elem;
            }
            x *= two_by_sqrt_pi;
        }
        else {
            x = float128::one() - erfc(x);
            /*
            //                     2
            //                   -x           -3       -5       -7
            //                   e         -1  x      3x      15x
            //  erf(x) = 1 - --------- * (x  - -   + ---- -  ----- - ...) 
            //               sqrt(pi)          2      4        8
            //                                    
            //
            // where each element is 
            //      n                -(2n+1)   -n
            //  (-1)  * (2n - 1)!! * x       * 2 
            
            // the left side of the equation
            float128 left_side = exp(-xx) / sqrt_pi;

            x = reciprocal(x); // first element
            float128 elem, elem1 = x, elem2;
            xx = x * x;
            expo = x.get_exponent();
            for (int i = 1, sign = 1; i < 25; ++i, sign = 1 - sign) {
                elem1 *= xx;
                elem2 = double_factorial(2 * i - 1);
                elem = (elem1 * elem2) >> i; // next element in the series
                // precision limit has been hit
                if (elem2.is_inf() || elem.get_exponent() + float128::FRAC_BITS < expo)
                    break;
                x += (sign) ? -elem : elem;
            }

            x *= left_side;
            x = float128::one() - x;
            */
        }


        x.set_sign(x_sign);
        return x;
    }
    /**
     * @brief Computes the complementary error function of a value.
     * The erfc functions return a value in the range 0 to 2. If x is too large for erfc, the errno variable is set to ERANGE.
     * @param x A floating-point value.
     * @return The erfc functions return the complementary Gauss error function of x.
    */
    friend FP128_INLINE float128 erfc(float128 x) noexcept {
        if (x.is_zero()) return 1;
        if (x.is_inf()) return (x.get_sign()) ? 2 : 0;
        if (x.is_nan())
            return nan();

        const float128 one = float128::one();
        if (fabs(x) < one)
            return one - erf(x);

        return (::erfc(static_cast<double>(x)));

        //auto x_sign = x.get_sign();
        //x.set_sign(0);

        //float128 t = reciprocal(one + x >> 1);
        //static const float128 a[13] = {
        //   0.00000000000000000,
        //    -0.0000000000000000087,
        //    0.0000000000000151128,
        //    -0.0000000000042211594,
        //    0.0000000005376941005,
        //    -0.0000000487263438743,
        //    0.0000033500543171203,
        //    -0.0001681266736668477,
        //    0.0061765455561937135,
        //    -0.1531173381477592232,
        //    2.2723844989269186145,
        //    -18.695741264372354978,
        //    78.905158514824205072
        //};
        //
        //float128 ans = 0.0;
        //for (int i = 12; i >= 0; i--) {
        //    ans = t * ans + a[i];
        //}
        //ans *= t * exp(-x * x);
        //return (x_sign) ? float128(2) - ans : ans;
    }

};

static_assert(sizeof(float128) == sizeof(uint64_t) * 2);

} //namespace fp128

#endif // FLOAT128_H