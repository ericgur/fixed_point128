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
#pragma warning(disable: 4996)  // This function or variable may be unsafe.Consider using strncpy_s instead.To disable deprecation, use _CRT_SECURE_NO_WARNINGS.See online help for details.fixed_point128	C : \GitHub\fixed_point128\inc\float128.h	255


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
float128 sqrt(const float128& x, uint32_t iterations = 3) noexcept;
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
float128 log(float128 x, int32_t f = 112) noexcept;
float128 log2(float128 x, int32_t f = 112) noexcept;
float128 log10(float128 x, int32_t f = 112) noexcept;
float128 logb(float128 x, int32_t f = 112) noexcept;
float128 log1p(float128 x, int32_t f = 112) noexcept;
// non CRT function
float128 reciprocal(const float128& x) noexcept;
void fact_reciprocal(int x, float128& res) noexcept;

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

class __declspec(align(16)) float128
{
    // build time validation of template parameters
    static_assert(sizeof(void*) == 8, "float128 is supported in 64 bit builds only!");
    friend class fp128_gtest;
    
    static constexpr int32_t EXP_BITS = 15;
    static constexpr int32_t EXP_BIAS = 0x3FFF;
    static constexpr int32_t ZERO_EXP_BIASED = EXP_BIAS;
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
            constexpr int32_t digit_group = 8;
            static_assert(digit_group <= 9); // must fit in 32 bit
            int32_t i = 0;
            const float128 group_base = exp10(-digit_group);
            float128 current_base = 1;

            for (; i + digit_group <= digits; i += digit_group) {
                // convert to an int
                uint32_t val = 0;
                for (auto j = 0; j < digit_group; ++j) {
                    val *= 10;
                    val += frac_start[i + j] - '0';
                }
                
                current_base *= group_base;
                frac_part += current_base * val;
            }
            // handle the remaining digits one at a time
            for (; i < digits; ++i) {
                // convert to an int
                uint32_t val = frac_start[i] - '0';
                current_base *= tenth();
                frac_part += current_base * val;
            }

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
        auto expo = get_exponent();
        if (expo > 63) return (get_sign()) ? static_cast<uint64_t>(INT64_MIN) : UINT64_MAX;
        if (expo < 0) return 0;
        
        auto shift = static_cast<int>(FRAC_BITS) - expo;
        uint64_t res = shift_right128(low, high_bits.f, shift);
        // Add msb
        res |= FP128_ONE_SHIFT(expo);
        return res;
    }
    /**
     * @brief operator int64_t converts to a int64_t
    */
    FP128_INLINE operator int64_t() const noexcept {
        auto expo = get_exponent();
        if (expo > 62) return (get_sign()) ? INT64_MIN : INT64_MAX;
        if (expo < 0) return 0;

        auto shift = static_cast<int>(FRAC_BITS) - expo;
        int64_t res = static_cast<int64_t>(shift_right128(low, high_bits.f, shift));
        // Add msb
        res |= FP128_ONE_SHIFT(expo);
        return  (get_sign()) ? -res : res;
    }
    /**
     * @brief operator uint32_t converts to a uint32_t
    */
    FP128_INLINE operator uint32_t() const noexcept {
        auto expo = get_exponent();
        if (expo > 31) return (get_sign()) ? static_cast<uint32_t>(INT32_MIN) : UINT32_MAX;
        if (expo < 0) return 0;

        auto shift = static_cast<int>(FRAC_BITS) - expo;
        uint64_t res = high_bits.f >> shift;
        // Add msb
        res |= FP128_ONE_SHIFT(expo);
        return static_cast<uint32_t>(res);
    }
    /**
     * @brief operator int32_t converts to a int32_t
    */
    FP128_INLINE operator int32_t() const noexcept {
        auto expo = get_exponent();
        if (expo > 30) return (get_sign()) ? INT32_MIN : INT32_MAX;
        if (expo < 0) return 0;

        auto shift = static_cast<int>(FRAC_BITS) - expo;
        int32_t res = static_cast<int32_t>(high_bits.f >> shift);
        // Add msb
        res |= FP128_ONE_SHIFT(expo);
        return  (get_sign()) ? -res : res;
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

        // use uint128_t to produce the integer part
        int32_t frac_bits = (expo >= FRAC_BITS) ? 0 : min(FRAC_BITS - expo, FRAC_BITS);
        // copy all the fraction bits + the unity bit to a 128 bit integer
        uint128_t int_part(low, high_bits.f | FRAC_UNITY);
        uint128_t frac_part = int_part;
        if (expo >= 0) {
            int32_t shift = expo - FRAC_BITS;
            if (shift >= 0)
                int_part <<= shift;
            else
                int_part >>= -shift;
            strcpy(p, (char*)int_part);
            p += strlen(p);
        }
        else {
            *p++ = '0';
        }
        // there are fraction bits left to print
        if (frac_bits > 0) {
            uint32_t digits = (uint32_t)((double)frac_bits / 3.29);
            *p++ = '.';
            uint128_t mask(UINT64_MAX, 0xFFFFFFFFFFFF); // keeps lower 112 bits
            frac_part <<= (128 - frac_bits); // erase the remaining integer bits
            frac_part >>= 16;
            while (frac_part != 0 && digits-- > 0) {
                frac_part *= 10;
                int32_t digit = frac_part >> FRAC_BITS;
                frac_part &= mask;
                if (digit > 0 || frac_part != 0)
                    *p++ = static_cast<char>(digit + '0');
            }
        }

        //auto msb = 

        *p = '\0';
        return str;
    }
    /**
     * @brief converts the stored value to a string with scientific notation
     * @param str Output buffer
     * @param buff_size Output buffer size in bytes
    */
    void to_e_format(char* str, int32_t buff_size) const {
        strcpy(str, "e format not supported yet");
    }
    //
    // math operators
    //
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
        int32_t shift = expo - other_expo;

        if (shift > 0) {
            // exponents are too far apart, result will stay the same
            if (shift >= FRAC_BITS)
                return *this;
            shift_right128_inplace_safe(l2, h2, shift);
        }
        else if (shift < 0) {
            shift = -shift;
            // exponents are too far apart, use the other value
            if (shift >= FRAC_BITS) {
                *this = other;
                return *this;
            }
            shift_right128_inplace_safe(l1, h1, shift);
            // result base exponent comes from the other value
            expo = other_expo;
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

        // fix the exponent
        auto msb = 127 - static_cast<int32_t>(lzcnt128(l1, h1));
        // all zeros
        if (msb < 0) {
            low = high = 0;
            return *this;
        }

        // if the msb is exactly msb == FRAC_BITS the exponent stays the same
        shift = msb - FRAC_BITS;

        expo += shift;
        if (shift > 0) {
            shift_right128_inplace_safe(l1, h1, shift);
        }
        else if (shift < 0) {
            shift_left128_inplace_safe(l1, h1, -shift);
        }
        
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
        constexpr int32_t lsb = FRAC_BITS & 63;            // bit within the 64bit data pointed by res[index]
        constexpr uint64_t half = 1ull << (lsb - 1);       // used for rounding
        const bool need_rounding = (res[index] & half) != 0;

        l1 = shift_right128(res[index], res[index + 1], lsb); // custom function is 20% faster in Mandelbrot than the intrinsic
        h1 = shift_right128(res[index + 1], res[index + 2], lsb);

        if (need_rounding) {
            ++l1; // low will wrap around to zero if overflowed
            h1 += l1 == 0;
        }

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
            // need to shift right 16 bit (128 - 112)
            l1 = shift_right128_round(q[0], q[1], 128 - FRAC_BITS);
            h1 = __shiftright128(q[1], q[2], 128 - FRAC_BITS);
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
     * @brief Returns true if the value positive (incuding zero)
     * @return True when the the value positive
    */
    __forceinline bool is_positive() const noexcept
    {
        return high_bits.s == 0;
    }
    /**
     * @brief Returns true if the value negative (smaller than zero)
     * @return True when the the value negative
    */
    __forceinline bool is_negative() const noexcept
    {
        return high_bits.s == 1;
    }
    /**
     * @brief Returns true if the value is zero
     * @return Returns true if the value is zero
    */
    __forceinline bool is_zero() const noexcept
    {
        return 0 == low && 0 == high;
    }
    /**
     * @brief Tests if the value is subnormal
     * @return True when the value is subnormal
    */
    __forceinline bool is_subnormal() const noexcept
    {
        return high_bits.e == 0;
    }
    /**
     * @brief Tests if the value is normal
     * @return True when the value is normal. Return false for subnormal, inf and nan
    */
    __forceinline bool is_normal() const noexcept
    {
        return high_bits.e != 0 && high_bits.e != INF_EXP_BIASED;
    }
    /**
     * @brief Tests if this value is a NaN
     * @return True when the value is a NaN
    */
    __forceinline bool is_nan() const {
        // fraction is zero for +- INF, non-zero for NaN 
        return high_bits.e == INF_EXP_BIASED && high_bits.f != 0;
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
    __forceinline void set_components(uint64_t l, uint64_t h, int32_t e, uint32_t s) noexcept
    {
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
    __forceinline void norm_fraction(uint64_t& l, uint64_t& h, int32_t& e) const noexcept
    {
        // fix the exponent
        auto msb = static_cast<int32_t>(log2(l, h));
        // l and h are both zero
        if (msb < 0) {
            e = ZERO_EXP_UNBIASED;
            return;
        }

        // if the msb is exactly msb == FRAC_BITS the exponent stays the same
        auto shift = msb - FRAC_BITS;

        e += shift;
        if (shift > 0) {
            shift_right128_inplace(l, h, shift);
        }
        else
            shift_left128_inplace_safe(l, h, -shift);

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
        if (x.is_negative() || x.is_zero())
            return 0;

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
            root = (norm_x / root + root);
            root.set_exponent(root.get_exponent() - 1);
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
     * @brief Calculates the reciprocal of a value. y = 1 / x
     * Using newton iterations: Yn+1 = Yn(2 - x * Yn)
     * @param x Input value
     * @return 1 / x. Returns zero on overflow or division by zero
    */
    friend FP128_INLINE float128 reciprocal(const float128& x) noexcept {
        static const float128 one = 1, two = 2;
        constexpr int max_iterations = 3;
        constexpr int debug = false;
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
        y.set_sign(x.get_sign());
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

        if (x >= 0 && x < series_len) {
            res = c[x];
        }
        else {
            res = 0;
        }
    }
};

static_assert(sizeof(float128) == sizeof(uint64_t) * 2);

} //namespace fp128

#endif // FLOAT128_H