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
            if (shift < 64)
                // Note shift cannot be zero
                shift_left128_inplace(low, high, shift);
            else {
                high = low << (shift - 64);
                low = 0;
            }
            
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
    FP128_INLINE float128(T x) noexcept {
        if constexpr (std::is_floating_point_v<T>) {
            new(this) float128(static_cast<double>(x));
            return;
        }
        uint64_t sign = 0;
        if constexpr (std::is_signed_v<T>) {
        #pragma warning(push) 
        #pragma warning(disable: 4702) // static analysis bug in VS 2022 17.4. This code _is_ reachable.
            // alway do positive multiplication
            if (x < 0) {
                x = -x;
                sign = 1;
            }
        #pragma warning(pop) 
        }

        // integers: convert to uint64 for a simpler operation.
        high = 0;
        low = static_cast<uint64_t>(x);
        if (low == 0)
            return;
        
        auto expo = log2(low); // this is the index of the msb as well
        auto shift = static_cast<int32_t>(FRAC_BITS - expo);
        if (shift < 64)
            shift_left128_inplace(low, high, shift);
        else {
            high = low << (shift - 64);
            low = 0;
        }
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

        uint32_t sign = 0;
        uint32_t base = 10;
        int32_t expo = 0;

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
            
        // check for a hex string
        if (0 == strncmp(p, "0x", 2)) {
            base = 16;
            p += 2;
        }

        // skip leading zeros
        while (*p == '0') ++p;

        int32_t int_digits = 0, frac_digits = 0, exp_digits = 0;
        bool exp_sign = 0;
        char* int_start = p;
        char* frac_start = nullptr;
        char* exp_start = nullptr;

        if (*p == '\0') {
            set_sign(sign);
            return;
        }

        // count the integer digits
        while (isdigit(*p) || (base == 16 && *p >= 'a' && *p <= 'f')) {
            ++int_digits;
            ++p;
        }
        
        // got a hex unsigned int literal
        // every digit is 4 bits, need to keep at most 112 bits after the msb.
        if (base == 16) {
            constexpr uint64_t max_digits = (112 + 4) / 4; // 29 hex digits. 28 for the fraction (112 bit) and another for the unity
            uint64_t* frac_bits = &low; // fill the internal data structure directly
            // start at the leftmost digit and iterate right

            int32_t digits_consumed = min(int_digits, max_digits);
            char* cur_digit = int_start;
            char* const end_digit = int_start + digits_consumed;

            if (int_digits == 0) {
                set_sign(sign);
                return;
            }
            // fill the internal structure starting the the top bits of high
            while (cur_digit < end_digit) {
                uint64_t d = *cur_digit;
                if (d >= '0' && d <= '9')
                    d -= '0';
                else
                    d = 10ull + d - 'a';

                ++cur_digit;
                auto index = 1 - (expo >> 6); // fill high part first
                assert(index >= 0 && index <= 1);
                auto shift = 60 - (expo & 63);
                frac_bits[index] |= d << shift;
                expo += 4;
            }
            // shift the result into position and fix the exponent
            uint64_t left_digit = high >> 60;
            assert(left_digit != 0);
            static constexpr int32_t digit_msb_lut[16] = { 0, 0, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3 };
            const auto digit_msb = digit_msb_lut[left_digit];
            expo -= 4 - digit_msb;
            expo += 4 * (int_digits - digits_consumed); // account for digits which do not fit in the fraction.
            shift_right128_inplace(low, high, 124 + digit_msb - FRAC_BITS); // move digit's 2nd  msb to bit 111
            set_exponent(expo);
            set_sign(sign);
            return;
        }

        // fraction and exponent are valid only with base 10
        if (base == 10)
        {
            // check for the optional decimal point
            if (*p == '.') {
                *p = '\0';
                ++p;
                frac_start = (isdigit(*p)) ? p : nullptr;
            }

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

            // check for the optional exponent
            if (*p == 'e') {
                *p = '\0'; // terminate the fraction string
                ++p;
                if (*p == '-') {
                    exp_sign = 1;
                    ++p;
                }
                if (*p == '+')
                    ++p;

                exp_start = (isdigit(*p)) ? p : nullptr;
            }

            // count the fraction digits if they exist
            if (exp_start) {
                while (isdigit(*p)) {
                    ++exp_digits;
                    ++p;
                }

                // convert the exponent
                expo = strtol(exp_start, nullptr, 10);
            }
        }

        // under flow
        if (exp_sign && exp_digits > 4965) {
            return;
        }

        // infinity
        if (!exp_sign && exp_digits > 4932) {
            //TODO: check overflow only when int and fraction are non zero to avoid 0E99999 become inf instead of zero
            *this == inf();
        }


        //const uint64_t int_val = std::strtoull(p, nullptr, 10);

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
        static thread_local char str[128]; // need roughly a (meaningful) decimals digit per 3.2 bits
        FP128_NOT_IMPLEMENTED_EXCEPTION;
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
            if (shift <  64)
                shift_right128_inplace(l2, h2 , shift);
            else {
                l2 = h2 >> (shift - 64);
                h2 = 0;
            }
        }
        else if (shift < 0) {
            shift = -shift;
            // exponents are too far apart, use the other value
            if (shift >= FRAC_BITS) {
                *this = other;
                return *this;
            }
            if (shift < 64)
                shift_right128_inplace(l1, h1, shift);
            else {
                l1 = h1 >> (shift - 64);
                h1 = 0;
            }
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
            shift_right128_inplace(l1, h1, shift);
        }
        else if (shift < 0) {
            shift_left128_inplace(l1, h1, -shift);
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

        // copy block #1 (lowest)
        l1 = shift_right128_round(res[index], res[index + 1], lsb); // custom function is 20% faster in Mandelbrot than the intrinsic
        h1 = shift_right128(res[index + 1], res[index + 2], lsb);

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
            FP128_FLOAT_DIVIDE_BY_ZERO_EXCEPTION;
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
            FP128_FLOAT_DIVIDE_BY_ZERO_EXCEPTION;
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
            auto shift = 127 - FRAC_BITS - static_cast<int32_t>(lzcnt128(l, h));
            if (shift >= 64) {
                h = l << (shift - 64);
                l = 0;
            }
            else if (shift > 0){
                shift_left128_inplace(l, h, shift);
            }

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
            // some bits stay in the high QWORD
            else if (shift < FRAC_BITS - 64) {
                shift_right128_inplace(l, h, shift);
            }
            else {
                l = shift_right128(l, h, shift);
                h = 0;
            }
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
        auto msb = 127 - static_cast<int32_t>(lzcnt128(l, h));
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
        else if (shift < -63) {
            h = l << (-shift);
            l = 0;
        }
        else if (shift < 0) {
            shift_left128_inplace(l, h, -shift);
        }


    }
    /**
     * @brief Return the infinite constant
     * @return INF
    */
    static float128 inf() {
        static const float128 _inf(0, 0, INF_EXP_BIASED, 0);
        return _inf;
    }
    /**
     * @brief Return the quiet (non-signaling) NaN constant
     * @return NaN
    */
    static float128 nan() {
        static const float128 _nan(1, 0, INF_EXP_BIASED, 0);
        return _nan;
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
};

static_assert(sizeof(float128) == sizeof(uint64_t) * 2);

} //namespace fp128

#endif // FLOAT128_H