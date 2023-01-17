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
    
    static constexpr uint64_t INF_EXPONENT = 0x7FFF;
    static constexpr uint64_t EXPONENT_BIAS = 0x3FFF;
    static constexpr uint64_t FRAC_BITS = 112;
    static constexpr uint64_t EXP_BITS = 15;
    static constexpr uint64_t FRAC_MASK = FP128_MAX_VALUE_64(FRAC_BITS - 64);
    static constexpr uint64_t FRAC_UNITY = FP128_ONE_SHIFT(FRAC_BITS - 64);
    static constexpr uint64_t EXP_MASK = INF_EXPONENT;
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
            auto shift = FRAC_BITS - msb + 1;
            if (shift < 64)
                shift_left128_inplace(low, high, FRAC_BITS - msb);
            else {
                high = low << (shift - 64);
                low = 0;
            }
            
            set_exponent(expo);
        }
        // NaN & INF
        else if (d.e == 0x7FF)
        {
            high_bits.e = INF_EXPONENT;
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
        high_bits.s = d.s;
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
        if (high_bits.e == INF_EXPONENT) {
            // zero fraction means inf, otherwise nan
            if (low == 0 && high_bits.f == 0)
                return HUGE_VAL;
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
            // TODO: handle subnormal
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
     * @brief Operator double
     */
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
        static thread_local char str[128]; // need roughly a (meaningful) decimal digit per 3.2 bits
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
        if (high_bits.e == INF_EXPONENT || other.high_bits.e == INF_EXPONENT) {
            if (isnan() || other.isnan())
                *this = nan();
            // return inf with the right sign
            if (other.isinf())
                *this = other;
            return *this;
        }
        //if (is_zero()) {
        //    *this = other;
        //    return *this;
        //}
        //else if (other.is_zero())
        //    return *this;

        int32_t expo = get_exponent();
        int32_t other_expo = other.get_exponent();
        uint64_t l1 = low;
        uint64_t h1 = high_bits.f | FRAC_UNITY;
        uint64_t l2 = other.low;
        uint64_t h2 = other.high_bits.f | FRAC_UNITY;

        if (expo >= other_expo) {
            // exponents are too far apart, result will stay the same
            if (expo >= other_expo + FRAC_BITS)
                return *this;
            shift_right128_inplace(l2, h2 , expo - other_expo);
        }
        else {
            // exponents are too far apart, use the other value
            if (other_expo >= expo + FRAC_BITS) {
                *this = other;
                return *this;
            }
            shift_right128_inplace(l1, h1, other_expo - expo);
            // result base exponent comes from the other value
            expo = other_expo;
        }

        // same sign: the simple case
        if (other.high_bits.s == high_bits.s) {
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
            if (h1 & FP128_ONE_SHIFT(63)) {
                twos_complement128(l1, h1);
                invert_sign();
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
        int32_t shift = msb - FRAC_BITS;

        expo += shift;
        if (shift >= 0) {
            shift_right128_inplace(l1, h1, shift);
        }
        else {
            shift_left128_inplace(l1, h1, -shift);
        }
        low = l1;
        high_bits.f = h1;
        set_exponent(expo);
        return *this;
    }
    /**
     * @brief Add a value to this object
     * @param other Right hand side operand
     * @return This object.
    */
    template<typename T>
    FP128_INLINE float128& operator+=(const T& other) {
        return operator+=(fixed_point128(other));
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
        return operator+=(-fixed_point128(other));
    }
    /**
     * @brief Multiply a value to this object
     * @param other Right hand side operand
     * @return This object.
    */
    FP128_INLINE float128& operator*=(const float128& other) noexcept {
        FP128_NOT_IMPLEMENTED_EXCEPTION;
        return *this;
    }
    /**
     * @brief Multiply a value to this object
     * @param other Right hand side operand
     * @return This object.
    */
    template<typename T>
    FP128_INLINE float128& operator*=(const T& other) {
        return operator*=(fixed_point128(other));
    }
    /**
     * @brief Divide this object by a value
     * @param other Right hand side operand
     * @return This object.
    */
    FP128_INLINE float128& operator/=(const float128& other) noexcept {
        FP128_NOT_IMPLEMENTED_EXCEPTION;
        return *this;
    }
    /**
     * @brief Divide this object by a value
     * @param other Right hand side operand
     * @return This object.
    */
    template<typename T>
    FP128_INLINE float128& operator/=(const T& other) {
        return operator/=(fixed_point128(other));
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
     * @brief Returns the exponent of the object - like the base 2 exponent of a floating point
     * A value of 2.1 would return 1, values in the range [0.5,1.0) would return -1.
     * @return Exponent of the number
    */
    __forceinline int32_t get_exponent() const noexcept
    {
        return static_cast<int32_t>(high_bits.e) - EXPONENT_BIAS;
    }
    __forceinline void set_exponent(int32_t e) noexcept
    {
        assert(e + EXPONENT_BIAS > 0);
        assert(e <= static_cast<int32_t>(INF_EXPONENT));
        high_bits.e =  e + EXPONENT_BIAS;
    }
    /**
     * @brief Tests if this value is a NaN
     * @return True when the value is a NaN
    */
    __forceinline bool isnan() const {
        // fraction is zero for +- INF, non-zero for NaN 
        return high_bits.e == INF_EXPONENT && high_bits.f != 0;
    }
    /**
     * @brief Tests if this value is an Infinite (negative or positive)
     * @return True when the value is an Infinite
    */
    __forceinline bool isinf() const {
        // fraction is zero for +- INF, non-zero for NaN 
        return high_bits.e == INF_EXPONENT && high_bits.f == 0;
    }


    /**
     * @brief Return the infinite constant
     * @return INF
    */
    static float128 inf() {
        static const float128 _inf(0, 0, INF_EXPONENT, 0);
        return _inf;
    }
    /**
     * @brief Return the quiet (non-signaling) NaN constant
     * @return NaN
    */
    static float128 nan() {
        static const float128 _nan(1, 0, INF_EXPONENT, 0);
        return _nan;
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
        return x.isnan();
    }
    /**
     * @brief Tests if the value is an Infinite (negative or positive)
     * @param x Value to test
     * @return True when the value is an Infinite
    */
    friend bool isinf(const float128& x) {
        return x.isinf();
    }
};

static_assert(sizeof(float128) == sizeof(uint64_t) * 2);

} //namespace fp128

#endif // FLOAT128_H