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

************************************************************************************/

#ifndef UINT128_T_H
#define UINT128_T_H

// override some static analysis checks
#pragma warning(push)
#pragma warning(disable: 26472) // Don't use a static_cast for arithmetic conversions. Use brace initialization
#pragma warning(disable: 26485) // No array to pointer decay
#pragma warning(disable: 26481) // Don't use pointer arithmetic. Use span instead
#pragma warning(disable: 26446) // Prefer to use gsl::at() instead of unchecked subscript operator
#pragma warning(disable: 26482) // Only index into arrays using constant expressions
#pragma warning(disable: 26408) // Avoid malloc() and free(), prefer the nothrow version of new with delete

#include "fixed_point128_shared.h"

/***********************************************************************************
*                                  Build Options
************************************************************************************/
// Set to TRUE to disable function inlining - useful for profiling a specific function
#ifndef UINT128_T_DISABLE_INLINE
#define UINT128_T_DISABLE_INLINE FALSE
#endif 

#if UINT128_T_DISABLE_INLINE != FALSE
#define UINT128_T_INLINE __declspec(noinline)
#else
#define UINT128_T_INLINE __forceinline
#endif

namespace fp128 {

/***********************************************************************************
*                                  Forward declarations
************************************************************************************/
class fp128_gtest; // Google test class

/***********************************************************************************
*                                  Main Code
************************************************************************************/

/**
 * @brief 128 bit unsigned integer class.
 * 
 * This class implements the standard operators an unsigned integer type has.<BR>
 * All of uint128_t's methods are inline for maximum performance.
 *
 * <B>Implementation notes:</B>
 * <UL>
 * <LI>Overflow is handled silently, similar to builtin integer operations.</LI>
 * <LI>A uint128_t object is not thread safe. Accessing a const object from multiple threads is safe.</LI>
 * <LI>uint128_t is <B>conditionally safe</B>, 2 different non const objects can be accessed concurrently.</LI>
 * <LI>Only 64 bit builds are supported.</LI>
 * </UL>
*/

class __declspec(align(16)) uint128_t
{
    // build time validation of template parameters
    static_assert(sizeof(void*) == 8, "uint128_t is supported in 64 bit builds only!");
    template<int32_t I> friend class fixed_point128;
    friend class fp128_gtest;

private:
    //
    // members
    //
#pragma warning(push)
#pragma warning(disable: 4201)
    union {
        struct {
            uint64_t low;  // lower QWORD
            uint64_t high; // upper QWORD
        };
        uint64_t words[2];
    };
#pragma warning(pop)
    // useful const calculations
    static constexpr uint64_t dbl_f_msb = 64ull + dbl_frac_bits; // msb location of the double precision fraction
    static constexpr uint64_t flt_f_msb = 64ull + flt_frac_bits; // msb location of the single precision fraction

public:
    typedef uint128_t type;
    typedef uint128_t* ptr_type;
    typedef uint128_t& ref_type;

    //
    // Constructors
    //

    /**
     * @brief Default constructor, creates an instance with a value of zero.
    */
    uint128_t() noexcept :
        low(0), high(0) {}
    /**
     * @brief Copy constructor
     * @param rhs Object to copy from
    */
    uint128_t(const uint128_t& rhs) noexcept :
        low(rhs.low), high(rhs.high) {}
    /**
     * @brief Move constructor
     * Doesn't modify the right hand side object. Acts like a copy constructor.
     * @param rhs Object to copy from
    */
    uint128_t(const uint128_t&& rhs) noexcept :
        low(rhs.low), high(rhs.high) {}
    /**
     * @brief Constructor from the double type
     * Underflow goes to zero. Overflow, NaN and +-INF go to max supported positive value.
     * @param x Input value
    */
    uint128_t(double x) noexcept {
        Double d;
        d.val = x;
        // very common case
        if (x == 0) {
            low = high = 0;
            return;
        }

        const int32_t e = static_cast<int32_t>(d.e) - 1023;
        uint64_t f = d.f;
        
        // overflow which catches NaN and Inf
        if (e > 127) {
            high = low = UINT64_MAX;
            return;
        }

        // normal number, produces non zero value
        if (e >= 0) {
            // bit 52 in f is the unity value of the float. it needs to move to the unity position in fixed point
            f |= FP128_ONE_SHIFT(dbl_frac_bits);
            const int32_t bits_to_shift = e - dbl_frac_bits;
            low = f;
            high = 0;

            // f fits in high QWORD
            if (bits_to_shift > 0) {
                *this <<= bits_to_shift;
            }
            // shift right
            else {
                *this >>= -bits_to_shift;
            }
        }
        // too small to be represented, no need to bother.
        else {
            high = low = 0;
        }
    }
    /**
     * @brief Constructor from uint64_t type
     * @param x Input value
    */
    uint128_t(uint64_t x) noexcept :
        low(x), high(0) {}
    /**
     * @brief Constructor from int64_t type
     * @param x Input value
    */
    uint128_t(int64_t x) noexcept {
        low = static_cast<uint64_t>(x);
        high = (x < 0) ? UINT64_MAX : 0; // sign extend the value to the higher QWORD
    }
    /**
     * @brief Constructor from uint32_t type
     * @param x Input value
   */
    uint128_t(uint32_t x) noexcept :
        low(x), high(0) {}
    /**
     * @brief Constructor from int32_t type
     * @param x Input value
   */
    uint128_t(int32_t x) noexcept {
        low = static_cast<uint64_t>(x);
        high = (x < 0) ? UINT64_MAX : 0; // sign extend the value to the higher QWORD
    }
    /**
     * @brief Constructor from const char* (C string).
     * Allows creating 128 bit values from a string. Much slower than the other constructors.<BR>
     * Input string can be decimal or hex. Hex initialization is faster.
     * Throws an std::invalid_argument when encountering an illegal character.
     * @param x Input string
    */
    uint128_t(const char* x) noexcept {
        low = high = 0;

        if (x == nullptr) return;

        // convert the input string to lowercase for simpler processing.
        char* str = _strdup(x);
        _strlwr_s(str, 1 + strlen(x));
        x = str;

        // trim leading white space
        while (*x && isspace(*x))
            ++x;
        
        // base 10 or 16 are supported
        const uint32_t base = (0 == strncmp("0x", x, 2)) ? 16u : 10u;
        if (base == 16)
            x += 2;
        
        // trim leading zeros
        while (*x == '0')
            ++x;
        
        // convert one digit at a time
        while (*x && !isspace(*x) && *x != '.') {
            uint64_t d = *x;
            if (d >= '0' && d <= '9')
                d -= '0';
            else if (base == 16 && (d >= 'a' && d <= 'f'))
                d = 10ull + d - 'a';
            else if (d == ',') // allow string to contain commas
                continue;
            else
                break; // treat an unknown char as end of string
            
            // 4 bits per digit
            if (base == 16) {
                *this <<= 4;
                low |= d;
            }
            else {
                *this *= base;
                *this += d;
            }

            ++x;
        }

        free(str);
    }
    /**
     * @brief Constructor from std::string.
     * Allows creating very high precision values. Much slower than the other constructors.
     * @param x Input string
    */
    uint128_t(const std::string& x) noexcept {
        *this = x.c_str();
    }
    /**
     * @brief Constructor from the 2 uint128_t base elements, useful for creating special constants.
     * @param l Low QWORD
     * @param h High QWORD
    */
    uint128_t(uint64_t l, uint64_t h) noexcept:
        low(l), high(h) {}    
    /**
     * @brief Destructor
    */
    ~uint128_t() noexcept {}
    /**
     * @brief Assignment operator
     * @param rhs Object to copy from 
     * @return This object.
    */
    UINT128_T_INLINE uint128_t& operator=(const uint128_t& rhs) noexcept {
        high = rhs.high;
        low = rhs.low;
        return *this;
    }
    /**
     * @brief Assignment operator
     * @param rhs Value to copy from
     * @return This object.
    */
    template<typename T>
    UINT128_T_INLINE uint128_t& operator=(const T& rhs) noexcept {
        *this = uint128_t(rhs);
        return *this;
    }
    /**
     * @brief Move assignment operator
     * @param rhs Object to copy from
     * @return This object.
    */
    UINT128_T_INLINE uint128_t& operator=(const uint128_t&& rhs) noexcept {
        high = rhs.high;
        low = rhs.low;
        return *this;
    }

    //
    // conversion operators
    //
    /**
     * @brief operator uint64_t - converts to a uint64_t
     * @return Object value.
    */
    UINT128_T_INLINE operator uint64_t() const noexcept {
        return low;
    }
    /**
     * @brief operator int64_t - converts to a int64_t
     * @return Object value.
    */
    UINT128_T_INLINE operator int64_t() const noexcept {
        return static_cast<int64_t>(low);
    }
    /**
     * @brief operator uint32_t - converts to a uint32_t
     * @return Object value.
    */
    UINT128_T_INLINE operator uint32_t() const noexcept {
        return static_cast<uint32_t>(low);
    }
    /**
     * @brief operator int32_t - converts to a int32_t
     * @return Object value.
    */
    UINT128_T_INLINE operator int32_t() const noexcept {
        return static_cast<uint32_t>(low);
    }
    /**
     * @brief operator float - converts to a float
     * @return Object value.
    */
    UINT128_T_INLINE operator float() const noexcept {
        if (!*this)
            return 0;

        const uint64_t expo = log2(*this); // returns the bit location of the msb
        Float d;
        d.e = expo + 127;

        // move bits to the high QWORD so the msb goes to bit 23. bit [22:0] will contain the fraction.
        // double actually doesn't hold the msb, it's implicit
        if (expo < flt_f_msb) {
            d.f = shift_left128(low, high, flt_f_msb - static_cast<int32_t>(expo));
        }
        else {
            d.f = high >> (expo - flt_f_msb);
        }

        return d.val;
    }
    /**
     * @brief operator double - converts to a double
     * @return Object value.
    */
    UINT128_T_INLINE operator double() const noexcept {
        if (!*this)
            return 0;
        
        const uint64_t expo = log2(*this); // returns the bit location of the msb
        Double d;
        d.e = expo + 1023;
        
        // move bits to the high QWORD so the msb goes to bit 52. bit [51:0] will contain the fraction.
        // double actually doesn't hold the msb, it's implicit
        if (expo < dbl_f_msb) {
            d.f = shift_left128(low, high, dbl_f_msb - static_cast<int32_t>(expo));
        }
        else  {
            d.f = high >> (expo - dbl_f_msb);
        }

        return d.val;
    }
    /**
     * @brief operator long double - converts to a long double
     * @return Object value.
    */
    UINT128_T_INLINE operator long double() const noexcept {
        return operator double();
    }
    /**
     * @brief Converts to a std::string (slow) string holds all meaningful fraction bits.
     * @return object string representation
    */
    UINT128_T_INLINE operator std::string() const {
        return operator char*();
    }
    /**
     * @brief Converts to a C string (slow) string holds all meaningful fraction bits.
     * The returned string is a statically, thread-allocated buffer.
     * Additional calls to this function from the same thread, overwrite the previous result.
     * @return C string with describing the value of the object.
     * @return object string representation
    */
    explicit UINT128_T_INLINE operator char*() const noexcept {
        static thread_local char str[45];
        // small number, use the fast snprintf method
        if (high == 0) {
            snprintf(str, sizeof(str), "%llu", low);
            return str;
        }

        uint128_t temp = *this;
        str[32] = 0;
        char* p = str + 31; // writing the string in reverse
        uint128_t base = 10;
        uint64_t q[2]{};
        uint64_t digit{};
        // as long as the intermediate value is >64 bit, use the more expensive long division.
        while (temp.high) {
            --p;
            const uint64_t nom[2] = { temp.low, temp.high };
            div_64bit((uint64_t*)q, &digit, (uint64_t*)nom, 10ull, 2);
            temp.low = q[0];
            temp.high = q[1];
            *p = static_cast<char>(digit + '0');
        }
        // the intermediate value fits in 64 bit, use the much faster native-64 bit division 
        while (temp.low) {
            --p;
            uint64_t r = temp.low % 10ull;
            temp.low /= 10ull;
            *p = static_cast<char>(r + '0');
        }
        return p;
    }

    //
    // math operators
    //

    /**
     * @brief Performs right shift operation.
     * @param shift bits to shift
     * @return Temporary object with the result of the operation
    */
    UINT128_T_INLINE uint128_t operator>>(int32_t shift) const noexcept {
        uint128_t temp(*this);
        return temp >>= shift;
    }
    /**
     * @brief Performs left shift operation.
     * @param shift bits to shift
     * @return Temporary object with the result of the operation
    */
    UINT128_T_INLINE uint128_t operator<<(int32_t shift) const noexcept {
        uint128_t temp(*this);
        return temp <<= shift;
    }
    /**
     * @brief Add a value to this object
     * like other uint types, overflow will result in a small value
     * @param rhs Right hand side operand
     * @return This object.
    */
    UINT128_T_INLINE uint128_t& operator+=(const uint128_t& rhs) noexcept {
        const uint8_t carry = _addcarry_u64(0, low, rhs.low, &low);
        high += rhs.high + carry;
        return *this;
    }
    /**
     * @brief Add a value to this object
     * like other uint types, overflow will result in a small value
     * @param rhs Right hand side operand
     * @return This object.
    */
    template<typename T>
    UINT128_T_INLINE uint128_t& operator+=(const T& rhs) noexcept {
        return operator+=(uint128_t(rhs));
    }
    /**
     * @brief Subtract a value from this object
     * like other uint types, underflow will result in a large value
     * @param rhs Right hand side operand
     * @return This object.
    */
    UINT128_T_INLINE uint128_t& operator-=(const uint128_t& rhs) noexcept {
        uint128_t temp = rhs;
        twos_complement128(temp.low, temp.high);
        const uint8_t carry = _addcarry_u64(0, low, temp.low, &low);
        high += temp.high + carry;
        return *this;
    }
    /**
     * @brief Subtract a value to this object
     * like other uint types, underflow will result in a large value
     * @param rhs Right hand side operand
     * @return This object.
    */
    template<typename T>
    UINT128_T_INLINE uint128_t& operator-=(const T& rhs) noexcept {
        return operator-=(uint128_t(rhs));
    }
    /**
     * @brief Multiplies a value to this object
     * @param rhs Right hand side operand
     * @return This object.
    */
    UINT128_T_INLINE uint128_t& operator*=(const uint128_t& rhs) noexcept{
        uint128_t temp = *this;

        // multiply low QWORDs
        low = _umul128(temp.low, rhs.low, &high);

        // multiply low this and high rhs
        high += temp.low * rhs.high;

        // multiply high this and low rhs
        high += temp.high * rhs.low;
        return *this;
    }
    /**
     * @brief Multiplies a value to this object
     * @param x Right hand side operand
     * @return This object.
    */
    template<typename T>
    UINT128_T_INLINE uint128_t& operator*=(T x) noexcept {
        // check if the type is signed or not
        // for negative values, convert to uint128_t and multiply.
        if constexpr (std::is_signed<T>::value) {
            if (x < 0) return operator*=(uint128_t(x));
        }

        uint64_t temp;
        const uint64_t uval = static_cast<uint64_t>(x);
        // multiply low QWORDs
        low = _umul128(low, uval, &temp);
        high = high * uval + temp;

        return *this;
    }
    /**
     * @brief Divide this object by rhs.
     * @param rhs Right hand side operator (denominator)
     * @return this object.
    */
    inline uint128_t& operator/=(const uint128_t& rhs) {
        // check some trivial cases
        if (rhs.is_zero()) {
            FP128_INT_DIVIDE_BY_ZERO_EXCEPTION;
        }
        if (is_zero() || rhs > *this) {
            low = high = 0;
            return *this;
        }
        if (rhs == *this) {
            low = 1;
            high = 0;
            return *this;
        }
        
        // exponent of 2, convert to a much faster shift operation
        if (1 == popcnt128(rhs.low, rhs.high)) {
            return *this >>= static_cast<int32_t>(log2(rhs));
        }

        uint64_t q[2]{}; 
        const uint64_t nom[2] = { low, high };

        // optimization for when dividing by a small integer
        if (rhs.high == 0) {
            if (div_64bit((uint64_t*)q, nullptr, (uint64_t*)nom, rhs.low, 2)) {
                FP128_INT_DIVIDE_BY_ZERO_EXCEPTION;
            }
        }
        // divide by a 128 bit divisor
        else {
            const uint64_t denom[2] = {rhs.low, rhs.high};
            if (div_32bit((uint32_t*)q, nullptr, (uint32_t*)nom, (uint32_t*)denom, 2ll * array_length(nom), 2ll * array_length(denom))) {
                FP128_INT_DIVIDE_BY_ZERO_EXCEPTION;
            }
        }
        low = q[0];
        high = q[1];
        return *this;
    }
    /**
     * @brief Divide this object by x.
     * @param x Right hand side operator (denominator)
     * @return this object.
    */
    template<typename T>
    inline uint128_t& operator/=(T x) {
        if (x == 0) FP128_INT_DIVIDE_BY_ZERO_EXCEPTION;

        // check if the type is signed or not
        // for negative values only, convert to uint128 and divide.
        if constexpr (std::is_signed<T>::value) {
            if (x < 0) return operator/=(uint128_t(x));
        }

        uint64_t uval = static_cast<uint64_t>(x);

        // check some trivial cases
        if (is_zero() || *this < uval) {
            low = high = 0;
            return *this;
        }

        if (*this == uval) {
            low = 1;
            high = 0;
            return *this;
        }

        // exponent of 2, convert to a much faster shift operation
        if (1 == __popcnt64(uval)) {
            return *this >>= (int32_t)log2(uval);
        }

        uint64_t nom[2] = { low, high };

        if (div_64bit(words, nullptr, (uint64_t*)nom, uval, 2)) {
            FP128_INT_DIVIDE_BY_ZERO_EXCEPTION;
        }
        return *this;
    }
    /**
     * @brief %= operator
     * @param rhs Modulo operand.
     * @return This object.
    */
    inline uint128_t& operator%=(const uint128_t& rhs) {
        // check some trivial cases
        if (rhs.is_zero()) FP128_INT_DIVIDE_BY_ZERO_EXCEPTION;

        if (*this < rhs) return *this;

        if (*this == rhs) {
            low = 0; high = 0;
        }

        uint64_t q[2]{};
        const uint64_t nom[2] = { low, high };

        // optimization for when dividing by a small integer
        if (rhs.high == 0) {
            if (div_64bit((uint64_t*)q, &low, (uint64_t*)nom, rhs.low, 2)) {
                FP128_INT_DIVIDE_BY_ZERO_EXCEPTION;
            }
            high = 0;
        }
        else {
            const uint64_t denom[2] = { rhs.low, rhs.high };
            if (div_32bit((uint32_t*)q, (uint32_t*)words, (uint32_t*)nom, (uint32_t*)denom, 2ll * array_length(nom), 2ll * array_length(denom))) {
                FP128_INT_DIVIDE_BY_ZERO_EXCEPTION;
            }
        }
        return *this;
    }
    /**
     * @brief %= operator
     * @param x Modulo operand.
     * @return This object.
    */
    template<typename T>
    inline uint128_t& operator%=(T x) {
        return operator%=(uint128_t(x));
    }
    /**
     * @brief Shift right this object.
     * Shifting is done without rounding to match uint64_t behavior
     * @param shift Bits to shift. Negative or very high values cause undefined behavior.
     * @return This object.
    */
    UINT128_T_INLINE uint128_t& operator>>=(int32_t shift) noexcept {
        if (shift < 1)
            return *this;
        // 1-64 bit shift - most common
        if (shift < 64) {
            low = __shiftright128(low, high, static_cast<uint8_t>(shift));
            high >>= shift;
        }
        else {
            low = high >> (shift - 64);
            high = 0;
        }
        return *this;
    }
    /**
     * @brief Shift left this object.
     * @param shift Bits to shift. Negative or very high values cause undefined behavior. 
     * @return This object.
    */
    UINT128_T_INLINE uint128_t& operator<<=(int32_t shift) noexcept {
        if (shift < 1)
            return *this;
        if (shift < 64) {
            high = __shiftleft128(low, high, static_cast<uint8_t>(shift));
            low <<= shift;
        }
        else {
            high = low << (shift - 64);
            low = 0;
        }
        return *this;
    }
    /**
     * @brief Bitwise AND=
     * @param rhs AND mask.
     * @return This object.
    */
    UINT128_T_INLINE uint128_t& operator&=(const uint128_t& rhs) noexcept {
        low &= rhs.low;
        high &= rhs.high;
        return *this;
    }
    /**
     * @brief Bitwise AND=
     * @param rhs AND mask.
     * @return This object.
    */
    template<typename T>
    UINT128_T_INLINE uint128_t& operator&=(const T& rhs) noexcept {
        return operator&=(uint128_t(rhs));
    }
    /**
     * @brief Bitwise OR=
     * @param rhs OR mask.
     * @return This object.
    */
    UINT128_T_INLINE uint128_t& operator|=(const uint128_t& rhs) noexcept {
        low |= rhs.low;
        high |= rhs.high;
        return *this;
    }
    /**
     * @brief Bitwise OR=
     * @param rhs OR mask.
     * @return This object.
    */
    template<typename T>
    UINT128_T_INLINE uint128_t& operator|=(const T& rhs) noexcept {
        return operator|=(uint128_t(rhs));
    }
    /**
     * @brief Bitwise XOR=
     * @param rhs XOR mask.
     * @return This object.
    */
    UINT128_T_INLINE uint128_t& operator^=(const uint128_t& rhs) noexcept {
        low ^= rhs.low;
        high ^= rhs.high;
        return *this;
    }
    /**
     * @brief Bitwise XOR=
     * @param rhs XOR mask.
     * @return This object.
    */
    template<typename T>
    UINT128_T_INLINE uint128_t& operator^=(const T& rhs) noexcept {
        return operator^=(uint128_t(rhs));
    }
    /**
     * @brief Prefix ++ operation (++a)
     * @return This object.
    */
    UINT128_T_INLINE uint128_t& operator++() noexcept {
        *this += 1;
        return *this;
    }
    /**
     * @brief Postfix ++ operation (a++)
     * @return This object.
    */
    UINT128_T_INLINE uint128_t operator++(int32_t) noexcept {
        uint128_t temp(*this);
        ++*this; // call the prefix implementation
        return temp;
    }
    /**
     * @brief Prefix -- operation (--a)
     * @return This object.
    */
    UINT128_T_INLINE uint128_t& operator--() noexcept {
        *this -= 1;
        return *this;
    }
    /**
     * @brief Postfix -- operation (a--)
     * @return This object.
    */
    UINT128_T_INLINE uint128_t operator--(int32_t) noexcept {
        uint128_t temp(*this);
        --*this; // call the prefix implementation
        return temp;
    }
    
    //
    // unary operations
    //     
    /**
     * @brief Convert to bool
    */
    UINT128_T_INLINE operator bool() const noexcept {
        return high != 0 || low != 0;
    }
    /**
     * @brief Logical not (!). Opposite of operator bool.
    */
    UINT128_T_INLINE bool operator!() const noexcept {
        return high == 0 && low == 0;
    }
    /**
     * @brief Bitwise not (~).
    */
    UINT128_T_INLINE uint128_t operator~() const noexcept {
        uint128_t temp(*this);
        temp.high = ~high;
        temp.low = ~low;
        return temp;
    }
    /**
     * @brief Unary +. Returns a copy of the object.
    */
    UINT128_T_INLINE uint128_t operator+() const noexcept {
        uint128_t temp(*this);
        return temp;
    }
    /**
     * @brief Unary -. Returns a copy of the object with sign inverted.
     * Performs a 2's complement operation just like native unsigned types
    */
    UINT128_T_INLINE uint128_t operator-() const noexcept {
        uint128_t temp = *this;
        twos_complement128(temp.low, temp.high);
        return temp;
    }

    //
    // Comparison operators
    //
    /**
     * @brief Compare logical/bitwise equal.
     * @param rhs Righthand operand
     * @return True if this and rhs are equal.
    */
    UINT128_T_INLINE bool operator==(const uint128_t& rhs) const noexcept {
        return high == rhs.high && low == rhs.low;
    }
    /**
     * @brief Compare logical/bitwise equal.
     * @param rhs Righthand operand
     * @return True if this and rhs are equal.
    */
    template<typename T>
    UINT128_T_INLINE bool operator==(const T& rhs) const noexcept {
        return *this == uint128_t(rhs);
    }
    /**
     * @brief Return true when objects are not equal. Can be used as logical XOR.
     * @param rhs Righthand operand.
     * @return True of not equal.
    */
    UINT128_T_INLINE bool operator!=(const uint128_t& rhs) const noexcept {
        return low != rhs.low || high != rhs.high;
    }
    /**
     * @brief Return true when objects are not equal. Can be used as logical XOR.
     * @param rhs Righthand operand.
     * @return True of not equal.
    */
    template<typename T>
    UINT128_T_INLINE bool operator!=(const T& rhs) const noexcept {
        return *this != uint128_t(rhs);
    }
    /**
     * @brief Return true if this object is small than the rhs
     * @param rhs Righthand operand.
     * @return True when this object is smaller.
    */
    UINT128_T_INLINE bool operator<(const uint128_t& rhs) const noexcept {
        return high < rhs.high || (high == rhs.high && low < rhs.low);
    }
    /**
     * @brief Return true if this object is small than the rhs
     * @param rhs Righthand operand.
     * @return True when this object is smaller.
    */
    template<typename T>
    UINT128_T_INLINE bool operator<(const T& rhs) const noexcept {
        return *this < uint128_t(rhs);
    }
    /**
     * @brief Return true this object is small or equal than the rhs
     * @param rhs Righthand operand.
     * @return True when this object is smaller or equal.
    */
    UINT128_T_INLINE bool operator<=(const uint128_t& rhs) const noexcept {
        return !(*this > rhs);
    }
    /**
     * @brief Return true this object is small or equal than the rhs
     * @param rhs Righthand operand.
     * @return True when this object is smaller or equal.
    */
    template<typename T>
    UINT128_T_INLINE bool operator<=(const T& rhs) const noexcept {
        return !(*this > uint128_t(rhs));
    }
    /**
     * @brief Return true this object is larger than the rhs
     * @param rhs Righthand operand.
     * @return True when this objext is larger.
    */
    UINT128_T_INLINE bool operator>(const uint128_t& rhs) const noexcept {
        return high > rhs.high || (high == rhs.high && low > rhs.low);
    }
    /**
     * @brief Return true this object is larger than the rhs
     * @param rhs Righthand operand.
     * @return True when this objext is larger.
    */
    template<typename T>
    UINT128_T_INLINE bool operator>(const T& rhs) const noexcept {
        return *this > uint128_t(rhs);
    }
    /**
     * @brief Return true this object is larger or equal than the rhs
     * @param rhs Righthand operand.
     * @return True when this objext is larger or equal.
    */
    UINT128_T_INLINE bool operator>=(const uint128_t& rhs) const noexcept {
        return !(*this < rhs);
    }
    /**
     * @brief Return true this object is larger or equal than the rhs
     * @param rhs Righthand operand.
     * @return True when this objext is larger or equal.
    */
    template<typename T>
    UINT128_T_INLINE bool operator>=(const T& rhs) const noexcept {
        return !(*this < uint128_t(rhs));
    }

    //
    // useful public functions
    //
    /**
     * @brief Returns true if the value is an int (fraction is zero)
     * @return True when the fraction is zero.
    */
    UINT128_T_INLINE constexpr bool is_int() const noexcept {
        return true;
    }
    /**
     * @brief Returns true if the value positive (incuding zero)
     * @return True when the the value positive
    */
    UINT128_T_INLINE constexpr bool is_positive() const noexcept {
        return true;
    }
    /**
     * @brief Returns true if the value negative (smaller than zero)
     * @return True when the the value negative
    */
    UINT128_T_INLINE constexpr bool is_negative() const noexcept {
        return false;
    }
    /**
     * @brief Returns true if the value is zero
     * @return Returns true if the value is zero 
    */
    UINT128_T_INLINE bool is_zero() const noexcept {
        return 0 == low && 0 == high;
    }
    /**
     * @brief get a specific bit within the 128 data
     * @param bit bit to get [0,127]
     * @return 0 or 1. Undefined when bit > 127
    */
    UINT128_T_INLINE int32_t get_bit(unsigned bit) const noexcept {
        if (bit < 64) {
            return FP128_GET_BIT(low, bit);
        }
        return FP128_GET_BIT(high, bit-64);
    }
    /**
     * @brief Return an instance of uint128_t with the value of 1
     * @return 1
    */
    UINT128_T_INLINE static const uint128_t& one() noexcept {
        static const uint128_t one = 1;
        return one;
    }
    /**
     * @brief Converts this object to a hex C string.
     * The returned string is a statically thread-allocated buffer.
     * Additional calls to this function from the same thread overwrite the previous result.
     * @return C string containing the HEX value of the object.
    */
    UINT128_T_INLINE char* hex() const noexcept {
        constexpr int buff_size = 35;
        static thread_local char str[buff_size];
        snprintf(str, buff_size, "0x%llX%016llX", high, low);
        return str;
    }

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
    friend __forceinline uint128_t operator+(uint128_t lhs, const T& rhs) noexcept {
        return lhs += rhs;
    }
    /**
     * @brief Subtracts 2 values and returns the result.
     * @param lhs left hand side operand
     * @param rhs Right hand side operand
     * @return Result of the operation
    */
    template<typename T>
    friend __forceinline uint128_t operator-(uint128_t lhs, const T& rhs) noexcept {
        return lhs -= rhs;
    }
    /**
     * @brief Multiplies 2 values and returns the result.
     * @param lhs left hand side operand
     * @param rhs Right hand side operand
     * @return Result of the operation
    */
    template<typename T>
    friend __forceinline uint128_t operator*(uint128_t lhs, const T& rhs) noexcept {
        return lhs *= rhs;
    }
    /**
     * @brief Divides 2 values and returns the result.
     * @param lhs left hand side operand
     * @param rhs Right hand side operand
     * @return Result of the operation
    */
    template<typename T>
    friend __forceinline uint128_t operator/(uint128_t lhs, const T& rhs) {
        return lhs /= rhs;
    }
    /**
     * @brief Performs modulo and returns the result.
     * @param lhs left hand side operand
     * @param rhs Right hand side operand
     * @return Result of the operation
    */
    template<typename T>
    friend __forceinline uint128_t operator%(uint128_t lhs, const T& rhs) {
        return lhs %= rhs;
    }

    //
    // Binary math operators
    //
    /**
     * @brief Performs bitwise AND (&).
     * @param lhs left hand side operand
     * @param rhs Right hand side operand
     * @return Result of the operation
    */
    template<typename T>
    friend __forceinline uint128_t operator&(uint128_t lhs, const T& rhs) {
        return lhs &= rhs;
    }
    /**
     * @brief Performs bitwise OR (|).
     * @param lhs left hand side operand
     * @param rhs Right hand side operand
     * @return Result of the operation
    */
    template<typename T>
    friend __forceinline uint128_t operator|(uint128_t lhs, const T& rhs) {
        return lhs |= rhs;
    }
    /**
     * @brief Performs bitwise XOR (^).
     * @param lhs left hand side operand
     * @param rhs Right hand side operand
     * @return Result of the operation
    */
    template<typename T>
    friend __forceinline uint128_t operator^(uint128_t lhs, const T& rhs) {
        return lhs &= rhs;
    }

    //
    // Various math functions, implemented as friend functions
    //

    /**
     * @brief Calculates the left zero count of value x.
     * @param x input value.
     * @return lzc (uint32_t) of the result.
    */
    friend __forceinline uint64_t lzcnt128(const uint128_t& x) noexcept
    {
        return (x.high != 0) ? __lzcnt64(x.high) : 64 + __lzcnt64(x.low);
    }
    /**
     * @brief Calculates the square root using Newton's method.
     * Based on the book "Math toolkit for real time programming" by Jack W. Crenshaw 
     * @param x Value to calculate the root of
     * @return Square root of (x), zero when x <= 0.
    */
    friend UINT128_T_INLINE uint64_t sqrt(const uint128_t& x) noexcept
    {
        const auto expo = static_cast<uint32_t>(log2(x));
        if (expo == 0) {
            return (x == 1ull) ? 1 : 0;
        }

        uint128_t root = uint128_t::one();
        uint128_t e, temp;
        root <<= ((expo + 1) >> 1);
        
        // Newton iterations to reduce the error
        do {
            temp = x / root;
            // working with unsigned numbers, must keep the below positive at all times
            e = ((root > temp) ? (root - temp) : (temp - root)) >> 1; 
            root = (root + temp) >> 1;
        } while (e);

        return root;
    }
    /**
     * @brief Calculates the Log base 2 of x: log2(x)
     * Rounding is always towards zero so the maximu error is close to 1.
     * @param x The number to perform log2 on.
     * @return log2(x). Returns zero when x is zero.
    */
    friend __forceinline uint64_t log2(const uint128_t& x) noexcept
    {
        return (x) ? 127ull - lzcnt128(x) : 0;
    }
    /**
     * @brief Calculates the natural Log (base e) of x: log(x), rounded down to the nearest integer.
     * @param x The number to perform log on.
     * @return log(x)
    */
    friend UINT128_T_INLINE uint64_t log(const uint128_t& x) noexcept
    {
        // the table below holds the value of various powers of e (~2.718).
        // index i holds: ceil(pow(e, i)).
        static uint128_t lan_table[] = {
            "1",                                      // e^0
            "3",                                      // e^1
            "8",                                      // e^2
            "21",                                     // e^3
            "55",                                     // e^4
            "149",                                    // e^5
            "404",                                    // e^6
            "1097",                                   // e^7
            "2981",                                   // e^8
            "8104",                                   // e^9
            "22027",                                  // e^10
            "59875",                                  // e^11
            "162755",                                 // e^12
            "442414",                                 // e^13
            "1202605",                                // e^14
            "3269018",                                // e^15
            "8886111",                                // e^16
            "24154953",                               // e^17
            "65659970",                               // e^18
            "178482301",                              // e^19
            "485165196",                              // e^20
            "1318815735",                             // e^21
            "3584912847",                             // e^22
            "9744803447",                             // e^23
            "26489122130",                            // e^24
            "72004899338",                            // e^25
            "195729609429",                           // e^26
            "532048240602",                           // e^27
            "1446257064292",                          // e^28
            "3931334297145",                          // e^29
            "10686474581525",                         // e^30
            "29048849665248",                         // e^31
            "78962960182681",                         // e^32
            "214643579785917",                        // e^33
            "583461742527455",                        // e^34
            "1586013452313431",                       // e^35
            "4311231547115196",                       // e^36
            "11719142372802612",                      // e^37
            "31855931757113757",                      // e^38
            "86593400423993747",                      // e^39
            "235385266837019986",                     // e^40
            "639843493530054950",                     // e^41
            "1739274941520501048",                    // e^42
            "4727839468229346562",                    // e^43
            "12851600114359308276",                   // e^44
            "34934271057485095349",                   // e^45
            "94961194206024488746",                   // e^46
            "258131288619006739624",                  // e^47
            "701673591209763173866",                  // e^48
            "1907346572495099690526",                 // e^49
            "5184705528587072464088",                 // e^50
            "14093490824269387964493",                // e^51
            "38310080007165768493036",                // e^52
            "104137594330290877971835",               // e^53
            "283075330327469390044207",               // e^54
            "769478526514201713818275",               // e^55
            "2091659496012996153907072",              // e^56
            "5685719999335932222640349",              // e^57
            "15455389355901039303530767",             // e^58
            "42012104037905142549565935",             // e^59
            "114200738981568428366295719",            // e^60
            "310429793570191990870734215",            // e^61
            "843835666874145448907332949",            // e^62
            "2293783159469609879099352841",           // e^63
            "6235149080811616882909238709",           // e^64
            "16948892444103337141417836115",          // e^65
            "46071866343312915426773184429",          // e^66
            "125236317084221378051352196075",         // e^67
            "340427604993174052137690718701",         // e^68
            "925378172558778760024239791669",         // e^69
            "2515438670919167006265781174253",        // e^70
            "6837671229762743866755892826678",        // e^71
            "18586717452841279803403701812546",       // e^72
            "50523936302761041945570383321858",       // e^73
            "137338297954017618778418852980854",      // e^74
            "373324199679900164025490831726471",      // e^75
            "1014800388113888727832461784131717",     // e^76
            "2758513454523170206286469819902662",     // e^77
            "7498416996990120434675630591224061",     // e^78
            "20382810665126687668323137537172633",    // e^79
            "55406223843935100525711733958316613",    // e^80
            "150609731458503054835259413016767499",   // e^81
            "409399696212745469666091422932782905",   // e^82
            "1112863754791759412087071478183940806",  // e^83
            "3025077322201142338266566396443428743",  // e^84
            "8223012714622913510304328016407774696",  // e^85
            "22352466037347150474430657323327147400", // e^86
            "60760302250568721495223289381302760756", // e^87
            "165163625499399482076757134095306194944" // e^88
        };
        static constexpr uint64_t lan_table_len = array_length(lan_table);
        if (x < 3) return 0;
        if (x >= lan_table[lan_table_len - 1]) return lan_table_len - 1;

        // binary search the result
        uint64_t l = 0, h = lan_table_len - 1, res = (h + l) >> 1;

        while (l < h) {
            if (x < lan_table[res]) {
                h = res;
            }
            else if (l == res)
                break;
            else {
                l = res;
            }
            res = (h + l) >> 1;
        }
        return res;
    }
    /**
     * @brief Calculates Log base 10 of x: log10(x), rounded to the nearest integer.
     * @param x The number to perform log on.
     * @return log10(x)
    */
    friend UINT128_T_INLINE uint64_t log10(const uint128_t& x) noexcept
    {
        static constexpr uint64_t log10_max = 38;
        static uint128_t log10_table[log10_max + 1]; // holds all log10 values for multiples of 10
        // initialize the table
        if (!log10_table[log10_max]) {
            log10_table[0] = 1ull;
            for (uint32_t i = 1u; i <= log10_max; ++i) {
                log10_table[i] = log10_table[i - 1] * 10ull;
            }
        }

        if (!x) return 0;

        if (x >= log10_table[log10_max])
            return log10_max;

        // binary search the result
        uint64_t l = 0, h = 38, res = (h + l) >> 1;

        while (l < h) {
            if (x < log10_table[res]) {
                h = res;
            }
            else if (l == res)
                break;
            else {
                l = res;
            }
            res = (h + l) >> 1;
        }
        return res;
    }
    /**
     * @brief Calculates x to the power of y.
     * @param x Base value
     * @param y Exponent 
     * @return x to the power of y
    */
    friend UINT128_T_INLINE uint128_t pow(const uint128_t& x, uint32_t y) noexcept
    {
        // zero power always yields 1
        // Even if x is zero! Same behavior as pow(double, double) but different from Python which returns zero.
        if (y == 0) return 1;

        // special case where base is zero
        if (!x) return x;

        uint128_t res = 1;
        uint128_t b = x;
        while (y > 0) {
            if (y & 1)
                res *= b;
            y >>= 1;
            b *= b;
        }
        return res;
    }
}; //class uint128_t

} //namespace fp128

#pragma warning(pop)

#endif // #ifndef UINT128_T_H
