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

class uint128_t
{
    // build time validation of template parameters
    static_assert(sizeof(void*) == 8, "uint128_t is supported in 64 bit builds only!");

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
    static constexpr int32_t dbl_exp_bits = 11;                         // exponent bit count of a double variable
    static constexpr int32_t dbl_frac_bits = 52;                        // mantisa bit count of a double variable
    static constexpr double upper_base = ((double)(1ull << 63)) * 2.0;
public:
    typedef uint128_t type;
    typedef uint128_t* ptr_type;
    typedef uint128_t& ref_type;

    //
    // ctors
    //

    /**
     * @brief Default constructor, creates an instance with a value of zero.
    */
    uint128_t() noexcept :
        low(0), high(0) {}
    /**
     * @brief Copy constructor
     * @param other Object to copy from
    */
    uint128_t(const uint128_t& other) noexcept :
        low(other.low), high(other.high) {}
    /**
     * @brief Move constructor
     * Doesn't modify the right hand side object. Acts like a copy constructor.
     * @param other Object to copy from
    */
    uint128_t(const uint128_t&& other) noexcept :
        low(other.low), high(other.high) {}
    /**
     * @brief Constructor from the double type
     * Underflow goes to zero. Overflow, NaN and +-INF go to max supported positive value.
     * @param x Input value
    */
    uint128_t(double x) noexcept {
        // brute convert to uint64_t for easy bit manipulation
        const uint64_t i = *reinterpret_cast<uint64_t*>(&x);
        // very common case
        if (i == 0) {
            low = high = 0;
            return;
        }

        const int32_t e = FP128_GET_BITS(i, dbl_frac_bits, dbl_exp_bits) - 1023;
        uint64_t f = (i & FP128_MAX_VALUE_64(dbl_frac_bits));
        
        // overflow which catches NaN and Inf
        if (e > 127) {
            high = low = UINT64_MAX;
            return;
        }

        // normal number, produces non zero value
        if (e > 0) {
            // bit 52 in f is the unity value of the float. it needs to move to the unity position in fixed point
            f |= FP128_ONE_SHIFT(dbl_frac_bits);
            int32_t bits_to_shift = e - dbl_frac_bits;

            // f fits in high QWORD
            if (bits_to_shift >= 0) {
                if (bits_to_shift > 128 - dbl_frac_bits - 1) {
                }

                high = f << bits_to_shift;
                low = 0;
            }
            // shift right
            else {
                bits_to_shift = -bits_to_shift;
                high = 0;
                low = f >> bits_to_shift;
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
    uint128_t(uint64_t x) noexcept {
        low = x;
        high = 0;
    }
    /**
     * @brief Constructor from int64_t type
     * @param x Input value
    */
    uint128_t(int64_t x) noexcept {
        low = static_cast<uint64_t>(x);
        high = (x < 0) ? UINT64_MAX : 0;
    }
    /**
     * @brief Constructor from uint32_t type
     * @param x Input value
   */
    uint128_t(uint32_t x) noexcept {
        low = x;
        high = 0;
    }
    /**
     * @brief Constructor from int32_t type
     * @param x Input value
   */
    uint128_t(int32_t x) noexcept {
        low = static_cast<uint64_t>(x);
        high = (x < 0) ? UINT64_MAX : 0;
    }
    /**
     * @brief Constructor from const char* (C string).
     * Allows creating 128 bit values. Much slower than the other constructors.<BR>
     * Input string can be decimal or hex. Hex initialization is faster.
     * Throws an std::invalid_argument when encountering an illegal character.
     * @param x Input string
    */
    uint128_t(const char* x) {
        low = high = 0;
        char* str = _strdup(x);
        _strlwr_s(str, 1 + strlen(x));
        x = str;

        // trim leading white space
        while (*x && isspace(*x))
            ++x;
        
        int base = (0 == strncmp("0x", x, 2)) ? 16 : 10;
        if (base == 16)
            x += 2;
        
        // trim leading zeros
        while (*x == '0')
            ++x;
        
        // convert one digit at a time
        while (*x && !isspace(*x)) {
            uint64_t d = *x;
            if (d >= '0' && d <= '9')
                d -= '0';
            else if (base == 16 && (d >= 'a' && d <= 'f'))
                d = 10ull + d - 'a';
            else {
                throw std::invalid_argument("uint128_t: Invalid characters used to create an object");
            }
            
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
        uint128_t temp = x.c_str();
        *this = temp;
    }
    /**
     * @brief Constructor from the 2 uint128_t base elements, useful for creating special constants.
     * @param l Low QWORD
     * @param h High QWORD
    */
    uint128_t(uint64_t l, uint64_t h) noexcept:
        low(l), high(h) {
    }
    
    /**
     * @brief Destructor
    */
    ~uint128_t() noexcept {}
    /**
     * @brief Assignment operator
     * @param other Object to copy from 
     * @return This object.
    */
    UINT128_T_INLINE uint128_t& operator=(const uint128_t& other) noexcept {
        high = other.high;
        low = other.low;
        return *this;
    }
    /**
     * @brief Move assignment operator
     * @param other Object to copy from
     * @return This object.
    */
    UINT128_T_INLINE uint128_t& operator=(const uint128_t&& other) noexcept {
        high = other.high;
        low = other.low;
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
        double res = upper_base * high;
        res += low;
        return (float)res;
    }
    /**
     * @brief operator double - converts to a double
     * @return Object value.
    */
    UINT128_T_INLINE operator double() const noexcept {
        double res = upper_base * high;
        res += low;
        return res;
    }
    /**
     * @brief operator long double - converts to a long double
     * @return Object value.
    */
    UINT128_T_INLINE operator long double() const noexcept {
        double res = upper_base * high;
        res += low;
        return res;
    }
    /**
     * @brief Converts to a std::string (slow) string holds all meaningful fraction bits.
     * @return object string representation
    */
    UINT128_T_INLINE operator std::string() const {
        return uint128tostr();
    }
    /**
     * @brief Converts to a C string (slow) string holds all meaningful fraction bits.
     * @return object string representation
    */
    explicit UINT128_T_INLINE operator char*() const {
        return uint128tostr();
    }
    //
    // math operators
    //
    /**
     * @brief Adds the right hand side operand to this object to and returns the result.
     * @param other Right hand side operand
     * @return Temporary object with the result of the operation
    */
    UINT128_T_INLINE uint128_t operator+(const uint128_t& other) const {
        uint128_t temp(*this);
        return temp += other;
    }
    /**
     * @brief subtracts the right hand side operand to this object to and returns the result.
     * @param other Right hand side operand
     * @return Temporary object with the result of the operation
    */
    UINT128_T_INLINE uint128_t operator-(const uint128_t& other) const {
        uint128_t temp(*this);
        return temp -= other;
    }
    /**
     * @brief Multiplies the right hand side operand with this object to and returns the result.
     * @param other Right hand side operand
     * @return Temporary object with the result of the operation
    */
    UINT128_T_INLINE uint128_t operator*(const uint128_t& other) const {
        uint128_t temp(*this);
        return temp *= other;
    }
    /**
     * @brief Multiplies the right hand side operand with this object to and returns the result.
     * @param x Right hand side operand
     * @return Temporary object with the result of the operation
    */
    UINT128_T_INLINE uint128_t operator*(double x) const {
        // TODO: check what is the correct behavior
        uint128_t temp(*this);
        return temp *= uint128_t(x);
    }
    /**
     * @brief Multiplies the right hand side operand with this object to and returns the result.
     * @param x Right hand side operand
     * @return Temporary object with the result of the operation
    */
    UINT128_T_INLINE uint128_t operator*(int64_t x) const {
        uint128_t temp(*this);
        return temp *= x;
    }
    /**
     * @brief Multiplies the right hand side operand with this object to and returns the result.
     * @param x Right hand side operand
     * @return Temporary object with the result of the operation
    */
    UINT128_T_INLINE uint128_t operator*(int32_t x) const {
        uint128_t temp(*this);
        return temp *= x;
    }
    /**
     * @brief Divides this object by the right hand side operand and returns the result.
     * @param other Right hand side operand (denominator)
     * @return Temporary object with the result of the operation
    */
    UINT128_T_INLINE uint128_t operator/(const uint128_t& other) const {
        uint128_t temp(*this);
        return temp /= other;
    }
    /**
     * @brief Divides this object by the right hand side operand and returns the result.
     * @param x Right hand side operand (denominator)
     * @return Temporary object with the result of the operation
    */
    UINT128_T_INLINE uint128_t operator/(double x) const {
        if (x == 0)
            FP128_FLOAT_DIVIDE_BY_ZERO_EXCEPTION;

        uint128_t temp(*this);
        temp /= x;
        return temp;
    }
    /**
     * @brief Calculates modulo.
     * @param other Right hand side operand (denominator)
     * @return Temporary object with the result of the operation
    */
    UINT128_T_INLINE uint128_t operator%(const uint128_t& other) const {
        uint128_t temp(*this);
        return temp %= other;
    }
    /**
     * @brief Performs right shift operation.
     * @param shift bits to shift
     * @return Temporary object with the result of the operation
    */
    UINT128_T_INLINE uint128_t operator>>(int32_t shift) const {
        uint128_t temp(*this);
        return temp >>= shift;
    }
    /**
     * @brief Performs left shift operation.
     * @param shift bits to shift
     * @return Temporary object with the result of the operation
    */
    UINT128_T_INLINE uint128_t operator<<(int32_t shift) const {
        uint128_t temp(*this);
        return temp <<= shift;
    }
    /**
     * @brief Performs bitwise AND (&)
     * @param other Right hand side operand
     * @return Temporary object with the result of the operation
    */
    UINT128_T_INLINE uint128_t operator&(const uint128_t& other) const {
        uint128_t temp(*this);
        return temp &= other;
    }
    /**
     * @brief Performs bitwise OR (|)
     * @param other Right hand side operand
     * @return Temporary object with the result of the operation
    */
    UINT128_T_INLINE uint128_t operator|(const uint128_t& other) const {
        uint128_t temp(*this);
        return temp |= other;
    }
    /**
     * @brief Performs bitwise XOR (^)
     * @param other Right hand side operand
     * @return Temporary object with the result of the operation
    */
    UINT128_T_INLINE uint128_t operator^(const uint128_t& other) const {
        uint128_t temp(*this);
        return temp ^= other;
    }
    /**
     * @brief Add a value to this object
     * @param other Right hand side operand
     * @return This object.
    */
    UINT128_T_INLINE uint128_t& operator+=(const uint128_t& other) {
        unsigned char carry = _addcarry_u64(0, low, other.low, &low);
        high += other.high + carry;
        return *this;
    }
    /**
     * @brief Subtract a value to this object
     * @param other Right hand side operand
     * @return This object.
    */
    UINT128_T_INLINE uint128_t& operator-=(const uint128_t& other) {
        uint128_t temp = other;
        twos_complement128(temp.low, temp.high);
        unsigned char carry = _addcarry_u64(0, low, temp.low, &low);
        high += temp.high + carry;
        return *this;
    }
    /**
     * @brief Multiplies a value to this object
     * @param other Right hand side operand
     * @return This object.
    */
    UINT128_T_INLINE uint128_t& operator*=(const uint128_t& other) noexcept{
        uint128_t temp = *this;

        // multiply low QWORDs
        low = _umul128(temp.low, other.low, &high);

        // multiply low this and high other
        high += temp.low * other.high;

        // multiply high this and low other
        high += temp.high * other.low;
        return *this;
    }
    /**
     * @brief Multiplies a value to this object
     * @param x Right hand side operand
     * @return This object.
    */
    UINT128_T_INLINE uint128_t& operator*=(int32_t x) {
        bool sign = (x < 0);
        // alway do positive multiplication
        uint64_t uval = (uint64_t)((sign) ? -x : x);
        uint64_t temp;

        // multiply low QWORDs
        low = _umul128(low, uval, &temp);
        high = high * uval + temp;
        if (sign) {
            twos_complement128(low, high);
        }
        return *this;
    }
    /**
     * @brief Multiplies a value to this object
     * @param x Right hand side operand
     * @return This object.
    */
    UINT128_T_INLINE uint128_t& operator*=(uint32_t x) {
        uint64_t temp;

        // multiply low QWORDs
        low = _umul128(low, x, &temp);
        high = high * x + temp;
        return *this;
    }
    /**
     * @brief Multiplies a value to this object
     * @param x Right hand side operand
     * @return This object.
    */
    UINT128_T_INLINE uint128_t& operator*=(int64_t x) {
        bool sign = (x < 0);
        // alway do positive multiplication
        uint64_t uval = (uint64_t)((sign) ? -x : x);
        uint64_t temp;

        // multiply low QWORDs
        low = _umul128(low, uval, &temp);
        high = high * uval + temp;
        return *this;
    }
    /**
     * @brief Multiplies a value to this object
     * @param x Right hand side operand
     * @return This object.
    */
    UINT128_T_INLINE uint128_t& operator*=(uint64_t x) {
        uint64_t temp;

        // multiply low QWORDs
        low = _umul128(low, x, &temp);
        high = high * x + temp;
        return *this;
    }
    /**
     * @brief Divide this object by x.
     * @param other Right hand side operator (denominator)
     * @return this object.
    */
    inline uint128_t& operator/=(const uint128_t& other) {
        // check some trivial cases
        if (other.is_zero()) {
            FP128_INT_DIVIDE_BY_ZERO_EXCEPTION;
        }
        if (is_zero() || other > *this) {
            low = high = 0;
            return *this;
        }
        if (other == *this) {
            low = 1;
            high = 0;
            return *this;
        }

        uint64_t q[2]{}; 
        uint64_t nom[2] = { low, high };

        // optimization for when dividing by a small integer
        if (other.high == 0) {
            if (div_64bit((uint64_t*)q, nullptr, (uint64_t*)nom, other.low, 2)) {
                FP128_INT_DIVIDE_BY_ZERO_EXCEPTION;
            }
        }
        else {
            uint64_t denom[2] = {other.low, other.high};
            if (div_32bit((uint32_t*)q, nullptr, (uint32_t*)nom, (uint32_t*)denom, 2ll * array_length(nom), 2ll * array_length(denom))) {
                FP128_INT_DIVIDE_BY_ZERO_EXCEPTION;
            }
        }
        low = q[0];
        high = q[1];
        return *this;
    }
    /**
     * @brief %= operator
     * @param other Modulo operand.
     * @return This object.
    */
    inline uint128_t& operator%=(const uint128_t& other) {
        // check some trivial cases
        if (other.is_zero()) {
            FP128_INT_DIVIDE_BY_ZERO_EXCEPTION;
        }
        if (*this < other) {
            return *this;
        }
        if (*this == other) {
            low = 0; high = 0;
        }

        uint64_t q[2]{};
        uint64_t nom[2] = { low, high };

        // optimization for when dividing by a small integer
        if (other.high == 0) {
            if (div_64bit((uint64_t*)q, &low, (uint64_t*)nom, other.low, 2)) {
                FP128_INT_DIVIDE_BY_ZERO_EXCEPTION;
            }
            high = 0;
        }
        else {
            uint64_t denom[2] = { other.low, other.high };
            if (div_32bit((uint32_t*)q, (uint32_t*)words, (uint32_t*)nom, (uint32_t*)denom, 2ll * array_length(nom), 2ll * array_length(denom))) {
                FP128_INT_DIVIDE_BY_ZERO_EXCEPTION;
            }
        }
        return *this;
    }
    /**
     * @brief Shift right this object.
     * Shifting is done without rounding to match uint64_t behavior
     * @param shift Bits to shift. Negative or very high values cause undefined behavior.
     * @return This object.
    */
    UINT128_T_INLINE uint128_t& operator>>=(int32_t shift) noexcept {
        // 0-64 bit shift - most common
        if (shift < 64) {
            low = shift_right128(low, high, (uint8_t)shift);
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
        if (shift < 64) {
            high = shift_left128(low, high, (unsigned char)shift);
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
     * @param other AND mask.
     * @return This object.
    */
    UINT128_T_INLINE uint128_t& operator&=(const uint128_t& other) noexcept {
        low &= other.low;
        high &= other.high;
        return *this;
    }
    /**
     * @brief Bitwise OR=
     * @param other OR mask.
     * @return This object.
    */
    UINT128_T_INLINE uint128_t& operator|=(const uint128_t& other) noexcept {
        low |= other.low;
        high |= other.high;
        return *this;
    }
    /**
     * @brief Bitwise XOR=
     * @param other XOR mask.
     * @return This object.
    */
    UINT128_T_INLINE uint128_t& operator^=(const uint128_t& other) noexcept {
        low ^= other.low;
        high ^= other.high;
        return *this;
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
    */
    UINT128_T_INLINE uint128_t operator-() const {
        uint128_t temp = *this;
        twos_complement128(temp.low, temp.high);
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
    UINT128_T_INLINE bool operator==(const uint128_t& other) const noexcept {
        return high == other.high && low == other.low;
    }
    /**
     * @brief Return true when objects are not equal. Can be used as logical XOR.
     * @param other Righthand operand.
     * @return True of not equal.
    */
    UINT128_T_INLINE bool operator!=(const uint128_t& other) const noexcept {
        return low != other.low || high != other.high;
    }
    /**
     * @brief Return true if this object is small than the other
     * @param other Righthand operand.
     * @return True when this object is smaller.
    */
    UINT128_T_INLINE bool operator<(const uint128_t& other) const noexcept {
        // MSB is the same, check the LSB
        if (high == other.high)
            return low < other.low;

        return high < other.high;
    }
    /**
     * @brief Return true this object is small or equal than the other
     * @param other Righthand operand.
     * @return True when this object is smaller or equal.
    */
    UINT128_T_INLINE bool operator<=(const uint128_t& other) const noexcept {
        return !(*this > other);
    }
    /**
     * @brief Return true this object is larger than the other
     * @param other Righthand operand.
     * @return True when this objext is larger.
    */
    UINT128_T_INLINE bool operator>(const uint128_t& other) const noexcept {
        // MSB is the same, check the LSB
        if (high == other.high)
            return low > other.low;

        return high > other.high;
    }
    /**
     * @brief Return true this object is larger or equal than the other
     * @param other Righthand operand.
     * @return True when this objext is larger or equal.
    */
    UINT128_T_INLINE bool operator>=(const uint128_t& other) const noexcept {
        return !(*this < other);
    }
    //
    // useful public functions
    //
    /**
     * @brief Returns true if the value is an int (fraction is zero)
     * @return True when the fraction is zero.
    */
    UINT128_T_INLINE constexpr bool is_int() const noexcept
    {
        return true;
    }
    /**
     * @brief Returns true if the value positive (incuding zero)
     * @return True when the the value positive
    */
    UINT128_T_INLINE constexpr bool is_positive() const noexcept
    {
        return true;
    }
    /**
     * @brief Returns true if the value negative (smaller than zero)
     * @return True when the the value negative
    */
    UINT128_T_INLINE constexpr bool is_negative() const noexcept
    {
        return false;
    }
    /**
     * @brief Returns true if the value is zero
     * @return Returns true if the value is zero 
    */
    UINT128_T_INLINE bool is_zero() const noexcept
    {
        return 0 == low && 0 == high;
    }
    /**
     * @brief get a specific bit within the 128 data
     * @param bit bit to get [0,127]
     * @return 0 or 1. Undefined when bit > 127
    */
    UINT128_T_INLINE int32_t get_bit(unsigned bit) const noexcept
    {
        if (bit < 64) {
            return FP128_GET_BIT(low, bit);
        }
        return FP128_GET_BIT(high, bit-64);
    }
    /**
     * @brief Returns the exponent of the object - like the base 2 exponent of a floating point
     * A value of 2.1 would return 1, [0.5,1.0) would return -1.
     * @return Exponent of the number
    */
    UINT128_T_INLINE int32_t get_exponent() const
    {
        int32_t s = lzcnt128(*this) - 1;
        return s;
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
     * Additional calls to this function from the same thread, overwrite the previous result.
     * @return C string with describing the value of the object.
    */
    UINT128_T_INLINE char* hex() const {
        constexpr int buff_size = 35;
        static thread_local char str[buff_size];
        
        if (high) {
            snprintf(str, buff_size, "0x%llX%016llX", high, low);
        }
        else {
            snprintf(str, buff_size, "0x%llX", low);
        }
        
        return str;
    }
private:
    /**
     * @brief Converts this object to a C string.
     * The returned string is a statically thread-allocated buffer.
     * Additional calls to this function from the same thread, overwrite the previous result.
     * @return C string with describing the value of the object.
    */
    UINT128_T_INLINE char* uint128tostr() const {
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
        while (temp.high) {
            --p;
            uint64_t nom[2] = { temp.low, temp.high };
            div_64bit((uint64_t*)q, &digit, (uint64_t*)nom, 10ull, 2);
            temp.low = q[0];
            temp.high = q[1];
            *p = static_cast<char>(digit + '0');
        }

        while (temp.low) {
            --p;
            lldiv_t res = lldiv(temp.low, 10ull);
            temp.low = res.quot;
            *p = static_cast<char>(res.rem + '0');
        }
        return p;
    }
public:
    /**
     * @brief Calculates the left zero count of value x, ignoring the sign.
     * @param x input value.
     * @return lzc (uint32_t) of the result.
    */
    friend UINT128_T_INLINE uint32_t lzcnt128(const uint128_t& x) noexcept
    {
        return (int32_t)(__lzcnt64(x.high) + __lzcnt64(x.low));
    }
    /**
     * @brief Calculates the square root using Newton's method.
     * Based on the book "Math toolkit for real time programming" by Jack W. Crenshaw 
     * @param x Value to calculate the root of
     * @param iterations how many iterations to perform (more is more accurate). Sensible values are 0-5.
     * @return Square root of (x), zero when x <= 0.
    */
    friend UINT128_T_INLINE uint128_t sqrt(const uint128_t& x, uint32_t iterations = 3) noexcept
    {
        FP128_NOT_IMPLEMENTED_EXCEPTION;
    }
    /**
     * @brief Calculates the Log base 2 of x: log2(x)
     * @param x The number to perform log2 on.
     * @return log2(x)
    */
    friend UINT128_T_INLINE uint128_t log2(uint128_t x) noexcept
    {
        uint128_t y = x.get_exponent();
        return y;
    }
    /**
     * @brief Calculates the natural Log (base e) of x: log(x)
     * @param x The number to perform log on.
     * @return log2(x)
    */
    friend UINT128_T_INLINE uint128_t log(uint128_t x) noexcept
    {
        FP128_NOT_IMPLEMENTED_EXCEPTION;
    }
    /**
     * @brief Calculates Log base 10 of x: log10(x)
     * @param x The number to perform log on.
     * @return log10(x)
    */
    friend UINT128_T_INLINE uint128_t log10(uint128_t x) noexcept
    {
        FP128_NOT_IMPLEMENTED_EXCEPTION;
    }
}; //class uint128_t

} //namespace fp128

#endif // #ifndef FIXED_POINT128_H
