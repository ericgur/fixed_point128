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

#ifndef FIXED_POINT128_H
#define FIXED_POINT128_H

// override some static analysis checks
#pragma warning(push)
#pragma warning(disable: 26472) // Don't use a static_cast for arithmetic conversions. Use brace initialization
#pragma warning(disable: 26485) // No array to pointer decay
#pragma warning(disable: 26481) // Don't use pointer arithmetic. Use span instead
#pragma warning(disable: 26446) // Prefer to use gsl::at() instead of unchecked subscript operator
#pragma warning(disable: 26482) // Only index into arrays using constant expressions
#pragma warning(disable: 26408) // Avoid malloc() and free(), prefer the nothrow version of new with delete


#include "fixed_point128_shared.h"
//#include "uint128_t.h"

namespace fp128 {


    
/***********************************************************************************
*                                  Forward declarations
************************************************************************************/
class fp128_gtest; // Google test class
template<int32_t I> class fixed_point128;
// Note: Release builds will fail without these forward declarations. Hints towards compiler a bug (VS2022 v17.4)
// The compiler and Inteliisense don't match these functions in some cases and try to use the CRT versions which 
// a causes a compilations error.
// CRT style function
template<int32_t I> fixed_point128<I> fabs(const fixed_point128<I>& x) noexcept;
template<int32_t I> fixed_point128<I> floor(const fixed_point128<I>& x) noexcept;
template<int32_t I> fixed_point128<I> ceil(const fixed_point128<I>& x) noexcept;
template<int32_t I> fixed_point128<I> trunc(const fixed_point128<I>& x) noexcept;
template<int32_t I> fixed_point128<I> round(const fixed_point128<I>& x) noexcept;
template<int32_t I> fixed_point128<I> fmod(const fixed_point128<I>& x, const fixed_point128<I>& y);
template<int32_t I> fixed_point128<I> modf(const fixed_point128<I>& x, fixed_point128<I>* iptr) noexcept;
template<int32_t I> fixed_point128<I> sqrt(const fixed_point128<I>& x, uint32_t iterations = 3) noexcept;
template<int32_t I> fixed_point128<I> sin(fixed_point128<I> x) noexcept;
template<int32_t I> fixed_point128<I> asin(fixed_point128<I> x) noexcept;
template<int32_t I> fixed_point128<I> cos(fixed_point128<I> x) noexcept;
template<int32_t I> fixed_point128<I> acos(fixed_point128<I> x) noexcept;
template<int32_t I> fixed_point128<I> tan(fixed_point128<I> x) noexcept;
template<int32_t I> fixed_point128<I> atan(fixed_point128<I> x) noexcept;
template<int32_t I> fixed_point128<I> exp(const fixed_point128<I>& x) noexcept;
template<int32_t I> fixed_point128<I> exp2(const fixed_point128<I>& x) noexcept;
template<int32_t I> fixed_point128<I> expm1(const fixed_point128<I>& x) noexcept;
template<int32_t I> fixed_point128<I> pow(const fixed_point128<I>& x, const fixed_point128<I>& y) noexcept;
template<int32_t I> fixed_point128<I> log(fixed_point128<I> x) noexcept;
template<int32_t I> fixed_point128<I> log2(fixed_point128<I> x) noexcept;
template<int32_t I> fixed_point128<I> log10(fixed_point128<I> x) noexcept;
template<int32_t I> fixed_point128<I> logb(fixed_point128<I> x) noexcept;
template<int32_t I> fixed_point128<I> log1p(fixed_point128<I> x) noexcept;
// non CRT function
template<int32_t I> uint64_t lzcnt128(const fixed_point128<I>& x) noexcept;
template<int32_t I> fixed_point128<I> reciprocal(const fixed_point128<I>& x) noexcept;
template<int32_t I> void fact_reciprocal(int x, fixed_point128<I>& res) noexcept;


/***********************************************************************************
*                                  Main Code
************************************************************************************/

/**
 * @brief 128 bit fixed point class template.
 * 
 * This template provides a floating point-like type that performs various math operations quickly compared to traditional high precision libraries.
 * The only template parameter <B>I</B> is the bit count of the integer part. The fraction part complements <B>I</B> to 128 bit.<BR>
 * <B>I</B> is limited to the range [1,64] in order to simplify the implementation and increase preformance.<BR>
 * This restriction is enforced at compile time.
 * All of fixed_point128's methods are inline for maximum performance.
 *
 * <B>Implementation notes:</B>
 * <UL>
 * <LI>Overflow is handled silently, similar to builtin integer operations.</LI>
 * <LI>Unlike integer operations, the sign is held in a separate data member which simplifies the implementation.</LI>
 * <LI>A fixed_point128 object is not thread safe. Accessing a const object from multiple threads is safe.</LI>
 * <LI>fixed_point128 is <B>conditionally safe</B>, 2 different non const objects can be accessed concurrently.</LI>
 * <LI>Only 64 bit builds are supported.</LI>
 * </UL>
*/
template<int32_t I>
class fixed_point128
{
    // build time validation of template parameters
    static_assert(1 <= I && I <= 64, "Template parameter <I> must be in the [1,64] range!");
    static_assert(sizeof(void*) == 8, "fixed_point128 is supported in 64 bit builds only!");

    // friends
    friend class fixed_point128; // this class is a friend of all its template instances. Avoids awkward getter/setter functions.
    friend class fp128_gtest;
private:
    //
    // members
    //
    uint64_t low;  // lower QWORD
    uint64_t high; // upper QWORD
    uint32_t sign; // 0 = positive, 1 negative

    // useful const calculations
    static constexpr int32_t F = 128 - I;                               // fraction bit count
    static constexpr int32_t upper_frac_bits = F - 64;                  // how many bits of the fraction exist in the upper QWORD
    static constexpr uint64_t unity = 1ull << upper_frac_bits;          // upper QWORD value equal to '1'
    static inline const double upper_unity = ::pow(2, 64 - F);     // convert upper QWORD to floating point
    static inline const double lower_unity = ::pow(2, -F);         // convert lower QWORD to floating point
    static constexpr uint64_t int_mask = UINT64_MAX << upper_frac_bits; // mask of the integer bits in the upper QWORD
    static constexpr int32_t max_frac_digits = (int)(F / 3.3);          // meaningful base 10 digits for the fraction
public:
    static constexpr uint64_t max_int_value = int_mask >> upper_frac_bits;
    typedef fixed_point128<I> type;
    typedef fixed_point128<I>* ptr_type;
    typedef fixed_point128<I>& ref_type;

    //
    // ctors
    //

    /**
     * @brief Default constructor, creates an instance with a value of zero.
    */
    __forceinline constexpr fixed_point128() noexcept :
        low(0), high(0), sign(0) {}
    /**
     * @brief Copy constructor
     * @param other Object to copy from
    */
    __forceinline fixed_point128(const fixed_point128& other) noexcept :
        low(other.low), high(other.high), sign(other.sign) {}
    /**
     * @brief cross-template Copy constructor, can be used between two different fixed_point128 templates
     * @param other fixed_point128 instance with from a different template instance.
     * @return This object.
    */
    template<int32_t I2>
    __forceinline fixed_point128(const fixed_point128<I2>& other) noexcept
    {
        sign = other.sign;
        if constexpr (I == I2) {
            high = other.high;
            low = other.low;
        }
        // other has less integer bits and more fraction bits
        else if constexpr (I < I2) {
            // shift left by I2 - I bits
            constexpr int shift = I2 - I;
            low = other.low << shift;
            high = shift_left128(other.low, other.high, shift);
        }
        // other has more integer bits and less fraction bits
        else { // I > I2
            // shift right by I - I2 bits
            constexpr int shift = I - I2;
            low = shift_right128_round(other.low, other.high, shift);
            high = other.high >> shift;
        }
    }

    /**
     * @brief Move constructor
     * Doesn't modify the right hand side object. Acts like a copy constructor.
     * @param other Object to copy from
    */
    __forceinline fixed_point128(const fixed_point128&& other) noexcept :
        low(other.low), high(other.high), sign(other.sign) {}
    /**
     * @brief Constructor from the double type
     * Underflow goes to zero. Overflow, NaN and +-INF go to max supported positive value.
     * @param x Input value
    */
    FP128_INLINE fixed_point128(double x) noexcept {
        // very common case
        if (x == 0) {
            low = high = 0;
            sign = 0;
            return;
        }

        // hack the double bit fields
        const Double d(x);

        sign = d.s;
        const int32_t e = static_cast<int32_t>(d.e) - 1023;
        uint64_t f = d.f;
        
        // overflow which also catches NaN and Inf
        if (e >= I) {
            high = low = UINT64_MAX;
            sign = 0;
            return;
        }

        // normal number, produces non zero value
        if (e >= -F) {
            // bit 52 in f is the unity value of the float. it needs to move to the unity position in fixed point
            f |= FP128_ONE_SHIFT(dbl_frac_bits);
            int32_t bits_to_shift = 64 - I - dbl_frac_bits + e;

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
                    // f has bit 52 high, shift it left to moves to bit 63
                    f <<= 63 - dbl_frac_bits;
                    bits_to_shift -= 64 - (63 - dbl_frac_bits);
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
     * @brief Constructor from uint64_t type
     * @param x Input value
    */
    __forceinline fixed_point128(uint64_t x) noexcept {
        low = 0;
        sign = 0;
        high = x << upper_frac_bits;
    }
    /**
     * @brief Constructor from int64_t type
     * @param x Input value
    */
    __forceinline fixed_point128(int64_t x) noexcept {
        low = 0;
        sign = FP128_GET_BIT(x, 63);
        high = ((sign != 0) ? -x : x) << upper_frac_bits;
    }
    /**
     * @brief Constructor from uint32_t type
     * @param x Input value
   */
    __forceinline fixed_point128(uint32_t x) noexcept {
        low = 0;
        sign = 0;
        high = (uint64_t)x << upper_frac_bits;
    }
    /**
     * @brief Constructor from int32_t type
     * @param x Input value
   */
    __forceinline fixed_point128(int32_t x) noexcept {
        low = 0;
        sign = FP128_GET_BIT(x, 31);
        high = static_cast<uint64_t>((sign != 0) ? -x : x) << upper_frac_bits;
    }
    /**
     * @brief Constructor from const char* (C string).
     * Accurate up to 37 digits after the decimal point.
     * Allows creating very high precision values. Much slower than the other constructors.
     * @param x Input string
    */
    fixed_point128(const char* x) noexcept {
        static fixed_point128<1> base10_table[] = {
        //  {low QWORD,             high QWORD,         sign}
            {0x0000000000000000ull, 0x8000000000000000ull, 0},  // 10^0, not used, makes the code simpler
            {0xCCCCCCCCCCCCCCCDull, 0x0CCCCCCCCCCCCCCCull, 0},  // 10^-1
            {0x147AE147AE147AE1ull, 0x0147AE147AE147AEull, 0},  // 10^-2
            {0xCED916872B020C4Aull, 0x0020C49BA5E353F7ull, 0},  // 10^-3
            {0x94AF4F0D844D013Bull, 0x000346DC5D638865ull, 0},  // 10^-4
            {0xC21187E7C06E19B9ull, 0x000053E2D6238DA3ull, 0},  // 10^-5
            {0xC69B5A63F9A49C2Cull, 0x000008637BD05AF6ull, 0},  // 10^-6
            {0x7A42BC3D32907604ull, 0x000000D6BF94D5E5ull, 0},  // 10^-7
            {0x8C39DF9FB841A567ull, 0x00000015798EE230ull, 0},  // 10^-8
            {0xDAD2965CC5A02A24ull, 0x0000000225C17D04ull, 0},  // 10^-9
            {0xAF7B756FAD5CD103ull, 0x0000000036F9BFB3ull, 0},  // 10^-10
            {0x5E592557F7BC7B4Dull, 0x00000000057F5FF8ull, 0},  // 10^-11
            {0x096F5088CBF93F88ull, 0x00000000008CBCCCull, 0},  // 10^-12
            {0x3424BB40E132865Aull, 0x00000000000E12E1ull, 0},  // 10^-13
            {0xB86A12B9B01EA709ull, 0x0000000000016849ull, 0},  // 10^-14
            {0x5F3DCEAC2B3643E7ull, 0x0000000000002407ull, 0},  // 10^-15
            {0x5652FB1137856D31ull, 0x000000000000039Aull, 0},  // 10^-16
            {0x3BD5191B525A2485ull, 0x000000000000005Cull, 0},  // 10^-17
            {0x392EE8E921D5D074ull, 0x0000000000000009ull, 0},  // 10^-18
            {0xEC1E4A7DB69561A5ull, 0x0000000000000000ull, 0},  // 10^-19
            {0x179CA10C9242235Dull, 0x0000000000000000ull, 0},  // 10^-20
            {0x025C768141D369F0ull, 0x0000000000000000ull, 0},  // 10^-21
            {0x003C7240202EBDCBull, 0x0000000000000000ull, 0},  // 10^-22
            {0x00060B6CD004AC94ull, 0x0000000000000000ull, 0},  // 10^-23
            {0x00009ABE14CD4475ull, 0x0000000000000000ull, 0},  // 10^-24
            {0x00000F79687AED3Full, 0x0000000000000000ull, 0},  // 10^-25
            {0x0000018C240C4AEDull, 0x0000000000000000ull, 0},  // 10^-26
            {0x000000279D346DE4ull, 0x0000000000000000ull, 0},  // 10^-27
            {0x00000003F61ED7CAull, 0x0000000000000000ull, 0},  // 10^-28
            {0x0000000065697BFBull, 0x0000000000000000ull, 0},  // 10^-29
            {0x000000000A2425FFull, 0x0000000000000000ull, 0},  // 10^-30
            {0x0000000001039D66ull, 0x0000000000000000ull, 0},  // 10^-31
            {0x000000000019F624ull, 0x0000000000000000ull, 0},  // 10^-32
            {0x000000000002989Dull, 0x0000000000000000ull, 0},  // 10^-33
            {0x0000000000004276ull, 0x0000000000000000ull, 0},  // 10^-34
            {0x00000000000006A5ull, 0x0000000000000000ull, 0},  // 10^-35
            {0x00000000000000AAull, 0x0000000000000000ull, 0},  // 10^-36
            {0x0000000000000011ull, 0x0000000000000000ull, 0},  // 10^-37
            {0x0000000000000002ull, 0x0000000000000000ull, 0}   // 10^-38    
        };
        constexpr uint64_t base10_table_size = array_length(base10_table);
        static_assert(max_frac_digits < base10_table_size);

        low = high = 0;
        sign = 0;
        if (x == nullptr) return;
        
        char* str = _strdup(x);
        char* p = str;
        if (p == nullptr) return;

        // skip white space
        while (*p == ' ') ++p;

        if (*p == '\0') {
            *this = fixed_point128();
            free(str);
            return;
        }

        // set negative sign if needed
        if (*p == '-') {
            sign = 1;
            ++p;
        }
        char* dec = strchr(p, '.');
        // number is an integer
        if (dec == nullptr) {
            high = std::strtoull(p, nullptr, 10) << upper_frac_bits;
            free(str);
            return;
        }

        // number is a float, get the integer part using strtoull()
        *dec = '\0';
        const uint64_t int_val = std::strtoull(p, nullptr, 10) << upper_frac_bits;

        p = dec + 1;
        int32_t digits = 0;
        fixed_point128<1> frac;
        // multiply each digits by 10^-n
        while (++digits < base10_table_size && isdigit(*p)) {
            uint32_t d = static_cast<uint64_t>(p[0] - '0');
            frac += base10_table[digits] * d;
            ++p;
        }

        frac >>= (I - 1);
        low = frac.low;
        high = frac.high + int_val;
        free(str);
    }
    /**
     * @brief Constructor from std::string.
     * Accurate to 37 digits after the decimal point.
     * Allows creating very high precision values. Much slower than the other constructors.
     * @param x Input string
    */
    FP128_INLINE fixed_point128(const std::string& x) noexcept {
        fixed_point128 temp = x.c_str();
        *this = temp;
    }
    /**
     * @brief Constructor from the 3 fixed_point128 base elements, useful for creating very small constants.
     * @param l Low QWORD
     * @param h High QWORD
     * @param s Sign - zero for positive, 1 for negative.
    */
    FP128_INLINE fixed_point128(uint64_t l, uint64_t h, uint32_t s) noexcept:
        low(l), high(h) ,sign(s) {
        sign = (sign != 0);
    }
    
    /**
     * @brief Destructor
    */
    ~fixed_point128() noexcept {}
    /**
     * @brief Assignment operator
     * @param other Object to copy from 
     * @return This object.
    */
    FP128_INLINE fixed_point128& operator=(const fixed_point128& other) noexcept {
        high = other.high;
        low = other.low;
        sign = other.sign;
        return *this;
    }
    /**
     * @brief Move assignment operator
     * @param other Object to copy from
     * @return This object.
    */
    FP128_INLINE fixed_point128& operator=(const fixed_point128&& other) noexcept {
        high = other.high;
        low = other.low;
        sign = other.sign;
        return *this;
    }
    /**
     * @brief cross-template assignment operator, can be used between two different fixed_point128 templates
     * @param other fixed_point128 instance with from a different template instance.
     * @return This object.
    */
    template<int32_t I2>
    FP128_INLINE fixed_point128<I>& operator=(const fixed_point128<I2>& other) noexcept
    {
        sign = other.sign;
        if constexpr (I == I2) {
            high = other.high;
            low = other.low;
        }
        // other has less integer bits and more fraction bits
        else if constexpr (I < I2) {
            // shift left by I2 - I bits
            constexpr int shift = I2 - I;
            low = other.low << shift;
            high = shift_left128(other.low, other.high, shift);
        }
        // other has more integer bits and less fraction bits
        else { // I > I2
            // shift right by I - I2 bits
            constexpr int shift = I - I2;
            low = shift_right128_round(other.low, other.high, shift);
            high = other.high >> shift;
        }

        return *this;
    }

    //
    // conversion operators
    //
    /**
     * @brief operator uint64_t - converts to a uint64_t
     * @return Object value.
    */
    FP128_INLINE operator uint64_t() const noexcept {
        return (high >> upper_frac_bits) & UINT64_MAX;
    }
    /**
     * @brief operator int64_t - converts to a int64_t
     * @return Object value.
    */
    FP128_INLINE operator int64_t() const noexcept {
        const int64_t res = (sign) ? -1ll : 1ll;
        return res * ((high >> upper_frac_bits) & UINT64_MAX);
    }
    /**
     * @brief operator uint32_t - converts to a uint32_t
     * @return Object value.
    */
    FP128_INLINE operator uint32_t() const noexcept {
        return (high >> upper_frac_bits) & UINT32_MAX;
    }
    /**
     * @brief operator int32_t - converts to a int32_t
     * @return Object value.
    */
    FP128_INLINE operator int32_t() const noexcept {
        int32_t res = (sign) ? -1 : 1;
        return res * ((int32_t)((int64_t)high >> upper_frac_bits) & (UINT32_MAX));
    }
    /**
     * @brief operator float - converts to a float
     * @return Object value.
    */
    FP128_INLINE operator float() const noexcept {
        if (!*this)
            return 0.0f;
        constexpr uint64_t f_mask = FP128_MAX_VALUE_64(flt_frac_bits);

        Float res;
        res.s = sign;
        int32_t s = (int32_t)lzcnt128(*this);
        int32_t msb = 127 - s;
        auto expo = I - 1 - s;

        res.e = 127u + expo;
        // get the 23 bits right of the msb.
        fixed_point128 temp(*this);
        int shift = msb - flt_frac_bits;
        if (shift > 0) {
            res.f = f_mask & shift_right128(temp.low, temp.high, shift);
        }
        else {
            res.f = f_mask & shift_left128(temp.low, temp.high, -shift);
        }
        // round
        if (res.f == f_mask) {
            res.f = 0;
            ++res.e;
        }
        return res.val;
    }
    /**
     * @brief operator double - converts to a double
     * @return Object value.
    */
    FP128_INLINE operator double() const noexcept {
        if (!*this)
            return 0.0;
        constexpr uint64_t f_mask = FP128_MAX_VALUE_64(dbl_frac_bits);
        Double res;
        res.s = sign;
        const int32_t s = static_cast<int32_t>(lzcnt128(*this));
        const int32_t msb = 127 - s;
        const auto expo = I - 1 - s;

        res.e = 1023ull + expo;
        // get the 52 bits right of the msb.
        fixed_point128 temp(*this);
        const int shift = msb - dbl_frac_bits;
        if (shift > 0) {
            res.f = f_mask & shift_right128(temp.low, temp.high, shift);
        }
        else {
            res.f = f_mask & (temp.low << -shift);
        }
        // round
        if (res.f == f_mask) {
            res.f = 0;
            ++res.e;
        }

        return res.val;
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
        return operator char*();
    }
    /**
     * @brief Converts to a C string (slow) string holds all meaningful fraction bits.
     * @return object string representation
    */
    explicit FP128_INLINE operator char*() const noexcept {
        static thread_local char str[128]; // need roughly a (meaningful) decimal digit per 3.2 bits

        char* p = &str[0];
        fixed_point128 temp = *this;

        //number is negative
        if (temp.sign)
            *p++ = '-';

        uint64_t integer = FP128_GET_BITS(temp.high, upper_frac_bits, 63);
        p += snprintf(p, sizeof(str) + p - str, "%llu", integer);
        temp.high &= ~int_mask; // remove the integer part
        // check if temp has additional digits (not zero)
        if (temp) {
            *p++ = '.';
        }
        // the faster way, requires temp *= 10 not overflowing
        int digits = 0;
        uint64_t res[2]{};
        while (digits++ < max_frac_digits && temp) {
            if constexpr (I < 4) {

                res[0] = _umul128(high, 10ull, &res[1]); // multiply by 10
                // extract the integer part
                integer = shift_right128_round(res[0], res[1], upper_frac_bits);
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
     * @brief Performs right shift operation.
     * @param shift bits to shift
     * @return Temporary object with the result of the operation
    */
    template<typename T>
    __forceinline fixed_point128 operator>>(T shift) const noexcept {
        fixed_point128 temp(*this);
        return temp >>= static_cast<int32_t>(shift);
    }
    /**
     * @brief Performs left shift operation.
     * @param shift bits to shift
     * @return Temporary object with the result of the operation
    */
    template<typename T>
    __forceinline fixed_point128 operator<<(T shift) const noexcept {
        fixed_point128 temp(*this);
        return temp <<= static_cast<int32_t>(shift);
    }
    /**
     * @brief Add a value to this object
     * @param other Right hand side operand
     * @return This object.
    */
    FP128_INLINE fixed_point128& operator+=(const fixed_point128& other) noexcept {
        bool result_has_different_sign = false;
        fixed_point128 temp = other;

        // different sign: invert the sign for other and subtract
        if (other.sign != sign) {
            temp.sign ^= 1;
            result_has_different_sign = (sign) ? temp < *this : temp > *this;
            twos_complement128(temp.low, temp.high);
        }
        //add the other value
        const uint8_t carry = _addcarry_u64(0, low, temp.low, &low);
        high += temp.high + carry;

        // if result is with a different sign, invert it along with the sign.
        if (result_has_different_sign) {
            sign ^= 1;
            twos_complement128(low, high);
        }
        
        reset_sign_for_zero();
        return *this;
    }
    /**
     * @brief Add a value to this object
     * @param other Right hand side operand
     * @return This object.
    */
    template<typename T>
    FP128_INLINE fixed_point128& operator+=(const T& other) {
        return operator+=(fixed_point128(other));
    }
    /**
     * @brief Subtract a value to this object
     * @param other Right hand side operand
     * @return This object.
    */
    FP128_INLINE fixed_point128& operator-=(const fixed_point128& other) noexcept {
        bool result_has_different_sign = false;
        fixed_point128 temp = other;

        // same sign: invert the sign for other and subtract
        if (other.sign == sign) {
            result_has_different_sign = (sign) ? temp < *this : temp > *this;
            twos_complement128(temp.low, temp.high);
        }
        //add the other value
        const uint8_t carry = _addcarry_u64(0, low, temp.low, &low);
        high += temp.high + carry;

        // if result is with a different sign, invert it along with the sign.
        if (result_has_different_sign) {
            sign ^= 1;
            twos_complement128(low, high);
        }

        reset_sign_for_zero();
        return *this;
    }
    /**
     * @brief Subtract a value to this object
     * @param other Right hand side operand
     * @return This object.
    */
    template<typename T>
    __forceinline fixed_point128& operator-=(const T& other) noexcept {
        return operator-=(fixed_point128(other));
    }
    /**
     * @brief Multiplies a value to this object
     * @param other Right hand side operand
     * @return This object.
    */
    FP128_INLINE fixed_point128& operator*=(const fixed_point128& other) noexcept{
        // Temporary arrays to store the result. They are uninitialzied to get 10-50% extra performance.
        // Zero initialization is a 10% penalty and using a thread_local static varible lowers
        //  performance by >50%.

        uint64_t res[4]; // 256 bit of result
        uint64_t temp1[2], temp2[2];

        // multiply low QWORDs
        res[0] = _umul128(low, other.low, &res[1]);

        // multiply high QWORDs (overflow can happen)
        res[2] = _umul128(high, other.high, &res[3]);

        // multiply low this and high other
        temp1[0] = _umul128(low, other.high, &temp1[1]);
        uint8_t carry = _addcarry_u64(0, res[1], temp1[0], &res[1]);
        res[3] += _addcarry_u64(carry, res[2], temp1[1], &res[2]);

        // multiply high this and low other
        temp2[0] = _umul128(high, other.low, &temp2[1]);
        carry = _addcarry_u64(0, res[1], temp2[0], &res[1]);
        res[3] += _addcarry_u64(carry, res[2], temp2[1], &res[2]);

        // extract the bits from res[] keeping the precision the same as this object
        // shift result by F
        static constexpr int32_t index = F >> 6; // / 64;
        static constexpr int32_t lsb = F & FP128_MAX_VALUE_64(6); // bit within the 64bit data pointed by res[index]
        static constexpr uint64_t half = 1ull << (lsb - 1);       // used for rounding
        const bool need_rounding = (res[index] & half) != 0;

        // copy block #1 (lowest)
        low = __shiftright128(res[index], res[index + 1], lsb); // intrinsic is faster when shift is < 64
        high = __shiftright128(res[index+1], res[index + 2], lsb);

        if (need_rounding) {
            ++low; // low will wrap around to zero if overflowed
            high += low == 0;
        }
        // set the sign
        sign ^= other.sign;
        reset_sign_for_zero();
        return *this;
    }
    /**
     * @brief Multiplies a 64 bit value to this object
     * @param x Right hand side operand
     * @return This object.
    */
    FP128_INLINE fixed_point128& operator*=(uint64_t x) noexcept{
        uint64_t temp;

        // multiply low QWORDs
        low = _umul128(low, x, &temp);
        high = high * x + temp;
        reset_sign_for_zero();
        return *this;
    }
    /**
     * @brief Multiplies a value to this object
     * @param x Right hand side operand
     * @return This object.
    */
    template<typename T>
    __forceinline fixed_point128& operator*=(T x) noexcept {
        // floating point
        if constexpr (std::is_floating_point_v<T>) {
            return operator*=(fixed_point128(x));
        }
        // integers: convert to uint64 for a simpler operation.
        if constexpr (std::is_signed_v<T>) {
            // alway do positive multiplication
            if (x < 0) {
                x = -x;
                sign ^= 1;
            }
        }
        
        return operator*=(static_cast<uint64_t>(x));
    }
    /**
     * @brief Divide this object by x.
     * @param other Right hand side operator (denominator)
     * @return this object.
    */
    FP128_INLINE fixed_point128& operator/=(const fixed_point128& other) {
        bool need_rounding = false;
        // trivial case, this object is zero
        if (!*this)
            return *this;

        // exponent of 2, convert to a much faster shift operation
        if (1 == popcnt128(other.low, other.high)) {
            const auto expo = other.get_exponent();
            if (expo > 0)
                *this >>= (int32_t)expo;
            else if (expo < 0)
                *this <<= (int32_t)-expo;
        }
        // optimization for when dividing by an integer
        else if (other.is_int() && (uint64_t)other <= UINT64_MAX) {
            //*this *= fabs(reciprocal(other));
            uint64_t q[2]{};
            const uint64_t nom[2] = { low, high };
            const uint64_t denom = (uint64_t)other;
            uint64_t r;
            if (0 == div_64bit((uint64_t*)q, &r, (uint64_t*)nom, denom, 2)) {
                need_rounding = r > (denom >> 1);
                high = q[1];
                low = q[0];
            }
            else { // error
                FP128_FLOAT_DIVIDE_BY_ZERO_EXCEPTION;
            }
        }
        else {
            *this *= fabs(reciprocal(other));
        }
        //else {
        //    uint64_t q[4]{};
        //    const uint64_t nom[4] = {0, 0, low, high};
        //    const uint64_t denom[2] = {other.low, other.high};

        //    if (0 == div_32bit((uint32_t*)q, nullptr, (uint32_t*)nom, (uint32_t*)denom, 2ll * array_length(nom), 2ll * array_length(denom))) {
        //        static constexpr uint64_t half = 1ull << (I - 1);  // used for rounding
        //        need_rounding = (q[0] & half) != 0;
        //        // result in q needs to shifted left by F (F bits were added to the right)
        //        // shifting right by 128-F is simpler.
        //        if constexpr (I == 64) {
        //            high = q[1];
        //            low = q[0];
        //        }
        //        else if constexpr (I < 64) {
        //            high = __shiftright128(q[1], q[2], I);
        //            low = __shiftright128(q[0], q[1], I);
        //        }
        //    }
        //    else { // error
        //        FP128_FLOAT_DIVIDE_BY_ZERO_EXCEPTION;
        //    }
        //}

        if (need_rounding) {
            ++low;
            high += low == 0;
        }
        sign ^= other.sign;
        reset_sign_for_zero();
        return *this;
    }
    /**
     * @brief Divide this object by x.
     * @param x Denominator.
     * @return This object.
    */
    template<typename T>
    __forceinline fixed_point128& operator/=(T x) {
        if constexpr (std::is_floating_point<T>::value) {
            return operator/=(static_cast<double>(x));
        }

        return operator/=(fixed_point128(x));
    }
    /**
     * @brief Divide this object by x.
     * @param x Denominator.
     * @return This object.
    */
    template<>
    FP128_INLINE fixed_point128& operator/=<double>(double x) {
        if (0 == x) FP128_FLOAT_DIVIDE_BY_ZERO_EXCEPTION;

        // Simple and common case, the value is an exponent of 2
        // Convert to a much faster shift operation
        const Double d(x);
        if (0 == d.f) {
            sign ^= d.s;
            const int32_t e = static_cast<int32_t>(d.e) - 1023;
            return (e >= 0) ? *this >>= e : *this <<= -e;
        }

        // normal division
        return *this /= fixed_point128(x);
    }
    /**
     * @brief %= operator
     * @param other Modulo operand.
     * @return This object.
    */
    FP128_INLINE fixed_point128& operator%=(const fixed_point128& other) {
        // trivial case, this object is zero
        if (!*this)
            return *this;

        // simple case, both are integers (fractions is zero)
        if (is_int() && other.is_int()) {
            if (!other)
                FP128_FLOAT_DIVIDE_BY_ZERO_EXCEPTION;
            return operator=(static_cast<int64_t>(*this) % static_cast<int64_t>(other));
        }
        // num or denom are fractions
        // x mod y =  x - y * floor(x/y)
        // do the division in with positive numbers
        const uint64_t nom[4] = { 0, 0, low, high };
        const uint64_t denom[2] = { other.low, other.high };
        uint64_t q[4]{};
        uint64_t r[2]{}; // same size as the denominator
        if (0 == div_32bit((uint32_t*)q, (uint32_t*)r, (uint32_t*)nom, (uint32_t*)denom, 2ll * array_length(nom), 2ll * array_length(denom))) {
            fixed_point128 x_div_y; // x/root. 
            x_div_y.high = (q[2] << upper_frac_bits) | (q[1] >> I);
            x_div_y.low = (q[1] << upper_frac_bits) | (q[0] >> I);
            x_div_y.sign = sign ^ other.sign;
            // Integer result - remainder is zero.
            if (x_div_y.is_int()) {
                *this = 0;
            }
            // Fraction result - remainder is non zero.
            else {
                if constexpr (FP128_CPP_STYLE_MODULO) {
                    const bool this_was_positive = is_positive();
                    *this -= other * floor(x_div_y);
                    // this was positive, return positive modulo
                    if (this_was_positive) {
                        if (is_negative())
                            *this += other.is_positive() ? other : -other;
                    }
                    // this was negative, return negative modulo
                    else {
                        if (is_positive())
                            *this += other.is_positive() ? -other : other;
                    }
                }
                else {
                    *this -= other * floor(x_div_y);
                    // common case (fractions and integers) where one of the values is negative
                    if (sign != other.sign) {
                        // the remainder + denominator
                        *this += other;
                    }
                }
            }

            // Note if signs are the same, for nom/denom, the result keeps the sign.
            reset_sign_for_zero();
        }
        else { // error
            FP128_FLOAT_DIVIDE_BY_ZERO_EXCEPTION;
        }
        return *this;
    }
    template<typename T>
    __forceinline fixed_point128& operator%=(T other) {
        return operator%=(fixed_point128(other));
    }
    /**
     * @brief Shift right this object.
     * @param shift Bits to shift. Negative or very high values cause undefined behavior.
     * @return This object.
    */
    __forceinline fixed_point128& operator>>=(int32_t shift) noexcept {
        if (shift < 1)
            return *this;
        // 0-64 bit shift - most common
        if (shift <= 64) {
            low = shift_right128_round(low, high, shift);
            high >>= shift;
        }
        else if (shift >= 128) {
            low = high = 0;
        }
        else if (shift >= 64) {
            low = shift_right64_round(high, shift - 64);
            high = 0;
        }
        reset_sign_for_zero();
        return *this;
    }
    /**
     * @brief Shift left this object.
     * @param shift Bits to shift. Negative or very high values cause undefined behavior. 
     * @return This object.
    */
    __forceinline fixed_point128& operator<<=(int32_t shift) noexcept {
        if (shift < 1)
            return *this;
        if (shift <= 64) {
            high = shift_left128(low, high, shift);
            low <<= shift;
        }
        else if (shift < 128) {
            high = low << (shift - 64);
            low = 0;
        }
        else {
            low = high = 0;
        }
        reset_sign_for_zero();
        return *this;
    }
    /**
     * @brief Bitwise AND=
     * @param rhs AND mask.
     * @return This object.
    */
    __forceinline fixed_point128& operator&=(const fixed_point128& rhs) noexcept {
        low &= rhs.low;
        high &= rhs.high;
        sign &= rhs.sign;
        reset_sign_for_zero();
        return *this;
    }
    /**
     * @brief Bitwise AND=
     * @param rhs Right hand side operand
     * @return This object.
    */
    template<typename T>
    FP128_INLINE fixed_point128& operator&=(const T& rhs) {
        return operator&=(fixed_point128(rhs));
    }
    /**
     * @brief Bitwise OR=
     * @param other OR mask.
     * @return This object.
    */
    __forceinline fixed_point128& operator|=(const fixed_point128& other) noexcept {
        low |= other.low;
        high |= other.high;
        sign |= other.sign;
        return *this;
    }
    /**
     * @brief Bitwise OR=
     * @param rhs Right hand side operand
     * @return This object.
    */
    template<typename T>
    FP128_INLINE fixed_point128& operator|= (const T & rhs){
        return operator|=(fixed_point128(rhs));
    }
    /**
     * @brief Bitwise XOR=
     * @param other XOR mask.
     * @return This object.
    */
    __forceinline fixed_point128& operator^=(const fixed_point128& other) noexcept {
        low ^= other.low;
        high ^= other.high;
        sign ^= other.sign;
        return *this;
    }
    /**
     * @brief Bitwise XOR=
     * @param rhs Right hand side operand
     * @return This object.
    */
    template<typename T>
    FP128_INLINE fixed_point128& operator^= (const T & rhs){
        return operator^=(fixed_point128(rhs));
    }
    /**
     * @brief Prefix ++ operation (++a)
     * @return This object.
    */
    __forceinline fixed_point128& operator++() noexcept {
        *this += one();
        return *this;
    }
    /**
     * @brief Postfix ++ operation (a++)
     * @return This object.
    */
    __forceinline fixed_point128 operator++(int32_t) noexcept {
        fixed_point128 temp(*this);
        ++*this; // call the prefix implementation
        return temp;
    }
    /**
     * @brief Prefix -- operation (--a)
     * @return This object.
    */
    __forceinline fixed_point128& operator--() noexcept {
        *this -= one();
        return *this;
    }
    /**
     * @brief Postfix -- operation (a--)
     * @return This object.
    */
    __forceinline fixed_point128 operator--(int32_t) noexcept {
        fixed_point128 temp(*this);
        --*this; // call the prefix implementation
        return temp;
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
     * @brief Bitwise not (~).
    */
    __forceinline fixed_point128 operator~() const noexcept {
        fixed_point128 temp(*this);
        temp.high = ~high;
        temp.low = ~low;
        temp.reset_sign_for_zero();
        return temp;
    }
    /**
     * @brief Unary +. Returns a copy of the object.
    */
    __forceinline fixed_point128 operator+() const noexcept {
        fixed_point128 temp(*this);
        return temp;
    }
    /**
     * @brief Unary -. Returns a copy of the object with sign inverted.
    */
    __forceinline fixed_point128 operator-() const noexcept {
        fixed_point128 temp(*this);
        temp.sign ^= 1;
        temp.reset_sign_for_zero();
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
    __forceinline bool operator==(const fixed_point128& other) const noexcept {
        return sign == other.sign && high == other.high && low == other.low;
    }
    /**
     * @brief Compare logical/bitwise equal.
     * @param other Righthand operand
     * @return True if this and other are equal.
    */
    template<typename T>
    __forceinline bool operator==(T other) const noexcept {
        return *this == fixed_point128(other);
    }
    /**
     * @brief Return true when objects are not equal. Can be used as logical XOR.
     * @param other Righthand operand.
     * @return True of not equal.
    */
    __forceinline bool operator!=(const fixed_point128& other) const noexcept {
        return sign != other.sign || high != other.high || low != other.low;
    }
    /**
     * @brief Return true when objects are not equal. Can be used as logical XOR.
     * @param other Righthand operand.
     * @return True of not equal.
    */
    template<typename T>
    __forceinline bool operator!=(T other) const noexcept {
        return *this != fixed_point128(other);
    }
    /**
     * @brief Return true if this object is small than the other
     * @param other Righthand operand.
     * @return True when this object is smaller.
    */
    __forceinline bool operator<(const fixed_point128& other) const noexcept {
        // signs are different
        if (sign != other.sign)
            return sign > other.sign; // true when sign is 1 and other.sign is 0

        // MSB is the same, check the LSB
        if (high == other.high)
            return (sign) ? low > other.low : low < other.low;

        return (sign) ? high > other.high : high < other.high;
    }
    /**
     * @brief Return true if this object is small than the other
     * @param other Righthand operand.
     * @return True when this object is smaller.
    */
    template<typename T>
    __forceinline bool operator<(T other) const noexcept {
        return operator<(fixed_point128(other));
    }
    /**
     * @brief Return true this object is small or equal than the other
     * @param other Righthand operand.
     * @return True when this object is smaller or equal.
    */
    __forceinline bool operator<=(const fixed_point128& other) const noexcept {
        return !(*this > other);
    }
    /**
     * @brief Return true this object is small or equal than the other
     * @param other Righthand operand.
     * @return True when this object is smaller or equal.
    */
    template<typename T>
    __forceinline bool operator<=(T other) const noexcept {
        return !(*this > fixed_point128(other));
    }
    /**
     * @brief Return true this object is larger than the other
     * @param other Righthand operand.
     * @return True when this objext is larger.
    */
    __forceinline bool operator>(const fixed_point128& other) const noexcept {
        // signs are different
        if (sign != other.sign)
            return sign < other.sign; // true when sign is 0 and other.sign is 1

        // MSB is the same, check the LSB
        if (high == other.high)
            return (sign) ? low < other.low : low > other.low;

        return (sign) ? high < other.high : high > other.high;
    }
    /**
     * @brief Return true this object is larger than the other
     * @param other Righthand operand.
     * @return True when this objext is larger.
    */
    template<typename T>
    __forceinline bool operator>(T other) const noexcept {
        return *this > fixed_point128(other);
    }
    /**
     * @brief Return true this object is larger or equal than the other
     * @param other Righthand operand.
     * @return True when this objext is larger or equal.
    */
    __forceinline bool operator>=(const fixed_point128& other) const noexcept {
        return !(*this < other);
    }
    /**
     * @brief Return true this object is larger or equal than the other
     * @param other Righthand operand.
     * @return True when this objext is larger or equal.
    */
    template<typename T>
    __forceinline bool operator>=(T other) const noexcept {
        return !(*this < fixed_point128(other));
    }
    //
    // useful public functions
    //
    /**
     * @brief Returns true if the value is an int (fraction is zero)
     * @return True when the fraction is zero.
    */
    __forceinline bool is_int() const noexcept
    {
        return 0 == low && 0 == (high << I);
    }
    /**
     * @brief Returns true if the value positive (incuding zero)
     * @return True when the the value positive
    */
    __forceinline bool is_positive() const noexcept
    {
        return 0 == sign;
    }
    /**
     * @brief Returns true if the value negative (smaller than zero)
     * @return True when the the value negative
    */
    __forceinline bool is_negative() const noexcept
    {
        return 1 == sign;
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
     * @brief get a specific bit within the 128 fixed point data
     * @param bit bit to get [0,127]
     * @return 0 or 1. Undefined when bit > 127
    */
    __forceinline int32_t get_bit(uint32_t bit) const noexcept
    {
        if (bit < 64) {
            return FP128_GET_BIT(low, bit);
        }
        return FP128_GET_BIT(high, bit-64);
    }
    /**
     * @brief Returns the exponent of the object - like the base 2 exponent of a floating point
     * A value of 2.1 would return 1, values in the range [0.5,1.0) would return -1.
     * @return Exponent of the number
    */
    __forceinline int32_t get_exponent() const noexcept
    {
        const int32_t s = static_cast<int32_t>(lzcnt128(*this));
        return I - 1 - s;
    }
    /**
     * @brief Returns an instance of fixed_point128 with the value of pi
     * @return pi
    */
    __forceinline static const fixed_point128& pi() noexcept {
        static const fixed_point128 pi = "3.14159265358979323846264338327950288419716939937510"; // 50 first digits of pi
        return pi;
    }
    /**
     * @brief Returns an instance of fixed_point128 with the value of pi * 2
     * @return pi * 2
    */
    __forceinline static const fixed_point128& pi2() noexcept {
        static const fixed_point128 pi2 = "6.28318530717958647692528676655900576839433879875021"; // 50 first digits of pi * 2
        return pi2;
    }
    /**
     * @brief Returns an instance of fixed_point128 with the value of pi / 2
     * @return pi / 2
    */
    __forceinline static const fixed_point128& half_pi() noexcept {
        static const fixed_point128 half_pi = "1.57079632679489661923132169163975144209858469968755"; // 50 first digits of pi / 2
        return half_pi;
    }    
    /**
     * @brief Returns an instance of fixed_point128 with the value of the golden ratio
     * @return golden ratio constant
    */
    __forceinline static const fixed_point128& golden_ratio() noexcept {
        static const fixed_point128 golden_ratio = "1.6180339887498948482045868343656381177203"; // 40 first digits of the golden ratio
        return golden_ratio;
    }
    /**
     * @brief Returns an instance of fixed_point128 with the value of e
     * @return e
    */
    __forceinline static const fixed_point128& e() noexcept {
        static const fixed_point128 e = "2.71828182845904523536028747135266249775724709369"; // 50 first digits of e
        return e;
    }
    /**
     * @brief Returns an instance of fixed_point128 with the value of sqrt(2)
     * @return e
    */
    __forceinline static const fixed_point128& sqrt_2() noexcept {
        static const fixed_point128 sqrt_2 = "1.41421356237309504880168872420969807856967187537"; // 50 first digits of e
        return sqrt_2;
    }
    /**
     * @brief Return an instance of fixed_point128 with the value of 1
     * @return 1
    */
    __forceinline static const fixed_point128& one() noexcept {
        static const fixed_point128 one = 1;
        return one;
    }
    /**
     * @brief Return an instance of fixed_point128 with the value of 0.5
     * @return 1
    */
    __forceinline static const fixed_point128& half() noexcept
    {
        static const fixed_point128 half = 0.5;
        return half;
    }
    /**
     * @brief Return an instance of fixed_point128 with the smallest positive value possible
     * @return 1
    */
    __forceinline static const fixed_point128& epsilon() noexcept {
        static const fixed_point128 epsilon(1, 0, 0);
        return epsilon;
    }
private:
    /**
     * @brief Set the sign to 0 when both low and high are zero, i.e. avoid having negative zero value
    */
    __forceinline void reset_sign_for_zero() {
        sign &= (0 != low || 0 != high);
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
    friend __forceinline fixed_point128 operator+(fixed_point128 lhs, const T& rhs) noexcept {
        return lhs += rhs;
    }
    /**
     * @brief subtracts the right hand side operand to this object to and returns the result.
     * @param lhs left hand side operand
     * @param rhs Right hand side operand
     * @return The fixed_point128 result
    */
    template<typename T>
    friend __forceinline fixed_point128 operator-(fixed_point128 lhs, const T& rhs) noexcept {
        return lhs -= rhs;
    }
    /**
     * @brief Multiplies the right hand side operand with this object to and returns the result.
     * @param lhs left hand side operand
     * @param rhs Right hand side operand
     * @return The fixed_point128 result
    */
    template<typename T>
    friend __forceinline fixed_point128 operator*(fixed_point128 lhs, const T& rhs) noexcept {
        return lhs *= rhs;
    }
    /**
     * @brief Divides this object by the right hand side operand and returns the result.
     * @param lhs left hand side operand
     * @param rhs Right hand side operand
     * @return The fixed_point128 result
    */
    template<typename T>
    friend __forceinline fixed_point128 operator/(fixed_point128 lhs, const T& rhs) {
        return lhs /= rhs;
    }
    /**
     * @brief Calculates modulo.
     * @param lhs left hand side operand
     * @param rhs Right hand side operand
     * @return The fixed_point128 result
    */
    template<typename T>
    friend __forceinline fixed_point128 operator%(fixed_point128 lhs, const T& rhs) {
        return lhs %= rhs;
    }
    
    //
    // Binary math operators
    //

    /**
     * @brief Performs bitwise AND (&)
     * @param lhs left hand side operand
     * @param rhs Right hand side operand
     * @return The fixed_point128 result
    */
    template<typename T>
    friend __forceinline fixed_point128 operator&(fixed_point128 lhs, const T& rhs) {
        return lhs &= rhs;
    }
    /**
     * @brief Performs bitwise OR (|)
     * @param lhs left hand side operand
     * @param rhs Right hand side operand
     * @return The fixed_point128 result
    */
    template<typename T>
    friend __forceinline fixed_point128 operator|(fixed_point128 lhs, const T& rhs) {
        return lhs |= rhs;
    }
    /**
     * @brief Performs bitwise XOR (^)
     * @param lhs left hand side operand
     * @param rhs Right hand side operand
     * @return The fixed_point128 result
    */
    template<typename T>
    friend __forceinline fixed_point128 operator^(fixed_point128 lhs, const T& rhs) {
        return lhs ^= rhs;
    }

    //
    // Floating point style functions, implemented as friend functions
    //

    /**
     * @brief Returns the absolute value
     * @param x Fixed_point128 object
     * @return A copy of x with sign removed
    */
    friend __forceinline fixed_point128 fabs(const fixed_point128& x) noexcept
    {
        fixed_point128 temp = x;
        temp.sign = 0;
        return temp;
    }
    /**
     * @brief Performs the floor() function, similar to libc's floor(), rounds down towards -infinity.
     * @param x Input value
     * @return A fixed_point128 holding the integer value. Overflow is not reported.
    */
    friend FP128_INLINE fixed_point128 floor(const fixed_point128& x) noexcept
    {
        if (x.is_int()) return x;

        fixed_point128 res;
        res.high = x.high & x.int_mask;
        res.high += (uint64_t)x.is_negative() << upper_frac_bits;
        res.low = 0;
        res.sign = x.sign;
        return res;
    }

    /**
     * @brief Performs the ceil() function, similar to libc's ceil(), rounds up towards infinity.
     * @param x Input value
     * @return A fixed_point128 holding the integer value. Overflow is not reported.
    */
    friend FP128_INLINE fixed_point128 ceil(const fixed_point128& x) noexcept
    {
        if (x.is_int()) return x;

        fixed_point128 res;
        res.high = x.high & x.int_mask;
        res.high += (uint64_t)x.is_positive() << upper_frac_bits;
        res.low = 0;
        res.sign = x.sign;
        return res;
    }
    /**
     * @brief Rounds towards zero
     * @param x Value to truncate
     * @return Integer value, rounded towards zero.
    */
    friend __forceinline fixed_point128 trunc(const fixed_point128& x) noexcept
    {
        return fixed_point128(0, x.high & x.int_mask, x.sign);
    }
    /**
     * @brief Rounds towards the nearest integer.
     * The halfway value (0.5) is rounded away from zero.
     * @param x Value to round
     * @return Integer value, rounded towards the nearest integer.
    */
    friend FP128_INLINE fixed_point128 round(const fixed_point128& x) noexcept
    {
        // save the sign
        auto sign = x.sign;
        fixed_point128 res = floor(fabs(x) + fixed_point128::half());

        // restore the sign
        res.sign = sign;
        return res;
    }
    /**
     * @brief Performs the fmod() function, similar to libc's fmod(), returns the remainder of a division x/root.
     * @param x Numerator
     * @param y Denominator
     * @return A fixed_point128 holding the modulo value.
    */
    friend __forceinline fixed_point128 fmod(const fixed_point128& x, const fixed_point128& y)
    {
        return x % y;
    }

    /**
     * @brief Split into integer and fraction parts.
     * Both results carry the sign of the input variable.
     * @param x Input value
     * @param iptr Pointer to fixed_point128 holding the integer part of x.
     * @return The fraction part of x. Undefined when iptr is nullptr.
    */
    friend FP128_INLINE fixed_point128 modf(const fixed_point128& x, fixed_point128* iptr) noexcept
    {
        if (iptr == nullptr)
            return 0;

        iptr->high = x.high & x.int_mask; // lose the fraction
        iptr->low = 0;
        iptr->sign = x.sign;

        fixed_point128 res = x;
        res.high &= ~x.int_mask; // lose the integer part
        return res;
    }
    /**
     * @brief Calculates the left zero count of value x, ignoring the sign.
     * @param x input value.
     * @return lzc (uint32_t) of th result.
    */
    friend __forceinline uint64_t lzcnt128(const fixed_point128& x) noexcept
    {
        return (x.high != 0) ? __lzcnt64(x.high) : 64 + __lzcnt64(x.low);
    }
    /**
     * @brief Calculates the square root using Newton's method.
     * Based on the book "Math toolkit for real time programming" by Jack W. Crenshaw 
     * @param x Value to calculate the root of
     * @param iterations how many iterations to perform (more is more accurate). Sensible values are 0-5.
     * @return Square root of (x), zero when x <= 0.
    */
    friend FP128_INLINE fixed_point128 sqrt(const fixed_point128& x, uint32_t iterations = 3) noexcept
    {
        static const fixed_point128 factor = "0.70710678118654752440084436210484903928483593768847403658833981"; // sqrt(2) / 2
        if (x.sign || !x)
            return 0;

        // normalize the input to the range [0.5, 1)
        int32_t expo = x.get_exponent() + 1;
        fixed_point128 norm_x = (expo > 0) ? x >> expo : x << -expo;
        
        // use existing HW to provide an excellent first estimate.
        fixed_point128 root = ::sqrt(static_cast<double>(norm_x));

        // iterate several times via Newton's method
        //                  X
        //   Xn+1 = 0.5 * (---- + Xn )
        //                  Xn
        for (auto i = iterations; i != 0; --i) {
            root = (norm_x * reciprocal(root) + root) >> 1;
        }

        if (expo & 1) {
            root *= factor;
            ++expo;
        }

        // Denormalize the result
        if (expo > 0)
            root <<= (expo + 1) / 2;
        else if (expo < 0)
            root >>= -expo / 2;

        return root;
    }
    /**
     * @brief Factorial reciprocal (inverse). Calculates 1 / x!
     * Maximum value of x that may produce non zero values is 34. 
     * This value depends on the amount of fraction bits.
     * @param x Input value
     * @param res Result of the function
    */
    friend FP128_INLINE void fact_reciprocal(int x, fixed_point128& res) noexcept
    {
        static const fixed_point128 c[] = {
            "1",                                         // 1 /  0!
            "1",                                         // 1 /  1!
            "0.5",                                       // 1 /  2!
            "0.166666666666666666666666666666666666667", // 1 /  3!
            "0.041666666666666666666666666666666666667", // 1 /  4!
            "0.008333333333333333333333333333333333333", // 1 /  5!
            "0.001388888888888888888888888888888888889", // 1 /  6!
            "0.000198412698412698412698412698412698413", // 1 /  7!
            "0.000024801587301587301587301587301587302", // 1 /  8!
            "0.000002755731922398589065255731922398589", // 1 /  9!
            "0.000000275573192239858906525573192239859", // 1 / 10!
            "0.000000025052108385441718775052108385442", // 1 / 11!
            "0.000000002087675698786809897921009032120", // 1 / 12!
            "0.000000000160590438368216145993923771702", // 1 / 13!
            "0.000000000011470745597729724713851697979", // 1 / 14!
            "0.000000000000764716373181981647590113199", // 1 / 15!
            "0.000000000000047794773323873852974382075", // 1 / 16!
            "0.000000000000002811457254345520763198946", // 1 / 17!
            "0.000000000000000156192069685862264622164", // 1 / 18!
            "0.000000000000000008220635246624329716956", // 1 / 19!
            "0.000000000000000000411031762331216485848", // 1 / 20!
            "0.000000000000000000019572941063391261231", // 1 / 21!
            "0.000000000000000000000889679139245057329", // 1 / 22!
            "0.000000000000000000000038681701706306840", // 1 / 23!
            "0.000000000000000000000001611737571096118", // 1 / 24!
            "0.000000000000000000000000064469502843845", // 1 / 25!
            "0.000000000000000000000000002479596263225", // 1 / 26!
            "0.000000000000000000000000000091836898638", // 1 / 27!
            "0.000000000000000000000000000003279889237", // 1 / 28!
            "0.000000000000000000000000000000113099628", // 1 / 29!
            "0.000000000000000000000000000000003769988", // 1 / 30!
            "0.000000000000000000000000000000000121613", // 1 / 31!
            "0.000000000000000000000000000000000003800", // 1 / 32!
            "0.000000000000000000000000000000000000115", // 1 / 33!
            "0.000000000000000000000000000000000000003"  // 1 / 34!
        };
        constexpr int series_len = array_length(c);

        if (x >= 0 && x < series_len) {
            res = c[x];
        }
        else {
            res = 0;
        }
    }
    /**
     * @brief Calculates the reciprocal of a value. y = 1 / x
     * Using newton iterations: Yn+1 = Yn(2 - x * Yn)
     * @param x Input value
     * @return 1 / x. Returns zero on overflow or division by zero
    */
    friend FP128_INLINE fixed_point128 reciprocal(const fixed_point128& x) noexcept
    {
        static const fixed_point128 one = 1, two = 2;
        static const fixed_point128 xy_max = one + (fixed_point128::epsilon() << 2);
        static const fixed_point128 xy_min = one - (fixed_point128::epsilon() << 2);
        constexpr int max_iterations = 3;
        constexpr int debug = false;
        fixed_point128 y = 1.0 / static_cast<double>(x);

        if (!y)
            return y; 
        
        fixed_point128 xy, y_prev;
        // Newton iterations:
        int i = 0;
        for (; i < max_iterations && (y_prev != y) && (xy < xy_min || xy > xy_max); ++i) {
            y_prev = y;
            xy = x * y;
            y = y * (two - xy);
        }

        if constexpr (debug) {
            static int debug_max_iter = 0;
            if (i > debug_max_iter || i == max_iterations) {
                debug_max_iter = i;
                printf("reciprocal took %i iterations for %.10lf\n", i, static_cast<double>(x));
            }
        }

        return y;
    }
    /**
         * @brief Calculate the sine function over a limited range [-0.5pi, 0.5pi]
         * Using the Maclaurin series expansion, the formula is:
         *              x^3   x^5   x^7
         * sin(x) = x - --- + --- - --- + ...
         *               3!    5!    7!
         * @param x value in Radians in the range [-0.5pi, 0.5pi]
         * @return Sine of x
        */
    friend FP128_INLINE fixed_point128 sin1(fixed_point128 x) noexcept
    {
        static_assert(I >= 4, "fixed_point128 must have at least 4 integer bits to use sin1()!");
        assert(fabs(x) <= fixed_point128::half_pi());

        // first part of the series is just 'x'
        const fixed_point128 xx = x * x;
        fixed_point128 elem_denom, elem_nom = x;

        // compute the rest of the series, starting with: -(x^3 / 3!)
        for (int i = 3, sign = 1; ; i += 2, sign = 1 - sign) {
            elem_nom *= xx;
            fact_reciprocal(i, elem_denom);
            fixed_point128 elem = elem_nom * elem_denom; // next element in the series
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
    friend __forceinline fixed_point128 cos1(const fixed_point128& x) noexcept
    {
        static_assert(I >= 4, "fixed_point128 must have at least 4 integer bits to use sin1()!");
        static const fixed_point128& half_pi = fixed_point128::half_pi();
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
    friend fixed_point128 sin(fixed_point128 x) noexcept
    {
        static_assert(I >= 4, "fixed_point128 must have at least 4 integer bits to use sin()!");
        static const fixed_point128& pi = fixed_point128::pi();
        static const fixed_point128& pi2 = fixed_point128::pi2(); // 2 * pi
        static const fixed_point128& half_pi = fixed_point128::half_pi(); // pi / 2
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
     * @param x value in the range [-1,1]
     * @return Inverse sine of x
    */
    friend fixed_point128 asin(fixed_point128 x) noexcept
    {
        static const fixed_point128 eps = fixed_point128::epsilon() << 1;
        constexpr int max_iterations = 6;
        if (x < -1 || x > 1) return 0;

        //              sin(Xn) - a
        // Xn+1 = Xn - -------------
        //                cos(Xn)
        // where 'a' is the argument, each iteration will converge on the result if the initial
        //  estimate is close enough.
        auto sign = x.sign;
        x.sign = 0;
        fixed_point128 res = ::asin(static_cast<double>(x));
        for (int i = 0; i < max_iterations; ++i) {
            fixed_point128 e = (sin(res) - x) / cos(res);
            res -= e;
            if (fabs(e) <= eps)
                break;
        }

        res.sign = sign;
        return res;
    }
    /**
     * @brief Calculate the cosine function
     * Ultimately uses sin1() with a reduced range of [-pi/4, pi/4]
     * Sine's Maclaurin series converges faster than Cosine's.
     * @param x value in Radians
     * @return Cosine of x
    */
    friend fixed_point128 cos(fixed_point128 x) noexcept
    {
        static_assert(I >= 4, "fixed_point128 must have at least 4 integer bits to use cos()!");
        static const fixed_point128& pi = fixed_point128::pi();
        static const fixed_point128& pi2 = fixed_point128::pi2(); // 2 * pi
        static const fixed_point128& half_pi = fixed_point128::half_pi(); // pi / 2
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
     * @param x value in the range [-1,1]
     * @return Inverse cosine of x
    */
    friend fixed_point128 acos(fixed_point128 x) noexcept
    {
        static const fixed_point128 eps = fixed_point128::epsilon() << 1;
        constexpr int max_iterations = 6;
        if (x < -1 || x > 1) return 0;
        //              cos(Xn) - a           a - cos(Xn)
        // Xn+1 = Xn - ------------- = Xn -  ------------
        //                -sin(Xn)              sin(Xn)
        // where 'a' is the argument, each iteration will converge on the result if the initial
        //  estimate is close enough.
        fixed_point128 res = ::acos(static_cast<double>(x));
        for (int i = 0; i < max_iterations; ++i) {
            fixed_point128 cos_xn = cos(res);
            fixed_point128 sin_xn = sin(res);
            fixed_point128 e = (x - cos_xn) / sin_xn;
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
    friend FP128_INLINE fixed_point128 tan(fixed_point128 x) noexcept
    {
        static_assert(I >= 4, "fixed_point128 must have at least 4 integer bits to use tan()!");
        return sin(x) / cos(x);
    }
    /**
     * @brief Calculate the inverse tangent function
     * @param x value
     * @return Arctangent of x
    */
    friend fixed_point128 atan(fixed_point128 x) noexcept
    {
        // constants for segmentation
        static const fixed_point128& half_pi = fixed_point128::half_pi(); // pi / 2
        static const fixed_point128 eps = fixed_point128::epsilon() << 1;
        bool comp = false;
        constexpr int max_iterations = 6;

        // make argument positive, save the sign
        auto sign = x.sign;
        x.sign = 0;

        // limit argument to 0..1
        if (x > 1) {
            comp = true;
            x = reciprocal(x);
        }

        // initial step uses the CRT function.
        fixed_point128 res = ::atan(static_cast<double>(x));
        //
        // Xn+1 =  Xn - cos(Xn) * ( sin(Xn) - a * cos(Xn))
        // 
        // where 'a' is the argument, each iteration will converge on the result if the initial
        //  estimate is close enough.
        for (int i = 0; i < max_iterations; ++i) {
            fixed_point128 cos_xn = cos(res);
            fixed_point128 sin_xn = sin(res);
            fixed_point128 e = cos_xn * (sin_xn - x * cos_xn); // this is the iteration estimated error
            res -= e;
            if (fabs(e) <= eps)
                break;
        }

        // restore complement if needed
        if (comp)
            res = half_pi - res;
        // restore sign if needed
        res.sign = sign;
        return res;
    }
    /**
     * @brief Calculates the exponent of x: e^x
     * Using the Maclaurin series expansion, the formula is:
     *                x^1     x^2     x^3
     * exp(x) = 1  +  ---  +  ---  +  --- + ...
     *                 1!      2!      3!
     *
     * The Maclaurin series will quickly overflow as x's power increases rapidly.
     * Using the equality e^x = e^ix * e^fx
     * Where ix is the integer part of x and fx is the fraction part.
     * ix is computed via multiplications which won't overflow if the result value can be held.
     * fx is computed via Maclaurin series expansion, but since fx < 1, it won't overflow.
     * @param x A number specifying a power. 
     * @return Exponent of x
    */
    friend FP128_INLINE fixed_point128 exp(const fixed_point128& x) noexcept
    {
        static const fixed_point128 e = fixed_point128::e();
        fixed_point128 _ix, exp_ix; // integer part of x
        fixed_point128 fx = modf(fabs(x), &_ix);
        uint64_t ix = static_cast<uint64_t>(_ix); // 64 bit is an overkill to hold the exponent
        
        // compute e^ix (integer part of x)
        if (ix > 0) {
            exp_ix = 1;      // result
            fixed_point128 b = e; // value of e^1
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

        // compute e^fx (fraction part of x)
        // first and second elements of the series
        fixed_point128 exp_fx = fixed_point128::one() + fx;
        fixed_point128 elem_denom, elem_nom = fx;

        for (int i = 2; ; ++i) {
            elem_nom *= fx;
            fact_reciprocal(i, elem_denom);
            fixed_point128 elem = elem_nom * elem_denom;
            if (!elem)
                break;
            exp_fx += elem; // next element in the series
        }

        fixed_point128 res = exp_ix * exp_fx;
        return (x.is_positive()) ? res : reciprocal(res);
    }
    /**
     * @brief Calculates the exponent of x and reduces 1 from the result: (e^x) - 1
     * @param x A number specifying a power.
     * @return Exponent of x
    */
    friend FP128_INLINE fixed_point128 expm1(const fixed_point128& x) noexcept
    {
        return exp(x) - fixed_point128::one();
    }
    /**
     * @brief Computes 2 to the power of x
     * @param x Exponent value
     * @return 2^x
    */
    friend FP128_INLINE fixed_point128 exp2(const fixed_point128& x) noexcept
    {
        //
        // Based on exponent law: (x^n)^m = x^(m*n)
        // Convert the exponent x (function parameter) to produce an exponent that will work with exp()
        // y = log(2) 
        // 2^x = e^(y*x) = exp(y*x)
        //
        static const fixed_point128 lan2 = "0.693147180559945309417232121458176575";
        return exp(x * lan2);
    }
    /**
     * @brief Computes x to the power of y
     * @param x Base value, must be positive
     * @param y Exponent value
     * @return x^y
    */
    friend FP128_INLINE fixed_point128 pow(const fixed_point128& x, const fixed_point128& y) noexcept
    {
        //
        // Based on exponent law: (x^n)^m = x^(m * n)
        // Convert the exponent y (function parameter) to produce an exponent that will work with exp()
        // z = log(x) 
        // pow(x, y) = x^y = e^(y * z) = exp(y * z)
        //
        if (x.is_negative()) return 0;

        fixed_point128 lan_x = log(x);
        if (!lan_x) 
            return lan_x;

        return exp(y * lan_x);
    }
    /**
     * @brief Calculates the Log base 2 of x: log2(x)
     * @param x The number to perform log2 on.
     * @return log2(x)
    */
    friend FP128_INLINE fixed_point128 log2(fixed_point128 x) noexcept
    {
        static const fixed_point128 two(2); 
        fixed_point128 b = fixed_point128::one() >> 1;
        fixed_point128 y = 0;

        if (!x.is_positive() || x.is_zero()) {
            return fixed_point128(UINT64_MAX, UINT64_MAX, 1); // represents negative infinity
        }

        // bring x to the range [1,2)
        while (x < fixed_point128::one()) {
            x <<= 1;
            --y;
        }

        while (x >= two) {
            x >>= 1;
            ++y;
        }

        // x is an exponent of 2.
        if (x == fixed_point128::one())
            return y;

        fixed_point128 z = x;
        for (size_t i = 0; i < fixed_point128::F; ++i) {
            z *= z;
            if (z >= two) {
                z >>= 1;
                y += b;
            }
            b >>= 1;
        }

        return y;
    }
    /**
     * @brief Calculates the natural Log (base e) of x: log(x)
     * @param x The number to perform log on.
     * @return log(x)
    */
    friend FP128_INLINE fixed_point128 log(fixed_point128 x) noexcept
    {
        static const fixed_point128 inv_log2_e = "0.693147180559945309417232121458176575";
        fixed_point128 y = log2(x);
        return y * inv_log2_e;
    }
    /**
     * @brief Calculates the natural Log (base e) of 1 + x: log(1 + x)
     * @param x The number to perform log on.
     * @return log1p(x)
    */
    friend FP128_INLINE fixed_point128 log1p(fixed_point128 x) noexcept
    {
        return log(fixed_point128::one() + x);
    }
    /**
     * @brief Calculates Log base 10 of x: log10(x)
     * @param x The number to perform log on.
     * @return log10(x)
    */
    friend FP128_INLINE fixed_point128 log10(fixed_point128 x) noexcept
    {
        static const fixed_point128 inv_log2_10 = "0.301029995663981195213738894724493068";
        fixed_point128 y = log2(x);
        return y * inv_log2_10;
    }
    /**
     * @brief Calculates Log base 2 of x as an integer ignoring the sign of x.
     * Similar to: floor(log2(fabs(x)))
     * @param x The number to perform log on.
     * @return logb(x)
    */
    friend FP128_INLINE fixed_point128 logb(fixed_point128 x) noexcept
    {
        return x.get_exponent();
    }
    /*
    static int div_64bit_test(uint64_t* q, uint64_t* r, const uint64_t* u, const uint64_t* v, int m, int n) noexcept
    {
        constexpr uint128_t b(0, 1); // Number base (64 bits).
        constexpr uint128_t mask(UINT64_MAX, 0);   // 64 bit mask
        uint64_t* un, * vn;                // Normalized form of u, v.
        uint128_t qhat;                     // Estimated quotient digit.
        uint128_t rhat;                     // A remainder.
        uint128_t p;                        // Product of two digits.
        //int128_t t, k;                      // Temporary variables
        int64_t t, k;                      // Temporary variables
        int32_t i, j;                      // Indexes
        // disable various warnings, some are bogus in VS2022.
        // the below code relies on the implied truncation (to 32 bit) of several expressions.
    #pragma warning(push)
    #pragma warning(disable: 6255)
    #pragma warning(disable: 4244)
    #pragma warning(disable: 6297)
    #pragma warning(disable: 6385)
    #pragma warning(disable: 6386)
    #pragma warning(disable: 26451)

    // shrink the arrays to avoid extra work on small numbers
        while (m > 0 && u[m - 1] == 0) --m;
        while (n > 0 && v[n - 1] == 0) --n;

        if (m < n || n <= 0 || v[n - 1] == 0)
            return 1; // Return if invalid param.

        // Take care of the case of a single-digit divisor here.
        if (n == 1)
            return div_64bit(q, r, u, v[0], m);

        // Normalize by shifting v left just enough so that its high-order
        // bit is on, and shift u left the same amount. We may have to append a
        // high-order digit on the dividend; we do that unconditionally.

        const int32_t s = __lzcnt64(v[n - 1]);             // 0 <= s <= 64.
        const int32_t s_comp = 64 - s;
        vn = (uint64_t*)_alloca(sizeof(uint64_t) * n);
        for (i = n - 1; i > 0; --i) {
            vn[i] = shift_left128(v[i - 1], v[i], s);
            //vn[i] = (v[i] << s) | ((uint64_t)v[i - 1] >> s_comp);
        }
        vn[0] = v[0] << s;

        un = (uint64_t*)_alloca(sizeof(uint64_t) * (m + 1));
        un[m] = (uint128_t)u[m - 1] >> s_comp;
        for (i = m - 1; i > 0; --i)
            un[i] = (u[i] << s) | (uint64_t)((uint128_t)u[i - 1] >> s_comp);
        un[0] = u[0] << s;

        for (j = m - n; j >= 0; --j) {       // Main loop.
            // Compute estimate qhat of q[j].
            uint128_t d(un[j + n - 1], un[j + n]);
            qhat = d / vn[n - 1];
            rhat = d - qhat * vn[n - 1];
            //qhat = (un[j + n] * b + un[j + n - 1]) / vn[n - 1];
            //rhat = (un[j + n] * b + un[j + n - 1]) - qhat * vn[n - 1];
        again:
            if (qhat >= b || qhat * vn[n - 2] > (rhat << 64) + un[j + n - 2]) {
                --qhat;
                rhat += vn[n - 1];
                if (rhat < b) goto again;
            }

            // Multiply and subtract.
            k = 0;
            for (i = 0; i < n; ++i) {
                p = qhat * vn[i];
                t = (uint128_t)un[i + j] - k - (p & mask);
                un[i + j] = t;
                k = (p.high) - (t.high);
            }
            t = un[j + n] - k;
            un[j + n] = t;

            q[j] = qhat;          // Store quotient digit.
            if (t < 0) {          // If we subtracted too
                q[j] = q[j] - 1;  // much, add back.
                k = 0;
                for (i = 0; i < n; ++i) {
                    t = (uint64_t)un[i + j] + vn[i] + k;
                    un[i + j] = t;
                    k = t >> 32;
                }
                un[j + n] = un[j + n] + k;
            }
        } // End j.
        // If the caller wants the remainder, unnormalize
        // it and pass it back.
        if (r != NULL) {
            for (i = 0; i < n - 1; ++i)
                r[i] = (un[i] >> s) | ((uint64_t)un[i + 1] << s_comp);

            r[n - 1] = un[n - 1] >> s;
        }
        return 0;
    #pragma warning(pop)
    }
    */

}; //class fixed_point128


} //namespace fp128

#pragma warning(pop)

#endif // #ifndef FIXED_POINT128_H
