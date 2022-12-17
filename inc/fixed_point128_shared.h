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
#ifndef FIXED_POINT128_SHARED_H
#define FIXED_POINT128_SHARED_H

#include <intrin.h>
#include <string>
#include <cstdint>
#include <cstdlib>
#include <cassert>
#include <stdexcept>
#include <type_traits>

/***********************************************************************************
*                                  Build Options
************************************************************************************/
// Set to TRUE to disable function inlining - useful for profiling a specific function
#ifndef FP128_DISABLE_INLINE
#define FP128_DISABLE_INLINE FALSE
#endif 

#ifndef FP128_CPP_STYLE_MODULO 
#define FP128_CPP_STYLE_MODULO TRUE
#endif 

/***********************************************************************************
*                                  Macros
************************************************************************************/
#define FP128_ONE_SHIFT(x)          (1ull << (x))
#define FP128_MAX_VALUE_64(x)       (((uint64_t)-1ll) >> (64 - x))
#define FP128_GET_BIT(x, n)         (((x) >> (n)) & 1)
#define FP128_GET_BITS(x, b, count) (((x) >> (b)) & FP128_MAX_VALUE_64(count))
#define FP128_INT_DIVIDE_BY_ZERO_EXCEPTION   throw std::logic_error("Integer divide by zero!")
#define FP128_FLOAT_DIVIDE_BY_ZERO_EXCEPTION throw std::logic_error("Floating point divide by zero!")
#define FP128_NOT_IMPLEMENTED_EXCEPTION throw std::exception("Not implemented!")
#if defined _DEBUG || defined DEBUG
#define FP128_ASSERT assert
#define FP128_THROW_ONLY_IN_DEBUG
#else
#define FP128_ASSERT(x)
#define FP128_THROW_ONLY_IN_DEBUG noexcept
#endif // _DEBUG

#if FP128_DISABLE_INLINE != FALSE
#define FP128_INLINE __declspec(noinline)
#else
#define FP128_INLINE __forceinline
#endif

namespace fp128 {

/***********************************************************************************
*                                  Constants
************************************************************************************/
// useful const calculations
static constexpr int32_t flt_frac_bits = 23;  // mantisa bit count of a float variable
static constexpr int32_t flt_exp_bits = 8;    // exponent bit count of a float variable
static constexpr int32_t dbl_frac_bits = 52;  // mantisa bit count of a double variable
static constexpr int32_t dbl_exp_bits = 11;   // exponent bit count of a double variable

/***********************************************************************************
*                                  Containers
************************************************************************************/
#pragma warning(push)
#pragma warning(disable: 4201) // nameless union/structs
struct Double {
    Double(double v = 0) : val(v) {}
    union {
        struct {
            uint64_t f : dbl_frac_bits; // mantisa/fraction
            uint64_t e : dbl_exp_bits; // exponent 
            uint64_t s : 1;  // sign
        };
        double val;
    };
};
struct Float {
    Float(float v = 0) : val(v) {}
    union {
        struct {
            uint32_t f : flt_frac_bits; // mantisa/fraction
            uint32_t e : flt_exp_bits;  // exponent 
            uint32_t s : 1;  // sign
        };
        float val;
    };
};
static_assert(sizeof(Double) == sizeof(double), "The Double union should have the same size as a double variable!");
static_assert(sizeof(Float) == sizeof(float), "The Float union should have the same size as a float variable!");

#pragma warning(pop)
/***********************************************************************************
*                                  Functions
************************************************************************************/

/**
    * @brief Calculates the element count of a C style array at build time
    * Example:
    * int a[5];
    * constexpr int len = array_length(a); // returns 5 at build time
    * @tparam T array type
    * @param a array instance
    * @return Element count in array
*/
template<typename T>
constexpr uint32_t array_length(const T& a) {
    static_assert(sizeof(a[0]) != 0, "Requires an array of non-zero sized elements!");
    return sizeof(a) / sizeof(a[0]);
}
/**
    * @brief shift right by 'shift' bits
    * Undefined behavior when shift is outside the range [0, 31]
    * @param x value to shift
    * @param shift how many bits to shift
    * @return result of 'x' the combined 64 bit element right shifted by 'shift' bits.
*/
__forceinline uint32_t shift_right64(uint32_t l, uint32_t h, int shift) 
{
    FP128_ASSERT(shift >= 0 && shift < 32);
    return (shift > 0) ? (l >> shift) | (h << (32 - shift)) : l;
}
/**
    * @brief shift left by 'shift' bits
    * Undefined behavior when shift is outside the range [0, 31]
    * @param x value to shift
    * @param shift how many bits to shift
    * @return result of the combined 64 bit element left shifted by 'shift' bits.
*/
__forceinline uint32_t shift_left64(uint32_t l, uint32_t h, int shift) 
{
    FP128_ASSERT(shift >= 0 && shift < 32);
    return (shift > 0) ?  (h << shift) | (l >> (32 - shift)) : h;
}
/**
    * @brief shift right 'x' by 'shift' bits with rounding
    * Undefined behavior when shift is outside the range [1, 63]
    * @param x value to shift
    * @param shift how many bits to shift
    * @return result of 'x' right shifted by 'shift' bits.
*/
__forceinline uint64_t shift_right64_round(uint64_t x, int shift) 
{
    FP128_ASSERT(shift >= 0 && shift < 64);
    x += 1ull << (shift - 1);
    return x >> shift;
}
/**
    * @brief Right shift a 128 bit integer.
    * @param l Low QWORD
    * @param h High QWORD
    * @param shift Bits to shift, between 0-127
    * @return Lower 64 bit of the result
*/
__forceinline uint64_t shift_right128(uint64_t l, uint64_t h, int shift)
{
    if (shift == 0) return l;
    if (shift < 64) return (l >> shift) | (h << (64 - shift));
    if (shift < 128) return h >> (shift ^ 64);
    return 0;
}
/**
    * @brief Right shift a 128 bit integer with rounding.
    * @param l Low QWORD
    * @param h High QWORD
    * @param shift Bits to shift, between 0-127
    * @return Lower 64 bit of the result
*/
__forceinline uint64_t shift_right128_round(uint64_t l, uint64_t h, int shift)
{
    if (shift == 0) return l;
    if (shift < 64) {
        const bool need_rounding = (l & 1ull << (shift - 1)) != 0;
        return need_rounding + ((l >> shift) | (h << (64 - shift)));
    }
    if (shift < 128) {
        shift ^= 64;
        const bool need_rounding = (h & 1ull << (shift - 1)) != 0;
        return need_rounding + (h >> shift);
    }
    return 0;
}
/**
    * @brief Left shift a 128 bit integer.
    * @param l Low QWORD
    * @param h High QWORD
    * @param shift Bits to shift, between 0-127
    * @return Upper 64 bit of the result
*/
__forceinline uint64_t shift_left128(uint64_t l, uint64_t h, int shift)
{
    if (shift == 0) return h;
    if (shift < 64) return (h << shift) | (l >> (64 - shift));
    if (shift < 128) return l << (shift ^ 64);
    return 0;
}
/**
    * @brief converts a 128 integer to negavite via 2's complement.
    * @param l Low QWORD (ref)
    * @param h High QWORD (ref)
    * @return void
*/
__forceinline void twos_complement128(uint64_t& l, uint64_t& h) noexcept
{
    l = ~l + 1ull;
    h = ~h + (l == 0);
}
/**
    * @brief 32 bit words unsigned divide function. Variation of the code from the book Hacker's Delight.
    * @param q (output) Pointer to receive the quote
    * @param r (output, optional) Pointer to receive the remainder. Can be nullptr
    * @param u Pointer Numerator, an array of uint32_t
    * @param v denominator (uint32_t)
    * @param m Count of elements in u
    * @return 0 for success
*/
FP128_INLINE static int32_t div_32bit(uint32_t* q, uint32_t* r, const uint32_t* u, uint32_t v, int64_t m) noexcept
{
    if (u == nullptr || q == nullptr || v == 0)
        return 1;

    while (m > 0 && u[m - 1] == 0) --m;

    uint64_t k = 0;
    for (auto j = m - 1; j >= 0; --j) {
        k = (k << 32) + u[j];
        q[j] = (uint32_t)(k / v);
        k -= (uint64_t)q[j] * v;
    }

    if (r != nullptr)
        *r = (uint32_t)k;
    return 0;
}
/**
    * @brief 32 bit words unsigned divide function. Variation of the code from the book Hacker's Delight.
    * @param q (output) Pointer to receive the quote
    * @param r (output, optional) Pointer to receive the remainder. Can be nullptr
    * @param u Pointer numerator, an array of uint32_t
    * @param v Pointer denominator, an array of uint32_t
    * @param m Count of elements in u
    * @param n Count of elements in v
    * @return 0 for success
*/
static int div_32bit(uint32_t* q, uint32_t* r, const uint32_t* u, const uint32_t* v, int m, int n) noexcept
{
    constexpr uint64_t b = 1ull << 32; // Number base (32 bits).
    constexpr uint64_t mask = b - 1;   // 32 bit mask
    uint32_t *un, *vn;                 // Normalized form of u, v.
    uint64_t qhat;                     // Estimated quotient digit.
    uint64_t rhat;                     // A remainder.
    uint64_t p;                        // Product of two digits.
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
        return div_32bit(q, r, u, v[0], m);

    /* Normalize by shifting v left just enough so that its high-order
    bit is on, and shift u left the same amount. We may have to append a
    high-order digit on the dividend; we do that unconditionally. */

    const int32_t s = __lzcnt(v[n - 1]);             // 0 <= s <= 31.
    const int32_t s_comp = 32 - s;
    vn = (uint32_t*)_alloca(sizeof(uint32_t) * n);
    for (i = n - 1; i > 0; --i) {
        //vn[i] = shift_left64(v[i - 1], v[i], s);
        vn[i] = (v[i] << s) | ((uint64_t)v[i - 1] >> s_comp);
    }
    vn[0] = v[0] << s;

    un = (uint32_t*)_alloca(sizeof(uint32_t) * (m + 1));
    un[m] = (uint64_t)u[m - 1] >> s_comp;
    for (i = m - 1; i > 0; --i)
        un[i] = (u[i] << s) | ((uint64_t)u[i - 1] >> s_comp);
    un[0] = u[0] << s;

    for (j = m - n; j >= 0; --j) {       // Main loop.
        // Compute estimate qhat of q[j].
        qhat = _udiv128(0, un[j + n] * b + un[j + n - 1], vn[n - 1], &rhat);
        //qhat = (un[j + n] * b + un[j + n - 1]) / vn[n - 1];
        //rhat = (un[j + n] * b + un[j + n - 1]) - qhat * vn[n - 1];
    again:
        if (qhat >= b || qhat * vn[n - 2] > (rhat << 32) + un[j + n - 2]) {
            --qhat;
            rhat += vn[n - 1];
            if (rhat < b) goto again;
        }

        // Multiply and subtract.
        k = 0;
        for (i = 0; i < n; ++i) {
            p = qhat * vn[i];
            t = un[i + j] - k - (p & mask);
            un[i + j] = t;
            k = (p >> 32) - (t >> 32);
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
/**
    * @brief 64 bit words unsigned divide function. Variation of the code from the book Hacker's Delight.
    * @param q (output) Pointer to receive the quote. Expected to be initialized to zero
    * @param r (output, optional) Pointer to receive the remainder. Can be nullptr
    * @param u Pointer to Numerator, an array of uint64_t
    * @param v denominator (uint64_t)
    * @param m Count of elements in u
    * @return 0 for success
*/
FP128_INLINE static int32_t div_64bit(uint64_t* q, uint64_t* r, const uint64_t* u, uint64_t v, int64_t m) noexcept
{
    if (u == nullptr || q == nullptr)
        return 1;
    uint64_t dummy_reminder{};
    if (r == nullptr) {
        r = &dummy_reminder;
    }

    if (v == 0) // error case
        return 1;

    while (m > 0 && u[m - 1] == 0) --m;

    // Trivial cases
    if (m < 2) {
        if (m == 0 || u[0] < v) {
            *r = v;
            return 0;
        }
        if (u[0] == v) {
            *q = 1;
            *r = 0;
            return 0;
        }
    }

    uint64_t k[2] = {};
    for (auto j = m - 1; j >= 0; --j) {
        k[1] = k[0];
        k[0] = u[j];
        q[j] = _udiv128(k[1], k[0], v, &k[0]);
    }

    if (r != nullptr)
        *r = k[0];
    return 0;
}
/**
 * @brief Calculates the Log base 2 of x: log2(x)
 * Rounding is always towards zero so the maximum error is close to 1.
 * @param x The number to perform log2 on.
 * @return log2(x). Returns zero when x is zero.
*/
__forceinline uint64_t log2(uint64_t x)
{
    return (x) ? 63ull - __lzcnt64(x) : 0;
}
/**
 * @brief Calculates the Log base 2 of x: log2(x)
 * Rounding is always towards zero so the maximum error is close to 1.
 * @param x The number to perform log2 on.
 * @return log2(x). Returns zero when x is zero.
*/
__forceinline uint32_t log2(uint32_t x)
{
    return (x) ? 31ull - __lzcnt(x) : 0;
}

} //namespace fp128 {
#endif // FIXED_POINT128_SHARED_H
