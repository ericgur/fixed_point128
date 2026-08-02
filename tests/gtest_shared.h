#ifndef GTEST_SHARED_H
#define GTEST_SHARED_H

#ifndef UNREFERENCED_PARAMETER
#define UNREFERENCED_PARAMETER(P) (P)
#endif

// Quoted includes for the library headers, angle brackets for system and third-party ones. The
// distinction is what lets MSVC treat GoogleTest and the standard library as external headers and
// hold only our own code to /W4; the library headers follow the same convention among themselves.
#include "fixed_point128.h"
#include "uint128_t.h"
#include "int128_t.h"
#include "float128.h"

/*************************************************
 * Fixed point 128 tests
 **************************************************/

using namespace fp128;

static constexpr int RANDOM_TEST_COUNT = 1 << 16;
static constexpr int RANDOM_SEED = 0x12345678;  // must have a repeatable seed for debugging

// These are defined below with internal linkage, and the definitions have to be preceded by a
// declaration because get_double_random() calls them. The 'static' has to appear here as well:
// a first declaration without it gives the name external linkage, which the later definition
// keeps, and every translation unit including this header then emits the same external symbol.
//
// [[maybe_unused]] is required on every helper in this header: internal linkage means each test
// translation unit gets its own copy, and no single one of them uses the whole set, which is
// exactly what MSVC reports as C4505 at /W4.
[[maybe_unused]] static uint64_t get_uint64_random();
[[maybe_unused]] static int64_t get_int64_random();
[[maybe_unused]] static uint32_t get_uint32_random();
[[maybe_unused]] static int32_t get_int32_random();

// friend class to all containers to simplify test cases
namespace fp128
{
class fp128_gtest
{
    template <int32_t I> inline static void get_fixed_point128_members(const fixed_point128<I>& obj, uint64_t& l, uint64_t& h, uint32_t& s)
    {
        l = obj.low;
        h = obj.high;
        s = obj.sign;
    }
    inline static void get_uint128_t_members(const uint128_t& obj, uint64_t& l, uint64_t& h)
    {
        l = obj.low;
        h = obj.high;
    }
};
}  // namespace fp128

// FP128_FORCE_INLINE rather than __forceinline: the latter only exists on the MSVC frontend and does
// not compile under Clang on macOS or Linux.
FP128_FORCE_INLINE int32_t get_random_sign()
{
    return (rand() & 1) ? 1 : -1;
}

// returns a random number
[[maybe_unused]] double static get_double_random(int32_t min_exponent = -10, int32_t max_exponent = 63)
{
    Double res;
    int expo = (get_uint32_random() % (max_exponent - min_exponent)) + min_exponent;
    res.e = (uint64_t)expo + 1023;
    res.f = get_uint64_random();
    res.s = get_random_sign() == 1;
    return res.val;
}

// returns a positive random number
uint64_t static get_uint64_random()
{
    return (((uint64_t)rand()) << 60) + (((uint64_t)rand()) << 45) + (((uint64_t)rand()) << 30) + (((uint64_t)rand()) << 15) + (uint64_t)rand();
}

// returns a random number
int64_t static get_int64_random()
{
    return (int64_t)get_uint64_random() * (int64_t)get_random_sign();
}

// returns a positive random number
uint32_t static get_uint32_random()
{
    return (((uint32_t)rand()) << 30) + (((uint32_t)rand()) << 15) + (uint32_t)rand();
}

// returns a random number
int32_t static get_int32_random()
{
    return (int32_t)get_uint32_random() * get_random_sign();
}

[[maybe_unused]] char static get_digit_random()
{
    return (char)(rand() % 10) + '0';
}
// return true on overflow
// 'd' is only present so that I can be deduced from the call site; its value is never read.
template <typename T, int I> bool check_overflow(T value, [[maybe_unused]] const fixed_point128<I>& d)
{
    if constexpr (std::is_floating_point<T>::value) {
        value = fabs(value);
        if (value <= 1.0)
            return false;
    } else if constexpr (std::is_signed<T>::value) {
        value = abs(value);
        if (value <= 1)
            return false;
    }

    return floor(log2(value)) >= I;
}

[[maybe_unused]] bool static check_overflow_uint128(double value)
{
    return floor(log2(abs(value))) > 127;
}

[[maybe_unused]] bool static check_overflow_int128(double value)
{
    return floor(log2(abs(value))) > 126;
}

// Exact check that r is the integer square root of x, i.e. sqrt(x) rounded down.
//
// floor(sqrt(double)) cannot serve as the reference once x goes above 2^106: the root then needs
// more than the 53 bits a double's mantissa holds, so sqrt() hands back a rounded value and the
// reference is the thing that is wrong. The defining property r^2 <= x < (r+1)^2 has no such
// ceiling and only needs a 128 bit multiply to verify.
//
// @tparam T uint128_t or int128_t
// @param r Candidate root
// @param x Value the root was taken of, must not be negative
// @return True when r is exactly floor(sqrt(x)).
template <typename T> bool static is_exact_isqrt(uint64_t r, const T& x)
{
    uint64_t l = 0, h = 0;
    x.get_components(l, h);
    const uint128_t value(l, h);

    uint64_t sqr_high = 0;
    const uint64_t sqr_low = mulx_u64(r, r, &sqr_high);
    if (uint128_t(sqr_low, sqr_high) > value)
        return false;

    // (r+1)^2 is 2^128 here, which is above every 128 bit value
    if (r == UINT64_MAX)
        return true;

    uint64_t next_high = 0;
    const uint64_t next_low = mulx_u64(r + 1, r + 1, &next_high);
    return uint128_t(next_low, next_high) > value;
}

// Reference implementation of the truncated 128 bit product, for use as a test oracle.
//
// Multiplies with the schoolbook algorithm over 32 bit limbs, using nothing but native 64 bit
// arithmetic. No partial product can overflow: the largest intermediate is
// (2^32-1)^2 + 2*(2^32-1), which is 2^64-1. The carry out of the top limb is dropped, which is
// the truncation the operators under test perform.
//
// The same computation serves both instantiations. A truncated 128 bit product is a multiplication
// modulo 2^128, and in that ring the two's complement bit pattern of a negative value is its value,
// so the signed product is the same bit pattern as the unsigned one.
//
// Deliberately shares no code with int128_base: it reaches for neither mulx_u64 nor any other
// intrinsic the implementation uses, so a defect in those cannot cancel out against it.
//
// @tparam T uint128_t or int128_t
// @param a Left hand side operand
// @param b Right hand side operand
// @return The low 128 bit of a * b.
template <typename T> T static ReferenceMultiply(const T& a, const T& b)
{
    uint64_t al = 0, ah = 0, bl = 0, bh = 0;
    a.get_components(al, ah);
    b.get_components(bl, bh);

    const uint32_t x[4] = {static_cast<uint32_t>(al), static_cast<uint32_t>(al >> 32), static_cast<uint32_t>(ah), static_cast<uint32_t>(ah >> 32)};
    const uint32_t y[4] = {static_cast<uint32_t>(bl), static_cast<uint32_t>(bl >> 32), static_cast<uint32_t>(bh), static_cast<uint32_t>(bh >> 32)};
    uint32_t r[4] {};

    for (auto i = 0u; i < 4u; ++i) {
        uint64_t carry = 0;
        for (auto j = 0u; i + j < 4u; ++j) {
            const uint64_t t = static_cast<uint64_t>(x[i]) * y[j] + r[i + j] + carry;
            r[i + j] = static_cast<uint32_t>(t);
            carry = t >> 32;
        }
    }

    return T((static_cast<uint64_t>(r[1]) << 32) | r[0], (static_cast<uint64_t>(r[3]) << 32) | r[2]);
}

// The high QWORD produced when a signed 64 bit value is sign extended to 128 bit. Mirrors what the
// integral constructor does, so a reference value can be assembled from the two QWORDs without
// going through the type under test.
[[maybe_unused]] uint64_t static SignExtension(int64_t x)
{
    return (x < 0) ? UINT64_MAX : 0ull;
}

// Absolute tolerance for comparing a fixed_point128<I> division result against a double reference.
//
// Every other fixed_point128 operation reproduces the double reference bit for bit, so only
// division needs a tolerance. Two independent error sources contribute to it, and neither one can
// stand in for the other:
//
//   - The type quantizes absolutely. fixed_point128<I> holds 128-I fraction bits and the division
//     leaves the result within a few units of that last place whatever its magnitude. This term
//     dominates for a small quotient, which keeps few fraction bits: a result near 2^-50 has only
//     38 of them left, a relative error near 1e-12 that no sane relative bound would allow.
//   - The double reference is rounded to 53 significant bits, and so is the conversion of the
//     result back to double. This term is relative and dominates for a large quotient.
//
// Both were measured at exactly 8 units of their respective error, so the 16 below is a factor of
// two of margin rather than a fitted constant. The resulting bound is around 28000 times tighter
// than the DOUBLE_REL_EPS comparison it replaced.
//
// @tparam I Integer bit count of the fixed_point128 type under test
// @param reference The double the result is being compared against
// @return Largest absolute difference to accept.
template <int32_t I> double static FixedPointDivisionTolerance(double reference)
{
    return 16.0 * (ldexp(1.0, I - 128) + fabs(reference) * ldexp(1.0, -52));
}

#endif  // #ifndef GTEST_SHARED_H
