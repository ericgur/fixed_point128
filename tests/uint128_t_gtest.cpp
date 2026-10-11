// remove warnings from gtest itself
#if defined(_MSC_VER)
#pragma warning(push)
#pragma warning(disable : 26439)
#pragma warning(disable : 26495)
#endif
#include <gtest/gtest.h>
#if defined(_MSC_VER)
#pragma warning(pop)
#endif
#include <bit>
#include <ostream>
#include <ctime>
#include <cfloat>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <string>
#include <type_traits>
#include "gtest_shared.h"

namespace
{
/**
 * @brief Builds the exact uint128_t value of a non negative, integral double below 2^128.
 *
 * The value is assembled from the double's own significand and exponent with native 64 bit shifts,
 * so the reference shares no code with the constructor or the shift operators under test.
 *
 * @param value Integral value in [0, 2^128).
 * @return The same value as a uint128_t.
 */
[[nodiscard]] uint128_t ExactFromDouble(double value)
{
    if (value == 0) {
        return uint128_t();
    }
    int exponent = 0;
    const double fraction = std::frexp(value, &exponent);  // value = fraction * 2^exponent, fraction in [0.5, 1)
    const uint64_t significand = static_cast<uint64_t>(std::ldexp(fraction, 53));
    const int32_t shift = exponent - 53;  // value = significand * 2^shift
    if (shift <= 0) {
        return uint128_t(significand >> -shift, 0ull);  // the bits shifted out are zero, the value is integral
    }
    if (shift < 64) {
        return uint128_t(significand << shift, significand >> (64 - shift));
    }
    return uint128_t(0ull, significand << (shift - 64));
}

/**
 * @brief Checks that q and r are exactly the quotient and remainder of a / b.
 *
 * q * b + r == a evaluated in 128 bit arithmetic only pins q down modulo 2^128 / gcd(b, 2^128): for an
 * even b a quotient off by 2^127 satisfies it as well, and so does garbage in the quotient's high
 * QWORD whenever the low QWORD of b is zero. The product is therefore taken in full, 256 bits wide,
 * and the sum is checked for a carry out of the top. Everything past the multiplication is native 64
 * bit arithmetic, so no operator under test takes part.
 *
 * @param a Dividend
 * @param b Divisor, non zero
 * @param q Quotient to check
 * @param r Remainder to check
 * @return Success when q * b + r == a holds over the integers and r < b.
 */
[[nodiscard]] testing::AssertionResult IsExactDivision(const uint128_t& a, const uint128_t& b, const uint128_t& q, const uint128_t& r)
{
    uint128_t high;
    const uint128_t low = ReferenceMultiplyWide(q, b, high);
    if (!high.is_zero()) {
        return testing::AssertionFailure() << "q * b overflows 128 bits, q=" << static_cast<std::string>(q);
    }

    uint64_t productLow = 0, productHigh = 0, rLow = 0, rHigh = 0, aLow = 0, aHigh = 0, bLow = 0, bHigh = 0;
    low.get_components(productLow, productHigh);
    r.get_components(rLow, rHigh);
    a.get_components(aLow, aHigh);
    b.get_components(bLow, bHigh);

    // q * b + r, with the carries made explicit
    const uint64_t sumLow = productLow + rLow;
    const uint64_t carry = (sumLow < productLow) ? 1ull : 0ull;
    const uint64_t partial = productHigh + rHigh;
    const uint64_t sumHigh = partial + carry;
    const bool carryOut = (partial < productHigh) || (sumHigh < partial);
    if (carryOut || sumLow != aLow || sumHigh != aHigh) {
        return testing::AssertionFailure() << "q * b + r != a, q=" << static_cast<std::string>(q) << ", r=" << static_cast<std::string>(r);
    }
    if (!(rHigh < bHigh || (rHigh == bHigh && rLow < bLow))) {
        return testing::AssertionFailure() << "remainder not below the divisor, r=" << static_cast<std::string>(r);
    }
    return testing::AssertionSuccess();
}

/**
 * @brief Reference shift of a 128 bit value, built one bit at a time from its QWORDs.
 *
 * Follows the documented contract of the shift operators: a count of zero or below leaves the value
 * unchanged, and a count of 128 or more shifts every bit out.
 *
 * @param x Value to shift
 * @param shift Bit count
 * @param left True for a left shift, false for a logical right shift
 * @return The shifted value.
 */
[[nodiscard]] uint128_t ReferenceShift(const uint128_t& x, int32_t shift, bool left)
{
    uint64_t in[2] {};
    x.get_components(in[0], in[1]);
    if (shift <= 0) {
        return uint128_t(in[0], in[1]);
    }

    uint64_t out[2] {};
    for (int32_t bit = 0; bit < 128; ++bit) {
        const int32_t source = left ? bit - shift : bit + shift;
        if (source < 0 || source > 127) {
            continue;
        }
        if (((in[source >> 6] >> (source & 63)) & 1ull) != 0) {
            out[bit >> 6] |= 1ull << (bit & 63);
        }
    }
    return uint128_t(out[0], out[1]);
}

/// @brief Satisfied when a * b compiles. A concept rather than a bare requires-expression, which
///        MSVC rejects outside a template.
template <class A, class B>
concept Multipliable = requires(A a, B b) { a * b; };

/// @brief Satisfied when a *= b compiles.
template <class A, class B>
concept MultiplyAssignable = requires(A& a, B b) { a *= b; };
}  // namespace

/**********************************************************************
 * uint128_t tests
 ***********************************************************************/
// Construct fixed_point128 and convert back to/from various elements.
TEST(uint128_t, DefaultConstructor)
{
    uint128_t i;
    EXPECT_EQ(static_cast<uint64_t>(i), 0ull);
}
// The round trip through operator double() used to be compared with EXPECT_DOUBLE_EQ. Its 4 ulp
// tolerance accepts a constructor that drops the two lowest significand bits, and a constructor and a
// conversion that err in opposite directions cancel out in a round trip. The object is now compared
// with a value assembled straight from the double's bits, and the round trip has to be exact.
//
// get_double_random's upper exponent is exclusive: 128 is what reaches the top binade [2^127, 2^128).
// Every value drawn is below 2^128, which is why no overflow guard is needed.
TEST(uint128_t, ConstructorFromDouble)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const double value = fabs(floor(get_double_random(1, 128)));
        const uint128_t f = value;
        EXPECT_TRUE(f == ExactFromDouble(value)) << "value=" << value << ", f=" << (std::string)f;
        EXPECT_EQ(static_cast<double>(f), value) << "value=" << value;
    }
}
// The float is built from random bits rather than rounded from a random double: a double in
// [2^128 - 2^103, 2^128) rounds to +inf as a float, and no uint128_t holds that.
TEST(uint128_t, ConstructorFromFloat)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        // 24 significant bits scaled to an exponent of 1 to 127, which spans [2, FLT_MAX]
        const float significand = static_cast<float>((1u << 23) | (get_uint32_random() & 0x7FFFFFu));
        const int32_t exponent = 1 + static_cast<int32_t>(get_uint32_random() % 127u);
        const float value = std::floor(std::ldexp(significand, exponent - 23));
        const uint128_t f = value;
        EXPECT_TRUE(f == ExactFromDouble(value)) << "value=" << value << ", f=" << (std::string)f;
        EXPECT_EQ(static_cast<float>(f), value) << "value=" << value;
    }
}
// Each integral constructor is also checked over the full 128 bits. Casting the object back to the
// source type only reads the low QWORD, so it cannot tell sign extension from zero extension.
TEST(uint128_t, ConstructorFromInt32)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        int32_t value = get_int32_random();
        uint128_t f = value;
        EXPECT_EQ(static_cast<int32_t>(f), value);
        // a negative value is sign extended through both QWORDs, wrapping like a builtin unsigned type
        EXPECT_TRUE(f == uint128_t(static_cast<uint64_t>(static_cast<int64_t>(value)), SignExtension(value))) << "value=" << value;
    }
}
TEST(uint128_t, ConstructorFromUnsignedInt32)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        uint32_t value = get_uint32_random();
        uint128_t f = value;
        EXPECT_EQ(static_cast<uint32_t>(f), value);
        EXPECT_TRUE(f == uint128_t(static_cast<uint64_t>(value), 0ull)) << "value=" << value;
    }
}
TEST(uint128_t, ConstructorFromInt64)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        int64_t value = get_int64_random();
        uint128_t f = value;
        EXPECT_EQ(static_cast<int64_t>(f), value);
        EXPECT_TRUE(f == uint128_t(static_cast<uint64_t>(value), SignExtension(value))) << "value=" << value;
    }
}
TEST(uint128_t, ConstructorFromUnsignedInt64)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        uint64_t value = get_uint64_random();
        uint128_t f = value;
        EXPECT_EQ(static_cast<uint64_t>(f), value);
        // zero extended even when the msb is set, which a signed reading would extend with ones
        EXPECT_TRUE(f == uint128_t(value, 0ull)) << "value=" << value;
    }
}
// Compared exactly. The previous version went through strtod() and EXPECT_DOUBLE_EQ, and the 57 bit
// value below differs from its neighbours in digits a double cannot hold.
TEST(uint128_t, ConstructorFromString)
{
    const struct {
        const char* text;
        uint64_t value;
    } values[] = {{"0xDEADBEAF", 0xDEADBEAFull}, {"1234", 1234ull}, {"123456789012345678", 123456789012345678ull}};
    for (const auto& v : values) {
        const uint128_t f = v.text;
        EXPECT_TRUE(f == uint128_t(v.value, 0ull)) << "text=" << v.text;
        EXPECT_STREQ(static_cast<char*>(f), std::to_string(v.value).c_str()) << "text=" << v.text;
    }
}
// The copy and move members are checked QWORD by QWORD on full width operands. Comparing the two
// objects through operator double() saw only their top 53 bits.
TEST(uint128_t, CopyConstructor)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const uint64_t low = get_uint64_random(), high = get_uint64_random();
        const uint128_t f1(low, high);
        const uint128_t f2 = f1;
        uint64_t l = 0, h = 0;
        f2.get_components(l, h);
        EXPECT_EQ(l, low);
        EXPECT_EQ(h, high);
    }
}
TEST(uint128_t, MoveConstructor)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const uint64_t low = get_uint64_random(), high = get_uint64_random();
        uint128_t f1(low, high);
        const uint128_t f2(std::move(f1));
        uint64_t l = 0, h = 0;
        f2.get_components(l, h);
        EXPECT_EQ(l, low);
        EXPECT_EQ(h, high);
        // documented to leave the source untouched, acting like a copy
        EXPECT_TRUE(f1 == f2);
    }
}
TEST(uint128_t, AssignmentOperator)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const uint64_t low = get_uint64_random(), high = get_uint64_random();
        const uint128_t f1(low, high);
        uint128_t f2;
        f2 = f1;
        uint64_t l = 0, h = 0;
        f2.get_components(l, h);
        EXPECT_EQ(l, low);
        EXPECT_EQ(h, high);
    }
}
TEST(uint128_t, MoveAssignmentOperator)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const uint64_t low = get_uint64_random(), high = get_uint64_random();
        uint128_t f1(low, high);
        uint128_t f2;
        f2 = std::move(f1);
        uint64_t l = 0, h = 0;
        f2.get_components(l, h);
        EXPECT_EQ(l, low);
        EXPECT_EQ(h, high);
        EXPECT_TRUE(f1 == f2);
    }
}
// Operands are raw QWORD pairs rather than values converted from a double. A double carries 53
// bits, so the previous version of this test could only ever see the top 53 bits of a 128 bit
// result, and it compared them through is_similar_double(), whose 1e-10 tolerance every iteration
// satisfied. The EXPECT below it therefore never executed once in 65536 iterations.
TEST(uint128_t, Add)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const uint64_t al = get_uint64_random(), ah = get_uint64_random();
        const uint64_t bl = get_uint64_random(), bh = get_uint64_random();
        const uint128_t a(al, ah), b(bl, bh);

        // the reference carry out of the low QWORD: an unsigned sum that wrapped is smaller than
        // either operand. Overflow out of the high QWORD is dropped, as the operator does.
        const uint64_t low = al + bl;
        const uint64_t high = ah + bh + ((low < al) ? 1ull : 0ull);
        const uint128_t expected(low, high);

        EXPECT_TRUE(a + b == expected) << "a=" << (std::string)a << ", b=" << (std::string)b;

        // the compound form has to agree with the binary one, including when rhs aliases *this
        uint128_t c = a;
        c += b;
        EXPECT_TRUE(c == expected) << "a=" << (std::string)a << ", b=" << (std::string)b;
        uint128_t d = a;
        d += d;
        EXPECT_TRUE(d == a + a) << "a=" << (std::string)a;

        // identities that hold in modulo 2^128 arithmetic whatever the operands wrapped through
        EXPECT_TRUE(a + b == b + a) << "a=" << (std::string)a << ", b=" << (std::string)b;
        EXPECT_TRUE((a + b) - b == a) << "a=" << (std::string)a << ", b=" << (std::string)b;
        EXPECT_TRUE(a + uint128_t() == a) << "a=" << (std::string)a;
    }
}
TEST(uint128_t, AddEdgeCases)
{
    const uint128_t zero, one(1ull), max(~0ull, ~0ull);

    // the carry has to cross the QWORD boundary
    EXPECT_TRUE(uint128_t(~0ull, 0ull) + one == uint128_t(0ull, 1ull));
    EXPECT_TRUE(uint128_t(~0ull, 0ull) + uint128_t(~0ull, 0ull) == uint128_t(~0ull - 1ull, 1ull));
    // and overflow out of the top has to wrap silently
    EXPECT_TRUE(max + one == zero);
    EXPECT_TRUE(max + max == uint128_t(~0ull - 1ull, ~0ull));
    EXPECT_TRUE(zero + zero == zero);
    // 2^64 + 2^64 == 2^65, entirely within the high QWORD
    EXPECT_TRUE(uint128_t(0ull, 1ull) + uint128_t(0ull, 1ull) == uint128_t(0ull, 2ull));
}
// The name is historical: both operands are uint64_t, the second holding the bit pattern of a negated
// value, and each is zero extended. Their sum carries into the high QWORD when it passes 2^64.
TEST(uint128_t, AddDifferentSign)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        uint64_t value1 = get_uint64_random();
        // two's complement negation written as a subtraction from zero: same bit pattern as unary
        // minus, without MSVC's C4146 complaint about negating an unsigned type
        uint64_t value2 = 0ull - get_uint64_random();
        const uint64_t low = value1 + value2;
        const uint128_t expected(low, (low < value1) ? 1ull : 0ull);
        uint128_t f1 = value1;
        uint128_t f2 = value2;
        uint128_t f3 = f1 + f2;
        EXPECT_TRUE(f3 == expected) << "value1=" << value1 << ", value2=" << value2 << ", f3=" << (std::string)f3;
    }
}
// The reference used to be computed as value1 + value2 in int64_t, which overflows on roughly a
// quarter of the operand pairs (16511 of 65536 measured) and is undefined behaviour rather than
// the wrap the test was counting on. It is now built in unsigned arithmetic, where wrapping is
// what the language actually promises.
//
// Checking the full 128 bit result rather than a cast to uint64_t is what the operation deserves:
// both operands are sign extended to 128 bit before being added, so their sum is exact and the
// high QWORD carries real information. The old cast discarded it.
TEST(uint128_t, AddInt64)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const int64_t value1 = get_int64_random();
        const int64_t value2 = get_int64_random();

        const uint64_t u1 = static_cast<uint64_t>(value1), u2 = static_cast<uint64_t>(value2);
        const uint64_t low = u1 + u2;
        const uint64_t high = SignExtension(value1) + SignExtension(value2) + ((low < u1) ? 1ull : 0ull);
        const uint128_t expected(low, high);

        const uint128_t f1 = value1;
        EXPECT_TRUE(f1 + value2 == expected) << "value1=" << value1 << ", value2=" << value2;
        // adding the scalar has to match promoting it to uint128_t first
        EXPECT_TRUE(f1 + value2 == f1 + uint128_t(value2)) << "value1=" << value1 << ", value2=" << value2;

        uint128_t f2 = f1;
        f2 += value2;
        EXPECT_TRUE(f2 == expected) << "value1=" << value1 << ", value2=" << value2;
    }
}
// About half of these sums pass 2^64, so the carry into the high QWORD is checked as well as the low
// QWORD. The scalar is zero extended even when its msb is set.
TEST(uint128_t, AddUnsignedInt64)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        auto value1 = get_uint64_random();
        auto value2 = get_uint64_random();
        const uint64_t low = value1 + value2;
        const uint128_t expected(low, (low < value1) ? 1ull : 0ull);
        uint128_t f1 = value1;
        uint128_t f3 = f1 + value2;
        EXPECT_TRUE(f3 == expected) << "value1=" << value1 << ", value2=" << value2 << ", f3=" << (std::string)f3;
        f1 += value2;
        EXPECT_TRUE(f1 == expected) << "value1=" << value1 << ", value2=" << value2;
    }
}

// The previous version derived the subtrahend as floor(value1 / 2.5), so the result was never
// negative and the borrow out of the top QWORD went untested. Independent random operands reach
// it roughly half the time.
TEST(uint128_t, Subtract)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const uint64_t al = get_uint64_random(), ah = get_uint64_random();
        const uint64_t bl = get_uint64_random(), bh = get_uint64_random();
        const uint128_t a(al, ah), b(bl, bh);

        // the reference borrow out of the low QWORD. Underflow out of the high QWORD wraps.
        const uint64_t low = al - bl;
        const uint64_t high = ah - bh - ((al < bl) ? 1ull : 0ull);
        const uint128_t expected(low, high);

        EXPECT_TRUE(a - b == expected) << "a=" << (std::string)a << ", b=" << (std::string)b;

        uint128_t c = a;
        c -= b;
        EXPECT_TRUE(c == expected) << "a=" << (std::string)a << ", b=" << (std::string)b;
        // a subtrahend that aliases the target always yields zero
        uint128_t d = a;
        d -= d;
        EXPECT_TRUE(d.is_zero()) << "a=" << (std::string)a;

        EXPECT_TRUE((a - b) + b == a) << "a=" << (std::string)a << ", b=" << (std::string)b;
        EXPECT_TRUE(a - b == -(b - a)) << "a=" << (std::string)a << ", b=" << (std::string)b;
        EXPECT_TRUE(a - uint128_t() == a) << "a=" << (std::string)a;
        EXPECT_TRUE(uint128_t() - a == -a) << "a=" << (std::string)a;
    }
}
TEST(uint128_t, SubtractEdgeCases)
{
    const uint128_t zero, one(1ull), max(~0ull, ~0ull);

    // the borrow has to cross the QWORD boundary
    EXPECT_TRUE(uint128_t(0ull, 1ull) - one == uint128_t(~0ull, 0ull));
    // and underflow past zero has to wrap silently
    EXPECT_TRUE(zero - one == max);
    EXPECT_TRUE(zero - max == one);
    EXPECT_TRUE(zero - zero == zero);
    EXPECT_TRUE(max - max == zero);
    EXPECT_TRUE(one - uint128_t(2ull) == max);
}
// The name is historical: both operands are uint64_t and zero extended, see AddDifferentSign. A
// difference below zero borrows through the high QWORD, which becomes all ones.
TEST(uint128_t, SubtractDifferentSign)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        uint64_t value1 = get_uint64_random();
        // two's complement negation written as a subtraction from zero: same bit pattern as unary
        // minus, without MSVC's C4146 complaint about negating an unsigned type
        uint64_t value2 = 0ull - get_uint64_random();
        const uint128_t expected(value1 - value2, (value1 < value2) ? ~0ull : 0ull);
        uint128_t f1 = value1;
        uint128_t f2 = value2;
        uint128_t f3 = f1 - f2;
        EXPECT_TRUE(f3 == expected) << "value1=" << value1 << ", value2=" << value2 << ", f3=" << (std::string)f3;
    }
}
// value1 - value2 in int64_t overflows on roughly a quarter of the operand pairs (16429 of 65536
// measured), which is undefined behaviour. See AddInt64 for the rest of the reasoning.
TEST(uint128_t, SubtractInt64)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const int64_t value1 = get_int64_random();
        const int64_t value2 = get_int64_random();

        const uint64_t u1 = static_cast<uint64_t>(value1), u2 = static_cast<uint64_t>(value2);
        const uint64_t low = u1 - u2;
        const uint64_t high = SignExtension(value1) - SignExtension(value2) - ((u1 < u2) ? 1ull : 0ull);
        const uint128_t expected(low, high);

        const uint128_t f1 = value1;
        EXPECT_TRUE(f1 - value2 == expected) << "value1=" << value1 << ", value2=" << value2;
        EXPECT_TRUE(f1 - value2 == f1 - uint128_t(value2)) << "value1=" << value1 << ", value2=" << value2;

        uint128_t f2 = f1;
        f2 -= value2;
        EXPECT_TRUE(f2 == expected) << "value1=" << value1 << ", value2=" << value2;
    }
}
// About half of these differences are negative, which borrows through the high QWORD.
TEST(uint128_t, SubtractUnsignedInt64)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        auto value1 = get_uint64_random();
        auto value2 = get_uint64_random();
        const uint128_t expected(value1 - value2, (value1 < value2) ? ~0ull : 0ull);
        uint128_t f1 = value1;
        uint128_t f3 = f1 - value2;
        EXPECT_TRUE(f3 == expected) << "value1=" << value1 << ", value2=" << value2 << ", f3=" << (std::string)f3;
        f1 -= value2;
        EXPECT_TRUE(f1 == expected) << "value1=" << value1 << ", value2=" << value2;
    }
}
// Checked against ReferenceMultiply() rather than a double product. Besides losing everything
// below the top 53 bits, the double oracle forced the old version to skip every operand pair whose
// product exceeded 2^128, which is where truncation happens and therefore the half of the input
// space most worth testing. It skipped 33508 of 65536 iterations on that count alone, and the
// 1e-10 tolerance swallowed the remaining 32028.
TEST(uint128_t, MultiplyByUint128)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        // shifting the high QWORDs by independent random amounts mixes operands that overflow
        // with operands whose product stays within 128 bit
        const uint128_t a(get_uint64_random(), get_uint64_random() >> (get_uint32_random() % 64));
        const uint128_t b(get_uint64_random(), get_uint64_random() >> (get_uint32_random() % 64));
        const uint128_t expected = ReferenceMultiply(a, b);

        EXPECT_TRUE(a * b == expected) << "a=" << (std::string)a << ", b=" << (std::string)b;

        uint128_t c = a;
        c *= b;
        EXPECT_TRUE(c == expected) << "a=" << (std::string)a << ", b=" << (std::string)b;

        // identities that hold in modulo 2^128 arithmetic, so they survive truncation
        EXPECT_TRUE(a * b == b * a) << "a=" << (std::string)a << ", b=" << (std::string)b;
        EXPECT_TRUE(a * uint128_t(1ull) == a) << "a=" << (std::string)a;
        EXPECT_TRUE((a * uint128_t()).is_zero()) << "a=" << (std::string)a;
        EXPECT_TRUE(a * uint128_t(2ull) == a + a) << "a=" << (std::string)a;
        EXPECT_TRUE(a * (b + uint128_t(1ull)) == a * b + a) << "a=" << (std::string)a << ", b=" << (std::string)b;
    }
}
TEST(uint128_t, MultiplyByUint128EdgeCases)
{
    const uint128_t zero, one(1ull), max(~0ull, ~0ull);
    const uint128_t two_64(0ull, 1ull);  // 2^64

    EXPECT_TRUE(zero * zero == zero);
    EXPECT_TRUE(max * one == max);
    // 2^64 * 2^64 is 2^128, every bit of which truncates away
    EXPECT_TRUE(two_64 * two_64 == zero);
    // (2^128-1)^2 is 2^256 - 2^129 + 1, whose low 128 bit are 1
    EXPECT_TRUE(max * max == one);
    // 2^64 * (2^64 + 1) keeps only the 2^64 term
    EXPECT_TRUE(two_64 * uint128_t(1ull, 1ull) == two_64);
    // the largest product that does not overflow at all
    EXPECT_TRUE(uint128_t(~0ull, 0ull) * uint128_t(1ull) == uint128_t(~0ull, 0ull));

    // the reference has to agree on all of the above
    const uint128_t values[] = {zero, one, two_64, max, uint128_t(~0ull, 0ull), uint128_t(1ull, 1ull), uint128_t(0ull, ~0ull)};
    for (const auto& x : values) {
        for (const auto& y : values) {
            EXPECT_TRUE(x * y == ReferenceMultiply(x, y)) << "x=" << (std::string)x << ", y=" << (std::string)y;
        }
    }
}
// sqr(x) must be bit identical to x * x, it computes the same truncated product with one
// less multiply. Raw bit patterns are used so the whole 128 bit range is covered, including
// the values that overflow, which double cannot represent.
TEST(uint128_t, SqrMatchesMultiply)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const uint64_t low = get_uint64_random();
        const uint64_t high = get_uint64_random();
        const uint128_t x(low, high);
        EXPECT_TRUE(sqr(x) == x * x) << "low=" << low << ", high=" << high;
    }
}
TEST(uint128_t, SqrEdgeCases)
{
    const uint128_t values[] = {uint128_t(0ull, 0ull),   uint128_t(1ull, 0ull),   uint128_t(0ull, 1ull),
                                uint128_t(~0ull, 0ull),  uint128_t(0ull, ~0ull),  uint128_t(~0ull, ~0ull),
                                uint128_t(2ull, 0ull),   uint128_t(0ull, 2ull),   uint128_t(1ull, 1ull)};
    for (const auto& x : values) {
        EXPECT_TRUE(sqr(x) == x * x) << "x=" << (std::string)x;
    }
    // small values whose square is exact
    for (uint64_t v = 0; v < 1000; ++v) {
        const uint128_t x(v, 0);
        EXPECT_TRUE(sqr(x) == uint128_t(v * v, 0)) << "v=" << v;
    }
}
// floor(value1 / value2) computed in double is not a usable reference: the quotient of two 128 bit
// values needs more than the 53 bits a double holds, so the old oracle was wrong rather than the
// implementation, and is_similar_double() hid the difference on all 65536 iterations.
//
// The defining property of integer division needs no reference implementation and no rounding
// argument: q*b + r == a with r < b identifies q and r uniquely, but only over the integers.
// Evaluated in 128 bit arithmetic it wraps, and then a quotient off by 2^127 passes for every even
// divisor. IsExactDivision() therefore takes the product in full and rejects a carry out of the top.
TEST(uint128_t, DivideByUint128)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        // One iteration in four keeps both operands inside 64 bit. That is what makes the native
        // reference at the bottom of the loop reachable: with a full width random pair the high
        // QWORD is set essentially always, and the check ran 0 times out of 65536.
        const bool narrow = (get_uint32_random() & 3) == 0;
        const uint128_t a = narrow ? uint128_t(get_uint64_random(), 0ull)
                                   : uint128_t(get_uint64_random(), get_uint64_random());

        // The divisor's magnitude is spread over the whole range so that all three paths in
        // operator/= are reached: the power of two shift, the 64 bit divisor shortcut when the
        // high QWORD shifts out entirely, and the full 128 bit long division.
        // the shift counts are int32_t: an unsigned one leaves the shift operators ambiguous
        // against the builtin shifts reachable through operator uint64_t()
        const int32_t width = narrow ? 64 : 128;
        const int32_t shift = static_cast<int32_t>(get_uint32_random() % static_cast<uint32_t>(width));
        uint128_t b = (narrow ? uint128_t(get_uint64_random(), 0ull)
                              : uint128_t(get_uint64_random(), get_uint64_random())) >> shift;
        if ((get_uint32_random() & 7) == 0)
            b = uint128_t(1ull) << static_cast<int32_t>(get_uint32_random() % static_cast<uint32_t>(width));  // exercise the shift path
        if (b.is_zero())
            continue;

        const uint128_t q = a / b;
        const uint128_t r = a % b;

        EXPECT_TRUE(IsExactDivision(a, b, q, r)) << "a=" << (std::string)a << ", b=" << (std::string)b;

        // the compound forms have to agree with the binary ones
        uint128_t c = a;
        c /= b;
        EXPECT_TRUE(c == q) << "a=" << (std::string)a << ", b=" << (std::string)b;
        uint128_t d = a;
        d %= b;
        EXPECT_TRUE(d == r) << "a=" << (std::string)a << ", b=" << (std::string)b;

        // when both operands fit in 64 bit the native operators are an exact reference
        uint64_t al = 0, ah = 0, bl = 0, bh = 0;
        a.get_components(al, ah);
        b.get_components(bl, bh);
        if (ah == 0 && bh == 0) {
            EXPECT_TRUE(q == uint128_t(al / bl, 0ull)) << "a=" << al << ", b=" << bl;
            EXPECT_TRUE(r == uint128_t(al % bl, 0ull)) << "a=" << al << ", b=" << bl;
        }
    }
}
TEST(uint128_t, DivideByUint128EdgeCases)
{
    const uint128_t zero, one(1ull), max(~0ull, ~0ull);

    EXPECT_TRUE(max / max == one);
    EXPECT_TRUE(max / one == max);
    EXPECT_TRUE((max % max).is_zero());
    EXPECT_TRUE((max % one).is_zero());
    EXPECT_TRUE((zero / max).is_zero());
    EXPECT_TRUE((zero % max).is_zero());
    // a divisor larger than the dividend truncates to zero and leaves the dividend as remainder
    EXPECT_TRUE((one / max).is_zero());
    EXPECT_TRUE(one % max == one);
    // a quotient that spans both QWORDs
    EXPECT_TRUE(max / uint128_t(3ull) == uint128_t(0x5555555555555555ull, 0x5555555555555555ull));
    EXPECT_TRUE(max % uint128_t(3ull) == zero);
    // 2^128-1 is 255 modulo 2^8, the divisor being a power of two takes the shift path
    EXPECT_TRUE(max % uint128_t(256ull) == uint128_t(255ull));

    // division by zero is reported, not silently wrong
    EXPECT_THROW((void)(max / zero), std::logic_error);
    EXPECT_THROW((void)(max % zero), std::logic_error);
    EXPECT_THROW((void)(max / 0ull), std::logic_error);
}
TEST(uint128_t, ModuloByUint128)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        auto value1 = get_uint64_random();
        auto value2 = get_uint32_random();
        if (value2 == 0)
            continue;
        uint64_t res = value1 % value2;
        uint128_t f1 = value1;
        uint128_t f2 = value2;
        uint128_t f3 = f1 % f2;
        // the whole 128 bit remainder, not just its low QWORD
        EXPECT_TRUE(f3 == uint128_t(res, 0ull)) << "value1=" << value1 << ", value2=" << value2 << ", f3=" << (std::string)f3;
    }
}
// Integral doubles convert exactly, so comparing the doubles is an exact reference. The upper exponent
// of get_double_random is exclusive: 128 reaches the top binade and keeps every operand below 2^128.
// CompareRawQwords covers operands that differ only in their low bits, which doubles cannot express.
TEST(uint128_t, CompareUint128)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = fabs(floor(get_double_random(1, 128)));
        double value2 = fabs(floor(get_double_random(1, 128)));
        bool res = value1 > value2;
        uint128_t f1 = value1;
        uint128_t f2 = value2;

        bool uint128_res = f1 > f2;
        EXPECT_TRUE(uint128_res == res) << "operator>: "
                                        << "value1=" << value1 << ", value2=" << value2;

        res = value1 >= value2;
        uint128_res = f1 >= f2;
        EXPECT_TRUE(uint128_res == res) << "operator>=: "
                                        << "value1=" << value1 << ", value2=" << value2;

        res = value1 < value2;
        uint128_res = f1 < f2;
        EXPECT_TRUE(uint128_res == res) << "operator<: "
                                        << "value1=" << value1 << ", value2=" << value2;

        res = value1 <= value2;
        uint128_res = f1 <= f2;
        EXPECT_TRUE(uint128_res == res) << "operator<=: "
                                        << "value1=" << value1 << ", value2=" << value2;
    }
}
// Full width operands whose low QWORD is all ones one iteration in eight, so the carry into the high
// QWORD is reached, and the whole result is checked. The postfix forms have to hand back the value
// from before the step, and the prefix forms the object itself.
TEST(uint128_t, OperatorPlusPlus)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const uint64_t low = ((i & 7u) == 0u) ? ~0ull : get_uint64_random();
        const uint64_t high = get_uint64_random();
        const uint128_t f1(low, high);
        // a low QWORD that wraps to zero carries into the high QWORD, which itself wraps at the top
        const uint128_t expected(low + 1ull, (low == ~0ull) ? high + 1ull : high);

        uint128_t f2 = f1;
        const uint128_t old = f2++;
        EXPECT_TRUE(f2 == expected) << "operator++(int) a=" << (std::string)f1;
        EXPECT_TRUE(old == f1) << "operator++(int) returned " << (std::string)old << ", a=" << (std::string)f1;

        f2 = f1;
        EXPECT_EQ(&(++f2), &f2) << "operator++() has to return the object itself";
        EXPECT_TRUE(f2 == expected) << "operator++() a=" << (std::string)f1;
    }
}
TEST(uint128_t, OperatorMinusMinus)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const uint64_t low = ((i & 7u) == 0u) ? 0ull : get_uint64_random();
        const uint64_t high = get_uint64_random();
        const uint128_t f1(low, high);
        // a low QWORD of zero borrows from the high QWORD, which itself wraps at the bottom
        const uint128_t expected(low - 1ull, (low == 0ull) ? high - 1ull : high);

        uint128_t f2 = f1;
        const uint128_t old = f2--;
        EXPECT_TRUE(f2 == expected) << "operator--(int) a=" << (std::string)f1;
        EXPECT_TRUE(old == f1) << "operator--(int) returned " << (std::string)old << ", a=" << (std::string)f1;

        f2 = f1;
        EXPECT_EQ(&(--f2), &f2) << "operator--() has to return the object itself";
        EXPECT_TRUE(f2 == expected) << "operator--() a=" << (std::string)f1;
    }
}
// The previous version asserted f1 == f1, which is true of any implementation that manages to read
// its own members twice, and never compared two distinct objects. In particular it could not tell
// a comparison of the whole 128 bit value from one that only looks at a single QWORD. Operands
// that differ in exactly one QWORD are what pins that down.
TEST(uint128_t, OperatorEqual)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const uint64_t al = get_uint64_random(), ah = get_uint64_random();
        const uint64_t bl = get_uint64_random(), bh = get_uint64_random();
        const uint128_t a(al, ah), b(bl, bh);

        // equality is the conjunction over both QWORDs, so independent random operands are
        // almost always unequal. The cases below supply the equal ones.
        EXPECT_EQ(a == b, (al == bl) && (ah == bh)) << "a=" << (std::string)a << ", b=" << (std::string)b;

        // a copy compares equal
        const uint128_t copy = a;
        EXPECT_TRUE(a == copy) << "a=" << (std::string)a;
        // and a value that differs in exactly one QWORD does not, in either position
        EXPECT_FALSE(a == uint128_t(al ^ 1ull, ah)) << "a=" << (std::string)a;
        EXPECT_FALSE(a == uint128_t(al, ah ^ 1ull)) << "a=" << (std::string)a;
        // including when the difference is confined to the msb of either QWORD
        EXPECT_FALSE(a == uint128_t(al ^ (1ull << 63), ah)) << "a=" << (std::string)a;
        EXPECT_FALSE(a == uint128_t(al, ah ^ (1ull << 63))) << "a=" << (std::string)a;

        // the templated overload promotes the scalar, so it has to weigh the high QWORD too.
        // Comparing only the low QWORD would report equality for every value of ah.
        EXPECT_EQ(a == al, ah == 0) << "a=" << (std::string)a;
    }
}
// operator!= is declared separately rather than being the compiler's rewrite of operator==, so it
// can disagree with it. Every case below checks both that it is correct and that the two agree.
TEST(uint128_t, OperatorNotEqual)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const uint64_t al = get_uint64_random(), ah = get_uint64_random();
        const uint64_t bl = get_uint64_random(), bh = get_uint64_random();
        const uint128_t a(al, ah), b(bl, bh);

        EXPECT_EQ(a != b, (al != bl) || (ah != bh)) << "a=" << (std::string)a << ", b=" << (std::string)b;
        EXPECT_EQ(a != b, !(a == b)) << "a=" << (std::string)a << ", b=" << (std::string)b;

        const uint128_t copy = a;
        EXPECT_FALSE(a != copy) << "a=" << (std::string)a;
        EXPECT_TRUE(a != uint128_t(al ^ 1ull, ah)) << "a=" << (std::string)a;
        EXPECT_TRUE(a != uint128_t(al, ah ^ 1ull)) << "a=" << (std::string)a;

        EXPECT_EQ(a != al, ah != 0) << "a=" << (std::string)a;
        EXPECT_EQ(a != al, !(a == al)) << "a=" << (std::string)a;
    }
}
TEST(uint128_t, OperatorEqualEdgeCases)
{
    const uint128_t zero, one(1ull), max(~0ull, ~0ull);
    const uint128_t two_64(0ull, 1ull);  // 2^64

    EXPECT_TRUE(zero == uint128_t());
    EXPECT_TRUE(max == uint128_t(~0ull, ~0ull));
    EXPECT_FALSE(zero == one);
    EXPECT_TRUE(zero != one);

    // the two values either side of the QWORD boundary share no bits yet both are "1" in one QWORD
    EXPECT_FALSE(one == two_64);
    EXPECT_TRUE(one != two_64);

    // a value above 2^64 must never compare equal to the scalar sitting in its low QWORD
    EXPECT_FALSE(two_64 == 0ull);
    EXPECT_TRUE(two_64 != 0ull);
    EXPECT_TRUE(uint128_t(5ull, 0ull) == 5ull);
    EXPECT_FALSE(uint128_t(5ull, 1ull) == 5ull);

    // default constructed objects are all equal, and equality is reflexive for the extremes
    EXPECT_TRUE(uint128_t() == uint128_t());
    EXPECT_FALSE(uint128_t() != uint128_t());
    EXPECT_TRUE(max == max);
    EXPECT_FALSE(max != max);
}
// The double reference is only trustworthy away from the powers of ten, which random operands
// essentially never land next to; Log10Boundaries checks the boundaries themselves exactly. The upper
// exponent of get_double_random is exclusive, so 128 reaches [2^127, 2^128) and nothing beyond it.
TEST(uint128_t, log10)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = floor(1.0 + fabs(get_double_random(1, 128)));
        uint64_t res = (uint64_t)floor(log10(value1));
        uint128_t i1 = value1;
        uint64_t i_res = log10(i1);
        EXPECT_EQ(i_res, res) << "double value1=" << value1;
    }
}
// Like log10, the double reference holds away from the table boundaries, which LogBoundaries covers.
TEST(uint128_t, log)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = 1.0 + floor(fabs(get_double_random(0, 128)));
        uint64_t res = (uint64_t)floor(log(value1));
        uint128_t i1 = value1;
        uint64_t i_res = log(i1);
        EXPECT_EQ(i_res, res) << "double value1=" << value1;
    }
}
// Operands that a double can hold; SqrtOfFullWidthValues covers the rest of the 128 bit range.
TEST(uint128_t, sqrt)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = floor(1.0 + fabs(get_double_random(1, 128)));
        uint128_t i1 = value1;
        uint64_t i_res = sqrt(i1);
        // floor(sqrt(value1)) is not a usable reference at this range, the root needs more than the
        // 53 bits a double holds. is_exact_isqrt() checks r^2 <= x < (r+1)^2 in 128 bit instead.
        EXPECT_TRUE(is_exact_isqrt(i_res, i1)) << "double value1=" << value1 << ", sqrt=" << i_res;
        // the conversion of the operand has to be exact, otherwise the check above passes vacuously
        EXPECT_EQ(static_cast<double>(i1), value1);
    }
}
// Small operands whose power fits in 64 bits (7^15 < 2^43), so repeated native multiplication is an
// exact reference. powEdgeCases covers results past 64 bits and the wrap around at 2^128.
TEST(uint128_t, pow)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const uint64_t base = get_uint32_random() % 8u;
        const uint32_t exponent = get_uint32_random() % 16u;
        uint64_t expected = 1;
        for (auto k = 0u; k < exponent; ++k) {
            expected *= base;
        }
        const uint128_t result = pow(uint128_t(base), exponent);
        EXPECT_TRUE(result == uint128_t(expected, 0ull)) << "base=" << base << ", exponent=" << exponent << ", result=" << (std::string)result;
    }
}

/**********************************************************************
 * uint128_t regression tests
 *
 * Each test below pins down a defect that was fixed. They deliberately use
 * operand ranges the original tests never reached.
 ***********************************************************************/

// udiv128() is the 128/64 division primitive underneath every integer divide in the library, and on
// x86-64 with a GCC style frontend it is now hand written assembly: Clang cannot prove the quotient
// fits in 64 bits, so it lowered the portable __uint128_t expression to a __udivti3 call - a full
// software 128/128 division that measured 2.7x slower than the DIV instruction. MSVC reaches the
// same instruction through the _udiv128 intrinsic, and the macro makes this test cover both.
//
// The division identity is checked rather than a second implementation, so the test says the same
// thing on every toolchain. Every case keeps hi < divisor, which is the precondition both spellings
// carry: DIV raises #DE instead of truncating when the quotient overflows 64 bits.
TEST(uint128_t, UDiv128SatisfiesTheDivisionIdentity)
{
    const auto check = [](uint64_t hi, uint64_t lo, uint64_t divisor) {
        ASSERT_NE(divisor, 0ull);
        ASSERT_LT(hi, divisor) << "the quotient has to fit in 64 bits";

        uint64_t remainder = ~0ull;
        const uint64_t quotient = udiv128(hi, lo, divisor, &remainder);

        EXPECT_LT(remainder, divisor) << "hi=" << hi << ", lo=" << lo << ", divisor=" << divisor;
        // quotient * divisor + remainder reconstructs the dividend exactly, over the integers
        EXPECT_TRUE(IsExactDivision(uint128_t(lo, hi), uint128_t(divisor, 0ull), uint128_t(quotient, 0ull), uint128_t(remainder, 0ull)))
            << "hi=" << hi << ", lo=" << lo << ", divisor=" << divisor << ", q=" << quotient << ", r=" << remainder;
        // the remainder pointer is optional, and dropping it must not change the quotient
        EXPECT_EQ(udiv128(hi, lo, divisor, nullptr), quotient);
    };

    // The interesting divisors are the extremes and the small values the library actually divides
    // by; the interesting high halves sit just below the divisor, where the quotient barely fits.
    const uint64_t divisors[] = {1ull, 2ull, 3ull, 5ull, 10ull, 0xFFFFFFFFull, 0x8000000000000000ull, ~0ull};
    const uint64_t lows[] = {0ull, 1ull, 0x123456789ABCDEFull, 0x8000000000000000ull, ~0ull};
    for (const uint64_t divisor : divisors) {
        for (const uint64_t lo : lows) {
            check(0ull, lo, divisor);
            check(divisor - 1, lo, divisor);
            check(divisor / 2, lo, divisor);
        }
    }

    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const uint64_t divisor = get_uint64_random() | 1ull;  // never zero
        check(get_uint64_random() % divisor, get_uint64_random(), divisor);
    }
}
// div_128bit() is the long division every 128 bit type reaches once the divisor needs more than 64
// bits, so it is checked directly rather than only through the operators that call it. The identity
// q * v + r == u with r < v is toolchain neutral and needs no reference implementation to compare
// against, which matters here: the routine answers an error rather than a quotient for the shapes
// it does not support, and a reference that made the same choices would hide that.
TEST(uint128_t, Div128BitSatisfiesTheDivisionIdentity)
{
    const auto check = [](const uint64_t* u, int64_t m, uint64_t vLow, uint64_t vHigh) {
        ASSERT_NE(vHigh, 0ull) << "div_128bit rejects a divisor that fits in one QWORD";

        const uint64_t v[2] = {vLow, vHigh};
        // The quotient has m - 1 words and only those may be written. The sentinel in the rest shows
        // they were left alone, which the 128 bit callers rely on to keep q[1] at zero.
        constexpr uint64_t sentinel = 0xA5A5A5A5A5A5A5A5ull;
        uint64_t q[4] = {sentinel, sentinel, sentinel, sentinel}, r[2] {};
        ASSERT_TRUE(div_128bit(q, r, u, v, m));
        for (int64_t w = m - 1; w < 4; ++w) {
            EXPECT_EQ(q[w], sentinel) << "quotient word " << w << " was written, m=" << m;
        }

        // Back-multiply the quotient by the divisor and add the remainder; the sum has to be the
        // dividend exactly, in every word.
        uint64_t product[5] {};
        for (int64_t i = 0; i < m - 1; ++i) {
            for (int64_t k = 0; k < 2; ++k) {
                uint64_t hi = 0;
                const uint64_t lo = mulx_u64(q[i], v[k], &hi);
                unsigned char carry = addcarryx_u64(0, product[i + k], lo, &product[i + k]);
                carry = addcarryx_u64(carry, product[i + k + 1], hi, &product[i + k + 1]);
                for (int64_t w = i + k + 2; carry != 0 && w < 5; ++w) {
                    carry = addcarryx_u64(0, product[w], carry, &product[w]);
                }
            }
        }
        unsigned char carry = addcarryx_u64(0, product[0], r[0], &product[0]);
        carry = addcarryx_u64(carry, product[1], r[1], &product[1]);
        for (int64_t w = 2; carry != 0 && w < 5; ++w) {
            carry = addcarryx_u64(0, product[w], carry, &product[w]);
        }

        for (int64_t i = 0; i < m; ++i) {
            EXPECT_EQ(product[i], u[i]) << "word " << i << " of q * v + r, m=" << m;
        }
        // Every word above the dividend has to be zero, not only the top one: a quotient off by a
        // multiple of 2^128 / v would otherwise still match the dividend in its low words.
        for (int64_t w = m; w < 5; ++w) {
            EXPECT_EQ(product[w], 0ull) << "q * v + r overflowed the dividend into word " << w << ", m=" << m;
        }
        EXPECT_TRUE(uint128_t(r[0], r[1]) < uint128_t(vLow, vHigh)) << "remainder not reduced, m=" << m;
    };

    // A divisor whose high QWORD equals the dividend's top word is what drives the quotient estimate
    // to its 2^64 - 1 cap, the one branch random operands almost never reach.
    const uint64_t edges[] = {0ull, 1ull, 0x123456789ABCDEFull, 0x8000000000000000ull, ~0ull};
    for (const uint64_t a : edges) {
        for (const uint64_t b : edges) {
            for (const uint64_t c : edges) {
                const uint64_t u4[4] = {a, b, c, a};
                const uint64_t u3[3] = {a, b, c};
                const uint64_t u2[2] = {a, b};
                check(u4, 4, c, (b == 0) ? 1ull : b);
                check(u4, 4, a, (a == 0) ? 0x8000000000000000ull : a);
                check(u3, 3, b, (c == 0) ? 1ull : c);
                check(u3, 3, a, (c == 0) ? 0x8000000000000000ull : c);
                check(u2, 2, c, (a == 0) ? 1ull : a);
            }
        }
    }

    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        // Sweep the divisor's magnitude so every normalization shift count is exercised.
        const uint64_t vHigh = (get_uint64_random() >> (i % 64)) | 1ull;
        const uint64_t u4[4] = {get_uint64_random(), get_uint64_random(), get_uint64_random(), get_uint64_random()};
        const uint64_t u3[3] = {get_uint64_random(), get_uint64_random(), get_uint64_random()};
        const uint64_t u2[2] = {get_uint64_random(), get_uint64_random()};
        check(u4, 4, get_uint64_random(), vHigh);
        check(u3, 3, get_uint64_random(), vHigh);
        check(u2, 2, get_uint64_random(), vHigh);
        // The shape the fractional types build: the value scaled by 2^128.
        const uint64_t scaled[4] = {0, 0, get_uint64_random(), get_uint64_random()};
        check(scaled, 4, get_uint64_random(), vHigh);
    }
}
// div_128bit rejects what it cannot divide rather than answering wrongly: a divisor that fits in one
// QWORD belongs in div_64bit, and the temporaries are sized for a 256 bit dividend.
TEST(uint128_t, Div128BitRejectsUnsupportedShapes)
{
    const uint64_t u[4] = {1, 2, 3, 4};
    const uint64_t v[2] = {5, 6};
    uint64_t q[4] {}, r[2] {};

    EXPECT_TRUE(div_128bit(q, r, u, v, 4));
    EXPECT_FALSE(div_128bit(nullptr, r, u, v, 4)) << "null quotient";
    EXPECT_FALSE(div_128bit(q, r, nullptr, v, 4)) << "null dividend";
    EXPECT_FALSE(div_128bit(q, r, u, nullptr, 4)) << "null divisor";
    EXPECT_FALSE(div_128bit(q, r, u, v, 1)) << "dividend narrower than the divisor";
    EXPECT_FALSE(div_128bit(q, r, u, v, 5)) << "dividend wider than 256 bits";

    const uint64_t narrow[2] = {7, 0};
    EXPECT_FALSE(div_128bit(q, r, u, narrow, 4)) << "divisor that fits in one QWORD";
}
// operator%= used to hand the long division this object's own QWORDs as the remainder buffer. The
// routine of the day shrank the denominator past its leading zero words and filled only that many,
// so a divisor whose high QWORD fit in 32 bits left the top 32 bits of the numerator untouched in
// the result. ModuloByUint128 misses this: it only uses a 64 bit numerator and a 32 bit divisor.
TEST(uint128_t, ModuloByWideUint128)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        // shifting the high QWORD by a random amount covers both the 3 and 4 word denominators
        const uint128_t a(get_uint64_random(), get_uint64_random());
        const uint128_t b(get_uint64_random(), get_uint64_random() >> (get_uint32_random() % 64));
        if (b.is_zero())
            continue;

        const uint128_t q = a / b;
        const uint128_t r = a % b;
        EXPECT_TRUE(IsExactDivision(a, b, q, r)) << "a=" << (std::string)a << ", b=" << (std::string)b;
    }
}
TEST(uint128_t, ModuloByWideUint128EdgeCases)
{
    // 2^64 = -1 modulo (2^64 + 1), so 2^127 = -2^63 and the remainder is 2^63 + 1.
    // The denominator shrinks to 3 words here, which is what used to leave a stale word behind.
    EXPECT_TRUE(uint128_t(0ull, 0x8000000000000000ull) % uint128_t(1ull, 1ull) == uint128_t(0x8000000000000001ull, 0ull));
    // 4 word denominator, the case that always worked
    EXPECT_TRUE(uint128_t(0ull, 0x8000000000000000ull) % uint128_t(1ull, 0x100000000ull) == uint128_t(0xFFFFFFFF80000001ull, 0xFFFFFFFFull));
    // denominator larger than the numerator, and equal to it
    EXPECT_TRUE(uint128_t(1ull, 1ull) % uint128_t(2ull, 1ull) == uint128_t(1ull, 1ull));
    EXPECT_TRUE((uint128_t(1ull, 1ull) % uint128_t(1ull, 1ull)).is_zero());
}
// operator*= read rhs only after it had overwritten this object's QWORDs, so it fed the partial
// result back into the cross products whenever rhs aliased *this. The friend operator* is immune
// because it copies the left operand, which is why MultiplyByUint128 never caught it.
TEST(uint128_t, MultiplyInPlaceAliasing)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const uint128_t x(get_uint64_random(), get_uint64_random());
        uint128_t y = x;
        y *= y;
        EXPECT_TRUE(y == x * x) << "x=" << (std::string)x;
        EXPECT_TRUE(y == sqr(x)) << "x=" << (std::string)x;
    }
}
TEST(uint128_t, MultiplyInPlaceAliasingEdgeCases)
{
    const uint128_t values[] = {uint128_t(0ull, 0ull),  uint128_t(1ull, 0ull),  uint128_t(0ull, 1ull),
                                uint128_t(1ull, 1ull),  uint128_t(~0ull, 0ull), uint128_t(0ull, ~0ull),
                                uint128_t(~0ull, ~0ull), uint128_t(2ull, 3ull)};
    for (const auto& x : values) {
        uint128_t y = x;
        y *= y;
        EXPECT_TRUE(y == x * x) << "x=" << (std::string)x;
    }
    // (2^64 + 1)^2 is 2^128 + 2^65 + 1, which truncates to 2^65 + 1
    uint128_t v(1ull, 1ull);
    v *= v;
    EXPECT_TRUE(v == uint128_t(1ull, 2ull)) << "v=" << (std::string)v;
}
// operator^ was implemented with &=, so it computed a bitwise AND.
TEST(uint128_t, BitwiseOperators)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const uint64_t al = get_uint64_random(), ah = get_uint64_random();
        const uint64_t bl = get_uint64_random(), bh = get_uint64_random();
        const uint128_t a(al, ah), b(bl, bh);

        EXPECT_TRUE((a ^ b) == uint128_t(al ^ bl, ah ^ bh)) << "a=" << (std::string)a << ", b=" << (std::string)b;
        EXPECT_TRUE((a & b) == uint128_t(al & bl, ah & bh)) << "a=" << (std::string)a << ", b=" << (std::string)b;
        EXPECT_TRUE((a | b) == uint128_t(al | bl, ah | bh)) << "a=" << (std::string)a << ", b=" << (std::string)b;

        // the compound form and the binary form must agree
        uint128_t c = a;
        c ^= b;
        EXPECT_TRUE(c == (a ^ b)) << "a=" << (std::string)a << ", b=" << (std::string)b;
        // x ^ x is always zero, x ^ 0 is always x
        EXPECT_TRUE((a ^ a).is_zero()) << "a=" << (std::string)a;
        EXPECT_TRUE((a ^ uint128_t()) == a) << "a=" << (std::string)a;
    }
}
// operator char*() put the terminator at index 32 of a 45 byte buffer and wrote the digits
// backwards from there, so any value of 32 digits or more wrote past the front of the buffer.
// 2^128-1 needs 39 digits and used to start 8 bytes below the buffer.
TEST(uint128_t, ToStringLargeValues)
{
    const char* values[] = {
        "18446744073709551616",                     // 2^64, the smallest value using the long division path
        "10000000000000000000000000000000",         // 10^31, 32 digits
        "99999999999999999999999999999999",         // 32 digits
        "170141183460469231731687303715884105728",  // 2^127, 39 digits
        "340282366920938463463374607431768211455"   // 2^128-1, 39 digits, the longest possible
    };
    // A value whose high QWORD is zero takes the snprintf path, which returns the start of the
    // static buffer. That yields its address without assuming where the terminator sits, so the
    // check below stays valid if the buffer layout changes. Comparing the strings alone is not
    // enough: the old code wrote out of bounds yet still returned the correct digits.
    const char* buffer_start = static_cast<char*>(uint128_t(1ull));
    for (auto i = 0u; i < array_length(values); ++i) {
        const uint128_t v = values[i];
        const char* p = static_cast<char*>(v);
        EXPECT_STREQ(p, values[i]);
        EXPECT_GE(p, buffer_start) << "wrote " << strlen(values[i]) << " digits before the start of the buffer";
    }

    // full round trip over the whole 128 bit range
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const uint128_t v(get_uint64_random(), get_uint64_random());
        const std::string s = (std::string)v;
        EXPECT_TRUE(uint128_t(s) == v) << "s=" << s;
    }
}
// A comma in the input used to continue the parse loop without advancing the read pointer,
// which spun forever. Commas are documented as allowed inside the number.
// The user defined literal in uint128_t.h was never exercised. It is a raw literal operator, so it
// receives the source characters of the token: the compiler does not strip the C++ digit separator,
// and the parser used to stop at the first apostrophe and return 1.
TEST(uint128_t, UserDefinedLiteral)
{
    EXPECT_TRUE(1234_uint128 == uint128_t(1234ull));
    EXPECT_TRUE(0xDEADBEAF_uint128 == uint128_t(0xDEADBEAFull));
    EXPECT_TRUE(340282366920938463463374607431768211455_uint128 == uint128_t(~0ull, ~0ull));
    EXPECT_TRUE(1'000'000_uint128 == uint128_t(1000000ull));
    EXPECT_TRUE(340'282'366'920'938'463'463'374'607'431'768'211'455_uint128 == uint128_t(~0ull, ~0ull));
}
TEST(uint128_t, ConstructorFromStringWithCommas)
{
    EXPECT_TRUE(uint128_t("1,000") == uint128_t(1000ull));
    // the C++ digit separator is accepted alongside the comma, and the two can be mixed
    EXPECT_TRUE(uint128_t("1'000'000") == uint128_t(1000000ull));
    EXPECT_TRUE(uint128_t("1'000,000") == uint128_t(1000000ull));
    EXPECT_TRUE(uint128_t("0xDE'AD'BE'AF") == uint128_t(0xDEADBEAFull));
    EXPECT_TRUE(uint128_t("1,000,000") == uint128_t(1000000ull));
    EXPECT_TRUE(uint128_t("1,") == uint128_t(1ull));
    EXPECT_TRUE(uint128_t(",,,42") == uint128_t(42ull));
    EXPECT_TRUE(uint128_t("0xDEAD,BEAF") == uint128_t(0xDEADBEAFull));
    EXPECT_TRUE(uint128_t("340,282,366,920,938,463,463,374,607,431,768,211,455") == uint128_t(~0ull, ~0ull));
}
// A leading sign used to stop the parse before it began, so "-5" produced zero while -5.0 and
// (int64_t)-5 both produced 2^128-5. A negative string now wraps around the same way.
TEST(uint128_t, ConstructorFromSignedString)
{
    EXPECT_TRUE(uint128_t("-5") == uint128_t(-5.0));
    EXPECT_TRUE(uint128_t("-5") == uint128_t(static_cast<int64_t>(-5)));
    EXPECT_TRUE(uint128_t("-1") == uint128_t(~0ull, ~0ull));
    EXPECT_TRUE(uint128_t("+5") == uint128_t(5ull));
    EXPECT_TRUE(uint128_t("-0").is_zero());
    EXPECT_TRUE(uint128_t("-").is_zero());       // a sign with no digits is still zero
    EXPECT_TRUE(uint128_t("   -5") == uint128_t(static_cast<int64_t>(-5)));

    // the sign is consumed before the base prefix, so hex still works
    EXPECT_TRUE(uint128_t("-0x10") == uint128_t(static_cast<int64_t>(-16)));
    // and it composes with digit grouping
    EXPECT_TRUE(uint128_t("-1,000") == uint128_t(static_cast<int64_t>(-1000)));

    // the whole negative int64_t range round trips through its decimal string
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const int64_t value = -llabs(get_int64_random());
        EXPECT_TRUE(uint128_t(std::to_string(value)) == uint128_t(value)) << "value=" << value;
    }
}
// The parser used to hand a raw char to isspace. Values above 0x7F become negative, which is
// outside the domain isspace accepts and trips an assert in the debug CRT.
TEST(uint128_t, ConstructorFromStringHighByte)
{
    const char high_byte[] = {'1', '2', static_cast<char>(0xB5), '3', 0};
    EXPECT_TRUE(uint128_t(high_byte) == uint128_t(12ull));
    const char leading[] = {static_cast<char>(0xB5), '9', 0};
    EXPECT_TRUE(uint128_t(leading).is_zero());
}
// The conversions used to truncate the bits that did not fit the mantissa. IEEE 754 and the
// builtin integer to floating point conversions round to nearest, ties to even. strtod and strtof
// are correctly rounded by the C standard, so the decimal string of a value is an exact reference.
TEST(uint128_t, ConversionToDoubleIsCorrectlyRounded)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        // shifting the high QWORD spreads the values over the whole exponent range
        const uint128_t v(get_uint64_random(), get_uint64_random() >> (get_uint32_random() % 64));
        if (v.is_zero())
            continue;
        const std::string s = (std::string)v;
        EXPECT_EQ(static_cast<double>(v), strtod(s.c_str(), nullptr)) << "value=" << s;
    }
}
TEST(uint128_t, ConversionToFloatIsCorrectlyRounded)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const uint128_t v(get_uint64_random(), get_uint64_random() >> (get_uint32_random() % 64));
        if (v.is_zero())
            continue;
        const std::string s = (std::string)v;
        EXPECT_EQ(static_cast<float>(v), strtof(s.c_str(), nullptr)) << "value=" << s;
    }
}
TEST(uint128_t, ConversionToFloatingPointEdgeCases)
{
    // exact powers of two need no rounding
    EXPECT_EQ(static_cast<double>(uint128_t(0ull, 1ull)), 18446744073709551616.0);
    // 2^65-1 is not a tie: all twelve dropped bits are set, which is above half, so it rounds up to
    // 2^65. ConversionRoundsTiesToEven covers the exact ties.
    EXPECT_EQ(static_cast<double>(uint128_t(~0ull, 1ull)), 36893488147419103232.0);
    // the largest value rounds up past the top of the 128 bit range, which a double can hold
    EXPECT_EQ(static_cast<double>(uint128_t(~0ull, ~0ull)), 340282366920938463463374607431768211456.0);
    // a float cannot, so it saturates to infinity exactly as the builtin conversions do
    EXPECT_TRUE(std::isinf(static_cast<float>(uint128_t(~0ull, ~0ull))));
    // FLT_MAX is (2^24-1) * 2^104 and must stay finite
    EXPECT_EQ(static_cast<float>(uint128_t(0xFFFFFFull) << 104), FLT_MAX);
}
// Negative values used to lose their sign, so -5.0 produced 5. They now wrap around via 2's
// complement, which is what the integer constructors and the builtin unsigned types do.
TEST(uint128_t, ConstructorFromNegativeDouble)
{
    EXPECT_TRUE(uint128_t(-5.0) == uint128_t(static_cast<int64_t>(-5)));
    EXPECT_TRUE(uint128_t(-1.0) == uint128_t(~0ull, ~0ull));
    EXPECT_TRUE(uint128_t(-0.0).is_zero());
    EXPECT_TRUE(uint128_t(-0.5).is_zero());  // truncates towards zero before wrapping

    // A double holds an int64_t exactly only up to 2^53, so the generator stays inside that range: an
    // arithmetic right shift by 11 to 63 bits spreads the magnitudes over every scale up to 2^52. The
    // previous version drew full width values and skipped the ones above 2^53, which left 69 of the
    // 65536 iterations asserting anything.
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const int64_t value = get_int64_random() >> (11u + get_uint32_random() % 53u);
        EXPECT_TRUE(uint128_t(static_cast<double>(value)) == uint128_t(value)) << "value=" << value;
    }
}
// The templated operators cast the operand straight to uint64_t, which is undefined behavior for
// a floating point value at or above 2^64. Floating point operands are now routed through the
// uint128_t constructor, which truncates them towards zero first. A divisor whose magnitude is below
// one therefore converts to zero and throws, exactly as dividing by the converted value does.
TEST(uint128_t, MultiplyAndDivideByFloatingPoint)
{
    EXPECT_THROW((void)(uint128_t(100ull) / 0.5), std::logic_error);
    EXPECT_THROW((void)(uint128_t(100ull) % 0.5), std::logic_error);
    EXPECT_THROW((void)(uint128_t(100ull) / -0.5f), std::logic_error);

    uint128_t big(2ull);
    big *= 1e30;
    EXPECT_TRUE(big == uint128_t(1e30) * uint128_t(2ull));

    uint128_t d(0ull, 1ull);
    d /= 1e19;
    EXPECT_TRUE(d == uint128_t(0ull, 1ull) / uint128_t(1e19));

    // a fractional operand truncates first, matching the uint128_t operand overloads
    uint128_t t(100ull);
    t *= 2.5;
    EXPECT_TRUE(t == uint128_t(200ull));

    uint128_t q(1000ull);
    q /= 4.0f;
    EXPECT_TRUE(q == uint128_t(250ull));
}
// Shifting by 128 or more used to be taken modulo 64 by the hardware and leave bits behind.
TEST(uint128_t, ShiftBeyondWidth)
{
    // volatile keeps the compiler from folding the shift and hiding the defect
    volatile int32_t shift = 128;
    const uint128_t all_ones(~0ull, ~0ull);
    EXPECT_TRUE((all_ones >> shift).is_zero());
    EXPECT_TRUE((all_ones << shift).is_zero());
    shift = 200;
    EXPECT_TRUE((all_ones >> shift).is_zero());
    EXPECT_TRUE((all_ones << shift).is_zero());

    // the boundary either side of it still behaves
    shift = 127;
    EXPECT_TRUE((all_ones >> shift) == uint128_t(1ull));
    EXPECT_TRUE((uint128_t(1ull) << shift) == uint128_t(0ull, 0x8000000000000000ull));
    shift = 0;
    EXPECT_TRUE((all_ones >> shift) == all_ones);
    EXPECT_TRUE((all_ones << shift) == all_ones);
}
// One overload per fixed width type left long and unsigned long ambiguous on MSVC, where both are
// distinct from int and long long. A constrained template covers every builtin integral type.
TEST(uint128_t, ConstructorFromAllIntegralTypes)
{
    EXPECT_TRUE(uint128_t(static_cast<unsigned long>(7)) == uint128_t(7ull));
    EXPECT_TRUE(uint128_t(static_cast<long>(-7)) == uint128_t(static_cast<int64_t>(-7)));
    EXPECT_TRUE(uint128_t(static_cast<unsigned short>(3)) == uint128_t(3ull));
    EXPECT_TRUE(uint128_t(static_cast<short>(-3)) == uint128_t(static_cast<int64_t>(-3)));
    EXPECT_TRUE(uint128_t(static_cast<char>(5)) == uint128_t(5ull));
    EXPECT_TRUE(uint128_t(static_cast<unsigned char>(5)) == uint128_t(5ull));
    EXPECT_TRUE(uint128_t(true) == uint128_t(1ull));
    EXPECT_TRUE(uint128_t(static_cast<size_t>(9)) == uint128_t(9ull));

    // the compound operators accept them too
    uint128_t x(10ull);
    x += static_cast<unsigned long>(5);
    x -= static_cast<long>(-5);
    EXPECT_TRUE(x == uint128_t(20ull));
}
// A scalar on the left used to be ambiguous: converting it to uint128_t and converting the
// uint128_t to a builtin type are both a single user defined conversion.
TEST(uint128_t, ScalarOnTheLeftHandSide)
{
    const uint128_t x(10ull);
    EXPECT_TRUE((1 + x) == uint128_t(11ull));
    EXPECT_TRUE((100 - x) == uint128_t(90ull));
    EXPECT_TRUE((3 * x) == uint128_t(30ull));
    EXPECT_TRUE((100 / x) == uint128_t(10ull));
    EXPECT_TRUE((105 % x) == uint128_t(5ull));
    EXPECT_TRUE((0xFF & x) == uint128_t(10ull));
    EXPECT_TRUE((1 | x) == uint128_t(11ull));
    EXPECT_TRUE((3 ^ x) == uint128_t(9ull));

    // the uint128_t on the left forms must still resolve, including both operands being uint128_t
    EXPECT_TRUE((x + 1) == uint128_t(11ull));
    EXPECT_TRUE((x + x) == uint128_t(20ull));
    EXPECT_TRUE((x * x) == uint128_t(100ull));

    // The left operand is widened rather than the object being narrowed, so the result is 128 bit.
    static_assert(std::is_same_v<decltype(1 + x), uint128_t>);
    static_assert(std::is_same_v<decltype(1 ^ x), uint128_t>);
}
// The comparisons only existed with the uint128_t on the left. A scalar there had nothing but the
// builtin comparisons to choose from, and the conversion operators of uint128_t make every one of
// them equally good, so the call was ambiguous.
TEST(uint128_t, ScalarOnTheLeftHandSideComparisons)
{
    const uint128_t x(10ull);
    EXPECT_TRUE(10 == x);
    EXPECT_FALSE(9 == x);
    EXPECT_TRUE(9 != x);
    EXPECT_FALSE(10 != x);
    EXPECT_TRUE(9 < x);
    EXPECT_FALSE(10 < x);
    EXPECT_TRUE(9 <= x);
    EXPECT_TRUE(10 <= x);
    EXPECT_TRUE(11 > x);
    EXPECT_FALSE(10 > x);
    EXPECT_TRUE(11 >= x);
    EXPECT_TRUE(10 >= x);

    // A value no builtin type can hold must not compare through a narrowed conversion.
    const uint128_t big = uint128_t(1ull) << 100;
    EXPECT_TRUE(1 < big);
    EXPECT_FALSE(1 > big);
    EXPECT_TRUE(1 != big);

    // the uint128_t on the left forms must still resolve
    EXPECT_TRUE(x == 10);
    EXPECT_TRUE(x > 9);
    EXPECT_TRUE(x == uint128_t(10ull));
}
// Widening the left operand is what makes (1 << 100) produce the expected power of two. The same
// expression on a builtin int is undefined behavior, the shift count being wider than the operand.
TEST(uint128_t, ScalarOnTheLeftHandSideShifts)
{
    const uint128_t ten(10ull);
    EXPECT_TRUE((1 << ten) == uint128_t(1024ull));
    EXPECT_TRUE((1024 >> ten) == uint128_t(1ull));
    EXPECT_TRUE((3 << uint128_t(2ull)) == uint128_t(12ull));

    const uint128_t hundred(100ull);
    EXPECT_TRUE((1 << hundred) == (uint128_t(1ull) << 100));
    EXPECT_TRUE((1 << hundred) != 0ull);
    static_assert(std::is_same_v<decltype(1 << hundred), uint128_t>);

    // the uint128_t on the left forms must still resolve
    EXPECT_TRUE((ten << 1) == uint128_t(20ull));
    EXPECT_TRUE((ten >> 1) == uint128_t(5ull));
}
// The documentation claimed these return zero for a zero input, they throw.
TEST(uint128_t, LogFunctionsRejectZero)
{
    const uint128_t zero;
    EXPECT_THROW((void)log2(zero), std::domain_error);
    EXPECT_THROW((void)log(zero), std::domain_error);
    EXPECT_THROW((void)log10(zero), std::domain_error);
    // sqrt is documented to return zero rather than throw
    EXPECT_EQ(sqrt(zero), 0ull);
    EXPECT_EQ(sqrt(uint128_t(1ull)), 1ull);
}
// The top three entries of the natural log table had been produced with a floating point pow(),
// which loses its low digits well before e^86. e^86 and e^87 came out one and three too large, so
// log() answered 85 and 86 at those boundaries, and e^88 came out 703476075845531179681763 too
// small, so log() answered 88 across that whole range of values whose log is really 87.
//
// Nothing in the random log test can see this: a uniform 128 bit value lands in the widest of
// those windows with a probability of 2e-15. The boundaries have to be named explicitly.
//
// Index k of the table holds ceil(e^k), the smallest integer whose log is k, which pins the entry
// exactly: log() has to answer k there and k-1 one below it.
TEST(uint128_t, LogBoundaries)
{
    // ceil(e^k) computed with exact arithmetic, paired with its exponent
    const struct {
        const char* boundary;
        uint64_t expo;
    } boundaries[] = {
        {"3", 1},
        {"8", 2},
        {"22027", 10},
        {"235385266837019986", 40},
        {"2293783159469609879099352841", 63},
        {"6235149080811616882909238709", 64},
        {"3025077322201142338266566396443428743", 84},
        {"8223012714622913510304328016407774696", 85},
        {"22352466037347150474430657323327147399", 86},
        {"60760302250568721495223289381302760753", 87},
        {"165163625499400185552832979626485876707", 88},
    };

    for (const auto& b : boundaries) {
        const uint128_t x(b.boundary);
        EXPECT_EQ(log(x), b.expo) << "at ceil(e^" << b.expo << ") = " << b.boundary;
        EXPECT_EQ(log(x - uint128_t(1ull)), b.expo - 1) << "just below ceil(e^" << b.expo << ") = " << b.boundary;
    }

    // the three values that used to be misreported, one from each broken entry
    EXPECT_EQ(log(uint128_t("22352466037347150474430657323327147399")), 86ull);
    EXPECT_EQ(log(uint128_t("60760302250568721495223289381302760755")), 87ull);
    EXPECT_EQ(log(uint128_t("165163625499399482076757134095306194944")), 87ull);

    // 88 is the largest exponent whose power of e fits, everything above the last boundary
    // saturates there rather than running off the end of the table
    EXPECT_EQ(log(uint128_t(~0ull, ~0ull)), 88ull);
}
// uint128_t and int128_t are now aliases of a single class template. Nothing but this test pins the
// properties the merge had to preserve: the storage is unchanged, the 16 byte alignment survived the
// template, and the two types still refuse to convert into one another.
//
// The type is deliberately *not* trivially copyable. The copy and move members were defaulted until
// the 128 bit multiplication benchmark showed what that costs: MSVC implements a trivial 16 byte
// copy as one `vmovups`, which overlaps the two QWORD stores every operator uses to build its
// result, fails store to load forwarding, and stalls. Writing the members out one at a time removes
// the stack temporary altogether and multiplies six times faster. The assertion below is what keeps
// the members from drifting back to `= default`.
TEST(uint128_t, LayoutAndTraits)
{
    static_assert(std::is_same_v<uint128_t, int128_base<false>>, "uint128_t is the unsigned instantiation");
    static_assert(uint128_t::is_signed == false, "uint128_t reports itself as unsigned");
    static_assert(sizeof(uint128_t) == 16, "the object holds exactly two QWORDs");
    static_assert(alignof(uint128_t) == 16, "FP128_ALIGN16 has to survive the template");
    static_assert(std::is_standard_layout_v<uint128_t>, "the division helpers alias the members as an array");
    static_assert(!std::is_trivially_copyable_v<uint128_t>, "the copy and move members are written out, not defaulted");
    // getting from one signedness to the other needs two user defined conversions, which the
    // language never performs implicitly, and leaves the explicit form ambiguous
    static_assert(!std::is_convertible_v<int128_t, uint128_t>, "no silent signed to unsigned conversion");
    static_assert(!std::is_constructible_v<uint128_t, int128_t>, "and no unambiguous explicit one either");
    SUCCEED();
}
// Every math function on these types is a hidden friend, found by argument dependent lookup. A value
// whose low QWORD is zero catches the failure mode where a call falls through to the narrower
// fp128::log2(uint64_t) overload through the implicit conversion to uint64_t: that overload would
// see a zero and answer zero.
TEST(uint128_t, MathFunctionsSeeTheWholeValue)
{
    const uint128_t x(0ull, 1ull << 40);  // 2^104, the low QWORD is zero

    EXPECT_EQ(log2(x), 104ull);
    EXPECT_EQ(lzcnt128(x), 23ull);
    EXPECT_EQ(log10(x), 31ull);
    EXPECT_EQ(log(x), 72ull);
    EXPECT_EQ(sqrt(x), 1ull << 52);
    EXPECT_TRUE(pow(uint128_t(2ull), 104u) == x);
    EXPECT_TRUE(sqr(uint128_t(1ull << 52)) == x);
}
// to_string() was declared on int128_t only.
TEST(uint128_t, ToStringMatchesTheStringConversion)
{
    const uint128_t x("340282366920938463463374607431768211455");  // 2^128 - 1
    const std::string via_operator = static_cast<std::string>(x);
    EXPECT_STREQ(x.to_string(), "340282366920938463463374607431768211455");
    EXPECT_STREQ(x.to_string(), via_operator.c_str());
}
// operator/= detects a power of two divisor and shifts instead of dividing. The shift is a logical
// one applied to the magnitude, so it has to be correct on both sides of the 64 bit boundary and
// for a divisor of exactly one, which shifts by zero.
TEST(uint128_t, DivideByPowersOfTwo)
{
    const uint128_t x("340282366920938463463374607431768211455");  // 2^128 - 1

    EXPECT_TRUE((x / uint128_t(1ull)) == x);
    EXPECT_TRUE((x / uint128_t(2ull)) == (x >> 1));
    EXPECT_TRUE((x / uint128_t(1ull << 63)) == (x >> 63));
    EXPECT_TRUE((x / uint128_t(0ull, 1ull)) == (x >> 64));         // 2^64
    EXPECT_TRUE((x / uint128_t(0ull, 1ull << 36)) == (x >> 100));  // 2^100
    EXPECT_TRUE((x / uint128_t(0ull, 1ull << 63)) == uint128_t(1ull));
}

/**********************************************************************
 * uint128_t paths the tests above never reached
 *
 * Found with a coverage run of the whole suite: each test below drives code
 * that executed zero times, or edge cases nothing above could construct.
 * Fixed expected values come from Python's arbitrary precision integers
 * (and its correctly rounded int to float conversion), never from the
 * library itself.
 ***********************************************************************/

// Division by a builtin integer has a path of its own for the unsigned type. It writes the quotient
// over the dividend in place, and has shortcuts for a dividend below or equal to the divisor and for
// a power of two. None of it ran anywhere in the suite: the only integral divisor was a zero that
// throws first.
TEST(uint128_t, DivideByScalar)
{
    const struct {
        uint64_t aLow, aHigh, divisor, qLow, qHigh, remainder;
    } rows[] = {
        {0xFFFFFFFFFFFFFFFFull, 0xFFFFFFFFFFFFFFFFull, 0x0000000000000003ull, 0x5555555555555555ull, 0x5555555555555555ull, 0x0000000000000000ull},
        {0xFFFFFFFFFFFFFFFFull, 0xFFFFFFFFFFFFFFFFull, 0x000000000000000Aull, 0x9999999999999999ull, 0x1999999999999999ull, 0x0000000000000005ull},
        {0xFFFFFFFFFFFFFFFFull, 0xFFFFFFFFFFFFFFFFull, 0xFFFFFFFFFFFFFFFFull, 0x0000000000000001ull, 0x0000000000000001ull, 0x0000000000000000ull},
        {0xFFFFFFFFFFFFFFFFull, 0xFFFFFFFFFFFFFFFFull, 0x8000000000000001ull, 0xFFFFFFFFFFFFFFFCull, 0x0000000000000001ull, 0x0000000000000003ull},
        {0xFFFFFFFFFFFFFFFFull, 0xFFFFFFFFFFFFFFFFull, 0x0000000000000001ull, 0xFFFFFFFFFFFFFFFFull, 0xFFFFFFFFFFFFFFFFull, 0x0000000000000000ull},
        {0xFFFFFFFFFFFFFFFFull, 0xFFFFFFFFFFFFFFFFull, 0x0000000000000002ull, 0xFFFFFFFFFFFFFFFFull, 0x7FFFFFFFFFFFFFFFull, 0x0000000000000001ull},
        {0xFFFFFFFFFFFFFFFFull, 0xFFFFFFFFFFFFFFFFull, 0x8000000000000000ull, 0xFFFFFFFFFFFFFFFFull, 0x0000000000000001ull, 0x7FFFFFFFFFFFFFFFull},
        {0x0123456789ABCDEFull, 0xFEDCBA9876543210ull, 0x0000000000000007ull, 0x494E2E7C8161AFB4ull, 0x2468ACF13579BE02ull, 0x0000000000000003ull},
        {0x0123456789ABCDEFull, 0xFEDCBA9876543210ull, 0x00000000FFFFFFFFull, 0x7530ECA976543211ull, 0x00000000FEDCBA99ull, 0x0000000000000000ull},
        {0x0123456789ABCDEFull, 0xFEDCBA9876543210ull, 0x0000000100000001ull, 0x7777777889ABCDEFull, 0x00000000FEDCBA97ull, 0x0000000000000000ull},
        {0x0123456789ABCDEFull, 0xFEDCBA9876543210ull, 0xFFFFFFFFFFFFFFFFull, 0xFEDCBA9876543211ull, 0x0000000000000000ull, 0x0000000000000000ull},
        {0x0000000000000000ull, 0x0000000000000001ull, 0x0000000000000003ull, 0x5555555555555555ull, 0x0000000000000000ull, 0x0000000000000001ull},
        {0x0000000000000005ull, 0x0000000000000001ull, 0xFFFFFFFFFFFFFFFFull, 0x0000000000000001ull, 0x0000000000000000ull, 0x0000000000000006ull},
        {0x0000000000003039ull, 0x0000000000000000ull, 0x0000000000003039ull, 0x0000000000000001ull, 0x0000000000000000ull, 0x0000000000000000ull},  // equal
        {0x0000000000003038ull, 0x0000000000000000ull, 0x0000000000003039ull, 0x0000000000000000ull, 0x0000000000000000ull, 0x0000000000003038ull},  // smaller
        {0x0000000000000000ull, 0x0000000000000000ull, 0x0000000000000007ull, 0x0000000000000000ull, 0x0000000000000000ull, 0x0000000000000000ull},  // zero
        {0xFFFFFFFFFFFFFFFFull, 0x0000000000000000ull, 0xFFFFFFFFFFFFFFFEull, 0x0000000000000001ull, 0x0000000000000000ull, 0x0000000000000001ull},
    };
    for (const auto& row : rows) {
        const uint128_t a(row.aLow, row.aHigh);
        const uint128_t expected(row.qLow, row.qHigh);
        EXPECT_TRUE(a / row.divisor == expected) << "a=" << (std::string)a << ", divisor=" << row.divisor;
        EXPECT_TRUE(a % row.divisor == uint128_t(row.remainder, 0ull)) << "a=" << (std::string)a << ", divisor=" << row.divisor;
        uint128_t c = a;
        c /= row.divisor;
        EXPECT_TRUE(c == expected) << "a=" << (std::string)a << ", divisor=" << row.divisor;
        // the narrower builtin types reach the same path through their own instantiations
        if (row.divisor <= UINT32_MAX) {
            EXPECT_TRUE(a / static_cast<uint32_t>(row.divisor) == expected) << "a=" << (std::string)a << ", divisor=" << row.divisor;
        }
        if (row.divisor <= static_cast<uint64_t>(INT32_MAX)) {
            EXPECT_TRUE(a / static_cast<int32_t>(row.divisor) == expected) << "a=" << (std::string)a << ", divisor=" << row.divisor;
        }
    }

    // A negative divisor is sign extended, so for the unsigned type it wraps around to 2^128 - |d|,
    // just as a builtin unsigned division by a negative int converts the int first.
    const struct {
        uint64_t aLow, aHigh;
        int64_t divisor;
        uint64_t qLow, qHigh, rLow, rHigh;
    } negative[] = {
        {0x0000000000000064ull, 0x0000000000000000ull, -1ll, 0x0000000000000000ull, 0x0000000000000000ull, 0x0000000000000064ull, 0x0000000000000000ull},
        {0xFFFFFFFFFFFFFFFFull, 0xFFFFFFFFFFFFFFFFull, -1ll, 0x0000000000000001ull, 0x0000000000000000ull, 0x0000000000000000ull, 0x0000000000000000ull},
        {0xFFFFFFFFFFFFFFFFull, 0xFFFFFFFFFFFFFFFFull, -2ll, 0x0000000000000001ull, 0x0000000000000000ull, 0x0000000000000001ull, 0x0000000000000000ull},
        {0xFFFFFFFFFFFFFFFEull, 0xFFFFFFFFFFFFFFFFull, -1ll, 0x0000000000000000ull, 0x0000000000000000ull, 0xFFFFFFFFFFFFFFFEull, 0xFFFFFFFFFFFFFFFFull},
        {0x0123456789ABCDEFull, 0xFEDCBA9876543210ull, -7ll, 0x0000000000000000ull, 0x0000000000000000ull, 0x0123456789ABCDEFull, 0xFEDCBA9876543210ull},
        {0x0000000000000000ull, 0x8000000000000000ull, INT64_MIN, 0x0000000000000000ull, 0x0000000000000000ull, 0x0000000000000000ull, 0x8000000000000000ull},
    };
    for (const auto& row : negative) {
        const uint128_t a(row.aLow, row.aHigh);
        EXPECT_TRUE(a / row.divisor == uint128_t(row.qLow, row.qHigh)) << "a=" << (std::string)a << ", divisor=" << row.divisor;
        EXPECT_TRUE(a % row.divisor == uint128_t(row.rLow, row.rHigh)) << "a=" << (std::string)a << ", divisor=" << row.divisor;
    }

    const uint128_t max(~0ull, ~0ull);
    EXPECT_THROW((void)(max / 0), std::logic_error);
    EXPECT_THROW((void)(max % 0ull), std::logic_error);

    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        // one dividend in four fits in 64 bits, which lets the shortcuts for a small dividend fire
        const uint64_t aHigh = ((i & 3u) == 0u) ? 0ull : get_uint64_random() >> (get_uint32_random() % 64u);
        const uint128_t a(get_uint64_random(), aHigh);
        // the forced top bit keeps the divisor non zero whatever the shift; one in eight is a power of two
        uint64_t divisor = (get_uint64_random() | (1ull << 63)) >> (get_uint32_random() % 64u);
        if ((i & 7u) == 1u) {
            divisor = 1ull << (get_uint32_random() % 64u);
        }

        const uint128_t q = a / divisor;
        const uint128_t r = a % divisor;
        EXPECT_TRUE(IsExactDivision(a, uint128_t(divisor, 0ull), q, r)) << "a=" << (std::string)a << ", divisor=" << divisor;
        // the in place path and the int128_base path are separate code and have to agree
        EXPECT_TRUE(q == a / uint128_t(divisor, 0ull)) << "a=" << (std::string)a << ", divisor=" << divisor;
        uint128_t c = a;
        c /= divisor;
        EXPECT_TRUE(c == q) << "a=" << (std::string)a << ", divisor=" << divisor;
    }
}
// Anything at or beyond 2^128 saturates to the maximum, NaN and both infinities included, whatever
// the sign: the branch that does it never ran. Within range a negative value wraps around via 2's
// complement. The expected values are exact powers of two and sums of them.
TEST(uint128_t, ConstructorFromDoubleSaturatesAndWraps)
{
    const uint128_t max(~0ull, ~0ull);
    const double inf = std::numeric_limits<double>::infinity();
    const double saturating[] = {inf, -inf, std::numeric_limits<double>::quiet_NaN(), 1e40, -1e40, std::ldexp(1.0, 128)};
    for (const double value : saturating) {
        EXPECT_TRUE(uint128_t(value) == max) << "value=" << value;
    }
    EXPECT_TRUE(uint128_t(std::numeric_limits<float>::infinity()) == max);

    // the largest double below 2^128 is (2^53 - 1) * 2^75, which converts exactly
    const double largest = std::nextafter(std::ldexp(1.0, 128), 0.0);
    EXPECT_TRUE(uint128_t(largest) == uint128_t(0ull, 0xFFFFFFFFFFFFF800ull));
    EXPECT_TRUE(uint128_t(-largest) == uint128_t(0ull, 0x800ull));  // 2^128 minus it is 2^75
    EXPECT_TRUE(uint128_t(std::ldexp(1.0, 127)) == uint128_t(0ull, 1ull << 63));
    EXPECT_TRUE(uint128_t(-std::ldexp(1.0, 127)) == uint128_t(0ull, 1ull << 63));
    EXPECT_TRUE(uint128_t(-std::ldexp(1.0, 64)) == uint128_t(0ull, ~0ull));

    // magnitudes below one truncate to zero
    EXPECT_TRUE(uint128_t(DBL_MIN).is_zero());
    EXPECT_TRUE(uint128_t(std::numeric_limits<double>::denorm_min()).is_zero());
    EXPECT_TRUE(uint128_t(0.9999).is_zero());
}
// The conversions round to nearest with ties to even, and a tie is the one case random operands never
// produce: every dropped bit but the highest has to be zero. Each row is an exact tie, or one away from
// it, with the deciding bits placed in each of the branches RoundedMantissa() has: a single dropped bit,
// dropped bits within the low QWORD, the round bit at bit 64, and dropped bits spanning both QWORDs.
TEST(uint128_t, ConversionRoundsTiesToEven)
{
    const struct {
        uint64_t low, high, doubleBits;
        uint32_t floatBits;
    } rows[] = {
        {0x0020000000000001ull, 0x0000000000000000ull, 0x4340000000000000ull, 0x5A000000u},  // 2^53+1: one dropped bit, a tie, rounds down to even
        {0x0020000000000003ull, 0x0000000000000000ull, 0x4340000000000002ull, 0x5A000000u},  // 2^53+3: one dropped bit, a tie, rounds up to even
        {0x0040000000000002ull, 0x0000000000000000ull, 0x4350000000000000ull, 0x5A800000u},  // 2^54+2: two dropped bits, a tie, rounds down
        {0x0040000000000006ull, 0x0000000000000000ull, 0x4350000000000002ull, 0x5A800000u},  // 2^54+6: two dropped bits, a tie, rounds up
        {0x0000000001000001ull, 0x0000000000000000ull, 0x4170000010000000ull, 0x4B800000u},  // 2^24+1: float tie on one dropped bit
        {0x0000000001000003ull, 0x0000000000000000ull, 0x4170000030000000ull, 0x4B800002u},  // 2^24+3: float tie, rounds up
        {0x0000000000000800ull, 0x0000000000000001ull, 0x43F0000000000000ull, 0x5F800000u},  // 2^64+2^11: double tie within the low QWORD
        {0x0000000000000801ull, 0x0000000000000001ull, 0x43F0000000000001ull, 0x5F800000u},  // one above it: the sticky bit decides
        {0x0000000000000000ull, 0x0020000000000001ull, 0x4740000000000000ull, 0x7A000000u},  // 2^117+2^64: round bit is bit 64
        {0x0000000000000001ull, 0x0020000000000001ull, 0x4740000000000001ull, 0x7A000000u},  // one above it: the sticky bit is the low QWORD
        {0x0000000000000000ull, 0x0020000000000003ull, 0x4740000000000002ull, 0x7A000000u},  // round bit at bit 64, rounds up to even
        {0x0000000000000000ull, 0x8000000000000400ull, 0x47E0000000000000ull, 0x7F000000u},  // 2^127+2^74: tie spanning both QWORDs
        {0x0000000000000001ull, 0x8000000000000400ull, 0x47E0000000000001ull, 0x7F000000u},  // sticky bit only in the low QWORD
        {0x0000000000000000ull, 0x8000000000000401ull, 0x47E0000000000001ull, 0x7F000000u},  // sticky bit only in the high QWORD
        {0x0000000000000000ull, 0x8000008000000000ull, 0x47E0000010000000ull, 0x7F000000u},  // 2^127+2^103: float tie near the top
        {0x0000000000000000ull, 0x8000018000000000ull, 0x47E0000030000000ull, 0x7F000002u},  // 2^127+3*2^103: float tie, rounds up
        {0x0000000000000000ull, 0xFFFFFF8000000000ull, 0x47EFFFFFF0000000ull, 0x7F800000u},  // 2^128-2^103: float tie past FLT_MAX, +inf
        {0xFFFFFFFFFFFFFFFFull, 0xFFFFFF7FFFFFFFFFull, 0x47EFFFFFF0000000ull, 0x7F7FFFFFu},  // one below that tie: FLT_MAX
        {0xFFFFFFFFFFFFFFFFull, 0x0000000000000001ull, 0x4400000000000000ull, 0x60000000u},  // 2^65-1: not a tie, rounds up
    };
    for (const auto& row : rows) {
        const uint128_t x(row.low, row.high);
        EXPECT_EQ(std::bit_cast<uint64_t>(static_cast<double>(x)), row.doubleBits) << "x=" << (std::string)x;
        EXPECT_EQ(std::bit_cast<uint32_t>(static_cast<float>(x)), row.floatBits) << "x=" << (std::string)x;
    }
}
// The shifts are compared with a reference built one bit at a time, over every count from -2 to 130:
// zero and negative counts leave the value unchanged, 128 and above clear it. Nothing checked counts
// 1 to 127 at runtime against anything but other shifts, or negative counts on the unsigned type.
TEST(uint128_t, ShiftsMatchABitReference)
{
    const uint128_t values[] = {uint128_t(~0ull, ~0ull), uint128_t(0x0123456789ABCDEFull, 0xFEDCBA9876543210ull), uint128_t(1ull),
                                uint128_t(0ull, 1ull << 63), uint128_t(1ull << 63, 1ull)};
    for (const auto& x : values) {
        for (int32_t s = -2; s <= 130; ++s) {
            // volatile keeps the compiler from folding the count and testing its own arithmetic
            volatile int32_t count = s;
            const uint128_t left = ReferenceShift(x, s, true);
            const uint128_t right = ReferenceShift(x, s, false);
            EXPECT_TRUE((x << count) == left) << "x=" << (std::string)x << ", shift=" << s;
            EXPECT_TRUE((x >> count) == right) << "x=" << (std::string)x << ", shift=" << s;
            uint128_t c = x;
            c <<= count;
            EXPECT_TRUE(c == left) << "x=" << (std::string)x << ", shift=" << s;
            c = x;
            c >>= count;
            EXPECT_TRUE(c == right) << "x=" << (std::string)x << ", shift=" << s;
        }
    }

    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const uint128_t x(get_uint64_random(), get_uint64_random());
        volatile int32_t count = static_cast<int32_t>(get_uint32_random() % 133u) - 2;
        const int32_t s = count;
        EXPECT_TRUE((x << count) == ReferenceShift(x, s, true)) << "x=" << (std::string)x << ", shift=" << s;
        EXPECT_TRUE((x >> count) == ReferenceShift(x, s, false)) << "x=" << (std::string)x << ", shift=" << s;
    }
}
// Multiplication by a builtin integer skips a partial product for a non negative scalar and widens a
// negative one first. The decimal parser was its only caller, always with a factor of ten. The
// reference widens the scalar from its QWORDs, so it does not depend on the integral constructor.
TEST(uint128_t, MultiplyByScalar)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const uint128_t a(get_uint64_random(), get_uint64_random());
        const uint64_t u = get_uint64_random() >> (get_uint32_random() % 64u);
        // both signs, including magnitudes whose top bit is set
        const int64_t s = ((i & 1u) == 0u) ? static_cast<int64_t>(u) : static_cast<int64_t>(0ull - u);
        const uint32_t w = static_cast<uint32_t>(u);

        const uint128_t expectedU = ReferenceMultiply(a, uint128_t(u, 0ull));
        EXPECT_TRUE(a * u == expectedU) << "a=" << (std::string)a << ", u=" << u;
        EXPECT_TRUE(u * a == expectedU) << "a=" << (std::string)a << ", u=" << u;
        uint128_t c = a;
        c *= u;
        EXPECT_TRUE(c == expectedU) << "a=" << (std::string)a << ", u=" << u;

        const uint128_t expectedS = ReferenceMultiply(a, uint128_t(static_cast<uint64_t>(s), SignExtension(s)));
        EXPECT_TRUE(a * s == expectedS) << "a=" << (std::string)a << ", s=" << s;
        c = a;
        c *= s;
        EXPECT_TRUE(c == expectedS) << "a=" << (std::string)a << ", s=" << s;

        EXPECT_TRUE(a * w == ReferenceMultiply(a, uint128_t(static_cast<uint64_t>(w), 0ull))) << "a=" << (std::string)a << ", w=" << w;
    }

    // -1 negates; INT64_MIN times 3 is -(2^64 + 2^63), i.e. 2^128 - 2^64 - 2^63
    EXPECT_TRUE(uint128_t(5ull) * -1 == uint128_t(~0ull - 4ull, ~0ull));
    EXPECT_TRUE(uint128_t(3ull) * INT64_MIN == uint128_t(0x8000000000000000ull, 0xFFFFFFFFFFFFFFFEull));
    // (2^128 - 1)(2^64 - 1) is -(2^64 - 1) modulo 2^128, i.e. 2^128 - 2^64 + 1
    EXPECT_TRUE(uint128_t(~0ull, ~0ull) * ~0ull == uint128_t(1ull, ~0ull));
}
// The carry and borrow across the QWORD boundary and the wrap at both ends of the range, for the
// prefix and postfix forms. The postfix forms return the value from before the step.
TEST(uint128_t, IncrementDecrementBoundaries)
{
    const uint128_t max(~0ull, ~0ull);

    uint128_t a(~0ull, 0ull);  // 2^64 - 1
    EXPECT_TRUE(a++ == uint128_t(~0ull, 0ull));
    EXPECT_TRUE(a == uint128_t(0ull, 1ull));
    EXPECT_TRUE(a-- == uint128_t(0ull, 1ull));
    EXPECT_TRUE(a == uint128_t(~0ull, 0ull));
    EXPECT_TRUE(++a == uint128_t(0ull, 1ull));
    EXPECT_TRUE(--a == uint128_t(~0ull, 0ull));

    uint128_t b = max;
    EXPECT_TRUE(b++ == max);
    EXPECT_TRUE(b.is_zero());
    EXPECT_TRUE(b-- == uint128_t());
    EXPECT_TRUE(b == max);
    EXPECT_TRUE((++b).is_zero());
    EXPECT_TRUE(--b == max);
}
// pow() past 64 bits and through the wrap at 2^128, which the random pow test never reaches. Each
// expected value is pow(base, exponent, 2**128) in Python.
TEST(uint128_t, powEdgeCases)
{
    const struct {
        uint64_t baseLow, baseHigh;
        uint32_t exponent;
        uint64_t resultLow, resultHigh;
    } rows[] = {
        {0x0000000000000002ull, 0x0000000000000000ull, 127u, 0x0000000000000000ull, 0x8000000000000000ull},
        {0x0000000000000002ull, 0x0000000000000000ull, 128u, 0x0000000000000000ull, 0x0000000000000000ull},
        {0x0000000000000003ull, 0x0000000000000000ull, 80u, 0x3CEA59789C79D441ull, 0x6F32F1EF8B18A2BCull},
        {0xFFFFFFFFFFFFFFFFull, 0xFFFFFFFFFFFFFFFFull, 2u, 0x0000000000000001ull, 0x0000000000000000ull},
        {0xFFFFFFFFFFFFFFFFull, 0xFFFFFFFFFFFFFFFFull, 3u, 0xFFFFFFFFFFFFFFFFull, 0xFFFFFFFFFFFFFFFFull},
        {0x000000000000000Aull, 0x0000000000000000ull, 38u, 0x098A224000000000ull, 0x4B3B4CA85A86C47Aull},
        {0x000000000000000Aull, 0x0000000000000000ull, 39u, 0x5F65568000000000ull, 0xF050FE938943ACC4ull},
        {0x0000000000000000ull, 0x0000000000000001ull, 2u, 0x0000000000000000ull, 0x0000000000000000ull},
        {0xFFFFFFFFFFFFFFFFull, 0x0000000000000000ull, 2u, 0x0000000000000001ull, 0xFFFFFFFFFFFFFFFEull},
        {0x0123456789ABCDEFull, 0x0000000000000000ull, 5u, 0x8FC8CE089782DBAFull, 0x94FECC0A73DCE9B0ull},
        {0x0000000000000007ull, 0x0000000000000000ull, 4294967295u, 0x66D8E0D5B6DB6DB7ull, 0x42A5292C21D74697ull},
        {0x0000000000000001ull, 0x0000000000000000ull, 4294967295u, 0x0000000000000001ull, 0x0000000000000000ull},
        {0x0000000000000000ull, 0x0000000000000000ull, 1u, 0x0000000000000000ull, 0x0000000000000000ull},
        {0x0000000000000000ull, 0x0000000000000000ull, 0u, 0x0000000000000001ull, 0x0000000000000000ull},
        {0x0000000000000005ull, 0x0000000000000000ull, 1u, 0x0000000000000005ull, 0x0000000000000000ull},
    };
    for (const auto& row : rows) {
        const uint128_t base(row.baseLow, row.baseHigh);
        const uint128_t result = pow(base, row.exponent);
        EXPECT_TRUE(result == uint128_t(row.resultLow, row.resultHigh))
            << "base=" << (std::string)base << ", exponent=" << row.exponent << ", result=" << (std::string)result;
    }
}
// The log10 counterpart of LogBoundaries: log10() has to answer k at 10^k and k - 1 just below it, for
// every power of ten that fits. The random log10 test compares against a double and cannot get this
// close to a boundary.
TEST(uint128_t, Log10Boundaries)
{
    // 10^k for k = 0 to 38, as {low, high}
    const uint64_t powers[][2] = {
        {0x0000000000000001ull, 0x0000000000000000ull}, {0x000000000000000Aull, 0x0000000000000000ull},
        {0x0000000000000064ull, 0x0000000000000000ull}, {0x00000000000003E8ull, 0x0000000000000000ull},
        {0x0000000000002710ull, 0x0000000000000000ull}, {0x00000000000186A0ull, 0x0000000000000000ull},
        {0x00000000000F4240ull, 0x0000000000000000ull}, {0x0000000000989680ull, 0x0000000000000000ull},
        {0x0000000005F5E100ull, 0x0000000000000000ull}, {0x000000003B9ACA00ull, 0x0000000000000000ull},
        {0x00000002540BE400ull, 0x0000000000000000ull}, {0x000000174876E800ull, 0x0000000000000000ull},
        {0x000000E8D4A51000ull, 0x0000000000000000ull}, {0x000009184E72A000ull, 0x0000000000000000ull},
        {0x00005AF3107A4000ull, 0x0000000000000000ull}, {0x00038D7EA4C68000ull, 0x0000000000000000ull},
        {0x002386F26FC10000ull, 0x0000000000000000ull}, {0x016345785D8A0000ull, 0x0000000000000000ull},
        {0x0DE0B6B3A7640000ull, 0x0000000000000000ull}, {0x8AC7230489E80000ull, 0x0000000000000000ull},
        {0x6BC75E2D63100000ull, 0x0000000000000005ull}, {0x35C9ADC5DEA00000ull, 0x0000000000000036ull},
        {0x19E0C9BAB2400000ull, 0x000000000000021Eull}, {0x02C7E14AF6800000ull, 0x000000000000152Dull},
        {0x1BCECCEDA1000000ull, 0x000000000000D3C2ull}, {0x161401484A000000ull, 0x0000000000084595ull},
        {0xDCC80CD2E4000000ull, 0x000000000052B7D2ull}, {0x9FD0803CE8000000ull, 0x00000000033B2E3Cull},
        {0x3E25026110000000ull, 0x00000000204FCE5Eull}, {0x6D7217CAA0000000ull, 0x00000001431E0FAEull},
        {0x4674EDEA40000000ull, 0x0000000C9F2C9CD0ull}, {0xC0914B2680000000ull, 0x0000007E37BE2022ull},
        {0x85ACEF8100000000ull, 0x000004EE2D6D415Bull}, {0x38C15B0A00000000ull, 0x0000314DC6448D93ull},
        {0x378D8E6400000000ull, 0x0001ED09BEAD87C0ull}, {0x2B878FE800000000ull, 0x0013426172C74D82ull},
        {0xB34B9F1000000000ull, 0x00C097CE7BC90715ull}, {0x00F436A000000000ull, 0x0785EE10D5DA46D9ull},
        {0x098A224000000000ull, 0x4B3B4CA85A86C47Aull},
    };
    static_assert(std::extent_v<decltype(powers)> == 39, "10^0 to 10^38");
    for (uint64_t k = 0; k < std::extent_v<decltype(powers)>; ++k) {
        const uint128_t power(powers[k][0], powers[k][1]);
        EXPECT_EQ(log10(power), k) << "at 10^" << k;
        if (k > 0) {
            EXPECT_EQ(log10(power - uint128_t(1ull)), k - 1) << "just below 10^" << k;
        }
    }
    // everything from 10^38 up to the top of the range answers 38
    EXPECT_EQ(log10(uint128_t(~0ull, ~0ull)), 38ull);
}
// Parsing corners no other test reaches: hex digits across the QWORD boundary (the longest hex string
// elsewhere has eight digits), overflow wrapping past 2^128 in both bases, the upper case prefix, and
// each of the characters that end the number.
TEST(uint128_t, ConstructorFromStringCorners)
{
    const struct {
        const char* text;
        uint64_t low, high;
    } rows[] = {
        {"0x0123456789ABCDEF0123456789ABCDEF", 0x0123456789ABCDEFull, 0x0123456789ABCDEFull},
        {"0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF", 0xFFFFFFFFFFFFFFFFull, 0xFFFFFFFFFFFFFFFFull},
        {"0x10000000000000000", 0x0000000000000000ull, 0x0000000000000001ull},
        {"0X1F", 0x000000000000001Full, 0x0000000000000000ull},
        {"0x1FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF", 0xFFFFFFFFFFFFFFFFull, 0xFFFFFFFFFFFFFFFFull},  // 33 digits, the top one shifts out
        {"340282366920938463463374607431768211456", 0x0000000000000000ull, 0x0000000000000000ull},  // 2^128 wraps to zero
        {"340282366920938463463374607431768211457", 0x0000000000000001ull, 0x0000000000000000ull},
        {"000000000000000000000000000000000000000042", 0x000000000000002Aull, 0x0000000000000000ull},
        {"123.999", 0x000000000000007Bull, 0x0000000000000000ull},
        {"42 17", 0x000000000000002Aull, 0x0000000000000000ull},
        {"\t\n 42", 0x000000000000002Aull, 0x0000000000000000ull},
        {"12abc", 0x000000000000000Cull, 0x0000000000000000ull},
        {"abc", 0x0000000000000000ull, 0x0000000000000000ull},
        {"", 0x0000000000000000ull, 0x0000000000000000ull},
        {"0x", 0x0000000000000000ull, 0x0000000000000000ull},
        {"-340282366920938463463374607431768211455", 0x0000000000000001ull, 0x0000000000000000ull},
    };
    for (const auto& row : rows) {
        const uint128_t expected(row.low, row.high);
        EXPECT_TRUE(uint128_t(row.text) == expected) << "text=\"" << row.text << "\"";
        EXPECT_TRUE(uint128_t(std::string(row.text)) == expected) << "text=\"" << row.text << "\"";
    }
    EXPECT_TRUE(uint128_t(static_cast<const char*>(nullptr)).is_zero());
}
// The random sqrt test only sees values a double can hold, 53 significant bits at most. Perfect squares
// and their neighbours with all 128 bits in play are where an integer Newton iteration goes wrong.
TEST(uint128_t, SqrtOfFullWidthValues)
{
    const auto check = [](uint64_t r) {
        const uint128_t square = ReferenceMultiply(uint128_t(r, 0ull), uint128_t(r, 0ull));
        EXPECT_EQ(sqrt(square), r) << "r=" << r;
        if (r > 0) {
            EXPECT_EQ(sqrt(square - uint128_t(1ull)), r - 1) << "r=" << r;
        }
        // (r + 1)^2 - 1 is the largest value whose root is still r; for r = 2^64 - 1 it is 2^128 - 1
        const uint128_t top = (r == UINT64_MAX) ? uint128_t(~0ull, ~0ull) : ReferenceMultiply(uint128_t(r + 1, 0ull), uint128_t(r + 1, 0ull)) - uint128_t(1ull);
        EXPECT_EQ(sqrt(top), r) << "r=" << r;
    };

    const uint64_t roots[] = {1ull, 2ull, 3ull, 0xFFFFFFFFull, 0x100000000ull, 0x7FFFFFFFFFFFFFFFull, 0x8000000000000000ull, ~0ull};
    for (const uint64_t r : roots) {
        check(r);
    }

    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        check(get_uint64_random() >> (get_uint32_random() % 64u));
        const uint128_t x(get_uint64_random(), get_uint64_random() >> (get_uint32_random() % 64u));
        const uint64_t root = sqrt(x);
        EXPECT_TRUE(is_exact_isqrt(root, x)) << "x=" << (std::string)x << ", sqrt=" << root;
    }
}
// Comparisons between operands that differ in a single QWORD by one, or only in the msb of the low
// QWORD, or not at all. CompareUint128 builds its operands from doubles, which cannot express any of
// these. The reference compares the QWORD pairs natively.
TEST(uint128_t, CompareRawQwords)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const uint64_t al = get_uint64_random(), ah = get_uint64_random();
        uint64_t bl = al, bh = ah;
        switch (i % 6u) {
        case 0:
            break;  // equal
        case 1:
            bl = al + 1ull;
            break;
        case 2:
            bl = al - 1ull;
            break;
        case 3:
            bh = ah + 1ull;
            break;
        case 4:
            bh = ah - 1ull;
            break;
        default:
            bl = al ^ (1ull << 63);  // the low QWORD has to be compared as unsigned
            break;
        }
        const uint128_t a(al, ah), b(bl, bh);
        const bool less = ah < bh || (ah == bh && al < bl);
        const bool equal = al == bl && ah == bh;

        EXPECT_EQ(a < b, less) << "a=" << (std::string)a << ", b=" << (std::string)b;
        EXPECT_EQ(a <= b, less || equal) << "a=" << (std::string)a << ", b=" << (std::string)b;
        EXPECT_EQ(a > b, !less && !equal) << "a=" << (std::string)a << ", b=" << (std::string)b;
        EXPECT_EQ(a >= b, !less) << "a=" << (std::string)a << ", b=" << (std::string)b;
        EXPECT_EQ(a == b, equal) << "a=" << (std::string)a << ", b=" << (std::string)b;
        EXPECT_EQ(a != b, !equal) << "a=" << (std::string)a << ", b=" << (std::string)b;

        // The templated forms widen a uint64_t scalar to 128 bits. Against a narrow operand the low
        // QWORDs decide; against a wide one the high QWORD does.
        const uint128_t narrow(al, 0ull);
        EXPECT_EQ(narrow < bl, al < bl) << "a=" << al << ", b=" << bl;
        EXPECT_EQ(narrow <= bl, al <= bl) << "a=" << al << ", b=" << bl;
        EXPECT_EQ(narrow > bl, al > bl) << "a=" << al << ", b=" << bl;
        EXPECT_EQ(narrow >= bl, al >= bl) << "a=" << al << ", b=" << bl;
        EXPECT_EQ(a <= bl, ah == 0 && al <= bl) << "a=" << (std::string)a << ", b=" << bl;
        EXPECT_EQ(a >= bl, ah != 0 || al >= bl) << "a=" << (std::string)a << ", b=" << bl;
    }

    // a negative scalar converts to the unsigned type first, so -1 is the maximum, as for the builtins
    const uint128_t five(5ull);
    EXPECT_TRUE(five < -1);
    EXPECT_TRUE(five <= -1);
    EXPECT_FALSE(five >= -1);
    EXPECT_TRUE(uint128_t(~0ull, ~0ull) == -1);
}
// The unsigned counterpart of the int128_t test of the same name. None of these ran on the unsigned
// type except inside a static_assert, where the portable constant evaluated code stands in.
TEST(uint128_t, UnaryOperatorsAndAccessors)
{
    const uint128_t x(0x0123456789ABCDEFull, 0xFEDCBA9876543210ull);
    const uint128_t max(~0ull, ~0ull);

    // unary plus is a copy, unary minus a two's complement negation, ~x the complement of each QWORD
    EXPECT_TRUE(+x == x);
    EXPECT_TRUE(-x == uint128_t(~0x0123456789ABCDEFull + 1ull, ~0xFEDCBA9876543210ull));
    EXPECT_TRUE(x + (-x) == uint128_t());
    EXPECT_TRUE(-uint128_t(1ull) == max);
    EXPECT_TRUE(-uint128_t() == uint128_t());
    EXPECT_TRUE(~x == uint128_t(~0x0123456789ABCDEFull, ~0xFEDCBA9876543210ull));
    EXPECT_TRUE(~x == max - x);
    EXPECT_TRUE(~(~x) == x);

    // operator bool and operator!, including a value whose low QWORD is zero
    EXPECT_TRUE(static_cast<bool>(uint128_t(0ull, 1ull)));
    EXPECT_FALSE(!uint128_t(0ull, 1ull));
    EXPECT_FALSE(static_cast<bool>(uint128_t()));
    EXPECT_TRUE(!uint128_t());

    // an unsigned value is never negative, the msb set included
    EXPECT_TRUE(max.is_positive());
    EXPECT_FALSE(max.is_negative());
    EXPECT_TRUE(uint128_t().is_positive());
    EXPECT_TRUE(x.is_int());

    // get_bit over both QWORDs
    EXPECT_EQ(x.get_bit(0), 1);
    EXPECT_EQ(x.get_bit(63), 0);
    EXPECT_EQ(x.get_bit(64), 0);
    EXPECT_EQ(x.get_bit(68), 1);
    EXPECT_EQ(x.get_bit(127), 1);

    EXPECT_TRUE(uint128_t::one() == uint128_t(1ull));
    uint64_t low = 0, high = 0;
    x.get_components(low, high);
    EXPECT_EQ(low, 0x0123456789ABCDEFull);
    EXPECT_EQ(high, 0xFEDCBA9876543210ull);

    // hex() prints both QWORDs, high first, zero padded
    EXPECT_STREQ(x.hex(), "0xFEDCBA98765432100123456789ABCDEF");
    EXPECT_STREQ(uint128_t().hex(), "0x00000000000000000000000000000000");
    EXPECT_STREQ(max.hex(), "0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF");

    // long double goes through the double conversion; 2^128 - 1 rounds to 2^128, which is exact
    EXPECT_EQ(static_cast<long double>(max), 340282366920938463463374607431768211456.0L);
    // zero takes an early return in both floating point conversions
    EXPECT_EQ(std::bit_cast<uint32_t>(static_cast<float>(uint128_t())), 0u);
    EXPECT_EQ(std::bit_cast<uint64_t>(static_cast<double>(uint128_t())), 0ull);
}
// div_64bit() is the multi word division by a 64 bit divisor that every integer type reaches. It was
// only ever called through operators that filter out its trivial cases first, so its own contract
// went untested: the error returns, the shortcut for a single word equal to the divisor, numerators
// wider than two words, and the in place use the uint128_t scalar division depends on. Expected
// values from Python's arbitrary precision integers.
TEST(uint128_t, Div64BitContract)
{
    {
        const uint64_t u[3] = {0x0123456789ABCDEFull, 0xFEDCBA9876543210ull, 0x1ull};
        uint64_t q[3] {}, r = 99;
        EXPECT_TRUE(div_64bit(q, &r, u, 10, 3));
        EXPECT_EQ(q[0], 0x99B6BA23F42AC7CBull);
        EXPECT_EQ(q[1], 0x331612A8D8A20501ull);
        EXPECT_EQ(q[2], 0ull);
        EXPECT_EQ(r, 1ull);
    }
    {
        const uint64_t u[4] = {1ull, 0ull, 0ull, 0x8000000000000000ull};
        uint64_t q[4] {}, r = 99;
        EXPECT_TRUE(div_64bit(q, &r, u, 3, 4));
        EXPECT_EQ(q[0], 0xAAAAAAAAAAAAAAABull);
        EXPECT_EQ(q[1], 0xAAAAAAAAAAAAAAAAull);
        EXPECT_EQ(q[2], 0xAAAAAAAAAAAAAAAAull);
        EXPECT_EQ(q[3], 0x2AAAAAAAAAAAAAAAull);
        EXPECT_EQ(r, 0ull);
    }
    {
        // the quotient written over the numerator, as uint128_t::operator/=(T) does
        uint64_t u[4] = {~0ull, ~0ull, ~0ull, ~0ull};
        EXPECT_TRUE(div_64bit(u, nullptr, u, ~0ull, 4));
        EXPECT_EQ(u[0], 1ull);
        EXPECT_EQ(u[1], 1ull);
        EXPECT_EQ(u[2], 1ull);
        EXPECT_EQ(u[3], 1ull);
    }

    // the trivial cases: a single word equal to, or below, the divisor, and a zero numerator
    const uint64_t u[2] = {42ull, 0ull};
    uint64_t q[2] {}, r = 0;
    EXPECT_TRUE(div_64bit(q, &r, u, 42, 2));
    EXPECT_EQ(q[0], 1ull);
    EXPECT_EQ(r, 0ull);
    q[0] = 0;
    EXPECT_TRUE(div_64bit(q, &r, u, 43, 2));
    EXPECT_EQ(q[0], 0ull);
    EXPECT_EQ(r, 42ull);
    const uint64_t zero[2] {};
    r = 5;
    EXPECT_TRUE(div_64bit(q, &r, zero, 7, 2));
    EXPECT_EQ(r, 0ull);

    // the error returns: a zero divisor and the missing pointers
    EXPECT_FALSE(div_64bit(q, &r, u, 0, 2));
    EXPECT_FALSE(div_64bit(nullptr, &r, u, 3, 2));
    EXPECT_FALSE(div_64bit(q, &r, nullptr, 3, 2));
}
// Whether an operation may throw is part of the interface: generic code branches on it, and a
// container only moves elements whose move constructor is noexcept. Division, modulo and the
// logarithms throw on a zero or out of domain operand and must not claim otherwise.
TEST(uint128_t, NoexceptContracts)
{
    uint128_t a, b;

    static_assert(noexcept(a + b) && noexcept(a - b) && noexcept(a * b));
    static_assert(noexcept(a + 1) && noexcept(1 + a) && noexcept(a * 1u) && noexcept(2 * a) && noexcept(a - 1) && noexcept(1 - a));
    static_assert(noexcept(a += b) && noexcept(a -= b) && noexcept(a *= b) && noexcept(a &= b) && noexcept(a |= b) && noexcept(a ^= b));
    static_assert(noexcept(a & b) && noexcept(a | b) && noexcept(a ^ b) && noexcept(a & 1) && noexcept(a | 1u) && noexcept(a ^ 1ull));
    static_assert(noexcept(a << 3) && noexcept(a >> 3) && noexcept(a <<= 3) && noexcept(a >>= 3) && noexcept(1 << a));
    static_assert(noexcept(a == b) && noexcept(a != b) && noexcept(a < b) && noexcept(a <= b) && noexcept(a > b) && noexcept(a >= b));
    static_assert(noexcept(a == 1) && noexcept(1 == a) && noexcept(a < 1) && noexcept(1 < a));
    static_assert(noexcept(~a) && noexcept(-a) && noexcept(+a) && noexcept(!a));
    static_assert(noexcept(++a) && noexcept(a++) && noexcept(--a) && noexcept(a--));
    static_assert(noexcept(sqr(a)) && noexcept(sqrt(a)) && noexcept(pow(a, 2u)) && noexcept(lzcnt128(a)));
    static_assert(noexcept(static_cast<double>(a)) && noexcept(static_cast<float>(a)) && noexcept(static_cast<uint64_t>(a)) && noexcept(a.hex()));
    static_assert(noexcept(uint128_t(1)) && noexcept(uint128_t(1.5)) && noexcept(uint128_t("1")));
    static_assert(std::is_nothrow_copy_constructible_v<uint128_t> && std::is_nothrow_move_constructible_v<uint128_t>);
    static_assert(std::is_nothrow_copy_assignable_v<uint128_t> && std::is_nothrow_move_assignable_v<uint128_t>);

    // the operations that throw
    static_assert(!noexcept(a / b) && !noexcept(a % b) && !noexcept(a / 2) && !noexcept(a % 2) && !noexcept(a /= b) && !noexcept(a %= b));
    static_assert(!noexcept(log2(a)) && !noexcept(log(a)) && !noexcept(log10(a)));
    SUCCEED();
}
// Comparisons with a floating point scalar compare real numbers. They used to convert the scalar to
// uint128_t first, which truncated it (2 == 2.5 held) and saturated NaN, the infinities and 1e39 to
// the maximum, which then compared equal to all three. Each row gives the order of one value against
// every scalar column, from Python's exact comparison of an int with a float: '<' the value is below
// the scalar, '=' equal, '>' above, 'u' unordered. Integral scalars are still widened instead, see
// CompareRawQwords.
TEST(uint128_t, CompareWithFloatingPointIsExact)
{
    const double inf = std::numeric_limits<double>::infinity();
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double doubles[] = {0.0,     -0.0,    0.5,     1.0,    1.5,    2.0,    2.5,    -0.5,    -1.0,
                              0x1p53,  0x1.0000000000001p53,   0x1.fffffffffffffp63,   0x1p64,  0x1p127,
                              0x1.fffffffffffffp127,           0x1p128, 1e39,   DBL_MAX, DBL_MIN, inf, -inf, nan};
    const float finf = std::numeric_limits<float>::infinity();
    const float fnan = std::numeric_limits<float>::quiet_NaN();
    const float floats[] = {0.0f, -0.0f, 2.5f, -0.5f, 0x1p24f, 0x1.000002p24f, 0x1p64f, 0x1.fffffep127f, 0x1p127f, finf, -finf, fnan};

    const struct {
        uint64_t low, high;
        const char* vsDouble;
        const char* vsFloat;
    } rows[] = {
        {0x0000000000000000ull, 0x0000000000000000ull, "==<<<<<>><<<<<<<<<<<>u", "==<><<<<<<>u"},  // 0
        {0x0000000000000001ull, 0x0000000000000000ull, ">>>=<<<>><<<<<<<<<><>u", ">><><<<<<<>u"},  // 1
        {0x0000000000000002ull, 0x0000000000000000ull, ">>>>>=<>><<<<<<<<<><>u", ">><><<<<<<>u"},  // 2
        {0x0000000000000003ull, 0x0000000000000000ull, ">>>>>>>>><<<<<<<<<><>u", ">>>><<<<<<>u"},  // 3
        {0x0020000000000000ull, 0x0000000000000000ull, ">>>>>>>>>=<<<<<<<<><>u", ">>>>>><<<<>u"},  // 2^53
        {0x0020000000000001ull, 0x0000000000000000ull, ">>>>>>>>>><<<<<<<<><>u", ">>>>>><<<<>u"},  // 2^53+1
        {0xFFFFFFFFFFFFFFFFull, 0x0000000000000000ull, ">>>>>>>>>>>><<<<<<><>u", ">>>>>><<<<>u"},  // 2^64-1
        {0x0000000000000000ull, 0x0000000000000001ull, ">>>>>>>>>>>>=<<<<<><>u", ">>>>>>=<<<>u"},  // 2^64
        {0x0000000000000001ull, 0x0000000000000001ull, ">>>>>>>>>>>>><<<<<><>u", ">>>>>>><<<>u"},  // 2^64+1
        {0x0000000000000000ull, 0x8000000000000000ull, ">>>>>>>>>>>>>=<<<<><>u", ">>>>>>><=<>u"},  // 2^127
        {0x0000000000000000ull, 0xFFFFFFFFFFFFF800ull, ">>>>>>>>>>>>>>=<<<><>u", ">>>>>>>>><>u"},  // (2^53-1)*2^75
        {0x0000000000000001ull, 0xFFFFFFFFFFFFF800ull, ">>>>>>>>>>>>>>><<<><>u", ">>>>>>>>><>u"},  // (2^53-1)*2^75+1
        {0xFFFFFFFFFFFFFFFFull, 0xFFFFFFFFFFFFF7FFull, ">>>>>>>>>>>>>><<<<><>u", ">>>>>>>>><>u"},  // (2^53-1)*2^75-1
        {0x0000000000000000ull, 0xFFFFFF0000000000ull, ">>>>>>>>>>>>>><<<<><>u", ">>>>>>>=><>u"},  // (2^24-1)*2^104
        {0xFFFFFFFFFFFFFFFFull, 0xFFFFFFFFFFFFFFFFull, ">>>>>>>>>>>>>>><<<><>u", ">>>>>>>>><>u"},  // 2^128-1
    };

    // all six operators, with the scalar on either side, against the expected relation
    const auto check = [](const uint128_t& v, auto scalar, char relation) {
        const bool less = relation == '<', equal = relation == '=', greater = relation == '>';
        EXPECT_EQ(v == scalar, equal) << "v=" << (std::string)v << ", scalar=" << std::hexfloat << scalar;
        EXPECT_EQ(v != scalar, !equal) << "v=" << (std::string)v << ", scalar=" << std::hexfloat << scalar;
        EXPECT_EQ(v < scalar, less) << "v=" << (std::string)v << ", scalar=" << std::hexfloat << scalar;
        EXPECT_EQ(v <= scalar, less || equal) << "v=" << (std::string)v << ", scalar=" << std::hexfloat << scalar;
        EXPECT_EQ(v > scalar, greater) << "v=" << (std::string)v << ", scalar=" << std::hexfloat << scalar;
        EXPECT_EQ(v >= scalar, greater || equal) << "v=" << (std::string)v << ", scalar=" << std::hexfloat << scalar;
        EXPECT_EQ(scalar == v, equal) << "v=" << (std::string)v << ", scalar=" << std::hexfloat << scalar;
        EXPECT_EQ(scalar != v, !equal) << "v=" << (std::string)v << ", scalar=" << std::hexfloat << scalar;
        EXPECT_EQ(scalar < v, greater) << "v=" << (std::string)v << ", scalar=" << std::hexfloat << scalar;
        EXPECT_EQ(scalar <= v, greater || equal) << "v=" << (std::string)v << ", scalar=" << std::hexfloat << scalar;
        EXPECT_EQ(scalar > v, less) << "v=" << (std::string)v << ", scalar=" << std::hexfloat << scalar;
        EXPECT_EQ(scalar >= v, less || equal) << "v=" << (std::string)v << ", scalar=" << std::hexfloat << scalar;
    };

    for (const auto& row : rows) {
        const uint128_t v(row.low, row.high);
        ASSERT_EQ(strlen(row.vsDouble), std::extent_v<decltype(doubles)>);
        ASSERT_EQ(strlen(row.vsFloat), std::extent_v<decltype(floats)>);
        for (size_t j = 0; j < std::extent_v<decltype(doubles)>; ++j) {
            check(v, doubles[j], row.vsDouble[j]);
        }
        for (size_t j = 0; j < std::extent_v<decltype(floats)>; ++j) {
            check(v, floats[j], row.vsFloat[j]);
        }
    }
    // a long double goes through double, which on every supported compiler holds these exactly
    check(uint128_t(2ull), 2.5L, '<');
    check(uint128_t(3ull), 2.5L, '>');

    // the same in constant evaluation
    static_assert(uint128_t(2ull) != 2.5 && uint128_t(2ull) < 2.5 && 2.5 > uint128_t(2ull) && uint128_t(3ull) > 2.5f);
    static_assert(!(uint128_t() == std::numeric_limits<double>::quiet_NaN()) && uint128_t() != std::numeric_limits<double>::quiet_NaN());
    static_assert(!(uint128_t() < std::numeric_limits<double>::quiet_NaN()) && !(std::numeric_limits<double>::quiet_NaN() >= uint128_t()));
    static_assert(uint128_t(~0ull, ~0ull) < std::numeric_limits<double>::infinity() && uint128_t(~0ull, ~0ull) < 1e39);
    static_assert(uint128_t() > -0.5 && uint128_t() > -std::numeric_limits<double>::infinity() && uint128_t() == -0.0);
}
// operator*=(T) used to accept any operand that converted to uint64_t, so the other signedness slipped
// through a conversion operator and was truncated to its low QWORD. Mixing the two now fails to compile,
// as it already did for the other arithmetic operators; the builtin operands keep working.
TEST(uint128_t, MultiplyRejectsTheOtherSignedness)
{
    static_assert(!Multipliable<uint128_t, int128_t>, "uint128_t * int128_t must not compile");
    static_assert(!MultiplyAssignable<uint128_t, int128_t>, "uint128_t *= int128_t must not compile");
    static_assert(Multipliable<uint128_t, uint128_t> && Multipliable<uint128_t, int> && Multipliable<uint128_t, uint64_t>);
    static_assert(Multipliable<uint128_t, double> && Multipliable<int, uint128_t> && MultiplyAssignable<uint128_t, int64_t>);
    SUCCEED();
}
// shift_right128() and shift_left128() take a count from 0 to 127. At 0 the bits crossing from the
// other QWORD moved by 64, which is undefined behavior and on x64 ORed the whole other QWORD in. Each
// is compared with the matching QWORD of the bit by bit reference, at every count, with the count
// hidden from the optimizer, and once in constant evaluation at the count that used to fail.
TEST(uint128_t, Shift128HelpersAcceptAZeroCount)
{
    static_assert(shift_right128(5ull, 0x7000000000000000ull, 0) == 5ull);
    static_assert(shift_left128(5ull, 0x7000000000000000ull, 0) == 0x7000000000000000ull);

    const uint128_t values[] = {uint128_t(5ull, 0x7000000000000000ull), uint128_t(0x0123456789ABCDEFull, 0xFEDCBA9876543210ull), uint128_t(~0ull, ~0ull)};
    for (const auto& x : values) {
        uint64_t low = 0, high = 0;
        x.get_components(low, high);
        for (int32_t s = 0; s < 128; ++s) {
            volatile int32_t count = s;
            uint64_t expectedLow = 0, unused = 0, expectedHigh = 0;
            ReferenceShift(x, s, false).get_components(expectedLow, unused);
            ReferenceShift(x, s, true).get_components(unused, expectedHigh);
            EXPECT_EQ(shift_right128(low, high, count), expectedLow) << "x=" << (std::string)x << ", shift=" << s;
            EXPECT_EQ(shift_left128(low, high, count), expectedHigh) << "x=" << (std::string)x << ", shift=" << s;
        }
    }
}

/**********************************************************************
 * Compile time (constexpr) evaluation
 *
 * Each test below pairs static_asserts, which fail the build the moment one of these operations
 * stops being usable in a constant expression, with a runtime check of the same expression built
 * from opaque() values the optimizer cannot fold. See the matching block in
 * fixed_point128_gtest.cpp for why the runtime half is not a restatement of the compile time one.
 ***********************************************************************/
TEST(uint128_t, ConstexprConstructionAndConversion)
{
    constexpr uint128_t zero;
    constexpr uint128_t seven(7);
    constexpr uint128_t copied(seven);
    constexpr uint128_t assigned = [] { uint128_t a; a = 5; return a; }();
    constexpr uint128_t from_bool(true);
    constexpr uint128_t from_negative(-5);  // wraps around via 2's complement, like the builtin types

    static_assert(zero.is_zero());
    static_assert(static_cast<uint64_t>(seven) == 7ull);
    static_assert(static_cast<uint32_t>(assigned) == 5u);
    static_assert(static_cast<uint64_t>(from_bool) == 1ull);
    static_assert(static_cast<uint64_t>(uint128_t(1ull, 2ull)) == 1ull);  // the high QWORD is dropped
    static_assert(copied == seven);
    static_assert(from_negative == uint128_t(0ull) - uint128_t(5ull));
    static_assert(seven.is_positive() && !seven.is_negative());  // always so for the unsigned type

    static_assert(seven.get_bit(0) == 1 && seven.get_bit(3) == 0);
    static_assert(uint128_t(0ull, 1ull).get_bit(64) == 1);
    static_assert(uint128_t(1ull, 2ull).get_bit(65) == 1);      // a bit of the high QWORD
    static_assert((~uint128_t(3ull, 0ull)).get_bit(0) == 0);
    static_assert(uint128_t(-1).get_bit(127) == 1);             // the integral constructor sign extends
    static_assert(static_cast<uint64_t>(uint128_t::one()) == 1ull);
    // constructed inside the lambda: a captureless lambda cannot odr-use a local, and passing an
    // object to get_components() by reference is exactly that
    constexpr uint64_t high_qword = [] {
        uint64_t l = 0, h = 0;
        uint128_t(0x00000000FFFFFFFFull, 0x1234ull).get_components(l, h);
        return h;
    }();
    static_assert(high_qword == 0x1234ull);

    EXPECT_TRUE(uint128_t(opaque(7ull)) == seven);
    EXPECT_TRUE(uint128_t(opaque(0ull)) - uint128_t(5ull) == from_negative);
}
TEST(uint128_t, ConstexprShiftsBitwiseAndUnary)
{
    constexpr uint128_t all_ones(~0ull, ~0ull);

    static_assert(static_cast<uint64_t>(uint128_t(1ull) << 3) == 8ull);
    static_assert(static_cast<uint64_t>(uint128_t(64ull) >> 3) == 8ull);
    static_assert((uint128_t(1ull) << 64) == uint128_t(0ull, 1ull));  // crosses the QWORD boundary
    static_assert((uint128_t(0ull, 1ull) >> 1) == uint128_t(1ull << 63, 0ull));
    static_assert((all_ones >> 64) == uint128_t(~0ull));
    static_assert((uint128_t(1ull) << 128).is_zero());  // shifting the whole value out
    static_assert((all_ones >> 128).is_zero());
    static_assert((uint128_t(1ull) << 0) == uint128_t(1ull));  // a zero count leaves the value alone
    constexpr uint128_t roundtrip = [] { uint128_t a(1ull); a <<= 100; a >>= 100; return a; }();
    static_assert(static_cast<uint64_t>(roundtrip) == 1ull);

    static_assert(static_cast<uint64_t>(uint128_t(12ull) & uint128_t(10ull)) == 8ull);
    static_assert(static_cast<uint64_t>(uint128_t(12ull) | uint128_t(3ull)) == 15ull);
    static_assert(static_cast<uint64_t>(uint128_t(12ull) ^ uint128_t(10ull)) == 6ull);
    static_assert(static_cast<uint64_t>(12ull & uint128_t(10ull)) == 8ull);  // scalar on the left
    static_assert(~uint128_t(0ull) == all_ones);
    static_assert(!uint128_t(0ull) && static_cast<bool>(uint128_t(1ull)));
    static_assert(+uint128_t(5ull) == uint128_t(5ull));

    EXPECT_TRUE(((uint128_t(opaque(1ull)) << 100) >> 100) == uint128_t(1ull));
    EXPECT_TRUE((uint128_t(0ull, opaque(1ull)) >> 1) == uint128_t(1ull << 63, 0ull));
    EXPECT_TRUE((uint128_t(opaque(~0ull), opaque(~0ull)) >> 64) == uint128_t(~0ull));
}
TEST(uint128_t, ConstexprArithmetic)
{
    constexpr uint128_t all_ones(~0ull, ~0ull);

    static_assert(static_cast<uint64_t>(uint128_t(3ull) + uint128_t(4ull)) == 7ull);
    static_assert(uint128_t(~0ull) + uint128_t(1ull) == uint128_t(0ull, 1ull));  // carry into high
    static_assert(uint128_t(0ull) - uint128_t(1ull) == all_ones);                // underflow wraps
    static_assert(all_ones + uint128_t(1ull) == uint128_t(0ull));                // overflow wraps
    static_assert(static_cast<uint64_t>(uint128_t(3ull) + 4) == 7ull);  // generic right hand side
    static_assert(static_cast<uint64_t>(4 + uint128_t(3ull)) == 7ull);  // scalar on the left
    constexpr uint128_t summed = [] { uint128_t a; for (int i = 1; i <= 10; ++i) a += i; return a; }();
    static_assert(static_cast<uint64_t>(summed) == 55ull);

    static_assert(static_cast<uint64_t>(uint128_t(6ull) * uint128_t(7ull)) == 42ull);
    static_assert(uint128_t(1ull << 32) * uint128_t(1ull << 32) == uint128_t(0ull, 1ull));
    static_assert(uint128_t(~0ull) * uint128_t(~0ull) == uint128_t(1ull, 0xFFFFFFFFFFFFFFFEull));
    static_assert(static_cast<uint64_t>(uint128_t(6ull) * 7) == 42ull);  // generic right hand side
    static_assert(static_cast<uint64_t>(7 * uint128_t(6ull)) == 42ull);  // scalar on the left
    // sqr is documented to be bit identical to x * x
    static_assert(sqr(uint128_t(~0ull)) == uint128_t(~0ull) * uint128_t(~0ull));
    constexpr uint128_t squared = [] { uint128_t a(5ull); a.square(); return a; }();
    static_assert(static_cast<uint64_t>(squared) == 25ull);

    constexpr uint128_t stepped = [] { uint128_t a(5ull); ++a; a++; --a; return a; }();
    static_assert(static_cast<uint64_t>(stepped) == 6ull);
    constexpr uint128_t decremented_zero = [] { uint128_t a(0ull); --a; return a; }();
    static_assert(decremented_zero == all_ones);

    static_assert(uint128_t(1ull) < uint128_t(2ull) && uint128_t(2ull) > uint128_t(1ull));
    static_assert(uint128_t(0ull, 1ull) > uint128_t(~0ull));  // compares across the QWORD boundary
    static_assert(all_ones > uint128_t(0ull));                // the all ones pattern is the largest
    static_assert(uint128_t(1ull) <= 1 && uint128_t(1ull) >= 1);
    static_assert(uint128_t(1ull) == 1 && uint128_t(1ull) != 2);

    EXPECT_TRUE(uint128_t(opaque(6ull)) * uint128_t(opaque(7ull)) == uint128_t(42ull));
    EXPECT_TRUE(uint128_t(opaque(~0ull)) * uint128_t(opaque(~0ull)) == uint128_t(1ull, 0xFFFFFFFFFFFFFFFEull));
    EXPECT_TRUE(uint128_t(opaque(~0ull)) + uint128_t(1ull) == uint128_t(0ull, 1ull));
    EXPECT_TRUE(sqr(uint128_t(opaque(~0ull))) == sqr(uint128_t(~0ull)));
}
TEST(uint128_t, ConstexprMathFunctions)
{
    static_assert(lzcnt128(uint128_t(1ull)) == 127);
    static_assert(lzcnt128(uint128_t(0ull, 1ull)) == 63);
    static_assert(lzcnt128(uint128_t(~0ull, ~0ull)) == 0);
    static_assert(log2(uint128_t(8ull)) == 3);
    static_assert(log2(uint128_t(0ull, 1ull)) == 64);  // 2^64
    static_assert(log2(uint128_t(~0ull, ~0ull)) == 127);
    static_assert(pow(uint128_t(2ull), 0u) == uint128_t(1ull));
    static_assert(pow(uint128_t(0ull), 0u) == uint128_t(1ull));  // matches pow(double, double)
    static_assert(pow(uint128_t(2ull), 10u) == uint128_t(1024ull));
    static_assert(pow(uint128_t(2ull), 100u) == uint128_t(0ull, 1ull << 36));
    // 10^38, the largest power of ten that fits. Spelled as its two QWORDs because the string
    // constructor allocates and is therefore not available in a constant expression.
    static_assert(pow(uint128_t(10ull), 38u) == uint128_t(0x098A224000000000ull, 0x4B3B4CA85A86C47Aull));

    EXPECT_TRUE(lzcnt128(uint128_t(0ull, opaque(1ull))) == 63);
    EXPECT_TRUE(log2(uint128_t(0ull, opaque(1ull))) == 64);
    EXPECT_TRUE(pow(uint128_t(opaque(2ull)), 100u) == pow(uint128_t(2ull), 100u));
}

// The double and float conversions used to be the one thing that could not happen at compile time:
// they read the inactive member of a union, which constant evaluation rejects. Now that Double and
// Float convert with std::bit_cast, a floating point literal crosses the boundary either way.
TEST(uint128_t, ConstexprFloatingPointConversion)
{
    constexpr uint128_t from_double = 42.0;
    constexpr uint128_t big = 1e18;

    static_assert(static_cast<uint64_t>(from_double) == 42ull);
    static_assert(static_cast<double>(from_double) == 42.0);
    static_assert(static_cast<float>(from_double) == 42.0f);
    static_assert(static_cast<double>(big) == 1e18);
    static_assert(uint128_t(0.0).is_zero());
    static_assert(uint128_t(3.99) == uint128_t(3ull));  // truncates towards zero
    // a negative double wraps around via 2's complement, like the integer constructors
    static_assert(uint128_t(-5.0) == uint128_t(0ull) - uint128_t(5ull));

    EXPECT_TRUE(uint128_t(static_cast<double>(opaque(42ull))) == from_double);
    EXPECT_DOUBLE_EQ(static_cast<double>(uint128_t(opaque(42ull))), 42.0);
}
