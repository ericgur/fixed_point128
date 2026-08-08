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
#include <ostream>
#include <ctime>
#include <cfloat>
#include <cstdlib>
#include "gtest_shared.h"

/**********************************************************************
 * int128_t tests
 ***********************************************************************/
// Construct fixed_point128 and convert back to/from various elements.
TEST(int128_t, DefaultConstructor)
{
    int128_t i;
    EXPECT_EQ(static_cast<int64_t>(i), 0ull);
}
TEST(int128_t, ConstructorFromDouble)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = floor(get_double_random(1, 127));
        int128_t f = value;
        if (check_overflow_int128(value)) {
            continue;
        }
        double f_value = static_cast<double>(f);
        EXPECT_DOUBLE_EQ(f_value, value) << "value=" << value;
    }
}
TEST(int128_t, ConstructorFromFloat)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        float value = (float)floor(get_double_random(1, 127));
        int128_t f = value;
        if (check_overflow_int128(value)) {
            continue;
        }
        EXPECT_FLOAT_EQ(static_cast<float>(f), value) << "value=" << value;
    }
}
TEST(int128_t, ConstructorFromInt32)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        int32_t value = get_int32_random();
        int128_t f = value;
        EXPECT_EQ(static_cast<int32_t>(f), value);
    }
}
TEST(int128_t, ConstructorFromUnsignedInt32)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        uint32_t value = get_uint32_random();
        int128_t f = value;
        EXPECT_EQ(static_cast<uint32_t>(f), value);
    }
}
TEST(int128_t, ConstructorFromInt64)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        int64_t value = get_int64_random();
        int128_t f = value;
        EXPECT_EQ(static_cast<int64_t>(f), value);
    }
}
// The previous version used get_int64_random() and int64_t throughout, making it a verbatim copy
// of ConstructorFromInt64 that never once fed the constructor an unsigned operand.
TEST(int128_t, ConstructorFromUnsignedInt64)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const uint64_t value = get_uint64_random();
        const int128_t f = value;

        EXPECT_EQ(static_cast<uint64_t>(f), value);
        // an unsigned operand is zero extended rather than sign extended, so a value with its msb
        // set becomes a positive 128 bit number rather than a negative one
        EXPECT_TRUE(f == int128_t(value, 0ull)) << "value=" << value;
        EXPECT_TRUE(f.is_positive()) << "value=" << value;
    }
}
// The previous version parsed the input, printed it back, and compared the two through strtod. That
// is a round trip, so a parse and a print that were wrong in opposite directions cancelled out, and
// the double funnel discarded everything below 53 bits: "123456789012345678" was compared as
// 123456789012345680. Values are now checked against what they are supposed to be.
TEST(int128_t, ConstructorFromString)
{
    EXPECT_TRUE(int128_t("0xDEADBEAF") == int128_t(0xDEADBEAFll));
    EXPECT_TRUE(int128_t("1234") == int128_t(1234ll));
    EXPECT_TRUE(int128_t("-1234") == int128_t(-1234ll));
    EXPECT_TRUE(int128_t("123456789012345678") == int128_t(123456789012345678ll));
    EXPECT_TRUE(int128_t("-123456789012345678") == int128_t(-123456789012345678ll));
    // values needing more than 64 bit, where a double reference has no chance
    EXPECT_TRUE(int128_t("18446744073709551616") == int128_t(0ull, 1ull));                   // 2^64
    EXPECT_TRUE(int128_t("170141183460469231731687303715884105727") == int128_t(~0ull, 0x7FFFFFFFFFFFFFFFull));
    EXPECT_TRUE(int128_t("-170141183460469231731687303715884105728") == int128_t(0ull, 0x8000000000000000ull));
    EXPECT_TRUE(int128_t("0x7FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF") == int128_t(~0ull, 0x7FFFFFFFFFFFFFFFull));
    // a nullptr, an empty string and a leading illegal character all produce zero
    EXPECT_TRUE(int128_t(static_cast<const char*>(nullptr)).is_zero());
    EXPECT_TRUE(int128_t("").is_zero());
    EXPECT_TRUE(int128_t("abc").is_zero());
    // parsing stops at the first character that is not a digit of the detected base
    EXPECT_TRUE(int128_t("12abc") == int128_t(12ll));
    EXPECT_TRUE(int128_t("  -42  ") == int128_t(-42ll));

    // every value round trips through its own decimal string, both signs, full width
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const int128_t v(get_uint64_random(), get_uint64_random());
        const std::string s = (std::string)v;
        EXPECT_TRUE(int128_t(s) == v) << "s=" << s;
    }
}
TEST(int128_t, CopyConstructor)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random(1, 127);
        int128_t f1 = value;
        int128_t f2 = f1;
        EXPECT_DOUBLE_EQ(static_cast<double>(f1), static_cast<double>(f2));
    }
}
TEST(int128_t, MoveConstructor)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random(1, 127);
        int128_t f1 = value;
        int128_t f2(std::move(f1));
        EXPECT_DOUBLE_EQ(static_cast<double>(f1), static_cast<double>(f2));
    }
}
TEST(int128_t, AssignmentOperator)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random(1, 127);
        int128_t f1 = value;
        int128_t f2;
        f2 = f1;
        EXPECT_DOUBLE_EQ(static_cast<double>(f1), static_cast<double>(f2));
    }
}
TEST(int128_t, MoveAssignmentOperator)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random(1, 127);
        int128_t f1 = value;
        int128_t f2;
        f2 = std::move(f1);
        EXPECT_DOUBLE_EQ(static_cast<double>(f1), static_cast<double>(f2));
    }
}
// Operands are raw QWORD pairs rather than values converted from a double. A double carries 53
// bits, so the previous version could only ever see the top 53 bits of a 128 bit result, and it
// compared them through is_similar_double(), whose 1e-10 tolerance every iteration satisfied. The
// EXPECT below it never executed once in 65536 iterations.
//
// Addition is bit identical for both instantiations: a two's complement sum is a sum modulo 2^128,
// so the reference is the same unsigned limb arithmetic the uint128_t suite uses.
TEST(int128_t, Add)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const uint64_t al = get_uint64_random(), ah = get_uint64_random();
        const uint64_t bl = get_uint64_random(), bh = get_uint64_random();
        const int128_t a(al, ah), b(bl, bh);

        // the reference carry out of the low QWORD: an unsigned sum that wrapped is smaller than
        // either operand. Overflow out of the high QWORD is dropped, as the operator does.
        const uint64_t low = al + bl;
        const uint64_t high = ah + bh + ((low < al) ? 1ull : 0ull);
        const int128_t expected(low, high);

        EXPECT_TRUE(a + b == expected) << "a=" << (std::string)a << ", b=" << (std::string)b;

        int128_t c = a;
        c += b;
        EXPECT_TRUE(c == expected) << "a=" << (std::string)a << ", b=" << (std::string)b;
        int128_t d = a;
        d += d;
        EXPECT_TRUE(d == a + a) << "a=" << (std::string)a;

        EXPECT_TRUE(a + b == b + a) << "a=" << (std::string)a << ", b=" << (std::string)b;
        EXPECT_TRUE((a + b) - b == a) << "a=" << (std::string)a << ", b=" << (std::string)b;
        EXPECT_TRUE(a + int128_t() == a) << "a=" << (std::string)a;
        EXPECT_TRUE(a + (-a) == int128_t()) << "a=" << (std::string)a;
    }
}
TEST(int128_t, AddEdgeCases)
{
    const int128_t zero, one(1ll), minus_one(-1ll);
    const int128_t max(~0ull, 0x7FFFFFFFFFFFFFFFull);   // 2^127 - 1
    const int128_t min(0ull, 0x8000000000000000ull);    // -2^127

    // the carry has to cross the QWORD boundary in both directions
    EXPECT_TRUE(int128_t(~0ull, 0ull) + one == int128_t(0ull, 1ull));
    EXPECT_TRUE(minus_one + one == zero);
    EXPECT_TRUE(zero + minus_one == minus_one);
    // signed overflow wraps silently, matching the builtin types
    EXPECT_TRUE(max + one == min);
    EXPECT_TRUE(min + minus_one == max);
    EXPECT_TRUE(max + max == int128_t(~0ull - 1ull, 0xFFFFFFFFFFFFFFFFull));  // -2
    // adding the most negative value to itself wraps to zero
    EXPECT_TRUE(min + min == zero);
}
// The reference int64_t res = value1 + value2 overflows on a quarter of the operand pairs, and
// -get_int64_random() is itself undefined when the draw comes back INT64_MIN. The operands are now
// full width 128 bit values of opposite sign, with magnitudes kept below 2^127 so that negating
// them is exact.
TEST(int128_t, AddDifferentSign)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const int128_t positive(get_uint64_random(), get_uint64_random() >> 1);
        const int128_t magnitude(get_uint64_random(), get_uint64_random() >> 1);
        const int128_t negative = -magnitude;

        uint64_t pl = 0, ph = 0, nl = 0, nh = 0;
        positive.get_components(pl, ph);
        negative.get_components(nl, nh);
        const uint64_t low = pl + nl;
        const uint64_t high = ph + nh + ((low < pl) ? 1ull : 0ull);

        const int128_t sum = positive + negative;
        EXPECT_TRUE(sum == int128_t(low, high)) << "p=" << (std::string)positive << ", n=" << (std::string)negative;
        // adding a negative value is subtracting its magnitude
        EXPECT_TRUE(sum == positive - magnitude) << "p=" << (std::string)positive << ", n=" << (std::string)negative;
        // and the result is negative exactly when the negative operand has the larger magnitude
        EXPECT_EQ(sum.is_negative(), magnitude > positive) << "p=" << (std::string)positive << ", n=" << (std::string)negative;
    }
}
// The reference used to be computed as value1 + value2 in int64_t, which overflows on roughly a
// quarter of the operand pairs (16281 of 65536 measured) and is undefined behaviour rather than the
// wrap the test was counting on. Both operands are sign extended to 128 bit before being added, so
// the exact result needs no wrapping at all and the full width is worth checking.
TEST(int128_t, AddInt64)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const int64_t value1 = get_int64_random();
        const int64_t value2 = get_int64_random();

        const uint64_t u1 = static_cast<uint64_t>(value1), u2 = static_cast<uint64_t>(value2);
        const uint64_t low = u1 + u2;
        const uint64_t high = SignExtension(value1) + SignExtension(value2) + ((low < u1) ? 1ull : 0ull);
        const int128_t expected(low, high);

        const int128_t f1 = value1;
        EXPECT_TRUE(f1 + value2 == expected) << "value1=" << value1 << ", value2=" << value2;
        EXPECT_TRUE(f1 + value2 == f1 + int128_t(value2)) << "value1=" << value1 << ", value2=" << value2;

        int128_t f2 = f1;
        f2 += value2;
        EXPECT_TRUE(f2 == expected) << "value1=" << value1 << ", value2=" << value2;

        // The sum of two int64_t always fits in 128 bit, so no wrap can occur and the sign of the
        // result is the mathematical one. Derived without overflowing: only a mixed sign pair can
        // change sign, and then the larger magnitude decides. Magnitudes are taken in unsigned
        // arithmetic, where negating INT64_MIN is well defined.
        bool expected_negative;
        if (value1 < 0 && value2 < 0) {
            expected_negative = true;
        } else if (value1 >= 0 && value2 >= 0) {
            expected_negative = false;
        } else {
            const uint64_t m1 = (value1 < 0) ? (0ull - u1) : u1;
            const uint64_t m2 = (value2 < 0) ? (0ull - u2) : u2;
            expected_negative = (value1 < 0) ? (m1 > m2) : (m2 > m1);
        }
        EXPECT_EQ(expected.is_negative(), expected_negative) << "value1=" << value1 << ", value2=" << value2;
    }
}
// The previous version used get_int64_random() and int64_t throughout, so it was a verbatim copy of
// AddInt64 and never fed the operator an unsigned scalar. An unsigned operand is zero extended, so
// a value with its msb set adds a large positive number rather than a small negative one.
TEST(int128_t, AddUnsignedInt64)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const int64_t value1 = get_int64_random();
        const uint64_t value2 = get_uint64_random();

        const uint64_t u1 = static_cast<uint64_t>(value1);
        const uint64_t low = u1 + value2;
        const uint64_t high = SignExtension(value1) + ((low < u1) ? 1ull : 0ull);  // value2 zero extends
        const int128_t expected(low, high);

        const int128_t f1 = value1;
        EXPECT_TRUE(f1 + value2 == expected) << "value1=" << value1 << ", value2=" << value2;
        EXPECT_TRUE(f1 + value2 == f1 + int128_t(value2)) << "value1=" << value1 << ", value2=" << value2;

        int128_t f2 = f1;
        f2 += value2;
        EXPECT_TRUE(f2 == expected) << "value1=" << value1 << ", value2=" << value2;
    }
}

// The previous version derived the subtrahend as floor(value1 / 2.5), so the two operands always
// shared a sign and the result never crossed zero. See Add for the rest of the reasoning.
TEST(int128_t, Subtract)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const uint64_t al = get_uint64_random(), ah = get_uint64_random();
        const uint64_t bl = get_uint64_random(), bh = get_uint64_random();
        const int128_t a(al, ah), b(bl, bh);

        const uint64_t low = al - bl;
        const uint64_t high = ah - bh - ((al < bl) ? 1ull : 0ull);
        const int128_t expected(low, high);

        EXPECT_TRUE(a - b == expected) << "a=" << (std::string)a << ", b=" << (std::string)b;

        int128_t c = a;
        c -= b;
        EXPECT_TRUE(c == expected) << "a=" << (std::string)a << ", b=" << (std::string)b;
        int128_t d = a;
        d -= d;
        EXPECT_TRUE(d.is_zero()) << "a=" << (std::string)a;

        EXPECT_TRUE((a - b) + b == a) << "a=" << (std::string)a << ", b=" << (std::string)b;
        EXPECT_TRUE(a - b == -(b - a)) << "a=" << (std::string)a << ", b=" << (std::string)b;
        EXPECT_TRUE(a - int128_t() == a) << "a=" << (std::string)a;
        EXPECT_TRUE(int128_t() - a == -a) << "a=" << (std::string)a;
        // subtracting is adding the negation, which is what the implementation actually does
        EXPECT_TRUE(a - b == a + (-b)) << "a=" << (std::string)a << ", b=" << (std::string)b;
    }
}
TEST(int128_t, SubtractEdgeCases)
{
    const int128_t zero, one(1ll), minus_one(-1ll);
    const int128_t max(~0ull, 0x7FFFFFFFFFFFFFFFull);   // 2^127 - 1
    const int128_t min(0ull, 0x8000000000000000ull);    // -2^127

    EXPECT_TRUE(zero - one == minus_one);
    EXPECT_TRUE(zero - minus_one == one);
    // the borrow has to cross the QWORD boundary
    EXPECT_TRUE(int128_t(0ull, 1ull) - one == int128_t(~0ull, 0ull));
    // signed underflow wraps silently
    EXPECT_TRUE(min - one == max);
    EXPECT_TRUE(max - minus_one == min);
    // the most negative value is its own negation, so subtracting it is the same as adding it
    EXPECT_TRUE(zero - min == min);
    EXPECT_TRUE(min - min == zero);
}
// Same undefined signed overflow as AddDifferentSign, plus the undefined negation of INT64_MIN.
TEST(int128_t, SubtractDifferentSign)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const int128_t positive(get_uint64_random(), get_uint64_random() >> 1);
        const int128_t magnitude(get_uint64_random(), get_uint64_random() >> 1);
        const int128_t negative = -magnitude;

        uint64_t pl = 0, ph = 0, nl = 0, nh = 0;
        positive.get_components(pl, ph);
        negative.get_components(nl, nh);
        const uint64_t low = pl - nl;
        const uint64_t high = ph - nh - ((pl < nl) ? 1ull : 0ull);

        const int128_t difference = positive - negative;
        EXPECT_TRUE(difference == int128_t(low, high)) << "p=" << (std::string)positive << ", n=" << (std::string)negative;
        // subtracting a negative value adds its magnitude, so the result cannot be smaller
        EXPECT_TRUE(difference == positive + magnitude) << "p=" << (std::string)positive << ", n=" << (std::string)negative;
    }
}
// The reference value1 - value2 in int64_t overflows on roughly a quarter of the operand pairs
// (16427 of 65536 measured), which is undefined behaviour. See AddInt64.
TEST(int128_t, SubtractInt64)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const int64_t value1 = get_int64_random();
        const int64_t value2 = get_int64_random();

        const uint64_t u1 = static_cast<uint64_t>(value1), u2 = static_cast<uint64_t>(value2);
        const uint64_t low = u1 - u2;
        const uint64_t high = SignExtension(value1) - SignExtension(value2) - ((u1 < u2) ? 1ull : 0ull);
        const int128_t expected(low, high);

        const int128_t f1 = value1;
        EXPECT_TRUE(f1 - value2 == expected) << "value1=" << value1 << ", value2=" << value2;
        EXPECT_TRUE(f1 - value2 == f1 - int128_t(value2)) << "value1=" << value1 << ", value2=" << value2;

        int128_t f2 = f1;
        f2 -= value2;
        EXPECT_TRUE(f2 == expected) << "value1=" << value1 << ", value2=" << value2;
    }
}
// The previous version used get_int64_random() and int64_t throughout, so it was a verbatim copy of
// SubtractInt64 and never fed the operator an unsigned scalar.
TEST(int128_t, SubtractUnsignedInt64)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const int64_t value1 = get_int64_random();
        const uint64_t value2 = get_uint64_random();

        const uint64_t u1 = static_cast<uint64_t>(value1);
        const uint64_t low = u1 - value2;
        const uint64_t high = SignExtension(value1) - ((u1 < value2) ? 1ull : 0ull);  // value2 zero extends
        const int128_t expected(low, high);

        const int128_t f1 = value1;
        EXPECT_TRUE(f1 - value2 == expected) << "value1=" << value1 << ", value2=" << value2;
        EXPECT_TRUE(f1 - value2 == f1 - int128_t(value2)) << "value1=" << value1 << ", value2=" << value2;

        int128_t f2 = f1;
        f2 -= value2;
        EXPECT_TRUE(f2 == expected) << "value1=" << value1 << ", value2=" << value2;
    }
}
// Checked against ReferenceMultiply() rather than a double product. Besides losing everything below
// the top 53 bits, the old oracle had a guard that did not work: log2(res) is NaN for a negative
// product and NaN > 127 is false, so half the iterations (32768 of 65536 measured) sailed past the
// overflow check it was supposed to enforce. Every iteration was then skipped by is_similar_double()
// anyway, so the assertion never ran.
TEST(int128_t, MultiplyByint128)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        // shifting the high QWORDs by independent random amounts mixes products that overflow with
        // products that stay in range, and both operands take either sign
        int128_t a(get_uint64_random(), get_uint64_random() >> (get_uint32_random() % 64));
        int128_t b(get_uint64_random(), get_uint64_random() >> (get_uint32_random() % 64));
        if (get_uint32_random() & 1)
            a = -a;
        if (get_uint32_random() & 1)
            b = -b;

        const int128_t expected = ReferenceMultiply(a, b);
        EXPECT_TRUE(a * b == expected) << "a=" << (std::string)a << ", b=" << (std::string)b;

        int128_t c = a;
        c *= b;
        EXPECT_TRUE(c == expected) << "a=" << (std::string)a << ", b=" << (std::string)b;

        // identities that hold in modulo 2^128 arithmetic, so they survive truncation
        EXPECT_TRUE(a * b == b * a) << "a=" << (std::string)a << ", b=" << (std::string)b;
        EXPECT_TRUE(a * int128_t(1ll) == a) << "a=" << (std::string)a;
        EXPECT_TRUE((a * int128_t()).is_zero()) << "a=" << (std::string)a;
        EXPECT_TRUE(a * int128_t(-1ll) == -a) << "a=" << (std::string)a;
        EXPECT_TRUE(a * int128_t(2ll) == a + a) << "a=" << (std::string)a;
        // the sign rules, which the magnitude free implementation has to reproduce
        EXPECT_TRUE((-a) * b == -(a * b)) << "a=" << (std::string)a << ", b=" << (std::string)b;
        EXPECT_TRUE((-a) * (-b) == a * b) << "a=" << (std::string)a << ", b=" << (std::string)b;
        EXPECT_TRUE(a * (b + int128_t(1ll)) == a * b + a) << "a=" << (std::string)a << ", b=" << (std::string)b;
    }
}
TEST(int128_t, MultiplyByint128EdgeCases)
{
    const int128_t zero, one(1ll), minus_one(-1ll);
    const int128_t max(~0ull, 0x7FFFFFFFFFFFFFFFull);   // 2^127 - 1
    const int128_t min(0ull, 0x8000000000000000ull);    // -2^127

    EXPECT_TRUE(zero * zero == zero);
    EXPECT_TRUE(max * one == max);
    EXPECT_TRUE(max * minus_one == -max);
    // the most negative value has no positive counterpart, so negating it wraps back to itself
    EXPECT_TRUE(min * one == min);
    EXPECT_TRUE(min * minus_one == min);
    EXPECT_TRUE((min * int128_t(2ll)).is_zero());
    // (2^127-1)^2 is 2^254 - 2^128 + 1, whose low 128 bit are 1
    EXPECT_TRUE(max * max == one);
    EXPECT_TRUE(minus_one * minus_one == one);

    // the reference has to agree on every pairing of the extremes
    const int128_t values[] = {zero, one, minus_one, max, min, int128_t(0ull, 1ull), int128_t(~0ull, 0ull)};
    for (const auto& x : values) {
        for (const auto& y : values) {
            EXPECT_TRUE(x * y == ReferenceMultiply(x, y)) << "x=" << (std::string)x << ", y=" << (std::string)y;
        }
    }
}
// sqr(x) must be bit identical to x * x. square() skips the absolute value conversions that
// operator*= performs, which is valid because (-x)^2 == x^2 modulo 2^128. Raw bit patterns
// are used so both signs and the overflowing values are covered.
TEST(int128_t, SqrMatchesMultiply)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const uint64_t low = get_uint64_random();
        const uint64_t high = get_uint64_random();
        const int128_t x(low, high);
        EXPECT_TRUE(sqr(x) == x * x) << "low=" << low << ", high=" << high;
    }
}
TEST(int128_t, SqrEdgeCases)
{
    const int128_t values[] = {int128_t(0ull, 0ull),  int128_t(1ull, 0ull),   int128_t(0ull, 1ull),
                               int128_t(~0ull, ~0ull),  // -1
                               int128_t(0ull, 0x8000000000000000ull),  // most negative value
                               int128_t(~0ull, 0x7FFFFFFFFFFFFFFFull)};  // most positive value
    for (const auto& x : values) {
        EXPECT_TRUE(sqr(x) == x * x) << "x=" << (std::string)x;
    }
    // small values of both signs whose square is exact
    for (int64_t v = -1000; v <= 1000; ++v) {
        const int128_t x(v);
        EXPECT_TRUE(sqr(x) == int128_t(v * v)) << "v=" << v;
        EXPECT_TRUE(sqr(x) == x * x) << "v=" << v;
    }
}
// trunc(value1 / value2) computed in double is not a usable reference: the quotient of two 128 bit
// values needs more than the 53 bits a double holds, so the old oracle was wrong rather than the
// implementation, and is_similar_double() hid the difference on all 65536 iterations.
//
// The defining property of truncating integer division needs no reference implementation: with
// q*b + r == a, |r| < |b| and r taking the sign of a, the pair is unique. Magnitudes stay below
// 2^127 so that negation and abs() are exact.
TEST(int128_t, DivideByint128)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        int128_t a(get_uint64_random(), get_uint64_random() >> 1);

        // spread the divisor over the whole magnitude range so all three paths in operator/= are
        // reached: the power of two shift, the 64 bit divisor shortcut and the 128 bit long division
        const int32_t shift = static_cast<int32_t>(get_uint32_random() % 127);
        int128_t b = int128_t(get_uint64_random(), get_uint64_random() >> 1) >> shift;
        if ((get_uint32_random() & 7) == 0)
            b = int128_t(1ll) << static_cast<int32_t>(get_uint32_random() % 126);
        if (b.is_zero())
            continue;
        if (get_uint32_random() & 1)
            a = -a;
        if (get_uint32_random() & 1)
            b = -b;

        const int128_t q = a / b;
        const int128_t r = a % b;

        EXPECT_TRUE(q * b + r == a) << "a=" << (std::string)a << ", b=" << (std::string)b;
        EXPECT_TRUE(abs(r) < abs(b)) << "a=" << (std::string)a << ", b=" << (std::string)b;
        // truncation towards zero puts the remainder's sign on the dividend's side
        if (!r.is_zero())
            EXPECT_EQ(r.is_negative(), a.is_negative()) << "a=" << (std::string)a << ", b=" << (std::string)b;
        // and makes the quotient's sign the product of the operand signs
        if (!q.is_zero())
            EXPECT_EQ(q.is_negative(), a.is_negative() != b.is_negative()) << "a=" << (std::string)a << ", b=" << (std::string)b;

        int128_t c = a;
        c /= b;
        EXPECT_TRUE(c == q) << "a=" << (std::string)a << ", b=" << (std::string)b;
        int128_t d = a;
        d %= b;
        EXPECT_TRUE(d == r) << "a=" << (std::string)a << ", b=" << (std::string)b;

        // the magnitude of the quotient does not depend on the signs
        EXPECT_TRUE(abs(q) == abs((-a) / b)) << "a=" << (std::string)a << ", b=" << (std::string)b;
    }
}
TEST(int128_t, DivideByint128EdgeCases)
{
    const int128_t one(1ll), minus_one(-1ll);
    const int128_t max(~0ull, 0x7FFFFFFFFFFFFFFFull);   // 2^127 - 1
    const int128_t min(0ull, 0x8000000000000000ull);    // -2^127

    EXPECT_TRUE(max / max == one);
    EXPECT_TRUE(max / one == max);
    EXPECT_TRUE(max / minus_one == -max);
    EXPECT_TRUE(min / one == min);
    EXPECT_TRUE(min / min == one);

    // -2^127 / -1 is 2^127, which is one past the top of the range. The builtin types make this
    // undefined; here it wraps back to the most negative value, and the remainder stays zero.
    EXPECT_TRUE(min / minus_one == min);
    EXPECT_TRUE((min % minus_one).is_zero());
    // the division identity survives it
    EXPECT_TRUE((min / minus_one) * minus_one + (min % minus_one) == min);

    EXPECT_TRUE(min / int128_t(2ll) == int128_t(0ull, 0xC000000000000000ull));  // -2^126
    EXPECT_TRUE((int128_t() / max).is_zero());
    EXPECT_TRUE((one / max).is_zero());
    EXPECT_TRUE(one % max == one);

    // division by zero is reported, not silently wrong
    EXPECT_THROW((void)(max / int128_t()), std::logic_error);
    EXPECT_THROW((void)(max % int128_t()), std::logic_error);
    EXPECT_THROW((void)(max / 0ll), std::logic_error);
}
// The assertion was guarded by "if (int128_res == res) continue;", so it only ever ran when it was
// about to fail and the test executed zero assertions on a passing run. The operands were also a
// 64 bit numerator and a 32 bit divisor, and value1 % value2 is undefined for INT64_MIN % -1.
TEST(int128_t, ModuloByint128)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const int64_t value1 = get_int64_random();
        const int32_t value2 = get_int32_random();
        if (value2 == 0)
            continue;
        // INT64_MIN % -1 is undefined for the builtin types, so it cannot serve as a reference.
        // The 128 bit behaviour at that point is pinned by DivideByint128EdgeCases instead.
        if (value1 == INT64_MIN && value2 == -1)
            continue;

        const int128_t f1 = value1, f2 = value2;
        EXPECT_EQ(static_cast<int64_t>(f1 % f2), value1 % value2) << "value1=" << value1 << ", value2=" << value2;
        EXPECT_EQ(static_cast<int64_t>(f1 / f2), value1 / value2) << "value1=" << value1 << ", value2=" << value2;
        // the 128 bit result has to be the sign extension of the 64 bit one, not merely equal in
        // its low QWORD
        EXPECT_TRUE(f1 % f2 == int128_t(value1 % value2)) << "value1=" << value1 << ", value2=" << value2;
        EXPECT_TRUE((f1 / f2) * f2 + (f1 % f2) == f1) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(int128_t, Compareint128)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = floor(get_double_random(1, 127));
        double value2 = floor(get_double_random(1, 127));
        bool res = value1 > value2;
        int128_t f1 = value1;
        int128_t f2 = value2;
        if (check_overflow_int128(value1) || check_overflow_int128(value2))
            continue;

        bool int128_res = f1 > f2;
        EXPECT_TRUE(int128_res == res) << "operator>: "
                                       << "value1=" << value1 << ", value2=" << value2;

        res = value1 >= value2;
        int128_res = f1 >= f2;
        EXPECT_TRUE(int128_res == res) << "operator>=: "
                                       << "value1=" << value1 << ", value2=" << value2;

        res = value1 < value2;
        int128_res = f1 < f2;
        EXPECT_TRUE(int128_res == res) << "operator<: "
                                       << "value1=" << value1 << ", value2=" << value2;

        res = value1 <= value2;
        int128_res = f1 <= f2;
        EXPECT_TRUE(int128_res == res) << "operator<=: "
                                       << "value1=" << value1 << ", value2=" << value2;
    }
}
// Operands were 64 bit values, so no iteration ever carried into the high QWORD, and the postfix
// form's return value was never looked at. value1 + 1 is also undefined at INT64_MAX.
TEST(int128_t, OperatorPlusPlus)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const int128_t f1(get_uint64_random(), get_uint64_random());
        const int128_t expected = f1 + int128_t(1ll);

        int128_t f2 = f1;
        // the postfix form returns the value from before the increment
        EXPECT_TRUE(f2++ == f1) << "f1=" << (std::string)f1;
        EXPECT_TRUE(f2 == expected) << "f1=" << (std::string)f1;

        f2 = f1;
        // the prefix form returns the value after it, as a reference to the object itself
        EXPECT_TRUE(++f2 == expected) << "f1=" << (std::string)f1;
        EXPECT_TRUE(f2 == expected) << "f1=" << (std::string)f1;
    }
}
TEST(int128_t, OperatorMinusMinus)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const int128_t f1(get_uint64_random(), get_uint64_random());
        const int128_t expected = f1 - int128_t(1ll);

        int128_t f2 = f1;
        EXPECT_TRUE(f2-- == f1) << "f1=" << (std::string)f1;
        EXPECT_TRUE(f2 == expected) << "f1=" << (std::string)f1;

        f2 = f1;
        EXPECT_TRUE(--f2 == expected) << "f1=" << (std::string)f1;
        EXPECT_TRUE(f2 == expected) << "f1=" << (std::string)f1;
    }
}
TEST(int128_t, IncrementDecrementBoundaries)
{
    const int128_t max(~0ull, 0x7FFFFFFFFFFFFFFFull);   // 2^127 - 1
    const int128_t min(0ull, 0x8000000000000000ull);    // -2^127

    // the carry has to cross the QWORD boundary
    int128_t a(~0ull, 0ull);  // 2^64 - 1
    ++a;
    EXPECT_TRUE(a == int128_t(0ull, 1ull));
    int128_t b(0ull, 1ull);   // 2^64
    --b;
    EXPECT_TRUE(b == int128_t(~0ull, 0ull));

    // and crossing zero flips every bit of the high QWORD
    int128_t c(-1ll);
    ++c;
    EXPECT_TRUE(c.is_zero());
    int128_t d;
    --d;
    EXPECT_TRUE(d == int128_t(-1ll));

    // overflow wraps silently, as it does for the builtin types
    int128_t e = max;
    ++e;
    EXPECT_TRUE(e == min);
    int128_t f = min;
    --f;
    EXPECT_TRUE(f == max);
}
// The previous version asserted f1 == f1, which is true of any implementation that manages to read
// its own members twice, and never compared two distinct objects. In particular it could not tell a
// comparison of the whole 128 bit value from one that only looks at a single QWORD.
TEST(int128_t, OperatorEqual)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const uint64_t al = get_uint64_random(), ah = get_uint64_random();
        const uint64_t bl = get_uint64_random(), bh = get_uint64_random();
        const int128_t a(al, ah), b(bl, bh);

        // equality is bit identical for both instantiations, no sign test is involved
        EXPECT_EQ(a == b, (al == bl) && (ah == bh)) << "a=" << (std::string)a << ", b=" << (std::string)b;

        const int128_t copy = a;
        EXPECT_TRUE(a == copy) << "a=" << (std::string)a;
        // a value differing in exactly one QWORD must not compare equal, in either position
        EXPECT_FALSE(a == int128_t(al ^ 1ull, ah)) << "a=" << (std::string)a;
        EXPECT_FALSE(a == int128_t(al, ah ^ 1ull)) << "a=" << (std::string)a;
        // including when the difference is the sign bit alone
        EXPECT_FALSE(a == int128_t(al, ah ^ (1ull << 63))) << "a=" << (std::string)a;

        // the templated overload sign extends the scalar, so a negative int64_t compares equal only
        // when the high QWORD is all ones, and a positive one only when it is zero
        const int64_t scalar = static_cast<int64_t>(al);
        EXPECT_EQ(a == scalar, ah == SignExtension(scalar)) << "a=" << (std::string)a << ", scalar=" << scalar;
    }
}
// operator!= is declared separately rather than being the compiler's rewrite of operator==, so it
// can disagree with it. Every case below checks both that it is correct and that the two agree.
TEST(int128_t, OperatorNotEqual)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const uint64_t al = get_uint64_random(), ah = get_uint64_random();
        const uint64_t bl = get_uint64_random(), bh = get_uint64_random();
        const int128_t a(al, ah), b(bl, bh);

        EXPECT_EQ(a != b, (al != bl) || (ah != bh)) << "a=" << (std::string)a << ", b=" << (std::string)b;
        EXPECT_EQ(a != b, !(a == b)) << "a=" << (std::string)a << ", b=" << (std::string)b;

        const int128_t copy = a;
        EXPECT_FALSE(a != copy) << "a=" << (std::string)a;
        EXPECT_TRUE(a != int128_t(al ^ 1ull, ah)) << "a=" << (std::string)a;
        EXPECT_TRUE(a != int128_t(al, ah ^ 1ull)) << "a=" << (std::string)a;

        const int64_t scalar = static_cast<int64_t>(al);
        EXPECT_EQ(a != scalar, !(a == scalar)) << "a=" << (std::string)a << ", scalar=" << scalar;
    }
}
TEST(int128_t, OperatorEqualEdgeCases)
{
    const int128_t zero, one(1ll), minus_one(-1ll);
    const int128_t max(~0ull, 0x7FFFFFFFFFFFFFFFull);
    const int128_t min(0ull, 0x8000000000000000ull);

    EXPECT_TRUE(zero == int128_t());
    EXPECT_FALSE(zero == one);
    EXPECT_TRUE(zero != one);
    // max and min differ in the sign bit alone
    EXPECT_FALSE(max == min);
    EXPECT_TRUE(max != min);
    // -1 is all ones in both QWORDs, 2^64-1 fills only the low one
    EXPECT_FALSE(minus_one == int128_t(~0ull, 0ull));
    EXPECT_TRUE(minus_one == int128_t(~0ull, ~0ull));

    // a scalar is sign extended before the comparison
    EXPECT_TRUE(minus_one == -1ll);
    EXPECT_FALSE(int128_t(~0ull, 0ull) == -1ll);  // 2^64-1, positive, not -1
    EXPECT_TRUE(int128_t(5ull, 0ull) == 5ll);
    EXPECT_FALSE(int128_t(5ull, 1ull) == 5ll);
}
TEST(int128_t, log10)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = floor(1.0 + fabs(get_double_random(1, 127)));
        int64_t res = (int64_t)floor(log10(value1));  // double doesn't lose any bits with this operation!
        if (check_overflow_int128(value1)) {
            continue;
        }
        int128_t i1 = value1;
        int64_t i_res = log10(i1);
        EXPECT_EQ(i_res, res) << "double value1=" << value1;
    }
}
TEST(int128_t, log)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = 1.0 + floor(fabs(get_double_random(0, 127)));
        int64_t res = (int64_t)floor(log(value1));  // double doesn't lose any bits with this operation!
        if (check_overflow_int128(value1)) {
            continue;
        }
        int128_t i1 = value1;
        int64_t i_res = log(i1);
        EXPECT_EQ(i_res, res) << "double value1=" << value1;
    }
}
TEST(int128_t, sqrt)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = floor(1.0 + fabs(get_double_random(1, 127)));
        if (check_overflow_int128(value1)) {
            continue;
        }
        int128_t i1 = value1;
        uint64_t i_res = sqrt(i1);
        // floor(sqrt(value1)) is not a usable reference at this range, the root needs more than the
        // 53 bits a double holds. is_exact_isqrt() checks r^2 <= x < (r+1)^2 in 128 bit instead.
        EXPECT_TRUE(is_exact_isqrt(i_res, i1)) << "double value1=" << value1 << ", sqrt=" << i_res;
        // the conversion of the operand has to be exact, otherwise the check above passes vacuously
        EXPECT_EQ(static_cast<double>(i1), value1);
    }
}
// The old version declared "int128_t int128_res = static_cast<double>(i2);" and handed that to
// EXPECT_DOUBLE_EQ, which converted it straight back to a double: a round trip through two
// conversions that happened to work only because the operands were tiny. The base never exceeded 7
// and the exponent 15, so the 128 bit range was never approached and negative bases never appeared.
//
// Repeated multiplication is an exact reference at any size, and it stays valid past the point
// where the result wraps because pow() and operator*= truncate identically.
TEST(int128_t, pow)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const int64_t base = static_cast<int64_t>(get_uint32_random() % 200) - 100;  // [-100, 99]
        const uint32_t expo = get_uint32_random() % 40;
        const int128_t x(base);

        int128_t reference(1ll);
        for (auto k = 0u; k < expo; ++k) {
            reference *= x;
        }

        EXPECT_TRUE(pow(x, expo) == reference) << "base=" << base << ", expo=" << expo;
        // an even power of a negative base is positive and matches the power of its magnitude
        if (base != INT64_MIN && (expo & 1) == 0)
            EXPECT_TRUE(pow(x, expo) == pow(-x, expo)) << "base=" << base << ", expo=" << expo;
        else if (base != INT64_MIN)
            EXPECT_TRUE(pow(x, expo) == -pow(-x, expo)) << "base=" << base << ", expo=" << expo;
    }
}
TEST(int128_t, powEdgeCases)
{
    // any base to the zero power is one, including zero, matching pow(double, double)
    EXPECT_TRUE(pow(int128_t(), 0u) == int128_t(1ll));
    EXPECT_TRUE(pow(int128_t(5ll), 0u) == int128_t(1ll));
    EXPECT_TRUE(pow(int128_t(-5ll), 0u) == int128_t(1ll));
    // zero to any positive power is zero
    EXPECT_TRUE(pow(int128_t(), 7u).is_zero());
    // one and minus one alternate
    EXPECT_TRUE(pow(int128_t(-1ll), 3u) == int128_t(-1ll));
    EXPECT_TRUE(pow(int128_t(-1ll), 4u) == int128_t(1ll));
    // 2^126 is the largest power of two that stays positive
    EXPECT_TRUE(pow(int128_t(2ll), 126u) == int128_t(0ull, 0x4000000000000000ull));
    // 2^127 wraps to the most negative value, and 2^128 to zero
    EXPECT_TRUE(pow(int128_t(2ll), 127u) == int128_t(0ull, 0x8000000000000000ull));
    EXPECT_TRUE(pow(int128_t(2ll), 128u).is_zero());
    // 3^80 overflows, and must truncate the same way repeated multiplication does
    int128_t reference(1ll);
    for (auto i = 0u; i < 80u; ++i) {
        reference *= int128_t(3ll);
    }
    EXPECT_TRUE(pow(int128_t(3ll), 80u) == reference);
}

/**********************************************************************
 * int128_t regression tests
 *
 * Each test below pins down a defect that was fixed. They deliberately use
 * operand ranges and sign combinations the original tests never reached.
 ***********************************************************************/

// Division reduced both operands to magnitudes, divided, negated, and then subtracted one whenever
// the remainder was non zero, which floors. The builtin integer types truncate towards zero, and
// operator%= was already truncating, so the two disagreed and the division identity did not hold.
TEST(int128_t, DivideTruncatesTowardsZero)
{
    struct {
        int64_t a, b, q, r;
    } cases[] = {{-7, 2, -3, -1}, {7, -2, -3, 1}, {-7, -2, 3, -1}, {7, 2, 3, 1},
                 {-100, 7, -14, -2}, {-1, 2, 0, -1}, {-9, 3, -3, 0}, {-10, 4, -2, -2}};
    for (const auto& c : cases) {
        const int128_t a(c.a), b(c.b);
        EXPECT_EQ(static_cast<int64_t>(a / b), c.q) << c.a << " / " << c.b;
        EXPECT_EQ(static_cast<int64_t>(a % b), c.r) << c.a << " % " << c.b;
        // the quotient and the remainder must agree with each other
        EXPECT_TRUE((a / b) * b + (a % b) == a) << c.a << ", " << c.b;
    }

    // every sign combination, cross checked against int64_t
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const int64_t a = get_int64_random() >> 8;
        const int64_t b = get_int32_random();
        if (b == 0)
            continue;
        const int128_t f1(a), f2(b);
        EXPECT_EQ(static_cast<int64_t>(f1 / f2), a / b) << "a=" << a << ", b=" << b;
        EXPECT_EQ(static_cast<int64_t>(f1 % f2), a % b) << "a=" << a << ", b=" << b;
        EXPECT_TRUE((f1 / f2) * f2 + (f1 % f2) == f1) << "a=" << a << ", b=" << b;
    }
}
// The 64 bit fast path shortcut small numerators with a signed comparison, so any negative value
// compared less than a positive divisor and the result was zeroed before the sign was applied.
TEST(int128_t, DivideByScalarWithNegativeNumerator)
{
    const int128_t a(static_cast<int64_t>(-100));
    EXPECT_EQ(static_cast<int64_t>(a / 7), -14);
    EXPECT_EQ(static_cast<int64_t>(a / 7ull), -14);
    EXPECT_EQ(static_cast<int64_t>(a / 7u), -14);
    EXPECT_EQ(static_cast<int64_t>(a % 7), -2);
    // a power of two divisor takes a different branch and must truncate too, not shift
    EXPECT_EQ(static_cast<int64_t>(int128_t(static_cast<int64_t>(-7)) / 2), -3);
    EXPECT_EQ(static_cast<int64_t>(int128_t(static_cast<int64_t>(-8)) / 2), -4);
}
// A divisor whose magnitude needs more than 64 bit exercises the 128 bit division path.
TEST(int128_t, DivideByWideInt128)
{
    // 2^66 / (3 * 2^64) truncates to 1, and the remainder lives entirely in the upper QWORD,
    // which the old rounding correction could not see
    const int128_t big(0ull, 4ull);          // 2^66
    const int128_t divisor(0ull, 3ull);      // 3 * 2^64
    EXPECT_EQ(static_cast<int64_t>(big / divisor), 1);
    EXPECT_EQ(static_cast<int64_t>((-big) / divisor), -1);
    EXPECT_TRUE((big / divisor) * divisor + (big % divisor) == big);
    EXPECT_TRUE(((-big) / divisor) * divisor + ((-big) % divisor) == -big);

    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        int128_t a(get_uint64_random(), get_uint64_random() >> 1);
        int128_t b(get_uint64_random(), get_uint64_random() >> (get_uint32_random() % 63));
        if (b.is_zero())
            continue;
        if (get_uint32_random() & 1)
            a = -a;
        if (get_uint32_random() & 1)
            b = -b;
        const int128_t q = a / b, r = a % b;
        EXPECT_TRUE(q * b + r == a) << "a=" << (std::string)a << ", b=" << (std::string)b;
        // the remainder takes the sign of the dividend and stays smaller in magnitude
        EXPECT_TRUE(abs(r) < abs(b)) << "a=" << (std::string)a << ", b=" << (std::string)b;
        if (!r.is_zero())
            EXPECT_EQ(r.is_negative(), a.is_negative()) << "a=" << (std::string)a;
    }
}
// sqrt only rejected negatives, then called log2, which throws on zero. sqrt is noexcept, so a
// zero input terminated the process.
TEST(int128_t, SqrtOfZeroDoesNotThrow)
{
    EXPECT_EQ(sqrt(int128_t()), 0ull);
    EXPECT_EQ(sqrt(int128_t(0ull)), 0ull);
    EXPECT_EQ(sqrt(int128_t(1ull)), 1ull);
    EXPECT_EQ(sqrt(int128_t(static_cast<int64_t>(-1))), 0ull);
    // the log functions do throw, that part is intentional
    EXPECT_THROW((void)log2(int128_t()), std::domain_error);
    EXPECT_THROW((void)log(int128_t()), std::domain_error);
    EXPECT_THROW((void)log10(int128_t()), std::domain_error);
    EXPECT_THROW((void)log2(int128_t(static_cast<int64_t>(-1))), std::domain_error);
}
// operator^ was implemented with &=, so it computed a bitwise AND.
TEST(int128_t, BitwiseOperators)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const uint64_t al = get_uint64_random(), ah = get_uint64_random();
        const uint64_t bl = get_uint64_random(), bh = get_uint64_random();
        const int128_t a(al, ah), b(bl, bh);
        EXPECT_TRUE((a ^ b) == int128_t(al ^ bl, ah ^ bh));
        EXPECT_TRUE((a & b) == int128_t(al & bl, ah & bh));
        EXPECT_TRUE((a | b) == int128_t(al | bl, ah | bh));
        EXPECT_TRUE((a ^ a).is_zero());
    }
}
// The terminator sat at index 32 of a 45 byte buffer with the digits written backwards from there,
// so a long value ran off the front. A negative value needs one more character than an unsigned
// one of the same width. Separately, the fast path tested high == 0 and formatted with %lld, which
// printed every value in [2^63, 2^64) as a negative number.
TEST(int128_t, ToStringLargeValues)
{
    const char* values[] = {"9223372036854775808",                       // 2^63, the %lld trap
                            "18446744073709551615",                      // 2^64-1, still positive
                            "18446744073709551616",                      // 2^64
                            "-9223372036854775809",                      // just past int64_t's range
                            "170141183460469231731687303715884105727",   // INT128_MAX, 39 digits
                            "-170141183460469231731687303715884105728"}; // INT128_MIN, 39 digits and a sign

    // a value that fits in an int64_t takes the snprintf path, which returns the buffer start
    const char* buffer_start = static_cast<char*>(int128_t(1ull));
    for (auto i = 0u; i < array_length(values); ++i) {
        const int128_t v = values[i];
        const char* p = static_cast<char*>(v);
        EXPECT_STREQ(p, values[i]);
        EXPECT_GE(p, buffer_start) << "wrote before the start of the buffer: " << values[i];
    }

    // full round trip over the whole 128 bit range, both signs
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const int128_t v(get_uint64_random(), get_uint64_random());
        const std::string s = (std::string)v;
        EXPECT_TRUE(int128_t(s) == v) << "s=" << s;
    }
}
// The conversions used to truncate the bits that did not fit the mantissa. strtod and strtof are
// correctly rounded by the C standard, so the decimal string of a value is an exact reference.
TEST(int128_t, ConversionToFloatingPointIsCorrectlyRounded)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        int128_t v(get_uint64_random(), get_uint64_random() >> (get_uint32_random() % 64));
        if (v.is_zero())
            continue;
        if (get_uint32_random() & 1)
            v = -v;
        const std::string s = (std::string)v;
        EXPECT_EQ(static_cast<double>(v), strtod(s.c_str(), nullptr)) << "value=" << s;
        EXPECT_EQ(static_cast<float>(v), strtof(s.c_str(), nullptr)) << "value=" << s;
    }
}
TEST(int128_t, ConversionToFloatingPointEdgeCases)
{
    // the most negative value: its absolute value keeps the sign bit set
    const int128_t min_val(0ull, 0x8000000000000000ull);
    EXPECT_EQ(static_cast<double>(min_val), -170141183460469231731687303715884105728.0);
    EXPECT_EQ(static_cast<float>(min_val), -170141183460469231731687303715884105728.0f);
    // INT128_MAX rounds up to 2^127, which both types can hold
    const int128_t max_val(~0ull, 0x7FFFFFFFFFFFFFFFull);
    EXPECT_EQ(static_cast<double>(max_val), 170141183460469231731687303715884105728.0);
    EXPECT_FALSE(std::isinf(static_cast<float>(max_val)));
    // exact values need no rounding
    EXPECT_EQ(static_cast<double>(int128_t(0ull, 1ull)), 18446744073709551616.0);
    EXPECT_EQ(static_cast<double>(int128_t(static_cast<int64_t>(-1))), -1.0);
}
// Shifting by 128 or more used to be taken modulo 64 by the hardware and leave bits behind.
// An arithmetic right shift has to leave copies of the sign bit, not zeros.
TEST(int128_t, ShiftBeyondWidth)
{
    volatile int32_t shift = 128;
    const int128_t positive(0ull, 0x4000000000000000ull);
    const int128_t negative(static_cast<int64_t>(-1));
    EXPECT_TRUE((positive >> shift).is_zero());
    EXPECT_EQ(static_cast<int64_t>(negative >> shift), -1);
    EXPECT_TRUE((positive << shift).is_zero());
    shift = 200;
    EXPECT_TRUE((positive >> shift).is_zero());
    EXPECT_EQ(static_cast<int64_t>(negative >> shift), -1);
    EXPECT_TRUE((positive << shift).is_zero());
    // the boundary either side still behaves
    shift = 127;
    EXPECT_EQ(static_cast<int64_t>(negative >> shift), -1);
    EXPECT_TRUE((int128_t(1ull) << shift) == int128_t(0ull, 0x8000000000000000ull));
}
// A comma used to continue the parse loop without advancing the read pointer, spinning forever.
// A high byte handed straight to isspace is outside the domain it accepts.
TEST(int128_t, ConstructorFromStringEdgeCases)
{
    EXPECT_TRUE(int128_t("1,000") == int128_t(1000ull));
    EXPECT_TRUE(int128_t("-1,000") == int128_t(static_cast<int64_t>(-1000)));
    EXPECT_TRUE(int128_t("+5") == int128_t(5ull));
    EXPECT_TRUE(int128_t("-0x10") == int128_t(static_cast<int64_t>(-16)));
    EXPECT_TRUE(int128_t("-0").is_zero());
    const char high_byte[] = {'1', '2', static_cast<char>(0xB5), '3', 0};
    EXPECT_TRUE(int128_t(high_byte) == int128_t(12ull));
}
// One overload per fixed width type left long and unsigned long ambiguous on MSVC, and a scalar on
// the left hand side of a binary operator could not resolve at all.
TEST(int128_t, IntegralTypesAndScalarOnTheLeft)
{
    EXPECT_TRUE(int128_t(static_cast<unsigned long>(7)) == int128_t(7ull));
    EXPECT_TRUE(int128_t(static_cast<long>(-7)) == int128_t(static_cast<int64_t>(-7)));
    EXPECT_TRUE(int128_t(static_cast<short>(-3)) == int128_t(static_cast<int64_t>(-3)));
    EXPECT_TRUE(int128_t(true) == int128_t(1ull));
    EXPECT_TRUE(int128_t(static_cast<size_t>(9)) == int128_t(9ull));

    const int128_t x(10ull);
    EXPECT_TRUE((1 + x) == int128_t(11ull));
    EXPECT_TRUE((100 - x) == int128_t(90ull));
    EXPECT_TRUE((3 * x) == int128_t(30ull));
    EXPECT_TRUE((100 / x) == int128_t(10ull));
    EXPECT_TRUE((105 % x) == int128_t(5ull));
    EXPECT_TRUE((3 ^ x) == int128_t(9ull));
    EXPECT_TRUE((-100 / x) == int128_t(static_cast<int64_t>(-10)));
    // the int128_t on the left forms must still resolve
    EXPECT_TRUE((x + 1) == int128_t(11ull));
    EXPECT_TRUE((x + x) == int128_t(20ull));
}
// int128_t and uint128_t are now aliases of a single class template. Nothing but this test pins the
// properties the merge had to preserve: the storage is unchanged, the 16 byte alignment survived the
// template, and the two types still refuse to convert into one another. Trivial copyability is new,
// a consequence of defaulting the copy and move members.
TEST(int128_t, LayoutAndTraits)
{
    static_assert(std::is_same_v<int128_t, int128_base<true>>, "int128_t is the signed instantiation");
    static_assert(int128_t::is_signed, "int128_t reports itself as signed");
    static_assert(sizeof(int128_t) == 16, "the object holds exactly two QWORDs");
    static_assert(alignof(int128_t) == 16, "FP128_ALIGN16 has to survive the template");
    static_assert(std::is_standard_layout_v<int128_t>, "the division helpers alias the members as an array");
    static_assert(std::is_trivially_copyable_v<int128_t>, "the copy and move members are defaulted");
    // getting from one signedness to the other needs two user defined conversions, which the
    // language never performs implicitly, and leaves the explicit form ambiguous
    static_assert(!std::is_convertible_v<uint128_t, int128_t>, "no silent unsigned to signed conversion");
    static_assert(!std::is_constructible_v<int128_t, uint128_t>, "and no unambiguous explicit one either");
    SUCCEED();
}
// Every math function on these types is a hidden friend, found by argument dependent lookup. A value
// whose low QWORD is zero catches the failure mode where a call falls through to the narrower
// fp128::log2(uint64_t) overload through the implicit conversion to uint64_t: that overload would
// see a zero and answer zero.
TEST(int128_t, MathFunctionsSeeTheWholeValue)
{
    const int128_t x(0ull, 1ull << 40);  // 2^104, the low QWORD is zero

    EXPECT_EQ(log2(x), 104ull);
    EXPECT_EQ(lzcnt128(x), 23ull);
    EXPECT_EQ(log10(x), 31ull);
    EXPECT_EQ(log(x), 72ull);
    EXPECT_EQ(sqrt(x), 1ull << 52);
    EXPECT_TRUE(pow(int128_t(2ll), 104u) == x);
    EXPECT_TRUE(sqr(int128_t(1ull << 52)) == x);
}
// abs() is a hidden friend constrained to the signed instantiation, so it is never found by ordinary
// unqualified lookup and can never hide ::abs() for the builtin types inside namespace fp128.
TEST(int128_t, AbsResolvesToTheHiddenFriend)
{
    // the result type has to be int128_t, not the int that ::abs() would return
    static_assert(std::is_same_v<decltype(abs(int128_t(1ll))), int128_t>, "abs must not fall through to ::abs");

    const int128_t x("-170141183460469231731687303715884105727");  // -(2^127 - 1)
    EXPECT_TRUE(abs(x) == int128_t("170141183460469231731687303715884105727"));
    EXPECT_TRUE(abs(int128_t(-5ll)) == int128_t(5ll));
    EXPECT_TRUE(abs(int128_t(5ll)) == int128_t(5ll));
    EXPECT_TRUE(abs(int128_t()).is_zero());
}
// operator*= multiplies the two's complement bit patterns directly instead of converting both
// operands to magnitudes and reapplying the sign afterwards. A truncated 128 bit product is a
// multiplication modulo 2^128, and in that ring the bit pattern of a negative value is its value,
// so the result is bit identical to what the magnitude based version produced.
TEST(int128_t, MultiplySignCombinations)
{
    const int128_t a("123456789012345678901234567890");
    const int128_t b(-987654321ll);

    EXPECT_TRUE((a * b) == -(a * -b));
    EXPECT_TRUE((-a * b) == (a * -b));
    EXPECT_TRUE((-a * -b) == (a * b));
    EXPECT_TRUE((a * b).is_negative());
    EXPECT_TRUE((-a * b).is_positive());

    // the most negative value is its own negation, which a magnitude based implementation has to
    // treat as a special case
    const int128_t most_negative(0ull, 0x8000000000000000ull);
    EXPECT_TRUE((most_negative * int128_t(1ll)) == most_negative);
    EXPECT_TRUE((most_negative * int128_t(-1ll)) == most_negative);
    EXPECT_TRUE((most_negative * int128_t(2ll)).is_zero());
}
// The right hand side is allowed to alias this object, which the implementation guards against by
// snapshotting both operands before it writes anything.
TEST(int128_t, MultiplyInPlaceAliasing)
{
    int128_t a(-123456789ll);
    const int128_t expected = a * a;
    a *= a;
    EXPECT_TRUE(a == expected);
    EXPECT_TRUE(a == int128_t(15241578750190521ll));
}
// The scalar overload used to convert to int128_t and take the full path. It now uses the same 64
// bit shortcut the unsigned type does, which is valid for a negative object for the same modulo
// 2^128 reason operator*=(const int128_t&) is.
TEST(int128_t, MultiplyByScalarKeepsTheSign)
{
    int128_t a(-1000ll);
    a *= 3u;
    EXPECT_TRUE(a == int128_t(-3000ll));

    int128_t b("-123456789012345678901234567890");
    b *= 7ull;
    EXPECT_TRUE(b == int128_t("-864197523086419752308641975230"));

    // a negative scalar still goes the long way round
    int128_t c(-1000ll);
    c *= -3;
    EXPECT_TRUE(c == int128_t(3000ll));

    // and so does a floating point one
    int128_t d(-1000ll);
    d *= 2.5;
    EXPECT_TRUE(d == int128_t(-2000ll));
}
// operator/= detects a power of two divisor and shifts instead of dividing. The shift is a logical
// one applied to the magnitude, so it has to be correct on both sides of the 64 bit boundary and
// for a divisor of exactly one, which shifts by zero.
TEST(int128_t, DivideByPowersOfTwo)
{
    const int128_t x("170141183460469231731687303715884105727");  // 2^127 - 1

    EXPECT_TRUE((x / int128_t(1ll)) == x);
    EXPECT_TRUE((x / int128_t(2ll)) == (x >> 1));
    EXPECT_TRUE((x / int128_t(1ull << 63)) == (x >> 63));
    EXPECT_TRUE((x / int128_t(0ull, 1ull)) == (x >> 64));         // 2^64
    EXPECT_TRUE((x / int128_t(0ull, 1ull << 36)) == (x >> 100));  // 2^100

    // a negative dividend truncates towards zero, which an arithmetic shift of a negative value
    // does not do, so the magnitude has to be shifted and negated afterwards
    EXPECT_TRUE((-x / int128_t(2ll)) == -(x / int128_t(2ll)));
    EXPECT_TRUE((-x / int128_t(0ull, 1ull)) == -(x >> 64));

    // the magnitude of the most negative value keeps its sign bit set, which is why the shift works
    // on the raw QWORDs instead of going through the signed operators
    const int128_t most_negative(0ull, 0x8000000000000000ull);  // -2^127
    EXPECT_TRUE((most_negative / int128_t(1ll)) == most_negative);
    EXPECT_TRUE((most_negative / int128_t(2ll)) == int128_t(0ull, 0xC000000000000000ull));  // -2^126
    EXPECT_TRUE((most_negative / int128_t(0ull, 1ull)) == -int128_t(1ull << 63));           // -2^63
}
// Public members that no other test in this file reaches. All of them work today, so this pins
// current behaviour rather than a past defect.
TEST(int128_t, UnaryOperatorsAndAccessors)
{
    const int128_t x(0x0123456789ABCDEFull, 0xFEDCBA9876543210ull);
    const int128_t min(0ull, 0x8000000000000000ull);

    // unary plus is a copy, unary minus is a two's complement negation
    EXPECT_TRUE(+x == x);
    EXPECT_TRUE(-(-x) == x);
    EXPECT_TRUE(x + (-x) == int128_t());
    // ~x == -x - 1 is the defining relation of two's complement
    EXPECT_TRUE(~x == -x - int128_t(1ll));
    EXPECT_TRUE(~(~x) == x);
    EXPECT_TRUE(~int128_t() == int128_t(-1ll));
    // the most negative value is its own negation, having no positive counterpart
    EXPECT_TRUE(-min == min);

    // operator bool and operator!, which are each other's opposite
    EXPECT_TRUE(static_cast<bool>(x));
    EXPECT_FALSE(static_cast<bool>(int128_t()));
    EXPECT_FALSE(!x);
    EXPECT_TRUE(!int128_t());
    // a value whose low QWORD is zero is still non zero
    EXPECT_TRUE(static_cast<bool>(int128_t(0ull, 1ull)));

    // the sign predicates, including the boundary between them
    EXPECT_TRUE(int128_t(1ll).is_positive());
    EXPECT_TRUE(int128_t().is_positive());  // zero counts as positive
    EXPECT_FALSE(int128_t().is_negative());
    EXPECT_TRUE(int128_t(-1ll).is_negative());
    EXPECT_TRUE(min.is_negative());
    EXPECT_TRUE(int128_t(~0ull, 0x7FFFFFFFFFFFFFFFull).is_positive());
    // an integer type has no fraction
    EXPECT_TRUE(x.is_int());

    // get_bit over both QWORDs, and the sign bit at the top
    EXPECT_EQ(int128_t(1ll).get_bit(0), 1);
    EXPECT_EQ(int128_t(1ll).get_bit(1), 0);
    EXPECT_EQ(int128_t(0ull, 1ull).get_bit(64), 1);
    EXPECT_EQ(min.get_bit(127), 1);
    EXPECT_EQ(int128_t(-1ll).get_bit(127), 1);
    EXPECT_EQ(int128_t(1ll).get_bit(127), 0);

    // one() and get_components()
    EXPECT_TRUE(int128_t::one() == int128_t(1ll));
    uint64_t low = 0, high = 0;
    x.get_components(low, high);
    EXPECT_EQ(low, 0x0123456789ABCDEFull);
    EXPECT_EQ(high, 0xFEDCBA9876543210ull);

    // hex() prints both QWORDs, high first, zero padded
    EXPECT_STREQ(x.hex(), "0xFEDCBA98765432100123456789ABCDEF");
    EXPECT_STREQ(int128_t().hex(), "0x00000000000000000000000000000000");
    EXPECT_STREQ(int128_t(-1ll).hex(), "0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF");

    // to_string() has to agree with the std::string conversion
    const int128_t big("-170141183460469231731687303715884105728");
    EXPECT_STREQ(big.to_string(), "-170141183460469231731687303715884105728");
    EXPECT_STREQ(big.to_string(), static_cast<std::string>(big).c_str());

    // long double goes through the double conversion
    EXPECT_EQ(static_cast<long double>(int128_t(-1ll)), -1.0L);
    EXPECT_EQ(static_cast<long double>(min), static_cast<long double>(static_cast<double>(min)));
}
// The user defined literal is declared at namespace scope in int128_t.h and was never exercised.
TEST(int128_t, UserDefinedLiteral)
{
    EXPECT_TRUE(1234_int128 == int128_t(1234ll));
    EXPECT_TRUE(-1234_int128 == int128_t(-1234ll));
    EXPECT_TRUE(0xDEADBEAF_int128 == int128_t(0xDEADBEAFll));
    EXPECT_TRUE(170141183460469231731687303715884105727_int128 == int128_t(~0ull, 0x7FFFFFFFFFFFFFFFull));
    // the operator takes the raw characters, so digit grouping survives
    EXPECT_TRUE(1'000'000_int128 == int128_t(1000000ll));
}
// The compound bitwise and shift assignments, and the scalar forms of the bitwise operators.
TEST(int128_t, CompoundBitwiseAndShiftOperators)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const uint64_t al = get_uint64_random(), ah = get_uint64_random();
        const uint64_t bl = get_uint64_random(), bh = get_uint64_random();
        const int128_t a(al, ah), b(bl, bh);

        int128_t c = a;
        c &= b;
        EXPECT_TRUE(c == (a & b)) << "a=" << (std::string)a << ", b=" << (std::string)b;
        c = a;
        c |= b;
        EXPECT_TRUE(c == (a | b)) << "a=" << (std::string)a << ", b=" << (std::string)b;
        c = a;
        c ^= b;
        EXPECT_TRUE(c == (a ^ b)) << "a=" << (std::string)a << ", b=" << (std::string)b;

        // the compound shifts have to match the binary ones
        const int32_t shift = static_cast<int32_t>(get_uint32_random() % 128);
        c = a;
        c <<= shift;
        EXPECT_TRUE(c == (a << shift)) << "a=" << (std::string)a << ", shift=" << shift;
        c = a;
        c >>= shift;
        EXPECT_TRUE(c == (a >> shift)) << "a=" << (std::string)a << ", shift=" << shift;
    }

    // a zero or negative shift count leaves the value alone, as documented
    volatile int32_t shift = 0;
    const int128_t x(0x0123456789ABCDEFull, 0xFEDCBA9876543210ull);
    EXPECT_TRUE((x >> shift) == x);
    EXPECT_TRUE((x << shift) == x);
    shift = -1;
    EXPECT_TRUE((x >> shift) == x);
    EXPECT_TRUE((x << shift) == x);

    // the right shift is arithmetic: it replicates the sign bit rather than shifting in zeros
    EXPECT_TRUE((int128_t(-1ll) >> 1) == int128_t(-1ll));
    EXPECT_TRUE((int128_t(-2ll) >> 1) == int128_t(-1ll));
    EXPECT_TRUE((int128_t(0ull, 0x8000000000000000ull) >> 127) == int128_t(-1ll));
    EXPECT_TRUE((int128_t(0ull, 0x4000000000000000ull) >> 126) == int128_t(1ll));
}

// A scalar on the left used to be ambiguous: converting it to int128_t and converting the int128_t
// to a builtin type are both a single user defined conversion.
TEST(int128_t, ScalarOnTheLeftHandSide)
{
    const int128_t x(10);
    EXPECT_TRUE((1 + x) == int128_t(11));
    EXPECT_TRUE((100 - x) == int128_t(90));
    EXPECT_TRUE((3 * x) == int128_t(30));
    EXPECT_TRUE((100 / x) == int128_t(10));
    EXPECT_TRUE((105 % x) == int128_t(5));
    EXPECT_TRUE((0xFF & x) == int128_t(10));
    EXPECT_TRUE((1 | x) == int128_t(11));
    EXPECT_TRUE((3 ^ x) == int128_t(9));

    // subtracting past zero must keep the sign rather than wrap
    EXPECT_TRUE((1 - x) == int128_t(-9));
    EXPECT_TRUE((-1 * x) == int128_t(-10));

    // The left operand is widened rather than the object being narrowed, so the result is 128 bit.
    static_assert(std::is_same_v<decltype(1 + x), int128_t>);
    static_assert(std::is_same_v<decltype(1 ^ x), int128_t>);
}
// The comparisons only existed with the int128_t on the left. A scalar there had nothing but the
// builtin comparisons to choose from, and the conversion operators of int128_t make every one of
// them equally good, so the call was ambiguous. Signedness is the point here: resolving to any of
// the unsigned builtin conversions would order the negative values above the positive ones.
TEST(int128_t, ScalarOnTheLeftHandSideComparisons)
{
    const int128_t x(10);
    EXPECT_TRUE(10 == x);
    EXPECT_TRUE(9 != x);
    EXPECT_TRUE(9 < x);
    EXPECT_FALSE(10 < x);
    EXPECT_TRUE(9 <= x);
    EXPECT_TRUE(10 <= x);
    EXPECT_TRUE(11 > x);
    EXPECT_FALSE(10 > x);
    EXPECT_TRUE(11 >= x);
    EXPECT_TRUE(10 >= x);

    // across zero, where a narrowed unsigned comparison would give the opposite answer
    const int128_t negative(-3);
    EXPECT_TRUE(-3 == negative);
    EXPECT_TRUE(-4 < negative);
    EXPECT_TRUE(-2 > negative);
    EXPECT_TRUE(0 > negative);
    EXPECT_TRUE(1 > negative);
    EXPECT_FALSE(0 < negative);
    EXPECT_TRUE(0 > int128_t(-1));
    EXPECT_TRUE(-1 < int128_t(0));

    // a value no builtin type can hold must not compare through a narrowed conversion
    const int128_t big = int128_t(1) << 100;
    EXPECT_TRUE(1 < big);
    EXPECT_FALSE(1 > big);
    EXPECT_TRUE(1 != big);
    EXPECT_TRUE(1 > -big);

    // the int128_t on the left forms must still resolve
    EXPECT_TRUE(x == 10);
    EXPECT_TRUE(x > 9);
    EXPECT_TRUE(x == int128_t(10));
}
// Widening the left operand is what makes (1 << 100) produce the expected power of two. The same
// expression on a builtin int is undefined behavior, the shift count being wider than the operand.
TEST(int128_t, ScalarOnTheLeftHandSideShifts)
{
    const int128_t ten(10);
    EXPECT_TRUE((1 << ten) == int128_t(1024));
    EXPECT_TRUE((1024 >> ten) == int128_t(1));
    EXPECT_TRUE((3 << int128_t(2)) == int128_t(12));

    const int128_t hundred(100);
    EXPECT_TRUE((1 << hundred) == (int128_t(1) << 100));
    EXPECT_TRUE((1 << hundred) != 0);
    static_assert(std::is_same_v<decltype(1 << hundred), int128_t>);

    // the right shift is arithmetic, as it is for the int128_t on the left form
    EXPECT_TRUE((-1024 >> ten) == int128_t(-1));
    EXPECT_TRUE((-8 >> int128_t(1)) == int128_t(-4));

    // the int128_t on the left forms must still resolve
    EXPECT_TRUE((ten << 1) == int128_t(20));
    EXPECT_TRUE((ten >> 1) == int128_t(5));
}

/**********************************************************************
 * Compile time (constexpr) evaluation
 *
 * Each test below pairs static_asserts, which fail the build the moment one of these operations
 * stops being usable in a constant expression, with a runtime check of the same expression built
 * from opaque() values the optimizer cannot fold. See the matching block in
 * fixed_point128_gtest.cpp for why the runtime half is not a restatement of the compile time one.
 *
 * The two instantiations share most of their code, so these concentrate on what the signedness
 * changes: sign extension, the arithmetic right shift, negation and the signed comparisons.
 * uint128_t_gtest.cpp covers the shared operations.
 ***********************************************************************/
TEST(int128_t, ConstexprConstructionAndConversion)
{
    constexpr int128_t zero;
    constexpr int128_t negative(-7);
    constexpr int128_t positive(7);
    constexpr int128_t copied(negative);
    constexpr int128_t assigned = [] { int128_t a; a = -5; return a; }();
    constexpr int128_t from_char(static_cast<char>(-3));
    constexpr int128_t most_negative(0ull, 0x8000000000000000ull);

    static_assert(zero.is_zero());
    static_assert(static_cast<int64_t>(negative) == -7ll);
    static_assert(static_cast<int64_t>(assigned) == -5ll);
    static_assert(static_cast<int64_t>(from_char) == -3ll);
    static_assert(copied == negative);

    // A negative value is sign extended into the high QWORD, a non negative one is zero extended.
    // The values are constructed inside the lambdas: a captureless lambda cannot odr-use a local,
    // and passing an object to get_components() by reference is exactly that.
    constexpr uint64_t neg_high = [] {
        uint64_t l = 0, h = 0;
        int128_t(-7).get_components(l, h);
        return h;
    }();
    constexpr uint64_t pos_high = [] {
        uint64_t l = 0, h = 0;
        int128_t(7).get_components(l, h);
        return h;
    }();
    static_assert(neg_high == UINT64_MAX);
    static_assert(pos_high == 0ull);
    static_assert(int128_t(-1).get_bit(127) == 1);  // the sign bit reaches the top of the high QWORD
    static_assert(static_cast<uint64_t>(int128_t(1ull, 2ull)) == 1ull);

    static_assert(negative.is_negative() && !negative.is_positive());
    static_assert(positive.is_positive() && !positive.is_negative());
    static_assert(zero.is_positive());  // zero counts as positive
    static_assert(most_negative.is_negative());
    static_assert(int128_t(~0ull, ~0ull) == int128_t(-1));  // the all ones pattern is -1

    EXPECT_TRUE(int128_t(static_cast<int64_t>(opaque(-7ll))) == negative);
}
TEST(int128_t, ConstexprShiftsAndUnary)
{
    constexpr int128_t minus_one(-1);
    constexpr int128_t most_negative(0ull, 0x8000000000000000ull);

    // the right shift is arithmetic: it replicates the sign bit rather than shifting in zeros
    static_assert((minus_one >> 1) == minus_one);
    static_assert((minus_one >> 200) == minus_one);  // shifting past the width leaves all sign bits
    static_assert((int128_t(-16) >> 2) == int128_t(-4));
    static_assert((int128_t(-2) >> 1) == minus_one);
    static_assert((most_negative >> 127) == minus_one);
    static_assert((int128_t(0ull, 0x4000000000000000ull) >> 126) == int128_t(1));
    static_assert((int128_t(16) >> 2) == int128_t(4));  // a positive value shifts in zeros
    static_assert((int128_t(1) << 127) == most_negative);
    static_assert((int128_t(1) << 200).is_zero());

    static_assert(-int128_t(5) == int128_t(-5));
    static_assert(-int128_t(-5) == int128_t(5));
    static_assert(-int128_t(0) == int128_t(0));
    static_assert(+int128_t(-5) == int128_t(-5));
    static_assert(-most_negative == most_negative);  // its magnitude has no representation, like INT64_MIN
    static_assert(~int128_t(0) == minus_one);

    static_assert(abs(int128_t(-5)) == int128_t(5));
    static_assert(abs(int128_t(5)) == int128_t(5));
    static_assert(abs(int128_t(0)) == int128_t(0));

    EXPECT_TRUE((int128_t(static_cast<int64_t>(opaque(-1ll))) >> 200) == minus_one);
    EXPECT_TRUE((int128_t(static_cast<int64_t>(opaque(-16ll))) >> 2) == int128_t(-4));
    EXPECT_TRUE(abs(int128_t(static_cast<int64_t>(opaque(-5ll)))) == int128_t(5));
}
TEST(int128_t, ConstexprArithmeticAndComparisons)
{
    static_assert(static_cast<int64_t>(int128_t(3) + int128_t(4)) == 7ll);
    static_assert(static_cast<int64_t>(int128_t(3) - int128_t(4)) == -1ll);
    static_assert(static_cast<int64_t>(int128_t(-3) + int128_t(-4)) == -7ll);
    static_assert(int128_t(3) - int128_t(3) == int128_t(0));
    static_assert(static_cast<int64_t>(4 - int128_t(3)) == 1ll);  // scalar on the left

    static_assert(static_cast<int64_t>(int128_t(-6) * int128_t(7)) == -42ll);
    static_assert(static_cast<int64_t>(int128_t(-6) * int128_t(-7)) == 42ll);
    static_assert(static_cast<int64_t>(int128_t(-6) * -7) == 42ll);   // generic right hand side
    static_assert(static_cast<int64_t>(-7 * int128_t(-6)) == 42ll);   // scalar on the left
    static_assert(sqr(int128_t(-9)) == int128_t(81));
    static_assert(sqr(int128_t(-9)) == int128_t(-9) * int128_t(-9));  // sqr is bit identical to x * x
    static_assert(pow(int128_t(-3), 3u) == int128_t(-27));
    static_assert(pow(int128_t(-3), 2u) == int128_t(9));

    constexpr int128_t stepped = [] { int128_t a(5); ++a; a++; --a; return a; }();
    static_assert(static_cast<int64_t>(stepped) == 6ll);
    constexpr int128_t crossed_zero = [] { int128_t a(1); --a; --a; return a; }();
    static_assert(crossed_zero == int128_t(-1));

    // the comparisons are signed: a negative value is smaller whatever the magnitude of its bits
    static_assert(int128_t(-2) < int128_t(-1));
    static_assert(int128_t(-1) < int128_t(1));
    static_assert(int128_t(1) > int128_t(-1));
    static_assert(int128_t(0ull, 0x8000000000000000ull) < int128_t(0));  // the most negative value
    static_assert(int128_t(3ull, 0ull) < int128_t(1ull, 2ull));  // two positives, decided by the high QWORD
    static_assert(int128_t(-1) <= -1 && int128_t(-1) >= -1);
    static_assert(int128_t(-1) == -1 && int128_t(-1) != 1);

    static_assert(lzcnt128(int128_t(1)) == 127);
    static_assert(log2(int128_t(8)) == 3);

    EXPECT_TRUE(int128_t(static_cast<int64_t>(opaque(-6ll))) * int128_t(7) == int128_t(-42));
    EXPECT_TRUE(sqr(int128_t(static_cast<int64_t>(opaque(-9ll)))) == int128_t(81));
    EXPECT_TRUE(pow(int128_t(static_cast<int64_t>(opaque(-3ll))), 3u) == int128_t(-27));
    EXPECT_TRUE(int128_t(static_cast<int64_t>(opaque(-2ll))) < int128_t(-1));
}

// The double and float conversions used to be the one thing that could not happen at compile time:
// they read the inactive member of a union, which constant evaluation rejects. Now that Double and
// Float convert with std::bit_cast, a floating point literal crosses the boundary either way.
TEST(int128_t, ConstexprFloatingPointConversion)
{
    constexpr int128_t positive = 42.0;
    constexpr int128_t negative = -42.0;

    static_assert(static_cast<int64_t>(positive) == 42ll);
    static_assert(static_cast<int64_t>(negative) == -42ll);
    static_assert(static_cast<double>(positive) == 42.0);
    static_assert(static_cast<double>(negative) == -42.0);
    static_assert(static_cast<float>(negative) == -42.0f);
    static_assert(int128_t(0.0).is_zero());
    static_assert(negative.is_negative() && positive.is_positive());
    static_assert(int128_t(-3.99) == int128_t(-3));  // truncates towards zero
    static_assert(static_cast<double>(int128_t(-1e18)) == -1e18);

    EXPECT_TRUE(int128_t(static_cast<double>(opaque(-42ll))) == negative);
    EXPECT_DOUBLE_EQ(static_cast<double>(int128_t(opaque(-42ll))), -42.0);
}
