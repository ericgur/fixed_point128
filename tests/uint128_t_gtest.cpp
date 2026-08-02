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
 * uint128_t tests
 ***********************************************************************/
// Construct fixed_point128 and convert back to/from various elements.
TEST(uint128_t, DefaultConstructor)
{
    uint128_t i;
    EXPECT_EQ(static_cast<uint64_t>(i), 0ull);
}
TEST(uint128_t, ConstructorFromDouble)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = fabs(floor(get_double_random(1, 127)));
        uint128_t f = value;
        if (check_overflow_uint128(value)) {
            continue;
        }
        double f_value = static_cast<double>(f);
        EXPECT_DOUBLE_EQ(f_value, value) << "value=" << value;
    }
}
TEST(uint128_t, ConstructorFromFloat)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        float value = (float)fabs(floor(get_double_random(1, 127)));
        uint128_t f = value;
        if (check_overflow_uint128(value)) {
            continue;
        }
        EXPECT_FLOAT_EQ(static_cast<float>(f), value) << "value=" << value;
    }
}
TEST(uint128_t, ConstructorFromInt32)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        int32_t value = get_int32_random();
        uint128_t f = value;
        EXPECT_EQ(static_cast<int32_t>(f), value);
    }
}
TEST(uint128_t, ConstructorFromUnsignedInt32)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        uint32_t value = get_uint32_random();
        uint128_t f = value;
        EXPECT_EQ(static_cast<uint32_t>(f), value);
    }
}
TEST(uint128_t, ConstructorFromInt64)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        int64_t value = get_int64_random();
        uint128_t f = value;
        EXPECT_EQ(static_cast<int64_t>(f), value);
    }
}
TEST(uint128_t, ConstructorFromUnsignedInt64)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        uint64_t value = get_uint64_random();
        uint128_t f = value;
        EXPECT_EQ(static_cast<uint64_t>(f), value);
    }
}
TEST(uint128_t, ConstructorFromString)
{
    const char* values[] = {"0xDEADBEAF", "1234", "123456789012345678"};
    for (auto i = 0u; i < array_length(values); ++i) {
        uint128_t f = values[i];
        double d1 = strtod(values[i], nullptr);
        double d2 = strtod(static_cast<char*>(f), nullptr);
        EXPECT_DOUBLE_EQ(d1, d2);
    }
}
TEST(uint128_t, CopyConstructor)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random(1, 127);
        uint128_t f1 = value;
        uint128_t f2 = f1;
        EXPECT_DOUBLE_EQ(static_cast<double>(f1), static_cast<double>(f2));
    }
}
TEST(uint128_t, MoveConstructor)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random(1, 127);
        uint128_t f1 = value;
        uint128_t f2(std::move(f1));
        EXPECT_DOUBLE_EQ(static_cast<double>(f1), static_cast<double>(f2));
    }
}
TEST(uint128_t, AssignmentOperator)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random(1, 127);
        uint128_t f1 = value;
        uint128_t f2;
        f2 = f1;
        EXPECT_DOUBLE_EQ(static_cast<double>(f1), static_cast<double>(f2));
    }
}
TEST(uint128_t, MoveAssignmentOperator)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random(1, 127);
        uint128_t f1 = value;
        uint128_t f2;
        f2 = std::move(f1);
        EXPECT_DOUBLE_EQ(static_cast<double>(f1), static_cast<double>(f2));
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
TEST(uint128_t, AddDifferentSign)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        uint64_t value1 = get_uint64_random();
        // two's complement negation written as a subtraction from zero: same bit pattern as unary
        // minus, without MSVC's C4146 complaint about negating an unsigned type
        uint64_t value2 = 0ull - get_uint64_random();
        uint64_t res = value1 + value2;
        uint128_t f1 = value1;
        uint128_t f2 = value2;
        uint128_t f3 = f1 + f2;
        uint64_t uint128_res = static_cast<uint64_t>(f3);
        EXPECT_EQ(uint128_res, res) << "value1=" << value1 << ", value2=" << value2;
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
TEST(uint128_t, AddUnsignedInt64)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        auto value1 = get_uint64_random();
        auto value2 = get_uint64_random();
        auto res = value1 + value2;
        uint128_t f1 = value1;
        uint128_t f3 = f1 + value2;
        uint64_t uint128_res = static_cast<uint64_t>(f3);
        EXPECT_EQ(uint128_res, res) << "value1=" << value1 << ", value2=" << value2;
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
TEST(uint128_t, SubtractDifferentSign)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        uint64_t value1 = get_uint64_random();
        // two's complement negation written as a subtraction from zero: same bit pattern as unary
        // minus, without MSVC's C4146 complaint about negating an unsigned type
        uint64_t value2 = 0ull - get_uint64_random();
        uint64_t res = value1 - value2;
        uint128_t f1 = value1;
        uint128_t f2 = value2;
        uint128_t f3 = f1 - f2;
        uint64_t uint128_res = static_cast<uint64_t>(f3);
        EXPECT_EQ(uint128_res, res) << "value1=" << value1 << ", value2=" << value2;
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
TEST(uint128_t, SubtractUnsignedInt64)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        auto value1 = get_uint64_random();
        auto value2 = get_uint64_random();
        auto res = value1 - value2;
        uint128_t f1 = value1;
        uint128_t f3 = f1 - value2;
        uint64_t uint128_res = static_cast<uint64_t>(f3);
        EXPECT_EQ(uint128_res, res) << "value1=" << value1 << ", value2=" << value2;
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
// argument: q*b + r == a with r < b identifies q and r uniquely. It is checked in full 128 bit
// arithmetic, and q*b + r cannot overflow because q*b <= a.
TEST(uint128_t, DivideByUint128)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const uint128_t a(get_uint64_random(), get_uint64_random());

        // The divisor's magnitude is spread over the whole range so that all three paths in
        // operator/= are reached: the power of two shift, the 64 bit divisor shortcut when the
        // high QWORD shifts out entirely, and the full 128 bit long division.
        // the shift counts are int32_t: an unsigned one leaves the shift operators ambiguous
        // against the builtin shifts reachable through operator uint64_t()
        const int32_t shift = static_cast<int32_t>(get_uint32_random() % 128);
        uint128_t b = uint128_t(get_uint64_random(), get_uint64_random()) >> shift;
        if ((get_uint32_random() & 7) == 0)
            b = uint128_t(1ull) << static_cast<int32_t>(get_uint32_random() % 128);  // exercise the shift path
        if (b.is_zero())
            continue;

        const uint128_t q = a / b;
        const uint128_t r = a % b;

        EXPECT_TRUE(r < b) << "a=" << (std::string)a << ", b=" << (std::string)b;
        EXPECT_TRUE(q * b + r == a) << "a=" << (std::string)a << ", b=" << (std::string)b;

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
        uint64_t uint128_res = static_cast<uint64_t>(f3);
        EXPECT_EQ(uint128_res, res) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(uint128_t, CompareUint128)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = fabs(floor(get_double_random(1, 127)));
        double value2 = fabs(floor(get_double_random(1, 127)));
        bool res = value1 > value2;
        uint128_t f1 = value1;
        uint128_t f2 = value2;
        if (check_overflow_uint128(value1) || check_overflow_uint128(value2))
            continue;

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
TEST(uint128_t, OperatorPlusPlus)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        uint64_t value1 = get_uint64_random();
        uint64_t res = value1 + 1;
        uint128_t f1 = value1;
        uint128_t f2 = f1;
        f2++;

        EXPECT_EQ(static_cast<uint64_t>(f2), res) << "operator++(int)"
                                                  << "value1=" << value1;

        f2 = f1;
        ++f2;

        EXPECT_EQ(static_cast<uint64_t>(f2), res) << "operator++()"
                                                  << "value1=" << value1;
    }
}
TEST(uint128_t, OperatorMinusMinus)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        uint64_t value1 = get_uint64_random();
        uint64_t res = value1 - 1;
        uint128_t f1 = value1;
        uint128_t f2 = f1;
        f2--;

        EXPECT_EQ(static_cast<uint64_t>(f2), res) << "operator++(int)"
                                                  << "value1=" << value1;

        f2 = f1;
        --f2;

        EXPECT_EQ(static_cast<uint64_t>(f2), res) << "operator++()"
                                                  << "value1=" << value1;
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
TEST(uint128_t, log10)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = floor(1.0 + fabs(get_double_random(1, 127)));
        uint64_t res = (uint64_t)floor(log10(value1));  // double doesn't lose any bits with this operation!
        if (check_overflow_uint128(value1)) {
            continue;
        }
        uint128_t i1 = value1;
        uint64_t i_res = log10(i1);
        EXPECT_EQ(i_res, res) << "double value1=" << value1;
    }
}
TEST(uint128_t, log)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = 1.0 + floor(fabs(get_double_random(0, 127)));
        uint64_t res = (uint64_t)floor(log(value1));  // double doesn't lose any bits with this operation!
        if (check_overflow_uint128(value1)) {
            continue;
        }
        uint128_t i1 = value1;
        uint64_t i_res = log(i1);
        EXPECT_EQ(i_res, res) << "double value1=" << value1;
    }
}
TEST(uint128_t, sqrt)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = floor(1.0 + fabs(get_double_random(1, 127)));
        if (check_overflow_uint128(value1)) {
            continue;
        }
        uint128_t i1 = value1;
        uint64_t i_res = sqrt(i1);
        // floor(sqrt(value1)) is not a usable reference at this range, the root needs more than the
        // 53 bits a double holds. is_exact_isqrt() checks r^2 <= x < (r+1)^2 in 128 bit instead.
        EXPECT_TRUE(is_exact_isqrt(i_res, i1)) << "double value1=" << value1 << ", sqrt=" << i_res;
        // the conversion of the operand has to be exact, otherwise the check above passes vacuously
        EXPECT_EQ(static_cast<double>(i1), value1);
    }
}
TEST(uint128_t, pow)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = get_uint32_random() % 8;
        double value2 = get_uint32_random() % 16;
        double res = pow(value1, value2);  // double doesn't lose any bits with this operation!
        if (check_overflow_uint128(value1) || check_overflow_uint128(res)) {
            continue;
        }
        uint128_t i1 = value1;
        uint128_t i2 = pow(i1, (uint32_t)value2);
        uint128_t uint128_res = static_cast<double>(i2);
        EXPECT_DOUBLE_EQ(uint128_res, res) << "value1=" << value1 << ", value2=" << value2;
    }
}

/**********************************************************************
 * uint128_t regression tests
 *
 * Each test below pins down a defect that was fixed. They deliberately use
 * operand ranges the original tests never reached.
 ***********************************************************************/

// operator%= used to hand div_32bit this object's own QWORDs as the remainder buffer. div_32bit
// shrinks the denominator past its leading zero words and fills only that many words, so a divisor
// whose high QWORD fits in 32 bits left the top 32 bits of the numerator untouched in the result.
// ModuloByUint128 misses this: it only uses a 64 bit numerator and a 32 bit divisor.
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
        EXPECT_TRUE(r < b) << "a=" << (std::string)a << ", b=" << (std::string)b;
        EXPECT_TRUE(q * b + r == a) << "a=" << (std::string)a << ", b=" << (std::string)b;
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
    // 2^65-1 is a tie that rounds up to 2^65
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

    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const int64_t value = get_int64_random();
        // a double holds any int64_t exactly only up to 2^53, keep the comparison meaningful
        if (value > (1ll << 53) || value < -(1ll << 53))
            continue;
        EXPECT_TRUE(uint128_t(static_cast<double>(value)) == uint128_t(value)) << "value=" << value;
    }
}
// The templated operators cast the operand straight to uint64_t, which is undefined behavior for
// a floating point value at or above 2^64 and turned a divisor between zero and one into a bogus
// divide by zero. Floating point operands are now routed through the uint128_t constructor.
TEST(uint128_t, MultiplyAndDivideByFloatingPoint)
{
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
// template, and the two types still refuse to convert into one another. Trivial copyability is new,
// a consequence of defaulting the copy and move members.
TEST(uint128_t, LayoutAndTraits)
{
    static_assert(std::is_same_v<uint128_t, int128_base<false>>, "uint128_t is the unsigned instantiation");
    static_assert(uint128_t::is_signed == false, "uint128_t reports itself as unsigned");
    static_assert(sizeof(uint128_t) == 16, "the object holds exactly two QWORDs");
    static_assert(alignof(uint128_t) == 16, "FP128_ALIGN16 has to survive the template");
    static_assert(std::is_standard_layout_v<uint128_t>, "the division helpers alias the members as an array");
    static_assert(std::is_trivially_copyable_v<uint128_t>, "the copy and move members are defaulted");
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
