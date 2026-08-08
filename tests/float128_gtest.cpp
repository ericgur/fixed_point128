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
#include <format>
#include "gtest_shared.h"

/**********************************************************************
 * fixed_point128 tests
 ***********************************************************************/
// Construct fixed_point128 and convert back to/from various elements.
TEST(float128, DefaultConstructor)
{
    float128 f;
    EXPECT_DOUBLE_EQ(static_cast<double>(f), 0.0);
}
TEST(float128, ConstructorFromDouble)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random(-1023, 1022);
        float128 f = value;
        double f_value = static_cast<double>(f);
        EXPECT_DOUBLE_EQ(f_value, value) << "value=" << value;
    }
}
TEST(float128, ConstructorFromFloat)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        float value = (float)get_double_random(-1023, 1022);
        float128 f = value;
        float f_value = static_cast<float>(f);
        EXPECT_FLOAT_EQ(f_value, value) << "value=" << value;
    }
}
TEST(float128, ConstructorFromInt32)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        int32_t value = get_int32_random();
        float128 f = value;
        EXPECT_EQ(static_cast<int32_t>(f), value);
    }
}
TEST(float128, ConstructorFromUnsignedInt32)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        uint32_t value = get_uint32_random();
        float128 f = value;
        EXPECT_EQ(static_cast<uint32_t>(f), value);
    }
}
TEST(float128, ConstructorFromInt64)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        int64_t value = get_int64_random();
        float128 f = value;
        EXPECT_EQ(static_cast<int64_t>(f), value);
    }
}
TEST(float128, ConstructorFromUnsignedInt64)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        uint64_t value = get_uint64_random();
        float128 f = value;
        EXPECT_EQ(static_cast<uint64_t>(f), value);
    }
}
TEST(float128, ConstructorFromString)
{
    char str[128];
    char str2[128];
    srand(RANDOM_SEED);
    int j = 0;
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        try {
            // construct random number string
            str[0] = get_digit_random();
            str[1] = '.';
            for (j = 2; j < 36; ++j) {
                str[j] = get_digit_random();
            }
            // fraction length is based on the value of the first digit
            switch (str[0]) {
            case 0:
                break;
            case 1:
                j -= 1;
                break;
            case 2:
            case 3:
                j -= 2;
                break;
            case 4:
            case 5:
            case 6:
            case 7:
                j -= 3;
                break;
            default:
                j -= 4;
                break;
            }
            str[j] = '\0';
            strcpy(str2, str);
            const float128 f = str;

            // The parse is exact, so the value has to agree with what the C library makes of the
            // same text to every bit a double holds.
            const double reference = strtod(str, nullptr);
            EXPECT_NEAR(static_cast<double>(f), reference, fabs(reference) * DBL_EPSILON * 2)
                << "parsed " << str << " as " << static_cast<double>(f);

            // And the shortest text the value produces has to read back as the same value. This is
            // what the previous check was reaching for by comparing the last digits of the two
            // strings, which only worked while both were laid out the same way: the shortest form
            // switches to scientific notation once that is shorter, and the comparison then read
            // an exponent as if it were digits.
            const std::string text = to_string(f);
            float128 round_tripped;
            fp128::from_chars(text.data(), text.data() + text.size(), round_tripped);
            EXPECT_TRUE(round_tripped == f) << "round trip of " << str << " through " << text;
        } catch (...) {
            // EXPECT_NO_THROW(i) used to sit here, which would have passed even on this path:
            // evaluating an integer cannot throw, so reaching the handler reported nothing.
            ADD_FAILURE() << "unexpected exception at iteration " << i;
        }
    }
}
TEST(float128, CopyConstructor)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random();
        float128 f1 = value;
        float128 f2 = f1;
        EXPECT_DOUBLE_EQ(static_cast<double>(f1), static_cast<double>(f2)) << "value=" << value;
        ;
    }
}
TEST(float128, MoveConstructor)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random();
        float128 f1 = value;
        float128 f2(std::move(f1));
        EXPECT_DOUBLE_EQ(static_cast<double>(f1), static_cast<double>(f2)) << "value=" << value;
        ;
    }
}
TEST(float128, AssignmentOperator)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random();
        float128 f1 = value;
        float128 f2;
        f2 = f1;
        EXPECT_DOUBLE_EQ(static_cast<double>(f1), static_cast<double>(f2)) << "value=" << value;
        ;
    }
}
TEST(float128, MoveAssignmentOperator)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random(-1023, 1022);
        float128 f1 = value;
        float128 f2;
        f2 = std::move(f1);
        EXPECT_DOUBLE_EQ(static_cast<double>(f1), static_cast<double>(f2)) << "value=" << value;
        ;
    }
}
TEST(float128, AddSameSign)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = fabs(get_double_random());
        double value2 = value1 * 2.5;
        double res = value1 + value2;
        float128 f1 = value1;
        float128 f2 = value2;
        float128 f3 = f1 + f2;
        double float128_res = static_cast<double>(f3);
        EXPECT_DOUBLE_EQ(float128_res, res) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(float128, AddDifferentSign)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = fabs(get_double_random());
        double value2 = value1 * -2.5;
        double res = value1 + value2;
        float128 f1 = value1;
        float128 f2 = value2;
        float128 f3 = f1 + f2;
        double float128_res = static_cast<double>(f3);
        EXPECT_DOUBLE_EQ(float128_res, res) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(float128, AddDouble)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = get_double_random();
        double value2 = get_double_random();
        double res = value1 + value2;
        float128 f1 = value1;
        float128 f3 = f1 + value2;
        double float128_res = static_cast<double>(f3);
        EXPECT_DOUBLE_EQ(float128_res, res) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(float128, AddFloat)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        float value1 = (float)get_double_random();
        float value2 = (float)get_double_random();
        float res = value1 + value2;
        float128 f1 = value1;
        float128 f3 = f1 + value2;
        float float128_res = static_cast<float>(f3);
        EXPECT_FLOAT_EQ(float128_res, res) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(float128, AddInt32)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const int32_t value1 = abs(get_int32_random());
        const int32_t value2 = abs(get_int32_random());
        // The sum of two positive 32 bit values can exceed the type. Detecting that after the
        // fact by testing for a negative result relies on signed overflow, which is undefined
        // behavior: the compiler may conclude the test can never be true and drop it, which is
        // exactly what Clang does. Compute the reference in a wider type instead.
        const int64_t res = static_cast<int64_t>(value1) + value2;

        float128 f1 = value1;
        float128 f3 = f1 + value2;
        auto float128_res = static_cast<int64_t>(f3);
        EXPECT_EQ(float128_res, res) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(float128, AddUnsignedInt32)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        auto value1 = get_uint32_random();
        auto value2 = get_uint32_random();
        auto res = value1 + value2;
        if (res < value1 || res < value2)
            continue;  // wrap around

        float128 f1 = value1;
        float128 f3 = f1 + value2;
        auto float128_res = static_cast<uint64_t>(f3);
        EXPECT_EQ(float128_res, res) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(float128, AddInt64)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const int64_t value1 = abs(get_int64_random());
        const int64_t value2 = abs(get_int64_random());
        // There is no wider type to compute the reference in, so test for the overflow before it
        // happens. Inspecting the result afterwards would rely on signed overflow, which is
        // undefined behavior and gets optimized away.
        if (value1 > INT64_MAX - value2)
            continue;
        const int64_t res = value1 + value2;
        float128 f1 = value1;
        float128 f3 = f1 + value2;
        auto float128_res = static_cast<int64_t>(f3);
        EXPECT_EQ(float128_res, res) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(float128, AddUnsignedInt64)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        auto value1 = get_uint64_random();
        auto value2 = get_uint64_random();
        auto res = value1 + value2;
        if (res < value1 || res < value2)
            continue;  // wrap around
        float128 f1 = value1;
        float128 f3 = f1 + value2;
        auto float128_res = static_cast<uint64_t>(f3);
        EXPECT_EQ(float128_res, res) << "value1=" << value1 << ", value2=" << value2;
    }
}
// IEEE 754 special value rules for addition: a NaN operand propagates, adding infinities of
// opposite signs is the invalid operation and produces a NaN, and any other infinite operand
// carries its own sign through to the result.
TEST(float128, AddSpecialValues)
{
    const float128 inf = float128::inf();
    const float128 nan = float128::nan();
    const float128 zero = float128(0);
    const float128 three = float128(3);

    // a finite operand does not change an infinity
    EXPECT_TRUE(inf + three == inf);
    EXPECT_TRUE(three + inf == inf);
    EXPECT_TRUE(-inf + three == -inf);
    EXPECT_TRUE(three + -inf == -inf);
    EXPECT_TRUE(inf + -three == inf);
    EXPECT_TRUE(-inf + -three == -inf);
    EXPECT_TRUE(inf + zero == inf);
    EXPECT_TRUE(-inf + zero == -inf);

    // like signed infinities add to the same infinity
    EXPECT_TRUE(inf + inf == inf);
    EXPECT_TRUE(-inf + -inf == -inf);

    // opposite signed infinities are the invalid operation
    EXPECT_TRUE(isnan(inf + -inf));
    EXPECT_TRUE(isnan(-inf + inf));
    EXPECT_TRUE(isnan(inf - inf));
    EXPECT_TRUE(isnan(-inf - -inf));

    // subtraction goes through operator+= and follows the same rules
    EXPECT_TRUE(inf - three == inf);
    EXPECT_TRUE(three - inf == -inf);
    EXPECT_TRUE(inf - -inf == inf);

    // a NaN operand always wins, including against an infinity
    EXPECT_TRUE(isnan(nan + three));
    EXPECT_TRUE(isnan(three + nan));
    EXPECT_TRUE(isnan(nan + inf));
    EXPECT_TRUE(isnan(inf + nan));
    EXPECT_TRUE(isnan(nan + -inf));
    EXPECT_TRUE(isnan(-inf + nan));
    EXPECT_TRUE(isnan(nan + nan));
    EXPECT_TRUE(isnan(nan - inf));

    // finite operands are untouched by the special value handling
    EXPECT_TRUE(three + three == float128(6));
    EXPECT_TRUE(three - three == zero);
}
TEST(float128, SubtractSameSign)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = fabs(get_double_random());
        double value2 = value1 * 2.5;
        double res = value1 - value2;
        float128 f1 = value1;
        float128 f2 = value2;
        float128 f3 = f1 - f2;
        double float128_res = static_cast<double>(f3);
        EXPECT_DOUBLE_EQ(float128_res, res) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(float128, SubtractDifferentSign)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = fabs(get_double_random());
        double value2 = value1 * -2.5;
        double res = value1 - value2;
        float128 f1 = value1;
        float128 f2 = value2;
        float128 f3 = f1 - f2;
        double float128_res = static_cast<double>(f3);
        EXPECT_DOUBLE_EQ(float128_res, res) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(float128, SubtractInt32)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        auto value1 = abs(get_int32_random());
        auto value2 = abs(get_int32_random());
        if (value2 > value1)
            value2 = value1 / 3;

        auto res = value1 - value2;

        float128 f1 = value1;
        float128 f3 = f1 - value2;
        auto float128_res = static_cast<int32_t>(f3);
        EXPECT_EQ(float128_res, res) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(float128, SubtractUnsignedInt32)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        auto value1 = get_uint32_random();
        auto value2 = get_uint32_random();
        if (value1 < value2) {
            std::swap(value1, value2);
        }
        auto res = value1 - value2;
        float128 f1 = value1;
        float128 f3 = f1 - value2;
        auto float128_res = static_cast<uint32_t>(f3);
        EXPECT_EQ(float128_res, res) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(float128, SubtractInt64)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        auto value1 = get_int64_random();
        auto value2 = get_int64_random();
        if ((value1 < 0 && value2 > 0) || (value1 > 0 && value2 < 0))
            value2 = -value2;
        auto res = value1 - value2;
        float128 f1 = value1;
        float128 f3 = f1 - value2;
        auto float128_res = static_cast<int64_t>(f3);
        EXPECT_EQ(float128_res, res) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(float128, SubtractUnsignedInt64)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        auto value1 = get_uint64_random();
        auto value2 = get_uint64_random();
        if (value1 < value2) {
            std::swap(value1, value2);
        }
        auto res = value1 - value2;
        float128 f1 = value1;
        float128 f3 = f1 - value2;
        auto float128_res = static_cast<uint64_t>(f3);
        EXPECT_EQ(float128_res, res) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(float128, MultiplyByFloat128)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = get_double_random();
        double value2 = get_double_random();
        double res = value1 * value2;
        float128 f1 = value1;
        float128 f2 = value2;
        float128 f3 = f1 * f2;

        double float128_res = static_cast<double>(f3);
        EXPECT_DOUBLE_EQ(float128_res, res) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(float128, MultiplyByExponentOf2)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = get_double_random();
        double value2 = ::pow(2.0, get_int32_random() % 100);
        double res = value1 * value2;
        float128 f1 = value1;
        float128 f2 = value2;
        float128 f3 = f1 * f2;
        double float128_res = static_cast<double>(f3);
        EXPECT_DOUBLE_EQ(float128_res, res) << "value1=" << value1 << ", value2=" << value2;

        f3 = f2 * f1;
        float128_res = static_cast<double>(f3);
        EXPECT_DOUBLE_EQ(float128_res, res) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(float128, MultiplyByDouble)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        // The two operands used to be overwritten here with a hardcoded pair, left behind from
        // debugging, so all 65536 iterations multiplied 0.80270613357871456 by 0.70710678118654757.
        double value1 = get_double_random();
        double value2 = get_double_random();
        double res = value1 * value2;
        float128 f1 = value1;
        float128 f3 = f1 * value2;
        double float128_res = static_cast<double>(f3);
        EXPECT_DOUBLE_EQ(float128_res, res) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(float128, MultiplyByFloat)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        float value1 = (float)get_double_random();
        float value2 = (float)get_double_random();
        float res = value1 * value2;
        float128 f1 = value1;
        float128 f3 = f1 * value2;
        float float128_res = static_cast<float>(f3);
        EXPECT_FLOAT_EQ(float128_res, res) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(float128, MultiplyByInt32)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        // A full width random pair overflows 31 bits every time, so the guard below used to reject
        // all 65536 iterations and the test ran no assertion at all. Shifting by at least 17 bounds
        // each operand to 2^14, which keeps every product inside the range the result is cast back
        // to. It also keeps value1 away from INT32_MIN, whose abs() below would be undefined.
        auto value1 = get_int32_random() >> (17 + static_cast<int32_t>(get_uint32_random() % 15));
        auto value2 = get_int32_random() >> (17 + static_cast<int32_t>(get_uint32_random() % 15));
        // check result overflow
        if (log2(abs(value1)) + log2(abs(value2)) > 31)
            continue;
        int64_t res = (int64_t)value1 * (int64_t)value2;
        float128 f1 = value1;
        float128 f3 = f1 * value2;
        auto float128_res = static_cast<int32_t>(f3);
        EXPECT_EQ(float128_res, res) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(float128, MultiplyByUnsignedInt32)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        // bounded to 2^15 for the reason given in MultiplyByInt32
        auto value1 = get_uint32_random() >> (17 + static_cast<int32_t>(get_uint32_random() % 15));
        auto value2 = get_uint32_random() >> (17 + static_cast<int32_t>(get_uint32_random() % 15));
        // check result overflow
        if (log2(value1) + log2(value2) > 31)
            continue;
        uint64_t res = (uint64_t)value1 * (uint64_t)value2;
        float128 f1 = value1;
        float128 f3 = f1 * value2;
        auto float128_res = static_cast<uint32_t>(f3);
        EXPECT_EQ(float128_res, res) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(float128, MultiplyByInt64)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        // bounded to 2^30 for the reason given in MultiplyByInt32
        auto value1 = get_int64_random() >> (33 + static_cast<int32_t>(get_uint32_random() % 31));
        auto value2 = get_int64_random() >> (33 + static_cast<int32_t>(get_uint32_random() % 31));
        // check result overflow
        if (log2(abs(value1)) + log2(abs(value2)) > 63)
            continue;
        int64_t res = value1 * value2;
        float128 f1 = value1;
        float128 f3 = f1 * value2;
        auto float128_res = static_cast<int64_t>(f3);
        EXPECT_EQ(float128_res, res) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(float128, MultiplyByUnsignedInt64)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        // bounded to 2^31 for the reason given in MultiplyByInt32
        auto value1 = get_uint64_random() >> (33 + static_cast<int32_t>(get_uint32_random() % 31));
        auto value2 = get_uint64_random() >> (33 + static_cast<int32_t>(get_uint32_random() % 31));
        // check result overflow
        if (log2(value1) + log2(value2) > 63)
            continue;
        uint64_t res = value1 * value2;
        float128 f1 = value1;
        float128 f3 = f1 * value2;
        auto float128_res = static_cast<uint64_t>(f3);
        EXPECT_EQ(float128_res, res) << "value1=" << value1 << ", value2=" << value2;
    }
}
// sqr(x) must be bit identical to x * x for every finite value. operator== on float128 is a
// raw bit comparison, so this is an exact check.
TEST(float128, SqrMatchesMultiply)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const double value = get_double_random();
        const float128 x = value;
        EXPECT_TRUE(sqr(x) == x * x) << "value=" << value;
    }
}
TEST(float128, SqrExponentOf2)
{
    srand(RANDOM_SEED);
    // exponents of 2 take a dedicated code path that skips the fraction multiply
    for (auto i = 0; i < 1000; ++i) {
        const double value = ::pow(2.0, get_int32_random() % 200);
        const float128 x = value;
        EXPECT_TRUE(sqr(x) == x * x) << "value=" << value;
        EXPECT_TRUE(sqr(-x) == (-x) * (-x)) << "value=" << -value;
    }
}
TEST(float128, SqrSpecialValues)
{
    // zero
    EXPECT_TRUE(sqr(float128(0)) == float128(0) * float128(0));
    EXPECT_TRUE(sqr(float128(0)).is_zero());
    // NaN stays a NaN
    EXPECT_TRUE(isnan(sqr(float128::nan())));
    // a square is never negative, and this agrees with x * x for the infinities as well
    EXPECT_TRUE(sqr(float128::inf()) == float128::inf());
    EXPECT_TRUE(sqr(-float128::inf()) == float128::inf());
    EXPECT_TRUE(sqr(-float128::inf()) == (-float128::inf()) * (-float128::inf()));
    // negating the operand must not change the result
    EXPECT_TRUE(sqr(-float128::pi()) == sqr(float128::pi()));
    EXPECT_TRUE(sqr(-float128::e()) == float128::e() * float128::e());
}
// IEEE 754 special value rules for multiplication: a NaN operand propagates, inf * zero is
// the invalid operation and produces a NaN, and any other infinite operand produces an
// infinity whose sign is the exclusive or of the operand signs.
TEST(float128, MultiplySpecialValues)
{
    const float128 inf = float128::inf();
    const float128 zero = float128(0);
    const float128 three = float128(3);

    // sign of an infinite result is the combination of both operand signs
    EXPECT_TRUE(inf * three == inf);
    EXPECT_TRUE(inf * -three == -inf);
    EXPECT_TRUE(-inf * three == -inf);
    EXPECT_TRUE(-inf * -three == inf);
    EXPECT_TRUE(three * inf == inf);
    EXPECT_TRUE(-three * inf == -inf);
    EXPECT_TRUE(three * -inf == -inf);
    EXPECT_TRUE(-three * -inf == inf);
    EXPECT_TRUE(inf * inf == inf);
    EXPECT_TRUE(inf * -inf == -inf);
    EXPECT_TRUE(-inf * inf == -inf);
    EXPECT_TRUE(-inf * -inf == inf);

    // inf * zero is invalid and yields a NaN, whatever the signs are
    EXPECT_TRUE(isnan(inf * zero));
    EXPECT_TRUE(isnan(zero * inf));
    EXPECT_TRUE(isnan(-inf * zero));
    EXPECT_TRUE(isnan(zero * -inf));
    EXPECT_TRUE(isnan(inf * -zero));

    // a NaN operand always wins, including against an infinity
    const float128 nan = float128::nan();
    EXPECT_TRUE(isnan(nan * three));
    EXPECT_TRUE(isnan(three * nan));
    EXPECT_TRUE(isnan(nan * inf));
    EXPECT_TRUE(isnan(inf * nan));
    EXPECT_TRUE(isnan(nan * -inf));
    EXPECT_TRUE(isnan(-inf * nan));
    EXPECT_TRUE(isnan(nan * zero));
    EXPECT_TRUE(isnan(nan * nan));

    // finite operands are untouched by the special value handling
    EXPECT_TRUE(three * three == float128(9));
    EXPECT_TRUE(-three * three == float128(-9));
}
TEST(float128, DivideByFloat128)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = get_double_random();
        double value2 = get_double_random();
        if (value2 == 0)
            continue;
        double res = value1 / value2;
        float128 f1 = value1;
        float128 f2 = value2;
        float128 f3 = f1 / f2;
        double float128_res = static_cast<double>(f3);
        EXPECT_DOUBLE_EQ(float128_res, res) << "value1=" << value1 << ", value2=" << value2;
    }
}
// IEEE 754 special value rules for division: a NaN operand propagates, zero / zero and
// inf / inf are the invalid operations, an infinite dividend or a zero divisor produce an
// infinity, an infinite divisor or a zero dividend produce a zero, and every one of those
// results carries the exclusive or of the operand signs.
TEST(float128, DivideSpecialValues)
{
    const float128 inf = float128::inf();
    const float128 nan = float128::nan();
    const float128 zero = float128(0);
    const float128 three = float128(3);

    // division by zero produces an infinity with the combined sign
    EXPECT_TRUE(three / zero == inf);
    EXPECT_TRUE(-three / zero == -inf);
    EXPECT_TRUE(three / -zero == -inf);
    EXPECT_TRUE(-three / -zero == inf);

    // an infinite dividend produces an infinity with the combined sign
    EXPECT_TRUE(inf / three == inf);
    EXPECT_TRUE(inf / -three == -inf);
    EXPECT_TRUE(-inf / three == -inf);
    EXPECT_TRUE(-inf / -three == inf);
    EXPECT_TRUE(inf / zero == inf);
    EXPECT_TRUE(-inf / zero == -inf);

    // an infinite divisor produces a zero with the combined sign
    EXPECT_TRUE((three / inf).is_zero());
    EXPECT_TRUE((three / -inf).is_zero());
    EXPECT_TRUE(three / inf == zero);
    EXPECT_TRUE(three / -inf == -zero);
    EXPECT_TRUE(-three / inf == -zero);
    EXPECT_TRUE(-three / -inf == zero);

    // a zero dividend produces a zero with the combined sign
    EXPECT_TRUE(zero / three == zero);
    EXPECT_TRUE(zero / -three == -zero);
    EXPECT_TRUE(-zero / three == -zero);
    EXPECT_TRUE(-zero / -three == zero);
    EXPECT_TRUE(zero / inf == zero);
    EXPECT_TRUE(zero / -inf == -zero);

    // the two invalid operations
    EXPECT_TRUE(isnan(zero / zero));
    EXPECT_TRUE(isnan(-zero / zero));
    EXPECT_TRUE(isnan(inf / inf));
    EXPECT_TRUE(isnan(inf / -inf));
    EXPECT_TRUE(isnan(-inf / inf));
    EXPECT_TRUE(isnan(-inf / -inf));

    // a NaN operand always wins, including against an infinity and a zero
    EXPECT_TRUE(isnan(nan / three));
    EXPECT_TRUE(isnan(three / nan));
    EXPECT_TRUE(isnan(nan / inf));
    EXPECT_TRUE(isnan(inf / nan));
    EXPECT_TRUE(isnan(nan / zero));
    EXPECT_TRUE(isnan(zero / nan));
    EXPECT_TRUE(isnan(nan / nan));

    // finite operands are untouched by the special value handling
    EXPECT_TRUE(float128(9) / three == three);
    EXPECT_TRUE(float128(-9) / three == -three);
}
TEST(float128, DivideByExponentOf2)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = get_double_random();
        double value2 = ::pow(2.0, get_int32_random() % 100);
        double res = value1 / value2;
        float128 f1 = value1;
        float128 f2 = value2;
        float128 f3 = f1 / f2;
        double float128_res = static_cast<double>(f3);
        EXPECT_DOUBLE_EQ(float128_res, res) << "value1=" << value1 << ", value2=" << value2;

        // do this in reverse
        res = value2 / value1;
        f3 = f2 / f1;
        float128_res = static_cast<double>(f3);
        EXPECT_DOUBLE_EQ(float128_res, res) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(float128, DivideByDouble)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = get_double_random();
        double value2 = get_double_random();
        double res = value1 / value2;
        if (value2 == 0)
            continue;
        float128 f1 = value1;
        float128 f3 = f1 / value2;
        double float128_res = static_cast<double>(f3);
        EXPECT_DOUBLE_EQ(float128_res, res) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(float128, DivideByFloat)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        float value1 = (float)get_double_random();
        float value2 = (float)get_double_random();
        if (value2 == 0)
            continue;
        float res = value1 / value2;
        float128 f1 = value1;
        float128 f3 = f1 / value2;
        float float128_res = static_cast<float>(f3);
        EXPECT_FLOAT_EQ(float128_res, res) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(float128, DivideByInt32)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = get_double_random(-1, 39);
        int32_t value2 = get_int32_random();
        double res = value1 / value2;
        if (value2 == 0)
            continue;
        float128 f1 = value1;
        float128 f3 = f1 / value2;
        double float128_res = static_cast<double>(f3);
        EXPECT_DOUBLE_EQ(float128_res, res) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(float128, DivideByUnsignedInt32)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = get_double_random(-1, 39);
        uint32_t value2 = get_uint32_random();
        double res = value1 / value2;
        if (value2 == 0)
            continue;
        float128 f1 = value1;
        float128 f3 = f1 / value2;
        double float128_res = static_cast<double>(f3);
        EXPECT_DOUBLE_EQ(float128_res, res) << "value1=" << value1 << ", value2=" << value2;
    }
}

TEST(float128, DivideByInt64)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = get_double_random(-1, 39);
        int64_t value2 = get_int32_random();  // on purpose
        double res = value1 / value2;
        if (value2 == 0)
            continue;
        float128 f1 = value1;
        float128 f3 = f1 / value2;
        double float128_res = static_cast<double>(f3);
        EXPECT_DOUBLE_EQ(float128_res, res) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(float128, DivideByUnsignedInt64)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = get_double_random(-1, 39);
        uint64_t value2 = get_uint32_random();  // on purpose
        double res = value1 / value2;
        if (value2 == 0)
            continue;
        float128 f1 = value1;
        float128 f3 = f1 / value2;
        double float128_res = static_cast<double>(f3);
        EXPECT_DOUBLE_EQ(float128_res, res) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(float128, CompareFloat128)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = get_double_random(-1, 39);
        double value2 = get_double_random(-1, 39);
        bool res = value1 > value2;
        float128 f1 = value1;
        float128 f2 = value2;

        bool float128_res = f1 > f2;
        EXPECT_TRUE(float128_res == res) << "operator>: "
                                         << "value1=" << value1 << ", value2=" << value2;

        res = value1 >= value2;
        float128_res = f1 >= f2;
        EXPECT_TRUE(float128_res == res) << "operator>=: "
                                         << "value1=" << value1 << ", value2=" << value2;

        res = value1 < value2;
        float128_res = f1 < f2;
        EXPECT_TRUE(float128_res == res) << "operator<: "
                                         << "value1=" << value1 << ", value2=" << value2;

        res = value1 <= value2;
        float128_res = f1 <= f2;
        EXPECT_TRUE(float128_res == res) << "operator<=: "
                                         << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(float128, CompareDouble)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = get_double_random();
        double value2 = get_double_random();
        bool res = value1 > value2;
        float128 f1 = value1;

        bool float128_res = f1 > value2;
        EXPECT_TRUE(float128_res == res) << "operator>: "
                                         << "value1=" << value1 << ", value2=" << value2;

        res = value1 >= value2;
        float128_res = f1 >= value2;
        EXPECT_TRUE(float128_res == res) << "operator>=: "
                                         << "value1=" << value1 << ", value2=" << value2;

        res = value1 < value2;
        float128_res = f1 < value2;
        EXPECT_TRUE(float128_res == res) << "operator<: "
                                         << "value1=" << value1 << ", value2=" << value2;

        res = value1 <= value2;
        float128_res = f1 <= value2;
        EXPECT_TRUE(float128_res == res) << "operator<=: "
                                         << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(float128, CompareFloat)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        float value1 = (float)get_double_random();
        float value2 = (float)get_double_random();
        bool res = value1 > value2;
        float128 f1 = value1;

        bool float128_res = f1 > value2;
        EXPECT_TRUE(float128_res == res) << "operator>: "
                                         << "value1=" << value1 << ", value2=" << value2;

        res = value1 >= value2;
        float128_res = f1 >= value2;
        EXPECT_TRUE(float128_res == res) << "operator>=: "
                                         << "value1=" << value1 << ", value2=" << value2;

        res = value1 < value2;
        float128_res = f1 < value2;
        EXPECT_TRUE(float128_res == res) << "operator<: "
                                         << "value1=" << value1 << ", value2=" << value2;

        res = value1 <= value2;
        float128_res = f1 <= value2;
        EXPECT_TRUE(float128_res == res) << "operator<=: "
                                         << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(float128, CompareInt32)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        int32_t value1 = get_int32_random();
        int32_t value2 = get_int32_random();
        bool res = value1 > value2;
        float128 f1 = value1;

        bool float128_res = f1 > value2;
        EXPECT_TRUE(float128_res == res) << "operator>: "
                                         << "value1=" << value1 << ", value2=" << value2;

        res = value1 >= value2;
        float128_res = f1 >= value2;
        EXPECT_TRUE(float128_res == res) << "operator>=: "
                                         << "value1=" << value1 << ", value2=" << value2;

        res = value1 < value2;
        float128_res = f1 < value2;
        EXPECT_TRUE(float128_res == res) << "operator<: "
                                         << "value1=" << value1 << ", value2=" << value2;

        res = value1 <= value2;
        float128_res = f1 <= value2;
        EXPECT_TRUE(float128_res == res) << "operator<=: "
                                         << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(float128, CompareUnsignedInt32)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        uint32_t value1 = get_uint32_random();
        uint32_t value2 = get_uint32_random();
        bool res = value1 > value2;
        float128 f1 = value1;

        bool float128_res = f1 > value2;
        EXPECT_TRUE(float128_res == res) << "operator>: "
                                         << "value1=" << value1 << ", value2=" << value2;

        res = value1 >= value2;
        float128_res = f1 >= value2;
        EXPECT_TRUE(float128_res == res) << "operator>=: "
                                         << "value1=" << value1 << ", value2=" << value2;

        res = value1 < value2;
        float128_res = f1 < value2;
        EXPECT_TRUE(float128_res == res) << "operator<: "
                                         << "value1=" << value1 << ", value2=" << value2;

        res = value1 <= value2;
        float128_res = f1 <= value2;
        EXPECT_TRUE(float128_res == res) << "operator<=: "
                                         << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(float128, CompareInt64)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        int64_t value1 = get_int64_random();
        int64_t value2 = get_int64_random();
        bool res = value1 > value2;
        float128 f1 = value1;

        bool float128_res = f1 > value2;
        EXPECT_TRUE(float128_res == res) << "operator>: "
                                         << "value1=" << value1 << ", value2=" << value2;

        res = value1 >= value2;
        float128_res = f1 >= value2;
        EXPECT_TRUE(float128_res == res) << "operator>=: "
                                         << "value1=" << value1 << ", value2=" << value2;

        res = value1 < value2;
        float128_res = f1 < value2;
        EXPECT_TRUE(float128_res == res) << "operator<: "
                                         << "value1=" << value1 << ", value2=" << value2;

        res = value1 <= value2;
        float128_res = f1 <= value2;
        EXPECT_TRUE(float128_res == res) << "operator<=: "
                                         << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(float128, CompareUnsignedInt64)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        uint64_t value1 = get_uint64_random();
        uint64_t value2 = get_uint64_random();
        bool res = value1 > value2;
        float128 f1 = value1;

        bool float128_res = f1 > value2;
        EXPECT_TRUE(float128_res == res) << "operator>: "
                                         << "value1=" << value1 << ", value2=" << value2;

        res = value1 >= value2;
        float128_res = f1 >= value2;
        EXPECT_TRUE(float128_res == res) << "operator>=: "
                                         << "value1=" << value1 << ", value2=" << value2;

        res = value1 < value2;
        float128_res = f1 < value2;
        EXPECT_TRUE(float128_res == res) << "operator<: "
                                         << "value1=" << value1 << ", value2=" << value2;

        res = value1 <= value2;
        float128_res = f1 <= value2;
        EXPECT_TRUE(float128_res == res) << "operator<=: "
                                         << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(float128, OperatorEqual)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = get_double_random();
        float128 f1 = value1;
        bool float128_res = f1 == f1;
        EXPECT_TRUE(float128_res == true) << "operator==: "
                                          << "value1=" << value1;
    }
}
TEST(float128, OperatorNotEqual)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = get_double_random();
        float128 f1 = value1;
        bool float128_res = f1 != f1;
        EXPECT_TRUE(float128_res == false) << "operator!=: "
                                           << "value1=" << value1;
    }
}
TEST(float128, TemplateOperatorNotEqual)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = get_double_random();
        float128 f1 = value1;
        bool float128_res = f1 != value1;
        EXPECT_TRUE(float128_res == false) << "operator!=<double>: "
                                           << "value1=" << value1;

        f1 = static_cast<float>(value1);
        float128_res = f1 != static_cast<float>(value1);
        EXPECT_TRUE(float128_res == false) << "operator!=<float>: "
                                           << "value1=" << value1;

        f1 = static_cast<uint64_t>(value1);
        float128_res = f1 != static_cast<uint64_t>(value1);
        EXPECT_TRUE(float128_res == false) << "operator!=<uint64_t>: "
                                           << "value1=" << value1;

        f1 = static_cast<int64_t>(value1);
        float128_res = f1 != static_cast<int64_t>(value1);
        EXPECT_TRUE(float128_res == false) << "operator!=<int64_t>: "
                                           << "value1=" << value1;

        f1 = static_cast<uint32_t>(value1);
        float128_res = f1 != static_cast<uint32_t>(value1);
        EXPECT_TRUE(float128_res == false) << "operator!=<uint32_t>: "
                                           << "value1=" << value1;

        f1 = static_cast<int32_t>(value1);
        float128_res = f1 != static_cast<int32_t>(value1);
        EXPECT_TRUE(float128_res == false) << "operator!=<int32_t>: "
                                           << "value1=" << value1;
    }
}
TEST(float128, floor)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random();
        double res = ::floor(value);
        float128 f1 = value;
        float128 float128_res = floor(f1);
        EXPECT_DOUBLE_EQ(float128_res, res) << "floor: "
                                            << "value=" << value;
    }
}
TEST(float128, ceil)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random();
        double res = ::ceil(value);
        float128 f1 = value;
        float128 float128_res = ceil(f1);
        EXPECT_DOUBLE_EQ(float128_res, res) << "ceil: "
                                            << "value=" << value;
    }
}
TEST(float128, trunc)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random();
        double res = trunc(value);
        float128 f1 = value;
        float128 float128_res = trunc(f1);
        EXPECT_DOUBLE_EQ(float128_res, res) << "trunc: "
                                            << "value=" << value;
    }
}
TEST(float128, round)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random();
        double res = round(value);
        float128 f1 = value;
        float128 float128_res = round(f1);
        EXPECT_DOUBLE_EQ(float128_res, res) << "round: "
                                            << "value=" << value;
    }
}
TEST(float128, copysign)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = get_double_random();
        double value2 = get_double_random();
        double res = copysign(value1, value2);
        float128 f1 = value1;
        float128 f2 = value2;
        float128 float128_res = copysign(f1, f2);
        EXPECT_DOUBLE_EQ(float128_res, res) << "copysign: "
                                            << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(float128, fmod_fraction)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = get_double_random();
        double value2 = get_double_random();
        // check if double precision can produce accurate results
        if (abs(ilogb(value1) - ilogb(value2)) > 52)
            continue;
        double res = ::fmod(value1, value2);
        float128 f1 = value1;
        float128 f2 = value2;
        float128 float128_res = fmod(f1, f2);
        EXPECT_DOUBLE_EQ(float128_res, res) << "fmod: "
                                            << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(float128, fmod_integer)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = get_int32_random();
        double value2 = get_int32_random();
        double res = ::fmod(value1, value2);
        float128 f1 = value1;
        float128 f2 = value2;
        float128 float128_res = fmod(f1, f2);
        EXPECT_DOUBLE_EQ(float128_res, res) << "fmod: "
                                            << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(float128, modf)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random();
        double res_int = 0;
        double res = ::modf(value, &res_int);
        float128 f1 = value;
        float128 float128_res_int;
        float128 float128_res_frac = modf(f1, &float128_res_int);
        EXPECT_DOUBLE_EQ(float128_res_frac, res) << "modf fraction: "
                                                 << "value=" << value;
        EXPECT_DOUBLE_EQ(float128_res_int, res_int) << "modf integer: "
                                                    << "value=" << value;
    }
}
TEST(float128, sqrt)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = fabs(get_double_random());
        double res = ::sqrt(value);
        float128 f1 = value;
        float128 float128_res = sqrt(f1);
        EXPECT_DOUBLE_EQ(float128_res, res) << "sqrt: "
                                            << "value=" << value;
    }
}
TEST(float128, cbrt)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = fabs(get_double_random());
        double res = ::cbrt(value);
        float128 f1 = value;
        float128 float128_res = cbrt(f1);
        EXPECT_DOUBLE_EQ(float128_res, res) << "cbrt: "
                                            << "value=" << value;
    }
}
TEST(float128, frexp)
{
    int res_exp, float128_exp;
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random();
        double res = ::frexp(value, &res_exp);
        float128 f1 = value;
        float128 float128_res = frexp(f1, &float128_exp);
        EXPECT_DOUBLE_EQ(float128_res, res) << "frexp mantissa: "
                                            << "value=" << value;
        EXPECT_EQ(float128_exp, res_exp) << "frexp exp: "
                                         << "value=" << value;
    }
}
TEST(float128, ldexp)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random();
        int value_exp = get_int32_random() % 300;
        double res = ::ldexp(value, value_exp);
        float128 f1 = value;
        float128 float128_res = ldexp(f1, value_exp);
        EXPECT_DOUBLE_EQ(float128_res, res) << "ldexp: "
                                            << "value=" << value;
    }
}
TEST(float128, erf)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random(-60, 5);
        double res = ::erf(value);
        float128 f1 = value;
        float128 float128_res = erf(f1);
        EXPECT_DOUBLE_EQ(float128_res, res) << "erf: "
                                            << "value=" << value;
    }
}
TEST(float128, erfc)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random(-60, 4);
        double res = ::erfc(value);
        float128 f1 = value;
        float128 float128_res = erfc(f1);
        EXPECT_DOUBLE_EQ(float128_res, res) << "erfc: "
                                            << "value=" << value;
    }
}
TEST(float128, reciprocal)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random();
        double res = 1.0 / value;
        float128 f1 = value;
        float128 float128_res = reciprocal(f1);
        EXPECT_DOUBLE_EQ(float128_res, res) << "reciprocal: "
                                            << "value=" << value;
    }
}
TEST(float128, hypot)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = get_double_random();
        double value2 = get_double_random();
        double res = hypot(value1, value2);
        float128 f1 = value1;
        float128 f2 = value2;
        // float128 t = f1 * f1;
        // float128 t2 = sqrt(t);

        float128 float128_res = hypot(f1, f2);
        EXPECT_DOUBLE_EQ(float128_res, res) << "hypot: "
                                            << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(float128, ilogb)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random();
        int32_t res = ::ilogb(value);
        float128 f1 = value;
        int32_t float128_res = ilogb(f1);
        EXPECT_EQ(float128_res, res) << "ilogb: "
                                     << "value=" << value;
    }
}
TEST(float128, log)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = fabs(get_double_random());
        double res = log(value);
        float128 f1 = value;
        float128 float128_res = log(f1);
        EXPECT_DOUBLE_EQ(float128_res, res) << "log: "
                                            << "value=" << value;
    }
}
TEST(float128, log2)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = fabs(get_double_random());
        double res = log2(value);
        float128 f1 = value;
        float128 float128_res = log2(f1);
        EXPECT_DOUBLE_EQ(float128_res, res) << "log2: "
                                            << "value=" << value;
    }
}
TEST(float128, log10)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = fabs(get_double_random());
        double res = log10(value);
        float128 f1 = value;
        float128 float128_res = log10(f1);
        EXPECT_DOUBLE_EQ(float128_res, res) << "log10: "
                                            << "value=" << value;
    }
}
TEST(float128, log1p)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = fabs(get_double_random());
        double res = log1p(value);
        float128 f1 = value;
        float128 float128_res = log1p(f1);
        EXPECT_DOUBLE_EQ(float128_res, res) << "log1p: "
                                            << "value=" << value;
    }
}
TEST(float128, logb)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = fabs(get_double_random(-40, 15));  // lower exponent results in lost bits
        double res = logb(value);
        float128 f1 = value;
        float128 float128_res = logb(f1);
        EXPECT_DOUBLE_EQ(float128_res, res) << "logb: "
                                            << "value=" << value;
    }
}
TEST(float128, exp)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random(-16, 16);
        double res = exp(value);
        float128 f1 = value;
        float128 float128_res = exp(f1);
        EXPECT_DOUBLE_EQ(float128_res, res) << "exp: "
                                            << "value=" << value;
    }
}
TEST(float128, exp2)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random(-16, 16);  // lower exponent results in lost bits
        double res = exp2(value);
        float128 f1 = value;
        float128 float128_res = exp2(f1);
        EXPECT_DOUBLE_EQ(float128_res, res) << "exp2: "
                                            << "value=" << value;
    }
}
TEST(float128, expm1)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random(-16, 16);  // lower exponent results in lost bits
        double res = expm1(value);
        float128 f1 = value;
        float128 float128_res_ = expm1(f1);
        double float128_res = (double)float128_res_;
        EXPECT_DOUBLE_EQ(float128_res, res) << "expm1: "
                                            << "value=" << value;
    }
}
TEST(float128, pow)
{
    // test special values
    double values[][2] = {{0, 1}, {1, INFINITY}, {1.1, INFINITY}, {INFINITY, 1}, {NAN, 1}, {1, NAN}, {-2, 5}, {-2, 5.5}};
    constexpr auto value_count = array_length(values);
    for (auto i = 0u; i < value_count; ++i) {
        double value1 = values[i][0];
        double value2 = values[i][1];
        double res = pow(value1, value2);
        float128 f1 = value1;
        float128 f2 = value2;
        float128 float128_res = pow(f1, f2);
        // since NaNs fail to compare, skip the exception
        if (isnan(res) && isnan(float128_res))
            continue;

        EXPECT_DOUBLE_EQ(float128_res, res) << "pow: "
                                            << " value1=" << value1 << ", value2=" << value2;
    }

    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        // double value1 = fabs(get_double_random());
        double value1 = fabs(get_double_random());
        double value2 = get_double_random(-16, 16);  // lower exponent results in lost bits
        double res = pow(value1, value2);
        float128 f1 = value1;
        float128 f2 = value2;
        float128 float128_res = pow(f1, f2);
        EXPECT_DOUBLE_EQ(float128_res, res) << "pow: "
                                            << " value1=" << value1 << ", value2=" << value2;
    }
}

TEST(float128, sin)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random(-60, 3);
        double res = ::sin(value);
        float128 f1 = value;
        float128 float128_res = sin(f1);
        EXPECT_DOUBLE_EQ(float128_res, res) << "sin: "
                                            << "value=" << value;
    }
}
TEST(float128, cos)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random(-60, 3);
        double res = ::cos(value);
        float128 f1 = value;
        float128 float128_res = cos(f1);
        EXPECT_DOUBLE_EQ(float128_res, res) << "cos: "
                                            << "value=" << value;
    }
}
TEST(float128, tan)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random(-60, 3);
        double res = ::tan(value);
        float128 f1 = value;
        float128 float128_res = tan(f1);
        EXPECT_DOUBLE_EQ(float128_res, res) << "tan: "
                                            << "value=" << value;
    }
}
TEST(float128, asin)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random(-60, -1);
        double res = ::asin(value);
        float128 f1 = value;
        float128 float128_res = asin(f1);
        EXPECT_DOUBLE_EQ(float128_res, res) << "asin: "
                                            << "value=" << value;
    }
}
TEST(float128, acos)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random(-60, -1);
        double res = acos(value);
        float128 f1 = value;
        float128 float128_res = acos(f1);
        EXPECT_DOUBLE_EQ(float128_res, res) << "acos: "
                                            << "value=" << value;
    }
}
TEST(float128, atan)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random(-60, 14);
        double res = atan(value);
        float128 f1 = value;
        float128 float128_res = atan(f1);
        EXPECT_DOUBLE_EQ(float128_res, res) << "atan: "
                                            << "value=" << value;
    }
}
TEST(float128, atan2)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = get_double_random(-10, 14);  // lower exponent results in lost bits
        double value2 = get_double_random(-10, 14);  // lower exponent results in lost bits
        double res = atan2(value1, value2);
        float128 f1 = value1;
        float128 f2 = value2;
        float128 float128_res = atan2(f1, f2);
        EXPECT_DOUBLE_EQ(float128_res, res) << "atan2: "
                                            << " value1=" << value1 << ", value2=" << value2;
    }
}
TEST(float128, sinh)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random(-60, 2);
        double res = ::sinh(value);
        float128 f1 = value;
        float128 float128_res = sinh(f1);
        EXPECT_DOUBLE_EQ(float128_res, res) << "sinh: "
                                            << "value=" << value;
    }
}
TEST(float128, asinh)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random(-60, 2);
        double res = ::asinh(value);
        float128 f1 = value;
        float128 float128_res = asinh(f1);
        EXPECT_DOUBLE_EQ(float128_res, res) << "asinh: "
                                            << "value=" << value;
    }
}
TEST(float128, cosh)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random(-60, 2);
        double res = ::cosh(value);
        float128 f1 = value;
        float128 float128_res = cosh(f1);
        EXPECT_DOUBLE_EQ(float128_res, res) << "cosh: "
                                            << "value=" << value;
    }
}
TEST(float128, acosh)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = 1.0 + fabs(get_double_random(-60, 2));  // must be >= 1
        double res = ::acosh(value);
        float128 f1 = value;
        float128 float128_res = acosh(f1);
        EXPECT_DOUBLE_EQ(float128_res, res) << "acosh: "
                                            << "value=" << value;
    }
}
TEST(float128, tanh)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random(-60, 2);
        double res = ::tanh(value);
        float128 f1 = value;
        float128 float128_res = tanh(f1);
        EXPECT_DOUBLE_EQ(float128_res, res) << "tanh: "
                                            << "value=" << value;
    }
}
TEST(float128, atanh)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random(-60, -1);  // must be: -1 < value < 1
        double res = ::atanh(value);
        float128 f1 = value;
        float128 float128_res = atanh(f1);
        EXPECT_DOUBLE_EQ(float128_res, res) << "atanh: "
                                            << "value=" << value;
    }
}

/**********************************************************************
 * float128 regression tests
 *
 * Each test below pins down a defect that was fixed.
 ***********************************************************************/

// The comparison operators compared the raw 128 bits, which is not how IEEE 754 orders floats:
// a NaN has to compare unordered with everything including itself, and the two zeros have to
// compare equal despite their differing sign bit.
TEST(float128, IEEEComparisonSemantics)
{
    const float128 zero = 0;
    const float128 neg_zero = -float128(0);
    const float128 qnan = float128::nan();
    const float128 one = 1;

    // a NaN equals nothing, not even itself
    EXPECT_FALSE(qnan == qnan);
    EXPECT_TRUE(qnan != qnan);
    EXPECT_FALSE(qnan == one);
    EXPECT_TRUE(qnan != one);

    // every relational test against a NaN is false, in both directions
    EXPECT_FALSE(qnan < one);
    EXPECT_FALSE(qnan > one);
    EXPECT_FALSE(qnan <= one);
    EXPECT_FALSE(qnan >= one);
    EXPECT_FALSE(one < qnan);
    EXPECT_FALSE(one > qnan);
    EXPECT_FALSE(one <= qnan);
    EXPECT_FALSE(one >= qnan);
    EXPECT_FALSE(qnan <= qnan);
    EXPECT_FALSE(qnan >= qnan);

    // the two zeros are numerically equal
    EXPECT_TRUE(zero == neg_zero);
    EXPECT_FALSE(zero != neg_zero);
    EXPECT_FALSE(neg_zero < zero);
    EXPECT_FALSE(zero > neg_zero);
    EXPECT_TRUE(neg_zero <= zero);
    EXPECT_TRUE(zero >= neg_zero);

    // ordinary ordering still works
    EXPECT_TRUE(float128(-1) < zero);
    EXPECT_TRUE(one > zero);
    EXPECT_TRUE(-float128::inf() < float128::inf());
    EXPECT_TRUE(float128::inf() > one);
}
// operator bool tested the raw words, so negative zero reported as a non zero value.
TEST(float128, NegativeZeroIsFalsy)
{
    const float128 neg_zero = -float128(0);
    EXPECT_TRUE(neg_zero.is_zero());
    EXPECT_FALSE(static_cast<bool>(neg_zero));
    EXPECT_TRUE(!neg_zero);
    EXPECT_TRUE(neg_zero.is_negative());  // the sign itself is preserved

    // the sign of a zero product is the combination of the operand signs
    const float128 p = float128(0) * float128(-5);
    EXPECT_TRUE(p.is_zero());
    EXPECT_TRUE(p.is_negative());
}
// fmod returned an infinity for a zero divisor, the CRT returns a NaN.
TEST(float128, FmodSpecialCases)
{
    EXPECT_TRUE(fmod(float128(5), float128(0)).is_nan());
    EXPECT_TRUE(fmod(float128(-5), float128(0)).is_nan());
    EXPECT_TRUE(fmod(float128::inf(), float128(5)).is_nan());
    EXPECT_TRUE(fmod(float128::nan(), float128(5)).is_nan());
    EXPECT_TRUE(fmod(float128(5), float128::nan()).is_nan());
    // an infinite divisor leaves the dividend alone
    EXPECT_TRUE(fmod(float128(5), float128::inf()) == float128(5));
    // ordinary behaviour is unchanged
    EXPECT_TRUE(fmod(float128(7), float128(3)) == float128(1));
}
// The integer digits were extracted with a plain divide, which leaves the discarded digit behind
// as a fraction and corrupts every digit produced after it: 123456789 printed as 123467899.
TEST(float128, IntegerToStringIsExact)
{
    const uint64_t values[] = {0ull, 1ull, 9ull, 10ull, 99ull, 100ull, 12345ull, 1000000ull,
                               123456789ull, 999999999999999999ull, 12345678901234567890ull};
    // The fixed layout is asked for explicitly. The default conversion produces the shortest text
    // that reads back as the same value, and for a round number that is the scientific form:
    // 1000000 comes out as 1e+06, which is what std::format gives a double as well.
    char buf[64];
    for (auto v : values) {
        snprintf(buf, sizeof(buf), "%llu", v);
        EXPECT_EQ(std::format("{:.0f}", float128(v)), buf);
    }

    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const uint64_t v = get_uint64_random() >> (get_uint32_random() % 64);
        snprintf(buf, sizeof(buf), "%llu", v);
        EXPECT_EQ(std::format("{:.0f}", float128(v)), buf) << "value=" << v;
    }
}
// to_e_format fed a negative value straight to log10, producing a NaN exponent of INT32_MIN and a
// mantissa of "inf". It also padded the mantissa out to the whole buffer with '0' characters.
TEST(float128, ScientificNotationFormat)
{
    // these exponents force the e-format path
    const std::string big = (std::string)float128("1e100");
    const std::string small = (std::string)float128("1e-100");
    const std::string neg_small = (std::string)float128("-1e-100");

    EXPECT_EQ(big.find("inf"), std::string::npos) << big;
    EXPECT_EQ(small.find("inf"), std::string::npos) << small;
    EXPECT_EQ(neg_small.find("inf"), std::string::npos) << neg_small;
    EXPECT_EQ(neg_small[0], '-') << neg_small;
    // no run of padding zeros, the mantissa is short
    EXPECT_LT(big.size(), 45u) << big;
    EXPECT_LT(neg_small.size(), 50u) << neg_small;
    // and the value survives the trip
    EXPECT_NEAR(static_cast<double>(float128(big.c_str())), 1e100, 1e86);
    EXPECT_NEAR(static_cast<double>(float128(neg_small.c_str())), -1e-100, 1e-114);
}
// The fraction was accumulated by multiplying a stored, and therefore inexact, 10^-9 constant
// once per group of 9 digits, so the error compounded. Values that are exactly representable in
// binary did not survive: "0.5" came out one unit in the last place low.
TEST(float128, ExactDecimalsParseExactly)
{
    struct {
        const char* text;
        double value;
    } cases[] = {{"0.5", 0.5},   {"0.25", 0.25},     {"0.125", 0.125},   {"2.5", 2.5},
                 {"-0.75", -0.75}, {"0.0625", 0.0625}, {"-8.5", -8.5},   {"1024.0", 1024.0}};
    for (const auto& c : cases) {
        EXPECT_TRUE(float128(c.text) == float128(c.value)) << c.text << " parsed as " << (std::string)float128(c.text);
    }
    EXPECT_TRUE(float128("0.5") == float128::half());
}

// The fraction digits used to be produced by repeatedly multiplying a float128 by 100000, which
// rounds once per group and compounds. They now come from an exact fixed point expansion of the
// encoding, so a value whose decimal expansion terminates prints all of its digits and no more.
TEST(float128, FractionDigitsAreExact)
{
    // every negative power of two terminates in decimal and must print exactly
    struct {
        const char* text;
        const char* expected;
    } cases[] = {{"0.5", "0.5"},         {"0.25", "0.25"},       {"0.125", "0.125"},
                 {"0.0625", "0.0625"},   {"0.03125", "0.03125"}, {"2.5", "2.5"},
                 {"-0.75", "-0.75"},     {"1.5", "1.5"},         {"-8.25", "-8.25"}};
    for (const auto& c : cases) {
        EXPECT_STREQ(static_cast<char*>(float128(c.text)), c.expected);
    }

    // a power of two reached by shifting has an exact, terminating expansion too
    float128 v = 1;
    std::string expected = "0.5";
    for (int i = 1; i <= 20; ++i) {
        v >>= 1;
        EXPECT_STREQ(static_cast<char*>(v), expected.c_str()) << "2^-" << i;
        // 2^-(i+1) is the previous string with its last digit halved, build it the easy way
        expected = (std::string)(v >> 1);
    }
}
// Rounding the multiply and the divide from a single guard bit, with no knowledge of the
// discarded low words, made ties resolve arbitrarily. With a real sticky bit the operations
// round half to even, which tightens the decimal round trip.
TEST(float128, RoundTripAccuracy)
{
    srand(RANDOM_SEED);
    int exact = 0, total = 0;
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const float128 v(get_uint64_random(), get_uint64_random() & 0xFFFFFFFFFFFFull,
                         0x3FFF + static_cast<int>(get_uint32_random() % 60) - 30, 0);
        if (v.is_zero() || v.is_special())
            continue;
        ++total;
        const float128 back = ((std::string)v).c_str();
        if (back == v)
            ++exact;
        // whatever else happens, the value must stay in the same neighbourhood
        EXPECT_NEAR(static_cast<double>(back), static_cast<double>(v), fabs(static_cast<double>(v)) * 1e-30);
    }
    // measured at 69% when this was written, well above the 58% the old code managed. The bound
    // is deliberately loose: it guards against a regression, it is not a precision promise.
    EXPECT_GT(exact * 100 / total, 60) << exact << " of " << total << " round tripped exactly";
}

// A float128 shift scales by a power of two. With a scalar on the left the left operand is widened,
// so the result is a float128 rather than the builtin type of the scalar. The count is the right
// operand converted to int32_t, which for float128 rounds to nearest.
TEST(float128, ScalarOnTheLeftHandSideShifts)
{
    EXPECT_DOUBLE_EQ(static_cast<double>(1 << float128(10)), 1024.0);
    EXPECT_DOUBLE_EQ(static_cast<double>(1024 >> float128(10)), 1.0);
    EXPECT_DOUBLE_EQ(static_cast<double>(3 << float128(2)), 12.0);
    EXPECT_DOUBLE_EQ(static_cast<double>(1 >> float128(2)), 0.25);
    EXPECT_DOUBLE_EQ(static_cast<double>(-1 << float128(3)), -8.0);
    EXPECT_DOUBLE_EQ(static_cast<double>(2.5 << float128(2)), 10.0);

    // A fractional shift count goes through operator int32_t, which rounds to nearest here rather
    // than truncating, so 3.75 shifts by 4. The same overload picked with the float128 on the left
    // has to agree, both reaching the count the same way.
    EXPECT_TRUE((1 << float128(3.75)) == (1 << float128(4)));
    EXPECT_TRUE((1 << float128(3.25)) == (1 << float128(3)));
    EXPECT_TRUE((1 << float128(3.75)) == (float128(1) << float128(3.75)));

    static_assert(std::is_same_v<decltype(1 << float128(3)), float128>);

    // the float128 on the left forms must still resolve
    EXPECT_TRUE((float128(4) << 1) == float128(8));
    EXPECT_TRUE((float128(4) >> 1) == float128(2));

    // the arithmetic scalar on the left forms that already existed must still resolve
    EXPECT_TRUE((1 + float128(2)) == float128(3));
    EXPECT_TRUE((1 - float128(2)) == float128(-1));
    static_assert(std::is_same_v<decltype(1 + float128(2)), float128>);
}

/**********************************************************************
 * Compile time (constexpr) evaluation
 *
 * Each test below pairs static_asserts, which fail the build the moment one of these operations
 * stops being usable in a constant expression, with a runtime check of the same expression built
 * from opaque() values the optimizer cannot fold. See the matching block in
 * fixed_point128_gtest.cpp for why the runtime half is not a restatement of the compile time one.
 *
 * float128 needs one thing beyond the shared intrinsic wrappers: the fields of the high QWORD are
 * read through shift and mask accessors rather than the _float128_bits view, because reading the
 * inactive member of a union is not allowed while the compiler is evaluating.
 ***********************************************************************/
TEST(float128, ConstexprConstructionAndConversion)
{
    constexpr float128 zero;
    constexpr float128 one = float128::one();
    constexpr float128 half = float128::half();
    constexpr float128 from_int(42);
    constexpr float128 from_neg(-42);
    constexpr float128 copied(from_int);
    constexpr float128 assigned = [] { float128 a; a = float128::one(); return a; }();

    static_assert(zero.is_zero());
    static_assert(static_cast<int32_t>(from_int) == 42);
    static_assert(static_cast<int64_t>(from_neg) == -42);
    static_assert(static_cast<uint64_t>(from_int) == 42ull);
    static_assert(static_cast<uint32_t>(from_int) == 42u);
    static_assert(copied == from_int && assigned == one);
    static_assert(from_int.is_positive() && from_neg.is_negative());
    static_assert(from_int.is_normal() && from_int.is_finite() && !from_int.is_special());
    static_assert(from_int.is_int() && !half.is_int());
    static_assert(one.get_exponent() == 0 && half.get_exponent() == -1);
    static_assert(one.is_exponent_of_2() && !float128(3).is_exponent_of_2());
    static_assert(one.get_bit(0) == 0);          // 1.0 stores a zero fraction, the unity bit is implicit
    static_assert(one.get_bit(127) == 0);        // bit 127 is the sign
    static_assert((-one).get_bit(127) == 1);

    // The most negative integers, whose magnitude has no positive counterpart. The integral
    // constructor used to reach these by negating the signed value, which is undefined.
    static_assert(float128(INT64_MIN) == -(float128::one() << 63));
    static_assert(float128(INT32_MIN) == -(float128::one() << 31));
    static_assert(static_cast<int64_t>(float128(INT64_MIN)) == INT64_MIN);
    static_assert(float128(INT64_MIN) < float128(INT64_MIN + 1));

    EXPECT_TRUE(float128(opaque(42)) == from_int);
    EXPECT_TRUE(float128(static_cast<int64_t>(opaque(-1)) * (INT64_MAX / 2 + 1) * 2) == float128(INT64_MIN));
}
TEST(float128, ConstexprSpecialValues)
{
    constexpr float128 one = float128::one();

    static_assert(float128::inf().is_inf() && float128::inf().is_special());
    static_assert(!float128::inf().is_finite() && !float128::inf().is_normal());
    static_assert(float128::nan().is_nan() && float128::nan().is_special());
    static_assert(isinf(float128::inf()) && isnan(float128::nan()));
    static_assert(!isfinite(float128::inf()) && isfinite(one));
    static_assert((-float128::inf()).is_negative() && (-float128::inf()).is_inf());
    static_assert(float128::nan() != float128::nan());  // a NaN compares unequal to everything
    static_assert(float128().is_zero() && (-float128()).is_zero());
    static_assert(float128() == -float128());  // the two zeros compare equal despite differing bits

    // infinity arithmetic follows IEEE 754
    static_assert((float128::inf() + one).is_inf());
    static_assert((float128::inf() - float128::inf()).is_nan());  // the invalid operation
    static_assert((float128::nan() + one).is_nan());

    // nextUp and nextDown walk one representable value at a time
    static_assert(float128::nextUp(one) > one);
    static_assert(float128::nextDown(one) < one);
    static_assert(float128::nextDown(float128::nextUp(one)) == one);

    EXPECT_TRUE(isnan(float128::nan() + float128(opaque(1))));
    EXPECT_TRUE(float128::nextUp(float128(opaque(1))) == float128::nextUp(one));
}
TEST(float128, ConstexprArithmetic)
{
    constexpr float128 one = float128::one();
    constexpr float128 half = float128::half();

    static_assert(one + one == float128(2));
    static_assert(float128(10) - float128(3) == float128(7));
    static_assert(half + half == one);
    static_assert(one - one == float128());
    static_assert(float128(6) * float128(7) == float128(42));
    static_assert(float128(-6) * float128(7) == float128(-42));
    static_assert(half * float128(8) == float128(4));
    static_assert(sqr(float128(-9)) == float128(81));
    static_assert(sqr(half) == half * half);  // sqr is bit identical to x * x
    constexpr float128 squared = [] { float128 a(5); a.square(); return a; }();
    static_assert(squared == float128(25));
    constexpr float128 summed = [] { float128 a; for (int i = 1; i <= 10; ++i) a += float128(i); return a; }();
    static_assert(summed == float128(55));

    // the scalar may sit on either side
    static_assert(float128(3) + 4 == float128(7));
    static_assert(4 + float128(3) == float128(7));
    static_assert(10 - float128(3) == float128(7));
    static_assert(7 * float128(6) == float128(42));

    // the shifts scale by a power of two
    static_assert((one << 3) == float128(8));
    static_assert((float128(8) >> 3) == one);
    static_assert((one >> 1) == half);
    constexpr float128 shifted = [] { float128 a(16); a >>= 2; a <<= 1; return a; }();
    static_assert(shifted == float128(8));

    static_assert(-one == float128(-1) && +one == one);
    static_assert(!float128() && static_cast<bool>(one));
    static_assert(float128(1) < float128(2) && float128(2) > float128(1));
    static_assert(float128(-2) < float128(-1) && float128(-1) < float128(1));
    static_assert(float128(1) <= 1 && float128(1) >= 1 && float128(1) == 1 && float128(1) != 2);

    EXPECT_TRUE(float128(opaque(6)) * float128(opaque(7)) == float128(42));
    EXPECT_TRUE(float128(opaque(10)) - float128(opaque(3)) == float128(7));
    EXPECT_TRUE(opaque(4) + float128(3) == float128(7));
    EXPECT_TRUE(sqr(float128(opaque(-9))) == float128(81));
    // the loop below runs the full add path, exponent alignment included
    EXPECT_TRUE(summed == [] { float128 a; for (int i = 1; i <= 10; ++i) a += float128(opaque(i)); return a; }());
}
TEST(float128, ConstexprMathFunctions)
{
    constexpr float128 half = float128::half();
    constexpr float128 two_and_a_half = float128(5) >> 1;

    static_assert(fabs(float128(-5)) == float128(5));
    static_assert(floor(two_and_a_half) == float128(2));
    static_assert(floor(-two_and_a_half) == float128(-3));
    static_assert(ceil(two_and_a_half) == float128(3));
    static_assert(trunc(-two_and_a_half) == float128(-2));
    static_assert(round(two_and_a_half) == float128(3));  // ties round away from zero
    static_assert(round(-two_and_a_half) == float128(-3));
    static_assert(llround(two_and_a_half) == 3);
    static_assert(lround(-two_and_a_half) == -3);
    static_assert(copysign(float128(5), float128(-1)).is_negative());
    static_assert(fmin(float128(1), float128(2)) == float128(1));
    static_assert(fmax(float128(1), float128(2)) == float128(2));
    static_assert(fdim(float128(5), float128(3)) == float128(2));
    static_assert(fdim(float128(3), float128(5)).is_zero());
    static_assert(ilogb(float128(8)) == 3);
    static_assert(fma(float128(6), float128(7), float128(1)) == float128(43));
    static_assert(ldexp(float128(3), 4) == float128(48));
    static_assert(ldexp(float128(48), -4) == float128(3));

    // modf and frexp report through a pointer, so they run inside a lambda
    constexpr float128 modf_frac = [] { float128 ip; return modf(float128(5) >> 1, &ip); }();
    constexpr float128 modf_int = [] { float128 ip; (void)modf(float128(5) >> 1, &ip); return ip; }();
    static_assert(modf_frac == half && modf_int == float128(2));
    constexpr int frexp_expo = [] { int e = 0; (void)frexp(float128(8), &e); return e; }();
    constexpr float128 frexp_mant = [] { int e = 0; return frexp(float128(8), &e); }();
    static_assert(frexp_expo == 4 && frexp_mant == half);

    static_assert(float128::exp10(0) == float128::one());
    static_assert(float128::exp10(3) == float128(1000));
    static_assert(float128::exp10(-1) == float128::tenth());

    EXPECT_TRUE(round(-(float128(opaque(5)) >> 1)) == float128(-3));
    EXPECT_TRUE(floor(-(float128(opaque(5)) >> 1)) == float128(-3));
    EXPECT_TRUE(ilogb(float128(opaque(8))) == 3);
    EXPECT_TRUE(float128::exp10(opaque(3)) == float128::exp10(3));
}

// The double and float conversions used to be the one thing that could not happen at compile time:
// they read the inactive member of a union, which constant evaluation rejects. Now that Double and
// Float convert with std::bit_cast, a floating point literal crosses the boundary either way. The
// float branch of the template constructor also had to stop using placement new, which is not
// allowed in a constant expression until C++26.
TEST(float128, ConstexprFloatingPointConversion)
{
    constexpr float128 from_double = 3.25;
    constexpr float128 from_float = 2.5f;
    constexpr float128 negative = -2.5;

    static_assert(static_cast<double>(from_double) == 3.25);
    static_assert(static_cast<float>(from_double) == 3.25f);
    static_assert(static_cast<long double>(from_double) == 3.25L);
    static_assert(from_float == float128(2.5));
    static_assert(static_cast<double>(negative) == -2.5);
    static_assert(float128(0.5) == float128::half());
    static_assert(float128(1.0) == float128::one());
    static_assert(float128(0.0).is_zero() && float128(-0.0).is_zero());
    static_assert(from_double + negative == float128(0.75));

    // the whole double range round trips, including the extremes
    static_assert(static_cast<double>(float128(1e300)) == 1e300);
    static_assert(static_cast<double>(float128(-1e-300)) == -1e-300);
    constexpr double tenth = static_cast<double>(float128(0.1));
    static_assert(tenth == 0.1);

    // the special values reach a double as HUGE_VAL and NAN
    static_assert(static_cast<double>(float128::inf()) > 1e308);
    static_assert(static_cast<double>(-float128::inf()) < -1e308);
    constexpr double quiet = static_cast<double>(float128::nan());
    static_assert(quiet != quiet);  // only a NaN compares unequal to itself

    EXPECT_TRUE(float128(static_cast<double>(opaque(13)) / 4.0) == from_double);
    EXPECT_DOUBLE_EQ(static_cast<double>(float128(opaque(1)) / 4), 0.25);
    EXPECT_TRUE(float128(0.1) == float128(tenth));
}

/**********************************************************************
 * Standard library surface
 *
 * numeric_limits and the <cmath> functions that had no float128 counterpart. How accurately the
 * ones that compute something do so is checked against correctly rounded references in
 * float128_accuracy_gtest.cpp. What is checked here is that they exist, that the lookup a caller
 * would use finds them, and that they behave as specified at the arguments where a formula on its
 * own gives the wrong answer.
 ***********************************************************************/
TEST(float128, NumericLimits)
{
    using limits = std::numeric_limits<float128>;

    static_assert(limits::is_specialized && limits::is_signed && limits::is_iec559);
    static_assert(!limits::is_integer && !limits::is_exact && !limits::is_modulo);
    static_assert(limits::digits == 113 && limits::digits10 == 33 && limits::max_digits10 == 36);
    static_assert(limits::radix == 2);
    static_assert(limits::min_exponent == -16381 && limits::max_exponent == 16384);
    static_assert(limits::min_exponent10 == -4931 && limits::max_exponent10 == 4932);
    static_assert(limits::has_infinity && limits::has_quiet_NaN && limits::has_signaling_NaN);
    static_assert(limits::round_style == std::round_to_nearest);
    // a cv qualified type has to answer the same way
    static_assert(std::numeric_limits<const float128>::digits == limits::digits);

    // epsilon is the gap above one, so adding it is visible and adding half of it is not
    const float128 one = float128::one();
    EXPECT_TRUE(one + limits::epsilon() > one);
    EXPECT_TRUE(one + (limits::epsilon() >> 1) == one);

    // the extremes are finite, and stepping past either one leaves the finite range
    EXPECT_TRUE(isfinite(limits::max()) && isfinite(limits::lowest()));
    EXPECT_TRUE(limits::lowest() == -limits::max());
    EXPECT_TRUE(isinf(nextafter(limits::max(), limits::infinity())));

    // min is the smallest normal, denorm_min the smallest value of any kind
    EXPECT_TRUE(isnormal(limits::min()));
    EXPECT_FALSE(isnormal(nextafter(limits::min(), float128())));
    EXPECT_TRUE(limits::denorm_min().is_subnormal());
    EXPECT_TRUE(nextafter(float128(), limits::infinity()) == limits::denorm_min());

    EXPECT_TRUE(isnan(limits::quiet_NaN()) && !limits::quiet_NaN().is_signaling());
    EXPECT_TRUE(isnan(limits::signaling_NaN()) && limits::signaling_NaN().is_signaling());
    EXPECT_TRUE(isinf(limits::infinity()) && limits::infinity().is_positive());
    EXPECT_TRUE(limits::round_error() == float128::half());
}
TEST(float128, Classification)
{
    const float128 one = float128::one();
    const float128 denorm = std::numeric_limits<float128>::denorm_min();

    EXPECT_EQ(fpclassify(float128()), FP_ZERO);
    EXPECT_EQ(fpclassify(-float128()), FP_ZERO);
    EXPECT_EQ(fpclassify(one), FP_NORMAL);
    EXPECT_EQ(fpclassify(denorm), FP_SUBNORMAL);
    EXPECT_EQ(fpclassify(float128::inf()), FP_INFINITE);
    EXPECT_EQ(fpclassify(float128::nan()), FP_NAN);

    // signbit reads the bit rather than comparing, so it sees a negative zero and a negative NaN
    EXPECT_FALSE(signbit(one));
    EXPECT_TRUE(signbit(-one));
    EXPECT_FALSE(signbit(float128()));
    EXPECT_TRUE(signbit(-float128()));
    EXPECT_TRUE(signbit(-float128::nan()));

    EXPECT_TRUE(isnormal(one));
    EXPECT_FALSE(isnormal(float128()));
    EXPECT_FALSE(isnormal(denorm));
    EXPECT_FALSE(isnormal(float128::inf()));
}
TEST(float128, ComparisonPredicates)
{
    const float128 one = float128::one();
    const float128 two(2);
    const float128 quiet = float128::nan();

    EXPECT_TRUE(isgreater(two, one) && !isgreater(one, two) && !isgreater(one, one));
    EXPECT_TRUE(isgreaterequal(one, one) && isgreaterequal(two, one));
    EXPECT_TRUE(isless(one, two) && !isless(two, one));
    EXPECT_TRUE(islessequal(one, one) && islessequal(one, two));
    EXPECT_TRUE(islessgreater(one, two) && !islessgreater(one, one));

    // every one of them is false against a NaN, and only isunordered says so
    EXPECT_TRUE(isunordered(quiet, one) && isunordered(one, quiet) && isunordered(quiet, quiet));
    EXPECT_FALSE(isunordered(one, two));
    EXPECT_FALSE(isgreater(quiet, one));
    EXPECT_FALSE(isgreaterequal(quiet, quiet));
    EXPECT_FALSE(isless(quiet, one));
    EXPECT_FALSE(islessequal(quiet, quiet));
    EXPECT_FALSE(islessgreater(quiet, one));
}
TEST(float128, MinMaxDim)
{
    const float128 one = float128::one();
    const float128 two(2);
    const float128 quiet = float128::nan();

    EXPECT_TRUE(fmin(one, two) == one && fmax(one, two) == two);
    // a NaN operand is missing data rather than a value, so the other operand is the answer
    EXPECT_TRUE(fmin(quiet, two) == two && fmin(two, quiet) == two);
    EXPECT_TRUE(fmax(quiet, two) == two && fmax(two, quiet) == two);
    EXPECT_TRUE(isnan(fmin(quiet, quiet)) && isnan(fmax(quiet, quiet)));
    // the zeros compare equal, so which one comes back has to be decided on the sign
    EXPECT_TRUE(signbit(fmin(float128(), -float128())));
    EXPECT_FALSE(signbit(fmax(float128(), -float128())));

    EXPECT_TRUE(fdim(two, one) == one);
    EXPECT_TRUE(fdim(one, two).is_zero());
    EXPECT_TRUE(isnan(fdim(quiet, one)) && isnan(fdim(one, quiet)));
}
TEST(float128, NanPayload)
{
    // Qualified, unlike every other function here: argument dependent lookup cannot reach this
    // one from a character pointer, so with <cmath> also in scope the unqualified name is
    // ambiguous against the ::nan that returns a double.
    EXPECT_TRUE(isnan(fp128::nan("")));
    EXPECT_FALSE(fp128::nan("").is_signaling());
    EXPECT_TRUE(isnan(fp128::nan("7")));
    EXPECT_TRUE(isnan(fp128::nan("0x1234")));

    // the payload lands in the fraction, so two different ones are different encodings
    uint64_t low_a = 0, high_a = 0, low_b = 0, high_b = 0;
    fp128::nan("7").get_bits(low_a, high_a);
    fp128::nan("9").get_bits(low_b, high_b);
    EXPECT_NE(low_a, low_b);
    EXPECT_EQ(high_a, high_b);
}
TEST(float128, NextAfter)
{
    const float128 one = float128::one();
    const float128 up = nextafter(one, float128::inf());
    const float128 down = nextafter(one, float128());

    EXPECT_TRUE(up > one && down < one);
    EXPECT_TRUE(up - one == std::numeric_limits<float128>::epsilon());
    // the gap below one is half the gap above it, one binade down
    EXPECT_TRUE(one - down == (std::numeric_limits<float128>::epsilon() >> 1));
    EXPECT_TRUE(nextafter(one, one) == one);
    EXPECT_TRUE(nextafter(down, float128::inf()) == one);
    EXPECT_TRUE(isnan(nextafter(float128::nan(), one)));
    EXPECT_TRUE(nexttoward(one, float128::inf()) == up);

    // zero steps to the smallest subnormal in either direction
    EXPECT_TRUE(nextafter(float128(), one) == std::numeric_limits<float128>::denorm_min());
    EXPECT_TRUE(nextafter(float128(), -one) == -std::numeric_limits<float128>::denorm_min());
}
TEST(float128, RintTiesToEven)
{
    const float128 half = float128::half();

    // round() takes a tie away from zero, rint() takes it to the even neighbour
    EXPECT_TRUE(rint(float128(2) + half) == float128(2));
    EXPECT_TRUE(rint(float128(3) + half) == float128(4));
    EXPECT_TRUE(rint(-(float128(2) + half)) == float128(-2));
    EXPECT_TRUE(rint(-(float128(3) + half)) == float128(-4));
    EXPECT_TRUE(round(float128(2) + half) == float128(3));

    EXPECT_TRUE(rint(float128(2)) == float128(2));
    EXPECT_TRUE(rint(float128("2.4")) == float128(2));
    EXPECT_TRUE(rint(float128("2.6")) == float128(3));
    EXPECT_TRUE(nearbyint(float128(3) + half) == float128(4));
    EXPECT_TRUE(isnan(rint(float128::nan())) && isinf(rint(float128::inf())));

    EXPECT_EQ(llrint(float128(2) + half), 2);
    EXPECT_EQ(llrint(float128(3) + half), 4);
    EXPECT_EQ(lrint(float128(3) + half), 4);
    EXPECT_EQ(llround(float128(2) + half), 3);
}
TEST(float128, ScalbnAndIlogb)
{
    const float128 one = float128::one();

    EXPECT_TRUE(scalbn(one, 3) == float128(8));
    EXPECT_TRUE(scalbn(float128(8), -3) == one);
    EXPECT_TRUE(scalbln(one, 3L) == float128(8));
    // an exponent far outside the range saturates rather than wrapping
    EXPECT_TRUE(isinf(scalbln(one, 1000000L)));
    EXPECT_TRUE(scalbln(one, -1000000L).is_zero());

    EXPECT_EQ(ilogb(one), 0);
    EXPECT_EQ(ilogb(float128(8)), 3);
    EXPECT_EQ(ilogb(float128::half()), -1);
    // a subnormal has an exponent, just not the one its exponent field holds
    EXPECT_EQ(ilogb(std::numeric_limits<float128>::denorm_min()), -16494);
    EXPECT_EQ(ilogb(float128()), FP_ILOGB0);
    EXPECT_EQ(ilogb(float128::nan()), FP_ILOGBNAN);
    EXPECT_EQ(ilogb(float128::inf()), INT_MAX);

    EXPECT_TRUE(logb(float128(8)) == float128(3));
    EXPECT_TRUE(logb(std::numeric_limits<float128>::denorm_min()) == float128(-16494));
    EXPECT_TRUE(isinf(logb(float128())) && logb(float128()).is_negative());
}
TEST(float128, RemainderAndRemquo)
{
    const float128 five(5), three(3), two(2), one = float128::one();

    // fmod keeps the sign of the dividend, remainder lands in [-|y|/2, |y|/2]
    EXPECT_TRUE(fmod(five, three) == two);
    EXPECT_TRUE(remainder(five, three) == -one);
    EXPECT_TRUE(fmod(-five, three) == -two);
    EXPECT_TRUE(remainder(-five, three) == one);

    int quotient = 0;
    EXPECT_TRUE(remquo(five, three, &quotient) == -one);
    EXPECT_EQ(quotient, 2);
    EXPECT_TRUE(remquo(-five, three, &quotient) == one);
    EXPECT_EQ(quotient, -2);

    // a tie in the quotient goes to the even one: 3/2 rounds to 2, 1/2 rounds to 0
    EXPECT_TRUE(remainder(three, two) == -one);
    EXPECT_TRUE(remainder(one, two) == one);

    EXPECT_TRUE(isnan(remainder(five, float128())));
    EXPECT_TRUE(isnan(remainder(float128::inf(), three)));
    EXPECT_TRUE(remainder(five, float128::inf()) == five);
}
TEST(float128, Hypot3)
{
    EXPECT_TRUE(hypot(float128(3), float128(4)) == float128(5));
    EXPECT_TRUE(hypot(float128(2), float128(3), float128(6)) == float128(7));
    // squaring these would overflow, the scaling inside keeps them in range
    EXPECT_TRUE(hypot(ldexp(float128(3), 16000), ldexp(float128(4), 16000), float128()) == ldexp(float128(5), 16000));
    EXPECT_TRUE(isinf(hypot(float128::inf(), float128::nan(), float128())));
    EXPECT_TRUE(isnan(hypot(float128::nan(), float128(1), float128(1))));
}
TEST(float128, Gamma)
{
    // the factorials are exact up to 34!
    EXPECT_TRUE(tgamma(float128(1)) == float128::one());
    EXPECT_TRUE(tgamma(float128(5)) == float128(24));
    EXPECT_TRUE(tgamma(float128(11)) == float128(3628800));
    // lgamma(1) and lgamma(2) are exactly zero rather than merely small
    EXPECT_TRUE(lgamma(float128::one()).is_zero());
    EXPECT_TRUE(lgamma(float128(2)).is_zero());

    // gamma(1/2) is sqrt(pi)
    const float128 root_pi = tgamma(float128::half());
    EXPECT_TRUE(fabs(sqr(root_pi) - float128::pi()) < ldexp(float128::one(), -100));

    // the poles, and the alternating sign between them
    EXPECT_TRUE(isinf(tgamma(float128())) && tgamma(float128()).is_positive());
    EXPECT_TRUE(isnan(tgamma(float128(-3))));
    EXPECT_TRUE(tgamma(float128("-0.5")).is_negative());
    EXPECT_TRUE(tgamma(float128("-1.5")).is_positive());
    EXPECT_TRUE(tgamma(float128("-2.5")).is_negative());
    EXPECT_TRUE(isinf(lgamma(float128())) && isinf(lgamma(float128(-2))));
    EXPECT_TRUE(isinf(tgamma(float128(2000))));
}
TEST(float128, SubnormalArithmetic)
{
    using limits = std::numeric_limits<float128>;
    const float128 denorm_min = limits::denorm_min();
    const float128 one = float128::one();

    // Scaling a subnormal all the way back up has to land on exactly one. get_components() used to
    // read a subnormal's exponent as if the stored field meant something, which put the smallest
    // one 223 binades off and made every operation that touched it return an unrelated value.
    EXPECT_TRUE(ldexp(denorm_min, 16494) == one);
    EXPECT_TRUE(ldexp(limits::min(), 16382) == one);
    EXPECT_TRUE((denorm_min << 16494) == one);

    // arithmetic inside the subnormal range stays exact, it is a fixed point grid
    EXPECT_TRUE(denorm_min + denorm_min == float128(2, 0));
    EXPECT_TRUE(float128(3, 0) - denorm_min == float128(2, 0));
    EXPECT_TRUE((denorm_min << 1) == float128(2, 0));

    // and the boundary between the two ranges is continuous in both directions
    EXPECT_TRUE(nextafter(limits::min(), float128()) == float128(UINT64_MAX, 0x0000FFFFFFFFFFFFull));
    EXPECT_TRUE(nextafter(nextafter(limits::min(), float128()), one) == limits::min());
    EXPECT_TRUE(limits::min() - denorm_min < limits::min());
    EXPECT_TRUE((limits::min() >> 1) == float128(0, 0x0000800000000000ull));

    // a subnormal converts to double as zero, and to itself through the components
    uint64_t low = 0, high = 0;
    int32_t expo = 0;
    uint32_t sign = 0;
    denorm_min.get_components(low, high, expo, sign);
    EXPECT_EQ(expo, -16494);
    float128 rebuilt;
    rebuilt.set_components(low, high, expo, sign);
    EXPECT_TRUE(rebuilt == denorm_min);
}
TEST(float128, AbsAlias)
{
    EXPECT_TRUE(abs(float128(-3)) == float128(3));
    EXPECT_TRUE(abs(float128(3)) == float128(3));
    EXPECT_FALSE(signbit(abs(-float128())));
}
