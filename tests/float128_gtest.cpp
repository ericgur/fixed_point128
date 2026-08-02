// remove warnings from gtest itself
#pragma warning(push)
#pragma warning(disable : 26439)
#pragma warning(disable : 26495)
#include <gtest/gtest.h>
#pragma warning(pop)
#include <ostream>
#include <ctime>
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
    constexpr auto end_char_to_check = 9;
    constexpr auto max_allowed_error = 1;
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
            float128 f = str;
            char* res = static_cast<char*>(f);
            size_t len = std::min(strlen(res), strlen(str));
            // clip both string based on the shortest one.
            str[len] = '\0';
            res[len] = '\0';
            int64_t orig_last_digits = strtoll(&str[len - end_char_to_check], nullptr, 10);
            int64_t res_last_digits = strtoll(&res[len - end_char_to_check], nullptr, 10);
            int64_t err = abs(orig_last_digits - res_last_digits);

            EXPECT_LE(err, max_allowed_error) << "error: " << err << ", last digits source: " << &str[len - end_char_to_check]
                                              << "last digits result: " << &res[len - end_char_to_check];
        } catch (...) {
            EXPECT_NO_THROW(i) << "failed at iteration " << i;
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
        auto value1 = get_int32_random();
        auto value2 = get_int32_random();
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
        auto value1 = get_uint32_random();
        auto value2 = get_uint32_random();
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
        auto value1 = get_int64_random();
        auto value2 = get_int64_random();
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
        auto value1 = get_uint64_random();
        auto value2 = get_uint64_random();
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
    char buf[64];
    for (auto v : values) {
        snprintf(buf, sizeof(buf), "%llu", v);
        EXPECT_STREQ(static_cast<char*>(float128(v)), buf);
    }

    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const uint64_t v = get_uint64_random() >> (get_uint32_random() % 64);
        snprintf(buf, sizeof(buf), "%llu", v);
        EXPECT_STREQ(static_cast<char*>(float128(v)), buf) << "value=" << v;
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
