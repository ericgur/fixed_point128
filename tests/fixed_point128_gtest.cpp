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
#include "gtest_shared.h"

/**********************************************************************
 * fixed_point128 tests
 ***********************************************************************/
// Construct fixed_point128 and convert back to/from various elements.
TEST(fixed_point128, DefaultConstructor)
{
    fixed_point128<20> f;
    EXPECT_EQ(static_cast<uint64_t>(f), 0ull);
}
TEST(fixed_point128, ConstructorFromDouble)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random(-43, 31);
        fixed_point128<32> f = value;
        double f_value = static_cast<double>(f);
        if (fabs(value / f_value) - 1.0 > 0.0001) {
            f_value = value;
        }
        EXPECT_DOUBLE_EQ(f_value, value) << "value=" << value;
    }
}
TEST(fixed_point128, ConstructorFromFloat)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        float value = (float)get_double_random(-43, 31);
        fixed_point128<32> f = value;
        EXPECT_FLOAT_EQ(static_cast<float>(f), value);
    }
}
TEST(fixed_point128, ConstructorFromInt32)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        int32_t value = get_int32_random();
        fixed_point128<32> f = value;
        if (check_overflow(value, f)) {
            continue;
        }
        EXPECT_EQ(static_cast<int32_t>(f), value);
    }
}
TEST(fixed_point128, ConstructorFromUnsignedInt32)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        uint32_t value = get_uint32_random();
        fixed_point128<32> f = value;
        if (check_overflow(value, f)) {
            continue;
        }
        EXPECT_EQ(static_cast<uint32_t>(f), value);
    }
}
TEST(fixed_point128, ConstructorFromInt64)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        int64_t value = get_int64_random();
        fixed_point128<32> f = value;
        if (check_overflow(value, f)) {
            continue;
        }
        EXPECT_EQ(static_cast<int64_t>(f), value);
    }
}
TEST(fixed_point128, ConstructorFromUnsignedInt64)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        uint64_t value = get_uint64_random();
        fixed_point128<32> f = value;
        if (check_overflow(value, f)) {
            continue;
        }
        EXPECT_EQ(static_cast<uint64_t>(f), value);
    }
}
TEST(fixed_point128, ConstructorFromString)
{
    const char* values[] = {"0", "12.34", "-3.1435", "12345.6789"};
    for (auto i = 0u; i < array_length(values); ++i) {
        fixed_point128<20> f = values[i];
        double d1 = strtod(values[i], nullptr);
        double d2 = strtod(static_cast<char*>(f), nullptr);
        EXPECT_DOUBLE_EQ(d1, d2) << "value=" << values[i];
    }
}
TEST(fixed_point128, CopyConstructor)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random();
        fixed_point128<20> f1 = value;
        fixed_point128<20> f2 = f1;
        EXPECT_DOUBLE_EQ(static_cast<double>(f1), static_cast<double>(f2)) << "value=" << value;
        ;
    }
}
TEST(fixed_point128, MoveConstructor)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random();
        fixed_point128<20> f1 = value;
        fixed_point128<20> f2(std::move(f1));
        EXPECT_DOUBLE_EQ(static_cast<double>(f1), static_cast<double>(f2)) << "value=" << value;
        ;
    }
}
TEST(fixed_point128, AssignmentOperator)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random();
        fixed_point128<20> f1 = value;
        fixed_point128<20> f2;
        f2 = f1;
        EXPECT_DOUBLE_EQ(static_cast<double>(f1), static_cast<double>(f2)) << "value=" << value;
        ;
    }
}
TEST(fixed_point128, AssignmentOperatorOtherType)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random(-50, 18);
        fixed_point128<20> f1 = value;
        fixed_point128<22> f2;
        f2 = f1;
        EXPECT_DOUBLE_EQ(static_cast<double>(f1), static_cast<double>(f2)) << "value=" << value;
        ;
        f1 = f2;
        EXPECT_DOUBLE_EQ(static_cast<double>(f1), static_cast<double>(f2)) << "value=" << value;
        ;
    }
}
TEST(fixed_point128, MoveAssignmentOperator)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random(-50, 18);
        fixed_point128<20> f1 = value;
        fixed_point128<20> f2;
        f2 = std::move(f1);
        EXPECT_DOUBLE_EQ(static_cast<double>(f1), static_cast<double>(f2)) << "value=" << value;
        ;
    }
}
TEST(fixed_point128, CopyConstructorOtherType)
{
    double value = 0;
    bool common_value = false;
    constexpr int f1_prec = 64, f2_prec= 32; 
    constexpr int min_int_precision = std::min(f1_prec, f2_prec);
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        while (!common_value) {
            value = get_double_random(); 
            common_value = fabs(value) < fixed_point128<min_int_precision>::max_int_value;
        }
        fixed_point128<f1_prec> f1 = value;
        fixed_point128<f2_prec> f2 = f1;
        EXPECT_DOUBLE_EQ(static_cast<double>(f1), static_cast<double>(f2)) << "value=" << value;
        fixed_point128<f1_prec> f3 = f2;
        EXPECT_DOUBLE_EQ(static_cast<double>(f2), static_cast<double>(f3)) << "value=" << value;
    }
}
TEST(fixed_point128, AddSameSign)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = fabs(get_double_random());
        double value2 = value1 * 2.5;
        double res = value1 + value2;
        fixed_point128<40> f1 = value1;
        fixed_point128<40> f2 = value2;
        fixed_point128<40> f3 = f1 + f2;
        if (check_overflow(value1, f1) || check_overflow(value2, f2) || check_overflow(res, f3)) {
            continue;
        }
        // The operands are exactly representable and so is their sum, so the result matches the
        // double reference bit for bit. The is_similar_double() guard that used to sit here was
        // satisfied by every iteration, which meant the assertion below never ran at all.
        EXPECT_DOUBLE_EQ(static_cast<double>(f3), res) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(fixed_point128, AddDifferentSign)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = fabs(get_double_random(-10, 35));
        double value2 = value1 * -2.5;
        double res = value1 + value2;
        fixed_point128<40> f1 = value1;
        fixed_point128<40> f2 = value2;
        fixed_point128<40> f3 = f1 + f2;
        if (check_overflow(value1, f1) || check_overflow(value2, f2) || check_overflow(res, f3)) {
            continue;
        }
        // The operands are exactly representable and so is their sum, so the result matches the
        // double reference bit for bit. The is_similar_double() guard that used to sit here was
        // satisfied by every iteration, which meant the assertion below never ran at all.
        EXPECT_DOUBLE_EQ(static_cast<double>(f3), res) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(fixed_point128, AddDouble)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = get_double_random();
        double value2 = get_double_random();
        double res = value1 + value2;
        fixed_point128<40> f1 = value1;
        fixed_point128<40> f3 = f1 + value2;
        if (check_overflow(value1, f1) || check_overflow(value2, f1) || check_overflow(res, f3)) {
            continue;
        }
        // exactly representable operands and an exactly representable sum, see AddSameSign
        EXPECT_DOUBLE_EQ(static_cast<double>(f3), res) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(fixed_point128, AddFloat)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        float value1 = (float)get_double_random();
        float value2 = (float)get_double_random();
        float res = value1 + value2;
        fixed_point128<40> f1 = value1;
        fixed_point128<40> f3 = f1 + value2;
        if (check_overflow(value1, f1) || check_overflow(value2, f1) || check_overflow(res, f3)) {
            continue;
        }
        EXPECT_FLOAT_EQ(static_cast<float>(f3), res) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(fixed_point128, AddInt32)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const int32_t value1 = abs(get_int32_random());
        const int32_t value2 = abs(get_int32_random());
        // The sum of two positive 32 bit values can exceed the type. Detecting that after the
        // fact by testing for a negative result relies on signed overflow, which is undefined
        // behavior and may be optimized away. Compute the reference in a wider type instead.
        const int64_t res = static_cast<int64_t>(value1) + value2;

        fixed_point128<40> f1 = value1;
        fixed_point128<40> f3 = f1 + value2;
        if (check_overflow(value1, f1) || check_overflow(value2, f1) || check_overflow(res, f3)) {
            continue;
        }
        EXPECT_EQ(static_cast<int64_t>(f3), res) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(fixed_point128, AddUnsignedInt32)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        auto value1 = get_uint32_random();
        auto value2 = get_uint32_random();
        auto res = value1 + value2;
        if (res < value1 || res < value2)
            continue;  // wrap around won't happen in uint128_t

        fixed_point128<40> f1 = value1;
        fixed_point128<40> f3 = f1 + value2;
        if (check_overflow(value1, f1) || check_overflow(value2, f1) || check_overflow(res, f3)) {
            continue;
        }
        EXPECT_EQ(static_cast<uint64_t>(f3), res) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(fixed_point128, AddInt64)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        auto value1 = get_int64_random();
        auto value2 = get_int64_random();
        auto res = value1 + value2;
        fixed_point128<40> f1 = value1;
        fixed_point128<40> f3 = f1 + value2;
        if (check_overflow(value1, f1) || check_overflow(value2, f1) || check_overflow(res, f3)) {
            continue;
        }
        EXPECT_EQ(static_cast<int64_t>(f3), res) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(fixed_point128, AddUnsignedInt64)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        auto value1 = get_uint64_random();
        auto value2 = get_uint64_random();
        auto res = value1 + value2;
        fixed_point128<40> f1 = value1;
        fixed_point128<40> f3 = f1 + value2;
        if (check_overflow(value1, f1) || check_overflow(value2, f1) || check_overflow(res, f3)) {
            continue;
        }
        EXPECT_EQ(static_cast<uint64_t>(f3), res) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(fixed_point128, SubtractSameSign)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = fabs(get_double_random());
        double value2 = value1 * 2.5;
        double res = value1 - value2;
        fixed_point128<40> f1 = value1;
        fixed_point128<40> f2 = value2;
        fixed_point128<40> f3 = f1 - f2;
        if (check_overflow(value1, f1) || check_overflow(value2, f1) || check_overflow(res, f3)) {
            continue;
        }
        // exactly representable operands and an exactly representable difference, see AddSameSign
        EXPECT_DOUBLE_EQ(static_cast<double>(f3), res) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(fixed_point128, SubtractDifferentSign)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = fabs(get_double_random());
        double value2 = value1 * -2.5;
        double res = value1 - value2;
        fixed_point128<40> f1 = value1;
        fixed_point128<40> f2 = value2;
        fixed_point128<40> f3 = f1 - f2;
        if (check_overflow(value1, f1) || check_overflow(value2, f1) || check_overflow(res, f3)) {
            continue;
        }
        // exactly representable operands and an exactly representable difference, see AddSameSign
        EXPECT_DOUBLE_EQ(static_cast<double>(f3), res) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(fixed_point128, SubtractInt32)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        auto value1 = abs(get_int32_random());
        auto value2 = abs(get_int32_random());
        if (value2 > value1)
            value2 = value1 / 3;

        auto res = value1 - value2;

        fixed_point128<40> f1 = value1;
        fixed_point128<40> f3 = f1 - value2;
        if (check_overflow(value1, f1) || check_overflow(value2, f1) || check_overflow(res, f3)) {
            continue;
        }
        EXPECT_EQ(static_cast<int64_t>(f3), res) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(fixed_point128, SubtractUnsignedInt32)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        auto value1 = get_uint32_random();
        auto value2 = get_uint32_random();
        if (value1 < value2) {
            std::swap(value1, value2);
        }
        auto res = value1 - value2;
        fixed_point128<40> f1 = value1;
        fixed_point128<40> f3 = f1 - value2;
        if (check_overflow(value1, f1) || check_overflow(value2, f1) || check_overflow(res, f3)) {
            continue;
        }
        EXPECT_EQ(static_cast<uint64_t>(f3), res) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(fixed_point128, SubtractInt64)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        auto value1 = get_int64_random();
        auto value2 = get_int64_random();
        auto res = value1 - value2;
        fixed_point128<40> f1 = value1;
        fixed_point128<40> f3 = f1 - value2;
        if (check_overflow(value1, f1) || check_overflow(value2, f1) || check_overflow(res, f3)) {
            continue;
        }
        EXPECT_EQ(static_cast<int64_t>(f3), res) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(fixed_point128, SubtractUnsignedInt64)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        auto value1 = get_uint64_random();
        auto value2 = get_uint64_random();
        if (value1 < value2) {
            std::swap(value1, value2);
        }
        auto res = value1 - value2;
        fixed_point128<40> f1 = value1;
        fixed_point128<40> f3 = f1 - value2;
        if (check_overflow(value1, f1) || check_overflow(value2, f1) || check_overflow(res, f3)) {
            continue;
        }
        EXPECT_EQ(static_cast<uint64_t>(f3), res) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(fixed_point128, MultiplyByFP128)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = get_double_random();
        double value2 = get_double_random();
        double res = value1 * value2;
        fixed_point128<40> f1 = value1;
        fixed_point128<40> f2 = value2;
        fixed_point128<40> f3 = f1 * f2;
        if (check_overflow(value1, f1) || check_overflow(value2, f1) || check_overflow(res, f3)) {
            continue;
        }
        // The product is truncated to 88 fraction bits, which is far more than the 53 significant
        // bits a double keeps, so rounding the result back to double lands on the same value the
        // double multiply produced. Measured exact on every one of the 20955 iterations that reach
        // this line.
        EXPECT_DOUBLE_EQ(static_cast<double>(f3), res) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(fixed_point128, MultiplyByDouble)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = get_double_random();
        double value2 = get_double_random();
        double res = value1 * value2;
        fixed_point128<40> f1 = value1;
        fixed_point128<40> f3 = f1 * value2;
        if (check_overflow(value1, f1) || check_overflow(value2, f1) || check_overflow(res, f3)) {
            continue;
        }
        // see MultiplyByFP128, the product carries more bits than a double can hold
        EXPECT_DOUBLE_EQ(static_cast<double>(f3), res) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(fixed_point128, MultiplyByFloat)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        float value1 = (float)get_double_random();
        float value2 = (float)get_double_random();
        float res = value1 * value2;
        fixed_point128<40> f1 = value1;
        fixed_point128<40> f3 = f1 * value2;
        if (check_overflow(value1, f1) || check_overflow(value2, f1) || check_overflow(res, f3)) {
            continue;
        }
        // see MultiplyByFP128, the product carries far more bits than a float can hold
        EXPECT_FLOAT_EQ(static_cast<float>(f3), res) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(fixed_point128, MultiplyByInt32)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        auto value1 = get_int32_random();
        auto value2 = get_int32_random();
        int64_t res = (int64_t)value1 * (int64_t)value2;
        fixed_point128<40> f1 = value1;
        fixed_point128<40> f3 = f1 * value2;
        if (check_overflow(value1, f1) || check_overflow(value2, f1) || check_overflow(res, f3)) {
            continue;
        }
        EXPECT_EQ(static_cast<int64_t>(f3), res) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(fixed_point128, MultiplyByUnsignedInt32)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        auto value1 = get_uint32_random();
        auto value2 = get_uint32_random();
        uint64_t res = (uint64_t)value1 * (uint64_t)value2;
        fixed_point128<40> f1 = value1;
        fixed_point128<40> f3 = f1 * value2;
        if (check_overflow(value1, f1) || check_overflow(value2, f1) || check_overflow(res, f3)) {
            continue;
        }
        EXPECT_EQ(static_cast<uint64_t>(f3), res) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(fixed_point128, MultiplyByInt64)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        auto value1 = get_int64_random();
        auto value2 = get_int64_random();
        int64_t res = value1 * value2;
        fixed_point128<40> f1 = value1;
        fixed_point128<40> f3 = f1 * value2;
        if (check_overflow(value1, f1) || check_overflow(value2, f1) || check_overflow(res, f3)) {
            continue;
        }
        EXPECT_EQ(static_cast<int64_t>(f3), res) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(fixed_point128, MultiplyByUnsignedInt64)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        auto value1 = get_uint64_random();
        auto value2 = get_uint64_random();
        uint64_t res = value1 * value2;
        fixed_point128<40> f1 = value1;
        fixed_point128<40> f3 = f1 * value2;
        if (check_overflow(value1, f1) || check_overflow(value2, f1) || check_overflow(res, f3)) {
            continue;
        }
        EXPECT_EQ(static_cast<uint64_t>(f3), res) << "value1=" << value1 << ", value2=" << value2;
    }
}
/**
 * @brief Checks that sqr(x) is bit identical to x * x for a single template parameter.
 *
 * square() accumulates the same 256 bit product as operator*=, it only skips the
 * redundant second cross multiply, so the results must match exactly - including
 * rounding and sign. Raw bit patterns are used to cover the full value range instead
 * of going through double, which cannot represent every fixed_point128 value.
 *
 * @tparam I Number of integer bits passed to fixed_point128.
 */
template <int32_t I> static void CheckSqrMatchesMultiply()
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const uint64_t low = get_uint64_random();
        const uint64_t high = get_uint64_random();
        const uint32_t sign = (uint32_t)(rand() & 1);
        const fixed_point128<I> x(low, high, sign);
        const fixed_point128<I> viaMultiply = x * x;
        const fixed_point128<I> viaSqr = sqr(x);

        // operator== compares sign, high and low, so this is a bit exact comparison
        EXPECT_TRUE(viaSqr == viaMultiply) << "I=" << I << ", low=" << low << ", high=" << high << ", sign=" << sign
                                           << ", x*x=" << (std::string)viaMultiply << ", sqr(x)=" << (std::string)viaSqr;
        EXPECT_FALSE(viaSqr.is_negative()) << "I=" << I << ", a square must never be negative";
    }
}
TEST(fixed_point128, SqrMatchesMultiply)
{
    CheckSqrMatchesMultiply<1>();
    CheckSqrMatchesMultiply<10>();
    CheckSqrMatchesMultiply<20>();
    CheckSqrMatchesMultiply<40>();
    CheckSqrMatchesMultiply<64>();
}
TEST(fixed_point128, SqrEdgeCases)
{
    using fp = fixed_point128<20>;
    const fp values[] = {fp(0), fp::epsilon(), fp::one(), -fp::one(), fp::half(), fp::pi(), -fp::pi(),
                         fp::e(), fp::golden_ratio(), fp(0ull, ~0ull, 0), fp(~0ull, ~0ull, 0), fp(~0ull, ~0ull, 1)};
    for (const auto& x : values) {
        EXPECT_TRUE(sqr(x) == x * x) << "x=" << (std::string)x;
    }
    // squaring zero stays zero and keeps a positive sign
    EXPECT_TRUE(sqr(fp(0)) == fp(0));
    EXPECT_FALSE(sqr(-fp::one()).is_negative());
}
TEST(fixed_point128, DivideByFP128)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = get_double_random();
        double value2 = get_double_random();
        if (value2 == 0)
            continue;
        double res = value1 / value2;
        fixed_point128<40> f1 = value1;
        fixed_point128<40> f2 = value2;
        fixed_point128<40> f3 = f1 / f2;
        if (check_overflow(value1, f1) || check_overflow(value2, f1) || check_overflow(res, f3)) {
            continue;
        }
        // Division is the one operation whose result does not match the double reference exactly:
        // the quotient is quantized to the type's last fraction bit, which for a small quotient
        // leaves fewer significant bits than a double carries. See FixedPointDivisionTolerance().
        const double fp128_res = static_cast<double>(f3);
        EXPECT_NEAR(fp128_res, res, FixedPointDivisionTolerance<40>(res)) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(fixed_point128, DivideByDouble)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = get_double_random();
        double value2 = get_double_random();
        double res = value1 / value2;
        if (value2 == 0)
            continue;
        fixed_point128<40> f1 = value1;
        fixed_point128<40> f3 = f1 / value2;
        if (check_overflow(value1, f1) || check_overflow(value2, f1) || check_overflow(res, f3)) {
            continue;
        }
        // see DivideByFP128 for why division needs a tolerance and the others do not
        EXPECT_NEAR(static_cast<double>(f3), res, FixedPointDivisionTolerance<40>(res)) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(fixed_point128, DivideByFloat)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        float value1 = (float)get_double_random();
        float value2 = (float)get_double_random();
        if (value2 == 0)
            continue;
        float res = value1 / value2;
        fixed_point128<40> f1 = value1;
        fixed_point128<40> f3 = f1 / value2;
        if (check_overflow(value1, f1) || check_overflow(value2, f1) || check_overflow(res, f3))
            continue;
        EXPECT_FLOAT_EQ(static_cast<float>(f3), res) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(fixed_point128, DivideByInt32)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = get_double_random(-1, 39);
        int32_t value2 = get_int32_random();
        double res = value1 / value2;
        if (value2 == 0)
            continue;
        fixed_point128<40> f1 = value1;
        fixed_point128<40> f3 = f1 / value2;
        if (check_overflow(value1, f1) || check_overflow(value2, f1) || check_overflow(res, f3)) {
            continue;
        }
        EXPECT_DOUBLE_EQ(static_cast<double>(f3), res) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(fixed_point128, DivideByUnsignedInt32)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = get_double_random(-1, 39);
        uint32_t value2 = get_uint32_random();
        double res = value1 / value2;
        if (value2 == 0)
            continue;
        fixed_point128<40> f1 = value1;
        fixed_point128<40> f3 = f1 / value2;
        if (check_overflow(value1, f1) || check_overflow(value2, f1) || check_overflow(res, f3)) {
            continue;
        }
        EXPECT_DOUBLE_EQ(static_cast<double>(f3), res) << "value1=" << value1 << ", value2=" << value2;
    }
}

TEST(fixed_point128, DivideByInt64)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = get_double_random(-1, 39);
        int64_t value2 = get_int32_random();  // on purpose
        double res = value1 / value2;
        if (value2 == 0)
            continue;
        fixed_point128<40> f1 = value1;
        fixed_point128<40> f3 = f1 / value2;
        if (check_overflow(value1, f1) || check_overflow(value2, f1) || check_overflow(res, f3)) {
            continue;
        }
        EXPECT_DOUBLE_EQ(static_cast<double>(f3), res) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(fixed_point128, DivideByUnsignedInt64)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = get_double_random(-1, 39);
        uint64_t value2 = get_uint32_random();  // on purpose
        double res = value1 / value2;
        if (value2 == 0)
            continue;
        fixed_point128<40> f1 = value1;
        fixed_point128<40> f3 = f1 / value2;
        if (check_overflow(value1, f1) || check_overflow(value2, f1) || check_overflow(res, f3)) {
            continue;
        }
        EXPECT_DOUBLE_EQ(static_cast<double>(f3), res) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(fixed_point128, ModuloByFP128)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = get_double_random(-1, 39);
        double value2 = get_double_random(-1, 39);
        double res = fmod(value1, value2);
        if (value2 == 0)
            continue;
        fixed_point128<40> f1 = value1;
        fixed_point128<40> f2 = value2;
        fixed_point128<40> f3 = f1 % f2;
        // The modulo needs the quotient to be representable, not merely the two operands: it is
        // computed by dividing first, and a quotient beyond 2^40 overflows the integer part.
        if (check_overflow(value1, f1) || check_overflow(value2, f2) || check_overflow(value1 / value2, f3)) {
            continue;
        }
        // the remainder is exact, it is a difference of representable values rather than a quotient
        EXPECT_DOUBLE_EQ(static_cast<double>(f3), res) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(fixed_point128, ModuloByDouble)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = get_double_random(-1, 39);
        double value2 = get_double_random(-1, 39);
        double res = fmod(value1, value2);
        if (value2 == 0)
            continue;
        fixed_point128<40> f1 = value1;
        fixed_point128<40> f3 = f1 % value2;
        if (check_overflow(value1, f1) || check_overflow(value2, f1) || check_overflow(value1 / value2, f3)) {
            continue;
        }
        // the remainder is exact, see ModuloByFP128
        EXPECT_DOUBLE_EQ(static_cast<double>(f3), res) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(fixed_point128, ModuloByFloat)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        float value1 = (float)get_double_random(-1, 39);
        float value2 = (float)get_double_random(-1, 39);
        if (value2 == 0)
            continue;
        float res = fmodf(value1, value2);
        fixed_point128<40> f1 = value1;
        fixed_point128<40> f3 = f1 % value2;
        if (check_overflow(value1, f1) || check_overflow(value2, f1) || check_overflow(value1 / value2, f3)) {
            continue;
        }
        // the remainder is exact, see ModuloByFP128
        EXPECT_FLOAT_EQ(static_cast<float>(f3), res) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(fixed_point128, ModuloByInt32)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = get_double_random(-1, 39);
        int32_t value2 = get_int32_random();
        double res = fmod(value1, value2);
        if (value2 == 0)
            continue;
        fixed_point128<40> f1 = value1;
        fixed_point128<40> f3 = f1 % value2;
        if (check_overflow(value1, f1) || check_overflow(value2, f1) || check_overflow(res, f3)) {
            continue;
        }
        EXPECT_DOUBLE_EQ(static_cast<double>(f3), res) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(fixed_point128, ModuloByUnsignedInt32)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = get_double_random(-1, 39);
        uint32_t value2 = get_uint32_random();
        double res = fmod(value1, value2);
        if (value2 == 0)
            continue;
        fixed_point128<40> f1 = value1;
        fixed_point128<40> f3 = f1 % value2;
        if (check_overflow(value1, f1) || check_overflow(value2, f1) || check_overflow(res, f3)) {
            continue;
        }
        EXPECT_DOUBLE_EQ(static_cast<double>(f3), res) << "value1=" << value1 << ", value2=" << value2;
    }
}

TEST(fixed_point128, ModuloByInt64)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = get_double_random(-1, 39);
        int64_t value2 = get_int32_random();  // on purpose
        double res = fmod(value1, value2);
        if (value2 == 0)
            continue;
        fixed_point128<40> f1 = value1;
        fixed_point128<40> f3 = f1 % value2;
        if (check_overflow(value1, f1) || check_overflow(value2, f1) || check_overflow(res, f3)) {
            continue;
        }
        EXPECT_DOUBLE_EQ(static_cast<double>(f3), res) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(fixed_point128, ModuloByUnsignedInt64)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = get_double_random(-1, 39);
        uint64_t value2 = get_uint32_random();  // on purpose
        double res = fmod(value1, value2);
        if (value2 == 0)
            continue;
        fixed_point128<40> f1 = value1;
        fixed_point128<40> f3 = f1 % value2;
        if (check_overflow(value1, f1) || check_overflow(value2, f1) || check_overflow(res, f3)) {
            continue;
        }
        EXPECT_DOUBLE_EQ(static_cast<double>(f3), res) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(fixed_point128, CompareFP128)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = get_double_random(-1, 39);
        double value2 = get_double_random(-1, 39);
        bool res = value1 > value2;
        fixed_point128<40> f1 = value1;
        fixed_point128<40> f2 = value2;
        if (check_overflow(value1, f1) || check_overflow(value2, f1))
            continue;

        bool fp128_res = f1 > f2;
        EXPECT_TRUE(fp128_res == res) << "operator>: "
                                      << "value1=" << value1 << ", value2=" << value2;

        res = value1 >= value2;
        fp128_res = f1 >= f2;
        EXPECT_TRUE(fp128_res == res) << "operator>=: "
                                      << "value1=" << value1 << ", value2=" << value2;

        res = value1 < value2;
        fp128_res = f1 < f2;
        EXPECT_TRUE(fp128_res == res) << "operator<: "
                                      << "value1=" << value1 << ", value2=" << value2;

        res = value1 <= value2;
        fp128_res = f1 <= f2;
        EXPECT_TRUE(fp128_res == res) << "operator<=: "
                                      << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(fixed_point128, CompareDouble)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = get_double_random();
        double value2 = get_double_random();
        bool res = value1 > value2;
        fixed_point128<40> f1 = value1;
        if (check_overflow(value1, f1) || check_overflow(value2, f1))
            continue;

        bool fp128_res = f1 > value2;
        EXPECT_TRUE(fp128_res == res) << "operator>: "
                                      << "value1=" << value1 << ", value2=" << value2;

        res = value1 >= value2;
        fp128_res = f1 >= value2;
        EXPECT_TRUE(fp128_res == res) << "operator>=: "
                                      << "value1=" << value1 << ", value2=" << value2;

        res = value1 < value2;
        fp128_res = f1 < value2;
        EXPECT_TRUE(fp128_res == res) << "operator<: "
                                      << "value1=" << value1 << ", value2=" << value2;

        res = value1 <= value2;
        fp128_res = f1 <= value2;
        EXPECT_TRUE(fp128_res == res) << "operator<=: "
                                      << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(fixed_point128, CompareFloat)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        float value1 = (float)get_double_random();
        float value2 = (float)get_double_random();
        bool res = value1 > value2;
        fixed_point128<40> f1 = value1;
        if (check_overflow(value1, f1) || check_overflow(value2, f1))
            continue;

        bool fp128_res = f1 > value2;
        EXPECT_TRUE(fp128_res == res) << "operator>: "
                                      << "value1=" << value1 << ", value2=" << value2;

        res = value1 >= value2;
        fp128_res = f1 >= value2;
        EXPECT_TRUE(fp128_res == res) << "operator>=: "
                                      << "value1=" << value1 << ", value2=" << value2;

        res = value1 < value2;
        fp128_res = f1 < value2;
        EXPECT_TRUE(fp128_res == res) << "operator<: "
                                      << "value1=" << value1 << ", value2=" << value2;

        res = value1 <= value2;
        fp128_res = f1 <= value2;
        EXPECT_TRUE(fp128_res == res) << "operator<=: "
                                      << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(fixed_point128, CompareInt32)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        int32_t value1 = get_int32_random();
        int32_t value2 = get_int32_random();
        bool res = value1 > value2;
        fixed_point128<40> f1 = value1;
        if (check_overflow(value1, f1) || check_overflow(value2, f1))
            continue;

        bool fp128_res = f1 > value2;
        EXPECT_TRUE(fp128_res == res) << "operator>: "
                                      << "value1=" << value1 << ", value2=" << value2;

        res = value1 >= value2;
        fp128_res = f1 >= value2;
        EXPECT_TRUE(fp128_res == res) << "operator>=: "
                                      << "value1=" << value1 << ", value2=" << value2;

        res = value1 < value2;
        fp128_res = f1 < value2;
        EXPECT_TRUE(fp128_res == res) << "operator<: "
                                      << "value1=" << value1 << ", value2=" << value2;

        res = value1 <= value2;
        fp128_res = f1 <= value2;
        EXPECT_TRUE(fp128_res == res) << "operator<=: "
                                      << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(fixed_point128, CompareUnsignedInt32)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        uint32_t value1 = get_uint32_random();
        uint32_t value2 = get_uint32_random();
        bool res = value1 > value2;
        fixed_point128<40> f1 = value1;
        if (check_overflow(value1, f1) || check_overflow(value2, f1))
            continue;

        bool fp128_res = f1 > value2;
        EXPECT_TRUE(fp128_res == res) << "operator>: "
                                      << "value1=" << value1 << ", value2=" << value2;

        res = value1 >= value2;
        fp128_res = f1 >= value2;
        EXPECT_TRUE(fp128_res == res) << "operator>=: "
                                      << "value1=" << value1 << ", value2=" << value2;

        res = value1 < value2;
        fp128_res = f1 < value2;
        EXPECT_TRUE(fp128_res == res) << "operator<: "
                                      << "value1=" << value1 << ", value2=" << value2;

        res = value1 <= value2;
        fp128_res = f1 <= value2;
        EXPECT_TRUE(fp128_res == res) << "operator<=: "
                                      << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(fixed_point128, CompareInt64)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        int64_t value1 = get_int64_random();
        int64_t value2 = get_int64_random();
        bool res = value1 > value2;
        fixed_point128<40> f1 = value1;
        if (check_overflow(value1, f1) || check_overflow(value2, f1))
            continue;

        bool fp128_res = f1 > value2;
        EXPECT_TRUE(fp128_res == res) << "operator>: "
                                      << "value1=" << value1 << ", value2=" << value2;

        res = value1 >= value2;
        fp128_res = f1 >= value2;
        EXPECT_TRUE(fp128_res == res) << "operator>=: "
                                      << "value1=" << value1 << ", value2=" << value2;

        res = value1 < value2;
        fp128_res = f1 < value2;
        EXPECT_TRUE(fp128_res == res) << "operator<: "
                                      << "value1=" << value1 << ", value2=" << value2;

        res = value1 <= value2;
        fp128_res = f1 <= value2;
        EXPECT_TRUE(fp128_res == res) << "operator<=: "
                                      << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(fixed_point128, CompareUnsignedInt64)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        uint64_t value1 = get_uint64_random();
        uint64_t value2 = get_uint64_random();
        bool res = value1 > value2;
        fixed_point128<40> f1 = value1;
        if (check_overflow(value1, f1) || check_overflow(value2, f1))
            continue;

        bool fp128_res = f1 > value2;
        EXPECT_TRUE(fp128_res == res) << "operator>: "
                                      << "value1=" << value1 << ", value2=" << value2;

        res = value1 >= value2;
        fp128_res = f1 >= value2;
        EXPECT_TRUE(fp128_res == res) << "operator>=: "
                                      << "value1=" << value1 << ", value2=" << value2;

        res = value1 < value2;
        fp128_res = f1 < value2;
        EXPECT_TRUE(fp128_res == res) << "operator<: "
                                      << "value1=" << value1 << ", value2=" << value2;

        res = value1 <= value2;
        fp128_res = f1 <= value2;
        EXPECT_TRUE(fp128_res == res) << "operator<=: "
                                      << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(fixed_point128, OperatorPlusPlus)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = get_double_random();
        double res = value1 + 1;
        fixed_point128<40> f1 = value1;
        fixed_point128<40> f2 = f1;
        f2++;

        if (check_overflow(value1 + 2, f1))
            continue;

        EXPECT_DOUBLE_EQ(static_cast<double>(f2), res) << "operator++(int)"
                                                       << "value1=" << value1;

        f2 = f1;
        ++f2;

        EXPECT_DOUBLE_EQ(static_cast<double>(f2), res) << "operator++()"
                                                       << "value1=" << value1;
    }
}
TEST(fixed_point128, OperatorMinusMinus)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = get_double_random();
        double res = value1 - 1;
        fixed_point128<40> f1 = value1;
        fixed_point128<40> f2 = f1;
        f2--;

        if (check_overflow(value1 + 2, f1) || check_overflow(value1 - 2, f1))
            continue;

        EXPECT_DOUBLE_EQ(static_cast<double>(f2), res) << "operator--(int)"
                                                       << "value1=" << value1;

        f2 = f1;
        --f2;

        EXPECT_DOUBLE_EQ(static_cast<double>(f2), res) << "operator--()"
                                                       << "value1=" << value1;
    }
}
TEST(fixed_point128, OperatorEqual)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = get_double_random();
        fixed_point128<40> f1 = value1;
        bool fp128_res = f1 == f1;

        if (check_overflow(value1, f1))
            continue;

        EXPECT_TRUE(fp128_res == true) << "operator==: "
                                       << "value1=" << value1;
    }
}
TEST(fixed_point128, TemplateOperatorEqual)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = get_double_random();
        fixed_point128<40> f1 = value1;
        bool fp128_res = f1 == value1;

        if (check_overflow(value1, f1))
            continue;

        EXPECT_TRUE(fp128_res == true) << "operator==<double>: "
                                       << "value1=" << value1;

        f1 = static_cast<float>(value1);
        fp128_res = f1 == static_cast<float>(value1);

        EXPECT_TRUE(fp128_res == true) << "operator==<float>: "
                                       << "value1=" << value1;

        f1 = static_cast<uint64_t>(value1);
        fp128_res = f1 == static_cast<uint64_t>(value1);

        EXPECT_TRUE(fp128_res == true) << "operator==<uint64_t>: "
                                       << "value1=" << value1;

        f1 = static_cast<int64_t>(value1);
        fp128_res = f1 == static_cast<int64_t>(value1);

        EXPECT_TRUE(fp128_res == true) << "operator==<int64_t>: "
                                       << "value1=" << value1;

        f1 = static_cast<uint32_t>(value1);
        fp128_res = f1 == static_cast<uint32_t>(value1);

        EXPECT_TRUE(fp128_res == true) << "operator==<uint32_t>: "
                                       << "value1=" << value1;

        f1 = static_cast<int32_t>(value1);
        fp128_res = f1 == static_cast<int32_t>(value1);

        EXPECT_TRUE(fp128_res == true) << "operator==<int32_t>: "
                                       << "value1=" << value1;
    }
}
TEST(fixed_point128, OperatorNotEqual)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = get_double_random();
        fixed_point128<40> f1 = value1;
        bool fp128_res = f1 != f1;

        if (check_overflow(value1, f1))
            continue;

        EXPECT_TRUE(fp128_res == false) << "operator!=: "
                                        << "value1=" << value1;
    }
}
/*
TEST(fixed_point128, OperatorCString)
{
    srand(RANDOM_SEED);
    constexpr uint64_t MAX_TEST_STR_LEN = 37;
    // Builds a random signed decimal string with fracDigits fractional digits.
    // The last fractional digit is always non-zero to match operator char*() which strips trailing zeros.
    // maxIntVal limits the integer part to avoid overflow for the given template parameter.
    auto genStr = [](char* buf, size_t bufSize, uint32_t maxIntVal, int fracDigits) {
        const char* sign = get_random_sign() < 0 ? "-" : "";
        uint32_t intPart = get_uint32_random() % (maxIntVal + 1);
        char frac[40] = {};
        for (int j = 0; j < fracDigits - 1; ++j) {
            frac[j] = get_digit_random();
        }
        frac[fracDigits - 1] = (char)('1' + rand() % 9);  // non-zero to avoid trailing zero mismatch
        frac[fracDigits] = '\0';
        snprintf(buf, bufSize, "%s%u.%s", sign, intPart, frac);
    };

    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        char buf[128];

        // fixed_point128<40>: I=40, F=88, max_frac_digits=26; integer part bounded to 999999
        genStr(buf, sizeof(buf), 999999u, 10);
        {
            fixed_point128<40> f40 = buf;
            const char* out = static_cast<char*>(f40);
            size_t cmpLen = std::min(strlen(out), (size_t)MAX_TEST_STR_LEN);
            EXPECT_EQ(strncmp(buf, out, cmpLen), 0)
                << "fixed_point128<40>: input=" << buf << " output=" << out;
        }

        // fixed_point128<1>: I=1, F=127, max_frac_digits=38; integer part is 0 or 1
        genStr(buf, sizeof(buf), 1u, 10);
        {
            fixed_point128<1> f1 = buf;
            const char* out = static_cast<char*>(f1);
            size_t cmpLen = std::min(strlen(out), (size_t)MAX_TEST_STR_LEN);
            EXPECT_EQ(strncmp(buf, out, cmpLen), 0)
                << "fixed_point128<1>: input=" << buf << " output=" << out;
        }

        // fixed_point128<64>: I=64, F=64, max_frac_digits=19; integer part bounded to 999999
        genStr(buf, sizeof(buf), 999999u, 10);
        {
            fixed_point128<64> f64 = buf;
            const char* out = static_cast<char*>(f64);
            size_t cmpLen = std::min(strlen(out), (size_t)MAX_TEST_STR_LEN);
            EXPECT_EQ(strncmp(buf, out, cmpLen), 0)
                << "fixed_point128<64>: input=" << buf << " output=" << out;
        }
    }
}
*/
TEST(fixed_point128, TemplateOperatorNotEqual)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = get_double_random();
        fixed_point128<40> f1 = value1;
        bool fp128_res = f1 != value1;

        if (check_overflow(value1, f1))
            continue;

        EXPECT_TRUE(fp128_res == false) << "operator!=<double>: "
                                        << "value1=" << value1;

        f1 = static_cast<float>(value1);
        fp128_res = f1 != static_cast<float>(value1);

        EXPECT_TRUE(fp128_res == false) << "operator!=<float>: "
                                        << "value1=" << value1;

        f1 = static_cast<uint64_t>(value1);
        fp128_res = f1 != static_cast<uint64_t>(value1);

        EXPECT_TRUE(fp128_res == false) << "operator!=<uint64_t>: "
                                        << "value1=" << value1;

        f1 = static_cast<int64_t>(value1);
        fp128_res = f1 != static_cast<int64_t>(value1);

        EXPECT_TRUE(fp128_res == false) << "operator!=<int64_t>: "
                                        << "value1=" << value1;

        f1 = static_cast<uint32_t>(value1);
        fp128_res = f1 != static_cast<uint32_t>(value1);

        EXPECT_TRUE(fp128_res == false) << "operator!=<uint32_t>: "
                                        << "value1=" << value1;

        f1 = static_cast<int32_t>(value1);
        fp128_res = f1 != static_cast<int32_t>(value1);

        EXPECT_TRUE(fp128_res == false) << "operator!=<int32_t>: "
                                        << "value1=" << value1;
    }
}
TEST(fixed_point128, ShiftRight)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = get_uint32_random();
        uint32_t shift = get_uint32_random() % 40u;
        fixed_point128<40> f1 = value1;

        double expo = pow(2, -static_cast<int32_t>(shift));
        double res = value1 * expo;  // double doesn't lose any bits with this operation!
        fixed_point128 f3 = f1 >> shift;

        if (check_overflow(value1, f1) || check_overflow(res, f3)) {
            continue;
        }
        EXPECT_DOUBLE_EQ(static_cast<double>(f3), res) << "double value1=" << value1 << "; uint32_t shift=" << shift << ";";
    }
}
TEST(fixed_point128, ShiftLeft)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = 1.0 / get_int32_random();
        uint32_t shift = get_uint32_random() % 40u;
        fixed_point128<40> f1 = value1;
        double expo = pow(2, static_cast<int32_t>(shift));
        double res = value1 * expo;  // double doesn't lose any bits with this operation!
        fixed_point128 f3 = f1 << shift;
        if (check_overflow(value1, f1) || check_overflow(res, f3))
            continue;
        EXPECT_DOUBLE_EQ(static_cast<double>(f3), res) << "double value1=" << value1 << "; uint32_t shift=" << shift << ";";
    }
}
TEST(fixed_point128, floor)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random(-15, 15);  // can't exceed this range to avoid overflow
        double res = floor(value);
        fixed_point128<16> f1 = value;
        fixed_point128<16> fp128_res = floor(f1);
        EXPECT_DOUBLE_EQ(fp128_res, res) << "floor: "
                                         << "value=" << value;
    }
}
TEST(fixed_point128, ceil)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random(-15, 15);  // can't exceed this range to avoid overflow
        double res = ceil(value);
        fixed_point128<16> f1 = value;
        fixed_point128<16> fp128_res = ceil(f1);
        EXPECT_DOUBLE_EQ(fp128_res, res) << "ceil: "
                                         << "value=" << value;
    }
}
TEST(fixed_point128, trunc)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random(-15, 15);  // can't exceed this range to avoid overflow
        double res = trunc(value);
        fixed_point128<16> f1 = value;
        fixed_point128<16> fp128_res = trunc(f1);
        EXPECT_DOUBLE_EQ(fp128_res, res) << "trunc: "
                                         << "value=" << value;
    }
}
TEST(fixed_point128, round)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random(-15, 15);  // can't exceed this range to avoid overflow
        double res = round(value);
        fixed_point128<16> f1 = value;
        fixed_point128<16> fp128_res = round(f1);
        EXPECT_DOUBLE_EQ(fp128_res, res) << "round: "
                                         << "value=" << value;
    }
}
TEST(fixed_point128, copysign)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = get_double_random(-15, 15);
        double value2 = get_double_random(-15, 15);
        double res = copysign(value1, value2);
        fixed_point128<16> f1 = value1;
        fixed_point128<16> f2 = value2;
        fixed_point128<16> fp128_res = copysign(f1, f2);
        EXPECT_DOUBLE_EQ(fp128_res, res) << "copysign: "
                                         << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(fixed_point128, reciprocal)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random(-15, 15);  // can't exceed this range to avoid overflow
        double res = 1.0 / value;
        fixed_point128<16> f1 = value;
        fixed_point128<16> fp128_res = reciprocal(f1);
        EXPECT_DOUBLE_EQ(fp128_res, res) << "reciprocal: "
                                         << "value=" << value;
    }
}
TEST(fixed_point128, sqrt)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = fabs(get_double_random(-30, 15));  // can't exceed this range to avoid overflow
        double res = ::sqrt(value);
        fixed_point128<16> f1 = value;
        fixed_point128<16> fp128_res = sqrt(f1);
        EXPECT_DOUBLE_EQ(fp128_res, res) << "sqrt: "
                                         << "value=" << value;
    }
}
TEST(fixed_point128, hypot)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = get_double_random(-15, 3);
        double value2 = get_double_random(-15, 3);
        double res = hypot(value1, value2);
        fixed_point128<16> f1 = value1;
        fixed_point128<16> f2 = value2;
        fixed_point128<16> fp128_res = hypot(f1, f2);
        EXPECT_DOUBLE_EQ(fp128_res, res) << "hypot: "
                                         << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(fixed_point128, ilogb)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random(-60, 15);  // lower exponent results in lost bits
        int32_t res = ::ilogb(value);
        fixed_point128<16> f1 = value;
        int32_t fp128_res = ilogb(f1);
        EXPECT_EQ(fp128_res, res) << "ilogb: "
                                  << "value=" << value;
    }
}
TEST(fixed_point128, log)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = fabs(get_double_random(-40, 15));  // lower exponent results in lost bits
        double res = log(value);
        fixed_point128<16> f1 = value;
        fixed_point128<16> fp128_res = log(f1);
        EXPECT_DOUBLE_EQ(fp128_res, res) << "log: "
                                         << "value=" << value;
    }
}
TEST(fixed_point128, log2)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = fabs(get_double_random(-40, 15));  // lower exponent results in lost bits
        double res = log2(value);
        fixed_point128<16> f1 = value;
        fixed_point128<16> fp128_res = log2(f1);
        EXPECT_DOUBLE_EQ(fp128_res, res) << "log2: "
                                         << "value=" << value;
    }
}
TEST(fixed_point128, log10)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = fabs(get_double_random(-40, 15));  // lower exponent results in lost bits
        double res = log10(value);
        fixed_point128<16> f1 = value;
        fixed_point128<16> fp128_res = log10(f1);
        EXPECT_DOUBLE_EQ(fp128_res, res) << "log10: "
                                         << "value=" << value;
    }
}
TEST(fixed_point128, log1p)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = fabs(get_double_random(-40, 15));  // lower exponent results in lost bits
        double res = log1p(value);
        fixed_point128<16> f1 = value;
        fixed_point128<16> fp128_res = log1p(f1);
        EXPECT_DOUBLE_EQ(fp128_res, res) << "log1p: "
                                         << "value=" << value;
    }
}
TEST(fixed_point128, logb)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = fabs(get_double_random(-40, 15));  // lower exponent results in lost bits
        double res = logb(value);
        fixed_point128<16> f1 = value;
        fixed_point128<16> fp128_res = logb(f1);
        EXPECT_DOUBLE_EQ(fp128_res, res) << "logb: "
                                         << "value=" << value;
    }
}
TEST(fixed_point128, sin)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random(-60, 2);
        double res = ::sin(value);
        fixed_point128<16> f1 = value;
        fixed_point128<16> fp128_res = sin(f1);
        EXPECT_DOUBLE_EQ(fp128_res, res) << "sin: "
                                         << "value=" << value;
    }
}
TEST(fixed_point128, cos)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random(-60, 2);
        double res = ::cos(value);
        fixed_point128<16> f1 = value;
        fixed_point128<16> fp128_res = cos(f1);
        EXPECT_DOUBLE_EQ(fp128_res, res) << "cos: "
                                         << "value=" << value;
    }
}
TEST(fixed_point128, tan)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random(-60, 2);
        double res = ::tan(value);
        fixed_point128<16> f1 = value;
        fixed_point128<16> fp128_res = tan(f1);
        EXPECT_DOUBLE_EQ(fp128_res, res) << "tan: "
                                         << "value=" << value;
    }
}
TEST(fixed_point128, asin)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random(-60, -1);
        double res = ::asin(value);
        fixed_point128<16> f1 = value;
        fixed_point128<16> fp128_res = asin(f1);
        EXPECT_DOUBLE_EQ(fp128_res, res) << "asin: "
                                         << "value=" << value;
    }
}
TEST(fixed_point128, acos)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random(-60, -1);
        double res = acos(value);
        fixed_point128<16> f1 = value;
        fixed_point128<16> fp128_res = acos(f1);
        EXPECT_DOUBLE_EQ(fp128_res, res) << "acos: "
                                         << "value=" << value;
    }
}
TEST(fixed_point128, atan)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random(-60, 14);  // lower exponent results in lost bits
        double res = atan(value);
        fixed_point128<16> f1 = value;
        fixed_point128<16> fp128_res = atan(f1);
        EXPECT_DOUBLE_EQ(fp128_res, res) << "atan: "
                                         << "value=" << value;
    }
}
TEST(fixed_point128, atan2)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = get_double_random(-10, 14);  // lower exponent results in lost bits
        double value2 = get_double_random(-10, 14);  // lower exponent results in lost bits
        double res = atan2(value1, value2);
        fixed_point128<16> f1 = value1;
        fixed_point128<16> f2 = value2;
        fixed_point128<16> fp128_res = atan2(f1, f2);
        EXPECT_DOUBLE_EQ(fp128_res, res) << "atan2: "
                                         << " value1=" << value1 << ", value2=" << value2;
    }
}
TEST(fixed_point128, exp)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random(-3, 3);  // lower exponent results in lost bits
        double res = exp(value);
        fixed_point128<16> f1 = value;
        fixed_point128<16> fp128_res = exp(f1);
        EXPECT_DOUBLE_EQ(fp128_res, res) << "exp: "
                                         << "value=" << value;
    }
}
TEST(fixed_point128, exp2)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random(-4, 4);  // lower exponent results in lost bits
        double res = exp2(value);
        fixed_point128<16> f1 = value;
        fixed_point128<16> fp128_res = exp2(f1);
        EXPECT_DOUBLE_EQ(fp128_res, res) << "exp2: "
                                         << "value=" << value;
    }
}
TEST(fixed_point128, expm1)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random(-3, 3);  // lower exponent results in lost bits
        double res = expm1(value);
        fixed_point128<16> f1 = value;
        fixed_point128<16> fp128_res = expm1(f1);
        EXPECT_DOUBLE_EQ(fp128_res, res) << "expm1: "
                                         << "value=" << value;
    }
}
TEST(fixed_point128, pow)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = fabs(get_double_random(-4, 4));  // lower exponent results in lost bits
        double value2 = get_double_random(-4, 2);        // lower exponent results in lost bits
        double res = pow(value1, value2);
        fixed_point128<16> f1 = value1;
        fixed_point128<16> f2 = value2;
        fixed_point128<16> fp128_res = pow(f1, f2);
        EXPECT_DOUBLE_EQ(fp128_res, res) << "pow: "
                                         << " value1=" << value1 << ", value2=" << value2;
    }
}
TEST(fixed_point128, sinh)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random(-60, 2);
        double res = ::sinh(value);
        fixed_point128<16> f1 = value;
        fixed_point128<16> fp128_res = sinh(f1);
        EXPECT_DOUBLE_EQ(fp128_res, res) << "sinh: "
                                         << "value=" << value;
    }
}
TEST(fixed_point128, asinh)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random(-60, 2);
        double res = ::asinh(value);
        fixed_point128<16> f1 = value;
        fixed_point128<16> fp128_res = asinh(f1);
        EXPECT_DOUBLE_EQ(fp128_res, res) << "asinh: "
                                         << "value=" << value;
    }
}
TEST(fixed_point128, cosh)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random(-60, 2);
        double res = ::cosh(value);
        fixed_point128<16> f1 = value;
        fixed_point128<16> fp128_res = cosh(f1);
        EXPECT_DOUBLE_EQ(fp128_res, res) << "cosh: "
                                         << "value=" << value;
    }
}
TEST(fixed_point128, acosh)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = 1.0 + fabs(get_double_random(-60, 2));  // must be >= 1
        double res = ::acosh(value);
        fixed_point128<16> f1 = value;
        fixed_point128<16> fp128_res = acosh(f1);
        EXPECT_DOUBLE_EQ(fp128_res, res) << "acosh: "
                                         << "value=" << value;
    }
}
TEST(fixed_point128, tanh)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random(-60, 2);
        double res = ::tanh(value);
        fixed_point128<16> f1 = value;
        fixed_point128<16> fp128_res = tanh(f1);
        EXPECT_DOUBLE_EQ(fp128_res, res) << "tanh: "
                                         << "value=" << value;
    }
}

TEST(fixed_point128, atanh)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random(-60, -1);  // must be: -1 < value < 1
        double res = ::atanh(value);
        fixed_point128<16> f1 = value;
        fixed_point128<16> fp128_res = atanh(f1);
        EXPECT_DOUBLE_EQ(fp128_res, res) << "atanh: "
                                         << "value=" << value;
    }
}

/**********************************************************************
 * fixed_point128 regression tests
 *
 * Each test below pins down a defect that was fixed.
 ***********************************************************************/

// The sign lives in its own field, so a result that lands on zero while carrying a sign is a
// distinct bit pattern from a plain zero. Every comparison operator tests the sign field first,
// so such a value compared unequal to zero and smaller than it. trunc, ceil, round, modf and
// copysign all produced one.
TEST(fixed_point128, NoNegativeZero)
{
    typedef fixed_point128<32> fp;
    const fp zero(0);

    EXPECT_TRUE(trunc(fp(-0.5)) == zero);
    EXPECT_FALSE(trunc(fp(-0.5)) < zero);
    EXPECT_FALSE(trunc(fp(-0.5)).is_negative());

    EXPECT_TRUE(ceil(fp(-0.5)) == zero);
    EXPECT_FALSE(ceil(fp(-0.5)).is_negative());

    EXPECT_TRUE(round(fp(-0.4)) == zero);
    EXPECT_FALSE(round(fp(-0.4)).is_negative());

    EXPECT_TRUE(copysign(zero, fp(-1)) == zero);
    EXPECT_FALSE(copysign(zero, fp(-1)).is_negative());

    fp int_part;
    const fp frac = modf(fp(-0.25), &int_part);
    EXPECT_TRUE(int_part == zero);
    EXPECT_FALSE(int_part.is_negative());
    EXPECT_TRUE(frac.is_negative());  // the fraction keeps the sign, it is non zero

    // floor was already correct, keep it covered
    EXPECT_TRUE(floor(fp(-0.5)) == fp(-1));
    EXPECT_TRUE(trunc(fp(-1.5)) == fp(-1));
}
// The conversions truncated the bits that did not fit the mantissa, and a value whose surviving
// fraction bits happened to be all ones was rounded up even when it was exactly representable,
// so (2^24-1)/2^23 came out as 2. Converting the 128 bit magnitude and scaling by a power of two
// with ldexp is exact, which makes it a correctly rounded reference.
TEST(fixed_point128, ConversionToFloatingPointIsCorrectlyRounded)
{
    typedef fixed_point128<32> fp;
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const uint64_t l = get_uint64_random();
        const uint64_t h = get_uint64_random() >> (get_uint32_random() % 40);
        const uint32_t s = get_uint32_random() & 1;
        const fp v(l, h, s);
        if (!v)
            continue;
        const double mag = static_cast<double>(uint128_t(l, h));
        const double ref = s ? -ldexp(mag, -fp::F) : ldexp(mag, -fp::F);
        EXPECT_EQ(static_cast<double>(v), ref) << "low=" << l << ", high=" << h;
        EXPECT_EQ(static_cast<float>(v), static_cast<float>(ref)) << "low=" << l << ", high=" << h;
    }
}
TEST(fixed_point128, ConversionToFloatingPointEdgeCases)
{
    typedef fixed_point128<32> fp;
    // (2^24-1)/2^23 is exactly representable as a float and must not be rounded up to 2
    const fp exact_flt = fp(0xFFFFFFull) / fp(0x800000ull);
    EXPECT_LT(static_cast<float>(exact_flt), 2.0f);
    EXPECT_EQ(static_cast<float>(exact_flt), 16777215.0f / 8388608.0f);

    // the same trap one mantissa wider, for double
    typedef fixed_point128<64> fp64;
    const fp64 exact_dbl = fp64(0x1FFFFFFFFFFFFFull) / fp64(0x10000000000000ull);
    EXPECT_LT(static_cast<double>(exact_dbl), 2.0);

    // powers of two convert exactly in both directions
    for (int e = -60; e <= 20; ++e) {
        fp v(1);
        if (e > 0)
            v <<= e;
        else if (e < 0)
            v >>= -e;
        EXPECT_EQ(static_cast<double>(v), ::pow(2.0, e)) << "e=" << e;
        EXPECT_EQ(static_cast<float>(v), static_cast<float>(::pow(2.0, e))) << "e=" << e;
    }
    // the smallest value of the widest fraction is a denormal float, it must not become zero
    const fixed_point128<1> tiny(1, 0, 0);
    EXPECT_GT(static_cast<double>(tiny), 0.0);
    EXPECT_GT(static_cast<float>(tiny), 0.0f);
}
// div_64bit shortcuts a numerator smaller than the divisor and returns without writing the
// quotient, so dividing by a scalar left the value completely unchanged instead of producing zero.
TEST(fixed_point128, DivideSmallValueByScalar)
{
    typedef fixed_point128<32> fp;
    // the raw 128 bit form of 2^-90 is 64, which is smaller than the divisor below
    const fp tiny = fp(1) >> 90;
    EXPECT_TRUE(tiny);  // the value itself is representable

    fp v = tiny;
    v /= 1000ull;
    EXPECT_TRUE(v == fp(0)) << "expected zero, got " << (std::string)v;

    // the fixed_point128 operand overload has always produced zero here, the two must agree
    fp w = tiny;
    w /= fp(1000);
    EXPECT_TRUE(v == w);

    // a normal division still works
    fp u = fp(1) >> 40;
    u /= 1000ull;
    EXPECT_NEAR(static_cast<double>(u), ::pow(2.0, -40) / 1000.0, ::pow(2.0, -40) / 1e9);
}
// These are noexcept and used to reach a throwing path, which terminates the process instead.
TEST(fixed_point128, NoexceptFunctionsDoNotTerminate)
{
    typedef fixed_point128<16> fp;

    // pow forwards to log, which rejects zero
    EXPECT_TRUE(pow(fp(0), fp(2)) == fp(0));
    EXPECT_TRUE(pow(fp(0), fp(0)) == fp(1));  // the CRT defines pow(0, 0) as 1
    EXPECT_TRUE(pow(fp(-1), fp(2)) == fp(0));

    // acos at the ends of its domain starts where the derivative vanishes and sin(res) is zero
    EXPECT_TRUE(acos(fp(1)) == fp(0));
    EXPECT_NEAR(static_cast<double>(acos(fp(-1))), ::acos(-1.0), 1e-12);
    EXPECT_NEAR(static_cast<double>(asin(fp(1))), ::asin(1.0), 1e-12);

    // tan's range reduction lands exactly on a zero cosine at the poles
    EXPECT_TRUE(tan(fp::half_pi()) == fp(0));

    // reciprocal is documented to return zero rather than saturate
    EXPECT_TRUE(reciprocal(fp(0)) == fp(0));
}
// log1p forwards to log, which throws on a non positive argument. It used to be noexcept, which
// turned that into a call to std::terminate.
TEST(fixed_point128, Log1pPropagatesDomainError)
{
    typedef fixed_point128<16> fp;
    EXPECT_THROW((void)log1p(fp(-1)), std::domain_error);
    EXPECT_THROW((void)log1p(fp(-2)), std::domain_error);
    EXPECT_NEAR(static_cast<double>(log1p(fp(1))), ::log(2.0), 1e-12);
    // the log family itself keeps throwing
    EXPECT_THROW((void)log(fp(0)), std::domain_error);
    EXPECT_THROW((void)log2(fp(-1)), std::domain_error);
    EXPECT_THROW((void)log10(fp(0)), std::domain_error);
}
// A leading sign or a high byte handed straight to the ctype functions is outside the domain they
// accept, and only a literal space was skipped as leading white space.
TEST(fixed_point128, ConstructorFromStringEdgeCases)
{
    typedef fixed_point128<32> fp;
    EXPECT_TRUE(fp("\t 2.5") == fp(2.5));
    EXPECT_TRUE(fp("  -2.5") == fp(-2.5));
    EXPECT_TRUE(fp("+2.5") == fp(2.5));
    const char high_byte[] = {'1', '2', static_cast<char>(0xB5), 0};
    EXPECT_TRUE(fp(high_byte) == fp(12));
}

// A scalar on the left used to be ambiguous: converting it to fixed_point128 and converting the
// fixed_point128 to a builtin type are both a single user defined conversion. The scalar is the one
// that gets widened, so the fraction of the object survives instead of being truncated away.
TEST(fixed_point128, ScalarOnTheLeftHandSide)
{
    typedef fixed_point128<20> fp;
    const fp quarter = fp(1) / 4;  // 0.25

    EXPECT_DOUBLE_EQ(static_cast<double>(1 + quarter), 1.25);
    EXPECT_DOUBLE_EQ(static_cast<double>(1 - quarter), 0.75);
    EXPECT_DOUBLE_EQ(static_cast<double>(3 * quarter), 0.75);
    EXPECT_DOUBLE_EQ(static_cast<double>(1 / quarter), 4.0);
    EXPECT_DOUBLE_EQ(static_cast<double>(5 % fp(3)), 2.0);
    EXPECT_DOUBLE_EQ(static_cast<double>(2.5 + quarter), 2.75);
    EXPECT_DOUBLE_EQ(static_cast<double>(2.5 - quarter), 2.25);

    // subtracting past zero must keep the sign
    EXPECT_DOUBLE_EQ(static_cast<double>(1 - fp(3)), -2.0);
    EXPECT_DOUBLE_EQ(static_cast<double>(-1 + quarter), -0.75);

    // the bitwise forms agree with the fixed_point128 on the left ones
    EXPECT_TRUE((7 & fp(3)) == (fp(7) & 3));
    EXPECT_TRUE((7 | fp(8)) == (fp(7) | 8));
    EXPECT_TRUE((7 ^ fp(3)) == (fp(7) ^ 3));

    // the fixed_point128 on the left forms must still resolve
    EXPECT_TRUE((quarter + 1) == fp(1.25));
    EXPECT_TRUE((quarter + quarter) == fp(0.5));
    EXPECT_TRUE((quarter * quarter) == fp(0.0625));

    // The left operand is widened rather than the object being narrowed, so the result is 128 bit.
    // Narrowing instead would make (1 - quarter) an int and lose the fraction entirely.
    static_assert(std::is_same_v<decltype(1 + quarter), fp>);
    static_assert(std::is_same_v<decltype(1 ^ quarter), fp>);
}
// Widening the left operand keeps the shift at 128 bits. The count is the right operand converted
// to int32_t, which for fixed_point128 is a shift and so truncates toward zero.
TEST(fixed_point128, ScalarOnTheLeftHandSideShifts)
{
    typedef fixed_point128<20> fp;

    EXPECT_DOUBLE_EQ(static_cast<double>(1 << fp(3)), 8.0);
    EXPECT_DOUBLE_EQ(static_cast<double>(8 >> fp(3)), 1.0);
    EXPECT_DOUBLE_EQ(static_cast<double>(3 << fp(2)), 12.0);
    EXPECT_DOUBLE_EQ(static_cast<double>(1 >> fp(2)), 0.25);

    // A fractional shift count goes through operator int32_t, which truncates toward zero here, so
    // 3.75 shifts by 3. The same overload picked with the fixed_point128 on the left has to agree,
    // both reaching the count the same way. float128 rounds instead, see its matching test.
    EXPECT_TRUE((1 << fp(3.75)) == (1 << fp(3)));
    EXPECT_TRUE((1 << fp(3.75)) == (fp(1) << fp(3.75)));

    static_assert(std::is_same_v<decltype(1 << fp(3)), fp>);

    // the fixed_point128 on the left forms must still resolve
    EXPECT_TRUE((fp(4) << 1) == fp(8));
    EXPECT_TRUE((fp(4) >> 1) == fp(2));
}

/**********************************************************************
 * Compile time (constexpr) evaluation
 *
 * Each test below pairs static_asserts, which fail the build the moment one of these operations
 * stops being usable in a constant expression, with a runtime check of the same expression.
 *
 * The runtime half earns its place. A constant evaluated call does not run the same code a runtime
 * call does: the bit counting and extended arithmetic intrinsics (mulx_u64, addcarryx_u64,
 * lzcnt64) are not constant expressions, so fp128_shared.h substitutes a portable
 * implementation of each one while the compiler is evaluating. Only comparing the two results
 * shows that the substitutes agree with the hardware. opaque(), declared in gtest_shared.h, is
 * what keeps the runtime half from being constant folded back into the compile time one.
 ***********************************************************************/
TEST(fixed_point128, ConstexprConstructionAndConversion)
{
    typedef fixed_point128<32> fp;

    constexpr fp zero;
    constexpr fp from_u64(static_cast<uint64_t>(7));
    constexpr fp from_i64(static_cast<int64_t>(-7));
    constexpr fp from_u32(static_cast<uint32_t>(7));
    constexpr fp from_i32(static_cast<int32_t>(-7));
    constexpr fp from_bits(1, 0, 0);  // the low/high/sign constructor
    constexpr fp copied(from_i32);
    constexpr fp assigned = [] { fp a; a = fp(5); return a; }();

    static_assert(zero.is_zero());
    static_assert(static_cast<uint64_t>(from_u64) == 7ull);
    static_assert(static_cast<int64_t>(from_i64) == -7ll);
    static_assert(static_cast<uint32_t>(from_u32) == 7u);
    static_assert(static_cast<int32_t>(from_i32) == -7);
    static_assert(from_bits == fp::epsilon());
    static_assert(copied == from_i32);
    static_assert(static_cast<int32_t>(assigned) == 5);

    // the cross template conversions, in both directions
    constexpr fixed_point128<16> fewer_int_bits(from_u32);  // I < I2, shifts left
    constexpr fixed_point128<48> more_int_bits(from_u32);   // I > I2, shifts right
    static_assert(static_cast<int32_t>(fewer_int_bits) == 7);
    static_assert(static_cast<int32_t>(more_int_bits) == 7);

    // the exact constants, which are compile time values rather than lazily initialized statics
    static_assert(static_cast<int32_t>(fp::one()) == 1);
    static_assert(fp::half() + fp::half() == fp::one());
    static_assert(fp::epsilon() > fp(0));
    static_assert(fixed_point128<64>::half() < fixed_point128<64>::one());  // 0.5 lands in the low QWORD
    static_assert(fixed_point128<1>::half() < fixed_point128<1>::one());

    EXPECT_TRUE(fp(opaque(-7)) == from_i32);
    EXPECT_TRUE(fixed_point128<16>(fp(opaque(7))) == fewer_int_bits);
    EXPECT_TRUE(fixed_point128<48>(fp(opaque(7))) == more_int_bits);
}
TEST(fixed_point128, ConstexprShiftsBitwiseAndUnary)
{
    typedef fixed_point128<32> fp;

    static_assert(static_cast<int32_t>(fp(8) >> 3) == 1);
    static_assert(static_cast<int32_t>(fp(1) << 3) == 8);
    static_assert(static_cast<int32_t>(fp(1) << 70) == 0);  // overflow is silent, the value wraps
    constexpr fp shifted = [] { fp a(16); a >>= 2; a <<= 1; return a; }();
    static_assert(static_cast<int32_t>(shifted) == 8);

    static_assert(static_cast<int32_t>(fp(12) & fp(10)) == 8);
    static_assert(static_cast<int32_t>(fp(12) | fp(3)) == 15);
    static_assert(static_cast<int32_t>(fp(12) ^ fp(10)) == 6);
    static_assert(!(~fp(0)).is_zero());

    static_assert((-fp(5)).is_negative());
    static_assert((+fp(5)).is_positive());
    static_assert((-fp(0)).is_positive());  // negating zero must not produce a negative zero
    static_assert(!fp(0));
    static_assert(static_cast<bool>(fp(1)));

    EXPECT_TRUE((fp(opaque(16)) >> 2 << 1) == shifted);
    EXPECT_TRUE((fp(opaque(12)) & fp(opaque(10))) == fp(8));
    EXPECT_TRUE(-fp(opaque(0)) == fp(0));
}
TEST(fixed_point128, ConstexprComparisonsAndQueries)
{
    typedef fixed_point128<32> fp;

    static_assert(fp(1) < fp(2));
    static_assert(fp(2) > fp(1));
    static_assert(fp(-2) < fp(-1));
    static_assert(fp(1) <= 1 && fp(1) >= 1);
    static_assert(fp(1) == 1 && fp(1) != 2);
    static_assert(1 == fp(1) && 2 > fp(1));  // the overloads taking the fixed_point128 on the right

    static_assert(fp(3).is_int());
    static_assert(!(fp(1) >> 1).is_int());
    static_assert(fixed_point128<64>(3).is_int());  // the I == 64 branch, where there are no fraction bits in high
    static_assert(fp(0).is_zero());
    static_assert(fp(-1).is_negative() && fp(1).is_positive());
    static_assert(fp(1).get_bit(fp::F) == 1);  // 1.0 has its single set bit at the radix point
    static_assert(fp(4).get_exponent() == 2);
    static_assert((fp::one() >> 1).get_exponent() == -1);

    EXPECT_TRUE(fp(opaque(4)).get_exponent() == 2);
    EXPECT_TRUE(fp(opaque(3)).is_int());
    EXPECT_TRUE(fp(opaque(1)) < fp(opaque(2)));
}
TEST(fixed_point128, ConstexprArithmetic)
{
    typedef fixed_point128<32> fp;

    static_assert(static_cast<int32_t>(fp(3) + fp(4)) == 7);
    static_assert(static_cast<int32_t>(fp(3) - fp(4)) == -1);
    static_assert(static_cast<int32_t>(fp(-3) + fp(-4)) == -7);
    static_assert((fp(3) - fp(3)).is_zero());
    static_assert((fp(3) - fp(3)).is_positive());         // cancelling to zero must drop the sign
    static_assert(static_cast<int32_t>(fp(3) + 4) == 7);  // the generic right hand side overload

    static_assert(static_cast<int32_t>(fp(6) * fp(7)) == 42);
    static_assert(static_cast<int32_t>(fp(-6) * fp(7)) == -42);
    static_assert(fp::half() * fp(8) == fp(4));
    static_assert(static_cast<int32_t>(fp(6) * static_cast<uint64_t>(7)) == 42);  // the uint64_t specialization
    static_assert(static_cast<int32_t>(fp(-6) * -7) == 42);
    static_assert(static_cast<int32_t>(fixed_point128<64>(6) * fixed_point128<64>(7)) == 42);  // the F == 64 branch
    static_assert(fixed_point128<1>::half() * fixed_point128<1>::half() == fixed_point128<1>::half() >> 1);

    // sqr is documented to be bit identical to x * x
    static_assert(sqr(fp::half()) == fp::half() * fp::half());
    static_assert(static_cast<int32_t>(sqr(fp(-9))) == 81);
    static_assert(sqr(fp(-9)).is_positive());

    constexpr fp stepped = [] { fp a(5); ++a; a++; --a; return a; }();
    static_assert(static_cast<int32_t>(stepped) == 6);
    constexpr fp accumulated = [] { fp a; for (int32_t i = 1; i <= 10; ++i) a += fp(i); return a; }();
    static_assert(static_cast<int32_t>(accumulated) == 55);

    EXPECT_TRUE(fp(opaque(3)) + fp(opaque(4)) == fp(7));
    EXPECT_TRUE(fp(opaque(3)) - fp(opaque(4)) == fp(-1));
    EXPECT_TRUE(fp(opaque(6)) * fp(opaque(7)) == fp(42));
    EXPECT_TRUE(sqr(fp(opaque(-9))) == fp(81));
    EXPECT_TRUE(fixed_point128<64>(opaque(6)) * fixed_point128<64>(opaque(7)) == fixed_point128<64>(42));
}
TEST(fixed_point128, ConstexprMathFunctions)
{
    typedef fixed_point128<32> fp;
    constexpr fp two_and_a_half = fp(5) >> 1;

    static_assert(fabs(fp(-5)) == fp(5));
    static_assert(floor(two_and_a_half) == fp(2));
    static_assert(floor(-two_and_a_half) == fp(-3));
    static_assert(ceil(two_and_a_half) == fp(3));
    static_assert(ceil(-two_and_a_half) == fp(-2));
    static_assert(trunc(-two_and_a_half) == fp(-2));
    static_assert(round(fp::half()) == fp(1));  // the halfway value rounds away from zero
    static_assert(round(-fp::half()) == fp(-1));
    static_assert(round(fp::half() >> 1).is_positive());  // -0.25 rounds to +0, never to -0
    static_assert(copysign(fp(5), fp(-1)).is_negative());
    static_assert(fmin(fp(1), fp(2)) == fp(1));
    static_assert(fmax(fp(1), fp(2)) == fp(2));
    static_assert(fdim(fp(5), fp(3)) == fp(2));
    static_assert(fdim(fp(3), fp(5)) == fp(0));
    static_assert(ilogb(fp(8)) == 3);
    static_assert(lzcnt128(fp::one()) == 31);  // 1.0 sets bit F, leaving the I-1 integer bits above it clear

    // modf splits into an integer and a fraction part, both carrying the sign of the input
    constexpr fp modf_int = [] { fp ip; (void)modf(-(fp(5) >> 1), &ip); return ip; }();
    constexpr fp modf_frac = [] { fp ip; return modf(-(fp(5) >> 1), &ip); }();
    static_assert(modf_int == fp(-2));
    static_assert(modf_frac == -fp::half());

    static_assert(log2(fp(8)) == fp(3));  // an exact power of two takes the shortcut
    static_assert(logb(fp(9)) == fp(3));
    constexpr fp log2_of_10 = log2(fp(10));  // the full fraction loop, one squaring per fraction bit
    static_assert(log2_of_10 > fp(3) && log2_of_10 < fp(4));

    EXPECT_TRUE(floor(-(fp(opaque(5)) >> 1)) == fp(-3));
    EXPECT_TRUE(round(-fp::half() * fp(opaque(1))) == fp(-1));
    EXPECT_TRUE(ilogb(fp(opaque(8))) == 3);
    // the loop above runs ~100 squarings, every one of them through mulx_u64
    EXPECT_TRUE(log2(fp(opaque(10))) == log2_of_10);
}

// The double and float conversions used to be the one thing that could not happen at compile time:
// they read the inactive member of a union, which constant evaluation rejects. Now that Double and
// Float convert with std::bit_cast, a floating point literal crosses the boundary either way.
TEST(fixed_point128, ConstexprFloatingPointConversion)
{
    typedef fixed_point128<32> fp;

    constexpr fp from_double = 3.25;
    constexpr fp negative = -2.5;
    static_assert(static_cast<double>(from_double) == 3.25);
    static_assert(static_cast<float>(from_double) == 3.25f);
    static_assert(static_cast<double>(negative) == -2.5);
    static_assert(static_cast<long double>(from_double) == 3.25L);
    static_assert(fp(0.5) == fp::half());
    static_assert(fp(1.0) == fp::one());
    static_assert(fp(0.0).is_zero());
    static_assert(fp(-0.0).is_zero() && fp(-0.0).is_positive());  // no negative zero
    static_assert(from_double + negative == fp(0.75));
    static_assert(static_cast<int32_t>(fp(3.99)) == 3);  // truncates towards zero

    // a value that needs the full fraction, not just a few bits
    constexpr fp tenth = 0.1;
    static_assert(tenth > fp(0.09) && tenth < fp(0.11));

    EXPECT_TRUE(fp(static_cast<double>(opaque(13)) / 4.0) == fp(3.25));
    EXPECT_DOUBLE_EQ(static_cast<double>(fp(opaque(13)) / 4), 3.25);
    EXPECT_TRUE(fp(0.1) == tenth);  // the runtime conversion agrees with the constant evaluated one
}
