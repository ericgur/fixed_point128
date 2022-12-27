// remove warnings from gtest itself
#pragma warning(push)
#pragma warning(disable: 26439) 
#pragma warning(disable: 26495) 
#include <gtest/gtest.h>
#pragma warning(pop)
#include <ostream>
#include <ctime>
#include "gtest_shared.h"

/**********************************************************************
* uint128_t tests
***********************************************************************/
// Construct fixed_point128 and convert back to/from various elements.
TEST(uint128_t, DefaultConstructor) {
    uint128_t i;
    EXPECT_EQ(static_cast<uint64_t>(i), 0ull);
}
TEST(uint128_t, ConstructorFromDouble) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = fabs(floor(get_double_random()));
        uint128_t f = value;
        if (check_overflow_uint128(value)) {
            continue;
        }
        double f_value = static_cast<double>(f);
        EXPECT_DOUBLE_EQ(f_value, value) << "value=" << value;
    }
}
TEST(uint128_t, ConstructorFromFloat) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        float value = (float)fabs(floor(get_double_random()));
        uint128_t f = value;
        if (check_overflow_uint128(value)) {
            continue;
        }
        EXPECT_FLOAT_EQ(static_cast<float>(f), value) << "value=" << value;
    }
}
TEST(uint128_t, ConstructorFromInt32) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        int32_t value = get_int32_random();
        uint128_t f = value;
        EXPECT_EQ(static_cast<int32_t>(f), value);
    }
}
TEST(uint128_t, ConstructorFromUnsignedInt32) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        uint32_t value = get_uint32_random();
        uint128_t f = value;
        EXPECT_EQ(static_cast<uint32_t>(f), value);
    }
}
TEST(uint128_t, ConstructorFromInt64) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        int64_t value = get_int64_random();
        uint128_t f = value;
        EXPECT_EQ(static_cast<int64_t>(f), value);
    }
}
TEST(uint128_t, ConstructorFromUnsignedInt64) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        uint64_t value = get_uint64_random();
        uint128_t f = value;
        EXPECT_EQ(static_cast<uint64_t>(f), value);
    }
}
TEST(uint128_t, ConstructorFromString) {
    const char* values[] = { "0xDEADBEAF", "1234", "123456789012345678" };
    for (auto i = 0u; i < array_length(values); ++i) {
        uint128_t f = values[i];
        double d1 = strtod(values[i], nullptr);
        double d2 = strtod(static_cast<char*>(f), nullptr);
        EXPECT_DOUBLE_EQ(d1, d2);
    }
}
TEST(uint128_t, CopyConstructor) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random();
        uint128_t f1 = value;
        uint128_t f2 = f1;
        EXPECT_DOUBLE_EQ(static_cast<double>(f1), static_cast<double>(f2));
    }
}
TEST(uint128_t, MoveConstructor) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random();
        uint128_t f1 = value;
        uint128_t f2(std::move(f1));
        EXPECT_DOUBLE_EQ(static_cast<double>(f1), static_cast<double>(f2));
    }
}
TEST(uint128_t, AssignmentOperator) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random();
        uint128_t f1 = value;
        uint128_t f2;
        f2 = f1;
        EXPECT_DOUBLE_EQ(static_cast<double>(f1), static_cast<double>(f2));
    }
}
TEST(uint128_t, MoveAssignmentOperator) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random();
        uint128_t f1 = value;
        uint128_t f2;
        f2 = std::move(f1);
        EXPECT_DOUBLE_EQ(static_cast<double>(f1), static_cast<double>(f2));
    }
}
TEST(uint128_t, Add) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = fabs(floor(get_double_random()));
        double value2 = floor(value1 * 2.5);
        double res = value1 + value2;
        uint128_t f1 = value1;
        uint128_t f2 = value2;
        uint128_t f3 = f1 + f2;
        if (check_overflow_uint128(value1) || check_overflow_uint128(value2) || check_overflow_uint128(res))
            continue;
        double uint128_res = static_cast<double>(f3);
        // note that uint128_t is more precise than double, this can lead to issues when param1 and param2 are far apart.
        if (is_similar_double(uint128_res, res))
            continue;

        EXPECT_DOUBLE_EQ(uint128_res, res) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(uint128_t, AddDifferentSign) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        uint64_t value1 = get_uint64_random();
    #pragma warning(push)
    #pragma warning(disable: 4146)
        uint64_t value2 = -get_uint64_random();
    #pragma warning(pop)
        uint64_t res = value1 + value2;
        uint128_t f1 = value1;
        uint128_t f2 = value2;
        uint128_t f3 = f1 + f2;
        uint64_t uint128_res = static_cast<uint64_t>(f3);
        EXPECT_EQ(uint128_res, res) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(uint128_t, AddInt64) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        auto value1 = get_int64_random();
        auto value2 = get_int64_random();
        auto res = value1 + value2;
        uint128_t f1 = value1;
        uint128_t f3 = f1 + value2;
        uint64_t uint128_res = static_cast<uint64_t>(f3);
        EXPECT_EQ(uint128_res, res) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(uint128_t, AddUnsignedInt64) {
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

TEST(uint128_t, Subtract) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = fabs(floor(get_double_random()));
        double value2 = floor(value1 / 2.5);
        double res = value1 - value2;
        uint128_t f1 = value1;
        uint128_t f2 = value2;
        uint128_t f3 = f1 - f2;
        if (check_overflow_uint128(value1) || check_overflow_uint128(value2) || check_overflow_uint128(res))
            continue;
        double uint128_res = static_cast<double>(f3);
        // note that fp128 is more precise than double, this can lead to issues when param1 and param2 are far apart.
        if (is_similar_double(uint128_res, res))
            continue;

        EXPECT_DOUBLE_EQ(uint128_res, res) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(uint128_t, SubtractDifferentSign) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        uint64_t value1 = get_uint64_random();
    #pragma warning(push)
    #pragma warning(disable: 4146)
        uint64_t value2 = -get_uint64_random();
    #pragma warning(pop)
        uint64_t res = value1 - value2;
        uint128_t f1 = value1;
        uint128_t f2 = value2;
        uint128_t f3 = f1 - f2;
        uint64_t uint128_res = static_cast<uint64_t>(f3);
        EXPECT_EQ(uint128_res, res) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(uint128_t, SubtractInt64) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        auto value1 = get_int64_random();
        auto value2 = get_int64_random();
        auto res = value1 - value2;
        uint128_t f1 = value1;
        uint128_t f3 = f1 - value2;
        uint64_t uint128_res = static_cast<uint64_t>(f3);
        EXPECT_EQ(uint128_res, res) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(uint128_t, SubtractUnsignedInt64) {
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
TEST(uint128_t, MultiplyByUint128) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = fabs(floor(get_double_random()));
        double value2 = fabs(floor(get_double_random()));
        double res = value1 * value2;
        uint128_t f1 = value1;
        uint128_t f2 = value2;
        uint128_t f3 = f1 * f2;
        if (check_overflow_uint128(value1) || check_overflow_uint128(value2) || check_overflow_uint128(res))
            continue;

        double uint128_res = static_cast<double>(f3);
        // note that uint128_res is more precise than double, this can lead to issues when param1 and param2 are far apart.
        if (is_similar_double(uint128_res, res))
            continue;

        EXPECT_EQ(uint128_res, res) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(uint128_t, DivideByUint128) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = fabs(floor(get_double_random()));
        double value2 = fabs(floor(get_double_random()));
        if (value2 == 0)
            continue;
        double res = floor(value1 / value2);
        uint128_t f1 = value1;
        uint128_t f2 = value2;
        uint128_t f3 = f1 / f2;
        if (check_overflow_uint128(value1) || check_overflow_uint128(value2) || check_overflow_uint128(res))
            continue;
        double uint128_res = static_cast<double>(f3);
        // note that uint128_res is more precise than double, this can lead to issues when param1 and param2 are far apart.
        if (is_similar_double(uint128_res, res))
            continue;

        EXPECT_DOUBLE_EQ(uint128_res, res) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(uint128_t, ModuloByUint128) {
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
TEST(uint128_t, CompareUint128) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = fabs(floor(get_double_random()));
        double value2 = fabs(floor(get_double_random()));
        bool res = value1 > value2;
        uint128_t f1 = value1;
        uint128_t f2 = value2;
        if (check_overflow_uint128(value1) || check_overflow_uint128(value2))
            continue;

        bool uint128_res = f1 > f2;
        EXPECT_TRUE(uint128_res == res) << "operator>: " << "value1=" << value1 << ", value2=" << value2;

        res = value1 >= value2;
        uint128_res = f1 >= f2;
        EXPECT_TRUE(uint128_res == res) << "operator>=: " << "value1=" << value1 << ", value2=" << value2;

        res = value1 < value2;
        uint128_res = f1 < f2;
        EXPECT_TRUE(uint128_res == res) << "operator<: " << "value1=" << value1 << ", value2=" << value2;

        res = value1 <= value2;
        uint128_res = f1 <= f2;
        EXPECT_TRUE(uint128_res == res) << "operator<=: " << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(uint128_t, OperatorPlusPlus) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        uint64_t value1 = get_uint64_random();
        uint64_t res = value1 + 1;
        uint128_t f1 = value1;
        uint128_t f2 = f1;
        f2++;

        EXPECT_EQ(static_cast<uint64_t>(f2), res) << "operator++(int)" << "value1=" << value1;

        f2 = f1;
        ++f2;

        EXPECT_EQ(static_cast<uint64_t>(f2), res) << "operator++()" << "value1=" << value1;
    }
}
TEST(uint128_t, OperatorMinusMinus) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        uint64_t value1 = get_uint64_random();
        uint64_t res = value1 - 1;
        uint128_t f1 = value1;
        uint128_t f2 = f1;
        f2--;

        EXPECT_EQ(static_cast<uint64_t>(f2), res) << "operator++(int)" << "value1=" << value1;

        f2 = f1;
        --f2;

        EXPECT_EQ(static_cast<uint64_t>(f2), res) << "operator++()" << "value1=" << value1;
    }
}
TEST(uint128_t, OperatorEqual) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = fabs(floor(get_double_random()));
        uint128_t f1 = value1;
        bool fp128_res = f1 == f1;

        if (check_overflow_uint128(value1))
            continue;

        EXPECT_TRUE(fp128_res == true) << "operator==: " << "value1=" << value1;
    }
}
TEST(uint128_t, OperatorNotEqual) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = fabs(floor(get_double_random()));
        uint128_t f1 = value1;
        bool fp128_res = f1 != f1;

        if (check_overflow_uint128(value1))
            continue;

        EXPECT_TRUE(fp128_res == false) << "operator==: " << "value1=" << value1;
    }
}
TEST(uint128_t, log10) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = floor(1.0 + fabs(get_double_random()));
        uint64_t res = (uint64_t)floor(log10(value1)); // double doesn't lose any bits with this operation!
        if (check_overflow_uint128(value1)) {
            continue;
        }
        uint128_t i1 = value1;
        uint64_t  i_res = log10(i1);
        EXPECT_EQ(i_res, res) << "double value1=" << value1;
    }
}
TEST(uint128_t, log) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = 1.0 + floor(fabs(get_double_random())) * get_uint64_random();
        uint64_t res = (uint64_t)floor(log(value1)); // double doesn't lose any bits with this operation!
        if (check_overflow_uint128(value1)) {
            continue;
        }
        uint128_t i1 = value1;
        uint64_t  i_res = log(i1);
        EXPECT_EQ(i_res, res) << "double value1=" << value1;
    }
}
TEST(uint128_t, sqrt) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = floor(1.0 + fabs(get_double_random()));
        uint64_t res = (uint64_t)floor(sqrt(value1)); // double doesn't lose any bits with this operation!
        if (check_overflow_uint128(value1)) {
            continue;
        }
        uint128_t i1 = value1;
        uint64_t  i_res = sqrt(i1);
        EXPECT_EQ(i_res, res) << "double value1=" << value1;
    }
}
TEST(uint128_t, pow) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = get_uint32_random() % 8;
        double value2 = get_uint32_random() % 16;
        double res = pow(value1, value2); // double doesn't lose any bits with this operation!
        if (check_overflow_uint128(value1) || check_overflow_uint128(res)) {
            continue;
        }
        uint128_t i1 = value1;
        uint128_t i2 = pow(i1, (uint32_t)value2);
        uint128_t uint128_res = static_cast<double>(i2);
        EXPECT_DOUBLE_EQ(uint128_res, res) << "value1=" << value1 << ", value2=" << value2;
    }
}
