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
* int128_t tests
***********************************************************************/
// Construct fixed_point128 and convert back to/from various elements.
TEST(int128_t, DefaultConstructor) {
    int128_t i;
    EXPECT_EQ(static_cast<int64_t>(i), 0ull);
}
TEST(int128_t, ConstructorFromDouble) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = floor(get_double_random());
        int128_t f = value;
        if (check_overflow_int128(value)) {
            continue;
        }
        double f_value = static_cast<double>(f);
        EXPECT_DOUBLE_EQ(f_value, value) << "value=" << value;
    }
}
TEST(int128_t, ConstructorFromFloat) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        float value = (float)floor(get_double_random());
        int128_t f = value;
        if (check_overflow_int128(value)) {
            continue;
        }
        EXPECT_FLOAT_EQ(static_cast<float>(f), value) << "value=" << value;
    }
}
TEST(int128_t, ConstructorFromInt32) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        int32_t value = get_int32_random();
        int128_t f = value;
        EXPECT_EQ(static_cast<int32_t>(f), value);
    }
}
TEST(int128_t, ConstructorFromUnsignedInt32) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        uint32_t value = get_uint32_random();
        int128_t f = value;
        EXPECT_EQ(static_cast<uint32_t>(f), value);
    }
}
TEST(int128_t, ConstructorFromInt64) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        int64_t value = get_int64_random();
        int128_t f = value;
        EXPECT_EQ(static_cast<int64_t>(f), value);
    }
}
TEST(int128_t, ConstructorFromUnsignedInt64) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        int64_t value = get_int64_random();
        int128_t f = value;
        EXPECT_EQ(static_cast<int64_t>(f), value);
    }
}
TEST(int128_t, ConstructorFromString) {
    const char* values[] = { "0xDEADBEAF", "1234", "123456789012345678" };
    for (auto i = 0u; i < array_length(values); ++i) {
        int128_t f = values[i];
        double d1 = strtod(values[i], nullptr);
        double d2 = strtod(static_cast<char*>(f), nullptr);
        EXPECT_DOUBLE_EQ(d1, d2);
    }
}
TEST(int128_t, CopyConstructor) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random();
        int128_t f1 = value;
        int128_t f2 = f1;
        EXPECT_DOUBLE_EQ(static_cast<double>(f1), static_cast<double>(f2));
    }
}
TEST(int128_t, MoveConstructor) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random();
        int128_t f1 = value;
        int128_t f2(std::move(f1));
        EXPECT_DOUBLE_EQ(static_cast<double>(f1), static_cast<double>(f2));
    }
}
TEST(int128_t, AssignmentOperator) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random();
        int128_t f1 = value;
        int128_t f2;
        f2 = f1;
        EXPECT_DOUBLE_EQ(static_cast<double>(f1), static_cast<double>(f2));
    }
}
TEST(int128_t, MoveAssignmentOperator) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random();
        int128_t f1 = value;
        int128_t f2;
        f2 = std::move(f1);
        EXPECT_DOUBLE_EQ(static_cast<double>(f1), static_cast<double>(f2));
    }
}
TEST(int128_t, Add) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = floor(get_double_random());
        double value2 = floor(value1 * 2.5);
        double res = value1 + value2;
        int128_t f1 = value1;
        int128_t f2 = value2;
        int128_t f3 = f1 + f2;
        if (check_overflow_int128(value1) || check_overflow_int128(value2) || check_overflow_int128(res))
            continue;
        double int128_res = static_cast<double>(f3);
        // note that int128_t is more precise than double, this can lead to issues when param1 and param2 are far apart.
        if (is_similar_double(int128_res, res))
            continue;

        EXPECT_DOUBLE_EQ(int128_res, res) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(int128_t, AddDifferentSign) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        int64_t value1 = get_int64_random();
    #pragma warning(push)
    #pragma warning(disable: 4146)
        int64_t value2 = -get_int64_random();
    #pragma warning(pop)
        int64_t res = value1 + value2;
        int128_t f1 = value1;
        int128_t f2 = value2;
        int128_t f3 = f1 + f2;
        int64_t int128_res = static_cast<int64_t>(f3);
        EXPECT_EQ(int128_res, res) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(int128_t, AddInt64) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        auto value1 = get_int64_random();
        auto value2 = get_int64_random();
        auto res = value1 + value2;
        int128_t f1 = value1;
        int128_t f3 = f1 + value2;
        int64_t int128_res = static_cast<int64_t>(f3);
        EXPECT_EQ(int128_res, res) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(int128_t, AddUnsignedInt64) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        auto value1 = get_int64_random();
        auto value2 = get_int64_random();
        auto res = value1 + value2;
        int128_t f1 = value1;
        int128_t f3 = f1 + value2;
        int64_t int128_res = static_cast<int64_t>(f3);
        EXPECT_EQ(int128_res, res) << "value1=" << value1 << ", value2=" << value2;
    }
}

TEST(int128_t, Subtract) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = floor(get_double_random());
        double value2 = floor(value1 / 2.5);
        double res = value1 - value2;
        int128_t f1 = value1;
        int128_t f2 = value2;
        int128_t f3 = f1 - f2;
        if (check_overflow_int128(value1) || check_overflow_int128(value2) || check_overflow_int128(res))
            continue;
        double int128_res = static_cast<double>(f3);
        // note that fp128 is more precise than double, this can lead to issues when param1 and param2 are far apart.
        if (is_similar_double(int128_res, res))
            continue;

        EXPECT_DOUBLE_EQ(int128_res, res) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(int128_t, SubtractDifferentSign) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        int64_t value1 = get_int64_random();
    #pragma warning(push)
    #pragma warning(disable: 4146)
        int64_t value2 = -get_int64_random();
    #pragma warning(pop)
        int64_t res = value1 - value2;
        int128_t f1 = value1;
        int128_t f2 = value2;
        int128_t f3 = f1 - f2;
        int64_t int128_res = static_cast<int64_t>(f3);
        EXPECT_EQ(int128_res, res) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(int128_t, SubtractInt64) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        auto value1 = get_int64_random();
        auto value2 = get_int64_random();
        auto res = value1 - value2;
        int128_t f1 = value1;
        int128_t f3 = f1 - value2;
        int64_t int128_res = static_cast<int64_t>(f3);
        EXPECT_EQ(int128_res, res) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(int128_t, SubtractUnsignedInt64) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        auto value1 = get_int64_random();
        auto value2 = get_int64_random();
        auto res = value1 - value2;
        int128_t f1 = value1;
        int128_t f3 = f1 - value2;
        int64_t int128_res = static_cast<int64_t>(f3);
        EXPECT_EQ(int128_res, res) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(int128_t, MultiplyByint128) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = floor(get_double_random());
        double value2 = floor(get_double_random());
        double res = value1 * value2;
        int128_t f1 = value1;
        int128_t f2 = value2;
        int128_t f3 = f1 * f2;
        if (check_overflow_int128(value1) || check_overflow_int128(value2) || check_overflow_int128(res))
            continue;

        double int128_res = static_cast<double>(f3);
        // note that int128_res is more precise than double, this can lead to issues when param1 and param2 are far apart.
        if (is_similar_double(int128_res, res))
            continue;

        EXPECT_EQ(int128_res, res) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(int128_t, DivideByint128) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = floor(get_double_random());
        double value2 = floor(get_double_random());
        if (value2 == 0)
            continue;
        double res = floor(value1 / value2);
        int128_t f1 = value1;
        int128_t f2 = value2;
        int128_t f3 = f1 / f2;
        if (check_overflow_int128(value1) || check_overflow_int128(value2) || check_overflow_int128(res))
            continue;
        double int128_res = static_cast<double>(f3);
        // note that int128_res is more precise than double, this can lead to issues when param1 and param2 are far apart.
        if (is_similar_double(int128_res, res))
            continue;

        EXPECT_DOUBLE_EQ(int128_res, res) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(int128_t, ModuloByint128) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        auto value1 = get_int64_random();
        auto value2 = get_int32_random();

        if (value2 == 0)
            continue;
        int64_t res = value1 % value2;
        int128_t f1 = value1;
        int128_t f2 = value2;
        int128_t f3 = f1 % f2;
        int64_t int128_res = static_cast<int64_t>(f3);
        if (int128_res == res) continue;

        EXPECT_EQ(int128_res, res) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(int128_t, Compareint128) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = floor(get_double_random());
        double value2 = floor(get_double_random());
        bool res = value1 > value2;
        int128_t f1 = value1;
        int128_t f2 = value2;
        if (check_overflow_int128(value1) || check_overflow_int128(value2))
            continue;

        bool int128_res = f1 > f2;
        EXPECT_TRUE(int128_res == res) << "operator>: " << "value1=" << value1 << ", value2=" << value2;

        res = value1 >= value2;
        int128_res = f1 >= f2;
        EXPECT_TRUE(int128_res == res) << "operator>=: " << "value1=" << value1 << ", value2=" << value2;

        res = value1 < value2;
        int128_res = f1 < f2;
        EXPECT_TRUE(int128_res == res) << "operator<: " << "value1=" << value1 << ", value2=" << value2;

        res = value1 <= value2;
        int128_res = f1 <= f2;
        EXPECT_TRUE(int128_res == res) << "operator<=: " << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(int128_t, OperatorPlusPlus) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        int64_t value1 = get_int64_random();
        int64_t res = value1 + 1;
        int128_t f1 = value1;
        int128_t f2 = f1;
        f2++;

        EXPECT_EQ(static_cast<int64_t>(f2), res) << "operator++(int)" << "value1=" << value1;

        f2 = f1;
        ++f2;

        EXPECT_EQ(static_cast<int64_t>(f2), res) << "operator++()" << "value1=" << value1;
    }
}
TEST(int128_t, OperatorMinusMinus) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        int64_t value1 = get_int64_random();
        int64_t res = value1 - 1;
        int128_t f1 = value1;
        int128_t f2 = f1;
        f2--;

        EXPECT_EQ(static_cast<int64_t>(f2), res) << "operator++(int)" << "value1=" << value1;

        f2 = f1;
        --f2;

        EXPECT_EQ(static_cast<int64_t>(f2), res) << "operator++()" << "value1=" << value1;
    }
}
TEST(int128_t, OperatorEqual) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = floor(get_double_random());
        int128_t f1 = value1;
        bool fp128_res = f1 == f1;

        if (check_overflow_int128(value1))
            continue;

        EXPECT_TRUE(fp128_res == true) << "operator==: " << "value1=" << value1;
    }
}
TEST(int128_t, OperatorNotEqual) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = floor(get_double_random());
        int128_t f1 = value1;
        bool fp128_res = f1 != f1;

        if (check_overflow_int128(value1))
            continue;

        EXPECT_TRUE(fp128_res == false) << "operator==: " << "value1=" << value1;
    }
}
TEST(int128_t, log10) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = floor(1.0 + fabs(get_double_random()));
        int64_t res = (int64_t)floor(log10(value1)); // double doesn't lose any bits with this operation!
        if (check_overflow_int128(value1)) {
            continue;
        }
        int128_t i1 = value1;
        int64_t  i_res = log10(i1);
        EXPECT_EQ(i_res, res) << "double value1=" << value1;
    }
}
TEST(int128_t, log) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = 1.0 + floor(fabs(get_double_random(0, 127)));
        int64_t res = (int64_t)floor(log(value1)); // double doesn't lose any bits with this operation!
        if (check_overflow_int128(value1)) {
            continue;
        }
        int128_t i1 = value1;
        int64_t  i_res = log(i1);
        EXPECT_EQ(i_res, res) << "double value1=" << value1;
    }
}
TEST(int128_t, sqrt) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = floor(1.0 + fabs(get_double_random()));
        int64_t res = (int64_t)floor(sqrt(value1)); // double doesn't lose any bits with this operation!
        if (check_overflow_int128(value1)) {
            continue;
        }
        int128_t i1 = value1;
        int64_t  i_res = sqrt(i1);
        EXPECT_EQ(i_res, res) << "double value1=" << value1;
    }
}
TEST(int128_t, pow) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = get_uint32_random() % 8;
        double value2 = get_uint32_random() % 16;
        double res = pow(value1, value2); // double doesn't lose any bits with this operation!
        if (check_overflow_int128(value1) || check_overflow_int128(res)) {
            continue;
        }
        int128_t i1 = value1;
        int128_t i2 = pow(i1, (uint32_t)value2);
        int128_t int128_res = static_cast<double>(i2);
        EXPECT_DOUBLE_EQ(int128_res, res) << "value1=" << value1 << ", value2=" << value2;
    }
}
