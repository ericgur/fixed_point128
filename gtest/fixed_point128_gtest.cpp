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
* fixed_point128 tests
***********************************************************************/
// Construct fixed_point128 and convert back to/from various elements.
TEST(fixed_point128, DefaultConstructor) {
    fixed_point128<20> f;
    EXPECT_EQ(static_cast<uint64_t>(f), 0ull);
}
TEST(fixed_point128, ConstructorFromDouble) {
    srand(RANDOM_SEED); 
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random(-43, 31);
        fixed_point128<32> f = value;
        double f_value = static_cast<double>(f);
        if (fabs(value/ f_value) - 1.0 > 0.0001) {
            f_value = value;
        }
        EXPECT_DOUBLE_EQ(f_value, value) << "value=" << value;
    }
}
TEST(fixed_point128, ConstructorFromFloat) {
    srand(RANDOM_SEED); 
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        float value = (float)get_double_random(-43, 31);
        fixed_point128<32> f = value;
        EXPECT_FLOAT_EQ(static_cast<float>(f), value);
    }
}
TEST(fixed_point128, ConstructorFromInt32) {
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
TEST(fixed_point128, ConstructorFromUnsignedInt32) {
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
TEST(fixed_point128, ConstructorFromInt64) {
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
TEST(fixed_point128, ConstructorFromUnsignedInt64) {
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
TEST(fixed_point128, ConstructorFromString) {
    const char* values[] = { "0", "12.34", "-3.1435", "12345.6789"};
    for (auto i = 0u; i < array_length(values); ++i) {
        fixed_point128<20> f = values[i];
        double d1 = strtod(values[i], nullptr);
        double d2 = strtod(static_cast<char*>(f), nullptr);
        EXPECT_DOUBLE_EQ(d1, d2) << "value=" << values[i];
    }
}
TEST(fixed_point128, CopyConstructor) {
    srand(RANDOM_SEED); 
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random();
        fixed_point128<20> f1 = value;
        fixed_point128<20> f2 = f1;
        EXPECT_DOUBLE_EQ(static_cast<double>(f1), static_cast<double>(f2)) << "value=" << value;;
    }
}
TEST(fixed_point128, MoveConstructor) {
    srand(RANDOM_SEED); 
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random();
        fixed_point128<20> f1 = value;
        fixed_point128<20> f2(std::move(f1));
        EXPECT_DOUBLE_EQ(static_cast<double>(f1), static_cast<double>(f2)) << "value=" << value;;
    }
}
TEST(fixed_point128, AssignmentOperator) {
    srand(RANDOM_SEED); 
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random();
        fixed_point128<20> f1 = value;
        fixed_point128<20> f2;
        f2 = f1;
        EXPECT_DOUBLE_EQ(static_cast<double>(f1), static_cast<double>(f2)) << "value=" << value;;
    }
}
TEST(fixed_point128, AssignmentOperatorOtherType) {
    srand(RANDOM_SEED); 
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random(-50, 18);
        fixed_point128<20> f1 = value;
        fixed_point128<22> f2;
        f2 = f1;
        EXPECT_DOUBLE_EQ(static_cast<double>(f1), static_cast<double>(f2)) << "value=" << value;;
        f1 = f2;
        EXPECT_DOUBLE_EQ(static_cast<double>(f1), static_cast<double>(f2)) << "value=" << value;;
    }
}
TEST(fixed_point128, MoveAssignmentOperator) {
    srand(RANDOM_SEED); 
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random(-50, 18);
        fixed_point128<20> f1 = value;
        fixed_point128<20> f2;
        f2 = std::move(f1);
        EXPECT_DOUBLE_EQ(static_cast<double>(f1), static_cast<double>(f2)) << "value=" << value;;
    }
}
TEST(fixed_point128, CopyConstructorOtherType) {
    srand(RANDOM_SEED); 
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random();
        fixed_point128<20> f1 = value;
        fixed_point128<22> f2 = f1;
        EXPECT_DOUBLE_EQ(static_cast<double>(f1), static_cast<double>(f2)) << "value=" << value;
        fixed_point128<20> f3 = f2;
        EXPECT_DOUBLE_EQ(static_cast<double>(f2), static_cast<double>(f3)) << "value=" << value;
    }
}
TEST(fixed_point128, AddSameSign) {
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
        double fp128_res = static_cast<double>(f3);
        // note that fp128 is more precise than double, this can lead to issues when param1 and param2 are far apart.
        if (is_similar_double(fp128_res, res))
            continue;

        EXPECT_DOUBLE_EQ(static_cast<double>(f3), res) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(fixed_point128, AddDifferentSign) {
    srand(RANDOM_SEED); 
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = fabs(get_double_random());
        double value2 = value1 * -2.5;
        double res = value1 + value2;
        fixed_point128<40> f1 = value1;
        fixed_point128<40> f2 = value2;
        fixed_point128<40> f3 = f1 + f2;
        if (check_overflow(value1, f1) || check_overflow(value2, f2) || check_overflow(res, f3)) {
            continue;
        }
        double fp128_res = static_cast<double>(f3);
        // note that fp128 is more precise than double, this can lead to issues when param1 and param2 are far apart.
        if (is_similar_double(fp128_res, res))
            continue;

        EXPECT_DOUBLE_EQ(static_cast<double>(f3), res) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(fixed_point128, AddDouble) {
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

        double fp128_res = static_cast<double>(f3);
        // note that fp128 is more precise than double, this can lead to issues when param1 and param2 are far apart.
        if (is_similar_double(fp128_res, res))
            continue;

        EXPECT_DOUBLE_EQ(static_cast<double>(f3), res) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(fixed_point128, AddFloat) {
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
TEST(fixed_point128, AddInt32) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        auto value1 = abs(get_int32_random());
        auto value2 = abs(get_int32_random());
        auto res = value1 + value2;
        if (res < 0) continue; //wrap around won't happen in uint128_t

        fixed_point128<40> f1 = value1;
        fixed_point128<40> f3 = f1 + value2;
        if (check_overflow(value1, f1) || check_overflow(value2, f1) || check_overflow(res, f3)) {
            continue;
        }
        EXPECT_EQ(static_cast<int64_t>(f3), res) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(fixed_point128, AddUnsignedInt32) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        auto value1 = get_uint32_random();
        auto value2 = get_uint32_random();
        auto res = value1 + value2;
        if (res < value1 || res < value2) continue; //wrap around won't happen in uint128_t

        fixed_point128<40> f1 = value1;
        fixed_point128<40> f3 = f1 + value2;
        if (check_overflow(value1, f1) || check_overflow(value2, f1) || check_overflow(res, f3)) {
            continue;
        }
        EXPECT_EQ(static_cast<uint64_t>(f3), res) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(fixed_point128, AddInt64) {
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
TEST(fixed_point128, AddUnsignedInt64) {
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
TEST(fixed_point128, SubtractSameSign) {
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
        double fp128_res = static_cast<double>(f3);
        // note that fp128 is more precise than double, this can lead to issues when param1 and param2 are far apart.
        if (is_similar_double(fp128_res, res))
            continue;

        EXPECT_DOUBLE_EQ(static_cast<double>(f3), res) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(fixed_point128, SubtractDifferentSign) {
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
        double fp128_res = static_cast<double>(f3);
        // note that fp128 is more precise than double, this can lead to issues when param1 and param2 are far apart.
        if (is_similar_double(fp128_res, res))
            continue;

        EXPECT_DOUBLE_EQ(static_cast<double>(f3), res) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(fixed_point128, SubtractInt32) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        auto value1 = abs(get_int32_random());
        auto value2 = abs(get_int32_random());
        if (value2 > value1) value2 = value1 / 3;

        auto res = value1 - value2;

        fixed_point128<40> f1 = value1;
        fixed_point128<40> f3 = f1 - value2;
        if (check_overflow(value1, f1) || check_overflow(value2, f1) || check_overflow(res, f3)) {
            continue;
        }
        EXPECT_EQ(static_cast<int64_t>(f3), res) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(fixed_point128, SubtractUnsignedInt32) {
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
TEST(fixed_point128, SubtractInt64) {
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
TEST(fixed_point128, SubtractUnsignedInt64) {
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
TEST(fixed_point128, MultiplyByFP128) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = get_double_random();
        double value2 = get_double_random();
        double res = value1 * value2;
        fixed_point128<40> f1 = value1;
        fixed_point128<40> f2 = value2;
        fixed_point128<40> f3 = f1 * f2;
        //assert((double)f1 == value1);
        //assert((double)f2 == value2);
        if (check_overflow(value1, f1) || check_overflow(value2, f1) || check_overflow(res, f3)) {
            continue;
        }

        double fp128_res = static_cast<double>(f3);
        // note that fp128 is more precise than double, this can lead to issues when param1 and param2 are far apart.
        if (is_similar_double(fp128_res, res))
            continue;

        EXPECT_DOUBLE_EQ(static_cast<double>(f3), res) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(fixed_point128, MultiplyByDouble) {
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
        double fp128_res = static_cast<double>(f3);
        // note that fp128 is more precise than double, this can lead to issues when param1 and param2 are far apart.
        if (is_similar_double(fp128_res, res))
            continue;

        EXPECT_DOUBLE_EQ(static_cast<double>(f3), res) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(fixed_point128, MultiplyByFloat) {
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
        double fp128_res = static_cast<double>(f3);
        // note that fp128 is more precise than double, this can lead to issues when param1 and param2 are far apart.
        if (is_similar_double(fp128_res, res))
            continue;

        EXPECT_FLOAT_EQ(static_cast<float>(f3), res) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(fixed_point128, MultiplyByInt32) {
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
TEST(fixed_point128, MultiplyByUnsignedInt32) {
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
TEST(fixed_point128, MultiplyByInt64) {
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
TEST(fixed_point128, MultiplyByUnsignedInt64) {
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
TEST(fixed_point128, DivideByFP128) {
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
        double fp128_res = static_cast<double>(f3);
        // note that fp128 is more precise than double, this can lead to issues when param1 and param2 are far apart.
        if (is_similar_double(fp128_res, res))
            continue;
        EXPECT_DOUBLE_EQ(fp128_res, res) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(fixed_point128, DivideByDouble) {
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
        double fp128_res = static_cast<double>(f3);
        // note that fp128 is more precise than double, this can lead to issues when param1 and param2 are far apart.
        if (is_similar_double(fp128_res, res))
            continue;

        EXPECT_DOUBLE_EQ(static_cast<double>(f3), res) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(fixed_point128, DivideByFloat) {
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
TEST(fixed_point128, DivideByInt32) {
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
TEST(fixed_point128, DivideByUnsignedInt32) {
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

TEST(fixed_point128, DivideByInt64) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = get_double_random(-1, 39);
        int64_t value2 = get_int32_random(); // on purpose
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
TEST(fixed_point128, DivideByUnsignedInt64) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = get_double_random(-1, 39);
        uint64_t value2 = get_uint32_random(); // on purpose
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
TEST(fixed_point128, ModuloByFP128) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = get_double_random(-1, 39);
        double value2 = get_double_random(-1, 39);
        //printf("%u\n", i);
        double res = value1 / value2;
        if (value2 == 0)
            continue;
        fixed_point128<40> f1 = value1;
        fixed_point128<40> f2 = value2;
        fixed_point128<40> f3 = f1 / f2;
        if (check_overflow(value1, f1) || check_overflow(value2, f1) || check_overflow(res, f3)) {
            continue;
        }
        double fp128_res = static_cast<double>(f3);
        // note that fp128 is more precise than double, this can lead to issues when param1 and param2 are far apart.
        if (is_similar_double(fp128_res, res))
            continue;

        EXPECT_DOUBLE_EQ(static_cast<double>(f3), res) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(fixed_point128, ModuloByDouble) {
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
        double fp128_res = static_cast<double>(f3);
        // note that fp128 is more precise than double, this can lead to issues when param1 and param2 are far apart.
        if (is_similar_double(fp128_res, res))
            continue;

        EXPECT_DOUBLE_EQ(static_cast<double>(f3), res) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(fixed_point128, ModuloByFloat) {
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
        double fp128_res = static_cast<double>(f3);
        // note that fp128 is more precise than double, this can lead to issues when param1 and param2 are far apart.
        if (is_similar_double(fp128_res, res))
            continue;
        EXPECT_FLOAT_EQ(static_cast<float>(f3), res) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(fixed_point128, ModuloByInt32) {
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
TEST(fixed_point128, ModuloByUnsignedInt32) {
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

TEST(fixed_point128, ModuloByInt64) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = get_double_random(-1, 39);
        int64_t value2 = get_int32_random(); // on purpose
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
TEST(fixed_point128, ModuloByUnsignedInt64) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = get_double_random(-1, 39);
        uint64_t value2 = get_uint32_random(); // on purpose
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
TEST(fixed_point128, CompareFP128) {
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
        EXPECT_TRUE(fp128_res == res) << "operator>: " << "value1=" << value1 << ", value2=" << value2;
        
        res = value1 >= value2;
        fp128_res = f1 >= f2;
        EXPECT_TRUE(fp128_res == res) << "operator>=: " << "value1=" << value1 << ", value2=" << value2;

        res = value1 < value2;
        fp128_res = f1 < f2;
        EXPECT_TRUE(fp128_res == res) << "operator<: " << "value1=" << value1 << ", value2=" << value2;

        res = value1 <= value2;
        fp128_res = f1 <= f2;
        EXPECT_TRUE(fp128_res == res) << "operator<=: " << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(fixed_point128, CompareDouble) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = get_double_random();
        double value2 = get_double_random();
        bool res = value1 > value2;
        fixed_point128<40> f1 = value1;
        if (check_overflow(value1, f1) || check_overflow(value2, f1))
            continue;

        bool fp128_res = f1 > value2;
        EXPECT_TRUE(fp128_res == res) << "operator>: " << "value1=" << value1 << ", value2=" << value2;

        res = value1 >= value2;
        fp128_res = f1 >= value2;
        EXPECT_TRUE(fp128_res == res) << "operator>=: " << "value1=" << value1 << ", value2=" << value2;

        res = value1 < value2;
        fp128_res = f1 < value2;
        EXPECT_TRUE(fp128_res == res) << "operator<: " << "value1=" << value1 << ", value2=" << value2;

        res = value1 <= value2;
        fp128_res = f1 <= value2;
        EXPECT_TRUE(fp128_res == res) << "operator<=: " << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(fixed_point128, CompareFloat) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        float value1 = (float)get_double_random();
        float value2 = (float)get_double_random();
        bool res = value1 > value2;
        fixed_point128<40> f1 = value1;
        if (check_overflow(value1, f1) || check_overflow(value2, f1))
            continue;

        bool fp128_res = f1 > value2;
        EXPECT_TRUE(fp128_res == res) << "operator>: " << "value1=" << value1 << ", value2=" << value2;

        res = value1 >= value2;
        fp128_res = f1 >= value2;
        EXPECT_TRUE(fp128_res == res) << "operator>=: " << "value1=" << value1 << ", value2=" << value2;

        res = value1 < value2;
        fp128_res = f1 < value2;
        EXPECT_TRUE(fp128_res == res) << "operator<: " << "value1=" << value1 << ", value2=" << value2;

        res = value1 <= value2;
        fp128_res = f1 <= value2;
        EXPECT_TRUE(fp128_res == res) << "operator<=: " << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(fixed_point128, CompareInt32) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        int32_t value1 = get_int32_random();
        int32_t value2 = get_int32_random();
        bool res = value1 > value2;
        fixed_point128<40> f1 = value1;
        if (check_overflow(value1, f1) || check_overflow(value2, f1))
            continue;

        bool fp128_res = f1 > value2;
        EXPECT_TRUE(fp128_res == res) << "operator>: " << "value1=" << value1 << ", value2=" << value2;

        res = value1 >= value2;
        fp128_res = f1 >= value2;
        EXPECT_TRUE(fp128_res == res) << "operator>=: " << "value1=" << value1 << ", value2=" << value2;

        res = value1 < value2;
        fp128_res = f1 < value2;
        EXPECT_TRUE(fp128_res == res) << "operator<: " << "value1=" << value1 << ", value2=" << value2;

        res = value1 <= value2;
        fp128_res = f1 <= value2;
        EXPECT_TRUE(fp128_res == res) << "operator<=: " << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(fixed_point128, CompareUnsignedInt32) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        uint32_t value1 = get_uint32_random();
        uint32_t value2 = get_uint32_random();
        bool res = value1 > value2;
        fixed_point128<40> f1 = value1;
        if (check_overflow(value1, f1) || check_overflow(value2, f1))
            continue;

        bool fp128_res = f1 > value2;
        EXPECT_TRUE(fp128_res == res) << "operator>: " << "value1=" << value1 << ", value2=" << value2;

        res = value1 >= value2;
        fp128_res = f1 >= value2;
        EXPECT_TRUE(fp128_res == res) << "operator>=: " << "value1=" << value1 << ", value2=" << value2;

        res = value1 < value2;
        fp128_res = f1 < value2;
        EXPECT_TRUE(fp128_res == res) << "operator<: " << "value1=" << value1 << ", value2=" << value2;

        res = value1 <= value2;
        fp128_res = f1 <= value2;
        EXPECT_TRUE(fp128_res == res) << "operator<=: " << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(fixed_point128, CompareInt64) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        int64_t value1 = get_int64_random();
        int64_t value2 = get_int64_random();
        bool res = value1 > value2;
        fixed_point128<40> f1 = value1;
        if (check_overflow(value1, f1) || check_overflow(value2, f1))
            continue;

        bool fp128_res = f1 > value2;
        EXPECT_TRUE(fp128_res == res) << "operator>: " << "value1=" << value1 << ", value2=" << value2;

        res = value1 >= value2;
        fp128_res = f1 >= value2;
        EXPECT_TRUE(fp128_res == res) << "operator>=: " << "value1=" << value1 << ", value2=" << value2;

        res = value1 < value2;
        fp128_res = f1 < value2;
        EXPECT_TRUE(fp128_res == res) << "operator<: " << "value1=" << value1 << ", value2=" << value2;

        res = value1 <= value2;
        fp128_res = f1 <= value2;
        EXPECT_TRUE(fp128_res == res) << "operator<=: " << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(fixed_point128, CompareUnsignedInt64) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        uint64_t value1 = get_uint64_random();
        uint64_t value2 = get_uint64_random();
        bool res = value1 > value2;
        fixed_point128<40> f1 = value1;
        if (check_overflow(value1, f1) || check_overflow(value2, f1))
            continue;

        bool fp128_res = f1 > value2;
        EXPECT_TRUE(fp128_res == res) << "operator>: " << "value1=" << value1 << ", value2=" << value2;

        res = value1 >= value2;
        fp128_res = f1 >= value2;
        EXPECT_TRUE(fp128_res == res) << "operator>=: " << "value1=" << value1 << ", value2=" << value2;

        res = value1 < value2;
        fp128_res = f1 < value2;
        EXPECT_TRUE(fp128_res == res) << "operator<: " << "value1=" << value1 << ", value2=" << value2;

        res = value1 <= value2;
        fp128_res = f1 <= value2;
        EXPECT_TRUE(fp128_res == res) << "operator<=: " << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(fixed_point128, OperatorPlusPlus) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = get_double_random();
        double res = value1 + 1;
        fixed_point128<40> f1 = value1;
        fixed_point128<40> f2 = f1;
        f2++;

        if (check_overflow(value1+2, f1))
            continue;

        EXPECT_DOUBLE_EQ(static_cast<double>(f2), res) << "operator++(int)" << "value1=" << value1;

        f2 = f1;
        ++f2;

        EXPECT_DOUBLE_EQ(static_cast<double>(f2), res) << "operator++()" << "value1=" << value1;
    }
}
TEST(fixed_point128, OperatorMinusMinus) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = get_double_random();
        double res = value1 - 1;
        fixed_point128<40> f1 = value1;
        fixed_point128<40> f2 = f1;
        f2--;

        if (check_overflow(value1 + 2, f1) || check_overflow(value1 - 2, f1))
            continue;

        EXPECT_DOUBLE_EQ(static_cast<double>(f2), res) << "operator--(int)" << "value1=" << value1;

        f2 = f1;
        --f2;

        EXPECT_DOUBLE_EQ(static_cast<double>(f2), res) << "operator--()" << "value1=" << value1;
    }
}
TEST(fixed_point128, OperatorEqual) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = get_double_random();
        fixed_point128<40> f1 = value1;
        bool fp128_res = f1 == f1;

        if (check_overflow(value1, f1))
            continue;

        EXPECT_TRUE(fp128_res == true) << "operator==: " << "value1=" << value1;
    }
}
TEST(fixed_point128, TemplateOperatorEqual) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = get_double_random();
        fixed_point128<40> f1 = value1;
        bool fp128_res = f1 == value1;

        if (check_overflow(value1, f1))
            continue;

        EXPECT_TRUE(fp128_res == true) << "operator==<double>: " << "value1=" << value1;

        f1 = static_cast<float>(value1);
        fp128_res = f1 == static_cast<float>(value1);

        EXPECT_TRUE(fp128_res == true) << "operator==<float>: " << "value1=" << value1;

        f1 = static_cast<uint64_t>(value1);
        fp128_res = f1 == static_cast<uint64_t>(value1);

        EXPECT_TRUE(fp128_res == true) << "operator==<uint64_t>: " << "value1=" << value1;

        f1 = static_cast<int64_t>(value1);
        fp128_res = f1 == static_cast<int64_t>(value1);

        EXPECT_TRUE(fp128_res == true) << "operator==<int64_t>: " << "value1=" << value1;

        f1 = static_cast<uint32_t>(value1);
        fp128_res = f1 == static_cast<uint32_t>(value1);

        EXPECT_TRUE(fp128_res == true) << "operator==<uint32_t>: " << "value1=" << value1;

        f1 = static_cast<int32_t>(value1);
        fp128_res = f1 == static_cast<int32_t>(value1);

        EXPECT_TRUE(fp128_res == true) << "operator==<int32_t>: " << "value1=" << value1;
    }
}
TEST(fixed_point128, OperatorNotEqual) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = get_double_random();
        fixed_point128<40> f1 = value1;
        bool fp128_res = f1 != f1;

        if (check_overflow(value1, f1))
            continue;

        EXPECT_TRUE(fp128_res == false) << "operator!=: " << "value1=" << value1;
    }
}
TEST(fixed_point128, TemplateOperatorNotEqual) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = get_double_random();
        fixed_point128<40> f1 = value1;
        bool fp128_res = f1 != value1;

        if (check_overflow(value1, f1))
            continue;

        EXPECT_TRUE(fp128_res == false) << "operator!=<double>: " << "value1=" << value1;

        f1 = static_cast<float>(value1);
        fp128_res = f1 != static_cast<float>(value1);

        EXPECT_TRUE(fp128_res == false) << "operator!=<float>: " << "value1=" << value1;

        f1 = static_cast<uint64_t>(value1);
        fp128_res = f1 != static_cast<uint64_t>(value1);

        EXPECT_TRUE(fp128_res == false) << "operator!=<uint64_t>: " << "value1=" << value1;

        f1 = static_cast<int64_t>(value1);
        fp128_res = f1 != static_cast<int64_t>(value1);

        EXPECT_TRUE(fp128_res == false) << "operator!=<int64_t>: " << "value1=" << value1;

        f1 = static_cast<uint32_t>(value1);
        fp128_res = f1 != static_cast<uint32_t>(value1);

        EXPECT_TRUE(fp128_res == false) << "operator!=<uint32_t>: " << "value1=" << value1;

        f1 = static_cast<int32_t>(value1);
        fp128_res = f1 != static_cast<int32_t>(value1);

        EXPECT_TRUE(fp128_res == false) << "operator!=<int32_t>: " << "value1=" << value1;
    }
}
TEST(fixed_point128, ShiftRight) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = get_uint32_random();
        uint32_t shift = get_uint32_random() % 40u;
        fixed_point128<40> f1 = value1;

        double expo = pow(2, -static_cast<int32_t>(shift));
        double res = value1 * expo; // double doesn't lose any bits with this operation!
        fixed_point128 f3 = f1 >> shift;

        if (check_overflow(value1, f1) || check_overflow(res, f3)) {
            continue;
        }
        EXPECT_DOUBLE_EQ(static_cast<double>(f3), res) << "double value1=" << value1 << "; uint32_t shift=" << shift << ";";
    }
}
TEST(fixed_point128, ShiftLeft) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = 1.0 / get_int32_random();
        uint32_t shift = get_uint32_random() % 40u;
        fixed_point128<40> f1 = value1;
        double expo = pow(2, static_cast<int32_t>(shift));
        double res = value1 * expo; // double doesn't lose any bits with this operation!
        fixed_point128 f3 = f1 << shift;
        if (check_overflow(value1, f1) || check_overflow(res, f3))
            continue;
        EXPECT_DOUBLE_EQ(static_cast<double>(f3), res) << "double value1=" << value1 << "; uint32_t shift=" << shift << ";";
    }
}
TEST(fixed_point128, floor) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random(-15, 15); // can't exceed this range to avoid overflow
        double res = floor(value);
        fixed_point128<16> f1 = value;
        fixed_point128<16> fp128_res = floor(f1);
        EXPECT_DOUBLE_EQ(fp128_res, res) << "floor: " << "value=" << value;
    }
}
TEST(fixed_point128, ceil) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random(-15, 15); // can't exceed this range to avoid overflow
        double res = ceil(value);
        fixed_point128<16> f1 = value;
        fixed_point128<16> fp128_res = ceil(f1);
        EXPECT_DOUBLE_EQ(fp128_res, res) << "ceil: " << "value=" << value;
    }
}
TEST(fixed_point128, trunc) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random(-15, 15); // can't exceed this range to avoid overflow
        double res = trunc(value);
        fixed_point128<16> f1 = value;
        fixed_point128<16> fp128_res = trunc(f1);
        EXPECT_DOUBLE_EQ(fp128_res, res) << "trunc: " << "value=" << value;
    }
}
TEST(fixed_point128, round) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random(-15, 15); // can't exceed this range to avoid overflow
        double res = round(value);
        fixed_point128<16> f1 = value;
        fixed_point128<16> fp128_res = round(f1);
        EXPECT_DOUBLE_EQ(fp128_res, res) << "round: " << "value=" << value;
    }
}
TEST(fixed_point128, copysign) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = get_double_random(-15, 15); // can't exceed this range to avoid overflow
        double value2 = get_double_random(-15, 15); // can't exceed this range to avoid overflow
        double res = copysign(value1, value2);
        fixed_point128<16> f1 = value1;
        fixed_point128<16> f2 = value2;
        fixed_point128<16> fp128_res = copysign(f1, f2);
        EXPECT_DOUBLE_EQ(fp128_res, res) << "copysign: " << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(fixed_point128, reciprocal) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random(-15, 15); // can't exceed this range to avoid overflow
        double res = 1.0 / value;
        fixed_point128<16> f1 = value;
        fixed_point128<16> fp128_res = reciprocal(f1);
        EXPECT_DOUBLE_EQ(fp128_res, res) << "reciprocal: " << "value=" << value;
    }
}
TEST(fixed_point128, log) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = fabs(get_double_random(-40, 15)); // lower exponent results in lost bits
        double res = log(value);
        fixed_point128<16> f1 = value;
        fixed_point128<16> fp128_res = log(f1);
        EXPECT_DOUBLE_EQ(fp128_res, res) << "log: " << "value=" << value;
    }
}
TEST(fixed_point128, log2) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = fabs(get_double_random(-40, 15)); // lower exponent results in lost bits
        double res = log2(value);
        fixed_point128<16> f1 = value;
        fixed_point128<16> fp128_res = log2(f1);
        EXPECT_DOUBLE_EQ(fp128_res, res) << "log2: " << "value=" << value;
    }
}
TEST(fixed_point128, log10) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = fabs(get_double_random(-40, 15)); // lower exponent results in lost bits
        double res = log10(value);
        fixed_point128<16> f1 = value;
        fixed_point128<16> fp128_res = log10(f1);
        EXPECT_DOUBLE_EQ(fp128_res, res) << "log10: " << "value=" << value;
    }
}
TEST(fixed_point128, log1p) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = fabs(get_double_random(-40, 15)); // lower exponent results in lost bits
        double res = log1p(value);
        fixed_point128<16> f1 = value;
        fixed_point128<16> fp128_res = log1p(f1);
        EXPECT_DOUBLE_EQ(fp128_res, res) << "log1p: " << "value=" << value;
    }
}
TEST(fixed_point128, logb) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = fabs(get_double_random(-40, 15)); // lower exponent results in lost bits
        double res = logb(value);
        fixed_point128<16> f1 = value;
        fixed_point128<16> fp128_res = logb(f1);
        EXPECT_DOUBLE_EQ(fp128_res, res) << "logb: " << "value=" << value;
    }
}
TEST(fixed_point128, sin) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random(-60, 2);
        double res = ::sin(value);
        fixed_point128<16> f1 = value;
        fixed_point128<16> fp128_res = sin(f1);
        EXPECT_DOUBLE_EQ(fp128_res, res) << "sin: " << "value=" << value;
    }
}
TEST(fixed_point128, cos) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random(-60, 2);
        double res = ::cos(value);
        fixed_point128<16> f1 = value;
        fixed_point128<16> fp128_res = cos(f1);
        EXPECT_DOUBLE_EQ(fp128_res, res) << "cos: " << "value=" << value;
    }
}
TEST(fixed_point128, tan) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random(-60, 2);
        double res = ::tan(value);
        fixed_point128<16> f1 = value;
        fixed_point128<16> fp128_res = tan(f1);
        EXPECT_DOUBLE_EQ(fp128_res, res) << "tan: " << "value=" << value;
    }
}
TEST(fixed_point128, asin) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random(-60, -1);
        double res = ::asin(value);
        fixed_point128<16> f1 = value;
        fixed_point128<16> fp128_res = asin(f1);
        EXPECT_DOUBLE_EQ(fp128_res, res) << "asin: " << "value=" << value;
    }
}
TEST(fixed_point128, acos) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random(-60, -1);
        double res = acos(value);
        fixed_point128<16> f1 = value;
        fixed_point128<16> fp128_res = acos(f1);
        EXPECT_DOUBLE_EQ(fp128_res, res) << "acos: " << "value=" << value;
    }
}
TEST(fixed_point128, atan) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random(-60, 14); //lower exponent results in lost bits
        double res = atan(value);
        fixed_point128<16> f1 = value;
        fixed_point128<16> fp128_res = atan(f1);
        EXPECT_DOUBLE_EQ(fp128_res, res) << "atan: " << "value=" << value;
    }
}
TEST(fixed_point128, atan2) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = get_double_random(-10, 14); //lower exponent results in lost bits
        double value2 = get_double_random(-10, 14); //lower exponent results in lost bits
        double res = atan2(value1, value2);
        fixed_point128<16> f1 = value1;
        fixed_point128<16> f2 = value2;
        fixed_point128<16> fp128_res = atan2(f1, f2);
        EXPECT_DOUBLE_EQ(fp128_res, res) << "atan2: " << " value1=" << value1 << ", value2=" << value2;
    }
}
TEST(fixed_point128, exp) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random(-3, 3); // lower exponent results in lost bits
        double res = exp(value);
        fixed_point128<16> f1 = value;
        fixed_point128<16> fp128_res = exp(f1);
        EXPECT_DOUBLE_EQ(fp128_res, res) << "exp: " << "value=" << value;
    }
}
TEST(fixed_point128, exp2) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random(-4, 4); // lower exponent results in lost bits
        double res = exp2(value);
        fixed_point128<16> f1 = value;
        fixed_point128<16> fp128_res = exp2(f1);
        EXPECT_DOUBLE_EQ(fp128_res, res) << "exp2: " << "value=" << value;
    }
}
TEST(fixed_point128, expm1) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random(-3, 3); // lower exponent results in lost bits
        double res = expm1(value);
        fixed_point128<16> f1 = value;
        fixed_point128<16> fp128_res = expm1(f1);
        EXPECT_DOUBLE_EQ(fp128_res, res) << "expm1: " << "value=" << value;
    }
}
TEST(fixed_point128, pow) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = fabs(get_double_random(-4, 4)); // lower exponent results in lost bits
        double value2 = get_double_random(-4, 4); // lower exponent results in lost bits
        value1 = 2.3505350060521568, value2 = 0.20224790564651351; //TODO: remove
        double res = pow(value1, value2);
        fixed_point128<16> f1 = value1;
        fixed_point128<16> f2 = value2;
        fixed_point128<16> fp128_res = pow(f1, f2);
        EXPECT_DOUBLE_EQ(fp128_res, res) << "pow: " << " value1=" << value1 << ", value2=" << value2;
    }
}
TEST(fixed_point128, sinh) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random(-60, 2);
        double res = ::sinh(value);
        fixed_point128<16> f1 = value;
        fixed_point128<16> fp128_res = sinh(f1);
        EXPECT_DOUBLE_EQ(fp128_res, res) << "sinh: " << "value=" << value;
    }
}
TEST(fixed_point128, asinh) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random(-60, 2);
        double res = ::asinh(value);
        fixed_point128<16> f1 = value;
        fixed_point128<16> fp128_res = asinh(f1);
        EXPECT_DOUBLE_EQ(fp128_res, res) << "asinh: " << "value=" << value;
    }
}
TEST(fixed_point128, cosh) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random(-60, 2);
        double res = ::cosh(value);
        fixed_point128<16> f1 = value;
        fixed_point128<16> fp128_res = cosh(f1);
        EXPECT_DOUBLE_EQ(fp128_res, res) << "cosh: " << "value=" << value;
    }
}
TEST(fixed_point128, acosh) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = 1.0 + fabs(get_double_random(-60, 2)); // must be >= 1
        double res = ::acosh(value);
        fixed_point128<16> f1 = value;
        fixed_point128<16> fp128_res = acosh(f1);
        EXPECT_DOUBLE_EQ(fp128_res, res) << "acosh: " << "value=" << value;
    }
}
TEST(fixed_point128, tanh) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random(-60, 2);
        double res = ::tanh(value);
        fixed_point128<16> f1 = value;
        fixed_point128<16> fp128_res = tanh(f1);
        EXPECT_DOUBLE_EQ(fp128_res, res) << "tanh: " << "value=" << value;
    }
}
