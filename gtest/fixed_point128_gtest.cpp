#include <gtest/gtest.h>
#include <ostream>
#include <ctime>
#include "..\inc\fixed_point128.h"

/*************************************************
* Fixed point 128 tests
**************************************************/
using namespace fp128;

static constexpr int RANDOM_TEST_COUNT = 32768;
static constexpr int RANDOM_SEED = 0x12345678;

__forceinline int32_t get_random_sign()
{
    return (rand() > (RAND_MAX >> 1)) ? -1 : 1;
}

// returns a random number
double get_double_random()
{
    return (double)rand() * (double)rand() / (double)(rand() + 1) * get_random_sign();
}

// returns a positive random number
uint64_t get_uint64_random()
{
    return (uint64_t)rand() * (uint64_t)rand();
}

// returns a positive random number
uint64_t get_int64_random()
{
    return (int64_t)rand() * (int64_t)rand() * (int64_t)get_random_sign();
}

// returns a positive random number
uint32_t get_uint32_random()
{
    return (uint32_t)rand();
}

// returns a random number
uint32_t get_int32_random()
{
    return rand() * get_random_sign();
}

// Construct fixed_point128 and convert back to/from various elements.
TEST(fixed_point128, DefaultConstructor) {
    fixed_point128<20> f;
    EXPECT_EQ(static_cast<uint64_t>(f), 0ull);
}
TEST(fixed_point128, ConstructorFromDouble) {
    srand(RANDOM_SEED); // must have a repeatable seed for debugging
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random();
        fixed_point128<20> f = value;
        if (fabs(value) > f.max_int_value) {
            printf("Warning: test used a value which is too high: %lf\n", value);
            continue;
        }
        double f_value = static_cast<double>(f);
        if (fabs(value/ f_value) - 1.0 > 0.0001) {
            f_value = value;
        }
        EXPECT_DOUBLE_EQ(f_value, value);
    }
}
TEST(fixed_point128, ConstructorFromFloat) {
    srand(RANDOM_SEED); // must have a repeatable seed for debugging
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        float value = (float)get_double_random();
        fixed_point128<20> f = value;
        if (fabs(value) > f.max_int_value) {
            printf("Warning: test used a value which is too high: %f\n", value);
            continue;
        }
        EXPECT_FLOAT_EQ(static_cast<float>(f), value);
    }
}
TEST(fixed_point128, ConstructorFromInt32) {
    srand(RANDOM_SEED); // must have a repeatable seed for debugging
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        int32_t value = get_int32_random();
        fixed_point128<32> f = value;
        if (abs(value) > f.max_int_value) {
            printf("Warning: test used a value which is too high: %d\n", value);
            continue;
        }
        EXPECT_EQ(static_cast<int32_t>(f), value);
    }
}
TEST(fixed_point128, ConstructorFromUnsignedInt32) {
    srand(RANDOM_SEED); // must have a repeatable seed for debugging
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        uint32_t value = get_uint32_random();
        fixed_point128<32> f = value;
        if (value > f.max_int_value) {
            printf("Warning: test used a value which is too high: %u\n", value);
            continue;
        }
        EXPECT_EQ(static_cast<uint32_t>(f), value);
    }
}
TEST(fixed_point128, ConstructorFromInt64) {
    srand(RANDOM_SEED); // must have a repeatable seed for debugging
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        int64_t value = get_int64_random();
        fixed_point128<32> f = value;
        if (abs(value) > f.max_int_value) {
            printf("Warning: test used a value which is too high: %lld\n", value);
            continue;
        }
        EXPECT_EQ(static_cast<int64_t>(f), value);
    }
}
TEST(fixed_point128, ConstructorFromUnsignedInt64) {
    srand(RANDOM_SEED); // must have a repeatable seed for debugging
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        uint64_t value = get_uint64_random();
        fixed_point128<32> f = value;
        if (value > f.max_int_value) {
            printf("Warning: test used a value which is too high: %llu\n", value);
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
        EXPECT_DOUBLE_EQ(d1, d2);
    }
}
TEST(fixed_point128, CopyConstructor) {

    srand(RANDOM_SEED); // must have a repeatable seed for debugging
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random();
        fixed_point128<20> f1 = value;
        fixed_point128<20> f2 = f1;
        EXPECT_DOUBLE_EQ(static_cast<double>(f1), static_cast<double>(f2));
    }
}
TEST(fixed_point128, MoveConstructor) {
    srand(RANDOM_SEED); // must have a repeatable seed for debugging
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random();
        fixed_point128<20> f1 = value;
        fixed_point128<20> f2(std::move(f1));
        EXPECT_DOUBLE_EQ(static_cast<double>(f1), static_cast<double>(f2));
    }
}
TEST(fixed_point128, AssignmentOperator) {
    srand(RANDOM_SEED); // must have a repeatable seed for debugging
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random();
        fixed_point128<20> f1 = value;
        fixed_point128<20> f2;
        f2 = f1;
        EXPECT_DOUBLE_EQ(static_cast<double>(f1), static_cast<double>(f2));
    }
}
TEST(fixed_point128, AssignmentOperatorOtherType) {

    srand(RANDOM_SEED); // must have a repeatable seed for debugging
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random();
        fixed_point128<20> f1 = value;
        fixed_point128<22> f2;
        f2 = f1;
        EXPECT_DOUBLE_EQ(static_cast<double>(f1), static_cast<double>(f2));
    }
}
TEST(fixed_point128, MoveAssignmentOperator) {
    srand(RANDOM_SEED); // must have a repeatable seed for debugging
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random();
        fixed_point128<20> f1 = value;
        fixed_point128<20> f2;
        f2 = std::move(f1);
        EXPECT_DOUBLE_EQ(static_cast<double>(f1), static_cast<double>(f2));
    }
}
TEST(fixed_point128, CopyConstructorOtherType) {

    srand(RANDOM_SEED); // must have a repeatable seed for debugging
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random();
        fixed_point128<20> f1 = value;
        fixed_point128<22> f2 = f1;
        EXPECT_DOUBLE_EQ(static_cast<double>(f1), static_cast<double>(f2));
    }
}
