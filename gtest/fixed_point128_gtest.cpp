#include <gtest/gtest.h>
#include <ostream>
//#include <gmock/gmock.h>
#include "..\inc\fixed_point128.h"

/*************************************************
* Fixed point 128 tests
**************************************************/
using namespace fp128;
//using ::testing::Gt;
//using ::testing::Lt;
//using ::testing::MatchesRegex;
//using ::testing::StartsWith;
//using ::testing::DoubleEq;


// Construct fixed_point128 and convert back to/from various elements.
TEST(fixed_point128, ConstructorFromDouble) {
    double values[] = { 0.0, 1.0, 2.5, 0.001, 255.000001, 1e38 };
    for (auto i = 0u; i < array_length(values); ++i) {
        fixed_point128<20> f = values[i];
        if (values[i] > f.max_int_value) {
            printf("Warning: test used a value which is too high: %lf\n", values[i]);
            continue;
        }

        EXPECT_DOUBLE_EQ(static_cast<double>(f), values[i]);
    }
}
TEST(fixed_point128, ConstructorFromFloat) {
    float values[] = { 0.0f, 1.0f, 2.5f, 0.001f, 255.000001f };
    for (auto i = 0u; i < array_length(values); ++i) {
        fixed_point128<20> f = values[i];
        if (values[i] > f.max_int_value) {
            printf("Warning: test used a value which is too high: %f\n", values[i]);
            continue;
        }
        EXPECT_FLOAT_EQ(static_cast<float>(f), values[i]);
    }
}
TEST(fixed_point128, ConstructorFromInt32) {
    int32_t values[] = { 0, 10, -1543 };
    for (auto i = 0u; i < array_length(values); ++i) {
        fixed_point128<20> f = values[i];
        if (abs(values[i]) > f.max_int_value) {
            printf("Warning: test used a value which is too high: %d\n", values[i]);
            continue;
        }
        EXPECT_EQ(static_cast<int32_t>(f), values[i]);
    }
}
TEST(fixed_point128, ConstructorFromUnsignedInt32) {
    uint32_t values[] = { 0u, 10u, 1543u };
    for (auto i = 0u; i < array_length(values); ++i) {
        fixed_point128<20> f = values[i];
        if (values[i] > f.max_int_value) {
            printf("Warning: test used a value which is too high: %u\n", values[i]);
            continue;
        }
        EXPECT_EQ(static_cast<uint32_t>(f), values[i]);
    }
}
TEST(fixed_point128, ConstructorFromInt64) {
    int64_t values[] = { 0ll, 10ll, 1543ll };
    for (auto i = 0u; i < array_length(values); ++i) {
        fixed_point128<20> f = values[i];
        if (abs(values[i]) > f.max_int_value) {
            printf("Warning: test used a value which is too high: %lld\n", values[i]);
            continue;
        }
        EXPECT_EQ(static_cast<int64_t>(f), values[i]);
    }
}
TEST(fixed_point128, ConstructorFromUnsignedInt64) {
    uint64_t values[] = { 0ull, 10ull, 1543ull };
    for (auto i = 0u; i < array_length(values); ++i) {
        fixed_point128<20> f = values[i];
        if (values[i] > f.max_int_value) {
            printf("Warning: test used a value which is too high: %llu\n", values[i]);
            continue;
        }
        EXPECT_EQ(static_cast<uint64_t>(f), values[i]);
    }
}
TEST(fixed_point128, ConstructorFromString) {
    const char* values[] = { "0", "12.34"};
    for (auto i = 0u; i < array_length(values); ++i) {
        fixed_point128<20> f = values[i];
        double d1 = strtod(values[i], nullptr);
        double d2 = strtod(static_cast<char*>(f), nullptr);
        EXPECT_DOUBLE_EQ(d1, d2);
    }
}
