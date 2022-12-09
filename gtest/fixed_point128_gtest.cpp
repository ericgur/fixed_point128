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


// Construct fixed_point128 from various elements.
TEST(FixedPoint128Construction, Double) {
    double values[] = { 0.0, 1.0, 2.5, 0.001, 255.000001 };
    for (auto i = 0u; i < array_length(values); ++i) {
        fixed_point128<10> f = values[i];
        EXPECT_DOUBLE_EQ(static_cast<double>(f), values[i]);
    }
}
TEST(FixedPoint128Construction, Float) {
    float values[] = { 0.0f, 1.0f, 2.5f, 0.001f, 255.000001f };
    for (auto i = 0u; i < array_length(values); ++i) {
        fixed_point128<10> f = values[i];
        EXPECT_FLOAT_EQ(static_cast<float>(f), values[i]);
    }
}
TEST(FixedPoint128Construction, Int32) {
    int32_t values[] = { 0, 10, -1543 };
    for (auto i = 0u; i < array_length(values); ++i) {
        fixed_point128<20> f = values[i];
        EXPECT_EQ(static_cast<int32_t>(f), values[i]);
    }
}
TEST(FixedPoint128Construction, UnsignedInt32) {
    uint32_t values[] = { 0u, 10u, 1543u };
    for (auto i = 0u; i < array_length(values); ++i) {
        fixed_point128<20> f = values[i];
        EXPECT_EQ(static_cast<uint32_t>(f), values[i]);
    }
}
TEST(FixedPoint128Construction, Int64) {
    uint64_t values[] = { 0ll, 10ll, 1543ll };
    for (auto i = 0u; i < array_length(values); ++i) {
        fixed_point128<20> f = values[i];
        EXPECT_EQ(static_cast<int64_t>(f), values[i]);
    }
}

TEST(FixedPoint128Construction, UnsignedInt64) {
    uint64_t values[] = { 0ull, 10ull, 1543ull };
    for (auto i = 0u; i < array_length(values); ++i) {
        fixed_point128<20> f = values[i];
        EXPECT_EQ(static_cast<uint64_t>(f), values[i]);
    }
}

TEST(FixedPoint128Construction, String) {
    const char* values[] = { "0", "12.34"};
    for (auto i = 0u; i < array_length(values); ++i) {
        fixed_point128<20> f = values[i];
        double d1 = strtod(values[i], nullptr);
        double d2 = strtod(static_cast<char*>(f), nullptr);
        EXPECT_DOUBLE_EQ(d1, d2);
    }
}
