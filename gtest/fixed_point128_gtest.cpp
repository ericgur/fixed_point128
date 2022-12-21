// remove warnings from gtest itself
#pragma warning(push)
#pragma warning(disable: 26439) 
#pragma warning(disable: 26495) 
#include <gtest/gtest.h>
#pragma warning(pop)
#include <ostream>
#include <ctime>
#include "..\inc\fixed_point128.h"
#include "..\inc\uint128_t.h"

/*************************************************
* Fixed point 128 tests
**************************************************/
using namespace fp128;

static constexpr int RANDOM_TEST_COUNT = 32768;
static constexpr int RANDOM_SEED = 0x12345678; // must have a repeatable seed for debugging
static constexpr double DOUBLE_REL_EPS = 1.0e-10;

uint64_t get_uint64_random();
int64_t get_int64_random();
uint32_t get_uint32_random();
int32_t get_int32_random();

// friend class to all containers to simplify test cases
namespace fp128 {
    class fp128_gtest
    {
        template<int32_t I>
        inline static void get_fixed_point128_members(const fixed_point128<I>& obj, uint64_t& l, uint64_t& h, uint32_t& s) {
            l = obj.low;
            h = obj.high;
            s = obj.sign;
        }
        inline static void get_uint128_t_members(const uint128_t& obj, uint64_t& l, uint64_t& h) {
            l = obj.low;
            h = obj.high;
        }
    };
}

__forceinline int32_t get_random_sign()
{
    return (rand() > (RAND_MAX >> 1)) ? -1 : 1;
}

// returns a random number
double get_double_random()
{
    Double res;
    res.e = ((uint64_t)get_uint32_random() % 73) - 10 + 1023; // exponents [-10,63]
    res.f = get_uint64_random();
    res.s = get_random_sign();
    return res.val;
}

// returns a positive random number
uint64_t get_uint64_random()
{
    return ((uint64_t)rand() + 1) * (uint64_t)rand();
}

// returns a positive random number
int64_t get_int64_random()
{
    return (int64_t)(rand() + 1) * (int64_t)rand() * (int64_t)get_random_sign();
}

// returns a positive random number
uint32_t get_uint32_random()
{
    return (uint32_t)rand();
}

// returns a random number
int32_t get_int32_random()
{
    return rand() * get_random_sign();
}

// return true on overflow
template<typename T, int I>
bool check_overflow(T value, const fixed_point128<I>& d)
{
    if constexpr (std::is_floating_point<T>::value) {
        value = fabs(value);
        if (value <= 1.0) return false;
    }
    else if constexpr (std::is_signed<T>::value) {
        value = abs(value);
        if (value <= 1) return false;
    }

    return floor(log2(value)) >= I;
}

bool check_overflow_uint128(double value)
{
    return floor(log2(value)) > 127;
}

bool is_similar_double(double v1, double v2)
{
    if (v1 == 0.0 || v2 == 0.0) {
        v1 += 1.0e-30;
        v2 += 1.0e-30;
    }
    double ratio = fabs(v1 / v2 - 1.0);
    return ratio < DOUBLE_REL_EPS;
        
}
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
        double value = get_double_random();
        //double value = -0.00048828129930906378;
        
        fixed_point128<20> f = value;
        if (check_overflow(value, f)) {
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
    srand(RANDOM_SEED); 
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        float value = (float)get_double_random();
        fixed_point128<20> f = value;
        if (check_overflow(value, f)) {
            continue;
        }
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
        EXPECT_DOUBLE_EQ(d1, d2);
    }
}
TEST(fixed_point128, CopyConstructor) {
    srand(RANDOM_SEED); 
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random();
        fixed_point128<20> f1 = value;
        fixed_point128<20> f2 = f1;
        EXPECT_DOUBLE_EQ(static_cast<double>(f1), static_cast<double>(f2));
    }
}
TEST(fixed_point128, MoveConstructor) {
    srand(RANDOM_SEED); 
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random();
        fixed_point128<20> f1 = value;
        fixed_point128<20> f2(std::move(f1));
        EXPECT_DOUBLE_EQ(static_cast<double>(f1), static_cast<double>(f2));
    }
}
TEST(fixed_point128, AssignmentOperator) {
    srand(RANDOM_SEED); 
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random();
        fixed_point128<20> f1 = value;
        fixed_point128<20> f2;
        f2 = f1;
        EXPECT_DOUBLE_EQ(static_cast<double>(f1), static_cast<double>(f2));
    }
}
TEST(fixed_point128, AssignmentOperatorOtherType) {
    srand(RANDOM_SEED); 
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random();
        fixed_point128<20> f1 = value;
        fixed_point128<22> f2;
        f2 = f1;
        EXPECT_DOUBLE_EQ(static_cast<double>(f1), static_cast<double>(f2));
    }
}
TEST(fixed_point128, MoveAssignmentOperator) {
    srand(RANDOM_SEED); 
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random();
        fixed_point128<20> f1 = value;
        fixed_point128<20> f2;
        f2 = std::move(f1);
        EXPECT_DOUBLE_EQ(static_cast<double>(f1), static_cast<double>(f2));
    }
}
TEST(fixed_point128, CopyConstructorOtherType) {
    srand(RANDOM_SEED); 
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random();
        fixed_point128<20> f1 = value;
        fixed_point128<22> f2 = f1;
        EXPECT_DOUBLE_EQ(static_cast<double>(f1), static_cast<double>(f2));
    }
}
TEST(fixed_point128, AddSameSign) {
    srand(RANDOM_SEED); 
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = fabs(get_double_random());
        double value2 = value1 * 2.5;
        //double value1 = 3.6379791744728668e-12, value2 = 9.0949479361821669e-12;
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
        //value1 = 0.0078125007889450204, value2 = -0.019531251972362551;
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
        //float value1 = -4096.00048828125, value2 = -1099511627776;
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
        auto value1 = get_int32_random();
        auto value2 = get_int32_random();
        //value1 = 23446, value2 = -6193;
        auto res = value1 + value2;
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
        auto value1 = get_int32_random();
        auto value2 = get_int32_random();
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
        //double value1 = -9.5367434602401654e-07, value2 = -3.8146973747715616e-06;
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
        //double value1 = -9.5367438958301344e-07, value2=-549755865672.46912;
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
        //float value1 = -4194304, value2 = -0.03125;
        float res = value1 / value2;
        if (value2 == 0)
            continue;
        fixed_point128<40> f1 = value1;
        fixed_point128<40> f3 = f1 / value2;
        if (check_overflow(value1, f1) || check_overflow(value2, f1) || check_overflow(res, f3)) {
            continue;
        }
        EXPECT_FLOAT_EQ(static_cast<float>(f3), res) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(fixed_point128, DivideByInt32) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = get_double_random();
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
        double value1 = get_double_random();
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
        double value1 = get_double_random();
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
        double value1 = get_double_random();
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
        double value1 = get_double_random();
        double value2 = get_double_random();
        //double value1 = 48034.270022883298, value2 = 168.09205560447856;
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
        double value1 = get_double_random();
        double value2 = get_double_random();
        //value1=-15943.387313816127, value2=-53860.693251533739;
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
        float value1 = (float)get_double_random();
        float value2 = (float)get_double_random();
        //float value1=-274877917612.19153, value2=-0.062500000199658484;
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
        double value1 = get_double_random();
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
        double value1 = get_double_random();
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
        double value1 = get_double_random();
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
        double value1 = get_double_random();
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
        double value1 = get_double_random();
        double value2 = get_double_random();
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
        //double value1 = -549755845143.65332;

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
        //double value1 = 284; uint32_t shift = 125;
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
        //double value1 = 284; uint32_t shift = 125;
        fixed_point128<40> f1 = value1;
        double expo = pow(2, static_cast<int32_t>(shift));
        double res = value1 * expo; // double doesn't lose any bits with this operation!
        fixed_point128 f3 = f1 << shift;

        if (check_overflow(value1, f1) || check_overflow(res, f3)) {
            continue;
        }
        EXPECT_DOUBLE_EQ(static_cast<double>(f3), res) << "double value1=" << value1 << "; uint32_t shift=" << shift << ";";
    }
}
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
    const char* values[] = { "0xDEADBEAF", "1234", "123456789012345678"};
    for (auto i = 0u; i < array_length(values); ++i) {
        uint128_t f = values[i];
        double d1 = strtod(values[i], nullptr);
        double d2 = strtod(static_cast<char*>(f), nullptr);
        EXPECT_DOUBLE_EQ(d1, d2);
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
