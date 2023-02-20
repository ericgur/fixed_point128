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
TEST(float128, DefaultConstructor) {
    float128 f;
    EXPECT_DOUBLE_EQ(static_cast<double>(f), 0.0);
}
TEST(float128, ConstructorFromDouble) {
    srand(RANDOM_SEED); 
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random(-1023, 1022);
        float128 f = value;
        double f_value = static_cast<double>(f);
        EXPECT_DOUBLE_EQ(f_value, value) << "value=" << value;
    }
}
TEST(float128, ConstructorFromFloat) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        float value = (float)get_double_random(-1023, 1022);
        float128 f = value;
        float f_value = static_cast<float>(f);
        EXPECT_FLOAT_EQ(f_value, value) << "value=" << value;
    }
}
TEST(float128, ConstructorFromInt32) {
    srand(RANDOM_SEED); 
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        int32_t value = get_int32_random();
        float128 f = value;
        EXPECT_EQ(static_cast<int32_t>(f), value);
    }
}
TEST(float128, ConstructorFromUnsignedInt32) {
    srand(RANDOM_SEED); 
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        uint32_t value = get_uint32_random();
        float128 f = value;
        EXPECT_EQ(static_cast<uint32_t>(f), value);
    }
}
TEST(float128, ConstructorFromInt64) {
    srand(RANDOM_SEED); 
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        int64_t value = get_int64_random();
        float128 f = value;
        EXPECT_EQ(static_cast<int64_t>(f), value);
    }
}
TEST(float128, ConstructorFromUnsignedInt64) {
    srand(RANDOM_SEED); 
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        uint64_t value = get_uint64_random();
        float128 f = value;
        EXPECT_EQ(static_cast<uint64_t>(f), value);
    }
}
TEST(float128, ConstructorFromString) {
    char str[128];
    char str2[128];
    srand(RANDOM_SEED);
    int j = 0;
    constexpr auto end_char_to_check = 9;
    constexpr auto max_allowed_error = 1;
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        // construct random number string
        str[0] = get_digit_random();
        str[1] = '.';
        for (j = 2; j < 36; ++j) {
            str[j] = get_digit_random();
        }
        // fraction length is based on the value of the first digit
        switch (str[0])
        {
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
        size_t len = min(strlen(res), strlen(str));
        // clip both string based on the shortest one.
        str[len] = '\0';
        res[len] = '\0';
        int64_t orig_last_digits = strtoll(&str[len - end_char_to_check], nullptr, 10);
        int64_t res_last_digits = strtoll(&res[len - end_char_to_check], nullptr, 10);
        int64_t err = abs(orig_last_digits - res_last_digits);

        EXPECT_LE(err, max_allowed_error)<< "error: " << err << ", last digits source: " << &str[len - end_char_to_check] << "last digits result: " << &res[len - end_char_to_check];
    }
}
TEST(float128, CopyConstructor) {
    srand(RANDOM_SEED); 
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random();
        float128 f1 = value;
        float128 f2 = f1;
        EXPECT_DOUBLE_EQ(static_cast<double>(f1), static_cast<double>(f2)) << "value=" << value;;
    }
}
TEST(float128, MoveConstructor) {
    srand(RANDOM_SEED); 
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random();
        float128 f1 = value;
        float128 f2(std::move(f1));
        EXPECT_DOUBLE_EQ(static_cast<double>(f1), static_cast<double>(f2)) << "value=" << value;;
    }
}
TEST(float128, AssignmentOperator) {
    srand(RANDOM_SEED); 
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random();
        float128 f1 = value;
        float128 f2;
        f2 = f1;
        EXPECT_DOUBLE_EQ(static_cast<double>(f1), static_cast<double>(f2)) << "value=" << value;;
    }
}
TEST(float128, MoveAssignmentOperator) {
    srand(RANDOM_SEED); 
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random(-1023, 1022);
        float128 f1 = value;
        float128 f2;
        f2 = std::move(f1);
        EXPECT_DOUBLE_EQ(static_cast<double>(f1), static_cast<double>(f2)) << "value=" << value;;
    }
}
TEST(float128, AddSameSign) {
    srand(RANDOM_SEED); 
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = fabs(get_double_random());
        double value2 = value1 * 2.5;
        double res = value1 + value2;
        float128 f1 = value1;
        float128 f2 = value2;
        float128 f3 = f1 + f2;
        double float128_res = static_cast<double>(f3);
        // note that float128 is more precise than double, this can lead to issues when param1 and param2 are far apart.
        //if (is_similar_double(float128_res, res))
        //    continue;

        EXPECT_DOUBLE_EQ(float128_res, res) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(float128, AddDifferentSign) {
    srand(RANDOM_SEED); 
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = fabs(get_double_random());
        double value2 = value1 * -2.5;
        double res = value1 + value2;
        float128 f1 = value1;
        float128 f2 = value2;
        float128 f3 = f1 + f2;
        double float128_res = static_cast<double>(f3);
        // note that float128 is more precise than double, this can lead to issues when param1 and param2 are far apart.
        //if (is_similar_double(float128_res, res))
        //    continue;

        EXPECT_DOUBLE_EQ(float128_res, res) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(float128, AddDouble) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = get_double_random();
        double value2 = get_double_random();
        double res = value1 + value2;
        float128 f1 = value1;
        float128 f3 = f1 + value2;
        double float128_res = static_cast<double>(f3);
        // note that float128 is more precise than double, this can lead to issues when param1 and param2 are far apart.
        //if (is_similar_double(float128_res, res))
        //    continue;

        EXPECT_DOUBLE_EQ(float128_res, res) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(float128, AddFloat) {
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
TEST(float128, AddInt32) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        auto value1 = abs(get_int32_random());
        auto value2 = abs(get_int32_random());
        auto res = value1 + value2;
        if (res < 0) continue; //wrap around won't happen in uint128_t

        float128 f1 = value1;
        float128 f3 = f1 + value2;
        auto float128_res = static_cast<int64_t>(f3);
        EXPECT_EQ(float128_res, res) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(float128, AddUnsignedInt32) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        auto value1 = get_uint32_random();
        auto value2 = get_uint32_random();
        auto res = value1 + value2;
        if (res < value1 || res < value2) continue; //wrap around

        float128 f1 = value1;
        float128 f3 = f1 + value2;
        auto float128_res = static_cast<uint64_t>(f3);
        EXPECT_EQ(float128_res, res) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(float128, AddInt64) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        auto value1 = abs(get_int64_random());
        auto value2 = abs(get_int64_random());
        auto res = value1 + value2;
        if (res < value1 || res < value2) continue; //wrap around
        float128 f1 = value1;
        float128 f3 = f1 + value2;
        auto float128_res = static_cast<int64_t>(f3);
        EXPECT_EQ(float128_res, res) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(float128, AddUnsignedInt64) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        auto value1 = get_uint64_random();
        auto value2 = get_uint64_random();
        auto res = value1 + value2;
        if (res < value1 || res < value2) continue; //wrap around
        float128 f1 = value1;
        float128 f3 = f1 + value2;
        auto float128_res = static_cast<uint64_t>(f3);
        EXPECT_EQ(float128_res, res) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(float128, SubtractSameSign) {
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
TEST(float128, SubtractDifferentSign) {
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
TEST(float128, SubtractInt32) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        auto value1 = abs(get_int32_random());
        auto value2 = abs(get_int32_random());
        if (value2 > value1) value2 = value1 / 3;

        auto res = value1 - value2;

        float128 f1 = value1;
        float128 f3 = f1 - value2;
        auto float128_res = static_cast<int32_t>(f3);
        EXPECT_EQ(float128_res, res) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(float128, SubtractUnsignedInt32) {
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
TEST(float128, SubtractInt64) {
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
TEST(float128, SubtractUnsignedInt64) {
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
TEST(float128, MultiplyByFloat128 ) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = get_double_random();
        double value2 = get_double_random();
        double res = value1 * value2;
        float128 f1 = value1;
        float128 f2 = value2;
        float128 f3 = f1 * f2;

        double float128_res = static_cast<double>(f3);
        // note that float128 is more precise than double, this can lead to issues when param1 and param2 are far apart.
        //if (is_similar_double(float128_res, res))
        //    continue;

        EXPECT_DOUBLE_EQ(float128_res, res) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(float128, MultiplyByExponentOf2) {
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
TEST(float128, MultiplyByDouble) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = get_double_random();
        double value2 = get_double_random();
        value1 = 0.80270613357871456, value2 = 0.70710678118654757;
        double res = value1 * value2;
        float128 f1 = value1;
        float128 f3 = f1 * value2;
        double float128_res = static_cast<double>(f3);
        // note that float128 is more precise than double, this can lead to issues when param1 and param2 are far apart.
        //if (is_similar_double(float128_res, res))
        //    continue;

        EXPECT_DOUBLE_EQ(float128_res, res) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(float128, MultiplyByFloat) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        float value1 = (float)get_double_random();
        float value2 = (float)get_double_random();
        float res = value1 * value2;
        float128 f1 = value1;
        float128 f3 = f1 * value2;
        float float128_res = static_cast<float>(f3);
        // note that float128 is more precise than double, this can lead to issues when param1 and param2 are far apart.
        //if (is_similar_double(float128_res, res))
        //    continue;

        EXPECT_FLOAT_EQ(float128_res, res) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(float128, MultiplyByInt32) {
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
TEST(float128, MultiplyByUnsignedInt32) {
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
TEST(float128, MultiplyByInt64) {
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
TEST(float128, MultiplyByUnsignedInt64) {
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
TEST(float128, DivideByFloat128) {
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
TEST(float128, DivideByExponentOf2) {
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
TEST(float128, DivideByDouble) {
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
TEST(float128, DivideByFloat) {
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
TEST(float128, DivideByInt32) {
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
TEST(float128, DivideByUnsignedInt32) {
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

TEST(float128, DivideByInt64) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = get_double_random(-1, 39);
        int64_t value2 = get_int32_random(); // on purpose
        double res = value1 / value2;
        if (value2 == 0)
            continue;
        float128 f1 = value1;
        float128 f3 = f1 / value2;
        double float128_res = static_cast<double>(f3);
        EXPECT_DOUBLE_EQ(float128_res, res) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(float128, DivideByUnsignedInt64) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = get_double_random(-1, 39);
        uint64_t value2 = get_uint32_random(); // on purpose
        double res = value1 / value2;
        if (value2 == 0)
            continue;
        float128 f1 = value1;
        float128 f3 = f1 / value2;
        double float128_res = static_cast<double>(f3);
        EXPECT_DOUBLE_EQ(float128_res, res) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(float128, CompareFloat128) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = get_double_random(-1, 39);
        double value2 = get_double_random(-1, 39);
        bool res = value1 > value2;
        float128 f1 = value1;
        float128 f2 = value2;

        bool float128_res = f1 > f2;
        EXPECT_TRUE(float128_res == res) << "operator>: " << "value1=" << value1 << ", value2=" << value2;
        
        res = value1 >= value2;
        float128_res = f1 >= f2;
        EXPECT_TRUE(float128_res == res) << "operator>=: " << "value1=" << value1 << ", value2=" << value2;

        res = value1 < value2;
        float128_res = f1 < f2;
        EXPECT_TRUE(float128_res == res) << "operator<: " << "value1=" << value1 << ", value2=" << value2;

        res = value1 <= value2;
        float128_res = f1 <= f2;
        EXPECT_TRUE(float128_res == res) << "operator<=: " << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(float128, CompareDouble) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = get_double_random();
        double value2 = get_double_random();
        bool res = value1 > value2;
        float128 f1 = value1;

        bool float128_res = f1 > value2;
        EXPECT_TRUE(float128_res == res) << "operator>: " << "value1=" << value1 << ", value2=" << value2;

        res = value1 >= value2;
        float128_res = f1 >= value2;
        EXPECT_TRUE(float128_res == res) << "operator>=: " << "value1=" << value1 << ", value2=" << value2;

        res = value1 < value2;
        float128_res = f1 < value2;
        EXPECT_TRUE(float128_res == res) << "operator<: " << "value1=" << value1 << ", value2=" << value2;

        res = value1 <= value2;
        float128_res = f1 <= value2;
        EXPECT_TRUE(float128_res == res) << "operator<=: " << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(float128, CompareFloat) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        float value1 = (float)get_double_random();
        float value2 = (float)get_double_random();
        bool res = value1 > value2;
        float128 f1 = value1;

        bool float128_res = f1 > value2;
        EXPECT_TRUE(float128_res == res) << "operator>: " << "value1=" << value1 << ", value2=" << value2;

        res = value1 >= value2;
        float128_res = f1 >= value2;
        EXPECT_TRUE(float128_res == res) << "operator>=: " << "value1=" << value1 << ", value2=" << value2;

        res = value1 < value2;
        float128_res = f1 < value2;
        EXPECT_TRUE(float128_res == res) << "operator<: " << "value1=" << value1 << ", value2=" << value2;

        res = value1 <= value2;
        float128_res = f1 <= value2;
        EXPECT_TRUE(float128_res == res) << "operator<=: " << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(float128, CompareInt32) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        int32_t value1 = get_int32_random();
        int32_t value2 = get_int32_random();
        bool res = value1 > value2;
        float128 f1 = value1;

        bool float128_res = f1 > value2;
        EXPECT_TRUE(float128_res == res) << "operator>: " << "value1=" << value1 << ", value2=" << value2;

        res = value1 >= value2;
        float128_res = f1 >= value2;
        EXPECT_TRUE(float128_res == res) << "operator>=: " << "value1=" << value1 << ", value2=" << value2;

        res = value1 < value2;
        float128_res = f1 < value2;
        EXPECT_TRUE(float128_res == res) << "operator<: " << "value1=" << value1 << ", value2=" << value2;

        res = value1 <= value2;
        float128_res = f1 <= value2;
        EXPECT_TRUE(float128_res == res) << "operator<=: " << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(float128, CompareUnsignedInt32) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        uint32_t value1 = get_uint32_random();
        uint32_t value2 = get_uint32_random();
        bool res = value1 > value2;
        float128 f1 = value1;

        bool float128_res = f1 > value2;
        EXPECT_TRUE(float128_res == res) << "operator>: " << "value1=" << value1 << ", value2=" << value2;

        res = value1 >= value2;
        float128_res = f1 >= value2;
        EXPECT_TRUE(float128_res == res) << "operator>=: " << "value1=" << value1 << ", value2=" << value2;

        res = value1 < value2;
        float128_res = f1 < value2;
        EXPECT_TRUE(float128_res == res) << "operator<: " << "value1=" << value1 << ", value2=" << value2;

        res = value1 <= value2;
        float128_res = f1 <= value2;
        EXPECT_TRUE(float128_res == res) << "operator<=: " << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(float128, CompareInt64) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        int64_t value1 = get_int64_random();
        int64_t value2 = get_int64_random();
        bool res = value1 > value2;
        float128 f1 = value1;

        bool float128_res = f1 > value2;
        EXPECT_TRUE(float128_res == res) << "operator>: " << "value1=" << value1 << ", value2=" << value2;

        res = value1 >= value2;
        float128_res = f1 >= value2;
        EXPECT_TRUE(float128_res == res) << "operator>=: " << "value1=" << value1 << ", value2=" << value2;

        res = value1 < value2;
        float128_res = f1 < value2;
        EXPECT_TRUE(float128_res == res) << "operator<: " << "value1=" << value1 << ", value2=" << value2;

        res = value1 <= value2;
        float128_res = f1 <= value2;
        EXPECT_TRUE(float128_res == res) << "operator<=: " << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(float128, CompareUnsignedInt64) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        uint64_t value1 = get_uint64_random();
        uint64_t value2 = get_uint64_random();
        bool res = value1 > value2;
        float128 f1 = value1;

        bool float128_res = f1 > value2;
        EXPECT_TRUE(float128_res == res) << "operator>: " << "value1=" << value1 << ", value2=" << value2;

        res = value1 >= value2;
        float128_res = f1 >= value2;
        EXPECT_TRUE(float128_res == res) << "operator>=: " << "value1=" << value1 << ", value2=" << value2;

        res = value1 < value2;
        float128_res = f1 < value2;
        EXPECT_TRUE(float128_res == res) << "operator<: " << "value1=" << value1 << ", value2=" << value2;

        res = value1 <= value2;
        float128_res = f1 <= value2;
        EXPECT_TRUE(float128_res == res) << "operator<=: " << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(float128, OperatorEqual) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = get_double_random();
        float128 f1 = value1;
        bool float128_res = f1 == f1;
        EXPECT_TRUE(float128_res == true) << "operator==: " << "value1=" << value1;
    }
}
TEST(float128, OperatorNotEqual) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = get_double_random();
        float128 f1 = value1;
        bool float128_res = f1 != f1;
        EXPECT_TRUE(float128_res == false) << "operator!=: " << "value1=" << value1;
    }
}
TEST(float128, TemplateOperatorNotEqual) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = get_double_random();
        float128 f1 = value1;
        bool float128_res = f1 != value1;
        EXPECT_TRUE(float128_res == false) << "operator!=<double>: " << "value1=" << value1;

        f1 = static_cast<float>(value1);
        float128_res = f1 != static_cast<float>(value1);
        EXPECT_TRUE(float128_res == false) << "operator!=<float>: " << "value1=" << value1;

        f1 = static_cast<uint64_t>(value1);
        float128_res = f1 != static_cast<uint64_t>(value1);
        EXPECT_TRUE(float128_res == false) << "operator!=<uint64_t>: " << "value1=" << value1;

        f1 = static_cast<int64_t>(value1);
        float128_res = f1 != static_cast<int64_t>(value1);
        EXPECT_TRUE(float128_res == false) << "operator!=<int64_t>: " << "value1=" << value1;

        f1 = static_cast<uint32_t>(value1);
        float128_res = f1 != static_cast<uint32_t>(value1);
        EXPECT_TRUE(float128_res == false) << "operator!=<uint32_t>: " << "value1=" << value1;

        f1 = static_cast<int32_t>(value1);
        float128_res = f1 != static_cast<int32_t>(value1);
        EXPECT_TRUE(float128_res == false) << "operator!=<int32_t>: " << "value1=" << value1;
    }
}
TEST(float128, floor) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random();
        double res = ::floor(value);
        float128 f1 = value;
        float128 float128_res = floor(f1);
        EXPECT_DOUBLE_EQ(float128_res, res) << "floor: " << "value=" << value;
    }
}
TEST(float128, ceil) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random();
        double res = ::ceil(value);
        float128 f1 = value;
        float128 float128_res = ceil(f1);
        EXPECT_DOUBLE_EQ(float128_res, res) << "ceil: " << "value=" << value;
    }
}
TEST(float128, trunc) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random();
        double res = trunc(value);
        float128 f1 = value;
        float128 float128_res = trunc(f1);
        EXPECT_DOUBLE_EQ(float128_res, res) << "trunc: " << "value=" << value;
    }
}
TEST(float128, round) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random();
        double res = round(value);
        float128 f1 = value;
        float128 float128_res = round(f1);
        EXPECT_DOUBLE_EQ(float128_res, res) << "round: " << "value=" << value;
    }
}
TEST(float128, copysign) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = get_double_random();
        double value2 = get_double_random();
        double res = copysign(value1, value2);
        float128 f1 = value1;
        float128 f2 = value2;
        float128 float128_res = copysign(f1, f2);
        EXPECT_DOUBLE_EQ(float128_res, res) << "copysign: " << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(float128, fmod_fraction) {
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
        EXPECT_DOUBLE_EQ(float128_res, res) << "fmod: " << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(float128, fmod_integer) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = get_int32_random();
        double value2 = get_int32_random();
        double res = ::fmod(value1, value2);
        float128 f1 = value1;
        float128 f2 = value2;
        float128 float128_res = fmod(f1, f2);
        EXPECT_DOUBLE_EQ(float128_res, res) << "fmod: " << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(float128, modf) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random();
        double res_int = 0;
        double res = ::modf(value, &res_int);
        float128 f1 = value;
        float128 float128_res_int;
        float128 float128_res_frac = modf(f1, &float128_res_int);
        EXPECT_DOUBLE_EQ(float128_res_frac, res) << "modf fraction: " << "value=" << value;
        EXPECT_DOUBLE_EQ(float128_res_int, res_int) << "modf integer: " << "value=" << value;
    }
}
TEST(float128, sqrt) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = fabs(get_double_random());
        double res = ::sqrt(value);
        float128 f1 = value;
        float128 float128_res = sqrt(f1);
        EXPECT_DOUBLE_EQ(float128_res, res) << "sqrt: " << "value=" << value;
    }
}
TEST(float128, reciprocal) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random();
        double res = 1.0 / value;
        float128 f1 = value;
        float128 float128_res = reciprocal(f1);
        EXPECT_DOUBLE_EQ(float128_res, res) << "reciprocal: " << "value=" << value;
    }
}
TEST(float128, hypot) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = get_double_random();
        double value2 = get_double_random();
        double res = hypot(value1, value2);
        float128 f1 = value1;
        float128 f2 = value2;
        //float128 t = f1 * f1;
        //float128 t2 = sqrt(t);

        float128 float128_res = hypot(f1, f2);
        EXPECT_DOUBLE_EQ(float128_res, res) << "hypot: " << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(float128, ilogb) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random();
        int32_t res = ::ilogb(value);
        float128 f1 = value;
        int32_t float128_res = ilogb(f1);
        EXPECT_EQ(float128_res, res) << "ilogb: " << "value=" << value;
    }
}
TEST(float128, log) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = fabs(get_double_random()); 
        double res = log(value);
        float128 f1 = value;
        float128 float128_res = log(f1);
        EXPECT_DOUBLE_EQ(float128_res, res) << "log: " << "value=" << value;
    }
}
TEST(float128, log2) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = fabs(get_double_random());
        double res = log2(value);
        float128 f1 = value;
        float128 float128_res = log2(f1);
        EXPECT_DOUBLE_EQ(float128_res, res) << "log2: " << "value=" << value;
    }
}
TEST(float128, log10) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = fabs(get_double_random());
        double res = log10(value);
        float128 f1 = value;
        float128 float128_res = log10(f1);
        EXPECT_DOUBLE_EQ(float128_res, res) << "log10: " << "value=" << value;
    }
}
TEST(float128, log1p) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = fabs(get_double_random());
        double res = log1p(value);
        float128 f1 = value;
        float128 float128_res = log1p(f1);
        EXPECT_DOUBLE_EQ(float128_res, res) << "log1p: " << "value=" << value;
    }
}
TEST(float128, logb) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = fabs(get_double_random(-40, 15)); // lower exponent results in lost bits
        double res = logb(value);
        float128 f1 = value;
        float128 float128_res = logb(f1);
        EXPECT_DOUBLE_EQ(float128_res, res) << "logb: " << "value=" << value;
    }
}
TEST(float128, exp) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random(-16, 16);
        double res = exp(value);
        float128 f1 = value;
        float128 float128_res = exp(f1);
        double float128_res_ = (double)float128_res;
        EXPECT_DOUBLE_EQ(float128_res, res) << "exp: " << "value=" << value;
    }
}
TEST(float128, exp2) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random(-16, 16); // lower exponent results in lost bits
        double res = exp2(value);
        float128 f1 = value;
        float128 float128_res = exp2(f1);
        EXPECT_DOUBLE_EQ(float128_res, res) << "exp2: " << "value=" << value;
    }
}
TEST(float128, expm1) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random(-16, 16); // lower exponent results in lost bits
        double res = expm1(value);
        float128 f1 = value;
        float128 float128_res_ = expm1(f1);
        double float128_res = (double)float128_res_;
        EXPECT_DOUBLE_EQ(float128_res, res) << "expm1: " << "value=" << value;
    }
}
TEST(float128, pow) {
    // test special values
    double values[][2] = {
        {0, 1},
        {1, INFINITY},
        {1.1, INFINITY},
        {INFINITY, 1},
        {NAN, 1},
        {1, NAN},
        {-2, 5},
        {-2, 5.5}
    };
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

        EXPECT_DOUBLE_EQ(float128_res, res) << "pow: " << " value1=" << value1 << ", value2=" << value2;
    }

    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        //double value1 = fabs(get_double_random());
        double value1 = fabs(get_double_random());
        double value2 = get_double_random(-16, 16); // lower exponent results in lost bits
        double res = pow(value1, value2);
        float128 f1 = value1;
        float128 f2 = value2;
        float128 float128_res = pow(f1, f2);
        EXPECT_DOUBLE_EQ(float128_res, res) << "pow: " << " value1=" << value1 << ", value2=" << value2;
    }
}

TEST(float128, sin) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random(-60, 3);
        double res = ::sin(value);
        float128 f1 = value;
        float128 float128_res = sin(f1);
        EXPECT_DOUBLE_EQ(float128_res, res) << "sin: " << "value=" << value;
    }
}
TEST(float128, cos) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random(-60, 3);
        double res = ::cos(value);
        float128 f1 = value;
        float128 float128_res = cos(f1);
        EXPECT_DOUBLE_EQ(float128_res, res) << "cos: " << "value=" << value;
    }
}
TEST(float128, tan) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random(-60, 3);
        double res = ::tan(value);
        float128 f1 = value;
        float128 float128_res = tan(f1);
        EXPECT_DOUBLE_EQ(float128_res, res) << "tan: " << "value=" << value;
    }
}
TEST(float128, asin) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random(-60, -1);
        double res = ::asin(value);
        float128 f1 = value;
        float128 float128_res = asin(f1);
        EXPECT_DOUBLE_EQ(float128_res, res) << "asin: " << "value=" << value;
    }
}
TEST(float128, acos) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random(-60, -1);
        double res = acos(value);
        float128 f1 = value;
        float128 float128_res = acos(f1);
        EXPECT_DOUBLE_EQ(float128_res, res) << "acos: " << "value=" << value;
    }
}
TEST(float128, atan) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random(-60, 14);
        double res = atan(value);
        float128 f1 = value;
        float128 float128_res = atan(f1);
        EXPECT_DOUBLE_EQ(float128_res, res) << "atan: " << "value=" << value;
    }
}
TEST(float128, atan2) {
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = get_double_random(-10, 14); //lower exponent results in lost bits
        double value2 = get_double_random(-10, 14); //lower exponent results in lost bits
        double res = atan2(value1, value2);
        float128 f1 = value1;
        float128 f2 = value2;
        float128 float128_res = atan2(f1, f2);
        EXPECT_DOUBLE_EQ(float128_res, res) << "atan2: " << " value1=" << value1 << ", value2=" << value2;
    }
}
//TEST(float128, sinh) {
//    srand(RANDOM_SEED);
//    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
//        double value = get_double_random(-60, 2);
//        double res = ::sinh(value);
//        float128 f1 = value;
//        float128 float128_res = sinh(f1);
//        EXPECT_DOUBLE_EQ(float128_res, res) << "sinh: " << "value=" << value;
//    }
//}
//TEST(float128, asinh) {
//    srand(RANDOM_SEED);
//    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
//        double value = get_double_random(-60, 2);
//        double res = ::asinh(value);
//        float128 f1 = value;
//        float128 float128_res = asinh(f1);
//        EXPECT_DOUBLE_EQ(float128_res, res) << "asinh: " << "value=" << value;
//    }
//}
//TEST(float128, cosh) {
//    srand(RANDOM_SEED);
//    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
//        double value = get_double_random(-60, 2);
//        double res = ::cosh(value);
//        float128 f1 = value;
//        float128 float128_res = cosh(f1);
//        EXPECT_DOUBLE_EQ(float128_res, res) << "cosh: " << "value=" << value;
//    }
//}
//TEST(float128, acosh) {
//    srand(RANDOM_SEED);
//    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
//        double value = 1.0 + fabs(get_double_random(-60, 2)); // must be >= 1
//        double res = ::acosh(value);
//        float128 f1 = value;
//        float128 float128_res = acosh(f1);
//        EXPECT_DOUBLE_EQ(float128_res, res) << "acosh: " << "value=" << value;
//    }
//}
//TEST(float128, tanh) {
//    srand(RANDOM_SEED);
//    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
//        double value = get_double_random(-60, 2);
//        double res = ::tanh(value);
//        float128 f1 = value;
//        float128 float128_res = tanh(f1);
//        EXPECT_DOUBLE_EQ(float128_res, res) << "tanh: " << "value=" << value;
//    }
//}
//
//TEST(float128, atanh) {
//    srand(RANDOM_SEED);
//    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
//        double value = get_double_random(-60, -1); // must be: -1 < value < 1
//        double res = ::atanh(value);
//        float128 f1 = value;
//        float128 float128_res = atanh(f1);
//        EXPECT_DOUBLE_EQ(float128_res, res) << "atanh: " << "value=" << value;
//    }
//}
