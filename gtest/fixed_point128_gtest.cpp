#include <gtest/gtest.h>
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
TEST(FixedPoint128Construction, FromDouble) {
    double doubleValues[] = { 0.0, 1.0, 2.5, 0.001, 255.000001 };

    for (auto i = 0; i < array_length(doubleValues); ++i) {
        double dval = doubleValues[i];
        fixed_point128<10> f = dval;
        EXPECT_DOUBLE_EQ(static_cast<double>(f), dval);
    }
    
    // Expect two strings not to be equal.
    //EXPECT_STRNE("hello", "world");
    // Expect equality.
    //EXPECT_EQ(7 * 6, 42);
}
