#ifndef UNREFERENCED_PARAMETER
#define UNREFERENCED_PARAMETER(P)          (P)
#endif 

#include "..\inc\fixed_point128.h"
#include "..\inc\uint128_t.h"
#include "..\inc\float128.h"
#ifndef GTEST_SHARED_H
#define GTEST_SHARED_H

/*************************************************
* Fixed point 128 tests
**************************************************/

using namespace fp128;

static constexpr int RANDOM_TEST_COUNT = 1 << 16;
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
    return (rand() & 1) ? 1 : -1;
}

// returns a random number
double static get_double_random(int32_t min_exponent = -10, int32_t max_exponent = 63)
{
    Double res;
    int expo = (get_uint32_random() % (max_exponent - min_exponent)) + min_exponent;
    res.e = (uint64_t)expo + 1023; 
    res.f = get_uint64_random();
    res.s = get_random_sign() == 1;
    return res.val;
}

// returns a positive random number
uint64_t static get_uint64_random()
{
    return (((uint64_t)rand()) << 60) + (((uint64_t)rand()) << 45) + (((uint64_t)rand()) << 30) + (((uint64_t)rand()) << 15) + (uint64_t)rand();
}

// returns a positive random number
int64_t static get_int64_random()
{
    return (int64_t)get_uint64_random() * (int64_t)get_random_sign();
}

// returns a positive random number
uint32_t static get_uint32_random()
{
    return (((uint32_t)rand()) << 30) + (((uint32_t)rand()) << 15) + (uint32_t)rand();
}

// returns a random number
int32_t static get_int32_random()
{
    return (int32_t)get_uint32_random() * get_random_sign();
}

char static get_digit_random()
{
    return (char)(rand() % 10) + '0';
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

bool static check_overflow_uint128(double value)
{
    return floor(log2(value)) > 127;
}

bool static is_similar_double(double v1, double v2)
{
    if (v1 == 0.0 || v2 == 0.0) {
        v1 += 1.0e-30;
        v2 += 1.0e-30;
    }
    double ratio = fabs(v1 / v2 - 1.0);
    return ratio < DOUBLE_REL_EPS;

}


#endif //#ifndef GTEST_SHARED_H
