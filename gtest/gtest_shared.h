#include "..\inc\fixed_point128.h"
#include "..\inc\uint128_t.h"
#ifndef GTEST_SHARED_H
#define GTEST_SHARED_H

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
double static get_double_random()
{
    Double res;
    res.e = ((uint64_t)get_uint32_random() % 73) - 10 + 1023; // exponents [-10,63]
    res.f = get_uint64_random();
    res.s = get_random_sign();
    return res.val;
}

// returns a positive random number
uint64_t static get_uint64_random()
{
    return ((uint64_t)rand() + 1) * (uint64_t)rand();
}

// returns a positive random number
int64_t static get_int64_random()
{
    return (int64_t)(rand() + 1) * (int64_t)rand() * (int64_t)get_random_sign();
}

// returns a positive random number
uint32_t static get_uint32_random()
{
    return (uint32_t)rand();
}

// returns a random number
int32_t static get_int32_random()
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
