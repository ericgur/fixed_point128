/***********************************************************************************
    MIT License

    Copyright (c) 2025 Eric Gur (ericgur@iname.com)

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.
************************************************************************************/

// fixed_point128.cpp : test executable for the fixed_point128 class template
//

//#define FP128_DISABLE_INLINE TRUE

#include <windows.h>
#include <profileapi.h>
#include <cstdio>
#include <cassert>
#include <cmath>
#include "../inc/fixed_point128.h" 
#include "../inc/uint128_t.h" 
#include "../inc/int128_t.h" 
#include "../inc/float128.h" 

#pragma warning(disable: 26493) // Don't use C-style casts.
#pragma warning(disable: 26467) 
#pragma warning(disable: 26485) 
#pragma warning(disable: 26440) 
#pragma warning(disable: 26482) 
#pragma warning(disable: 26446) 
#pragma warning(disable: 26496) 

void bench();

using namespace fp128;

void test_conversion()
{
    printf("\nTest Conversion\n");
    float128 _pi = float128::pi();
    float128 e = float128::e();
    float128 sqrt2 = float128::sqrt_2();
    //double e_0p5 = exp(0.5);
    float128 e2 = exp(float128(0.5));

    float128 e3 = e2 * e2;
    float128 f128_zero;
    double d1 = 0.51;
    float128 df128_v1 = d1;
    //double f128_v1 = df128_v1;
    d1 = pow(2, -1025)/3.0;
    float128 df128_v2 = d1;
    double d2 = df128_v2;
    df128_v1 = "10.1234567890";
    char* str = (char*)df128_v1;
    d2 = 1e100;
    //int32_t i = (int32_t)d2;
    df128_v1 = "3.14159265358979323846264338327950288419716939937510"; // 50 first digits of pi
    df128_v1 *= float128::exp10(96);
    char buff[128];
    df128_v1.to_e_format(buff, 40);
    printf(buff);
    printf("\n");

    str = (char*)df128_v1;
    df128_v1 = float128::pi();
    str = (char*)df128_v1;
    df128_v1 = "2.71828182845904523536028747135266249775724709369"; // 50 first digits of e
    str = (char*)df128_v1;
    df128_v1 >>= 1;
    df128_v1 = 1;
    df128_v1 /= 3;

    str = (char*)df128_v1;
    df128_v1 = "31415.9265358979323846264338327950288419716939937510E-4"; // 50 first digits of pi
    float128 factor = "0.70710678118654752440084436210484903928483593768847403658833981";
    df128_v1 = "10";
    df128_v1 = "100";
    df128_v1 = "1000000000000000000000000000000000000000000";
    df128_v1 = "0x1";
    df128_v1 = "0x6";
    df128_v1 = "0xF";
    df128_v1 = "0xFF";
    df128_v1 = "0xDEAD";
    df128_v1 = "1234567890123456789012345678901234";
    df128_v1 = "0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF";
    df128_v1 = "0xF123456789ABCDEFFEDCBA9876543210";
    df128_v1 = df128_v1.inf();
    df128_v1 = -df128_v1.inf();
    df128_v1 = df128_v1.nan();
    df128_v1 = -3.1415926535897932384626433832795;
    df128_v2 = 9.1415926535897932384626433832795;
    df128_v2 = -9.123E+0019;
    df128_v1 += df128_v2;
    d1 = HUGE_VAL;
    d1 = -d1;
    d1 = nan("");
    df128_v1 = -3.1415926535897932384626433832795;
    float128 float128_res_int;
    float128 float128_res_frac = modf(df128_v1, &float128_res_int);

    uint128_t i1 = 1;
    uint128_t i2 = UINT64_MAX;
    uint128_t i3 = "0xDEADBEAFDEADBEAF";
    uint128_t i4 = "0xF123456789ABCDEFFEDCBA9876543210";
    d1 = 1e38;
    uint128_t i5a = d1;  // highest supported base 10 exponent
    uint128_t i5b = 1.0; // lowest supported
    assert(d1 == (double)i5a);
    d1 = pow(2.1, 80);
    uint128_t i6 = d1;
    // assert(fabs((d1 / (double)i6) - 1.0) < 0.00001);
    assert(d1 == (double)i6);
    uint128_t i7 = -2;
    assert((int32_t)i7 == -2);
    //uint64_t i64 = (uint64_t)d1;
    printf("uint128_t: 1=%llu, UINT64_MAX=0x%llX, 0xDEADBEAFDEADBEAF=0x%llX\n", (uint64_t)i1, (uint64_t)i2, (uint64_t)i3);
    printf("uint128_t: 0xF123456789ABCDEFFEDCBA9876543210=%s\n", (char*)i4);
    printf("uint128_t: 0xF123456789ABCDEFFEDCBA9876543210=%s\n", i4.hex());
    double a1 = 5.5;
    double a2 = 1.0 / (1ull << 22);
    double a3 = 3.000003;
    fixed_point128<4> f1 = a1;
    fixed_point128<4> f2 = a2;
    fixed_point128<32> f3 = a3;
    double r1 = f1;
    double r2 = f2;
    double r3 = f3;

    printf("a1: %0.15lf --> %0.15lf\n", a1, r1);
    printf("a2: %0.15lf --> %0.15lf\n", a2, r2);
    printf("a3: %0.15lf --> %0.15lf\n", a3, r3);
    f3 = fabs(-f3);
    printf("Pi (double): %0.40lf\n", 3.14159265358979323846264338327950288419716939937510);
    printf("f1 = \"3.14159265358979323846264338327950288419716939937510\"\n"),
        f1 = "3.14159265358979323846264338327950288419716939937510"; // 50 first digits of pi

    printf("f1: %s\n", (char*)f1);
    f1 = f1.pi();
    printf("fixed_point128<4>.pi(): %s\n", (char*)f1);
    f1 = f1.e();
    printf("fixed_point128<4>.e(): %s\n", (char*)f1);
    f1 = f1.one();
    printf("fixed_point128<4>.one(): %s\n", (char*)f1);
    f1 = f1.pi();
    printf("f3 = f1\n");
    f3 = f1;
    printf("f1: %s\n", (char*)f1);
    printf("f3: %s\n", (char*)f3);
    printf("f1 = f3\n");
    f1 = f3;
    printf("f1: %s\n", (char*)f1);
    printf("f3: %s\n", (char*)f3);
    printf("\n");
}

void test_addition()
{
    printf("\nTest Addition\n");
    
    srand(0xDEADBEEF);
    for (auto i = 0ull; i < 32768; ++i) {
        uint64_t val1 = i * rand();
        uint64_t val2 = i * rand();
        uint128_t i128a = val1;
        uint128_t i128b = val2;
        assert((uint64_t)i128a == val1);
        assert((uint64_t)i128b == val2);
        assert((uint64_t)(i128a + i128b) == (val1 + val2));
        assert((uint64_t)(i128a - i128b) == (val1 - val2));
        assert((uint64_t)(i128a + val2) == (val1 + val2));
        assert((uint64_t)(i128a - val2) == (val1 - val2));
    }

    printf("uint128_t addition and subtraction passed!\n");
    
    fixed_point128<4> f1 = 5.123456789012345;
    fixed_point128<4> f2 =  7.0 - 5.123456789012345;
    printf("f1: %0.15lf\n", (double)f1);
    printf("f2: %0.15lf\n", (double)f2);
    ++f1;
    printf("++f1 --> %0.15lf\n", (double)f1);
    f1++;
    printf("f1++ --> %0.15lf\n", (double)f1);
    fixed_point128<4> f3 = f1 + f2;
    printf("f1+f2 --> %0.15lf\n", (double)f3);
    f3 = f2 + f1;
    printf("f2+f1 --> %0.15lf\n", (double)f3);
    f1 = -f1;
    printf("f1 = -f1 --> %0.15lf\n", (double)f1);
    f3 = f1 + f2;
    printf("f1+f2 --> %0.15lf\n", (double)f3);
    f3 = f1 - f2;
    printf("f1-f2 --> %0.15lf\n", (double)f3);
    f1 = 2.25; f2 = 4.54;
    printf("f1: %0.15lf\n", (double)f1);
    printf("f2: %0.15lf\n", (double)f2);
    f3 = f1 - f2;
    printf("f1-f2 --> %0.15lf\n", (double)f3);
}

void test_shift()
{
    printf("\nTest Shift\n");

    srand(0xDEADBEEF);
    for (auto i = 0ull; i < 32768; ++i) {
        int32_t shift = 1 + (rand() >> 10);
        uint64_t val = i * rand();
        uint128_t i128 = val;
        uint128_t res = i128 >> shift;
        assert((uint64_t)(res) == (val >> shift));
        assert(((i128 << shift) >> shift) == i128);
        assert((uint64_t)(i128 << 1) == val * 2);
    }
    printf("uint128_t bit shifts passed!\n");

    fixed_point128<4> f1 = 1.0;
    fixed_point128<4> f2 = 5.123456789012345;
    printf("f1: %0.15lf\n", (double)f1);
    f1 >>= 10;
    printf("f1 >>= 10 --> %0.15lf\n", (double)f1);
    f1 <<= 11;
    printf("f1 <<= 11 --> %0.15lf\n", (double)f1);
    printf("f2: %0.15lf\n", (double)f2);
    f2 >>= 65;
    printf("f2 >>= 65 --> %0.15lf\n", (double)f2);
    f2 <<= 66;
    printf("f2 <<= 66 --> %0.15lf\n", (double)f2);
}

void test_multiplication()
{
    printf("\nTest Multiplication\n");

    srand(0xDEADBEEF);
    for (auto i = 0ull; i < 32768; ++i) {
        uint64_t val1 = 1ull + rand();
        uint64_t val2 = 1ull + rand();
        uint128_t i128a = val1;
        uint128_t i128b = val2;
        assert((uint64_t)i128a == val1);
        assert((uint64_t)i128b == val2);
        assert((uint64_t)(i128a * i128b) == (val1 * val2));
        assert((i128a * i128b) == (i128b * i128a));
        assert(((i128a * i128b) << 51) == i128b * i128a * uint128_t(1ull << 51u));
    }

    printf("uint128_t multiplication passed!\n");

    fixed_point128<16> f1 = 1.0;
    fixed_point128<16> f2 = 1.4142135623730951; // == sqrt of 2
    fixed_point128<16> f1_sq = f1 * f1;
    fixed_point128<16> f2_sq = f2 * f2;
    printf("f1: %0.15lf\n", (double)f1);
    printf("f1 * f1: %0.15lf\n", (double)f1_sq);
    printf("f2: %0.15lf\n", (double)f2);
    printf("f2 * f2: %0.15lf\n", (double)f2_sq);
    f1 = 0.00001;
    printf("f1: %0.15lf\n", (double)f1);
    f1 *= -100;
    printf("f1 *= -100: %0.15lf\n", (double)f1);
    f1 = f1 * 100;
}

void test_division()
{
    printf("\nTest Division\n");

    srand(0xDEADBEEF);
    for (auto i = 0ull; i < 32768; ++i) {
        uint64_t val1 = i * rand();
        uint64_t val2 = 1ull + rand();
        uint128_t i128a = val1;
        uint128_t i128b = val2;
        assert((uint64_t)i128a == val1);
        assert((uint64_t)i128b == val2);
        assert((uint64_t)(i128a / i128b) == (val1 / val2));
        assert((i128a * i128b * i128b / i128b / i128b) == i128a);
    }

    printf("uint128_t division passed!\n");

    fixed_point128<20> f1 = 1.0 / 3.0;
    printf("f1: %0.15lf\n", (double)f1);
    fixed_point128<20> f3 = f1 / 64.0;
    printf("f3 = f1 / 64.0: %0.15lf\n", (double)f3);
    fixed_point128<20> f2 = f1 / 2.0;
    printf("f2 = f1 / 2.0: %0.15lf\n", (double)f2);
    f2 = f1 / 32.0;
    printf("f2 = f1 / 32.0: %0.15lf\n", (double)f2);
    f2 = f1 / 3.0;
    printf("f2 = f1 / 3.0: %0.15lf\n", (double)f2);

    f1 = 0.01 / 3.0;
    printf("f1 = 0.01 / 3.0: %0.15lf\n", (double)f1);
    f2 = 1.0;
    printf("f2: %0.15lf\n", (double)f2);
    f2 /= f1;
    printf("f2 /= f1: %0.15lf\n", (double)f2);

    f1 = 8;
    printf("f1: %0.15lf\n", (double)f1);
    f1 %= 5;
    printf("f1 %%= 5: %0.15lf\n", (double)f1);
    f1 = 8;
    printf("f1: %0.15lf\n", (double)f1);
    f1 %= -5;
    printf("f1 %%= -5: %0.15lf\n", (double)f1);
    f1 = -8;
    printf("f1: %0.15lf\n", (double)f1);
    f1 %= 5;
    printf("f1 %%= 5: %0.15lf\n", (double)f1);
    f1 = 11;
    printf("f1: %0.15lf\n", (double)f1);
    f1 %= 0.7;
    printf("f1 %%= 0.7: %0.15lf\n", (double)f1);
    f1 = 11; 
    f1 %= -0.7;
    printf("f1 %%= -0.7: %0.15lf\n", (double)f1);
}


void test_precision()
{
    printf("\nTest Precision\n");
    fixed_point128<10> f = 1;
    for (int i = 1; i <= 118; ++i)         {
        f /= 2;
        printf("exp = -%d, f: %0.40lf\n", i, (double)f);
    }
}


void test_string()
{
    char digit2char[16] = { '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'A', 'B', 'C', 'D', 'E', 'F' };
    char str[64]{};
    srand(0xDEADBEEF);
    for (auto i = 0ull; i < 32768; ++i) {
        int j;
        for (j = 0; j < 35; ++j) {
            str[j] = digit2char[rand() % 10];
        }
        if (str[0] == '0')
            str[0] = '1';
        str[j] = '\0';
        uint128_t i128 = str;
        assert(0 == strcmp((char*)i128, str));
    }

    printf("uint128_t string operations passed!\n");

    printf("\nTest string\n");
    double val = 1.0 / (3.0 * (1 << 21));
    fixed_point128<5> f1 = val;

    printf("double val:\n%0.38lf\n", val);
    printf("f1:\n%s\n", (char*)f1);
    val = -5.12345678901234567890;
    f1 = val;
    printf("double val:\n%0.38lf\n", val);
    printf("f1:\n%s\n", (char*)f1);
}


void test_functions()
{
    printf("\nTest functions\n");
    uint128_t i128a = 1ull;
    for (auto i = 1; i < 127; ++i) {
        uint128_t i128b = i128a << i;
        assert(log2(i128b) == (uint32_t)i);
    }
    printf("uint128_t log2 passed!\n");

    i128a = 1ull;
    for (auto i = 0; i < 39; ++i) { // higher than 10^38 will overflow
        assert(log10(i128a) == (uint32_t)i);
        i128a *= 10ull;
    }
    printf("uint128_t log10 passed!\n");

    srand(0xDEADBEEF);
    for (auto i = 0ull; i < 256; ++i) {
        uint128_t i1 = (1ull + rand()) * (1ull + rand());
        uint128_t i2 = i1 * i1;
        assert(i1 == sqrt(i2));
    }
    printf("uint128_t sqrt passed!\n");


    const double values[] = {3.45, -7.5, 0.27};
    std::string s;
    fixed_point128<10> f1;
    int len = sizeof(values) / sizeof(double);
    for (int i = 0; i < len; ++i) {
        f1 = values[i];
        printf("f1: %0.15lf\n", (double)f1);
        printf("floor(f1): %0.15lf\n", (double)floor(f1));
        printf("ceil(f1): %0.15lf\n", (double)ceil(f1));
        printf("fabs(f1): %0.15lf\n", (double)fabs(f1));
        printf("reciprocal(double): %0.15lf\n", 1.0 / values[i]);
        auto y = reciprocal(f1);
        printf("reciprocal(f1): %0.15lf\n", (double)y);
        fixed_point128<10> int_part;
        fixed_point128<10> f2 = modf(f1, &int_part);
        printf("modf(f1): int: %lf fraction %0.15lf\n", (double)int_part, (double)f2);
       
        f2 = 0.25;
        printf("f2: %s\n", (char*)f2);
        printf("fmod(f1, f2): %s\n", (char*)fmod(f1, f2));
        printf("sqrt(f1): %s\n", (char*)sqrt(f1));
        printf("\n");
    }

    fixed_point128<4> fvalues[] = {
        0,
        fixed_point128<4>::half_pi(),
        -fixed_point128<4>::half_pi(),
        fixed_point128<4>::half_pi() >> 1,
        (-fixed_point128<4>::half_pi()) >> 1,
        fixed_point128<4>::pi() * 1.5,
        (fixed_point128<4>::pi()) * -1.5
    };

    fixed_point128<16> f4;
    len = sizeof(fvalues) / sizeof(fvalues[0]);
    for (int i = 0; i < len; ++i) {
        f4 = fvalues[i];
        s = f4;
        printf("f4 (value): %s\n", s.c_str());
        fixed_point128<16> res = sin(f4);
        printf("sin(f4): %s\n", (char*)res);
        s = asin(res);
        printf("asin(%s): %s\n", (char*)res, s.c_str());

        res = cos(f4);
        printf("cos(f4): %s\n", (char*)res);
        s = acos(res);
        printf("acos(%s): %s\n", (char*)res, s.c_str());
    }
    
    fixed_point128<16> fvalues2[] = {
        fixed_point128<16>::one(),
        fixed_point128<16>(2),
        fixed_point128<16>::one() / fixed_point128<16>(10), // 0.1 with high precision
        fixed_point128<16>::e()
    };
    len = sizeof(fvalues2) / sizeof(fvalues2[0]);
    for (int i = 0; i < len; ++i) {
        fixed_point128<16> f5 = fvalues2[i];
        s = f5;
        printf("f5: %s\n", s.c_str());
        s = exp(f5);
        printf("exp(f5): %s\n", s.c_str());
    }

    for (int i = 0; i < len; ++i) {
        fixed_point128<16> f5 = fvalues2[i];
        s = f5;
        printf("\n");
        printf("f5: %s\n", s.c_str());
        s = log2(f5);
        printf("log2(f5): %s\n", s.c_str());
        s = log(f5);
        printf("log(f5): %s\n", s.c_str());
        s = log10(f5);
        printf("log10(f5): %s\n", s.c_str());
    }

    fixed_point128<16> f6;
    f6 = asin(fixed_point128<16>::one());
    f6 = asin(-fixed_point128<16>::one());
    f6 = asin(fixed_point128<16>(0));

    //f6 = fixed_point128<16>::pi() / 5;
    //fixed_point128<16> f6_sin, f6_cos;
    //_sincos_cordic(f6, f6_sin, f6_cos, true);
    //_sincos_cordic(-f6, f6_sin, f6_cos, true);
    printf("\n");
}

void test_comparison()
{
    printf("\nTest comparison\n");

    fixed_point128<10> t = fixed_point128<10>::one();
    assert(1 == t);
    assert(t == 1);
    assert(2 > t);
    assert(!(t > 2));
    assert(t >= 1.0);
    assert(1.0 >= t);
    assert(t > 0.0);
    assert(!(0.0 > t));

    double values[][2] = {
            {2, 1},
            {-2, 1},
            {1.0000000001, 1.0000000002},
            {-1.0000000001, 1.0000000002},
            {-1.00000000000001, 1.00000000000002},
            { 1.00000000000002, -1.00000000000001},
            { 1.0000000001, 1.0000000001 }
    };

    constexpr int len = sizeof(values) / sizeof(double) / 2;
    for (int i = 0; i < len; ++i) {
        double d0 = values[i][0];
        double d1 = values[i][1];
        fixed_point128<16> f1 = d0;
        fixed_point128<16> f2 = d1;
        printf("f1: %0.15lf\n", (double)f1);
        printf("f2: %0.15lf\n", (double)f2);
        if (d0 < d1 != f1 < f2) {
            FP128_ASSERT(d0 < d1 == f1 < f2);
            printf("Bug: f1 < f2: %s\n", (f1 < f2) ? "True" : "False");
        }
        if (d0 <= d1 != f1 <= f2) {
            FP128_ASSERT(d0 <= d1 == f1 <= f2);
            printf("Bug: f1 <= f2: %s\n", (f1 <= f2) ? "True" : "False");
        }
        if (d0 > d1 != f1 > f2) {
            FP128_ASSERT(d0 > d1 == f1 > f2);
            printf("Bug: f1 > f2: %s\n", (f1 > f2) ? "True" : "False");
        }
        if (d0 >= d1 != f1 >= f2) {
            FP128_ASSERT(d0 >= d1 == f1 >= f2);
            printf("Bug: f1 >= f2: %s\n", (f1 >= f2) ? "True" : "False");
        }
        if ((d0 == d1) != (f1 == f2)) {
            FP128_ASSERT((d0 == d1) == (f1 == f2));
            printf("Bug: f1 == f2: %s\n", (f1 == f2) ? "True" : "False");
        }
        if ((d0 != d1) != (f1 != f2)) {
            FP128_ASSERT((d0 != d1) == (f1 != f2));
            printf("Bug: f1 != f2: %s\n", (f1 != f2) ? "True" : "False");
        }
        printf("\n");
    }
}

void test_literals()
{
    printf("\nTest Literals\n");
    printf("float128:\n");

    float128 pi1 = 3.14159265358979323846264338327950288419716939937510_f128;
    float128 pi2 = float128::pi();
    auto diff = pi2 - pi1;
    printf("True Pi vs. literal Pi: %s\n", (char*)diff);

    printf("uint128_t:\n");
    uint128_t u1 = 0x1234567890abcdef1234567890abcdef_uint128;
    uint128_t u2 = "0x1234567890abcdef1234567890abcdef";
    printf("uint128_t hex string vs. literal: %s\n", (u2 - u1).hex());
    u1 = 24197857200151252728969465429440056815_uint128;
    u2 = "24197857200151252728969465429440056815";
    printf("uint128_t decimal string vs. literal: %s\n", (char*)(u2 - u1));

    printf("int128_t:\n");
    int128_t i1 = 0x1234567890abcdef1234567890abcdef_int128;
    int128_t i2 = "0x1234567890abcdef1234567890abcdef";
    printf("int128_t hex string vs. literal: %s\n", (i2 - i1).hex());
    i1 = -24197857200151252728969465429440056815_int128;
    i2 = "-24197857200151252728969465429440056815";
    printf("int128_t decimal string vs. literal: %s\n", (i2 - i1).to_string());
    i1 = -700;

}

void test()
{
    test_precision();
    test_conversion();
    test_addition();
    test_shift();
    test_multiplication();
    test_division();
    test_functions();
    test_string();
    test_comparison();
    test_literals();
}

int main()
{
    test();
    bench();
}
