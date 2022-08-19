/***********************************************************************************
    MIT License

    Copyright (c) 2022 Eric Gur (ericgur@iname.com)

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

//#define FP128_DISABLE_INLINE

#include <windows.h>
#include <stdio.h>
#include <profileapi.h>
#include "../inc/fixed_point128.h" 

using namespace fp128;

void test_conversion()
{
    printf("\nTest Conversion\n");
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
    fixed_point128<20> f1 = 1.0 / 3.0;
    printf("f1: %0.15lf\n", (double)f1);
    fixed_point128<20> f2 = f1 / 2.0;
    printf("f2 = f1 / 2.0: %0.15lf\n", (double)f2);
    f2 /= 32.0;
    printf("f2 = f1 / 32.0: %0.15lf\n", (double)f2);
    f1 = 0.01 / 3.0;
    printf("f1 = 0.01 / 3.0 %0.15lf\n", (double)f1);
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
    double values[] = {3.45, -7.5, 0.27};
    std::string s;
    fixed_point128<10> f1;
    int len = sizeof(values) / sizeof(double);
    for (int i = 0; i < len; ++i) {
        f1 = values[i];
        printf("f1: %0.15lf\n", (double)f1);
        printf("floor(f1): %0.15lf\n", (double)floor(f1));
        printf("ciel(f1): %0.15lf\n", (double)ciel(f1));
        printf("fabs(f1): %0.15lf\n", (double)fabs(f1));
        fixed_point128<10> int_part;
        fixed_point128<10> f2 = modf(f1, &int_part);
        printf("modf(f1): int: %lf fraction %0.15lf\n", (double)int_part, (double)f2);
       
        f2 = 0.25;
        printf("f2: %0.15lf\n", (double)f2);
        printf("fmod(f1, f2): %0.15lf\n", (double)fmod(f1, f2));
        printf("sqrt(f1): %0.15lf\n", (double)sqrt(f1));

        printf("\n");
    }

    fixed_point128<4> fvalues[] = {
        fixed_point128<4>::pi(),
        fixed_point128<4>::pi() >> 1,
        fixed_point128<4>::pi() >> 2,
        (-fixed_point128<4>::pi()) >> 2
    };

    fixed_point128<4> f4;
    len = sizeof(fvalues) / sizeof(fvalues[0]);
    for (int i = 0; i < len; ++i) {
        f4 = fvalues[i];
        s = f4;
        printf("f4: %s\n", s.c_str());
        s = sin(f4);
        printf("sin(f4): %s\n", s.c_str());
        s = cos(f4);
        printf("cos(f4): %s\n", s.c_str());
    }
    
    printf("\n");
}

void test_comparison()
{
    printf("\nTest comparison\n");
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

void print_ips(const char* name, int64 ips)
{
    if (ips < 1000) {
        printf("%s: %lld per second\n", name, ips);
    }
    else if (ips < 1000000) {
        double dips = ips / 1000.0;
        printf("%s: %0.3lfK per second\n", name, dips);
    }
    else if (ips < 1000000000) {
        double dips = ips / 1000000.0;
        printf("%s: %0.3lfM per second\n", name, dips);
    }
    else {
        double dips = ips / 1000000000.0;
        printf("%s: %0.3lfG per second\n", name, dips);
    }
}

void bench()
{
    LARGE_INTEGER li, time_start, time_end;
    QueryPerformanceFrequency(&li);
    double totalTime, frequency = double(li.QuadPart);
    int iterations = 1000000000;
    uint64 ips = 0;
    fixed_point128<10> f1 = fixed_point128<10>::pi();
    fixed_point128<10> f2 = fixed_point128<10>::e();
    fixed_point128<10> f3;

    QueryPerformanceCounter(&time_start);
    for (int i = 0; i < iterations; ++i)
        f3 = f1 + f2;
    QueryPerformanceCounter(&time_end);
    totalTime = (time_end.QuadPart - time_start.QuadPart) / frequency;
    ips = (uint64)(iterations / totalTime);
    print_ips("Addition", ips);

    QueryPerformanceCounter(&time_start);
    for (int i = 0; i < iterations; ++i)
        f3 = f2 - f1;
    QueryPerformanceCounter(&time_end);
    totalTime = (time_end.QuadPart - time_start.QuadPart) / frequency;
    ips = (uint64)(iterations / totalTime);
    print_ips("Subtraction", ips);

    QueryPerformanceCounter(&time_start);
    for (int i = 0; i < iterations; ++i)
        f3 = f1 * f2;
    QueryPerformanceCounter(&time_end);
    totalTime = (time_end.QuadPart - time_start.QuadPart) / frequency;
    ips = (uint64)(iterations / totalTime);
    print_ips("Multiplication", ips);


    // slower functions
    iterations /= 10;
    QueryPerformanceCounter(&time_start);
    for (int i = 0; i < iterations; ++i)
        f3 = f1 / f2;
    QueryPerformanceCounter(&time_end);
    totalTime = (time_end.QuadPart - time_start.QuadPart) / frequency;
    ips = (uint64)(iterations / totalTime);
    print_ips("Division", ips);
}

int main()
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
    bench();
}

