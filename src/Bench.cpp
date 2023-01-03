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

// bench.cpp : benchamrk and profile for fixed_point128 and uint128_t classes
//

//#define FP128_DISABLE_INLINE TRUE

#include <windows.h>
#include <profileapi.h>
#include <cstdio>
#include <cassert>
#include "../inc/fixed_point128.h" 
#include "../inc/uint128_t.h" 

#pragma warning(disable: 26493) // Don't use C-style casts.
#pragma warning(disable: 26467) 
#pragma warning(disable: 26485) 
#pragma warning(disable: 26440) 
#pragma warning(disable: 26482) 
#pragma warning(disable: 26446) 
#pragma warning(disable: 26496) 

using namespace fp128;

void print_ips(const char* name, int64_t ips)
{
    if (ips < 1000) {
        printf("%s: %lld/s\n", name, ips);
    }
    else if (ips < 1000000) {
        const double dips = ips / 1000.0;
        printf("%s: %0.3lfK/s\n", name, dips);
    }
    else if (ips < 1000000000) {
        const double dips = ips / 1000000.0;
        printf("%s: %0.3lfM/s\n", name, dips);
    }
    else {
        const double dips = ips / 1000000000.0;
        printf("%s: %0.3lfG/s\n", name, dips);
    }
}

void bench()
{
    LARGE_INTEGER li, time_start{}, time_end{};
    QueryPerformanceFrequency(&li);
    double totalTime = 0;
    const double frequency = double(li.QuadPart);
    int iterations = 2000000000;
    uint64_t ips = 0;
    fixed_point128<10> f1 = fabs(fixed_point128<10>::pi());
    fixed_point128<10> f2 = fixed_point128<10>::e();
    fixed_point128<10> f3;

    QueryPerformanceCounter(&time_start);
    for (int i = 0; i < iterations; ++i) {
        ips += (f1 > f2) || (f1 >= f2) || (f1 < f2) || (f1 <= f2);
    }
    QueryPerformanceCounter(&time_end);
    totalTime = (time_end.QuadPart - time_start.QuadPart) / frequency;
    if (ips > 5) { // trick compiler to not optimzie out the above loop
        printf("");
    }
    ips = (uint64_t)(iterations / totalTime);
    print_ips("Operators >, >=, <, <=", ips);


    QueryPerformanceCounter(&time_start);
    for (int i = 0; i < iterations; ++i) {
        f3 = f1 + f2;
    }
    QueryPerformanceCounter(&time_end);
    totalTime = (time_end.QuadPart - time_start.QuadPart) / frequency;
    ips = (uint64_t)(iterations / totalTime);
    print_ips("Addition", ips);


    QueryPerformanceCounter(&time_start);
    for (int i = 0; i < iterations; ++i) {
        f3 = f2 - f1;
    }
    QueryPerformanceCounter(&time_end);
    totalTime = (time_end.QuadPart - time_start.QuadPart) / frequency;
    ips = (uint64_t)(iterations / totalTime);
    print_ips("Subtraction", ips);


    QueryPerformanceCounter(&time_start);
    for (int i = 0; i < iterations; ++i)
        f3 = f1 * f2;
    QueryPerformanceCounter(&time_end);
    totalTime = (time_end.QuadPart - time_start.QuadPart) / frequency;
    ips = (uint64_t)(iterations / totalTime);
    print_ips("Multiplication by fixed_point128", ips);

    fixed_point128<32> f10;
    uint32_t int_val = 123456789;
    QueryPerformanceCounter(&time_start);
    for (int i = 0; i < iterations; ++i)
        f10 = f10 * int_val;
    QueryPerformanceCounter(&time_end);
    totalTime = (time_end.QuadPart - time_start.QuadPart) / frequency;
    ips = (uint64_t)(iterations / totalTime);
    print_ips("Multiplication by int", ips);


    iterations /= 2;
    double values[2] = { 64, 64 };
    QueryPerformanceCounter(&time_start);
    for (int i = 0; i < iterations; ++i) {
        f3 = f1 / values[0]; // trick the compiler to not optimize away this code to nothing
        f3 = f1 / values[1];
    }
    iterations *= 2;
    QueryPerformanceCounter(&time_end);
    totalTime = (time_end.QuadPart - time_start.QuadPart) / frequency;
    ips = (uint64_t)(iterations / totalTime);
    print_ips("Division by float (exponent of 2)", ips);


    // slower functions
    //iterations /= 30;
    iterations /= 10;
    fixed_point128<10> f4 = 5;

    QueryPerformanceCounter(&time_start);
    for (int i = 0; i < iterations; ++i)
        f3 = f1 / f4;
    QueryPerformanceCounter(&time_end);
    totalTime = (time_end.QuadPart - time_start.QuadPart) / frequency;
    ips = (uint64_t)(iterations / totalTime);
    print_ips("Division by int", ips);

    // even slower
    iterations /= 5;

    QueryPerformanceCounter(&time_start);
    for (int i = 0; i < iterations; ++i)
        f3 = f1 / f2;
    QueryPerformanceCounter(&time_end);
    totalTime = (time_end.QuadPart - time_start.QuadPart) / frequency;
    ips = (uint64_t)(iterations / totalTime);
    print_ips("Division by fixed_point128", ips);

    QueryPerformanceCounter(&time_start);
    for (int i = 0; i < iterations; ++i)
        f3 = reciprocal(f2);
    QueryPerformanceCounter(&time_end);
    totalTime = (time_end.QuadPart - time_start.QuadPart) / frequency;
    ips = (uint64_t)(iterations / totalTime);
    print_ips("reciprocal", ips);

    // even slower
    iterations /= 5;

    QueryPerformanceCounter(&time_start);
    for (int i = 0; i < iterations; ++i)
        f3 = sqrt(f1);
    QueryPerformanceCounter(&time_end);
    totalTime = (time_end.QuadPart - time_start.QuadPart) / frequency;
    ips = (uint64_t)(iterations / totalTime);
    print_ips("sqrt", ips);


    QueryPerformanceCounter(&time_start);
    for (int i = 0; i < iterations; ++i)
        f3 = exp(f1);
    QueryPerformanceCounter(&time_end);
    totalTime = (time_end.QuadPart - time_start.QuadPart) / frequency;
    ips = (uint64_t)(iterations / totalTime);
    print_ips("exp", ips);


    QueryPerformanceCounter(&time_start);
    for (int i = 0; i < iterations; ++i)
        f3 = log(f1);
    QueryPerformanceCounter(&time_end);
    totalTime = (time_end.QuadPart - time_start.QuadPart) / frequency;
    ips = (uint64_t)(iterations / totalTime);
    print_ips("log", ips);

    
    QueryPerformanceCounter(&time_start);
    for (int i = 0; i < iterations; ++i)
        f3 = pow(f1, f2);
    QueryPerformanceCounter(&time_end);
    totalTime = (time_end.QuadPart - time_start.QuadPart) / frequency;
    ips = (uint64_t)(iterations / totalTime);
    print_ips("pow", ips);
}
