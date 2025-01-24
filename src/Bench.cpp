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
#include <chrono>
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
constexpr uint64_t BENCH_ITERATIONS = 5000;
constexpr double TIME_PER_FUNCTION = 1.0;

struct Duration
{
    Duration() : t1(), t2() {
        if (frequency == 0) {
            LARGE_INTEGER li;
            QueryPerformanceFrequency(&li);
            frequency = double(li.QuadPart);
        }
    }
    void start(){ QueryPerformanceCounter(&t1); }
    double cur_duration() { 
        QueryPerformanceCounter(&t2);
        return (t2.QuadPart - t1.QuadPart) / frequency;
    }
    double duration() {
        return (t2.QuadPart - t1.QuadPart) / frequency;
    }
    void clear() {
        t2.QuadPart = t1.QuadPart = 0;
    }
    LARGE_INTEGER t1, t2;
    inline static double frequency;
};

// get a value that makes the complier not optimize away certain expressions.
template<typename T> __declspec(noinline) T get_const(T val) {
    return val;
}

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

/**
 * @brief Benches all comparison functions
 * @param time_per_function Time spent in each sub-test
*/
void bench_comparison_operators(double time_per_function = 1.0)
{
    printf("\n");
    printf("-----------------------------\n");
    printf("Comparison operator benchmark\n");
    printf("-----------------------------\n");

    Duration dur;
    uint64_t total_iterations = 0;
    // setup
    fixed_point128<10> f1 = fixed_point128<10>::e();
    fixed_point128<10> f2 = fixed_point128<10>::golden_ratio();
    int64_t dummy = 0;
    // start the clock
    dur.start();
    while (dur.cur_duration() < time_per_function) {
        for (uint64_t i = BENCH_ITERATIONS; i != 0; --i) {
            dummy += (f1 > f2);
            dummy += (f1 >= f2);
            dummy += (f1 < f2);
            dummy += (f1 <= f2);
        }
        total_iterations += 4 * BENCH_ITERATIONS;
    }
    if (dummy > 5) { // trick the compiler to not optimize out the above loop
        printf("");
    }

    print_ips("Operators >, >=, <, <= (average of all 4)", (uint64_t)(total_iterations / dur.duration()));

    //dur.start();
    //while (dur.cur_duration() < time_per_function) {
    //    for (uint64_t i = BENCH_ITERATIONS; i != 0; --i) {
    //        dummy -= (f1 == f2) || (f1 != f2);
    //    }
    //    total_iterations += BENCH_ITERATIONS;
    //}
    //if (dummy <  5) { // trick the compiler to not optimize out the above loop
    //    printf("");
    //}

    //print_ips("Operators ==, !=> (both in one line)", (uint64_t)(total_iterations / dur.duration()));
}

void bench_addition(double time_per_function = 1.0)
{
    Duration dur;
    uint64_t total_iterations = 0;
    // setup
    fixed_point128<10> f1 = fabs(fixed_point128<10>::pi());
    fixed_point128<10> f2 = fixed_point128<10>::e();
    fixed_point128<10> f3;
    // start the clock
    dur.start();
    while (dur.cur_duration() < time_per_function) {
        for (uint64_t i = BENCH_ITERATIONS; i != 0; --i) {
            f3 = f2 + f1;
        }
        total_iterations += BENCH_ITERATIONS;
    }
    if (f3) { f3++; } // fool the complier into not optimizing away the benchmark
    print_ips("Addition", (uint64_t)(total_iterations / dur.duration()));
}

void bench_subtraction(double time_per_function = 1.0)
{
    Duration dur;
    uint64_t total_iterations = 0;
    // setup
    fixed_point128<10> f1 = fabs(fixed_point128<10>::pi());
    fixed_point128<10> f2 = fixed_point128<10>::e();
    fixed_point128<10> f3;
    // start the clock
    dur.start();
    while (dur.cur_duration() < time_per_function) {
        for (uint64_t i = BENCH_ITERATIONS; i != 0; --i) {
            f3 = f1 - f2;
        }
        total_iterations += BENCH_ITERATIONS;
    }
    if (f3) { f3++; } // fool the complier into not optimizing away the benchmark

    print_ips("Subtraction", (uint64_t)(total_iterations / dur.duration()));
}

void bench_multiplication(double time_per_function = 1.0)
{
    Duration dur;
    uint64_t total_iterations = 0;
    srand(0x12345678);
    // setup
    fixed_point128<8> f1 = fabs(fixed_point128<10>::pi());
    fixed_point128<8> f2 = (double)rand() / 1.0101010101010101;
    fixed_point128<8> f3;

    // start the clock
    total_iterations = 0;
    dur.start();
    while (dur.cur_duration() < time_per_function) {
        for (uint64_t i = BENCH_ITERATIONS; i != 0; --i) {
            f3 = f1 * f2;
        }
        total_iterations += BENCH_ITERATIONS;
    }
    if (f3) { f3++; } // fool the complier into not optimizing away the benchmark
    print_ips("Multiplication by fixed_point128", (uint64_t)(total_iterations / dur.duration()));

    fixed_point128<32> f10;
    uint32_t int_val = 123456789;
    total_iterations = 0;
    dur.start();
    while (dur.cur_duration() < time_per_function) {
        for (uint64_t i = BENCH_ITERATIONS; i != 0; --i) {
            f10 = f10 * int_val;
        }
        total_iterations += BENCH_ITERATIONS;
    }
    if (f10) { f10++; } // fool the complier into not optimizing away the benchmark
    print_ips("Multiplication by int32_t", (uint64_t)(total_iterations / dur.duration()));
}

void bench_division(double time_per_function = 1.0)
{
    Duration dur;
    uint64_t total_iterations = 0;
    // setup
    fixed_point128<10> f1 = fabs(fixed_point128<10>::pi());
    fixed_point128<10> f2 = fixed_point128<10>::e();
    fixed_point128<10> f3;
    total_iterations = 0;

    // start the clock
    double dval = get_const(64.0);
    dur.start();
    while (dur.cur_duration() < time_per_function) {
        for (uint64_t i = BENCH_ITERATIONS; i != 0; --i) {
            f3 = f1 / dval; // fix 
        }
        total_iterations += 2 * BENCH_ITERATIONS;
    }
    if (f3) { f3++; } // fool the complier into not optimizing away the benchmark
    print_ips("Division by double (exponent of 2)", (uint64_t)(total_iterations / dur.duration()));

    total_iterations = 0;
    dur.start();
    while (dur.cur_duration() < time_per_function) {
        for (uint64_t i = BENCH_ITERATIONS; i != 0; --i) {
            f3 = f1 / 5ll;
        }
        total_iterations += BENCH_ITERATIONS;
    }
    if (f3) { f3++; } // fool the complier into not optimizing away the benchmark
    print_ips("Division by int64", (uint64_t)(total_iterations / dur.duration()));


    fixed_point128<10> f4 = 5;
    total_iterations = 0;
    dur.start();
    while (dur.cur_duration() < time_per_function) {
        for (uint64_t i = BENCH_ITERATIONS; i != 0; --i) {
            f3 = f1 / f4;
        }
        total_iterations += BENCH_ITERATIONS;
    }
    if (f3) { f3++; } // fool the complier into not optimizing away the benchmark
    print_ips("Division by fixed_point128 (int)", (uint64_t)(total_iterations / dur.duration()));

    total_iterations = 0;
    dur.start();
    while (dur.cur_duration() < time_per_function) {
        for (uint64_t i = BENCH_ITERATIONS; i != 0; --i) {
            f3 = f1 / f2;
        }
        total_iterations += BENCH_ITERATIONS;
    }
    if (f3) { f3++; } // fool the complier into not optimizing away the benchmark
    print_ips("Division by fixed_point128 (float)", (uint64_t)(total_iterations / dur.duration()));
}

void bench_reciprocal(double time_per_function = 1.0)
{
    Duration dur;
    uint64_t total_iterations = 0;
    // setup
    fixed_point128<10> f1 = fixed_point128<10>::e();
    fixed_point128<10> f2;
    // start the clock
    dur.start();
    while (dur.cur_duration() < time_per_function) {
        for (uint64_t i = BENCH_ITERATIONS; i != 0; --i) {
            f2 = reciprocal(f1);
        }
        total_iterations += BENCH_ITERATIONS;
    }
    if (f2) { f2++; } // fool the complier into not optimizing away the benchmark

    print_ips("reciprocal", (uint64_t)(total_iterations / dur.duration()));
}

void bench_sqrt(double time_per_function = 1.0)
{
    Duration dur;
    uint64_t total_iterations = 0;
    // setup
    fixed_point128<10> f1 = fixed_point128<10>::e();
    fixed_point128<10> f2;
    // start the clock
    dur.start();
    while (dur.cur_duration() < time_per_function) {
        for (uint64_t i = BENCH_ITERATIONS; i != 0; --i) {
            f2 = sqrt(f1);
        }
        total_iterations += BENCH_ITERATIONS;
    }
    if (f2) { f2++; } // fool the complier into not optimizing away the benchmark

    print_ips("sqrt", (uint64_t)(total_iterations / dur.duration()));
}

void bench_exp(double time_per_function = 1.0)
{
    Duration dur;
    uint64_t total_iterations = 0;
    // setup
    fixed_point128<10> f1 = fixed_point128<10>::e();
    fixed_point128<10> f2;
    // start the clock
    dur.start();
    while (dur.cur_duration() < time_per_function) {
        for (uint64_t i = BENCH_ITERATIONS; i != 0; --i) {
            f2 = exp(f1);
        }
        total_iterations += BENCH_ITERATIONS;
    }
    if (f2) { f2++; } // fool the complier into not optimizing away the benchmark

    print_ips("exp", (uint64_t)(total_iterations / dur.duration()));
}

void bench_exp2(double time_per_function = 1.0)
{
    Duration dur;
    uint64_t total_iterations = 0;
    // setup
    fixed_point128<10> f1 = fixed_point128<10>::e();
    fixed_point128<10> f2;
    // start the clock
    dur.start();
    while (dur.cur_duration() < time_per_function) {
        for (uint64_t i = BENCH_ITERATIONS; i != 0; --i) {
            f2 = exp2(f1);
        }
        total_iterations += BENCH_ITERATIONS;
    }
    if (f2) { f2++; } // fool the complier into not optimizing away the benchmark

    print_ips("exp2", (uint64_t)(total_iterations / dur.duration()));
}

void bench_expm1(double time_per_function = 1.0)
{
    Duration dur;
    uint64_t total_iterations = 0;
    // setup
    fixed_point128<10> f1 = fixed_point128<10>::e();
    fixed_point128<10> f2;
    // start the clock
    dur.start();
    while (dur.cur_duration() < time_per_function) {
        for (uint64_t i = BENCH_ITERATIONS; i != 0; --i) {
            f2 = expm1(f1);
        }
        total_iterations += BENCH_ITERATIONS;
    }
    if (f2) { f2++; } // fool the complier into not optimizing away the benchmark

    print_ips("expm1", (uint64_t)(total_iterations / dur.duration()));
}

void bench_pow(double time_per_function = 1.0)
{
    Duration dur;
    uint64_t total_iterations = 0;
    // setup
    fixed_point128<10> f1 = fixed_point128<10>::e();
    fixed_point128<10> f2 = fixed_point128<10>::golden_ratio();
    fixed_point128<10> f3; 
    // start the clock
    dur.start();
    while (dur.cur_duration() < time_per_function) {
        for (uint64_t i = BENCH_ITERATIONS; i != 0; --i) {
            f3 = pow(f1, f2);
        }
        total_iterations += BENCH_ITERATIONS;
    }
    if (f3) { f3++; } // fool the complier into not optimizing away the benchmark

    print_ips("pow", (uint64_t)(total_iterations / dur.duration()));
}


void bench_log(double time_per_function = 1.0)
{
    Duration dur;
    uint64_t total_iterations = 0;
    // setup
    fixed_point128<10> f1 = fixed_point128<10>::e();
    fixed_point128<10> f2;
    // start the clock
    dur.start();
    while (dur.cur_duration() < time_per_function) {
        for (uint64_t i = BENCH_ITERATIONS; i != 0; --i) {
            f2 = log(f1);
        }
        total_iterations += BENCH_ITERATIONS;
    }
    if (f2) { f2++; } // fool the complier into not optimizing away the benchmark

    print_ips("log", (uint64_t)(total_iterations / dur.duration()));
}

void bench_log2(double time_per_function = 1.0)
{
    Duration dur;
    uint64_t total_iterations = 0;
    // setup
    fixed_point128<10> f1 = fixed_point128<10>::e();
    fixed_point128<10> f2;
    // start the clock
    dur.start();
    while (dur.cur_duration() < time_per_function) {
        for (uint64_t i = BENCH_ITERATIONS; i != 0; --i) {
            f2 = log2(f1);
        }
        total_iterations += BENCH_ITERATIONS;
    }
    if (f2) { f2++; } // fool the complier into not optimizing away the benchmark

    print_ips("log2", (uint64_t)(total_iterations / dur.duration()));
}

void bench_log10(double time_per_function = 1.0)
{
    Duration dur;
    uint64_t total_iterations = 0;
    // setup
    fixed_point128<10> f1 = fixed_point128<10>::e();
    fixed_point128<10> f2;
    // start the clock
    dur.start();
    while (dur.cur_duration() < time_per_function) {
        for (uint64_t i = BENCH_ITERATIONS; i != 0; --i) {
            f2 = log10(f1);
        }
        total_iterations += BENCH_ITERATIONS;
    }
    if (f2) { f2++; } // fool the complier into not optimizing away the benchmark

    print_ips("log10", (uint64_t)(total_iterations / dur.duration()));
}

void bench_log1p(double time_per_function = 1.0)
{
    Duration dur;
    uint64_t total_iterations = 0;
    // setup
    fixed_point128<10> f1 = fixed_point128<10>::e();
    fixed_point128<10> f2;
    // start the clock
    dur.start();
    while (dur.cur_duration() < time_per_function) {
        for (uint64_t i = BENCH_ITERATIONS; i != 0; --i) {
            f2 = log1p(f1);
        }
        total_iterations += BENCH_ITERATIONS;
    }
    if (f2) { f2++; } // fool the complier into not optimizing away the benchmark

    print_ips("log1p", (uint64_t)(total_iterations / dur.duration()));
}

void bench_sin(double time_per_function = 1.0)
{
    Duration dur;
    uint64_t total_iterations = 0;
    // setup
    fixed_point128<10> f1 = fixed_point128<10>::e() / 2;
    fixed_point128<10> f2;
    // start the clock
    dur.start();
    while (dur.cur_duration() < time_per_function) {
        for (uint64_t i = BENCH_ITERATIONS; i != 0; --i) {
            f2 = sin(f1);
        }
        total_iterations += BENCH_ITERATIONS;
    }
    if (f2) { f2++; } // fool the complier into not optimizing away the benchmark

    print_ips("sin", (uint64_t)(total_iterations / dur.duration()));
}

void bench_asin(double time_per_function = 1.0)
{
    Duration dur;
    uint64_t total_iterations = 0;
    // setup
    fixed_point128<10> f1 = fixed_point128<10>::pi() / 5;
    fixed_point128<10> f2;
    // start the clock
    dur.start();
    while (dur.cur_duration() < time_per_function) {
        for (uint64_t i = BENCH_ITERATIONS; i != 0; --i) {
            f2 = asin(f1);
        }
        total_iterations += BENCH_ITERATIONS;
    }
    if (f2) { f2++; } // fool the complier into not optimizing away the benchmark

    print_ips("asin", (uint64_t)(total_iterations / dur.duration()));
}

void bench_cos(double time_per_function = 1.0)
{
    Duration dur;
    uint64_t total_iterations = 0;
    // setup
    fixed_point128<10> f1 = fixed_point128<10>::e() / 2;
    fixed_point128<10> f2;
    // start the clock
    dur.start();
    while (dur.cur_duration() < time_per_function) {
        for (uint64_t i = BENCH_ITERATIONS; i != 0; --i) {
            f2 = cos(f1);
        }
        total_iterations += BENCH_ITERATIONS;
    }
    if (f2) { f2++; } // fool the complier into not optimizing away the benchmark

    print_ips("cos", (uint64_t)(total_iterations / dur.duration()));
}

void bench_acos(double time_per_function = 1.0)
{
    Duration dur;
    uint64_t total_iterations = 0;
    // setup
    fixed_point128<10> f1 = fixed_point128<10>::pi() / 5;
    fixed_point128<10> f2;
    // start the clock
    dur.start();
    while (dur.cur_duration() < time_per_function) {
        for (uint64_t i = BENCH_ITERATIONS; i != 0; --i) {
            f2 = acos(f1);
        }
        total_iterations += BENCH_ITERATIONS;
    }
    if (f2) { f2++; } // fool the complier into not optimizing away the benchmark

    print_ips("acos", (uint64_t)(total_iterations / dur.duration()));
}

void bench_tan(double time_per_function = 1.0)
{
    Duration dur;
    uint64_t total_iterations = 0;
    // setup
    fixed_point128<10> f1 = fixed_point128<10>::e() / 2;
    fixed_point128<10> f2;
    // start the clock
    dur.start();
    while (dur.cur_duration() < time_per_function) {
        for (uint64_t i = BENCH_ITERATIONS; i != 0; --i) {
            f2 = tan(f1);
        }
        total_iterations += BENCH_ITERATIONS;
    }
    if (f2) { f2++; } // fool the complier into not optimizing away the benchmark

    print_ips("tan", (uint64_t)(total_iterations / dur.duration()));
}

void bench_atan(double time_per_function = 1.0)
{
    Duration dur;
    uint64_t total_iterations = 0;
    // setup
    fixed_point128<10> f1 = fixed_point128<10>::pi() / 5;
    fixed_point128<10> f2;
    // start the clock
    dur.start();
    while (dur.cur_duration() < time_per_function) {
        for (uint64_t i = BENCH_ITERATIONS; i != 0; --i) {
            f2 = atan(f1);
        }
        total_iterations += BENCH_ITERATIONS;
    }
    if (f2) { f2++; } // fool the complier into not optimizing away the benchmark

    print_ips("atan", (uint64_t)(total_iterations / dur.duration()));
}

void bench_sinh(double time_per_function = 1.0)
{
    Duration dur;
    uint64_t total_iterations = 0;
    // setup
    fixed_point128<10> f1 = fixed_point128<10>::e() / 2;
    fixed_point128<10> f2;
    // start the clock
    dur.start();
    while (dur.cur_duration() < time_per_function) {
        for (uint64_t i = BENCH_ITERATIONS; i != 0; --i) {
            f2 = sinh(f1);
        }
        total_iterations += BENCH_ITERATIONS;
    }
    if (f2) { f2++; } // fool the complier into not optimizing away the benchmark

    print_ips("sinh", (uint64_t)(total_iterations / dur.duration()));
}

void bench_asinh(double time_per_function = 1.0)
{
    Duration dur;
    uint64_t total_iterations = 0;
    // setup
    fixed_point128<10> f1 = fixed_point128<10>::pi() / 5;
    fixed_point128<10> f2;
    // start the clock
    dur.start();
    while (dur.cur_duration() < time_per_function) {
        for (uint64_t i = BENCH_ITERATIONS; i != 0; --i) {
            f2 = asinh(f1);
        }
        total_iterations += BENCH_ITERATIONS;
    }
    if (f2) { f2++; } // fool the complier into not optimizing away the benchmark

    print_ips("asinh", (uint64_t)(total_iterations / dur.duration()));
}

void bench_cosh(double time_per_function = 1.0)
{
    Duration dur;
    uint64_t total_iterations = 0;
    // setup
    fixed_point128<10> f1 = fixed_point128<10>::e() / 2;
    fixed_point128<10> f2;
    // start the clock
    dur.start();
    while (dur.cur_duration() < time_per_function) {
        for (uint64_t i = BENCH_ITERATIONS; i != 0; --i) {
            f2 = cosh(f1);
        }
        total_iterations += BENCH_ITERATIONS;
    }
    if (f2) { f2++; } // fool the complier into not optimizing away the benchmark

    print_ips("cosh", (uint64_t)(total_iterations / dur.duration()));
}

void bench_acosh(double time_per_function = 1.0)
{
    Duration dur;
    uint64_t total_iterations = 0;
    // setup
    fixed_point128<10> f1 = fixed_point128<10>::e() / 2;
    fixed_point128<10> f2;
    // start the clock
    dur.start();
    while (dur.cur_duration() < time_per_function) {
        for (uint64_t i = BENCH_ITERATIONS; i != 0; --i) {
            f2 = acosh(f1);
        }
        total_iterations += BENCH_ITERATIONS;
    }
    if (f2) { f2++; } // fool the complier into not optimizing away the benchmark

    print_ips("acosh", (uint64_t)(total_iterations / dur.duration()));
}

void bench_tanh(double time_per_function = 1.0)
{
    Duration dur;
    uint64_t total_iterations = 0;
    // setup
    fixed_point128<10> f1 = fixed_point128<10>::e() / 2;
    fixed_point128<10> f2;
    // start the clock
    dur.start();
    while (dur.cur_duration() < time_per_function) {
        for (uint64_t i = BENCH_ITERATIONS; i != 0; --i) {
            f2 = tanh(f1);
        }
        total_iterations += BENCH_ITERATIONS;
    }
    if (f2) { f2++; } // fool the complier into not optimizing away the benchmark

    print_ips("tanh", (uint64_t)(total_iterations / dur.duration()));
}

void bench_atanh(double time_per_function = 1.0)
{
    Duration dur;
    uint64_t total_iterations = 0;
    // setup
    fixed_point128<10> f1 = fixed_point128<10>::e() / 4; // abs(v) < 1
    fixed_point128<10> f2;
    // start the clock
    dur.start();
    while (dur.cur_duration() < time_per_function) {
        for (uint64_t i = BENCH_ITERATIONS; i != 0; --i) {
            f2 = tanh(f1);
        }
        total_iterations += BENCH_ITERATIONS;
    }
    if (f2) { f2++; } // fool the complier into not optimizing away the benchmark

    print_ips("atanh", (uint64_t)(total_iterations / dur.duration()));
}

void bench_mandelbrot(double time_per_function = 1.0)
{
    Duration dur;
    uint64_t total_iterations = 0;
    // setup
    fixed_point128<10> usq, vsq, tmp, modulus, u, v;
    // a point that doesn't diverge quickly
    fixed_point128<10> x = -0.7294734415;
    fixed_point128<10> y = 0.242809;
    // start the clock
    dur.start();
    while (dur.cur_duration() < time_per_function) {
        for (uint64_t i = BENCH_ITERATIONS; i != 0; --i) {
            // real
            tmp = usq - vsq + x;

            // imaginary
            // v = 2.0 * (u * v) + y;
            v = ((u * v) << 1) + y;
            u = tmp;
            usq = u * u;
            vsq = v * v;
            // check uv vector amplitude is smaller than 2
            modulus = usq + vsq;

        }
        total_iterations += BENCH_ITERATIONS;
    }
    if (modulus) { modulus++; } // fool the complier into not optimizing away the benchmark

    print_ips("Mandelbrot", (uint64_t)(total_iterations / dur.duration()));
}
/**
 * @brief Benches all simple arithmatic functions
 * @param time_per_function Time spent in each sub-test
*/

void bench_arithmatic(double time_per_function = 1.0)
{
    printf("\n");
    printf("--------------------\n");
    printf("Arithmatic benchmark\n");
    printf("--------------------\n");

    bench_addition(time_per_function);
    bench_subtraction(time_per_function);
    bench_multiplication(time_per_function);
    bench_division(time_per_function);
    bench_reciprocal(time_per_function);
}

/**
 * @brief Benches all exponent functions
 * @param time_per_function Time spent in each sub-test
*/

void bench_exponents(double time_per_function = 1.0)
{
    printf("\n");
    printf("-------------------\n");
    printf("Exponents benchmark\n");
    printf("-------------------\n");
    bench_sqrt(time_per_function);
    bench_exp(time_per_function);
    bench_exp2(time_per_function);
    bench_pow(time_per_function);
    bench_expm1(time_per_function);
}

/**
 * @brief Benches all log functions
 * @param time_per_function Time spent in each sub-test
*/
void bench_log_functions(double time_per_function = 1.0)
{
    printf("\n");
    printf("---------------------\n");
    printf("Logarithmic benchmark\n");
    printf("---------------------\n");

    bench_log(time_per_function);
    bench_log2(time_per_function);
    bench_log10(time_per_function);
    bench_log1p(time_per_function);
}

/**
 * @brief Benches all trig functions
 * @param time_per_function Time spent in each sub-test
*/
void bench_trig_functions(double time_per_function = 1.0)
{
    printf("\n");
    printf("----------------------\n");
    printf("Trigonometic benchmark\n");
    printf("----------------------\n");

    bench_sin(time_per_function);
    bench_asin(time_per_function);
    bench_cos(time_per_function);
    bench_acos(time_per_function);
    bench_tan(time_per_function);
    bench_atan(time_per_function);
}

/**
 * @brief Benches all special functions
 * @param time_per_function Time spent in each sub-test
*/
void bench_special_functions(double time_per_function = 1.0)
{
    printf("\n");
    printf("--------------------------\n");
    printf("Special function benchmark\n");
    printf("--------------------------\n");

    bench_mandelbrot(time_per_function);
}

void bench_hyperbolic_trig_functions(double time_per_function = 1.0)
{
    printf("\n");
    printf("---------------------------------\n");
    printf("Hyperbolic trigonometic benchmark\n");
    printf("---------------------------------\n");

    bench_sinh(time_per_function);
    bench_asinh(time_per_function);
    bench_cosh(time_per_function);
    bench_acosh(time_per_function);
    bench_tanh(time_per_function);
    bench_atanh(time_per_function);
}
/**
 * @brief Main benchmark function
*/
void bench() 
{
    printf("=========================\n");
    printf("Single threaded benchmark\n");
    printf("=========================\n");

    // run the function groups
    bench_comparison_operators(TIME_PER_FUNCTION);
    bench_arithmatic(TIME_PER_FUNCTION);
    bench_exponents(TIME_PER_FUNCTION);
    bench_log_functions(TIME_PER_FUNCTION);
    bench_trig_functions(TIME_PER_FUNCTION);
    bench_hyperbolic_trig_functions(TIME_PER_FUNCTION);
    bench_special_functions(TIME_PER_FUNCTION);
}
