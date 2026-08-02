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

// bench.cpp : benchmark and profile for fixed_point128 and uint128_t classes
//

// #define FP128_DISABLE_INLINE TRUE

#include <cstdio>
#include <chrono>
#include <format>
#include "fixed_point128.h"
#include "uint128_t.h"

using namespace std;
using namespace fp128;
constexpr uint64_t BENCH_ITERATIONS = 5000;
constexpr double TIME_PER_FUNCTION = 0.5;  // in seconds

struct Duration {
    using time_point = std::chrono::high_resolution_clock::time_point;
    Duration() = default;

    void start() { t1 = std::chrono::high_resolution_clock::now(); }
    double cur_duration()
    {
        t2 = std::chrono::high_resolution_clock::now();
        return std::chrono::duration<double>(t2 - t1).count();
    }
    double duration() { return std::chrono::duration<double>(t2 - t1).count(); }
    void clear() { t1 = t2 = time_point {}; }
    time_point t1 {}, t2 {};
    inline static double frequency;
};

// get a value that makes the complier not optimize away certain expressions.
template <typename T> __declspec(noinline) T get_const(T val)
{
    return val;
}

void print_ips(const char* name, int64_t ips)
{
    if (ips < 1000) {
        printf("%s: %lld/s\n", name, ips);
    } else if (ips < 1000000) {
        const double dips = ips / 1000.0;
        printf("%s: %0.3lfK/s\n", name, dips);
    } else if (ips < 1000000000) {
        const double dips = ips / 1000000.0;
        printf("%s: %0.3lfM/s\n", name, dips);
    } else {
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
    // get_const on both operands each iteration prevents LICM from hoisting the
    // loop-invariant comparisons. get_const(dummy) as a sink forces dummy to be
    // live, preventing dead-code elimination of the loop body.
    dur.start();
    while (dur.cur_duration() < time_per_function) {
        for (uint64_t i = BENCH_ITERATIONS; i != 0; --i) {
            const auto v1 = get_const(f1);
            const auto v2 = get_const(f2);
            dummy += (v1 > v2);
            dummy += (v1 >= v2);
            dummy += (v1 < v2);
            dummy += (v1 <= v2);
        }
        total_iterations += 4 * BENCH_ITERATIONS;
    }
    dummy = get_const(dummy);  // opaque sink: forces the loop to actually execute

    print_ips("Operators >, >=, <, <= (average of all 4)", (uint64_t)(total_iterations / dur.duration()));
}

void bench_addition(double time_per_function = 1.0)
{
    Duration dur;
    uint64_t total_iterations = 0;
    // setup
    fixed_point128<10> f1 = fabs(fixed_point128<10>::pi());
    fixed_point128<10> f2 = fixed_point128<10>::e();
    // f3 accumulates each iteration, creating a loop-carried dependency that
    // prevents LICM from hoisting the addition out of the loop. get_const(f3)
    // after the loop is an opaque noinline call that forces f3 to be live,
    // preventing dead-code elimination of the entire loop body.
    fixed_point128<10> f3 = f2;
    // start the clock
    dur.start();
    while (dur.cur_duration() < time_per_function) {
        for (uint64_t i = BENCH_ITERATIONS; i != 0; --i) {
            f3 = f3 + f1;
        }
        total_iterations += BENCH_ITERATIONS;
    }
    f3 = get_const(f3);  // opaque sink: forces the loop to actually execute
    print_ips("Addition", (uint64_t)(total_iterations / dur.duration()));
}

void bench_subtraction(double time_per_function = 1.0)
{
    Duration dur;
    uint64_t total_iterations = 0;
    // setup
    fixed_point128<10> f1 = fabs(fixed_point128<10>::pi());
    fixed_point128<10> f2 = fixed_point128<10>::e();
    // f3 accumulates each iteration, creating a loop-carried dependency that
    // prevents LICM from hoisting the subtraction out of the loop. get_const(f3)
    // after the loop is an opaque noinline call that forces f3 to be live,
    // preventing dead-code elimination of the entire loop body.
    fixed_point128<10> f3 = f1;
    // start the clock
    dur.start();
    while (dur.cur_duration() < time_per_function) {
        for (uint64_t i = BENCH_ITERATIONS; i != 0; --i) {
            f3 = f3 - f2;
        }
        total_iterations += BENCH_ITERATIONS;
    }
    f3 = get_const(f3);  // opaque sink: forces the loop to actually execute
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
            f3 = get_const(f1) * f2;
        }
        total_iterations += BENCH_ITERATIONS;
    }
    f3 = get_const(f3);  // opaque sink: forces the loop to actually execute
    print_ips("Multiplication by fixed_point128", (uint64_t)(total_iterations / dur.duration()));

    // Initialize to 1 so the compiler cannot prove f10 is always 0 and eliminate the loop.
    fixed_point128<32> f10 = 1;
    uint32_t int_val = 123456789;
    total_iterations = 0;
    dur.start();
    while (dur.cur_duration() < time_per_function) {
        for (uint64_t i = BENCH_ITERATIONS; i != 0; --i) {
            f10 = f10 * int_val;
        }
        total_iterations += BENCH_ITERATIONS;
    }
    f10 = get_const(f10);  // opaque sink: forces the loop to actually execute
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
    // get_const on the dividend each iteration prevents LICM from hoisting the
    // loop-invariant division out of the loop, without causing value accumulation
    // that would produce degenerate (zero/overflow) inputs.
    double dval = get_const(64.0);
    dur.start();
    while (dur.cur_duration() < time_per_function) {
        for (uint64_t i = BENCH_ITERATIONS; i != 0; --i) {
            f3 = get_const(f1) / dval;
        }
        total_iterations += 2 * BENCH_ITERATIONS;
    }
    f3 = get_const(f3);  // opaque sink: forces the loop to actually execute
    print_ips("Division by double (exponent of 2)", (uint64_t)(total_iterations / dur.duration()));

    total_iterations = 0;
    dur.start();
    while (dur.cur_duration() < time_per_function) {
        for (uint64_t i = BENCH_ITERATIONS; i != 0; --i) {
            f3 = get_const(f1) / 5ll;
        }
        total_iterations += BENCH_ITERATIONS;
    }
    f3 = get_const(f3);  // opaque sink: forces the loop to actually execute
    print_ips("Division by int64", (uint64_t)(total_iterations / dur.duration()));

    fixed_point128<10> f4 = 5;
    total_iterations = 0;
    dur.start();
    while (dur.cur_duration() < time_per_function) {
        for (uint64_t i = BENCH_ITERATIONS; i != 0; --i) {
            f3 = get_const(f1) / f4;
        }
        total_iterations += BENCH_ITERATIONS;
    }
    f3 = get_const(f3);  // opaque sink: forces the loop to actually execute
    print_ips("Division by fixed_point128 (int)", (uint64_t)(total_iterations / dur.duration()));

    total_iterations = 0;
    dur.start();
    while (dur.cur_duration() < time_per_function) {
        for (uint64_t i = BENCH_ITERATIONS; i != 0; --i) {
            f3 = get_const(f1) / f2;
        }
        total_iterations += BENCH_ITERATIONS;
    }
    f3 = get_const(f3);  // opaque sink: forces the loop to actually execute
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
    if (f2) {
        f2++;
    }  // fool the complier into not optimizing away the benchmark

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
    if (f2) {
        f2++;
    }  // fool the complier into not optimizing away the benchmark

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
    if (f2) {
        f2++;
    }  // fool the complier into not optimizing away the benchmark

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
    if (f2) {
        f2++;
    }  // fool the complier into not optimizing away the benchmark

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
    if (f2) {
        f2++;
    }  // fool the complier into not optimizing away the benchmark

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
    if (f3) {
        f3++;
    }  // fool the complier into not optimizing away the benchmark

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
    if (f2) {
        f2++;
    }  // fool the complier into not optimizing away the benchmark

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
    if (f2) {
        f2++;
    }  // fool the complier into not optimizing away the benchmark

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
    if (f2) {
        f2++;
    }  // fool the complier into not optimizing away the benchmark

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
    if (f2) {
        f2++;
    }  // fool the complier into not optimizing away the benchmark

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
    if (f2) {
        f2++;
    }  // fool the complier into not optimizing away the benchmark

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
    if (f2) {
        f2++;
    }  // fool the complier into not optimizing away the benchmark

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
    if (f2) {
        f2++;
    }  // fool the complier into not optimizing away the benchmark

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
    if (f2) {
        f2++;
    }  // fool the complier into not optimizing away the benchmark

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
    if (f2) {
        f2++;
    }  // fool the complier into not optimizing away the benchmark

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
    if (f2) {
        f2++;
    }  // fool the complier into not optimizing away the benchmark

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
    if (f2) {
        f2++;
    }  // fool the complier into not optimizing away the benchmark

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
    if (f2) {
        f2++;
    }  // fool the complier into not optimizing away the benchmark

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
    if (f2) {
        f2++;
    }  // fool the complier into not optimizing away the benchmark

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
    if (f2) {
        f2++;
    }  // fool the complier into not optimizing away the benchmark

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
    if (f2) {
        f2++;
    }  // fool the complier into not optimizing away the benchmark

    print_ips("tanh", (uint64_t)(total_iterations / dur.duration()));
}

void bench_atanh(double time_per_function = 1.0)
{
    Duration dur;
    uint64_t total_iterations = 0;
    // setup
    fixed_point128<10> f1 = fixed_point128<10>::e() / 4;  // abs(v) < 1
    fixed_point128<10> f2;
    // start the clock
    dur.start();
    while (dur.cur_duration() < time_per_function) {
        for (uint64_t i = BENCH_ITERATIONS; i != 0; --i) {
            f2 = tanh(f1);
        }
        total_iterations += BENCH_ITERATIONS;
    }
    if (f2) {
        f2++;
    }  // fool the complier into not optimizing away the benchmark

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
            usq = sqr(u);
            vsq = sqr(v);
            // check uv vector amplitude is smaller than 2
            modulus = usq + vsq;
        }
        total_iterations += BENCH_ITERATIONS;
    }
    if (modulus) {
        modulus++;
    }  // fool the complier into not optimizing away the benchmark

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
#ifdef __clang__
    auto compiler = format("Clang {}", __clang_version__);
#elif defined(__GNUC__) || defined(__GNUG__)
    auto compiler = "GCC";
#elif defined(_MSC_VER)
    auto compiler = format("MSVC {}.{}", _MSC_VER / 100, _MSC_VER - 100 * (_MSC_VER / 100));
#else
    auto compiler = "an unknown compiler";
#endif
    printf("Compiled with %s\n", compiler.c_str());
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

/**
 * @brief Forces the compiler to instantiate all public methods and friend functions of fixed_point128<I>.
 *
 * Calls every constructor, assignment operator, conversion operator, arithmetic operator,
 * comparison operator, query method, static constant accessor, and CRT-style friend math
 * function for the given template parameter I.  Intended to catch compilation errors across
 * the full range of supported I values (1, 40, 64).
 *
 * @tparam I Number of integer bits passed to fixed_point128.
 */
template <int32_t I> FP128_NO_INLINE void force_instantiation()
{
    using fp = fixed_point128<I>;

    // --- Constructors ---
    fp def;                                               // default
    fp from_double(1.5);                                  // double
    fp copy_ctor(from_double);                            // copy
    fp move_ctor(std::move(fp(2.0)));                     // move
    fp from_u64((uint64_t)1);                             // uint64_t
    fp from_i64((int64_t)1);                              // int64_t
    fp from_u32((uint32_t)1);                             // uint32_t
    fp from_i32((int32_t)1);                              // int32_t
    fp from_cstr("1.5");                                  // const char*
    fp from_str(std::string("1.5"));                      // std::string
    fp from_raw(0ull, 1ull, 0u);                          // raw (low, high, sign)

    // cross-template copy constructor (I2 = 10)
    fixed_point128<10> f10(1.5);
    fp cross_ctor(f10);

    // --- Assignment operators ---
    def = copy_ctor;                                      // copy assign
    def = std::move(fp(3.0));                             // move assign
    def = f10;                                            // cross-template assign

    // --- Conversion operators ---
    (void)(uint64_t)from_double;
    (void)(int64_t)from_double;
    (void)(uint32_t)from_double;
    (void)(int32_t)from_double;
    (void)(float)from_double;
    (void)(double)from_double;
    (void)(long double)from_double;
    (void)(std::string)from_double;
    (void)(char*)from_double;
    (void)(bool)from_double;

    // --- Arithmetic compound-assignment operators ---
    fp a = fp::e();
    fp b = fp::golden_ratio();
    fp c;

    c += a;
    c -= a;
    c *= a;
    c /= b;
    c %= b;
    c >>= 1;
    c <<= 1;
    c &= a;
    c |= a;
    c ^= a;
    c.square();

    // compound-assignment with scalar types (exercises template overloads)
    c += 1.5;
    c -= 1.5;
    c *= 2.0;
    c /= 2.0;
    c %= 1.5;
    c *= (uint64_t)2;                                     // operator*=<uint64_t> specialization
    c /= (uint64_t)2;                                     // operator/=<uint64_t> specialization
    c /= (double)2.0;                                     // operator/=<double> specialization

    // --- Binary arithmetic operators (friend) ---
    c = a + b;
    c = a - b;
    c = a * b;
    c = a / b;
    c = a % b;
    c = a >> 1;
    c = a << 1;
    c = a & b;
    c = a | b;
    c = a ^ b;

    // --- Unary operators ---
    c = -a;
    c = +a;
    c = ~a;
    (void)!a;
    ++c;
    --c;
    c++;
    c--;

    // --- Comparison operators ---
    (void)(a == b);
    (void)(a != b);
    (void)(a < b);
    (void)(a <= b);
    (void)(a > b);
    (void)(a >= b);

    // comparison overloads with scalar type T
    (void)(a == 1.5);
    (void)(1.5 == a);
    (void)(a != 1.5);
    (void)(1.5 != a);
    (void)(a < 1.5);
    (void)(1.5 < a);
    (void)(a <= 1.5);
    (void)(1.5 <= a);
    (void)(a > 1.5);
    (void)(1.5 > a);
    (void)(a >= 1.5);
    (void)(1.5 >= a);

    // --- Query methods ---
    (void)a.is_int();
    (void)a.is_positive();
    (void)a.is_negative();
    (void)a.is_zero();
    (void)a.get_bit(0);
    (void)a.get_exponent();

    // --- Static constant accessors ---
    (void)fp::pi();
    (void)fp::pi2();
    (void)fp::half_pi();
    (void)fp::golden_ratio();
    (void)fp::e();
    (void)fp::sqrt_2();
    (void)fp::one();
    (void)fp::half();
    (void)fp::epsilon();

    // --- Friend math functions (CRT-style) ---
    fp val = fp::e();
    fp half_val = fp::half();

    (void)fabs(val);
    (void)floor(val);
    (void)ceil(val);
    (void)trunc(val);
    (void)round(val);
    (void)ilogb(val);
    (void)copysign(val, val);
    (void)fmod(val, val);
    fp iptr;
    (void)modf(val, &iptr);
    (void)fdim(val, val);
    (void)fmin(val, val);
    (void)fmax(val, val);
    (void)hypot(val, val);
    (void)sqr(val);
    (void)sqrt(val);
    (void)exp(val);
    (void)exp2(val);
    (void)expm1(val);
    (void)pow(val, val);
    (void)log(val);
    (void)log2(val);
    (void)log10(val);
    (void)logb(val);
    (void)log1p(val);
    (void)lzcnt128(val);
    (void)reciprocal(val);
    fp fact_res;
    fact_reciprocal(5, fact_res);
    // friend trigonometric functions (require minimum template parameter I >= 4)
    if constexpr (I >= 4) {
        (void)sin(val);
        (void)asin(half_val);                                 // |x| <= 1 required
        (void)cos(val);
        (void)acos(half_val);                                 // |x| <= 1 required
        (void)tan(val);
        (void)atan(val);
        (void)atan2(val, val);
        (void)sinh(val);
        (void)asinh(val);
        (void)cosh(val);
        (void)acosh(fp::e());  // x >= 1 required
        (void)tanh(val);
        (void)atanh(half_val);  // |x| < 1 required
    }
    
    // suppress unused-variable warnings
    (void)def; (void)from_double; (void)copy_ctor; (void)move_ctor;
    (void)from_u64; (void)from_i64; (void)from_u32; (void)from_i32;
    (void)from_cstr; (void)from_str; (void)from_raw; (void)cross_ctor;
    (void)a; (void)b; (void)c; (void)f10; (void)fact_res;
}

int main()
{
    // Force instantiation of all public methods and friend functions for I=1, 40 and 64.
    force_instantiation<1>();
    force_instantiation<40>();
    force_instantiation<64>();

    bench();
    return 0;
}
