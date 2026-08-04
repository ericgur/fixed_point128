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
#include <cstring>
#include <chrono>
#include <format>
#include <fstream>
#include <string>
#include <vector>
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

/**
 * @brief Volatile sink written by DoNotOptimize().
 *
 * A store to a volatile object is an observable side effect, so no optimizer is permitted to remove
 * it - nor the computation that produced the stored value.
 */
static volatile unsigned char benchSink = 0;

/**
 * @brief Makes @p value observable so the computation that produced it cannot be eliminated.
 *
 * Every byte of the object is stored to a volatile sink, which forces the complete value to be
 * materialized. Call this after a timed loop on whatever the loop produced; without it the compiler
 * is free to delete the entire loop body as dead code, and both MSVC and Clang do exactly that.
 *
 * The argument is taken by value and the function is deliberately not inlined, so that the address
 * of the *caller's* variable is never taken. Taking it costs real throughput: a by reference,
 * inlined version of this sink pins the accumulator of the cheap operators to a second register
 * pair and MSVC then loses 24% on the addition benchmark, measuring the harness instead of the add.
 *
 * This is deliberately not cheap - it costs sizeof(T) volatile stores - so call it outside the timed
 * loop, never inside it.
 *
 * @tparam T Type of the consumed value.
 * @param value Value to make observable.
 */
template <typename T> FP128_NO_INLINE void DoNotOptimize(T value) noexcept
{
    const unsigned char* bytes = reinterpret_cast<const unsigned char*>(&value);
    for (size_t i = 0; i < sizeof(T); ++i) {
        benchSink = bytes[i];
    }
}

/**
 * @brief Identity function used as the target of MakeOpaque()'s indirect call.
 * @tparam T Type of the value passed through.
 * @param value Value to return unchanged.
 * @return @p value.
 */
template <typename T> [[nodiscard]] FP128_NO_INLINE T Identity(T value) noexcept
{
    return value;
}

/** @brief Pointer to function type used by MakeOpaque(). */
template <typename T> using IdentityFunc = T (*)(T) noexcept;

/**
 * @brief Volatile pointer to Identity<T>(), reloaded from memory at every call.
 *
 * Being volatile is the whole point: the optimizer may not assume which function the pointer
 * designates, so it can neither devirtualize the call nor propagate the argument through it.
 */
template <typename T> volatile IdentityFunc<T> identityPtr = &Identity<T>;

/**
 * @brief Returns @p value through a call the optimizer cannot see through.
 *
 * Used on the *input* of a timed loop. Because the result is opaque, loop invariant code motion
 * cannot hoist an expression computed from it out of the loop, which would otherwise turn a timed
 * loop over a pure function into a single evaluation plus an empty loop.
 *
 * The cost is one volatile load plus an indirect call: measured at roughly 5 cycles under MSVC and
 * 18 under Clang, so it is noise for anything from sqrt() upwards, and a fixed floor under the cheap
 * arithmetic operators - which is why the loop carried benchmarks (addition, subtraction, Mandelbrot)
 * do not use it, being unhoistable on their own.
 *
 * @tparam T Type of the value passed through.
 * @param value Value to hide from the optimizer.
 * @return @p value, unchanged but opaque.
 */
template <typename T> [[nodiscard]] FP128_INLINE T MakeOpaque(T value) noexcept
{
    return identityPtr<T>(value);
}

/**
 * @brief A single benchmark measurement, recorded by print_ips().
 *
 * Measurements are always collected, whether or not a JSON report was requested: recording happens
 * outside every timed loop, so it cannot perturb a result, and collecting unconditionally keeps the
 * benchmark functions free of any knowledge of the output format.
 */
struct BenchResult {
    std::string group;  ///< Title of the group the measurement belongs to.
    std::string name;   ///< Name of the measured operation, as printed to the terminal.
    int64_t ips;        ///< Measured iterations per second.
};

/** @brief Every measurement taken so far, in execution order. */
static std::vector<BenchResult> benchResults;

/** @brief Title of the group currently being benchmarked; tags each new result. */
static std::string currentGroup;

/**
 * @brief Prints a group banner and tags all subsequent results as belonging to that group.
 *
 * The rules above and below the title are sized to the title, which is why the banner is built here
 * rather than being spelled out at each call site.
 *
 * @param title Group title, e.g. "Arithmetic benchmark".
 */
void PrintGroupHeader(const char* title)
{
    currentGroup = title;
    const std::string rule(std::strlen(title), '-');
    printf("\n%s\n%s\n%s\n", rule.c_str(), title, rule.c_str());
}

void print_ips(const char* name, int64_t ips)
{
    benchResults.push_back(BenchResult { currentGroup, name, ips });

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
    PrintGroupHeader("Comparison operator benchmark");

    Duration dur;
    uint64_t total_iterations = 0;
    // setup
    fixed_point128<10> f1 = fixed_point128<10>::e();
    fixed_point128<10> f2 = fixed_point128<10>::golden_ratio();
    int64_t dummy = 0;
    // start the clock
    // Hiding one operand is enough to stop LICM from hoisting the loop invariant comparisons out of
    // the loop, and costs half of what hiding both would. DoNotOptimize(dummy) after the loop keeps
    // the accumulated result observable; the previous sink went through an int64_t overload of the
    // opaque helper that MSVC proved pure and deleted, taking the whole loop body with it.
    dur.start();
    while (dur.cur_duration() < time_per_function) {
        for (uint64_t i = BENCH_ITERATIONS; i != 0; --i) {
            const auto v1 = MakeOpaque(f1);
            dummy += (v1 > f2);
            dummy += (v1 >= f2);
            dummy += (v1 < f2);
            dummy += (v1 <= f2);
        }
        total_iterations += 4 * BENCH_ITERATIONS;
    }
    DoNotOptimize(dummy);

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
    // prevents LICM from hoisting the addition out of the loop. DoNotOptimize(f3)
    // after the loop stores f3 to a volatile sink, forcing it to be live and
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
    DoNotOptimize(f3);
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
    // prevents LICM from hoisting the subtraction out of the loop. DoNotOptimize(f3)
    // after the loop stores f3 to a volatile sink, forcing it to be live and
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
    DoNotOptimize(f3);
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
            f3 = MakeOpaque(f1) * f2;
        }
        total_iterations += BENCH_ITERATIONS;
    }
    DoNotOptimize(f3);
    print_ips("Multiplication by fixed_point128", (uint64_t)(total_iterations / dur.duration()));

    // The multiplicand is hidden behind MakeOpaque() rather than accumulated into a single value.
    // The accumulating form (f10 = f10 * int_val) was degenerate: the low QWORD of the product does
    // not depend on the high QWORD, so with nothing observing the result both compilers proved the
    // high half dead. MSVC collapsed the whole operation into a single 64 bit mulx and Clang went
    // further, vectorizing the remaining chain 4 wide - neither was timing a 128 bit multiply.
    // pi * 123456789 needs 29 integer bits, so it fits fixed_point128<32> without overflowing.
    fixed_point128<32> f10 = fixed_point128<32>::pi();
    fixed_point128<32> f11;
    const uint32_t int_val = 123456789;
    total_iterations = 0;
    dur.start();
    while (dur.cur_duration() < time_per_function) {
        for (uint64_t i = BENCH_ITERATIONS; i != 0; --i) {
            f11 = MakeOpaque(f10) * int_val;
        }
        total_iterations += BENCH_ITERATIONS;
    }
    DoNotOptimize(f11);
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
    // MakeOpaque on the dividend each iteration prevents LICM from hoisting the
    // loop-invariant division out of the loop, without causing value accumulation
    // that would produce degenerate (zero/overflow) inputs.
    double dval = MakeOpaque(64.0);
    dur.start();
    while (dur.cur_duration() < time_per_function) {
        for (uint64_t i = BENCH_ITERATIONS; i != 0; --i) {
            f3 = MakeOpaque(f1) / dval;
        }
        total_iterations += BENCH_ITERATIONS;
    }
    DoNotOptimize(f3);
    print_ips("Division by double (exponent of 2)", (uint64_t)(total_iterations / dur.duration()));

    total_iterations = 0;
    dur.start();
    while (dur.cur_duration() < time_per_function) {
        for (uint64_t i = BENCH_ITERATIONS; i != 0; --i) {
            f3 = MakeOpaque(f1) / 5ll;
        }
        total_iterations += BENCH_ITERATIONS;
    }
    DoNotOptimize(f3);
    print_ips("Division by int64", (uint64_t)(total_iterations / dur.duration()));

    fixed_point128<10> f4 = 5;
    total_iterations = 0;
    dur.start();
    while (dur.cur_duration() < time_per_function) {
        for (uint64_t i = BENCH_ITERATIONS; i != 0; --i) {
            f3 = MakeOpaque(f1) / f4;
        }
        total_iterations += BENCH_ITERATIONS;
    }
    DoNotOptimize(f3);
    print_ips("Division by fixed_point128 (int)", (uint64_t)(total_iterations / dur.duration()));

    total_iterations = 0;
    dur.start();
    while (dur.cur_duration() < time_per_function) {
        for (uint64_t i = BENCH_ITERATIONS; i != 0; --i) {
            f3 = MakeOpaque(f1) / f2;
        }
        total_iterations += BENCH_ITERATIONS;
    }
    DoNotOptimize(f3);
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
            f2 = reciprocal(MakeOpaque(f1));
        }
        total_iterations += BENCH_ITERATIONS;
    }
    DoNotOptimize(f2);

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
            f2 = sqrt(MakeOpaque(f1));
        }
        total_iterations += BENCH_ITERATIONS;
    }
    DoNotOptimize(f2);

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
            f2 = exp(MakeOpaque(f1));
        }
        total_iterations += BENCH_ITERATIONS;
    }
    DoNotOptimize(f2);

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
            f2 = exp2(MakeOpaque(f1));
        }
        total_iterations += BENCH_ITERATIONS;
    }
    DoNotOptimize(f2);

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
            f2 = expm1(MakeOpaque(f1));
        }
        total_iterations += BENCH_ITERATIONS;
    }
    DoNotOptimize(f2);

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
            f3 = pow(MakeOpaque(f1), f2);
        }
        total_iterations += BENCH_ITERATIONS;
    }
    DoNotOptimize(f3);

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
            f2 = log(MakeOpaque(f1));
        }
        total_iterations += BENCH_ITERATIONS;
    }
    DoNotOptimize(f2);

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
            f2 = log2(MakeOpaque(f1));
        }
        total_iterations += BENCH_ITERATIONS;
    }
    DoNotOptimize(f2);

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
            f2 = log10(MakeOpaque(f1));
        }
        total_iterations += BENCH_ITERATIONS;
    }
    DoNotOptimize(f2);

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
            f2 = log1p(MakeOpaque(f1));
        }
        total_iterations += BENCH_ITERATIONS;
    }
    DoNotOptimize(f2);

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
            f2 = sin(MakeOpaque(f1));
        }
        total_iterations += BENCH_ITERATIONS;
    }
    DoNotOptimize(f2);

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
            f2 = asin(MakeOpaque(f1));
        }
        total_iterations += BENCH_ITERATIONS;
    }
    DoNotOptimize(f2);

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
            f2 = cos(MakeOpaque(f1));
        }
        total_iterations += BENCH_ITERATIONS;
    }
    DoNotOptimize(f2);

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
            f2 = acos(MakeOpaque(f1));
        }
        total_iterations += BENCH_ITERATIONS;
    }
    DoNotOptimize(f2);

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
            f2 = tan(MakeOpaque(f1));
        }
        total_iterations += BENCH_ITERATIONS;
    }
    DoNotOptimize(f2);

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
            f2 = atan(MakeOpaque(f1));
        }
        total_iterations += BENCH_ITERATIONS;
    }
    DoNotOptimize(f2);

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
            f2 = sinh(MakeOpaque(f1));
        }
        total_iterations += BENCH_ITERATIONS;
    }
    DoNotOptimize(f2);

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
            f2 = asinh(MakeOpaque(f1));
        }
        total_iterations += BENCH_ITERATIONS;
    }
    DoNotOptimize(f2);

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
            f2 = cosh(MakeOpaque(f1));
        }
        total_iterations += BENCH_ITERATIONS;
    }
    DoNotOptimize(f2);

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
            f2 = acosh(MakeOpaque(f1));
        }
        total_iterations += BENCH_ITERATIONS;
    }
    DoNotOptimize(f2);

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
            f2 = tanh(MakeOpaque(f1));
        }
        total_iterations += BENCH_ITERATIONS;
    }
    DoNotOptimize(f2);

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
            f2 = atanh(MakeOpaque(f1));
        }
        total_iterations += BENCH_ITERATIONS;
    }
    DoNotOptimize(f2);

    print_ips("atanh", (uint64_t)(total_iterations / dur.duration()));
}

void bench_mandelbrot(double time_per_function = 1.0)
{
    Duration dur;
    uint64_t total_iterations = 0;
    // setup
    fixed_point128<10> usq, vsq, tmp, modulus, u, v;
    // A point that doesn't diverge quickly. Both coordinates are hidden from the optimizer so the
    // iteration cannot be constant folded; the orbit itself is loop carried, so nothing inside the
    // loop is hoistable and no per-iteration barrier is needed.
    fixed_point128<10> x = MakeOpaque(fixed_point128<10>(-0.7294734415));
    fixed_point128<10> y = MakeOpaque(fixed_point128<10>(0.242809));
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
    DoNotOptimize(modulus);

    print_ips("Mandelbrot", (uint64_t)(total_iterations / dur.duration()));
}
/**
 * @brief Benches all simple arithmetic functions
 * @param time_per_function Time spent in each sub-test
 */

void bench_arithmetic(double time_per_function = 1.0)
{
    PrintGroupHeader("Arithmetic benchmark");

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
    PrintGroupHeader("Exponents benchmark");
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
    PrintGroupHeader("Logarithmic benchmark");

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
    PrintGroupHeader("Trigonometric benchmark");

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
    PrintGroupHeader("Special function benchmark");

    bench_mandelbrot(time_per_function);
}

void bench_hyperbolic_trig_functions(double time_per_function = 1.0)
{
    PrintGroupHeader("Hyperbolic trigonometric benchmark");

    bench_sinh(time_per_function);
    bench_asinh(time_per_function);
    bench_cosh(time_per_function);
    bench_acosh(time_per_function);
    bench_tanh(time_per_function);
    bench_atanh(time_per_function);
}
/**
 * @brief Returns the human readable name and version of the compiler that built this binary.
 *
 * clang-cl defines both `__clang__` and `_MSC_VER`; `__clang__` is therefore tested first, so a
 * clang-cl build is reported as Clang - the code generator, which is what a benchmark result depends
 * on, rather than the flag syntax.
 *
 * @return Compiler description, e.g. "MSVC 19.44".
 */
[[nodiscard]] std::string CompilerName()
{
#ifdef __clang__
    return format("Clang {}", __clang_version__);
#elif defined(__GNUC__) || defined(__GNUG__)
    return format("GCC {}.{}", __GNUC__, __GNUC_MINOR__);
#elif defined(_MSC_VER)
    return format("MSVC {}.{}", _MSC_VER / 100, _MSC_VER - 100 * (_MSC_VER / 100));
#else
    return "an unknown compiler";
#endif
}

/**
 * @brief Returns a file name safe identifier of the compiler and its version.
 *
 * Unlike CompilerName(), this never embeds `__clang_version__`, which carries the vendor string and
 * repository URL and so contains spaces and slashes.
 *
 * @return Compiler tag, e.g. "clang-19.1.5" or "msvc-19.44".
 */
[[nodiscard]] std::string CompilerTag()
{
#ifdef __clang__
    return format("clang-{}.{}.{}", __clang_major__, __clang_minor__, __clang_patchlevel__);
#elif defined(__GNUC__) || defined(__GNUG__)
    return format("gcc-{}.{}", __GNUC__, __GNUC_MINOR__);
#elif defined(_MSC_VER)
    return format("msvc-{}.{}", _MSC_VER / 100, _MSC_VER - 100 * (_MSC_VER / 100));
#else
    return "unknown";
#endif
}

/**
 * @brief Returns the build type this binary was compiled with.
 *
 * Uses the same condition as the library itself (fixed_point128_shared.h), so the reported build
 * type cannot disagree with the one the measured code was compiled under.
 *
 * @return "debug" or "release".
 */
[[nodiscard]] constexpr const char* BuildType() noexcept
{
#if defined _DEBUG || defined DEBUG
    return "debug";
#else
    return "release";
#endif
}

/**
 * @brief Returns the name of the JSON report for this binary.
 *
 * The name encodes the compiler and the build type so that reports from different toolchains, or
 * from a debug and a release build of the same toolchain, never overwrite one another.
 *
 * @return File name, e.g. "bench_msvc-19.44_release.json".
 */
[[nodiscard]] std::string JsonFileName()
{
    return format("bench_{}_{}.json", CompilerTag(), BuildType());
}

/**
 * @brief Returns the current UTC time as an ISO 8601 timestamp, e.g. "2026-08-04T09:15:42Z".
 */
[[nodiscard]] std::string Timestamp()
{
    const auto now = std::chrono::floor<std::chrono::seconds>(std::chrono::system_clock::now());
    return format("{:%FT%TZ}", now);
}

/**
 * @brief Escapes a string for embedding in a JSON string literal.
 *
 * Benchmark names are ASCII and currently contain nothing that needs escaping, but they are free
 * text; escaping here means adding a name later can never silently produce invalid JSON.
 *
 * @param text String to escape.
 * @return @p text with all JSON control characters escaped.
 */
[[nodiscard]] std::string EscapeJson(const std::string& text)
{
    std::string res;
    res.reserve(text.size());
    for (const char c : text) {
        switch (c) {
        case '"':  res += "\\\""; break;
        case '\\': res += "\\\\"; break;
        case '\b': res += "\\b";  break;
        case '\f': res += "\\f";  break;
        case '\n': res += "\\n";  break;
        case '\r': res += "\\r";  break;
        case '\t': res += "\\t";  break;
        default:
            if (static_cast<unsigned char>(c) < 0x20) {
                res += format("\\u{:04x}", static_cast<unsigned>(static_cast<unsigned char>(c)));
            } else {
                res += c;
            }
            break;
        }
    }

    return res;
}

/**
 * @brief Writes every recorded measurement to @p path as JSON.
 *
 * The report holds the build metadata needed to compare two runs meaningfully - compiler, build type
 * and the timing parameters - followed by the results in execution order, each tagged with the group
 * it was printed under.
 *
 * @param path Destination file, overwritten if it exists.
 * @return true on success, false if the file could not be written.
 */
bool WriteJsonReport(const std::string& path)
{
    std::ofstream file(path, std::ios::binary | std::ios::trunc);
    if (!file) {
        return false;
    }

    file << "{\n";
    file << format("  \"compiler\": \"{}\",\n", EscapeJson(CompilerName()));
    file << format("  \"compilerTag\": \"{}\",\n", CompilerTag());
    file << format("  \"build\": \"{}\",\n", BuildType());
    file << format("  \"timestamp\": \"{}\",\n", Timestamp());
    file << format("  \"timePerFunction\": {},\n", TIME_PER_FUNCTION);
    file << format("  \"benchIterations\": {},\n", BENCH_ITERATIONS);
    file << "  \"results\": [\n";
    for (size_t i = 0; i < benchResults.size(); ++i) {
        const auto& result = benchResults[i];
        const char* separator = (i + 1 < benchResults.size()) ? "," : "";
        file << format("    {{ \"group\": \"{}\", \"name\": \"{}\", \"iterationsPerSecond\": {} }}{}\n",
                       EscapeJson(result.group), EscapeJson(result.name), result.ips, separator);
    }
    file << "  ]\n";
    file << "}\n";
    file.flush();

    return file.good();
}

/**
 * @brief Main benchmark function
 */
void bench()
{
    printf("Compiled with %s\n", CompilerName().c_str());
    printf("=========================\n");
    printf("Single threaded benchmark\n");
    printf("=========================\n");

    // run the function groups
    bench_comparison_operators(TIME_PER_FUNCTION);
    bench_arithmetic(TIME_PER_FUNCTION);
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

/**
 * @brief Prints the command line usage of the benchmark.
 * @param exeName Name the executable was invoked with.
 */
void PrintUsage(const char* exeName)
{
    printf("Usage: %s [-j|--json] [-h|--help]\n", exeName);
    printf("  -j, --json  Write the results to a JSON file in addition to the terminal.\n");
    printf("              The file is named after the compiler and build type, e.g. %s\n", JsonFileName().c_str());
    printf("  -h, --help  Print this help text and exit.\n");
}

int main(int argc, char* argv[])
{
    bool jsonOutput = false;
    for (int i = 1; i < argc; ++i) {
        if (std::strcmp(argv[i], "-j") == 0 || std::strcmp(argv[i], "--json") == 0) {
            jsonOutput = true;
        } else if (std::strcmp(argv[i], "-h") == 0 || std::strcmp(argv[i], "--help") == 0) {
            PrintUsage(argv[0]);
            return 0;
        } else {
            printf("Unknown option: %s\n", argv[i]);
            PrintUsage(argv[0]);
            return 1;
        }
    }

    // Force instantiation of all public methods and friend functions for I=1, 40 and 64.
    force_instantiation<1>();
    force_instantiation<40>();
    force_instantiation<64>();

    bench();

    if (jsonOutput) {
        const std::string path = JsonFileName();
        if (!WriteJsonReport(path)) {
            fprintf(stderr, "\nFailed to write the JSON report to '%s'\n", path.c_str());
            return 1;
        }
        printf("\nJSON report written to %s\n", path.c_str());
    }

    return 0;
}
