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

// bench.cpp : benchmark and profile for the fixed_point128, float128, int128_t and uint128_t classes
//

// #define FP128_DISABLE_INLINE TRUE

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <algorithm>
#include <chrono>
#include <format>
#include <fstream>
#include <string>
#include <vector>
#include "fixed_point128.h"
#include "float128.h"
#include "int128_t.h"
#include "uint128_t.h"

using namespace std;
using namespace fp128;
constexpr uint64_t BENCH_ITERATIONS = 5000;
constexpr double TIME_PER_FUNCTION = 0.5;  // in seconds

/**
 * @brief Number of integer bits of the benchmarked fixed_point128 instantiation.
 *
 * Kept at the value the single type benchmark has always used, so the fixed_point128 numbers stay
 * comparable with earlier runs: the count of fraction bits sets how many iterations the series based
 * functions need, and changing it would move every transcendental result.
 */
constexpr int32_t BENCH_FP_INT_BITS = 10;

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
    std::string type;   ///< Type the measurement was taken on, e.g. "float128".
    std::string group;  ///< Title of the group the measurement belongs to.
    std::string name;   ///< Name of the measured operation, as printed to the terminal.
    int64_t ips;        ///< Measured iterations per second.
};

/** @brief Every measurement taken so far, in execution order. */
static std::vector<BenchResult> benchResults;

/** @brief Name of the type currently being benchmarked; tags each new result. */
static std::string currentType;

/** @brief Title of the group currently being benchmarked; tags each new result. */
static std::string currentGroup;

/**
 * @brief Prints a type banner and tags all subsequent results as measurements of that type.
 * @param name Type name, e.g. "uint128_t".
 */
void PrintTypeHeader(const std::string& name)
{
    currentType = name;
    const std::string title = format("Benchmarking {}", name);
    const std::string rule(title.size(), '=');
    printf("\n%s\n%s\n%s\n", rule.c_str(), title.c_str(), rule.c_str());
}

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
    benchResults.push_back(BenchResult { currentType, currentGroup, name, ips });

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
 * @brief Per type benchmark configuration: display name, capabilities and operands.
 *
 * Capabilities are stated explicitly instead of being detected with a `requires` expression. All
 * four types convert *implicitly* to double, so a detection based check would silently succeed for
 * a function the type does not have - `requires { exp(v); }` is satisfied for uint128_t by the CRT's
 * exp(double) - and the benchmark would then time the CRT while labelling the result uint128_t.
 *
 * The three general operands are ordered by magnitude and map onto the constants the single type
 * benchmark has always used for the fractional types, which keeps those results comparable with
 * earlier runs:
 *   operandA() - |pi|, the first operand of add/sub/div and the argument of the inverse trig tests
 *   operandB() - e, the second operand of add/sub/div and the argument of most unary functions
 *   operandC() - the golden ratio (sqrt(2) for float128, which has no golden_ratio()), used as the
 *                comparison right hand side and the pow() exponent
 * The integer types have no such constants, so they supply three positive values in the same order.
 *
 * @tparam T Benchmarked type.
 */
template <typename T> struct BenchTraits;

/**
 * @brief Shared configuration of the two 128 bit integer types.
 *
 * Neither type has a fraction or any transcendental function, so both skip the same benchmarks.
 * Operands are built from the (low, high) constructor rather than parsed from strings, and the top
 * bit is left clear so the same values are positive when reused for the signed type.
 *
 * @tparam T uint128_t or int128_t.
 */
template <typename T> struct IntegerBenchTraits {
    static constexpr bool isFractional = false;      ///< Cannot represent values between 0 and 1.
    static constexpr bool hasTranscendental = false; ///< No exp/trig/reciprocal family.

    using MulType = T;     ///< Type used by the 128 bit multiplication benchmark.
    using IntMulType = T;  ///< Type used by the multiply by int32_t benchmark.

    [[nodiscard]] static T operandA() noexcept { return T(0xFEDCBA9876543210ull, 0x0123456789ABCDEFull); }
    [[nodiscard]] static T operandB() noexcept { return T(0x0F0F0F0F0F0F0F0Full, 0x00FEDCBA98765432ull); }
    [[nodiscard]] static T operandC() noexcept { return T(0x9E3779B97F4A7C15ull, 0ull); }
    [[nodiscard]] static MulType mulA() noexcept { return operandA(); }
    [[nodiscard]] static MulType mulB() noexcept { return operandC(); }
    /** @brief 57 bit value, so multiplying it by a 27 bit int cannot overflow 128 bits. */
    [[nodiscard]] static IntMulType intMulA() noexcept { return T(0x0123456789ABCDEFull, 0ull); }
};

/** @brief Benchmark configuration of uint128_t. */
template <> struct BenchTraits<uint128_t> : IntegerBenchTraits<uint128_t> {
    [[nodiscard]] static std::string name() { return "uint128_t"; }
};

/** @brief Benchmark configuration of int128_t. */
template <> struct BenchTraits<int128_t> : IntegerBenchTraits<int128_t> {
    [[nodiscard]] static std::string name() { return "int128_t"; }
};

/**
 * @brief Benchmark configuration of fixed_point128.
 *
 * The multiplication benchmarks keep the instantiations the single type benchmark used: I=8 for the
 * 128 bit multiply, and I=32 for the multiply by int32_t, where pi * 123456789 needs 29 integer bits
 * and would otherwise overflow.
 *
 * @tparam I Number of integer bits.
 */
template <int32_t I> struct BenchTraits<fixed_point128<I>> {
    using T = fixed_point128<I>;
    static constexpr bool isFractional = true;
    static constexpr bool hasTranscendental = true;

    using MulType = fixed_point128<8>;
    using IntMulType = fixed_point128<32>;

    [[nodiscard]] static std::string name() { return format("fixed_point128<{}>", I); }
    [[nodiscard]] static T operandA() noexcept { return fabs(T::pi()); }
    [[nodiscard]] static T operandB() noexcept { return T::e(); }
    [[nodiscard]] static T operandC() noexcept { return T::golden_ratio(); }
    [[nodiscard]] static MulType mulA() noexcept { return MulType(fabs(T::pi())); }
    [[nodiscard]] static MulType mulB() noexcept
    {
        srand(0x12345678);
        return MulType((double)rand() / 1.0101010101010101);
    }
    [[nodiscard]] static IntMulType intMulA() noexcept { return IntMulType::pi(); }
};

/**
 * @brief Benchmark configuration of float128.
 *
 * float128 has no golden_ratio(), so operandC() is sqrt(2): the benchmarks only need a constant
 * larger than one that differs from the other two operands.
 */
template <> struct BenchTraits<float128> {
    using T = float128;
    static constexpr bool isFractional = true;
    static constexpr bool hasTranscendental = true;

    using MulType = float128;
    using IntMulType = float128;

    [[nodiscard]] static std::string name() { return "float128"; }
    [[nodiscard]] static T operandA() noexcept { return fabs(T::pi()); }
    [[nodiscard]] static T operandB() noexcept { return T::e(); }
    [[nodiscard]] static T operandC() noexcept { return T::sqrt_2(); }
    [[nodiscard]] static MulType mulA() noexcept { return fabs(T::pi()); }
    [[nodiscard]] static MulType mulB() noexcept
    {
        srand(0x12345678);
        return MulType((double)rand() / 1.0101010101010101);
    }
    [[nodiscard]] static IntMulType intMulA() noexcept { return T::pi(); }
};

/**
 * @brief Benches all comparison functions
 * @tparam T Benchmarked type.
 * @param time_per_function Time spent in each sub-test
 */
template <typename T> void bench_comparison_operators(double time_per_function = 1.0)
{
    using Traits = BenchTraits<T>;
    PrintGroupHeader("Comparison operator benchmark");

    Duration dur;
    uint64_t total_iterations = 0;
    // setup
    T f1 = Traits::operandB();
    T f2 = Traits::operandC();
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

template <typename T> void bench_addition(double time_per_function = 1.0)
{
    using Traits = BenchTraits<T>;
    Duration dur;
    uint64_t total_iterations = 0;
    // setup
    T f1 = Traits::operandA();
    T f2 = Traits::operandB();
    // f3 accumulates each iteration, creating a loop-carried dependency that
    // prevents LICM from hoisting the addition out of the loop. DoNotOptimize(f3)
    // after the loop stores f3 to a volatile sink, forcing it to be live and
    // preventing dead-code elimination of the entire loop body.
    T f3 = f2;
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

template <typename T> void bench_subtraction(double time_per_function = 1.0)
{
    using Traits = BenchTraits<T>;
    Duration dur;
    uint64_t total_iterations = 0;
    // setup
    T f1 = Traits::operandA();
    T f2 = Traits::operandB();
    // f3 accumulates each iteration, creating a loop-carried dependency that
    // prevents LICM from hoisting the subtraction out of the loop. DoNotOptimize(f3)
    // after the loop stores f3 to a volatile sink, forcing it to be live and
    // preventing dead-code elimination of the entire loop body.
    T f3 = f1;
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

template <typename T> void bench_multiplication(double time_per_function = 1.0)
{
    using Traits = BenchTraits<T>;
    using MulType = typename Traits::MulType;
    using IntMulType = typename Traits::IntMulType;
    Duration dur;
    uint64_t total_iterations = 0;
    // setup
    MulType f1 = Traits::mulA();
    MulType f2 = Traits::mulB();
    MulType f3;

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
    print_ips("Multiplication by 128-bit value", (uint64_t)(total_iterations / dur.duration()));

    // The multiplicand is hidden behind MakeOpaque() rather than accumulated into a single value.
    // The accumulating form (f10 = f10 * int_val) was degenerate: the low QWORD of the product does
    // not depend on the high QWORD, so with nothing observing the result both compilers proved the
    // high half dead. MSVC collapsed the whole operation into a single 64 bit mulx and Clang went
    // further, vectorizing the remaining chain 4 wide - neither was timing a 128 bit multiply.
    IntMulType f10 = Traits::intMulA();
    IntMulType f11;
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

template <typename T> void bench_division(double time_per_function = 1.0)
{
    using Traits = BenchTraits<T>;
    Duration dur;
    uint64_t total_iterations = 0;
    // setup
    T f1 = Traits::operandA();
    T f3;
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

    T f4 = 5;
    total_iterations = 0;
    dur.start();
    while (dur.cur_duration() < time_per_function) {
        for (uint64_t i = BENCH_ITERATIONS; i != 0; --i) {
            f3 = MakeOpaque(f1) / f4;
        }
        total_iterations += BENCH_ITERATIONS;
    }
    DoNotOptimize(f3);
    print_ips("Division by 128-bit integer value", (uint64_t)(total_iterations / dur.duration()));

    // Only meaningful where the divisor can hold a fraction. For the integer types this loop would
    // be the previous one with a different divisor.
    if constexpr (Traits::isFractional) {
        T f5 = Traits::operandB();
        total_iterations = 0;
        dur.start();
        while (dur.cur_duration() < time_per_function) {
            for (uint64_t i = BENCH_ITERATIONS; i != 0; --i) {
                f3 = MakeOpaque(f1) / f5;
            }
            total_iterations += BENCH_ITERATIONS;
        }
        DoNotOptimize(f3);
        print_ips("Division by 128-bit fractional value", (uint64_t)(total_iterations / dur.duration()));
    }
}

template <typename T> void bench_reciprocal(double time_per_function = 1.0)
{
    Duration dur;
    uint64_t total_iterations = 0;
    // setup
    T f1 = BenchTraits<T>::operandB();
    T f2;
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

/**
 * @brief Benches sqrt().
 *
 * The result is held in `auto`: the integer types return the truncated integer root as a uint64_t
 * rather than a value of their own type.
 */
template <typename T> void bench_sqrt(double time_per_function = 1.0)
{
    Duration dur;
    uint64_t total_iterations = 0;
    // setup
    T f1 = BenchTraits<T>::operandB();
    auto f2 = sqrt(f1);
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

template <typename T> void bench_exp(double time_per_function = 1.0)
{
    Duration dur;
    uint64_t total_iterations = 0;
    // setup
    T f1 = BenchTraits<T>::operandB();
    T f2;
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

template <typename T> void bench_exp2(double time_per_function = 1.0)
{
    Duration dur;
    uint64_t total_iterations = 0;
    // setup
    T f1 = BenchTraits<T>::operandB();
    T f2;
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

template <typename T> void bench_expm1(double time_per_function = 1.0)
{
    Duration dur;
    uint64_t total_iterations = 0;
    // setup
    T f1 = BenchTraits<T>::operandB();
    T f2;
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

/**
 * @brief Benches pow().
 *
 * The integer types only offer pow(x, uint32_t) - raising to a real power needs exp() and log(),
 * which they do not have - so they are timed on a small base and integer exponent instead. The two
 * measurements are not comparable across the categories, which is what the type field of each result
 * is there to make clear.
 */
template <typename T> void bench_pow(double time_per_function = 1.0)
{
    using Traits = BenchTraits<T>;
    Duration dur;
    uint64_t total_iterations = 0;
    T f3;
    // start the clock
    if constexpr (Traits::hasTranscendental) {
        T f1 = Traits::operandB();
        T f2 = Traits::operandC();
        dur.start();
        while (dur.cur_duration() < time_per_function) {
            for (uint64_t i = BENCH_ITERATIONS; i != 0; --i) {
                f3 = pow(MakeOpaque(f1), f2);
            }
            total_iterations += BENCH_ITERATIONS;
        }
    } else {
        // A small base keeps the result inside 128 bits; the cost of the binary exponentiation
        // depends on the exponent, not on the base.
        T f1 = T(7);
        dur.start();
        while (dur.cur_duration() < time_per_function) {
            for (uint64_t i = BENCH_ITERATIONS; i != 0; --i) {
                f3 = pow(MakeOpaque(f1), 5u);
            }
            total_iterations += BENCH_ITERATIONS;
        }
    }
    DoNotOptimize(f3);

    print_ips((Traits::hasTranscendental) ? "pow" : "pow (integer exponent)",
              (uint64_t)(total_iterations / dur.duration()));
}

/**
 * @brief Benches log().
 *
 * As with sqrt(), the integer types return a uint64_t: the integer part of the logarithm.
 */
template <typename T> void bench_log(double time_per_function = 1.0)
{
    Duration dur;
    uint64_t total_iterations = 0;
    // setup
    T f1 = BenchTraits<T>::operandB();
    auto f2 = log(f1);
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

template <typename T> void bench_log2(double time_per_function = 1.0)
{
    Duration dur;
    uint64_t total_iterations = 0;
    // setup
    T f1 = BenchTraits<T>::operandB();
    auto f2 = log2(f1);
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

template <typename T> void bench_log10(double time_per_function = 1.0)
{
    Duration dur;
    uint64_t total_iterations = 0;
    // setup
    T f1 = BenchTraits<T>::operandB();
    auto f2 = log10(f1);
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

template <typename T> void bench_log1p(double time_per_function = 1.0)
{
    Duration dur;
    uint64_t total_iterations = 0;
    // setup
    T f1 = BenchTraits<T>::operandB();
    T f2;
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

template <typename T> void bench_sin(double time_per_function = 1.0)
{
    Duration dur;
    uint64_t total_iterations = 0;
    // setup
    T f1 = BenchTraits<T>::operandB() / 2;
    T f2;
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

template <typename T> void bench_asin(double time_per_function = 1.0)
{
    Duration dur;
    uint64_t total_iterations = 0;
    // setup
    T f1 = BenchTraits<T>::operandA() / 5;
    T f2;
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

template <typename T> void bench_cos(double time_per_function = 1.0)
{
    Duration dur;
    uint64_t total_iterations = 0;
    // setup
    T f1 = BenchTraits<T>::operandB() / 2;
    T f2;
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

template <typename T> void bench_acos(double time_per_function = 1.0)
{
    Duration dur;
    uint64_t total_iterations = 0;
    // setup
    T f1 = BenchTraits<T>::operandA() / 5;
    T f2;
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

template <typename T> void bench_tan(double time_per_function = 1.0)
{
    Duration dur;
    uint64_t total_iterations = 0;
    // setup
    T f1 = BenchTraits<T>::operandB() / 2;
    T f2;
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

template <typename T> void bench_atan(double time_per_function = 1.0)
{
    Duration dur;
    uint64_t total_iterations = 0;
    // setup
    T f1 = BenchTraits<T>::operandA() / 5;
    T f2;
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

template <typename T> void bench_sinh(double time_per_function = 1.0)
{
    Duration dur;
    uint64_t total_iterations = 0;
    // setup
    T f1 = BenchTraits<T>::operandB() / 2;
    T f2;
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

template <typename T> void bench_asinh(double time_per_function = 1.0)
{
    Duration dur;
    uint64_t total_iterations = 0;
    // setup
    T f1 = BenchTraits<T>::operandA() / 5;
    T f2;
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

template <typename T> void bench_cosh(double time_per_function = 1.0)
{
    Duration dur;
    uint64_t total_iterations = 0;
    // setup
    T f1 = BenchTraits<T>::operandB() / 2;
    T f2;
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

template <typename T> void bench_acosh(double time_per_function = 1.0)
{
    Duration dur;
    uint64_t total_iterations = 0;
    // setup
    T f1 = BenchTraits<T>::operandB() / 2;  // x >= 1
    T f2;
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

template <typename T> void bench_tanh(double time_per_function = 1.0)
{
    Duration dur;
    uint64_t total_iterations = 0;
    // setup
    T f1 = BenchTraits<T>::operandB() / 2;
    T f2;
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

template <typename T> void bench_atanh(double time_per_function = 1.0)
{
    Duration dur;
    uint64_t total_iterations = 0;
    // setup
    T f1 = BenchTraits<T>::operandB() / 4;  // abs(v) < 1
    T f2;
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

template <typename T> void bench_mandelbrot(double time_per_function = 1.0)
{
    Duration dur;
    uint64_t total_iterations = 0;
    // setup
    T usq, vsq, tmp, modulus, u, v;
    // A point that doesn't diverge quickly. Both coordinates are hidden from the optimizer so the
    // iteration cannot be constant folded; the orbit itself is loop carried, so nothing inside the
    // loop is hoistable and no per-iteration barrier is needed.
    T x = MakeOpaque(T(-0.7294734415));
    T y = MakeOpaque(T(0.242809));
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
 * @tparam T Benchmarked type.
 * @param time_per_function Time spent in each sub-test
 */
template <typename T> void bench_arithmetic(double time_per_function = 1.0)
{
    PrintGroupHeader("Arithmetic benchmark");

    bench_addition<T>(time_per_function);
    bench_subtraction<T>(time_per_function);
    bench_multiplication<T>(time_per_function);
    bench_division<T>(time_per_function);
    // reciprocal() returns 1/x, which the integer types cannot represent.
    if constexpr (BenchTraits<T>::hasTranscendental) {
        bench_reciprocal<T>(time_per_function);
    }
}

/**
 * @brief Benches all exponent functions
 * @tparam T Benchmarked type.
 * @param time_per_function Time spent in each sub-test
 */
template <typename T> void bench_exponents(double time_per_function = 1.0)
{
    PrintGroupHeader("Exponents benchmark");

    bench_sqrt<T>(time_per_function);
    if constexpr (BenchTraits<T>::hasTranscendental) {
        bench_exp<T>(time_per_function);
        bench_exp2<T>(time_per_function);
    }
    bench_pow<T>(time_per_function);
    if constexpr (BenchTraits<T>::hasTranscendental) {
        bench_expm1<T>(time_per_function);
    }
}

/**
 * @brief Benches all log functions
 * @tparam T Benchmarked type.
 * @param time_per_function Time spent in each sub-test
 */
template <typename T> void bench_log_functions(double time_per_function = 1.0)
{
    PrintGroupHeader("Logarithmic benchmark");

    bench_log<T>(time_per_function);
    bench_log2<T>(time_per_function);
    bench_log10<T>(time_per_function);
    // log1p() is log(1 + x), which is only interesting where x can be a fraction.
    if constexpr (BenchTraits<T>::hasTranscendental) {
        bench_log1p<T>(time_per_function);
    }
}

/**
 * @brief Benches all trig functions
 * @tparam T Benchmarked type.
 * @param time_per_function Time spent in each sub-test
 */
template <typename T> void bench_trig_functions(double time_per_function = 1.0)
{
    PrintGroupHeader("Trigonometric benchmark");

    bench_sin<T>(time_per_function);
    bench_asin<T>(time_per_function);
    bench_cos<T>(time_per_function);
    bench_acos<T>(time_per_function);
    bench_tan<T>(time_per_function);
    bench_atan<T>(time_per_function);
}

/**
 * @brief Benches all hyperbolic trig functions
 * @tparam T Benchmarked type.
 * @param time_per_function Time spent in each sub-test
 */
template <typename T> void bench_hyperbolic_trig_functions(double time_per_function = 1.0)
{
    PrintGroupHeader("Hyperbolic trigonometric benchmark");

    bench_sinh<T>(time_per_function);
    bench_asinh<T>(time_per_function);
    bench_cosh<T>(time_per_function);
    bench_acosh<T>(time_per_function);
    bench_tanh<T>(time_per_function);
    bench_atanh<T>(time_per_function);
}

/**
 * @brief Benches all special functions
 * @tparam T Benchmarked type.
 * @param time_per_function Time spent in each sub-test
 */
template <typename T> void bench_special_functions(double time_per_function = 1.0)
{
    PrintGroupHeader("Special function benchmark");

    bench_mandelbrot<T>(time_per_function);
}

/**
 * @brief Runs every benchmark that is meaningful for @p T.
 *
 * Groups that would end up empty are skipped whole, banner included, so the output never advertises
 * a section with nothing under it.
 *
 * @tparam T Benchmarked type.
 * @param time_per_function Time spent in each sub-test
 */
template <typename T> void BenchOneType(double time_per_function)
{
    PrintTypeHeader(BenchTraits<T>::name());

    bench_comparison_operators<T>(time_per_function);
    bench_arithmetic<T>(time_per_function);
    bench_exponents<T>(time_per_function);
    bench_log_functions<T>(time_per_function);
    if constexpr (BenchTraits<T>::hasTranscendental) {
        bench_trig_functions<T>(time_per_function);
        bench_hyperbolic_trig_functions<T>(time_per_function);
    }
    // The Mandelbrot orbit lives inside the unit disc, so it needs a type with a fraction.
    if constexpr (BenchTraits<T>::isFractional) {
        bench_special_functions<T>(time_per_function);
    }
}

/**
 * @brief The types the benchmark can measure, as selected on the command line.
 *
 * The members are prefixed to keep them from shadowing the type names themselves, which are all
 * visible here through `using namespace fp128`.
 */
struct TypeSelection {
    bool runUint128 = false;
    bool runInt128 = false;
    bool runFixedPoint128 = false;
    bool runFloat128 = false;

    /** @brief Returns true when at least one type is selected. */
    [[nodiscard]] bool any() const noexcept { return runUint128 || runInt128 || runFixedPoint128 || runFloat128; }
    /** @brief Selects every type. */
    void selectAll() noexcept { runUint128 = runInt128 = runFixedPoint128 = runFloat128 = true; }
};

/**
 * @brief Adds the type named by @p name to @p selection.
 * @param name Type name as given on the command line.
 * @param selection Selection to update.
 * @return true if @p name was recognized.
 */
bool SelectType(const char* name, TypeSelection& selection)
{
    if (std::strcmp(name, "all") == 0) {
        selection.selectAll();
    } else if (std::strcmp(name, "uint128") == 0 || std::strcmp(name, "uint128_t") == 0) {
        selection.runUint128 = true;
    } else if (std::strcmp(name, "int128") == 0 || std::strcmp(name, "int128_t") == 0) {
        selection.runInt128 = true;
    } else if (std::strcmp(name, "fixed_point128") == 0 || std::strcmp(name, "fp128") == 0) {
        selection.runFixedPoint128 = true;
    } else if (std::strcmp(name, "float128") == 0 || std::strcmp(name, "f128") == 0) {
        selection.runFloat128 = true;
    } else {
        return false;
    }

    return true;
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
 * from a debug and a release build of the same toolchain, never overwrite one another. It does not
 * encode the selected types: which types a report covers is recorded inside it, in the "types"
 * array and in the type field of every result.
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
 * @brief Returns the distinct type names present in the results, in the order they were measured.
 */
[[nodiscard]] std::vector<std::string> MeasuredTypes()
{
    std::vector<std::string> types;
    for (const auto& result : benchResults) {
        if (std::find(types.begin(), types.end(), result.type) == types.end()) {
            types.push_back(result.type);
        }
    }

    return types;
}

/**
 * @brief Writes every recorded measurement to @p path as JSON.
 *
 * The report holds the build metadata needed to compare two runs meaningfully - compiler, build type
 * and the timing parameters - followed by the results in execution order, each tagged with the type
 * it was measured on and the group it was printed under.
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

    const auto types = MeasuredTypes();
    std::string typeList;
    for (size_t i = 0; i < types.size(); ++i) {
        typeList += format("{}\"{}\"", (i == 0) ? "" : ", ", EscapeJson(types[i]));
    }

    file << "{\n";
    file << format("  \"compiler\": \"{}\",\n", EscapeJson(CompilerName()));
    file << format("  \"compilerTag\": \"{}\",\n", CompilerTag());
    file << format("  \"build\": \"{}\",\n", BuildType());
    file << format("  \"timestamp\": \"{}\",\n", Timestamp());
    file << format("  \"timePerFunction\": {},\n", TIME_PER_FUNCTION);
    file << format("  \"benchIterations\": {},\n", BENCH_ITERATIONS);
    file << format("  \"types\": [{}],\n", typeList);
    file << "  \"results\": [\n";
    for (size_t i = 0; i < benchResults.size(); ++i) {
        const auto& result = benchResults[i];
        const char* separator = (i + 1 < benchResults.size()) ? "," : "";
        file << format("    {{ \"type\": \"{}\", \"group\": \"{}\", \"name\": \"{}\", \"iterationsPerSecond\": {} }}{}\n",
                       EscapeJson(result.type), EscapeJson(result.group), EscapeJson(result.name),
                       result.ips, separator);
    }
    file << "  ]\n";
    file << "}\n";
    file.flush();

    return file.good();
}

/**
 * @brief Main benchmark function
 * @param types Types to measure.
 */
void bench(const TypeSelection& types)
{
    printf("Compiled with %s\n", CompilerName().c_str());
    printf("=========================\n");
    printf("Single threaded benchmark\n");
    printf("=========================\n");

    // run the selected types
    if (types.runUint128) {
        BenchOneType<uint128_t>(TIME_PER_FUNCTION);
    }
    if (types.runInt128) {
        BenchOneType<int128_t>(TIME_PER_FUNCTION);
    }
    if (types.runFixedPoint128) {
        BenchOneType<fixed_point128<BENCH_FP_INT_BITS>>(TIME_PER_FUNCTION);
    }
    if (types.runFloat128) {
        BenchOneType<float128>(TIME_PER_FUNCTION);
    }
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
    printf("Usage: %s [-j|--json] [-t|--type <type>] [-h|--help]\n", exeName);
    printf("  -j, --json         Write the results to a JSON file in addition to the terminal.\n");
    printf("                     The file is named after the compiler and build type, e.g. %s\n", JsonFileName().c_str());
    printf("  -t, --type <type>  Measure only the given type. May be repeated; all types are\n");
    printf("                     measured when the option is absent. One of:\n");
    printf("                       uint128 (uint128_t)          128 bit unsigned integer\n");
    printf("                       int128 (int128_t)            128 bit signed integer\n");
    printf("                       fixed_point128 (fp128)       fixed point, %d integer bits\n", BENCH_FP_INT_BITS);
    printf("                       float128 (f128)              quadruple precision float\n");
    printf("                       all                          every type above\n");
    printf("  -h, --help         Print this help text and exit.\n");
}

int main(int argc, char* argv[])
{
    bool jsonOutput = false;
    TypeSelection types;

    for (int i = 1; i < argc; ++i) {
        if (std::strcmp(argv[i], "-j") == 0 || std::strcmp(argv[i], "--json") == 0) {
            jsonOutput = true;
        } else if (std::strcmp(argv[i], "-t") == 0 || std::strcmp(argv[i], "--type") == 0) {
            if (++i == argc) {
                printf("Missing type name after %s\n", argv[i - 1]);
                PrintUsage(argv[0]);
                return 1;
            }
            if (!SelectType(argv[i], types)) {
                printf("Unknown type: %s\n", argv[i]);
                PrintUsage(argv[0]);
                return 1;
            }
        } else if (std::strcmp(argv[i], "-h") == 0 || std::strcmp(argv[i], "--help") == 0) {
            PrintUsage(argv[0]);
            return 0;
        } else {
            printf("Unknown option: %s\n", argv[i]);
            PrintUsage(argv[0]);
            return 1;
        }
    }

    // no explicit selection means measure everything
    if (!types.any()) {
        types.selectAll();
    }

    // Force instantiation of all public methods and friend functions for I=1, 40 and 64.
    force_instantiation<1>();
    force_instantiation<40>();
    force_instantiation<64>();

    bench(types);

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
