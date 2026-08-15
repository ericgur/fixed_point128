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
#include <atomic>
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

/** @brief Wall clock time each measured operation is given, in seconds. */
constexpr double TIME_PER_FUNCTION = 0.5;

/**
 * @brief Shortest a timed batch is allowed to be, in seconds.
 *
 * Every batch is bracketed by two clock reads, and a clock read costs on the order of 20ns here. At
 * two milliseconds a batch that is the only thing standing between them, they contribute about one
 * part in a hundred thousand - far below the run to run scatter - which is what lets the harness
 * measure a one cycle operation without the timer showing up in the result.
 */
constexpr double BENCH_MIN_BATCH_TIME = 0.002;

/** @brief Iteration count the batch sizing starts from. */
constexpr uint64_t BENCH_MIN_BATCH_ITERATIONS = 256;

/**
 * @brief Largest batch the sizing loop will produce.
 *
 * Only reached by a loop that costs nothing per iteration, which in practice means one the optimizer
 * managed to delete. The cap keeps such a benchmark from running away; the absurd rate it then
 * reports is the signal that the loop needs looking at.
 */
constexpr uint64_t BENCH_MAX_BATCH_ITERATIONS = 1ull << 30;

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
 * @brief Volatile sink that swallows the addresses handed to Escape().
 *
 * Storing a pointer here is an observable side effect, so the address genuinely leaves the function
 * and the optimizer has to assume that whatever it points at is reachable from outside.
 */
static void* volatile escapeSink = nullptr;

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
 * loop, never inside it. Escape() and Barrier() are what a timed loop uses.
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
 * @brief Publishes the address of @p object, so the optimizer must treat it as externally reachable.
 *
 * On its own this only says that the object could be read or written from elsewhere; it does not by
 * itself stop anything being hoisted. Paired with Barrier() inside the loop it does, because the
 * barrier is a point at which that elsewhere could have run. The two are always used together:
 * Escape() once during setup, Barrier() once per iteration.
 *
 * Escaping an object costs it its registers - it has to live in memory - but nothing per iteration
 * beyond the loads and stores that implies, and on x86 those fold into the instructions that use
 * them. This is the whole reason the harness is built on it rather than on an opaque call.
 *
 * @tparam T Type of the object.
 * @param object Object whose address is published. Must outlive the loop that relies on it.
 */
template <typename T> FP128_INLINE void Escape(T& object) noexcept
{
    escapeSink = static_cast<void*>(&object);
}

/**
 * @brief Compiler only memory barrier: emits no instructions at all.
 *
 * Placed at the top of a timed loop it does two jobs at once, both of which the loops here need:
 * - An escaped operand read after it cannot have been read before it, so an expression computed
 *   from that operand cannot be hoisted out of the loop, nor rewritten in closed form. Without this
 *   Clang 17 turns the 128 bit addition benchmark into an AVX2 computation of the sum of an
 *   arithmetic progression and reports over 10G/s - a rate no chain of 64 bit adds can reach.
 * - A store to an escaped result made before it cannot be dropped in favour of the next iteration's
 *   store, so the operation that produced the result stays alive without the loop having to feed it
 *   back into itself.
 *
 * std::atomic_signal_fence() is the portable spelling: MSVC expands it to _ReadWriteBarrier() and
 * Clang and GCC to a single thread fence, none of which generate code.
 */
FP128_INLINE void Barrier() noexcept
{
    std::atomic_signal_fence(std::memory_order_acq_rel);
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
 * @brief Runs @p body in timed batches and returns the fastest per iteration rate observed.
 *
 * @p body runs a whole batch: it takes the iteration count, sets up its own operands, loops, and
 * consumes whatever it produced. Owning its state is what lets its accumulators stay in registers -
 * state passed in from here would have to be addressable, and the loop would then measure the
 * memory traffic instead of the operation. It also bounds how far a benchmark that accumulates into
 * itself can drift, since every batch starts from the same operands. bench_mandelbrot() is the one
 * benchmark that deliberately does the opposite; the reason is on that function.
 *
 * Two things make this more accurate than timing one long run:
 * - The batch is sized so it lasts at least BENCH_MIN_BATCH_TIME, which pushes the cost of the two
 *   clock reads bracketing it below the noise floor however cheap a single iteration is.
 * - The *fastest* batch is reported rather than the average over the whole budget. Every source of
 *   error left - a scheduler tick, a migration to another core, a frequency excursion - can only
 *   make a batch slower, so the fastest one is the sample least contaminated by them, and it is far
 *   more repeatable from run to run than the mean. Measured over three consecutive runs, the spread
 *   of a result is typically under half a percent.
 *
 * @tparam Body Callable invoked as body(uint64_t iterations).
 * @param time_budget Time to spend measuring, in seconds. Excludes the sizing phase.
 * @param body Batch to run.
 * @return Iterations per second of the fastest batch, or 0 if no batch could be timed.
 */
template <typename Body> [[nodiscard]] int64_t MeasureRate(double time_budget, Body body)
{
    Duration dur;

    // Sizing phase. It doubles as the warm up: by the time the last batch runs, the caches and the
    // branch predictors have seen the loop and the core has had time to clock up.
    uint64_t batch = BENCH_MIN_BATCH_ITERATIONS;
    double elapsed = 0.0;
    for (;;) {
        dur.start();
        body(batch);
        elapsed = dur.cur_duration();
        if (elapsed >= BENCH_MIN_BATCH_TIME || batch >= BENCH_MAX_BATCH_ITERATIONS) {
            break;
        }

        // Aim at twice the target, so a batch that already nearly reaches it does not need several
        // more rounds. The step is clamped because the first batches are too short to extrapolate
        // from: below, so a batch delayed by an interrupt cannot stall the search, and above, so a
        // clock that reported zero cannot overshoot the cap in one go.
        const double growth = (elapsed > 0.0) ? (2.0 * BENCH_MIN_BATCH_TIME / elapsed) : 100.0;
        batch = (uint64_t)((double)batch * std::clamp(growth, 2.0, 100.0));
        batch = std::min(batch, BENCH_MAX_BATCH_ITERATIONS);
    }

    // Measurement phase.
    double best_seconds = elapsed / (double)batch;
    Duration budget;
    budget.start();
    while (budget.cur_duration() < time_budget) {
        dur.start();
        body(batch);
        const double per_iteration = dur.cur_duration() / (double)batch;
        best_seconds = std::min(best_seconds, per_iteration);
    }

    return (best_seconds > 0.0) ? (int64_t)(1.0 / best_seconds) : 0;
}

/**
 * @brief Number of distinct arguments a timed loop cycles through. Must be a power of two.
 *
 * Sixty four arguments of at most 24 bytes stay inside L1, so cycling through them costs an L1 load
 * per iteration and nothing else.
 */
constexpr uint64_t BENCH_ARG_COUNT = 64;

/**
 * @brief Fills @p args with values that differ in their low order bits but not in their magnitude.
 *
 * Calling a function on one argument for half a second measures something no caller ever sees. The
 * branch predictor learns the whole sequence of data dependent branches inside the function and
 * then gets it right every time: fixed_point128<10>::log2() runs 118 iterations with one such
 * branch each, and timing it on a single argument overstates it by about 20%.
 *
 * Worse, a fixed argument can invert a comparison rather than just shift it. Making the rounding
 * step inside square() branchless measures 10% *slower* on a fixed argument and 16% faster on
 * varying ones, on the same compiler - so the benchmark would have rejected a change that is worth
 * having.
 *
 * The perturbation is deliberately small, at most a part in 2^10. The arguments have to stay inside
 * the domain their caller picked them for - asin() needs |x| <= 1, acosh() needs x >= 1 - and they
 * have to keep the magnitude they were chosen with, since the series based functions iterate a
 * number of times that depends on it and the results would otherwise stop being comparable with
 * earlier runs.
 *
 * @tparam T Argument type.
 * @param args Array of BENCH_ARG_COUNT elements to fill.
 * @param base Argument the benchmark would otherwise have used on its own. Becomes args[0].
 */
template <typename T> void BuildArgs(T* args, const T& base)
{
    const T step = base >> 16;
    for (uint64_t i = 0; i < BENCH_ARG_COUNT; ++i) {
        args[i] = base + step * static_cast<uint32_t>(i);
    }
}

/**
 * @brief Times a unary function over a rotating set of operands.
 *
 * The operands are escaped and a barrier runs at the top of every iteration, so the call cannot be
 * hoisted out of the loop; the result is escaped too, so the store the loop makes to it on every
 * iteration keeps the call from being deleted as dead. Operands and result therefore live in
 * memory, which on x86 costs the loads and stores only - they fold into the instructions around
 * them - and nothing at all in extra instructions.
 *
 * See BuildArgs() for why the loop cycles through a set of arguments rather than reusing one.
 *
 * Assigning the returned result is what BenchBinary() had to stop doing, because MSVC widens the
 * copy of a trivially copyable 128 bit return value to 16 bytes and then stalls on the store
 * forwarding. Nothing measured here reaches that: every unary function the two integer types have -
 * sqrt() and the log() family - returns a uint64_t, which is copied with a single store, and the
 * two fractional types assign their QWORDs one at a time and so never see the wide copy. A unary
 * function returning uint128_t or int128_t would reintroduce it, and would need the treatment
 * BenchBinary() documents.
 *
 * @tparam T Operand type.
 * @tparam Func Callable invoked as func(const T&).
 * @param name Name to print the measurement under.
 * @param time_per_function Time to spend measuring, in seconds.
 * @param argument Value the function is applied to, and the base the rest are derived from.
 * @param func Function to measure.
 */
template <typename T, typename Func> void BenchUnary(const char* name, double time_per_function, T argument, Func func)
{
    const int64_t ips = MeasureRate(time_per_function, [argument, func](uint64_t count) {
        T args[BENCH_ARG_COUNT];
        BuildArgs(args, argument);
        auto result = func(args[0]);
        Escape(args[0]);
        Escape(result);
        for (uint64_t i = count; i != 0; --i) {
            Barrier();
            result = func(args[i & (BENCH_ARG_COUNT - 1)]);
        }
        DoNotOptimize(result);
    });

    print_ips(name, ips);
}

/**
 * @brief Times a binary operation over a rotating set of left hand operands.
 *
 * Both operands are escaped, not just the one that would be enough to defeat hoisting: leaving the
 * right hand side visible would let the compiler specialize the operation for it, and a division by
 * a literal 5 rewritten as a multiplication by its reciprocal is not the division this is meant to
 * be timing.
 *
 * Only the left hand side rotates. That is enough to keep the operation's data dependent branches
 * from repeating - see BuildArgs() - and it leaves the right hand side free to be a type that has
 * no arithmetic of its own, such as the uint32_t exponent of the integer pow().
 *
 * @p op has to apply the operation to the left hand operand in place rather than return the result,
 * for the reason spelled out on BenchAccumulate(): written as `result = lhs op rhs`, MSVC builds the
 * operator's return value in a stack temporary with two QWORD stores and then copies it into the
 * escaped result with a single 16 byte load. That load overlaps both stores and cannot be forwarded,
 * so it waits for them to reach L1 - about 15 cycles, on every iteration.
 *
 * The stall only lands on the types whose copy assignment is the compiler generated one, because
 * that is the copy MSVC is free to widen to 16 bytes; fixed_point128 and float128 assign their two
 * QWORDs individually and never see it. Timing `lhs op rhs` therefore charged uint128_t and int128_t
 * for a stall the fractional types were not paying, which put fixed_point128 ahead of uint128_t on
 * the 128 bit multiplication - 491M/s against 258M/s, for an operation that is three multiplies
 * against four plus a shift and a rounding step. In place, the two come out at 494M/s and 1.73G/s,
 * and clang-cl - which never emitted the wide copy and so always measured the multiply itself -
 * agrees with both figures.
 *
 * Nothing is lost by measuring the compound assignment: the binary operator is defined as `lhs op=
 * rhs` on a copy of the left operand, and the copy is exactly what the loop makes when it reloads
 * the next operand into the result.
 *
 * @tparam T Type of the left hand operand.
 * @tparam U Type of the right hand operand.
 * @tparam Op Callable invoked as op(T& lhs, const U& rhs), which must update lhs in place.
 * @param name Name to print the measurement under.
 * @param time_per_function Time to spend measuring, in seconds.
 * @param left Left hand operand, and the base the rest of the rotating set is derived from.
 * @param right Right hand operand, the same on every iteration.
 * @param op Operation to measure.
 */
template <typename T, typename U, typename Op> void BenchBinary(const char* name, double time_per_function, T left, U right, Op op)
{
    const int64_t ips = MeasureRate(time_per_function, [left, right, op](uint64_t count) {
        T args[BENCH_ARG_COUNT];
        BuildArgs(args, left);
        U rhs = right;
        T result = args[0];
        Escape(args[0]);
        Escape(rhs);
        Escape(result);
        for (uint64_t i = count; i != 0; --i) {
            Barrier();
            result = args[i & (BENCH_ARG_COUNT - 1)];
            op(result, rhs);
        }
        DoNotOptimize(result);
    });

    print_ips(name, ips);
}

/**
 * @brief Times an operation that folds the right hand side into an accumulator, in place.
 *
 * This measures the latency of the operation rather than its throughput, which is the honest thing
 * to report for operators cheap enough that a dependent chain of them is what limits real code.
 * Only the right hand side is escaped: the accumulator stays in registers, so the chain measured is
 * the operation itself and not a round trip through the stack.
 *
 * @p func has to update the accumulator in place rather than return the new value. Written the
 * other way round - acc = acc + rhs - MSVC materializes the operator's return value in a stack
 * temporary and reads the accumulator back out of it on every iteration. The 16 byte store and the
 * 8 byte reload that overlaps it do not forward, and the resulting stall costs a factor of 25 on
 * the uint128_t addition benchmark. Nothing is lost by measuring the compound assignment instead:
 * operator+ is defined as `lhs += rhs` on a copy, so this is the same operation with a copy removed
 * that the measurement has no business including.
 *
 * @tparam T Operand type.
 * @tparam Func Callable invoked as func(T& acc, const T& rhs), which must update acc in place.
 * @param name Name to print the measurement under.
 * @param time_per_function Time to spend measuring, in seconds.
 * @param seed Initial value of the accumulator.
 * @param operand Right hand operand, the same on every iteration.
 * @param func Operation to measure.
 */
template <typename T, typename Func> void BenchAccumulate(const char* name, double time_per_function, T seed, T operand, Func func)
{
    const int64_t ips = MeasureRate(time_per_function, [seed, operand, func](uint64_t count) {
        T rhs = operand;
        T acc = seed;
        Escape(rhs);
        for (uint64_t i = count; i != 0; --i) {
            Barrier();
            func(acc, rhs);
        }
        DoNotOptimize(acc);
    });

    print_ips(name, ips);
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

    // The four comparisons accumulate into a plain integer. That is enough to keep them alive and,
    // unlike escaping their results, it leaves the accumulator in a register.
    //
    // The right hand side is taken from the middle of the rotating set rather than from
    // operandC(). Rotating the left hand side alone would not have helped: every value BuildArgs()
    // produces is within a part in 2^10 of the base, so all of them would still land on the same
    // side of an unrelated constant and every comparison would take the branch it took last time.
    // Sitting the right hand side inside the set splits the outcomes evenly, which is what a
    // comparison in real code does.
    const int64_t ips = MeasureRate(time_per_function, [](uint64_t count) {
        T args[BENCH_ARG_COUNT];
        BuildArgs(args, Traits::operandB());
        T f2 = args[BENCH_ARG_COUNT / 2];
        int64_t matches = 0;
        Escape(args[0]);
        Escape(f2);
        for (uint64_t i = count; i != 0; --i) {
            Barrier();
            const T& f1 = args[i & (BENCH_ARG_COUNT - 1)];
            matches += (f1 > f2);
            matches += (f1 >= f2);
            matches += (f1 < f2);
            matches += (f1 <= f2);
        }
        DoNotOptimize(matches);
    });

    // One iteration is four comparisons, and the reported figure is the average cost of one.
    print_ips("Operators >, >=, <, <= (average of all 4)", 4 * ips);
}

template <typename T> void bench_addition(double time_per_function = 1.0)
{
    using Traits = BenchTraits<T>;
    BenchAccumulate<T>("Addition", time_per_function, Traits::operandB(), Traits::operandA(), [](T& acc, const T& rhs) { acc += rhs; });
}

template <typename T> void bench_subtraction(double time_per_function = 1.0)
{
    using Traits = BenchTraits<T>;
    BenchAccumulate<T>("Subtraction", time_per_function, Traits::operandA(), Traits::operandB(), [](T& acc, const T& rhs) { acc -= rhs; });
}

template <typename T> void bench_multiplication(double time_per_function = 1.0)
{
    using Traits = BenchTraits<T>;
    using MulType = typename Traits::MulType;
    using IntMulType = typename Traits::IntMulType;

    BenchBinary<MulType, MulType>("Multiplication by 128-bit value", time_per_function, Traits::mulA(), Traits::mulB(),
                                  [](MulType& lhs, const MulType& rhs) { lhs *= rhs; });

    // The result is escaped rather than accumulated. The accumulating form (f10 = f10 * int_val) was
    // degenerate: the low QWORD of the product does not depend on the high QWORD, so with nothing
    // observing the result both compilers proved the high half dead. MSVC collapsed the whole
    // operation into a single 64 bit mulx and Clang went further, vectorizing the remaining chain
    // 4 wide - neither was timing a 128 bit multiply.
    BenchBinary<IntMulType, uint32_t>("Multiplication by int32_t", time_per_function, Traits::intMulA(), 123456789u,
                                      [](IntMulType& lhs, const uint32_t& rhs) { lhs *= rhs; });
}

template <typename T> void bench_division(double time_per_function = 1.0)
{
    using Traits = BenchTraits<T>;

    BenchBinary<T, double>("Division by double (exponent of 2)", time_per_function, Traits::operandA(), 64.0,
                           [](T& lhs, const double& rhs) { lhs /= rhs; });

    BenchBinary<T, int64_t>("Division by int64", time_per_function, Traits::operandA(), 5ll, [](T& lhs, const int64_t& rhs) { lhs /= rhs; });

    BenchBinary<T, T>("Division by 128-bit integer value", time_per_function, Traits::operandA(), T(5), [](T& lhs, const T& rhs) { lhs /= rhs; });

    // Only meaningful where the divisor can hold a fraction. For the integer types this loop would
    // be the previous one with a different divisor.
    if constexpr (Traits::isFractional) {
        BenchBinary<T, T>("Division by 128-bit fractional value", time_per_function, Traits::operandA(), Traits::operandB(),
                          [](T& lhs, const T& rhs) { lhs /= rhs; });
    }
}

template <typename T> void bench_reciprocal(double time_per_function = 1.0)
{
    BenchUnary<T>("reciprocal", time_per_function, BenchTraits<T>::operandB(), [](const T& v) { return reciprocal(v); });
}

/**
 * @brief Benches sqrt().
 *
 * BenchUnary() holds the result in `auto`, which matters here: the integer types return the
 * truncated integer root as a uint64_t rather than a value of their own type.
 */
template <typename T> void bench_sqrt(double time_per_function = 1.0)
{
    BenchUnary<T>("sqrt", time_per_function, BenchTraits<T>::operandB(), [](const T& v) { return sqrt(v); });
}

template <typename T> void bench_exp(double time_per_function = 1.0)
{
    BenchUnary<T>("exp", time_per_function, BenchTraits<T>::operandB(), [](const T& v) { return exp(v); });
}

template <typename T> void bench_exp2(double time_per_function = 1.0)
{
    BenchUnary<T>("exp2", time_per_function, BenchTraits<T>::operandB(), [](const T& v) { return exp2(v); });
}

template <typename T> void bench_expm1(double time_per_function = 1.0)
{
    BenchUnary<T>("expm1", time_per_function, BenchTraits<T>::operandB(), [](const T& v) { return expm1(v); });
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

    if constexpr (Traits::hasTranscendental) {
        BenchBinary<T, T>("pow", time_per_function, Traits::operandB(), Traits::operandC(),
                          [](T& base, const T& exponent) { base = pow(base, exponent); });
    } else {
        // A small base keeps the result inside 128 bits; the cost of the binary exponentiation
        // depends on the exponent, not on the base.
        BenchBinary<T, uint32_t>("pow (integer exponent)", time_per_function, T(7), 5u,
                                 [](T& base, const uint32_t& exponent) { base = pow(base, exponent); });
    }
}

/**
 * @brief Benches log().
 *
 * As with sqrt(), the integer types return a uint64_t: the integer part of the logarithm.
 */
template <typename T> void bench_log(double time_per_function = 1.0)
{
    BenchUnary<T>("log", time_per_function, BenchTraits<T>::operandB(), [](const T& v) { return log(v); });
}

template <typename T> void bench_log2(double time_per_function = 1.0)
{
    BenchUnary<T>("log2", time_per_function, BenchTraits<T>::operandB(), [](const T& v) { return log2(v); });
}

template <typename T> void bench_log10(double time_per_function = 1.0)
{
    BenchUnary<T>("log10", time_per_function, BenchTraits<T>::operandB(), [](const T& v) { return log10(v); });
}

template <typename T> void bench_log1p(double time_per_function = 1.0)
{
    BenchUnary<T>("log1p", time_per_function, BenchTraits<T>::operandB(), [](const T& v) { return log1p(v); });
}

template <typename T> void bench_sin(double time_per_function = 1.0)
{
    BenchUnary<T>("sin", time_per_function, BenchTraits<T>::operandB() / 2, [](const T& v) { return sin(v); });
}

template <typename T> void bench_asin(double time_per_function = 1.0)
{
    BenchUnary<T>("asin", time_per_function, BenchTraits<T>::operandA() / 5, [](const T& v) { return asin(v); });
}

template <typename T> void bench_cos(double time_per_function = 1.0)
{
    BenchUnary<T>("cos", time_per_function, BenchTraits<T>::operandB() / 2, [](const T& v) { return cos(v); });
}

template <typename T> void bench_acos(double time_per_function = 1.0)
{
    BenchUnary<T>("acos", time_per_function, BenchTraits<T>::operandA() / 5, [](const T& v) { return acos(v); });
}

template <typename T> void bench_tan(double time_per_function = 1.0)
{
    BenchUnary<T>("tan", time_per_function, BenchTraits<T>::operandB() / 2, [](const T& v) { return tan(v); });
}

template <typename T> void bench_atan(double time_per_function = 1.0)
{
    BenchUnary<T>("atan", time_per_function, BenchTraits<T>::operandA() / 5, [](const T& v) { return atan(v); });
}

template <typename T> void bench_sinh(double time_per_function = 1.0)
{
    BenchUnary<T>("sinh", time_per_function, BenchTraits<T>::operandB() / 2, [](const T& v) { return sinh(v); });
}

template <typename T> void bench_asinh(double time_per_function = 1.0)
{
    BenchUnary<T>("asinh", time_per_function, BenchTraits<T>::operandA() / 5, [](const T& v) { return asinh(v); });
}

template <typename T> void bench_cosh(double time_per_function = 1.0)
{
    BenchUnary<T>("cosh", time_per_function, BenchTraits<T>::operandB() / 2, [](const T& v) { return cosh(v); });
}

template <typename T> void bench_acosh(double time_per_function = 1.0)
{
    // x >= 1
    BenchUnary<T>("acosh", time_per_function, BenchTraits<T>::operandB() / 2, [](const T& v) { return acosh(v); });
}

template <typename T> void bench_tanh(double time_per_function = 1.0)
{
    BenchUnary<T>("tanh", time_per_function, BenchTraits<T>::operandB() / 2, [](const T& v) { return tanh(v); });
}

template <typename T> void bench_atanh(double time_per_function = 1.0)
{
    // abs(v) < 1
    BenchUnary<T>("atanh", time_per_function, BenchTraits<T>::operandB() / 4, [](const T& v) { return atanh(v); });
}

/**
 * @brief Benches one Mandelbrot iteration.
 *
 * Unlike every other benchmark here, the state lives outside the batch and the orbit runs on across
 * batch boundaries. It has to: the seed coordinates come from doubles, so the low two thirds of
 * their fraction bits are zero, and it takes on the order of a million iterations for the products
 * to fill them in. MSVC is roughly twice as slow over that transient as it is afterwards, so a
 * benchmark that restarted the orbit for every batch would report a figure that depended on the
 * batch size the sizing phase happened to pick - 68M/s at one batch length against 152M/s at
 * another, for the same loop. Running the orbit on makes the transient a one time cost.
 *
 * The point is inside the set, so the orbit stays bounded for as long as it is iterated and no
 * value ever leaves the range the type can hold.
 */
template <typename T> void bench_mandelbrot(double time_per_function = 1.0)
{
    // Both coordinates are escaped so the iteration cannot be constant folded; the orbit itself is
    // loop carried, so nothing else here is hoistable and no result needs escaping to keep the body
    // alive.
    T x = T(-0.7294734415);
    T y = T(0.242809);
    T usq, vsq, tmp, modulus, u, v;
    Escape(x);
    Escape(y);

    const int64_t ips = MeasureRate(time_per_function, [&](uint64_t count) {
        for (uint64_t i = count; i != 0; --i) {
            Barrier();
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
        DoNotOptimize(modulus);
    });

    print_ips("Mandelbrot", ips);
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
 * Uses the same condition as the library itself (fp128_shared.h), so the reported build
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
        case '"':
            res += "\\\"";
            break;
        case '\\':
            res += "\\\\";
            break;
        case '\b':
            res += "\\b";
            break;
        case '\f':
            res += "\\f";
            break;
        case '\n':
            res += "\\n";
            break;
        case '\r':
            res += "\\r";
            break;
        case '\t':
            res += "\\t";
            break;
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
 * The report holds the build metadata needed to compare two runs meaningfully - the library version,
 * compiler, build type and the timing parameters - followed by the results in execution order, each
 * tagged with the type it was measured on and the group it was printed under.
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
    file << format("  \"libraryVersion\": \"{}\",\n", version_string);
    file << format("  \"compiler\": \"{}\",\n", EscapeJson(CompilerName()));
    file << format("  \"compilerTag\": \"{}\",\n", CompilerTag());
    file << format("  \"build\": \"{}\",\n", BuildType());
    file << format("  \"timestamp\": \"{}\",\n", Timestamp());
    file << format("  \"timePerFunction\": {},\n", TIME_PER_FUNCTION);
    file << format("  \"minBatchTime\": {},\n", BENCH_MIN_BATCH_TIME);
    file << "  \"sampling\": \"fastest batch\",\n";
    file << format("  \"types\": [{}],\n", typeList);
    file << "  \"results\": [\n";
    for (size_t i = 0; i < benchResults.size(); ++i) {
        const auto& result = benchResults[i];
        const char* separator = (i + 1 < benchResults.size()) ? "," : "";
        file << format("    {{ \"type\": \"{}\", \"group\": \"{}\", \"name\": \"{}\", \"iterationsPerSecond\": {} }}{}\n", EscapeJson(result.type),
                       EscapeJson(result.group), EscapeJson(result.name), result.ips, separator);
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
    printf("fixed_point128 %s\n", version_string);
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
    fp def;                            // default
    fp from_double(1.5);               // double
    fp copy_ctor(from_double);         // copy
    fp move_ctor(std::move(fp(2.0)));  // move
    fp from_u64((uint64_t)1);          // uint64_t
    fp from_i64((int64_t)1);           // int64_t
    fp from_u32((uint32_t)1);          // uint32_t
    fp from_i32((int32_t)1);           // int32_t
    fp from_cstr("1.5");               // const char*
    fp from_str(std::string("1.5"));   // std::string
    fp from_raw(0ull, 1ull);           // raw (low, high)

    // cross-template copy constructor (I2 = 10)
    fixed_point128<10> f10(1.5);
    fp cross_ctor(f10);

    // --- Assignment operators ---
    def = copy_ctor;           // copy assign
    def = std::move(fp(3.0));  // move assign
    def = f10;                 // cross-template assign

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
    // The operands are kept small enough that nothing below leaves the range of any instantiation
    // this is called with, fixed_point128<1> included, whose range stops at 2. An overflowed value
    // lands wherever the wrap takes it, possibly negative, and log() rejects a non positive
    // argument by throwing - which is a crash rather than a measurement.
    fp a = fp::half();
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
    c *= (uint64_t)2;  // operator*=<uint64_t> specialization
    c /= (uint64_t)2;  // operator/=<uint64_t> specialization
    c /= (double)2.0;  // operator/=<double> specialization

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
    // pi, pi2 and e do not fit in the smallest instantiations and say so with a static_assert
    if constexpr (I >= 2) {
        (void)fp::pi();
        (void)fp::e();
    }
    if constexpr (I >= 3) {
        (void)fp::pi2();
    }
    (void)fp::half_pi();
    (void)fp::golden_ratio();
    (void)fp::sqrt_2();
    (void)fp::one();
    (void)fp::half();
    (void)fp::epsilon();

    // --- Friend math functions (CRT-style) ---
    // below one, so log1p() and exp() stay in range as well - see the note at 'a' above
    fp val = fp::half();
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
    // the exponential family multiplies by e(), which needs 2 integer bits
    if constexpr (I >= 2) {
        (void)exp(val);
        (void)exp2(val);
        (void)expm1(val);
        (void)pow(val, val);
    }
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
        (void)asin(half_val);  // |x| <= 1 required
        (void)cos(val);
        (void)acos(half_val);  // |x| <= 1 required
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
    (void)def;
    (void)from_double;
    (void)copy_ctor;
    (void)move_ctor;
    (void)from_u64;
    (void)from_i64;
    (void)from_u32;
    (void)from_i32;
    (void)from_cstr;
    (void)from_str;
    (void)from_raw;
    (void)cross_ctor;
    (void)a;
    (void)b;
    (void)c;
    (void)f10;
    (void)fact_res;
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

    // Force instantiation of all public methods and friend functions for I=1, 40 and 63.
    force_instantiation<1>();
    force_instantiation<40>();
    force_instantiation<63>();

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
