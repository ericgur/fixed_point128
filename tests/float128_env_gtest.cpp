// remove warnings from gtest itself
#if defined(_MSC_VER)
#pragma warning(push)
#pragma warning(disable : 26439)
#pragma warning(disable : 26495)
#endif
#include <gtest/gtest.h>
#if defined(_MSC_VER)
#pragma warning(pop)
#endif
#include <bit>
#include <format>
#include <thread>
#include "float128_ref_check.h"

/**********************************************************************
 * float128 with FP128_IEEE_ENV: rounding directions and exception flags
 *
 * This file is the whole of its own test program, built with FP128_IEEE_ENV defined. The macro
 * changes the definitions of inline functions, so it cannot share an executable with the rest of
 * the suite, which is built without it.
 *
 * The correctly rounded operations are run over their exact reference tables here as well: the
 * environment must not change what rounding to nearest produces.
 ***********************************************************************/

#ifndef FP128_IEEE_ENV
#error "float128_env_gtest.cpp has to be built with FP128_IEEE_ENV defined"
#endif

using namespace f128_ref;

namespace
{
/// @brief A value from its encoding, high QWORD first.
[[nodiscard]] float128 Encoded(uint64_t high, uint64_t low)
{
    return float128(low, high);
}

[[nodiscard]] float128 Pow2(int e)
{
    return ldexp(float128::one(), e);
}

#define EXPECT_SAME(actual, expected) \
    EXPECT_TRUE(SameResult((actual), (expected))) << #actual << ": expected " << Bits(expected) << ", got " << Bits(actual)

/// @brief Selects a rounding direction for the length of a scope and puts nearest back after.
class RoundingScope
{
public:
    explicit RoundingScope(int mode) { EXPECT_EQ(fp128::fesetround(mode), 0); }
    ~RoundingScope() { fp128::fesetround(FE_TONEAREST); }
    RoundingScope(const RoundingScope&) = delete;
    RoundingScope& operator=(const RoundingScope&) = delete;
};

/// @brief The flags an expression raised, starting from none.
template <typename Fn> [[nodiscard]] int FlagsOf(Fn fn)
{
    fp128::feclearexcept(FE_ALL_EXCEPT);
    // Not const: a const integral local initialized from a constant expression is evaluated at
    // compile time, where no flag can be raised.
    [[maybe_unused]] auto result = fn();
    const int flags = fp128::fetestexcept(FE_ALL_EXCEPT);
    fp128::feclearexcept(FE_ALL_EXCEPT);
    return flags;
}

constexpr int modes[] = {FE_TONEAREST, FE_TOWARDZERO, FE_UPWARD, FE_DOWNWARD};
}  // namespace

TEST(float128_env, Conformance)
{
    static_assert(std::numeric_limits<float128>::is_iec559);
    static_assert(float128::is754version2008());
    EXPECT_EQ(fp128::fegetround(), FE_TONEAREST);
    for (const int mode : modes) {
        EXPECT_EQ(fp128::fesetround(mode), 0);
        EXPECT_EQ(fp128::fegetround(), mode);
    }
    EXPECT_NE(fp128::fesetround(-12345), 0);
    EXPECT_EQ(fp128::fegetround(), FE_DOWNWARD);
    fp128::fesetround(FE_TONEAREST);
}

TEST(float128_env, NearestIsUnchanged)
{
    CheckExactBinary("add", add_exact_ref, [](const float128& x, const float128& y) { return x + y; });
    CheckExactBinary("mul", mul_exact_ref, [](const float128& x, const float128& y) { return x * y; });
    CheckExactBinary("div", div_exact_ref, [](const float128& x, const float128& y) { return x / y; });
    CheckExactUnary("sqrt", sqrt_exact_ref, [](const float128& x) { return sqrt(x); });
    CheckExactTernary("fma", fma_exact_ref, [](const float128& x, const float128& y, const float128& z) { return fma(x, y, z); });
    CheckExactUnary("to_double", to_double_exact_ref, [](const float128& x) { return float128(static_cast<double>(x)); });
}

TEST(float128_env, DirectedRoundingOfTheArithmetic)
{
    // 1/3 is 0x3FFD5555...5555 and a third of an ulp more; -1/3 the mirror image.
    const float128 down = Encoded(0x3ffd555555555555, 0x5555555555555555);
    const float128 up = float128::nextUp(down);
    const float128 third_nearest = float128(1) / float128(3);
    EXPECT_SAME(third_nearest, down);
    {
        RoundingScope scope(FE_UPWARD);
        EXPECT_SAME(float128(1) / float128(3), up);
        EXPECT_SAME(float128(-1) / float128(3), -down);
        EXPECT_SAME(sqrt(float128(2)), float128::nextUp(float128::sqrt_2()));
        EXPECT_SAME(float128::one() + Pow2(-200), float128::nextUp(float128::one()));
    }
    {
        RoundingScope scope(FE_DOWNWARD);
        EXPECT_SAME(float128(1) / float128(3), down);
        EXPECT_SAME(float128(-1) / float128(3), -up);
        EXPECT_SAME(sqrt(float128(2)), float128::sqrt_2());
        EXPECT_SAME(float128::one() - Pow2(-200), float128::nextDown(float128::one()));
    }
    {
        RoundingScope scope(FE_TOWARDZERO);
        EXPECT_SAME(float128(-1) / float128(3), -down);
        EXPECT_SAME(float128(2) / float128(3), Encoded(0x3ffe555555555555, 0x5555555555555555));
        // 3 * (1/3 rounded down) is exactly 1 - 2^-114, which fma reaches without a rounding
        EXPECT_SAME(fma(float128(3), float128(1) / float128(3), -float128::one()), -Pow2(-114));
    }
}

TEST(float128_env, ExactZeroSign)
{
    const float128 x(1.5);
    for (const int mode : modes) {
        RoundingScope scope(mode);
        const uint32_t want = (mode == FE_DOWNWARD) ? 1u : 0u;
        EXPECT_EQ((x - x).get_sign(), want) << mode;
        EXPECT_EQ((x + (-x)).get_sign(), want) << mode;
        EXPECT_EQ((float128() + (-float128())).get_sign(), want) << mode;
        EXPECT_EQ(fma(x, float128(2), float128(-3)).get_sign(), want) << mode;
        // like signs keep their sign whatever the direction
        EXPECT_EQ((-float128() + -float128()).get_sign(), 1u) << mode;
    }
}

TEST(float128_env, OverflowDependsOnTheDirection)
{
    const float128 max = std::numeric_limits<float128>::max();
    const float128 inf = float128::inf();
    // the aligned members first, so that no padding is inserted between them and the int (C4324)
    const struct {
        float128 positive, negative;
        int mode;
    } cases[] = {{inf, -inf, FE_TONEAREST}, {max, -max, FE_TOWARDZERO}, {inf, -max, FE_UPWARD}, {max, -inf, FE_DOWNWARD}};
    for (const auto& c : cases) {
        RoundingScope scope(c.mode);
        EXPECT_SAME(max * float128(2), c.positive);
        EXPECT_SAME(-max * float128(2), c.negative);
        EXPECT_SAME(max + max, c.positive);
        EXPECT_SAME(float128("1e5000"), c.positive);
        EXPECT_SAME(ldexp(max, 10), c.positive);
    }
}

TEST(float128_env, UnderflowDependsOnTheDirection)
{
    const float128 tiny = std::numeric_limits<float128>::denorm_min();
    {
        RoundingScope scope(FE_TONEAREST);
        EXPECT_SAME(tiny / float128(2), float128());  // a tie, and zero is even
        EXPECT_SAME(tiny * float128(0.75), tiny);
    }
    {
        RoundingScope scope(FE_UPWARD);
        EXPECT_SAME(tiny / float128(4), tiny);
        EXPECT_SAME(-tiny / float128(4), -float128());
    }
    {
        RoundingScope scope(FE_DOWNWARD);
        EXPECT_SAME(tiny * float128(0.75), float128());
        EXPECT_SAME(-tiny / float128(4), -tiny);
    }
}

TEST(float128_env, Flags)
{
    const float128 one = float128::one();
    const float128 max = std::numeric_limits<float128>::max();
    const float128 min = std::numeric_limits<float128>::min();

    EXPECT_EQ(FlagsOf([&] { return one / float128(4); }), 0);
    EXPECT_EQ(FlagsOf([&] { return one / float128(3); }), FE_INEXACT);
    EXPECT_EQ(FlagsOf([&] { return max * float128(2); }), FE_OVERFLOW | FE_INEXACT);
    EXPECT_EQ(FlagsOf([&] { return min / float128(3); }), FE_UNDERFLOW | FE_INEXACT);
    // An exact subnormal result is not an underflow.
    EXPECT_EQ(FlagsOf([&] { return min / float128(4); }), 0);
    EXPECT_EQ(FlagsOf([&] { return one / float128(); }), FE_DIVBYZERO);
    EXPECT_EQ(FlagsOf([&] { return float128::inf() / float128(); }), 0);

    EXPECT_EQ(FlagsOf([&] { return float128::inf() - float128::inf(); }), FE_INVALID);
    EXPECT_EQ(FlagsOf([&] { return float128::inf() * float128(); }), FE_INVALID);
    EXPECT_EQ(FlagsOf([&] { return float128() / float128(); }), FE_INVALID);
    EXPECT_EQ(FlagsOf([&] { return sqrt(-one); }), FE_INVALID);
    EXPECT_EQ(FlagsOf([&] { return fmod(one, float128()); }), FE_INVALID);
    EXPECT_EQ(FlagsOf([&] { return float128::signaling_nan() + one; }), FE_INVALID);
    EXPECT_EQ(FlagsOf([&] { return float128::nan() + one; }), 0);

    // Relational comparisons signal on any NaN, equality and the <cmath> predicates only on a
    // signaling one.
    EXPECT_EQ(FlagsOf([&] { return float128::nan() < one; }), FE_INVALID);
    EXPECT_EQ(FlagsOf([&] { return float128::nan() >= one; }), FE_INVALID);
    EXPECT_EQ(FlagsOf([&] { return float128::nan() == one; }), 0);
    EXPECT_EQ(FlagsOf([&] { return float128::signaling_nan() == one; }), FE_INVALID);
    EXPECT_EQ(FlagsOf([&] { return isless(float128::nan(), one); }), 0);

    // Integer conversions: out of range and NaN are invalid, an inexact one raises nothing
    EXPECT_EQ(FlagsOf([&] { return static_cast<int32_t>(float128(1.5)); }), 0);
    EXPECT_EQ(FlagsOf([&] { return static_cast<int32_t>(float128(1e10)); }), FE_INVALID);
    EXPECT_EQ(FlagsOf([&] { return static_cast<int64_t>(float128::nan()); }), FE_INVALID);
    EXPECT_EQ(FlagsOf([&] { return static_cast<uint64_t>(float128(-1)); }), FE_INVALID);
    EXPECT_EQ(FlagsOf([&] { return static_cast<uint64_t>(float128(-0.5)); }), 0);

    // roundToIntegralExact raises inexact, the other integral roundings do not
    EXPECT_EQ(FlagsOf([&] { return rint(float128(1.5)); }), FE_INEXACT);
    EXPECT_EQ(FlagsOf([&] { return rint(float128(2)); }), 0);
    EXPECT_EQ(FlagsOf([&] { return nearbyint(float128(1.5)); }), 0);
    EXPECT_EQ(FlagsOf([&] { return floor(float128(1.5)); }), 0);

    // Conversions
    EXPECT_EQ(FlagsOf([&] { return static_cast<double>(one / float128(3)); }), FE_INEXACT);
    EXPECT_EQ(FlagsOf([&] { return static_cast<double>(float128(0.5)); }), 0);
    EXPECT_EQ(FlagsOf([&] { return static_cast<double>(max); }), FE_OVERFLOW | FE_INEXACT);
    EXPECT_EQ(FlagsOf([&] { return float128(std::numeric_limits<double>::signaling_NaN()); }), FE_INVALID);
    EXPECT_EQ(FlagsOf([&] { return float128("0.5"); }), 0);
    EXPECT_EQ(FlagsOf([&] { return float128("0.1"); }), FE_INEXACT);
    EXPECT_EQ(FlagsOf([&] { return float128("1e5000"); }), FE_OVERFLOW | FE_INEXACT);
}

TEST(float128_env, FlagsOfTheMathFunctions)
{
    // Only what the result justifies: the many roundings inside a function raise nothing of their
    // own beyond the inexact flag.
    const float128 one = float128::one();
    EXPECT_EQ(FlagsOf([&] { return exp(float128(20000)); }), FE_OVERFLOW | FE_INEXACT);
    EXPECT_EQ(FlagsOf([&] { return exp(float128(-20000)); }), FE_UNDERFLOW | FE_INEXACT);
    EXPECT_EQ(FlagsOf([&] { return exp(one); }), FE_INEXACT);
    EXPECT_EQ(FlagsOf([&] { return sin(Pow2(-13000)); }), FE_INEXACT);
    EXPECT_EQ(FlagsOf([&] { return log(one); }), 0);
    EXPECT_EQ(FlagsOf([&] { return log(float128()); }), FE_DIVBYZERO);
    EXPECT_EQ(FlagsOf([&] { return log(-one); }), FE_INVALID);
    EXPECT_EQ(FlagsOf([&] { return atanh(one); }), FE_DIVBYZERO);
    EXPECT_EQ(FlagsOf([&] { return pow(float128(2), float128(10)); }), 0);
    EXPECT_EQ(FlagsOf([&] { return pow(float128(), float128(-1)); }), FE_DIVBYZERO);
    EXPECT_EQ(FlagsOf([&] { return sqrt(float128(4)); }), 0);
    EXPECT_EQ(FlagsOf([&] { return sqrt(float128(2)); }), FE_INEXACT);
    EXPECT_EQ(FlagsOf([&] { return cosh(float128(11356.6)); }), FE_INEXACT);
    EXPECT_EQ(FlagsOf([&] { return hypot(float128(3), float128(4)); }), 0);
}

TEST(float128_env, FlagsStayRaisedUntilLowered)
{
    fp128::feclearexcept(FE_ALL_EXCEPT);
    [[maybe_unused]] const float128 a = float128(1) / float128(3);
    [[maybe_unused]] const float128 b = float128(1) / float128(4);
    EXPECT_EQ(fp128::fetestexcept(FE_ALL_EXCEPT), FE_INEXACT);

    std::fexcept_t saved {};
    fp128::fegetexceptflag(&saved, FE_ALL_EXCEPT);
    fp128::feclearexcept(FE_INEXACT);
    EXPECT_EQ(fp128::fetestexcept(FE_ALL_EXCEPT), 0);
    fp128::feraiseexcept(FE_OVERFLOW);
    fp128::fesetexceptflag(&saved, FE_INEXACT);
    EXPECT_EQ(fp128::fetestexcept(FE_ALL_EXCEPT), FE_INEXACT | FE_OVERFLOW);
    fp128::fesetexceptflag(&saved, FE_OVERFLOW);
    EXPECT_EQ(fp128::fetestexcept(FE_ALL_EXCEPT), FE_INEXACT);
    fp128::feclearexcept(FE_ALL_EXCEPT);

    // The float128 environment is separate from the hardware one.
    EXPECT_EQ(fp128::fetestexcept(FE_ALL_EXCEPT), 0);
}

TEST(float128_env, TheEnvironmentIsPerThread)
{
    fp128::feclearexcept(FE_ALL_EXCEPT);
    fp128::fesetround(FE_UPWARD);
    int other_round = 0;
    int other_flags = 0;
    std::thread other([&] {
        other_round = fp128::fegetround();
        [[maybe_unused]] const float128 x = float128(1) / float128(3);
        other_flags = fp128::fetestexcept(FE_ALL_EXCEPT);
    });
    other.join();
    EXPECT_EQ(other_round, FE_TONEAREST);
    EXPECT_EQ(other_flags, FE_INEXACT);
    EXPECT_EQ(fp128::fetestexcept(FE_ALL_EXCEPT), 0);
    fp128::fesetround(FE_TONEAREST);
}

TEST(float128_env, ConversionsFollowTheDirection)
{
    const float128 third = float128(1) / float128(3);
    {
        RoundingScope scope(FE_UPWARD);
        // to nearest, 1/3 rounds down as a double and up as a float
        EXPECT_EQ(std::bit_cast<uint64_t>(static_cast<double>(third)), 0x3FD5555555555556ull);
        EXPECT_EQ(std::bit_cast<uint32_t>(static_cast<float>(third)), 0x3EAAAAABu);
        EXPECT_EQ(std::format("{:.3e}", third), "3.334e-01");
        EXPECT_EQ(std::format("{:.3e}", -third), "-3.333e-01");
        EXPECT_EQ(std::format("{:.1a}", float128("0x1.01p0")), "0x1.1p+00");
        // to nearest, 0.1 already rounds up
        EXPECT_SAME(float128("0.1"), float128::tenth());
    }
    {
        RoundingScope scope(FE_DOWNWARD);
        EXPECT_EQ(std::bit_cast<uint64_t>(static_cast<double>(third)), 0x3FD5555555555555ull);
        EXPECT_EQ(std::bit_cast<uint32_t>(static_cast<float>(third)), 0x3EAAAAAAu);
        EXPECT_EQ(std::format("{:.3e}", third), "3.333e-01");
        EXPECT_EQ(std::format("{:.3e}", -third), "-3.334e-01");
        EXPECT_EQ(std::format("{:.2f}", float128(0.001)), "0.00");
        EXPECT_EQ(std::format("{:.2f}", float128(-0.001)), "-0.01");
        EXPECT_SAME(float128("0.1"), float128::nextDown(float128::tenth()));
    }
    {
        RoundingScope scope(FE_TOWARDZERO);
        EXPECT_EQ(std::format("{:.3e}", -third), "-3.333e-01");
        EXPECT_EQ(static_cast<int32_t>(float128(2.9)), 2);  // a C++ conversion truncates in any direction
        EXPECT_SAME(rint(float128(2.9)), float128(2));
        EXPECT_SAME(round(float128(2.5)), float128(3));  // round() ignores the direction
    }
}

TEST(float128_env, ConstantEvaluationRoundsToNearest)
{
    // A constant cannot depend on the thread that compiled it, so constant evaluation always
    // rounds to nearest; the same expression at run time follows the direction.
    constexpr float128 sum = float128::one() + ldexp(float128::one(), -200);
    static_assert(sum == float128::one());
    RoundingScope scope(FE_UPWARD);
    EXPECT_SAME(sum, float128::one());
    EXPECT_SAME(float128::one() + Pow2(-200), float128::nextUp(float128::one()));
}
