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
#include <vector>
#include "float128_ref_check.h"

/**********************************************************************
 * float128 IEEE 754-2008 conformance
 *
 * The operations clause 5 requires to be correctly rounded are checked bit for bit against
 * references computed on exact rationals, the sign of a zero included. Most of the rest of the
 * suite compares against a double, which only pins the top 53 of the 113 mantissa bits, and so
 * could not see that division truncated close to half of its results, that addition misjudged the
 * rounding whenever the bits that settled it lay below a three bit window, or that every result in
 * the subnormal range was rounded twice.
 *
 * Each of those defects also has a test of its own below, named for what it guards, alongside the
 * special values of clauses 6 and 9 and the operations clause 5 lists that the library lacked.
 ***********************************************************************/

using namespace f128_ref;

namespace
{
/// @brief A value from its encoding, written high QWORD first, the way the bits read.
[[nodiscard]] float128 Encoded(uint64_t high, uint64_t low)
{
    return float128(low, high);
}

/// @brief 2^e, exactly.
[[nodiscard]] float128 Pow2(int e)
{
    return ldexp(float128::one(), e);
}

/// @brief EXPECT that two values have the same encoding, or are both quiet NaNs.
#define EXPECT_SAME(actual, expected) \
    EXPECT_TRUE(SameResult((actual), (expected))) << #actual << ": expected " << Bits(expected) << ", got " << Bits(actual)

/// @brief True for a quiet NaN.
[[nodiscard]] bool IsQuietNan(const float128& x)
{
    return isnan(x) && !x.is_signaling();
}
}  // namespace

/**********************************************************************
 * The correctly rounded operations, against their exact reference tables
 ***********************************************************************/

TEST(float128_ieee, AdditionIsCorrectlyRounded)
{
    CheckExactBinary("add", add_exact_ref, [](const float128& x, const float128& y) { return x + y; });
    CheckExactBinary("sub", sub_exact_ref, [](const float128& x, const float128& y) { return x - y; });
}
TEST(float128_ieee, MultiplicationIsCorrectlyRounded)
{
    CheckExactBinary("mul", mul_exact_ref, [](const float128& x, const float128& y) { return x * y; });
}
TEST(float128_ieee, DivisionIsCorrectlyRounded)
{
    CheckExactBinary("div", div_exact_ref, [](const float128& x, const float128& y) { return x / y; });
}
TEST(float128_ieee, SquareRootIsCorrectlyRounded)
{
    CheckExactUnary("sqrt", sqrt_exact_ref, [](const float128& x) { return sqrt(x); });
    CheckExactUnary("sqrt", sqrt_ref, [](const float128& x) { return sqrt(x); });
}
TEST(float128_ieee, FusedMultiplyAddIsCorrectlyRounded)
{
    CheckExactTernary("fma", fma_exact_ref, [](const float128& x, const float128& y, const float128& z) { return fma(x, y, z); });
}
TEST(float128_ieee, NarrowingConversionsAreCorrectlyRounded)
{
    // A double and a float convert back exactly, so the round trip shows the narrowed value.
    CheckExactUnary("to_double", to_double_exact_ref, [](const float128& x) { return float128(static_cast<double>(x)); });
    CheckExactUnary("to_float", to_float_exact_ref, [](const float128& x) { return float128(static_cast<float>(x)); });
}
TEST(float128_ieee, RoundingToIntegralIsExact)
{
    CheckExactUnary("floor", floor_exact_ref, [](const float128& x) { return floor(x); });
    CheckExactUnary("ceil", ceil_exact_ref, [](const float128& x) { return ceil(x); });
    CheckExactUnary("trunc", trunc_exact_ref, [](const float128& x) { return trunc(x); });
    CheckExactUnary("round", round_exact_ref, [](const float128& x) { return round(x); });
    CheckExactUnary("rint", rint_exact_ref, [](const float128& x) { return rint(x); });
    CheckExactUnary("nearbyint", rint_exact_ref, [](const float128& x) { return nearbyint(x); });
}

/**********************************************************************
 * The defects the tables above found, one at a time
 ***********************************************************************/

TEST(float128_ieee, AdditionKeepsTheBitsBelowTheRoundingWindow)
{
    const float128 one = float128::one();
    const float128 next = float128::nextUp(one);
    // 1 + 2^-113 * (1 + 2^-12) is above the halfway point between 1 and its neighbour, by a bit
    // that sits twelve places under the guard bit. Deciding the rounding from a three bit window,
    // as addition used to, saw a tie and kept 1.
    EXPECT_SAME(one + Pow2(-113) * (one + Pow2(-12)), next);
    EXPECT_SAME(one + Pow2(-113) * next, next);
    // the tie itself goes to the even neighbour
    EXPECT_SAME(one + Pow2(-113), one);
    EXPECT_SAME(next + Pow2(-113), float128::nextUp(next));
    EXPECT_SAME(-one - Pow2(-113) * (one + Pow2(-12)), -next);
}

TEST(float128_ieee, DivisionRoundsWhenTheDividendMantissaIsTheSmaller)
{
    // The quotient of a smaller mantissa by a larger one leads at bit 127 rather than 128, and a
    // fixed normalization shift took that for nothing to round: such quotients were truncated.
    const float128 x = Encoded(0x3fe9fb7c7824cb69, 0xba1540834a61b1b2);
    const float128 y = Encoded(0x406effffffffffff, 0xffffffffffffffff);
    EXPECT_SAME(x / y, Encoded(0x3f79fb7c7824cb69, 0xba1540834a61b1b3));
    EXPECT_SAME(float128(2) / float128(3), Encoded(0x3ffe555555555555, 0x5555555555555555));
    EXPECT_SAME(float128(1) / float128(3), Encoded(0x3ffd555555555555, 0x5555555555555555));
    EXPECT_SAME(float128(5) / float128(7), Encoded(0x3ffe6db6db6db6db, 0x6db6db6db6db6db7));
}

TEST(float128_ieee, SubnormalResultsAreRoundedOnce)
{
    // Products, quotients and scalings into the subnormal range used to be rounded to 113 bits and
    // then again to the subnormal width, which rounds the tie the first rounding created.
    EXPECT_SAME(Encoded(0xa1f41d27cfff6699, 0x26bba92dcacd475f) * Encoded(0x9dbe7091cd370d57, 0xe2bd7958f67c9b9a),
                Encoded(0x0000000000000000, 0x000000066a2e8c89));
    EXPECT_SAME(ldexp(Encoded(0x005678a8c2245154, 0x2d35b5dabd4c8b09), -117), Encoded(0x00000000000178a8, 0xc22451542d35b5db));

    const float128 denorm_min = std::numeric_limits<float128>::denorm_min();
    // half of the smallest subnormal is a tie with zero, which is even; anything above it is not
    EXPECT_SAME(denorm_min / float128(2), float128());
    EXPECT_SAME(denorm_min * float128(0.75), denorm_min);
    EXPECT_SAME((denorm_min * float128(3)) / float128(2), denorm_min * float128(2));
    EXPECT_SAME(float128("3.3e-4966"), denorm_min);
    EXPECT_SAME(float128("4e-4966"), denorm_min);
    EXPECT_SAME(float128("3.2e-4966"), float128());
}

TEST(float128_ieee, SquareRoot)
{
    EXPECT_SAME(sqrt(float128(2)), float128::sqrt_2());
    EXPECT_SAME(sqrt(float128(4)), float128(2));
    EXPECT_SAME(sqrt(float128(3)), float128::sqrt_3());
    EXPECT_SAME(sqrt(std::numeric_limits<float128>::denorm_min()), Pow2(-8247));
    EXPECT_SAME(sqrt(Pow2(-16493)), float128::sqrt_2() * Pow2(-8247));
    EXPECT_SAME(sqrt(-float128()), -float128());
    EXPECT_TRUE(IsQuietNan(sqrt(float128(-1))));
    EXPECT_SAME(sqrt(float128::inf()), float128::inf());
    // The iteration count no longer matters: the last bit is settled exactly.
    EXPECT_SAME(sqrt(float128(2), 0), float128::sqrt_2());
}

/**********************************************************************
 * Conversions
 ***********************************************************************/

TEST(float128_ieee, ConversionFromDouble)
{
    static_assert(float128(-0.0).is_negative() && float128(-0.0).is_zero());
    EXPECT_TRUE(float128(-0.0).is_negative());

    // A quiet NaN stays quiet and keeps its payload through the round trip; a signaling one is
    // quieted, which is what an operation on it has to do.
    const double quiet = std::bit_cast<double>(0x7FF8000000000123ull);
    const float128 q = quiet;
    EXPECT_TRUE(IsQuietNan(q));
    EXPECT_EQ(std::bit_cast<uint64_t>(static_cast<double>(q)), 0x7FF8000000000123ull);
    EXPECT_EQ(std::bit_cast<uint64_t>(static_cast<double>(-q)), 0xFFF8000000000123ull);

    const float128 s = std::bit_cast<double>(0x7FF0000000000001ull);
    EXPECT_TRUE(IsQuietNan(s));
}

TEST(float128_ieee, ConversionToDouble)
{
    // Overflow and underflow keep the sign.
    EXPECT_EQ(static_cast<double>(-Pow2(2000)), -HUGE_VAL);
    EXPECT_TRUE(std::signbit(static_cast<double>(-Pow2(-2000))));
    // The smallest subnormal double converts as itself, and half of it ties to zero.
    EXPECT_EQ(std::bit_cast<uint64_t>(static_cast<double>(Pow2(-1074))), 1u);
    EXPECT_EQ(std::bit_cast<uint64_t>(static_cast<double>(Pow2(-1075))), 0u);
    EXPECT_EQ(std::bit_cast<uint64_t>(static_cast<double>(Pow2(-1075) * float128::nextUp(float128::one()))), 1u);
    // Bits below the window decide: 1 + 2^-53 + 2^-100 is above the tie.
    const float128 above_tie = float128::one() + Pow2(-53) + Pow2(-100);
    EXPECT_EQ(static_cast<double>(above_tie), 1.0 + 0x1p-52);
    EXPECT_EQ(static_cast<double>(-above_tie), -1.0 - 0x1p-52);
}

TEST(float128_ieee, ConversionToFloatRoundsOnce)
{
    // Through double this rounded twice: to 1 + 2^-24, a tie, and from there to 1.
    const float128 v = float128::one() + Pow2(-24) + Pow2(-82);
    EXPECT_EQ(static_cast<float>(v), 1.0f + 0x1p-23f);
    EXPECT_EQ(static_cast<float>(-v), -1.0f - 0x1p-23f);
    EXPECT_EQ(static_cast<float>(float128::one() + Pow2(-24)), 1.0f);
}

TEST(float128_ieee, LongDoubleConversions)
{
    // Every long double is a float128, so the round trip is exact - including the bits an x87
    // long double has beyond a double, which a conversion through double used to drop.
    using ld_limits = std::numeric_limits<long double>;
    for (const long double v : {1.0L / 3.0L, -ld_limits::max(), ld_limits::min(), ld_limits::denorm_min(), 0.1L, 1.0L + ld_limits::epsilon()}) {
        EXPECT_EQ(static_cast<long double>(float128(v)), v);
    }
    EXPECT_SAME(float128(1.0L + ld_limits::epsilon()), float128::one() + Pow2(1 - LDBL_MANT_DIG));
    EXPECT_TRUE(std::signbit(static_cast<long double>(float128(-0.0L))));
    EXPECT_TRUE(isnan(float128(ld_limits::quiet_NaN())));
    EXPECT_TRUE(std::isnan(static_cast<long double>(float128::nan())));
    EXPECT_EQ(static_cast<long double>(-float128::inf()), -HUGE_VALL);

    // The conversion to long double rounds once, to its own precision.
    if constexpr (LDBL_MANT_DIG < 100) {
        const float128 above_tie = float128::one() + Pow2(-LDBL_MANT_DIG) + Pow2(-110);
        EXPECT_EQ(static_cast<long double>(above_tie), 1.0L + ld_limits::epsilon());
        EXPECT_EQ(static_cast<long double>(float128::one() + Pow2(-LDBL_MANT_DIG)), 1.0L);
    }
}

TEST(float128_ieee, IntegerConversionsTruncate)
{
    static_assert(static_cast<int32_t>(float128(1.75)) == 1);
    static_assert(static_cast<int32_t>(float128(-1.75)) == -1);
    EXPECT_EQ(static_cast<int32_t>(float128(1.75)), 1);
    EXPECT_EQ(static_cast<int32_t>(float128(-1.75)), -1);
    EXPECT_EQ(static_cast<int64_t>(float128(0.75)), 0);
    EXPECT_EQ(static_cast<int64_t>(float128(-0.99)), 0);
    EXPECT_EQ(static_cast<int32_t>(float128(2147483647.5)), INT32_MAX);
    EXPECT_EQ(static_cast<uint64_t>(float128(2.9)), 2u);

    // Out of range saturates, a NaN converts to zero.
    EXPECT_EQ(static_cast<int32_t>(float128(1e10)), INT32_MAX);
    EXPECT_EQ(static_cast<int32_t>(float128(-1e10)), INT32_MIN);
    EXPECT_EQ(static_cast<int64_t>(float128::inf()), INT64_MAX);
    EXPECT_EQ(static_cast<int64_t>(-float128::inf()), INT64_MIN);
    EXPECT_EQ(static_cast<int64_t>(-Pow2(63)), INT64_MIN);
    EXPECT_EQ(static_cast<int64_t>(float128::nan()), 0);
    EXPECT_EQ(static_cast<uint64_t>(float128(-5)), 0u);
    EXPECT_EQ(static_cast<uint64_t>(float128(-0.5)), 0u);
    EXPECT_EQ(static_cast<uint64_t>(float128(1e30)), UINT64_MAX);
    EXPECT_EQ(static_cast<uint32_t>(float128(-1)), 0u);
}

TEST(float128_ieee, RoundingToIntegralKeepsTheSignOfZero)
{
    static_assert(ceil(float128(-0.5)).is_negative());
    EXPECT_TRUE(ceil(float128(-0.5)).is_negative() && ceil(float128(-0.5)).is_zero());
    EXPECT_TRUE(trunc(float128(-0.75)).is_negative() && trunc(float128(-0.75)).is_zero());
    EXPECT_TRUE(round(float128(-0.25)).is_negative() && round(float128(-0.25)).is_zero());
    EXPECT_TRUE(rint(float128(-0.25)).is_negative() && rint(float128(-0.25)).is_zero());
    EXPECT_TRUE(nearbyint(float128(-0.5)).is_negative() && nearbyint(float128(-0.5)).is_zero());
    EXPECT_TRUE(floor(float128(0.25)).is_positive() && floor(float128(0.25)).is_zero());

    float128 integral;
    const float128 fraction = modf(float128(-3), &integral);
    EXPECT_TRUE(fraction.is_zero() && fraction.is_negative());
    EXPECT_SAME(integral, float128(-3));
    const float128 inf_fraction = modf(-float128::inf(), &integral);
    EXPECT_TRUE(inf_fraction.is_zero() && inf_fraction.is_negative());
    EXPECT_SAME(integral, -float128::inf());
    EXPECT_TRUE(IsQuietNan(modf(float128::signaling_nan(), &integral)) && IsQuietNan(integral));
}

TEST(float128_ieee, RoundBreaksTiesAwayFromZero)
{
    // round() used to be trunc(x + 0.5), and that addition rounds: the largest value below one half
    // came out as one, and 2^112 + 1 as 2^112 + 2.
    EXPECT_SAME(round(float128::half() - Pow2(-114)), float128());
    EXPECT_SAME(round(Pow2(112) + float128(1)), Pow2(112) + float128(1));
    EXPECT_SAME(round(float128(2.5)), float128(3));
    EXPECT_SAME(round(float128(-2.5)), float128(-3));
    EXPECT_SAME(round(float128(0.5)), float128(1));
    EXPECT_SAME(rint(float128(2.5)), float128(2));
    EXPECT_SAME(rint(float128(3.5)), float128(4));
    EXPECT_EQ(llround(float128(2.5)), 3);
    EXPECT_EQ(llrint(float128(2.5)), 2);
    EXPECT_EQ(llround(float128::nan()), 0);
    EXPECT_EQ(llrint(float128(1e30)), 0);
}

/**********************************************************************
 * NaNs and the operations of clause 5.7
 ***********************************************************************/

TEST(float128_ieee, SignalingNaNsAreQuieted)
{
    const float128 s = float128::signaling_nan();
    const float128 one = float128::one();
    EXPECT_TRUE(IsQuietNan(s + one));
    EXPECT_TRUE(IsQuietNan(s * one));
    EXPECT_TRUE(IsQuietNan(s / one));
    EXPECT_TRUE(IsQuietNan(sqrt(s)));
    EXPECT_TRUE(IsQuietNan(fma(one, one, s)));
    EXPECT_TRUE(IsQuietNan(exp(s)));
    EXPECT_TRUE(IsQuietNan(sin(s)));
    EXPECT_TRUE(IsQuietNan(pow(s, float128(2))));
    EXPECT_TRUE(IsQuietNan(floor(s)));
    EXPECT_TRUE(IsQuietNan(rint(s)));
    EXPECT_TRUE(IsQuietNan(ldexp(s, 1)));
    EXPECT_TRUE(IsQuietNan(s << 1));
    EXPECT_TRUE(IsQuietNan(nextafter(s, one)));
    EXPECT_TRUE(IsQuietNan(logb(s)));
    // minNum treats a quiet NaN as missing, but a signaling one is invalid
    EXPECT_TRUE(IsQuietNan(fmin(s, one)));
    EXPECT_TRUE(IsQuietNan(fmax(one, s)));
    EXPECT_SAME(fmin(float128::nan(), one), one);
    // an infinity beats a quiet NaN in hypot, not a signaling one
    EXPECT_SAME(hypot(float128::inf(), float128::nan()), float128::inf());
    EXPECT_TRUE(IsQuietNan(hypot(float128::inf(), s)));
    // copying, negating and taking the magnitude are not operations on the value
    EXPECT_TRUE((-s).is_signaling() && fabs(s).is_signaling() && copysign(s, -one).is_signaling());
}

TEST(float128_ieee, NaNPayloadsPropagate)
{
    const float128 payload = fp128::nan("123");
    uint64_t want_low = 0, want_high = 0;
    payload.get_bits(want_low, want_high);

    for (const float128& result : {payload + float128(1), float128(1) * payload, sqrt(payload), exp(payload), -(-payload - float128(2))}) {
        uint64_t low = 0, high = 0;
        result.get_bits(low, high);
        EXPECT_EQ(low, want_low);
        EXPECT_EQ(high & ~(1ull << 63), want_high);
    }
    // the first NaN operand is the one delivered
    EXPECT_EQ(Bits(fp128::nan("7") + fp128::nan("9")), Bits(fp128::nan("7")));
}

TEST(float128_ieee, ClassDistinguishesSignalingNaNs)
{
    EXPECT_NE(float128::signaling_nan().get_class(), float128::nan().get_class());
    EXPECT_TRUE(issignaling(float128::signaling_nan()));
    EXPECT_FALSE(issignaling(float128::nan()));
}

TEST(float128_ieee, TotalOrder)
{
    const float128 one = float128::one();
    const float128 max = std::numeric_limits<float128>::max();
    const float128 tiny = std::numeric_limits<float128>::denorm_min();
    const float128 quiet = fp128::nan("5");
    const float128 signaling = float128::signaling_nan();
    const std::vector<float128> ascending = {-quiet, -signaling, -float128::inf(), -max, -one, -tiny, -float128(), float128(), tiny, one, max,
                                             float128::inf(), signaling, quiet};
    for (size_t i = 0; i < ascending.size(); ++i) {
        for (size_t j = 0; j < ascending.size(); ++j)
            EXPECT_EQ(totalorder(ascending[i], ascending[j]), i <= j) << i << ", " << j;
    }
    EXPECT_TRUE(totalordermag(-one, float128(2)));
    EXPECT_FALSE(totalordermag(float128(-2), one));
    EXPECT_TRUE(totalordermag(-float128(), float128()) && totalordermag(float128(), -float128()));
}

TEST(float128_ieee, MagnitudeMinimumAndMaximum)
{
    EXPECT_SAME(fmaxmag(float128(-3), float128(2)), float128(-3));
    EXPECT_SAME(fminmag(float128(-3), float128(2)), float128(2));
    // equal magnitudes fall back on fmax and fmin
    EXPECT_SAME(fmaxmag(float128(-2), float128(2)), float128(2));
    EXPECT_SAME(fminmag(float128(-2), float128(2)), float128(-2));
    EXPECT_SAME(fmaxmag(float128::nan(), float128(-2)), float128(-2));
    EXPECT_TRUE(IsQuietNan(fminmag(float128::signaling_nan(), float128(2))));
}

TEST(float128_ieee, Predicates)
{
    EXPECT_TRUE(iscanonical(float128::nan()) && iscanonical(float128(3)));
    EXPECT_TRUE(issubnormal(std::numeric_limits<float128>::denorm_min()));
    EXPECT_FALSE(issubnormal(float128()));
    EXPECT_FALSE(issubnormal(std::numeric_limits<float128>::min()));
    EXPECT_TRUE(iszero(-float128()));
    EXPECT_FALSE(iszero(std::numeric_limits<float128>::denorm_min()));
}

/**********************************************************************
 * Special operands and the edges of the range
 *
 * The exact tables above hold finite, non zero operands. The cases below are the ones they cannot:
 * infinities, zeros and the sign rules that decide them, NaNs, and the ends of the range. Every
 * expected encoding comes from exact rationals (reviewer D's encodings_for_proposals.py and
 * fma_tie_case.txt), not from the library.
 ***********************************************************************/

TEST(float128_ieee, FmaSpecialOperands)
{
    // IEEE 754-2008 7.2 and 6.3: the product's own special cases come first, an infinite addend
    // wins over a finite product, and a zero sum follows the sign rule of an exact zero.
    const float128 inf = float128::inf();
    const float128 zero;
    const float128 one = float128::one();
    const float128 denorm = std::numeric_limits<float128>::denorm_min();
    const float128 max = std::numeric_limits<float128>::max();
    EXPECT_TRUE(IsQuietNan(fma(inf, zero, one)));
    EXPECT_TRUE(IsQuietNan(fma(zero, -inf, float128::nan())));
    EXPECT_TRUE(IsQuietNan(fma(zero, inf, float128::signaling_nan())));
    EXPECT_TRUE(IsQuietNan(fma(inf, float128(2), -inf)));
    EXPECT_SAME(fma(inf, float128(-2), -inf), -inf);
    EXPECT_SAME(fma(-inf, -inf, inf), inf);
    EXPECT_SAME(fma(float128(2), float128(3), inf), inf);
    EXPECT_SAME(fma(float128(2), float128(3), -inf), -inf);
    EXPECT_SAME(fma(zero, float128(5), float128(7)), float128(7));
    EXPECT_SAME(fma(zero, float128(5), -zero), zero);
    EXPECT_SAME(fma(-zero, float128(5), -zero), -zero);
    EXPECT_SAME(fma(zero, float128(-5), zero), zero);
    EXPECT_SAME(fma(float128(2), float128(3), zero), float128(6));
    EXPECT_SAME(fma(float128(2), float128(3), -zero), float128(6));
    // a zero addend leaves the product, which underflows to a zero of its own sign
    EXPECT_SAME(fma(denorm, denorm, -zero), zero);
    EXPECT_SAME(fma(-denorm, denorm, zero), -zero);
    // the product is past the range, the sum is not
    EXPECT_SAME(fma(max, float128(2), -max), max);

    // Exact cancellation: (1 + 2^-60)^2 - (1 + 2^-59) is 2^-120, and 3 * nearest(1/3) - 1 is -2^-114.
    EXPECT_SAME(fma(one + Pow2(-60), one + Pow2(-60), -(one + Pow2(-59))), Pow2(-120));
    EXPECT_SAME(fma(Encoded(0x3FFD555555555555, 0x5555555555555555), float128(3), -one), -Pow2(-114));

    // y is (2^113 + 7) / 3 * 2^-112, so 3y = 2 + 7 * 2^-112 is an exact tie above the odd
    // significand of 2 + 3 * 2^-111: the product rounds it to even, up. An addend of the smallest
    // subnormal, far below the product, falls off the accumulator and only its sign survives - as
    // a borrow, which puts the sum below the tie, so fma() rounds down.
    const float128 y = Encoded(0x3FFE555555555555, 0x555555555555555A);
    EXPECT_SAME(float128(3) * y, Encoded(0x4000000000000000, 4));
    EXPECT_SAME(fma(float128(3), y, -denorm), Encoded(0x4000000000000000, 3));
    EXPECT_SAME(fma(float128(3), y, denorm), Encoded(0x4000000000000000, 4));
}

TEST(float128_ieee, ExactRemaindersKeepTheSignOfZero)
{
    const float128 zero;
    const float128 denorm = std::numeric_limits<float128>::denorm_min();
    const float128 max = std::numeric_limits<float128>::max();

    // A zero remainder takes the sign of x (IEEE 754-2008 5.3.1, and C's fmod).
    EXPECT_SAME(fmod(float128(6), float128(3)), zero);
    EXPECT_SAME(fmod(float128(-6), float128(3)), -zero);
    EXPECT_SAME(fmod(-zero, float128(3)), -zero);
    EXPECT_SAME(fmod(zero, float128(-3)), zero);
    EXPECT_SAME(fmod(max, denorm), zero);
    EXPECT_SAME(fmod(max, float128(3)), float128(2));
    EXPECT_SAME(fmod(float128(-5), float128::inf()), float128(-5));
    EXPECT_SAME(remainder(float128(6), float128(3)), zero);
    EXPECT_SAME(remainder(float128(-6), float128(3)), -zero);
    EXPECT_SAME(remainder(-zero, float128(3)), -zero);

    // 3d / 2d is a tie, and the even quotient is 2; 3 * 2^16000 / 7 is a quotient of 4800 digits.
    EXPECT_SAME(remainder(denorm * float128(3), denorm * float128(2)), -denorm);
    EXPECT_SAME(remainder(Pow2(16000) * float128(3), float128(7)), float128(-1));

    // The quotient's low bits carry the sign of x / y.
    int quotient = 0;
    EXPECT_SAME(remquo(float128(-7), float128(2), &quotient), float128(1));
    EXPECT_EQ(quotient, -4);
    EXPECT_SAME(remquo(float128(1000), float128(3), &quotient), float128(1));
    EXPECT_EQ(quotient, 333 & 0x7F);
    EXPECT_SAME(remquo(float128(-1000), float128(-3), &quotient), float128(-1));
    EXPECT_EQ(quotient, 333 & 0x7F);

    // At the top of the range. x / max is 0.75, so the remainder is x - max, exactly; doubling the
    // partial remainder in floating point, as remquo() used to, overflowed here.
    const float128 x = Encoded(0x7FFE800000000000, 0);  // 1.5 * 2^16383
    EXPECT_SAME(remainder(x, max), Encoded(0xFFFCFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFC));
    EXPECT_SAME(remainder(-x, max), Encoded(0x7FFCFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFC));
    EXPECT_SAME(remainder(x, -max), Encoded(0xFFFCFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFC));
    EXPECT_SAME(remquo(x, max, &quotient), Encoded(0xFFFCFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFC));
    EXPECT_EQ(quotient, 1);
    // half the largest value is a tie whose even quotient is 0, the next value up is past it
    const float128 half_max = Encoded(0x7FFDFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF);
    EXPECT_SAME(remainder(half_max, max), half_max);
    EXPECT_SAME(remainder(float128::nextUp(half_max), max), Encoded(0xFFFDFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFE));

    EXPECT_TRUE(IsQuietNan(remainder(float128::nan(), float128(1))));
    EXPECT_TRUE(IsQuietNan(remquo(float128(1), float128::signaling_nan(), &quotient)));
    EXPECT_EQ(quotient, 0);
}

TEST(float128_ieee, NeighboursInEveryDirection)
{
    const float128 inf = float128::inf();
    const float128 zero;
    const float128 one = float128::one();
    const float128 denorm = std::numeric_limits<float128>::denorm_min();
    const float128 min = std::numeric_limits<float128>::min();
    const float128 max = std::numeric_limits<float128>::max();
    const float128 largest_subnormal = Encoded(0x0000FFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF);

    // IEEE 754-2008 5.3.1 nextUp and nextDown, at every boundary a step can cross
    EXPECT_SAME(float128::nextUp(-inf), -max);
    EXPECT_SAME(float128::nextUp(inf), inf);
    EXPECT_SAME(float128::nextDown(inf), max);
    EXPECT_SAME(float128::nextDown(-inf), -inf);
    EXPECT_SAME(float128::nextUp(max), inf);
    EXPECT_SAME(float128::nextDown(-max), -inf);
    EXPECT_SAME(float128::nextUp(-denorm), -zero);
    EXPECT_SAME(float128::nextDown(denorm), zero);
    EXPECT_SAME(float128::nextUp(-zero), denorm);
    EXPECT_SAME(float128::nextDown(zero), -denorm);
    EXPECT_SAME(float128::nextUp(-min), -largest_subnormal);
    EXPECT_SAME(float128::nextDown(-largest_subnormal), -min);
    EXPECT_SAME(float128::nextUp(-one), Encoded(0xBFFEFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF));
    EXPECT_SAME(float128::nextDown(-one), Encoded(0xBFFF000000000000, 1));
    EXPECT_TRUE(IsQuietNan(float128::nextUp(float128::signaling_nan())));
    EXPECT_TRUE(IsQuietNan(float128::nextDown(float128::nan())));

    // nextafter returns y itself when the two are equal, so the sign of a zero comes from y
    EXPECT_SAME(nextafter(zero, -zero), -zero);
    EXPECT_SAME(nextafter(-zero, zero), zero);
    EXPECT_SAME(nextafter(-one, -one), -one);
    EXPECT_SAME(nextafter(inf, zero), max);
    EXPECT_SAME(nextafter(-inf, zero), -max);
    EXPECT_SAME(nextafter(-denorm, float128(-1)), Encoded(0x8000000000000000, 2));
    EXPECT_SAME(nexttoward(-one, -inf), Encoded(0xBFFF000000000000, 1));
}

TEST(float128_ieee, NarrowingSpecialsAndEdges)
{
    const float128 one = float128::one();
    // A float NaN keeps its sign and the top of its payload, and a signaling one is quieted.
    EXPECT_EQ(std::bit_cast<uint32_t>(static_cast<float>(-float128(std::bit_cast<float>(0x7FC00123u)))), 0xFFC00123u);
    EXPECT_EQ(std::bit_cast<uint32_t>(static_cast<float>(float128::signaling_nan())), 0x7FC00000u);
    EXPECT_EQ(std::bit_cast<uint32_t>(static_cast<float>(-float128::inf())), 0xFF800000u);
    EXPECT_EQ(std::bit_cast<uint32_t>(static_cast<float>(-float128())), 0x80000000u);
    EXPECT_EQ(std::bit_cast<uint32_t>(static_cast<float>(std::numeric_limits<float128>::max())), 0x7F800000u);
    // the overflow threshold of float is a tie, which goes to the even neighbour: up, to infinity
    EXPECT_EQ(std::bit_cast<uint32_t>(static_cast<float>(float128(FLT_MAX) + Pow2(103))), 0x7F800000u);
    EXPECT_EQ(std::bit_cast<uint32_t>(static_cast<float>(float128(FLT_MAX) + Pow2(102))), 0x7F7FFFFFu);
    // half the smallest subnormal ties with zero; one and a half of it rounds to two
    EXPECT_EQ(std::bit_cast<uint32_t>(static_cast<float>(Pow2(-150))), 0u);
    EXPECT_EQ(std::bit_cast<uint32_t>(static_cast<float>(-Pow2(-149) * float128(1.5))), 0x80000002u);

    // the same edges for double
    EXPECT_EQ(std::bit_cast<uint64_t>(static_cast<double>(float128(DBL_MAX) + Pow2(970))), 0x7FF0000000000000ull);
    EXPECT_EQ(std::bit_cast<uint64_t>(static_cast<double>(float128(DBL_MAX) + Pow2(963))), 0x7FEFFFFFFFFFFFFFull);
    EXPECT_EQ(std::bit_cast<uint64_t>(static_cast<double>(float128(DBL_MIN) * (one - Pow2(-60)))), 0x0010000000000000ull);
    EXPECT_EQ(std::bit_cast<uint64_t>(static_cast<double>(-Pow2(-1080))), 0x8000000000000000ull);
    EXPECT_EQ(std::bit_cast<uint64_t>(static_cast<double>(std::numeric_limits<float128>::denorm_min())), 0u);
    EXPECT_EQ(std::bit_cast<uint64_t>(static_cast<double>(float128::signaling_nan())), 0x7FF8000000000000ull);
}

TEST(float128_ieee, FrexpLogbIlogbOfSpecials)
{
    const float128 inf = float128::inf();
    const float128 zero;
    const float128 denorm = std::numeric_limits<float128>::denorm_min();
    const float128 max = std::numeric_limits<float128>::max();

    int e = 99;
    EXPECT_SAME(frexp(-zero, &e), -zero);
    EXPECT_EQ(e, 0);
    e = 99;
    EXPECT_SAME(frexp(-inf, &e), -inf);
    EXPECT_EQ(e, 0);
    e = 99;
    // a signaling NaN is quieted, as every other operation does to it
    EXPECT_TRUE(IsQuietNan(frexp(float128::signaling_nan(), &e)));
    EXPECT_EQ(e, 0);
    EXPECT_SAME(frexp(-denorm, &e), -float128::half());
    EXPECT_EQ(e, -16493);
    EXPECT_SAME(frexp(max, &e), Encoded(0x3FFEFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF));
    EXPECT_EQ(e, 16384);

    EXPECT_SAME(logb(-inf), inf);
    EXPECT_SAME(logb(-zero), -inf);
    EXPECT_TRUE(IsQuietNan(logb(float128::signaling_nan())));
    EXPECT_SAME(logb(denorm), float128(-16494));
    EXPECT_SAME(logb(-max), float128(16383));
    EXPECT_EQ(ilogb(max), 16383);
    EXPECT_EQ(ilogb(-denorm), -16494);
    EXPECT_EQ(ilogb(-zero), FP_ILOGB0);
    EXPECT_EQ(ilogb(-inf), INT_MAX);
}

TEST(float128_ieee, ModfFdimMagnitudeSigns)
{
    const float128 inf = float128::inf();
    const float128 zero;
    const float128 one = float128::one();
    const float128 denorm = std::numeric_limits<float128>::denorm_min();
    const float128 max = std::numeric_limits<float128>::max();

    float128 integral;
    EXPECT_SAME(modf(float128(3), &integral), zero);
    EXPECT_SAME(integral, float128(3));
    EXPECT_SAME(modf(float128(-2.5), &integral), -float128::half());
    EXPECT_SAME(integral, float128(-2));
    EXPECT_SAME(modf(-denorm, &integral), -denorm);
    EXPECT_SAME(integral, -zero);

    // fdim is +0 unless x > y, whatever the signs of the zeros
    EXPECT_SAME(fdim(one, one), zero);
    EXPECT_SAME(fdim(-zero, zero), zero);
    EXPECT_SAME(fdim(-inf, -inf), zero);
    EXPECT_SAME(fdim(inf, -inf), inf);
    EXPECT_SAME(fdim(max, -max), inf);
    EXPECT_TRUE(IsQuietNan(fdim(one, float128::signaling_nan())));

    // the zeros compare equal, so which one comes back is decided on the sign, in either order
    EXPECT_SAME(fmax(-zero, zero), zero);
    EXPECT_SAME(fmax(zero, -zero), zero);
    EXPECT_SAME(fmin(zero, -zero), -zero);
    EXPECT_SAME(fmin(-zero, zero), -zero);
    EXPECT_SAME(fmaxmag(-zero, zero), zero);
    EXPECT_SAME(fminmag(-zero, zero), -zero);
    EXPECT_SAME(fmaxmag(one, -one), one);
    EXPECT_SAME(fminmag(one, -one), -one);
    EXPECT_SAME(fminmag(-inf, one), one);
    EXPECT_SAME(fmaxmag(-inf, one), -inf);
    EXPECT_SAME(fmaxmag(float128(2), float128(-3)), float128(-3));
    EXPECT_SAME(fminmag(float128(-3), float128(2)), float128(2));
    EXPECT_TRUE(IsQuietNan(fmaxmag(float128::nan(), fp128::nan("7"))));
}

/**********************************************************************
 * Special values of the recommended functions (IEEE 754-2008 9.2.1)
 ***********************************************************************/

TEST(float128_ieee, LogarithmSpecialValues)
{
    const float128 inf = float128::inf();
    for (auto fn : {+[](const float128& x) { return log(x); }, +[](const float128& x) { return log2(x); }, +[](const float128& x) { return log10(x); }}) {
        EXPECT_TRUE(IsQuietNan(fn(float128(-1))));
        EXPECT_TRUE(IsQuietNan(fn(-inf)));
        EXPECT_TRUE(IsQuietNan(fn(float128::nan())));
        EXPECT_SAME(fn(inf), inf);
        EXPECT_SAME(fn(float128()), -inf);
        EXPECT_SAME(fn(-float128()), -inf);
        EXPECT_SAME(fn(float128::one()), float128());
    }
    // a subnormal argument has an exponent of its own, not the one its field holds
    EXPECT_SAME(log2(std::numeric_limits<float128>::denorm_min()), float128(-16494));
    EXPECT_SAME(log1p(float128(-1)), -inf);
    EXPECT_TRUE(IsQuietNan(log1p(float128(-2))));
    EXPECT_SAME(log1p(-float128()), -float128());
    EXPECT_SAME(logb(float128()), -inf);
}

TEST(float128_ieee, PowSpecialValues)
{
    const float128 inf = float128::inf();
    const float128 nan = float128::nan();
    const float128 one = float128::one();
    const float128 half = float128::half();
    const float128 zero;
    struct Case {
        float128 x, y, want;
    };
    const Case cases[] = {
        {nan, zero, one},         {nan, -zero, one},       {one, nan, one},          {one, inf, one},          {-one, inf, one},
        {-one, -inf, one},        {half, inf, zero},       {float128(2), inf, inf},  {half, -inf, inf},        {float128(2), -inf, zero},
        {-zero, half, zero},      {-zero, -half, inf},     {-zero, float128(3), -zero}, {-zero, float128(-3), -inf}, {zero, float128(-3), inf},
        {-inf, half, inf},        {-inf, float128(3), -inf}, {-inf, float128(-3), -zero}, {-inf, -half, zero},  {inf, -half, zero},
        {float128(2), float128(-20000), zero}, {half, float128(20000), zero}, {half, float128(-20000), inf},
        {-one, Pow2(100), one},   {float128(-2), Pow2(100), inf}, {half, Pow2(100), zero}, {float128(2), float128(10), float128(1024)},
        {float128(-2), float128(3), float128(-8)},
    };
    for (const Case& c : cases)
        EXPECT_SAME(pow(c.x, c.y), c.want) << Bits(c.x) << " ^ " << Bits(c.y);
    EXPECT_TRUE(IsQuietNan(pow(float128(2), nan)));
    EXPECT_TRUE(IsQuietNan(pow(nan, one)));
    EXPECT_TRUE(IsQuietNan(pow(float128(-8), one / float128(3))));

    // Integer powers that fit are exact.
    EXPECT_SAME(pow(float128(10), float128(30)), float128("1e30"));
    EXPECT_SAME(pow(float128(3), 70), float128("2503155504993241601315571986085849"));
    EXPECT_SAME(pow(float128(2), -3), float128(0.125));

    // A large exponent on a base near one is an ordinary finite value; it used to overflow. The
    // references are correctly rounded (mpmath), and the bound is the 4 + 2.5 |y log2 x| ulp pow()
    // documents, which for these is 4: |y log2(1 + 2^-40)| is at most 0.004.
    const float128 base = one + Pow2(-40);
    const struct {
        float128 y, want;
    } large[] = {
        {float128(12000), Encoded(0x3FFF0000002EE000, 0x044A8A904309D8B1)},
        {float128(-12000), Encoded(0x3FFEFFFFFFA24000, 0x089572DF79DB2423)},
        {float128(3000000000), Encoded(0x3FFF00B30EDFC302, 0x839BC35E3C7E0DFD)},
    };
    for (const auto& c : large) {
        EXPECT_LE(Saturate64(UlpDistance(pow(base, c.y), c.want)), 4u) << Bits(c.y);
    }
}

TEST(float128_ieee, Atan2SpecialValues)
{
    const float128 pi = float128::pi();
    const float128 half_pi = float128::half_pi();
    const float128 quarter_pi = float128::quarter_pi();
    const float128 inf = float128::inf();
    const float128 zero;
    const float128 one = float128::one();
    struct Case {
        float128 y, x, want;
    };
    const Case cases[] = {
        {zero, zero, zero},     {-zero, zero, -zero},   {zero, -zero, pi},       {-zero, -zero, -pi},     {zero, one, zero},
        {-zero, one, -zero},    {zero, -one, pi},       {-zero, -one, -pi},      {one, zero, half_pi},    {-one, -zero, -half_pi},
        {inf, inf, quarter_pi}, {-inf, inf, -quarter_pi}, {one, inf, zero},      {-one, inf, -zero},      {one, -inf, pi},
        {-one, -inf, -pi},      {inf, one, half_pi},    {-inf, -one, -half_pi},
    };
    for (const Case& c : cases)
        EXPECT_SAME(atan2(c.y, c.x), c.want) << Bits(c.y) << ", " << Bits(c.x);
    // 3pi/4, correctly rounded
    EXPECT_SAME(atan2(inf, -inf), Encoded(0x40002d97c7f3321d, 0x234f272993d1414a));
    EXPECT_SAME(atan2(-inf, -inf), -Encoded(0x40002d97c7f3321d, 0x234f272993d1414a));
}

TEST(float128_ieee, InverseTrigonometricDomain)
{
    const float128 one = float128::one();
    EXPECT_TRUE(IsQuietNan(asin(float128(2))));
    EXPECT_TRUE(IsQuietNan(asin(float128(-2))));
    EXPECT_TRUE(IsQuietNan(acos(float128(2))));
    EXPECT_TRUE(IsQuietNan(acos(float128::inf())));
    EXPECT_SAME(asin(one), float128::half_pi());
    EXPECT_SAME(asin(-one), -float128::half_pi());
    EXPECT_SAME(acos(one), float128());
    EXPECT_SAME(acos(-one), float128::pi());
    EXPECT_SAME(asin(-float128()), -float128());

    // acos(1 - d) is about sqrt(2d); the Newton iteration this used divided by sin(acos(x)) and
    // returned a NaN here. The reference is correctly rounded (mpmath), and the bound the 3 ulp of
    // the acos accuracy test: comparing against sqrt(2d) allowed 2^17 ulp.
    EXPECT_LE(Saturate64(UlpDistance(acos(one - Pow2(-100)), Encoded(0x3FCD6A09E667F3BC, 0xC908B2FB1366EC78))), 3u);
}

TEST(float128_ieee, LargeTrigonometricArguments)
{
    // Above 2^321 the series sin() evaluated never terminated; between 2^62 and that it returned
    // values such as 2^3240. The reduction now reads 2/pi as deep as the argument needs.
    const float128 huge = Pow2(16000) * float128(1.375);
    const float128 s = sin(huge);
    const float128 c = cos(huge);
    EXPECT_TRUE(isfinite(s) && fabs(s) <= float128::one());
    EXPECT_LT(fabs(sqr(s) + sqr(c) - float128::one()), Pow2(-108));
    EXPECT_SAME(sin(-huge), -s);
    EXPECT_SAME(cos(-huge), c);
    EXPECT_TRUE(isfinite(tan(std::numeric_limits<float128>::max())));
    EXPECT_TRUE(isfinite(sin(Pow2(70) * float128(1.1))) && fabs(sin(Pow2(70) * float128(1.1))) <= float128::one());
}

TEST(float128_ieee, HyperbolicFunctionsNearOverflow)
{
    // exp(x) overflows a little before cosh(x) = exp(x)/2 does.
    EXPECT_TRUE(isfinite(cosh(float128(11356.6))));
    EXPECT_TRUE(isfinite(sinh(float128(-11356.6))) && sinh(float128(-11356.6)).is_negative());
    EXPECT_SAME(cosh(float128(11357.3)), float128::inf());
    EXPECT_SAME(cosh(-float128::inf()), float128::inf());
}

TEST(float128_ieee, HypotScales)
{
    EXPECT_SAME(hypot(float128::inf(), float128::nan()), float128::inf());
    EXPECT_SAME(hypot(float128::nan(), -float128::inf()), float128::inf());
    EXPECT_SAME(hypot(Pow2(9000), Pow2(9000)), float128::sqrt_2() * Pow2(9000));
    EXPECT_SAME(hypot(Pow2(-9000), Pow2(-9000)), float128::sqrt_2() * Pow2(-9000));
    EXPECT_SAME(hypot(float128(3), float128(-4)), float128(5));
}

/**********************************************************************
 * Text conversions (IEEE 754-2008 5.12)
 ***********************************************************************/

TEST(float128_ieee, HexadecimalInput)
{
    EXPECT_SAME(float128("0x1.8p1"), float128(3));
    EXPECT_SAME(float128("-0x1.8p1"), float128(-3));
    EXPECT_SAME(float128("0X1P-16494"), std::numeric_limits<float128>::denorm_min());
    EXPECT_SAME(float128("0x.8"), float128::half());
    EXPECT_SAME(float128("0x10"), float128(16));
    // Correctly rounded: 125 one bits round up to 2^125, the tie above 2 - 2^-112 to 2, the tie
    // above 1 to 1, and the one above 1 + 2^-112 to 1 + 2^-111.
    EXPECT_SAME(float128("0x1fffffffffffffffffffffffffffffff"), Pow2(125));
    EXPECT_SAME(float128("0x1.ffffffffffffffffffffffffffff8p0"), float128(2));
    EXPECT_SAME(float128("0x1.00000000000000000000000000008p0"), float128::one());
    EXPECT_SAME(float128("0x1.00000000000000000000000000018p0"), float128::one() + Pow2(-111));

    float128 value;
    const char text[] = "1.8p1";
    auto result = fp128::from_chars(text, text + 5, value, std::chars_format::hex);
    EXPECT_TRUE(result.ec == std::errc {} && result.ptr == text + 5);
    EXPECT_SAME(value, float128(3));
    const char partial[] = "1p";
    result = fp128::from_chars(partial, partial + 2, value, std::chars_format::hex);
    EXPECT_TRUE(result.ec == std::errc {} && result.ptr == partial + 1);
    EXPECT_SAME(value, float128::one());
    const char huge[] = "1p100000";
    result = fp128::from_chars(huge, huge + 8, value, std::chars_format::hex);
    EXPECT_TRUE(result.ec == std::errc::result_out_of_range);
    EXPECT_SAME(value, float128::inf());
}

TEST(float128_ieee, SpecialValueInput)
{
    EXPECT_TRUE(float128("snan").is_signaling());
    EXPECT_TRUE(float128("-SNaN").is_signaling() && float128("-SNaN").is_negative());
    EXPECT_TRUE(IsQuietNan(float128("nan")));
    EXPECT_EQ(Bits(float128("nan(0x7b)")), Bits(fp128::nan("0x7b")));
    EXPECT_EQ(Bits(float128("nan(123)")), Bits(fp128::nan("123")));
    EXPECT_SAME(float128("-infinity"), -float128::inf());

    // An unclosed payload is not part of the number.
    float128 value;
    const char text[] = "nan(12";
    const auto result = fp128::from_chars(text, text + 6, value);
    EXPECT_TRUE(result.ptr == text + 3 && IsQuietNan(value));
}

TEST(float128_ieee, HexadecimalOutputTiesToEven)
{
    // The first dropped digit decides; a tie goes to the even digit, it used to go up.
    EXPECT_EQ(std::format("{:.0a}", float128(1.5)), "0x2p+00");
    EXPECT_EQ(std::format("{:.0a}", float128(2.5)), "0x1p+01");
    EXPECT_EQ(std::format("{:.1a}", float128("0x1.18p0")), "0x1.2p+00");
    EXPECT_EQ(std::format("{:.1a}", float128("0x1.28p0")), "0x1.2p+00");
    EXPECT_EQ(std::format("{:.1a}", float128("0x1.2801p0")), "0x1.3p+00");
}

TEST(float128_ieee, DecimalOutputOfSubnormalsIsShortest)
{
    // The shortest form is the one that reads back to the same value at the value's own precision,
    // which for a subnormal is fewer than 113 bits.
    EXPECT_EQ(to_string(std::numeric_limits<float128>::denorm_min()), "6e-4966");
    EXPECT_SAME(float128(to_string(std::numeric_limits<float128>::denorm_min() * float128(3))), std::numeric_limits<float128>::denorm_min() * float128(3));
}

/**********************************************************************
 * The floating point environment, as seen without FP128_IEEE_ENV
 ***********************************************************************/

#ifndef FP128_IEEE_ENV
TEST(float128_ieee, EnvironmentIsFixedWithoutTheMacro)
{
    EXPECT_EQ(fp128::fegetround(), FE_TONEAREST);
    EXPECT_NE(fp128::fesetround(FE_UPWARD), 0);
    EXPECT_EQ(fp128::fegetround(), FE_TONEAREST);
    EXPECT_EQ(fp128::fesetround(FE_TONEAREST), 0);
    EXPECT_NE(fp128::feraiseexcept(FE_INEXACT), 0);
    EXPECT_EQ(fp128::fetestexcept(FE_ALL_EXCEPT), 0);
    static_assert(!std::numeric_limits<float128>::is_iec559);
    static_assert(!float128::is754version2008());
}
#endif
