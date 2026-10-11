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
#include <array>
#include <bit>
#include <ostream>
#include <ctime>
#include <cfloat>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <string>
#include <type_traits>
#include "gtest_shared.h"

namespace
{
/**
 * @brief Two's complement negation of a QWORD pair in native 64 bit arithmetic.
 * @param l Low QWORD, negated in place
 * @param h High QWORD, negated in place
 */
void NegateQwords(uint64_t& l, uint64_t& h)
{
    l = ~l + 1ull;
    h = ~h + ((l == 0) ? 1ull : 0ull);
}

/**
 * @brief Reads the sign and the magnitude of a value straight from its QWORDs.
 *
 * The magnitude of -2^127 is 2^127, which needs the unsigned type: that is the one value whose
 * negation int128_t cannot represent.
 *
 * @param x Value to split
 * @param negative Receives true when x is below zero
 * @return |x|
 */
[[nodiscard]] uint128_t MagnitudeOf(const int128_t& x, bool& negative)
{
    uint64_t l = 0, h = 0;
    x.get_components(l, h);
    negative = (h >> 63) != 0;
    if (negative) {
        NegateQwords(l, h);
    }
    return uint128_t(l, h);
}

/**
 * @brief Builds an int128_t from a sign and a magnitude. A magnitude of 2^127 with a positive sign
 *        wraps to -2^127, as the type itself does.
 * @param negative True for a negative result
 * @param l Low QWORD of the magnitude
 * @param h High QWORD of the magnitude
 * @return The signed value.
 */
[[nodiscard]] int128_t FromSignAndMagnitude(bool negative, uint64_t l, uint64_t h)
{
    if (negative) {
        NegateQwords(l, h);
    }
    return int128_t(l, h);
}

/**
 * @brief The 128 bit value a signed 64 bit integer sign extends to, assembled from QWORDs.
 * @param x Value to extend
 * @return x as an int128_t.
 */
[[nodiscard]] int128_t FromInt64(int64_t x)
{
    return int128_t(static_cast<uint64_t>(x), SignExtension(x));
}

/**
 * @brief Reference truncating division: restoring shift and subtract over the magnitudes.
 *
 * One quotient bit per step, in native 64 bit arithmetic on the QWORDs, so it shares nothing with
 * div_64bit(), div_128bit() or the power of two shortcut of the operator under test. The signs are
 * applied the way the builtin types do: the quotient is negative when exactly one operand is, and the
 * remainder takes the sign of the dividend. -2^127 / -1 wraps back to -2^127.
 *
 * @param a Dividend
 * @param b Divisor, non zero
 * @param q Receives the quotient
 * @param r Receives the remainder
 */
void ReferenceDivide(const int128_t& a, const int128_t& b, int128_t& q, int128_t& r)
{
    bool aNegative = false, bNegative = false;
    uint64_t al = 0, ah = 0, bl = 0, bh = 0;
    MagnitudeOf(a, aNegative).get_components(al, ah);
    MagnitudeOf(b, bNegative).get_components(bl, bh);

    uint64_t ql = 0, qh = 0, rl = 0, rh = 0;
    for (int32_t bit = 127; bit >= 0; --bit) {
        // The partial remainder is below |b| <= 2^127 before the shift, so it still fits afterwards.
        const uint64_t next = (bit >= 64) ? ((ah >> (bit - 64)) & 1ull) : ((al >> bit) & 1ull);
        rh = (rh << 1) | (rl >> 63);
        rl = (rl << 1) | next;
        if (rh > bh || (rh == bh && rl >= bl)) {
            const uint64_t borrow = (rl < bl) ? 1ull : 0ull;
            rl -= bl;
            rh -= bh + borrow;
            if (bit >= 64) {
                qh |= 1ull << (bit - 64);
            } else {
                ql |= 1ull << bit;
            }
        }
    }
    q = FromSignAndMagnitude(aNegative != bNegative, ql, qh);
    r = FromSignAndMagnitude(aNegative, rl, rh);
}

/**
 * @brief Checks that q and r are exactly the truncated quotient and remainder of a / b.
 *
 * q * b + r == a in 128 bit arithmetic only pins q down modulo 2^128 / gcd(b, 2^128), and for this
 * type it holds for any q at all, because operator%= derives the remainder from the quotient. The
 * identity is therefore checked on the magnitudes over the integers: |q| * |b| is taken 256 bits
 * wide and must not reach past 128 bits, and adding |r| must give |a| without a carry out. The signs
 * follow truncation towards zero. Not valid for -2^127 / -1, whose quotient wraps.
 *
 * @param a Dividend
 * @param b Divisor, non zero
 * @param q Quotient to check
 * @param r Remainder to check
 * @return Success when |a| = |q| * |b| + |r| with |r| < |b| and the signs are right.
 */
[[nodiscard]] testing::AssertionResult IsExactSignedDivision(const int128_t& a, const int128_t& b, const int128_t& q, const int128_t& r)
{
    bool aNegative = false, bNegative = false, qNegative = false, rNegative = false;
    const uint128_t ma = MagnitudeOf(a, aNegative), mb = MagnitudeOf(b, bNegative);
    const uint128_t mq = MagnitudeOf(q, qNegative), mr = MagnitudeOf(r, rNegative);

    uint128_t high;
    const uint128_t low = ReferenceMultiplyWide(mq, mb, high);
    if (!high.is_zero()) {
        return testing::AssertionFailure() << "|q| * |b| overflows 128 bits, q=" << static_cast<std::string>(q);
    }

    uint64_t productLow = 0, productHigh = 0, rLow = 0, rHigh = 0, aLow = 0, aHigh = 0, bLow = 0, bHigh = 0;
    low.get_components(productLow, productHigh);
    mr.get_components(rLow, rHigh);
    ma.get_components(aLow, aHigh);
    mb.get_components(bLow, bHigh);

    // |q| * |b| + |r|, with the carries made explicit
    const uint64_t sumLow = productLow + rLow;
    const uint64_t carry = (sumLow < productLow) ? 1ull : 0ull;
    const uint64_t partial = productHigh + rHigh;
    const uint64_t sumHigh = partial + carry;
    const bool carryOut = (partial < productHigh) || (sumHigh < partial);
    if (carryOut || sumLow != aLow || sumHigh != aHigh) {
        return testing::AssertionFailure() << "|q| * |b| + |r| != |a|, q=" << static_cast<std::string>(q) << ", r=" << static_cast<std::string>(r);
    }
    if (!(rHigh < bHigh || (rHigh == bHigh && rLow < bLow))) {
        return testing::AssertionFailure() << "|r| not below |b|, r=" << static_cast<std::string>(r);
    }
    if (!mq.is_zero() && qNegative != (aNegative != bNegative)) {
        return testing::AssertionFailure() << "quotient has the wrong sign, q=" << static_cast<std::string>(q);
    }
    if (!mr.is_zero() && rNegative != aNegative) {
        return testing::AssertionFailure() << "remainder does not take the sign of the dividend, r=" << static_cast<std::string>(r);
    }
    return testing::AssertionSuccess();
}

/**
 * @brief Reference shift of a 128 bit value, built one bit at a time from its QWORDs.
 *
 * Follows the documented contract of the shift operators: a count of zero or below leaves the value
 * unchanged, the right shift is arithmetic and fills with copies of the sign bit, and a count of 128
 * or more shifts every bit out (leaving zero, or -1 for a negative value shifted right).
 *
 * @param x Value to shift
 * @param shift Bit count
 * @param left True for a left shift, false for an arithmetic right shift
 * @return The shifted value.
 */
[[nodiscard]] int128_t ReferenceShift(const int128_t& x, int32_t shift, bool left)
{
    uint64_t in[2] {};
    x.get_components(in[0], in[1]);
    if (shift <= 0) {
        return int128_t(in[0], in[1]);
    }

    const uint64_t sign = in[1] >> 63;
    uint64_t out[2] {};
    for (int32_t bit = 0; bit < 128; ++bit) {
        // 64 bit arithmetic, so that a huge count cannot overflow
        const int64_t source = left ? static_cast<int64_t>(bit) - shift : static_cast<int64_t>(bit) + shift;
        uint64_t value = 0;
        if (source > 127) {
            value = sign;
        } else if (source >= 0) {
            value = (in[source >> 6] >> (source & 63)) & 1ull;
        }
        if (value != 0) {
            out[bit >> 6] |= 1ull << (bit & 63);
        }
    }
    return int128_t(out[0], out[1]);
}

/**
 * @brief The value the double constructor has to produce, built from the IEEE fields.
 *
 * The significand and exponent come from frexp() and are placed with native 64 bit shifts, so the
 * reference shares no code with the constructor or the shift operators. Fractions truncate towards
 * zero; NaN, the infinities and every magnitude of 2^127 or more saturate to the end of the range
 * with the sign of the input, as documented.
 *
 * @param value Any double
 * @return The expected int128_t.
 */
[[nodiscard]] int128_t ExpectedFromDouble(double value)
{
    const int128_t max(~0ull, 0x7FFFFFFFFFFFFFFFull);
    const int128_t min(0ull, 0x8000000000000000ull);
    if (std::isnan(value)) {
        return std::signbit(value) ? min : max;
    }
    if (value >= 0x1p127) {
        return max;
    }
    if (value <= -0x1p127) {
        return min;
    }

    const double magnitude = std::fabs(value);
    if (magnitude < 1.0) {
        return int128_t();
    }
    int exponent = 0;
    const double fraction = std::frexp(magnitude, &exponent);  // magnitude = fraction * 2^exponent, fraction in [0.5, 1)
    const uint64_t significand = static_cast<uint64_t>(std::ldexp(fraction, 53));
    const int32_t shift = exponent - 53;  // magnitude = significand * 2^shift, and shift >= -52
    uint64_t l = 0, h = 0;
    if (shift <= 0) {
        l = significand >> -shift;  // drops the fraction, truncating the magnitude
    } else if (shift < 64) {
        l = significand << shift;
        h = significand >> (64 - shift);
    } else {
        h = significand << (shift - 64);
    }
    if (std::signbit(value)) {
        NegateQwords(l, h);
    }
    return int128_t(l, h);
}

/**
 * @brief The integer powers of ten that fit in an int128_t, 10^0 to 10^38.
 *
 * Built with ReferenceMultiply() rather than the operator under test, so log10() is checked against
 * values it did not compute itself.
 *
 * @return The table, index k holding 10^k.
 */
[[nodiscard]] std::array<int128_t, 39> PowersOfTen()
{
    std::array<int128_t, 39> table;
    table[0] = int128_t(1ull, 0ull);
    for (size_t k = 1; k < table.size(); ++k) {
        table[k] = ReferenceMultiply(table[k - 1], int128_t(10ull, 0ull));
    }
    return table;
}

/**
 * @brief Orders two non negative values by their QWORDs, independently of the comparison operators.
 * @param a Left operand, not negative
 * @param b Right operand, not negative
 * @return True when a <= b.
 */
[[nodiscard]] bool NonNegativeLessOrEqual(const int128_t& a, const int128_t& b)
{
    uint64_t al = 0, ah = 0, bl = 0, bh = 0;
    a.get_components(al, ah);
    b.get_components(bl, bh);
    return ah < bh || (ah == bh && al <= bl);
}

/// @brief Satisfied when a * b compiles. A concept rather than a bare requires-expression, which
///        MSVC rejects outside a template.
template <class A, class B>
concept Multipliable = requires(A a, B b) { a * b; };

/// @brief Satisfied when a *= b compiles.
template <class A, class B>
concept MultiplyAssignable = requires(A& a, B b) { a *= b; };
}  // namespace

/**********************************************************************
 * int128_t tests
 ***********************************************************************/
// Construct fixed_point128 and convert back to/from various elements.
// Both QWORDs are checked: reading the value back through int64_t only ever looked at the low one.
TEST(int128_t, DefaultConstructor)
{
    const int128_t i;
    uint64_t low = 1, high = 1;
    i.get_components(low, high);
    EXPECT_EQ(low, 0ull);
    EXPECT_EQ(high, 0ull);
    EXPECT_TRUE(i.is_zero());
}
// The previous version converted the result back to a double and compared with EXPECT_DOUBLE_EQ, a
// round trip with a 4 ulp tolerance that cannot see any bit below the 53 a double holds: a
// constructor leaving junk in bit 0 of every value above 2^64 passed it. Its operands were also
// floor()ed doubles in [2, 2^127), so fractions, magnitudes below one, zero and the saturation were
// never reached. The result is now compared bit for bit with a value built from the IEEE fields, on
// operands that keep their fraction.
TEST(int128_t, ConstructorFromDouble)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const double value = get_double_random(-3, 127);
        const int128_t f = value;
        EXPECT_TRUE(f == ExpectedFromDouble(value)) << "value=" << std::hexfloat << value << ", f=" << (std::string)f;
    }
}
TEST(int128_t, ConstructorFromFloat)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        // rounding to float can carry a value just below 2^127 up to 2^127, which saturates
        const float value = static_cast<float>(get_double_random(-3, 127));
        const int128_t f = value;
        EXPECT_TRUE(f == ExpectedFromDouble(value)) << "value=" << std::hexfloat << value << ", f=" << (std::string)f;
    }
}
// Reading the value back through int32_t only checked the low 32 bits, so a constructor that zero
// extended a negative operand passed. The whole 128 bit pattern is compared now.
TEST(int128_t, ConstructorFromInt32)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const int32_t value = get_int32_random();
        const int128_t f = value;
        EXPECT_EQ(static_cast<int32_t>(f), value);
        EXPECT_TRUE(f == FromInt64(value)) << "value=" << value;
        EXPECT_EQ(f.is_negative(), value < 0) << "value=" << value;
    }
}
TEST(int128_t, ConstructorFromUnsignedInt32)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const uint32_t value = get_uint32_random();
        const int128_t f = value;
        EXPECT_EQ(static_cast<uint32_t>(f), value);
        // zero extended, a value with its msb set stays positive
        EXPECT_TRUE(f == int128_t(static_cast<uint64_t>(value), 0ull)) << "value=" << value;
    }
}
TEST(int128_t, ConstructorFromInt64)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const int64_t value = get_int64_random();
        const int128_t f = value;
        EXPECT_EQ(static_cast<int64_t>(f), value);
        EXPECT_TRUE(f == FromInt64(value)) << "value=" << value;
        EXPECT_EQ(f.is_negative(), value < 0) << "value=" << value;
    }
}
// The previous version used get_int64_random() and int64_t throughout, making it a verbatim copy
// of ConstructorFromInt64 that never once fed the constructor an unsigned operand.
TEST(int128_t, ConstructorFromUnsignedInt64)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const uint64_t value = get_uint64_random();
        const int128_t f = value;

        EXPECT_EQ(static_cast<uint64_t>(f), value);
        // an unsigned operand is zero extended rather than sign extended, so a value with its msb
        // set becomes a positive 128 bit number rather than a negative one
        EXPECT_TRUE(f == int128_t(value, 0ull)) << "value=" << value;
        EXPECT_TRUE(f.is_positive()) << "value=" << value;
    }
}
// The previous version parsed the input, printed it back, and compared the two through strtod. That
// is a round trip, so a parse and a print that were wrong in opposite directions cancelled out, and
// the double funnel discarded everything below 53 bits: "123456789012345678" was compared as
// 123456789012345680. Values are now checked against what they are supposed to be.
TEST(int128_t, ConstructorFromString)
{
    EXPECT_TRUE(int128_t("0xDEADBEAF") == int128_t(0xDEADBEAFll));
    EXPECT_TRUE(int128_t("1234") == int128_t(1234ll));
    EXPECT_TRUE(int128_t("-1234") == int128_t(-1234ll));
    EXPECT_TRUE(int128_t("123456789012345678") == int128_t(123456789012345678ll));
    EXPECT_TRUE(int128_t("-123456789012345678") == int128_t(-123456789012345678ll));
    // values needing more than 64 bit, where a double reference has no chance
    EXPECT_TRUE(int128_t("18446744073709551616") == int128_t(0ull, 1ull));                   // 2^64
    EXPECT_TRUE(int128_t("170141183460469231731687303715884105727") == int128_t(~0ull, 0x7FFFFFFFFFFFFFFFull));
    EXPECT_TRUE(int128_t("-170141183460469231731687303715884105728") == int128_t(0ull, 0x8000000000000000ull));
    EXPECT_TRUE(int128_t("0x7FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF") == int128_t(~0ull, 0x7FFFFFFFFFFFFFFFull));
    // a nullptr, an empty string and a leading illegal character all produce zero
    EXPECT_TRUE(int128_t(static_cast<const char*>(nullptr)).is_zero());
    EXPECT_TRUE(int128_t("").is_zero());
    EXPECT_TRUE(int128_t("abc").is_zero());
    // parsing stops at the first character that is not a digit of the detected base
    EXPECT_TRUE(int128_t("12abc") == int128_t(12ll));
    EXPECT_TRUE(int128_t("  -42  ") == int128_t(-42ll));

    // every value round trips through its own decimal string, both signs, full width
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const uint64_t low = get_uint64_random();
        const uint64_t high = get_uint64_random();
        const int128_t v(low, high);
        const std::string s = (std::string)v;
        EXPECT_TRUE(int128_t(s) == v) << "s=" << s;
    }
}
// The copy and move tests used to compare the two objects through EXPECT_DOUBLE_EQ of their double
// conversions, which cannot see anything below the top 53 bits. Every QWORD is compared now, against
// the pattern the source was built from rather than against the source alone.
TEST(int128_t, CopyConstructor)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const uint64_t low = get_uint64_random();
        const uint64_t high = get_uint64_random();
        const int128_t f1(low, high);
        const int128_t f2 = f1;
        EXPECT_TRUE(f2 == int128_t(low, high)) << "f1=" << (std::string)f1;
    }
}
TEST(int128_t, MoveConstructor)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const uint64_t low = get_uint64_random();
        const uint64_t high = get_uint64_random();
        int128_t f1(low, high);
        const int128_t f2(std::move(f1));
        EXPECT_TRUE(f2 == int128_t(low, high)) << "f2=" << (std::string)f2;
        // documented to leave the source untouched
        EXPECT_TRUE(f1 == int128_t(low, high)) << "f1=" << (std::string)f1;
    }
}
TEST(int128_t, AssignmentOperator)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const uint64_t low = get_uint64_random();
        const uint64_t high = get_uint64_random();
        const int128_t f1(low, high);
        int128_t f2(~low, ~high);
        f2 = f1;
        EXPECT_TRUE(f2 == int128_t(low, high)) << "f1=" << (std::string)f1;
    }
}
TEST(int128_t, MoveAssignmentOperator)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const uint64_t low = get_uint64_random();
        const uint64_t high = get_uint64_random();
        int128_t f1(low, high);
        int128_t f2(~low, ~high);
        f2 = std::move(f1);
        EXPECT_TRUE(f2 == int128_t(low, high)) << "f2=" << (std::string)f2;
        EXPECT_TRUE(f1 == int128_t(low, high)) << "f1=" << (std::string)f1;
    }
}
// Operands are raw QWORD pairs rather than values converted from a double. A double carries 53
// bits, so the previous version could only ever see the top 53 bits of a 128 bit result, and it
// compared them through is_similar_double(), whose 1e-10 tolerance every iteration satisfied. The
// EXPECT below it never executed once in 65536 iterations.
//
// Addition is bit identical for both instantiations: a two's complement sum is a sum modulo 2^128,
// so the reference is the same unsigned limb arithmetic the uint128_t suite uses.
TEST(int128_t, Add)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const uint64_t al = get_uint64_random(), ah = get_uint64_random();
        const uint64_t bl = get_uint64_random(), bh = get_uint64_random();
        const int128_t a(al, ah), b(bl, bh);

        // the reference carry out of the low QWORD: an unsigned sum that wrapped is smaller than
        // either operand. Overflow out of the high QWORD is dropped, as the operator does.
        const uint64_t low = al + bl;
        const uint64_t high = ah + bh + ((low < al) ? 1ull : 0ull);
        const int128_t expected(low, high);

        EXPECT_TRUE(a + b == expected) << "a=" << (std::string)a << ", b=" << (std::string)b;

        int128_t c = a;
        c += b;
        EXPECT_TRUE(c == expected) << "a=" << (std::string)a << ", b=" << (std::string)b;
        int128_t d = a;
        d += d;
        EXPECT_TRUE(d == a + a) << "a=" << (std::string)a;

        EXPECT_TRUE(a + b == b + a) << "a=" << (std::string)a << ", b=" << (std::string)b;
        EXPECT_TRUE((a + b) - b == a) << "a=" << (std::string)a << ", b=" << (std::string)b;
        EXPECT_TRUE(a + int128_t() == a) << "a=" << (std::string)a;
        EXPECT_TRUE(a + (-a) == int128_t()) << "a=" << (std::string)a;
    }
}
TEST(int128_t, AddEdgeCases)
{
    const int128_t zero, one(1ll), minus_one(-1ll);
    const int128_t max(~0ull, 0x7FFFFFFFFFFFFFFFull);   // 2^127 - 1
    const int128_t min(0ull, 0x8000000000000000ull);    // -2^127

    // the carry has to cross the QWORD boundary in both directions
    EXPECT_TRUE(int128_t(~0ull, 0ull) + one == int128_t(0ull, 1ull));
    EXPECT_TRUE(minus_one + one == zero);
    EXPECT_TRUE(zero + minus_one == minus_one);
    // signed overflow wraps silently, matching the builtin types
    EXPECT_TRUE(max + one == min);
    EXPECT_TRUE(min + minus_one == max);
    EXPECT_TRUE(max + max == int128_t(~0ull - 1ull, 0xFFFFFFFFFFFFFFFFull));  // -2
    // adding the most negative value to itself wraps to zero
    EXPECT_TRUE(min + min == zero);
}
// The reference int64_t res = value1 + value2 overflows on a quarter of the operand pairs, and
// -get_int64_random() is itself undefined when the draw comes back INT64_MIN. The operands are now
// full width 128 bit values of opposite sign, with magnitudes kept below 2^127 so that negating
// them is exact.
TEST(int128_t, AddDifferentSign)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        // one draw per statement: the order in which a call's arguments are evaluated is unspecified,
        // and MSVC and clang-cl would otherwise test different operands
        const uint64_t positiveLow = get_uint64_random();
        const uint64_t positiveHigh = get_uint64_random() >> 1;
        const uint64_t magnitudeLow = get_uint64_random();
        const uint64_t magnitudeHigh = get_uint64_random() >> 1;
        const int128_t positive(positiveLow, positiveHigh);
        const int128_t magnitude(magnitudeLow, magnitudeHigh);
        const int128_t negative = -magnitude;

        uint64_t pl = 0, ph = 0, nl = 0, nh = 0;
        positive.get_components(pl, ph);
        negative.get_components(nl, nh);
        const uint64_t low = pl + nl;
        const uint64_t high = ph + nh + ((low < pl) ? 1ull : 0ull);

        const int128_t sum = positive + negative;
        EXPECT_TRUE(sum == int128_t(low, high)) << "p=" << (std::string)positive << ", n=" << (std::string)negative;
        // adding a negative value is subtracting its magnitude
        EXPECT_TRUE(sum == positive - magnitude) << "p=" << (std::string)positive << ", n=" << (std::string)negative;
        // and the result is negative exactly when the negative operand has the larger magnitude
        EXPECT_EQ(sum.is_negative(), magnitude > positive) << "p=" << (std::string)positive << ", n=" << (std::string)negative;
    }
}
// The reference used to be computed as value1 + value2 in int64_t, which overflows on roughly a
// quarter of the operand pairs (16281 of 65536 measured) and is undefined behaviour rather than the
// wrap the test was counting on. Both operands are sign extended to 128 bit before being added, so
// the exact result needs no wrapping at all and the full width is worth checking.
TEST(int128_t, AddInt64)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const int64_t value1 = get_int64_random();
        const int64_t value2 = get_int64_random();

        const uint64_t u1 = static_cast<uint64_t>(value1), u2 = static_cast<uint64_t>(value2);
        const uint64_t low = u1 + u2;
        const uint64_t high = SignExtension(value1) + SignExtension(value2) + ((low < u1) ? 1ull : 0ull);
        const int128_t expected(low, high);

        const int128_t f1 = value1;
        EXPECT_TRUE(f1 + value2 == expected) << "value1=" << value1 << ", value2=" << value2;
        EXPECT_TRUE(f1 + value2 == f1 + int128_t(value2)) << "value1=" << value1 << ", value2=" << value2;

        int128_t f2 = f1;
        f2 += value2;
        EXPECT_TRUE(f2 == expected) << "value1=" << value1 << ", value2=" << value2;

        // The sum of two int64_t always fits in 128 bit, so no wrap can occur and the sign of the
        // result is the mathematical one. Derived without overflowing: only a mixed sign pair can
        // change sign, and then the larger magnitude decides. Magnitudes are taken in unsigned
        // arithmetic, where negating INT64_MIN is well defined.
        bool expected_negative;
        if (value1 < 0 && value2 < 0) {
            expected_negative = true;
        } else if (value1 >= 0 && value2 >= 0) {
            expected_negative = false;
        } else {
            const uint64_t m1 = (value1 < 0) ? (0ull - u1) : u1;
            const uint64_t m2 = (value2 < 0) ? (0ull - u2) : u2;
            expected_negative = (value1 < 0) ? (m1 > m2) : (m2 > m1);
        }
        EXPECT_EQ(expected.is_negative(), expected_negative) << "value1=" << value1 << ", value2=" << value2;
    }
}
// The previous version used get_int64_random() and int64_t throughout, so it was a verbatim copy of
// AddInt64 and never fed the operator an unsigned scalar. An unsigned operand is zero extended, so
// a value with its msb set adds a large positive number rather than a small negative one.
TEST(int128_t, AddUnsignedInt64)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const int64_t value1 = get_int64_random();
        const uint64_t value2 = get_uint64_random();

        const uint64_t u1 = static_cast<uint64_t>(value1);
        const uint64_t low = u1 + value2;
        const uint64_t high = SignExtension(value1) + ((low < u1) ? 1ull : 0ull);  // value2 zero extends
        const int128_t expected(low, high);

        const int128_t f1 = value1;
        EXPECT_TRUE(f1 + value2 == expected) << "value1=" << value1 << ", value2=" << value2;
        EXPECT_TRUE(f1 + value2 == f1 + int128_t(value2)) << "value1=" << value1 << ", value2=" << value2;

        int128_t f2 = f1;
        f2 += value2;
        EXPECT_TRUE(f2 == expected) << "value1=" << value1 << ", value2=" << value2;
    }
}

// The previous version derived the subtrahend as floor(value1 / 2.5), so the two operands always
// shared a sign and the result never crossed zero. See Add for the rest of the reasoning.
TEST(int128_t, Subtract)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const uint64_t al = get_uint64_random(), ah = get_uint64_random();
        const uint64_t bl = get_uint64_random(), bh = get_uint64_random();
        const int128_t a(al, ah), b(bl, bh);

        const uint64_t low = al - bl;
        const uint64_t high = ah - bh - ((al < bl) ? 1ull : 0ull);
        const int128_t expected(low, high);

        EXPECT_TRUE(a - b == expected) << "a=" << (std::string)a << ", b=" << (std::string)b;

        int128_t c = a;
        c -= b;
        EXPECT_TRUE(c == expected) << "a=" << (std::string)a << ", b=" << (std::string)b;
        int128_t d = a;
        d -= d;
        EXPECT_TRUE(d.is_zero()) << "a=" << (std::string)a;

        EXPECT_TRUE((a - b) + b == a) << "a=" << (std::string)a << ", b=" << (std::string)b;
        EXPECT_TRUE(a - b == -(b - a)) << "a=" << (std::string)a << ", b=" << (std::string)b;
        EXPECT_TRUE(a - int128_t() == a) << "a=" << (std::string)a;
        EXPECT_TRUE(int128_t() - a == -a) << "a=" << (std::string)a;
        // subtracting is adding the negation, which is what the implementation actually does
        EXPECT_TRUE(a - b == a + (-b)) << "a=" << (std::string)a << ", b=" << (std::string)b;
    }
}
TEST(int128_t, SubtractEdgeCases)
{
    const int128_t zero, one(1ll), minus_one(-1ll);
    const int128_t max(~0ull, 0x7FFFFFFFFFFFFFFFull);   // 2^127 - 1
    const int128_t min(0ull, 0x8000000000000000ull);    // -2^127

    EXPECT_TRUE(zero - one == minus_one);
    EXPECT_TRUE(zero - minus_one == one);
    // the borrow has to cross the QWORD boundary
    EXPECT_TRUE(int128_t(0ull, 1ull) - one == int128_t(~0ull, 0ull));
    // signed underflow wraps silently
    EXPECT_TRUE(min - one == max);
    EXPECT_TRUE(max - minus_one == min);
    // the most negative value is its own negation, so subtracting it is the same as adding it
    EXPECT_TRUE(zero - min == min);
    EXPECT_TRUE(min - min == zero);
}
// Same undefined signed overflow as AddDifferentSign, plus the undefined negation of INT64_MIN.
TEST(int128_t, SubtractDifferentSign)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        // one draw per statement: the order in which a call's arguments are evaluated is unspecified,
        // and MSVC and clang-cl would otherwise test different operands
        const uint64_t positiveLow = get_uint64_random();
        const uint64_t positiveHigh = get_uint64_random() >> 1;
        const uint64_t magnitudeLow = get_uint64_random();
        const uint64_t magnitudeHigh = get_uint64_random() >> 1;
        const int128_t positive(positiveLow, positiveHigh);
        const int128_t magnitude(magnitudeLow, magnitudeHigh);
        const int128_t negative = -magnitude;

        uint64_t pl = 0, ph = 0, nl = 0, nh = 0;
        positive.get_components(pl, ph);
        negative.get_components(nl, nh);
        const uint64_t low = pl - nl;
        const uint64_t high = ph - nh - ((pl < nl) ? 1ull : 0ull);

        const int128_t difference = positive - negative;
        EXPECT_TRUE(difference == int128_t(low, high)) << "p=" << (std::string)positive << ", n=" << (std::string)negative;
        // subtracting a negative value adds its magnitude, so the result cannot be smaller
        EXPECT_TRUE(difference == positive + magnitude) << "p=" << (std::string)positive << ", n=" << (std::string)negative;
    }
}
// The reference value1 - value2 in int64_t overflows on roughly a quarter of the operand pairs
// (16427 of 65536 measured), which is undefined behaviour. See AddInt64.
TEST(int128_t, SubtractInt64)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const int64_t value1 = get_int64_random();
        const int64_t value2 = get_int64_random();

        const uint64_t u1 = static_cast<uint64_t>(value1), u2 = static_cast<uint64_t>(value2);
        const uint64_t low = u1 - u2;
        const uint64_t high = SignExtension(value1) - SignExtension(value2) - ((u1 < u2) ? 1ull : 0ull);
        const int128_t expected(low, high);

        const int128_t f1 = value1;
        EXPECT_TRUE(f1 - value2 == expected) << "value1=" << value1 << ", value2=" << value2;
        EXPECT_TRUE(f1 - value2 == f1 - int128_t(value2)) << "value1=" << value1 << ", value2=" << value2;

        int128_t f2 = f1;
        f2 -= value2;
        EXPECT_TRUE(f2 == expected) << "value1=" << value1 << ", value2=" << value2;
    }
}
// The previous version used get_int64_random() and int64_t throughout, so it was a verbatim copy of
// SubtractInt64 and never fed the operator an unsigned scalar.
TEST(int128_t, SubtractUnsignedInt64)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const int64_t value1 = get_int64_random();
        const uint64_t value2 = get_uint64_random();

        const uint64_t u1 = static_cast<uint64_t>(value1);
        const uint64_t low = u1 - value2;
        const uint64_t high = SignExtension(value1) - ((u1 < value2) ? 1ull : 0ull);  // value2 zero extends
        const int128_t expected(low, high);

        const int128_t f1 = value1;
        EXPECT_TRUE(f1 - value2 == expected) << "value1=" << value1 << ", value2=" << value2;
        EXPECT_TRUE(f1 - value2 == f1 - int128_t(value2)) << "value1=" << value1 << ", value2=" << value2;

        int128_t f2 = f1;
        f2 -= value2;
        EXPECT_TRUE(f2 == expected) << "value1=" << value1 << ", value2=" << value2;
    }
}
// Checked against ReferenceMultiply() rather than a double product. Besides losing everything below
// the top 53 bits, the old oracle had a guard that did not work: log2(res) is NaN for a negative
// product and NaN > 127 is false, so half the iterations (32768 of 65536 measured) sailed past the
// overflow check it was supposed to enforce. Every iteration was then skipped by is_similar_double()
// anyway, so the assertion never ran.
TEST(int128_t, MultiplyByint128)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        // shifting the high QWORDs by independent random amounts mixes products that overflow with
        // products that stay in range, and both operands take either sign
        const uint64_t aLow = get_uint64_random();
        const uint64_t aHigh = get_uint64_random();
        const uint32_t aShift = get_uint32_random() % 64;
        const uint64_t bLow = get_uint64_random();
        const uint64_t bHigh = get_uint64_random();
        const uint32_t bShift = get_uint32_random() % 64;
        int128_t a(aLow, aHigh >> aShift);
        int128_t b(bLow, bHigh >> bShift);
        if (get_uint32_random() & 1)
            a = -a;
        if (get_uint32_random() & 1)
            b = -b;

        const int128_t expected = ReferenceMultiply(a, b);
        EXPECT_TRUE(a * b == expected) << "a=" << (std::string)a << ", b=" << (std::string)b;

        int128_t c = a;
        c *= b;
        EXPECT_TRUE(c == expected) << "a=" << (std::string)a << ", b=" << (std::string)b;

        // identities that hold in modulo 2^128 arithmetic, so they survive truncation
        EXPECT_TRUE(a * b == b * a) << "a=" << (std::string)a << ", b=" << (std::string)b;
        EXPECT_TRUE(a * int128_t(1ll) == a) << "a=" << (std::string)a;
        EXPECT_TRUE((a * int128_t()).is_zero()) << "a=" << (std::string)a;
        EXPECT_TRUE(a * int128_t(-1ll) == -a) << "a=" << (std::string)a;
        EXPECT_TRUE(a * int128_t(2ll) == a + a) << "a=" << (std::string)a;
        // the sign rules, which the magnitude free implementation has to reproduce
        EXPECT_TRUE((-a) * b == -(a * b)) << "a=" << (std::string)a << ", b=" << (std::string)b;
        EXPECT_TRUE((-a) * (-b) == a * b) << "a=" << (std::string)a << ", b=" << (std::string)b;
        EXPECT_TRUE(a * (b + int128_t(1ll)) == a * b + a) << "a=" << (std::string)a << ", b=" << (std::string)b;
    }
}
TEST(int128_t, MultiplyByint128EdgeCases)
{
    const int128_t zero, one(1ll), minus_one(-1ll);
    const int128_t max(~0ull, 0x7FFFFFFFFFFFFFFFull);   // 2^127 - 1
    const int128_t min(0ull, 0x8000000000000000ull);    // -2^127

    EXPECT_TRUE(zero * zero == zero);
    EXPECT_TRUE(max * one == max);
    EXPECT_TRUE(max * minus_one == -max);
    // the most negative value has no positive counterpart, so negating it wraps back to itself
    EXPECT_TRUE(min * one == min);
    EXPECT_TRUE(min * minus_one == min);
    EXPECT_TRUE((min * int128_t(2ll)).is_zero());
    // (2^127-1)^2 is 2^254 - 2^128 + 1, whose low 128 bit are 1
    EXPECT_TRUE(max * max == one);
    EXPECT_TRUE(minus_one * minus_one == one);

    // the reference has to agree on every pairing of the extremes
    const int128_t values[] = {zero, one, minus_one, max, min, int128_t(0ull, 1ull), int128_t(~0ull, 0ull)};
    for (const auto& x : values) {
        for (const auto& y : values) {
            EXPECT_TRUE(x * y == ReferenceMultiply(x, y)) << "x=" << (std::string)x << ", y=" << (std::string)y;
        }
    }
}
// sqr(x) must be bit identical to x * x. square() skips the absolute value conversions that
// operator*= performs, which is valid because (-x)^2 == x^2 modulo 2^128. Raw bit patterns
// are used so both signs and the overflowing values are covered.
TEST(int128_t, SqrMatchesMultiply)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const uint64_t low = get_uint64_random();
        const uint64_t high = get_uint64_random();
        const int128_t x(low, high);
        EXPECT_TRUE(sqr(x) == x * x) << "low=" << low << ", high=" << high;
    }
}
TEST(int128_t, SqrEdgeCases)
{
    const int128_t values[] = {int128_t(0ull, 0ull),  int128_t(1ull, 0ull),   int128_t(0ull, 1ull),
                               int128_t(~0ull, ~0ull),  // -1
                               int128_t(0ull, 0x8000000000000000ull),  // most negative value
                               int128_t(~0ull, 0x7FFFFFFFFFFFFFFFull)};  // most positive value
    for (const auto& x : values) {
        EXPECT_TRUE(sqr(x) == x * x) << "x=" << (std::string)x;
    }
    // small values of both signs whose square is exact
    for (int64_t v = -1000; v <= 1000; ++v) {
        const int128_t x(v);
        EXPECT_TRUE(sqr(x) == int128_t(v * v)) << "v=" << v;
        EXPECT_TRUE(sqr(x) == x * x) << "v=" << v;
    }
}
// trunc(value1 / value2) computed in double is not a usable reference: the quotient of two 128 bit
// values needs more than the 53 bits a double holds, so the old oracle was wrong rather than the
// implementation, and is_similar_double() hid the difference on all 65536 iterations.
//
// With q*b + r == a, |r| < |b| and r taking the sign of a, the pair is unique - but only when the
// identity holds over the integers. Evaluated in 128 bit arithmetic it holds for ANY quotient here,
// because operator%= derives the remainder from the quotient, and a quotient off by 2^126 with a
// divisor divisible by four kept every other check green as well. IsExactSignedDivision() checks the
// identity on the magnitudes without the wrap, and ReferenceDivide() supplies the exact pair.
// Magnitudes stay below 2^127 so that negation and abs() are exact.
TEST(int128_t, DivideByint128)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const uint64_t aLow = get_uint64_random();
        const uint64_t aHigh = get_uint64_random() >> 1;
        int128_t a(aLow, aHigh);

        // spread the divisor over the whole magnitude range so all three paths in operator/= are
        // reached: the power of two shift, the 64 bit divisor shortcut and the 128 bit long division
        const int32_t shift = static_cast<int32_t>(get_uint32_random() % 127);
        const uint64_t bLow = get_uint64_random();
        const uint64_t bHigh = get_uint64_random() >> 1;
        int128_t b = int128_t(bLow, bHigh) >> shift;
        if ((get_uint32_random() & 7) == 0)
            b = int128_t(1ll) << static_cast<int32_t>(get_uint32_random() % 126);
        if (b.is_zero())
            continue;
        if (get_uint32_random() & 1)
            a = -a;
        if (get_uint32_random() & 1)
            b = -b;

        const int128_t q = a / b;
        const int128_t r = a % b;

        EXPECT_TRUE(IsExactSignedDivision(a, b, q, r)) << "a=" << (std::string)a << ", b=" << (std::string)b;
        int128_t expectedQ, expectedR;
        ReferenceDivide(a, b, expectedQ, expectedR);
        EXPECT_TRUE(q == expectedQ) << "a=" << (std::string)a << ", b=" << (std::string)b << ", q=" << (std::string)q;
        EXPECT_TRUE(r == expectedR) << "a=" << (std::string)a << ", b=" << (std::string)b << ", r=" << (std::string)r;
        EXPECT_TRUE(abs(r) < abs(b)) << "a=" << (std::string)a << ", b=" << (std::string)b;

        int128_t c = a;
        c /= b;
        EXPECT_TRUE(c == q) << "a=" << (std::string)a << ", b=" << (std::string)b;
        int128_t d = a;
        d %= b;
        EXPECT_TRUE(d == r) << "a=" << (std::string)a << ", b=" << (std::string)b;

        // the magnitude of the quotient does not depend on the signs
        EXPECT_TRUE(abs(q) == abs((-a) / b)) << "a=" << (std::string)a << ", b=" << (std::string)b;
    }
}
TEST(int128_t, DivideByint128EdgeCases)
{
    const int128_t one(1ll), minus_one(-1ll);
    const int128_t max(~0ull, 0x7FFFFFFFFFFFFFFFull);   // 2^127 - 1
    const int128_t min(0ull, 0x8000000000000000ull);    // -2^127

    EXPECT_TRUE(max / max == one);
    EXPECT_TRUE(max / one == max);
    EXPECT_TRUE(max / minus_one == -max);
    EXPECT_TRUE(min / one == min);
    EXPECT_TRUE(min / min == one);

    // -2^127 / -1 is 2^127, which is one past the top of the range. The builtin types make this
    // undefined; here it wraps back to the most negative value, and the remainder stays zero.
    EXPECT_TRUE(min / minus_one == min);
    EXPECT_TRUE((min % minus_one).is_zero());
    // the division identity survives it
    EXPECT_TRUE((min / minus_one) * minus_one + (min % minus_one) == min);

    EXPECT_TRUE(min / int128_t(2ll) == int128_t(0ull, 0xC000000000000000ull));  // -2^126
    EXPECT_TRUE((int128_t() / max).is_zero());
    EXPECT_TRUE((one / max).is_zero());
    EXPECT_TRUE(one % max == one);

    // division by zero is reported, not silently wrong
    EXPECT_THROW((void)(max / int128_t()), std::logic_error);
    EXPECT_THROW((void)(max % int128_t()), std::logic_error);
    EXPECT_THROW((void)(max / 0ll), std::logic_error);
}
// The assertion was guarded by "if (int128_res == res) continue;", so it only ever ran when it was
// about to fail and the test executed zero assertions on a passing run. The operands were also a
// 64 bit numerator and a 32 bit divisor, and value1 % value2 is undefined for INT64_MIN % -1.
TEST(int128_t, ModuloByint128)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const int64_t value1 = get_int64_random();
        const int32_t value2 = get_int32_random();
        if (value2 == 0)
            continue;
        // INT64_MIN % -1 is undefined for the builtin types, so it cannot serve as a reference.
        // The 128 bit behaviour at that point is pinned by DivideByint128EdgeCases instead.
        if (value1 == INT64_MIN && value2 == -1)
            continue;

        const int128_t f1 = value1, f2 = value2;
        EXPECT_EQ(static_cast<int64_t>(f1 % f2), value1 % value2) << "value1=" << value1 << ", value2=" << value2;
        EXPECT_EQ(static_cast<int64_t>(f1 / f2), value1 / value2) << "value1=" << value1 << ", value2=" << value2;
        // the 128 bit results have to be the sign extension of the 64 bit ones, not merely equal in
        // their low QWORD: an error in the high QWORD of the quotient is invisible to the two above
        EXPECT_TRUE(f1 % f2 == FromInt64(value1 % value2)) << "value1=" << value1 << ", value2=" << value2;
        EXPECT_TRUE(f1 / f2 == FromInt64(value1 / value2)) << "value1=" << value1 << ", value2=" << value2;
        EXPECT_TRUE(IsExactSignedDivision(f1, f2, f1 / f2, f1 % f2)) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(int128_t, Compareint128)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = floor(get_double_random(1, 127));
        double value2 = floor(get_double_random(1, 127));
        bool res = value1 > value2;
        int128_t f1 = value1;
        int128_t f2 = value2;
        if (check_overflow_int128(value1) || check_overflow_int128(value2))
            continue;

        bool int128_res = f1 > f2;
        EXPECT_TRUE(int128_res == res) << "operator>: "
                                       << "value1=" << value1 << ", value2=" << value2;

        res = value1 >= value2;
        int128_res = f1 >= f2;
        EXPECT_TRUE(int128_res == res) << "operator>=: "
                                       << "value1=" << value1 << ", value2=" << value2;

        res = value1 < value2;
        int128_res = f1 < f2;
        EXPECT_TRUE(int128_res == res) << "operator<: "
                                       << "value1=" << value1 << ", value2=" << value2;

        res = value1 <= value2;
        int128_res = f1 <= f2;
        EXPECT_TRUE(int128_res == res) << "operator<=: "
                                       << "value1=" << value1 << ", value2=" << value2;
    }
}
// Operands were 64 bit values, so no iteration ever carried into the high QWORD, and the postfix
// form's return value was never looked at. value1 + 1 is also undefined at INT64_MAX.
TEST(int128_t, OperatorPlusPlus)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const uint64_t low = get_uint64_random();
        const uint64_t high = get_uint64_random();
        const int128_t f1(low, high);
        const int128_t expected = f1 + int128_t(1ll);

        int128_t f2 = f1;
        // the postfix form returns the value from before the increment
        EXPECT_TRUE(f2++ == f1) << "f1=" << (std::string)f1;
        EXPECT_TRUE(f2 == expected) << "f1=" << (std::string)f1;

        f2 = f1;
        // the prefix form returns the value after it, as a reference to the object itself
        EXPECT_TRUE(++f2 == expected) << "f1=" << (std::string)f1;
        EXPECT_TRUE(f2 == expected) << "f1=" << (std::string)f1;
    }
}
TEST(int128_t, OperatorMinusMinus)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const uint64_t low = get_uint64_random();
        const uint64_t high = get_uint64_random();
        const int128_t f1(low, high);
        const int128_t expected = f1 - int128_t(1ll);

        int128_t f2 = f1;
        EXPECT_TRUE(f2-- == f1) << "f1=" << (std::string)f1;
        EXPECT_TRUE(f2 == expected) << "f1=" << (std::string)f1;

        f2 = f1;
        EXPECT_TRUE(--f2 == expected) << "f1=" << (std::string)f1;
        EXPECT_TRUE(f2 == expected) << "f1=" << (std::string)f1;
    }
}
TEST(int128_t, IncrementDecrementBoundaries)
{
    const int128_t max(~0ull, 0x7FFFFFFFFFFFFFFFull);   // 2^127 - 1
    const int128_t min(0ull, 0x8000000000000000ull);    // -2^127

    // the carry has to cross the QWORD boundary
    int128_t a(~0ull, 0ull);  // 2^64 - 1
    ++a;
    EXPECT_TRUE(a == int128_t(0ull, 1ull));
    int128_t b(0ull, 1ull);   // 2^64
    --b;
    EXPECT_TRUE(b == int128_t(~0ull, 0ull));

    // and crossing zero flips every bit of the high QWORD
    int128_t c(-1ll);
    ++c;
    EXPECT_TRUE(c.is_zero());
    int128_t d;
    --d;
    EXPECT_TRUE(d == int128_t(-1ll));

    // overflow wraps silently, as it does for the builtin types
    int128_t e = max;
    ++e;
    EXPECT_TRUE(e == min);
    int128_t f = min;
    --f;
    EXPECT_TRUE(f == max);
}
// The previous version asserted f1 == f1, which is true of any implementation that manages to read
// its own members twice, and never compared two distinct objects. In particular it could not tell a
// comparison of the whole 128 bit value from one that only looks at a single QWORD.
TEST(int128_t, OperatorEqual)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const uint64_t al = get_uint64_random(), ah = get_uint64_random();
        const uint64_t bl = get_uint64_random(), bh = get_uint64_random();
        const int128_t a(al, ah), b(bl, bh);

        // equality is bit identical for both instantiations, no sign test is involved
        EXPECT_EQ(a == b, (al == bl) && (ah == bh)) << "a=" << (std::string)a << ", b=" << (std::string)b;

        const int128_t copy = a;
        EXPECT_TRUE(a == copy) << "a=" << (std::string)a;
        // a value differing in exactly one QWORD must not compare equal, in either position
        EXPECT_FALSE(a == int128_t(al ^ 1ull, ah)) << "a=" << (std::string)a;
        EXPECT_FALSE(a == int128_t(al, ah ^ 1ull)) << "a=" << (std::string)a;
        // including when the difference is the sign bit alone
        EXPECT_FALSE(a == int128_t(al, ah ^ (1ull << 63))) << "a=" << (std::string)a;

        // the templated overload sign extends the scalar, so a negative int64_t compares equal only
        // when the high QWORD is all ones, and a positive one only when it is zero
        const int64_t scalar = static_cast<int64_t>(al);
        EXPECT_EQ(a == scalar, ah == SignExtension(scalar)) << "a=" << (std::string)a << ", scalar=" << scalar;
    }
}
// operator!= is declared separately rather than being the compiler's rewrite of operator==, so it
// can disagree with it. Every case below checks both that it is correct and that the two agree.
TEST(int128_t, OperatorNotEqual)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const uint64_t al = get_uint64_random(), ah = get_uint64_random();
        const uint64_t bl = get_uint64_random(), bh = get_uint64_random();
        const int128_t a(al, ah), b(bl, bh);

        EXPECT_EQ(a != b, (al != bl) || (ah != bh)) << "a=" << (std::string)a << ", b=" << (std::string)b;
        EXPECT_EQ(a != b, !(a == b)) << "a=" << (std::string)a << ", b=" << (std::string)b;

        const int128_t copy = a;
        EXPECT_FALSE(a != copy) << "a=" << (std::string)a;
        EXPECT_TRUE(a != int128_t(al ^ 1ull, ah)) << "a=" << (std::string)a;
        EXPECT_TRUE(a != int128_t(al, ah ^ 1ull)) << "a=" << (std::string)a;

        const int64_t scalar = static_cast<int64_t>(al);
        EXPECT_EQ(a != scalar, !(a == scalar)) << "a=" << (std::string)a << ", scalar=" << scalar;
    }
}
TEST(int128_t, OperatorEqualEdgeCases)
{
    const int128_t zero, one(1ll), minus_one(-1ll);
    const int128_t max(~0ull, 0x7FFFFFFFFFFFFFFFull);
    const int128_t min(0ull, 0x8000000000000000ull);

    EXPECT_TRUE(zero == int128_t());
    EXPECT_FALSE(zero == one);
    EXPECT_TRUE(zero != one);
    // max and min differ in the sign bit alone
    EXPECT_FALSE(max == min);
    EXPECT_TRUE(max != min);
    // -1 is all ones in both QWORDs, 2^64-1 fills only the low one
    EXPECT_FALSE(minus_one == int128_t(~0ull, 0ull));
    EXPECT_TRUE(minus_one == int128_t(~0ull, ~0ull));

    // a scalar is sign extended before the comparison
    EXPECT_TRUE(minus_one == -1ll);
    EXPECT_FALSE(int128_t(~0ull, 0ull) == -1ll);  // 2^64-1, positive, not -1
    EXPECT_TRUE(int128_t(5ull, 0ull) == 5ll);
    EXPECT_FALSE(int128_t(5ull, 1ull) == 5ll);
}
// The reference used to be floor(log10(value1)) in double, which is not exact: for k = 23, 24, 28, 29,
// 31, 33, 34, 35, 37 and 38 the double nearest 10^k lies below 10^k, yet log10() of it returns k. No
// random draw happened to land there. The reference is now the exact power of ten table.
TEST(int128_t, log10)
{
    const auto powers = PowersOfTen();
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const double value1 = floor(1.0 + fabs(get_double_random(1, 127)));
        const int128_t i1 = value1;

        uint64_t expected = 0;
        while (expected + 1 < powers.size() && NonNegativeLessOrEqual(powers[expected + 1], i1)) {
            ++expected;
        }
        EXPECT_EQ(log10(i1), expected) << "value=" << (std::string)i1;
    }
}
TEST(int128_t, log)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = 1.0 + floor(fabs(get_double_random(0, 127)));
        int64_t res = (int64_t)floor(log(value1));  // double doesn't lose any bits with this operation!
        if (check_overflow_int128(value1)) {
            continue;
        }
        int128_t i1 = value1;
        int64_t i_res = log(i1);
        EXPECT_EQ(i_res, res) << "double value1=" << value1;
    }
}
TEST(int128_t, sqrt)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = floor(1.0 + fabs(get_double_random(1, 127)));
        if (check_overflow_int128(value1)) {
            continue;
        }
        int128_t i1 = value1;
        uint64_t i_res = sqrt(i1);
        // floor(sqrt(value1)) is not a usable reference at this range, the root needs more than the
        // 53 bits a double holds. is_exact_isqrt() checks r^2 <= x < (r+1)^2 in 128 bit instead.
        EXPECT_TRUE(is_exact_isqrt(i_res, i1)) << "double value1=" << value1 << ", sqrt=" << i_res;
        // the conversion of the operand has to be exact, otherwise the check above passes vacuously
        EXPECT_EQ(static_cast<double>(i1), value1);
    }
}
// The old version declared "int128_t int128_res = static_cast<double>(i2);" and handed that to
// EXPECT_DOUBLE_EQ, which converted it straight back to a double: a round trip through two
// conversions that happened to work only because the operands were tiny. The base never exceeded 7
// and the exponent 15, so the 128 bit range was never approached and negative bases never appeared.
//
// Repeated multiplication is an exact reference at any size, and it stays valid past the point
// where the result wraps because pow() and operator*= truncate identically.
TEST(int128_t, pow)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const int64_t base = static_cast<int64_t>(get_uint32_random() % 200) - 100;  // [-100, 99]
        const uint32_t expo = get_uint32_random() % 40;
        const int128_t x(base);

        int128_t reference(1ll);
        for (auto k = 0u; k < expo; ++k) {
            reference *= x;
        }

        EXPECT_TRUE(pow(x, expo) == reference) << "base=" << base << ", expo=" << expo;
        // an even power of a negative base is positive and matches the power of its magnitude
        if ((expo & 1) == 0) {
            EXPECT_TRUE(pow(x, expo) == pow(-x, expo)) << "base=" << base << ", expo=" << expo;
        } else {
            EXPECT_TRUE(pow(x, expo) == -pow(-x, expo)) << "base=" << base << ", expo=" << expo;
        }
    }
}
TEST(int128_t, powEdgeCases)
{
    // any base to the zero power is one, including zero, matching pow(double, double)
    EXPECT_TRUE(pow(int128_t(), 0u) == int128_t(1ll));
    EXPECT_TRUE(pow(int128_t(5ll), 0u) == int128_t(1ll));
    EXPECT_TRUE(pow(int128_t(-5ll), 0u) == int128_t(1ll));
    // zero to any positive power is zero
    EXPECT_TRUE(pow(int128_t(), 7u).is_zero());
    // one and minus one alternate
    EXPECT_TRUE(pow(int128_t(-1ll), 3u) == int128_t(-1ll));
    EXPECT_TRUE(pow(int128_t(-1ll), 4u) == int128_t(1ll));
    // 2^126 is the largest power of two that stays positive
    EXPECT_TRUE(pow(int128_t(2ll), 126u) == int128_t(0ull, 0x4000000000000000ull));
    // 2^127 wraps to the most negative value, and 2^128 to zero
    EXPECT_TRUE(pow(int128_t(2ll), 127u) == int128_t(0ull, 0x8000000000000000ull));
    EXPECT_TRUE(pow(int128_t(2ll), 128u).is_zero());
    // 3^80 overflows, and must truncate the same way repeated multiplication does
    int128_t reference(1ll);
    for (auto i = 0u; i < 80u; ++i) {
        reference *= int128_t(3ll);
    }
    EXPECT_TRUE(pow(int128_t(3ll), 80u) == reference);
}

/**********************************************************************
 * int128_t regression tests
 *
 * Each test below pins down a defect that was fixed. They deliberately use
 * operand ranges and sign combinations the original tests never reached.
 ***********************************************************************/

// Division reduced both operands to magnitudes, divided, negated, and then subtracted one whenever
// the remainder was non zero, which floors. The builtin integer types truncate towards zero, and
// operator%= was already truncating, so the two disagreed and the division identity did not hold.
TEST(int128_t, DivideTruncatesTowardsZero)
{
    struct {
        int64_t a, b, q, r;
    } cases[] = {{-7, 2, -3, -1}, {7, -2, -3, 1}, {-7, -2, 3, -1}, {7, 2, 3, 1},
                 {-100, 7, -14, -2}, {-1, 2, 0, -1}, {-9, 3, -3, 0}, {-10, 4, -2, -2}};
    // The results are compared at full width: static_cast<int64_t> reads the low QWORD only, and
    // (a / b) * b + (a % b) == a holds for any quotient, the remainder being derived from it.
    for (const auto& c : cases) {
        const int128_t a(c.a), b(c.b);
        EXPECT_TRUE(a / b == FromInt64(c.q)) << c.a << " / " << c.b;
        EXPECT_TRUE(a % b == FromInt64(c.r)) << c.a << " % " << c.b;
        EXPECT_TRUE(IsExactSignedDivision(a, b, a / b, a % b)) << c.a << ", " << c.b;
    }

    // every sign combination, cross checked against int64_t
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const int64_t a = get_int64_random() >> 8;
        const int64_t b = get_int32_random();
        if (b == 0)
            continue;
        const int128_t f1(a), f2(b);
        EXPECT_TRUE(f1 / f2 == FromInt64(a / b)) << "a=" << a << ", b=" << b;
        EXPECT_TRUE(f1 % f2 == FromInt64(a % b)) << "a=" << a << ", b=" << b;
    }
}
// The 64 bit fast path shortcut small numerators with a signed comparison, so any negative value
// compared less than a positive divisor and the result was zeroed before the sign was applied.
TEST(int128_t, DivideByScalarWithNegativeNumerator)
{
    const int128_t a(static_cast<int64_t>(-100));
    EXPECT_TRUE(a / 7 == FromInt64(-14));
    EXPECT_TRUE(a / 7ull == FromInt64(-14));
    EXPECT_TRUE(a / 7u == FromInt64(-14));
    EXPECT_TRUE(a % 7 == FromInt64(-2));
    // a power of two divisor takes a different branch and must truncate too, not shift
    EXPECT_TRUE(int128_t(static_cast<int64_t>(-7)) / 2 == FromInt64(-3));
    EXPECT_TRUE(int128_t(static_cast<int64_t>(-8)) / 2 == FromInt64(-4));
}
// A divisor whose magnitude needs more than 64 bit exercises the 128 bit division path.
TEST(int128_t, DivideByWideInt128)
{
    // 2^66 / (3 * 2^64) truncates to 1, and the remainder lives entirely in the upper QWORD,
    // which the old rounding correction could not see
    const int128_t big(0ull, 4ull);          // 2^66
    const int128_t divisor(0ull, 3ull);      // 3 * 2^64
    EXPECT_TRUE(big / divisor == int128_t(1ll));
    EXPECT_TRUE((-big) / divisor == int128_t(-1ll));
    EXPECT_TRUE(big % divisor == int128_t(0ull, 1ull));      // 2^64
    EXPECT_TRUE((-big) % divisor == -int128_t(0ull, 1ull));  // -2^64

    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const uint64_t aLow = get_uint64_random();
        const uint64_t aHigh = get_uint64_random() >> 1;
        const uint64_t bLow = get_uint64_random();
        const uint64_t bHigh = get_uint64_random();
        const uint32_t bShift = get_uint32_random() % 63;
        int128_t a(aLow, aHigh);
        int128_t b(bLow, bHigh >> bShift);
        if (b.is_zero())
            continue;
        if (get_uint32_random() & 1)
            a = -a;
        if (get_uint32_random() & 1)
            b = -b;
        const int128_t q = a / b, r = a % b;
        EXPECT_TRUE(IsExactSignedDivision(a, b, q, r)) << "a=" << (std::string)a << ", b=" << (std::string)b;
        int128_t expectedQ, expectedR;
        ReferenceDivide(a, b, expectedQ, expectedR);
        EXPECT_TRUE(q == expectedQ && r == expectedR) << "a=" << (std::string)a << ", b=" << (std::string)b;
    }
}
// sqrt only rejected negatives, then called log2, which throws on zero. sqrt is noexcept, so a
// zero input terminated the process.
TEST(int128_t, SqrtOfZeroDoesNotThrow)
{
    EXPECT_EQ(sqrt(int128_t()), 0ull);
    EXPECT_EQ(sqrt(int128_t(0ull)), 0ull);
    EXPECT_EQ(sqrt(int128_t(1ull)), 1ull);
    EXPECT_EQ(sqrt(int128_t(static_cast<int64_t>(-1))), 0ull);
    // the log functions do throw, that part is intentional
    EXPECT_THROW((void)log2(int128_t()), std::domain_error);
    EXPECT_THROW((void)log(int128_t()), std::domain_error);
    EXPECT_THROW((void)log10(int128_t()), std::domain_error);
    EXPECT_THROW((void)log2(int128_t(static_cast<int64_t>(-1))), std::domain_error);
}
// operator^ was implemented with &=, so it computed a bitwise AND.
TEST(int128_t, BitwiseOperators)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const uint64_t al = get_uint64_random(), ah = get_uint64_random();
        const uint64_t bl = get_uint64_random(), bh = get_uint64_random();
        const int128_t a(al, ah), b(bl, bh);
        EXPECT_TRUE((a ^ b) == int128_t(al ^ bl, ah ^ bh));
        EXPECT_TRUE((a & b) == int128_t(al & bl, ah & bh));
        EXPECT_TRUE((a | b) == int128_t(al | bl, ah | bh));
        EXPECT_TRUE((a ^ a).is_zero());
    }
}
// The terminator sat at index 32 of a 45 byte buffer with the digits written backwards from there,
// so a long value ran off the front. A negative value needs one more character than an unsigned
// one of the same width. Separately, the fast path tested high == 0 and formatted with %lld, which
// printed every value in [2^63, 2^64) as a negative number.
TEST(int128_t, ToStringLargeValues)
{
    const char* values[] = {"9223372036854775808",                       // 2^63, the %lld trap
                            "18446744073709551615",                      // 2^64-1, still positive
                            "18446744073709551616",                      // 2^64
                            "-9223372036854775809",                      // just past int64_t's range
                            "170141183460469231731687303715884105727",   // INT128_MAX, 39 digits
                            "-170141183460469231731687303715884105728"}; // INT128_MIN, 39 digits and a sign

    // a value that fits in an int64_t takes the snprintf path, which returns the buffer start
    const char* buffer_start = static_cast<char*>(int128_t(1ull));
    for (auto i = 0u; i < array_length(values); ++i) {
        const int128_t v = values[i];
        const char* p = static_cast<char*>(v);
        EXPECT_STREQ(p, values[i]);
        EXPECT_GE(p, buffer_start) << "wrote before the start of the buffer: " << values[i];
    }

    // full round trip over the whole 128 bit range, both signs
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const uint64_t low = get_uint64_random();
        const uint64_t high = get_uint64_random();
        const int128_t v(low, high);
        const std::string s = (std::string)v;
        EXPECT_TRUE(int128_t(s) == v) << "s=" << s;
    }
}
// The conversions used to truncate the bits that did not fit the mantissa. strtod and strtof are
// correctly rounded by the C standard, so the decimal string of a value is an exact reference.
TEST(int128_t, ConversionToFloatingPointIsCorrectlyRounded)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const uint64_t low = get_uint64_random();
        const uint64_t high = get_uint64_random();
        const uint32_t shift = get_uint32_random() % 64;
        int128_t v(low, high >> shift);
        if (v.is_zero())
            continue;
        if (get_uint32_random() & 1)
            v = -v;
        const std::string s = (std::string)v;
        EXPECT_EQ(static_cast<double>(v), strtod(s.c_str(), nullptr)) << "value=" << s;
        EXPECT_EQ(static_cast<float>(v), strtof(s.c_str(), nullptr)) << "value=" << s;
    }
}
TEST(int128_t, ConversionToFloatingPointEdgeCases)
{
    // the most negative value: its absolute value keeps the sign bit set
    const int128_t min_val(0ull, 0x8000000000000000ull);
    EXPECT_EQ(static_cast<double>(min_val), -170141183460469231731687303715884105728.0);
    EXPECT_EQ(static_cast<float>(min_val), -170141183460469231731687303715884105728.0f);
    // INT128_MAX rounds up to 2^127, which both types can hold
    const int128_t max_val(~0ull, 0x7FFFFFFFFFFFFFFFull);
    EXPECT_EQ(static_cast<double>(max_val), 170141183460469231731687303715884105728.0);
    EXPECT_FALSE(std::isinf(static_cast<float>(max_val)));
    // exact values need no rounding
    EXPECT_EQ(static_cast<double>(int128_t(0ull, 1ull)), 18446744073709551616.0);
    EXPECT_EQ(static_cast<double>(int128_t(static_cast<int64_t>(-1))), -1.0);
}
// Shifting by 128 or more used to be taken modulo 64 by the hardware and leave bits behind.
// An arithmetic right shift has to leave copies of the sign bit, not zeros.
TEST(int128_t, ShiftBeyondWidth)
{
    volatile int32_t shift = 128;
    const int128_t positive(0ull, 0x4000000000000000ull);
    const int128_t negative(static_cast<int64_t>(-1));
    EXPECT_TRUE((positive >> shift).is_zero());
    EXPECT_EQ(static_cast<int64_t>(negative >> shift), -1);
    EXPECT_TRUE((positive << shift).is_zero());
    shift = 200;
    EXPECT_TRUE((positive >> shift).is_zero());
    EXPECT_EQ(static_cast<int64_t>(negative >> shift), -1);
    EXPECT_TRUE((positive << shift).is_zero());
    // the boundary either side still behaves
    shift = 127;
    EXPECT_EQ(static_cast<int64_t>(negative >> shift), -1);
    EXPECT_TRUE((int128_t(1ull) << shift) == int128_t(0ull, 0x8000000000000000ull));
}
// A comma used to continue the parse loop without advancing the read pointer, spinning forever.
// A high byte handed straight to isspace is outside the domain it accepts.
TEST(int128_t, ConstructorFromStringEdgeCases)
{
    EXPECT_TRUE(int128_t("1,000") == int128_t(1000ull));
    EXPECT_TRUE(int128_t("-1,000") == int128_t(static_cast<int64_t>(-1000)));
    EXPECT_TRUE(int128_t("+5") == int128_t(5ull));
    EXPECT_TRUE(int128_t("-0x10") == int128_t(static_cast<int64_t>(-16)));
    EXPECT_TRUE(int128_t("-0").is_zero());
    const char high_byte[] = {'1', '2', static_cast<char>(0xB5), '3', 0};
    EXPECT_TRUE(int128_t(high_byte) == int128_t(12ull));
}
// One overload per fixed width type left long and unsigned long ambiguous on MSVC, and a scalar on
// the left hand side of a binary operator could not resolve at all.
TEST(int128_t, IntegralTypesAndScalarOnTheLeft)
{
    EXPECT_TRUE(int128_t(static_cast<unsigned long>(7)) == int128_t(7ull));
    EXPECT_TRUE(int128_t(static_cast<long>(-7)) == int128_t(static_cast<int64_t>(-7)));
    EXPECT_TRUE(int128_t(static_cast<short>(-3)) == int128_t(static_cast<int64_t>(-3)));
    EXPECT_TRUE(int128_t(true) == int128_t(1ull));
    EXPECT_TRUE(int128_t(static_cast<size_t>(9)) == int128_t(9ull));

    const int128_t x(10ull);
    EXPECT_TRUE((1 + x) == int128_t(11ull));
    EXPECT_TRUE((100 - x) == int128_t(90ull));
    EXPECT_TRUE((3 * x) == int128_t(30ull));
    EXPECT_TRUE((100 / x) == int128_t(10ull));
    EXPECT_TRUE((105 % x) == int128_t(5ull));
    EXPECT_TRUE((3 ^ x) == int128_t(9ull));
    EXPECT_TRUE((-100 / x) == int128_t(static_cast<int64_t>(-10)));
    // the int128_t on the left forms must still resolve
    EXPECT_TRUE((x + 1) == int128_t(11ull));
    EXPECT_TRUE((x + x) == int128_t(20ull));
}
// int128_t and uint128_t are now aliases of a single class template. Nothing but this test pins the
// properties the merge had to preserve: the storage is unchanged, the 16 byte alignment survived the
// template, and the two types still refuse to convert into one another.
//
// The type is deliberately *not* trivially copyable; see the uint128_t version of this test for the
// store forwarding stall that a defaulted 16 byte copy costs on MSVC.
TEST(int128_t, LayoutAndTraits)
{
    static_assert(std::is_same_v<int128_t, int128_base<true>>, "int128_t is the signed instantiation");
    static_assert(int128_t::is_signed, "int128_t reports itself as signed");
    static_assert(sizeof(int128_t) == 16, "the object holds exactly two QWORDs");
    static_assert(alignof(int128_t) == 16, "FP128_ALIGN16 has to survive the template");
    static_assert(std::is_standard_layout_v<int128_t>, "the division helpers alias the members as an array");
    static_assert(!std::is_trivially_copyable_v<int128_t>, "the copy and move members are written out, not defaulted");
    // getting from one signedness to the other needs two user defined conversions, which the
    // language never performs implicitly, and leaves the explicit form ambiguous
    static_assert(!std::is_convertible_v<uint128_t, int128_t>, "no silent unsigned to signed conversion");
    static_assert(!std::is_constructible_v<int128_t, uint128_t>, "and no unambiguous explicit one either");
    SUCCEED();
}
// Every math function on these types is a hidden friend, found by argument dependent lookup. A value
// whose low QWORD is zero catches the failure mode where a call falls through to the narrower
// fp128::log2(uint64_t) overload through the implicit conversion to uint64_t: that overload would
// see a zero and answer zero.
TEST(int128_t, MathFunctionsSeeTheWholeValue)
{
    const int128_t x(0ull, 1ull << 40);  // 2^104, the low QWORD is zero

    EXPECT_EQ(log2(x), 104ull);
    EXPECT_EQ(lzcnt128(x), 23ull);
    EXPECT_EQ(log10(x), 31ull);
    EXPECT_EQ(log(x), 72ull);
    EXPECT_EQ(sqrt(x), 1ull << 52);
    EXPECT_TRUE(pow(int128_t(2ll), 104u) == x);
    EXPECT_TRUE(sqr(int128_t(1ull << 52)) == x);
}
// abs() is a hidden friend constrained to the signed instantiation, so it is never found by ordinary
// unqualified lookup and can never hide ::abs() for the builtin types inside namespace fp128.
TEST(int128_t, AbsResolvesToTheHiddenFriend)
{
    // the result type has to be int128_t, not the int that ::abs() would return
    static_assert(std::is_same_v<decltype(abs(int128_t(1ll))), int128_t>, "abs must not fall through to ::abs");

    const int128_t x("-170141183460469231731687303715884105727");  // -(2^127 - 1)
    EXPECT_TRUE(abs(x) == int128_t("170141183460469231731687303715884105727"));
    EXPECT_TRUE(abs(int128_t(-5ll)) == int128_t(5ll));
    EXPECT_TRUE(abs(int128_t(5ll)) == int128_t(5ll));
    EXPECT_TRUE(abs(int128_t()).is_zero());
}
// operator*= multiplies the two's complement bit patterns directly instead of converting both
// operands to magnitudes and reapplying the sign afterwards. A truncated 128 bit product is a
// multiplication modulo 2^128, and in that ring the bit pattern of a negative value is its value,
// so the result is bit identical to what the magnitude based version produced.
TEST(int128_t, MultiplySignCombinations)
{
    const int128_t a("123456789012345678901234567890");
    const int128_t b(-987654321ll);

    EXPECT_TRUE((a * b) == -(a * -b));
    EXPECT_TRUE((-a * b) == (a * -b));
    EXPECT_TRUE((-a * -b) == (a * b));
    EXPECT_TRUE((a * b).is_negative());
    EXPECT_TRUE((-a * b).is_positive());

    // the most negative value is its own negation, which a magnitude based implementation has to
    // treat as a special case
    const int128_t most_negative(0ull, 0x8000000000000000ull);
    EXPECT_TRUE((most_negative * int128_t(1ll)) == most_negative);
    EXPECT_TRUE((most_negative * int128_t(-1ll)) == most_negative);
    EXPECT_TRUE((most_negative * int128_t(2ll)).is_zero());
}
// The right hand side is allowed to alias this object, which the implementation guards against by
// snapshotting both operands before it writes anything.
TEST(int128_t, MultiplyInPlaceAliasing)
{
    int128_t a(-123456789ll);
    const int128_t expected = a * a;
    a *= a;
    EXPECT_TRUE(a == expected);
    EXPECT_TRUE(a == int128_t(15241578750190521ll));
}
// The scalar overload used to convert to int128_t and take the full path. It now uses the same 64
// bit shortcut the unsigned type does, which is valid for a negative object for the same modulo
// 2^128 reason operator*=(const int128_t&) is.
TEST(int128_t, MultiplyByScalarKeepsTheSign)
{
    int128_t a(-1000ll);
    a *= 3u;
    EXPECT_TRUE(a == int128_t(-3000ll));

    int128_t b("-123456789012345678901234567890");
    b *= 7ull;
    EXPECT_TRUE(b == int128_t("-864197523086419752308641975230"));

    // a negative scalar still goes the long way round
    int128_t c(-1000ll);
    c *= -3;
    EXPECT_TRUE(c == int128_t(3000ll));

    // and so does a floating point one
    int128_t d(-1000ll);
    d *= 2.5;
    EXPECT_TRUE(d == int128_t(-2000ll));
}
// operator/= detects a power of two divisor and shifts instead of dividing. The shift is a logical
// one applied to the magnitude, so it has to be correct on both sides of the 64 bit boundary and
// for a divisor of exactly one, which shifts by zero.
TEST(int128_t, DivideByPowersOfTwo)
{
    const int128_t x("170141183460469231731687303715884105727");  // 2^127 - 1

    EXPECT_TRUE((x / int128_t(1ll)) == x);
    EXPECT_TRUE((x / int128_t(2ll)) == (x >> 1));
    EXPECT_TRUE((x / int128_t(1ull << 63)) == (x >> 63));
    EXPECT_TRUE((x / int128_t(0ull, 1ull)) == (x >> 64));         // 2^64
    EXPECT_TRUE((x / int128_t(0ull, 1ull << 36)) == (x >> 100));  // 2^100

    // a negative dividend truncates towards zero, which an arithmetic shift of a negative value
    // does not do, so the magnitude has to be shifted and negated afterwards
    EXPECT_TRUE((-x / int128_t(2ll)) == -(x / int128_t(2ll)));
    EXPECT_TRUE((-x / int128_t(0ull, 1ull)) == -(x >> 64));

    // the magnitude of the most negative value keeps its sign bit set, which is why the shift works
    // on the raw QWORDs instead of going through the signed operators
    const int128_t most_negative(0ull, 0x8000000000000000ull);  // -2^127
    EXPECT_TRUE((most_negative / int128_t(1ll)) == most_negative);
    EXPECT_TRUE((most_negative / int128_t(2ll)) == int128_t(0ull, 0xC000000000000000ull));  // -2^126
    EXPECT_TRUE((most_negative / int128_t(0ull, 1ull)) == -int128_t(1ull << 63));           // -2^63
}
// Public members that no other test in this file reaches. All of them work today, so this pins
// current behaviour rather than a past defect.
TEST(int128_t, UnaryOperatorsAndAccessors)
{
    const int128_t x(0x0123456789ABCDEFull, 0xFEDCBA9876543210ull);
    const int128_t min(0ull, 0x8000000000000000ull);

    // unary plus is a copy, unary minus is a two's complement negation
    EXPECT_TRUE(+x == x);
    EXPECT_TRUE(-(-x) == x);
    EXPECT_TRUE(x + (-x) == int128_t());
    // ~x == -x - 1 is the defining relation of two's complement
    EXPECT_TRUE(~x == -x - int128_t(1ll));
    EXPECT_TRUE(~(~x) == x);
    EXPECT_TRUE(~int128_t() == int128_t(-1ll));
    // the most negative value is its own negation, having no positive counterpart
    EXPECT_TRUE(-min == min);

    // operator bool and operator!, which are each other's opposite
    EXPECT_TRUE(static_cast<bool>(x));
    EXPECT_FALSE(static_cast<bool>(int128_t()));
    EXPECT_FALSE(!x);
    EXPECT_TRUE(!int128_t());
    // a value whose low QWORD is zero is still non zero
    EXPECT_TRUE(static_cast<bool>(int128_t(0ull, 1ull)));

    // the sign predicates, including the boundary between them
    EXPECT_TRUE(int128_t(1ll).is_positive());
    EXPECT_TRUE(int128_t().is_positive());  // zero counts as positive
    EXPECT_FALSE(int128_t().is_negative());
    EXPECT_TRUE(int128_t(-1ll).is_negative());
    EXPECT_TRUE(min.is_negative());
    EXPECT_TRUE(int128_t(~0ull, 0x7FFFFFFFFFFFFFFFull).is_positive());
    // an integer type has no fraction
    EXPECT_TRUE(x.is_int());

    // get_bit over both QWORDs, and the sign bit at the top
    EXPECT_EQ(int128_t(1ll).get_bit(0), 1);
    EXPECT_EQ(int128_t(1ll).get_bit(1), 0);
    EXPECT_EQ(int128_t(0ull, 1ull).get_bit(64), 1);
    EXPECT_EQ(min.get_bit(127), 1);
    EXPECT_EQ(int128_t(-1ll).get_bit(127), 1);
    EXPECT_EQ(int128_t(1ll).get_bit(127), 0);

    // one() and get_components()
    EXPECT_TRUE(int128_t::one() == int128_t(1ll));
    uint64_t low = 0, high = 0;
    x.get_components(low, high);
    EXPECT_EQ(low, 0x0123456789ABCDEFull);
    EXPECT_EQ(high, 0xFEDCBA9876543210ull);

    // hex() prints both QWORDs, high first, zero padded
    EXPECT_STREQ(x.hex(), "0xFEDCBA98765432100123456789ABCDEF");
    EXPECT_STREQ(int128_t().hex(), "0x00000000000000000000000000000000");
    EXPECT_STREQ(int128_t(-1ll).hex(), "0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF");

    // to_string() has to agree with the std::string conversion
    const int128_t big("-170141183460469231731687303715884105728");
    EXPECT_STREQ(big.to_string(), "-170141183460469231731687303715884105728");
    EXPECT_STREQ(big.to_string(), static_cast<std::string>(big).c_str());

    // long double goes through the double conversion
    EXPECT_EQ(static_cast<long double>(int128_t(-1ll)), -1.0L);
    EXPECT_EQ(static_cast<long double>(min), static_cast<long double>(static_cast<double>(min)));
}
// The user defined literal is declared at namespace scope in int128_t.h and was never exercised.
TEST(int128_t, UserDefinedLiteral)
{
    EXPECT_TRUE(1234_int128 == int128_t(1234ll));
    EXPECT_TRUE(-1234_int128 == int128_t(-1234ll));
    EXPECT_TRUE(0xDEADBEAF_int128 == int128_t(0xDEADBEAFll));
    EXPECT_TRUE(170141183460469231731687303715884105727_int128 == int128_t(~0ull, 0x7FFFFFFFFFFFFFFFull));
    // the operator takes the raw characters, so digit grouping survives
    EXPECT_TRUE(1'000'000_int128 == int128_t(1000000ll));
}
// The compound bitwise and shift assignments. They used to be compared with the binary operators,
// which are implemented by calling them, so the random half checked a function against itself: an
// arithmetic right shift that went logical for counts 11 to 62 passed the whole file. Both forms are
// now compared with references built from the QWORDs. The scalar operands are covered by
// BitwiseWithScalarOperandsExtendTheScalar.
TEST(int128_t, CompoundBitwiseAndShiftOperators)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const uint64_t al = get_uint64_random(), ah = get_uint64_random();
        const uint64_t bl = get_uint64_random(), bh = get_uint64_random();
        const int128_t a(al, ah), b(bl, bh);

        int128_t c = a;
        c &= b;
        EXPECT_TRUE(c == int128_t(al & bl, ah & bh)) << "a=" << (std::string)a << ", b=" << (std::string)b;
        c = a;
        c |= b;
        EXPECT_TRUE(c == int128_t(al | bl, ah | bh)) << "a=" << (std::string)a << ", b=" << (std::string)b;
        c = a;
        c ^= b;
        EXPECT_TRUE(c == int128_t(al ^ bl, ah ^ bh)) << "a=" << (std::string)a << ", b=" << (std::string)b;

        // both shift forms against the bit by bit reference
        const int32_t shift = static_cast<int32_t>(get_uint32_random() % 128);
        const int128_t expectedLeft = ReferenceShift(a, shift, true);
        const int128_t expectedRight = ReferenceShift(a, shift, false);
        c = a;
        c <<= shift;
        EXPECT_TRUE(c == expectedLeft) << "a=" << (std::string)a << ", shift=" << shift;
        EXPECT_TRUE((a << shift) == expectedLeft) << "a=" << (std::string)a << ", shift=" << shift;
        c = a;
        c >>= shift;
        EXPECT_TRUE(c == expectedRight) << "a=" << (std::string)a << ", shift=" << shift;
        EXPECT_TRUE((a >> shift) == expectedRight) << "a=" << (std::string)a << ", shift=" << shift;
    }

    // a zero or negative shift count leaves the value alone, as documented
    volatile int32_t shift = 0;
    const int128_t x(0x0123456789ABCDEFull, 0xFEDCBA9876543210ull);
    EXPECT_TRUE((x >> shift) == x);
    EXPECT_TRUE((x << shift) == x);
    shift = -1;
    EXPECT_TRUE((x >> shift) == x);
    EXPECT_TRUE((x << shift) == x);

    // the right shift is arithmetic: it replicates the sign bit rather than shifting in zeros
    EXPECT_TRUE((int128_t(-1ll) >> 1) == int128_t(-1ll));
    EXPECT_TRUE((int128_t(-2ll) >> 1) == int128_t(-1ll));
    EXPECT_TRUE((int128_t(0ull, 0x8000000000000000ull) >> 127) == int128_t(-1ll));
    EXPECT_TRUE((int128_t(0ull, 0x4000000000000000ull) >> 126) == int128_t(1ll));
}

// A scalar on the left used to be ambiguous: converting it to int128_t and converting the int128_t
// to a builtin type are both a single user defined conversion.
TEST(int128_t, ScalarOnTheLeftHandSide)
{
    const int128_t x(10);
    EXPECT_TRUE((1 + x) == int128_t(11));
    EXPECT_TRUE((100 - x) == int128_t(90));
    EXPECT_TRUE((3 * x) == int128_t(30));
    EXPECT_TRUE((100 / x) == int128_t(10));
    EXPECT_TRUE((105 % x) == int128_t(5));
    EXPECT_TRUE((0xFF & x) == int128_t(10));
    EXPECT_TRUE((1 | x) == int128_t(11));
    EXPECT_TRUE((3 ^ x) == int128_t(9));

    // subtracting past zero must keep the sign rather than wrap
    EXPECT_TRUE((1 - x) == int128_t(-9));
    EXPECT_TRUE((-1 * x) == int128_t(-10));

    // The left operand is widened rather than the object being narrowed, so the result is 128 bit.
    static_assert(std::is_same_v<decltype(1 + x), int128_t>);
    static_assert(std::is_same_v<decltype(1 ^ x), int128_t>);
}
// The comparisons only existed with the int128_t on the left. A scalar there had nothing but the
// builtin comparisons to choose from, and the conversion operators of int128_t make every one of
// them equally good, so the call was ambiguous. Signedness is the point here: resolving to any of
// the unsigned builtin conversions would order the negative values above the positive ones.
TEST(int128_t, ScalarOnTheLeftHandSideComparisons)
{
    const int128_t x(10);
    EXPECT_TRUE(10 == x);
    EXPECT_TRUE(9 != x);
    EXPECT_TRUE(9 < x);
    EXPECT_FALSE(10 < x);
    EXPECT_TRUE(9 <= x);
    EXPECT_TRUE(10 <= x);
    EXPECT_TRUE(11 > x);
    EXPECT_FALSE(10 > x);
    EXPECT_TRUE(11 >= x);
    EXPECT_TRUE(10 >= x);

    // across zero, where a narrowed unsigned comparison would give the opposite answer
    const int128_t negative(-3);
    EXPECT_TRUE(-3 == negative);
    EXPECT_TRUE(-4 < negative);
    EXPECT_TRUE(-2 > negative);
    EXPECT_TRUE(0 > negative);
    EXPECT_TRUE(1 > negative);
    EXPECT_FALSE(0 < negative);
    EXPECT_TRUE(0 > int128_t(-1));
    EXPECT_TRUE(-1 < int128_t(0));

    // a value no builtin type can hold must not compare through a narrowed conversion
    const int128_t big = int128_t(1) << 100;
    EXPECT_TRUE(1 < big);
    EXPECT_FALSE(1 > big);
    EXPECT_TRUE(1 != big);
    EXPECT_TRUE(1 > -big);

    // the int128_t on the left forms must still resolve
    EXPECT_TRUE(x == 10);
    EXPECT_TRUE(x > 9);
    EXPECT_TRUE(x == int128_t(10));
}
// Widening the left operand is what makes (1 << 100) produce the expected power of two. The same
// expression on a builtin int is undefined behavior, the shift count being wider than the operand.
TEST(int128_t, ScalarOnTheLeftHandSideShifts)
{
    const int128_t ten(10);
    EXPECT_TRUE((1 << ten) == int128_t(1024));
    EXPECT_TRUE((1024 >> ten) == int128_t(1));
    EXPECT_TRUE((3 << int128_t(2)) == int128_t(12));

    const int128_t hundred(100);
    EXPECT_TRUE((1 << hundred) == (int128_t(1) << 100));
    EXPECT_TRUE((1 << hundred) != 0);
    static_assert(std::is_same_v<decltype(1 << hundred), int128_t>);

    // the right shift is arithmetic, as it is for the int128_t on the left form
    EXPECT_TRUE((-1024 >> ten) == int128_t(-1));
    EXPECT_TRUE((-8 >> int128_t(1)) == int128_t(-4));

    // the int128_t on the left forms must still resolve
    EXPECT_TRUE((ten << 1) == int128_t(20));
    EXPECT_TRUE((ten >> 1) == int128_t(5));
}

/**********************************************************************
 * int128_t edge cases and contracts
 *
 * Paths the coverage report showed as never executed for the signed type, and behavior a mutation
 * of the library survived. Expected values come from Python's arbitrary precision integers
 * (truncating division written out, since Python's // floors) or are built from QWORDs here.
 ***********************************************************************/

// The zero, saturation and below-one branches of the double constructor never ran: the random operands
// were floor()ed doubles in [2, 2^127). NaN, the infinities and every magnitude of 2^127 or more saturate
// to the end of the range with the sign of the input; a fraction truncates towards zero.
TEST(int128_t, ConstructorFromDoubleSaturatesAndTruncates)
{
    const int128_t max(~0ull, 0x7FFFFFFFFFFFFFFFull);
    const int128_t min(0ull, 0x8000000000000000ull);
    const double inf = std::numeric_limits<double>::infinity();
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double denorm = std::numeric_limits<double>::denorm_min();

    // the expected value as QWORDs: an int128_t member would pad the struct to its 16 byte alignment
    const struct {
        double value;
        uint64_t low, high;
    } cases[] = {
        {0.0, 0ull, 0ull},
        {-0.0, 0ull, 0ull},
        {0.5, 0ull, 0ull},
        {-0.5, 0ull, 0ull},
        {-0.9999999, 0ull, 0ull},
        {denorm, 0ull, 0ull},
        {-denorm, 0ull, 0ull},
        {-DBL_MIN, 0ull, 0ull},
        {3.99, 3ull, 0ull},
        {-3.99, ~0ull - 2ull, ~0ull},
        {-1.0, ~0ull, ~0ull},
        {-0x1p53, 0xFFE0000000000000ull, ~0ull},
        {0x1p126, 0ull, 0x4000000000000000ull},
        {-0x1p126, 0ull, 0xC000000000000000ull},
        {0x1.fffffffffffffp126, 0ull, 0x7FFFFFFFFFFFFC00ull},  // the largest double below 2^127
        {-0x1.fffffffffffffp126, 0ull, 0x8000000000000400ull},
        {0x1p127, ~0ull, 0x7FFFFFFFFFFFFFFFull},
        {-0x1p127, 0ull, 0x8000000000000000ull},  // exactly the most negative value, no saturation needed
        {-0x1.0000000000001p127, 0ull, 0x8000000000000000ull},
        {1e39, ~0ull, 0x7FFFFFFFFFFFFFFFull},
        {-1e39, 0ull, 0x8000000000000000ull},
        {DBL_MAX, ~0ull, 0x7FFFFFFFFFFFFFFFull},
        {-DBL_MAX, 0ull, 0x8000000000000000ull},
        {inf, ~0ull, 0x7FFFFFFFFFFFFFFFull},
        {-inf, 0ull, 0x8000000000000000ull},
        {nan, ~0ull, 0x7FFFFFFFFFFFFFFFull},
        {std::copysign(nan, -1.0), 0ull, 0x8000000000000000ull},
    };
    for (const auto& c : cases) {
        const double value = opaque(c.value);  // the runtime path, not the constant evaluated one
        const int128_t expected(c.low, c.high);
        EXPECT_TRUE(int128_t(value) == expected) << "value=" << std::hexfloat << value;
        // and the reference the random test relies on agrees with every row
        EXPECT_TRUE(ExpectedFromDouble(value) == expected) << "value=" << std::hexfloat << value;
    }

    EXPECT_TRUE(int128_t(opaque(-1.5f)) == int128_t(-1ll));
    EXPECT_TRUE(int128_t(opaque(3.4e38f)) == max);  // a float above 2^127
    EXPECT_TRUE(int128_t(opaque(-0x1p127f)) == min);
    EXPECT_TRUE(int128_t(opaque(-std::numeric_limits<float>::infinity())) == min);
}
// Round to nearest, ties to even, at the places random operands never reach: an exact tie, and the
// sticky bit living in the high QWORD alone while the low QWORD is zero. That last branch was never
// evaluated for the signed type, and ignoring it went unnoticed by both integer suites. Each row is
// checked for the value and its negation, whose result differs in the sign bit only. Bit patterns
// from Python's exact integer rounding and its struct module.
TEST(int128_t, ConversionRoundsTiesToEven)
{
    const struct {
        uint64_t low, high, doubleBits;
        uint32_t floatBits;
    } rows[] = {
        {0x0020000000000001ull, 0x0000000000000000ull, 0x4340000000000000ull, 0x5A000000u},  // 2^53+1: one dropped bit, a tie, rounds down to even
        {0x0020000000000003ull, 0x0000000000000000ull, 0x4340000000000002ull, 0x5A000000u},  // 2^53+3: one dropped bit, a tie, rounds up to even
        {0x0040000000000002ull, 0x0000000000000000ull, 0x4350000000000000ull, 0x5A800000u},  // 2^54+2: two dropped bits, a tie, rounds down
        {0x0040000000000006ull, 0x0000000000000000ull, 0x4350000000000002ull, 0x5A800000u},  // 2^54+6: two dropped bits, a tie, rounds up
        {0x0000000001000001ull, 0x0000000000000000ull, 0x4170000010000000ull, 0x4B800000u},  // 2^24+1: float tie on one dropped bit
        {0x0000000001000003ull, 0x0000000000000000ull, 0x4170000030000000ull, 0x4B800002u},  // 2^24+3: float tie, rounds up
        {0x0000000000000000ull, 0x0100000000000008ull, 0x4770000000000000ull, 0x7B800000u},  // 2^120+2^67: low QWORD zero, a double tie
        {0x0000000000000000ull, 0x0100000000000009ull, 0x4770000000000001ull, 0x7B800000u},  // 2^120+2^67+2^64: sticky bit in the high QWORD only
        {0x0000000000000000ull, 0x0100000100000000ull, 0x4770000010000000ull, 0x7B800000u},  // 2^120+2^96: low QWORD zero, a float tie
        {0x0000000000000000ull, 0x0100000100000001ull, 0x4770000010000000ull, 0x7B800001u},  // 2^120+2^96+2^64: float sticky bit in the high QWORD only
        {0xFFFFFFFFFFFFFFFFull, 0x7FFFFFFFFFFFFFFFull, 0x47E0000000000000ull, 0x7F000000u},  // 2^127-1: rounds up to 2^127
    };
    for (const auto& row : rows) {
        const int128_t v(row.low, row.high);
        EXPECT_EQ(std::bit_cast<uint64_t>(static_cast<double>(v)), row.doubleBits) << "v=" << (std::string)v;
        EXPECT_EQ(std::bit_cast<uint32_t>(static_cast<float>(v)), row.floatBits) << "v=" << (std::string)v;
        const int128_t negated = -v;
        EXPECT_EQ(std::bit_cast<uint64_t>(static_cast<double>(negated)), row.doubleBits | 0x8000000000000000ull) << "v=" << (std::string)negated;
        EXPECT_EQ(std::bit_cast<uint32_t>(static_cast<float>(negated)), row.floatBits | 0x80000000u) << "v=" << (std::string)negated;
    }

    // zero converts to +0.0, never -0.0
    const int128_t zero(opaque(0ll));
    EXPECT_EQ(std::bit_cast<uint64_t>(static_cast<double>(zero)), 0ull);
    EXPECT_EQ(std::bit_cast<uint32_t>(static_cast<float>(zero)), 0u);
}
// The magnitude of -2^127 is the only one with bit 127 set, and the edge tests only ever divided it
// through the power of two and equality shortcuts. Dropping that bit in the long division paths passed
// the whole file. Quotients and remainders from Python, truncated towards zero.
TEST(int128_t, DivideMostNegativeThroughLongDivision)
{
    const int128_t max(~0ull, 0x7FFFFFFFFFFFFFFFull);
    const int128_t min(0ull, 0x8000000000000000ull);
    const int128_t minusOne(-1ll);

    const struct {
        int128_t divisor, quotient, remainder;
    } cases[] = {
        {max, minusOne, minusOne},
        {-max, int128_t(1ll), minusOne},
        {int128_t(3ll), int128_t(0x5555555555555556ull, 0xD555555555555555ull), int128_t(-2ll)},
        {int128_t(-3ll), int128_t(0xAAAAAAAAAAAAAAAAull, 0x2AAAAAAAAAAAAAAAull), int128_t(-2ll)},
        {int128_t(7ll), int128_t(0xDB6DB6DB6DB6DB6Eull, 0xEDB6DB6DB6DB6DB6ull), int128_t(-2ll)},
        {int128_t(10ll), int128_t(0x3333333333333334ull, 0xF333333333333333ull), int128_t(-8ll)},
        {FromInt64(INT64_MIN), int128_t(0ull, 1ull), int128_t()},
        {int128_t(UINT64_MAX, 0ull), FromInt64(INT64_MIN), FromInt64(INT64_MIN)},
        {int128_t(1ull, 1ull), int128_t(0x8000000000000001ull, ~0ull), int128_t(0x7FFFFFFFFFFFFFFFull, ~0ull)},     // 2^64+1
        {-int128_t(1ull, 1ull), int128_t(0x7FFFFFFFFFFFFFFFull, 0ull), int128_t(0x7FFFFFFFFFFFFFFFull, ~0ull)},
        {int128_t(1ull, 0x4000000000000000ull), minusOne, int128_t(1ull, 0xC000000000000000ull)},                // 2^126+1
        {int128_t(12345ull, 1ull << 36), int128_t(0xFFFFFFFFF8000001ull, ~0ull), int128_t(0x00000181C7FFCFC7ull, 0xFFFFFFF000000000ull)},  // 2^100+12345
    };
    for (const auto& c : cases) {
        const int128_t q = min / c.divisor;
        const int128_t r = min % c.divisor;
        EXPECT_TRUE(q == c.quotient) << "divisor=" << (std::string)c.divisor << ", q=" << (std::string)q;
        EXPECT_TRUE(r == c.remainder) << "divisor=" << (std::string)c.divisor << ", r=" << (std::string)r;
        EXPECT_TRUE(IsExactSignedDivision(min, c.divisor, q, r)) << "divisor=" << (std::string)c.divisor;
        int128_t expectedQ, expectedR;
        ReferenceDivide(min, c.divisor, expectedQ, expectedR);
        EXPECT_TRUE(expectedQ == c.quotient && expectedR == c.remainder) << "the reference disagrees, divisor=" << (std::string)c.divisor;

        int128_t compound = min;
        compound /= c.divisor;
        EXPECT_TRUE(compound == c.quotient) << "divisor=" << (std::string)c.divisor;
        compound = min;
        compound %= c.divisor;
        EXPECT_TRUE(compound == c.remainder) << "divisor=" << (std::string)c.divisor;
    }

    // the same through the scalar overloads
    EXPECT_TRUE(min / 3 == int128_t(0x5555555555555556ull, 0xD555555555555555ull));
    EXPECT_TRUE(min % 3 == int128_t(-2ll));
    EXPECT_TRUE(min / -3 == int128_t(0xAAAAAAAAAAAAAAAAull, 0x2AAAAAAAAAAAAAAAull));
    EXPECT_TRUE(min / 10u == int128_t(0x3333333333333334ull, 0xF333333333333333ull));
    EXPECT_TRUE(min % 10u == int128_t(-8ll));
    EXPECT_TRUE(min / INT64_MIN == int128_t(0ull, 1ull));
    EXPECT_TRUE(min / UINT64_MAX == FromInt64(INT64_MIN));
    EXPECT_TRUE(min % UINT64_MAX == FromInt64(INT64_MIN));
    EXPECT_TRUE(min / INT64_MAX == int128_t(0xFFFFFFFFFFFFFFFEull, 0xFFFFFFFFFFFFFFFEull));  // -(2^64+2)
    EXPECT_TRUE(min % INT64_MAX == int128_t(-2ll));

    // and the other way round: a smaller magnitude truncates to zero and is its own remainder
    EXPECT_TRUE((max / min).is_zero());
    EXPECT_TRUE(max % min == max);
}
// Division by every 64 bit scalar type, of both signs and every size, against the reference division.
// The int64_t and uint64_t forms of the same bits differ exactly when bit 63 is set, which is where a
// sign or zero extension mistake shows. A scalar on the left is checked as well.
TEST(int128_t, DivideByScalarMatchesReferenceDivision)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const uint64_t aLow = get_uint64_random();
        const uint64_t aHigh = get_uint64_random();
        const int32_t aShift = static_cast<int32_t>(get_uint32_random() % 128);
        const uint64_t divisorBits = get_uint64_random();
        const uint32_t divisorShift = get_uint32_random() % 64;

        // full width with either sign, at every magnitude
        const int128_t a = ReferenceShift(int128_t(aLow, aHigh), aShift, false);
        const uint64_t u = divisorBits >> divisorShift;
        const int64_t s = static_cast<int64_t>(u);  // the same bits, negative when bit 63 is set
        if (u == 0)
            continue;

        int128_t q, r;
        ReferenceDivide(a, FromInt64(s), q, r);
        EXPECT_TRUE(a / s == q) << "a=" << (std::string)a << ", s=" << s;
        EXPECT_TRUE(a % s == r) << "a=" << (std::string)a << ", s=" << s;
        int128_t c = a;
        c /= s;
        EXPECT_TRUE(c == q) << "a=" << (std::string)a << ", s=" << s;

        ReferenceDivide(a, int128_t(u, 0ull), q, r);
        EXPECT_TRUE(a / u == q) << "a=" << (std::string)a << ", u=" << u;
        EXPECT_TRUE(a % u == r) << "a=" << (std::string)a << ", u=" << u;
        c = a;
        c %= u;
        EXPECT_TRUE(c == r) << "a=" << (std::string)a << ", u=" << u;

        if (!a.is_zero()) {
            ReferenceDivide(FromInt64(s), a, q, r);
            EXPECT_TRUE(s / a == q) << "a=" << (std::string)a << ", s=" << s;
            EXPECT_TRUE(s % a == r) << "a=" << (std::string)a << ", s=" << s;
        }
    }
}
// The shift operators against the bit by bit reference at every count from -2 to 130, on values of
// both signs. The arithmetic right shift fills with the sign bit; counts of zero or below leave the
// value alone. Hidden behind volatile so that the runtime code is what gets tested.
TEST(int128_t, ShiftsMatchABitReference)
{
    const int128_t values[] = {int128_t(0x0123456789ABCDEFull, 0xFEDCBA9876543210ull),
                               int128_t(0x0123456789ABCDEFull, 0x7EDCBA9876543210ull),
                               int128_t(-1ll),
                               int128_t(1ll),
                               int128_t(0ull, 0x8000000000000000ull),
                               int128_t(~0ull, 0x7FFFFFFFFFFFFFFFull)};
    for (const auto& x : values) {
        for (int32_t s = -2; s <= 130; ++s) {
            volatile int32_t count = s;
            EXPECT_TRUE((x >> count) == ReferenceShift(x, s, false)) << "x=" << (std::string)x << ", shift=" << s;
            EXPECT_TRUE((x << count) == ReferenceShift(x, s, true)) << "x=" << (std::string)x << ", shift=" << s;
        }
        for (const int32_t s : {200, INT32_MAX, INT32_MIN}) {
            volatile int32_t count = s;
            EXPECT_TRUE((x >> count) == ReferenceShift(x, s, false)) << "x=" << (std::string)x << ", shift=" << s;
            EXPECT_TRUE((x << count) == ReferenceShift(x, s, true)) << "x=" << (std::string)x << ", shift=" << s;
        }
    }

    // a few results spelled out, from Python's arithmetic >> on the same value
    const int128_t x(0x0123456789ABCDEFull, 0xFEDCBA9876543210ull);
    EXPECT_TRUE((x >> 63) == int128_t(0xFDB97530ECA86420ull, ~0ull));
    EXPECT_TRUE((x >> 64) == int128_t(0xFEDCBA9876543210ull, ~0ull));
    EXPECT_TRUE((x >> 65) == int128_t(0xFF6E5D4C3B2A1908ull, ~0ull));
    EXPECT_TRUE((x >> 100) == int128_t(-1193047ll));
    EXPECT_TRUE((x << 65) == int128_t(0ull, 0x02468ACF13579BDEull));

    // a scalar on the left is widened to 128 bit first, then shifted
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const int64_t v = static_cast<int64_t>(get_uint64_random());
        const int32_t n = static_cast<int32_t>(get_uint32_random() % 128);
        const int128_t count(n);
        EXPECT_TRUE((v << count) == ReferenceShift(FromInt64(v), n, true)) << "v=" << v << ", n=" << n;
        EXPECT_TRUE((v >> count) == ReferenceShift(FromInt64(v), n, false)) << "v=" << v << ", n=" << n;
    }
}
// The scalar overloads of operator*= skip one partial product for an operand that fits in 64 bits. The
// binary form with the scalar on the right had never run, and neither had int64_t operands or an
// unsigned one with bit 63 set; dropping the carry of the low product for such a scalar passed both
// integer suites. Every scalar type is checked against ReferenceMultiply() of its 128 bit extension.
TEST(int128_t, MultiplyByScalar)
{
    const auto check = [](const int128_t& a, auto scalar) {
        using S = decltype(scalar);
        int128_t wide;
        if constexpr (std::is_signed_v<S>) {
            wide = FromInt64(static_cast<int64_t>(scalar));
        } else {
            wide = int128_t(static_cast<uint64_t>(scalar), 0ull);
        }
        const int128_t expected = ReferenceMultiply(a, wide);
        EXPECT_TRUE(a * scalar == expected) << "a=" << (std::string)a << ", scalar=" << scalar;
        EXPECT_TRUE(scalar * a == expected) << "a=" << (std::string)a << ", scalar=" << scalar;
        int128_t c = a;
        c *= scalar;
        EXPECT_TRUE(c == expected) << "a=" << (std::string)a << ", scalar=" << scalar;
    };

    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const uint64_t aLow = get_uint64_random();
        const uint64_t aHigh = get_uint64_random();
        const uint64_t bits = get_uint64_random();
        const int128_t a(aLow, aHigh);
        check(a, static_cast<int64_t>(bits));
        check(a, bits);
        check(a, static_cast<int32_t>(bits));
        check(a, static_cast<uint32_t>(bits));
    }

    const int128_t values[] = {int128_t(-1ll), int128_t(0ull, 0x8000000000000000ull), int128_t(~0ull, 0x7FFFFFFFFFFFFFFFull),
                               int128_t(0x0123456789ABCDEFull, 0xFEDCBA9876543210ull)};
    const int64_t signedScalars[] = {0, 1, -1, 2, -2, INT32_MIN, INT32_MAX, INT64_MIN, INT64_MIN + 1, INT64_MAX};
    const uint64_t unsignedScalars[] = {0ull, 1ull, UINT32_MAX, 0x7FFFFFFFFFFFFFFFull, 0x8000000000000000ull, 0x8000000000000001ull, UINT64_MAX};
    for (const auto& a : values) {
        for (const int64_t s : signedScalars) {
            check(a, s);
        }
        for (const uint64_t u : unsignedScalars) {
            check(a, u);
        }
    }
}
// The scalar forms of the bitwise operators had never been instantiated. The scalar is widened first:
// sign extended when it is signed, zero extended when it is not, so -256 and 0xFF...00ull clear the
// same low bits but only the unsigned one clears the high QWORD.
TEST(int128_t, BitwiseWithScalarOperandsExtendTheScalar)
{
    const int128_t pattern(0xFFFFFFFFFFFFFFFFull, 0x1234ull);
    EXPECT_TRUE((pattern & -256) == int128_t(0xFFFFFFFFFFFFFF00ull, 0x1234ull));
    EXPECT_TRUE((pattern & 0xFFFFFFFFFFFFFF00ull) == int128_t(0xFFFFFFFFFFFFFF00ull, 0ull));
    EXPECT_TRUE((pattern | INT64_MIN) == int128_t(-1ll));
    EXPECT_TRUE((pattern ^ -1) == int128_t(0ull, 0xFFFFFFFFFFFFEDCBull));
    EXPECT_TRUE((-256 & pattern) == int128_t(0xFFFFFFFFFFFFFF00ull, 0x1234ull));
    EXPECT_TRUE((0xFFFFFFFFFFFFFF00ull & pattern) == int128_t(0xFFFFFFFFFFFFFF00ull, 0ull));
    int128_t y = pattern;
    y &= -256;
    EXPECT_TRUE(y == int128_t(0xFFFFFFFFFFFFFF00ull, 0x1234ull));
    y = pattern;
    y |= -1;
    EXPECT_TRUE(y == int128_t(-1ll));
    y = pattern;
    y ^= 1u;
    EXPECT_TRUE(y == int128_t(0xFFFFFFFFFFFFFFFEull, 0x1234ull));

    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const uint64_t al = get_uint64_random();
        const uint64_t ah = get_uint64_random();
        const uint64_t bits = get_uint64_random();
        const int128_t a(al, ah);
        const int64_t s = static_cast<int64_t>(bits);
        const uint64_t sh = SignExtension(s);

        EXPECT_TRUE((a & s) == int128_t(al & bits, ah & sh)) << "a=" << (std::string)a << ", s=" << s;
        EXPECT_TRUE((a | s) == int128_t(al | bits, ah | sh)) << "a=" << (std::string)a << ", s=" << s;
        EXPECT_TRUE((a ^ s) == int128_t(al ^ bits, ah ^ sh)) << "a=" << (std::string)a << ", s=" << s;
        EXPECT_TRUE((s & a) == int128_t(al & bits, ah & sh)) << "a=" << (std::string)a << ", s=" << s;
        EXPECT_TRUE((a & bits) == int128_t(al & bits, 0ull)) << "a=" << (std::string)a << ", u=" << bits;
        EXPECT_TRUE((a | bits) == int128_t(al | bits, ah)) << "a=" << (std::string)a << ", u=" << bits;
        EXPECT_TRUE((bits ^ a) == int128_t(al ^ bits, ah)) << "a=" << (std::string)a << ", u=" << bits;

        int128_t c = a;
        c &= s;
        EXPECT_TRUE(c == int128_t(al & bits, ah & sh)) << "a=" << (std::string)a << ", s=" << s;
        c = a;
        c |= bits;
        EXPECT_TRUE(c == int128_t(al | bits, ah)) << "a=" << (std::string)a << ", u=" << bits;
        c = a;
        c ^= s;
        EXPECT_TRUE(c == int128_t(al ^ bits, ah ^ sh)) << "a=" << (std::string)a << ", s=" << s;
    }
}
// Comparisons with an integral scalar widen the scalar, which makes them exact at 128 bits - unlike
// the builtin types, where -1 < UINT64_MAX is false. x <= scalar and x >= scalar had never run
// outside a constant expression, and computing <= as >= survived both integer suites. The reference
// orders the QWORDs directly: the high ones as signed values, then the low ones as unsigned. Every
// other iteration gives both sides the same high QWORD, so that the low one decides.
TEST(int128_t, CompareWithIntegralScalars)
{
    const int128_t max(~0ull, 0x7FFFFFFFFFFFFFFFull);
    const int128_t min(0ull, 0x8000000000000000ull);
    EXPECT_TRUE(int128_t(-2) <= -1);
    EXPECT_FALSE(int128_t(0) <= -1);
    EXPECT_FALSE(int128_t(-1) >= 0);
    EXPECT_TRUE(int128_t(0) >= -1);
    EXPECT_TRUE(int128_t(-1) < UINT64_MAX);
    EXPECT_FALSE(int128_t(-1) == UINT64_MAX);
    EXPECT_TRUE(UINT64_MAX > int128_t(-1));
    EXPECT_TRUE(int128_t(~0ull, 0ull) == UINT64_MAX);
    EXPECT_TRUE(int128_t(-1) <= 0u);
    EXPECT_FALSE(int128_t(-1) >= 0u);
    EXPECT_TRUE(min < INT64_MIN);
    EXPECT_TRUE(max > UINT64_MAX);
    EXPECT_TRUE(FromInt64(INT64_MIN) == INT64_MIN);
    EXPECT_TRUE(INT64_MIN >= FromInt64(INT64_MIN));

    const auto check = [](const int128_t& a, auto scalar, uint64_t sl, uint64_t sh) {
        uint64_t al = 0, ah = 0;
        a.get_components(al, ah);
        const bool less = static_cast<int64_t>(ah) < static_cast<int64_t>(sh) || (ah == sh && al < sl);
        const bool equal = ah == sh && al == sl;
        const bool greater = !less && !equal;
        EXPECT_EQ(a == scalar, equal) << "a=" << (std::string)a << ", scalar=" << scalar;
        EXPECT_EQ(a != scalar, !equal) << "a=" << (std::string)a << ", scalar=" << scalar;
        EXPECT_EQ(a < scalar, less) << "a=" << (std::string)a << ", scalar=" << scalar;
        EXPECT_EQ(a <= scalar, less || equal) << "a=" << (std::string)a << ", scalar=" << scalar;
        EXPECT_EQ(a > scalar, greater) << "a=" << (std::string)a << ", scalar=" << scalar;
        EXPECT_EQ(a >= scalar, greater || equal) << "a=" << (std::string)a << ", scalar=" << scalar;
        EXPECT_EQ(scalar == a, equal) << "a=" << (std::string)a << ", scalar=" << scalar;
        EXPECT_EQ(scalar != a, !equal) << "a=" << (std::string)a << ", scalar=" << scalar;
        EXPECT_EQ(scalar < a, greater) << "a=" << (std::string)a << ", scalar=" << scalar;
        EXPECT_EQ(scalar <= a, greater || equal) << "a=" << (std::string)a << ", scalar=" << scalar;
        EXPECT_EQ(scalar > a, less) << "a=" << (std::string)a << ", scalar=" << scalar;
        EXPECT_EQ(scalar >= a, less || equal) << "a=" << (std::string)a << ", scalar=" << scalar;
    };

    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const uint64_t al = get_uint64_random();
        const uint64_t ah = get_uint64_random();
        const uint64_t bits = get_uint64_random();
        const bool sameHigh = (i & 1) != 0;
        const bool sameLow = (i & 6) == 0;  // one iteration in four of those compares equal values

        const int64_t s = static_cast<int64_t>(bits);
        check(int128_t(sameLow ? bits : al, sameHigh ? SignExtension(s) : ah), s, bits, SignExtension(s));
        check(int128_t(sameLow ? bits : al, sameHigh ? 0ull : ah), bits, bits, 0ull);
        const int32_t s32 = static_cast<int32_t>(bits);
        const uint64_t s32Low = static_cast<uint64_t>(static_cast<int64_t>(s32));
        check(int128_t(sameLow ? s32Low : al, sameHigh ? SignExtension(s32) : ah), s32, s32Low, SignExtension(s32));
    }
}
// Comparisons with a floating point scalar compare real numbers. They used to convert the scalar to
// int128_t first, which truncated it towards zero: int128_t(-2) == -2.5 held, int128_t(-2) > -2.5 did
// not, and NaN, the infinities and 1e39 saturated and compared equal to the ends of the range. Each
// row gives the order of one value against every scalar column, from Python's exact comparison of an
// int with a float: '<' the value is below the scalar, '=' equal, '>' above, 'u' unordered. Integral
// scalars are still widened instead, see CompareWithIntegralScalars.
TEST(int128_t, CompareWithFloatingPointIsExact)
{
    const double inf = std::numeric_limits<double>::infinity();
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double doubles[] = {0.0,     -0.0,    0.5,     -0.5,     1.0,     -1.0,     2.5,     -2.5,     -2.0,
                              0x1p53,  -0x1p53, 0x1.0000000000001p53,   -0x1.0000000000001p53,   0x1p63,   -0x1p63,
                              0x1p64,  -0x1p64, 0x1.fffffffffffffp126,  -0x1.fffffffffffffp126,  0x1p127,  -0x1p127,
                              -0x1.0000000000001p127,  1e39,    -1e39,   DBL_MAX,  -DBL_MAX, DBL_MIN,  -DBL_MIN, inf,
                              -inf,    nan};
    const float finf = std::numeric_limits<float>::infinity();
    const float fnan = std::numeric_limits<float>::quiet_NaN();
    const float floats[] = {0.0f,    -0.0f,    2.5f,     -2.5f,    0x1p24f, -0x1p24f, 0x1.000002p24f, 0x1p64f,
                            -0x1p64f, 0x1.fffffep126f,  0x1p127f, -0x1p127f, finf,    -finf,    fnan};

    const struct {
        uint64_t low, high;
        const char* vsDouble;
        const char* vsFloat;
    } rows[] = {
        {0x0000000000000000ull, 0x0000000000000000ull, "==<><><>><><><><><><>><><><><>u", "==<><><<><<><>u"},  // 0
        {0x0000000000000001ull, 0x0000000000000000ull, ">>>>=><>><><><><><><>><><>>><>u", ">><><><<><<><>u"},  // 1
        {0xFFFFFFFFFFFFFFFFull, 0xFFFFFFFFFFFFFFFFull, "<<<<<=<>><><><><><><>><><><<<>u", "<<<><><<><<><>u"},  // -1
        {0x0000000000000002ull, 0x0000000000000000ull, ">>>>>><>><><><><><><>><><>>><>u", ">><><><<><<><>u"},  // 2
        {0xFFFFFFFFFFFFFFFEull, 0xFFFFFFFFFFFFFFFFull, "<<<<<<<>=<><><><><><>><><><<<>u", "<<<><><<><<><>u"},  // -2
        {0x0000000000000003ull, 0x0000000000000000ull, ">>>>>>>>><><><><><><>><><>>><>u", ">>>><><<><<><>u"},  // 3
        {0xFFFFFFFFFFFFFFFDull, 0xFFFFFFFFFFFFFFFFull, "<<<<<<<<<<><><><><><>><><><<<>u", "<<<<<><<><<><>u"},  // -3
        {0x0020000000000000ull, 0x0000000000000000ull, ">>>>>>>>>=><><><><><>><><>>><>u", ">>>>>>><><<><>u"},  // 2^53
        {0x0020000000000001ull, 0x0000000000000000ull, ">>>>>>>>>>><><><><><>><><>>><>u", ">>>>>>><><<><>u"},  // 2^53+1
        {0xFFE0000000000000ull, 0xFFFFFFFFFFFFFFFFull, "<<<<<<<<<<=<><><><><>><><><<<>u", "<<<<<<<<><<><>u"},  // -2^53
        {0xFFDFFFFFFFFFFFFFull, 0xFFFFFFFFFFFFFFFFull, "<<<<<<<<<<<<><><><><>><><><<<>u", "<<<<<<<<><<><>u"},  // -(2^53+1)
        {0x8000000000000000ull, 0x0000000000000000ull, ">>>>>>>>>>>>>=><><><>><><>>><>u", ">>>>>>><><<><>u"},  // 2^63
        {0x8000000000000000ull, 0xFFFFFFFFFFFFFFFFull, "<<<<<<<<<<<<<<=<><><>><><><<<>u", "<<<<<<<<><<><>u"},  // -2^63
        {0x0000000000000000ull, 0x0000000000000001ull, ">>>>>>>>>>>>>>>=><><>><><>>><>u", ">>>>>>>=><<><>u"},  // 2^64
        {0x0000000000000000ull, 0xFFFFFFFFFFFFFFFFull, "<<<<<<<<<<<<<<<<=<><>><><><<<>u", "<<<<<<<<=<<><>u"},  // -2^64
        {0xFFFFFFFFFFFFFFFFull, 0xFFFFFFFFFFFFFFFEull, "<<<<<<<<<<<<<<<<<<><>><><><<<>u", "<<<<<<<<<<<><>u"},  // -2^64-1
        {0x0000000000000000ull, 0x7FFFFFFFFFFFFC00ull, ">>>>>>>>>>>>>>>>>=><>><><>>><>u", ">>>>>>>>>><><>u"},  // 2^127-2^74
        {0x0000000000000001ull, 0x7FFFFFFFFFFFFC00ull, ">>>>>>>>>>>>>>>>>>><>><><>>><>u", ">>>>>>>>>><><>u"},  // 2^127-2^74+1
        {0x0000000000000000ull, 0x8000000000000400ull, "<<<<<<<<<<<<<<<<<<=<>><><><<<>u", "<<<<<<<<<<<><>u"},  // -(2^127-2^74)
        {0x0000000000000000ull, 0x7FFFFF8000000000ull, ">>>>>>>>>>>>>>>>><><>><><>>><>u", ">>>>>>>>>=<><>u"},  // (2^24-1)*2^103
        {0xFFFFFFFFFFFFFFFFull, 0x7FFFFFFFFFFFFFFFull, ">>>>>>>>>>>>>>>>>>><>><><>>><>u", ">>>>>>>>>><><>u"},  // MAX
        {0x0000000000000001ull, 0x8000000000000000ull, "<<<<<<<<<<<<<<<<<<<<>><><><<<>u", "<<<<<<<<<<<><>u"},  // MIN+1
        {0x0000000000000000ull, 0x8000000000000000ull, "<<<<<<<<<<<<<<<<<<<<=><><><<<>u", "<<<<<<<<<<<=<>u"},  // MIN
    };

    // all six operators, with the scalar on either side, against the expected relation
    const auto check = [](const int128_t& v, auto scalar, char relation) {
        const bool less = relation == '<', equal = relation == '=', greater = relation == '>';
        EXPECT_EQ(v == scalar, equal) << "v=" << (std::string)v << ", scalar=" << std::hexfloat << scalar;
        EXPECT_EQ(v != scalar, !equal) << "v=" << (std::string)v << ", scalar=" << std::hexfloat << scalar;
        EXPECT_EQ(v < scalar, less) << "v=" << (std::string)v << ", scalar=" << std::hexfloat << scalar;
        EXPECT_EQ(v <= scalar, less || equal) << "v=" << (std::string)v << ", scalar=" << std::hexfloat << scalar;
        EXPECT_EQ(v > scalar, greater) << "v=" << (std::string)v << ", scalar=" << std::hexfloat << scalar;
        EXPECT_EQ(v >= scalar, greater || equal) << "v=" << (std::string)v << ", scalar=" << std::hexfloat << scalar;
        EXPECT_EQ(scalar == v, equal) << "v=" << (std::string)v << ", scalar=" << std::hexfloat << scalar;
        EXPECT_EQ(scalar != v, !equal) << "v=" << (std::string)v << ", scalar=" << std::hexfloat << scalar;
        EXPECT_EQ(scalar < v, greater) << "v=" << (std::string)v << ", scalar=" << std::hexfloat << scalar;
        EXPECT_EQ(scalar <= v, greater || equal) << "v=" << (std::string)v << ", scalar=" << std::hexfloat << scalar;
        EXPECT_EQ(scalar > v, less) << "v=" << (std::string)v << ", scalar=" << std::hexfloat << scalar;
        EXPECT_EQ(scalar >= v, less || equal) << "v=" << (std::string)v << ", scalar=" << std::hexfloat << scalar;
    };

    for (const auto& row : rows) {
        const int128_t v(row.low, row.high);
        ASSERT_EQ(strlen(row.vsDouble), std::extent_v<decltype(doubles)>);
        ASSERT_EQ(strlen(row.vsFloat), std::extent_v<decltype(floats)>);
        for (size_t j = 0; j < std::extent_v<decltype(doubles)>; ++j) {
            check(v, doubles[j], row.vsDouble[j]);
        }
        for (size_t j = 0; j < std::extent_v<decltype(floats)>; ++j) {
            check(v, floats[j], row.vsFloat[j]);
        }
    }
    // a long double goes through double, which on every supported compiler holds these exactly
    check(int128_t(-3ll), -2.5L, '<');
    check(int128_t(-2ll), -2.5L, '>');

    // the same in constant evaluation
    constexpr int128_t max(~0ull, 0x7FFFFFFFFFFFFFFFull);
    constexpr int128_t min(0ull, 0x8000000000000000ull);
    constexpr double cnan = std::numeric_limits<double>::quiet_NaN();
    constexpr double cinf = std::numeric_limits<double>::infinity();
    static_assert(int128_t(-2) != -2.5 && int128_t(-2) > -2.5 && int128_t(-3) < -2.5 && -2.5 < int128_t(-2) && -2.5f > int128_t(-3));
    static_assert(!(int128_t(-1) == cnan) && int128_t(-1) != cnan && !(int128_t(-1) < cnan) && !(cnan >= int128_t(-1)));
    static_assert(max < cinf && max < 1e39 && max < 0x1p127 && min > -cinf && min > -1e39 && min == -0x1p127 && min > -0x1.0000000000001p127);
    static_assert(int128_t() == -0.0 && int128_t() > -0.5 && int128_t() < 0.5);
}
// operator*=(T) used to accept any operand that converted to uint64_t, so a uint128_t, a float128 or a
// fixed_point128 went through a conversion operator and was truncated: int128_t(3) * uint128_t(2^64+1)
// gave 3. Such operands now fail to compile, as they already did for the other arithmetic operators,
// and the builtin operands keep working.
TEST(int128_t, MultiplyRejectsTheOtherSignedness)
{
    static_assert(!Multipliable<int128_t, uint128_t>, "int128_t * uint128_t must not compile");
    static_assert(!MultiplyAssignable<int128_t, uint128_t>, "int128_t *= uint128_t must not compile");
    static_assert(!MultiplyAssignable<int128_t, float128>, "int128_t *= float128 must not compile");
    static_assert(!MultiplyAssignable<int128_t, fixed_point128<20>>, "int128_t *= fixed_point128 must not compile");
    static_assert(Multipliable<int128_t, int128_t> && Multipliable<int128_t, int> && Multipliable<int128_t, int64_t>);
    static_assert(Multipliable<int128_t, uint64_t> && Multipliable<int128_t, double> && Multipliable<int, int128_t>);
    static_assert(Multipliable<double, int128_t> && MultiplyAssignable<int128_t, int64_t> && MultiplyAssignable<int128_t, double>);
    static_assert(MultiplyAssignable<int128_t, int128_t> && MultiplyAssignable<int128_t, uint32_t>);

    // and the accepted forms compute the full product
    EXPECT_TRUE(int128_t(3) * int128_t(1ull, 1ull) == int128_t(3ull, 3ull));
    EXPECT_TRUE(int128_t(-1) * UINT64_MAX == int128_t(1ull, ~0ull));  // -(2^64-1)
    EXPECT_TRUE(-3 * int128_t(1ull, 1ull) == -int128_t(3ull, 3ull));
}
// Whether an operation may throw is part of the interface: generic code branches on it, and a
// container only moves elements whose move constructor is noexcept. Division, modulo and the
// logarithms throw on a zero or out of domain operand and must not claim otherwise.
TEST(int128_t, NoexceptContracts)
{
    int128_t a, b;

    static_assert(noexcept(a + b) && noexcept(a - b) && noexcept(a * b));
    static_assert(noexcept(a + 1) && noexcept(1 + a) && noexcept(a * -1) && noexcept(-2 * a) && noexcept(a - 1) && noexcept(1 - a));
    static_assert(noexcept(a += b) && noexcept(a -= b) && noexcept(a *= b) && noexcept(a &= b) && noexcept(a |= b) && noexcept(a ^= b));
    static_assert(noexcept(a & b) && noexcept(a | b) && noexcept(a ^ b) && noexcept(a & -1) && noexcept(a | 1u) && noexcept(a ^ 1ull));
    static_assert(noexcept(-1 & a) && noexcept(1u | a) && noexcept(1ull ^ a));
    static_assert(noexcept(a << 3) && noexcept(a >> 3) && noexcept(a <<= 3) && noexcept(a >>= 3) && noexcept(1 << a) && noexcept(-1 >> a));
    static_assert(noexcept(a == b) && noexcept(a != b) && noexcept(a < b) && noexcept(a <= b) && noexcept(a > b) && noexcept(a >= b));
    static_assert(noexcept(a == -1) && noexcept(-1 == a) && noexcept(a < -1) && noexcept(-1 < a) && noexcept(a <= 2.5) && noexcept(2.5 >= a));
    static_assert(noexcept(~a) && noexcept(-a) && noexcept(+a) && noexcept(!a) && noexcept(abs(a)));
    static_assert(noexcept(++a) && noexcept(a++) && noexcept(--a) && noexcept(a--));
    static_assert(noexcept(sqr(a)) && noexcept(sqrt(a)) && noexcept(pow(a, 2u)) && noexcept(lzcnt128(a)));
    static_assert(noexcept(static_cast<double>(a)) && noexcept(static_cast<float>(a)) && noexcept(static_cast<int64_t>(a)) && noexcept(a.hex()));
    static_assert(noexcept(int128_t(-1)) && noexcept(int128_t(-1.5)) && noexcept(int128_t("-1")));
    static_assert(std::is_nothrow_copy_constructible_v<int128_t> && std::is_nothrow_move_constructible_v<int128_t>);
    static_assert(std::is_nothrow_copy_assignable_v<int128_t> && std::is_nothrow_move_assignable_v<int128_t>);

    // the operations that throw
    static_assert(!noexcept(a / b) && !noexcept(a % b) && !noexcept(a / -2) && !noexcept(a % -2) && !noexcept(a /= b) && !noexcept(a %= b));
    static_assert(!noexcept(-2 / a) && !noexcept(-2 % a) && !noexcept(a / 2.0));
    static_assert(!noexcept(log2(a)) && !noexcept(log(a)) && !noexcept(log10(a)));
    SUCCEED();
}
// The negative branch of log() and log10() never ran. A negative argument is outside the domain, as
// zero is, and the most negative value is no exception.
TEST(int128_t, LogFunctionsRejectNegatives)
{
    const int128_t min(0ull, 0x8000000000000000ull);
    EXPECT_THROW((void)log(int128_t(-1ll)), std::domain_error);
    EXPECT_THROW((void)log10(int128_t(-1ll)), std::domain_error);
    EXPECT_THROW((void)log(min), std::domain_error);
    EXPECT_THROW((void)log10(min), std::domain_error);
    EXPECT_THROW((void)log2(min), std::domain_error);
    // sqrt() is noexcept and documented to answer zero instead
    EXPECT_EQ(sqrt(min), 0ull);
    EXPECT_EQ(sqrt(int128_t(-5ll)), 0ull);
}
// log10() at both sides of every power of ten that fits, and log() at both sides of ceil(e^k) - the
// smallest integer whose natural log is k - for a spread of k. The random tests cannot reach these: a
// uniform operand lands on a boundary with negligible probability. The ceil(e^k) values come from
// Python's decimal module at 80 significant digits.
TEST(int128_t, LogAndLog10Boundaries)
{
    const auto powers = PowersOfTen();
    for (uint64_t k = 1; k < powers.size(); ++k) {
        EXPECT_EQ(log10(powers[k]), k) << "10^" << k;
        EXPECT_EQ(log10(powers[k] - int128_t(1ll)), k - 1) << "10^" << k << " - 1";
    }

    const struct {
        const char* boundary;
        uint64_t expo;
    } boundaries[] = {
        {"3", 1},
        {"8", 2},
        {"22027", 10},
        {"235385266837019986", 40},
        {"4727839468229346562", 43},
        {"12851600114359308276", 44},  // above 2^63, still in the low QWORD
        {"2293783159469609879099352841", 63},
        {"6235149080811616882909238709", 64},
        {"22352466037347150474430657323327147399", 86},
        {"60760302250568721495223289381302760753", 87},
        {"165163625499400185552832979626485876707", 88},
    };
    for (const auto& b : boundaries) {
        const int128_t x(b.boundary);
        EXPECT_EQ(log(x), b.expo) << "at ceil(e^" << b.expo << ") = " << b.boundary;
        EXPECT_EQ(log(x - int128_t(1ll)), b.expo - 1) << "just below ceil(e^" << b.expo << ") = " << b.boundary;
    }

    // the top of the range: 10^38 <= max < 10^39, e^88 <= max < e^89, isqrt from Python's math.isqrt
    const int128_t max(~0ull, 0x7FFFFFFFFFFFFFFFull);
    EXPECT_EQ(log10(max), 38ull);
    EXPECT_EQ(log(max), 88ull);
    EXPECT_EQ(log2(max), 126ull);
    EXPECT_EQ(sqrt(max), 13043817825332782212ull);
}
// The parse's sign, prefix and wrap around corners, values from Python reduced modulo 2^128. Overflow
// wraps like the builtin types, so 2^127 parses as the most negative value.
TEST(int128_t, ConstructorFromStringSignsAndWrap)
{
    const int128_t max(~0ull, 0x7FFFFFFFFFFFFFFFull);
    const int128_t min(0ull, 0x8000000000000000ull);
    EXPECT_TRUE(int128_t("-0x80000000000000000000000000000000") == min);
    EXPECT_TRUE(int128_t("170141183460469231731687303715884105728") == min);    // 2^127 wraps
    EXPECT_TRUE(int128_t("-170141183460469231731687303715884105729") == max);   // -2^127-1 wraps
    EXPECT_TRUE(int128_t("340282366920938463463374607431768211457") == int128_t(1ll));  // 2^128+1
    // a sign without digits, or a second sign, parses as nothing
    for (const char* text : {"+-5", "--5", "-", "+", "-0x", "- 5", "-0x-5", "0x"}) {
        EXPECT_TRUE(int128_t(text).is_zero()) << "text=" << text;
    }
    // parsing stops at the first character that is not a digit of the detected base
    EXPECT_TRUE(int128_t("-1.5") == int128_t(-1ll));
    EXPECT_TRUE(int128_t("1e5") == int128_t(1ll));
    EXPECT_TRUE(int128_t("-0XFF") == int128_t(-255ll));  // the prefix is case insensitive
    EXPECT_TRUE(int128_t("\t-42") == int128_t(-42ll));
    EXPECT_TRUE(int128_t("-000123") == int128_t(-123ll));
    EXPECT_TRUE(int128_t("-'5") == int128_t(-5ll));     // a digit separator is skipped
    EXPECT_TRUE(int128_t(std::string("-0x8000000000000000")) == FromInt64(INT64_MIN));
    // the literal: the minus is applied to 2^127, which wraps to itself
    EXPECT_TRUE(-170141183460469231731687303715884105728_int128 == min);
}
// The generic assignment and the narrow integral constructors had only been evaluated at compile
// time. opaque() keeps the operands out of the optimizer's sight.
TEST(int128_t, AssignmentAndConstructionFromScalars)
{
    int128_t x;
    x = opaque(-5);
    EXPECT_TRUE(x == FromInt64(-5));
    x = opaque(5ull);
    EXPECT_TRUE(x == int128_t(5ull, 0ull));
    x = opaque(UINT64_MAX);
    EXPECT_TRUE(x == int128_t(UINT64_MAX, 0ull));  // zero extended
    x = opaque(INT64_MIN);
    EXPECT_TRUE(x == FromInt64(INT64_MIN));        // sign extended
    x = opaque(-2.5);
    EXPECT_TRUE(x == FromInt64(-2));
    x = opaque(0x1p100);
    EXPECT_TRUE(x == int128_t(0ull, 1ull << 36));
    const char* text = "-123";
    x = text;
    EXPECT_TRUE(x == FromInt64(-123));

    EXPECT_TRUE(int128_t(opaque(static_cast<signed char>(-3))) == FromInt64(-3));
    EXPECT_TRUE(int128_t(opaque(static_cast<char>(65))) == FromInt64(65));
    EXPECT_TRUE(int128_t(opaque(static_cast<unsigned char>(255))) == FromInt64(255));
    EXPECT_TRUE(int128_t(opaque(static_cast<short>(-32768))) == FromInt64(-32768));
    EXPECT_TRUE(int128_t(opaque(static_cast<unsigned short>(65535))) == FromInt64(65535));
    EXPECT_TRUE(int128_t(opaque(static_cast<char16_t>(0xFFFF))) == FromInt64(65535));
    EXPECT_TRUE(int128_t(opaque(static_cast<char32_t>(0xFFFFFFFF))) == int128_t(0xFFFFFFFFull, 0ull));
    EXPECT_TRUE(int128_t(opaque(true)) == FromInt64(1));
}
// A scalar on the left is widened before the operation, so the result is exact at 128 bits even where
// the same expression in int64_t would overflow: INT64_MIN / -1 is 2^63 here.
TEST(int128_t, ScalarOnTheLeftWith64BitExtremes)
{
    EXPECT_TRUE(UINT64_MAX + int128_t(1) == int128_t(0ull, 1ull));
    EXPECT_TRUE(0ull - int128_t(1) == int128_t(-1ll));
    EXPECT_TRUE(5.0 - int128_t(10) == int128_t(-5ll));
    EXPECT_TRUE(-7 / int128_t(2) == int128_t(-3ll));
    EXPECT_TRUE(-7 % int128_t(2) == int128_t(-1ll));
    EXPECT_TRUE(INT64_MIN / int128_t(-1) == int128_t(0x8000000000000000ull, 0ull));
    EXPECT_TRUE(INT64_MIN * int128_t(-1) == int128_t(0x8000000000000000ull, 0ull));
    EXPECT_TRUE(INT64_MIN - int128_t(1) == int128_t(0x7FFFFFFFFFFFFFFFull, ~0ull));
    EXPECT_TRUE((-1 << int128_t(64)) == int128_t(0ull, ~0ull));
    EXPECT_TRUE((-1ll >> int128_t(200)) == int128_t(-1ll));
    EXPECT_TRUE((UINT64_MAX & int128_t(-1)) == int128_t(UINT64_MAX, 0ull));
}
// The conversions to the builtin integer types keep the low bits, as a builtin narrowing does.
TEST(int128_t, NarrowingConversionsTruncate)
{
    const int128_t v(0x0000000500000007ull, 0x123ull);
    EXPECT_EQ(static_cast<int64_t>(v), 21474836487ll);
    EXPECT_EQ(static_cast<uint64_t>(v), 21474836487ull);
    EXPECT_EQ(static_cast<int32_t>(v), 7);
    EXPECT_EQ(static_cast<uint32_t>(v), 7u);

    const int128_t min(0ull, 0x8000000000000000ull);
    EXPECT_EQ(static_cast<int64_t>(min), 0ll);
    EXPECT_TRUE(static_cast<bool>(min));
    EXPECT_EQ(static_cast<int32_t>(int128_t(-1ll)), -1);
    EXPECT_EQ(static_cast<uint32_t>(int128_t(-1ll)), 4294967295u);
    EXPECT_EQ(static_cast<uint64_t>(int128_t(-1ll)), UINT64_MAX);
    EXPECT_EQ(static_cast<int64_t>(int128_t(0x8000000000000000ull)), INT64_MIN);  // 2^63, positive, reads back negative
    EXPECT_EQ(static_cast<int32_t>(int128_t(0x80000000ull)), INT32_MIN);
}
// A floating point divisor is converted first, truncating towards zero as the constructor does, which
// is the semantics the multiplication has too. A magnitude below one converts to zero, and dividing by
// it throws.
TEST(int128_t, DivideAndModuloByFloatingPoint)
{
    EXPECT_TRUE(int128_t(-10) / 2.0 == int128_t(-5ll));
    EXPECT_TRUE(int128_t(-10) % 3.0 == int128_t(-1ll));
    EXPECT_TRUE(int128_t(-10) / -2.5 == int128_t(5ll));  // -2.5 converts to -2
    int128_t c(-1000);
    c /= 4.0f;
    EXPECT_TRUE(c == int128_t(-250ll));
    int128_t d(-10);
    d *= -2.5;
    EXPECT_TRUE(d == int128_t(20ll));

    int128_t big = -int128_t(0ull, 1ull);  // -2^64
    big /= 1e19;
    EXPECT_TRUE(big == int128_t(-1ll));
    int128_t product(-3);
    product *= 1e30;  // the double nearest 1e30 is 1000000000000000019884624838656
    EXPECT_TRUE(product == ReferenceMultiply(int128_t(-3ll), ExpectedFromDouble(1e30)));

    EXPECT_THROW((void)(int128_t(10) / 0.5), std::logic_error);
    EXPECT_THROW((void)(int128_t(10) % 0.5), std::logic_error);
    EXPECT_THROW((void)(int128_t(-10) / -0.9), std::logic_error);
    int128_t e(10);
    EXPECT_THROW(e /= 0.25f, std::logic_error);
}
// The most negative value is its own negation and absolute value, and the compound operators accept an
// operand that aliases the object. pow() values from Python, reduced modulo 2^128.
TEST(int128_t, MostNegativeAndAliasing)
{
    const int128_t max(~0ull, 0x7FFFFFFFFFFFFFFFull);
    const int128_t min(0ull, 0x8000000000000000ull);
    EXPECT_TRUE(abs(min) == min);
    EXPECT_TRUE(abs(min + int128_t(1ll)) == max);
    EXPECT_TRUE(abs(int128_t(-1ll)) == int128_t(1ll));

    const int128_t values[] = {int128_t(-7ll), int128_t(7ll), min, max, int128_t(0x0123456789ABCDEFull, 0xFEDCBA9876543210ull)};
    for (const auto& v : values) {
        int128_t a = v;
        const int128_t& alias = a;
        a /= alias;
        EXPECT_TRUE(a == int128_t(1ll)) << "v=" << (std::string)v;
        a = v;
        a %= alias;
        EXPECT_TRUE(a.is_zero()) << "v=" << (std::string)v;
        a = v;
        a *= alias;
        EXPECT_TRUE(a == ReferenceMultiply(v, v)) << "v=" << (std::string)v;
        a = v;
        a /= -v;
        // -min is min, so that quotient is 1; every other value divided by its negation is -1
        EXPECT_TRUE(a == ((v == min) ? int128_t(1ll) : int128_t(-1ll))) << "v=" << (std::string)v;
    }

    EXPECT_TRUE(pow(int128_t(-2ll), 127u) == min);
    EXPECT_TRUE(pow(int128_t(-2ll), 128u).is_zero());
    EXPECT_TRUE(pow(int128_t(-2ll), 126u) == int128_t(0ull, 0x4000000000000000ull));
    EXPECT_TRUE(pow(min, 2u).is_zero());
    EXPECT_TRUE(pow(min, 1u) == min);
    EXPECT_TRUE(pow(int128_t(-1ll), UINT32_MAX) == int128_t(-1ll));
    EXPECT_TRUE(pow(int128_t(-3ll), 80u) == int128_t(0x3CEA59789C79D441ull, 0x6F32F1EF8B18A2BCull));
    EXPECT_TRUE(pow(int128_t(-7ll), 45u) == int128_t(0xA6EE598222C2CA19ull, 0xAF7F38482F1CE458ull));
}

/**********************************************************************
 * Compile time (constexpr) evaluation
 *
 * Each test below pairs static_asserts, which fail the build the moment one of these operations
 * stops being usable in a constant expression, with a runtime check of the same expression built
 * from opaque() values the optimizer cannot fold. See the matching block in
 * fixed_point128_gtest.cpp for why the runtime half is not a restatement of the compile time one.
 *
 * The two instantiations share most of their code, so these concentrate on what the signedness
 * changes: sign extension, the arithmetic right shift, negation and the signed comparisons.
 * uint128_t_gtest.cpp covers the shared operations.
 ***********************************************************************/
TEST(int128_t, ConstexprConstructionAndConversion)
{
    constexpr int128_t zero;
    constexpr int128_t negative(-7);
    constexpr int128_t positive(7);
    constexpr int128_t copied(negative);
    constexpr int128_t assigned = [] { int128_t a; a = -5; return a; }();
    constexpr int128_t from_char(static_cast<char>(-3));
    constexpr int128_t most_negative(0ull, 0x8000000000000000ull);

    static_assert(zero.is_zero());
    static_assert(static_cast<int64_t>(negative) == -7ll);
    static_assert(static_cast<int64_t>(assigned) == -5ll);
    static_assert(static_cast<int64_t>(from_char) == -3ll);
    static_assert(copied == negative);

    // A negative value is sign extended into the high QWORD, a non negative one is zero extended.
    // The values are constructed inside the lambdas: a captureless lambda cannot odr-use a local,
    // and passing an object to get_components() by reference is exactly that.
    constexpr uint64_t neg_high = [] {
        uint64_t l = 0, h = 0;
        int128_t(-7).get_components(l, h);
        return h;
    }();
    constexpr uint64_t pos_high = [] {
        uint64_t l = 0, h = 0;
        int128_t(7).get_components(l, h);
        return h;
    }();
    static_assert(neg_high == UINT64_MAX);
    static_assert(pos_high == 0ull);
    static_assert(int128_t(-1).get_bit(127) == 1);  // the sign bit reaches the top of the high QWORD
    static_assert(static_cast<uint64_t>(int128_t(1ull, 2ull)) == 1ull);

    static_assert(negative.is_negative() && !negative.is_positive());
    static_assert(positive.is_positive() && !positive.is_negative());
    static_assert(zero.is_positive());  // zero counts as positive
    static_assert(most_negative.is_negative());
    static_assert(int128_t(~0ull, ~0ull) == int128_t(-1));  // the all ones pattern is -1

    EXPECT_TRUE(int128_t(static_cast<int64_t>(opaque(-7ll))) == negative);
}
TEST(int128_t, ConstexprShiftsAndUnary)
{
    constexpr int128_t minus_one(-1);
    constexpr int128_t most_negative(0ull, 0x8000000000000000ull);

    // the right shift is arithmetic: it replicates the sign bit rather than shifting in zeros
    static_assert((minus_one >> 1) == minus_one);
    static_assert((minus_one >> 200) == minus_one);  // shifting past the width leaves all sign bits
    static_assert((int128_t(-16) >> 2) == int128_t(-4));
    static_assert((int128_t(-2) >> 1) == minus_one);
    static_assert((most_negative >> 127) == minus_one);
    static_assert((int128_t(0ull, 0x4000000000000000ull) >> 126) == int128_t(1));
    static_assert((int128_t(16) >> 2) == int128_t(4));  // a positive value shifts in zeros
    static_assert((int128_t(1) << 127) == most_negative);
    static_assert((int128_t(1) << 200).is_zero());

    static_assert(-int128_t(5) == int128_t(-5));
    static_assert(-int128_t(-5) == int128_t(5));
    static_assert(-int128_t(0) == int128_t(0));
    static_assert(+int128_t(-5) == int128_t(-5));
    static_assert(-most_negative == most_negative);  // its magnitude has no representation, like INT64_MIN
    static_assert(~int128_t(0) == minus_one);

    static_assert(abs(int128_t(-5)) == int128_t(5));
    static_assert(abs(int128_t(5)) == int128_t(5));
    static_assert(abs(int128_t(0)) == int128_t(0));

    EXPECT_TRUE((int128_t(static_cast<int64_t>(opaque(-1ll))) >> 200) == minus_one);
    EXPECT_TRUE((int128_t(static_cast<int64_t>(opaque(-16ll))) >> 2) == int128_t(-4));
    EXPECT_TRUE(abs(int128_t(static_cast<int64_t>(opaque(-5ll)))) == int128_t(5));
}
TEST(int128_t, ConstexprArithmeticAndComparisons)
{
    static_assert(static_cast<int64_t>(int128_t(3) + int128_t(4)) == 7ll);
    static_assert(static_cast<int64_t>(int128_t(3) - int128_t(4)) == -1ll);
    static_assert(static_cast<int64_t>(int128_t(-3) + int128_t(-4)) == -7ll);
    static_assert(int128_t(3) - int128_t(3) == int128_t(0));
    static_assert(static_cast<int64_t>(4 - int128_t(3)) == 1ll);  // scalar on the left

    static_assert(static_cast<int64_t>(int128_t(-6) * int128_t(7)) == -42ll);
    static_assert(static_cast<int64_t>(int128_t(-6) * int128_t(-7)) == 42ll);
    static_assert(static_cast<int64_t>(int128_t(-6) * -7) == 42ll);   // generic right hand side
    static_assert(static_cast<int64_t>(-7 * int128_t(-6)) == 42ll);   // scalar on the left
    static_assert(sqr(int128_t(-9)) == int128_t(81));
    static_assert(sqr(int128_t(-9)) == int128_t(-9) * int128_t(-9));  // sqr is bit identical to x * x
    static_assert(pow(int128_t(-3), 3u) == int128_t(-27));
    static_assert(pow(int128_t(-3), 2u) == int128_t(9));

    constexpr int128_t stepped = [] { int128_t a(5); ++a; a++; --a; return a; }();
    static_assert(static_cast<int64_t>(stepped) == 6ll);
    constexpr int128_t crossed_zero = [] { int128_t a(1); --a; --a; return a; }();
    static_assert(crossed_zero == int128_t(-1));

    // the comparisons are signed: a negative value is smaller whatever the magnitude of its bits
    static_assert(int128_t(-2) < int128_t(-1));
    static_assert(int128_t(-1) < int128_t(1));
    static_assert(int128_t(1) > int128_t(-1));
    static_assert(int128_t(0ull, 0x8000000000000000ull) < int128_t(0));  // the most negative value
    static_assert(int128_t(3ull, 0ull) < int128_t(1ull, 2ull));  // two positives, decided by the high QWORD
    static_assert(int128_t(-1) <= -1 && int128_t(-1) >= -1);
    static_assert(int128_t(-1) == -1 && int128_t(-1) != 1);

    static_assert(lzcnt128(int128_t(1)) == 127);
    static_assert(log2(int128_t(8)) == 3);

    EXPECT_TRUE(int128_t(static_cast<int64_t>(opaque(-6ll))) * int128_t(7) == int128_t(-42));
    EXPECT_TRUE(sqr(int128_t(static_cast<int64_t>(opaque(-9ll)))) == int128_t(81));
    EXPECT_TRUE(pow(int128_t(static_cast<int64_t>(opaque(-3ll))), 3u) == int128_t(-27));
    EXPECT_TRUE(int128_t(static_cast<int64_t>(opaque(-2ll))) < int128_t(-1));
}

// The double and float conversions used to be the one thing that could not happen at compile time:
// they read the inactive member of a union, which constant evaluation rejects. Now that Double and
// Float convert with std::bit_cast, a floating point literal crosses the boundary either way.
TEST(int128_t, ConstexprFloatingPointConversion)
{
    constexpr int128_t positive = 42.0;
    constexpr int128_t negative = -42.0;

    static_assert(static_cast<int64_t>(positive) == 42ll);
    static_assert(static_cast<int64_t>(negative) == -42ll);
    static_assert(static_cast<double>(positive) == 42.0);
    static_assert(static_cast<double>(negative) == -42.0);
    static_assert(static_cast<float>(negative) == -42.0f);
    static_assert(int128_t(0.0).is_zero());
    static_assert(negative.is_negative() && positive.is_positive());
    static_assert(int128_t(-3.99) == int128_t(-3));  // truncates towards zero
    static_assert(static_cast<double>(int128_t(-1e18)) == -1e18);

    EXPECT_TRUE(int128_t(static_cast<double>(opaque(-42ll))) == negative);
    EXPECT_DOUBLE_EQ(static_cast<double>(int128_t(opaque(-42ll))), -42.0);
}
