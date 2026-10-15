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
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <ctime>
#include <limits>
#include <ostream>
#include <string>
#include <type_traits>
#include "gtest_shared.h"

/**********************************************************************
 * Reference arithmetic
 *
 * Exact references for the operations under test, written against the raw 128 bit two's complement
 * pattern and sharing no code with fixed_point128. The 256 bit products come from
 * ReferenceMultiplyWide() in gtest_shared.h, a 32 bit limb schoolbook multiply; everything else is
 * plain 64 bit arithmetic below.
 *
 * They exist because a double reference sees 53 of the up to 127 significant bits of a result. A
 * multiply that truncated instead of rounding, a division that did the same, and a square root that
 * skipped its refinement all passed every test that compared against a double.
 ***********************************************************************/

/// @brief A raw 128 bit value, low QWORD first. Read as two's complement or as unsigned, as stated.
struct Raw128 {
    uint64_t low;
    uint64_t high;
};

/// @brief A 256 bit unsigned value, w[0] lowest.
struct U256 {
    uint64_t w[4];
};

/// @brief The raw bits of a fixed_point128.
template <int32_t I> [[nodiscard]] static Raw128 RawOf(const fixed_point128<I>& x)
{
    Raw128 r {};
    x.get_components(r.low, r.high);
    return r;
}

/// @brief A fixed_point128 made from raw bits.
template <int32_t I> [[nodiscard]] static fixed_point128<I> FromRaw(const Raw128& r)
{
    return fixed_point128<I>(r.low, r.high);
}

/// @brief The raw bits of a fixed_point128 as text, for failure messages.
template <int32_t I> [[nodiscard]] static std::string HexOf(const fixed_point128<I>& x)
{
    const Raw128 r = RawOf(x);
    char text[48];
    snprintf(text, sizeof(text), "0x%016llx_%016llx", static_cast<unsigned long long>(r.high), static_cast<unsigned long long>(r.low));
    return text;
}

/// @brief Two's complement negation of a raw value.
[[nodiscard]] static Raw128 RefNegate(Raw128 v)
{
    v.low = ~v.low + 1;
    v.high = ~v.high + ((v.low == 0) ? 1 : 0);
    return v;
}

/// @brief True when the raw value, read as two's complement, is negative.
[[nodiscard]] static bool RefIsNegative(const Raw128& v)
{
    return (v.high >> 63) != 0;
}

/// @brief The magnitude of a two's complement raw value, as unsigned: 2^127 for the most negative one.
[[nodiscard]] static Raw128 RefMagnitude(const Raw128& v)
{
    return RefIsNegative(v) ? RefNegate(v) : v;
}

/// @brief True when two raw values are the same bit pattern.
[[nodiscard]] static bool RefEqual(const Raw128& a, const Raw128& b)
{
    return a.low == b.low && a.high == b.high;
}

/// @brief True when a is below b, both read as unsigned.
[[nodiscard]] static bool RefLess(const Raw128& a, const Raw128& b)
{
    return (a.high != b.high) ? a.high < b.high : a.low < b.low;
}

/// @brief True when a is below b, both read as two's complement.
[[nodiscard]] static bool RefSignedLess(const Raw128& a, const Raw128& b)
{
    if (RefIsNegative(a) != RefIsNegative(b)) {
        return RefIsNegative(a);
    }
    return RefLess(a, b);
}

/// @brief a + b, modulo 2^128.
[[nodiscard]] static Raw128 RefAdd(const Raw128& a, const Raw128& b)
{
    Raw128 r {a.low + b.low, a.high + b.high};
    r.high += (r.low < a.low) ? 1 : 0;
    return r;
}

/// @brief a - b, modulo 2^128.
[[nodiscard]] static Raw128 RefSub(const Raw128& a, const Raw128& b)
{
    return RefAdd(a, RefNegate(b));
}

/// @brief An unsigned 128 bit value shifted right by n in [0, 127], vacated bits zero.
[[nodiscard]] static Raw128 RefShr(const Raw128& v, int32_t n)
{
    if (n == 0)
        return v;
    if (n < 64)
        return {(v.low >> n) | (v.high << (64 - n)), v.high >> n};
    return {v.high >> (n - 64), 0};
}

/// @brief An unsigned 128 bit value shifted left by n in [0, 127], bits shifted out lost.
[[nodiscard]] static Raw128 RefShl(const Raw128& v, int32_t n)
{
    if (n == 0)
        return v;
    if (n < 64)
        return {v.low << n, (v.high << n) | (v.low >> (64 - n))};
    return {0, v.low << (n - 64)};
}

/// @brief Bit n of an unsigned 128 bit value.
[[nodiscard]] static uint64_t RefBit(const Raw128& v, int32_t n)
{
    return (n < 64) ? (v.low >> n) & 1 : (v.high >> (n - 64)) & 1;
}

/// @brief Full 256 bit product of two unsigned 128 bit values.
[[nodiscard]] static U256 RefMul(const Raw128& a, const Raw128& b)
{
    uint128_t high;
    const uint128_t low = ReferenceMultiplyWide(uint128_t(a.low, a.high), uint128_t(b.low, b.high), high);
    U256 r {};
    low.get_components(r.w[0], r.w[1]);
    high.get_components(r.w[2], r.w[3]);
    return r;
}

/// @brief An unsigned 128 bit value times 2^n, n in [0, 255], as 256 bits; bits above 2^256 are lost.
[[nodiscard]] static U256 RefWiden(const Raw128& v, int32_t n)
{
    U256 r {};
    const int32_t words = n / 64, bits = n % 64;
    const uint64_t source[2] = {v.low, v.high};
    for (int32_t i = 0; i < 2; ++i) {
        if (i + words < 4) {
            r.w[i + words] |= (bits == 0) ? source[i] : (source[i] << bits);
        }
        if (bits != 0 && i + words + 1 < 4) {
            r.w[i + words + 1] |= source[i] >> (64 - bits);
        }
    }
    return r;
}

/// @brief A 256 bit value shifted right by n in [0, 255].
[[nodiscard]] static U256 RefShr256(const U256& v, int32_t n)
{
    U256 r {};
    const int32_t words = n / 64, bits = n % 64;
    for (int32_t i = 0; i + words < 4; ++i) {
        r.w[i] = (bits == 0) ? v.w[i + words] : (v.w[i + words] >> bits);
        if (bits != 0 && i + words + 1 < 4) {
            r.w[i] |= v.w[i + words + 1] << (64 - bits);
        }
    }
    return r;
}

/// @brief a + b on 256 bits, the carry out of the top dropped.
[[nodiscard]] static U256 RefAdd256(const U256& a, const U256& b)
{
    U256 r {};
    uint64_t carry = 0;
    for (int32_t i = 0; i < 4; ++i) {
        const uint64_t sum = a.w[i] + b.w[i];
        const uint64_t carry_out = (sum < a.w[i]) ? 1 : 0;
        r.w[i] = sum + carry;
        carry = carry_out | ((r.w[i] < sum) ? 1 : 0);
    }
    return r;
}

/// @brief Three way comparison of two 256 bit values.
[[nodiscard]] static int32_t RefCompare(const U256& a, const U256& b)
{
    for (int32_t i = 3; i >= 0; --i) {
        if (a.w[i] != b.w[i]) {
            return (a.w[i] < b.w[i]) ? -1 : 1;
        }
    }
    return 0;
}

/**
 * @brief Unsigned division n / d by restoring long division, one bit at a time. d must not be zero.
 * @param n Numerator
 * @param d Denominator, at most 2^127
 * @param quotient Receives the quotient
 * @param remainder Receives the remainder
 */
static void RefDivide(const Raw128& n, const Raw128& d, Raw128& quotient, Raw128& remainder)
{
    quotient = {0, 0};
    remainder = {0, 0};
    for (int32_t i = 127; i >= 0; --i) {
        remainder = {(remainder.low << 1) | RefBit(n, i), (remainder.high << 1) | (remainder.low >> 63)};
        if (!RefLess(remainder, d)) {
            remainder = RefSub(remainder, d);
            if (i >= 64) {
                quotient.high |= 1ull << (i - 64);
            } else {
                quotient.low |= 1ull << i;
            }
        }
    }
}

/**
 * @brief The product a * b on the grid of fixed_point128<I>, from the exact 256 bit product, rounded to
 *        nearest with a tie going away from zero (the magnitude rounded half up).
 * @param a Raw left operand
 * @param b Raw right operand
 * @param res Receives the raw product
 * @return false when the rounded product is 2^I or more in magnitude, where the operator wraps.
 */
template <int32_t I> [[nodiscard]] static bool RefProduct(const Raw128& a, const Raw128& b, Raw128& res)
{
    constexpr int32_t F = 127 - I;
    const U256 p = RefMul(RefMagnitude(a), RefMagnitude(b));
    U256 q = RefShr256(p, F);
    if ((p.w[(F - 1) / 64] >> ((F - 1) % 64)) & 1) {
        for (int32_t i = 0; i < 4; ++i) {
            if (++q.w[i] != 0)
                break;
        }
    }
    if (q.w[2] != 0 || q.w[3] != 0 || (q.w[1] >> 63) != 0)
        return false;

    res = {q.w[0], q.w[1]};
    if (RefIsNegative(a) != RefIsNegative(b)) {
        res = RefNegate(res);
    }
    return true;
}

/**
 * @brief Checks that q is n / d rounded to nearest with a tie going up, all three unsigned.
 *
 * The rounded quotient q satisfies (2q - 1) d <= 2n < (2q + 1) d, the lower bound dropping out for
 * q = 0, which needs no division to check: the numerator is passed doubled, so that the half step
 * stays an integer.
 *
 * @param twiceNumerator 2n, 256 bits
 * @param d Denominator, not zero
 * @param q Candidate quotient, below 2^127
 * @return True when q is the rounded quotient.
 */
[[nodiscard]] static bool RefIsRoundedQuotient(const U256& twiceNumerator, const Raw128& d, const Raw128& q)
{
    const Raw128 twiceQ = RefShl(q, 1);
    const Raw128 above = {twiceQ.low | 1, twiceQ.high};  // 2q + 1, bit 0 of 2q being clear
    if (RefCompare(twiceNumerator, RefMul(above, d)) >= 0)
        return false;
    if ((q.low | q.high) == 0)
        return true;
    const Raw128 below = RefSub(twiceQ, Raw128 {1, 0});  // 2q - 1
    return RefCompare(RefMul(below, d), twiceNumerator) <= 0;
}

/**
 * @brief Checks a fixed_point128<I> quotient a / b against the exact one: rounded to nearest, a tie away
 *        from zero, with the sign of a * b.
 * @return True when q is the correctly rounded quotient.
 */
template <int32_t I> [[nodiscard]] static bool RefIsQuotient(const Raw128& a, const Raw128& b, const Raw128& q)
{
    constexpr int32_t F = 127 - I;
    const Raw128 magnitude = RefMagnitude(q);
    const bool negative = RefIsNegative(a) != RefIsNegative(b);
    if ((magnitude.low | magnitude.high) != 0 && RefIsNegative(q) != negative)
        return false;
    // a / b scaled by 2^F is |a| * 2^F / |b|; the numerator goes in doubled
    return RefIsRoundedQuotient(RefWiden(RefMagnitude(a), F + 1), RefMagnitude(b), magnitude);
}

/// @brief True when a / b, exact, has a magnitude of 2^I or more and so does not fit the range.
template <int32_t I> [[nodiscard]] static bool RefQuotientOverflows(const Raw128& a, const Raw128& b)
{
    constexpr int32_t F = 127 - I;
    // |a| * 2^F >= 2^127 * |b|
    return RefCompare(RefWiden(RefMagnitude(a), F), RefWiden(RefMagnitude(b), 127)) >= 0;
}

/// @brief Checks a raw quotient a / k by an integer of magnitude k: rounded to nearest, a tie away from zero.
[[nodiscard]] static bool RefIsIntegerQuotient(const Raw128& a, uint64_t k, bool kNegative, const Raw128& q)
{
    const Raw128 magnitude = RefMagnitude(q);
    const bool negative = RefIsNegative(a) != kNegative;
    if ((magnitude.low | magnitude.high) != 0 && RefIsNegative(q) != negative)
        return false;
    return RefIsRoundedQuotient(RefWiden(RefMagnitude(a), 1), Raw128 {k, 0}, magnitude);
}

/**
 * @brief v / 2^s rounded to nearest with a tie going away from zero: the magnitude rounded half up.
 *
 * Rounding the magnitude half up needs nothing but the highest dropped bit, so this shares nothing with
 * the two's complement form the library computes.
 */
[[nodiscard]] static Raw128 RefShiftRightRounded(const Raw128& v, int32_t s)
{
    if (s == 0)
        return v;
    const Raw128 m = RefMagnitude(v);
    Raw128 q = (s < 128) ? RefShr(m, s) : Raw128 {0, 0};
    const uint64_t roundBit = (s <= 128) ? RefBit(m, s - 1) : 0;
    q = RefAdd(q, Raw128 {roundBit, 0});
    return RefIsNegative(v) ? RefNegate(q) : q;
}

/**
 * @brief The raw value of a double that lies on the grid of fixed_point128<I>, built from frexp() rather
 *        than through the constructor under test.
 */
template <int32_t I> [[nodiscard]] static Raw128 RefFromDouble(double d)
{
    constexpr int32_t F = 127 - I;
    if (d == 0)
        return {0, 0};
    int exponent = 0;
    const double mantissa = std::frexp(d, &exponent);
    const int64_t m = static_cast<int64_t>(std::ldexp(mantissa, 53));  // exact, below 2^53 in magnitude
    uint64_t magnitude = static_cast<uint64_t>((m < 0) ? -m : m);
    int32_t shift = exponent - 53 + F;  // d = m * 2^(exponent - 53)
    // trailing zero bits of the mantissa may lie below the grid; only set bits have to be on it
    while (shift < 0 && (magnitude & 1) == 0) {
        magnitude >>= 1;
        ++shift;
    }
    EXPECT_GE(shift, 0) << "the reference needs a double on the grid, d=" << d;
    const Raw128 raw = RefShl(Raw128 {magnitude, 0}, std::max(shift, 0));
    return (m < 0) ? RefNegate(raw) : raw;
}

/// @brief The raw value of the integer k on the grid of fixed_point128<I>, spelled without the type.
template <int32_t I> [[nodiscard]] static Raw128 RefFromInt(int64_t k)
{
    // k * 2^F: the integer bits start at bit F, which is bit 63 - I of the upper QWORD
    return Raw128 {0, static_cast<uint64_t>(k) << (63 - I)};
}

/// @brief The distance between two values in units of the last place.
template <int32_t I> [[nodiscard]] static double UlpDistance(const fixed_point128<I>& a, const fixed_point128<I>& b)
{
    const Raw128 d = RefMagnitude(RefSub(RawOf(a), RawOf(b)));
    return static_cast<double>(d.high) * 18446744073709551616.0 + static_cast<double>(d.low);
}

/// @brief A random raw value: a full random pattern with a random count of leading bits cleared, either sign.
[[nodiscard]] static Raw128 RandomRaw()
{
    const int32_t shift = 1 + static_cast<int32_t>(get_uint32_random() % 127);
    const Raw128 r = RefShr(Raw128 {get_uint64_random(), get_uint64_random()}, shift);
    return (rand() & 1) ? RefNegate(r) : r;
}

/// @brief A random signed integer of at most bits significant bits, sign included, bits in [1, 64].
[[nodiscard]] static int64_t RandomBoundedInt(int32_t bits)
{
    return get_int64_random() >> (64 - bits);
}

/// @brief A random unsigned integer below 2^bits, bits in [1, 64].
[[nodiscard]] static uint64_t RandomBoundedUInt(int32_t bits)
{
    return get_uint64_random() >> (64 - bits);
}

/**********************************************************************
 * fixed_point128 tests
 ***********************************************************************/
// Construct fixed_point128 and convert back to/from various elements.
TEST(fixed_point128, DefaultConstructor)
{
    fixed_point128<20> f;
    EXPECT_EQ(static_cast<uint64_t>(f), 0ull);
    EXPECT_TRUE(RefEqual(RawOf(f), Raw128 {0, 0}));
}
/**
 * @brief Every double on the grid of fixed_point128<I> converts exactly, and back again.
 * @param minExponent Smallest exponent whose 53 bit mantissa still ends on the grid, 52 - F.
 */
template <int32_t I> static void CheckConstructorFromDouble(int32_t minExponent)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const double value = get_double_random(minExponent, I);
        const fixed_point128<I> f = value;
        EXPECT_TRUE(RefEqual(RawOf(f), RefFromDouble<I>(value))) << "I=" << I << ", value=" << value << ", f=" << HexOf(f);
        EXPECT_EQ(static_cast<double>(f), value) << "I=" << I << ", value=" << value;
    }
}
// A guard used to replace any result more than 0.01% off with the expected value before asserting it.
// Every operand is exactly representable, so the conversion has to be exact, both ways and bit for bit.
TEST(fixed_point128, ConstructorFromDouble)
{
    CheckConstructorFromDouble<1>(-74);
    CheckConstructorFromDouble<32>(-43);
    CheckConstructorFromDouble<63>(-12);
}
TEST(fixed_point128, ConstructorFromFloat)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        float value = (float)get_double_random(-43, 31);
        fixed_point128<32> f = value;
        EXPECT_FLOAT_EQ(static_cast<float>(f), value);
    }
}
TEST(fixed_point128, ConstructorFromInt32)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        int32_t value = get_int32_random();
        fixed_point128<32> f = value;
        if (check_overflow(value, f)) {
            continue;
        }
        EXPECT_EQ(static_cast<int32_t>(f), value);
    }
}
TEST(fixed_point128, ConstructorFromUnsignedInt32)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        uint32_t value = get_uint32_random();
        fixed_point128<32> f = value;
        if (check_overflow(value, f)) {
            continue;
        }
        EXPECT_EQ(static_cast<uint32_t>(f), value);
    }
}
// The operands used to be full width 64 bit values, which fixed_point128<32> cannot hold: the guard
// rejected all 65536 iterations and the test asserted nothing. They are bounded to 31 bits now, and the
// guard has become an assertion, so a generator that drifts out of range fails loudly.
TEST(fixed_point128, ConstructorFromInt64)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const int64_t value = get_int64_random() >> (32 + get_uint32_random() % 32);
        const fixed_point128<32> f = value;
        ASSERT_FALSE(check_overflow(value, f)) << "value=" << value;
        EXPECT_EQ(static_cast<int64_t>(f), value);
        EXPECT_TRUE(RefEqual(RawOf(f), RefFromInt<32>(value))) << "value=" << value;
    }
}
TEST(fixed_point128, ConstructorFromUnsignedInt64)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const uint64_t value = get_uint64_random() >> (32 + get_uint32_random() % 32);
        const fixed_point128<32> f = value;
        ASSERT_FALSE(check_overflow(value, f)) << "value=" << value;
        EXPECT_EQ(static_cast<uint64_t>(f), value);
    }
}
// fixed_point128<63> holds every int64_t, INT64_MIN included, and a uint64_t above INT64_MAX wraps into
// the sign and back out of it unchanged.
TEST(fixed_point128, ConstructorFromInt64AtFullWidth)
{
    using fp = fixed_point128<63>;
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const int64_t value = get_int64_random();
        EXPECT_EQ(static_cast<int64_t>(fp(value)), value);
        EXPECT_TRUE(RefEqual(RawOf(fp(value)), RefFromInt<63>(value))) << "value=" << value;
        const uint64_t unsigned_value = get_uint64_random();
        EXPECT_EQ(static_cast<uint64_t>(fp(unsigned_value)), unsigned_value);
    }
    EXPECT_EQ(static_cast<int64_t>(fp(std::numeric_limits<int64_t>::min())), std::numeric_limits<int64_t>::min());
    EXPECT_EQ(static_cast<int64_t>(fp(std::numeric_limits<int64_t>::max())), std::numeric_limits<int64_t>::max());
    EXPECT_TRUE(fp(std::numeric_limits<int64_t>::min()) == std::numeric_limits<fp>::lowest());
}
TEST(fixed_point128, ConstructorFromString)
{
    const char* values[] = {"0", "12.34", "-3.1435", "12345.6789"};
    for (auto i = 0u; i < array_length(values); ++i) {
        fixed_point128<20> f = values[i];
        double d1 = strtod(values[i], nullptr);
        double d2 = strtod(static_cast<char*>(f), nullptr);
        EXPECT_DOUBLE_EQ(d1, d2) << "value=" << values[i];
    }
}
// The early returns of the string constructor: no string, nothing but white space, and an integer
// without a decimal point, whose sign takes a path of its own.
TEST(fixed_point128, ConstructorFromStringWithoutAFraction)
{
    using fp = fixed_point128<32>;
    EXPECT_TRUE(fp(static_cast<const char*>(nullptr)) == fp(0));
    EXPECT_TRUE(fp("") == fp(0));
    EXPECT_TRUE(fp(" \t ") == fp(0));
    EXPECT_TRUE(RefEqual(RawOf(fp("-5")), RefFromInt<32>(-5)));
    EXPECT_TRUE(RefEqual(RawOf(fp("+5")), RefFromInt<32>(5)));
    EXPECT_TRUE(RefEqual(RawOf(fp("7")), RefFromInt<32>(7)));
    EXPECT_TRUE(RefEqual(RawOf(fixed_point128<63>("-9223372036854775807")), RefFromInt<63>(-9223372036854775807ll)));
    EXPECT_TRUE(fp(std::string("-5")) == fp(-5));
}
/// @brief A decimal string and its exact value rounded to the grid, from Python fractions.
struct StringVector {
    int32_t i;
    const char* text;
    uint64_t expectedHigh, expectedLow;
};
static constexpr StringVector kStringVectors[] = {
    {1, "0.7101901153000597972121087318719189", 0x2D73C13DC9FF9B92ull, 0x2FC58C8431E02318ull},
    {1, "0.814754337808105276", 0x3424EF60CA3C9891ull, 0x1CA86E0D16D2D043ull},
    {1, "0.19721832051", 0x0C9F399730CCBC31ull, 0xD02BD9D208C0D7C3ull},
    {1, "0.4706402424270965759770792787441", 0x1E1EF84059FA202Aull, 0x84A4D0A1119BA1C3ull},
    {1, "0.5", 0x2000000000000000ull, 0x0000000000000000ull},
    {2, "1.18176435129425711630403", 0x25D103790C665090ull, 0x09F78D9EBE017CB7ull},
    {2, "1.7194163163783538629", 0x3705755DE1B8E2A1ull, 0x79EDA1DAE14084EEull},
    {2, "1.58278752209", 0x32A632047C96F569ull, 0x94EA4079BA7A9386ull},
    {2, "-0.815200123066536484659281446694816434", 0xE5E9E16E777D579Dull, 0x093659FF17E984AFull},
    {4, "-0.737264017646726724959", 0xFA1A15529D84975Dull, 0xFAF9CFA9C1330833ull},
    {4, "-3.897707105712955", 0xE0D17EEFDC9BBCB3ull, 0x570C3C97E4CDE59Dull},
    {4, "4.9594938418935855741899813452060", 0x27AD0B1B7D2BB6CCull, 0x4BB848A3D8F1F9C1ull},
    {4, "3.2585681447935810959327728991283", 0x1A118C2CED681623ull, 0xB69DF682863E4BA3ull},
    {16, "-25765.5289899705156099153935869700714652730", 0xCDAD3C4A0E805C9Aull, 0x0C276730BBC7CF6Aull},
    {16, "-24340.21033998282409759546653", 0xD075E51394565D59ull, 0xAC52441818A59D5Full},
    {16, "-573.59418903899582545427208633742596879977", 0xFEE133F19D12EF85ull, 0x7AF06ADD7967C4AFull},
    {16, "24446.18692177188243156485683977578", 0x2FBF17ED0D78929Cull, 0xFD3793F0D349B5BCull},
    {63, "-485700.132952643402", 0xFFFFFFFFFFF896BBull, 0xDDF6D0C8ABF1E5DDull},
    {63, "-104595.04655496926", 0xFFFFFFFFFFFE676Cull, 0xF414F9398FE031AAull},
    {63, "-9281.9945932397976858466", 0xFFFFFFFFFFFFDBBEull, 0x016256623F0A7A62ull},
    {63, "517524.34226203985836708611", 0x000000000007E594ull, 0x579E7C2BDA9B7DBDull},
};
/// @brief The bound the string constructor documents, plus the half ulp the rounded reference adds.
template <int32_t I> [[nodiscard]] static double StringConstructorBound()
{
    return 1.0 + 23.0 / std::ldexp(1.0, I - 1);
}
/// @brief Checks the string vectors of one I.
template <int32_t I> static void CheckStringVectors()
{
    int32_t checked = 0;
    for (const auto& v : kStringVectors) {
        if (v.i != I)
            continue;
        const fixed_point128<I> f(v.text);
        EXPECT_LE(UlpDistance(f, fixed_point128<I>(v.expectedLow, v.expectedHigh)), StringConstructorBound<I>())
            << "I=" << I << ", text=" << v.text << ", result=" << HexOf(f);
        ++checked;
    }
    EXPECT_GT(checked, 0) << "I=" << I;
}
// The string constructor holds each power of ten to 126 bits and sums the digits with no guard bits,
// so a narrow instantiation is a few ulp off: "0.5" itself is two ulp short at I = 1. The bound
// documented on the constructor is checked against exact values (Python fractions).
TEST(fixed_point128, ConstructorFromStringAccuracy)
{
    CheckStringVectors<1>();
    CheckStringVectors<2>();
    CheckStringVectors<4>();
    CheckStringVectors<16>();
    CheckStringVectors<63>();
}
/// @brief The transcendental constants per I, correctly rounded (mpmath 1.4.1), {high, low}; zero where
///        the constant needs more integer bits than I.
struct ConstantVector {
    int32_t i;
    uint64_t pi[2], e[2], halfPi[2], pi2[2], sqrt2[2], goldenRatio[2];
};
static constexpr ConstantVector kConstantVectors[] = {
    {1, {0ull, 0ull}, {0ull, 0ull}, {0x6487ED5110B4611Aull, 0x62633145C06E0E69ull}, {0ull, 0ull}, {0x5A827999FCEF3242ull, 0x2CBEC4D9BAA55F50ull}, {0x678DDE6E5FD29F05ull, 0x7CE73018173B720Dull}},
    {2, {0x6487ED5110B4611Aull, 0x62633145C06E0E69ull}, {0x56FC2A2C515DA54Dull, 0x57EE2B10139E9E79ull}, {0x3243F6A8885A308Dull, 0x313198A2E0370734ull}, {0ull, 0ull}, {0x2D413CCCFE779921ull, 0x165F626CDD52AFA8ull}, {0x33C6EF372FE94F82ull, 0xBE73980C0B9DB907ull}},
    {3, {0x3243F6A8885A308Dull, 0x313198A2E0370734ull}, {0x2B7E151628AED2A6ull, 0xABF7158809CF4F3Cull}, {0x1921FB54442D1846ull, 0x9898CC51701B839Aull}, {0x6487ED5110B4611Aull, 0x62633145C06E0E69ull}, {0x16A09E667F3BCC90ull, 0x8B2FB1366EA957D4ull}, {0x19E3779B97F4A7C1ull, 0x5F39CC0605CEDC83ull}},
    {4, {0x1921FB54442D1846ull, 0x9898CC51701B839Aull}, {0x15BF0A8B14576953ull, 0x55FB8AC404E7A79Eull}, {0x0C90FDAA22168C23ull, 0x4C4C6628B80DC1CDull}, {0x3243F6A8885A308Dull, 0x313198A2E0370734ull}, {0x0B504F333F9DE648ull, 0x4597D89B3754ABEAull}, {0x0CF1BBCDCBFA53E0ull, 0xAF9CE60302E76E42ull}},
    {16, {0x0001921FB54442D1ull, 0x8469898CC51701B8ull}, {0x00015BF0A8B14576ull, 0x95355FB8AC404E7Aull}, {0x0000C90FDAA22168ull, 0xC234C4C6628B80DCull}, {0x0003243F6A8885A3ull, 0x08D313198A2E0370ull}, {0x0000B504F333F9DEull, 0x6484597D89B3754Bull}, {0x0000CF1BBCDCBFA5ull, 0x3E0AF9CE60302E77ull}},
    {32, {0x00000001921FB544ull, 0x42D18469898CC517ull}, {0x000000015BF0A8B1ull, 0x457695355FB8AC40ull}, {0x00000000C90FDAA2ull, 0x2168C234C4C6628Cull}, {0x00000003243F6A88ull, 0x85A308D313198A2Eull}, {0x00000000B504F333ull, 0xF9DE6484597D89B3ull}, {0x00000000CF1BBCDCull, 0xBFA53E0AF9CE6030ull}},
    {63, {0x0000000000000003ull, 0x243F6A8885A308D3ull}, {0x0000000000000002ull, 0xB7E151628AED2A6Bull}, {0x0000000000000001ull, 0x921FB54442D1846Aull}, {0x0000000000000006ull, 0x487ED5110B4611A6ull}, {0x0000000000000001ull, 0x6A09E667F3BCC909ull}, {0x0000000000000001ull, 0x9E3779B97F4A7C16ull}},
};
/// @brief Checks the constants of one I against the table, within the string constructor's bound.
template <int32_t I> static void CheckConstants()
{
    using fp = fixed_point128<I>;
    const ConstantVector* v = nullptr;
    for (const auto& c : kConstantVectors) {
        if (c.i == I)
            v = &c;
    }
    ASSERT_NE(v, nullptr) << "I=" << I;
    const double bound = StringConstructorBound<I>();
    if constexpr (I >= 2) {
        EXPECT_LE(UlpDistance(fp::pi(), fp(v->pi[1], v->pi[0])), bound) << "pi, I=" << I;
        EXPECT_LE(UlpDistance(fp::e(), fp(v->e[1], v->e[0])), bound) << "e, I=" << I;
    }
    if constexpr (I >= 3) {
        EXPECT_LE(UlpDistance(fp::pi2(), fp(v->pi2[1], v->pi2[0])), bound) << "pi2, I=" << I;
    }
    EXPECT_LE(UlpDistance(fp::half_pi(), fp(v->halfPi[1], v->halfPi[0])), bound) << "half_pi, I=" << I;
    EXPECT_LE(UlpDistance(fp::sqrt_2(), fp(v->sqrt2[1], v->sqrt2[0])), bound) << "sqrt_2, I=" << I;
    EXPECT_LE(UlpDistance(fp::golden_ratio(), fp(v->goldenRatio[1], v->goldenRatio[0])), bound) << "golden_ratio, I=" << I;
}
// The constants are parsed from decimal strings and carry the string constructor's bound; at I = 1
// sqrt_2() is 17 ulp off. Checked against mpmath at every width the table covers.
TEST(fixed_point128, ConstantsAccuracy)
{
    CheckConstants<1>();
    CheckConstants<2>();
    CheckConstants<3>();
    CheckConstants<4>();
    CheckConstants<16>();
    CheckConstants<32>();
    CheckConstants<63>();
}
/// @brief Checks fact_reciprocal(k) against 2^F / k!, divided exactly here; k! fits 128 bits up to 34.
template <int32_t I> static void CheckFactorialReciprocal()
{
    using fp = fixed_point128<I>;
    constexpr int32_t F = 127 - I;
    Raw128 factorial {1, 0};
    for (int32_t k = 0; k <= 34; ++k) {
        if (k > 1) {
            // factorial *= k, which stays below 2^128 up to 34!
            uint128_t high;
            const uint128_t low = ReferenceMultiplyWide(uint128_t(factorial.low, factorial.high), uint128_t(static_cast<uint64_t>(k), 0ull), high);
            low.get_components(factorial.low, factorial.high);
        }
        Raw128 quotient {}, remainder {};
        RefDivide(RefShl(Raw128 {1, 0}, F), factorial, quotient, remainder);
        // rounded to nearest: up when twice the remainder reaches the divisor
        if (!RefLess(RefShl(remainder, 1), factorial)) {
            quotient = RefAdd(quotient, Raw128 {1, 0});
        }
        fp result;
        fact_reciprocal(k, result);
        EXPECT_LE(UlpDistance(result, FromRaw<I>(quotient)), StringConstructorBound<I>()) << "I=" << I << ", k=" << k;
    }
    for (int32_t k : {-1, 35, 40}) {
        fp result = fp::one();
        fact_reciprocal(k, result);
        EXPECT_TRUE(result == fp(0)) << "I=" << I << ", k=" << k;
    }
}
// fact_reciprocal() has no test of its own and feeds the exponential series; its table is parsed
// from strings, so it carries the string constructor's bound.
TEST(fixed_point128, FactorialReciprocalAccuracy)
{
    CheckFactorialReciprocal<1>();
    CheckFactorialReciprocal<2>();
    CheckFactorialReciprocal<4>();
    CheckFactorialReciprocal<16>();
    CheckFactorialReciprocal<32>();
    CheckFactorialReciprocal<63>();
}
TEST(fixed_point128, CopyConstructor)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const Raw128 raw = RandomRaw();
        const fixed_point128<20> f1 = FromRaw<20>(raw);
        const fixed_point128<20> f2 = f1;
        EXPECT_TRUE(RefEqual(RawOf(f2), raw));
    }
}
TEST(fixed_point128, MoveConstructor)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const Raw128 raw = RandomRaw();
        fixed_point128<20> f1 = FromRaw<20>(raw);
        fixed_point128<20> f2(std::move(f1));
        EXPECT_TRUE(RefEqual(RawOf(f2), raw));
        // documented to leave the source unchanged, it acts as a copy
        EXPECT_TRUE(RefEqual(RawOf(f1), raw));  // NOLINT(bugprone-use-after-move)
    }
}
TEST(fixed_point128, AssignmentOperator)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const Raw128 raw = RandomRaw();
        const fixed_point128<20> f1 = FromRaw<20>(raw);
        fixed_point128<20> f2;
        f2 = f1;
        EXPECT_TRUE(RefEqual(RawOf(f2), raw));
    }
}
TEST(fixed_point128, AssignmentOperatorOtherType)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value = get_double_random(-50, 18);
        fixed_point128<20> f1 = value;
        fixed_point128<22> f2;
        f2 = f1;
        // two more integer bits drop the two lowest fraction bits; the operand has none set there
        EXPECT_TRUE(RefEqual(RawOf(f2), RefFromDouble<22>(value))) << "value=" << value;
        f1 = f2;
        EXPECT_TRUE(RefEqual(RawOf(f1), RefFromDouble<20>(value))) << "value=" << value;
    }
}
TEST(fixed_point128, MoveAssignmentOperator)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const Raw128 raw = RandomRaw();
        fixed_point128<20> f1 = FromRaw<20>(raw);
        fixed_point128<20> f2;
        f2 = std::move(f1);
        EXPECT_TRUE(RefEqual(RawOf(f2), raw));
    }
}
// The value used to be picked once, by a loop whose flag was never reset, so all 65536 iterations
// converted the same number.
TEST(fixed_point128, CopyConstructorOtherType)
{
    constexpr int f1_prec = 63, f2_prec = 32;
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        // below 2^31, inside both ranges, and on the grid of both
        const double value = get_double_random(-12, 31);
        const fixed_point128<f1_prec> f1 = value;
        const fixed_point128<f2_prec> f2 = f1;  // more fraction bits: exact
        EXPECT_TRUE(RefEqual(RawOf(f2), RefFromDouble<f2_prec>(value))) << "value=" << value;
        const fixed_point128<f1_prec> f3 = f2;  // back: the dropped fraction bits are all zero
        EXPECT_TRUE(f3 == f1) << "value=" << value;
    }
}
/// @brief Checks the conversion from fixed_point128<From> to fixed_point128<To> against the reference.
template <int32_t From, int32_t To> static void CheckCrossWidthConversion()
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT / 4; ++i) {
        const Raw128 raw = RandomRaw();
        const fixed_point128<From> x = FromRaw<From>(raw);
        const fixed_point128<To> converted(x);
        fixed_point128<To> assigned;
        assigned = x;
        Raw128 expected {};
        if constexpr (To > From) {
            // fewer fraction bits: a right shift, rounded to nearest with a tie away from zero
            expected = RefShiftRightRounded(raw, To - From);
        } else {
            // more fraction bits: a left shift, exact as long as the value fits, and wrapping otherwise
            expected = RefShl(raw, From - To);
        }
        EXPECT_TRUE(RefEqual(RawOf(converted), expected)) << From << "->" << To << ", x=" << HexOf(x) << ", got " << HexOf(converted);
        EXPECT_TRUE(converted == assigned) << From << "->" << To << ", x=" << HexOf(x);
    }
}
// The widest width difference, 62 bits each way, at compile time and at run time. The narrowing
// direction drops up to 62 fraction bits and has to round them to nearest: 0.5 + 2^-62 of the new last
// place used to be rounded down, the sticky bit covering one dropped bit instead of all of them.
TEST(fixed_point128, CrossWidthConversionAtTheExtremes)
{
    constexpr fixed_point128<63> widened(fixed_point128<1>(0x8000000000000001ull, 0x2000000000000000ull));
    static_assert(widened == fixed_point128<63>(0x8000000000000002ull, 0));
    constexpr fixed_point128<1> narrowed(fixed_point128<63>(1));
    static_assert(narrowed == fixed_point128<1>::one());

    // raw 2^125 + 2^61 + 1 at I = 1 is 2^63 + 1/2 + 2^-62 at I = 63, which rounds up
    const fixed_point128<1> src(0x2000000000000001ull, 0x2000000000000000ull);
    EXPECT_TRUE(RefEqual(RawOf(fixed_point128<63>(src)), Raw128 {0x8000000000000001ull, 0}));
    // an exact tie goes away from zero: 2^63 + 1/2 and its negation
    const fixed_point128<1> tie(0x2000000000000000ull, 0x2000000000000000ull);
    EXPECT_TRUE(RefEqual(RawOf(fixed_point128<63>(tie)), Raw128 {0x8000000000000001ull, 0}));
    EXPECT_TRUE(RefEqual(RawOf(fixed_point128<63>(-tie)), RefNegate(Raw128 {0x8000000000000001ull, 0})));
    // a value too large for the narrower type wraps: 3.25 is -0.75 in fixed_point128<1>
    EXPECT_TRUE(fixed_point128<1>(fixed_point128<63>(3.25)) == fixed_point128<1>(-0.75));

    CheckCrossWidthConversion<1, 63>();
    CheckCrossWidthConversion<63, 1>();
    CheckCrossWidthConversion<20, 22>();
    CheckCrossWidthConversion<32, 16>();
}

/**********************************************************************
 * Addition and subtraction
 *
 * The operands are doubles on the grid of fixed_point128<40> and their sum or difference is checked
 * bit for bit against the sum of the raw values, built without the type. The double generators are
 * bounded so that every operand and every result fits: with the default exponent range a third to a
 * half of the iterations used to be skipped.
 ***********************************************************************/
TEST(fixed_point128, AddSameSign)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const double value1 = fabs(get_double_random(-10, 38));
        const double value2 = value1 * 2.5;
        const fixed_point128<40> f1 = value1;
        const fixed_point128<40> f2 = value2;
        const fixed_point128<40> f3 = f1 + f2;
        ASSERT_FALSE(check_overflow(value1 + value2, f3)) << "value1=" << value1;
        EXPECT_TRUE(RefEqual(RawOf(f3), RefAdd(RefFromDouble<40>(value1), RefFromDouble<40>(value2)))) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(fixed_point128, AddDifferentSign)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const double value1 = fabs(get_double_random(-10, 38));
        const double value2 = value1 * -2.5;
        const fixed_point128<40> f1 = value1;
        const fixed_point128<40> f2 = value2;
        const fixed_point128<40> f3 = f1 + f2;
        EXPECT_TRUE(RefEqual(RawOf(f3), RefAdd(RefFromDouble<40>(value1), RefFromDouble<40>(value2)))) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(fixed_point128, AddDouble)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const double value1 = get_double_random(-10, 39);
        const double value2 = get_double_random(-10, 39);
        const fixed_point128<40> f1 = value1;
        const fixed_point128<40> f3 = f1 + value2;
        EXPECT_TRUE(RefEqual(RawOf(f3), RefAdd(RefFromDouble<40>(value1), RefFromDouble<40>(value2)))) << "value1=" << value1 << ", value2=" << value2;
        // the scalar on the left takes the same path
        EXPECT_TRUE(value2 + f1 == f3) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(fixed_point128, AddFloat)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const float value1 = (float)get_double_random(-10, 39);
        const float value2 = (float)get_double_random(-10, 39);
        const fixed_point128<40> f1 = value1;
        const fixed_point128<40> f3 = f1 + value2;
        EXPECT_TRUE(RefEqual(RawOf(f3), RefAdd(RefFromDouble<40>(value1), RefFromDouble<40>(value2)))) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(fixed_point128, AddInt32)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const int32_t value1 = static_cast<int32_t>(RandomBoundedInt(32));
        const int32_t value2 = static_cast<int32_t>(RandomBoundedInt(32));
        // computed in a wider type, the sum of two 32 bit values can exceed 32 bits
        const int64_t res = static_cast<int64_t>(value1) + value2;
        const fixed_point128<40> f1 = value1;
        const fixed_point128<40> f3 = f1 + value2;
        EXPECT_EQ(static_cast<int64_t>(f3), res) << "value1=" << value1 << ", value2=" << value2;
        EXPECT_TRUE(RefEqual(RawOf(f3), RefFromInt<40>(res))) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(fixed_point128, AddUnsignedInt32)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        // 31 bit operands, whose sum cannot wrap the uint32_t reference (half used to be skipped)
        const uint32_t value1 = get_uint32_random() >> 1;
        const uint32_t value2 = get_uint32_random() >> 1;
        const uint32_t res = value1 + value2;
        const fixed_point128<40> f1 = value1;
        const fixed_point128<40> f3 = f1 + value2;
        EXPECT_EQ(static_cast<uint64_t>(f3), res) << "value1=" << value1 << ", value2=" << value2;
    }
}
// The 64 bit operands used to be full width, beyond the range of fixed_point128<40>: every iteration was
// skipped, and the reference sum overflowed int64_t, which is undefined behavior. Both are bounded to
// 39 bits now, so the sum fits the type and the reference.
TEST(fixed_point128, AddInt64)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const int64_t value1 = get_int64_random() >> (25 + get_uint32_random() % 39);
        const int64_t value2 = get_int64_random() >> (25 + get_uint32_random() % 39);
        const int64_t res = value1 + value2;
        const fixed_point128<40> f1 = value1;
        const fixed_point128<40> f3 = f1 + value2;
        ASSERT_FALSE(check_overflow(res, f3)) << "value1=" << value1 << ", value2=" << value2;
        EXPECT_EQ(static_cast<int64_t>(f3), res) << "value1=" << value1 << ", value2=" << value2;
        EXPECT_TRUE(RefEqual(RawOf(f3), RefFromInt<40>(res))) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(fixed_point128, AddUnsignedInt64)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const uint64_t value1 = get_uint64_random() >> (25 + get_uint32_random() % 39);
        const uint64_t value2 = get_uint64_random() >> (25 + get_uint32_random() % 39);
        const uint64_t res = value1 + value2;
        const fixed_point128<40> f1 = value1;
        const fixed_point128<40> f3 = f1 + value2;
        ASSERT_FALSE(check_overflow(res, f3)) << "value1=" << value1 << ", value2=" << value2;
        EXPECT_EQ(static_cast<uint64_t>(f3), res) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(fixed_point128, SubtractSameSign)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const double value1 = fabs(get_double_random(-10, 38));
        const double value2 = value1 * 2.5;
        const fixed_point128<40> f1 = value1;
        const fixed_point128<40> f2 = value2;
        const fixed_point128<40> f3 = f1 - f2;
        EXPECT_TRUE(RefEqual(RawOf(f3), RefSub(RefFromDouble<40>(value1), RefFromDouble<40>(value2)))) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(fixed_point128, SubtractDifferentSign)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const double value1 = fabs(get_double_random(-10, 38));
        const double value2 = value1 * -2.5;
        const fixed_point128<40> f1 = value1;
        const fixed_point128<40> f2 = value2;
        const fixed_point128<40> f3 = f1 - f2;
        EXPECT_TRUE(RefEqual(RawOf(f3), RefSub(RefFromDouble<40>(value1), RefFromDouble<40>(value2)))) << "value1=" << value1 << ", value2=" << value2;
        EXPECT_TRUE(value1 - f2 == f3) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(fixed_point128, SubtractInt32)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const int32_t value1 = static_cast<int32_t>(RandomBoundedInt(32));
        const int32_t value2 = static_cast<int32_t>(RandomBoundedInt(32));
        const int64_t res = static_cast<int64_t>(value1) - value2;
        const fixed_point128<40> f1 = value1;
        const fixed_point128<40> f3 = f1 - value2;
        EXPECT_EQ(static_cast<int64_t>(f3), res) << "value1=" << value1 << ", value2=" << value2;
        EXPECT_TRUE(RefEqual(RawOf(f3), RefFromInt<40>(res))) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(fixed_point128, SubtractUnsignedInt32)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        auto value1 = get_uint32_random();
        auto value2 = get_uint32_random();
        if (value1 < value2) {
            std::swap(value1, value2);
        }
        const uint32_t res = value1 - value2;
        const fixed_point128<40> f1 = value1;
        const fixed_point128<40> f3 = f1 - value2;
        EXPECT_EQ(static_cast<uint64_t>(f3), res) << "value1=" << value1 << ", value2=" << value2;
    }
}
// See AddInt64: the operands are bounded to 39 bits so the difference fits and does not overflow the
// reference.
TEST(fixed_point128, SubtractInt64)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const int64_t value1 = get_int64_random() >> (25 + get_uint32_random() % 39);
        const int64_t value2 = get_int64_random() >> (25 + get_uint32_random() % 39);
        const int64_t res = value1 - value2;
        const fixed_point128<40> f1 = value1;
        const fixed_point128<40> f3 = f1 - value2;
        ASSERT_FALSE(check_overflow(res, f3)) << "value1=" << value1 << ", value2=" << value2;
        EXPECT_EQ(static_cast<int64_t>(f3), res) << "value1=" << value1 << ", value2=" << value2;
        EXPECT_TRUE(RefEqual(RawOf(f3), RefFromInt<40>(res))) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(fixed_point128, SubtractUnsignedInt64)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        uint64_t value1 = get_uint64_random() >> (25 + get_uint32_random() % 39);
        uint64_t value2 = get_uint64_random() >> (25 + get_uint32_random() % 39);
        if (value1 < value2) {
            std::swap(value1, value2);
        }
        const uint64_t res = value1 - value2;
        const fixed_point128<40> f1 = value1;
        const fixed_point128<40> f3 = f1 - value2;
        ASSERT_FALSE(check_overflow(res, f3)) << "value1=" << value1 << ", value2=" << value2;
        EXPECT_EQ(static_cast<uint64_t>(f3), res) << "value1=" << value1 << ", value2=" << value2;
    }
}

/**********************************************************************
 * Multiplication
 *
 * Checked against the exact 256 bit product, rounded to nearest with a tie away from zero. The
 * double reference these tests used before could not see the last 30 or more bits of a product: a
 * multiply that truncated instead of rounding passed every one of them.
 ***********************************************************************/
TEST(fixed_point128, MultiplyByFP128)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        // the product of two operands below 2^20 and 2^19 stays below 2^39
        const double value1 = get_double_random(-10, 20);
        const double value2 = get_double_random(-10, 19);
        const fixed_point128<40> f1 = value1;
        const fixed_point128<40> f2 = value2;
        Raw128 expected {};
        ASSERT_TRUE(RefProduct<40>(RefFromDouble<40>(value1), RefFromDouble<40>(value2), expected)) << "value1=" << value1;
        EXPECT_TRUE(RefEqual(RawOf(f1 * f2), expected)) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(fixed_point128, MultiplyByDouble)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const double value1 = get_double_random(-10, 20);
        const double value2 = get_double_random(-10, 19);
        const fixed_point128<40> f1 = value1;
        Raw128 expected {};
        ASSERT_TRUE(RefProduct<40>(RefFromDouble<40>(value1), RefFromDouble<40>(value2), expected)) << "value1=" << value1;
        EXPECT_TRUE(RefEqual(RawOf(f1 * value2), expected)) << "value1=" << value1 << ", value2=" << value2;
        EXPECT_TRUE(RefEqual(RawOf(value2 * f1), expected)) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(fixed_point128, MultiplyByFloat)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const float value1 = (float)get_double_random(-10, 20);
        const float value2 = (float)get_double_random(-10, 19);
        const fixed_point128<40> f1 = value1;
        Raw128 expected {};
        ASSERT_TRUE(RefProduct<40>(RefFromDouble<40>(value1), RefFromDouble<40>(value2), expected)) << "value1=" << value1;
        EXPECT_TRUE(RefEqual(RawOf(f1 * value2), expected)) << "value1=" << value1 << ", value2=" << value2;
    }
}
/**
 * @brief Multiplies by an integer scalar of type T, the two operands sharing a budget of bits so that
 *        the exact product always fits fixed_point128<I>.
 *
 * The four tests below used two full width operands against fixed_point128<40>: the product left the
 * range on practically every iteration and the guard skipped all of them, while the int64_t reference
 * product overflowed. Splitting a bit budget between the operands keeps every product in range and
 * still reaches the full width of T for one operand when the other is small.
 *
 * @tparam I Integer bits of the fixed_point128
 * @tparam T The scalar type
 * @param budget Combined significant bits of the two operands, at most I - 1
 */
template <int32_t I, typename T> static void CheckMultiplyByScalar(int32_t budget)
{
    constexpr int32_t scalarBits = static_cast<int32_t>(sizeof(T) * 8);
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const int32_t bits1 = 1 + static_cast<int32_t>(get_uint32_random() % static_cast<uint32_t>(std::min(budget - 1, scalarBits)));
        const int32_t bits2 = 1 + static_cast<int32_t>(get_uint32_random() % static_cast<uint32_t>(std::min(budget - bits1, scalarBits)));
        T value1, value2;
        if constexpr (std::is_signed_v<T>) {
            value1 = static_cast<T>(RandomBoundedInt(bits1));
            value2 = static_cast<T>(RandomBoundedInt(bits2));
        } else {
            value1 = static_cast<T>(RandomBoundedUInt(bits1));
            value2 = static_cast<T>(RandomBoundedUInt(bits2));
        }
        const fixed_point128<I> f1 = value1;
        const fixed_point128<I> f3 = f1 * value2;
        if constexpr (std::is_signed_v<T>) {
            // the operands are below 2^(bits - 1) in magnitude, so the product fits int64_t
            const int64_t res = static_cast<int64_t>(value1) * static_cast<int64_t>(value2);
            ASSERT_FALSE(check_overflow(res, f3)) << "I=" << I << ", value1=" << value1 << ", value2=" << value2;
            EXPECT_EQ(static_cast<int64_t>(f3), res) << "I=" << I << ", value1=" << value1 << ", value2=" << value2;
            EXPECT_TRUE(RefEqual(RawOf(f3), RefFromInt<I>(res))) << "I=" << I << ", value1=" << value1 << ", value2=" << value2;
        } else {
            const uint64_t res = static_cast<uint64_t>(value1) * static_cast<uint64_t>(value2);
            ASSERT_FALSE(check_overflow(res, f3)) << "I=" << I << ", value1=" << value1 << ", value2=" << value2;
            EXPECT_EQ(static_cast<uint64_t>(f3), res) << "I=" << I << ", value1=" << value1 << ", value2=" << value2;
        }
    }
}
TEST(fixed_point128, MultiplyByInt32)
{
    CheckMultiplyByScalar<40, int32_t>(40);
}
TEST(fixed_point128, MultiplyByUnsignedInt32)
{
    CheckMultiplyByScalar<40, uint32_t>(39);
}
TEST(fixed_point128, MultiplyByInt64)
{
    CheckMultiplyByScalar<40, int64_t>(40);
    // fixed_point128<63> holds a product of up to 63 bits, which takes one operand to the full width
    CheckMultiplyByScalar<63, int64_t>(63);
}
TEST(fixed_point128, MultiplyByUnsignedInt64)
{
    CheckMultiplyByScalar<40, uint64_t>(39);
    CheckMultiplyByScalar<63, uint64_t>(62);
}
/// @brief Checks the multiply of random raw operands against the exact product, for one I.
template <int32_t I> static void CheckMultiplyAgainstReference()
{
    using fp = fixed_point128<I>;
    srand(RANDOM_SEED);
    int32_t checked = 0;
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const Raw128 a = RandomRaw(), b = RandomRaw();
        Raw128 expected {};
        if (!RefProduct<I>(a, b, expected))
            continue;  // beyond the range, where the product wraps by design
        ++checked;
        const fp x = FromRaw<I>(a), y = FromRaw<I>(b);
        EXPECT_TRUE(RefEqual(RawOf(x * y), expected)) << "I=" << I << ", x=" << HexOf(x) << ", y=" << HexOf(y) << ", got " << HexOf(x * y);
        if (RefEqual(a, b)) {
            EXPECT_TRUE(sqr(x) == x * y) << "I=" << I << ", x=" << HexOf(x);
        }
    }
    // the generator keeps most products in range; a tripwire against a loop that asserts nothing
    EXPECT_GT(checked, static_cast<int32_t>(RANDOM_TEST_COUNT) / 2) << "I=" << I;

    // exact ties: half an epsilon goes away from zero, 1.5 epsilon to two of them
    const fp eps = fp::epsilon(), half = fp::half();
    EXPECT_TRUE(eps * half == eps) << "I=" << I;
    EXPECT_TRUE(-eps * half == -eps) << "I=" << I;
    EXPECT_TRUE((eps * 3) * half == eps * 2) << "I=" << I;
    EXPECT_TRUE(-(eps * 3) * half == -(eps * 2)) << "I=" << I;
    // aliasing: an object multiplied by itself in place
    fp z = fp::one() + fp::half();
    z *= z;
    EXPECT_TRUE(z == fp::one() + fp::one() + (fp::one() >> 2)) << "I=" << I;  // 2.25 needs I >= 2
}
// Every product of random raw operands that stays in range, bit for bit against the exact 256 bit
// product rounded half away from zero, at both ends of the width range and in between: I = 63 reads the
// product from the lowest QWORDs, every other width from a shifted window.
TEST(fixed_point128, MultiplyMatchesExactReference)
{
    CheckMultiplyAgainstReference<2>();
    CheckMultiplyAgainstReference<32>();
    CheckMultiplyAgainstReference<62>();
    CheckMultiplyAgainstReference<63>();

    // fixed_point128<1> cannot hold the 2.25 of the aliasing check above
    using fp1 = fixed_point128<1>;
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const Raw128 a = RandomRaw(), b = RandomRaw();
        Raw128 expected {};
        if (!RefProduct<1>(a, b, expected))
            continue;
        EXPECT_TRUE(RefEqual(RawOf(FromRaw<1>(a) * FromRaw<1>(b)), expected)) << "x=" << HexOf(FromRaw<1>(a)) << ", y=" << HexOf(FromRaw<1>(b));
    }
    EXPECT_TRUE(fp1::epsilon() * fp1::half() == fp1::epsilon());
}

/**********************************************************************
 * Division
 *
 * A quotient is checked without dividing: q is a / b rounded to nearest, a tie away from zero, when
 * (2|q| - 1) |b| <= 2 |a| 2^F < (2|q| + 1) |b|, two products of the exact multiply. The double
 * comparison these tests made before accepted an error of 2^-48 relative, which a division that
 * truncated instead of rounding met with room to spare.
 ***********************************************************************/
TEST(fixed_point128, DivideByFP128)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        // operands from 2^-10 to 2^29 keep the quotient below 2^39
        const double value1 = get_double_random(-10, 29);
        const double value2 = get_double_random(-10, 29);
        const fixed_point128<40> f1 = value1;
        const fixed_point128<40> f2 = value2;
        const fixed_point128<40> f3 = f1 / f2;
        EXPECT_TRUE(RefIsQuotient<40>(RefFromDouble<40>(value1), RefFromDouble<40>(value2), RawOf(f3))) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(fixed_point128, DivideByDouble)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const double value1 = get_double_random(-10, 29);
        const double value2 = get_double_random(-10, 29);
        const fixed_point128<40> f1 = value1;
        const fixed_point128<40> f3 = f1 / value2;
        EXPECT_TRUE(RefIsQuotient<40>(RefFromDouble<40>(value1), RefFromDouble<40>(value2), RawOf(f3))) << "value1=" << value1 << ", value2=" << value2;
        EXPECT_TRUE(value1 / fixed_point128<40>(value2) == f3) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(fixed_point128, DivideByFloat)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const float value1 = (float)get_double_random(-10, 29);
        const float value2 = (float)get_double_random(-10, 29);
        const fixed_point128<40> f1 = value1;
        const fixed_point128<40> f3 = f1 / value2;
        EXPECT_TRUE(RefIsQuotient<40>(RefFromDouble<40>(value1), RefFromDouble<40>(value2), RawOf(f3))) << "value1=" << value1 << ", value2=" << value2;
    }
}
/// @brief Divides by integer scalars of type T and checks the rounded quotient against the reference.
template <typename T> static void CheckDivideByScalar()
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const double value1 = get_double_random(-1, 39);
        T value2 {};
        if constexpr (std::is_signed_v<T>) {
            value2 = static_cast<T>(RandomBoundedInt(static_cast<int32_t>(sizeof(T) * 8)));
        } else {
            value2 = static_cast<T>(RandomBoundedUInt(static_cast<int32_t>(sizeof(T) * 8)));
        }
        if (value2 == 0)
            continue;
        const fixed_point128<40> f1 = value1;
        const fixed_point128<40> f3 = f1 / value2;
        bool negative = false;
        uint64_t magnitude = static_cast<uint64_t>(value2);
        if constexpr (std::is_signed_v<T>) {
            negative = value2 < 0;
            if (negative)
                magnitude = 0ull - magnitude;
        }
        EXPECT_TRUE(RefIsIntegerQuotient(RefFromDouble<40>(value1), magnitude, negative, RawOf(f3))) << "value1=" << value1 << ", value2=" << value2 << ", got " << HexOf(f3);
    }
}
// Dividing by an integer used to truncate the quotient where dividing by the same value held in a
// fixed_point128 rounded it. Both round to nearest now, a tie away from zero.
TEST(fixed_point128, DivideByInt32)
{
    CheckDivideByScalar<int32_t>();
}
TEST(fixed_point128, DivideByUnsignedInt32)
{
    CheckDivideByScalar<uint32_t>();
}
// The divisors used to be 32 bit values "on purpose"; the full 64 bit range takes the 128 by 64 bit
// division through divisors that need all of its width.
TEST(fixed_point128, DivideByInt64)
{
    CheckDivideByScalar<int64_t>();
}
TEST(fixed_point128, DivideByUnsignedInt64)
{
    CheckDivideByScalar<uint64_t>();
}
/// @brief Checks the three division paths and the scalar overloads against the reference, for one I.
template <int32_t I> static void CheckDivisionAgainstReference()
{
    using fp = fixed_point128<I>;
    constexpr int32_t F = 127 - I;
    srand(RANDOM_SEED);
    int32_t checked = 0;
    for (auto i = 0u; i < RANDOM_TEST_COUNT / 2; ++i) {
        const Raw128 a = RandomRaw();
        const fp x = FromRaw<I>(a);

        // a divisor with a fraction, which takes the long division
        const Raw128 b = RandomRaw();
        if ((b.low | b.high) != 0 && !RefQuotientOverflows<I>(a, b)) {
            ++checked;
            const fp q = x / FromRaw<I>(b);
            EXPECT_TRUE(RefIsQuotient<I>(a, b, RawOf(q))) << "I=" << I << ", x=" << HexOf(x) << ", y=" << HexOf(FromRaw<I>(b)) << ", got " << HexOf(q);
        }

        // a divisor holding an integer, which takes div_64bit
        const int64_t k = 3 + static_cast<int64_t>(get_uint32_random() % 1000);
        if constexpr (I >= 11) {
            const Raw128 kr = RefFromInt<I>(get_random_sign() * k);
            if (!RefQuotientOverflows<I>(a, kr)) {
                const fp q = x / FromRaw<I>(kr);
                EXPECT_TRUE(RefIsQuotient<I>(a, kr, RawOf(q))) << "I=" << I << ", x=" << HexOf(x) << ", k=" << k;
            }
        }
        // the same integer as a scalar
        const fp qk = x / k;
        EXPECT_TRUE(RefIsIntegerQuotient(a, static_cast<uint64_t>(k), false, RawOf(qk))) << "I=" << I << ", x=" << HexOf(x) << ", k=" << k;
        EXPECT_TRUE(RefIsIntegerQuotient(a, static_cast<uint64_t>(k), true, RawOf(x / -k))) << "I=" << I << ", x=" << HexOf(x) << ", k=" << k;

        // a power of two divisor 2^j, which becomes a shift: rounded for j > 0, exact for j <= 0
        const int32_t j = static_cast<int32_t>(get_uint32_random() % static_cast<uint32_t>(I + F)) - F;  // [-F, I)
        const fp power = FromRaw<I>(RefShl(Raw128 {1, 0}, F + j));
        const Raw128 expected = (j > 0) ? RefShiftRightRounded(a, j) : RefShl(a, -j);
        if (j > 0 || !RefQuotientOverflows<I>(a, RawOf(power))) {
            EXPECT_TRUE(RefEqual(RawOf(x / power), expected)) << "I=" << I << ", x=" << HexOf(x) << ", j=" << j;
            EXPECT_TRUE(RefEqual(RawOf(x / std::ldexp(1.0, j)), expected)) << "I=" << I << ", x=" << HexOf(x) << ", j=" << j;
            EXPECT_TRUE(RefEqual(RawOf(x / -std::ldexp(1.0, j)), RefNegate(expected))) << "I=" << I << ", x=" << HexOf(x) << ", j=" << j;
        }
    }
    EXPECT_GT(checked, static_cast<int32_t>(RANDOM_TEST_COUNT) / 4) << "I=" << I;

    // an object divided by itself in place
    fp z = fp::one() + fp::half();
    z /= z;
    EXPECT_TRUE(z == fp::one()) << "I=" << I;
}
// Random raw operands at both ends of the width range, through every division path: the long
// division, an integer held in a fixed_point128, an integer scalar of either sign, and a power of two
// given as a fixed_point128 or as a double. I = 63 is the one width whose long division keeps its
// quotient in whole QWORDs, a branch of its own that no test reached.
TEST(fixed_point128, DivideIsCorrectlyRounded)
{
    CheckDivisionAgainstReference<1>();
    CheckDivisionAgainstReference<2>();
    CheckDivisionAgainstReference<32>();
    CheckDivisionAgainstReference<62>();
    CheckDivisionAgainstReference<63>();
}
// Values from Python fractions. The integer divisor truncated, the fixed_point128 one rounded; and an
// exact tie goes away from zero on every path, as it does in the multiply.
TEST(fixed_point128, IntegerScalarDivisionRounds)
{
    using fp = fixed_point128<32>;
    const Raw128 third {0xAAAAAAAAAAAAAAABull, 0x000000002AAAAAAAull};
    EXPECT_TRUE(RefEqual(RawOf(fp(1) / 3), third));
    EXPECT_TRUE(RefEqual(RawOf(fp(1) / 3ull), third));
    EXPECT_TRUE(RefEqual(RawOf(fp(1) / static_cast<int64_t>(3)), third));
    EXPECT_TRUE(RefEqual(RawOf(fp(1) / fp(3)), third));
    EXPECT_TRUE(RefEqual(RawOf(fp(1) / 3.0), third));
    EXPECT_TRUE(RefEqual(RawOf(fp(-1) / 3), RefNegate(third)));
    EXPECT_TRUE(RefEqual(RawOf(fp(1) / -3), RefNegate(third)));
    // ties: 1.5, 2.5 and -1.5 units of the last place
    EXPECT_TRUE(RefEqual(RawOf(fp(3, 0) / 2), Raw128 {2, 0}));
    EXPECT_TRUE(RefEqual(RawOf(fp(5, 0) / 2), Raw128 {3, 0}));
    EXPECT_TRUE(RefEqual(RawOf(-fp(3, 0) / 2), RefNegate(Raw128 {2, 0})));
    EXPECT_TRUE(RefEqual(RawOf(fp(3, 0) / fp(2)), Raw128 {2, 0}));
    EXPECT_TRUE(RefEqual(RawOf(fp(5, 0) / fp(6)), Raw128 {1, 0}));  // 0.83 rounds up
    EXPECT_TRUE(RefEqual(RawOf(fp(3, 0) / fp(6)), Raw128 {1, 0}));  // the tie 0.5 too
}
// Every division by zero throws, zero divided by zero included, which used to return zero: the
// numerator's shortcut ran before the divisor was looked at. numeric_limits::traps says it throws.
template <int32_t I> static void CheckDivisionByZeroThrows()
{
    using fp = fixed_point128<I>;
    for (const fp numerator : {fp(0), fp::one(), -fp::one(), fp::epsilon()}) {
        EXPECT_THROW((void)(numerator / fp(0)), std::logic_error) << "I=" << I;
        EXPECT_THROW((void)(numerator / 0), std::logic_error) << "I=" << I;
        EXPECT_THROW((void)(numerator / 0ull), std::logic_error) << "I=" << I;
        EXPECT_THROW((void)(numerator / 0.0), std::logic_error) << "I=" << I;
        EXPECT_THROW((void)(numerator / 0.0f), std::logic_error) << "I=" << I;
        EXPECT_THROW((void)(numerator % fp(0)), std::logic_error) << "I=" << I;
        EXPECT_THROW((void)(numerator % 0), std::logic_error) << "I=" << I;
        EXPECT_THROW((void)fmod(numerator, fp(0)), std::logic_error) << "I=" << I;
        fp compound = numerator;
        EXPECT_THROW(compound /= fp(0), std::logic_error) << "I=" << I;
    }
}
TEST(fixed_point128, DivisionByZeroThrows)
{
    CheckDivisionByZeroThrows<1>();
    CheckDivisionByZeroThrows<32>();
    CheckDivisionByZeroThrows<63>();
}

/**********************************************************************
 * Modulo
 ***********************************************************************/
TEST(fixed_point128, ModuloByFP128)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = get_double_random(-1, 39);
        double value2 = get_double_random(-1, 39);
        double res = fmod(value1, value2);
        fixed_point128<40> f1 = value1;
        fixed_point128<40> f2 = value2;
        fixed_point128<40> f3 = f1 % f2;
        // the remainder is exact, a difference of representable values that a double holds too
        EXPECT_EQ(static_cast<double>(f3), res) << "value1=" << value1 << ", value2=" << value2;
        EXPECT_TRUE(fmod(f1, f2) == f3) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(fixed_point128, ModuloByDouble)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = get_double_random(-1, 39);
        double value2 = get_double_random(-1, 39);
        double res = fmod(value1, value2);
        fixed_point128<40> f1 = value1;
        fixed_point128<40> f3 = f1 % value2;
        EXPECT_EQ(static_cast<double>(f3), res) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(fixed_point128, ModuloByFloat)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        float value1 = (float)get_double_random(-1, 39);
        float value2 = (float)get_double_random(-1, 39);
        if (value2 == 0)
            continue;
        float res = fmodf(value1, value2);
        fixed_point128<40> f1 = value1;
        fixed_point128<40> f3 = f1 % value2;
        EXPECT_EQ(static_cast<float>(f3), res) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(fixed_point128, ModuloByInt32)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = get_double_random(-1, 39);
        int32_t value2 = get_int32_random();
        if (value2 == 0)
            continue;
        double res = fmod(value1, value2);
        fixed_point128<40> f1 = value1;
        fixed_point128<40> f3 = f1 % value2;
        EXPECT_EQ(static_cast<double>(f3), res) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(fixed_point128, ModuloByUnsignedInt32)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = get_double_random(-1, 39);
        uint32_t value2 = get_uint32_random();
        if (value2 == 0)
            continue;
        double res = fmod(value1, value2);
        fixed_point128<40> f1 = value1;
        fixed_point128<40> f3 = f1 % value2;
        EXPECT_EQ(static_cast<double>(f3), res) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(fixed_point128, ModuloByInt64)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = get_double_random(-1, 39);
        int64_t value2 = get_int32_random();  // a divisor beyond 2^40 leaves the type
        if (value2 == 0)
            continue;
        double res = fmod(value1, static_cast<double>(value2));
        fixed_point128<40> f1 = value1;
        fixed_point128<40> f3 = f1 % value2;
        EXPECT_EQ(static_cast<double>(f3), res) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(fixed_point128, ModuloByUnsignedInt64)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = get_double_random(-1, 39);
        uint64_t value2 = get_uint32_random();  // see ModuloByInt64
        if (value2 == 0)
            continue;
        double res = fmod(value1, static_cast<double>(value2));
        fixed_point128<40> f1 = value1;
        fixed_point128<40> f3 = f1 % value2;
        EXPECT_EQ(static_cast<double>(f3), res) << "value1=" << value1 << ", value2=" << value2;
    }
}
/// @brief Checks the remainder of random raw operands against a bitwise long division, for one I.
template <int32_t I> static void CheckModuloAgainstReference()
{
    using fp = fixed_point128<I>;
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT / 16; ++i) {
        const Raw128 a = RandomRaw(), b = RandomRaw();
        if ((b.low | b.high) == 0)
            continue;
        Raw128 quotient {}, remainder {};
        RefDivide(RefMagnitude(a), RefMagnitude(b), quotient, remainder);
        const Raw128 expected = RefIsNegative(a) ? RefNegate(remainder) : remainder;
        const fp x = FromRaw<I>(a), y = FromRaw<I>(b);
        EXPECT_TRUE(RefEqual(RawOf(x % y), expected)) << "I=" << I << ", x=" << HexOf(x) << ", y=" << HexOf(y) << ", got " << HexOf(x % y);
        EXPECT_TRUE(fmod(x, y) == x % y) << "I=" << I;
    }
}
// x % y is the remainder of a truncating division and keeps the sign of x. It used to be computed from
// the rounded quotient: one that rounded up onto an integer made the remainder zero, a whole divisor
// short. Values from Python fractions.
TEST(fixed_point128, ModuloMatchesExactRemainder)
{
    using fp = fixed_point128<20>;
    const fp eps = fp::epsilon();
    EXPECT_TRUE(RefEqual(RawOf((fp(9) - eps) % fp(3)), Raw128 {0xFFFFFFFFFFFFFFFFull, 0x000017FFFFFFFFFFull}));
    EXPECT_TRUE(RefEqual(RawOf(fmod(fp(9) - eps, fp(3))), Raw128 {0xFFFFFFFFFFFFFFFFull, 0x000017FFFFFFFFFFull}));
    EXPECT_TRUE(RefEqual(RawOf((eps - fp(9)) % fp(3)), Raw128 {0x0000000000000001ull, 0xFFFFE80000000000ull}));
    EXPECT_TRUE(RefEqual(RawOf((fp(7.5) - eps) % fp(2.5)), Raw128 {0xFFFFFFFFFFFFFFFFull, 0x000013FFFFFFFFFFull}));
    // the sign follows the numerator only
    EXPECT_TRUE(fp(-7) % fp(2) == fp(-1));
    EXPECT_TRUE(fp(7) % fp(-2) == fp(1));
    EXPECT_TRUE(fp(-7.5) % fp(-2) == fp(-1.5));
    // zero modulo anything is zero, and anything modulo itself is zero
    EXPECT_TRUE(fp(0) % fp(3) == fp(0));
    EXPECT_TRUE(fp(2.5) % fp(2.5) == fp(0));
    // INT64_MIN % -1 sat in the old integer shortcut: lowest() of fixed_point128<63> % -1 killed the
    // process with an integer overflow exception
    EXPECT_TRUE(std::numeric_limits<fixed_point128<63>>::lowest() % fixed_point128<63>(-1) == fixed_point128<63>(0));
    EXPECT_TRUE(std::numeric_limits<fixed_point128<63>>::lowest() % fixed_point128<63>(3) == fixed_point128<63>(-2));

    CheckModuloAgainstReference<1>();
    CheckModuloAgainstReference<20>();
    CheckModuloAgainstReference<63>();
}

/**********************************************************************
 * Comparisons
 ***********************************************************************/
TEST(fixed_point128, CompareFP128)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = get_double_random(-1, 39);
        double value2 = get_double_random(-1, 39);
        fixed_point128<40> f1 = value1;
        fixed_point128<40> f2 = value2;
        EXPECT_EQ(f1 > f2, value1 > value2) << "value1=" << value1 << ", value2=" << value2;
        EXPECT_EQ(f1 >= f2, value1 >= value2) << "value1=" << value1 << ", value2=" << value2;
        EXPECT_EQ(f1 < f2, value1 < value2) << "value1=" << value1 << ", value2=" << value2;
        EXPECT_EQ(f1 <= f2, value1 <= value2) << "value1=" << value1 << ", value2=" << value2;
    }
}
/// @brief Compares a fixed_point128<40> with scalars of type T, on both sides, against the scalars' own order.
template <typename T, typename Gen> static void CheckCompareWithScalar(Gen generate)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const T value1 = generate();
        const T value2 = (i % 8 == 0) ? value1 : generate();  // equal operands an eighth of the time
        const fixed_point128<40> f1 = value1;
        EXPECT_EQ(f1 > value2, value1 > value2) << "value1=" << value1 << ", value2=" << value2;
        EXPECT_EQ(f1 >= value2, value1 >= value2) << "value1=" << value1 << ", value2=" << value2;
        EXPECT_EQ(f1 < value2, value1 < value2) << "value1=" << value1 << ", value2=" << value2;
        EXPECT_EQ(f1 <= value2, value1 <= value2) << "value1=" << value1 << ", value2=" << value2;
        EXPECT_EQ(f1 == value2, value1 == value2) << "value1=" << value1 << ", value2=" << value2;
        EXPECT_EQ(f1 != value2, value1 != value2) << "value1=" << value1 << ", value2=" << value2;
        // the scalar on the left: these overloads ran only in constant expressions
        EXPECT_EQ(value2 < f1, value2 < value1) << "value1=" << value1 << ", value2=" << value2;
        EXPECT_EQ(value2 <= f1, value2 <= value1) << "value1=" << value1 << ", value2=" << value2;
        EXPECT_EQ(value2 > f1, value2 > value1) << "value1=" << value1 << ", value2=" << value2;
        EXPECT_EQ(value2 >= f1, value2 >= value1) << "value1=" << value1 << ", value2=" << value2;
        EXPECT_EQ(value2 == f1, value2 == value1) << "value1=" << value1 << ", value2=" << value2;
        EXPECT_EQ(value2 != f1, value2 != value1) << "value1=" << value1 << ", value2=" << value2;
    }
}
TEST(fixed_point128, CompareDouble)
{
    CheckCompareWithScalar<double>([] { return get_double_random(-10, 39); });
}
TEST(fixed_point128, CompareFloat)
{
    CheckCompareWithScalar<float>([] { return static_cast<float>(get_double_random(-10, 39)); });
}
TEST(fixed_point128, CompareInt32)
{
    CheckCompareWithScalar<int32_t>([] { return get_int32_random(); });
}
TEST(fixed_point128, CompareUnsignedInt32)
{
    CheckCompareWithScalar<uint32_t>([] { return get_uint32_random(); });
}
// A full width 64 bit random is never representable in fixed_point128<40>: the guard used to reject all
// 65536 iterations. Shifting by an independent random amount bounds each operand to 2^39 while keeping
// a spread of magnitudes.
TEST(fixed_point128, CompareInt64)
{
    CheckCompareWithScalar<int64_t>([] { return get_int64_random() >> (24 + static_cast<int32_t>(get_uint32_random() % 40)); });
}
TEST(fixed_point128, CompareUnsignedInt64)
{
    CheckCompareWithScalar<uint64_t>([] { return get_uint64_random() >> (25 + static_cast<int32_t>(get_uint32_random() % 39)); });
}
// operator== and operator!= used to be tested on an object against itself, which a broken comparison
// passes as well as a correct one. Two values one bit apart have to compare unequal, for every bit,
// and order the way their two's complement patterns do.
template <int32_t I> static void CheckEqualityDistinguishesEveryBit()
{
    using fp = fixed_point128<I>;
    srand(RANDOM_SEED);
    for (auto i = 0u; i < 64; ++i) {
        const Raw128 a = {get_uint64_random(), get_uint64_random()};
        const fp x = FromRaw<I>(a);
        EXPECT_TRUE(x == FromRaw<I>(a)) << "I=" << I;
        EXPECT_FALSE(x != FromRaw<I>(a)) << "I=" << I;
        for (int32_t bit = 0; bit < 128; ++bit) {
            Raw128 b = a;
            if (bit < 64) {
                b.low ^= 1ull << bit;
            } else {
                b.high ^= 1ull << (bit - 64);
            }
            const fp y = FromRaw<I>(b);
            const bool less = RefSignedLess(a, b);
            EXPECT_FALSE(x == y) << "I=" << I << ", bit=" << bit;
            EXPECT_TRUE(x != y) << "I=" << I << ", bit=" << bit;
            EXPECT_EQ(x < y, less) << "I=" << I << ", bit=" << bit;
            EXPECT_EQ(x > y, !less) << "I=" << I << ", bit=" << bit;
            EXPECT_EQ(x <= y, less) << "I=" << I << ", bit=" << bit;
            EXPECT_EQ(x >= y, !less) << "I=" << I << ", bit=" << bit;
        }
    }
}
TEST(fixed_point128, OperatorEqual)
{
    CheckEqualityDistinguishesEveryBit<1>();
    CheckEqualityDistinguishesEveryBit<40>();
}
TEST(fixed_point128, OperatorNotEqual)
{
    CheckEqualityDistinguishesEveryBit<63>();
}
// The scalar forms of == and != used to compare a value with itself converted the same way on both
// sides, and cast doubles out of the range of the integer types, which is undefined behavior. The
// fixed_point128 side is spelled from raw bits here, so neither side goes through the other's conversion.
TEST(fixed_point128, TemplateOperatorEqual)
{
    using fp = fixed_point128<40>;
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const int32_t k = get_int32_random();
        const int32_t other = (k < 0) ? k + 1 : k - 1;
        const fp x = FromRaw<40>(RefFromInt<40>(k));
        EXPECT_TRUE(x == k) << "k=" << k;
        EXPECT_TRUE(k == x) << "k=" << k;
        EXPECT_TRUE(x == static_cast<int64_t>(k)) << "k=" << k;
        EXPECT_TRUE(x == static_cast<double>(k)) << "k=" << k;
        EXPECT_FALSE(x == other) << "k=" << k;
        EXPECT_FALSE(other == x) << "k=" << k;

        const uint32_t u = get_uint32_random();
        const fp xu = FromRaw<40>(RefFromInt<40>(static_cast<int64_t>(u)));
        EXPECT_TRUE(xu == u) << "u=" << u;
        EXPECT_TRUE(xu == static_cast<uint64_t>(u)) << "u=" << u;
        EXPECT_FALSE(xu == u + 1ull) << "u=" << u;

        const double d = get_double_random(-10, 39);
        const fp xd = FromRaw<40>(RefFromDouble<40>(d));
        EXPECT_TRUE(xd == d) << "d=" << d;
        EXPECT_TRUE(d == xd) << "d=" << d;
        EXPECT_FALSE(xd == std::nextafter(d, 0.0)) << "d=" << d;
        EXPECT_TRUE(xd == static_cast<float>(d) || static_cast<double>(static_cast<float>(d)) != d) << "d=" << d;
    }
}
TEST(fixed_point128, TemplateOperatorNotEqual)
{
    using fp = fixed_point128<40>;
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const int32_t k = get_int32_random();
        const int32_t other = (k < 0) ? k + 1 : k - 1;
        const fp x = FromRaw<40>(RefFromInt<40>(k));
        EXPECT_FALSE(x != k) << "k=" << k;
        EXPECT_FALSE(k != x) << "k=" << k;
        EXPECT_TRUE(x != other) << "k=" << k;
        EXPECT_TRUE(other != x) << "k=" << k;
        EXPECT_TRUE(x != static_cast<int64_t>(other)) << "k=" << k;

        const double d = get_double_random(-10, 39);
        const fp xd = FromRaw<40>(RefFromDouble<40>(d));
        EXPECT_FALSE(xd != d) << "d=" << d;
        EXPECT_FALSE(d != xd) << "d=" << d;
        EXPECT_TRUE(xd != std::nextafter(d, 0.0)) << "d=" << d;
        EXPECT_TRUE(std::nextafter(d, 0.0) != xd) << "d=" << d;
    }
}

/**********************************************************************
 * Increment, decrement and shifts
 ***********************************************************************/
TEST(fixed_point128, OperatorPlusPlus)
{
    using fp = fixed_point128<40>;
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const double value1 = get_double_random(-10, 38);
        const fp f1 = value1;
        const Raw128 expected = RefAdd(RefFromDouble<40>(value1), RefFromInt<40>(1));
        fp f2 = f1;
        // the postfix form hands back the value from before the increment
        EXPECT_TRUE(f2++ == f1) << "value1=" << value1;
        EXPECT_TRUE(RefEqual(RawOf(f2), expected)) << "value1=" << value1;
        f2 = f1;
        EXPECT_TRUE(RefEqual(RawOf(++f2), expected)) << "value1=" << value1;
        EXPECT_TRUE(RefEqual(RawOf(f2), expected)) << "value1=" << value1;
    }
}
TEST(fixed_point128, OperatorMinusMinus)
{
    using fp = fixed_point128<40>;
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const double value1 = get_double_random(-10, 38);
        const fp f1 = value1;
        const Raw128 expected = RefSub(RefFromDouble<40>(value1), RefFromInt<40>(1));
        fp f2 = f1;
        EXPECT_TRUE(f2-- == f1) << "value1=" << value1;
        EXPECT_TRUE(RefEqual(RawOf(f2), expected)) << "value1=" << value1;
        f2 = f1;
        EXPECT_TRUE(RefEqual(RawOf(--f2), expected)) << "value1=" << value1;
        EXPECT_TRUE(RefEqual(RawOf(f2), expected)) << "value1=" << value1;
    }
}
TEST(fixed_point128, ShiftRight)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = get_uint32_random();
        uint32_t shift = get_uint32_random() % 40u;
        fixed_point128<40> f1 = value1;

        double expo = pow(2, -static_cast<int32_t>(shift));
        double res = value1 * expo;  // no set bit is shifted out, so this one is exact
        fixed_point128 f3 = f1 >> shift;

        if (check_overflow(value1, f1) || check_overflow(res, f3)) {
            continue;
        }
        EXPECT_DOUBLE_EQ(static_cast<double>(f3), res) << "double value1=" << value1 << "; uint32_t shift=" << shift << ";";
    }
}
TEST(fixed_point128, ShiftLeft)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        double value1 = 1.0 / get_int32_random();
        uint32_t shift = get_uint32_random() % 40u;
        fixed_point128<40> f1 = value1;
        double expo = pow(2, static_cast<int32_t>(shift));
        double res = value1 * expo;  // double doesn't lose any bits with this operation!
        fixed_point128 f3 = f1 << shift;
        if (check_overflow(value1, f1) || check_overflow(res, f3)) {
            continue;
        }
        EXPECT_DOUBLE_EQ(static_cast<double>(f3), res) << "double value1=" << value1 << "; uint32_t shift=" << shift << ";";
    }
}
/// @brief Checks x >> s and x >>= s for every count from 0 to 130, on random raw values, for one I.
template <int32_t I> static void CheckShiftRightAgainstReference()
{
    using fp = fixed_point128<I>;
    srand(RANDOM_SEED);
    for (auto i = 0u; i < 512; ++i) {
        const Raw128 a = (i % 2 == 0) ? RandomRaw() : Raw128 {get_uint64_random(), get_uint64_random()};
        const fp x = FromRaw<I>(a);
        for (int32_t s = 0; s <= 130; ++s) {
            const Raw128 expected = RefShiftRightRounded(a, s);
            EXPECT_TRUE(RefEqual(RawOf(x >> s), expected)) << "I=" << I << ", x=" << HexOf(x) << ", s=" << s << ", got " << HexOf(x >> s);
            fp y = x;
            y >>= s;
            EXPECT_TRUE(RefEqual(RawOf(y), expected)) << "I=" << I << ", x=" << HexOf(x) << ", s=" << s;
        }
    }
}
// The shift rounded to nearest with only one of the dropped bits below the rounding bit as its sticky
// bit, so a dropped part above one half whose second highest bit was clear rounded down: raw 5 >> 3 is
// 0.625 and came out as 0. Every count is checked against the magnitude rounded half up - a tie goes
// away from zero, the rule of the multiply and the divisions.
TEST(fixed_point128, ShiftRightRoundsToNearest)
{
    using fp = fixed_point128<32>;
    EXPECT_TRUE(RefEqual(RawOf(fp(5, 0) >> 3), Raw128 {1, 0}));
    EXPECT_TRUE(RefEqual(RawOf(fp(21, 0) >> 3), Raw128 {3, 0}));
    EXPECT_TRUE(RefEqual(RawOf(fp(0x8000000000000001ull, 2) >> 66), Raw128 {1, 0}));
    EXPECT_TRUE(RefEqual(RawOf(-fp(21, 0) >> 3), RefNegate(Raw128 {3, 0})));
    // ties go away from zero
    EXPECT_TRUE(RefEqual(RawOf(fp(4, 0) >> 3), Raw128 {1, 0}));
    EXPECT_TRUE(RefEqual(RawOf(fp(12, 0) >> 3), Raw128 {2, 0}));
    EXPECT_TRUE(RefEqual(RawOf(fp(20, 0) >> 3), Raw128 {3, 0}));
    EXPECT_TRUE(RefEqual(RawOf(-fp(4, 0) >> 3), RefNegate(Raw128 {1, 0})));
    EXPECT_TRUE(RefEqual(RawOf(-fp(12, 0) >> 3), RefNegate(Raw128 {2, 0})));
    // which makes +-epsilon a fixed point of a shift by one, as documented
    EXPECT_TRUE((fp::epsilon() >> 1) == fp::epsilon());
    EXPECT_TRUE((-fp::epsilon() >> 1) == -fp::epsilon());
    // 128 and beyond: zero, but for lowest() >> 128, which is the tie -epsilon / 2
    EXPECT_TRUE((std::numeric_limits<fp>::lowest() >> 128) == -fp::epsilon());
    EXPECT_TRUE((std::numeric_limits<fp>::lowest() >> 129) == fp(0));
    EXPECT_TRUE((std::numeric_limits<fp>::max() >> 128) == fp(0));

    CheckShiftRightAgainstReference<1>();
    CheckShiftRightAgainstReference<32>();
    CheckShiftRightAgainstReference<63>();
}

/**********************************************************************
 * Rounding functions
 ***********************************************************************/
/**
 * @brief floor, ceil, trunc, round and copysign of random doubles on the grid, bit for bit against the
 *        same functions of the double: the results are integers, or the operand itself, which the
 *        double holds exactly.
 * @param minExponent Smallest exponent whose 53 bit mantissa still ends on the grid of fixed_point128<I>.
 */
template <int32_t I> static void CheckRoundingFunctions(int32_t minExponent)
{
    using fp = fixed_point128<I>;
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT / 4; ++i) {
        // below 2^(I-1), so that rounding up stays in range
        const double value = get_double_random(minExponent, I - 1);
        const double other = get_double_random(minExponent, I - 1);
        const fp f = value;
        EXPECT_TRUE(RefEqual(RawOf(floor(f)), RefFromDouble<I>(::floor(value)))) << "floor, I=" << I << ", value=" << value;
        EXPECT_TRUE(RefEqual(RawOf(ceil(f)), RefFromDouble<I>(::ceil(value)))) << "ceil, I=" << I << ", value=" << value;
        EXPECT_TRUE(RefEqual(RawOf(trunc(f)), RefFromDouble<I>(::trunc(value)))) << "trunc, I=" << I << ", value=" << value;
        EXPECT_TRUE(RefEqual(RawOf(round(f)), RefFromDouble<I>(::round(value)))) << "round, I=" << I << ", value=" << value;
        EXPECT_TRUE(RefEqual(RawOf(copysign(f, fp(other))), RefFromDouble<I>(::copysign(value, other)))) << "copysign, I=" << I << ", value=" << value;
        fp int_part;
        const fp frac = modf(f, &int_part);
        double int_reference = 0;
        const double frac_reference = ::modf(value, &int_reference);
        EXPECT_TRUE(RefEqual(RawOf(int_part), RefFromDouble<I>(int_reference))) << "modf, I=" << I << ", value=" << value;
        EXPECT_TRUE(RefEqual(RawOf(frac), RefFromDouble<I>(frac_reference))) << "modf, I=" << I << ", value=" << value;
    }
}
TEST(fixed_point128, floor)
{
    CheckRoundingFunctions<2>(-72);
    CheckRoundingFunctions<16>(-58);
    CheckRoundingFunctions<32>(-42);
    CheckRoundingFunctions<63>(-11);
}
// ceil() returned an integer argument unchanged through a branch no test reached: only fractions were
// ever passed to it. Integers are their own floor, ceiling, truncation and rounding at every width.
template <int32_t I> static void CheckRoundingOfIntegers()
{
    using fp = fixed_point128<I>;
    const int64_t top = static_cast<int64_t>(fp::max_int_value);
    for (const int64_t k : {-top - 1, -top, int64_t(-3), int64_t(-2), int64_t(-1), int64_t(0), int64_t(1), int64_t(2), int64_t(5), top}) {
        if (k < -top - 1 || k > top)
            continue;
        const fp x = FromRaw<I>(RefFromInt<I>(k));
        EXPECT_TRUE(ceil(x) == x) << "I=" << I << ", k=" << k;
        EXPECT_TRUE(floor(x) == x) << "I=" << I << ", k=" << k;
        EXPECT_TRUE(trunc(x) == x) << "I=" << I << ", k=" << k;
        if (k > -top - 1) {
            // round() works on fabs(), which lowest() survives only as itself
            EXPECT_TRUE(round(x) == x) << "I=" << I << ", k=" << k;
        }
        fp int_part;
        EXPECT_TRUE(modf(x, &int_part) == fp(0)) << "I=" << I << ", k=" << k;
        EXPECT_TRUE(int_part == x) << "I=" << I << ", k=" << k;
    }
}
TEST(fixed_point128, ceil)
{
    CheckRoundingOfIntegers<1>();
    CheckRoundingOfIntegers<32>();
    CheckRoundingOfIntegers<63>();
}
TEST(fixed_point128, trunc)
{
    CheckRoundingFunctions<1>(-73);
}
TEST(fixed_point128, round)
{
    using fp = fixed_point128<16>;
    // the halfway values round away from zero
    EXPECT_TRUE(round(fp(2.5)) == fp(3));
    EXPECT_TRUE(round(fp(-2.5)) == fp(-3));
    EXPECT_TRUE(round(fp(0.5)) == fp(1));
    EXPECT_TRUE(round(fp(-0.5)) == fp(-1));
    EXPECT_TRUE(round(fp(0.5) - fp::epsilon()) == fp(0));
}
TEST(fixed_point128, copysign)
{
    using fp = fixed_point128<16>;
    EXPECT_TRUE(copysign(fp(3), fp(-1)) == fp(-3));
    EXPECT_TRUE(copysign(fp(-3), fp(1)) == fp(3));
    EXPECT_TRUE(copysign(fp(-3), fp(0)) == fp(3));  // zero is positive, there is no other
    EXPECT_TRUE(copysign(fp(0), fp(-1)) == fp(0));
}

/**********************************************************************
 * Reciprocal, square root, hypot, ilogb and logb: exact references
 ***********************************************************************/
/// @brief Checks reciprocal() of random raw values against the exact quotient, for one I.
template <int32_t I> static void CheckReciprocal()
{
    using fp = fixed_point128<I>;
    constexpr int32_t F = 127 - I;
    const Raw128 one = RefShl(Raw128 {1, 0}, F);
    srand(RANDOM_SEED);
    int32_t checked = 0;
    for (auto i = 0u; i < RANDOM_TEST_COUNT / 2; ++i) {
        // magnitudes from about 2^-I up, where most reciprocals fit, and some below, where they do not
        const int32_t shift = 1 + static_cast<int32_t>(get_uint32_random() % static_cast<uint32_t>(2 * I + 4));
        Raw128 a = RefShr(Raw128 {get_uint64_random(), get_uint64_random()}, std::min(shift, 127));
        if (rand() & 1) {
            a = RefNegate(a);
        }
        if ((a.low | a.high) == 0)
            continue;
        const fp x = FromRaw<I>(a);
        const fp r = reciprocal(x);
        if (RefQuotientOverflows<I>(one, a)) {
            // 1 / x does not fit: documented to give zero, or lowest() for x = -2^-I exactly
            EXPECT_TRUE(r == fp(0) || (r == std::numeric_limits<fp>::lowest() && x == -(fp::one() >> I))) << "I=" << I << ", x=" << HexOf(x);
            continue;
        }
        ++checked;
        EXPECT_TRUE(RefIsQuotient<I>(one, a, RawOf(r))) << "I=" << I << ", x=" << HexOf(x) << ", got " << HexOf(r);
    }
    EXPECT_GT(checked, static_cast<int32_t>(RANDOM_TEST_COUNT) / 16) << "I=" << I;
}
// reciprocal() refined a double estimate with Newton steps: a test against the double reference could
// not tell the refined result from the bare estimate, and the refinement lost accuracy as 1 / x grew -
// tens of thousands of ulp at fixed_point128<16>. It is the correctly rounded quotient now.
TEST(fixed_point128, reciprocal)
{
    CheckReciprocal<1>();
    CheckReciprocal<16>();
    CheckReciprocal<32>();
    CheckReciprocal<63>();
}
// Documented to return zero where 1 / x leaves the range; it returned whatever the Newton steps made of
// a saturated estimate instead.
TEST(fixed_point128, ReciprocalOutOfRange)
{
    using fp = fixed_point128<16>;
    EXPECT_TRUE(reciprocal(fp(0)) == fp(0));
    EXPECT_TRUE(reciprocal(fp::epsilon()) == fp(0));
    EXPECT_TRUE(reciprocal(fp::one() >> 20) == fp(0));
    EXPECT_TRUE(reciprocal(fp::one() >> 16) == fp(0));  // 2^16 is one past max()
    EXPECT_TRUE(reciprocal(-(fp::one() >> 16)) == std::numeric_limits<fp>::lowest());
    EXPECT_TRUE(reciprocal(fp::one() >> 15) == fp(32768));
    EXPECT_TRUE(reciprocal(fp(-4)) == -(fp::one() >> 2));
}
/**
 * @brief Checks that r is sqrt(s) rounded to nearest: r^2 - r < s <= r^2 + r, all integers, the lower
 *        bound dropping out for r = 0.
 * @param s The radicand, 256 bits
 * @param r The candidate root, unsigned
 */
[[nodiscard]] static bool RefIsRoundedRoot(const U256& s, const Raw128& r)
{
    if ((r.low | r.high) == 0) {
        return (s.w[0] | s.w[1] | s.w[2] | s.w[3]) == 0;
    }
    const U256 square = RefMul(r, r);
    const U256 r256 = {{r.low, r.high, 0, 0}};
    const U256 above = RefAdd256(square, r256);
    // r^2 - r < s is r^2 < s + r
    return RefCompare(square, RefAdd256(s, r256)) < 0 && RefCompare(s, above) <= 0;
}
/// @brief Checks sqrt() of random positive raw values against the exact integer root, for one I.
template <int32_t I> static void CheckSqrt()
{
    using fp = fixed_point128<I>;
    constexpr int32_t F = 127 - I;
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT / 2; ++i) {
        const Raw128 a = RefMagnitude(RandomRaw());
        const fp x = FromRaw<I>(a);
        const fp root = sqrt(x);
        // sqrt(raw / 2^F) * 2^F is the root of raw * 2^F
        EXPECT_TRUE(RefIsRoundedRoot(RefWiden(a, F), RawOf(root))) << "I=" << I << ", x=" << HexOf(x) << ", got " << HexOf(root);
        // the step count is ignored: the result is already correctly rounded
        EXPECT_TRUE(sqrt(x, 0) == root) << "I=" << I << ", x=" << HexOf(x);
    }
    // perfect squares come out exact, and the documented zero for a non positive argument
    for (int64_t k = 0; k * k <= static_cast<int64_t>(fp::max_int_value) && k < 3000000000; k = 2 * k + 1) {
        EXPECT_TRUE(sqrt(fp(k * k)) == fp(k)) << "I=" << I << ", k=" << k;
    }
    EXPECT_TRUE(sqrt(fp::one() >> 2) == fp::half()) << "I=" << I;
    EXPECT_TRUE(sqrt(-fp::one()) == fp(0)) << "I=" << I;
    EXPECT_TRUE(sqrt(std::numeric_limits<fp>::lowest()) == fp(0)) << "I=" << I;
    EXPECT_TRUE(sqrt(fp::epsilon()) == FromRaw<I>(RefShl(Raw128 {1, 0}, F / 2)) || (F % 2) != 0) << "I=" << I;
}
// sqrt() refined a double estimate after normalizing the argument with a right shift, and scaled the
// root back with a left shift: up to half of the integer bits were lost, 2.9e9 ulp at I = 63. Its test
// compared with a double, which a root without any refinement passed as well. Checked here against the
// exact integer root at both ends of the width range.
TEST(fixed_point128, sqrt)
{
    CheckSqrt<1>();
    CheckSqrt<2>();
    CheckSqrt<16>();
    CheckSqrt<32>();
    CheckSqrt<63>();
}
/// @brief Checks hypot() of random raw values whose result fits against the exact root, for one I.
template <int32_t I> static void CheckHypot()
{
    using fp = fixed_point128<I>;
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT / 2; ++i) {
        // magnitudes below 2^126 each, so the root stays below 2^127, either sign
        Raw128 a = RefShr(RefMagnitude(RandomRaw()), 1), b = RefShr(RefMagnitude(RandomRaw()), 1);
        if (rand() & 1) {
            a = RefNegate(a);
        }
        if (rand() & 1) {
            b = RefNegate(b);
        }
        const fp x = FromRaw<I>(a), y = FromRaw<I>(b);
        const fp h = hypot(x, y);
        const U256 s = RefAdd256(RefMul(RefMagnitude(a), RefMagnitude(a)), RefMul(RefMagnitude(b), RefMagnitude(b)));
        EXPECT_TRUE(RefIsRoundedRoot(s, RawOf(h))) << "I=" << I << ", x=" << HexOf(x) << ", y=" << HexOf(y) << ", got " << HexOf(h);
    }
}
TEST(fixed_point128, hypot)
{
    CheckHypot<1>();
    CheckHypot<16>();
    CheckHypot<32>();
    CheckHypot<63>();
}
// hypot() squared its arguments in the type itself, which overflowed once one of them reached 2^(I/2):
// hypot(300, 400) was zero at fixed_point128<16>, though 500 is well inside the range.
TEST(fixed_point128, HypotHasNoIntermediateOverflow)
{
    using fp = fixed_point128<16>;
    EXPECT_TRUE(hypot(fp(300), fp(400)) == fp(500));
    EXPECT_TRUE(hypot(fp(-300), fp(400)) == fp(500));
    EXPECT_TRUE(hypot(fp(300), fp(0)) == fp(300));
    EXPECT_TRUE(hypot(fp(0), fp(-300)) == fp(300));
    EXPECT_TRUE(hypot(std::numeric_limits<fp>::max(), fp(0)) == std::numeric_limits<fp>::max());
    EXPECT_TRUE(hypot(fixed_point128<63>(3000000000ll), fixed_point128<63>(4000000000ll)) == fixed_point128<63>(5000000000ll));
}
/// @brief Checks ilogb() and logb() of random raw values against the position of the leading one.
template <int32_t I> static void CheckIlogb()
{
    using fp = fixed_point128<I>;
    constexpr int32_t F = 127 - I;
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT / 4; ++i) {
        const Raw128 a = RandomRaw();
        if ((a.low | a.high) == 0)
            continue;
        const Raw128 m = RefMagnitude(a);
        int32_t msb = 127;
        while (RefBit(m, msb) == 0) {
            --msb;
        }
        const fp x = FromRaw<I>(a);
        EXPECT_EQ(ilogb(x), msb - F) << "I=" << I << ", x=" << HexOf(x);
        if (static_cast<int64_t>(msb - F) >= -static_cast<int64_t>(fp::max_int_value) - 1) {
            EXPECT_TRUE(logb(x) == fp(msb - F)) << "I=" << I << ", x=" << HexOf(x);
        }
    }
}
// logb() is documented, like the C function, to ignore the sign of its argument, and threw for a
// negative one; ilogb(0) returned I - 128 instead of FP_ILOGB0.
TEST(fixed_point128, ilogb)
{
    CheckIlogb<1>();
    CheckIlogb<16>();
    CheckIlogb<63>();
    using fp = fixed_point128<16>;
    EXPECT_EQ(ilogb(fp(0)), FP_ILOGB0);
    EXPECT_EQ(ilogb(fp(-8)), 3);
    EXPECT_EQ(ilogb(fp(0.75)), -1);
}
TEST(fixed_point128, logb)
{
    using fp = fixed_point128<16>;
    EXPECT_TRUE(logb(fp(-8)) == fp(3));
    EXPECT_TRUE(logb(fp(9)) == fp(3));
    EXPECT_TRUE(logb(fp(-0.75)) == fp(-1));
    EXPECT_THROW((void)logb(fp(0)), std::domain_error);
}

/**********************************************************************
 * Math functions against mpmath
 *
 * Every vector holds an argument and the correctly rounded result, computed with mpmath 1.4.1 at 600
 * bits. Per function and width there are the two arguments with the largest error found among 2000
 * random ones, and one or four (at I = 16 and 63) random ones. The bound is 1.2 times the largest
 * error of those 2000 plus 0.6 ulp, which covers the rounding of the reference itself, and is the
 * accuracy each function documents.
 *
 * These replace tests that compared a fixed_point128<16> result with the double function. A double
 * holds 53 of the result's 111 fraction bits: cos() was 8e-21 off for every argument at I = 16 and
 * 7026 ulp off at I = 63, and both passed.
 ***********************************************************************/
/// @brief A reference vector of a math function: arguments and the correctly rounded result, raw.
struct MathVector {
    int32_t i;
    uint64_t argHigh, argLow;
    uint64_t arg2High, arg2Low;  ///< second argument for hypot, pow and atan2, zero otherwise
    uint64_t expectedHigh, expectedLow;
};
/// @brief The largest error accepted at one I: absUlps + relUlps * max(1, |result|) + perArgUlps * |x| ulp.
struct MathBound {
    int32_t i;
    double absUlps;
    double relUlps;
    double perArgUlps;
};

/// @brief log2: arguments and the correctly rounded result, from mpmath 1.4.1.
static constexpr MathVector kLog2Vectors[] = {
    {1, 0x243DA5C09274AFA1ull, 0x3F07DCA4005B8485ull, 0ull, 0ull, 0xCB7DA33BFE0DA2ECull, 0x63EDF6F422F79D51ull},
    {1, 0x3B9BD1528BCB4085ull, 0x5362F9CBFB204EECull, 0ull, 0ull, 0xF96FD28F27B88B78ull, 0xE12C27C89CA7FECFull},
    {1, 0x4D2BFC826A7A5F16ull, 0xCAA21DB5B8165231ull, 0ull, 0ull, 0x1147B78CA3225984ull, 0x620E4812F5FF06E2ull},
    {2, 0x39A2BD5A91EC2F13ull, 0xB9E3EB3C40E4F6EFull, 0ull, 0ull, 0x1B2A1D2244256783ull, 0x3145431E1DD9E037ull},
    {2, 0x6603CA4F5B8080B8ull, 0xAA60CB66BF1416ABull, 0ull, 0ull, 0x3586394A85B3F7D1ull, 0x9926DB3BCD806BDAull},
    {2, 0x41CDE31A980C1C3Aull, 0xF5FA5A2BADD43541ull, 0ull, 0ull, 0x21489237796B8CB7ull, 0x6556A44DDDBCE836ull},
    {4, 0x000021186B7C260Cull, 0x4DB83311C1614FFFull, 0ull, 0ull, 0x906372922E043F76ull, 0xE80B4378AF17800Bull},
    {4, 0x2BB1D2AF9CABB285ull, 0x3C243F25FA360374ull, 0ull, 0ull, 0x13985686B149E87Cull, 0x292EB2A830396295ull},
    {4, 0x055F8930C8DB55F3ull, 0x5420D110BE05C622ull, 0ull, 0ull, 0xFB67FE564BF9CDD1ull, 0xE7EC17DBB8C904D4ull},
    {16, 0x0000000000000000ull, 0x66F86C8351C087EEull, 0ull, 0ull, 0xFFE7D7D1AA1426BDull, 0x17ACDA69B65C417Full},
    {16, 0x0000000000000000ull, 0x00EF378991B57596ull, 0ull, 0ull, 0xFFE4737A7ED49ECEull, 0xD3410D0E2BD17DF7ull},
    {16, 0x0000000000000019ull, 0x7F2BB476838D96AEull, 0ull, 0ull, 0xFFEAD60C071D61E3ull, 0x41F39CFF8947A863ull},
    {16, 0x0000000000000000ull, 0x00000000000002EFull, 0ull, 0ull, 0xFFCD46BDDC69AAB4ull, 0x75780E47156A9C63ull},
    {16, 0x000000000000C686ull, 0xC6D68C8F49E0DA7Bull, 0ull, 0ull, 0xFFF0510C4A63FECFull, 0xC6E94FDEDE8243B6ull},
    {16, 0x0000000000000003ull, 0xFF06C9F37D185865ull, 0ull, 0ull, 0xFFE97FD30962FD03ull, 0x73FAF27D0326B47Full},
    {32, 0x0000000000000000ull, 0x0003E728A7ABFF1Aull, 0ull, 0ull, 0xFFFFFFE97B770938ull, 0xBB6F52DC2F2E672Aull},
    {32, 0x0000000000000000ull, 0x000000260EFDCD52ull, 0ull, 0ull, 0xFFFFFFE32004E290ull, 0xFD2CA02ED9E579F5ull},
    {32, 0x000000000242BF19ull, 0x230036D1D60661C3ull, 0ull, 0ull, 0xFFFFFFFD16A0FAB0ull, 0x1B0F851DF51BC0A9ull},
    {62, 0x0000000000000000ull, 0x4811D8DD6C56844Eull, 0ull, 0ull, 0xFFFFFFFFFFFFFFFAull, 0x57B76AF748383C1Bull},
    {62, 0x0000004F6BC6C0C4ull, 0x4B1F4ACD3F80E9E0ull, 0ull, 0ull, 0x000000000000004Aull, 0x9F76327A4AA21CCEull},
    {62, 0x00798B8C6D62B507ull, 0x44820D3F1A0226B1ull, 0ull, 0ull, 0x000000000000006Bull, 0xD9C7388B2EA9D773ull},
    {63, 0x0000000000000000ull, 0x0000000028CAF7B5ull, 0ull, 0ull, 0xFFFFFFFFFFFFFFDDull, 0x59A99952E9E4156Eull},
    {63, 0x000000000000A26Dull, 0x908FE1FABCFE631Cull, 0ull, 0ull, 0x000000000000000Full, 0x57F9DE19DCAD7F21ull},
    {63, 0x0000000000008495ull, 0xE2451128502F058Eull, 0ull, 0ull, 0x000000000000000Full, 0x0CFFDA96C8056C4Aull},
    {63, 0x00000003F7221569ull, 0xC02A61559C825615ull, 0ull, 0ull, 0x0000000000000021ull, 0xFCC9BC4DB77DB092ull},
    {63, 0x0000000006332631ull, 0x430F0013CA689F13ull, 0ull, 0ull, 0x000000000000001Aull, 0xA1D94D8E78B1B066ull},
    {63, 0x0000000000000000ull, 0x00005CE6915D4813ull, 0ull, 0ull, 0xFFFFFFFFFFFFFFEEull, 0x89A14264B0BBD67Eull},
};
/// @brief log2: largest error accepted, per I (1.2 times the largest of 2000 random samples + 0.6).
static constexpr MathBound kLog2Bounds[] = {
    {1, 3.0, 0.0, 0.0},
    {2, 2.7, 0.0, 0.0},
    {4, 2.0, 0.0, 0.0},
    {16, 1.8, 0.0, 0.0},
    {32, 1.8, 0.0, 0.0},
    {62, 1.8, 0.0, 0.0},
    {63, 1.8, 0.0, 0.0},
};
/// @brief log: arguments and the correctly rounded result, from mpmath 1.4.1.
static constexpr MathVector kLogVectors[] = {
    {1, 0x24895B4495F306CAull, 0x7A54AF9FA9A18DBDull, 0ull, 0ull, 0xDC1F9F127CC92ACAull, 0x8EEFC901CB6773F3ull},
    {1, 0x2EE16A19B6C69252ull, 0x1097FD57A1B6FAE6ull, 0ull, 0ull, 0xEC13F9F810973A6Full, 0xB72F8E60A7C5C51Cull},
    {1, 0x6BA5A406CCB9E5A2ull, 0xD4D6500DC7FE98CBull, 0ull, 0ull, 0x2147437F52888A22ull, 0x6DB09ACA912DD08Eull},
    {2, 0x39A2BD5A91EC2F13ull, 0xB9E3EB3C40E4F6EFull, 0ull, 0ull, 0x12D43970BEE3CB28ull, 0x948B0361CACF4345ull},
    {2, 0x1F7B1E9779958109ull, 0x64F262CEF5E132FEull, 0ull, 0ull, 0xFF7A07AD61759EDBull, 0xABE65B16B6E26618ull},
    {2, 0x218921B2A53A51B0ull, 0x017AE7DC589665C9ull, 0ull, 0ull, 0x017FFD659AFFCCCFull, 0x4388B43FC5DD4256ull},
    {4, 0x000085AD0ED29F51ull, 0x37C9E0104DC0F4A7ull, 0ull, 0ull, 0xBDCE1336F4126C1Cull, 0x5A6CC8EAB59DE02Eull},
    {4, 0x00855B5412975A9Dull, 0x590F9C9C7064EF22ull, 0ull, 0ull, 0xEA25B38FE7AEAC8Bull, 0x047E1578638E17B9ull},
    {4, 0x0000725BE96FB976ull, 0x7124E3DB13C9ED47ull, 0ull, 0ull, 0xBC8E6E5FFB8EB645ull, 0xC7D50026323C5244ull},
    {16, 0x0000000000000000ull, 0x0000000618EAA229ull, 0ull, 0ull, 0xFFE5864CF6AE6803ull, 0x3D287DAD19EC6C5Dull},
    {16, 0x0000000000000000ull, 0x000BDBBC5A133107ull, 0ull, 0ull, 0xFFEB67025A199C57ull, 0xDC417BF1C499CC19ull},
    {16, 0x0000000000000000ull, 0x42CA808838F9E5B4ull, 0ull, 0ull, 0xFFEF0A0B6DA87F07ull, 0x27FBE71126C2C9B4ull},
    {16, 0x00360739772D5645ull, 0xB00F1236C6990F85ull, 0ull, 0ull, 0x00025761323DF92Bull, 0xD0A5A9AB21632CC6ull},
    {16, 0x0000000000226CECull, 0xBE5D2066307366ECull, 0ull, 0ull, 0xFFF7068EDB06C20Cull, 0x76C5A503D601DF57ull},
    {16, 0x0000000000000000ull, 0x00000000000000DEull, 0ull, 0ull, 0xFFDC3B4EBC7B0101ull, 0x044C70B6E2FFAA45ull},
    {32, 0x0000000000000000ull, 0x0003E728A7ABFF1Aull, 0ull, 0ull, 0xFFFFFFF064542CADull, 0xB215BE77BB346293ull},
    {32, 0x000000000266FF78ull, 0xB05A55ADA3E5B101ull, 0ull, 0ull, 0xFFFFFFFE031FE801ull, 0x892428556E859168ull},
    {32, 0x0000000000000000ull, 0x00000026BB60E159ull, 0ull, 0ull, 0xFFFFFFEBFE83A44Cull, 0xA7A11E98DD6D87A3ull},
    {62, 0x000000000075C296ull, 0xDE666F9BDF0CDA46ull, 0ull, 0ull, 0x000000000000001Eull, 0x54EB40B6A927BD68ull},
    {62, 0x0000000000000000ull, 0x0000000000D56FAFull, 0ull, 0ull, 0xFFFFFFFFFFFFFFC6ull, 0xCC5A4C094CFD0A07ull},
    {62, 0x000005538CD929FAull, 0x570DF4376E210660ull, 0ull, 0ull, 0x0000000000000039ull, 0x692B72B6676CE03Full},
    {63, 0x000017E573E9E36Cull, 0xC5D97C5E3A8EC0D9ull, 0ull, 0ull, 0x000000000000001Eull, 0xE64CE5A3F7640F04ull},
    {63, 0x0000000000000000ull, 0x00000001057DF9A1ull, 0ull, 0ull, 0xFFFFFFFFFFFFFFE9ull, 0xD72C1C2050DD46EBull},
    {63, 0x00000E43533289CEull, 0x2C2CBBACECBE97EDull, 0ull, 0ull, 0x000000000000001Eull, 0x6230ABED72F414F3ull},
    {63, 0x000014664B1ED75Eull, 0x9BC9CB0FA12881CEull, 0ull, 0ull, 0x000000000000001Eull, 0xBDCC80C62BF0E985ull},
    {63, 0x000003F41018BD2Eull, 0x473F3DE4715FA295ull, 0ull, 0ull, 0x000000000000001Dull, 0x19B777F22FD7D8CCull},
    {63, 0x0000000000000000ull, 0x000000000001377Cull, 0ull, 0ull, 0xFFFFFFFFFFFFFFDEull, 0xECD3D10DE5E90FBEull},
};
/// @brief log: largest error accepted, per I (1.2 times the largest of 2000 random samples + 0.6).
static constexpr MathBound kLogBounds[] = {
    {1, 2.8, 0.0, 0.0},
    {2, 2.6, 0.0, 0.0},
    {4, 2.2, 0.0, 0.0},
    {16, 2.0, 0.0, 0.0},
    {32, 2.0, 0.0, 0.0},
    {62, 2.0, 0.0, 0.0},
    {63, 2.0, 0.0, 0.0},
};
/// @brief log10: arguments and the correctly rounded result, from mpmath 1.4.1.
static constexpr MathVector kLog10Vectors[] = {
    {1, 0x238ADAE15A81059Dull, 0xBB4E68302D43A693ull, 0ull, 0ull, 0xEFA6FB872DAC175Aull, 0x08B5C2023A6E22E0ull},
    {1, 0x31D03DACE0022247ull, 0xA4D0C3A34215DDFEull, 0ull, 0ull, 0xF908DF81FCE8B614ull, 0xBF8AEE1E64659B41ull},
    {1, 0x5DB2E10F31B09782ull, 0x17F41778C4D4B536ull, 0ull, 0ull, 0x0A986EFF515BFF18ull, 0x8884BD28D2BA2228ull},
    {2, 0x0738229524E38AD4ull, 0x2C4A4D350122D93Bull, 0ull, 0ull, 0xEB4E979774486892ull, 0x05E16B29C3B4C15Aull},
    {2, 0x31375273CAB23117ull, 0xC1584D40D554F5F1ull, 0ull, 0ull, 0x05FB8DBCF1B2D8CEull, 0xCFD1292AD75B0D75ull},
    {2, 0x5D94B7F322106EE6ull, 0x88A8B96857805A6Eull, 0ull, 0ull, 0x0EE9C73D0DDDBE15ull, 0xE83783A23CE3829Bull},
    {4, 0x0010A9773B3A0ABBull, 0xBB3D31D057E66CC5ull, 0ull, 0ull, 0xEF487DFB577D8B6Cull, 0x9190D876D61EF7BFull},
    {4, 0x0AFA55AAC764C307ull, 0xEEF5B4F601AFEAE3ull, 0ull, 0ull, 0x011973E52A67A53Bull, 0x10FE2EE5821141D9ull},
    {4, 0x0000219345965E3Dull, 0x260DE936B805B6FAull, 0ull, 0ull, 0xDE739BEAE4ECFE21ull, 0x2B24296682928499ull},
    {16, 0x0000000000000000ull, 0x00005BEDDDE0427Bull, 0ull, 0ull, 0xFFF64B8FB98A55D1ull, 0xBC340C81DEA646ECull},
    {16, 0x0000000000000000ull, 0x013F7804F3AE6A3Eull, 0ull, 0ull, 0xFFF7C50F9B250B42ull, 0xB401BCBA9F61F023ull},
    {16, 0x000000000000028Bull, 0x2352E9AA9A39CBBCull, 0ull, 0ull, 0xFFFA5527587D40D5ull, 0xF42E43C999337BC2ull},
    {16, 0x000000000ACD109Bull, 0x78370299CFFAEA48ull, 0ull, 0ull, 0xFFFD0E0D1F4BDD6Eull, 0x611DE7B422DD57EAull},
    {16, 0x0000000000000010ull, 0xA2939E016090F1F8ull, 0ull, 0ull, 0xFFF9894B704B6AE6ull, 0x0F6B3C4FD56B54F0ull},
    {16, 0x0000000000000002ull, 0xC4BF529885F15435ull, 0ull, 0ull, 0xFFF9259C9ABE5654ull, 0xC0C1863E437FED09ull},
    {32, 0x0000000000000000ull, 0x00000000000F1685ull, 0ull, 0ull, 0xFFFFFFF4B2D9735Aull, 0xA40D98311865C6ECull},
    {32, 0x000099B2CC6FFD82ull, 0xA6678F0A27F8598Eull, 0ull, 0ull, 0x0000000272AE2466ull, 0x6D544DE529C0C374ull},
    {32, 0x0000000000000000ull, 0x00000000000002ACull, 0ull, 0ull, 0xFFFFFFF31E5CC875ull, 0x3EBD3B216C6D6046ull},
    {62, 0x00000091788427E8ull, 0x7912B76B0544C390ull, 0ull, 0ull, 0x0000000000000016ull, 0xFD4A0D2D266CFD78ull},
    {62, 0x0000000000000000ull, 0x0000003870FE2F65ull, 0ull, 0ull, 0xFFFFFFFFFFFFFFEFull, 0xA29D7BD7EC42FE95ull},
    {62, 0x00000000000025F9ull, 0xDA93374D3E57901Dull, 0ull, 0ull, 0x0000000000000007ull, 0x5F99A50641421E0Bull},
    {63, 0x0000002FBDE2026Full, 0x4EC4CC70D6CD164Cull, 0ull, 0ull, 0x000000000000000Bull, 0x4FD5ED3CAD7BB5FAull},
    {63, 0x0000000000000067ull, 0xF70D306E9BCAD558ull, 0ull, 0ull, 0x0000000000000002ull, 0x0452BAA64C260582ull},
    {63, 0x000006102AB7DBEFull, 0xF394EF60143C4A53ull, 0ull, 0ull, 0x000000000000000Cull, 0xD2EAFF3D7E2D7B12ull},
    {63, 0x00096F2C64362B6Full, 0x227B6EA8A06EB764ull, 0ull, 0ull, 0x000000000000000Full, 0x6C94F907C7BF293Full},
    {63, 0x0000000000000000ull, 0x00000C515858365Bull, 0ull, 0ull, 0xFFFFFFFFFFFFFFF9ull, 0xDDA5D4E83DD66E7Cull},
    {63, 0x0000000000000000ull, 0x0000000000000002ull, 0ull, 0ull, 0xFFFFFFFFFFFFFFEDull, 0x08FCFCA303E908BAull},
};
/// @brief log10: largest error accepted, per I (1.2 times the largest of 2000 random samples + 0.6).
static constexpr MathBound kLog10Bounds[] = {
    {1, 1.8, 0.0, 0.0},
    {2, 1.8, 0.0, 0.0},
    {4, 1.6, 0.0, 0.0},
    {16, 1.6, 0.0, 0.0},
    {32, 1.5, 0.0, 0.0},
    {62, 1.5, 0.0, 0.0},
    {63, 1.5, 0.0, 0.0},
};
/// @brief log1p: arguments and the correctly rounded result, from mpmath 1.4.1.
static constexpr MathVector kLog1pVectors[] = {
    {2, 0xE88370B75D072927ull, 0xCD85CF66A20B2852ull, 0ull, 0ull, 0xD5A10E1304DDFAF5ull, 0xE1232676267E538Aull},
    {2, 0xE92CCA340F287D84ull, 0xA40D2E8A036B4B6Bull, 0ull, 0ull, 0xD8061472EEC093EAull, 0x8E4984C8750FD095ull},
    {2, 0x0000000000000000ull, 0x0000000000000034ull, 0ull, 0ull, 0x0000000000000000ull, 0x0000000000000034ull},
    {4, 0x005B2A94ADCF5A4Aull, 0x247CF1B764E89A39ull, 0ull, 0ull, 0x0059320B0F1FBF97ull, 0x70F946157729C0A9ull},
    {4, 0xFFFFFFFFFFFFFFFFull, 0xFFFFFFFFFFFFFBA8ull, 0ull, 0ull, 0xFFFFFFFFFFFFFFFFull, 0xFFFFFFFFFFFFFBA8ull},
    {4, 0x000000001DEB6944ull, 0x2866C9B2B255264Dull, 0ull, 0ull, 0x000000001DEB6943ull, 0xF073E47544AD49F7ull},
    {16, 0xFFFFFFFFFFFFFE51ull, 0x1D64762FDB4D67C8ull, 0ull, 0ull, 0xFFFFFFFFFFFFFE51ull, 0x1D64762D060F5CAEull},
    {16, 0xFFFFFFFFFFFFFFFFull, 0xFFFFFFFFFFE498B6ull, 0ull, 0ull, 0xFFFFFFFFFFFFFFFFull, 0xFFFFFFFFFFE498B6ull},
    {16, 0xFFFFFFFFFFFFFFFFull, 0x631A9811B85D87C4ull, 0ull, 0ull, 0xFFFFFFFFFFFFFFFFull, 0x631A9811B85D279Cull},
    {16, 0x0000000000000000ull, 0x000000C94BF788ACull, 0ull, 0ull, 0x0000000000000000ull, 0x000000C94BF788ACull},
    {16, 0x0000000000000000ull, 0x00009C5060A890E7ull, 0ull, 0ull, 0x0000000000000000ull, 0x00009C5060A890E7ull},
    {16, 0xFFFFFFFFFFFFFEB0ull, 0x624568C7D98289BAull, 0ull, 0ull, 0xFFFFFFFFFFFFFEB0ull, 0x624568C621845A34ull},
    {32, 0xFFFFFFFFFFFFFFFFull, 0xFFFFE20A0E532CCDull, 0ull, 0ull, 0xFFFFFFFFFFFFFFFFull, 0xFFFFE20A0E532CCDull},
    {32, 0xFFFFFFFFFD521A34ull, 0x3FF71FCC81D3189Full, 0ull, 0ull, 0xFFFFFFFFFD4AD26Eull, 0xD37CD97DB69DEE0Bull},
    {32, 0x0000000078E806FFull, 0x8691B387A57DD927ull, 0ull, 0ull, 0x0000000055203E42ull, 0xE9AAA082F5574998ull},
    {62, 0x0000000000000000ull, 0x5524ACFDDD9E5BB4ull, 0ull, 0ull, 0x0000000000000000ull, 0x4EC31FF0D1F76AFFull},
    {62, 0x0000000000000000ull, 0x69CF5E35B457CCC0ull, 0ull, 0ull, 0x0000000000000000ull, 0x602EB68B7955799Cull},
    {62, 0x0000000000000000ull, 0x001DFD1ACB3E1229ull, 0ull, 0ull, 0x0000000000000000ull, 0x001DFC39FF712319ull},
    {63, 0xFFFFFFFFFFFFFFFFull, 0xA2321C4808DE9328ull, 0ull, 0ull, 0xFFFFFFFFFFFFFFFFull, 0x8B2B1A81E23DA615ull},
    {63, 0xFFFFFFFFFFFFFFFFull, 0xFFF55A2A4884836Eull, 0ull, 0ull, 0xFFFFFFFFFFFFFFFFull, 0xFFF559F196E24D4Eull},
    {63, 0x0000000000000000ull, 0x0000000000000267ull, 0ull, 0ull, 0x0000000000000000ull, 0x0000000000000267ull},
    {63, 0x0000000000000000ull, 0x000000001AB6CD1Aull, 0ull, 0ull, 0x0000000000000000ull, 0x000000001AB6CD1Aull},
    {63, 0x0000000000000000ull, 0x000000000000A3CFull, 0ull, 0ull, 0x0000000000000000ull, 0x000000000000A3CFull},
    {63, 0xFFFFFFFFFFFFFFFFull, 0xFFFFFFFFF40F40E6ull, 0ull, 0ull, 0xFFFFFFFFFFFFFFFFull, 0xFFFFFFFFF40F40E6ull},
};
/// @brief log1p: largest error accepted, per I (1.2 times the largest of 2000 random samples + 0.6).
static constexpr MathBound kLog1pBounds[] = {
    {2, 2.2, 0.0, 0.0},
    {4, 1.9, 0.0, 0.0},
    {16, 1.9, 0.0, 0.0},
    {32, 1.8, 0.0, 0.0},
    {62, 2.0, 0.0, 0.0},
    {63, 2.0, 0.0, 0.0},
};
/// @brief exp: arguments and the correctly rounded result, from mpmath 1.4.1.
static constexpr MathVector kExpVectors[] = {
    {2, 0xF1D76BC8D999C6D5ull, 0xBEAF094F8B98EF47ull, 0ull, 0ull, 0x148F048C0658E89Dull, 0x50694CE1BDA2CAD6ull},
    {2, 0xF847BA0EF6B7A44Dull, 0x5CADD92F6D3CB05Dull, 0ull, 0ull, 0x19240ACB209235B3ull, 0xE1D5A6B4866F5CE5ull},
    {2, 0xFFF5B462AD0B927Cull, 0x1A21B81D0394866Dull, 0ull, 0ull, 0x1FF5B60A7A17597Eull, 0xF1518EFCA49F688Cull},
    {4, 0x07B0FFC64C947B29ull, 0x177244F3DD401BC8ull, 0ull, 0ull, 0x14EC62380C76936Aull, 0x396F729B82919B80ull},
    {4, 0x079173D7953498DAull, 0xD610A37347E7B9C9ull, 0ull, 0ull, 0x149A8218573072F2ull, 0x320B1FED40686B61ull},
    {4, 0xFFFFFFFFFFFFFFFFull, 0xFFFFFB431C80033Cull, 0ull, 0ull, 0x07FFFFFFFFFFFFFFull, 0xFFFFFB431C80033Cull},
    {16, 0xFFFFBF44F1712D1Eull, 0x60D3FE192FFB0488ull, 0ull, 0ull, 0x00004D31AA05B516ull, 0xCBE3D3576319CB84ull},
    {16, 0xFFFFEE2CFE12B50Aull, 0x5EF7A5429D32EDCBull, 0ull, 0ull, 0x00006F5C72E14E6Eull, 0x6DC46E7ABCA76D22ull},
    {16, 0xFFFFFFFF80E35619ull, 0x17140E6C60746CBDull, 0ull, 0ull, 0x00007FFF80E39536ull, 0x75F85648E02914F6ull},
    {16, 0x0000000000000000ull, 0x00000000000078CEull, 0ull, 0ull, 0x0000800000000000ull, 0x00000000000078CEull},
    {16, 0x0000000000000000ull, 0x000000000507784Cull, 0ull, 0ull, 0x0000800000000000ull, 0x000000000507784Cull},
    {16, 0x0000000000000000ull, 0x0000000000000068ull, 0ull, 0ull, 0x0000800000000000ull, 0x0000000000000068ull},
    {32, 0xFFFFFFFFA7D5E342ull, 0x37F309BD71D58654ull, 0ull, 0ull, 0x0000000040479F95ull, 0x132A49D507E07FB8ull},
    {32, 0xFFFFFFFFD7E503B6ull, 0x8015EC74A0FD7ECBull, 0ull, 0ull, 0x000000005D91D896ull, 0x718F9E9AA0907A11ull},
    {32, 0x0000000000000000ull, 0x00000000000000F1ull, 0ull, 0ull, 0x0000000080000000ull, 0x00000000000000F1ull},
    {62, 0x000000000000003Full, 0xF46FA49E4612ED1Eull, 0ull, 0ull, 0x00008C6CD322F783ull, 0x95DE425D8F3D32DCull},
    {62, 0x000000000000003Bull, 0xF8B6F5FB1DE8198Bull, 0ull, 0ull, 0x00001329F889ACB3ull, 0xF2F927C55ECD8534ull},
    {62, 0x0000000000000000ull, 0x0000004032AA0834ull, 0ull, 0ull, 0x0000000000000002ull, 0x0000004032AA0C3Aull},
    {63, 0x0000000000000012ull, 0x2D7F3026BB8AA85Eull, 0ull, 0ull, 0x0000000004ACC0B1ull, 0xF5BE3DA98E78E51Bull},
    {63, 0x000000000000000Full, 0xA35135AC71E6F4D0ull, 0ull, 0ull, 0x00000000005E67E5ull, 0xAA5B59058ECD1041ull},
    {63, 0xFFFFFFFFFFFFFFFFull, 0xFC039A61D1F67F6Full, 0ull, 0ull, 0x0000000000000000ull, 0xFC0B817B52E8D882ull},
    {63, 0xFFFFFFFFFFFFFFFFull, 0xFFFFFFFFFFFFFCCFull, 0ull, 0ull, 0x0000000000000000ull, 0xFFFFFFFFFFFFFCCFull},
    {63, 0xFFFFFFFFFFFFFFFFull, 0xFFFFFFFFFFFF5DDAull, 0ull, 0ull, 0x0000000000000000ull, 0xFFFFFFFFFFFF5DDAull},
    {63, 0xFFFFFFFFFFFFFFFFull, 0xFFFFFF823A144A4Aull, 0ull, 0ull, 0x0000000000000000ull, 0xFFFFFF823A14692Full},
};
/// @brief exp: largest error accepted, per I (1.2 times the largest of 2000 random samples + 0.6).
static constexpr MathBound kExpBounds[] = {
    {2, 0.6, 6.9, 0.0},
    {4, 0.6, 7.5, 0.0},
    {16, 0.6, 5.0, 0.0},
    {32, 0.6, 5.5, 0.0},
    {62, 0.6, 8.0, 0.0},
    {63, 0.6, 4.1, 0.0},
};
/// @brief exp2: arguments and the correctly rounded result, from mpmath 1.4.1.
static constexpr MathVector kExp2Vectors[] = {
    {2, 0xE069AEAE50A170ACull, 0x6B8263B93E26C2ADull, 0ull, 0ull, 0x1024CA7C4A5C4BCAull, 0x5D0A96A74AE40800ull},
    {2, 0x1CFCB6CFD15254F4ull, 0x3226771ADDD9B060ull, 0ull, 0ull, 0x3BF4E9C9D1F7CBC5ull, 0x6CDAB6BB28DB686Aull},
    {2, 0x0000000000000000ull, 0x000001098EBA6791ull, 0ull, 0ull, 0x2000000000000000ull, 0x000000B812094D82ull},
    {4, 0xFCE35982D5340A36ull, 0x5F0A99AD961C4262ull, 0ull, 0ull, 0x061BFC0D633DB791ull, 0x853013834090D90Cull},
    {4, 0x0CF825554267E8B6ull, 0xC5B7C494418C8162ull, 0ull, 0ull, 0x189C203EBE7BA7DCull, 0xFB21F7C228EF8003ull},
    {4, 0xFFFFFFFFE816326Dull, 0x493763066A309FACull, 0ull, 0ull, 0x07FFFFFFEF6CB082ull, 0xEDCF0827F8BE0FBAull},
    {16, 0xFFFF926FC3758D5Aull, 0x9722EE1ABA2CBA35ull, 0ull, 0ull, 0x000046B82FCF4E3Aull, 0x923B33A13FB7C484ull},
    {16, 0x000233651DA670BEull, 0xA8D5FB7159601EB7ull, 0ull, 0ull, 0x000A9135020B67FEull, 0x5985F95C29BFBCBCull},
    {16, 0xFFFFFFFFFFFFFFFFull, 0xFFFFFFAEE8CF3B7Eull, 0ull, 0ull, 0x00007FFFFFFFFFFFull, 0xFFFFFFC7CAD35CD5ull},
    {16, 0x0000000000000000ull, 0x00000000000000C7ull, 0ull, 0ull, 0x0000800000000000ull, 0x000000000000008Aull},
    {16, 0x00000000762CE3EBull, 0x6D156E5364873C05ull, 0ull, 0ull, 0x0000800051E9CED9ull, 0x7328134358A52447ull},
    {16, 0xFFFDB34DA73932CBull, 0x1759D0281C276527ull, 0ull, 0ull, 0x00000547EEB10D70ull, 0x59C137ABB878CA4Cull},
    {32, 0xFFFFFFFFDBB6D620ull, 0x2331C6ECAE4FB291ull, 0ull, 0ull, 0x00000000692A5BFCull, 0xFAB377DEA3E850C2ull},
    {32, 0xFFFFFFFFBDEB87C5ull, 0xFFF393A0CCC3FC9Aull, 0ull, 0ull, 0x00000000597EF5DDull, 0x9D73A8D80CA93CF0ull},
    {32, 0xFFFFFFFFFFFFEDA6ull, 0x9C3040D1188D7B8Dull, 0ull, 0ull, 0x000000007FFFF348ull, 0x0518E2D5AC49BB7Cull},
    {62, 0x000000000000003Full, 0xDAB6861B0F2FEB29ull, 0ull, 0ull, 0x00000001E6CBCD23ull, 0xC897818912A6B36Aull},
    {62, 0x000000000000003Cull, 0x5BF35F32737E29ACull, 0ull, 0ull, 0x0000000090F7D20Full, 0xB5F49815546D08D9ull},
    {62, 0x0000000000000000ull, 0x00000000004BBDD3ull, 0ull, 0ull, 0x0000000000000002ull, 0x0000000000348001ull},
    {63, 0x0000000000000014ull, 0x9A04BC0CD347E083ull, 0ull, 0ull, 0x0000000000184769ull, 0xF31726FF655EAE42ull},
    {63, 0x0000000000000006ull, 0x0D49B03D067A4A5Full, 0ull, 0ull, 0x0000000000000042ull, 0x583334655200A640ull},
    {63, 0x0000000000000017ull, 0xF235DD995FF82C5Dull, 0ull, 0ull, 0x0000000000F69E36ull, 0xBADA12D50DF08DC8ull},
    {63, 0xFFFFFFFFFFFFFFFFull, 0xE9DC4905CE42A54Cull, 0ull, 0ull, 0x0000000000000000ull, 0xF11ADED57A50817Eull},
    {63, 0x0000000000000000ull, 0x00032A9CE54B5D9Dull, 0ull, 0ull, 0x0000000000000001ull, 0x000231E228E8FEECull},
    {63, 0x0000000000000000ull, 0x0001DE0DFEA7038Full, 0ull, 0ull, 0x0000000000000001ull, 0x00014B5D9283BD7Full},
};
/// @brief exp2: largest error accepted, per I (1.2 times the largest of 2000 random samples + 0.6).
static constexpr MathBound kExp2Bounds[] = {
    {2, 0.6, 4.6, 0.0},
    {4, 0.6, 5.6, 0.0},
    {16, 0.6, 5.0, 0.0},
    {32, 0.6, 4.9, 0.0},
    {62, 0.6, 7.0, 0.0},
    {63, 0.6, 4.3, 0.0},
};
/// @brief expm1: arguments and the correctly rounded result, from mpmath 1.4.1.
static constexpr MathVector kExpm1Vectors[] = {
    {2, 0xF1D76BC8D999C6D5ull, 0xBEAF094F8B98EF47ull, 0ull, 0ull, 0xF48F048C0658E89Dull, 0x50694CE1BDA2CAD6ull},
    {2, 0xF847BA0EF6B7A44Dull, 0x5CADD92F6D3CB05Dull, 0ull, 0ull, 0xF9240ACB209235B3ull, 0xE1D5A6B4866F5CE5ull},
    {2, 0xFFFFFFFFFFFFFFFFull, 0xFFFFFFFF925AEBF4ull, 0ull, 0ull, 0xFFFFFFFFFFFFFFFFull, 0xFFFFFFFF925AEBF4ull},
    {4, 0x07B0FFC64C947B29ull, 0x177244F3DD401BC8ull, 0ull, 0ull, 0x0CEC62380C76936Aull, 0x396F729B82919B80ull},
    {4, 0x079173D7953498DAull, 0xD610A37347E7B9C9ull, 0ull, 0ull, 0x0C9A8218573072F2ull, 0x320B1FED40686B61ull},
    {4, 0x000C2B4255F10988ull, 0x29FDCC049126A175ull, 0ull, 0ull, 0x000C3488618BDD3Full, 0xC0C190A87DB86776ull},
    {16, 0xFFFFBF44F1712D1Eull, 0x60D3FE192FFB0488ull, 0ull, 0ull, 0xFFFFCD31AA05B516ull, 0xCBE3D3576319CB84ull},
    {16, 0xFFFFEE2CFE12B50Aull, 0x5EF7A5429D32EDCBull, 0ull, 0ull, 0xFFFFEF5C72E14E6Eull, 0x6DC46E7ABCA76D22ull},
    {16, 0x0000000000000000ull, 0x0000AC68A83CCE9Cull, 0ull, 0ull, 0x0000000000000000ull, 0x0000AC68A83CCE9Cull},
    {16, 0x0000000000001169ull, 0x61F21D968F41BB0Bull, 0ull, 0ull, 0x0000000000001169ull, 0x61F21EC5B9A5611Cull},
    {16, 0xFFFFFFFFFFFFFFFFull, 0xFFFEE190281A65D8ull, 0ull, 0ull, 0xFFFFFFFFFFFFFFFFull, 0xFFFEE190281A65D8ull},
    {16, 0x0002202D08835D50ull, 0x6D46E9283090743Bull, 0ull, 0ull, 0x002299D565A439ABull, 0xD32C1A71C9C78A91ull},
    {32, 0xFFFFFFFFA7D5E342ull, 0x37F309BD71D58654ull, 0ull, 0ull, 0xFFFFFFFFC0479F95ull, 0x132A49D507E07FB8ull},
    {32, 0xFFFFFFFFD7E503B6ull, 0x8015EC74A0FD7ECBull, 0ull, 0ull, 0xFFFFFFFFDD91D896ull, 0x718F9E9AA0907A11ull},
    {32, 0x0000000000000000ull, 0x000921E7293C47A7ull, 0ull, 0ull, 0x0000000000000000ull, 0x000921E7293C47FAull},
    {62, 0x000000000000003Full, 0xF46FA49E4612ED1Eull, 0ull, 0ull, 0x00008C6CD322F781ull, 0x95DE425D8F3D32DCull},
    {62, 0x000000000000003Bull, 0xF8B6F5FB1DE8198Bull, 0ull, 0ull, 0x00001329F889ACB1ull, 0xF2F927C55ECD8534ull},
    {62, 0xFFFFFFFFFFFFFFFFull, 0xFFFFFFFFC913C37Eull, 0ull, 0ull, 0xFFFFFFFFFFFFFFFFull, 0xFFFFFFFFC913C37Eull},
    {63, 0x0000000000000012ull, 0x2D7F3026BB8AA85Eull, 0ull, 0ull, 0x0000000004ACC0B0ull, 0xF5BE3DA98E78E51Bull},
    {63, 0x0000000000000000ull, 0x78E4B8116E17C57Dull, 0ull, 0ull, 0x0000000000000000ull, 0x9A845C3B11495368ull},
    {63, 0xFFFFFFFFFFFFFFFFull, 0xFFFFF01295E12BACull, 0ull, 0ull, 0xFFFFFFFFFFFFFFFFull, 0xFFFFF012966002F8ull},
    {63, 0x0000000000000000ull, 0x000000000169679Full, 0ull, 0ull, 0x0000000000000000ull, 0x000000000169679Full},
    {63, 0xFFFFFFFFFFFFFFFFull, 0xF6F63FC094AC88EEull, 0ull, 0ull, 0xFFFFFFFFFFFFFFFFull, 0xF71E9DB9D802C150ull},
    {63, 0x0000000000000000ull, 0x000001B4B2DA5382ull, 0ull, 0ull, 0x0000000000000000ull, 0x000001B4B2DBC7FBull},
};
/// @brief expm1: largest error accepted, per I (1.2 times the largest of 2000 random samples + 0.6).
static constexpr MathBound kExpm1Bounds[] = {
    {2, 0.6, 6.9, 0.0},
    {4, 0.6, 11.8, 0.0},
    {16, 0.6, 5.0, 0.0},
    {32, 0.6, 5.5, 0.0},
    {62, 0.6, 8.0, 0.0},
    {63, 0.6, 4.1, 0.0},
};
/// @brief pow: arguments and the correctly rounded result, from mpmath 1.4.1.
static constexpr MathVector kPowVectors[] = {
    {2, 0x3DACD8754E641599ull, 0xED22AA4735D16805ull, 0xF82D897A6B804C25ull, 0x112BF7388CFEA897ull, 0x1B420BDE397483F9ull, 0xD5E2467B884B36B4ull},
    {2, 0x165205DEE780DA9Cull, 0xA46B2EE1769F7D5Full, 0x0051F4C7D2A618E8ull, 0x7F3DEFAA4B90CCADull, 0x1FE28798C05C1277ull, 0x963AF53EA3D2556Eull},
    {2, 0x1774454939716B43ull, 0x6D9129B2A87A8A71ull, 0xFFFFFFFFFFFFFFFFull, 0xFFFFFFFFC42BD1DFull, 0x2000000000000000ull, 0x0000000012968A0Dull},
    {4, 0x02285E779AA95FE5ull, 0xA4683321E3C17574ull, 0xF40AF77941F76168ull, 0x2F5D08BB6118FF66ull, 0x38B6E46B0F419536ull, 0x3F604897244072BBull},
    {4, 0x0085DE967AA77F5Dull, 0x2717E27512B7D636ull, 0xFD2A702FEF518C91ull, 0x2367DC2D9F899E66ull, 0x150702213BE4B5E8ull, 0x1FA1D58084777C9Full},
    {4, 0x01F2ED7CE196A220ull, 0x78E729FE12D74BA3ull, 0x0000000000000000ull, 0x00000000000094FEull, 0x07FFFFFFFFFFFFFFull, 0xFFFFFFFFFFFF2D99ull},
    {16, 0x00013BE6D4500F91ull, 0x230126CA0FA94B98ull, 0xFFFFE3AF6DD3C67Dull, 0x6185E2C14B5E8CA0ull, 0x000068D074C53689ull, 0x277FDA1F73AD12CAull},
    {16, 0x00002B24711EC26Bull, 0x43844731B20818F0ull, 0x00000C0A40C1776Full, 0x1496B73E4369AAAEull, 0x0000738DB8CED524ull, 0x8AD0E8C602BD0D1Aull},
    {16, 0x0001F7DBFDB59FD3ull, 0x176873E6908967E9ull, 0xFFFFFFFFFFFFFFFFull, 0xFFFFFFFFFFFFFFFDull, 0x00007FFFFFFFFFFFull, 0xFFFFFFFFFFFFFFFCull},
    {16, 0x0000534A67694E47ull, 0x1F866B45BFEA10E9ull, 0x00000000000003E0ull, 0x57877C7A8185C326ull, 0x00007FFFFFFFFE55ull, 0x98B7CE34016E3B93ull},
    {16, 0x00000B36BD315221ull, 0xE4B85A56E743326Cull, 0xFFFFFFFFFFFFFFFFull, 0xFFFFFFFC51EC112Eull, 0x0000800000000000ull, 0x00000008F5D96C2Bull},
    {16, 0x00000B6773F30D35ull, 0xDB8E843656702A66ull, 0x0000000000000000ull, 0x01F7E98B82AE2A87ull, 0x00007FFFFFFFFFFFull, 0xFB3D830433E70B3Bull},
    {32, 0x00000000318EC385ull, 0x5E5B941D756CFD25ull, 0x000000001F3ECF0Full, 0x9E7EEDACD1FBAE10ull, 0x000000006588E732ull, 0x1BEB0D8C5A5512FCull},
    {32, 0x00000000288CB19Eull, 0x36E131283FDC02B9ull, 0x0000000034982230ull, 0xB689155D4468892Eull, 0x000000004FD08D20ull, 0x5D4B0C20A5A73F33ull},
    {32, 0x000000000A6C9C67ull, 0xC42052ACC33F8DD0ull, 0x000000000638C75Cull, 0x0242ED5ED1F79202ull, 0x00000000714F5A37ull, 0x1803DFC9C7D84492ull},
    {62, 0x0000000000000003ull, 0x81CB707FED78E4F3ull, 0xFFFFFFFFFFFFFFFEull, 0x9AB567C622717E4Eull, 0x0000000000000001ull, 0x59FD5B2B27D940CAull},
    {62, 0x0000000000000003ull, 0x5ABE901EDC3EBC2Aull, 0xFFFFFFFFFFFFFFFFull, 0xC2376DD981D9F07Aull, 0x0000000000000001ull, 0xE106811705DC6B37ull},
    {62, 0x0000000000000007ull, 0x247227D7F1BF9476ull, 0xFFFFFFFFFFFFFFFFull, 0x9B7C708B452254B9ull, 0x0000000000000001ull, 0x8ECA2BACB78ADB7Dull},
    {63, 0x0000000000000002ull, 0x371E691CCE8A9F27ull, 0xFFFFFFFFFFFFFFFEull, 0xCF7160709C6F80CCull, 0x0000000000000000ull, 0x63607999384E890Eull},
    {63, 0x0000000000000001ull, 0xC6E28ACDD71737ADull, 0xFFFFFFFFFFFFFFFFull, 0x7E8155142C9DDF34ull, 0x0000000000000000ull, 0xBF677778EFA03C6Eull},
    {63, 0x0000000000000000ull, 0xEE6533177E2E64A1ull, 0xFFFFFFFFFFFFFFFFull, 0xFFFFFFB075ECBA6Bull, 0x0000000000000001ull, 0x00000005AABD8C09ull},
    {63, 0x0000000000000000ull, 0x1717AA39DA3CB40Cull, 0x0000000000000000ull, 0x00000DFC387027A7ull, 0x0000000000000000ull, 0xFFFFDE5B2C124E47ull},
    {63, 0x0000000000000003ull, 0x32638E45A188639Aull, 0x0000000000000000ull, 0x000000000036C6D9ull, 0x0000000000000001ull, 0x00000000003FA8CCull},
    {63, 0x0000000000000000ull, 0x12D6285FCAE90CBEull, 0xFFFFFFFFFFFFFFFFull, 0xFFFFFFFFFFFFC068ull, 0x0000000000000001ull, 0x000000000000A5F1ull},
};
/// @brief pow: largest error accepted, per I (1.2 times the largest of 2000 random samples + 0.6).
static constexpr MathBound kPowBounds[] = {
    {2, 0.6, 4.1, 0.0},
    {4, 0.6, 9.4, 0.0},
    {16, 0.6, 5.8, 0.0},
    {32, 0.6, 4.6, 0.0},
    {62, 0.6, 4.6, 0.0},
    {63, 0.6, 4.6, 0.0},
};
/// @brief sin: arguments and the correctly rounded result, from mpmath 1.4.1.
static constexpr MathVector kSinVectors[] = {
    {4, 0x1CEF2B21DBABBD8Full, 0x67D62844BA1E40F7ull, 0ull, 0ull, 0xFC5706E6D92804DAull, 0xC2D0CA6EF8107907ull},
    {4, 0xE949F9E8055BC998ull, 0x860CE7B25FF0104Dull, 0ull, 0ull, 0xFD9D779662B455F4ull, 0x84F37E65FEC5EEF6ull},
    {4, 0x0000000000000000ull, 0x004E1376F7F807F1ull, 0ull, 0ull, 0x0000000000000000ull, 0x004E1376F7F807F1ull},
    {16, 0x0001655995B17360ull, 0x9801BA883E8D3414ull, 0ull, 0ull, 0x00002BDDCDA15E99ull, 0xF7EA25343218DBEAull},
    {16, 0xFFFE67D85EE0981Full, 0x3AF20281FE105137ull, 0ull, 0ull, 0x0000060759A204EBull, 0x771B196CE728B9B4ull},
    {16, 0x0000000000000000ull, 0x0000000000000057ull, 0ull, 0ull, 0x0000000000000000ull, 0x0000000000000057ull},
    {16, 0x00000000004EECA1ull, 0x5422E0CEA1AEB5A4ull, 0ull, 0ull, 0x00000000004EECA1ull, 0x5422E0C9A1699174ull},
    {16, 0xFFFFFFFFFFFFFFFFull, 0xFFFFFFFFFFFD8D96ull, 0ull, 0ull, 0xFFFFFFFFFFFFFFFFull, 0xFFFFFFFFFFFD8D96ull},
    {16, 0xFFFFFFFFFFF537BEull, 0x746FBD0CD5F4687Full, 0ull, 0ull, 0xFFFFFFFFFFF537BEull, 0x746FBD0CD9381512ull},
    {32, 0xFFFFFFFE4642F4D8ull, 0x48AD93F29C32F5B8ull, 0ull, 0ull, 0x0000000026FC36B8ull, 0x1EB314D347C66758ull},
    {32, 0xFFFFFFFE2AF5C6F4ull, 0x7B3F0A9F6B3D7CDFull, 0ull, 0ull, 0x000000003FE8CFA4ull, 0x8E7913943228629Eull},
    {32, 0x0000000000000000ull, 0x000000006690FD82ull, 0ull, 0ull, 0x0000000000000000ull, 0x000000006690FD82ull},
    {62, 0x0000000000000006ull, 0x68CA0DDCC284470Cull, 0ull, 0ull, 0xFFFFFFFFFFFFFFFFull, 0xDFBA42373EC3A7D1ull},
    {62, 0xFFFFFFFFFFFFFFF8ull, 0x9D7C4398D1D700F9ull, 0ull, 0ull, 0x0000000000000001ull, 0x0BF91DD316046612ull},
    {62, 0xFFFFFFFFFFFFFFFFull, 0xFFD67429515E5DA5ull, 0ull, 0ull, 0xFFFFFFFFFFFFFFFFull, 0xFFD674295D0A63C7ull},
    {63, 0xFFFFFFFFFFFFFFFCull, 0xE4376899651BAE25ull, 0ull, 0ull, 0xFFFFFFFFFFFFFFFFull, 0xF78991ED465A6A0Eull},
    {63, 0xFFFFFFFFFFFFFFFDull, 0x1AFAF9E78484A0ECull, 0ull, 0ull, 0xFFFFFFFFFFFFFFFFull, 0xC169AC0C62F06D8Full},
    {63, 0x0000000000000000ull, 0x00004B2CB4EDD95Dull, 0ull, 0ull, 0x0000000000000000ull, 0x00004B2CB4EDD848ull},
    {63, 0x0000000000000000ull, 0x000001DB502A49E4ull, 0ull, 0ull, 0x0000000000000000ull, 0x000001DB502A49E4ull},
    {63, 0xFFFFFFFFFFFFFFFFull, 0xFFFFFFFFFFFFFFFEull, 0ull, 0ull, 0xFFFFFFFFFFFFFFFFull, 0xFFFFFFFFFFFFFFFEull},
    {63, 0xFFFFFFFFFFFFFFFFull, 0xFFFF4F8794AC3C32ull, 0ull, 0ull, 0xFFFFFFFFFFFFFFFFull, 0xFFFF4F8794AC4A2Cull},
};
/// @brief sin: largest error accepted, per I (1.2 times the largest of 2000 random samples + 0.6).
static constexpr MathBound kSinBounds[] = {
    {4, 3.6, 0.0, 0.0},
    {16, 1.5, 0.0, 0.0},
    {32, 2.2, 0.0, 0.0},
    {62, 1.5, 0.0, 0.0},
    {63, 2.3, 0.0, 0.0},
};
/// @brief cos: arguments and the correctly rounded result, from mpmath 1.4.1.
static constexpr MathVector kCosVectors[] = {
    {4, 0x1DFA7B441DE48C45ull, 0x82CF7A284A3332E8ull, 0ull, 0ull, 0xF96C58ACA506F303ull, 0xCB366FDE773B0291ull},
    {4, 0x0A85E000B28EA194ull, 0x5BF136F76418A72Aull, 0ull, 0ull, 0x0205723330CCD313ull, 0x105528A9CDCE2F8Dull},
    {4, 0xFFFFFFFFFFFFFFFEull, 0x9CA8F6CBD51DD433ull, 0ull, 0ull, 0x07FFFFFFFFFFFFFFull, 0xFFFFFFFFFFFFFFE1ull},
    {16, 0x0001E9FCDED10173ull, 0x62B7B13AA5FA7077ull, 0ull, 0ull, 0xFFFF9CFDA68A967Eull, 0x1965FA35876B5F40ull},
    {16, 0xFFFFF9C3B8C4A60Cull, 0xC46432276F87714Dull, 0ull, 0ull, 0x00007FD9206F9EFCull, 0x257BE4F15D45BABAull},
    {16, 0xFFFFFFFFFFFFFFFFull, 0xFFFF8A3BE651A9AAull, 0ull, 0ull, 0x0000800000000000ull, 0x0000000000000000ull},
    {16, 0xFFFFFFFFFFFF4E8Cull, 0x1CA240284FC02C7Full, 0ull, 0ull, 0x00007FFFFFFFFFFFull, 0xFFFF84FE8B224C8Eull},
    {16, 0xFFFFFFFFFFFFFFFFull, 0xFFFFFFFFFFFC00F1ull, 0ull, 0ull, 0x0000800000000000ull, 0x0000000000000000ull},
    {16, 0xFFFFFFFFFFFFFFFFull, 0xF8FB6E8A7B7348AEull, 0ull, 0ull, 0x00007FFFFFFFFFFFull, 0xFFFFFFFFFFFFFFCFull},
    {32, 0xFFFFFFFF52DBF3FFull, 0xCFEFDAA9276F868Full, 0ull, 0ull, 0x000000001BB341C8ull, 0xB7E331A31B8E3EB4ull},
    {32, 0x00000001DE9C06D9ull, 0x6964A38C87167C3Aull, 0ull, 0ull, 0xFFFFFFFF962E0A80ull, 0x2515A318DC22596Cull},
    {32, 0xFFFFFFFFFFFFFFFFull, 0xFFF08A79754233E7ull, 0ull, 0ull, 0x000000007FFFFFFFull, 0xFFFFFFFFFFFFFF11ull},
    {62, 0x0000000000000000ull, 0x000223B7E52D88DBull, 0ull, 0ull, 0x0000000000000001ull, 0xFFFFFFFEDB0927A2ull},
    {62, 0xFFFFFFFFFFFFFFFFull, 0xCC3CD21E2C55C610ull, 0ull, 0ull, 0x0000000000000001ull, 0xFD62BBBBABA7F23Dull},
    {62, 0x0000000000000000ull, 0x000000000000001Cull, 0ull, 0ull, 0x0000000000000002ull, 0x0000000000000000ull},
    {63, 0x0000000000000003ull, 0xF8F404BC273FC04Cull, 0ull, 0ull, 0xFFFFFFFFFFFFFFFFull, 0x5365F4583694BB37ull},
    {63, 0x0000000000000003ull, 0xF792AE62B52A1C98ull, 0ull, 0ull, 0xFFFFFFFFFFFFFFFFull, 0x5261A63AA9D18639ull},
    {63, 0x0000000000000000ull, 0x0000000000000015ull, 0ull, 0ull, 0x0000000000000001ull, 0x0000000000000000ull},
    {63, 0x0000000000000000ull, 0x000000000006666Aull, 0ull, 0ull, 0x0000000000000001ull, 0x0000000000000000ull},
    {63, 0x0000000000000000ull, 0x00822B39EBD93AD6ull, 0ull, 0ull, 0x0000000000000000ull, 0xFFFFDEE809A68281ull},
    {63, 0x0000000000000003ull, 0x2782BFA0C523625Eull, 0ull, 0ull, 0xFFFFFFFFFFFFFFFFull, 0x000552D563046DFFull},
};
/// @brief cos: largest error accepted, per I (1.2 times the largest of 2000 random samples + 0.6).
static constexpr MathBound kCosBounds[] = {
    {4, 2.8, 0.0, 0.0},
    {16, 1.4, 0.0, 0.0},
    {32, 1.8, 0.0, 0.0},
    {62, 1.4, 0.0, 0.0},
    {63, 2.2, 0.0, 0.0},
};
/// @brief tan: arguments and the correctly rounded result, from mpmath 1.4.1.
static constexpr MathVector kTanVectors[] = {
    {4, 0xF4735FC7BDE534D7ull, 0x78D995C1298287DDull, 0ull, 0ull, 0xC1698F49650CBEA8ull, 0x5F6791237FA9EDFAull},
    {4, 0xF461EAC5FFF67DDFull, 0xA12E963A6BC8ADD4ull, 0ull, 0ull, 0xBCDDFEB71833D35Full, 0xB996DA36BEEDFFF1ull},
    {4, 0x00000000000002B0ull, 0x20A75105D9AB0586ull, 0ull, 0ull, 0x00000000000002B0ull, 0x20A75105D9AB0586ull},
    {16, 0x0000CB0DB139802Cull, 0xC0E031F339EE03B0ull, 0ull, 0ull, 0xFFDFDDEDDD9449F0ull, 0x3F6A64DD0CB79111ull},
    {16, 0xFFFF3BAB2115ED52ull, 0x64A5F8ED71BECDE8ull, 0ull, 0ull, 0xFFF27A064A86A1F8ull, 0x516BB2C681248855ull},
    {16, 0xFFFFFFFFFFFF2DE2ull, 0xB1EAE2DD4E59E43Cull, 0ull, 0ull, 0xFFFFFFFFFFFF2DE2ull, 0xB1EAE2DD4E592783ull},
    {16, 0xFFFFFFFFFFFFFFFEull, 0x69A898EFC9B20C67ull, 0ull, 0ull, 0xFFFFFFFFFFFFFFFEull, 0x69A898EFC9B20C67ull},
    {16, 0x00000001E5AF74A4ull, 0xC63922AE3EE90E16ull, 0ull, 0ull, 0x00000001E5AF74ADull, 0xE11FDD48611B1417ull},
    {16, 0xFFFFFFFFFFFF9F5Bull, 0x4E270B5A829FE074ull, 0ull, 0ull, 0xFFFFFFFFFFFF9F5Bull, 0x4E270B5A829FCE17ull},
    {32, 0xFFFFFFFF33A9A38Full, 0x539410AADE499545ull, 0ull, 0ull, 0x0000001389047B18ull, 0xE0E573F20303909Cull},
    {32, 0x00000000CCCB9316ull, 0xFFBD2B0317A09522ull, 0ull, 0ull, 0xFFFFFFEEDC9D4251ull, 0xC37EFC8CBB92F612ull},
    {32, 0x0000000000000000ull, 0x002D1718BF69F032ull, 0ull, 0ull, 0x0000000000000000ull, 0x002D1718BF69F032ull},
    {62, 0x0000000000000003ull, 0x51A246B32490E113ull, 0ull, 0ull, 0xFFFFFFFFFFFFFFE9ull, 0x7F47CFF21F95FAA8ull},
    {62, 0x0000000000000002ull, 0xD7FD84F26EFF6682ull, 0ull, 0ull, 0x000000000000000Dull, 0x5429CA79DC4528ACull},
    {62, 0x0000000000000000ull, 0x00000000045B3D00ull, 0ull, 0ull, 0x0000000000000000ull, 0x00000000045B3D00ull},
    {63, 0x0000000000000001ull, 0x8E9EBC8C4E698863ull, 0ull, 0ull, 0x0000000000000049ull, 0x0F1F35F60A0175ADull},
    {63, 0xFFFFFFFFFFFFFFFEull, 0x73AA9A575F9551A4ull, 0ull, 0ull, 0xFFFFFFFFFFFFFFD3ull, 0xCBA5818A4F720E77ull},
    {63, 0xFFFFFFFFFFFFFFFFull, 0xFFDBE2B0FAB00E35ull, 0ull, 0ull, 0xFFFFFFFFFFFFFFFFull, 0xFFDBE2B0BD5B34C0ull},
    {63, 0x0000000000000000ull, 0x0000000000000142ull, 0ull, 0ull, 0x0000000000000000ull, 0x0000000000000142ull},
    {63, 0x0000000000000000ull, 0x000000000000B9F0ull, 0ull, 0ull, 0x0000000000000000ull, 0x000000000000B9F0ull},
    {63, 0x0000000000000000ull, 0x000000AB3769B4F6ull, 0ull, 0ull, 0x0000000000000000ull, 0x000000AB3769B4F6ull},
};
/// @brief tan: largest error accepted, per I (1.2 times the largest of 2000 random samples + 0.6).
static constexpr MathBound kTanBounds[] = {
    {4, 0.6, 13.6, 0.0},
    {16, 0.6, 19.9, 0.0},
    {32, 0.6, 20.7, 0.0},
    {62, 0.6, 6.2, 0.0},
    {63, 0.6, 53.7, 0.0},
};
/// @brief asin: arguments and the correctly rounded result, from mpmath 1.4.1.
static constexpr MathVector kAsinVectors[] = {
    {4, 0x041291D208B1BBEBull, 0x9A3E07C79056658Aull, 0ull, 0ull, 0x0445D69657B2D559ull, 0xF32D45FE18023817ull},
    {4, 0x05FEEAB118EE622Eull, 0xBAEE150294E6BC6Cull, 0ull, 0ull, 0x06C731B5E3FE3422ull, 0x57F03721A609A58Dull},
    {4, 0x0000000000000000ull, 0x0006CFF93ABDC9BAull, 0ull, 0ull, 0x0000000000000000ull, 0x0006CFF93ABDC9BAull},
    {16, 0xFFFF9AD18BB7E8E7ull, 0xFD349A0CBAEDCAA4ull, 0ull, 0ull, 0xFFFF8B50EF48C0B1ull, 0xE74F82ED41DED51Eull},
    {16, 0xFFFFB33DF708C977ull, 0x5E06C74ACE0919BAull, 0ull, 0ull, 0xFFFFADAF35B6A4F8ull, 0x518516FF85F0E392ull},
    {16, 0x000000000000001Bull, 0x579ADE77E3132ADCull, 0ull, 0ull, 0x000000000000001Bull, 0x579ADE77E3132ADCull},
    {16, 0x0000000000000000ull, 0x00D6597696865269ull, 0ull, 0ull, 0x0000000000000000ull, 0x00D6597696865269ull},
    {16, 0x0000000009C38F17ull, 0x64AD2D3511820B4Aull, 0ull, 0ull, 0x0000000009C38F17ull, 0x64AF99C2ADDD11AAull},
    {16, 0x0000000000000000ull, 0x0000000000000E89ull, 0ull, 0ull, 0x0000000000000000ull, 0x0000000000000E89ull},
    {32, 0xFFFFFFFFA374BB32ull, 0xED53E2EEE5CCB90Dull, 0ull, 0ull, 0xFFFFFFFF988F0685ull, 0xFE5065745C742C3Full},
    {32, 0x000000007237EF26ull, 0x6D414876FA382586ull, 0ull, 0ull, 0x000000008D1E3DE0ull, 0x42473B70C4CFBD35ull},
    {32, 0xFFFFFFFFFEBAF5BCull, 0x7C0AC2F9079F7E03ull, 0ull, 0ull, 0xFFFFFFFFFEBAF45Full, 0x22D5B701CEE34C5Cull},
    {62, 0xFFFFFFFFFFFFFFFEull, 0x7B3441A630323AC3ull, 0ull, 0ull, 0xFFFFFFFFFFFFFFFEull, 0x467B6C0548988957ull},
    {62, 0xFFFFFFFFFFFFFFFEull, 0x98626DD548891DF6ull, 0ull, 0ull, 0xFFFFFFFFFFFFFFFEull, 0x714AA8E7868E2631ull},
    {62, 0x0000000000000000ull, 0x0000000000000004ull, 0ull, 0ull, 0x0000000000000000ull, 0x0000000000000004ull},
    {63, 0xFFFFFFFFFFFFFFFFull, 0x39DDD996ACB1DF0Cull, 0ull, 0ull, 0xFFFFFFFFFFFFFFFFull, 0x1D6C20B36C0E0564ull},
    {63, 0xFFFFFFFFFFFFFFFFull, 0x3D8C3F7F79745E52ull, 0ull, 0ull, 0xFFFFFFFFFFFFFFFFull, 0x23285FE220D13B06ull},
    {63, 0x0000000000000000ull, 0x0000000000000003ull, 0ull, 0ull, 0x0000000000000000ull, 0x0000000000000003ull},
    {63, 0xFFFFFFFFFFFFFFFFull, 0xFFD78D6F448A78D2ull, 0ull, 0ull, 0xFFFFFFFFFFFFFFFFull, 0xFFD78D6F1975C426ull},
    {63, 0x0000000000000000ull, 0x00000000000D2664ull, 0ull, 0ull, 0x0000000000000000ull, 0x00000000000D2664ull},
    {63, 0x0000000000000000ull, 0x00000000000359FEull, 0ull, 0ull, 0x0000000000000000ull, 0x00000000000359FEull},
};
/// @brief asin: largest error accepted, per I (1.2 times the largest of 2000 random samples + 0.6).
static constexpr MathBound kAsinBounds[] = {
    {4, 3.9, 0.0, 0.0},
    {16, 2.7, 0.0, 0.0},
    {32, 3.9, 0.0, 0.0},
    {62, 3.4, 0.0, 0.0},
    {63, 3.4, 0.0, 0.0},
};
/// @brief acos: arguments and the correctly rounded result, from mpmath 1.4.1.
static constexpr MathVector kAcosVectors[] = {
    {4, 0xF9EC6CC3ED0EB18Eull, 0xC1D48618121DB1BCull, 0ull, 0ull, 0x1377AA17CA92B337ull, 0x3DD64BE620FE7F2Bull},
    {4, 0xF8B43214D62F2011ull, 0x347F8ED9745E860Eull, 0ull, 0ull, 0x15C0701A81DB288Aull, 0xDC633604FFB707CFull},
    {4, 0x0000000000000000ull, 0x0000000000000001ull, 0ull, 0ull, 0x0C90FDAA22168C23ull, 0x4C4C6628B80DC1CCull},
    {16, 0xFFFFB33DF708C977ull, 0x5E06C74ACE0919BAull, 0ull, 0ull, 0x00011B60A4EB7C70ull, 0x70AFADC6DC9A9D4Aull},
    {16, 0xFFFFBEA9DAC8BA60ull, 0x889B585029559457ull, 0ull, 0ull, 0x00010DA19BCD02C0ull, 0x4F8A6B29155B6050ull},
    {16, 0xFFFFFFFFFFFFFFFFull, 0xFFFFFFFFFFFFD328ull, 0ull, 0ull, 0x0000C90FDAA22168ull, 0xC234C4C6628BADB4ull},
    {16, 0x0000000000000000ull, 0x000E59B348486A8Eull, 0ull, 0ull, 0x0000C90FDAA22168ull, 0xC2266B131A43164Eull},
    {16, 0xFFFFFFDDC8C620D0ull, 0xFB2375706B1885A5ull, 0ull, 0ull, 0x0000C93211DC68E8ull, 0x63975500190B9408ull},
    {16, 0xFFFFFFFF52EB4BC8ull, 0x5A5C12008245E171ull, 0ull, 0ull, 0x0000C91087B6D5A0ull, 0x9C9741614DFE4995ull},
    {32, 0xFFFFFFFFA374BB32ull, 0xED53E2EEE5CCB90Dull, 0ull, 0ull, 0x000000013080D41Cull, 0x23185CC06852364Dull},
    {32, 0x000000007237EF26ull, 0x6D414876FA382586ull, 0ull, 0ull, 0x000000003BF19CC1ull, 0xDF2186C3FFF6A556ull},
    {32, 0x0000000000000002ull, 0x5520094B2C8A77CAull, 0ull, 0ull, 0x00000000C90FDA9Full, 0xCC48B8E9983BEAB9ull},
    {62, 0xFFFFFFFFFFFFFFFEull, 0x7B3441A630323AC3ull, 0ull, 0ull, 0x0000000000000004ull, 0xDDC3FE833D0A7F7Cull},
    {62, 0xFFFFFFFFFFFFFFFEull, 0x98626DD548891DF6ull, 0ull, 0ull, 0x0000000000000004ull, 0xB2F4C1A0FF14E2A2ull},
    {62, 0xFFFFFFFFFFFFFFFFull, 0xFFFFFE7FAD1910C2ull, 0ull, 0ull, 0x0000000000000003ull, 0x243F6C08D889F811ull},
    {63, 0xFFFFFFFFFFFFFFFFull, 0x7D17FE64C4123E83ull, 0ull, 0ull, 0x0000000000000002ull, 0x1B88AECAD9873757ull},
    {63, 0x0000000000000000ull, 0xA0DDAD5952B049CEull, 0ull, 0ull, 0x0000000000000000ull, 0xE42DD197415D4274ull},
    {63, 0x0000000000000000ull, 0x0000000000000192ull, 0ull, 0ull, 0x0000000000000001ull, 0x921FB54442D182D8ull},
    {63, 0x0000000000000000ull, 0x0000000006F88123ull, 0ull, 0ull, 0x0000000000000001ull, 0x921FB5443BD90347ull},
    {63, 0xFFFFFFFFFFFFFFFFull, 0xFFFFFFFF6CBF977Eull, 0ull, 0ull, 0x0000000000000001ull, 0x921FB544D611ECECull},
    {63, 0x0000000000000000ull, 0x000000222E008026ull, 0ull, 0ull, 0x0000000000000001ull, 0x921FB52214D10444ull},
};
/// @brief acos: largest error accepted, per I (1.2 times the largest of 2000 random samples + 0.6).
static constexpr MathBound kAcosBounds[] = {
    {4, 4.8, 0.0, 0.0},
    {16, 2.8, 0.0, 0.0},
    {32, 3.3, 0.0, 0.0},
    {62, 3.5, 0.0, 0.0},
    {63, 3.4, 0.0, 0.0},
};
/// @brief atan: arguments and the correctly rounded result, from mpmath 1.4.1.
static constexpr MathVector kAtanVectors[] = {
    {4, 0xE00A5A639B1FF9EAull, 0x38763F7AD29CE9FDull, 0ull, 0ull, 0xF56555CCB5EAC205ull, 0x76C9BE41D4E7E266ull},
    {4, 0xD3D1BA6BA4C70E50ull, 0x10CAD983B72A9A35ull, 0ull, 0ull, 0xF4DDDFCF5B4BBC1Full, 0xD0DCAC0FED63D588ull},
    {4, 0xFFFFFFFFFFFFFFFFull, 0xFFFFFD3297874291ull, 0ull, 0ull, 0xFFFFFFFFFFFFFFFFull, 0xFFFFFD3297874291ull},
    {16, 0xF504688FB31128A8ull, 0x6C1F4B5F5150E35Cull, 0ull, 0ull, 0xFFFF36F5F928170Dull, 0xD59E9C7CF7FA5567ull},
    {16, 0xFFFD05E88743437Eull, 0x335D14360BC1332Cull, 0ull, 0ull, 0xFFFF4C3CECDCD440ull, 0xB89F8F82DFF69E02ull},
    {16, 0x0000000000000000ull, 0x0000000000005BA1ull, 0ull, 0ull, 0x0000000000000000ull, 0x0000000000005BA1ull},
    {16, 0x000000000005AA83ull, 0x3F89386D14254F6Cull, 0ull, 0ull, 0x000000000005AA83ull, 0x3F89386D1332C52Cull},
    {16, 0x0000000000000000ull, 0x00037617F902B4ACull, 0ull, 0ull, 0x0000000000000000ull, 0x00037617F902B4ACull},
    {16, 0x0000000000000000ull, 0x0000000000002F64ull, 0ull, 0ull, 0x0000000000000000ull, 0x0000000000002F64ull},
    {32, 0xFFFBD45126FAD476ull, 0xEA642114CF5C070Full, 0ull, 0ull, 0xFFFFFFFF36F034B6ull, 0x4955115F6B0CF5A3ull},
    {32, 0xFFFFFFEE9CA9C9F5ull, 0x2DC868DB40333987ull, 0ull, 0ull, 0xFFFFFFFF3A9E24C0ull, 0x31FB9CB7A76DD262ull},
    {32, 0xFFFFFFFFFFFFFFFFull, 0xFFFFFFFD7829E6FDull, 0ull, 0ull, 0xFFFFFFFFFFFFFFFFull, 0xFFFFFFFD7829E6FDull},
    {62, 0x0000000000000047ull, 0x9B654A48C890560Full, 0ull, 0ull, 0x0000000000000003ull, 0x15F37F4AB2BA8383ull},
    {62, 0x0000000000000008ull, 0xC452F3DCCFD3047Cull, 0ull, 0ull, 0x0000000000000002ull, 0xB168EF35EA42773Aull},
    {62, 0x000000CC3ABBF204ull, 0x569D19E8A4E41B46ull, 0ull, 0ull, 0x0000000000000003ull, 0x243F6A88809F755Dull},
    {63, 0xFFFFFFFFFFFFFF5Dull, 0x5128ADFC3523D4AEull, 0ull, 0ull, 0xFFFFFFFFFFFFFFFEull, 0x6F7321B82DF76EDEull},
    {63, 0xFFFFFFFFFFFFFC80ull, 0x708D3A83CE4350E3ull, 0ull, 0ull, 0xFFFFFFFFFFFFFFFEull, 0x6E29787D4957A5A5ull},
    {63, 0x0000000000000000ull, 0x1A565101C9C485B1ull, 0ull, 0ull, 0x0000000000000000ull, 0x1A3EADD3071973F7ull},
    {63, 0xFFFFFFFFC858EA9Full, 0xB95F5C6090B5898Full, 0ull, 0ull, 0xFFFFFFFFFFFFFFFEull, 0x6DE04AC056C5655Bull},
    {63, 0x0000000000000000ull, 0x0000000422B40170ull, 0ull, 0ull, 0x0000000000000000ull, 0x0000000422B40170ull},
    {63, 0x000014E6A38AC7F2ull, 0x271F1B6AC0E6D48Dull, 0ull, 0ull, 0x0000000000000001ull, 0x921FB54442C544DCull},
};
/// @brief atan: largest error accepted, per I (1.2 times the largest of 2000 random samples + 0.6).
static constexpr MathBound kAtanBounds[] = {
    {4, 3.3, 0.0, 0.0},
    {16, 2.2, 0.0, 0.0},
    {32, 3.0, 0.0, 0.0},
    {62, 2.3, 0.0, 0.0},
    {63, 2.7, 0.0, 0.0},
};
/// @brief atan2: arguments and the correctly rounded result, from mpmath 1.4.1.
static constexpr MathVector kAtan2Vectors[] = {
    {4, 0x0000000000000000ull, 0x0000000F20EA609Eull, 0xFFFFFFFFFFFFFFD8ull, 0x79965116C2A723CEull, 0x1921FB54411D3439ull, 0x3767F9C6F5A01F4Full},
    {4, 0x0000000000000000ull, 0x000000659E08CE1Bull, 0xFFFFFFFFFFFF7717ull, 0xFA7748D72E386991ull, 0x1921FB544427282Bull, 0x3548CB862615433Full},
    {4, 0x0000000000000000ull, 0x00000E8C8D868318ull, 0x0000000000000000ull, 0x00000A3F8287797Bull, 0x07A83676B3E8936Full, 0xF35EF4A972095EF5ull},
    {16, 0xFFFFFFFFFFFFF945ull, 0xDFF1933DDFEFE00Aull, 0x0000000000000001ull, 0xCB6E1C8854DF9CDCull, 0xFFFF37124B3D7BF3ull, 0x961FAD1F2F59DCA1ull},
    {16, 0x0000000000000000ull, 0x0000004AE82BE942ull, 0xFFFFFFFFFFFFFFFFull, 0xFFFFFFFFFE92CD29ull, 0x0000C9124AAE1913ull, 0x400449F91AD09DDDull},
    {16, 0x0000000000000000ull, 0x00000000000FDE15ull, 0x0000000000000000ull, 0x03862D1DD30E8E86ull, 0x0000000000000240ull, 0x5323430C2EC55EDBull},
    {16, 0xFFFFFFFFFFFE4EF5ull, 0xB8D4539FAFCC8967ull, 0xFFFFFFFFDFF19654ull, 0x00ABAE761B4D5E3Eull, 0xFFFE6DE70BDA0D05ull, 0x556B8801CD662C3Aull},
    {16, 0xFFFFFFFFFFFFFFFFull, 0xFFFDED3E1F7EA5BBull, 0x0000000000000000ull, 0x0000019B0F23BB2Dull, 0xFFFF37534756BB80ull, 0xE91EDD71F25E1A71ull},
    {16, 0xFFFFFFFFFFFFFFFFull, 0xFFFFFFFFFFFFF508ull, 0x00000005D7365D66ull, 0x465695896AD0F818ull, 0xFFFFFFFFFFFFFFFFull, 0xFFFFFFFFFF0F9DE4ull},
    {32, 0x0000000000000000ull, 0x000003779EE8BED8ull, 0xFFFFFFFFFFFFFFFFull, 0xFFFFFFF34395A7A4ull, 0x00000000CAE5FE0Cull, 0x0A5BCDA504D93250ull},
    {32, 0xFEFF073BAFECC9E7ull, 0x7633032633D6B645ull, 0x00000000297C038Full, 0x94FDABEE88F9A820ull, 0xFFFFFFFF36F02572ull, 0x88848F1E97E5F27Full},
    {32, 0xFFFFFFFFEF808F4Bull, 0xE7729DB8D5C60698ull, 0x0000000000000000ull, 0x03F611E444083E3Aull, 0xFFFFFFFF36F0255Dull, 0xFD53020B51811156ull},
    {62, 0x0000000000000000ull, 0x00003CAB56FCF09Eull, 0x0000000000000000ull, 0x000008A4C3FD76A0ull, 0x0000000000000002ull, 0xDBCA49A927BDBB56ull},
    {62, 0x0000000000000000ull, 0x0000000010C86AFCull, 0x0000000000000000ull, 0x000000000002E9C9ull, 0x0000000000000003ull, 0x23E68AAA14BF5DE3ull},
    {62, 0x0000000000000000ull, 0x00010120B84DE2CBull, 0x0000000000006DD5ull, 0x47DEA4DE3CC5A5AEull, 0x0000000000000000ull, 0x00000004AEA1A5D5ull},
    {63, 0xFFFFFFFFFFFFFFFFull, 0xFFFFFFFFFF18A141ull, 0x0000000000000000ull, 0x0000000000043789ull, 0xFFFFFFFFFFFFFFFEull, 0x728A9D2AFAAD1ED3ull},
    {63, 0x0000000000000000ull, 0x0000000000032578ull, 0xFFFFFFFFFFFFFFFFull, 0xFFFFFFFFFFFFFF03ull, 0x0000000000000001ull, 0x92701E4C5A69A0A0ull},
    {63, 0x0000000000000000ull, 0x000000005211E84Aull, 0x0000000000000000ull, 0x000772C2E5B03286ull, 0x0000000000000000ull, 0x00000B04C546DEBCull},
    {63, 0xFFFF75249B72D404ull, 0x4E507F58C85192D6ull, 0x0052D4A772C4DAAEull, 0x5F5BDA8FD145BE50ull, 0xFFFFFFFFFFFFFFFFull, 0xFE52D94BCC755CEAull},
    {63, 0x0000000000000000ull, 0x00024D6030E85E46ull, 0x0000000000000000ull, 0x000FB90A909EE514ull, 0x0000000000000000ull, 0x25388A97716029A1ull},
    {63, 0xFFFFFFFFFFFFEA83ull, 0xFD15BBE25615FDE0ull, 0xFFFFFFFFFFFFFFFFull, 0xFEE346A508DF419Cull, 0xFFFFFFFFFFFFFFFEull, 0x6DE03D7B13ED43B5ull},
};
/// @brief atan2: largest error accepted, per I (1.2 times the largest of 2000 random samples + 0.6).
static constexpr MathBound kAtan2Bounds[] = {
    {4, 4.6, 0.0, 0.0},
    {16, 2.5, 0.0, 0.0},
    {32, 2.7, 0.0, 0.0},
    {62, 2.2, 0.0, 0.0},
    {63, 2.8, 0.0, 0.0},
};
/// @brief sinh: arguments and the correctly rounded result, from mpmath 1.4.1.
static constexpr MathVector kSinhVectors[] = {
    {4, 0xF815DB3A2A217219ull, 0xF28B5F71A432C02Full, 0ull, 0ull, 0xF6BAC7146BB9B9C0ull, 0x45D3650837E8E87Eull},
    {4, 0x0FB6B0C22D1F17D0ull, 0x558CDD72996B99C7ull, 0ull, 0ull, 0x1BF4B4A46FB5E573ull, 0xACCDB617B1BD34FCull},
    {4, 0x0000000000000067ull, 0x2BCAF62BDE908362ull, 0ull, 0ull, 0x0000000000000067ull, 0x2BCAF62BDE908362ull},
    {16, 0xFFFFB0422A172455ull, 0x8A48EA3E043B6099ull, 0ull, 0ull, 0xFFFFAAFFD9B0069Full, 0x96CA5F2EBE2AB314ull},
    {16, 0x00003471A2E815FCull, 0x6C8FEF21C88377EEull, 0ull, 0ull, 0x000035EC6D094FAFull, 0x03A0C08495579959ull},
    {16, 0xFFFFFFFFFFFFFFFFull, 0xFFFFFFFFF0AE1650ull, 0ull, 0ull, 0xFFFFFFFFFFFFFFFFull, 0xFFFFFFFFF0AE1650ull},
    {16, 0x0000000000000000ull, 0x0000000000000001ull, 0ull, 0ull, 0x0000000000000000ull, 0x0000000000000001ull},
    {16, 0x000000000001C0BCull, 0xDF2E60023D858BCAull, 0ull, 0ull, 0x000000000001C0BCull, 0xDF2E60023D8922FBull},
    {16, 0xFFFFFFFFFFFFFEDEull, 0x1EBF48271C77CE36ull, 0ull, 0ull, 0xFFFFFFFFFFFFFEDEull, 0x1EBF48271C77CE36ull},
    {32, 0x000000009055F234ull, 0x1531DB250B13117Aull, 0ull, 0ull, 0x00000000B0EDA3D1ull, 0xD2F898B1D2BC8592ull},
    {32, 0x000000002D538CECull, 0x2DE8BB5F8575B58Eull, 0ull, 0ull, 0x000000002E4794BAull, 0x3CA6B0E7A50E46B6ull},
    {32, 0x0000000000000003ull, 0xC3E9B981D80A5FDFull, 0ull, 0ull, 0x0000000000000003ull, 0xC3E9B981D80A6003ull},
    {62, 0x0000000000000000ull, 0x768093DDD0DAD9FFull, 0ull, 0ull, 0x0000000000000000ull, 0x77902763F56855E9ull},
    {62, 0xFFFFFFFFFFFFFFFEull, 0x57C9BA025133B33Eull, 0ull, 0ull, 0xFFFFFFFFFFFFFFFEull, 0x258F2D8F34A47A8Cull},
    {62, 0xFFFFFFFFFFFFFFFFull, 0xFFFFFFFFE1F12D24ull, 0ull, 0ull, 0xFFFFFFFFFFFFFFFFull, 0xFFFFFFFFE1F12D24ull},
    {63, 0xFFFFFFFFFFFFFFFFull, 0x5D684639190C6A1Cull, 0ull, 0ull, 0xFFFFFFFFFFFFFFFFull, 0x5240DFDA70BBA448ull},
    {63, 0xFFFFFFFFFFFFFFFFull, 0x418BBE2458266E5Cull, 0ull, 0ull, 0xFFFFFFFFFFFFFFFFull, 0x2F7C0860EAB0457Aull},
    {63, 0xFFFFFFFFFFFFFFFFull, 0xFFFFFFFFFFFFFF29ull, 0ull, 0ull, 0xFFFFFFFFFFFFFFFFull, 0xFFFFFFFFFFFFFF29ull},
    {63, 0xFFFFFFFFFFFFFFFFull, 0xFFFFFFFFFFFFFFE5ull, 0ull, 0ull, 0xFFFFFFFFFFFFFFFFull, 0xFFFFFFFFFFFFFFE5ull},
    {63, 0xFFFFFFFFFFFFFFFFull, 0xFFFFFFFBCC7EB9D7ull, 0ull, 0ull, 0xFFFFFFFFFFFFFFFFull, 0xFFFFFFFBCC7EB9D7ull},
    {63, 0xFFFFFFFFFFFFFFFFull, 0xFFFFFFFFEC11E62Bull, 0ull, 0ull, 0xFFFFFFFFFFFFFFFFull, 0xFFFFFFFFEC11E62Bull},
};
/// @brief sinh: largest error accepted, per I (1.2 times the largest of 2000 random samples + 0.6).
static constexpr MathBound kSinhBounds[] = {
    {4, 0.6, 11.6, 0.0},
    {16, 0.6, 4.0, 0.0},
    {32, 0.6, 2.9, 0.0},
    {62, 0.6, 3.4, 0.0},
    {63, 0.6, 3.3, 0.0},
};
/// @brief cosh: arguments and the correctly rounded result, from mpmath 1.4.1.
static constexpr MathVector kCoshVectors[] = {
    {4, 0xF815DB3A2A217219ull, 0xF28B5F71A432C02Full, 0ull, 0ull, 0x0C3EB916A3160129ull, 0x137C36EF95D07621ull},
    {4, 0x07F691962A5EB8AAull, 0x49A8F2C955BF2B01ull, 0ull, 0ull, 0x0C4D2DCE8CCF89B2ull, 0xB7FFC050809AE621ull},
    {4, 0x000000AEAFB2B24Cull, 0xD6E565368A39C7D0ull, 0ull, 0ull, 0x0800000000077334ull, 0xF802D081CEC65671ull},
    {16, 0x0000010AF0B6D88Cull, 0x66F4272257A095C8ull, 0ull, 0ull, 0x0000800116598333ull, 0x4F502718B417DCF6ull},
    {16, 0x000045CA0F81B32Bull, 0x019DA21F864F14C4ull, 0ull, 0ull, 0x0000938067DC67A7ull, 0x66EC954F7FF3C563ull},
    {16, 0xFFFFFFFFFFFFFFFFull, 0xFFFFFFFFFFFF8EE0ull, 0ull, 0ull, 0x0000800000000000ull, 0x0000000000000000ull},
    {16, 0xFFFFFFFFFFFFFFFFull, 0xFFFFFFFFFFFFFF87ull, 0ull, 0ull, 0x0000800000000000ull, 0x0000000000000000ull},
    {16, 0x000001AA0FA9C859ull, 0xE4F4D2B73C6B28ACull, 0ull, 0ull, 0x00008002C51AB0B9ull, 0xF7F470C99C6DA9F4ull},
    {16, 0x00000000000009DAull, 0x457F86FC46CC9E67ull, 0ull, 0ull, 0x0000800000000000ull, 0x0000006112FD678Eull},
    {32, 0xFFFFFFFFE3BC51A9ull, 0x089898B787666884ull, 0ull, 0ull, 0x000000008322205Full, 0x9F6110E94C721C5Cull},
    {32, 0xFFFFFFFFFEA1E7D9ull, 0xC15F96776072F37Cull, 0ull, 0ull, 0x000000008001DEC7ull, 0x3552C86C0D27E21Cull},
    {32, 0xFFFFFFFFFFFFFFFFull, 0xFFFF8F701D3C0422ull, 0ull, 0ull, 0x0000000080000000ull, 0x0000000000000000ull},
    {62, 0xFFFFFFFFFFFFFFFEull, 0x37872F65DD8AF3FDull, 0ull, 0ull, 0x0000000000000002ull, 0xD952D4C4874E0BCEull},
    {62, 0xFFFFFFFFFFFFFFFFull, 0xA6DE46E2BC6EADF5ull, 0ull, 0ull, 0x0000000000000002ull, 0x07C72327AAC3AD4Dull},
    {62, 0xFFFFFFFFFFFFFFFFull, 0xFFFFFFFFFFFFFF7Bull, 0ull, 0ull, 0x0000000000000002ull, 0x0000000000000000ull},
    {63, 0xFFFFFFFFFFFFFFFFull, 0x52F51178BAEC548Cull, 0ull, 0ull, 0x0000000000000001ull, 0x3CBEB41E612CBB24ull},
    {63, 0xFFFFFFFFFFFFFFFFull, 0x8EC9741958D79F55ull, 0ull, 0ull, 0x0000000000000001ull, 0x1971BADDDDE58A49ull},
    {63, 0xFFFFFFFFFFFFFFFFull, 0xFFFFFFFFFFFFFFA8ull, 0ull, 0ull, 0x0000000000000001ull, 0x0000000000000000ull},
    {63, 0xFFFFFFFFFFFFFFFFull, 0xFFFFFFFFFFFFB9B2ull, 0ull, 0ull, 0x0000000000000001ull, 0x0000000000000000ull},
    {63, 0x0000000000000000ull, 0x000000000000049Cull, 0ull, 0ull, 0x0000000000000001ull, 0x0000000000000000ull},
    {63, 0xFFFFFFFFFFFFFFFFull, 0xFDA6786EFBE1A4E7ull, 0ull, 0ull, 0x0000000000000001ull, 0x0002C2B82C80ED22ull},
};
/// @brief cosh: largest error accepted, per I (1.2 times the largest of 2000 random samples + 0.6).
static constexpr MathBound kCoshBounds[] = {
    {4, 0.6, 12.5, 0.0},
    {16, 0.6, 3.1, 0.0},
    {32, 0.6, 3.4, 0.0},
    {62, 0.6, 2.8, 0.0},
    {63, 0.6, 3.6, 0.0},
};
/// @brief tanh: arguments and the correctly rounded result, from mpmath 1.4.1.
static constexpr MathVector kTanhVectors[] = {
    {2, 0xF48D61B45D06BF8Bull, 0x6ECAE7F9AE785C44ull, 0ull, 0ull, 0xF50450EB00260726ull, 0xF89937011ABD0C64ull},
    {2, 0xF424E0A151987C07ull, 0xE1A112558B3CB6EDull, 0ull, 0ull, 0xF4A88791ABA32B27ull, 0x7284F69153594C43ull},
    {2, 0xFFFFFFFFFFFFFFFFull, 0xFFFFFFFFFFF06D38ull, 0ull, 0ull, 0xFFFFFFFFFFFFFFFFull, 0xFFFFFFFFFFF06D38ull},
    {4, 0x02F56AE7A198D54Cull, 0x7A8CEE957D9CACAAull, 0ull, 0ull, 0x02D4ACEE6A22BF59ull, 0xDB8E33BB919C51ADull},
    {4, 0xFC2A01B8BB1C1B44ull, 0x36E821A9AB006E05ull, 0ull, 0ull, 0xFC6EEE7B033445E1ull, 0x466FBFC4E74DCE4Bull},
    {4, 0xFFF86E920FF99B1Bull, 0x337AF8D1B359C19Aull, 0ull, 0ull, 0xFFF86E9451EE2150ull, 0xA5BD7C0E2892EE0Aull},
    {16, 0xFFFFF52B2DAEE360ull, 0xFF4DDB3712F8DFE0ull, 0ull, 0ull, 0xFFFFF531C71F458Full, 0x56AFA9EF34137A74ull},
    {16, 0x00002B4275A189A8ull, 0xB4B2DD98565770FEull, 0ull, 0ull, 0x000029AF3A5AC2AAull, 0xCDBD6F1B39B584B3ull},
    {16, 0x0000000000000E2Dull, 0x0DBA0D3EC28097A4ull, 0ull, 0ull, 0x0000000000000E2Dull, 0x0DBA0D3EC2809795ull},
    {16, 0x0000000000000017ull, 0xA8B8BA8D3F3BE3F3ull, 0ull, 0ull, 0x0000000000000017ull, 0xA8B8BA8D3F3BE3F3ull},
    {16, 0x0000003BBA45369Cull, 0x6CE8E0AA91ECDACAull, 0ull, 0ull, 0x0000003BBA40E0DDull, 0x50E0A7EAC086DBDEull},
    {16, 0xFFFFFFFFFFFFFFFFull, 0xFFFFFF5249450721ull, 0ull, 0ull, 0xFFFFFFFFFFFFFFFFull, 0xFFFFFF5249450721ull},
    {32, 0x000000002C7E4957ull, 0x6907367D7613C2C6ull, 0ull, 0ull, 0x000000002AC8AB00ull, 0xDCB571156C3CAD5Aull},
    {32, 0x0000000031EE15AAull, 0xE0F343602664897Bull, 0ull, 0ull, 0x000000002F8AF15Dull, 0xB64FAB75DD09AD91ull},
    {32, 0x0000000000000000ull, 0x0000000000007971ull, 0ull, 0ull, 0x0000000000000000ull, 0x0000000000007971ull},
    {62, 0x0000000000000000ull, 0x960CF381D09DAEE8ull, 0ull, 0ull, 0x0000000000000000ull, 0x91E5B825149FFC63ull},
    {62, 0xFFFFFFFFFFFFFFFFull, 0x6362F1D9EA53437Bull, 0ull, 0ull, 0xFFFFFFFFFFFFFFFFull, 0x68184D749B1DDF62ull},
    {62, 0x0000000000000000ull, 0x00000000013BDBABull, 0ull, 0ull, 0x0000000000000000ull, 0x00000000013BDBABull},
    {63, 0x0000000000000000ull, 0xFBA75ED510C943F9ull, 0ull, 0ull, 0x0000000000000000ull, 0xC11E7BDD8BBEE453ull},
    {63, 0x0000000000000000ull, 0x72E597FC075ADA8Aull, 0ull, 0ull, 0x0000000000000000ull, 0x6BC1BFB5D556C3F2ull},
    {63, 0xFFFFFFFFFFFFFFFFull, 0xFF985AC7CBF583F0ull, 0ull, 0ull, 0xFFFFFFFFFFFFFFFFull, 0xFF985ACD75B16647ull},
    {63, 0xFFFFFFFFFFFFFFFFull, 0xFFFFFFFE674841CBull, 0ull, 0ull, 0xFFFFFFFFFFFFFFFFull, 0xFFFFFFFE674841CBull},
    {63, 0xFFFFFFFFFFFFFFFFull, 0xFFFFFFFF6AD69D96ull, 0ull, 0ull, 0xFFFFFFFFFFFFFFFFull, 0xFFFFFFFF6AD69D96ull},
    {63, 0x0000000000000000ull, 0x000000000003CF5Bull, 0ull, 0ull, 0x0000000000000000ull, 0x000000000003CF5Bull},
};
/// @brief tanh: largest error accepted, per I (1.2 times the largest of 2000 random samples + 0.6).
static constexpr MathBound kTanhBounds[] = {
    {2, 4.2, 0.0, 0.0},
    {4, 5.6, 0.0, 0.0},
    {16, 3.1, 0.0, 0.0},
    {32, 4.2, 0.0, 0.0},
    {62, 3.7, 0.0, 0.0},
    {63, 4.8, 0.0, 0.0},
};
/// @brief asinh: arguments and the correctly rounded result, from mpmath 1.4.1.
static constexpr MathVector kAsinhVectors[] = {
    {1, 0xCB4E7DDBA5134AD9ull, 0x38072ED81DD13379ull, 0ull, 0ull, 0xCFF3206E2DEDB281ull, 0x96B414B48284DC99ull},
    {1, 0x75000635413A0091ull, 0x3F5BA4501C645349ull, 0ull, 0ull, 0x574C17074973BC75ull, 0x68C5220A53228F2Dull},
    {1, 0x0000000076EFED86ull, 0xF6AD8C8313F4E277ull, 0ull, 0ull, 0x0000000076EFED86ull, 0xF6AD8C82CF7EE9D0ull},
    {2, 0x0FA3B3C34EAE650Dull, 0x381115B8D740CDE1ull, 0ull, 0ull, 0x0F1359208B2CE515ull, 0xD9AA7F23ACCAF370ull},
    {2, 0xFFFFFFFFF9BAB922ull, 0x23A8345235C13C5Full, 0ull, 0ull, 0xFFFFFFFFF9BAB922ull, 0x23A8345235CB8263ull},
    {2, 0xFD0B0FEDE7EB9050ull, 0xE7449CA8B07BAAB1ull, 0ull, 0ull, 0xFD0C229BD46B9D05ull, 0x3A72D6D4EFE93787ull},
    {4, 0xB095DF6BED209F0Bull, 0xD7868B4A4AF53F05ull, 0ull, 0ull, 0xE8129B48D35F2F1Cull, 0x5E8C9185641AE876ull},
    {4, 0x4EE5AD43E6A531B0ull, 0x95490BC9708111B3ull, 0ull, 0ull, 0x17E0132E90E6CC2Full, 0xC819931AF1DCDB04ull},
    {4, 0xFFFFFFFFFFFFFFFFull, 0xFFFFFFFFFFFCC84Aull, 0ull, 0ull, 0xFFFFFFFFFFFFFFFFull, 0xFFFFFFFFFFFCC84Aull},
    {16, 0x5BF00FD35359C041ull, 0x8302B7CE7B7352CDull, 0ull, 0ull, 0x0005B9EE40FC994Bull, 0x07CD98984F64D0B5ull},
    {16, 0x82036BB7284D5054ull, 0x97DFD5E63AD35904ull, 0ull, 0ull, 0xFFFA1DBDB8BA2C78ull, 0x141D421F6C5B3E9Cull},
    {16, 0x1409E362E4ADFB50ull, 0x072C37C387BD60BBull, 0ull, 0ull, 0x0004F6EDD8D2F597ull, 0x5B3948EFCAF6344Dull},
    {16, 0x0000000000000000ull, 0x0000010B829EAF83ull, 0ull, 0ull, 0x0000000000000000ull, 0x0000010B829EAF83ull},
    {16, 0xFFFFFFFFFFFFFFFFull, 0xFFFFFFBCBCEED0CBull, 0ull, 0ull, 0xFFFFFFFFFFFFFFFFull, 0xFFFFFFBCBCEED0CBull},
    {16, 0x0000000000BED4D6ull, 0x3856407F8E651FA8ull, 0ull, 0ull, 0x0000000000BED4D6ull, 0x38564038DCE8C8B6ull},
    {32, 0x0000000002D39D07ull, 0x564E99646C1E3B65ull, 0ull, 0ull, 0x0000000002D38DF9ull, 0xDAF6ADA074497A77ull},
    {32, 0x0000000004014E51ull, 0x9742BAB632291CCAull, 0ull, 0ull, 0x0000000004012381ull, 0xE8AE59E5FC8975A1ull},
    {32, 0xFFFFFFFFFFFFFFFFull, 0xFFFFFB3249810029ull, 0ull, 0ull, 0xFFFFFFFFFFFFFFFFull, 0xFFFFFB3249810029ull},
    {62, 0x869BF39CAFDF2D23ull, 0x09F3D6C4B693B917ull, 0ull, 0ull, 0xFFFFFFFFFFFFFFA8ull, 0xC4FCFE1702C13465ull},
    {62, 0x7242741D879BC097ull, 0x765AB50D4FCCEDFBull, 0ull, 0ull, 0x0000000000000057ull, 0x1C039BA0E6881AD2ull},
    {62, 0xFFFFFFFFFFFFFFFFull, 0xFE81CC832C50D78Bull, 0ull, 0ull, 0xFFFFFFFFFFFFFFFFull, 0xFE81CCA6AB5E21C5ull},
    {63, 0xFFFFFFFFFFFFFFFFull, 0xC1C11826F9801190ull, 0ull, 0ull, 0xFFFFFFFFFFFFFFFFull, 0xC25A1273CB2E0E85ull},
    {63, 0x0000000000000000ull, 0x00072C30F9ACE793ull, 0ull, 0ull, 0x0000000000000000ull, 0x00072C30F96F674Dull},
    {63, 0xFFFFFFFFFFFFFFFFull, 0xFFFFFFFFFFFFF094ull, 0ull, 0ull, 0xFFFFFFFFFFFFFFFFull, 0xFFFFFFFFFFFFF094ull},
    {63, 0x0000000000000000ull, 0x0000000000001EA5ull, 0ull, 0ull, 0x0000000000000000ull, 0x0000000000001EA5ull},
    {63, 0x0000000000000000ull, 0x00000036E4383201ull, 0ull, 0ull, 0x0000000000000000ull, 0x00000036E4383201ull},
    {63, 0xFFFFFFFFF7628674ull, 0xE7D21450A20E6AB8ull, 0ull, 0ull, 0xFFFFFFFFFFFFFFECull, 0x848E87DF43439A4Aull},
};
/// @brief asinh: largest error accepted, per I (1.2 times the largest of 2000 random samples + 0.6).
static constexpr MathBound kAsinhBounds[] = {
    {1, 3.9, 0.0, 0.0},
    {2, 2.4, 0.0, 0.0},
    {4, 2.8, 0.0, 0.0},
    {16, 2.6, 0.0, 0.0},
    {32, 2.2, 0.0, 0.0},
    {62, 2.5, 0.0, 0.0},
    {63, 2.2, 0.0, 0.0},
};
/// @brief acosh: arguments and the correctly rounded result, from mpmath 1.4.1.
static constexpr MathVector kAcoshVectors[] = {
    {1, 0x5A630A8278D56BA6ull, 0xA2DF464EBAD93686ull, 0ull, 0ull, 0x3848F2BA59371C79ull, 0x4DC63907CF9B05C0ull},
    {1, 0x5CD3483795E57AAEull, 0xBC559F92E9447DD9ull, 0ull, 0ull, 0x3AAAADDE3E3723E8ull, 0x8A7D318EBDF19346ull},
    {1, 0x5FFC8E1579470E5Aull, 0xE52D1FFFA8FBBDB8ull, 0ull, 0ull, 0x3D954474C2CAAE4Aull, 0x41E2AB27919AF23Dull},
    {2, 0x4D106068CC55B499ull, 0xE358088041D8E5FFull, 0ull, 0ull, 0x30D3C89CC4786904ull, 0xA7DA1B083257E9FBull},
    {2, 0x6BC3A823566233C6ull, 0x3D191B461B23F11Cull, 0ull, 0ull, 0x3C4E34FE8BBC0057ull, 0x2147C09B62E15760ull},
    {2, 0x3F8F2407FD70B1F7ull, 0x5CF401B4F82B969Aull, 0ull, 0ull, 0x29E30EF78849C51Bull, 0xFB69D8E342CF0AB7ull},
    {4, 0x4E67E2B937003951ull, 0xB9A1BA584D7939D1ull, 0ull, 0ull, 0x17C8AF54D09CCBA2ull, 0x1381FCDF6C77E7F3ull},
    {4, 0x6AE152E678EE962Dull, 0x4F7DCFFAEEA02E51ull, 0ull, 0ull, 0x1A45A7396CE92407ull, 0x42BDE1527AD2B327ull},
    {4, 0x0E786E6667062BC8ull, 0x6579C8269FFEE7FFull, 0ull, 0ull, 0x0997161FF88D6C48ull, 0xD865B895D92212CBull},
    {16, 0x64BECD98F40E175Aull, 0x6C8E15E8055B0BD2ull, 0ull, 0ull, 0x0005C5A401D25D13ull, 0xA8B69CDE5CCA3901ull},
    {16, 0x717EAD82055E0C3Eull, 0xD6C4D8512862BE1Dull, 0ull, 0ull, 0x0005D4E4B13A40C0ull, 0xC14F5854BE5FE323ull},
    {16, 0x06E82EB880E2D8DDull, 0x713CBB71DF5FD3B2ull, 0ull, 0ull, 0x00046E979C732532ull, 0xE305B64D5C6D3B6Dull},
    {16, 0x0005BAF3E5DE4114ull, 0x7FE6134695D08EFDull, 0ull, 0ull, 0x000190A8C68D5139ull, 0xF47F6BC94CEF1588ull},
    {16, 0x048DCA0E885907EEull, 0xF74166D5F7CAB26Eull, 0ull, 0ull, 0x00043945F971B0B1ull, 0xEDEEB0A69E94C4DAull},
    {16, 0x000B913114BED291ull, 0x189D0C90FDB71187ull, 0ull, 0ull, 0x0001EAC0661918ADull, 0x3AF0DC7A325751CCull},
    {32, 0x000000017FC924AFull, 0x60B31F766E809A9Bull, 0ull, 0ull, 0x00000000E18E4C6Cull, 0xE24C478A5D387A57ull},
    {32, 0x06B08CAA7AD3D3A5ull, 0x18B3802DA9C1E32Dull, 0ull, 0ull, 0x00000009F610CCEEull, 0xA76F236D9F34FFCDull},
    {32, 0x00CDF997E8847F7Full, 0xFA8F8465099433A8ull, 0ull, 0ull, 0x00000008E6F6B537ull, 0x8675ACD12FE74B2Full},
    {62, 0x44AA320CE0BB46EDull, 0xA3FE2848BD7265E0ull, 0ull, 0ull, 0x0000000000000056ull, 0x17490C54DF252F8Bull},
    {62, 0x71955B79393BBF58ull, 0xB4E18F493FC562F4ull, 0ull, 0ull, 0x0000000000000057ull, 0x18F9A7BC3B9D87E7ull},
    {62, 0x0000000000003B70ull, 0x5B202A6BBA747FB3ull, 0ull, 0ull, 0x0000000000000013ull, 0x429FE5DB0549DD89ull},
    {63, 0x4BFFEE6B5B81EA8Cull, 0x7BCADDCA3AC849F2ull, 0ull, 0ull, 0x000000000000002Bull, 0xD7120BA5DB3B2EBFull},
    {63, 0x000000000000036Cull, 0x09F5649B5F7690DFull, 0ull, 0ull, 0x0000000000000007ull, 0x77F35FE6A91467E8ull},
    {63, 0x00000B017009171Dull, 0xF4BD7375AC14B7C7ull, 0ull, 0ull, 0x000000000000001Eull, 0xD143C0A1A16D6B31ull},
    {63, 0x318B7AD1CFC3C671ull, 0x34C2E39A7934D6D9ull, 0ull, 0ull, 0x000000000000002Bull, 0x698A5C58769D51BFull},
    {63, 0x000000000091E357ull, 0x9D6D83BC6AB6D3B0ull, 0ull, 0ull, 0x0000000000000010ull, 0xC42EF779C2D64DE1ull},
    {63, 0x0014F52385BF6D00ull, 0x079E40CB1E1B6052ull, 0ull, 0ull, 0x0000000000000025ull, 0x01B7DF5B84AE8E89ull},
};
/// @brief acosh: largest error accepted, per I (1.2 times the largest of 2000 random samples + 0.6).
static constexpr MathBound kAcoshBounds[] = {
    {1, 4.1, 0.0, 0.0},
    {2, 3.1, 0.0, 0.0},
    {4, 3.0, 0.0, 0.0},
    {16, 3.1, 0.0, 0.0},
    {32, 2.1, 0.0, 0.0},
    {62, 2.6, 0.0, 0.0},
    {63, 2.0, 0.0, 0.0},
};
/// @brief atanh: arguments and the correctly rounded result, from mpmath 1.4.1.
static constexpr MathVector kAtanhVectors[] = {
    {1, 0x081AD8BD1D21CB1Eull, 0x90E33F8D81AD696Cull, 0ull, 0ull, 0x08260BD85FD49950ull, 0xD474B55E53A42014ull},
    {1, 0xFEE81339D9AC68FFull, 0xFE5BB437866ADCD4ull, 0ull, 0ull, 0xFEE80C40842FF884ull, 0x7C8D5558CC1085F8ull},
    {1, 0xFFFFFFFFFFFFFFFFull, 0xFFFFFE60ACE56FB9ull, 0ull, 0ull, 0xFFFFFFFFFFFFFFFFull, 0xFFFFFE60ACE56FB9ull},
    {2, 0x0000000051EAF77Dull, 0xFCE197B094CD2A94ull, 0ull, 0ull, 0x0000000051EAF77Dull, 0xFCE197B147BE8676ull},
    {2, 0x007AECA921F73EC0ull, 0xB37C6ABC8935850Bull, 0ull, 0ull, 0x007AEF05D994286Cull, 0x25AD46DDE019DBF5ull},
    {2, 0xFFFFFFFF7C6B9EDAull, 0x9B2182F4B51B6B4Bull, 0ull, 0ull, 0xFFFFFFFF7C6B9EDAull, 0x9B2182F1CF8CD2DDull},
    {4, 0x000000000092433Full, 0xE71752B2BB4C8AE0ull, 0ull, 0ull, 0x000000000092433Full, 0xE71752B2BB4CCA89ull},
    {4, 0xFFFC00F9C860B7BFull, 0x6E1CB504C3D6EAA2ull, 0ull, 0ull, 0xFFFC00F97349B88Aull, 0xD30FC3DA115C2334ull},
    {4, 0xFFFFFFFFFFFFFFFFull, 0xFFFFFFFFFFFFFF8Full, 0ull, 0ull, 0xFFFFFFFFFFFFFFFFull, 0xFFFFFFFFFFFFFF8Full},
    {16, 0x000000000220F6E9ull, 0xD8DEF764B2FC899Bull, 0ull, 0ull, 0x000000000220F6E9ull, 0xD8DF04417C1B33A4ull},
    {16, 0x0000000000007F8Eull, 0xC74DB4A53DDCF137ull, 0ull, 0ull, 0x0000000000007F8Eull, 0xC74DB4A53DDD1B71ull},
    {16, 0x0000000000000000ull, 0x000000010B7CBCE8ull, 0ull, 0ull, 0x0000000000000000ull, 0x000000010B7CBCE8ull},
    {16, 0x0000000000000000ull, 0x000000000001671Cull, 0ull, 0ull, 0x0000000000000000ull, 0x000000000001671Cull},
    {16, 0xFFFFFFFFFFFFFFFFull, 0x702233964ED5D856ull, 0ull, 0ull, 0xFFFFFFFFFFFFFFFFull, 0x702233964ED5D856ull},
    {16, 0x0000000000000000ull, 0x273DBA04F0C2D10Dull, 0ull, 0ull, 0x0000000000000000ull, 0x273DBA04F0C2D10Dull},
    {32, 0x0000000000000002ull, 0x11052D2AFB1C4E3Eull, 0ull, 0ull, 0x0000000000000002ull, 0x11052D2AFB1C4E4Aull},
    {32, 0x0000000000003E13ull, 0x018921BEF9A1B994ull, 0ull, 0ull, 0x0000000000003E13ull, 0x0189269CBB36B5B6ull},
    {32, 0xFFFFFFFFFFFFFFFFull, 0xFFFBCF2A4FC5AD1Dull, 0ull, 0ull, 0xFFFFFFFFFFFFFFFFull, 0xFFFBCF2A4FC5AD1Dull},
    {62, 0x0000000000000000ull, 0x366F2B53EAE0CCF3ull, 0ull, 0ull, 0x0000000000000000ull, 0x36A4085A2001C8AAull},
    {62, 0x0000000000000000ull, 0x00000D2BFE152BD9ull, 0ull, 0ull, 0x0000000000000000ull, 0x00000D2BFE152BDAull},
    {62, 0x0000000000000000ull, 0x0003D1507E3DB0F9ull, 0ull, 0ull, 0x0000000000000000ull, 0x0003D1507E4253F3ull},
    {63, 0x0000000000000000ull, 0x0E7FD508F072F4D5ull, 0ull, 0ull, 0x0000000000000000ull, 0x0E83CF10C8E41E49ull},
    {63, 0x0000000000000000ull, 0x042715B77ECC6253ull, 0ull, 0ull, 0x0000000000000000ull, 0x04272D9752DEE283ull},
    {63, 0x0000000000000000ull, 0x0000000000001E1Aull, 0ull, 0ull, 0x0000000000000000ull, 0x0000000000001E1Aull},
    {63, 0x0000000000000000ull, 0x0000000000DC73A7ull, 0ull, 0ull, 0x0000000000000000ull, 0x0000000000DC73A7ull},
    {63, 0x0000000000000000ull, 0x0000000000000B36ull, 0ull, 0ull, 0x0000000000000000ull, 0x0000000000000B36ull},
    {63, 0xFFFFFFFFFFFFFFFFull, 0xFFFFFFFFFD3704FCull, 0ull, 0ull, 0xFFFFFFFFFFFFFFFFull, 0xFFFFFFFFFD3704FCull},
};
/// @brief atanh: largest error accepted, per I (1.2 times the largest of 2000 random samples + 0.6).
static constexpr MathBound kAtanhBounds[] = {
    {1, 2.1, 0.0, 0.0},
    {2, 2.2, 0.0, 0.0},
    {4, 2.2, 0.0, 0.0},
    {16, 2.1, 0.0, 0.0},
    {32, 2.1, 0.0, 0.0},
    {62, 2.2, 0.0, 0.0},
    {63, 2.0, 0.0, 0.0},
};
/// @brief sin big: arguments and the correctly rounded result, from mpmath 1.4.1.
static constexpr MathVector kSinBigVectors[] = {
    {32, 0x000001F837CF4A72ull, 0xA68C979CFB2C8846ull, 0ull, 0ull, 0x0000000001F2FB35ull, 0x3134109C312AEEEDull},
    {32, 0xFFFFFE07C28CEA53ull, 0x110742EB8475E74Bull, 0ull, 0ull, 0x0000000003B0A992ull, 0x74AC7F012F923685ull},
    {32, 0xFFFFFF92380374AEull, 0x4F14EC1287472FFEull, 0ull, 0ull, 0x000000002BC818D6ull, 0xDDE1559C59D53203ull},
    {62, 0x00000000001F221Full, 0x05AEBB25F53C83E9ull, 0ull, 0ull, 0xFFFFFFFFFFFFFFFFull, 0xB1378DB55F0A5817ull},
    {62, 0x00000000001F7AACull, 0x66475213917AE3C9ull, 0ull, 0ull, 0xFFFFFFFFFFFFFFFFull, 0x5910A57640CCECB5ull},
    {62, 0x00000000000DAC75ull, 0x4E6ECF74A144AC31ull, 0ull, 0ull, 0xFFFFFFFFFFFFFFFEull, 0x0000C83EDE247FE9ull},
    {63, 0xFFFFFFFFFFF00B12ull, 0x24993A08AAAAE4FBull, 0ull, 0ull, 0x0000000000000000ull, 0x16E5C9E101E4206Eull},
    {63, 0xFFFFFFFFFFF03C96ull, 0x59E17C15D86D7277ull, 0ull, 0ull, 0x0000000000000000ull, 0x0755FFF56CFF7E31ull},
    {63, 0xFFFFFFFFFFFFBCDAull, 0x3FCC89D2DC0344E7ull, 0ull, 0ull, 0x0000000000000000ull, 0xDD51DCEAE96FA4B0ull},
    {63, 0xFFFFFFFFFFFFFE21ull, 0xA0AE9D8E9011F127ull, 0ull, 0ull, 0xFFFFFFFFFFFFFFFFull, 0x3FA12DDC9FED5BEFull},
    {63, 0xFFFFFFFFFFFFFE4Dull, 0x9A57E54851936FDCull, 0ull, 0ull, 0xFFFFFFFFFFFFFFFFull, 0x3E71393416F74B37ull},
    {63, 0x00000000000134E3ull, 0xCC7C6215CC9F7391ull, 0ull, 0ull, 0x0000000000000000ull, 0xF144FCAC35B22F56ull},
};
/// @brief sin big: largest error accepted, per I (1.2 times the largest of 2000 random samples + 0.6).
static constexpr MathBound kSinBigBounds[] = {
    {32, 2.5, 0.0, 1.0 / 3.0},
    {62, 2.5, 0.0, 1.0 / 3.0},
    {63, 2.5, 0.0, 1.0 / 3.0},
};
/// @brief cos big: arguments and the correctly rounded result, from mpmath 1.4.1.
static constexpr MathVector kCosBigVectors[] = {
    {32, 0xFFFFFE022AA4A6C7ull, 0xB59072EC05BCCD6Cull, 0ull, 0ull, 0xFFFFFFFFE40FF5B5ull, 0x887F0F6FE964CF12ull},
    {32, 0x000001EF8AD6C64Dull, 0xEDFDB5347795A218ull, 0ull, 0ull, 0xFFFFFFFFF4C692B0ull, 0x887D2F57788445CFull},
    {32, 0x000001E400E24081ull, 0xE6734EAEDDFF0F35ull, 0ull, 0ull, 0x0000000076137B4Aull, 0xCA6BBEA3D31645F7ull},
    {62, 0x00000000001F932Eull, 0xD425363DF3270059ull, 0ull, 0ull, 0xFFFFFFFFFFFFFFFFull, 0xFEDF5D1B385F0D00ull},
    {62, 0xFFFFFFFFFFE0D4E8ull, 0x5B43D7B6776E02C7ull, 0ull, 0ull, 0xFFFFFFFFFFFFFFFFull, 0xCEF29E86C72DF199ull},
    {62, 0xFFFFFFFFFFFE298Full, 0xCA096CEB84DAA630ull, 0ull, 0ull, 0xFFFFFFFFFFFFFFFFull, 0x45FF92D83825B8DAull},
    {63, 0x00000000000FD746ull, 0x74A4F01E4B9C4F40ull, 0ull, 0ull, 0x0000000000000000ull, 0x2ABA376166AD0B2Bull},
    {63, 0x00000000000FBD35ull, 0x5DC0299030BB54B0ull, 0ull, 0ull, 0xFFFFFFFFFFFFFFFFull, 0xD26F40561684112Bull},
    {63, 0xFFFFFFFFFFFFC2C7ull, 0xCF39C84D5DC34FCDull, 0ull, 0ull, 0xFFFFFFFFFFFFFFFFull, 0xA6E12A2038A107F4ull},
    {63, 0xFFFFFFFFFFFFFFA3ull, 0x898C352EFCEBFC96ull, 0ull, 0ull, 0xFFFFFFFFFFFFFFFFull, 0xC9903C09A3684F1Full},
    {63, 0xFFFFFFFFFFFFA6B4ull, 0x7C2E603AD779F32Dull, 0ull, 0ull, 0x0000000000000000ull, 0x47BCB30BCF0E4EA7ull},
    {63, 0x000000000000579Eull, 0x0DE343A5AC32A067ull, 0ull, 0ull, 0x0000000000000000ull, 0x9BA3C9FF4A68E418ull},
};
/// @brief cos big: largest error accepted, per I (1.2 times the largest of 2000 random samples + 0.6).
static constexpr MathBound kCosBigBounds[] = {
    {32, 2.5, 0.0, 1.0 / 3.0},
    {62, 2.5, 0.0, 1.0 / 3.0},
    {63, 2.5, 0.0, 1.0 / 3.0},
};
/// @brief cordic sin: arguments and the correctly rounded result, from mpmath 1.4.1.
static constexpr MathVector kCordicSinVectors[] = {
    {16, 0x000000004D270F6Bull, 0x4491D891CABE892Full, 0ull, 0ull, 0x000000004D270F6Bull, 0x3FE5E05DD593068Aull},
    {16, 0x0000000000000000ull, 0x001A8077DE054A88ull, 0ull, 0ull, 0x0000000000000000ull, 0x001A8077DE054A88ull},
    {16, 0xFFFFFFFFFFFFFFFFull, 0xFFFFFFEBB881C36Eull, 0ull, 0ull, 0xFFFFFFFFFFFFFFFFull, 0xFFFFFFEBB881C36Eull},
    {16, 0x00000000014B02B9ull, 0x09147AA200389BD7ull, 0ull, 0ull, 0x00000000014B02B9ull, 0x091479310FD14DC2ull},
    {16, 0x0000000000000001ull, 0xA4C1BF7826D57EE2ull, 0ull, 0ull, 0x0000000000000001ull, 0xA4C1BF7826D57EE2ull},
    {16, 0xFFFFFFFFFFFFFFFFull, 0xFFFFFFF983FAD7B0ull, 0ull, 0ull, 0xFFFFFFFFFFFFFFFFull, 0xFFFFFFF983FAD7B0ull},
    {32, 0xFFFFFFFFFFFFFFFFull, 0xFFFE5C697F653CA7ull, 0ull, 0ull, 0xFFFFFFFFFFFFFFFFull, 0xFFFE5C697F653CA7ull},
    {32, 0xFFFFFFFFFFFFFFF3ull, 0xEC37E8678054F62Bull, 0ull, 0ull, 0xFFFFFFFFFFFFFFF3ull, 0xEC37E8678054FAC1ull},
    {32, 0x0000000000000000ull, 0x00000001B51D7C01ull, 0ull, 0ull, 0x0000000000000000ull, 0x00000001B51D7C01ull},
    {62, 0xFFFFFFFFFFFFFFFFull, 0xFFFBB87CDE18CB35ull, 0ull, 0ull, 0xFFFFFFFFFFFFFFFFull, 0xFFFBB87CDE1C0F1Eull},
    {62, 0xFFFFFFFFFFFFFFFFull, 0xFFFFABAA4693AE4Bull, 0ull, 0ull, 0xFFFFFFFFFFFFFFFFull, 0xFFFFABAA4693AEADull},
    {62, 0xFFFFFFFFFFFFFFFFull, 0x8BC3BF9B404CA4D9ull, 0ull, 0ull, 0xFFFFFFFFFFFFFFFFull, 0x8CC2B13389F587EAull},
    {63, 0xFFFFFFFFFFFFFFFFull, 0xFF0A9012E78D5052ull, 0ull, 0ull, 0xFFFFFFFFFFFFFFFFull, 0xFF0A90388131A921ull},
    {63, 0xFFFFFFFFFFFFFFFFull, 0xFFFFFF69D39C2CCEull, 0ull, 0ull, 0xFFFFFFFFFFFFFFFFull, 0xFFFFFF69D39C2CCEull},
    {63, 0xFFFFFFFFFFFFFFFFull, 0xFFC2725D6C752003ull, 0ull, 0ull, 0xFFFFFFFFFFFFFFFFull, 0xFFC2725E0449FF85ull},
    {63, 0x0000000000000000ull, 0x0000000000007063ull, 0ull, 0ull, 0x0000000000000000ull, 0x0000000000007063ull},
    {63, 0xFFFFFFFFFFFFFFFFull, 0xFFFFFFF814F97159ull, 0ull, 0ull, 0xFFFFFFFFFFFFFFFFull, 0xFFFFFFF814F97159ull},
    {63, 0xFFFFFFFFFFFFFFFFull, 0xFFFFFD3392B4B4A6ull, 0ull, 0ull, 0xFFFFFFFFFFFFFFFFull, 0xFFFFFD3392B4B4A6ull},
};
/// @brief cordic sin: largest error accepted, per I (1.2 times the largest of 2000 random samples + 0.6).
static constexpr MathBound kCordicSinBounds[] = {
    {16, 337744410924777.6, 0.0, 0.0},
    {32, 5146167441.0, 0.0, 0.0},
    {62, 10.9, 0.0, 0.0},
    {63, 8.4, 0.0, 0.0},
};
/// @brief cordic cos: arguments and the correctly rounded result, from mpmath 1.4.1.
static constexpr MathVector kCordicCosVectors[] = {
    {16, 0xFFFF932166C5892Dull, 0xE0ED836DE49F90F7ull, 0ull, 0ull, 0x0000546CF3ECF0CBull, 0xA89352DFDD475D08ull},
    {16, 0xFFFF8D8A92978D19ull, 0xFB1E6F4C51266BE8ull, 0ull, 0ull, 0x000050253A1E4090ull, 0x6A713D23AA7C1D97ull},
    {16, 0xFFFFD2F86164EF37ull, 0x50ABB6AACBE1133Full, 0ull, 0ull, 0x00007829257911B4ull, 0xF5343F7212569CFEull},
    {16, 0x000000003E8F1055ull, 0x0B0B8E8E670EF96Bull, 0ull, 0ull, 0x00007FFFFFFFF0B6ull, 0x64238EE7398CADA7ull},
    {16, 0xFFFFFFFFFFFFFFF7ull, 0xA3FAABF2B9946BBDull, 0ull, 0ull, 0x00007FFFFFFFFFFFull, 0xFFFFFFFFFFBA1E97ull},
    {16, 0xFFFFFFFFFFFFFFFFull, 0xFFFFFFFFEE4A87F8ull, 0ull, 0ull, 0x0000800000000000ull, 0x0000000000000000ull},
    {32, 0x0000000068DE1A13ull, 0x1906F28F91D830F9ull, 0ull, 0ull, 0x00000000576450F5ull, 0x94AFD89875BB5FCAull},
    {32, 0x000000007FDAEAC5ull, 0x01CFA85CE3EF6853ull, 0ull, 0ull, 0x000000004547D19Eull, 0x8F014D71DF400038ull},
    {32, 0xFFFFFFFFFFFFFFFFull, 0xFFFFFFFFFFFFFCD4ull, 0ull, 0ull, 0x0000000080000000ull, 0x0000000000000000ull},
    {62, 0x0000000000000001ull, 0xAE282E1762C6F7D6ull, 0ull, 0ull, 0x0000000000000001ull, 0x55AF045B4C548C2Eull},
    {62, 0x0000000000000000ull, 0xB941C53F0A05E8A1ull, 0ull, 0ull, 0x0000000000000001ull, 0xDED92B843C00109Eull},
    {62, 0xFFFFFFFFFFFFFFFFull, 0xFFFFF9CBBD7088C9ull, 0ull, 0ull, 0x0000000000000001ull, 0xFFFFFFFFFFF6608Eull},
    {63, 0xFFFFFFFFFFFFFFFFull, 0xED76377F0407C4E8ull, 0ull, 0ull, 0x0000000000000000ull, 0xFF543E0B0B31AFD3ull},
    {63, 0x0000000000000000ull, 0x001662CA3D9DF5D3ull, 0ull, 0ull, 0x0000000000000000ull, 0xFFFFFF056F8F21D3ull},
    {63, 0x0000000000000000ull, 0x0000CD1FA2CC3D70ull, 0ull, 0ull, 0x0000000000000000ull, 0xFFFFFFFFADD228AEull},
    {63, 0x0000000000000000ull, 0x2F604F7A90513322ull, 0ull, 0ull, 0x0000000000000000ull, 0xFBA0F24FB3E9A548ull},
    {63, 0xFFFFFFFFFFFFFFFFull, 0xFFFFFFFFFFFFA8D3ull, 0ull, 0ull, 0x0000000000000001ull, 0x0000000000000000ull},
    {63, 0xFFFFFFFFFFFFFFFFull, 0x913A73450B08FA87ull, 0ull, 0ull, 0x0000000000000000ull, 0xE867F2634DDB0367ull},
};
/// @brief cordic cos: largest error accepted, per I (1.2 times the largest of 2000 random samples + 0.6).
static constexpr MathBound kCordicCosBounds[] = {
    {16, 250919841397305.6, 0.0, 0.0},
    {32, 3567210928.9, 0.0, 0.0},
    {62, 8.1, 0.0, 0.0},
    {63, 6.3, 0.0, 0.0},
};

/// @brief Checks the vectors of one width, when the function exists at that width.
template <int32_t I, int32_t MinI, typename Fn, size_t N, size_t M>
static void CheckMathVectorsAt(const char* name, const MathVector (&vectors)[N], const MathBound (&bounds)[M], Fn fn)
{
    if constexpr (I >= MinI) {
        using fp = fixed_point128<I>;
        const MathBound* bound = nullptr;
        for (const auto& b : bounds) {
            if (b.i == I)
                bound = &b;
        }
        if (bound == nullptr)
            return;
        int32_t checked = 0;
        for (const auto& v : vectors) {
            if (v.i != I)
                continue;
            const fp x(v.argLow, v.argHigh), y(v.arg2Low, v.arg2High), expected(v.expectedLow, v.expectedHigh);
            fp result;
            if constexpr (std::is_invocable_v<Fn, fp, fp>) {
                result = fn(x, y);
            } else {
                result = fn(x);
            }
            const double allowed = bound->absUlps + bound->relUlps * std::max(1.0, fabs(static_cast<double>(expected))) + bound->perArgUlps * fabs(static_cast<double>(x));
            EXPECT_LE(UlpDistance(result, expected), allowed)
                << name << ", I=" << I << ", x=" << HexOf(x) << ", y=" << HexOf(y) << ", expected " << HexOf(expected) << ", got " << HexOf(result);
            ++checked;
        }
        EXPECT_GT(checked, 0) << name << ", I=" << I << ": a bound without vectors";
    }
}
/// @brief Checks the vectors of a math function at every width it has a bound for.
template <int32_t MinI, typename Fn, size_t N, size_t M>
static void CheckMathVectors(const char* name, const MathVector (&vectors)[N], const MathBound (&bounds)[M], Fn fn)
{
    CheckMathVectorsAt<1, MinI>(name, vectors, bounds, fn);
    CheckMathVectorsAt<2, MinI>(name, vectors, bounds, fn);
    CheckMathVectorsAt<4, MinI>(name, vectors, bounds, fn);
    CheckMathVectorsAt<16, MinI>(name, vectors, bounds, fn);
    CheckMathVectorsAt<32, MinI>(name, vectors, bounds, fn);
    CheckMathVectorsAt<62, MinI>(name, vectors, bounds, fn);
    CheckMathVectorsAt<63, MinI>(name, vectors, bounds, fn);
}
TEST(fixed_point128, log)
{
    CheckMathVectors<1>("log", kLogVectors, kLogBounds, [](auto x) { return log(x); });
}
TEST(fixed_point128, log2)
{
    CheckMathVectors<1>("log2", kLog2Vectors, kLog2Bounds, [](auto x) { return log2(x); });
}
TEST(fixed_point128, log10)
{
    CheckMathVectors<1>("log10", kLog10Vectors, kLog10Bounds, [](auto x) { return log10(x); });
}
TEST(fixed_point128, log1p)
{
    CheckMathVectors<1>("log1p", kLog1pVectors, kLog1pBounds, [](auto x) { return log1p(x); });
}
TEST(fixed_point128, sin)
{
    CheckMathVectors<4>("sin", kSinVectors, kSinBounds, [](auto x) { return sin(x); });
}
TEST(fixed_point128, cos)
{
    CheckMathVectors<4>("cos", kCosVectors, kCosBounds, [](auto x) { return cos(x); });
}
// Beyond pi/4 the reduction x - n pi/2 carries n times the rounding error of pi/2: the bound grows by a
// third of an ulp per unit of |x|, as documented.
TEST(fixed_point128, SinAndCosOfLargeArguments)
{
    CheckMathVectors<4>("sin", kSinBigVectors, kSinBigBounds, [](auto x) { return sin(x); });
    CheckMathVectors<4>("cos", kCosBigVectors, kCosBigBounds, [](auto x) { return cos(x); });
}
TEST(fixed_point128, tan)
{
    CheckMathVectors<4>("tan", kTanVectors, kTanBounds, [](auto x) { return tan(x); });
}
TEST(fixed_point128, asin)
{
    CheckMathVectors<4>("asin", kAsinVectors, kAsinBounds, [](auto x) { return asin(x); });
}
TEST(fixed_point128, acos)
{
    CheckMathVectors<4>("acos", kAcosVectors, kAcosBounds, [](auto x) { return acos(x); });
}
TEST(fixed_point128, atan)
{
    CheckMathVectors<4>("atan", kAtanVectors, kAtanBounds, [](auto x) { return atan(x); });
}
TEST(fixed_point128, atan2)
{
    CheckMathVectors<4>("atan2", kAtan2Vectors, kAtan2Bounds, [](auto y, auto x) { return atan2(y, x); });
}
TEST(fixed_point128, exp)
{
    CheckMathVectors<2>("exp", kExpVectors, kExpBounds, [](auto x) { return exp(x); });
}
TEST(fixed_point128, exp2)
{
    CheckMathVectors<2>("exp2", kExp2Vectors, kExp2Bounds, [](auto x) { return exp2(x); });
}
TEST(fixed_point128, expm1)
{
    CheckMathVectors<2>("expm1", kExpm1Vectors, kExpm1Bounds, [](auto x) { return expm1(x); });
}
TEST(fixed_point128, pow)
{
    CheckMathVectors<2>("pow", kPowVectors, kPowBounds, [](auto x, auto y) { return pow(x, y); });
}
TEST(fixed_point128, sinh)
{
    CheckMathVectors<4>("sinh", kSinhVectors, kSinhBounds, [](auto x) { return sinh(x); });
}
TEST(fixed_point128, asinh)
{
    CheckMathVectors<1>("asinh", kAsinhVectors, kAsinhBounds, [](auto x) { return asinh(x); });
}
TEST(fixed_point128, cosh)
{
    CheckMathVectors<4>("cosh", kCoshVectors, kCoshBounds, [](auto x) { return cosh(x); });
}
TEST(fixed_point128, acosh)
{
    CheckMathVectors<1>("acosh", kAcoshVectors, kAcoshBounds, [](auto x) { return acosh(x); });
}
TEST(fixed_point128, tanh)
{
    CheckMathVectors<2>("tanh", kTanhVectors, kTanhBounds, [](auto x) { return tanh(x); });
}
TEST(fixed_point128, atanh)
{
    CheckMathVectors<1>("atanh", kAtanhVectors, kAtanhBounds, [](auto x) { return atanh(x); });
}
// _sincos_cordic() is public through argument dependent lookup and had no test. Its 64 steps reach an
// absolute accuracy of about 2^-63, which is what the bounds express at every width.
TEST(fixed_point128, CordicSinCos)
{
    CheckMathVectors<16>("cordic sin", kCordicSinVectors, kCordicSinBounds, [](auto x) {
        decltype(x) s, c;
        _sincos_cordic(x, s, c, true);
        return s;
    });
    CheckMathVectors<16>("cordic cos", kCordicCosVectors, kCordicCosBounds, [](auto x) {
        decltype(x) s, c;
        _sincos_cordic(x, s, c, true);
        return c;
    });
}

/// @brief Expects a result within ulps of a correctly rounded raw reference.
template <int32_t I> static void ExpectWithinUlps(const fixed_point128<I>& result, uint64_t expectedHigh, uint64_t expectedLow, double ulps, const char* label)
{
    const fixed_point128<I> expected(expectedLow, expectedHigh);
    EXPECT_LE(UlpDistance(result, expected), ulps) << label << ", I=" << I << ", expected " << HexOf(expected) << ", got " << HexOf(result);
}
/// @brief Exact results at every width a function exists in.
template <int32_t I> static void CheckExactMathResults()
{
    using fp = fixed_point128<I>;
    const fp one = fp::one();
    EXPECT_TRUE(sqrt(one) == one) << "I=" << I;
    EXPECT_TRUE(log(one) == fp(0)) << "I=" << I;
    EXPECT_TRUE(log2(one) == fp(0)) << "I=" << I;
    EXPECT_TRUE(log10(one) == fp(0)) << "I=" << I;
    EXPECT_TRUE(log2(one >> 1) == -one) << "I=" << I;
    EXPECT_TRUE(asinh(fp(0)) == fp(0)) << "I=" << I;
    EXPECT_TRUE(acosh(one) == fp(0)) << "I=" << I;
    EXPECT_TRUE(atanh(fp(0)) == fp(0)) << "I=" << I;
    if constexpr (I >= 2) {
        EXPECT_TRUE(exp(fp(0)) == one) << "I=" << I;
        EXPECT_TRUE(expm1(fp(0)) == fp(0)) << "I=" << I;
        EXPECT_TRUE(tanh(fp(0)) == fp(0)) << "I=" << I;
        EXPECT_TRUE(pow(fp(1.5), fp(0)) == one) << "I=" << I;
        // pow(1, y) used to return the log of 1, zero, for every y
        EXPECT_TRUE(pow(one, fp(1.5)) == one) << "I=" << I;
        EXPECT_TRUE(pow(one, -fp(1.5)) == one) << "I=" << I;
        EXPECT_TRUE(pow(one, fp(0)) == one) << "I=" << I;
        EXPECT_TRUE(log2(fp(2)) == one) << "I=" << I;
    }
    if constexpr (I >= 7) {
        // exp() of a negative argument whose positive counterpart is far beyond the range: from I = 7 on,
        // e^-(2^I) is below half the last place and the result underflows to zero
        EXPECT_TRUE(exp(-std::numeric_limits<fp>::max()) == fp(0)) << "I=" << I;
        EXPECT_TRUE(exp(std::numeric_limits<fp>::lowest()) == fp(0)) << "I=" << I;
    }
    if constexpr (I >= 4) {
        // cos(0) is the cosine series at zero, exactly one; it was 7026 ulp short at I = 63
        EXPECT_TRUE(cos(fp(0)) == one) << "I=" << I;
        EXPECT_TRUE(sin(fp(0)) == fp(0)) << "I=" << I;
        EXPECT_TRUE(cos1(fp(0)) == one) << "I=" << I;
        EXPECT_TRUE(sin1(fp(0)) == fp(0)) << "I=" << I;
        EXPECT_TRUE(atan(fp(0)) == fp(0)) << "I=" << I;
        EXPECT_TRUE(tan(fp(0)) == fp(0)) << "I=" << I;
        EXPECT_TRUE(asin(fp(0)) == fp(0)) << "I=" << I;
        EXPECT_TRUE(sinh(fp(0)) == fp(0)) << "I=" << I;
        EXPECT_TRUE(cosh(fp(0)) == one) << "I=" << I;
    }
}
TEST(fixed_point128, MathIdentities)
{
    CheckExactMathResults<1>();
    CheckExactMathResults<2>();
    CheckExactMathResults<4>();
    CheckExactMathResults<16>();
    CheckExactMathResults<32>();
    CheckExactMathResults<63>();
    // exact powers of two take a shortcut
    for (int32_t k = -100; k <= 15; ++k) {
        const fixed_point128<16> power = (k >= 0) ? fixed_point128<16>::one() << k : fixed_point128<16>::one() >> -k;
        EXPECT_TRUE(log2(power) == fixed_point128<16>(k)) << "k=" << k;
    }
    // the documented error of an exact power: not exact, but within pow()'s bound
    ExpectWithinUlps(pow(fixed_point128<16>(2), fixed_point128<16>(3)), 0x0004000000000000ull, 0ull, 0.6 + 5.8 * 8, "pow(2, 3)");
}
// The axes are handled before the ratio is formed. And steeper than the diagonal, the angle used to
// take its sign from atan(x / y), which lost it when x / y rounded to zero: atan2(1000, -epsilon) came
// out as 3 pi / 2. Values from mpmath.
TEST(fixed_point128, Atan2OnTheAxes)
{
    using fp = fixed_point128<16>;
    EXPECT_TRUE(atan2(fp(0), fp(-1)) == fp::pi());
    EXPECT_TRUE(atan2(fp(1), fp(0)) == fp::half_pi());
    EXPECT_TRUE(atan2(fp(-1), fp(0)) == -fp::half_pi());
    EXPECT_TRUE(atan2(fp(0), fp(0)) == fp(0));
    EXPECT_TRUE(atan2(fp(0), fp(1)) == fp(0));
    ExpectWithinUlps(atan2(fp(1000), -fp::epsilon()), 0x0000C90FDAA22168ull, 0xC234C4C6628B80DCull, 2.5, "atan2(1000, -eps)");
    ExpectWithinUlps(atan2(fp(-1000), fp::epsilon()), 0xFFFF36F0255DDE97ull, 0x3DCB3B399D747F24ull, 2.5, "atan2(-1000, eps)");
    ExpectWithinUlps(atan2(fp(1000), fp::epsilon()), 0x0000C90FDAA22168ull, 0xC234C4C6628B80DCull, 2.5, "atan2(1000, eps)");
}
// exp() of an integer argument skips the series, and a negative argument used to be computed as
// 1 / e^|x|, which overflowed long before e^-|x| left the grid: exp(-12) was garbage at
// fixed_point128<16>. tanh() divided e^x - e^-x by e^x + e^-x and overflowed with them. Values from
// mpmath, within the documented bounds.
TEST(fixed_point128, ExpOfIntegersAndNegativeArguments)
{
    using fp16 = fixed_point128<16>;
    ExpectWithinUlps(exp(fp16(1)), 0x00015BF0A8B14576ull, 0x95355FB8AC404E7Aull, 0.6 + 5.0 * 2.72, "exp(1)");
    ExpectWithinUlps(exp(fp16(-12)), 0x00000000338A9876ull, 0x877915E321D83680ull, 0.6 + 5.0, "exp(-12)");
    ExpectWithinUlps(exp(fp16(-40)), 0x0000000000000000ull, 0x00272F2493363448ull, 0.6 + 5.0, "exp(-40)");
    ExpectWithinUlps(exp(fixed_point128<4>(-5)), 0x000DCC9FF006A9BBull, 0x0E22E1248A4D2DCBull, 0.6 + 7.5, "exp(-5)");
    ExpectWithinUlps(tanh(fp16(12)), 0x00007FFFFFFFD67Dull, 0xF0E03475D6D756CBull, 3.1, "tanh(12)");
    ExpectWithinUlps(tanh(fp16(-12)), 0xFFFF800000002982ull, 0x0F1FCB8A2928A935ull, 3.1, "tanh(-12)");
    ExpectWithinUlps(tanh(fixed_point128<2>(3)), 0x1FD77D111A0AFFE2ull, 0x0AF0A8B0B7FB5E58ull, 4.2, "tanh(3)");
    EXPECT_TRUE(tanh(std::numeric_limits<fp16>::max()) == fp16::one());
    EXPECT_TRUE(tanh(std::numeric_limits<fp16>::lowest()) == -fp16::one());
}
// atanh() formed (1 + x) / (1 - x), which left the range as x approached one, and a wrapped, negative
// value made log() throw out of the noexcept function: the process terminated. asinh() and acosh()
// added x to a root of x^2 +- 1 and overflowed the same way for large arguments, and acosh() lost
// hundreds of millions of ulp just above one, where x^2 - 1 cancels. Values from mpmath.
TEST(fixed_point128, InverseHyperbolicFunctionsAtTheEdges)
{
    using fp16 = fixed_point128<16>;
    ExpectWithinUlps(atanh(fp16::one() - (fp16::one() >> 20)), 0x0003A396FBD50D79ull, 0x3EC5B940C3336A0Aull, 2.1, "atanh(1 - 2^-20)");
    ExpectWithinUlps(atanh(fixed_point128<4>(0.95)), 0x0EA77CB44406EA1Eull, 0x8D6647213D2D192Bull, 2.2, "atanh(0.95)");
    EXPECT_TRUE(atanh(std::numeric_limits<fp16>::lowest()) == fp16(0));
    EXPECT_TRUE(atanh(fp16(1)) == fp16(0));
    ExpectWithinUlps(asinh(fp16(40000)), 0x0005A51792E165E2ull, 0xD2C16D450E7792BCull, 2.6, "asinh(40000)");
    ExpectWithinUlps(acosh(fp16(40000)), 0x0005A51792E0BA16ull, 0x5BAFE8E33F799313ull, 3.1, "acosh(40000)");
    ExpectWithinUlps(asinh(std::numeric_limits<fp16>::lowest()), 0xFFFA1BB63445689Cull, 0x75CBD7F08973D7C4ull, 2.6, "asinh(lowest)");
    ExpectWithinUlps(acosh(fixed_point128<1>::one() + fixed_point128<1>(1ull << 26, 0)), 0x00000000000016A0ull, 0x9E667F3BCC908B30ull, 4.1, "acosh(1 + 2^-100)");
}

/**
 * @brief Checks that sqr(x) is bit identical to x * x for a single template parameter.
 *
 * square() accumulates the same 256 bit product as operator*=, it only skips the
 * redundant second cross multiply, so the results must match exactly - including
 * rounding and the two's complement correction the sign needs. Raw bit patterns are used
 * to cover the full value range instead of going through double, which cannot represent
 * every fixed_point128 value.
 *
 * The product of two arbitrary 128 bit patterns overflows the range far more often than not,
 * and an overflowed square can land anywhere, sign included. What is being checked here is that
 * the two spellings agree bit for bit, whatever they produce; SqrOfSmallValuesIsPositive covers
 * the sign for the arguments whose square actually fits.
 *
 * @tparam I Number of integer bits passed to fixed_point128.
 */
template <int32_t I> static void CheckSqrMatchesMultiply()
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const uint64_t low = get_uint64_random();
        const uint64_t high = get_uint64_random();
        const fixed_point128<I> x(low, high);
        const fixed_point128<I> viaMultiply = x * x;
        const fixed_point128<I> viaSqr = sqr(x);

        // operator== compares high and low, so this is a bit exact comparison
        EXPECT_TRUE(viaSqr == viaMultiply) << "I=" << I << ", low=" << low << ", high=" << high
                                           << ", x*x=" << (std::string)viaMultiply << ", sqr(x)=" << (std::string)viaSqr;
    }
}
TEST(fixed_point128, SqrMatchesMultiply)
{
    CheckSqrMatchesMultiply<1>();
    CheckSqrMatchesMultiply<10>();
    CheckSqrMatchesMultiply<20>();
    CheckSqrMatchesMultiply<40>();
    CheckSqrMatchesMultiply<63>();
}
/**
 * @brief A square that stays inside the range is never negative, whatever the sign of the argument.
 *
 * The arguments are kept below 2^(I/2) so the product cannot overflow, which is the condition the
 * claim rests on. Both signs are covered, the negative one being the case the two's complement
 * correction in square() exists for.
 */
TEST(fixed_point128, SqrOfSmallValuesIsPositive)
{
    using fp = fixed_point128<32>;
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        // |value| < 2^16, so the square stays below 2^32 and inside the range of fixed_point128<32>
        const double value = get_double_random(-20, 15);
        const fp x(value);
        EXPECT_FALSE(sqr(x).is_negative()) << "value=" << value;
        EXPECT_FALSE(sqr(-x).is_negative()) << "value=" << -value;
        EXPECT_TRUE(sqr(x) == sqr(-x)) << "value=" << value;
        Raw128 expected {};
        ASSERT_TRUE(RefProduct<32>(RawOf(x), RawOf(x), expected)) << "value=" << value;
        EXPECT_TRUE(RefEqual(RawOf(sqr(x)), expected)) << "value=" << value;
    }
}
TEST(fixed_point128, SqrEdgeCases)
{
    using fp = fixed_point128<20>;
    using lim = std::numeric_limits<fp>;
    const fp values[] = {fp(0),  fp::epsilon(), -fp::epsilon(), fp::one(),      -fp::one(),   fp::half(),
                         fp::pi(), -fp::pi(),   fp::e(),        fp::golden_ratio(), lim::max(), lim::lowest(),
                         fp(0ull, ~0ull), fp(~0ull, ~0ull)};
    for (const auto& x : values) {
        EXPECT_TRUE(sqr(x) == x * x) << "x=" << (std::string)x;
    }
    // squaring zero stays zero, and a square that fits in the range is never negative
    EXPECT_TRUE(sqr(fp(0)) == fp(0));
    EXPECT_FALSE(sqr(-fp::one()).is_negative());
    EXPECT_TRUE(sqr(-fp::one()) == fp::one());
}

/**********************************************************************
 * fixed_point128 regression tests
 *
 * Each test below pins down a defect that was fixed.
 ***********************************************************************/

// The sign used to live in its own field, which made a zero carrying a sign a distinct bit pattern
// from a plain zero; every comparison tested the sign field first, so such a value compared unequal
// to zero and smaller than it, and trunc, ceil, round, modf and copysign all produced one. The two's
// complement representation cannot express a signed zero at all, which is what these now pin down.
TEST(fixed_point128, NoNegativeZero)
{
    typedef fixed_point128<32> fp;
    const fp zero(0);

    EXPECT_TRUE(trunc(fp(-0.5)) == zero);
    EXPECT_FALSE(trunc(fp(-0.5)) < zero);
    EXPECT_FALSE(trunc(fp(-0.5)).is_negative());

    EXPECT_TRUE(ceil(fp(-0.5)) == zero);
    EXPECT_FALSE(ceil(fp(-0.5)).is_negative());

    EXPECT_TRUE(round(fp(-0.4)) == zero);
    EXPECT_FALSE(round(fp(-0.4)).is_negative());

    EXPECT_TRUE(copysign(zero, fp(-1)) == zero);
    EXPECT_FALSE(copysign(zero, fp(-1)).is_negative());

    fp int_part;
    const fp frac = modf(fp(-0.25), &int_part);
    EXPECT_TRUE(int_part == zero);
    EXPECT_FALSE(int_part.is_negative());
    EXPECT_TRUE(frac.is_negative());  // the fraction keeps the sign, it is non zero

    // floor was already correct, keep it covered
    EXPECT_TRUE(floor(fp(-0.5)) == fp(-1));
    EXPECT_TRUE(trunc(fp(-1.5)) == fp(-1));

    // There is one zero and it is its own negation. operator== compares the raw QWORDs, so these
    // are bit pattern comparisons rather than value ones.
    EXPECT_TRUE(-zero == zero);
    EXPECT_TRUE(zero.is_positive() && !zero.is_negative());
    EXPECT_TRUE((fp(3) - fp(3)) == zero);
    EXPECT_FALSE((fp(3) - fp(3)).is_negative());
    EXPECT_TRUE(fp(-0.5) + fp(0.5) == zero);
}
// The conversions truncated the bits that did not fit the mantissa, and a value whose surviving
// fraction bits happened to be all ones was rounded up even when it was exactly representable,
// so (2^24-1)/2^23 came out as 2. Converting the 128 bit magnitude and scaling by a power of two
// with ldexp is exact, which makes it a correctly rounded reference.
TEST(fixed_point128, ConversionToFloatingPointIsCorrectlyRounded)
{
    typedef fixed_point128<32> fp;
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const uint64_t l = get_uint64_random();
        // shifted by at least one, so the pair is a magnitude with the sign bit clear
        const uint64_t h = get_uint64_random() >> (1 + get_uint32_random() % 40);
        const bool negative = (get_uint32_random() & 1) != 0;
        const fp v = (negative) ? -fp(l, h) : fp(l, h);
        if (!v)
            continue;
        const double mag = static_cast<double>(uint128_t(l, h));
        const double ref = negative ? -ldexp(mag, -fp::F) : ldexp(mag, -fp::F);
        EXPECT_EQ(static_cast<double>(v), ref) << "low=" << l << ", high=" << h;
        EXPECT_EQ(static_cast<float>(v), static_cast<float>(ref)) << "low=" << l << ", high=" << h;
    }
}
TEST(fixed_point128, ConversionToFloatingPointEdgeCases)
{
    typedef fixed_point128<32> fp;
    // (2^24-1)/2^23 is exactly representable as a float and must not be rounded up to 2
    const fp exact_flt = fp(0xFFFFFFull) / fp(0x800000ull);
    EXPECT_LT(static_cast<float>(exact_flt), 2.0f);
    EXPECT_EQ(static_cast<float>(exact_flt), 16777215.0f / 8388608.0f);

    // the same trap one mantissa wider, for double
    typedef fixed_point128<63> fp63;
    const fp63 exact_dbl = fp63(0x1FFFFFFFFFFFFFull) / fp63(0x10000000000000ull);
    EXPECT_LT(static_cast<double>(exact_dbl), 2.0);

    // powers of two convert exactly in both directions
    for (int e = -60; e <= 20; ++e) {
        fp v(1);
        if (e > 0) {
            v <<= e;
        } else if (e < 0) {
            v >>= -e;
        }
        EXPECT_EQ(static_cast<double>(v), ::pow(2.0, e)) << "e=" << e;
        EXPECT_EQ(static_cast<float>(v), static_cast<float>(::pow(2.0, e))) << "e=" << e;
    }
    // The widest fraction is 126 bits, so its last place is 2^-126 - exactly the smallest normal
    // float. Nothing this type can hold goes denormal, and nothing must become zero either.
    const fixed_point128<1> tiny(1, 0);
    EXPECT_GT(static_cast<double>(tiny), 0.0);
    EXPECT_GT(static_cast<float>(tiny), 0.0f);
    EXPECT_EQ(static_cast<float>(tiny), ::ldexpf(1.0f, -126));
    EXPECT_EQ(static_cast<float>(-tiny), -::ldexpf(1.0f, -126));
    EXPECT_EQ(static_cast<float>(fp(0)), 0.0f);
    EXPECT_EQ(static_cast<long double>(fp(-2.5)), -2.5L);
}
// The double constructor truncates what falls below the grid, and saturates what is beyond the
// range - NaN and the infinities included - to the extreme of the sign of the input.
TEST(fixed_point128, ConstructorFromDoubleSaturatesAndTruncates)
{
    using fp = fixed_point128<63>;
    using lim = std::numeric_limits<fp>;
    EXPECT_TRUE(fp(std::numeric_limits<double>::quiet_NaN()) == lim::max());
    EXPECT_TRUE(fp(std::numeric_limits<double>::infinity()) == lim::max());
    EXPECT_TRUE(fp(-std::numeric_limits<double>::infinity()) == lim::lowest());
    EXPECT_TRUE(fp(std::ldexp(1.0, 63)) == lim::max());
    EXPECT_TRUE(fp(-std::ldexp(1.0, 63)) == lim::lowest());
    EXPECT_TRUE(fp(-0.0) == fp(0) && fp(-0.0).is_positive());
    // below the last place: three quarters of it, and 2^-60 plus a bit 52 places further down
    EXPECT_TRUE(fp(0.75 * std::ldexp(1.0, -64)) == fp(0));
    EXPECT_TRUE(fp(std::ldexp(1.0 + std::ldexp(1.0, -52), -60)) == fp(16, 0));
    EXPECT_TRUE(fp(-std::ldexp(1.0 + std::ldexp(1.0, -52), -60)) == -fp(16, 0));
    EXPECT_TRUE(fixed_point128<1>(std::ldexp(1.0, -127)) == fixed_point128<1>(0));
}
// div_64bit shortcuts a numerator smaller than the divisor and returns without writing the
// quotient, so dividing by a scalar left the value completely unchanged instead of producing zero.
TEST(fixed_point128, DivideSmallValueByScalar)
{
    typedef fixed_point128<32> fp;
    // the raw 128 bit form of 2^-90 is 64, which is smaller than the divisor below
    const fp tiny = fp(1) >> 90;
    EXPECT_TRUE(tiny);  // the value itself is representable

    fp v = tiny;
    v /= 1000ull;
    EXPECT_TRUE(v == fp(0)) << "expected zero, got " << (std::string)v;

    // the fixed_point128 operand overload has always produced zero here, the two must agree
    fp w = tiny;
    w /= fp(1000);
    EXPECT_TRUE(v == w);

    // a normal division still works, and rounds
    fp u = fp(1) >> 40;
    u /= 1000ull;
    EXPECT_TRUE(RefIsIntegerQuotient(RawOf(fp(1) >> 40), 1000, false, RawOf(u)));
    // 1000 raw units divided by 1999 is just over a half and rounds up to one
    EXPECT_TRUE(fp(1000, 0) / 1999ull == fp(1, 0));
}
// These are noexcept and used to reach a throwing path, which terminates the process instead.
// Values from mpmath, within the documented bounds.
TEST(fixed_point128, NoexceptFunctionsDoNotTerminate)
{
    typedef fixed_point128<16> fp;

    // pow forwards to log, which rejects zero
    EXPECT_TRUE(pow(fp(0), fp(2)) == fp(0));
    EXPECT_TRUE(pow(fp(0), fp(0)) == fp(1));  // the CRT defines pow(0, 0) as 1
    EXPECT_TRUE(pow(fp(-1), fp(2)) == fp(0));

    // acos at the ends of its domain starts where the derivative vanishes and sin(res) is zero
    EXPECT_TRUE(acos(fp(1)) == fp(0));
    ExpectWithinUlps(acos(fp(-1)), 0x0001921FB54442D1ull, 0x8469898CC51701B8ull, 2.8, "acos(-1)");
    ExpectWithinUlps(asin(fp(1)), 0x0000C90FDAA22168ull, 0xC234C4C6628B80DCull, 2.7, "asin(1)");
    EXPECT_TRUE(asin(fp(2)) == fp(0));   // outside the domain: documented zero
    EXPECT_TRUE(acos(fp(-2)) == fp(0));
    EXPECT_TRUE(acosh(fp(0.5)) == fp(0));

    // tan's range reduction lands exactly on a zero cosine at the poles
    EXPECT_TRUE(tan(fp::half_pi()) == fp(0));

    // reciprocal is documented to return zero rather than saturate
    EXPECT_TRUE(reciprocal(fp(0)) == fp(0));
}
// Newton's method on sin() stalls near pi/2, where the derivative it divides by vanishes: asin(1) and
// acos(-1) were good to double precision only, 4e16 and 8e16 ulp off at fixed_point128<16>, which the
// 1e-12 tolerance above used to let through. Arguments beyond one half now reduce through
// sqrt((1 - x) / 2). Values from mpmath, within the documented bounds.
TEST(fixed_point128, AsinAcosAtTheEndsOfTheDomain)
{
    using fp = fixed_point128<16>;
    ExpectWithinUlps(asin(fp(1) - fp::epsilon()), 0x0000C90FDAA22168ull, 0xC134C4C6628B80DCull, 2.7, "asin(1 - eps)");
    ExpectWithinUlps(asin(-fp(1)), 0xFFFF36F0255DDE97ull, 0x3DCB3B399D747F24ull, 2.7, "asin(-1)");
    ExpectWithinUlps(acos(-fp(1) + fp::epsilon()), 0x0001921FB54442D1ull, 0x8369898CC51701B8ull, 2.8, "acos(-1 + eps)");
    ExpectWithinUlps(acos(fp(1) - (fp(1) >> 100)), 0x0000000000000000ull, 0x2D413CCCFE779921ull, 2.8, "acos(1 - 2^-100)");
    EXPECT_TRUE(acos(fp(1)) == fp(0));
    using fp63 = fixed_point128<63>;
    ExpectWithinUlps(asin(fp63(-1) + (fp63(1) >> 60)), 0xFFFFFFFFFFFFFFFEull, 0x6DE04AC165561536ull, 3.4, "asin(-1 + 2^-60)");
}
// log1p forwards to log, which throws on a non positive argument. It used to be noexcept, which
// turned that into a call to std::terminate.
TEST(fixed_point128, Log1pPropagatesDomainError)
{
    typedef fixed_point128<16> fp;
    EXPECT_THROW((void)log1p(fp(-1)), std::domain_error);
    EXPECT_THROW((void)log1p(fp(-2)), std::domain_error);
    ExpectWithinUlps(log1p(fp(1)), 0x000058B90BFBE8E7ull, 0xBCD5E4F1D9CC01F9ull, 1.9, "log1p(1)");
    // the log family itself keeps throwing
    EXPECT_THROW((void)log(fp(0)), std::domain_error);
    EXPECT_THROW((void)log2(fp(-1)), std::domain_error);
    EXPECT_THROW((void)log10(fp(0)), std::domain_error);
}
// A leading sign or a high byte handed straight to the ctype functions is outside the domain they
// accept, and only a literal space was skipped as leading white space.
TEST(fixed_point128, ConstructorFromStringEdgeCases)
{
    typedef fixed_point128<32> fp;
    EXPECT_TRUE(fp("\t 2.5") == fp(2.5));
    EXPECT_TRUE(fp("  -2.5") == fp(-2.5));
    EXPECT_TRUE(fp("+2.5") == fp(2.5));
    const char high_byte[] = {'1', '2', static_cast<char>(0xB5), 0};
    EXPECT_TRUE(fp(high_byte) == fp(12));
}
// The decimal digits of fixed_point128<1> to <3> come from a product that needs four integer bits, and
// were rounded rather than truncated: 1.75 printed as 1.85, and 0.96 as "0.:6::...", the character after
// '9'. The digits are truncated from the full product now. Expected text from Python's decimal module.
TEST(fixed_point128, NarrowInstantiationsPrintTheirDigits)
{
    EXPECT_EQ(static_cast<std::string>(fixed_point128<1>(0.25)), "0.25");
    EXPECT_EQ(static_cast<std::string>(fixed_point128<1>(1.75)), "1.75");
    EXPECT_EQ(static_cast<std::string>(fixed_point128<1>(0.96)), "0.95999999999999996447286321199499070644");
    EXPECT_EQ(static_cast<std::string>(fixed_point128<1>(0.1)), "0.10000000000000000555111512312578270212");
    // 38 places for I = 2 and 3 as for I = 1: ceil(F * log10(2)) is 38 for F = 124 to 126, the fewest
    // that tell neighbouring values apart. F / 3.3 gave 37, one short.
    EXPECT_EQ(static_cast<std::string>(fixed_point128<2>(-0.96)), "-0.95999999999999996447286321199499070644");
    EXPECT_EQ(static_cast<std::string>(fixed_point128<3>(0.96)), "0.95999999999999996447286321199499070644");
}

/**********************************************************************
 * Two's complement representation
 *
 * The value is a 128 bit two's complement integer scaled by 2^-F, with the sign in the MSB of the
 * upper QWORD. These pin down what that costs and what it guarantees.
 ***********************************************************************/

/// @brief The raw constructor and get_components read the same 128 bits the object stores.
template <int32_t I> static void CheckRawRepresentation()
{
    using fp = fixed_point128<I>;

    // 1.0 is a single bit at the radix point, and its negation is the two's complement of that
    const fp one = fp::one();
    uint64_t low = 0, high = 0;
    one.get_components(low, high);
    EXPECT_EQ(low, 0ull) << "I=" << I;
    EXPECT_EQ(high, 1ull << (fp::F - 64)) << "I=" << I;

    (-one).get_components(low, high);
    EXPECT_EQ(low, 0ull) << "I=" << I;
    EXPECT_EQ(high, 0ull - (1ull << (fp::F - 64))) << "I=" << I;
    EXPECT_TRUE(fp(low, high) == -one) << "I=" << I;

    // a value round trips through its components
    srand(RANDOM_SEED);
    for (auto i = 0u; i < 100; ++i) {
        const fp x(get_uint64_random(), get_uint64_random());
        uint64_t l = 0, h = 0;
        x.get_components(l, h);
        EXPECT_TRUE(fp(l, h) == x) << "I=" << I;
        // the sign of the value is the top bit of the upper QWORD, nothing else
        EXPECT_EQ(x.is_negative(), (h >> 63) == 1) << "I=" << I;
        EXPECT_NE(x.is_negative(), x.is_positive()) << "I=" << I;
        // the bit queries read the same bits, at run time as well as in a constant expression
        for (uint32_t bit : {0u, 1u, 63u, 64u, 65u, 127u}) {
            const uint64_t expected = (bit < 64) ? (l >> bit) & 1 : (h >> (bit - 64)) & 1;
            EXPECT_EQ(static_cast<uint64_t>(x.get_bit(bit)), expected) << "I=" << I << ", bit=" << bit;
        }
        // ~x is -x - epsilon, and unary + is the value itself
        EXPECT_TRUE(~x == -x - fp::epsilon()) << "I=" << I;
        EXPECT_TRUE(+x == x) << "I=" << I;
    }
}
TEST(fixed_point128, RawRepresentation)
{
    CheckRawRepresentation<1>();
    CheckRawRepresentation<20>();
    CheckRawRepresentation<32>();
    CheckRawRepresentation<63>();
}
/// @brief The range is [-2^I, 2^I - epsilon], and it wraps at both ends.
template <int32_t I> static void CheckRangeIsAsymmetric()
{
    using fp = fixed_point128<I>;
    using lim = std::numeric_limits<fp>;

    EXPECT_DOUBLE_EQ(static_cast<double>(lim::lowest()), -::ldexp(1.0, I)) << "I=" << I;
    EXPECT_TRUE(RefEqual(RawOf(lim::max()), Raw128 {~0ull, ~0ull >> 1})) << "I=" << I;
    // max() is one step short of the magnitude of lowest()
    EXPECT_TRUE(lim::lowest() + lim::max() == -lim::epsilon()) << "I=" << I;

    // lowest() has no positive counterpart, so negating it and taking its magnitude both wrap
    EXPECT_TRUE(-lim::lowest() == lim::lowest()) << "I=" << I;
    EXPECT_TRUE(fabs(lim::lowest()) == lim::lowest()) << "I=" << I;
    EXPECT_TRUE(lim::lowest().is_negative()) << "I=" << I;

    // the ends of the range wrap into each other, the way the builtin integer types do
    EXPECT_TRUE(lim::max() + lim::epsilon() == lim::lowest()) << "I=" << I;
    EXPECT_TRUE(lim::lowest() - lim::epsilon() == lim::max()) << "I=" << I;

    // and they compare in the right order
    EXPECT_TRUE(lim::lowest() < fp(-1) && fp(-1) < fp(0) && fp(0) < fp(1) && fp(1) < lim::max()) << "I=" << I;
    EXPECT_TRUE(lim::lowest() < lim::max()) << "I=" << I;
}
TEST(fixed_point128, RangeIsAsymmetric)
{
    CheckRangeIsAsymmetric<1>();
    CheckRangeIsAsymmetric<16>();
    CheckRangeIsAsymmetric<32>();
    CheckRangeIsAsymmetric<63>();
}
// The multiply and the divide both work on magnitudes and put the sign back, which is the part of
// the representation with the most room to get a sign wrong. Every combination is checked against
// the exact references, and the identity x * -y == -(x * y) has to hold exactly.
TEST(fixed_point128, SignedMultiplyAndDivide)
{
    typedef fixed_point128<32> fp;
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const double a = get_double_random(-20, 10);
        const double b = get_double_random(-20, 10);

        const fp fa(a), fb(b);
        EXPECT_TRUE(fa * -fb == -(fa * fb)) << "a=" << a << ", b=" << b;
        EXPECT_TRUE(-fa * fb == -(fa * fb)) << "a=" << a << ", b=" << b;
        EXPECT_TRUE(-fa * -fb == fa * fb) << "a=" << a << ", b=" << b;
        Raw128 product {};
        ASSERT_TRUE(RefProduct<32>(RawOf(-fa), RawOf(fb), product)) << "a=" << a << ", b=" << b;
        EXPECT_TRUE(RefEqual(RawOf(-fa * fb), product)) << "a=" << a << ", b=" << b;

        EXPECT_TRUE(fa / -fb == -(fa / fb)) << "a=" << a << ", b=" << b;
        EXPECT_TRUE(-fa / -fb == fa / fb) << "a=" << a << ", b=" << b;
        EXPECT_TRUE(RefIsQuotient<32>(RawOf(-fa), RawOf(fb), RawOf(-fa / fb))) << "a=" << a << ", b=" << b;

        // the integer overloads take the same path with the sign handled separately
        const int32_t n = (abs(get_int32_random()) % 1000) + 1;
        EXPECT_TRUE(fa * -n == -(fa * n)) << "a=" << a << ", n=" << n;
        EXPECT_TRUE(-fa / n == -(fa / n)) << "a=" << a << ", n=" << n;
        EXPECT_TRUE(fa / -n == -(fa / n)) << "a=" << a << ", n=" << n;
    }
}
// lowest() is the one numerator whose magnitude is 2^127, which reads as a negative value when it
// is held in the type. The division has a path for it, and this is what makes sure it is taken.
TEST(fixed_point128, DivideTheMostNegativeValue)
{
    typedef fixed_point128<32> fp;
    using lim = std::numeric_limits<fp>;
    const Raw128 lowest = RawOf(lim::lowest());

    // a power of two divisor, which is turned into a shift
    EXPECT_TRUE(lim::lowest() / fp(4) == fp(-1073741824));
    EXPECT_TRUE(lim::lowest() / fp(-4) == fp(1073741824));
    // an integer divisor, which goes through div_64bit, and neither of those, which is the case that
    // cannot use the reciprocal
    for (const double divisor : {3.0, -3.0, 2.5, -2.5, 7.75}) {
        EXPECT_TRUE(RefIsQuotient<32>(lowest, RawOf(fp(divisor)), RawOf(lim::lowest() / fp(divisor)))) << "divisor=" << divisor;
    }
    // and as a divisor, where its magnitude is a power of two
    EXPECT_DOUBLE_EQ(static_cast<double>(fp(1) / lim::lowest()), -1.0 / ::ldexp(1.0, 32));
}
// Shifting right is arithmetic: the sign bit is replicated, so a negative value stays negative and
// the shift keeps dividing by a power of two.
TEST(fixed_point128, ArithmeticShiftRight)
{
    typedef fixed_point128<32> fp;

    EXPECT_DOUBLE_EQ(static_cast<double>(fp(-8) >> 1), -4.0);
    EXPECT_DOUBLE_EQ(static_cast<double>(fp(-8) >> 3), -1.0);
    EXPECT_DOUBLE_EQ(static_cast<double>(fp(-1) >> 1), -0.5);
    EXPECT_TRUE((fp(-8) >> 2).is_negative());

    // The shift rounds to nearest, so a value shifted past its last place lands on zero rather than
    // on the sign bits it is filled with. Only a value whose leading one survives the shift can
    // round to the last place, which is what the second pair below does.
    EXPECT_TRUE((fp(-8) >> 127) == fp(0));
    EXPECT_TRUE((fp(8) >> 127) == fp(0));
    EXPECT_TRUE((std::numeric_limits<fp>::lowest() >> 127) == -std::numeric_limits<fp>::epsilon());
    EXPECT_TRUE((std::numeric_limits<fp>::max() >> 127) == std::numeric_limits<fp>::epsilon());

    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT; ++i) {
        const double value = get_double_random(-20, 10);
        const fp x(value);
        for (int32_t shift = 1; shift <= 8; ++shift) {
            EXPECT_TRUE(RefEqual(RawOf(x >> shift), RefShiftRightRounded(RawOf(x), shift))) << "value=" << value << ", shift=" << shift;
            EXPECT_EQ((x >> shift).is_negative(), value < 0 && (x >> shift) != fp(0)) << "value=" << value;
        }
    }
}

/**
 * @brief The instantiations that are too narrow for a constant or a function reject it at compile
 *        time, through a static_assert in the body of each one.
 *
 * A violation is a compile error, so it cannot be exercised from here. What this pins down is the
 * other half of the claim - that at exactly the documented minimum the value is representable and
 * the function works - which is what makes the bound the right one rather than merely a safe one.
 * Expected values from mpmath, within each function's documented bound.
 */
TEST(fixed_point128, MinimumIntegerBitsAreTight)
{
    // pi and e need 2 integer bits, pi2 needs 3, and each is exact at its minimum
    EXPECT_DOUBLE_EQ(static_cast<double>(fixed_point128<2>::pi()), 3.14159265358979323846);
    EXPECT_DOUBLE_EQ(static_cast<double>(fixed_point128<2>::e()), 2.71828182845904523536);
    EXPECT_DOUBLE_EQ(static_cast<double>(fixed_point128<3>::pi2()), 6.28318530717958647692);

    // the constants below two are usable in the narrowest instantiation there is
    EXPECT_DOUBLE_EQ(static_cast<double>(fixed_point128<1>::half_pi()), 1.57079632679489661923);
    EXPECT_DOUBLE_EQ(static_cast<double>(fixed_point128<1>::sqrt_2()), 1.41421356237309504880);
    EXPECT_DOUBLE_EQ(static_cast<double>(fixed_point128<1>::golden_ratio()), 1.61803398874989484820);
    EXPECT_DOUBLE_EQ(static_cast<double>(fixed_point128<1>::one()), 1.0);
    EXPECT_DOUBLE_EQ(static_cast<double>(fixed_point128<1>::half()), 0.5);

    // the exponential family at its minimum of 2 integer bits
    const fixed_point128<2> half2 = fixed_point128<2>::half();
    ExpectWithinUlps(exp(half2), 0x34C2531C3C0D3792ull, 0xE5BFDF56DBE67F36ull, 0.6 + 6.9 * 1.65, "exp(0.5)");
    ExpectWithinUlps(exp2(half2), 0x2D413CCCFE779921ull, 0x165F626CDD52AFA8ull, 0.6 + 4.6 * 1.42, "exp2(0.5)");
    ExpectWithinUlps(expm1(half2), 0x14C2531C3C0D3792ull, 0xE5BFDF56DBE67F36ull, 0.6 + 6.9, "expm1(0.5)");
    ExpectWithinUlps(pow(fixed_point128<2>(1.5), half2), 0x27311C2812425CFAull, 0x06432286434AAC8Eull, 0.6 + 4.1 * 1.23, "pow(1.5, 0.5)");
    ExpectWithinUlps(tanh(half2), 0x0EC9A9EBAB4579B2ull, 0x8947710E32D8469Dull, 4.2, "tanh(0.5)");

    // The trigonometric group at its minimum of 4. The Maclaurin series used to carry powers of its
    // argument in the instantiation itself, which left the range of fixed_point128<4> long before the
    // series converged and held it to about three decimal digits; the series now runs in a wider
    // fraction and fixed_point128<4> is good to a few ulp.
    const fixed_point128<4> half4 = fixed_point128<4>::half();
    ExpectWithinUlps(sin(half4), 0x03D5DD0E8960BDFDull, 0x0EC9786C9FB06F6Dull, 3.6, "sin(0.5)");
    ExpectWithinUlps(cos(half4), 0x07054A0196DF53E7ull, 0x6DEEECED17D7D6CCull, 2.8, "cos(0.5)");
    ExpectWithinUlps(atan(half4), 0x03B58CE0AC3769EDull, 0x15BF9117B2F0EA34ull, 3.3, "atan(0.5)");
    ExpectWithinUlps(asin(half4), 0x0430548E0B5CD961ull, 0x196ECCB83D59EB44ull, 3.9, "asin(0.5)");
    ExpectWithinUlps(acos(half4), 0x0860A91C16B9B2C2ull, 0x32DD99707AB3D689ull, 4.8, "acos(0.5)");

    // and the functions that carry no bound still work in the narrowest one
    const fixed_point128<1> x = fixed_point128<1>::half();
    ExpectWithinUlps(sqrt(x), 0x2D413CCCFE779921ull, 0x165F626CDD52AFA8ull, 0.5, "sqrt(0.5)");
    ExpectWithinUlps(log(x), 0xD3A37A020B8C2195ull, 0x0D871319FF034254ull, 2.8, "log(0.5)");
    ExpectWithinUlps(reciprocal(fixed_point128<1>(0.75)), 0x5555555555555555ull, 0x5555555555555555ull, 0.5, "reciprocal(0.75)");
    ExpectWithinUlps(asinh(x), 0x1ECC2CAEC5160994ull, 0xBE04204A968C7021ull, 3.9, "asinh(0.5)");
}

// A scalar on the left used to be ambiguous: converting it to fixed_point128 and converting the
// fixed_point128 to a builtin type are both a single user defined conversion. The scalar is the one
// that gets widened, so the fraction of the object survives instead of being truncated away.
TEST(fixed_point128, ScalarOnTheLeftHandSide)
{
    typedef fixed_point128<20> fp;
    const fp quarter = fp(1) / 4;  // 0.25

    EXPECT_DOUBLE_EQ(static_cast<double>(1 + quarter), 1.25);
    EXPECT_DOUBLE_EQ(static_cast<double>(1 - quarter), 0.75);
    EXPECT_DOUBLE_EQ(static_cast<double>(3 * quarter), 0.75);
    EXPECT_DOUBLE_EQ(static_cast<double>(1 / quarter), 4.0);
    EXPECT_DOUBLE_EQ(static_cast<double>(5 % fp(3)), 2.0);
    EXPECT_DOUBLE_EQ(static_cast<double>(2.5 + quarter), 2.75);
    EXPECT_DOUBLE_EQ(static_cast<double>(2.5 - quarter), 2.25);

    // subtracting past zero must keep the sign
    EXPECT_DOUBLE_EQ(static_cast<double>(1 - fp(3)), -2.0);
    EXPECT_DOUBLE_EQ(static_cast<double>(-1 + quarter), -0.75);

    // the bitwise forms agree with the fixed_point128 on the left ones
    EXPECT_TRUE((7 & fp(3)) == (fp(7) & 3));
    EXPECT_TRUE((7 | fp(8)) == (fp(7) | 8));
    EXPECT_TRUE((7 ^ fp(3)) == (fp(7) ^ 3));

    // the fixed_point128 on the left forms must still resolve
    EXPECT_TRUE((quarter + 1) == fp(1.25));
    EXPECT_TRUE((quarter + quarter) == fp(0.5));
    EXPECT_TRUE((quarter * quarter) == fp(0.0625));

    // The left operand is widened rather than the object being narrowed, so the result is 128 bit.
    // Narrowing instead would make (1 - quarter) an int and lose the fraction entirely.
    static_assert(std::is_same_v<decltype(1 + quarter), fp>);
    static_assert(std::is_same_v<decltype(1 ^ quarter), fp>);
}
// Widening the left operand keeps the shift at 128 bits. The count is the right operand converted
// to int32_t, which for fixed_point128 is a shift and so truncates toward zero.
TEST(fixed_point128, ScalarOnTheLeftHandSideShifts)
{
    typedef fixed_point128<20> fp;

    EXPECT_DOUBLE_EQ(static_cast<double>(1 << fp(3)), 8.0);
    EXPECT_DOUBLE_EQ(static_cast<double>(8 >> fp(3)), 1.0);
    EXPECT_DOUBLE_EQ(static_cast<double>(3 << fp(2)), 12.0);
    EXPECT_DOUBLE_EQ(static_cast<double>(1 >> fp(2)), 0.25);

    // A fractional shift count goes through operator int32_t, which truncates toward zero here, so
    // 3.75 shifts by 3. The same overload picked with the fixed_point128 on the left has to agree,
    // both reaching the count the same way. float128 rounds instead, see its matching test.
    EXPECT_TRUE((1 << fp(3.75)) == (1 << fp(3)));
    EXPECT_TRUE((1 << fp(3.75)) == (fp(1) << fp(3.75)));

    static_assert(std::is_same_v<decltype(1 << fp(3)), fp>);

    // the fixed_point128 on the left forms must still resolve
    EXPECT_TRUE((fp(4) << 1) == fp(8));
    EXPECT_TRUE((fp(4) >> 1) == fp(2));
}

/**********************************************************************
 * Compile time (constexpr) evaluation
 *
 * Each test below pairs static_asserts, which fail the build the moment one of these operations
 * stops being usable in a constant expression, with a runtime check of the same expression.
 *
 * The runtime half earns its place. A constant evaluated call does not run the same code a runtime
 * call does: the bit counting and extended arithmetic intrinsics (mulx_u64, addcarryx_u64,
 * lzcnt64) are not constant expressions, so fp128_shared.h substitutes a portable
 * implementation of each one while the compiler is evaluating. Only comparing the two results
 * shows that the substitutes agree with the hardware. opaque(), declared in gtest_shared.h, is
 * what keeps the runtime half from being constant folded back into the compile time one.
 ***********************************************************************/
TEST(fixed_point128, ConstexprConstructionAndConversion)
{
    typedef fixed_point128<32> fp;

    constexpr fp zero;
    constexpr fp from_u64(static_cast<uint64_t>(7));
    constexpr fp from_i64(static_cast<int64_t>(-7));
    constexpr fp from_u32(static_cast<uint32_t>(7));
    constexpr fp from_i32(static_cast<int32_t>(-7));
    constexpr fp from_bits(1, 0);  // the raw low/high constructor
    constexpr fp copied(from_i32);
    constexpr fp assigned = [] { fp a; a = fp(5); return a; }();

    static_assert(zero.is_zero());
    static_assert(static_cast<uint64_t>(from_u64) == 7ull);
    static_assert(static_cast<int64_t>(from_i64) == -7ll);
    static_assert(static_cast<uint32_t>(from_u32) == 7u);
    static_assert(static_cast<int32_t>(from_i32) == -7);
    static_assert(from_bits == fp::epsilon());
    static_assert(copied == from_i32);
    static_assert(static_cast<int32_t>(assigned) == 5);

    // the cross template conversions, in both directions
    constexpr fixed_point128<16> fewer_int_bits(from_u32);  // I < I2, shifts left
    constexpr fixed_point128<48> more_int_bits(from_u32);   // I > I2, shifts right
    static_assert(static_cast<int32_t>(fewer_int_bits) == 7);
    static_assert(static_cast<int32_t>(more_int_bits) == 7);

    // the exact constants, which are compile time values rather than lazily initialized statics
    static_assert(static_cast<int32_t>(fp::one()) == 1);
    static_assert(fp::half() + fp::half() == fp::one());
    static_assert(fp::epsilon() > fp(0));
    static_assert(fixed_point128<63>::half() < fixed_point128<63>::one());  // 0.5 lands in the low QWORD
    static_assert(fixed_point128<1>::half() < fixed_point128<1>::one());

    // the raw constructor spells a negative value with the sign bit of the high QWORD
    constexpr fp minus_one_from_bits(0, 0xFFFFFFFF80000000ull);  // -2^95 scaled by 2^-95
    static_assert(minus_one_from_bits == fp(-1));
    static_assert(minus_one_from_bits.is_negative());

    EXPECT_TRUE(fp(opaque(-7)) == from_i32);
    EXPECT_TRUE(fixed_point128<16>(fp(opaque(7))) == fewer_int_bits);
    EXPECT_TRUE(fixed_point128<48>(fp(opaque(7))) == more_int_bits);
}
TEST(fixed_point128, ConstexprShiftsBitwiseAndUnary)
{
    typedef fixed_point128<32> fp;

    static_assert(static_cast<int32_t>(fp(8) >> 3) == 1);
    static_assert(static_cast<int32_t>(fp(1) << 3) == 8);
    static_assert(static_cast<int32_t>(fp(1) << 70) == 0);  // overflow is silent, the value wraps
    constexpr fp shifted = [] { fp a(16); a >>= 2; a <<= 1; return a; }();
    static_assert(static_cast<int32_t>(shifted) == 8);

    static_assert(static_cast<int32_t>(fp(12) & fp(10)) == 8);
    static_assert(static_cast<int32_t>(fp(12) | fp(3)) == 15);
    static_assert(static_cast<int32_t>(fp(12) ^ fp(10)) == 6);
    static_assert(!(~fp(0)).is_zero());

    static_assert((-fp(5)).is_negative());
    static_assert((+fp(5)).is_positive());
    static_assert((-fp(0)).is_positive());  // negating zero must not produce a negative zero
    static_assert(!fp(0));
    static_assert(static_cast<bool>(fp(1)));

    EXPECT_TRUE((fp(opaque(16)) >> 2 << 1) == shifted);
    EXPECT_TRUE((fp(opaque(12)) & fp(opaque(10))) == fp(8));
    EXPECT_TRUE(-fp(opaque(0)) == fp(0));
}
TEST(fixed_point128, ConstexprComparisonsAndQueries)
{
    typedef fixed_point128<32> fp;

    static_assert(fp(1) < fp(2));
    static_assert(fp(2) > fp(1));
    static_assert(fp(-2) < fp(-1));
    static_assert(fp(1) <= 1 && fp(1) >= 1);
    static_assert(fp(1) == 1 && fp(1) != 2);
    static_assert(1 == fp(1) && 2 > fp(1));  // the overloads taking the fixed_point128 on the right

    static_assert(fp(3).is_int());
    static_assert(!(fp(1) >> 1).is_int());
    static_assert(fixed_point128<63>(3).is_int());  // the I == 63 branch, where the fraction is exactly the low QWORD
    static_assert(!(fixed_point128<63>(1) >> 1).is_int());
    static_assert(fp(-3).is_int() && !fp(-2.5).is_int());  // a negative value has its fraction bits in the same place
    static_assert(fp(0).is_zero());
    static_assert(fp(-1).is_negative() && fp(1).is_positive());
    static_assert(fp(1).get_bit(fp::F) == 1);  // 1.0 has its single set bit at the radix point
    static_assert(fp(4).get_exponent() == 2);
    static_assert((fp::one() >> 1).get_exponent() == -1);

    EXPECT_TRUE(fp(opaque(4)).get_exponent() == 2);
    EXPECT_TRUE(fp(opaque(3)).is_int());
    EXPECT_TRUE(fp(opaque(1)) < fp(opaque(2)));
}
TEST(fixed_point128, ConstexprArithmetic)
{
    typedef fixed_point128<32> fp;

    static_assert(static_cast<int32_t>(fp(3) + fp(4)) == 7);
    static_assert(static_cast<int32_t>(fp(3) - fp(4)) == -1);
    static_assert(static_cast<int32_t>(fp(-3) + fp(-4)) == -7);
    static_assert((fp(3) - fp(3)).is_zero());
    static_assert((fp(3) - fp(3)).is_positive());         // cancelling leaves the one and only zero
    static_assert(static_cast<int32_t>(fp(3) + 4) == 7);  // the generic right hand side overload

    static_assert(static_cast<int32_t>(fp(6) * fp(7)) == 42);
    static_assert(static_cast<int32_t>(fp(-6) * fp(7)) == -42);
    static_assert(fp::half() * fp(8) == fp(4));
    static_assert(static_cast<int32_t>(fp(6) * static_cast<uint64_t>(7)) == 42);  // the uint64_t specialization
    static_assert(static_cast<int32_t>(fp(-6) * -7) == 42);
    static_assert(static_cast<int32_t>(fixed_point128<63>(6) * fixed_point128<63>(7)) == 42);  // the F == 64 branch
    static_assert(fixed_point128<1>::half() * fixed_point128<1>::half() == fixed_point128<1>::half() >> 1);

    // sqr is documented to be bit identical to x * x
    static_assert(sqr(fp::half()) == fp::half() * fp::half());
    static_assert(static_cast<int32_t>(sqr(fp(-9))) == 81);
    static_assert(sqr(fp(-9)).is_positive());

    constexpr fp stepped = [] { fp a(5); ++a; a++; --a; return a; }();
    static_assert(static_cast<int32_t>(stepped) == 6);
    constexpr fp accumulated = [] { fp a; for (int32_t i = 1; i <= 10; ++i) a += fp(i); return a; }();
    static_assert(static_cast<int32_t>(accumulated) == 55);

    EXPECT_TRUE(fp(opaque(3)) + fp(opaque(4)) == fp(7));
    EXPECT_TRUE(fp(opaque(3)) - fp(opaque(4)) == fp(-1));
    EXPECT_TRUE(fp(opaque(6)) * fp(opaque(7)) == fp(42));
    EXPECT_TRUE(sqr(fp(opaque(-9))) == fp(81));
    EXPECT_TRUE(fixed_point128<63>(opaque(6)) * fixed_point128<63>(opaque(7)) == fixed_point128<63>(42));
}
TEST(fixed_point128, ConstexprMathFunctions)
{
    typedef fixed_point128<32> fp;
    constexpr fp two_and_a_half = fp(5) >> 1;

    static_assert(fabs(fp(-5)) == fp(5));
    static_assert(floor(two_and_a_half) == fp(2));
    static_assert(floor(-two_and_a_half) == fp(-3));
    static_assert(ceil(two_and_a_half) == fp(3));
    static_assert(ceil(-two_and_a_half) == fp(-2));
    static_assert(trunc(-two_and_a_half) == fp(-2));
    static_assert(round(fp::half()) == fp(1));  // the halfway value rounds away from zero
    static_assert(round(-fp::half()) == fp(-1));
    static_assert(round(-(fp::half() >> 1)).is_positive());  // -0.25 rounds to +0, never to -0
    static_assert(round(-(fp::half() >> 1)).is_zero());
    static_assert(copysign(fp(5), fp(-1)).is_negative());
    static_assert(fmin(fp(1), fp(2)) == fp(1));
    static_assert(fmax(fp(1), fp(2)) == fp(2));
    static_assert(fdim(fp(5), fp(3)) == fp(2));
    static_assert(fdim(fp(3), fp(5)) == fp(0));
    static_assert(ilogb(fp(8)) == 3);
    static_assert(lzcnt128(fp::one()) == 32);  // 1.0 sets bit F, leaving the sign and the I-1 integer bits above it clear
    static_assert(lzcnt128(-fp::one()) == 32);  // counted on the magnitude, so the sign makes no difference

    // modf splits into an integer and a fraction part, both carrying the sign of the input
    constexpr fp modf_int = [] { fp ip; (void)modf(-(fp(5) >> 1), &ip); return ip; }();
    constexpr fp modf_frac = [] { fp ip; return modf(-(fp(5) >> 1), &ip); }();
    static_assert(modf_int == fp(-2));
    static_assert(modf_frac == -fp::half());

    static_assert(log2(fp(8)) == fp(3));  // an exact power of two takes the shortcut
    static_assert(logb(fp(9)) == fp(3));
    constexpr fp log2_of_10 = log2(fp(10));  // the full fraction loop, one squaring per fraction bit
    static_assert(log2_of_10 > fp(3) && log2_of_10 < fp(4));

    EXPECT_TRUE(floor(-(fp(opaque(5)) >> 1)) == fp(-3));
    EXPECT_TRUE(round(-fp::half() * fp(opaque(1))) == fp(-1));
    EXPECT_TRUE(ilogb(fp(opaque(8))) == 3);
    // The constant evaluated log2 takes a bitwise path that truncates, the runtime one rounds to
    // nearest: they used to agree on this value only by chance. Both against mpmath: the runtime
    // result is log2(10) rounded, the compile time one log2(10) rounded down.
    static_assert(log2_of_10 == fp(0x9A3715FC9257EDFEull, 0x00000001A934F097ull));
    EXPECT_TRUE(log2(fp(opaque(10))) == fp(0x9A3715FC9257EDFFull, 0x00000001A934F097ull));
}

// The double and float conversions used to be the one thing that could not happen at compile time:
// they read the inactive member of a union, which constant evaluation rejects. Now that Double and
// Float convert with std::bit_cast, a floating point literal crosses the boundary either way.
TEST(fixed_point128, ConstexprFloatingPointConversion)
{
    typedef fixed_point128<32> fp;

    constexpr fp from_double = 3.25;
    constexpr fp negative = -2.5;
    static_assert(static_cast<double>(from_double) == 3.25);
    static_assert(static_cast<float>(from_double) == 3.25f);
    static_assert(static_cast<double>(negative) == -2.5);
    static_assert(static_cast<long double>(from_double) == 3.25L);
    static_assert(fp(0.5) == fp::half());
    static_assert(fp(1.0) == fp::one());
    static_assert(fp(0.0).is_zero());
    static_assert(fp(-0.0).is_zero() && fp(-0.0).is_positive());  // no negative zero
    static_assert(from_double + negative == fp(0.75));
    static_assert(static_cast<int32_t>(fp(3.99)) == 3);  // truncates towards zero

    // a value that needs the full fraction, not just a few bits
    constexpr fp tenth = 0.1;
    static_assert(tenth > fp(0.09) && tenth < fp(0.11));

    EXPECT_TRUE(fp(static_cast<double>(opaque(13)) / 4.0) == fp(3.25));
    EXPECT_DOUBLE_EQ(static_cast<double>(fp(opaque(13)) / 4), 3.25);
    EXPECT_TRUE(fp(0.1) == tenth);  // the runtime conversion agrees with the constant evaluated one
}
