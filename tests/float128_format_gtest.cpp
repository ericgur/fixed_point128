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
#include <format>
#include <iomanip>
#include <sstream>
#include <unordered_map>
#include <unordered_set>
#include "gtest_shared.h"

/**********************************************************************
 * float128 text conversion
 *
 * The conversions in both directions are exact: the decimal expansion of a binary128 is finite and
 * is produced in full, so the digits are correctly rounded at any precision and a string reads
 * back to the nearest representable value. That is what makes a round trip meaningful, and most of
 * what is checked here.
 ***********************************************************************/

TEST(float128_format, ShortestRoundTrip)
{
    // A value written by a short decimal literal comes back as that literal rather than as the 35
    // digits of its binary expansion.
    EXPECT_EQ(to_string(float128("0.1")), "0.1");
    EXPECT_EQ(to_string(float128("-2.5")), "-2.5");
    EXPECT_EQ(to_string(float128(1)), "1");
    EXPECT_EQ(to_string(float128("1e300")), "1e+300");
    EXPECT_EQ(to_string(float128("1e-300")), "1e-300");
    EXPECT_EQ(to_string(float128()), "0");
    EXPECT_EQ(to_string(-float128()), "-0");
    EXPECT_EQ(to_string(float128::inf()), "inf");
    EXPECT_EQ(to_string(-float128::inf()), "-inf");
    EXPECT_EQ(to_string(float128::nan()), "nan");

    // A value that needs the whole mantissa gets all 34 digits of it.
    EXPECT_EQ(to_string(float128::pi()), "3.1415926535897932384626433832795028");
}
TEST(float128_format, RoundTripIsExact)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT / 8; ++i) {
        // A random encoding across the whole normal range, subnormals included at the bottom.
        const uint64_t low = get_uint64_random();
        const uint64_t high = (get_uint64_random() & 0x0000FFFFFFFFFFFFull) | (static_cast<uint64_t>(get_uint32_random() % 0x7FFF) << 48);
        const float128 value(low, high);
        if (value.is_nan() || value.is_inf())
            continue;

        const std::string text = to_string(value);
        float128 parsed;
        const std::from_chars_result result = fp128::from_chars(text.data(), text.data() + text.size(), parsed);
        EXPECT_EQ(result.ptr, text.data() + text.size()) << "did not consume " << text;
        EXPECT_TRUE(parsed == value) << "round trip of " << text;
    }
}
TEST(float128_format, ParseIsCorrectlyRounded)
{
    // 0.1 is not representable, and the nearest binary128 to it is above it. Its exact expansion
    // is what the fixed layout has to show, digit for digit.
    EXPECT_EQ(std::format("{:.40f}", float128("0.1")), "0.1000000000000000000000000000000000048148");

    // Values that are representable exactly must come back exactly.
    EXPECT_TRUE(float128("0.5") == float128::half());
    EXPECT_TRUE(float128("0.25") == (float128::half() >> 1));
    EXPECT_TRUE(float128("1024") == ldexp(float128::one(), 10));
    EXPECT_TRUE(float128("12345678901234567890") == float128(12345678901234567890ull));

    // The extremes of the format
    using limits = std::numeric_limits<float128>;
    EXPECT_TRUE(float128(to_string(limits::max()).c_str()) == limits::max());
    EXPECT_TRUE(float128(to_string(limits::min()).c_str()) == limits::min());
    EXPECT_TRUE(float128(to_string(limits::denorm_min()).c_str()) == limits::denorm_min());

    // Beyond them the parse saturates rather than wrapping
    EXPECT_TRUE(isinf(float128("1e5000")));
    EXPECT_TRUE(float128("1e-5000").is_zero());
}
TEST(float128_format, FormatSpecification)
{
    const float128 pi = float128::pi();

    EXPECT_EQ(std::format("{:.3e}", pi), "3.142e+00");
    EXPECT_EQ(std::format("{:.3E}", pi), "3.142E+00");
    EXPECT_EQ(std::format("{:.4f}", pi), "3.1416");
    EXPECT_EQ(std::format("{:.0f}", float128("1.5")), "2");
    // a tie at the last place goes to the even digit
    EXPECT_EQ(std::format("{:.0f}", float128("0.5")), "0");
    EXPECT_EQ(std::format("{:.0f}", float128("2.5")), "2");

    // %g rules: scientific once the exponent leaves the range the precision covers
    EXPECT_EQ(std::format("{:g}", float128(1200000)), "1.2e+06");
    EXPECT_EQ(std::format("{:g}", float128("0.00001")), "1e-05");
    EXPECT_EQ(std::format("{:g}", float128("123.5")), "123.5");

    // width, fill and alignment
    EXPECT_EQ(std::format("{:>20.4f}", pi), "              3.1416");
    EXPECT_EQ(std::format("{:<20.4f}", pi), "3.1416              ");
    EXPECT_EQ(std::format("{:^20.4f}", pi), "       3.1416       ");
    EXPECT_EQ(std::format("{:*^20.4f}", pi), "*******3.1416*******");
    // zero padding goes after the sign
    EXPECT_EQ(std::format("{:020.4f}", -pi), "-00000000000003.1416");
    EXPECT_EQ(std::format("{:+.4f}", pi), "+3.1416");
    EXPECT_EQ(std::format("{: .4f}", pi), " 3.1416");
    // the alternate form keeps the point
    EXPECT_EQ(std::format("{:#.0f}", float128(3)), "3.");
    // a special value is right aligned but never zero padded
    EXPECT_EQ(std::format("{:10}", float128::inf()), "       inf");
    EXPECT_EQ(std::format("{:010}", float128::nan()), "       nan");
}
TEST(float128_format, HexadecimalFormat)
{
    EXPECT_EQ(std::format("{:a}", float128::one()), "0x1p+00");
    EXPECT_EQ(std::format("{:a}", float128(2)), "0x1p+01");
    EXPECT_EQ(std::format("{:a}", float128::half()), "0x1p-01");
    EXPECT_EQ(std::format("{:a}", float128()), "0x0p+00");
    EXPECT_EQ(std::format("{:a}", float128::pi()), "0x1.921fb54442d18469898cc51701b8p+01");
    EXPECT_EQ(std::format("{:A}", float128::pi()), "0X1.921FB54442D18469898CC51701B8P+01");
    // the fraction is rounded to the requested width
    EXPECT_EQ(std::format("{:.3a}", float128::pi()), "0x1.922p+01");
    EXPECT_EQ(std::format("{:.0a}", float128::pi()), "0x2p+01");
}
TEST(float128_format, StreamInsertion)
{
    std::ostringstream stream;
    stream << float128::one() / float128(3);
    EXPECT_EQ(stream.str(), "0.333333");  // the default precision is six significant digits

    stream.str("");
    stream << std::fixed << std::setprecision(8) << float128::pi();
    EXPECT_EQ(stream.str(), "3.14159265");

    stream.str("");
    stream << std::scientific << std::setprecision(3) << float128::pi();
    EXPECT_EQ(stream.str(), "3.142e+00");

    stream.str("");
    stream << std::uppercase << float128::pi();
    EXPECT_EQ(stream.str(), "3.142E+00");

    stream.str("");
    stream << std::nouppercase << std::hexfloat << float128::one();
    EXPECT_EQ(stream.str(), "0x1p+00");

    // width, fill and adjustment come from the stream too
    stream.str("");
    stream << std::defaultfloat << std::setw(10) << std::setfill('.') << float128(1);
    EXPECT_EQ(stream.str(), ".........1");

    stream.str("");
    stream << std::setw(10) << std::internal << std::showpos << float128(1);
    EXPECT_EQ(stream.str(), "+........1");
}
TEST(float128_format, StreamExtraction)
{
    float128 value;
    std::istringstream stream("3.14159265358979323846264338327950288 -1.5e10 bad");

    stream >> value;
    EXPECT_TRUE(value == float128::pi());
    stream >> value;
    EXPECT_TRUE(value == float128("-1.5e10"));

    // a token that is not a number sets failbit and leaves the value alone
    const float128 before = value;
    stream >> value;
    EXPECT_TRUE(stream.fail());
    EXPECT_TRUE(value == before);
}
TEST(float128_format, ToCharsAndFromChars)
{
    char buffer[64];

    auto written = fp128::to_chars(buffer, buffer + sizeof(buffer), float128("0.1"));
    EXPECT_EQ(written.ec, std::errc {});
    EXPECT_EQ(std::string(buffer, written.ptr), "0.1");

    written = fp128::to_chars(buffer, buffer + sizeof(buffer), float128::pi(), std::chars_format::scientific, 4);
    EXPECT_EQ(std::string(buffer, written.ptr), "3.1416e+00");

    written = fp128::to_chars(buffer, buffer + sizeof(buffer), float128::pi(), std::chars_format::fixed, 2);
    EXPECT_EQ(std::string(buffer, written.ptr), "3.14");

    // too small a range reports the failure rather than truncating
    written = fp128::to_chars(buffer, buffer + 2, float128::pi());
    EXPECT_EQ(written.ec, std::errc::value_too_large);

    float128 value;
    const char text[] = "  ";
    auto read = fp128::from_chars(text, text + 2, value);
    EXPECT_EQ(read.ec, std::errc::invalid_argument);

    const char number[] = "12.5rest";
    read = fp128::from_chars(number, number + sizeof(number) - 1, value);
    EXPECT_EQ(read.ec, std::errc {});
    EXPECT_EQ(read.ptr, number + 4);
    EXPECT_TRUE(value == float128("12.5"));

    // an exponent marker with no digits after it is not part of the number
    const char partial[] = "1e";
    read = fp128::from_chars(partial, partial + 2, value);
    EXPECT_EQ(read.ptr, partial + 1);
    EXPECT_TRUE(value == float128::one());

    // the named values
    const char infinity[] = "-Infinity";
    read = fp128::from_chars(infinity, infinity + sizeof(infinity) - 1, value);
    EXPECT_EQ(read.ptr, infinity + 9);
    EXPECT_TRUE(isinf(value) && value.is_negative());
}
TEST(float128_format, Hash)
{
    std::unordered_map<float128, int> counts;
    counts[float128::one()] = 1;
    counts[float128(2)] = 2;
    // the two zeros compare equal, so they have to be the same key
    counts[float128()] = 3;
    counts[-float128()] = 4;

    EXPECT_EQ(counts.size(), 3u);
    EXPECT_EQ(counts[float128()], 4);
    EXPECT_EQ(counts[float128::one()], 1);
}
TEST(float128_format, CommonType)
{
    // A mixed expression widens to float128 rather than narrowing to the builtin type.
    static_assert(std::is_same_v<std::common_type_t<float128, double>, float128>);
    static_assert(std::is_same_v<std::common_type_t<double, float128>, float128>);
    static_assert(std::is_same_v<std::common_type_t<float128, int>, float128>);
    static_assert(std::is_same_v<std::common_type_t<float128, float128>, float128>);
}

/**********************************************************************
 * Standard library surface for the integer and fixed point types
 *
 * The same treatment float128 has: numeric_limits so that generic code compiles against these
 * types, a hash so they can be container keys, a formatter and the stream operators.
 ***********************************************************************/
TEST(int128_t, NumericLimits)
{
    using signed_limits = std::numeric_limits<int128_t>;
    using unsigned_limits = std::numeric_limits<uint128_t>;

    static_assert(signed_limits::is_specialized && signed_limits::is_integer && signed_limits::is_exact);
    static_assert(signed_limits::is_signed && !unsigned_limits::is_signed);
    static_assert(signed_limits::digits == 127 && unsigned_limits::digits == 128);
    static_assert(signed_limits::radix == 2 && !signed_limits::has_infinity);
    static_assert(unsigned_limits::is_modulo && !signed_limits::is_modulo);
    static_assert(std::numeric_limits<const int128_t>::digits == signed_limits::digits);

    // The extremes wrap into each other, which is what makes them the extremes.
    EXPECT_TRUE(signed_limits::max() + int128_t(1) == signed_limits::min());
    EXPECT_TRUE(signed_limits::min() - int128_t(1) == signed_limits::max());
    EXPECT_TRUE(signed_limits::lowest() == signed_limits::min());
    EXPECT_TRUE(unsigned_limits::min().is_zero());
    EXPECT_TRUE(unsigned_limits::max() + uint128_t(1) == uint128_t(0));

    EXPECT_EQ(static_cast<std::string>(signed_limits::max()), "170141183460469231731687303715884105727");
    EXPECT_EQ(static_cast<std::string>(signed_limits::min()), "-170141183460469231731687303715884105728");
    EXPECT_EQ(static_cast<std::string>(unsigned_limits::max()), "340282366920938463463374607431768211455");
}
TEST(int128_t, FormatAndStream)
{
    const int128_t value("-12345678901234567890123456789");

    EXPECT_EQ(std::format("{}", value), "-12345678901234567890123456789");
    EXPECT_EQ(std::format("{}", uint128_t(42)), "42");
    EXPECT_EQ(std::format("{:>15}", int128_t(42)), "             42");
    EXPECT_EQ(std::format("{:<15}", int128_t(42)), "42             ");
    EXPECT_EQ(std::format("{:*^10}", int128_t(42)), "****42****");
    EXPECT_EQ(std::format("{:08}", int128_t(-42)), "-0000042");
    EXPECT_EQ(std::format("{:+}", int128_t(42)), "+42");
    EXPECT_EQ(std::format("{:d}", int128_t(42)), "42");

    std::ostringstream stream;
    stream << value;
    EXPECT_EQ(stream.str(), "-12345678901234567890123456789");

    stream.str("");
    stream << std::setw(8) << std::setfill('.') << uint128_t(42);
    EXPECT_EQ(stream.str(), "......42");

    int128_t read_back;
    std::istringstream input("-98765432109876543210");
    input >> read_back;
    EXPECT_TRUE(read_back == int128_t("-98765432109876543210"));
}
TEST(int128_t, Hash)
{
    std::unordered_set<int128_t> keys;
    keys.insert(int128_t(1));
    keys.insert(int128_t(2));
    keys.insert(int128_t(1));
    EXPECT_EQ(keys.size(), 2u);

    std::unordered_set<uint128_t> unsigned_keys;
    unsigned_keys.insert(std::numeric_limits<uint128_t>::max());
    unsigned_keys.insert(uint128_t(0));
    EXPECT_EQ(unsigned_keys.size(), 2u);
}
TEST(fixed_point128, NumericLimits)
{
    using limits = std::numeric_limits<fixed_point128<32>>;

    static_assert(limits::is_specialized && limits::is_signed);
    // A fixed point value is exactly the number it stands for, unlike a floating point one.
    static_assert(limits::is_exact && !limits::is_integer);
    static_assert(limits::digits == 128 && limits::radix == 2);
    static_assert(!limits::has_infinity && !limits::has_quiet_NaN);
    static_assert(limits::min_exponent == 32 - 128 && limits::max_exponent == 32);

    // The grid spacing is the same everywhere, which is the whole point of the representation.
    const fixed_point128<32> one = fixed_point128<32>::one();
    EXPECT_TRUE(one + limits::epsilon() > one);
    EXPECT_TRUE(limits::min() == limits::epsilon());
    // The sign is a separate field, so the range is symmetric.
    EXPECT_TRUE(limits::lowest() == -limits::max());
}
TEST(fixed_point128, FormatAndStream)
{
    const fixed_point128<32> value("-1234.5");

    EXPECT_EQ(std::format("{}", value), static_cast<std::string>(value));
    EXPECT_EQ(std::format("{:>15}", fixed_point128<32>(42)), "             42");
    EXPECT_EQ(std::format("{:*^10}", fixed_point128<32>(42)), "****42****");
    EXPECT_EQ(std::format("{:+}", fixed_point128<32>(42)), "+42");

    std::ostringstream stream;
    stream << value;
    EXPECT_EQ(stream.str(), static_cast<std::string>(value));

    fixed_point128<32> read_back;
    std::istringstream input("2.5");
    input >> read_back;
    EXPECT_TRUE(read_back == fixed_point128<32>("2.5"));
}
TEST(fixed_point128, Hash)
{
    std::unordered_set<fixed_point128<32>> keys;
    keys.insert(fixed_point128<32>(1));
    keys.insert(fixed_point128<32>(2));
    keys.insert(fixed_point128<32>(1));
    EXPECT_EQ(keys.size(), 2u);

    // the two zeros compare equal and so have to be the same key
    keys.insert(fixed_point128<32>());
    keys.insert(-fixed_point128<32>());
    EXPECT_EQ(keys.size(), 3u);
}
