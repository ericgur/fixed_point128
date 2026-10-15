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
#include <charconv>
#include <format>
#include <iomanip>
#include <sstream>
#include <unordered_map>
#include <unordered_set>
#include <utility>  // std::integer_sequence
#include "gtest_shared.h"

namespace
{
/// @brief The binary128 encoding of a value, high QWORD first, the way the bits read.
[[nodiscard]] std::pair<uint64_t, uint64_t> EncodingOf(const float128& x)
{
    uint64_t low = 0, high = 0;
    x.get_bits(low, high);
    return {high, low};
}

/// @brief The high QWORD of a binary128 encoding, built from its fields as IEEE 754 lays them out.
/// @param sign Sign bit
/// @param biased Biased exponent, 15 bits
/// @param fraction Top 48 bits of the fraction
[[nodiscard]] constexpr uint64_t HighQword(uint64_t sign, uint64_t biased, uint64_t fraction)
{
    return (sign << 63) | (biased << 48) | fraction;
}

/// @brief Formats a float128 to text with fp128::to_chars in the given format, as std::to_chars would.
[[nodiscard]] std::string ToChars(const float128& x, std::chars_format fmt, int precision = -1)
{
    char buffer[128];
    const std::to_chars_result result = (precision < 0) ? fp128::to_chars(buffer, buffer + sizeof(buffer), x, fmt)
                                                        : fp128::to_chars(buffer, buffer + sizeof(buffer), x, fmt, precision);
    EXPECT_EQ(result.ec, std::errc {});
    return std::string(buffer, result.ptr);
}

/// @brief Number of significant digits in a decimal string: from its first non zero digit to its
///        last one before any exponent, the point not counting. The zeros that fill out the integer
///        part of a large value in fixed notation are not significant.
[[nodiscard]] int32_t SignificantDigits(const std::string& text)
{
    int32_t count = 0, last_nonzero = 0;
    bool started = false;
    for (char c : text) {
        if (c == 'e' || c == 'E')
            break;
        if (c < '0' || c > '9')
            continue;
        started = started || (c != '0');
        if (started)
            ++count;
        if (c != '0')
            last_nonzero = count;
    }
    return last_nonzero;
}
}  // namespace

/**********************************************************************
 * float128 text conversion
 *
 * The conversions in both directions are exact: the decimal expansion of a binary128 is finite and
 * is produced in full, so the digits are correctly rounded at any precision and a string reads
 * back to the nearest representable value. That is what makes a round trip meaningful, and most of
 * what is checked here. Expected values come from Python's fractions and decimal modules.
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
// Random encodings over the whole range, either sign, and a run of subnormals of their own: drawn
// with the rest, a subnormal turned up once in 32767 draws, a quarter of a time per run.
TEST(float128_format, RoundTripIsExact)
{
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT / 8; ++i) {
        const uint64_t low = get_uint64_random();
        const uint64_t fraction = get_uint64_random() & 0x0000FFFFFFFFFFFFull;
        const uint64_t biased = (i % 8 == 0) ? 0 : get_uint32_random() % 0x7FFF;  // 0 is a subnormal
        const uint64_t sign = get_uint32_random() & 1;
        const float128 value(low, HighQword(sign, biased, fraction));

        const std::string text = to_string(value);
        float128 parsed;
        const std::from_chars_result result = fp128::from_chars(text.data(), text.data() + text.size(), parsed);
        EXPECT_EQ(result.ptr, text.data() + text.size()) << "did not consume " << text;
        EXPECT_TRUE(EncodingOf(parsed) == EncodingOf(value)) << "round trip of " << text;
    }
}
// to_string() promises the fewest digits that read back. Rounded to one digit fewer, the value
// must therefore read back as something else - which a conversion that always wrote all 36 digits,
// and still round tripped, would fail.
TEST(float128_format, ShortestIsShortest)
{
    // A value far from one takes a few milliseconds to convert, its exact expansion running to
    // thousands of digits, so this samples a thousand of them.
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT / 64; ++i) {
        uint64_t low = get_uint64_random();
        uint64_t fraction = get_uint64_random() & 0x0000FFFFFFFFFFFFull;
        // keep only the leading bits of the fraction now and then, so that short forms occur
        const uint32_t keep = get_uint32_random() % 113;
        if (keep <= 48) {
            low = 0;
            fraction &= ~((1ull << (48 - keep)) - 1);
        } else if (keep < 112) {
            low &= ~((1ull << (112 - keep)) - 1);
        }
        const uint64_t biased = 1 + get_uint32_random() % 0x7FFD;
        const float128 value(low, HighQword(0, biased, fraction));

        const std::string text = to_string(value);
        const int32_t digits = SignificantDigits(text);
        if (digits <= 1)
            continue;  // nothing shorter exists
        const std::string shorter = ToChars(value, std::chars_format::scientific, digits - 2);
        float128 back;
        (void)fp128::from_chars(shorter.data(), shorter.data() + shorter.size(), back);
        EXPECT_FALSE(EncodingOf(back) == EncodingOf(value)) << text << " has a shorter form " << shorter;
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
/// @brief A decimal text and the binary128 it rounds to, high QWORD first.
struct ParseCase {
    const char* text;
    uint64_t high, low;
};
// Points exactly halfway between two binary128 values, written out in full, and the values a hair
// above and below them. A tie goes to the even neighbour, anything off it to the nearer one. All
// but the integer ones need more than 40 significant digits; past the fortieth the digits used to
// be dropped, which turned a value just above a tie into the tie: 2^113 + 1 + 10^-6, the last
// case, read as 2^113.
TEST(float128_format, ParseRoundsTiesBeyondFortyDigits)
{
    static constexpr ParseCase cases[] = {
        {"1.00000000000000000000000000000000067407548053553254856959227990472456148833557687538586833397857844829559326171875", 0x3FFF000000000000ull, 0x0000000000000004ull},
        {"1.0000000000000000000000000000000006740754805355325485695922799047245614883355768753858683339785784482955932617187500001", 0x3FFF000000000000ull, 0x0000000000000004ull},
        {"1.0000000000000000000000000000000006740754805355325485695922799047245614883355768753858683339785784482955932617187499999", 0x3FFF000000000000ull, 0x0000000000000003ull},
        {"10384593717069655257060992658440193", 0x4070000000000000ull, 0x0000000000000000ull},
        {"10384593717069655257060992658440193.00001", 0x4070000000000000ull, 0x0000000000000001ull},
        {"10384593717069655257060992658440192.99999", 0x4070000000000000ull, 0x0000000000000000ull},
        {"0.000000000000000000000000000001479114197289397135146991059905224256601804031907808113730211151490774421267726668959559730238239926632065845158489272264951454171308800744221063905303736873975140042603015899658203125", 0x3F9BE00000000000ull, 0x0000000000000000ull},
        {"0.00000000000000000000000000000147911419728939713514699105990522425660180403190780811373021115149077442126772666895955973023823992663206584515848927226495145417130880074422106390530373687397514004260301589965820312500001", 0x3F9BE00000000000ull, 0x0000000000000001ull},
        {"0.00000000000000000000000000000147911419728939713514699105990522425660180403190780811373021115149077442126772666895955973023823992663206584515848927226495145417130880074422106390530373687397514004260301589965820312499999", 0x3F9BE00000000000ull, 0x0000000000000000ull},
        {"0.3333333333333333333333333333333333413580414349468160543999080941038638272420902008974508134997449815273284912109375", 0x3FFD555555555555ull, 0x5555555555555556ull},
        {"0.333333333333333333333333333333333341358041434946816054399908094103863827242090200897450813499744981527328491210937500001", 0x3FFD555555555555ull, 0x5555555555555556ull},
        {"0.333333333333333333333333333333333341358041434946816054399908094103863827242090200897450813499744981527328491210937499999", 0x3FFD555555555555ull, 0x5555555555555555ull},
        {"10384593717069655257060992658440193.000001", 0x4070000000000000ull, 0x0000000000000001ull},
    };
    for (const ParseCase& c : cases) {
        const float128 parsed(c.text);
        EXPECT_EQ(EncodingOf(parsed), std::make_pair(c.high, c.low)) << c.text;
        EXPECT_EQ(EncodingOf(-float128(("-" + std::string(c.text)).c_str())), std::make_pair(c.high, c.low)) << "-" << c.text;
    }
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
    EXPECT_EQ(std::format("{:.3g}", float128("0.0001234")), "0.000123");
    EXPECT_EQ(std::format("{:.3g}", float128("0.00001234")), "1.23e-05");
    EXPECT_EQ(std::format("{:#g}", float128("1.5")), "1.50000");
    EXPECT_EQ(std::format("{:G}", float128("1e-10")), "1E-10");
    EXPECT_EQ(std::format("{:G}", float128::inf()), "INF");
    EXPECT_EQ(std::format("{:F}", float128::nan()), "NAN");
    EXPECT_EQ(std::format("{:+}", float128::nan()), "+nan");

    // width, fill and alignment
    EXPECT_EQ(std::format("{:>20.4f}", pi), "              3.1416");
    EXPECT_EQ(std::format("{:<20.4f}", pi), "3.1416              ");
    EXPECT_EQ(std::format("{:^20.4f}", pi), "       3.1416       ");
    EXPECT_EQ(std::format("{:*^20.4f}", pi), "*******3.1416*******");
    // zero padding goes after the sign
    EXPECT_EQ(std::format("{:020.4f}", -pi), "-00000000000003.1416");
    EXPECT_EQ(std::format("{:+.4f}", pi), "+3.1416");
    EXPECT_EQ(std::format("{: .4f}", pi), " 3.1416");
    EXPECT_EQ(std::format("{: .4f}", -pi), "-3.1416");
    // the alternate form keeps the point
    EXPECT_EQ(std::format("{:#.0f}", float128(3)), "3.");
    // a special value is right aligned but never zero padded
    EXPECT_EQ(std::format("{:10}", float128::inf()), "       inf");
    EXPECT_EQ(std::format("{:010}", float128::nan()), "       nan");
}
/// @brief A value, a precision and the fixed notation text it rounds to.
struct FixedCase {
    const char* text;
    int precision;
    const char* expected;
};
// The value lies below the last digit asked for, so every digit shown is zero but possibly the
// last. 0.7 with no decimals printed as 0 - the integer one that rounding up adds was dropped - and
// 0.046 with one as 0.1: the half place was judged on its first digit rounded, 4.6 to 5.
TEST(float128_format, FixedNotationBelowTheLastDigit)
{
    static constexpr FixedCase cases[] = {
        {"0.7", 0, "1"},           {"0.9", 0, "1"},           {"0.5", 0, "0"},           {"1.5", 0, "2"},
        {"0.5000001", 0, "1"},     {"0.49999999", 0, "0"},    {"0.046", 1, "0.0"},       {"0.05", 1, "0.1"},
        {"0.0451", 1, "0.0"},      {"0.0449", 1, "0.0"},      {"0.0046", 2, "0.00"},     {"0.0049", 2, "0.00"},
        {"0.006", 2, "0.01"},      {"0.0051", 2, "0.01"},     {"-0.0046", 2, "-0.00"},   {"-0.7", 0, "-1"},
        {"0.00095", 3, "0.001"},   {"0.000951", 3, "0.001"},  {"9.5e-7", 6, "0.000001"}, {"4.5e-7", 6, "0.000000"},
    };
    for (const FixedCase& c : cases) {
        const float128 value(c.text);
        const std::string specification = "{:." + std::to_string(c.precision) + "f}";
        EXPECT_EQ(std::vformat(specification, std::make_format_args(value)), c.expected) << c.text;
        EXPECT_EQ(ToChars(value, std::chars_format::fixed, c.precision), c.expected) << c.text;
        std::ostringstream stream;
        stream << std::fixed << std::setprecision(c.precision) << value;
        EXPECT_EQ(stream.str(), c.expected) << c.text;
    }
    EXPECT_EQ(std::format("{:#.0f}", float128("0.7")), "1.");
    EXPECT_EQ(std::format("{:#.0f}", float128("0.3")), "0.");
}
TEST(float128_format, ZeroAndDefaultFormats)
{
    EXPECT_EQ(std::format("{:e}", float128()), "0.000000e+00");
    EXPECT_EQ(std::format("{:.2e}", -float128()), "-0.00e+00");
    EXPECT_EQ(std::format("{:g}", float128()), "0");
    EXPECT_EQ(std::format("{:#g}", float128()), "0.00000");
    EXPECT_EQ(std::format("{:f}", float128()), "0.000000");
    EXPECT_EQ(std::format("{:.3a}", float128()), "0x0.000p+00");
    // An empty specification is the shortest form, as for a double.
    EXPECT_EQ(std::format("{}", float128("0.1")), "0.1");
    EXPECT_EQ(std::format("{}", std::numeric_limits<float128>::max()), "1.189731495357231765085759326628007e+4932");
    EXPECT_EQ(std::format("{}", float128::pi()), to_string(float128::pi()));
    EXPECT_EQ(std::format("{}", -float128::inf()), "-inf");
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
    // a carry runs through a fraction of f digits into the leading one
    EXPECT_EQ(std::format("{:.2a}", float128("0x1.fffp0")), "0x2.00p+00");
    // past the 28 digits of the fraction the precision is filled with zeros; it used to stop at 28
    EXPECT_EQ(std::format("{:.30a}", float128::pi()), "0x1.921fb54442d18469898cc51701b800p+01");
    EXPECT_EQ(std::format("{:a}", std::numeric_limits<float128>::denorm_min()), "0x1p-16494");
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

    stream.str("");
    stream << std::setw(6) << std::left << std::noshowpos << std::setfill('*') << float128(-1);
    EXPECT_EQ(stream.str(), "-1****");

    stream.str("");
    stream << std::right << std::setfill(' ') << std::showpoint << std::setprecision(4) << float128(2);
    EXPECT_EQ(stream.str(), "2.000");
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
// The extractor reads what the string constructor reads, and like a builtin one only as far as the
// number goes: "1.5,2.5" used to fail as one token, and a hexadecimal number to fail outright.
TEST(float128_format, StreamExtractionStopsWhereTheNumberEnds)
{
    float128 value;
    char separator = 0;
    std::istringstream csv("1.5,2.5");
    csv >> value >> separator;
    EXPECT_TRUE(value == float128("1.5"));
    EXPECT_EQ(separator, ',');
    csv >> value;
    EXPECT_TRUE(value == float128("2.5"));
    EXPECT_TRUE(csv.eof() && !csv.fail());

    std::istringstream mixed("0x1.8p1 -0X1P-2 inf -Infinity nan(7) 1.5e3x 7");
    mixed >> value;
    EXPECT_TRUE(value == float128(3));
    mixed >> value;
    EXPECT_TRUE(value == float128("-0.25"));
    mixed >> value;
    EXPECT_TRUE(isinf(value) && value.is_positive());
    mixed >> value;
    EXPECT_TRUE(isinf(value) && value.is_negative());
    mixed >> value;
    EXPECT_TRUE(isnan(value));
    mixed >> value >> separator;
    EXPECT_TRUE(value == float128(1500));
    EXPECT_EQ(separator, 'x');
    mixed >> value;
    EXPECT_TRUE(value == float128(7));

    // a number that stops halfway fails and leaves the value alone
    for (const char* text : {"1e", "0x", "abc", "-", ".", "1e+"}) {
        float128 untouched(42);
        std::istringstream bad(text);
        bad >> untouched;
        EXPECT_TRUE(bad.fail()) << text;
        EXPECT_TRUE(untouched == float128(42)) << text;
    }

    // with skipws off a leading space is not skipped
    float128 untouched(42);
    std::istringstream spaced(" 1");
    spaced >> std::noskipws >> untouched;
    EXPECT_TRUE(spaced.fail());
    EXPECT_TRUE(untouched == float128(42));
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

    written = fp128::to_chars(buffer, buffer + sizeof(buffer), float128::pi(), std::chars_format::general, 3);
    EXPECT_EQ(std::string(buffer, written.ptr), "3.14");

    // too small a range reports the failure rather than truncating
    written = fp128::to_chars(buffer, buffer + 2, float128::pi());
    EXPECT_EQ(written.ec, std::errc::value_too_large);
    written = fp128::to_chars(buffer, buffer + 2, float128::pi(), std::chars_format::fixed, 4);
    EXPECT_EQ(written.ec, std::errc::value_too_large);
    EXPECT_EQ(written.ptr, buffer + 2);

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
TEST(float128_format, FromCharsEdges)
{
    float128 value(7);
    const char* empty = "";
    auto read = fp128::from_chars(empty, empty, value);
    EXPECT_EQ(read.ec, std::errc::invalid_argument);
    EXPECT_TRUE(value == float128(7));  // untouched

    // the zeros, with their sign
    for (const char* zero : {"0", "-0.0", "0e5000", "-000.000e-7"}) {
        read = fp128::from_chars(zero, zero + strlen(zero), value);
        EXPECT_EQ(read.ec, std::errc {}) << zero;
        EXPECT_EQ(read.ptr, zero + strlen(zero)) << zero;
        EXPECT_TRUE(value.is_zero() && value.is_negative() == (zero[0] == '-')) << zero;
    }

    // out of range: the infinity or zero it rounds to, and the error that says so
    const char* huge = "-1e5000";
    read = fp128::from_chars(huge, huge + 7, value);
    EXPECT_EQ(read.ec, std::errc::result_out_of_range);
    EXPECT_TRUE(isinf(value) && value.is_negative());
    const char* tiny = "1e-5000";
    read = fp128::from_chars(tiny, tiny + 7, value);
    EXPECT_EQ(read.ec, std::errc::result_out_of_range);
    EXPECT_TRUE(value.is_zero());
    // the smallest subnormal is 6.48e-4966, so 4e-4966 rounds up to it and 3e-4966 down to zero
    const char* sub = "4e-4966";
    read = fp128::from_chars(sub, sub + 7, value);
    EXPECT_EQ(read.ec, std::errc {});
    EXPECT_TRUE(value == std::numeric_limits<float128>::denorm_min());
    const char* gone = "3e-4966";
    read = fp128::from_chars(gone, gone + 7, value);
    EXPECT_EQ(read.ec, std::errc::result_out_of_range);
    EXPECT_TRUE(value.is_zero());

    // hexadecimal: upper case digits, leading zeros after the point, a sign, a zero, no digits
    const auto hex = [&](const char* text, std::errc ec) {
        const std::from_chars_result result = fp128::from_chars(text, text + strlen(text), value, std::chars_format::hex);
        EXPECT_EQ(result.ec, ec) << text;
        return result.ptr - text;
    };
    EXPECT_EQ(hex("A.Bp1", std::errc {}), 5);
    EXPECT_TRUE(value == float128(21.375));
    EXPECT_EQ(hex("0.0010p4", std::errc {}), 8);
    EXPECT_TRUE(value == ldexp(float128::one(), -8));
    EXPECT_EQ(hex("-1.8", std::errc {}), 4);
    EXPECT_TRUE(value == float128(-1.5));
    EXPECT_EQ(hex("0p5", std::errc {}), 3);
    EXPECT_TRUE(value.is_zero() && value.is_positive());
    EXPECT_EQ(hex("-inf", std::errc {}), 4);
    EXPECT_TRUE(isinf(value) && value.is_negative());
    value = float128(7);
    EXPECT_EQ(hex("p5", std::errc::invalid_argument), 0);
    EXPECT_TRUE(value == float128(7));
    // more digits than the mantissa holds: everything past them is a sticky bit. 1 + 2^-113 is the
    // tie between 1 and its neighbour and goes to the even one; a trailing 1 far below it does not.
    EXPECT_EQ(hex("1.00000000000000000000000000008", std::errc {}), 31);
    EXPECT_TRUE(value == float128::one());
    EXPECT_EQ(hex("1.000000000000000000000000000080000001", std::errc {}), 38);
    EXPECT_TRUE(value == float128::nextUp(float128::one()));

    // every other format reads the decimal form
    const char* decimal = "2.5e1";
    read = fp128::from_chars(decimal, decimal + 5, value, std::chars_format::general);
    EXPECT_EQ(read.ptr, decimal + 5);
    EXPECT_TRUE(value == float128(25));
}
// What std::to_chars writes for a double: no 0x prefix, as few exponent digits as needed. That is
// what from_chars with chars_format::hex reads, so the two make a round trip; with the prefix and
// a padded exponent, as {:a} writes, it read the leading 0 and stopped.
TEST(float128_format, HexadecimalToCharsRoundTrips)
{
    EXPECT_EQ(ToChars(float128(3), std::chars_format::hex), "1.8p+1");
    EXPECT_EQ(ToChars(float128::one(), std::chars_format::hex), "1p+0");
    EXPECT_EQ(ToChars(float128("-0.5"), std::chars_format::hex), "-1p-1");
    EXPECT_EQ(ToChars(float128(), std::chars_format::hex), "0p+0");
    EXPECT_EQ(ToChars(-float128(), std::chars_format::hex), "-0p+0");
    EXPECT_EQ(ToChars(std::numeric_limits<float128>::denorm_min(), std::chars_format::hex), "1p-16494");
    EXPECT_EQ(ToChars(std::numeric_limits<float128>::max(), std::chars_format::hex), "1.ffffffffffffffffffffffffffffp+16383");
    EXPECT_EQ(ToChars(float128::pi(), std::chars_format::hex, 3), "1.922p+1");
    EXPECT_EQ(ToChars(float128(1.5), std::chars_format::hex, 0), "2p+0");
    EXPECT_EQ(ToChars(float128::inf(), std::chars_format::hex), "inf");

    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT / 16; ++i) {
        const uint64_t low = get_uint64_random();
        const uint64_t fraction = get_uint64_random() & 0x0000FFFFFFFFFFFFull;
        const uint64_t biased = get_uint32_random() % 0x7FFF;
        const uint64_t sign = get_uint32_random() & 1;
        const float128 value(low, HighQword(sign, biased, fraction));
        const std::string text = ToChars(value, std::chars_format::hex);
        float128 back;
        const std::from_chars_result result = fp128::from_chars(text.data(), text.data() + text.size(), back, std::chars_format::hex);
        EXPECT_EQ(result.ptr, text.data() + text.size()) << text;
        EXPECT_TRUE(EncodingOf(back) == EncodingOf(value)) << text;
    }
}
// The hexadecimal path of the string constructor applied the sign only when the value was in
// range: -0x1p100000 read as +inf.
TEST(float128_format, HexadecimalConstructorKeepsTheSignOutOfRange)
{
    EXPECT_EQ(EncodingOf(float128("-0x1p100000")), EncodingOf(-float128::inf()));
    EXPECT_EQ(EncodingOf(float128("-0x1p-20000")), EncodingOf(-float128()));
    EXPECT_EQ(EncodingOf(float128("0x1p100000")), EncodingOf(float128::inf()));
    EXPECT_EQ(EncodingOf(float128("-0x1p-16494")), EncodingOf(-std::numeric_limits<float128>::denorm_min()));
}
TEST(float128_format, InvalidSpecificationsThrow)
{
    const float128 one = float128::one();
    const int128_t signed_one(1);
    const fixed_point128<32> fixed_one(1);
    EXPECT_THROW((void)std::vformat("{:q}", std::make_format_args(one)), std::format_error);
    EXPECT_THROW((void)std::vformat("{:.3fx}", std::make_format_args(one)), std::format_error);
    EXPECT_THROW((void)std::vformat("{:x}", std::make_format_args(signed_one)), std::format_error);
    EXPECT_THROW((void)std::vformat("{:.3}", std::make_format_args(fixed_one)), std::format_error);
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
    static_assert(std::is_same_v<std::common_type_t<float128, long double>, float128>);
    static_assert(std::is_same_v<std::common_type_t<uint64_t, float128>, float128>);
    static_assert(std::is_same_v<std::common_type_t<float128, float128>, float128>);
}
// Every member against the binary128 format as IEEE 754 defines it: 113 bit significand, 15 bit
// exponent biased by 16383. The decimal members follow from those: digits10 = floor(112 log10 2),
// max_digits10 = ceil(1 + 113 log10 2), min_exponent10 = ceil(log10(2^-16382)) and
// max_exponent10 = floor(log10(max())).
TEST(float128_format, NumericLimitsEveryMember)
{
    using limits = std::numeric_limits<float128>;
    static_assert(limits::is_specialized && limits::is_signed && !limits::is_integer && !limits::is_exact);
    static_assert(limits::has_infinity && limits::has_quiet_NaN && limits::has_signaling_NaN);
    FP128_SUPPRESS_DEPRECATED_BEGIN
    static_assert(limits::has_denorm == std::denorm_present && !limits::has_denorm_loss);
    FP128_SUPPRESS_DEPRECATED_END
    static_assert(limits::round_style == std::round_to_nearest && limits::is_bounded && !limits::is_modulo);
    static_assert(!limits::traps && !limits::tinyness_before);
    static_assert(limits::digits == 113 && limits::digits10 == 33 && limits::max_digits10 == 36 && limits::radix == 2);
    static_assert(limits::min_exponent == -16381 && limits::max_exponent == 16384);
    static_assert(limits::min_exponent10 == -4931 && limits::max_exponent10 == 4932);

    // the encodings, from the fields
    using encoding = std::pair<uint64_t, uint64_t>;
    EXPECT_EQ(EncodingOf(limits::max()), encoding(HighQword(0, 0x7FFE, 0xFFFFFFFFFFFFull), UINT64_MAX));
    EXPECT_EQ(EncodingOf(limits::lowest()), encoding(HighQword(1, 0x7FFE, 0xFFFFFFFFFFFFull), UINT64_MAX));
    EXPECT_EQ(EncodingOf(limits::min()), encoding(HighQword(0, 1, 0), 0));
    EXPECT_EQ(EncodingOf(limits::epsilon()), encoding(HighQword(0, 16383 - 112, 0), 0));
    EXPECT_EQ(EncodingOf(limits::round_error()), encoding(HighQword(0, 16383 - 1, 0), 0));
    EXPECT_EQ(EncodingOf(limits::denorm_min()), encoding(0, 1));
    EXPECT_EQ(EncodingOf(limits::infinity()), encoding(HighQword(0, 0x7FFF, 0), 0));
    // a quiet NaN has the top fraction bit set, a signaling one has it clear and some other one set
    const auto [quiet_high, quiet_low] = EncodingOf(limits::quiet_NaN());
    EXPECT_EQ(quiet_high & HighQword(0, 0x7FFF, 0x800000000000ull), HighQword(0, 0x7FFF, 0x800000000000ull));
    const auto [signaling_high, signaling_low] = EncodingOf(limits::signaling_NaN());
    EXPECT_EQ(signaling_high & HighQword(0, 0x7FFF, 0x800000000000ull), HighQword(0, 0x7FFF, 0));
    EXPECT_TRUE((signaling_high & 0x7FFFFFFFFFFFull) != 0 || signaling_low != 0);
    EXPECT_TRUE(isnan(limits::signaling_NaN()));

    // and what the members mean: the decimal ones at the edge of their definitions
    const float128 ten(10);
    EXPECT_TRUE(pow(ten, float128(limits::max_exponent10)) <= limits::max());
    EXPECT_TRUE(isinf(pow(ten, float128(limits::max_exponent10 + 1))));
    EXPECT_TRUE(float128("1e-4931") >= limits::min());
    EXPECT_TRUE(float128("1e-4932") < limits::min());
    static_assert(std::numeric_limits<const volatile float128>::digits == 113);
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
    static_assert(signed_limits::digits10 == 38 && unsigned_limits::digits10 == 38);
    static_assert(std::numeric_limits<const int128_t>::digits == signed_limits::digits);
    static_assert(std::numeric_limits<const volatile uint128_t>::digits == unsigned_limits::digits);

    // Both types are modulo: the signed one wraps too, its overflow being defined rather than
    // undefined like a builtin's. The extremes wrap into each other under +, - and *.
    static_assert(unsigned_limits::is_modulo && signed_limits::is_modulo);
    EXPECT_TRUE(signed_limits::max() + int128_t(1) == signed_limits::min());
    EXPECT_TRUE(signed_limits::min() - int128_t(1) == signed_limits::max());
    EXPECT_TRUE(signed_limits::max() * int128_t(2) == int128_t(-2));
    EXPECT_TRUE(signed_limits::min() * int128_t(-1) == signed_limits::min());

    // A zero divisor throws, which is the types' trap.
    static_assert(signed_limits::traps && unsigned_limits::traps);
    EXPECT_THROW((void)(int128_t(1) / int128_t(0)), std::logic_error);
    EXPECT_THROW((void)(uint128_t(1) / uint128_t(0)), std::logic_error);
    EXPECT_TRUE(signed_limits::lowest() == signed_limits::min());
    EXPECT_TRUE(unsigned_limits::min().is_zero());
    EXPECT_TRUE(unsigned_limits::max() + uint128_t(1) == uint128_t(0));

    EXPECT_EQ(static_cast<std::string>(signed_limits::max()), "170141183460469231731687303715884105727");
    EXPECT_EQ(static_cast<std::string>(signed_limits::min()), "-170141183460469231731687303715884105728");
    EXPECT_EQ(static_cast<std::string>(unsigned_limits::max()), "340282366920938463463374607431768211455");
}
/// @brief The members numeric_limits gives an integer type, for either 128 bit one.
template <typename T> static void CheckIntegerLimitsEveryMember(bool is_signed)
{
    using limits = std::numeric_limits<T>;
    static_assert(limits::is_specialized && limits::is_integer && limits::is_exact && limits::is_bounded && limits::is_modulo);
    static_assert(!limits::has_infinity && !limits::has_quiet_NaN && !limits::has_signaling_NaN && !limits::is_iec559);
    FP128_SUPPRESS_DEPRECATED_BEGIN
    static_assert(limits::has_denorm == std::denorm_absent && !limits::has_denorm_loss);
    FP128_SUPPRESS_DEPRECATED_END
    // integer division truncates, as the builtin integers' does
    static_assert(limits::round_style == std::round_toward_zero && limits::traps && !limits::tinyness_before);
    static_assert(limits::radix == 2 && limits::max_digits10 == 0);
    static_assert(limits::min_exponent == 0 && limits::max_exponent == 0 && limits::min_exponent10 == 0 && limits::max_exponent10 == 0);
    EXPECT_EQ(limits::is_signed, is_signed);
    // digits are the value bits; digits10 is floor(digits * log10(2)), 38 for 127 and for 128
    EXPECT_EQ(limits::digits, is_signed ? 127 : 128);
    EXPECT_EQ(limits::digits10, 38);
    EXPECT_TRUE(T(7) / T(2) == T(3));
    if (is_signed) {
        EXPECT_TRUE(T(-7) / T(2) == T(-3));
    }

    uint64_t low = 0, high = 0;
    limits::max().get_components(low, high);
    EXPECT_EQ(low, UINT64_MAX);
    EXPECT_EQ(high, is_signed ? (UINT64_MAX >> 1) : UINT64_MAX);
    limits::min().get_components(low, high);
    EXPECT_EQ(low, 0u);
    EXPECT_EQ(high, is_signed ? (1ull << 63) : 0u);
    EXPECT_TRUE(limits::lowest() == limits::min());
    for (const T& zero : {limits::epsilon(), limits::round_error(), limits::infinity(), limits::quiet_NaN(), limits::signaling_NaN(), limits::denorm_min()}) {
        EXPECT_TRUE(zero.is_zero());
    }
}
TEST(int128_t, NumericLimitsEveryMember)
{
    CheckIntegerLimitsEveryMember<int128_t>(true);
    CheckIntegerLimitsEveryMember<uint128_t>(false);
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
// Values past 2^64, every sign option, zero padding with a sign and with the most negative value,
// and the stream flags, with literal expectations.
TEST(int128_t, StreamAndFormatFlags)
{
    const uint128_t big = std::numeric_limits<uint128_t>::max();
    const int128_t min = std::numeric_limits<int128_t>::min();
    EXPECT_EQ(std::format("{}", big), "340282366920938463463374607431768211455");
    EXPECT_EQ(std::format("{:+}", uint128_t(1, 1)), "+18446744073709551617");
    EXPECT_EQ(std::format("{:045}", min), "-00000170141183460469231731687303715884105728");
    EXPECT_EQ(std::format("{:>45}", min), "     -170141183460469231731687303715884105728");
    EXPECT_EQ(std::format("{:+06}", int128_t(42)), "+00042");
    EXPECT_EQ(std::format("{:-}", int128_t(42)), "42");
    EXPECT_EQ(std::format("{: }", int128_t(42)), " 42");
    EXPECT_EQ(std::format("{: }", int128_t(-42)), "-42");
    EXPECT_EQ(std::format("{:<06}", int128_t(42)), "42    ");  // an alignment overrides the zero padding
    EXPECT_EQ(std::format("{:^7}", uint128_t(42)), "  42   ");
    EXPECT_EQ(std::format("{:023}", uint128_t(0, 1)), "00018446744073709551616");

    std::ostringstream stream;
    stream << std::showpos << int128_t(42);
    EXPECT_EQ(stream.str(), "+42");
    stream.str("");
    stream << std::noshowpos << std::left << std::setw(6) << std::setfill('*') << int128_t(-42);
    EXPECT_EQ(stream.str(), "-42***");
    stream.str("");
    stream << std::internal << std::setw(6) << std::setfill('0') << int128_t(-42);
    EXPECT_EQ(stream.str(), "-00042");
    stream.str("");
    stream << std::right << std::setfill(' ') << std::setw(41) << big;
    EXPECT_EQ(stream.str(), "  340282366920938463463374607431768211455");
}
// The extractors read like a builtin integer's: the digits only, failure and overflow set failbit
// and store zero or the end of the range. They used to read a whole token through the string
// constructor, which never failed: "abc" read as zero and an overflow wrapped.
TEST(int128_t, StreamExtractionFailsLikeABuiltinInteger)
{
    int128_t value(7);
    char separator = 0;
    std::istringstream csv("12,-34");
    csv >> value >> separator;
    EXPECT_TRUE(value == int128_t(12));
    EXPECT_EQ(separator, ',');
    csv >> value;
    EXPECT_TRUE(value == int128_t(-34) && csv.eof() && !csv.fail());

    // 2^127 - 1 is the largest, -2^127 the smallest; one past either fails at that end
    std::istringstream edges("170141183460469231731687303715884105727 -170141183460469231731687303715884105728");
    edges >> value;
    EXPECT_TRUE(value == std::numeric_limits<int128_t>::max() && !edges.fail());
    edges >> value;
    EXPECT_TRUE(value == std::numeric_limits<int128_t>::min() && !edges.fail());
    std::istringstream above("170141183460469231731687303715884105728");
    above >> value;
    EXPECT_TRUE(above.fail() && value == std::numeric_limits<int128_t>::max());
    std::istringstream below("-170141183460469231731687303715884105729");
    below >> value;
    EXPECT_TRUE(below.fail() && value == std::numeric_limits<int128_t>::min());

    for (const char* text : {"abc", "-", "0x", "+x"}) {
        int128_t bad(7);
        std::istringstream stream(text);
        stream >> bad;
        EXPECT_TRUE(stream.fail()) << text;
        EXPECT_TRUE(bad.is_zero()) << text;
    }

    // hexadecimal after the prefix, or in the base the stream is set to
    std::istringstream hex("0x1F -0x10");
    hex >> value;
    EXPECT_TRUE(value == int128_t(31));
    hex >> value;
    EXPECT_TRUE(value == int128_t(-16));
    std::istringstream based("ff 17");
    based >> std::hex >> value;
    EXPECT_TRUE(value == int128_t(255));
    based >> std::oct >> value;
    EXPECT_TRUE(value == int128_t(15));

    uint128_t unsigned_value(7);
    std::istringstream unsigned_text("340282366920938463463374607431768211455 340282366920938463463374607431768211456");
    unsigned_text >> unsigned_value;
    EXPECT_TRUE(unsigned_value == std::numeric_limits<uint128_t>::max() && !unsigned_text.fail());
    unsigned_text >> unsigned_value;
    EXPECT_TRUE(unsigned_text.fail() && unsigned_value == std::numeric_limits<uint128_t>::max());
    // a minus sign wraps an unsigned value around, as strtoull does
    std::istringstream wrapped("-1");
    wrapped >> unsigned_value;
    EXPECT_TRUE(unsigned_value == std::numeric_limits<uint128_t>::max() && !wrapped.fail());
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
/// @brief Counts the distinct hashes of values that differ only in the high QWORD of their encoding.
template <typename T, typename Make> static size_t DistinctHashesOfHighWords(Make make)
{
    std::unordered_set<size_t> hashes;
    for (uint64_t i = 1; i <= 4096; ++i) {
        hashes.insert(std::hash<T> {}(make(i)));
    }
    return hashes.size();
}
// The hash tests above insert values that differ in their low QWORD, so a hash that ignored the
// high one, or a constant one, would pass them. The algorithm is fixed, so the count is exact.
TEST(int128_t, HashDistinguishesHighWord)
{
    EXPECT_EQ(DistinctHashesOfHighWords<uint128_t>([](uint64_t i) { return uint128_t(5, i); }), 4096u);
    EXPECT_EQ(DistinctHashesOfHighWords<int128_t>([](uint64_t i) { return int128_t(5, i << 40); }), 4096u);
    EXPECT_EQ(DistinctHashesOfHighWords<fixed_point128<32>>([](uint64_t i) { return fixed_point128<32>(5, i); }), 4096u);
    EXPECT_EQ(DistinctHashesOfHighWords<float128>([](uint64_t i) { return float128(5, HighQword(0, 0x3FFF, i)); }), 4096u);
    // and equal values hash equally
    EXPECT_EQ(std::hash<int128_t> {}(int128_t("-5")), std::hash<int128_t> {}(int128_t(-5)));
    EXPECT_EQ(std::hash<float128> {}(float128("0.5")), std::hash<float128> {}(float128::half()));
}
TEST(fixed_point128, NumericLimits)
{
    using fp = fixed_point128<32>;
    using limits = std::numeric_limits<fp>;

    static_assert(limits::is_specialized && limits::is_signed);
    // A fixed point value is exactly the number it stands for, unlike a floating point one.
    static_assert(limits::is_exact && !limits::is_integer);
    // one of the 128 bits is the sign, the other 127 carry magnitude
    static_assert(limits::digits == 127 && limits::digits10 == 38 && limits::radix == 2);
    static_assert(!limits::has_infinity && !limits::has_quiet_NaN && !limits::has_signaling_NaN);
    static_assert(limits::is_bounded && !limits::is_iec559 && !limits::tinyness_before);
    FP128_SUPPRESS_DEPRECATED_BEGIN
    static_assert(limits::has_denorm == std::denorm_absent && !limits::has_denorm_loss);
    FP128_SUPPRESS_DEPRECATED_END
    static_assert(limits::min_exponent == 32 - 126 && limits::max_exponent == 32);
    static_assert(limits::min_exponent10 == -28 && limits::max_exponent10 == 9);
    static_assert(limits::max_digits10 == 39 && std::numeric_limits<fixed_point128<10>>::max_digits10 == 40);
    static_assert(std::numeric_limits<const fp>::min_exponent == limits::min_exponent);

    // Every value is a constant expression, as the standard requires of a specialization.
    static_assert(limits::min() == fp::epsilon() && limits::epsilon() == fp::epsilon());
    static_assert(limits::lowest() < fp(0) && fp(0) < limits::max());
    static_assert(limits::round_error() == fp::half());
    static_assert(limits::infinity() == fp(0) && limits::quiet_NaN() == fp(0) && limits::signaling_NaN() == fp(0));
    // no value is subnormal, so the smallest one of any kind is the smallest positive one
    static_assert(limits::denorm_min() == limits::min());

    // Multiplication, division and the right shift round to nearest, a tie away from zero.
    // Truncating would leave the quotients at one epsilon, not two, and a tie at two, not three.
    static_assert(limits::round_style == std::round_to_nearest);
    const fp eps = limits::epsilon();
    EXPECT_TRUE((eps * 7) * (fp::half() >> 1) == eps * 2);  // 1.75 epsilon
    EXPECT_TRUE((eps * 5) / fp(3) == eps * 2);               // 1.67 epsilon
    EXPECT_TRUE(-(eps * 5) / fp(3) == -(eps * 2));
    EXPECT_TRUE((eps * 5) / fp(4) == eps);                   // 1.25 epsilon
    EXPECT_TRUE((eps * 7) / fp(4) == eps * 2);               // 1.75 epsilon
    EXPECT_TRUE((eps * 11) / fp(7) == eps * 2);              // 1.57 epsilon
    EXPECT_TRUE((eps * 10) / fp(7) == eps);                  // 1.43 epsilon
    EXPECT_TRUE((eps * 5) / fp(2) == eps * 3);               // a tie, away from zero
    EXPECT_TRUE(-(eps * 5) / fp(2) == -(eps * 3));
    EXPECT_TRUE((eps * 5) * fp::half() == eps * 3);
    EXPECT_TRUE(((eps * 5) >> 1) == eps * 3);
    EXPECT_TRUE(((-(eps * 5)) >> 1) == -(eps * 3));
    EXPECT_TRUE(((eps * 6) >> 2) == eps * 2);  // 1.5, a tie
    EXPECT_TRUE(((eps * 5) >> 2) == eps);      // 1.25

    // A zero divisor throws, which is the type's trap.
    static_assert(limits::traps);
    EXPECT_THROW((void)(fp::one() / fp(0)), std::logic_error);

    // +, - and * wrap around at the ends of the range, as is_modulo says.
    static_assert(limits::is_modulo);
    EXPECT_TRUE(limits::max() + eps == limits::lowest());
    EXPECT_TRUE(limits::lowest() - eps == limits::max());
    EXPECT_TRUE(limits::max() * fp(2) == -(eps * 2));

    // The grid spacing is the same everywhere, which is the whole point of the representation.
    const fixed_point128<32> one = fixed_point128<32>::one();
    EXPECT_TRUE(one + limits::epsilon() > one);
    EXPECT_TRUE(limits::min() == limits::epsilon());
    // Two's complement makes the range asymmetric: lowest() is one step below -max() and is its
    // own negation, having no positive counterpart.
    EXPECT_TRUE(limits::lowest() == -limits::max() - limits::epsilon());
    EXPECT_TRUE(-limits::lowest() == limits::lowest());
    EXPECT_TRUE(limits::lowest() < -limits::max());
}
/**
 * @brief Checks the exponent and decimal digit members of numeric_limits<fixed_point128<I>>
 *        against their definitions rather than against the formulas that compute them.
 *
 * float128 holds the powers of ten exactly: 10^k needs k * log2(5) bits of significand, about 91
 * for the largest k used here, well within its 113; scaling one by a power of two is exact too.
 */
template <int32_t I> static void CheckFixedPointLimitsAgainstDefinitions()
{
    using fp = fixed_point128<I>;
    using limits = std::numeric_limits<fp>;
    const fp one = fp::one();

    // 2^(min_exponent - 1) is the smallest positive value. The shift drops no set bit, so it is exact.
    EXPECT_TRUE((one >> (1 - limits::min_exponent)) == limits::min()) << "I=" << I;
    // 2^(max_exponent - 1) is in range and max() is one epsilon short of twice that.
    const fp top = one << (limits::max_exponent - 1);
    EXPECT_TRUE(limits::max() - top == top - limits::epsilon()) << "I=" << I;

    const auto ten_to = [](int32_t k) {
        float128 power(1);
        for (int32_t i = 0; i < k; ++i) {
            power *= float128(10);
        }
        return power;
    };
    // max_exponent10: 10^max_exponent10 <= max() < 10^(max_exponent10 + 1), max() being 2^I - 2^(I-127)
    const float128 max_value = ldexp(float128::one(), I) - ldexp(float128::one(), I - 127);
    EXPECT_TRUE(ten_to(limits::max_exponent10) <= max_value) << "I=" << I;
    EXPECT_TRUE(ten_to(limits::max_exponent10 + 1) > max_value) << "I=" << I;
    // min_exponent10: 10^min_exponent10 >= epsilon > 10^(min_exponent10 - 1), checked as
    // 1 >= epsilon * 10^-min_exponent10 and 10 > epsilon * 10^-min_exponent10, both exact.
    const float128 scaled_epsilon = ldexp(ten_to(-limits::min_exponent10), I - 127);
    EXPECT_TRUE(scaled_epsilon <= float128::one()) << "I=" << I;
    EXPECT_TRUE(scaled_epsilon * float128(10) > float128::one()) << "I=" << I;

    // max_digits10: in the top decade, from 10^max_exponent10 up, n significant digits are
    // 10^(max_exponent10 + 1 - n) apart. They tell neighbours apart only when that is less than the
    // grid spacing 2^(I - 127), that is when 10^(n - 1 - max_exponent10) exceeds 2^(127 - I). This
    // has to hold for max_digits10 and fail for one digit fewer. With max_exponent10 checked above,
    // this pins both terms of the definition.
    const float128 inverse_spacing = ldexp(float128(1), 127 - I);
    EXPECT_TRUE(ten_to(limits::max_digits10 - 2 - limits::max_exponent10) < inverse_spacing) << "I=" << I;
    EXPECT_TRUE(ten_to(limits::max_digits10 - 1 - limits::max_exponent10) > inverse_spacing) << "I=" << I;
    // every value's text form has at most max_digits10 significant digits
    EXPECT_LE(SignificantDigits(static_cast<std::string>(limits::max())), limits::max_digits10) << "I=" << I;

    static_assert(limits::is_modulo && limits::is_bounded && !limits::is_iec559 && limits::traps);
}
/// @brief Runs CheckFixedPointLimitsAgainstDefinitions for I = 1 + each of Offsets.
template <int32_t... Offsets> static void CheckFixedPointLimitsForEveryI(std::integer_sequence<int32_t, Offsets...>)
{
    (CheckFixedPointLimitsAgainstDefinitions<Offsets + 1>(), ...);
}
// Whether max_digits10 is 39 or 40 depends on I, so every one of the 63 instantiations is checked.
TEST(fixed_point128, NumericLimitsForEveryI)
{
    CheckFixedPointLimitsForEveryI(std::make_integer_sequence<int32_t, 63>{});
}
TEST(fixed_point128, FormatAndStream)
{
    const fixed_point128<32> value("-1234.5");

    EXPECT_EQ(std::format("{}", value), "-1234.5");
    EXPECT_EQ(std::format("{:>15}", fixed_point128<32>(42)), "             42");
    EXPECT_EQ(std::format("{:*^10}", fixed_point128<32>(42)), "****42****");
    EXPECT_EQ(std::format("{:+}", fixed_point128<32>(42)), "+42");
    EXPECT_EQ(std::format("{:08}", fixed_point128<32>("-1.5")), "-00001.5");
    EXPECT_EQ(std::format("{:<8}", fixed_point128<32>("-1.5")), "-1.5    ");
    EXPECT_EQ(std::format("{: }", fixed_point128<32>("1.5")), " 1.5");

    std::ostringstream stream;
    stream << value;
    EXPECT_EQ(stream.str(), "-1234.5");

    stream.str("");
    stream << std::setw(10) << std::left << fixed_point128<32>("1.5");
    EXPECT_EQ(stream.str(), "1.5       ");
    stream.str("");
    stream << std::setw(8) << std::internal << std::setfill('0') << std::showpos << fixed_point128<32>("1.5");
    EXPECT_EQ(stream.str(), "+00001.5");

    fixed_point128<32> read_back;
    std::istringstream input("2.5");
    input >> read_back;
    EXPECT_TRUE(read_back == fixed_point128<32>("2.5"));
}
// Like the integer extractors: the number only, and failbit for no digits or a value out of range,
// which stores the nearer end of the range. "xyz" used to read as zero, silently.
TEST(fixed_point128, StreamExtractionFailsOutOfRange)
{
    using fp = fixed_point128<32>;
    using limits = std::numeric_limits<fp>;
    fp value;
    char separator = 0;
    std::istringstream csv("2.5,-3");
    csv >> value >> separator;
    EXPECT_TRUE(value == fp("2.5"));
    EXPECT_EQ(separator, ',');
    csv >> value;
    EXPECT_TRUE(value == fp(-3) && !csv.fail());

    std::istringstream edges("-4294967296 4294967295.99999999999999999999999999997");
    edges >> value;
    EXPECT_TRUE(value == limits::lowest() && !edges.fail());
    edges >> value;
    EXPECT_TRUE(value == limits::max() && !edges.fail());
    // 2^32 is just out of range, and so is a value that rounds to it
    for (const char* text : {"4294967296", "4294967295.9999999999999999999999999999999", "99999999999999999999999"}) {
        std::istringstream stream(text);
        stream >> value;
        EXPECT_TRUE(stream.fail() && value == limits::max()) << text;
    }
    std::istringstream below("-4294967296.0000000000000000000000000001");
    below >> value;
    EXPECT_TRUE(below.fail() && value == limits::lowest());

    for (const char* text : {"xyz", "-", "."}) {
        fp bad(7);
        std::istringstream stream(text);
        stream >> bad;
        EXPECT_TRUE(stream.fail()) << text;
        EXPECT_TRUE(bad == fp(0)) << text;
    }
}
/// @brief Checks that random values and the extremes of fixed_point128<I> survive the trip through
///        their own text form, value -> string -> value, exactly.
template <int32_t I> static void CheckFixedPointStringRoundTrip()
{
    using fp = fixed_point128<I>;
    using limits = std::numeric_limits<fp>;
    srand(RANDOM_SEED);
    for (auto i = 0u; i < RANDOM_TEST_COUNT / 16; ++i) {
        const uint64_t low = get_uint64_random();
        const uint64_t high = get_uint64_random() >> (get_uint32_random() % 64);  // small magnitudes too
        const uint64_t sign = get_uint32_random() & 1;
        const fp value(low, sign ? ~high : high);
        const std::string text = static_cast<std::string>(value);
        EXPECT_TRUE(fp(text.c_str()) == value) << "I=" << I << ", text=" << text;
    }
    for (const fp& value : {limits::epsilon(), -limits::epsilon(), limits::max(), limits::lowest(), fp::one() - limits::epsilon()}) {
        const std::string text = static_cast<std::string>(value);
        EXPECT_TRUE(fp(text.c_str()) == value) << "I=" << I << ", text=" << text;
    }
}
// The text form had F / 3.3 decimal places, one short of what tells neighbours apart for most I,
// and the string constructor summed a series of rounded powers of ten: 1509 of 2012 random values
// of fixed_point128<32> did not come back. The expected text is from Python's decimal module, the
// exact value rounded half up to ceil(F log10 2) places.
TEST(fixed_point128, StringRoundTrip)
{
    CheckFixedPointStringRoundTrip<1>();
    CheckFixedPointStringRoundTrip<4>();
    CheckFixedPointStringRoundTrip<16>();
    CheckFixedPointStringRoundTrip<32>();
    CheckFixedPointStringRoundTrip<63>();

    using limits32 = std::numeric_limits<fixed_point128<32>>;
    EXPECT_EQ(static_cast<std::string>(limits32::epsilon()), "0.00000000000000000000000000003");
    EXPECT_EQ(static_cast<std::string>(-limits32::epsilon()), "-0.00000000000000000000000000003");
    EXPECT_EQ(static_cast<std::string>(limits32::max()), "4294967295.99999999999999999999999999997");
    EXPECT_EQ(static_cast<std::string>(limits32::lowest()), "-4294967296");
    EXPECT_EQ(static_cast<std::string>(fixed_point128<32>::one() - limits32::epsilon()), "0.99999999999999999999999999997");
    EXPECT_EQ(static_cast<std::string>(std::numeric_limits<fixed_point128<1>>::max()), "1.99999999999999999999999999999999999999");
    EXPECT_EQ(static_cast<std::string>(std::numeric_limits<fixed_point128<4>>::epsilon()), "0.00000000000000000000000000000000000009");
    EXPECT_EQ(static_cast<std::string>(std::numeric_limits<fixed_point128<16>>::max()), "65535.9999999999999999999999999999999996");
    EXPECT_EQ(static_cast<std::string>(std::numeric_limits<fixed_point128<63>>::max()), "9223372036854775807.99999999999999999995");
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
