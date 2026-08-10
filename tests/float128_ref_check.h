#ifndef FP128_FLOAT128_REF_CHECK_H
#define FP128_FLOAT128_REF_CHECK_H

// Checks float128 math results against the correctly rounded binary128 references in
// float128_ref_data.h, measured in units in the last place.
//
// The rest of the suite compares against the double the CRT returns, which can only pin down the
// top 53 of the 113 mantissa bits: a result that is wrong from bit 54 down passes it. These
// helpers close that gap. A bound is stated per function rather than shared, so a regression in
// one function cannot hide behind the loosest one in the set, and so the bound doubles as the
// documented accuracy of the implementation.
//
// Set the FP128_PRINT_ULP environment variable to have every check print the error it measured,
// which is how the bounds below were arrived at and how a change to an algorithm is evaluated.

#include <cstdlib>
#include <string>
#include "gtest_shared.h"
#include "float128_ref_data.h"

namespace f128_ref
{

/// @brief The binary128 encoding with the sign bit cleared, as a 128 bit integer.
[[nodiscard]] inline uint128_t Magnitude(const float128& x) noexcept
{
    uint64_t l = 0, h = 0;
    x.get_bits(l, h);
    return uint128_t(l, h & ~(1ull << 63));
}

/// @brief Number of representable values between a and b.
///
/// The encoding is sign and magnitude, and consecutive encodings of the same sign are consecutive
/// representable values, so a difference of magnitudes counts the steps directly. Across the sign
/// the two runs meet at zero, which makes the count the sum rather than the difference, and puts
/// +0 and -0 zero steps apart as IEEE 754 comparison requires.
///
/// @param a First value, must be finite or infinite but not NaN
/// @param b Second value, same
/// @return Steps from a to b, which cannot overflow: each magnitude is below 2^127.
[[nodiscard]] inline uint128_t UlpDistance(const float128& a, const float128& b) noexcept
{
    const uint128_t ma = Magnitude(a);
    const uint128_t mb = Magnitude(b);
    if (a.get_sign() == b.get_sign())
        return (ma > mb) ? ma - mb : mb - ma;
    return ma + mb;
}

/// @brief Saturates a 128 bit ulp count to 64 bit so a failure message can print it.
[[nodiscard]] inline uint64_t Saturate64(const uint128_t& v) noexcept
{
    uint64_t l = 0, h = 0;
    v.get_components(l, h);
    return (h != 0) ? UINT64_MAX : l;
}

/// @brief The raw 128 bit encoding as a hex string, for failure messages.
[[nodiscard]] inline std::string Bits(const float128& x)
{
    uint64_t l = 0, h = 0;
    x.get_bits(l, h);
    char buf[40];
    snprintf(buf, sizeof(buf), "%016llX%016llX", static_cast<unsigned long long>(h), static_cast<unsigned long long>(l));
    return buf;
}

/// @brief True when the harness was asked to report the error it measured.
[[nodiscard]] inline bool ShouldPrint()
{
    // getenv rather than _dupenv_s: the latter is MSVC only, and the string is neither stored nor
    // modified here. The CRT deprecation warning is turned off for the test project already.
    static const bool print = std::getenv("FP128_PRINT_ULP") != nullptr;
    return print;
}

/// @brief Accumulates the worst error seen over a table and reports it once.
class UlpTracker
{
public:
    explicit UlpTracker(const char* name) noexcept : funcName(name) {}

    /// @brief Compares one result against its reference, keeping the worst error so far.
    /// @param actual Value the implementation produced
    /// @param expected Correctly rounded reference
    /// @param index Position in the table, reported when this case turns out to be the worst
    /// @return True when the pair could be compared, false when a special value mismatched.
    bool Check(const float128& actual, const float128& expected, size_t index)
    {
        // A reference that is NaN or infinite is an exact requirement, not a rounding question.
        if (isnan(expected) || isnan(actual) || isinf(expected) || isinf(actual)) {
            if (isnan(expected) && isnan(actual))
                return true;
            if (expected == actual && !isnan(expected))
                return true;
            ADD_FAILURE() << funcName << "[" << index << "]: expected " << Bits(expected) << ", got " << Bits(actual);
            return false;
        }

        const uint128_t distance = UlpDistance(actual, expected);
        if (distance > worstDistance) {
            worstDistance = distance;
            worstIndex = index;
            worstActual = actual;
            worstExpected = expected;
        }
        return true;
    }

    /// @brief Fails the test when the worst error exceeded the bound, and reports where.
    /// @param bound Largest error to accept, in ulp
    void Expect(uint64_t bound) const
    {
        const uint64_t worst = Saturate64(worstDistance);
        if (ShouldPrint())
            printf("[ ULP      ] %-12s max %llu (bound %llu)\n", funcName, static_cast<unsigned long long>(worst),
                   static_cast<unsigned long long>(bound));

        EXPECT_LE(worst, bound) << funcName << ": worst case at index " << worstIndex << ", expected " << Bits(worstExpected) << ", got "
                                << Bits(worstActual);
    }

private:
    // The three 16 byte aligned members come first so that no padding has to be inserted between
    // them and the scalars, which is what MSVC reports as C4324.
    uint128_t worstDistance = 0;  ///< Largest error seen so far, in ulp.
    float128 worstActual {};      ///< Value the implementation produced there.
    float128 worstExpected {};    ///< Reference it was compared against.
    const char* funcName;         ///< Name of the function under test, for messages.
    size_t worstIndex = 0;        ///< Table index that produced it.
};

/// @brief Runs a single argument function over its reference table.
/// @tparam N Table size, deduced
/// @tparam Fn Callable taking a float128
/// @param name Function name, used in failure messages
/// @param table Reference cases
/// @param fn The function under test
/// @param bound Largest error to accept, in ulp
template <size_t N, typename Fn> inline void CheckUnary(const char* name, const ref_unary (&table)[N], Fn fn, uint64_t bound)
{
    UlpTracker tracker(name);
    for (size_t i = 0; i < N; ++i) {
        tracker.Check(fn(float128(table[i].xl, table[i].xh)), float128(table[i].rl, table[i].rh), i);
    }
    tracker.Expect(bound);
}

/// @brief Runs a two argument function over its reference table.
/// @copydetails CheckUnary
template <size_t N, typename Fn> inline void CheckBinary(const char* name, const ref_binary (&table)[N], Fn fn, uint64_t bound)
{
    UlpTracker tracker(name);
    for (size_t i = 0; i < N; ++i) {
        tracker.Check(fn(float128(table[i].xl, table[i].xh), float128(table[i].yl, table[i].yh)), float128(table[i].rl, table[i].rh), i);
    }
    tracker.Expect(bound);
}

/// @brief Runs a three argument function over its reference table.
/// @copydetails CheckUnary
template <size_t N, typename Fn> inline void CheckTernary(const char* name, const ref_ternary (&table)[N], Fn fn, uint64_t bound)
{
    UlpTracker tracker(name);
    for (size_t i = 0; i < N; ++i) {
        tracker.Check(fn(float128(table[i].xl, table[i].xh), float128(table[i].yl, table[i].yh), float128(table[i].zl, table[i].zh)),
                      float128(table[i].rl, table[i].rh), i);
    }
    tracker.Expect(bound);
}

}  // namespace f128_ref

#endif  // FP128_FLOAT128_REF_CHECK_H
