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

/***********************************************************************************
                                Acknologements
    The function div_128bit is derived from the book "Hacker's Delight" 2nd Edition
    by Henry S. Warren Jr. It implements Knuth's "Algorithm D", specialized to a 128 bit
    divisor and written in 64 bit limbs.

    The functions log, log2, log10 are derived from Dan Moulding's code:
    https://github.com/dmoulding/log2fix

    The function sqrt is based on the book "Math toolkit for real time programming"
    by Jack W. Crenshaw. The sin/cos/atan functions use some ideas from the book.

************************************************************************************/

/**
 * @file float128.h
 * @brief 128-bit IEEE 754-2008 quadruple-precision floating-point type.
 *
 * Provides the @ref fp128::float128 class with the same bit layout as
 * IEEE 754 binary128: 112-bit fraction, 15-bit exponent, and 1-bit sign.
 * Supports full arithmetic, comparison, and conversion operators, as well as
 * a comprehensive set of math functions (sqrt, cbrt, sin, cos, tan, asin,
 * acos, atan, atan2, sinh, cosh, tanh, asinh, acosh, atanh, exp, exp2, log,
 * log2, log10, pow, erf, erfc, etc.).
 * All methods are inline for maximum performance.
 *
 * The arithmetic, the conversions and the operations of IEEE 754-2008 clause 5 are correctly
 * rounded. The rounding directions and exception flags the standard also requires are opt in,
 * through FP128_IEEE_ENV; see the floating point environment section below.
 *
 * @see fp128_shared.h for supporting intrinsics and utilities.
 */

#ifndef FP128_FLOAT128_H
#define FP128_FLOAT128_H

#include <algorithm>
#include <array>
#include <cfenv>    // FE_ rounding directions and exception flags, which the float128 environment reuses
#include <cfloat>   // LDBL_MANT_DIG, which selects the long double conversion
#include <cmath>    // FP_NAN and friends, and the double seeds of sqrt/cbrt
#include <climits>  // INT_MAX, the ilogb answer for an infinity
#include <cstdlib>  // strtoull, for the nan() payload
#include <charconv>   // to_chars_result, from_chars_result, chars_format
#include <format>     // std::formatter
#include <functional> // std::hash
#include <istream>
#include <limits>
#include <ostream>
#include "fp128_shared.h"
#include "fp128_decimal.h"
#include "uint128_t.h"

namespace fp128
{
/***********************************************************************************
 *                                  Forward declarations
 ************************************************************************************/

class fp128_gtest;  // Google test class
class float128;

/// @brief Reads a decimal string into a float128. Defined below, declared here because the
///        constructor from a string delegates to it.
std::from_chars_result from_chars(const char* first, const char* last, float128& value);
/// @brief Reads a decimal or hexadecimal string into a float128. Defined below.
std::from_chars_result from_chars(const char* first, const char* last, float128& value, std::chars_format fmt);

namespace detail
{
/// @brief Renders a float128 as text. Defined below, declared here because the conversions to a
///        string delegate to it.
[[nodiscard]] std::string render(const float128& value, char style, int32_t precision, bool uppercase, bool alternate, char sign_char);
}  // namespace detail

/// @brief User-defined literal for constructing float128 from a string.
float128 operator""_f128(const char*);

/// @name CRT-Style Math Functions (Forward Declarations)
/// @brief Free functions providing standard math library equivalents for float128.
/// @{
constexpr float128 fabs(const float128& x) noexcept;
constexpr float128 floor(const float128& x) noexcept;
constexpr float128 ceil(const float128& x) noexcept;
constexpr float128 trunc(const float128& x) noexcept;
constexpr float128 round(const float128& x) noexcept;
constexpr int64_t llrint(const float128& x) noexcept;
constexpr int64_t llround(const float128& x) noexcept;
constexpr int32_t lrint(const float128& x) noexcept;
constexpr int32_t lround(const float128& x) noexcept;
constexpr int32_t ilogb(const float128& x) noexcept;
constexpr float128 copysign(const float128& x, const float128& y) noexcept;
float128 fmod(const float128& x, const float128& y) noexcept;
constexpr float128 modf(const float128& x, float128* iptr) noexcept;
constexpr float128 fdim(const float128& x, const float128& y) noexcept;
constexpr float128 fmin(const float128& x, const float128& y) noexcept;
constexpr float128 fmax(const float128& x, const float128& y) noexcept;
constexpr float128 fma(float128 x, float128 y, float128 z) noexcept;
float128 hypot(const float128& x, const float128& y) noexcept;
float128 cbrt(const float128 x, uint32_t iterations = 2) noexcept;
float128 sqrt(const float128& x, uint32_t iterations = 3) noexcept;
float128 erf(float128 x) noexcept;
float128 erfc(float128 x) noexcept;
float128 sin(float128 x) noexcept;
float128 asin(float128 x) noexcept;
float128 cos(float128 x) noexcept;
float128 acos(float128 x) noexcept;
float128 tan(float128 x) noexcept;
float128 atan(float128 x) noexcept;
float128 atan2(float128 y, float128 x) noexcept;
float128 sinh(float128 x) noexcept;
float128 asinh(float128 x) noexcept;
float128 cosh(float128 x) noexcept;
float128 acosh(float128 x) noexcept;
float128 tanh(float128 x) noexcept;
float128 atanh(float128 x) noexcept;
float128 exp(const float128& x) noexcept;
float128 exp2(const float128& x) noexcept;
float128 expm1(const float128& x) noexcept;
float128 pow(const float128& x, const float128& y) noexcept;
float128 pow(const float128& x, int32_t y) noexcept;
float128 log(float128 x) noexcept;
float128 log2(float128 x) noexcept;
float128 log10(float128 x) noexcept;
float128 logb(float128 x) noexcept;
float128 log1p(float128 x) noexcept;
constexpr float128 frexp(float128 x, int* expptr) noexcept;
constexpr float128 ldexp(float128 x, int exp) noexcept;
constexpr bool isfinite(const float128& x) noexcept;
constexpr bool signbit(const float128& x) noexcept;
constexpr bool isnormal(const float128& x) noexcept;
constexpr int fpclassify(const float128& x) noexcept;
constexpr bool isunordered(const float128& x, const float128& y) noexcept;
constexpr bool isgreater(const float128& x, const float128& y) noexcept;
constexpr bool isgreaterequal(const float128& x, const float128& y) noexcept;
constexpr bool isless(const float128& x, const float128& y) noexcept;
constexpr bool islessequal(const float128& x, const float128& y) noexcept;
constexpr bool islessgreater(const float128& x, const float128& y) noexcept;
constexpr float128 rint(const float128& x) noexcept;
constexpr float128 nearbyint(const float128& x) noexcept;
constexpr float128 scalbn(const float128& x, int n) noexcept;
constexpr float128 scalbln(const float128& x, long n) noexcept;
constexpr float128 nextafter(const float128& x, const float128& y) noexcept;
constexpr float128 nexttoward(const float128& x, const float128& y) noexcept;
float128 remainder(const float128& x, const float128& y) noexcept;
float128 remquo(const float128& x, const float128& y, int* quo) noexcept;
float128 hypot(const float128& x, const float128& y, const float128& z) noexcept;
float128 lgamma(float128 x) noexcept;
float128 tgamma(float128 x) noexcept;
constexpr float128 abs(const float128& x) noexcept;
constexpr bool totalorder(const float128& x, const float128& y) noexcept;
constexpr bool totalordermag(const float128& x, const float128& y) noexcept;
constexpr float128 fmaxmag(const float128& x, const float128& y) noexcept;
constexpr float128 fminmag(const float128& x, const float128& y) noexcept;
constexpr bool iscanonical(const float128& x) noexcept;
constexpr bool issignaling(const float128& x) noexcept;
constexpr bool issubnormal(const float128& x) noexcept;
constexpr bool iszero(const float128& x) noexcept;
/// The only one of these that a call cannot reach through argument dependent lookup: its argument
/// is a plain character pointer, so the name has to be visible at namespace scope - and a call has
/// to be qualified as fp128::nan(...) wherever `<cmath>` puts its own nan(const char*) in scope too.
float128 nan(const char* payload) noexcept;
/// @}

/// @name Non-CRT Utility Functions (Forward Declarations)
/// @{
float128 reciprocal(const float128& x) noexcept;
void fact_reciprocal(int x, float128& res) noexcept;
float128 double_factorial(int x) noexcept;
/// @}

/***********************************************************************************
 *                                  Floating point environment
 ************************************************************************************/

/**
 * @name Floating point environment
 *
 * IEEE 754-2008 asks two things of a binary format beyond its arithmetic. Clause 4 requires four
 * rounding-direction attributes the program can choose between: to nearest with ties to even, the
 * default, and towards zero, towards positive infinity and towards negative infinity. Clause 7
 * requires five exceptions - invalid operation, division by zero, overflow, underflow and inexact -
 * each reported by a status flag that stays raised until the program lowers it.
 *
 * float128 provides both, as a per thread environment of its own that is entirely separate from the
 * hardware one <cfenv> controls: setting the rounding direction here does not change how double
 * rounds, and an exception raised by a float128 operation shows up here and not in fetestexcept().
 *
 * <B>The environment is opt in.</B> It exists only when every translation unit that includes this
 * header is compiled with FP128_IEEE_ENV defined. Without it, which is the default, every operation
 * rounds to nearest with ties to even, no flag is ever raised, and the code that would maintain the
 * environment compiles away, leaving the arithmetic exactly as fast as it was before the
 * environment existed. The functions below are still declared, so the same source builds either
 * way: they report the fixed default and refuse to change it.
 *
 * The macro changes the definition of inline functions, so it must have the same value in every
 * translation unit of a program. Mixing the two violates the one definition rule.
 *
 * The rounding directions and exception flags are named with the FE_ constants of <cfenv>, so the
 * calls read the way their hardware counterparts do:
 * @code
 *     fp128::fesetround(FE_UPWARD);
 *     float128 bound = a / b;                       // rounded towards +inf
 *     if (fp128::fetestexcept(FE_INEXACT)) { ... }  // the quotient was not exact
 * @endcode
 * The calls have to be qualified: their arguments are plain integers, so argument dependent lookup
 * cannot tell them from the <cfenv> functions of the same name.
 * @{
 */

namespace detail
{
/// @brief The rounding-direction attributes IEEE 754-2008 requires of a binary format.
enum class rounding : uint8_t {
    nearest_even,  ///< roundTiesToEven, the default
    toward_zero,   ///< roundTowardZero
    upward,        ///< roundTowardPositive
    downward       ///< roundTowardNegative
};

/// @brief The per thread state behind the float128 environment.
struct float_env {
    rounding mode = rounding::nearest_even;  ///< Rounding direction in effect.
    int flags = 0;                           ///< Raised exceptions, as FE_ bits.
};

#ifdef FP128_IEEE_ENV
/// @brief The environment of the calling thread.
inline thread_local float_env env;
#endif

/**
 * @brief The rounding direction every float128 operation applies.
 *
 * Constant evaluation always rounds to nearest: a thread local cannot be read there, and a constant
 * whose value depended on the thread that happened to compile it would be meaningless anyway.
 *
 * Never store the result in a const local: `const rounding mode = current_rounding();` makes the
 * initializer manifestly constant evaluated, so it would always see nearest_even. Pass the call
 * as an argument, or keep the local non const.
 *
 * @return The current rounding direction, or nearest_even without FP128_IEEE_ENV.
 */
[[nodiscard]] FP128_FORCE_INLINE constexpr rounding current_rounding() noexcept
{
#ifdef FP128_IEEE_ENV
    if (!std::is_constant_evaluated())
        return env.mode;
#endif
    return rounding::nearest_even;
}

/**
 * @brief Raises exception flags. Does nothing without FP128_IEEE_ENV or during constant evaluation.
 * @param flags FE_ bits to raise
 */
FP128_FORCE_INLINE constexpr void raise_flags([[maybe_unused]] int flags) noexcept
{
#ifdef FP128_IEEE_ENV
    if (!std::is_constant_evaluated())
        env.flags |= flags;
#endif
}

/**
 * @brief Sets the exception flags aside for its lifetime, and discards whatever is raised meanwhile.
 *
 * For computations whose intermediate operations raise exceptions the result does not: the trial
 * parses of the shortest decimal search, or the Newton steps of sqrt(), which are inexact even when
 * the root they converge on is exact. Compiles to nothing without FP128_IEEE_ENV.
 */
class flag_quiet
{
public:
    flag_quiet() noexcept
    {
#ifdef FP128_IEEE_ENV
        saved = env.flags;
#endif
    }
    ~flag_quiet()
    {
#ifdef FP128_IEEE_ENV
        env.flags = saved;
#endif
    }
    flag_quiet(const flag_quiet&) = delete;
    flag_quiet& operator=(const flag_quiet&) = delete;

private:
    [[maybe_unused]] int saved = 0;  ///< Flags raised before the scope began.
};

/**
 * @brief Reduces what a math function's internal operations raise to what its result justifies.
 *
 * A function such as sin() is built from dozens of additions and multiplications, and those raise
 * exceptions of their own: a series term underflows, an intermediate overflows on a path whose
 * result does not. IEEE 754-2008 9.2 asks the function to signal as one operation would - inexact
 * when the result is, overflow when the result overflowed, underflow when it is tiny and inexact.
 * Constructed after a function has dealt with its special cases, which raise invalid and division
 * by zero themselves, it sets the caller's flags aside; operator() then keeps inexact if any step
 * raised it and derives the rest from the result. Compiles to nothing without FP128_IEEE_ENV.
 */
class flag_filter
{
public:
    flag_filter() noexcept
    {
#ifdef FP128_IEEE_ENV
        saved = env.flags;
        env.flags = 0;
#endif
    }
    ~flag_filter()
    {
#ifdef FP128_IEEE_ENV
        env.flags = saved | kept;
#endif
    }
    flag_filter(const flag_filter&) = delete;
    flag_filter& operator=(const flag_filter&) = delete;

    /**
     * @brief The same as operator(), for a result known to be inexact without an operation
     *        having said so: an overflow or underflow decided from the argument's range alone.
     * @tparam T float128
     * @param result The function's result
     * @return result
     */
    template <typename T> [[nodiscard]] T inexact(const T& result) noexcept
    {
        raise_flags(FE_INEXACT);
        return (*this)(result);
    }

    /**
     * @brief Records the flags a result justifies and passes the result through.
     * @tparam T float128, which is incomplete where this is declared
     * @param result The function's result, from finite operands
     * @return result
     */
    template <typename T> [[nodiscard]] T operator()(const T& result) noexcept
    {
#ifdef FP128_IEEE_ENV
        const int raised = env.flags;
        int justified = raised & FE_INEXACT;
        if (result.is_nan())
            justified |= FE_INVALID;
        else if (result.is_inf())
            justified |= FE_OVERFLOW | FE_INEXACT;
        else if ((result.is_zero() || result.is_subnormal()) && (raised & FE_INEXACT) != 0)
            justified |= FE_UNDERFLOW;
        kept |= justified;
#endif
        return result;
    }

private:
    [[maybe_unused]] int saved = 0;  ///< Flags raised before the function began.
    [[maybe_unused]] int kept = 0;   ///< Flags the result justified.
};

/**
 * @brief The digit rounding the current direction comes down to for a value of a given sign.
 * @param sign Sign of the value being converted
 * @return How to round its last digit.
 */
[[nodiscard]] inline digit_rounding digit_rounding_for(uint32_t sign) noexcept
{
    switch (current_rounding()) {
    case rounding::toward_zero: return digit_rounding::toward_zero;
    case rounding::upward:      return (sign != 0) ? digit_rounding::toward_zero : digit_rounding::away_from_zero;
    case rounding::downward:    return (sign != 0) ? digit_rounding::away_from_zero : digit_rounding::toward_zero;
    default:                    return digit_rounding::nearest_even;
    }
}
}  // namespace detail

/**
 * @brief The rounding direction float128 arithmetic uses on the calling thread.
 * @return One of FE_TONEAREST, FE_TOWARDZERO, FE_UPWARD or FE_DOWNWARD.
 */
[[nodiscard]] inline int fegetround() noexcept
{
    switch (detail::current_rounding()) {
    case detail::rounding::toward_zero: return FE_TOWARDZERO;
    case detail::rounding::upward:      return FE_UPWARD;
    case detail::rounding::downward:    return FE_DOWNWARD;
    default:                            return FE_TONEAREST;
    }
}

/**
 * @brief Selects the rounding direction float128 arithmetic uses on the calling thread.
 * @param round One of FE_TONEAREST, FE_TOWARDZERO, FE_UPWARD or FE_DOWNWARD
 * @return Zero on success. Non zero for an unknown direction, and for anything but FE_TONEAREST
 *         when the program was built without FP128_IEEE_ENV.
 */
inline int fesetround(int round) noexcept
{
    detail::rounding mode = detail::rounding::nearest_even;
    if (round == FE_TOWARDZERO)
        mode = detail::rounding::toward_zero;
    else if (round == FE_UPWARD)
        mode = detail::rounding::upward;
    else if (round == FE_DOWNWARD)
        mode = detail::rounding::downward;
    else if (round != FE_TONEAREST)
        return 1;

#ifdef FP128_IEEE_ENV
    detail::env.mode = mode;
    return 0;
#else
    return (mode == detail::rounding::nearest_even) ? 0 : 1;
#endif
}

/**
 * @brief Tests float128 exception flags (IEEE 754 testFlags).
 * @param excepts FE_ bits to test
 * @return The subset of @p excepts that is raised. Always zero without FP128_IEEE_ENV.
 */
[[nodiscard]] inline int fetestexcept([[maybe_unused]] int excepts) noexcept
{
#ifdef FP128_IEEE_ENV
    return detail::env.flags & excepts & FE_ALL_EXCEPT;
#else
    return 0;
#endif
}

/**
 * @brief Lowers float128 exception flags (IEEE 754 lowerFlags).
 * @param excepts FE_ bits to clear
 * @return Zero.
 */
inline int feclearexcept([[maybe_unused]] int excepts) noexcept
{
#ifdef FP128_IEEE_ENV
    detail::env.flags &= ~excepts;
#endif
    return 0;
}

/**
 * @brief Raises float128 exception flags (IEEE 754 raiseFlags).
 * @param excepts FE_ bits to raise
 * @return Zero on success; non zero without FP128_IEEE_ENV, when there is nowhere to raise them.
 */
inline int feraiseexcept(int excepts) noexcept
{
#ifdef FP128_IEEE_ENV
    detail::env.flags |= excepts & FE_ALL_EXCEPT;
    return 0;
#else
    return (excepts & FE_ALL_EXCEPT) != 0 ? 1 : 0;
#endif
}

/**
 * @brief Saves the state of float128 exception flags (IEEE 754 saveAllFlags).
 * @param flagp Receives the saved state
 * @param excepts FE_ bits to save
 * @return Zero.
 */
inline int fegetexceptflag(std::fexcept_t* flagp, int excepts) noexcept
{
    if (flagp != nullptr)
        *flagp = static_cast<std::fexcept_t>(fetestexcept(excepts));
    return 0;
}

/**
 * @brief Restores float128 exception flags saved by fegetexceptflag (IEEE 754 restoreFlags).
 *
 * Only the flags named in @p excepts are written, each set to the state it had when saved.
 *
 * @param flagp State saved by fegetexceptflag()
 * @param excepts FE_ bits to restore
 * @return Zero.
 */
inline int fesetexceptflag([[maybe_unused]] const std::fexcept_t* flagp, [[maybe_unused]] int excepts) noexcept
{
#ifdef FP128_IEEE_ENV
    if (flagp != nullptr) {
        const int mask = excepts & FE_ALL_EXCEPT;
        detail::env.flags = (detail::env.flags & ~mask) | (static_cast<int>(*flagp) & mask);
    }
#endif
    return 0;
}
/// @}

namespace detail
{
/**
 * @brief The binary expansion of 2/pi after the binary point, most significant word first.
 *
 * The trigonometric functions reduce a large argument by multiplying it with a window of these
 * bits chosen by its exponent (see float128::reduce_large()). An argument near the top of the
 * binary128 range reads about 16650 bits in. Generated by tools/gen_two_over_pi.py, which also
 * verifies it; do not edit by hand.
 */
inline constexpr uint64_t two_over_pi_bits[] = {
    0xA2F9836E4E441529, 0xFC2757D1F534DDC0, 0xDB6295993C439041, 0xFE5163ABDEBBC561,
    0xB7246E3A424DD2E0, 0x06492EEA09D1921C, 0xFE1DEB1CB129A73E, 0xE88235F52EBB4484,
    0xE99C7026B45F7E41, 0x3991D639835339F4, 0x9C845F8BBDF9283B, 0x1FF897FFDE05980F,
    0xEF2F118B5A0A6D1F, 0x6D367ECF27CB09B7, 0x4F463F669E5FEA2D, 0x7527BAC7EBE5F17B,
    0x3D0739F78A5292EA, 0x6BFB5FB11F8D5D08, 0x56033046FC7B6BAB, 0xF0CFBC209AF4361D,
    0xA9E391615EE61B08, 0x6599855F14A06840, 0x8DFFD8804D732731, 0x06061556CA73A8C9,
    0x60E27BC08C6B47C4, 0x19C367CDDCE8092A, 0x8359C4768B961CA6, 0xDDAF44D15719053E,
    0xA5FF07053F7E33E8, 0x32C2DE4F98327DBB, 0xC33D26EF6B1E5EF8, 0x9F3A1F35CAF27F1D,
    0x87F121907C7C246A, 0xFA6ED5772D30433B, 0x15C614B59D19C3C2, 0xC4AD414D2C5D000C,
    0x467D862D71E39AC6, 0x9B0062337CD2B497, 0xA7B4D55537F63ED7, 0x1810A3FC764D2A9D,
    0x64ABD770F87C6357, 0xB07AE715175649C0, 0xD9D63B3884A7CB23, 0x24778AD623545AB9,
    0x1F001B0AF1DFCE19, 0xFF319F6A1E666157, 0x9947FBACD87F7EB7, 0x652289E83260BFE6,
    0xCDC4EF09366CD43F, 0x5DD7DE16DE3B5892, 0x9BDE2822D2E88628, 0x4D58E232CAC616E3,
    0x08CB7DE050C017A7, 0x1DF35BE01834132E, 0x6212830148835B8E, 0xF57FB0ADF2E91E43,
    0x4A48D36710D8DDAA, 0x425FAECE616AA428, 0x0AB499D3F2A6067F, 0x775C83C2A3883C61,
    0x78738A5A8CAFBDD7, 0x6F63A62DCBBFF4EF, 0x818D67C12645CA55, 0x36D9CAD2A8288D61,
    0xC277C9121426049B, 0x4612C459C444C5C8, 0x91B24DF31700AD43, 0xD4E5492910D5FDFC,
    0xBE00CC941EEECE70, 0xF53E1380F1ECC3E7, 0xB328F8C79405933E, 0x71C1B3092EF3450B,
    0x9C12887B20AB9FB5, 0x2EC292472F327B6D, 0x550C90A7721FE76B, 0x96CB314A1679E279,
    0x4189DFF49794E884, 0xE6E29731996BED88, 0x365F5F0EFDBBB49A, 0x486CA46742727132,
    0x5D8DB8159F09E5BC, 0x25318D3974F71C05, 0x30010C0D68084B58, 0xEE2C90AA4702E774,
    0x24D6BDA67DF77248, 0x6EEF169FA6948EF6, 0x91B45153D1F20ACF, 0x3398207E4BF56863,
    0xB25F3EDD035D407F, 0x8985295255C06437, 0x10D86D324832754C, 0x5BD4714E6E5445C1,
    0x090B69F52AD56614, 0x9D072750045DDB3B, 0xB4C576EA17F9877D, 0x6B49BA271D296996,
    0xACCCC65414AD6AE2, 0x9089D98850722CBE, 0xA4049407777030F3, 0x27FC00A871EA49C2,
    0x663DE06483DD9797, 0x3FA3FD94438C860D, 0xDE41319D39928C70, 0xDDE7B7173BDF082B,
    0x3715A0805C93805A, 0x921110D8E80FAF80, 0x6C4BFFDB0F903876, 0x185915A562BBCB61,
    0xB989C7BD401004F2, 0xD2277549F6B6EBBB, 0x22DBAA140A2F2689, 0x768364333B091A94,
    0x0EAA3A51C2A31DAE, 0xEDAF12265C4DC26D, 0x9C7A2D9756C0833F, 0x03F6F0098C402B99,
    0x316D07B43915200C, 0x5BC3D8C492F54BAD, 0xC6A5CA4ECD37A736, 0xA9E69492AB6842DD,
    0xDE6319EF8C76528B, 0x6837DBFCABA1AE31, 0x15DFA1AE00DAFB0C, 0x664D64B705ED3065,
    0x29BF56573AFF47B9, 0xF96AF3BE75DF9328, 0x3080ABF68C6615CB, 0x040622FA1DE4D9A4,
    0xB33D8F1B5709CD36, 0xE9424EA4BE13B523, 0x331AAAF0A8654FA5, 0xC1D20F3F0BCD785B,
    0x76F923048B7B7217, 0x8953A6C6E26E6F00, 0xEBEF584A9BB7DAC4, 0xBA66AACFCF761D02,
    0xD12DF1B1C1998C77, 0xADC3DA4886A05DF7, 0xF480C62FF0AC9AEC, 0xDDBC5C3F6DDED01F,
    0xC790B6DB2A3A25A3, 0x9AAF009353AD0457, 0xB6B42D297E804BA7, 0x07DA0EAA76A1597B,
    0x2A12162DB7DCFDE5, 0xFAFEDB89FDBE896C, 0x76E4FCA90670803E, 0x156E85FF87FD073E,
    0x2833676186182AEA, 0xBD4DAFE7B36E6D8F, 0x3967955BBF3148D7, 0x8416DF30432DC735,
    0x6125CE70C9B8CB30, 0xFD6CBFA200A4E46C, 0x05A0DD5A476F21D2, 0x1262845CB9496170,
    0xE0566B0152993755, 0x50B7D51EC4F1335F, 0x6E13E4305DA92E85, 0xC3B21D3632A1A4B7,
    0x08D4B1EA21F716E4, 0x698F77FF2780030C, 0x2D408DA0CD4F99A5, 0x20D3A2B30A5D2F42,
    0xF9B4CBDA11D0BE7D, 0xC1DB9BBD17AB81A2, 0xCA5C6A0817552E55, 0x0027F0147F8607E1,
    0x640B148D4196DEBE, 0x872AFDDAB6256B34, 0x897BFEF3059EBFB9, 0x4F6A68A82A4A5AC4,
    0x4FBCF82D985AD795, 0xC7F48D4D0DA63A20, 0x5F57A4B13F149538, 0x800120CC86DD71B6,
    0xDEC9F560BF11654D, 0x6B0701ACB08CD0C0, 0xB24855510EFB1EC3, 0x72953B06A33540C0,
    0x7BDC06CC45E0FA29, 0x4EC8CAD641F3E8DE, 0x647CD8649B31BED9, 0xC397A4D45877C5E3,
    0x6913DAF03C3ABA46, 0x18465F7555F5BDD2, 0xC6926E5D2EACED44, 0x0E423E1C87C461E9,
    0xFD29F3D6E7CA7C22, 0x35916FC5E0088DD7, 0xFFE26A6EC6FDB0C1, 0x0893745D7CB2AD6B,
    0x9D6ECD7B723E6A11, 0xC6A9CFF7DF7329BA, 0xC9B55100B70DB2E2, 0x24BA74607DE58AD8,
    0x742C150D0C188194, 0x667E162901767A9F, 0xBEFDFDEF4556367E, 0xD913D9ECB9BA8BFC,
    0x97C427A831C36EF1, 0x36C59456A8D8B5A8, 0xB40ECCCF2D891234, 0x576F89562CE3CE99,
    0xB920D6AA5E6B9C2A, 0x3ECC5F114A0BFDFB, 0xF4E16D3B8E2C86E2, 0x84D4E9A9B4FCD1EE,
    0xEFC9352E61392F44, 0x2138C8D91B0AFC81, 0x6A4AFBD81C2F84B4, 0x538C994ECC2254DC,
    0x552AD6C6C096190B, 0xB8701A649569605A, 0x26EE523F0F117F11, 0xB5F4F5CBFC2DBC34,
    0xEEBC34CC5DE8605E, 0xDD9B8E67EF3392B8, 0x17C99B5861BC57E1, 0xC68351103ED84871,
    0xDDDD1C2DA118AF46, 0x2C21D7F359987AD9, 0xC0549EFA864FFC06, 0x56AE79E536228922,
    0xAD38DC9367AAE855, 0x3826829BE7CAA40D, 0x51B133990ED7A948, 0x0569F0B265A7887F,
    0x974C8836D1F9B392, 0x214A827B21CF98DC, 0x9F405547DC3A74E1, 0x42EB67DF9DFE5FD4,
    0x5EA4677B7AACBAA2, 0xF65523882B55BA41, 0x086E59862A218347, 0x39E6E389D49EE540,
    0xFB49E956FFCA0F1C, 0x8A59C52BFA94C5C1, 0xD3CFC50FAE5ADB86, 0xC5476243853B8621,
    0x94792C8761107B4C, 0x2A1A2C8012BF4390, 0x2688893C78E4C4A8, 0x7BDBE5C23AC4EAF4,
    0x268A67F7BF920D2B, 0xA365B1933D0B7CBD, 0xDC51A463DD27DDE1, 0x6919949A9529A828,
    0xCE68B4ED09209F44, 0xCA984E638270237C, 0x7E32B90F8EF5A7E7, 0x561408F1212A9DB5,
    0x4D7E6F5119A5ABF9, 0xB5D6DF8261DD9602, 0x36169F3AC4A1A283, 0x6DED727A8D39A9B8,
    0x825C326B5B2746ED, 0x34007700D255F4FC, 0x4D59018071E0E13F, 0x89B295F364A8F1AE,
    0xA74B38FC4CEAB2BB, 0x47270BABC3A734BA, 0x6052DD34F8563AEB, 0x7E8A31BB365895B7,
};
}  // namespace detail

/***********************************************************************************
 *                                  Main Code
 ************************************************************************************/

/**
 * @brief 128 bit floating point class.
 *
 * This class implements the standard operators a floating point data type.<BR>
 * All of float128's methods are inline for maximum performance.
 *
 * <B>IEEE 754-2008:</B>
 * <UL>
 * <LI>Addition, subtraction, multiplication, division, sqrt(), fma(), the remainders, the
 *     rounding to an integral value and the conversions to and from double, float, long double,
 *     the integer types and text are correctly rounded - each computes the exact result and rounds
 *     it once, to the destination's width, its subnormal range included (see round_pack()).</LI>
 * <LI>Every one of them follows the standard's rules for zeros, infinities and NaNs: the sign of
 *     an exact zero sum, the invalid operations, a signaling NaN quieted, the payload of the first
 *     NaN operand delivered.</LI>
 * <LI>Round to nearest, ties to even, is the only rounding direction and no exception flag is
 *     kept, unless FP128_IEEE_ENV is defined; with it fp128::fesetround() chooses among the four
 *     directions clause 4 requires, and the five flags of clause 7 are raised and tested through
 *     fp128::fetestexcept() and its siblings. Only with the macro does the type claim conformance
 *     (numeric_limits::is_iec559, is754version2008()).</LI>
 * <LI>The math functions of clause 9 follow its special values but are not correctly rounded,
 *     which the standard recommends rather than requires. Each one's error is measured and stated
 *     in tests/float128_accuracy_gtest.cpp.</LI>
 * </UL>
 *
 * <B>Implementation notes:</B>
 * <UL>
 * <LI>Same bit layout as binary128: 112 bit fraction, 15 bit exponent and 1 bit for the sign.</LI>
 * <LI>A float128 object is not thread safe. Accessing a const object from multiple threads is safe.</LI>
 * <LI>Only 64 bit builds are supported.</LI>
 * </UL>
 *
 * <B>Compile time evaluation:</B><BR>
 * Everything that stays within the 128 bit encoding is constexpr:
 * <UL>
 * <LI>Construction from any builtin arithmetic type and from the raw QWORDs, copy, move,
 *     assignment and the conversions to the integer and floating point types.</LI>
 * <LI>Addition, subtraction, multiplication, square(), the shifts, the unary operators and the
 *     comparisons.</LI>
 * <LI>The queries and the component accessors: is_zero(), is_normal(), is_nan(), is_inf(),
 *     is_int(), get_exponent(), get_class(), get_components(), set_components() and the rest,
 *     along with every named constant (one(), pi(), inf(), nan(), ...).</LI>
 * <LI>The math functions fabs, floor, ceil, trunc, round, llround, lround, llrint, lrint, ilogb,
 *     copysign, modf, fdim, fmin, fmax, fma, frexp, ldexp, sqr, isnan, isinf, isfinite, and the
 *     nextUp/nextDown/exp10 helpers.</LI>
 * </UL>
 *
 * Two things make that possible. The bit counting and extended arithmetic intrinsics are not
 * constant expressions, so fp128_shared.h wraps each one in a constexpr function that
 * serves a constant evaluated call from a portable implementation of the same operation. And the
 * fields of the high QWORD are read through the shift and mask accessors rather than the
 * _float128_bits view, because reading the inactive member of a union is not allowed during
 * constant evaluation. A runtime call is unaffected by either.
 *
 * The rest cannot be constexpr, for one of two reasons:
 * <UL>
 * <LI>Division and modulo, and everything built on them. Nothing in div_128bit itself is barred
 *     from a constant expression, so what keeps operator/=() out of one is the surrounding code
 *     rather than the long division.</LI>
 * <LI>The string conversions allocate, and the transcendental functions parse their constants
 *     from strings held in function local statics, which a constexpr function may not declare.</LI>
 * </UL>
 */

class FP128_ALIGN16 float128
{
    // build time validation of template parameters
    static_assert(sizeof(void*) == 8, "float128 is supported in 64 bit builds only!");
    friend class fp128_gtest;

    static constexpr int32_t EXP_BITS = 15;            ///< Number of exponent bits in binary128.
    static constexpr int32_t EXP_BIAS = 0x3FFF;         ///< Exponent bias (16383) for binary128.
    static constexpr int32_t ZERO_EXP_BIASED = -EXP_BIAS;        ///< Biased exponent value representing zero.
    static constexpr int32_t ZERO_EXP_UNBIASED = 0;              ///< Unbiased exponent value for zero encoding.
    static constexpr int32_t SUBNORM_EXP_BIASED = 0;             ///< Biased exponent value for subnormal numbers.
    static constexpr int32_t SUBNORM_EXP_UNBIASED = -EXP_BIAS;   ///< Unbiased exponent value for subnormal numbers.
    static constexpr int32_t INF_EXP_BIASED = 0x7FFF;            ///< Biased exponent value for infinity/NaN.
    static constexpr int32_t INF_EXP_UNBIASED = INF_EXP_BIASED - EXP_BIAS; ///< Unbiased exponent value for infinity/NaN.
    static constexpr uint64_t EXP_MASK = INF_EXP_BIASED;         ///< Bitmask for the exponent field.
    static constexpr int32_t FRAC_BITS = 112;                    ///< Number of fraction (mantissa) bits.
    static constexpr int32_t EXP_SHIFT = FRAC_BITS - 64;         ///< Bit position of the exponent field within the high QWORD.
    static constexpr uint64_t UPPER_FRAC_MASK = FP128_MAX_VALUE_64(FRAC_BITS - 64); ///< Bitmask for upper fraction bits within the high QWORD.
    static constexpr uint64_t FRAC_UNITY = FP128_ONE_SHIFT(FRAC_BITS - 64);         ///< The implicit unity bit position in the high QWORD.
    static constexpr uint64_t SIGN_MASK = 1ull << 63;            ///< Bitmask for the sign bit.
    static constexpr uint64_t QUIET_NAN_BIT = 1ull << (EXP_SHIFT - 1);  ///< Leading fraction bit, set on a quiet NaN.

    /// @brief Bit-field view of the upper 64 bits of a float128.
    ///
    /// Retained for the debugger visualizer (fixed_point128.natvis), which reads the three fields
    /// by name. The code itself reaches the same fields through the shift and mask accessors
    /// below: reading the inactive member of a union is not allowed during constant evaluation,
    /// and every operation on this type funnels through those accessors.
    struct _float128_bits {
        uint64_t f : 48; ///< Upper 48 bits of the 112-bit fraction.
        uint64_t e : 15; ///< 15-bit biased exponent.
        uint64_t s : 1;  ///< Sign bit (0 = positive, 1 = negative).
    };

    /// @brief IEEE 754 classification categories for float128 values.
    enum float128_class_t {
        signalingNaN,       ///< Signaling NaN.
        quietNaN,           ///< Quiet NaN.
        negativeInfinity,   ///< Negative infinity.
        negativeNormal,     ///< Negative normal number.
        negativeSubnormal,  ///< Negative subnormal (denormalized) number.
        negativeZero,       ///< Negative zero.
        positiveZero,       ///< Positive zero.
        positiveSubnormal,  ///< Positive subnormal (denormalized) number.
        positiveNormal,     ///< Positive normal number.
        positiveInfinity    ///< Positive infinity.
    };

    /// @brief Internal storage: 128-bit value split into low and high QWORDs.
    struct {
        uint64_t low;  ///< Lower 64 bits of the float128 encoding.
        union {
            uint64_t high;              ///< Upper 64 bits (raw).
            _float128_bits high_bits;   ///< Upper 64 bits (bit-field view).
        };
    };

public:
    /**
     * @brief Default constructor, creates an instance with a value of zero.
     */
    FP128_FORCE_INLINE constexpr float128() noexcept : low(0), high(0) {}
    /**
     * @brief Copy constructor
     * @param rhs Object to copy from
     */
    FP128_FORCE_INLINE constexpr float128(const float128& rhs) noexcept : low(rhs.low), high(rhs.high) {}
    /**
     * @brief Move constructor
     * Doesn't modify the right hand side object. Acts like a copy constructor.
     * @param rhs Object to copy from
     */
    FP128_FORCE_INLINE constexpr float128(float128&& rhs) noexcept : low(rhs.low), high(rhs.high) {}
    /**
     * @brief Low level constructor
     * @param l Low QWORD
     * @param h High QWORD
     */
    FP128_FORCE_INLINE constexpr float128(uint64_t l, uint64_t h) noexcept : low(l), high(h) {}
    /**
     * @brief Low level constructor
     * @param lf Low fraction part (bits 63:0)
     * @param hf High fraction part (bits 111:64)
     * @param e Exponent
     * @param s sign
     */
    FP128_FORCE_INLINE constexpr float128(uint64_t lf, uint64_t hf, uint32_t e, uint32_t s) noexcept :
        low(lf), high((hf & FP128_MAX_VALUE_64(48)) | (((0x7FFFull & e) << 48)) | ((1ull & s) << 63))
    {
    }
    /**
     * @brief Constructor from the double type
     * @param x Input value
     */
    FP128_INLINE constexpr float128(double x) noexcept
    {
        low = high = 0;
        // hack the double bit fields
        const Double d(x);

        // very common case. The sign is kept: -0.0 is a value of its own.
        if (x == 0) {
            set_sign(d.s());
            return;
        }

        // subnormal numbers
        if (d.e() == 0) {
            // the exponent is -1022 (1-1023)
            auto msb = 64 - static_cast<int32_t>(lzcnt64(d.f()));
            // exponent
            int32_t x_expo = static_cast<int32_t>(d.e()) - 1023;
            int32_t expo = x_expo + msb - dbl_frac_bits;
            // fraction
            low = d.f() & ~(1ull << (msb - 1));  // clear the msb
            auto shift = static_cast<int32_t>(FRAC_BITS - msb + 1);
            shift_left128_inplace_safe(low, high, shift);
            set_exponent(expo);
        }
        // NaN & INF
        else if (d.e() == 0x7FF) {
            set_exponent_bits(INF_EXP_BIASED);
            if (d.f() != 0) {
                // A NaN keeps its payload, moved to the top of the wider fraction where the
                // conversion back to double looks for it, and its quiet bit lands on the quiet bit.
                // A signaling NaN is quieted: the conversion is an operation, and a signaling
                // operand makes it the invalid one.
                low = d.f() << (FRAC_BITS - dbl_frac_bits);
                set_fraction_bits(d.f() >> (64 - (FRAC_BITS - dbl_frac_bits)));
                if (is_signaling()) {
                    detail::raise_flags(FE_INVALID);
                    high |= QUIET_NAN_BIT;
                }
            }
        }
        // normal numbers
        else {
            low = d.f() << 60;
            high = d.f() >> 4;
            set_exponent(static_cast<int32_t>(d.e()) - 1023);
        }

        // copy the sign
        set_sign(d.s());
    }
    /**
     * @brief Generic constructor for integral and floating-point types.
     *
     * Every value of every builtin arithmetic type is representable, so the conversion is exact.
     * float goes through double, which holds it exactly. A long double wider than a double - the
     * x87 extended format, or binary128 itself - is read from its encoding, since a double would
     * round it. Integral types are converted directly. Character pointer types delegate to the
     * const char* constructor.
     *
     * @tparam T Source type (must be arithmetic or a character pointer type)
     * @param x Input value
     */
    template <typename T> constexpr float128(T x) noexcept : low(0), high(0)
    {
        if constexpr (std::is_floating_point_v<T>) {
            if constexpr (std::is_same_v<T, long double> && LDBL_MANT_DIG != 53)
                *this = from_long_double(x);
            else
                *this = float128(static_cast<double>(x));
            return;
        } else if constexpr (std::is_same_v<char*, T> || std::is_same_v<unsigned char*, T> || std::is_same_v<const unsigned char*, T>) {
            *this = float128(static_cast<const char*>(x));
            return;
        } else if constexpr (std::is_integral_v<T>) {
            uint64_t sign = 0;
            // The magnitude is taken in the unsigned domain, where the conversion sign extends and
            // the negation wraps. Negating the signed value instead is undefined for the most
            // negative one, which has no positive counterpart, and produces the same bit pattern
            // for every other value.
            low = static_cast<uint64_t>(x);
            if constexpr (std::is_signed_v<T>) {
                // always do positive multiplication
                if (x < 0) {
                    low = 0ull - low;
                    sign = 1;
                }
            }

            if (low == 0)
                return;

            auto expo = log2(low);  // this is the index of the msb as well
            auto shift = static_cast<int32_t>(FRAC_BITS - expo);
            shift_left128_inplace_safe(low, high, shift);
            set_sign(sign);
            set_exponent(static_cast<int32_t>(expo));
            return;
        }
    }
    /**
     * @brief Construct from a string
     * Allows creating very high precision values, approximately 34 decimal digits.
     * Much slower than the other constructors.
     * @param x Input string
     */
    float128(const char* x) noexcept
    {
        low = high = 0;
        if (x == nullptr)
            return;

        // Everything is read by from_chars(), which builds the exact value and rounds it once:
        // a decimal number, a hexadecimal one (0x1.8p1, the form printf's %a writes), or one of
        // inf, infinity, nan, nan(payload) and snan. Leading white space is skipped, as strtod()
        // skips it.
        //
        // The decimal path used to live here: it accumulated the digits as a float128 and divided
        // by a power of ten held in the type. Negative powers of ten are not representable in
        // binary, so that division rounded, and the value came back about an ulp away from the
        // one the string named. The hexadecimal path that remained here read integers only, so
        // 0x1.8p1 came back as 1, and it truncated the digits past the 29th rather than rounding.
        while (*x != 0 && isspace(static_cast<unsigned char>(*x)))
            ++x;
        const char* const end = x + strlen(x);
        const char* body = x;
        if (*body == '-' || *body == '+')
            ++body;
        if (body[0] == '0' && (body[1] == 'x' || body[1] == 'X')) {
            // std::from_chars takes a hexadecimal number without its prefix.
            if (from_chars(body + 2, end, *this, std::chars_format::hex).ec == std::errc {} && *x == '-')
                invert_sign();
            return;
        }
        from_chars(x, end, *this);
    }
    /**
     * @brief Constructor from std::string.
     * Allows creating very high precision values. Much slower than the other constructors.
     * @param x Input string
     */
    FP128_INLINE float128(const std::string& x) noexcept
    {
        // delegate to the char* c'tor
        *this = float128(x.c_str());
    }
    /**
     * @brief Destructor
     */
    constexpr ~float128() = default;

    /**
     * @brief Assignment operator
     * @param rhs Object to copy from
     * @return This object.
     */
    constexpr FP128_FORCE_INLINE float128& operator=(const float128& rhs) noexcept
    {
        high = rhs.high;
        low = rhs.low;
        return *this;
    }
    /**
     * @brief Move assignment operator
     * @param rhs Object to copy from
     * @return This object.
     */
    constexpr FP128_FORCE_INLINE float128& operator=(float128&& rhs) noexcept
    {
        high = rhs.high;
        low = rhs.low;
        return *this;
    }

    //
    // conversion operators
    //

    /**
     * @name Conversion to narrower floating point formats
     *
     * IEEE 754-2008 5.4.2 convertFormat: the value is rounded once, in the current direction, to
     * the destination's precision and range - including its subnormal range, which keeps fewer
     * bits than its normal one. A NaN keeps its sign and as much of its payload as fits, and is
     * quieted.
     * @{
     */

    /// @brief A finite value rounded to a narrower binary format, see narrow().
    struct narrowed {
        uint64_t sig = 0;      ///< Significand, below 2^P. Zero for a zero result.
        int32_t lsb = 0;       ///< Weight of the significand's lowest bit, as a power of two.
        uint32_t sign = 0;     ///< Sign of the result.
        bool to_inf = false;   ///< The value overflowed to infinity.
    };

    /**
     * @brief Rounds this finite, non zero value to a narrower binary format.
     *
     * Overflow produces infinity, or the largest finite value in the two directions that round
     * towards it, and raises the overflow and inexact flags. A result that is tiny after rounding
     * and inexact raises underflow.
     *
     * @tparam P Precision of the destination, counting its leading bit whether stored or implicit.
     *         At most 64, so the significand fits a single word.
     * @tparam EMIN Exponent of the destination's smallest normal value
     * @tparam EMAX Exponent of its largest finite value
     * @return The rounded significand and the weight of its lowest bit.
     */
    template <int32_t P, int32_t EMIN, int32_t EMAX> [[nodiscard]] FP128_INLINE constexpr narrowed narrow() const noexcept
    {
        static_assert(P >= 2 && P <= 64, "the significand has to fit in a QWORD");
        constexpr uint64_t all_ones = (P == 64) ? UINT64_MAX : ((1ull << P) - 1);

        narrowed res;
        uint64_t l = 0, h = 0;
        int32_t e = 0;
        get_components(l, h, e, res.sign);
        // Not const: a const local of enumeration type initialized from a constant expression is
        // usable in constant expressions, which makes its initializer manifestly constant evaluated
        // - std::is_constant_evaluated() would answer true there, and fix the direction to nearest.
        detail::rounding mode = detail::current_rounding();

        // The lowest bit kept is P-1 places below the leading one, or the destination's subnormal
        // floor when that is higher. The 113 bit significand's own lowest bit weighs 2^(e-112).
        res.lsb = ((e < EMIN) ? EMIN : e) - (P - 1);
        const int32_t drop = res.lsb - (e - FRAC_BITS);

        // Tininess is detected after rounding: a value just below the smallest normal one is not
        // tiny when rounding it to P bits with an unbounded exponent carries it up to there.
        bool tiny = e < EMIN;
        if (e == EMIN - 1) {
            uint64_t tl = l, th = h, textra = 0;
            shift_right_jam128_extra(tl, th, textra, FRAC_BITS + 1 - P);
            if (tl == all_ones && round_increment(mode, res.sign, textra))
                tiny = false;
        }

        uint64_t extra = 0;
        shift_right_jam128_extra(l, h, extra, drop);
        res.sig = l;  // at most P bits are left, and P is at most 64
        if (extra != 0) {
            detail::raise_flags(tiny ? (FE_INEXACT | FE_UNDERFLOW) : FE_INEXACT);
            if (round_increment(mode, res.sign, extra)) {
                ++res.sig;
                // A carry out of the top lands on the next power of two. It has to be read before
                // the tie is settled: for a 64 bit significand the carry shows as a wrap to zero,
                // and a tie that rounded zero up to one and back to even leaves a zero as well.
                const bool carry = (P == 64) ? (res.sig == 0) : ((res.sig >> (P % 64)) != 0);
                if (mode == detail::rounding::nearest_even && (extra << 1) == 0)
                    res.sig &= ~1ull;  // an exact tie goes to the even neighbour
                if (carry) {
                    res.sig = 1ull << (P - 1);
                    ++res.lsb;
                }
            }
        }

        if (res.lsb + (P - 1) > EMAX) {
            detail::raise_flags(FE_OVERFLOW | FE_INEXACT);
            res.to_inf = mode == detail::rounding::nearest_even || mode == ((res.sign != 0) ? detail::rounding::downward : detail::rounding::upward);
            res.sig = all_ones;
            res.lsb = EMAX - (P - 1);
        }
        return res;
    }

    /**
     * @brief The top of a NaN's payload, for a narrower format to carry.
     *
     * The payload is the fraction below the quiet bit, bits 110:0. A narrower format keeps its top
     * bits, which is where the conversion from that format puts them, so a payload survives the
     * round trip.
     *
     * @param shift FRAC_BITS - 1 less the payload width of the destination, at least 47
     * @return Bits 110:shift of the fraction, right aligned.
     */
    [[nodiscard]] FP128_FORCE_INLINE constexpr uint64_t nan_payload_top(int32_t shift) const noexcept
    {
        const uint64_t frac_high = get_fraction_bits() & (UPPER_FRAC_MASK >> 1);
        return (shift >= 64) ? (frac_high >> (shift - 64)) : ((frac_high << (64 - shift)) | (low >> shift));
    }

    /**
     * @brief Operator double
     *
     * Correctly rounded in the current direction. Overflow keeps its sign, as does a value too small
     * for the subnormal range, and a NaN keeps its sign and the top 51 bits of its payload.
     */
    [[nodiscard]] FP128_INLINE constexpr operator double() const noexcept
    {
        const uint64_t sign = static_cast<uint64_t>(get_sign()) << 63;
        if (is_special()) {
            if (is_nan()) {
                if (is_signaling())
                    detail::raise_flags(FE_INVALID);
                return std::bit_cast<double>(sign | (0x7FFull << 52) | (1ull << 51) | nan_payload_top(FRAC_BITS - 1 - 51));
            }
            return std::bit_cast<double>(sign | (0x7FFull << 52));
        }
        if (is_zero())
            return std::bit_cast<double>(sign);

        const narrowed n = narrow<53, -1022, 1023>();
        if (n.to_inf)
            return std::bit_cast<double>(sign | (0x7FFull << 52));
        // A normal significand's leading one adds the last unit to the exponent field, which is
        // why the field is written one short; a subnormal one has no leading one and no exponent.
        if (n.sig >> 52)
            return std::bit_cast<double>(sign + (static_cast<uint64_t>(n.lsb + 52 + 1023 - 1) << 52) + n.sig);
        return std::bit_cast<double>(sign | n.sig);
    }
    /**
     * @brief operator float converts to a float
     *
     * Rounded once, straight to single precision. Going through double, as this used to, rounds
     * twice: 1 + 2^-24 + 2^-82 became 1 + 2^-24 in double, a tie, which then went to 1.0f.
     */
    [[nodiscard]] FP128_INLINE constexpr operator float() const noexcept
    {
        const uint32_t sign = get_sign() << 31;
        if (is_special()) {
            if (is_nan()) {
                if (is_signaling())
                    detail::raise_flags(FE_INVALID);
                return std::bit_cast<float>(sign | 0x7F800000u | 0x00400000u | static_cast<uint32_t>(nan_payload_top(FRAC_BITS - 1 - 22)));
            }
            return std::bit_cast<float>(sign | 0x7F800000u);
        }
        if (is_zero())
            return std::bit_cast<float>(sign);

        const narrowed n = narrow<24, -126, 127>();
        if (n.to_inf)
            return std::bit_cast<float>(sign | 0x7F800000u);
        const uint32_t sig = static_cast<uint32_t>(n.sig);
        if (sig >> 23)
            return std::bit_cast<float>(sign + (static_cast<uint32_t>(n.lsb + 23 + 127 - 1) << 23) + sig);
        return std::bit_cast<float>(sign | sig);
    }
    /// @}

    /**
     * @name Conversion to integers
     *
     * The C++ conversions truncate towards zero, which is IEEE 754-2008's
     * convertToIntegerTowardZero. They used to round to nearest for a magnitude of one or more and
     * truncate below it, which agreed with neither: 0.75 became 0 but 1.75 became 2.
     *
     * A value outside the range of the integer type, or a NaN, is the invalid operation. The result
     * is then defined rather than left undefined as it is for the builtin types: an out of range
     * value saturates at the nearest limit, and a NaN converts to zero.
     * @{
     */

    /**
     * @brief The integer part of |x|, truncated towards zero.
     * @param overflow Set when |x| is 2^64 or more, an infinity included
     * @return The truncated magnitude, or UINT64_MAX on overflow.
     */
    [[nodiscard]] FP128_INLINE constexpr uint64_t truncated_magnitude(bool& overflow) const noexcept
    {
        uint64_t l = 0, h = 0;
        int32_t e = 0;
        uint32_t s = 0;
        get_components(l, h, e, s);
        overflow = e > 63;
        if (overflow)
            return UINT64_MAX;
        if (e < 0)
            return 0;
        // The units bit is bit 112 - e of the significand, which is between 49 and 112.
        const int32_t shift = FRAC_BITS - e;
        return (shift >= 64) ? (h >> (shift - 64)) : ((l >> shift) | (h << (64 - shift)));
    }

    /**
     * @brief Converts to an integer type, truncating towards zero.
     * @tparam I The integer type
     * @return The converted value, saturated when out of range, zero for a NaN.
     */
    template <typename I> [[nodiscard]] FP128_INLINE constexpr I to_integer() const noexcept
    {
        if (is_nan()) {
            detail::raise_flags(FE_INVALID);
            return 0;
        }

        const bool negative = get_sign() != 0;
        bool overflow = false;
        const uint64_t magnitude = truncated_magnitude(overflow);
        if constexpr (std::is_signed_v<I>) {
            // The negative range reaches one further than the positive one.
            const uint64_t limit = static_cast<uint64_t>(std::numeric_limits<I>::max()) + (negative ? 1u : 0u);
            if (overflow || magnitude > limit) {
                detail::raise_flags(FE_INVALID);
                return negative ? std::numeric_limits<I>::min() : std::numeric_limits<I>::max();
            }
            // Negated in the unsigned domain, which also covers the most negative value.
            return static_cast<I>(negative ? (0ull - magnitude) : magnitude);
        } else {
            // A negative value that truncates to zero converts fine; anything below that does not.
            if (negative && magnitude != 0) {
                detail::raise_flags(FE_INVALID);
                return 0;
            }
            if (overflow || magnitude > std::numeric_limits<I>::max()) {
                detail::raise_flags(FE_INVALID);
                return std::numeric_limits<I>::max();
            }
            return static_cast<I>(magnitude);
        }
    }

    /**
     * @brief operator uint64_t converts to a uint64_t, truncating towards zero
     */
    [[nodiscard]] FP128_INLINE constexpr operator uint64_t() const noexcept { return to_integer<uint64_t>(); }
    /**
     * @brief operator int64_t converts to a int64_t, truncating towards zero
     */
    [[nodiscard]] FP128_INLINE constexpr operator int64_t() const noexcept { return to_integer<int64_t>(); }
    /**
     * @brief operator uint32_t converts to a uint32_t, truncating towards zero
     */
    [[nodiscard]] FP128_INLINE constexpr operator uint32_t() const noexcept { return to_integer<uint32_t>(); }
    /**
     * @brief operator int32_t converts to a int32_t, truncating towards zero
     */
    [[nodiscard]] FP128_INLINE constexpr operator int32_t() const noexcept { return to_integer<int32_t>(); }
    /// @}

    /**
     * @name long double
     *
     * What a long double is depends on the platform: a plain double under MSVC, the x87 80 bit
     * extended format on x86 under GCC and Clang, and binary128 itself on AArch64 Linux. The two
     * wider ones are converted through their encoding rather than through double, which would
     * round away the bits they have over it.
     * @{
     */
#if LDBL_MANT_DIG == 64
    /// @brief The x87 extended format: an explicit 64 bit significand, then sign and exponent.
    struct x87_bits {
        uint64_t mantissa;                                ///< Significand, its leading bit stored
        uint16_t sign_exponent;                           ///< Sign in bit 15, biased exponent below it
        unsigned char padding[sizeof(long double) - 10];  ///< Unused, outside the value representation
    };
    static_assert(sizeof(x87_bits) == sizeof(long double), "unexpected x87 long double layout");
#elif LDBL_MANT_DIG == 113
    /// @brief long double is binary128: the same encoding, in two QWORDs.
    struct binary128_bits {
        uint64_t low;   ///< Low QWORD of the encoding
        uint64_t high;  ///< High QWORD of the encoding
    };
    static_assert(sizeof(binary128_bits) == sizeof(long double), "unexpected binary128 long double layout");
#endif

    /**
     * @brief Converts a long double wider than a double, exactly.
     * @param x Value to convert
     * @return The same value. A signaling NaN is quieted, the conversion being an operation.
     */
    [[nodiscard]] static FP128_INLINE constexpr float128 from_long_double(long double x) noexcept
    {
#if LDBL_MANT_DIG == 64
        const x87_bits bits = std::bit_cast<x87_bits>(x);
        const uint32_t sign = static_cast<uint32_t>(bits.sign_exponent >> 15);
        const int32_t biased = bits.sign_exponent & 0x7FFF;
        if (biased == 0x7FFF) {
            // An infinity has nothing but its leading bit set. Anything else is a NaN, whose
            // payload and quiet bit move to the top of the fraction.
            const uint64_t fraction = bits.mantissa & ~(1ull << 63);
            float128 res(fraction << 49, (static_cast<uint64_t>(sign) << 63) | (static_cast<uint64_t>(INF_EXP_BIASED) << EXP_SHIFT) | (fraction >> 15));
            if (res.is_signaling()) {
                detail::raise_flags(FE_INVALID);
                res = quiet_nan(res);
            }
            return res;
        }
        // A denormal is scaled like the smallest normal value, it only lacks the leading bit.
        return norm_round_pack(sign, ((biased == 0) ? 1 : biased) - EXP_BIAS - 63 + FRAC_BITS, bits.mantissa, 0, false);
#elif LDBL_MANT_DIG == 113
        const binary128_bits bits = std::bit_cast<binary128_bits>(x);
        return float128(bits.low, bits.high);
#else
        // A format this code does not know, such as the double-double of PowerPC: only the double
        // nearest to it is converted.
        return float128(static_cast<double>(x));
#endif
    }

    /**
     * @brief operator long double
     *
     * The same as the conversion to double where long double is one. Elsewhere correctly rounded to
     * the wider format, or exact where long double is binary128.
     *
     * @return Object value.
     */
    [[nodiscard]] FP128_INLINE constexpr operator long double() const noexcept
    {
#if LDBL_MANT_DIG == 64
        x87_bits bits {};
        const uint16_t sign = static_cast<uint16_t>(get_sign() << 15);
        if (is_special()) {
            bits.sign_exponent = static_cast<uint16_t>(sign | 0x7FFF);
            bits.mantissa = 1ull << 63;
            if (is_nan()) {
                if (is_signaling())
                    detail::raise_flags(FE_INVALID);
                bits.mantissa |= (1ull << 62) | nan_payload_top(FRAC_BITS - 1 - 62);
            }
        } else if (is_zero()) {
            bits.sign_exponent = sign;
        } else {
            // The exponent range is binary128's own, so only rounding the largest values up can
            // overflow; the subnormal range is narrower, its last bit being worth 2^-16445.
            const narrowed n = narrow<64, -16382, 16383>();
            if (n.to_inf) {
                bits.sign_exponent = static_cast<uint16_t>(sign | 0x7FFF);
                bits.mantissa = 1ull << 63;
            } else {
                bits.mantissa = n.sig;
                const int32_t biased = (n.sig >> 63) ? (n.lsb + 63 + EXP_BIAS) : 0;
                bits.sign_exponent = static_cast<uint16_t>(sign | static_cast<uint16_t>(biased));
            }
        }
        return std::bit_cast<long double>(bits);
#elif LDBL_MANT_DIG == 113
        return std::bit_cast<long double>(binary128_bits {low, high});
#else
        return operator double();
#endif
    }
    /// @}
    /**
     * @brief Converts to a std::string (slow) string holds all meaningful fraction bits.
     * @return object string representation
     */
    [[nodiscard]] FP128_INLINE operator std::string() const noexcept { return operator char*(); }
    /**
     * @brief Converts to a C string (slow) string holds all meaningful fraction bits.
     * @return object string representation
     */
    [[nodiscard]] explicit FP128_INLINE operator char*() const noexcept
    {
        // The buffer is thread local so that the pointer stays valid after the call returns, which
        // is what an implicit conversion to char* has to provide. It is large enough for the
        // longest shortest-form string: a sign, 36 digits, a point and a five digit exponent.
        static thread_local char str[64];
        const std::string text = detail::render(*this, '\0', -1, false, false, '-');
        const size_t length = (text.size() < sizeof(str) - 1) ? text.size() : sizeof(str) - 1;
        for (size_t i = 0; i < length; ++i)
            str[i] = text[i];
        str[length] = '\0';
        return str;
    }
    /**
     * @brief Writes the value in scientific notation into a caller supplied buffer.
     *
     * Kept for the callers that predate std::format support. Everything it does is also reachable
     * as fp128::to_chars(..., std::chars_format::scientific) or as std::format("{:e}", value),
     * both of which say what they produce more clearly.
     *
     * @param str Output buffer
     * @param buff_size Output buffer size in bytes
     */
    void to_e_format(char* str, int32_t buff_size) const
    {
        if (str == nullptr || buff_size < 1)
            return;

        const std::string text = detail::render(*this, 'e', -1, false, false, '-');
        const size_t length = (text.size() < static_cast<size_t>(buff_size) - 1) ? text.size() : static_cast<size_t>(buff_size) - 1;
        for (size_t i = 0; i < length; ++i)
            str[i] = text[i];
        str[length] = '\0';
    }

    //
    // math operators
    //
    /// @brief Largest scaling a shift or ldexp() ever needs: enough to take the largest finite value
    ///        below half the smallest subnormal, or the smallest subnormal above the largest finite
    ///        value. Larger counts are clamped to it, which keeps the exponent arithmetic in range.
    static constexpr int32_t SCALE_LIMIT = 2 * (EXP_BIAS + FRAC_BITS + 1);

    /**
     * @brief Shift right this object, dividing it by a power of two.
     * @param shift Bits to shift. Values less than 1 do nothing, high values can cause the value to reach zero.
     * @return This object.
     */
    FP128_INLINE constexpr float128& operator>>=(int32_t shift) noexcept
    {
        if (shift < 1)
            return *this;
        if (is_special()) {
            if (is_nan())
                *this = propagate_nan(*this);
            return *this;
        }

        uint64_t l, h;
        int32_t e;
        uint32_t s;
        get_components(l, h, e, s);
        e -= (shift < SCALE_LIMIT) ? shift : SCALE_LIMIT;
        // A subnormal or zero result is rounded once, here, in the current direction.
        set_components(l, h, e, s);
        return *this;
    }
    /**
     * @brief Shift left this object, multiplying it by a power of two.
     * @param shift Bits to shift. Values less than 1 do nothing, high values can cause the value to reach infinity.
     * @return This object.
     */
    FP128_INLINE constexpr float128& operator<<=(int32_t shift) noexcept
    {
        if (shift < 1)
            return *this;
        if (is_special()) {
            if (is_nan())
                *this = propagate_nan(*this);
            return *this;
        }

        uint64_t l, h;
        int32_t e;
        uint32_t s;
        get_components(l, h, e, s);
        e += (shift < SCALE_LIMIT) ? shift : SCALE_LIMIT;
        set_components(l, h, e, s);
        return *this;
    }
    /**
     * @brief Performs right shift operation.
     *
     * Hinted rather than forced, unlike the other one line forwarders (see FP128_FORCE_INLINE).
     * What it forwards to is not cheap: operator>>= splits the value into components and
     * reassembles it, which brings the subnormal handling of get_components() and
     * set_components() along with it. Forcing this pair open cost ~9% of the float128 Mandelbrot
     * benchmark, which shifts inside a loop that already inlines an add, a multiply and two
     * squares. Left as a hint, the shift is still inlined wherever the caller has room for it.
     *
     * @param shift bits to shift
     * @return Temporary object with the result of the operation
     */
    template <typename T> FP128_INLINE constexpr float128 operator>>(T shift) const noexcept
    {
        float128 temp(*this);
        return temp >>= static_cast<int32_t>(shift);
    }
    /**
     * @brief Performs left shift operation.
     * Hinted rather than forced, see operator>> above.
     * @param shift bits to shift
     * @return Temporary object with the result of the operation
     */
    template <typename T> FP128_INLINE constexpr float128 operator<<(T shift) const noexcept
    {
        float128 temp(*this);
        return temp <<= static_cast<int32_t>(shift);
    }

    /**
     * @brief Add a value to this object
     * @param rhs Right hand side operand
     * @return This object.
     */
    FP128_INLINE constexpr float128& operator+=(const float128& rhs) noexcept
    {
        // check trivial cases
        if (is_special() || rhs.is_special()) {
            // A NaN operand propagates. Adding infinities of opposite signs is the invalid
            // operation and produces a NaN as well.
            if (is_nan() || rhs.is_nan()) {
                *this = propagate_nan(*this, rhs);
                return *this;
            }
            if (is_inf() && rhs.is_inf() && get_sign() != rhs.get_sign()) {
                *this = invalid_operation();
                return *this;
            }

            // Either a single operand is infinite, or both are infinite with the same sign.
            // In both cases the result is that infinity, keeping its own sign.
            if (rhs.is_inf())
                *this = rhs;
            return *this;
        }

        // A zero operand leaves the other one exactly. Two zeros of opposite sign make an exact
        // zero sum, which IEEE 754 makes positive in every rounding direction but downward.
        if (rhs.is_zero()) {
            if (is_zero() && get_sign() != rhs.get_sign())
                set_sign(exact_zero_sign());
            return *this;
        }
        if (is_zero()) {
            *this = rhs;
            return *this;
        }

        uint32_t sign, rhs_sign;
        int32_t expo, rhs_expo;
        uint64_t l1, h1, l2, h2;
        get_components(l1, h1, expo, sign);
        rhs.get_components(l2, h2, rhs_expo, rhs_sign);

        // Both mantissas are moved up so their leading one sits at bit 125. That leaves the
        // alignment thirteen bits to shift the smaller operand into before anything falls off,
        // a bit above for the carry of an addition, and the top bit clear.
        constexpr int32_t room = 127 - 2 - FRAC_BITS;

        // Exponents this far apart leave the smaller operand below a quarter of the larger one's
        // last place, so rounding to nearest returns the larger one unchanged. The other
        // directions can still move the result by one place, so they take the full path.
        if (detail::current_rounding() == detail::rounding::nearest_even) {
            if (expo - rhs_expo > FRAC_BITS + room) {
                detail::raise_flags(FE_INEXACT);
                return *this;
            }
            if (rhs_expo - expo > FRAC_BITS + room) {
                detail::raise_flags(FE_INEXACT);
                *this = rhs;
                return *this;
            }
        }

        // Make the first operand the one with the larger exponent, so only the second one is
        // ever shifted right.
        if (expo < rhs_expo) {
            std::swap(l1, l2);
            std::swap(h1, h2);
            std::swap(expo, rhs_expo);
            std::swap(sign, rhs_sign);
        }

        const int32_t diff = expo - rhs_expo;
        shift_left128_inplace_safe(l1, h1, room);
        if (diff <= room) {
            shift_left128_inplace_safe(l2, h2, room - diff);
        } else {
            // What is shifted out is jammed into the lowest bit rather than rounded away. Rounding
            // here and again after the sum, as this used to, decided the direction from a window
            // of three bits and misjudged every case the bits further down would have settled.
            // The jammed bit lands at least eleven places under the guard bit of the sum, which
            // is enough for it to stand in for everything it replaced.
            shift_right_jam128(l2, h2, diff - room);
        }

        if (sign == rhs_sign) {
            const uint8_t carry = addcarryx_u64(0, l1, l2, &l1);
            addcarryx_u64(carry, h1, h2, &h1);
        } else {
            uint8_t borrow = subborrow_u64(0, l1, l2, &l1);
            borrow = subborrow_u64(borrow, h1, h2, &h1);
            // Only operands with equal exponents can leave a negative difference, in which case
            // the second operand was the larger and its sign wins.
            if (borrow != 0) {
                twos_complement128(l1, h1);
                sign = rhs_sign;
            }
            if ((l1 | h1) == 0) {
                *this = float128(0, static_cast<uint64_t>(exact_zero_sign()) << 63);
                return *this;
            }
        }

        *this = norm_round_pack(sign, expo - room, l1, h1, false);
        return *this;
    }
    /**
     * @brief Add a value to this object
     * @param rhs Right hand side operand
     * @return This object.
     */
    template <typename T> FP128_FORCE_INLINE constexpr float128& operator+=(const T& rhs) { return operator+=(float128(rhs)); }
    /**
     * @brief Subtract a value from this object
     * @param rhs Right hand side operand
     * @return This object.
     */
    FP128_FORCE_INLINE constexpr float128& operator-=(const float128& rhs) noexcept { return *this += (-rhs); }
    /**
     * @brief Subtract a value from this object
     * @param rhs Right hand side operand
     * @return This object.
     */
    template <typename T> FP128_FORCE_INLINE constexpr float128& operator-=(const T& rhs) { return operator+=(-float128(rhs)); }
    /**
     * @brief Multiply a value to this object
     * @param rhs Right hand side operand
     * @return This object.
     */
    FP128_INLINE constexpr float128& operator*=(const float128& rhs) noexcept
    {
        // check trivial cases
        if (is_special() || rhs.is_special()) {
            // A NaN operand propagates, and inf * zero is the invalid operation.
            // Note the zero tests are only reachable for the operand that is not special.
            if (is_nan() || rhs.is_nan()) {
                *this = propagate_nan(*this, rhs);
                return *this;
            }
            if ((is_inf() && rhs.is_zero()) || (is_zero() && rhs.is_inf())) {
                *this = invalid_operation();
                return *this;
            }

            // at least one operand is infinite, the result is infinite with the combined sign
            const uint32_t res_sign = get_sign() ^ rhs.get_sign();
            *this = inf();
            set_sign(res_sign);
            return *this;
        } else if (is_zero() || rhs.is_zero()) {
            // the sign of a zero product is the combination of both operand signs
            const uint32_t res_sign = get_sign() ^ rhs.get_sign();
            *this = 0;
            set_sign(res_sign);
            return *this;
        }
        // extract fractions and exponents
        uint32_t sign, rhs_sign;
        int32_t expo, rhs_expo;
        uint64_t l1, h1, l2, h2;
        get_components(l1, h1, expo, sign);
        rhs.get_components(l2, h2, rhs_expo, rhs_sign);
        // Deliberately asked of the objects rather than derived from the mantissas above, even
        // though get_components() has just produced the bits this needs. Testing the raw value
        // leaves the question independent of the extraction, so it can be answered alongside it;
        // phrasing it as (h1 == FRAC_UNITY && l1 == 0) chains it behind instead and costs Clang
        // 15% of the multiply benchmark, against no gain on MSVC.
        bool is_exp2 = is_exponent_of_2();
        bool rhs_is_exp2 = rhs.is_exponent_of_2();

        // add the exponents
        expo += rhs_expo;

        // optimize for exponents of 2
        if (is_exp2 || rhs_is_exp2) {
            // copy the fraction as needed
            if (is_exp2) {
                l1 = l2;
                h1 = h2;
            }

            set_components(l1, h1, expo, sign ^ rhs_sign);
            return *this;
        }

        // multiply the fractions
        // the fractions are in u16.112 precision
        // the result will be u32.224 precision and will be shifted-right by 112 bit
        uint64_t res[4];  // 256 bit of result
        uint64_t temp1[2], temp2[2];

        // multiply low QWORDs
        res[0] = mulx_u64(l1, l2, &res[1]);

        // multiply high QWORDs (overflow can happen)
        res[2] = mulx_u64(h1, h2, &res[3]);

        // multiply low this and high rhs
        temp1[0] = mulx_u64(l1, h2, &temp1[1]);
        uint8_t carry = addcarryx_u64(0, res[1], temp1[0], &res[1]);
        res[3] += addcarryx_u64(carry, res[2], temp1[1], &res[2]);

        // multiply high this and low rhs
        temp2[0] = mulx_u64(h1, l2, &temp2[1]);
        carry = addcarryx_u64(0, res[1], temp2[0], &res[1]);
        res[3] += addcarryx_u64(carry, res[2], temp2[1], &res[2]);

        // extract the bits from res[] keeping the precision the same as this object
        // shift result by F
        constexpr int32_t index = 1;
        constexpr int32_t lsb = (FRAC_BITS & 63) - 1;  // bit within the 64bit data pointed by res[index], one short so the guard bit stays

        const uint64_t sticky_bits = res[0] | (res[1] & FP128_MAX_VALUE_64(lsb));
        l1 = shift_right128(res[index], res[index + 1], lsb);  // custom function is 20% faster in Mandelbrot than the intrinsic
        h1 = shift_right128(res[index + 1], res[index + 2], lsb);
        --expo;

        const uint64_t extra = norm_product(l1, h1, expo, sticky_bits != 0);
        *this = round_pack(sign ^ rhs_sign, expo, l1, h1, extra);
        return *this;
    }
    /**
     * @brief Squares this object in place.
     *
     * Cheaper than operator*=(*this) on several counts:
     * - a square is symmetric, so the two cross products of the fraction multiply are the
     *   same value and only three 64x64->128 bit multiplies are needed instead of four
     * - the result sign is always positive, so no sign combining is needed
     * - only one operand has to be tested for the special cases and for being an exponent
     *   of 2, and in the exponent of 2 case the fraction does not have to be copied
     * - the exponents are added to each other, which is a doubling
     *
     * The result is identical to (*this) * (*this) for every value, specials included.
     *
     * @return This object.
     */
    FP128_INLINE constexpr float128& square() noexcept
    {
        // check trivial cases
        if (is_special()) {
            // a NaN stays a NaN, +/-inf squared is +inf
            *this = is_nan() ? propagate_nan(*this) : inf();
            return *this;
        } else if (is_zero()) {
            *this = 0;
            return *this;
        }
        // extract the fraction and exponent, the sign of a square is always positive
        uint32_t sign;
        int32_t expo;
        uint64_t l, h;
        get_components(l, h, expo, sign);

        // add the exponent to itself
        expo += expo;

        // optimize for exponents of 2, the fraction is unchanged by the multiply
        if (is_exponent_of_2()) {
            set_components(l, h, expo, 0);
            return *this;
        }

        // square the fraction
        // the fraction is in u16.112 precision
        // the result will be u32.224 precision and will be shifted-right by 112 bit
        uint64_t res[4];  // 256 bit of result
        uint64_t cross[2];

        // multiply the low QWORD by itself
        res[0] = mulx_u64(l, l, &res[1]);

        // multiply the high QWORD by itself (overflow can happen)
        res[2] = mulx_u64(h, h, &res[3]);

        // the low * high cross product, which appears twice in the sum
        cross[0] = mulx_u64(l, h, &cross[1]);

        uint8_t carry = addcarryx_u64(0, res[1], cross[0], &res[1]);
        res[3] += addcarryx_u64(carry, res[2], cross[1], &res[2]);

        carry = addcarryx_u64(0, res[1], cross[0], &res[1]);
        res[3] += addcarryx_u64(carry, res[2], cross[1], &res[2]);

        // extract the bits from res[] keeping the precision the same as this object
        constexpr int32_t index = 1;
        constexpr int32_t lsb = (FRAC_BITS & 63) - 1;  // bit within the 64bit data pointed by res[index] minus 1 to improve rounding

        const uint64_t sticky_bits = res[0] | (res[1] & FP128_MAX_VALUE_64(lsb));
        l = shift_right128(res[index], res[index + 1], lsb);
        h = shift_right128(res[index + 1], res[index + 2], lsb);
        --expo;

        const uint64_t extra = norm_product(l, h, expo, sticky_bits != 0);
        *this = round_pack(0, expo, l, h, extra);
        return *this;
    }
    /**
     * @brief Multiply a value to this object
     * @param rhs Right hand side operand
     * @return This object.
     */
    template <typename T> FP128_FORCE_INLINE constexpr float128& operator*=(const T& rhs) { return operator*=(float128(rhs)); }
    /**
     * @brief Divide this object by a value
     * @param rhs Right hand side operand
     * @return This object.
     */
    FP128_INLINE float128& operator/=(const float128& rhs)
    {
        // check trivial cases
        // A NaN operand propagates, and zero / zero and inf / inf are the invalid operation.
        if (is_nan() || rhs.is_nan()) {
            *this = propagate_nan(*this, rhs);
            return *this;
        }
        if ((is_zero() && rhs.is_zero()) || (is_inf() && rhs.is_inf())) {
            *this = invalid_operation();
            return *this;
        }

        // the sign of any of the results below is the combination of both operand signs
        const uint32_t res_sign = get_sign() ^ rhs.get_sign();

        // An infinite dividend or a zero divisor produce an infinity, a zero dividend or an
        // infinite divisor produce a zero. The combinations where both apply, which are the
        // two invalid operations, were already handled above. Only a finite dividend over a zero
        // is the division by zero exception; an infinite one is exact.
        if (is_inf() || rhs.is_zero()) {
            if (!is_inf())
                detail::raise_flags(FE_DIVBYZERO);
            *this = inf();
            set_sign(res_sign);
            return *this;
        } else if (is_zero() || rhs.is_inf()) {
            *this = 0;
            set_sign(res_sign);
            return *this;
        }

        // extract fractions and exponents
        uint32_t sign, rhs_sign;
        int32_t expo, rhs_expo;
        uint64_t l1, h1, l2, h2;
        get_components(l1, h1, expo, sign);
        rhs.get_components(l2, h2, rhs_expo, rhs_sign);

        // subtract the exponents
        expo -= rhs_expo;

        // optimize for rhs value is an exponent of 2
        if (rhs.is_exponent_of_2()) {
            set_components(l1, h1, expo, sign ^ rhs_sign);
            return *this;
        }

        // divide the fractions
        uint64_t q[4] {};
        const uint64_t nom[4] = {0, 0, l1, h1};
        const uint64_t denom[2] = {l2, h2};

        uint64_t rem[2] {};
        // get_components() normalizes even a subnormal, so the divisor's fraction always carries the
        // unity bit at bit 112 and its high QWORD is never zero. div_128bit's precondition therefore
        // holds unconditionally here, unlike in fixed_point128 where a tiny divisor has to be routed
        // to div_64bit instead.
        if (div_128bit(q, rem, nom, denom, array_length(nom))) {
            // The dividend was scaled by 2^128 and both mantissas are in [2^112, 2^113), so the
            // quotient is in (2^127, 2^129): its leading one is at bit 128 when the dividend's
            // mantissa is the larger and at bit 127 otherwise. Keeping 113 bits means shifting out
            // 16 or 15 of them.
            //
            // The shift has to follow the quotient. A fixed shift of 15, as this used to have,
            // left a quotient that led at bit 127 already in place, and the normalization after it
            // took that as nothing to round: every such quotient was truncated, which made close
            // to half of all divisions with the smaller mantissa on top come out one place low.
            const int32_t top = static_cast<int32_t>(q[2] & 1);
            const int32_t shift = (127 - FRAC_BITS) + top;
            // What is shifted out leads the extra word; a non zero remainder means the quotient
            // continues below it and only has to register as a sticky bit.
            const uint64_t extra = (q[0] << (64 - shift)) | (((rem[0] | rem[1]) != 0) ? 1 : 0);
            l1 = shift_right128(q[0], q[1], shift);
            h1 = shift_right128(q[1], q[2], shift);
            expo += top - 1;
            *this = round_pack(sign ^ rhs_sign, expo, l1, h1, extra);
        } else {  // error
            *this = inf();
        }
        return *this;
    }
    /**
     * @brief Divide this object by a value
     * @param rhs Right hand side operand
     * @return This object.
     */
    template <typename T> FP128_FORCE_INLINE float128& operator/=(const T& rhs) { return operator/=(float128(rhs)); }

    //
    // unary operations
    //
    /**
     * @brief Convert to bool
     */
    [[nodiscard]] FP128_FORCE_INLINE constexpr operator bool() const noexcept { return !is_zero(); }
    /**
     * @brief Logical not (!). Opposite of operator bool.
     * Uses is_zero() rather than testing the raw words: negative zero has its sign bit set and
     * would otherwise be reported as a non zero value.
     */
    [[nodiscard]] FP128_FORCE_INLINE constexpr bool operator!() const noexcept { return is_zero(); }
    /**
     * @brief Unary +. Returns a copy of the object.
     */
    [[nodiscard]] FP128_FORCE_INLINE constexpr float128 operator+() const noexcept { return *this; }
    /**
     * @brief Unary -. Returns a copy of the object with sign inverted.
     */
    [[nodiscard]] FP128_FORCE_INLINE constexpr float128 operator-() const noexcept { return float128(low, high ^ SIGN_MASK); }

    //
    // useful public functions
    //
    /**
     * @brief Returns true if the value is positive (including zero and NaN)
     * @return True when the sign is 0
     */
    [[nodiscard]] FP128_FORCE_INLINE constexpr bool is_positive() const noexcept { return 0 == (high & SIGN_MASK); }
    /**
     * @brief Returns true if the value is negative (including zero and NaN).
     * @return True when the sign is 1
     */
    [[nodiscard]] FP128_FORCE_INLINE constexpr bool is_negative() const noexcept { return 0 != (high & SIGN_MASK); }
    /**
     * @brief Returns true if and only if the value is ±0.
     * @return Returns true if the value is zero
     */
    [[nodiscard]] FP128_FORCE_INLINE constexpr bool is_zero() const noexcept { return 0 == low && 0 == (high & ~SIGN_MASK); }
    /**
     * @brief Returns true if and only if x is zero, subnormal or normal (not infinite or NaN).
     * @return True if and only if x is zero, subnormal or normal (not infinite or NaN).
     */
    [[nodiscard]] FP128_FORCE_INLINE constexpr bool is_finite() const noexcept { return !is_special(); }
    /**
     * @brief Tests if the value is subnormal
     * @return True when the value is subnormal
     */
    [[nodiscard]] FP128_FORCE_INLINE constexpr bool is_subnormal() const noexcept { return get_exponent_bits() == 0; }
    /**
     * @brief Tests if the value is normal (not zero, subnormal, infinite, or NaN)
     * @return True if and only if the value is normal
     */
    [[nodiscard]] FP128_FORCE_INLINE constexpr bool is_normal() const noexcept { return get_exponent_bits() != 0 && get_exponent_bits() != INF_EXP_BIASED; }
    /**
     * @brief Tests if this value is a NaN
     * @return True when the value is a NaN
     */
    [[nodiscard]] FP128_FORCE_INLINE constexpr bool is_nan() const
    {
        // fraction is zero for +- INF, non-zero for NaN
        return get_exponent_bits() == INF_EXP_BIASED && (get_fraction_bits() != 0 || low != 0);
    }
    /**
     * @brief Tests if this value is a signaling NaN
     * @return True if this value is a signaling NaN
     */
    [[nodiscard]] FP128_FORCE_INLINE constexpr bool is_signaling() const
    {
        // A NaN is quiet when the leading fraction bit is set and signaling when it is clear,
        // which is what every binary interchange format in IEEE 754-2008 uses.
        return is_nan() && (high & QUIET_NAN_BIT) == 0;
    }
    /**
     * @brief Tests if this value is an Infinite (negative or positive)
     * @return True when the value is an Infinite
     */
    [[nodiscard]] FP128_FORCE_INLINE constexpr bool is_inf() const
    {
        // fraction is zero for +- INF, non-zero for NaN
        return get_exponent_bits() == INF_EXP_BIASED && get_fraction_bits() == 0;
    }
    /**
     * @brief Tests if the value is an exponent of 2 (fraction part is zero)
     * @return True when the value is an exponent of 2
     */
    [[nodiscard]] FP128_FORCE_INLINE constexpr bool is_exponent_of_2() const
    {
        // fraction is zero for +- INF, non-zero for NaN
        return get_fraction_bits() == 0 && low == 0;
    }
    /**
     * @brief return true when the value is either an inf or nan
     * @return true for inf and nan
     */
    [[nodiscard]] FP128_FORCE_INLINE constexpr bool is_special() const
    {
        // fraction is zero for +- INF, non-zero for NaN
        return get_exponent_bits() == INF_EXP_BIASED;
    }
    /**
     * @brief Returns if the value is an integer (fraction is zero).
     * @return True when the value is an integer.
     */
    [[nodiscard]] FP128_FORCE_INLINE constexpr bool is_int() const
    {
        // An infinity or a NaN has an exponent beyond every integer's, but is not one.
        if (is_special())
            return false;
        int32_t expo = get_exponent();
        if (expo < 0)
            return false;
        if (expo >= FRAC_BITS)
            return true;
        return get_fraction().is_zero();
    }
    /**
     * @brief get a specific bit within the float128 data
     * @param bit bit to get [0,127]
     * @return 0 or 1. Undefined when bit > 127
     */
    [[nodiscard]] FP128_FORCE_INLINE constexpr int32_t get_bit(uint32_t bit) const noexcept
    {
        if (bit < 64) {
            return FP128_GET_BIT(low, bit);
        }
        return FP128_GET_BIT(high, bit - 64);
    }
    /**
     * @brief Return the fraction part as a float128
     */
    [[nodiscard]] FP128_INLINE constexpr float128 get_fraction() const
    {
        auto expo = get_exponent();
        int32_t frac_bits = static_cast<int32_t>(FRAC_BITS) - expo;
        // all the bits are fraction
        if (frac_bits > FRAC_BITS)
            return *this;
        // the exponent is too large to hold a fraction
        if (frac_bits <= 0)
            return 0;

        uint64_t l = low, h = get_fraction_bits();
        if (frac_bits <= 64) {
            h = 0;
            l &= FP128_MAX_VALUE_64(frac_bits);
        } else {
            h &= FP128_MAX_VALUE_64(frac_bits - 64);
        }

        // no fraction bits are high
        if (l == 0 && h == 0) {
            return 0;
        }

        // find the msb and shift to bit 112
        int32_t msb = static_cast<int32_t>(log2(l, h));
        int32_t shift = FRAC_BITS - msb;  // how many bits to shift left
        shift_left128_inplace_safe(l, h, shift);
        expo -= shift;
        float128 res(l, h, expo + EXP_BIAS, get_sign());
        return res;
    }
    /**
     * @brief Reads the biased exponent field out of the high QWORD.
     *
     * This accessor and the three below it are the only places that know where the fields sit
     * within the high QWORD, which is the same layout the (l, h, e, s) constructor assembles.
     * They deliberately shift and mask rather than going through the _float128_bits view: that
     * would be a read of the inactive member of a union, which constant evaluation rejects.
     *
     * @return The biased exponent, in [0, 0x7FFF].
     */
    [[nodiscard]] FP128_FORCE_INLINE constexpr uint64_t get_exponent_bits() const noexcept { return (high >> EXP_SHIFT) & EXP_MASK; }
    /**
     * @brief Sets the biased exponent field, leaving the fraction and the sign alone.
     * @param e Biased exponent. Bits above the 15 the field holds are dropped.
     */
    FP128_FORCE_INLINE constexpr void set_exponent_bits(uint64_t e) noexcept { high = (high & ~(EXP_MASK << EXP_SHIFT)) | ((e & EXP_MASK) << EXP_SHIFT); }
    /**
     * @brief Reads the upper 48 bits of the fraction out of the high QWORD.
     * @return Bits [111:64] of the fraction.
     */
    [[nodiscard]] FP128_FORCE_INLINE constexpr uint64_t get_fraction_bits() const noexcept { return high & UPPER_FRAC_MASK; }
    /**
     * @brief Sets the upper 48 bits of the fraction, leaving the exponent and the sign alone.
     * @param f Fraction bits. Bits above the 48 the field holds are dropped.
     */
    FP128_FORCE_INLINE constexpr void set_fraction_bits(uint64_t f) noexcept { high = (high & ~UPPER_FRAC_MASK) | (f & UPPER_FRAC_MASK); }
    /**
     * @brief Inverts the sign
     */
    FP128_FORCE_INLINE constexpr void invert_sign() noexcept { high ^= SIGN_MASK; }
    /**
     * @brief Sets the sign
     */
    FP128_FORCE_INLINE constexpr void set_sign(uint64_t s) noexcept { high = (high & ~SIGN_MASK) | ((s & 1) << 63); }
    /**
     * @brief Gets the sign
     */
    [[nodiscard]] FP128_FORCE_INLINE constexpr uint32_t get_sign() const noexcept { return static_cast<uint32_t>(high >> 63); }
    /**
     * @brief Reads the raw binary128 encoding.
     *
     * The counterpart of the float128(low, high) constructor: the two together are what a
     * std::bit_cast to and from the storage would be, without needing the layout to be spelled
     * out at the call site. Nothing is interpreted, so a NaN keeps its payload and a subnormal is
     * not normalized - use get_components() for the value the encoding stands for.
     *
     * @param l Receives the low QWORD (fraction bits 63:0)
     * @param h Receives the high QWORD (sign, exponent, fraction bits 111:64)
     */
    FP128_FORCE_INLINE constexpr void get_bits(uint64_t& l, uint64_t& h) const noexcept
    {
        l = low;
        h = high;
    }
    /**
     * @brief Classifies the float128 value according to IEEE 754.
     * @return The classification category (normal, subnormal, zero, inf, or NaN).
     */
    [[nodiscard]] FP128_INLINE constexpr float128_class_t get_class() const noexcept
    {
        // inf and Nan
        if (get_exponent_bits() == INF_EXP_BIASED) {
            if (is_nan())
                return is_signaling() ? signalingNaN : quietNaN;
            return is_positive() ? positiveInfinity : negativeInfinity;
        }
        if (is_zero()) {
            return is_positive() ? positiveZero : negativeZero;
        }
        if (is_subnormal()) {
            return is_positive() ? positiveSubnormal : negativeSubnormal;
        }

        return is_positive() ? positiveNormal : negativeNormal;
    }
    /**
     * @brief Returns the exponent of the object - like the base 2 exponent of a floating point
     * A value of 2.1 would return 1, values in the range [0.5,1.0) would return -1.
     * @return Exponent of the number
     */
    [[nodiscard]] FP128_FORCE_INLINE constexpr int32_t get_exponent() const noexcept { return static_cast<int32_t>(get_exponent_bits()) - EXP_BIAS; }
    /**
     * @brief Set the exponent
     * @param e Exponent value
     */
    FP128_FORCE_INLINE constexpr void set_exponent(int32_t e) noexcept
    {
        e += EXP_BIAS;
        assert(e >= 0);
        assert(e <= INF_EXP_BIASED);
        set_exponent_bits(static_cast<uint64_t>(e));
    }
    /**
     * @brief break the float into its components.
     * Normalizes subnormal values
     *
     * Forced open rather than merely offered for inlining, which is an exception to the rule stated
     * on FP128_FORCE_INLINE: this has a body of its own, so it would otherwise be FP128_INLINE. The
     * body is a handful of shifts and masks off the object's own bits, and every arithmetic operator
     * and math function opens with two of these calls, so what the caller gains is not the call
     * itself but the constant folding across it - the exponent and sign arithmetic that follows
     * collapses only once the components are visible. Clang declined the invitation where MSVC
     * accepted it, which is what made this measurable rather than academic: on 2026-08-15 clang-cl
     * left this and the rounding step after the division (then norm_fraction_sticky(), now
     * round_pack(), which is forced open for the same reason) out of line in operator/=, and float128
     * division measured 13.0 M/s against MSVC's 20.9 M/s.
     *
     * @param l Reference to receive the low fraction
     * @param h Reference to receive the high fraction
     * @param e Reference to receive the unbiased exponent
     * @param s Reference to receive the sign
     */
    FP128_FORCE_INLINE constexpr void get_components(uint64_t& l, uint64_t& h, int32_t& e, uint32_t& s) const noexcept
    {
        l = low;
        h = get_fraction_bits();
        e = get_exponent();
        s = get_sign();

        if (is_subnormal()) {
            // shift the bits to the left so the msb is on bit 112
            auto shift = FRAC_BITS - static_cast<int32_t>(log2(l, h));
            shift_left128_inplace_safe(l, h, shift);

            // A subnormal has no implicit leading one, so its stored exponent field says nothing
            // about its magnitude - every subnormal holds the same one. What the value is worth is
            // decided entirely by where its leading bit sits: the fraction is scaled by 2^-16494,
            // and normalizing it left by `shift` puts the exponent that many places below the
            // smallest normal one. Adding the shift to the stored exponent instead, as this used
            // to, put the smallest subnormal 223 binades above where it belongs and made every
            // operation on a subnormal operand return an unrelated value. set_components() has
            // always used this convention, which is why the two did not round trip.
            e = SUBNORM_EXP_UNBIASED + 1 - shift;
        }
        // normal numbers
        else if (is_normal()) {
            // add the unity value
            h |= FRAC_UNITY;
        }
    }
    /**
     * @brief Number of Mercator series terms log2() needs.
     *
     * The reduction leaves |z| <= 2^-6, so term n is bounded by 2^(-6n), and the series is cut off
     * once that is below the last bit of a 113 bit mantissa with eight bits to spare.
     */
    static constexpr int32_t LOG2_TERMS = (FRAC_BITS + 1 + 8 + log2_reduction_bits - 1) / log2_reduction_bits;

    /**
     * @name Rounding
     *
     * Every operation whose exact result may not be representable ends in round_pack(), which is
     * the one place a result is rounded. IEEE 754 requires each operation to behave as if it first
     * computed the exact result and then rounded it once, to the destination format, in the
     * current rounding direction - and the destination format includes its subnormal range, where
     * fewer than 113 bits are kept. The operations used to round to 113 bits first and then let
     * set_components() round a second time when the result turned out to be subnormal, which gave
     * the wrong last bit for about one subnormal result in ten. They also decided the direction
     * from a window of three bits, which cannot tell a tie from a value just above one when the
     * bits that settle it lie further down.
     *
     * The shape follows Berkeley SoftFloat's roundPackToF128: the caller hands over a 113 bit
     * significand and a 64 bit word of the bits below it, left aligned, with anything further down
     * jammed into its lowest bit. The top bit of that word is the first bit dropped, and the word
     * being non zero means the result is inexact - which is all any rounding direction needs.
     * @{
     */

    /**
     * @brief Shifts a significand and its extra word right, keeping a sticky record of what fell off.
     *
     * The bits shifted out of the significand move into the top of @p extra, and whatever drops
     * out of @p extra is jammed into its lowest bit rather than lost, so the word still tells an
     * exact result, a tie and a value above the tie apart.
     *
     * @param l Low QWORD of the significand
     * @param h High QWORD of the significand
     * @param extra Bits below the significand, left aligned
     * @param dist Bits to shift, at least one
     */
    FP128_FORCE_INLINE static constexpr void shift_right_jam128_extra(uint64_t& l, uint64_t& h, uint64_t& extra, int32_t dist) noexcept
    {
        FP128_ASSERT(dist >= 1);
        uint64_t sticky = extra;
        if (dist < 64) {
            extra = l << (64 - dist);
            l = (l >> dist) | (h << (64 - dist));
            h >>= dist;
        } else {
            if (dist == 64) {
                extra = l;
                l = h;
            } else {
                sticky |= l;
                if (dist < 128) {
                    extra = h << (128 - dist);
                    l = h >> (dist - 64);
                } else {
                    // Everything is gone. At exactly 128 the top bit of h is the first bit dropped,
                    // further out it is below it and only counts towards the sticky bit.
                    extra = (dist == 128) ? h : ((h != 0) ? 1 : 0);
                    l = 0;
                }
            }
            h = 0;
        }
        extra |= (sticky != 0) ? 1 : 0;
    }

    /**
     * @brief Shifts a 128 bit value right, jamming any set bit that falls off into the lowest bit.
     *
     * The alignment step of an addition: the shifted operand only has to be accurate to the bit
     * below the guard bit of the sum, and a jammed lowest bit is enough to say there was more.
     *
     * Forced open, as are shift_right_jam128_extra() and norm_round_pack(), which addition reaches
     * through. Offered for inlining only, Clang kept all three out of line, the three outputs went
     * through memory on every call, and float128 addition measured 95 M/s against 174 M/s inlined
     * (P-core, 2026-10-05). MSVC inlined them either way.
     *
     * @param l Low QWORD
     * @param h High QWORD
     * @param dist Bits to shift, at least one
     */
    FP128_FORCE_INLINE static constexpr void shift_right_jam128(uint64_t& l, uint64_t& h, int32_t dist) noexcept
    {
        FP128_ASSERT(dist >= 1);
        if (dist < 64) {
            const uint64_t lost = l << (64 - dist);
            l = (l >> dist) | (h << (64 - dist)) | ((lost != 0) ? 1 : 0);
            h >>= dist;
        } else if (dist < 128) {
            const uint64_t lost = (dist == 64) ? l : (l | (h << (128 - dist)));
            l = ((dist == 64) ? h : (h >> (dist - 64))) | ((lost != 0) ? 1 : 0);
            h = 0;
        } else {
            l = ((l | h) != 0) ? 1 : 0;
            h = 0;
        }
    }

    /**
     * @brief Whether rounding should add one unit in the last place.
     * @param mode Rounding direction
     * @param sign Sign of the result
     * @param extra Bits below the significand, left aligned and jammed
     * @return True when the magnitude has to be rounded up.
     */
    [[nodiscard]] FP128_FORCE_INLINE static constexpr bool round_increment(detail::rounding mode, uint32_t sign, uint64_t extra) noexcept
    {
        if (mode == detail::rounding::nearest_even)
            return extra >= (1ull << 63);
        if (mode == detail::rounding::toward_zero)
            return false;
        return extra != 0 && mode == ((sign != 0) ? detail::rounding::downward : detail::rounding::upward);
    }

    /// @brief The largest finite magnitude's significand, high QWORD, with its leading one at bit 48.
    static constexpr uint64_t MAX_SIG_HIGH = FRAC_UNITY | UPPER_FRAC_MASK;

    /**
     * @brief The slow half of round_pack(): results that overflow or are subnormal.
     * @param sign Sign of the result
     * @param exp Biased exponent less one, as round_pack() computes it
     * @param l Low QWORD of the significand
     * @param h High QWORD of the significand
     * @param extra Bits below the significand, left aligned and jammed
     * @param mode Rounding direction
     * @return The rounded result.
     */
    [[nodiscard]] FP128_INLINE static constexpr float128 round_pack_edge(uint32_t sign, int32_t exp, uint64_t l, uint64_t h, uint64_t extra,
                                                                         detail::rounding mode) noexcept
    {
        bool increment = round_increment(mode, sign, extra);
        if (exp < 0) {
            // Tininess is detected after rounding, the way x86 does it for double: the result is
            // tiny unless rounding it to 113 bits with an unbounded exponent would have carried it
            // up to the smallest normal value.
            const bool tiny = (exp < -1) || !increment || h != MAX_SIG_HIGH || l != UINT64_MAX;
            shift_right_jam128_extra(l, h, extra, -exp);
            exp = 0;
            if (tiny && extra != 0)
                detail::raise_flags(FE_UNDERFLOW);
            increment = round_increment(mode, sign, extra);
        } else if (exp > INF_EXP_BIASED - 2 || (exp == INF_EXP_BIASED - 2 && increment && h == MAX_SIG_HIGH && l == UINT64_MAX)) {
            detail::raise_flags(FE_OVERFLOW | FE_INEXACT);
            // Rounding to nearest and rounding away from zero overflow to infinity; the other two
            // directions stop at the largest finite value.
            if (mode == detail::rounding::nearest_even || mode == ((sign != 0) ? detail::rounding::downward : detail::rounding::upward))
                return float128(0, (static_cast<uint64_t>(sign) << 63) | (static_cast<uint64_t>(INF_EXP_BIASED) << EXP_SHIFT));
            return float128(UINT64_MAX, (static_cast<uint64_t>(sign) << 63) | (static_cast<uint64_t>(INF_EXP_BIASED - 1) << EXP_SHIFT) | UPPER_FRAC_MASK);
        }

        if (extra != 0)
            detail::raise_flags(FE_INEXACT);
        if (increment) {
            ++l;
            h += (l == 0) ? 1 : 0;
            if (mode == detail::rounding::nearest_even && (extra << 1) == 0)
                l &= ~1ull;  // an exact tie goes to the even neighbour
        } else if ((l | h) == 0) {
            exp = 0;
        }
        return float128(l, (static_cast<uint64_t>(sign) << 63) + (static_cast<uint64_t>(exp) << EXP_SHIFT) + h);
    }

    /**
     * @brief Rounds an exact result once, in the current rounding direction, and encodes it.
     *
     * The significand is added to the exponent field rather than merged into it: its leading one
     * sits on the lowest exponent bit, which is why the exponent is passed in less one. That makes
     * the carry of a significand that rounds up to the next power of two land in the exponent on
     * its own, and the carry out of the largest finite value land on infinity.
     *
     * @param sign Sign of the result
     * @param e Unbiased exponent of bit 112 of the significand
     * @param l Low QWORD of the significand
     * @param h High QWORD of the significand. The leading one must be at bit 48: a zero significand
     *        is not accepted, the callers deal with an exact zero themselves.
     * @param extra Bits below the significand, left aligned, with anything further down jammed into
     *        bit 0. Zero when the significand is exact.
     * @return The correctly rounded result, including its subnormal and overflow cases.
     */
    [[nodiscard]] FP128_FORCE_INLINE static constexpr float128 round_pack(uint32_t sign, int32_t e, uint64_t l, uint64_t h, uint64_t extra) noexcept
    {
        return round_pack(sign, e, l, h, extra, detail::current_rounding());
    }
    /// @overload
    /// @param mode Rounding direction to use instead of the current one
    [[nodiscard]] FP128_FORCE_INLINE static constexpr float128 round_pack(uint32_t sign, int32_t e, uint64_t l, uint64_t h, uint64_t extra,
                                                                          detail::rounding mode) noexcept
    {
        const int32_t exp = e + (EXP_BIAS - 1);
        // One unsigned comparison catches both a subnormal result (negative) and one that may
        // overflow (the top two exponents).
        if (static_cast<uint32_t>(exp) >= static_cast<uint32_t>(INF_EXP_BIASED - 2))
            return round_pack_edge(sign, exp, l, h, extra, mode);

        if (extra != 0)
            detail::raise_flags(FE_INEXACT);
        // Without a branch: whether to round up depends on the bits being dropped, which is as good
        // as random, and a mispredicted branch here costs more than the arithmetic it would skip.
        const uint64_t increment = round_increment(mode, sign, extra) ? 1 : 0;
        h += addcarryx_u64(0, l, increment, &l);
        // an exact tie that was rounded up goes back down to the even neighbour
        if (mode == detail::rounding::nearest_even)
            l &= ~(((extra << 1) == 0) ? increment : 0);
        return float128(l, (static_cast<uint64_t>(sign) << 63) + (static_cast<uint64_t>(exp) << EXP_SHIFT) + h);
    }

    /**
     * @brief Normalizes an exact non negative 128 bit significand and rounds it once.
     *
     * For results whose leading one can be anywhere: a sum after cancellation, a remainder, a
     * fraction. The value is (h:l) * 2^(e - 112).
     *
     * @param sign Sign of the result
     * @param e Exponent the value would have if its leading one were at bit 112
     * @param l Low QWORD
     * @param h High QWORD
     * @param sticky True when set bits below (h:l) were already dropped
     * @return The correctly rounded result.
     */
    [[nodiscard]] FP128_FORCE_INLINE static constexpr float128 norm_round_pack(uint32_t sign, int32_t e, uint64_t l, uint64_t h, bool sticky) noexcept
    {
        if ((l | h) == 0)
            return float128(0, static_cast<uint64_t>(sign) << 63);

        const int32_t msb = static_cast<int32_t>(log2(l, h));
        uint64_t extra = sticky ? 1 : 0;
        if (msb > FRAC_BITS) {
            shift_right_jam128_extra(l, h, extra, msb - FRAC_BITS);
        } else if (msb < FRAC_BITS) {
            shift_left128_inplace_safe(l, h, FRAC_BITS - msb);
        }
        return round_pack(sign, e + msb - FRAC_BITS, l, h, extra);
    }

    /**
     * @brief Builds a float128 from a 128 bit fraction, a value in [0,1).
     *
     * The counterpart of reading a mantissa out with get_components(): log2() does its argument
     * reduction and series on plain 128 bit fractions, and this is how the two results it ends up
     * with re-enter floating point.
     *
     * @param low Low QWORD of the fraction.
     * @param high High QWORD of the fraction.
     * @return The value the fraction stands for, correctly rounded, or zero if it is zero.
     */
    [[nodiscard]] static FP128_INLINE constexpr float128 from_fraction128(uint64_t low, uint64_t high) noexcept
    {
        // Bit 127 of the fraction is worth 2^-1, so a leading one at bit 112 would be worth 2^-16.
        return norm_round_pack(0, -16, low, high, false);
    }

    /**
     * @brief Stores a significand and exponent, rounding when the value is subnormal.
     *
     * @param l Low part of the significand
     * @param h High part of the significand, which must have its leading one at bit 48 (bit 112 of
     *        the whole) unless the significand is zero
     * @param e Unbiased exponent, can be any value: one past the format's range produces an
     *        infinity, one below it a subnormal or a zero, rounded once in the current direction.
     * @param s Sign (1 is negative)
     */
    FP128_INLINE constexpr void set_components(uint64_t l, uint64_t h, int32_t e, uint32_t s) noexcept
    {
        if ((l | h) == 0) {
            *this = float128(0, static_cast<uint64_t>(s != 0) << 63);
            return;
        }
        *this = round_pack((s != 0) ? 1u : 0u, e, l, h, 0);
    }
    /// @}

    /**
     * @name NaN results
     * @{
     */
    /// @brief x with its quiet bit set. The payload, and the sign, are kept.
    [[nodiscard]] FP128_FORCE_INLINE static constexpr float128 quiet_nan(const float128& x) noexcept { return float128(x.low, x.high | QUIET_NAN_BIT); }

    /**
     * @brief The result of an operation with a NaN operand (IEEE 754-2008 6.2).
     *
     * A signaling NaN operand is the invalid operation. Either way the result is a quiet NaN, and
     * the standard recommends it be one of the inputs so that a payload survives the computation:
     * the first NaN operand is returned, quieted.
     *
     * @param a First operand
     * @param b Second operand
     * @return The quiet NaN to deliver.
     */
    [[nodiscard]] FP128_INLINE static constexpr float128 propagate_nan(const float128& a, const float128& b) noexcept
    {
        if (a.is_signaling() || b.is_signaling())
            detail::raise_flags(FE_INVALID);
        return quiet_nan(a.is_nan() ? a : b);
    }
    /// @overload
    [[nodiscard]] FP128_INLINE static constexpr float128 propagate_nan(const float128& a, const float128& b, const float128& c) noexcept
    {
        if (a.is_signaling() || b.is_signaling() || c.is_signaling())
            detail::raise_flags(FE_INVALID);
        return quiet_nan(a.is_nan() ? a : (b.is_nan() ? b : c));
    }
    /// @overload
    [[nodiscard]] FP128_INLINE static constexpr float128 propagate_nan(const float128& a) noexcept
    {
        if (a.is_signaling())
            detail::raise_flags(FE_INVALID);
        return quiet_nan(a);
    }

    /**
     * @brief The result of an invalid operation on non NaN operands, such as inf - inf or 0 / 0.
     * @return The default quiet NaN, after raising the invalid flag.
     */
    [[nodiscard]] FP128_FORCE_INLINE static constexpr float128 invalid_operation() noexcept
    {
        detail::raise_flags(FE_INVALID);
        return nan();
    }

    /**
     * @brief The sign of an exact zero sum of two operands of opposite sign (IEEE 754-2008 6.3).
     * @return 1 (negative) when rounding downward, 0 in every other direction.
     */
    [[nodiscard]] FP128_FORCE_INLINE static constexpr uint32_t exact_zero_sign() noexcept
    {
        return (detail::current_rounding() == detail::rounding::downward) ? 1u : 0u;
    }
    /// @}
    /**
     * @brief Normalize the fraction so the msb (unity bit) is on bit 112.
     * The fraction value must contain the unity value
     * The exponent component is adjusted based on the bit shift direction required.
     * @param l Low part of the fraction
     * @param h High part of the fraction
     * @param e Unbiased exponent, can be any value.
     */
    /// @brief Word count of the fixed point expansion used to print the fraction digits.
    static constexpr int32_t FIXED_WORDS = 4;

    /**
     * @brief Expands the fraction part of this value into an exact fixed point number.
     *
     * The result is the fraction scaled by 2^256 and held in four words, frac[0] being the least
     * significant. Nothing is lost: a fraction bit of a binary128 has weight 2^(e-112) with e no
     * smaller than -111 on the path that calls this, so scaling by 2^256 always lands on an
     * integer, and the fraction is below one so the product stays under 2^256.
     *
     * Printing the digits from this expansion is exact. Deriving them by repeatedly multiplying a
     * float128 by 100000 instead, as this used to, rounds once per group and the error compounds
     * over the eight or so groups that a full precision value needs.
     *
     * @param frac Output, four words holding the fraction scaled by 2^256
     * @return True when a fraction exists, false when the value is an integer
     */
    [[nodiscard]] FP128_INLINE constexpr bool fraction_to_fixed(uint64_t frac[FIXED_WORDS]) const noexcept
    {
        for (int32_t i = 0; i < FIXED_WORDS; ++i)
            frac[i] = 0;

        uint64_t l, h;
        int32_t e;
        uint32_t s;
        get_components(l, h, e, s);

        // no bits sit below the binary point
        if (e >= FRAC_BITS)
            return false;

        // keep only the bits whose weight is below one. A value under one contributes all of them.
        if (e >= 0) {
            const int32_t keep = FRAC_BITS - e;  // 1 - 112
            if (keep <= 64) {
                l &= FP128_MAX_VALUE_64(keep);
                h = 0;
            } else {
                h &= FP128_MAX_VALUE_64(keep - 64);
            }
        }

        if (l == 0 && h == 0)
            return false;

        // bit i of the surviving fraction has weight 2^(e-112+i), which scaled by 2^256 lands on
        // bit e+144+i. The shift is non negative for every exponent that reaches this function.
        const int32_t shift = e + 2 * 64 + (FIXED_WORDS * 64 - 2 * 64) - FRAC_BITS;  // e + 144
        FP128_ASSERT(shift >= 0 && shift < FIXED_WORDS * 64);
        const int32_t word = shift >> 6;
        const int32_t bit = shift & 63;

        if (bit == 0) {
            frac[word] = l;
            if (word + 1 < FIXED_WORDS)
                frac[word + 1] = h;
        } else {
            frac[word] = l << bit;
            if (word + 1 < FIXED_WORDS)
                frac[word + 1] = (l >> (64 - bit)) | (h << bit);
            if (word + 2 < FIXED_WORDS)
                frac[word + 2] = h >> (64 - bit);
        }
        return true;
    }
    /**
     * @brief Multiplies the fixed point fraction by a small value and returns the part that
     *        crossed the binary point.
     * @param frac Fraction scaled by 2^256, updated in place to hold the new fraction
     * @param mul Multiplier, small enough that the result stays below 2^64
     * @return The integer part produced by the multiply, which is the next group of digits
     */
    [[nodiscard]] FP128_INLINE static constexpr uint64_t fixed_mul_extract(uint64_t frac[FIXED_WORDS], uint64_t mul) noexcept
    {
        uint64_t carry = 0;
        for (int32_t i = 0; i < FIXED_WORDS; ++i) {
            uint64_t hi;
            const uint64_t lo = mulx_u64(frac[i], mul, &hi);
            const uint8_t c = addcarryx_u64(0, lo, carry, &frac[i]);
            carry = hi + c;
        }
        // the fraction was below one, so the overflow is below mul and fits in a single word
        return carry;
    }
    /// @brief True when every word of the fixed point fraction is zero.
    [[nodiscard]] FP128_INLINE static constexpr bool fixed_is_zero(const uint64_t frac[FIXED_WORDS]) noexcept
    {
        uint64_t acc = 0;
        for (int32_t i = 0; i < FIXED_WORDS; ++i)
            acc |= frac[i];
        return acc == 0;
    }

    /**
     * @name Wide accumulator
     *
     * A fixed 384 bit unsigned integer, used by fma() to hold the exact sum of a 226 bit product
     * and a 113 bit addend before a single rounding is applied to it. The width is what the two
     * together need in the worst case: whichever of them is larger is placed with its leading bit
     * at the top of the accumulator, which leaves 383-225 = 158 bits below the product and
     * 383-112 = 271 bits below the addend. An operand that does not fit under the other one is
     * more than 158 bits smaller than it, so it cannot influence anything but the sticky bit, and
     * is folded into that instead of being placed.
     *
     * These are implementation details of fma() rather than part of the interface.
     * @{
     */
    static constexpr int32_t WIDE_WORDS = 6;               ///< Words in the accumulator.
    static constexpr int32_t WIDE_BITS = WIDE_WORDS * 64;  ///< Bits in the accumulator.
    /// @brief Bit the leading one of the larger operand is placed on.
    ///
    /// One below the top, so that adding two operands of the same magnitude cannot carry out of
    /// the accumulator.
    static constexpr int32_t WIDE_TOP = WIDE_BITS - 2;

    /**
     * @brief Full 256 bit product of two 128 bit values.
     * @param al Low QWORD of the first operand
     * @param ah High QWORD of the first operand
     * @param bl Low QWORD of the second operand
     * @param bh High QWORD of the second operand
     * @param res Output, four words with res[0] the least significant
     */
    FP128_INLINE static constexpr void mul128to256(uint64_t al, uint64_t ah, uint64_t bl, uint64_t bh, uint64_t res[4]) noexcept
    {
        uint64_t hi = 0;
        res[0] = mulx_u64(al, bl, &hi);
        res[1] = hi;
        res[2] = 0;
        res[3] = 0;

        uint64_t lo = mulx_u64(al, bh, &hi);
        unsigned char c = addcarryx_u64(0, res[1], lo, &res[1]);
        c = addcarryx_u64(c, res[2], hi, &res[2]);
        res[3] += c;

        lo = mulx_u64(ah, bl, &hi);
        c = addcarryx_u64(0, res[1], lo, &res[1]);
        c = addcarryx_u64(c, res[2], hi, &res[2]);
        res[3] += c;

        lo = mulx_u64(ah, bh, &hi);
        c = addcarryx_u64(0, res[2], lo, &res[2]);
        res[3] += hi + c;
    }

    /**
     * @brief Index of the highest set bit of a multi word value, or -1 when it is zero.
     * @param a Words, a[0] the least significant
     * @param words Word count
     */
    [[nodiscard]] FP128_INLINE static constexpr int32_t wide_msb(const uint64_t* a, int32_t words) noexcept
    {
        for (int32_t i = words - 1; i >= 0; --i) {
            if (a[i] != 0)
                return i * 64 + 63 - static_cast<int32_t>(lzcnt64(a[i]));
        }
        return -1;
    }

    /**
     * @brief Writes a value into a zeroed accumulator, shifted to a given bit position.
     *
     * A negative offset puts part or all of the value below the accumulator's least significant
     * bit. Those bits are not representable there but still decide the rounding, so whether any of
     * them was set is reported back rather than dropped.
     *
     * @param acc Accumulator, must be zero on entry
     * @param src Value to place, src[0] the least significant word
     * @param words Word count of src
     * @param offset Bit position the least significant bit of src lands on, may be negative
     * @return True when a set bit fell below the accumulator.
     */
    [[nodiscard]] FP128_INLINE static constexpr bool wide_place(uint64_t* acc, const uint64_t* src, int32_t words, int32_t offset) noexcept
    {
        if (offset >= 0) {
            const int32_t word = offset >> 6;
            const int32_t bit = offset & 63;
            for (int32_t i = 0; i < words; ++i) {
                const int32_t dst = i + word;
                if (dst >= WIDE_WORDS)
                    break;
                acc[dst] |= src[i] << bit;
                if (bit != 0 && dst + 1 < WIDE_WORDS)
                    acc[dst + 1] |= src[i] >> (64 - bit);
            }
            return false;
        }

        const int32_t drop = -offset;
        if (drop >= words * 64) {
            uint64_t any = 0;
            for (int32_t i = 0; i < words; ++i)
                any |= src[i];
            return any != 0;
        }

        const int32_t word = drop >> 6;
        const int32_t bit = drop & 63;
        uint64_t lost = 0;
        for (int32_t i = 0; i < word; ++i)
            lost |= src[i];
        if (bit != 0)
            lost |= src[word] & FP128_MAX_VALUE_64(bit);

        for (int32_t i = 0; word + i < words; ++i) {
            uint64_t value = src[word + i] >> bit;
            if (bit != 0 && word + i + 1 < words)
                value |= src[word + i + 1] << (64 - bit);
            if (i < WIDE_WORDS)
                acc[i] |= value;
        }
        return lost != 0;
    }

    /// @brief Adds b into a. The accumulator is wide enough that the sum cannot carry out.
    FP128_INLINE static constexpr void wide_add(uint64_t* a, const uint64_t* b) noexcept
    {
        unsigned char c = 0;
        for (int32_t i = 0; i < WIDE_WORDS; ++i)
            c = addcarryx_u64(c, a[i], b[i], &a[i]);
    }

    /// @brief Subtracts b from a, which must not be smaller.
    FP128_INLINE static constexpr void wide_sub(uint64_t* a, const uint64_t* b) noexcept
    {
        unsigned char borrow = 0;
        for (int32_t i = 0; i < WIDE_WORDS; ++i) {
            const uint64_t rhs = b[i];
            const uint64_t diff = a[i] - rhs - borrow;
            borrow = (a[i] < rhs || (borrow && a[i] == rhs)) ? 1 : 0;
            a[i] = diff;
        }
    }

    /// @brief Subtracts one from the accumulator, which must not be zero.
    FP128_INLINE static constexpr void wide_dec(uint64_t* a) noexcept
    {
        for (int32_t i = 0; i < WIDE_WORDS; ++i) {
            if (a[i]-- != 0)
                break;
        }
    }

    /// @brief Three way comparison of two accumulators.
    [[nodiscard]] FP128_INLINE static constexpr int32_t wide_cmp(const uint64_t* a, const uint64_t* b) noexcept
    {
        for (int32_t i = WIDE_WORDS - 1; i >= 0; --i) {
            if (a[i] != b[i])
                return (a[i] > b[i]) ? 1 : -1;
        }
        return 0;
    }
    /// @}
    /**
     * @brief Normalizes the product of two mantissas to bit 112, handing back the bits it drops.
     *
     * The multiply and the square both arrive here with (h:l) equal to their 256 bit product shifted
     * right by 111. Both operands were normalized to [2^112, 2^113), so the product is in [2^224, 2^226)
     * and (h:l) is in [2^113, 2^115): its leading one is at bit 113 or bit 114 and nowhere else, which
     * makes the remaining shift 1 or 2 and the choice a single bit test.
     *
     * Nothing is rounded here. The bits shifted out, and the sticky bit for everything the caller
     * dropped before them, are returned in the form round_pack() takes, so that the result is
     * rounded once and to its final width - which for a subnormal product is narrower than 113 bits.
     *
     * @param l Low part of the shifted product, replaced by the normalized fraction
     * @param h High part of the shifted product, replaced by the normalized fraction
     * @param e Unbiased exponent, adjusted to match the normalized fraction
     * @param sticky True when the caller already dropped one or more set bits below l
     * @return The extra word for round_pack(): the dropped bits left aligned, the sticky bit in bit 0.
     */
    [[nodiscard]] FP128_FORCE_INLINE static constexpr uint64_t norm_product(uint64_t& l, uint64_t& h, int32_t& e, bool sticky) noexcept
    {
        // bit 114 of the product decides between the two, and it is bit 50 of the high QWORD
        const int32_t shift = 1 + static_cast<int32_t>(h >> 50);
        e += shift;

        const uint64_t extra = (l << (64 - shift)) | (sticky ? 1 : 0);
        l = shift_right128(l, h, shift);
        h >>= shift;
        return extra;
    }
    /**
     * @brief Produces the closest value larger than x
     * nextUp(x) is the least floating-point number in the format of x that compares greater than x.
     * If x is the negative number of least magnitude in x’s format, nextUp(x) is −0.
     * nextUp(±0) is the positive number of least magnitude in x’s format.
     * nextUp(+∞) is +∞, and nextUp(−∞) is the finite negative number largest in magnitude.
     * When x is NaN, then the result is according to 6.2. nextUp(x) is quiet except for sNaNs.
     * @param x Source value
     * @return Higher value closest to x
     */
    [[nodiscard]] FP128_INLINE static constexpr float128 nextUp(float128 x)
    {
        switch (x.get_class()) {
        case quietNaN:
        case signalingNaN:
            return propagate_nan(x);
        case positiveInfinity:
            break;
        case negativeInfinity:
            return float128(UINT64_MAX, UINT64_MAX, EXP_MASK - 1, 1);
        case negativeZero:
        case positiveZero:
            x.low = 1;
            x.set_sign(0);
            break;
        case positiveSubnormal:
        case positiveNormal:
            ++x.low;
            if (x.low == 0) {
                ++x.high;  // takes care of raising the exponent if needed as well
            }
            break;
        case negativeSubnormal:
        case negativeNormal:
            --x.low;
            if (x.low == UINT64_MAX) {
                --x.high;  // takes care of decreasing the exponent if needed as well
            }
            break;
        }
        return x;
    }
    /**
     * @brief Produces the closest value smaller than x
     * @param x Source value
     * @return Lower value closest to x
     */
    [[nodiscard]] FP128_INLINE static constexpr float128 nextDown(float128 x)
    {
        switch (x.get_class()) {
        case quietNaN:
        case signalingNaN:
            return propagate_nan(x);
        case negativeInfinity:
            break;
        case positiveInfinity:
            return float128(UINT64_MAX, UINT64_MAX, EXP_MASK - 1, 0);
        case negativeZero:
        case positiveZero:
            x.low = 1;
            x.set_sign(1);
            break;
        case negativeSubnormal:
        case negativeNormal:
            ++x.low;
            if (x.low == 0) {
                ++x.high;  // takes care of raising the exponent if needed as well
            }
            break;
        case positiveSubnormal:
        case positiveNormal:
            --x.low;
            if (x.low == UINT64_MAX) {
                --x.high;  // takes care of decreasing the exponent if needed as well
            }
            break;
        }
        return x;
    }

    /**
     * @brief Return the infinite constant
     * @return INF
     */
    [[nodiscard]] FP128_FORCE_INLINE static constexpr float128 inf() { return float128(0, 0, INF_EXP_BIASED, 0); }
    /**
     * @brief Return the quiet (non-signaling) NaN constant
     * @return NaN
     */
    [[nodiscard]] FP128_FORCE_INLINE static constexpr float128 nan() { return float128(0, QUIET_NAN_BIT, INF_EXP_BIASED, 0); }
    /**
     * @brief Return a signaling NaN
     *
     * Nothing in the library raises an exception on encountering one - there is no floating point
     * status to raise it in - but the encoding is distinguishable, which is what
     * numeric_limits::signaling_NaN() and is_signaling() need.
     * @return A signaling NaN
     */
    [[nodiscard]] FP128_FORCE_INLINE static constexpr float128 signaling_nan() { return float128(1, 0, INF_EXP_BIASED, 0); }
    /**
     * @brief Return the value of pi
     * @return pi
     */
    [[nodiscard]] FP128_FORCE_INLINE static constexpr float128 pi() { return float128(0x8469898CC51701B8, 0x921FB54442D1, 0x4000, 0); }
    /**
     * @brief Return the value of pi / 2
     * @return pi / 2
     */
    [[nodiscard]] FP128_FORCE_INLINE static constexpr float128 half_pi() { return float128(0x8469898CC51701B8, 0x921FB54442D1, 0x3FFF, 0); }
    /**
     * @brief Return the value of pi / 4
     * @return pi / 4
     */
    [[nodiscard]] FP128_FORCE_INLINE static constexpr float128 quarter_pi() { return float128(0x8469898CC51701B8, 0x921FB54442D1, 0x3FFE, 0); }
    /**
     * @brief Return the value of e
     * @return e
     */
    [[nodiscard]] FP128_FORCE_INLINE static constexpr float128 e()
    {
        // static const float128 e = "2.71828182845904523536028747135266249775724709369"; // 50 first digits of e
        // return e;
        return float128(0x95355FB8AC404E7A, 0x5BF0A8B14576, 0x4000, 0);
    }
    /**
     * @brief Returns a value of sqrt(2)
     * @return
     */
    [[nodiscard]] FP128_FORCE_INLINE static constexpr float128 sqrt_2() noexcept
    {
        // static const float128 sqrt_2 = "1.41421356237309504880168872420969807856967187537"; // 50 first digits of sqrt(2)
        // return sqrt_2;
        return float128(0xC908B2FB1366EA95, 0x00006A09E667F3BC, 0x3FFF, 0);
    }
    /**
     * @brief  Returns a value of 1
     * @return 1
     */
    [[nodiscard]] FP128_FORCE_INLINE static constexpr float128 one() noexcept { return float128(0, 0, 0x3FFF, 0); }
    /**
     * @brief  Returns a value of 0.5
     * @return 0.5
     */
    [[nodiscard]] FP128_FORCE_INLINE static constexpr float128 half() noexcept { return float128(0, 0, 0x3FFE, 0); }
    /**
     * @brief Return 0.1 using maximum precision
     * @return
     */
    [[nodiscard]] FP128_FORCE_INLINE static constexpr float128 tenth() noexcept
    {
        // 0.1 using maximum precision
        return float128(0x999999999999999A, 0x999999999999, EXP_BIAS - 4, 0u);
    }

    /**
     * @name Mathematical constants
     *
     * The binary128 nearest each constant, given as the encoding rather than parsed from a decimal
     * string: the string constructor is accurate to about an ulp, which is one bit too few to
     * define a constant the rest of the library computes from. The set mirrors `<numbers>`, so a
     * generic function template written against std::numbers has the same names available here.
     * @{
     */
    /// @brief 2*pi
    [[nodiscard]] FP128_FORCE_INLINE static constexpr float128 two_pi() noexcept { return float128(0x8469898CC51701B8, 0x921FB54442D1, 0x4001, 0); }
    /// @brief 1/pi
    [[nodiscard]] FP128_FORCE_INLINE static constexpr float128 inv_pi() noexcept { return float128(0x2A53F84EAFA3EA6A, 0x45F306DC9C88, 0x3FFD, 0); }
    /// @brief 1/sqrt(pi)
    [[nodiscard]] FP128_FORCE_INLINE static constexpr float128 inv_sqrt_pi() noexcept { return float128(0xD11AE3A914FED7FE, 0x20DD750429B6, 0x3FFE, 0); }
    /// @brief Natural logarithm of 2
    [[nodiscard]] FP128_FORCE_INLINE static constexpr float128 ln2() noexcept { return float128(0xF35793C7673007E6, 0x62E42FEFA39E, 0x3FFE, 0); }
    /// @brief Natural logarithm of 10
    [[nodiscard]] FP128_FORCE_INLINE static constexpr float128 ln10() noexcept { return float128(0x582DD4ADAC5705A6, 0x26BB1BBB5551, 0x4000, 0); }
    /// @brief Base 2 logarithm of e
    [[nodiscard]] FP128_FORCE_INLINE static constexpr float128 log2_e() noexcept { return float128(0xE1777D0FFDA0D23A, 0x71547652B82F, 0x3FFF, 0); }
    /// @brief Base 10 logarithm of e
    [[nodiscard]] FP128_FORCE_INLINE static constexpr float128 log10_e() noexcept { return float128(0xE32A6AB7555F5A68, 0xBCB7B1526E50, 0x3FFD, 0); }
    /// @brief Base 10 logarithm of 2
    [[nodiscard]] FP128_FORCE_INLINE static constexpr float128 log10_2() noexcept { return float128(0xEF311F12B35816F9, 0x34413509F79F, 0x3FFD, 0); }
    /// @brief Square root of 3
    [[nodiscard]] FP128_FORCE_INLINE static constexpr float128 sqrt_3() noexcept { return float128(0xA73B25742D7078B8, 0xBB67AE8584CA, 0x3FFF, 0); }
    /// @brief 1/sqrt(3)
    [[nodiscard]] FP128_FORCE_INLINE static constexpr float128 inv_sqrt_3() noexcept { return float128(0xC4D218F81E4AFB25, 0x279A74590331, 0x3FFE, 0); }
    /// @brief The Euler-Mascheroni constant
    [[nodiscard]] FP128_FORCE_INLINE static constexpr float128 egamma() noexcept { return float128(0x8F49A37C7F0202A6, 0x2788CFC6FB61, 0x3FFE, 0); }
    /// @brief The golden ratio
    [[nodiscard]] FP128_FORCE_INLINE static constexpr float128 phi() noexcept { return float128(0x7C15F39CC0605CEE, 0x9E3779B97F4A, 0x3FFF, 0); }
    /// @}
    /**
     * @brief calculates 10^e
     * @param e integer exponent, in the range
     * @return 10^e
     */
    [[nodiscard]] FP128_INLINE static constexpr float128 exp10(int32_t e) noexcept
    {
        // check the limits first
        if (e < -4965) {
            return 0;
        } else if (e > 4932) {
            return inf();
        }

        // calculate the exponent optimally
        float128 res = 1;
        float128 b;
        if (e > 0)
            b = 10;  // 10^1
        else if (e < 0) {
            b = tenth();
            e = -e;
        }

        while (e > 0) {
            if (e & 1)
                res *= b;
            e >>= 1;
            b.square();
        }
        return res;
    }

    /// @brief Returns false; this type does not conform to IEEE 754-1985.
    static constexpr bool is754version1985(void) { return false; }
    /**
     * @brief Whether the type conforms to IEEE 754-2008 (binary128).
     *
     * The arithmetic, the conversions and the operations of clause 5 conform in every build, but
     * the standard also requires the rounding-direction attributes of clause 4 and the status flags
     * of clause 7, which exist only when FP128_IEEE_ENV is defined.
     *
     * @return True when the program is built with FP128_IEEE_ENV, false otherwise.
     */
    static constexpr bool is754version2008(void)
    {
#ifdef FP128_IEEE_ENV
        return true;
#else
        return false;
#endif
    }
    //
    // End of class method implementation
    //

    //
    // Binary math operators
    //
    // Each of these applies the compound assignment to the by value left hand side and then returns
    // it on a line of its own, rather than the shorter `return lhs OP= rhs;`. The two are
    // equivalent, but the short form makes the return statement copy construct from an lvalue, which
    // MSVC compiles into a store forwarding stall. int128_base's binary operators carry the full
    // explanation and the measurement.
    //
    /**
     * @brief Adds 2 values and returns the result.
     * @param lhs left hand side operand
     * @param rhs Right hand side operand
     * @return Result of the operation
     */
    template <typename T> [[nodiscard]] friend FP128_FORCE_INLINE constexpr float128 operator+(float128 lhs, const T& rhs) noexcept
    {
        lhs += rhs;
        return lhs;
    }
    /**
     * @brief subtracts the right hand side operand to this object to and returns the result.
     * @param lhs left hand side operand
     * @param rhs Right hand side operand
     * @return The float128 result
     */
    template <typename T> [[nodiscard]] friend FP128_FORCE_INLINE constexpr float128 operator-(float128 lhs, const T& rhs) noexcept
    {
        lhs -= rhs;
        return lhs;
    }
    /**
     * @brief Multiplies the right hand side operand with this object to and returns the result.
     * @param lhs left hand side operand
     * @param rhs Right hand side operand
     * @return The float128 result
     */
    template <typename T> [[nodiscard]] friend FP128_FORCE_INLINE constexpr float128 operator*(float128 lhs, const T& rhs) noexcept
    {
        lhs *= rhs;
        return lhs;
    }
    /**
     * @brief Divides this object by the right hand side operand and returns the result.
     * @param lhs left hand side operand
     * @param rhs Right hand side operand
     * @return The float128 result
     */
    template <typename T> [[nodiscard]] friend FP128_FORCE_INLINE float128 operator/(float128 lhs, const T& rhs)
    {
        lhs /= rhs;
        return lhs;
    }

    //
    // Binary math operators with the scalar on the left hand side
    //
    // Without these, an expression like (1 + x) is ambiguous: converting the literal to float128
    // and converting x to a builtin type are both one user defined conversion, so neither
    // overload wins. Restricting the left operand to the arithmetic types keeps these from
    // competing with the float128 on the left versions above, which would otherwise be an equally
    // good match. The comparison operators already carry the same pair of overloads.
    //

    /// @brief Adds a scalar and a float128, in that order. @param lhs Left operand @param rhs Right operand @return The float128 result
    template <typename T>
        requires std::is_arithmetic_v<T>
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr float128 operator+(const T& lhs, const float128& rhs) noexcept
    {
        return float128(lhs) += rhs;
    }
    /// @brief Subtracts a float128 from a scalar. @param lhs Left operand @param rhs Right operand @return The float128 result
    template <typename T>
        requires std::is_arithmetic_v<T>
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr float128 operator-(const T& lhs, const float128& rhs) noexcept
    {
        return float128(lhs) -= rhs;
    }
    /// @brief Multiplies a scalar and a float128, in that order. @param lhs Left operand @param rhs Right operand @return The float128 result
    template <typename T>
        requires std::is_arithmetic_v<T>
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr float128 operator*(const T& lhs, const float128& rhs) noexcept
    {
        return float128(lhs) *= rhs;
    }
    /// @brief Divides a scalar by a float128. @param lhs Left operand @param rhs Right operand @return The float128 result
    template <typename T>
        requires std::is_arithmetic_v<T>
    [[nodiscard]] friend FP128_FORCE_INLINE float128 operator/(const T& lhs, const float128& rhs)
    {
        return float128(lhs) /= rhs;
    }
    /**
     * @brief Shifts a scalar left by a float128 shift count, scaling it by a power of two.
     *
     * The left operand is widened to float128 first, so the result is 128 bit rather than the
     * builtin type of lhs.
     *
     * The count is rhs converted to int32_t, the same conversion the float128 on the left overload
     * applies, and that conversion truncates, so a count of 3.75 shifts by 3 - as it does in
     * fixed_point128, whose conversion is a plain shift.
     *
     * @param lhs Left operand, the value being shifted
     * @param rhs Right operand, the shift count
     * @return The float128 result
     */
    template <typename T>
        requires std::is_arithmetic_v<T>
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr float128 operator<<(const T& lhs, const float128& rhs) noexcept
    {
        return float128(lhs) <<= static_cast<int32_t>(rhs);
    }
    /**
     * @brief Shifts a scalar right by a float128 shift count, scaling it by a power of two.
     * The left operand is widened to float128 first, see operator<< above.
     * @param lhs Left operand, the value being shifted
     * @param rhs Right operand, the shift count
     * @return The float128 result
     */
    template <typename T>
        requires std::is_arithmetic_v<T>
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr float128 operator>>(const T& lhs, const float128& rhs) noexcept
    {
        return float128(lhs) >>= static_cast<int32_t>(rhs);
    }

    //
    // Comparison operators
    //

    /**
     * @brief Compare logical/bitwise equal.
     * @param lhs left hand side operand
     * @param rhs Right hand side operand
     * @return True if lhs and rhs are equal.
     */
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr bool operator==(const float128& lhs, const float128& rhs) noexcept
    {
        // A NaN compares equal to nothing, not even to another NaN with the same bits. Equality is
        // a quiet comparison: only a signaling NaN is the invalid operation.
        if (lhs.is_nan() || rhs.is_nan()) {
            if (lhs.is_signaling() || rhs.is_signaling())
                detail::raise_flags(FE_INVALID);
            return false;
        }
        // Positive and negative zero are numerically equal even though their bits differ.
        if (lhs.is_zero() && rhs.is_zero())
            return true;
        return lhs.high == rhs.high && lhs.low == rhs.low;
    }
    /// @overload
    // Constrained to a builtin type. Unconstrained, T deduces to float128 for a comparison of two
    // of them, which makes this template an exact match alongside the non-template above - and
    // under the C++20 rewriting of == and != the reversed forms join in, which MSVC reports as an
    // ambiguity from inside any standard container that compares keys.
    template <typename T>
        requires std::is_arithmetic_v<T>
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr bool operator==(const float128& lhs, const T& rhs) noexcept { return lhs == float128(rhs); }
    /// @overload
    // Constrained to a builtin type, see operator== above.
    template <typename T>
        requires std::is_arithmetic_v<T>
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr bool operator==(const T& lhs, const float128& rhs) noexcept
    {
        return rhs == float128(lhs);
    }
    /**
     * @brief Return true when objects are not equal. Can be used as logical XOR.
     * @param lhs left hand side operand
     * @param rhs Right hand side operand
     * @return True if not equal.
     */
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr bool operator!=(const float128& lhs, const float128& rhs) noexcept { return !(lhs == rhs); }
    /// @overload
    // Constrained to a builtin type, see operator== above.
    template <typename T>
        requires std::is_arithmetic_v<T>
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr bool operator!=(const float128& lhs, const T& rhs) noexcept { return lhs != float128(rhs); }
    /// @overload
    /// @copydoc operator!=
    template <typename T>
        requires std::is_arithmetic_v<T>
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr bool operator!=(const T& lhs, const float128& rhs) noexcept { return rhs != float128(lhs); }
    /**
     * @brief Return true if this object is small than the other
     * @param lhs left hand side operand
     * @param rhs Right hand side operand
     * @return True when this object is smaller.
     */
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr bool operator<(const float128& lhs, const float128& rhs) noexcept
    {
        // A NaN is unordered with everything, so every relational test involving one is false. The
        // relational operators are IEEE 754's signaling comparisons, so any NaN is the invalid
        // operation; isless() and the other <cmath> predicates are the quiet ones.
        if (lhs.is_nan() || rhs.is_nan()) {
            detail::raise_flags(FE_INVALID);
            return false;
        }
        // The two zeros are numerically equal, so neither is smaller than the other. Without this
        // the differing sign bits would make -0 compare smaller than +0.
        if (lhs.is_zero() && rhs.is_zero())
            return false;

        auto rhs_sign = rhs.get_sign();
        auto lhs_sign = lhs.get_sign();

        // signs are different
        if (lhs_sign != rhs_sign)
            return lhs_sign > rhs_sign;  // true when lhs_sign is 1 and rhs.sign is 0

        // MSB is the same, check the LSB, implies the exponent is identical
        if (lhs.high == rhs.high)
            return (lhs_sign) ? lhs.low > rhs.low : lhs.low < rhs.low;

        return (lhs_sign) ? lhs.high > rhs.high : lhs.high < rhs.high;
    }
    /// @overload
    template <typename T>
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr bool operator<(const float128& lhs, const T& rhs) noexcept { return lhs < float128(rhs); }
    /// @overload
    template <typename T>
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr bool operator<(const T& lhs, const float128& rhs) noexcept { return float128(lhs) < rhs; }
    /**
     * @brief Return true this object is small or equal than the other
     * @param lhs left hand side operand
     * @param rhs Right hand side operand
     * @return True when this object is smaller or equal.
     */
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr bool operator<=(const float128& lhs, const float128& rhs) noexcept
    {
        // Not simply !(lhs > rhs): a NaN makes every relational test false, so negating the
        // opposite test would wrongly report that a NaN is less than or equal to everything.
        if (lhs.is_nan() || rhs.is_nan()) {
            detail::raise_flags(FE_INVALID);
            return false;
        }
        return !(lhs > rhs);
    }
    /// @overload
    template <typename T>
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr bool operator<=(const float128& lhs, const T& rhs) noexcept { return lhs <= float128(rhs); }
    /// @overload
    template <typename T>
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr bool operator<=(const T& lhs, const float128& rhs) noexcept { return float128(lhs) <= rhs; }
    /**
     * @brief Return true this object is larger than the other
     * @param lhs left hand side operand
     * @param rhs Right hand side operand
     * @return True when this object is larger.
     */
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr bool operator>(const float128& lhs, const float128& rhs) noexcept
    {
        // A NaN is unordered with everything, so every relational test involving one is false. The
        // relational operators are IEEE 754's signaling comparisons, so any NaN is the invalid
        // operation; isless() and the other <cmath> predicates are the quiet ones.
        if (lhs.is_nan() || rhs.is_nan()) {
            detail::raise_flags(FE_INVALID);
            return false;
        }
        // the two zeros are numerically equal, so neither is larger than the other
        if (lhs.is_zero() && rhs.is_zero())
            return false;

        auto rhs_sign = rhs.get_sign();
        auto lhs_sign = lhs.get_sign();

        // signs are different
        if (lhs_sign != rhs_sign)
            return lhs_sign < rhs_sign;  // true when lhs_sign is 1 and rhs.sign is 0

        // MSB is the same, check the LSB, implies the exponent is identical
        if (lhs.high == rhs.high)
            return (lhs_sign) ? lhs.low < rhs.low : lhs.low > rhs.low;

        return (lhs_sign) ? lhs.high < rhs.high : lhs.high > rhs.high;
    }
    /// @overload
    template <typename T>
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr bool operator>(const float128& lhs, const T& rhs) noexcept { return lhs > float128(rhs); }
    /// @overload
    template <typename T>
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr bool operator>(const T& lhs, const float128& rhs) noexcept { return float128(lhs) > rhs; }
    /**
     * @brief Return true this object is larger or equal than the other
     * @param lhs left hand side operand
     * @param rhs Right hand side operand
     * @return True when this objext is larger or equal.
     */
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr bool operator>=(const float128& lhs, const float128& rhs) noexcept
    {
        // see the note on operator<= about why this is not simply !(lhs < rhs)
        if (lhs.is_nan() || rhs.is_nan()) {
            detail::raise_flags(FE_INVALID);
            return false;
        }
        return !(lhs < rhs);
    }
    /// @overload
    template <typename T>
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr bool operator>=(const float128& lhs, const T& rhs) noexcept { return lhs >= float128(rhs); }
    /// @overload
    template <typename T>
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr bool operator>=(const T& lhs, const float128& rhs) noexcept { return float128(lhs) >= rhs; }

    /**
     * @brief Return the NaN constant
     * @return A quiet NaN. The argument only exists so that the call can be found by argument
     *         dependent lookup, its value is ignored.
     */
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr float128 nan(const float128&) { return float128::nan(); }
    /**
     * @brief Tests if the value is a NaN
     * @param x Value to test
     * @return True when the value is a NaN
     */
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr bool isnan(const float128& x)
    {
        // zero for +- INF, non-zero for NaN
        return x.is_nan();
    }
    /**
     * @brief Tests if the value is an Infinite (negative or positive)
     * @param x Value to test
     * @return True when the value is an Infinite
     */
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr bool isinf(const float128& x) { return x.is_inf(); }

    /**
     * @brief Returns the absolute value of x.
     * @param x Input value
     * @return |x|
     */
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr float128 fabs(const float128& x) noexcept
    {
        float128 temp = x;
        temp.set_sign(0);
        return temp;
    }
    /**
     * @brief Rounds to an integral value (IEEE 754-2008 roundToIntegral, 5.3.1).
     *
     * The shared engine of floor, ceil, trunc, round, rint and nearbyint, which differ only in the
     * direction they round. It works on the encoding: the fraction bits below the units bit are
     * cleared, and when the value has to move up one unit is added at the units bit, a carry out of
     * the fraction landing in the exponent on its own. Nothing is subtracted, so nothing rounds,
     * and the sign is left alone - a negative value that rounds to zero gives -0.
     *
     * The previous implementations subtracted the fraction from the value, which made the zero
     * results of a negative argument positive, and round() added one half and truncated, which
     * rounded the largest value below one half up to one and 2^112 + 1 to 2^112 + 2.
     *
     * @param x Value to round
     * @param mode Direction, ignored when @p ties_away is set
     * @param ties_away Round to nearest with ties away from zero instead, as round() does
     * @param signal_inexact Raise the inexact exception when the result differs from x, which
     *        distinguishes roundToIntegralExact (rint) from the other five
     * @return The integral value. A NaN comes back quieted, an infinity unchanged.
     */
    [[nodiscard]] FP128_INLINE static constexpr float128 round_integral(const float128& x, detail::rounding mode, bool ties_away,
                                                                        bool signal_inexact) noexcept
    {
        if (x.is_nan())
            return propagate_nan(x);
        const int32_t expo = x.get_exponent();
        // Infinities, and every finite value from 2^112 up, have no fraction bits.
        if (expo >= FRAC_BITS || x.is_zero())
            return x;

        const uint32_t sign = x.get_sign();
        const uint64_t sign_bit = static_cast<uint64_t>(sign) << 63;
        const bool round_up_away = (mode == ((sign != 0) ? detail::rounding::downward : detail::rounding::upward));
        if (expo < 0) {
            // |x| < 1, subnormals included: the result is zero or one, with the sign of x.
            bool one = false;
            if (ties_away)
                one = expo == -1;
            else if (mode == detail::rounding::nearest_even)
                one = expo == -1 && (x.get_fraction_bits() | x.low) != 0;  // one half itself ties to zero
            else
                one = round_up_away;
            if (signal_inexact)
                detail::raise_flags(FE_INEXACT);
            return float128(0, sign_bit | (one ? (static_cast<uint64_t>(EXP_BIAS) << EXP_SHIFT) : 0));
        }

        // The units bit is bit 112 - expo of the encoding. For expo == 0 that is the exponent
        // field's lowest bit, which is set for the exponent of [1, 2): the implicit one is odd.
        const int32_t units = FRAC_BITS - expo;
        const uint64_t mask_low = (units >= 64) ? UINT64_MAX : ((1ull << units) - 1);
        const uint64_t mask_high = (units >= 64) ? ((1ull << (units - 64)) - 1) : 0;
        uint64_t l = x.low, h = x.high;
        if (((l & mask_low) | (h & mask_high)) == 0)
            return x;

        const int32_t half_index = units - 1;
        const uint64_t half_low = (half_index < 64) ? (1ull << half_index) : 0;
        const uint64_t half_high = (half_index < 64) ? 0 : (1ull << (half_index - 64));
        const bool at_or_above_half = ((l & half_low) | (h & half_high)) != 0;
        const bool rest = ((l & mask_low & ~half_low) | (h & mask_high & ~half_high)) != 0;
        const bool odd = ((units < 64) ? (l >> units) : (h >> (units - 64))) & 1;

        bool up = false;
        if (ties_away)
            up = at_or_above_half;
        else if (mode == detail::rounding::nearest_even)
            up = at_or_above_half && (rest || odd);
        else
            up = round_up_away;

        if (signal_inexact)
            detail::raise_flags(FE_INEXACT);
        l &= ~mask_low;
        h &= ~mask_high;
        if (up) {
            if (units >= 64) {
                h += 1ull << (units - 64);
            } else {
                const uint64_t unit = 1ull << units;
                l += unit;
                h += (l < unit) ? 1 : 0;
            }
        }
        return float128(l, h);
    }

    /**
     * @brief Converts an integral value to an integer type, the way the C rounding functions do.
     * @tparam I The integer type
     * @param x An integral value, an infinity or a NaN
     * @return The value, or zero when it is a NaN or out of range for the type, which is the
     *         invalid operation.
     */
    template <typename I> [[nodiscard]] FP128_INLINE static constexpr I integral_to(const float128& x) noexcept
    {
        if (x.is_nan() || x > std::numeric_limits<I>::max() || x < std::numeric_limits<I>::min()) {
            detail::raise_flags(FE_INVALID);
            return 0;
        }
        return x.to_integer<I>();
    }

    /**
     * @brief Performs the floor() function, similar to libc's floor(), rounds down towards -infinity.
     * @param x Input value
     * @return A float128 holding the integer value.
     */
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr float128 floor(const float128& x) noexcept
    {
        return round_integral(x, detail::rounding::downward, false, false);
    }
    /**
     * @brief Performs the ceil() function, similar to libc's ceil(), rounds up towards infinity.
     * @param x Input value
     * @return A float128 holding the integer value. A negative value above -1 gives -0.
     */
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr float128 ceil(const float128& x) noexcept
    {
        return round_integral(x, detail::rounding::upward, false, false);
    }
    /**
     * @brief Rounds towards zero
     * @param x Value to truncate
     * @return Integer value, rounded towards zero. A negative value above -1 gives -0.
     */
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr float128 trunc(const float128& x) noexcept
    {
        return round_integral(x, detail::rounding::toward_zero, false, false);
    }
    /**
     * @brief Rounds towards the nearest integer.
     * The halfway value (0.5) is rounded away from zero.
     * @param x Value to round
     * @return Integer value, rounded towards the nearest integer.
     */
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr float128 round(const float128& x) noexcept
    {
        return round_integral(x, detail::rounding::nearest_even, true, false);
    }
    /**
     * @brief Rounds x to an integer in the current rounding direction and returns it as int64_t.
     * @param x Input value
     * @return Nearest integer as int64_t. Returns 0 on overflow or for a NaN, raising invalid.
     */
    // rint rounds ties to even, round() rounds them away from zero, so these cannot forward to
    // llround/lround the way they used to: llrint(2.5) is 2 and llround(2.5) is 3.
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr int64_t llrint(const float128& x) noexcept { return integral_to<int64_t>(rint(x)); }
    /**
     * @brief Rounds towards the nearest integer.
     * The halfway value (0.5) is rounded away from zero.
     * @param x Value to round
     * @return Integer value, rounded towards the nearest integer. Returns 0 on overflow or for a
     *         NaN, raising invalid.
     */
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr int64_t llround(const float128& x) noexcept { return integral_to<int64_t>(round(x)); }
    /**
     * @brief Rounds x to an integer in the current rounding direction and returns it as int32_t.
     * @param x Input value
     * @return Nearest integer as int32_t. Returns 0 on overflow or for a NaN, raising invalid.
     */
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr int32_t lrint(const float128& x) noexcept { return integral_to<int32_t>(rint(x)); }
    /**
     * @brief Rounds towards the nearest integer.
     * The halfway value (0.5) is rounded away from zero.
     * @param x Value to round
     * @return Integer value, rounded towards the nearest integer. Returns 0 on overflow or for a
     *         NaN, raising invalid.
     */
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr int32_t lround(const float128& x) noexcept { return integral_to<int32_t>(round(x)); }

    /**
     * @brief Retrieves an integer that represents the base-2 exponent of the specified value.
     * @param x The specified value.
     * @return Integer value, rounded towards the nearest integer.
     */
    [[nodiscard]] friend FP128_INLINE constexpr int32_t ilogb(const float128& x) noexcept
    {
        // The three answers <cmath> reserves for the arguments that have no exponent. A subnormal
        // does have one, but not the one its exponent field holds, so it goes through
        // get_components() with the rest.
        if (x.is_zero())
            return FP_ILOGB0;
        if (x.is_nan())
            return FP_ILOGBNAN;
        if (x.is_inf())
            return INT_MAX;

        uint64_t l = 0, h = 0;
        int32_t expo = 0;
        uint32_t sign = 0;
        x.get_components(l, h, expo, sign);
        return expo;
    }
    /**
     * @brief returns the value of x with the sign of y.
     * @param x The value that's returned as the magnitude of the result.
     * @param y The sign of the result.
     * @return The copysign functions return a floating-point value that combines the magnitude of x and the sign of y.
     */
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr float128 copysign(const float128& x, const float128& y) noexcept
    {
        float128 temp = x;
        temp.set_sign(y.get_sign());
        return temp;
    }
    /**
     * @brief Performs the fmod() function, similar to libc's fmod(), returns the remainder of a division x/root.
     * @param x Numerator
     * @param y Denominator
     * @return The modulo value.
     */
    [[nodiscard]] friend float128 fmod(const float128& x, const float128& y) noexcept
    {
        // a NaN operand propagates
        if (x.is_nan() || y.is_nan())
            return propagate_nan(x, y);

        // fmod(x, 0) is the invalid operation and produces a NaN, matching the CRT.
        // An infinite dividend is invalid for the same reason.
        if (y.is_zero() || x.is_inf())
            return invalid_operation();

        // trivial case, x is zero
        if (x.is_zero())
            return x;

        // an infinite divisor leaves the dividend untouched
        if (y.is_inf())
            return x;

        // fmod is exact: the result is x - n*y for an integer n, and it is representable. Reaching
        // it through trunc(x/y), as this used to, is only exact while the quotient is: once x/y
        // needs more than 113 bits the truncated quotient is a rounded value and the subtraction
        // returns something unrelated to the remainder. The mantissas are integers, so the whole
        // operation is a long division on them, done here one 64 bit block of the exponent
        // difference at a time.
        uint64_t lx = 0, hx = 0, ly = 0, hy = 0;
        int32_t ex = 0, ey = 0;
        uint32_t sx = 0, sy = 0;
        x.get_components(lx, hx, ex, sx);
        y.get_components(ly, hy, ey, sy);

        // |x| < |y| leaves the dividend as the remainder
        if (ex < ey || (ex == ey && uint128_t(lx, hx) < uint128_t(ly, hy)))
            return x;

        // Both mantissas hold 113 bits with the leading one at bit 112, so the values compared
        // here are x = rem * 2^(ex-112) and y = mod * 2^(ey-112).
        uint128_t rem(lx, hx);
        const uint128_t mod(ly, hy);
        rem %= mod;
        int32_t shift = ex - ey;

        // A remainder is below the divisor and so holds at most 113 bits; shifting 15 more in
        // keeps the dividend inside the 128 the modulo can take.
        constexpr int32_t block = 15;
        while (shift > 0 && !rem.is_zero()) {
            const int32_t step = (shift < block) ? shift : block;
            rem <<= step;
            shift -= step;
            rem %= mod;
        }

        if (rem.is_zero()) {
            float128 zero;
            zero.set_sign(sx);
            return zero;
        }

        // Rebuild the value: the remainder is an integer scaled by the divisor's last place.
        uint64_t rl = 0, rh = 0;
        rem.get_components(rl, rh);
        const int32_t msb = static_cast<int32_t>(log2(rl, rh));
        shift_left128_inplace_safe(rl, rh, FRAC_BITS - msb);

        float128 res;
        res.set_components(rl, rh, ey - FRAC_BITS + msb, sx);
        return res;
    }
    /**
     * @brief Split into integer and fraction parts.
     * Both results carry the sign of the input variable.
     * @param x Input value
     * @param iptr Pointer to float128 holding the integer part of x.
     * @return The fraction part of x. Undefined when iptr is nullptr.
     */
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr float128 modf(const float128& x, float128* iptr) noexcept
    {
        if (iptr == nullptr)
            return 0;

        // A NaN splits into two NaNs, an infinity into itself and a zero fraction.
        if (x.is_nan()) {
            *iptr = propagate_nan(x);
            return *iptr;
        }
        *iptr = trunc(x);
        if (x.is_inf())
            return float128(0, static_cast<uint64_t>(x.get_sign()) << 63);

        // The difference is exact, and a zero one takes the sign of x as C requires.
        float128 res = x - *iptr;
        res.set_sign(x.get_sign());
        return res;
    }
    /**
     * @brief Determines the positive difference between the first and second values.
     * @param x First value
     * @param y Second value
     * @return If x > y returns x - y. Otherwise zero.
     */
    [[nodiscard]] friend FP128_INLINE constexpr float128 fdim(const float128& x, const float128& y) noexcept
    {
        if (x.is_nan() || y.is_nan())
            return propagate_nan(x, y);
        return (x > y) ? x - y : float128();
    }
    /**
     * @brief Returns the minimum between x and y.
     * @param x First value
     * @param y Second value
     * @return If x < y returns x. Otherwise y.
     */
    [[nodiscard]] friend FP128_INLINE constexpr float128 fmin(const float128& x, const float128& y) noexcept
    {
        // A quiet NaN operand is treated as missing rather than as a value, so the other one wins.
        // Every comparison against a NaN is false, which made the plain conditional return
        // whichever operand happened to sit on the false branch. A signaling NaN is the invalid
        // operation instead, and produces a quiet NaN (IEEE 754-2008 minNum).
        if (x.is_signaling() || y.is_signaling())
            return propagate_nan(x, y);
        if (x.is_nan())
            return y;
        if (y.is_nan())
            return x;
        // The comparison operators treat the two zeros as equal, and the standard asks for the
        // negative one here.
        if (x.is_zero() && y.is_zero())
            return x.is_negative() ? x : y;
        return (x < y) ? x : y;
    }
    /**
     * @brief Returns the maximum between x and y.
     * @param x First value
     * @param y Second value
     * @return If x > y returns x. Otherwise y.
     */
    [[nodiscard]] friend FP128_INLINE constexpr float128 fmax(const float128& x, const float128& y) noexcept
    {
        if (x.is_signaling() || y.is_signaling())
            return propagate_nan(x, y);
        if (x.is_nan())
            return y;
        if (y.is_nan())
            return x;
        if (x.is_zero() && y.is_zero())
            return x.is_positive() ? x : y;
        return (x > y) ? x : y;
    }
    /**
     * @brief Calculates the hypotenuse. i.e. sqrt(x^2 + y^2)
     * @param x First value
     * @param y Second value
     * @return sqrt(x^2 + y^2).
     */
    [[nodiscard]] friend FP128_INLINE float128 hypot(const float128& x, const float128& y) noexcept
    {
        // An infinity wins even over a quiet NaN (IEEE 754-2008 9.2.1); a signaling NaN does not.
        if (x.is_signaling() || y.is_signaling())
            return propagate_nan(x, y);
        if (x.is_inf() || y.is_inf())
            return inf();
        if (x.is_nan() || y.is_nan())
            return propagate_nan(x, y);

        // Scaling by the larger side keeps both squares inside the format's range. Squaring the
        // values as they came, as this used to, overflowed for a side above 2^8192 and lost
        // everything below 2^-8247 to underflow.
        const float128 largest = fmax(fabs(x), fabs(y));
        if (largest.is_zero())
            return largest;
        // From here on only the flags the result justifies reach the caller.
        detail::flag_filter filter;
        const int32_t expo = ilogb(largest);
        const float128 a = ldexp(x, -expo);
        const float128 b = ldexp(y, -expo);
        return filter(ldexp(sqrt(sqr(a) + sqr(b)), expo));
    }
    /**
     * @brief Calculates the square of a value. i.e. x^2
     *
     * Faster than x * x and identical to it, see float128::square().
     *
     * @param x Value to square
     * @return x^2, which is never negative.
     */
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr float128 sqr(float128 x) noexcept { return x.square(); }
    /**
     * @brief Square of a 113 bit significand compared with a 256 bit value.
     * @param rl Low QWORD of the significand
     * @param rh High QWORD of the significand
     * @param n The value to compare against, n[0] the least significant word
     * @return Negative, zero or positive as r*r is below, equal to or above n.
     */
    [[nodiscard]] FP128_INLINE static constexpr int32_t square_compare(uint64_t rl, uint64_t rh, const uint64_t n[4]) noexcept
    {
        uint64_t sq[4] {};
        mul128to256(rl, rh, rl, rh, sq);
        for (int32_t i = 3; i >= 0; --i) {
            if (sq[i] != n[i])
                return (sq[i] > n[i]) ? 1 : -1;
        }
        return 0;
    }

    /**
     * @brief Calculates the square root, correctly rounded.
     *
     * IEEE 754 requires the square root to be correctly rounded, the same as the four arithmetic
     * operations, and Newton's method alone cannot promise that: every step rounds, and the last one
     * can land a unit either side of the answer. It used to, for a quarter of all arguments.
     *
     * So Newton's method only gets close, and the last bit is then settled exactly. The mantissa m
     * (doubled when the exponent is odd, so that what is left of it halves exactly) is in
     * [2^112, 2^114), which puts the square root of N = m * 2^112 in [2^112, 2^113) - exactly the
     * 113 bit significand the result needs. Squaring a candidate r in integer arithmetic and
     * comparing with N moves it onto floor(sqrt(N)), and the remainder N - r*r then says which way
     * to round: the root is past the halfway point r + 1/2 exactly when the remainder exceeds r,
     * and it can never sit on it, because (r + 1/2)^2 is not an integer.
     *
     * Based on the book "Math toolkit for real time programming" by Jack W. Crenshaw.
     *
     * @param x Value to calculate the root of
     * @param iterations Ignored. The result is correctly rounded whatever its value; the parameter
     *        is kept so that existing calls still compile.
     * @return Square root of x. sqrt(-0) is -0, and any other negative argument is the invalid
     *         operation and returns a NaN.
     */
    [[nodiscard]] friend float128 sqrt(const float128& x, [[maybe_unused]] uint32_t iterations) noexcept
    {
        if (x.is_nan())
            return propagate_nan(x);
        if (x.is_zero())
            return x;  // sqrt(-0) is -0
        if (x.is_negative())
            return invalid_operation();
        if (x.is_inf())
            return x;

        // Split off an even power of two, leaving a mantissa in [1, 4). Halving an even exponent is
        // exact. Reading the mantissa through get_components() also normalizes a subnormal
        // argument.
        uint64_t l = 0, h = 0;
        int32_t expo = 0;
        uint32_t sign = 0;
        x.get_components(l, h, expo, sign);
        const int32_t odd = expo & 1;
        const float128 norm_x(l, h, static_cast<uint32_t>(EXP_BIAS + odd), 0);
        const int32_t half_expo = (expo - odd) / 2;

        // The hardware double gives 53 correct bits to start from and every Newton step doubles
        // them, so two passes reach the last bit give or take the rounding of the steps themselves.
        //                  X
        //   Xn+1 = 0.5 * (---- + Xn )
        //                  Xn
        float128 root;
        {
            // The steps are inexact even when the root is not; only the final rounding may say so.
            const detail::flag_quiet quiet;
            root = ::sqrt(static_cast<double>(norm_x));
            for (int32_t i = 0; i < 2; ++i)
                root = (norm_x / root + root) >> 1;
        }

        // The candidate as a 113 bit integer. The root of a value below 4 is below 2, so an
        // exponent of one only arises from rounding up to 2 itself, and the largest significand
        // stands in for it.
        uint64_t rl = 0, rh = 0;
        int32_t root_expo = 0;
        uint32_t root_sign = 0;
        root.get_components(rl, rh, root_expo, root_sign);
        if (root_expo > 0) {
            rl = UINT64_MAX;
            rh = MAX_SIG_HIGH;
        }

        // N = m * 2^112, with m the mantissa shifted up by one for an odd exponent.
        if (odd != 0)
            shift_left128_inplace_safe(l, h, 1);
        const uint64_t n[4] = {0, l << 48, (l >> 16) | (h << 48), h >> 16};

        // Walk the candidate onto floor(sqrt(N)). It starts within a unit or two, so neither loop
        // runs more than twice.
        while (square_compare(rl, rh, n) > 0) {
            if (rl-- == 0)
                --rh;
        }
        for (;;) {
            uint64_t nl = rl + 1, nh = rh + ((rl == UINT64_MAX) ? 1 : 0);
            if (square_compare(nl, nh, n) > 0)
                break;
            rl = nl;
            rh = nh;
        }

        // The remainder N - r*r is below 2r + 1, so it fits in the low 128 bits.
        uint64_t sq[4] {};
        mul128to256(rl, rh, rl, rh, sq);
        uint64_t rem_l = 0, rem_h = 0;
        const uint8_t borrow = subborrow_u64(0, n[0], sq[0], &rem_l);
        subborrow_u64(borrow, n[1], sq[1], &rem_h);

        const bool exact = (rem_l | rem_h) == 0;
        const bool above_half = (rem_h > rh) || (rem_h == rh && rem_l > rl);
        const uint64_t extra = above_half ? ((1ull << 63) | 1) : (exact ? 0 : 1);
        return round_pack(0, half_expo, rl, rh, extra);
    }
    /**
     * @brief Calculates the cube root.
     * Uses the Halley's method.
     * @param x Floating point value
     * @param iterations how many Halley to perform, usually 1 is enough
     * @return cube root of x
     */
    [[nodiscard]] friend FP128_INLINE float128 cbrt(const float128 x, uint32_t iterations) noexcept
    {
        if (x.is_nan())
            return propagate_nan(x);
        // The real cube root is defined for a negative argument and is an odd function, which is
        // what the C library's cbrt computes. Returning a NaN, as this used to, made cbrt(-8)
        // undefined where the standard has it equal to -2.
        if (x.is_zero() || x.is_inf())
            return x;
        // From here on only the flags the result justifies reach the caller.
        detail::flag_filter filter;

        // Split off a power of two whose exponent is a multiple of three, leaving a mantissa in
        // [1, 8). Dividing that exponent by three is exact, so no correction factor is needed.
        //
        // The factors this used to multiply by, cbrt(2)/2 and cbrt(4)/2, were written out to the
        // 17 digits a double round trips through and no further, so two thirds of all arguments
        // came back with 53 correct bits instead of 113.
        uint64_t l = 0, h = 0;
        int32_t expo = 0;
        uint32_t sign = 0;
        x.get_components(l, h, expo, sign);
        int32_t rem = expo % 3;
        if (rem < 0)
            rem += 3;
        const float128 norm_x(l, h, static_cast<uint32_t>(EXP_BIAS + rem), 0);
        const int32_t third_expo = (expo - rem) / 3;

        // Halley's method triples the correct digits per pass, so two passes take the 53 bits the
        // hardware seed provides past the 113 the mantissa holds.
        float128 root = ::cbrt(static_cast<double>(norm_x));

        //                3
        //              Xn  + 2X
        //   Xn+1 = Xn ----------
        //                3
        //              2Xn  + X
        const auto x2 = norm_x << 1;
        for (auto i = iterations; i != 0; --i) {
            float128 r_cube = root * root * root;
            root = root * (r_cube + x2) / ((r_cube << 1) + norm_x);
        }

        root = ldexp(root, third_expo);
        root.set_sign(sign);
        return filter(root);
    }
    /**
     * @brief Calculates the reciprocal of a value. y = 1 / x
     * Using newton iterations: Yn+1 = Yn(2 - x * Yn)
     * @param x Input value
     * @return 1 / x. Returns zero on overflow or division by zero
     */
    [[nodiscard]] friend FP128_INLINE float128 reciprocal(const float128& x) noexcept
    {
        static const float128 two = 2;
        constexpr int max_iterations = 3;
        constexpr int debug = false;
        const auto x_sign = x.get_sign();
        if (x.is_nan())
            return propagate_nan(x);
        if (x.is_inf()) {
            float128 zero;
            zero.set_sign(x_sign);
            return zero;
        }
        if (x.is_zero()) {
            detail::raise_flags(FE_DIVBYZERO);
            return (x_sign) ? -inf() : inf();
        }

        // get_components() normalizes a subnormal, so the mantissa below is always in [1, 2) and
        // the Newton iteration sees the same problem whatever the magnitude of x was. Rebuilding
        // the value from the mantissa alone is also what keeps the scaling exact.
        uint64_t l = 0, h = 0;
        int32_t expo = 0;
        uint32_t sign = 0;
        x.get_components(l, h, expo, sign);
        const float128 norm_x(l, h, EXP_BIAS, 0);

        float128 y = 1.0 / static_cast<double>(norm_x);

        if (!y)
            return y;

        float128 xy, y_prev;
        // Newton iterations:
        int i = 0;
        for (; i < max_iterations && (y_prev != y); ++i) {
            y_prev = y;
            xy = norm_x * y;
            // y = y * (two - xy);
            y *= two - xy;
        }

        if constexpr (debug) {
            static int debug_max_iter = 0;
            if (i > debug_max_iter || i == max_iterations) {
                debug_max_iter = i;
                printf("reciprocal took %i iterations for %.10lf\n", i, static_cast<double>(x));
            }
        }

        // The mantissa was divided out above, so undoing the exponent is a scaling rather than a
        // replacement. Overwriting it, as this used to, was right only while y kept an exponent of
        // -1, which is every value except the one case y comes out exactly 1: reciprocal() of any
        // power of two was off by a factor of two, and exp() of a small negative argument with it.
        y = ldexp(y, -expo);
        y.set_sign(x_sign);
        return y;
    }
    /**
     * @brief Factorial reciprocal (inverse). Calculates 1 / x!
     * Maximum supported value of x is 50.
     * @param x Input value
     * @param res Result of the function
     */
    friend FP128_INLINE void fact_reciprocal(int x, float128& res) noexcept
    {
        // Given as encodings rather than decimal strings. Several of the strings this replaced
        // were the 17 digit expansion of a double - 1/6! among them - which capped the accuracy of
        // every series built on the table at 53 bits. The exponential's own reduction is good to
        // better than 2^-110, and a single term formed from a double sized constant was throwing
        // 40 of those bits away.
        static constexpr float128 c[] = {
            float128(0x0000000000000000, 0x000000000000, 0x3FFF, 0),  // 1 / 0!
            float128(0x0000000000000000, 0x000000000000, 0x3FFF, 0),  // 1 / 1!
            float128(0x0000000000000000, 0x000000000000, 0x3FFE, 0),  // 1 / 2!
            float128(0x5555555555555555, 0x555555555555, 0x3FFC, 0),  // 1 / 3!
            float128(0x5555555555555555, 0x555555555555, 0x3FFA, 0),  // 1 / 4!
            float128(0x1111111111111111, 0x111111111111, 0x3FF8, 0),  // 1 / 5!
            float128(0x6C16C16C16C16C17, 0x6C16C16C16C1, 0x3FF5, 0),  // 1 / 6!
            float128(0xA01A01A01A01A01A, 0xA01A01A01A01, 0x3FF2, 0),  // 1 / 7!
            float128(0xA01A01A01A01A01A, 0xA01A01A01A01, 0x3FEF, 0),  // 1 / 8!
            float128(0x38FAAC1C88E50017, 0x71DE3A556C73, 0x3FEC, 0),  // 1 / 9!
            float128(0xC72EF016D3EA6679, 0x27E4FB7789F5, 0x3FE9, 0),  // 1 / 10!
            float128(0x38FE747E4B837DC7, 0xAE64567F544E, 0x3FE5, 0),  // 1 / 11!
            float128(0x7B544DA987ACFE85, 0x1EED8EFF8D89, 0x3FE2, 0),  // 1 / 12!
            float128(0x97CA38331D23AF68, 0x6124613A86D0, 0x3FDE, 0),  // 1 / 13!
            float128(0xD20BADF145DFA3E5, 0x93974A8C07C9, 0x3FDA, 0),  // 1 / 14!
            float128(0xF11D8656B0EE8CB0, 0xAE7F3E733B81, 0x3FD6, 0),  // 1 / 15!
            float128(0xF11D8656B0EE8CB0, 0xAE7F3E733B81, 0x3FD2, 0),  // 1 / 16!
            float128(0xA6B2605197771B00, 0x952C77030AD4, 0x3FCE, 0),  // 1 / 17!
            float128(0x77BB004886A2C2AB, 0x6827863B97D9, 0x3FCA, 0),  // 1 / 18!
            float128(0x724CA1EC3B7B9675, 0x2F49B4681415, 0x3FC6, 0),  // 1 / 19!
            float128(0x507A9CAD2BF8F0BB, 0xE542BA402022, 0x3FC1, 0),  // 1 / 20!
            float128(0x18BEF146FCEE6E45, 0x71B8EF6DCF57, 0x3FBD, 0),  // 1 / 21!
            float128(0x29450C90B7F338EC, 0x0CE396DB7F85, 0x3FB9, 0),  // 1 / 22!
            float128(0x9D97B8704DD7F628, 0x761B41316381, 0x3FB4, 0),  // 1 / 23!
            float128(0x7CCA4B4067CA9D8A, 0xF2CF01972F57, 0x3FAF, 0),  // 1 / 24!
            float128(0x8D4E44A419776F11, 0x3F3CCDD165FA, 0x3FAB, 0),  // 1 / 25!
            float128(0x9A38F2050BA6B015, 0x88E85FC6A4E5, 0x3FA6, 0),  // 1 / 26!
            float128(0x320A9A18F15D4277, 0xD1AB1C2DCCEA, 0x3FA1, 0),  // 1 / 27!
            float128(0xD373C5C51C354A8D, 0x0A18A2635085, 0x3F9D, 0),  // 1 / 28!
            float128(0xD7ABE30E7766F129, 0x259F98B4358A, 0x3F98, 0),  // 1 / 29!
            float128(0xE60CADED4C2989C5, 0x3932C5047D60, 0x3F93, 0),  // 1 / 30!
            float128(0xC42E1EE46FA6BFC4, 0x434D2E783F5B, 0x3F8E, 0),  // 1 / 31!
            float128(0xC42E1EE46FA6BFC4, 0x434D2E783F5B, 0x3F89, 0),  // 1 / 32!
            float128(0x1B5382CDFFA97422, 0x3981254DD0D5, 0x3F84, 0),  // 1 / 33!
            float128(0xA13F8A2B4AF9D6B7, 0x2710231C0FD7, 0x3F7F, 0),  // 1 / 34!
            float128(0xF2833C7F5A7E0624, 0x0DC59C716D91, 0x3F7A, 0),  // 1 / 35!
            float128(0x92B06B8D12A7275C, 0xDF983290C2CA, 0x3F74, 0),  // 1 / 36!
            float128(0xAF4C78B15C3D89D3, 0x9EC8D1C94E85, 0x3F6F, 0),  // 1 / 37!
            float128(0xAE913D370A4EC4E8, 0x5D4ACB9C0C3A, 0x3F6A, 0),  // 1 / 38!
            float128(0xDE0104476AEB4C3B, 0x1E99449A4BAC, 0x3F65, 0),  // 1 / 39!
            float128(0x3001A07244ABAD2B, 0xCA8ED42A12AE, 0x3F5F, 0),  // 1 / 40!
            float128(0x0C7E25CFD1B1B2DD, 0x65E61C39D024, 0x3F5A, 0),  // 1 / 41!
            float128(0x836C4D9225DCB90A, 0x10AF527530DE, 0x3F55, 0),  // 1 / 42!
            float128(0x22DCBAE56DEF3720, 0x95DB45257E51, 0x3F4F, 0),  // 1 / 43!
            float128(0xA4FD9F327E7F6DE9, 0x272B1B03FEC6, 0x3F4A, 0),  // 1 / 44!
            float128(0x78E02C5E91C64CAC, 0xA3CB87222064, 0x3F44, 0),  // 1 / 45!
            float128(0x062CA46E4F25C608, 0x240804F65951, 0x3F3F, 0),  // 1 / 46!
            float128(0x9B78B454D8B628E5, 0x8DA8E0A127EB, 0x3F39, 0),  // 1 / 47!
            float128(0x67A5CD8DE5CEC5EE, 0x091B406B6FF2, 0x3F34, 0),  // 1 / 48!
            float128(0x5D94A3FD410E1232, 0x5A42F0DFEB08, 0x3F2E, 0),  // 1 / 49!
            float128(0x8205F0A053453602, 0xBB36F6E12CD7, 0x3F28, 0),  // 1 / 50!
        };
        constexpr int series_len = array_length(c);
        static_assert(series_len == 51);

        res = (x >= 0 && x < series_len) ? c[x] : float128();
    }
    /**
     * @brief returns the double factorial of a number
     * @param x The number to compute the double factorial of. Values below 100 come from a
     *          precomputed table.
     * @return x!!, or infinity for the values above the table.
     */
    [[nodiscard]] friend FP128_INLINE float128 double_factorial(int x) noexcept
    {
        constexpr int32_t arr_size = 100;
        static float128 c[arr_size];
        if (c[0].is_zero()) {
            c[0] = 1;
            c[1] = 1;

            // compute the odd and even double factorials
            for (int i = 2; i < arr_size - 1; i += 2) {
                c[i] = c[i - 2] * i;
                c[i + 1] = c[i - 1] * (i + 1);
            }
        }
        if (x < arr_size)
            return c[x];

        // TODO: compute the following members
        return float128::inf();
    }
    /**
     *                                       x
     * @brief Calculates the exponent of x: e
     * Using the Maclaurin series expansion, the formula is:
     *                  1       2       3
     *                 x       x       x
     * exp(x) = 1  +  ---  +  ---  +  --- + ...
     *                 1!      2!      3!
     *
     * The Maclaurin series will quickly overflow as x's power increases rapidly.
     *                     x   ix   fx
     * Using the equality e = e  * e
     * Where ix is the integer part of x and fx is the fraction part.
     * ix is computed via multiplications which won't overflow if the result value can be held.
     * fx is computed via Maclaurin series expansion, but since fx < 1, it won't overflow.
     * @param x A number specifying a power.
     * @return Exponent of x
     */
    /// @brief Terms of the Maclaurin series exp_reduced() runs. |r| <= ln2/2, so r^28/28! is past the last bit.
    static constexpr int32_t EXP_TERMS = 28;

    /**
     * @brief e to the power of a reduced argument, |r| <= ln2/2.
     * @param r Reduced argument
     * @return exp(r)
     */
    [[nodiscard]] FP128_INLINE static float128 exp_reduced(const float128& r) noexcept
    {
        float128 acc, factorial;
        fact_reciprocal(EXP_TERMS, acc);
        for (int32_t n = EXP_TERMS - 1; n >= 0; --n) {
            fact_reciprocal(n, factorial);
            acc = factorial + r * acc;
        }
        return acc;
    }

    /**
     * @brief Computes e to the power of x.
     *
     * x is written as k*ln2 + r with k an integer and |r| <= ln2/2, which turns the answer into
     * 2^k * exp(r): the power of two is an exact scaling and only the reduced argument goes
     * through a series.
     *
     * The reduction is where the accuracy of the whole function is decided. ln2 is irrational, so
     * k*ln2 cannot be represented and subtracting a rounded one leaves an error of k ulp - at the
     * top of the range k reaches 16000 and fourteen bits of the answer are gone before the series
     * even starts. Splitting ln2 into a high part whose low 14 bits are zero, so that k times it
     * is exact, and a low part carrying the rest, keeps r accurate to far below its last bit.
     *
     * The previous implementation raised e to the integer part by repeated squaring, which doubles
     * its own error at every step, and reached the negative half of the range through reciprocal().
     *
     * @param x A number specifying a power.
     * @return e^x
     */
    [[nodiscard]] friend FP128_INLINE float128 exp(const float128& x) noexcept
    {
        if (x.is_nan())
            return propagate_nan(x);
        if (x.is_inf())
            return x.is_negative() ? float128() : x;
        if (x.is_zero())
            return float128::one();
        // From here on only the flags the result justifies reach the caller.
        detail::flag_filter filter;

        // log of the largest finite value is 11356.52, log of the smallest subnormal is -11433.46
        if (x > float128(11357))
            return filter.inexact(inf());
        if (x < float128(-11434))
            return filter.inexact(float128());

        // ln2 to 321 bits, in three pieces. The first two have their low 16 mantissa bits zeroed,
        // so multiplying either by a k of up to 2^16 is exact and the products can be removed from
        // the argument without a rounding of their own.
        //
        // Two pieces would only carry ln2 to the 113 bits a float128 holds, which is not enough:
        // the reduction multiplies whatever error the constant has by k, and at the top of the
        // range k reaches 16000, so a constant good to 2^-114 leaves the reduced argument good to
        // only 2^-100 and the answer 8000 ulp wide.
        constexpr float128 ln2_hi(0xF35793C767300000, 0x62E42FEFA39E, 0x3FFE, 0);
        constexpr float128 ln2_mid(0x93394C5B16C50000, 0xF97B57A079A1, 0x3F98, 0);
        constexpr float128 ln2_lo(0x7CF70EC40DBD7593, 0xA2EB71755F45, 0x3F32, 0);

        const int32_t k = static_cast<int32_t>(llround(x * float128::log2_e()));
        const float128 kf(k);
        // The first subtraction is exact, cancelling two values within a factor of two of each
        // other; the others remove quantities far below the last place of the result.
        const float128 r = ((x - kf * ln2_hi) - kf * ln2_mid) - kf * ln2_lo;

        return filter(ldexp(exp_reduced(r), k));
    }
    /**
     * @brief Computes 2 to the power of x
     * @param x Exponent value
     * @return 2^x
     */
    [[nodiscard]] friend FP128_INLINE float128 exp2(const float128& x) noexcept
    {
        //
        // Based on exponent law: (x^n)^m = x^(m*n)
        // Convert the exponent x (function parameter) to produce an exponent that will work with exp()
        // y = log(2)
        // 2^x = e^(y*x) = exp(y*x)
        //
        if (x.is_nan())
            return propagate_nan(x);
        if (x.is_inf())
            return x.is_negative() ? float128() : x;
        if (x.is_zero())
            return float128::one();
        // From here on only the flags the result justifies reach the caller.
        detail::flag_filter filter;
        if (x > float128(16384))
            return filter.inexact(inf());
        if (x < float128(-16495))
            return filter.inexact(float128());

        // The integer part scales exactly, so only the fraction reaches exp() and the argument it
        // is handed stays small enough that the rounding of the multiply below is the only error.
        // That rounding is then recovered with an exact fma and folded back in: exp(p + e) is
        // exp(p) * (1 + e) to well past the last bit for an e this small.
        const int32_t k = static_cast<int32_t>(llround(x));
        const float128 f = x - float128(k);
        const float128 p = f * float128::ln2();
        const float128 correction = fma(f, float128::ln2(), -p);

        return filter(ldexp(exp(p) * (float128::one() + correction), k));
    }
    /**
     * @brief Calculates the exponent of x and reduces 1 from the result: (e^x) - 1
     * @param x A number specifying a power.
     * @return Exponent of x
     */
    [[nodiscard]] friend FP128_INLINE float128 expm1(const float128& x) noexcept
    {
        if (x.is_nan())
            return propagate_nan(x);
        if (x.is_zero())
            return x;
        if (x.is_inf())
            return x.is_negative() ? -float128::one() : x;
        // From here on only the flags the result justifies reach the caller.
        detail::flag_filter filter;

        // Half is the point where exp(x) - 1 stops throwing away significant bits: below it the
        // subtraction cancels most of what exp() produced, and at 2^-60 there is nothing left of
        // the answer at all. The Maclaurin series has no cancellation - every term is formed from
        // the argument itself - and needs 30 terms to reach the last bit at this magnitude.
        if (fabs(x) <= float128::half()) {
            constexpr int32_t terms = 30;
            float128 acc, factorial;
            fact_reciprocal(terms + 1, acc);
            for (int32_t n = terms - 1; n >= 0; --n) {
                fact_reciprocal(n + 1, factorial);
                acc = factorial + x * acc;
            }
            return filter(x * acc);
        }

        return filter(exp(x) - float128::one());
    }
    /**
     * @brief a^n by repeated squaring, for a positive finite a.
     *
     * A negative exponent inverts the power, unless the power overflowed or underflowed on its way
     * there while its reciprocal would not have: then exp() and log() reach it directly.
     *
     * @param a Base, positive and finite
     * @param n Exponent
     * @return a^n.
     */
    [[nodiscard]] FP128_INLINE static float128 pown_magnitude(const float128& a, int32_t n) noexcept
    {
        // The magnitude is taken in the unsigned domain: negating the most negative int32_t is
        // undefined.
        uint32_t expo = (n < 0) ? (0u - static_cast<uint32_t>(n)) : static_cast<uint32_t>(n);
        float128 res = float128::one();
        float128 b = a;
        while (expo > 0) {
            if (expo & 1)
                res *= b;
            expo >>= 1;
            if (expo > 0)
                b.square();
        }
        if (n >= 0)
            return res;
        if (res.is_inf() || res.is_zero() || res.is_subnormal())
            return exp(float128(n) * log(a));
        return float128::one() / res;
    }
    /**
     * @brief Computes x to the power of y (IEEE 754-2008 pown)
     * @param x Base value
     * @param y Exponent value (integer)
     * @return x^y
     */
    [[nodiscard]] friend FP128_INLINE float128 pow(const float128& x, int32_t y) noexcept
    {
        // Every int32_t is exact in a float128, so the general function answers this too, with
        // IEEE 754's pown special cases being the same as its pow ones for an integer exponent.
        return pow(x, float128(y));
    }
    /**
     * @brief Computes x to the power of y
     *
     * The special values follow IEEE 754-2008 9.2.1 throughout. Among those the previous version
     * got wrong: pow(0.5, inf) and pow(2, -inf) were infinite rather than zero, pow(2, NaN) was
     * infinite, pow(-0, 0.5) a NaN, and every integer exponent beyond 11355 in magnitude overflowed
     * whatever the base, so pow(1 + 2^-40, 12000), which is about 1.00000001, came out infinite.
     *
     * @param x Base value
     * @param y Exponent value
     * @return x^y
     */
    [[nodiscard]] friend FP128_INLINE float128 pow(const float128& x, const float128& y) noexcept
    {
        const float128 one_value = float128::one();

        // IEEE 754-2008 9.2.1, in the order its exceptions to the NaN rule need: x^0 is one and
        // so is 1^y, for any x and y, quiet NaNs included.
        if (y.is_zero() || x == one_value) {
            if (x.is_signaling() || y.is_signaling())
                return propagate_nan(x, y);
            return one_value;
        }
        if (x.is_nan() || y.is_nan())
            return propagate_nan(x, y);

        // The sign of a zero or infinite x survives only an odd integer exponent.
        const bool y_odd = is_odd_int(y);
        const bool negative_odd = x.is_negative() && y_odd;
        if (x.is_zero()) {
            if (y.is_negative()) {
                // a pole, except that 0^-inf is an exact infinity
                if (!y.is_inf())
                    detail::raise_flags(FE_DIVBYZERO);
                return negative_odd ? -inf() : inf();
            }
            return negative_odd ? x : float128();
        }
        if (y.is_inf()) {
            const float128 a = fabs(x);
            if (a == one_value)
                return one_value;  // (-1)^+-inf
            return ((a > one_value) != y.is_negative()) ? inf() : float128();
        }
        if (x.is_inf()) {
            if (y.is_negative())
                return negative_odd ? -float128() : float128();
            return negative_odd ? -inf() : inf();
        }

        // Both finite and non zero. A negative base has a real power only for an integer exponent.
        if (x.is_negative() && !y.is_int())
            return invalid_operation();
        // From here on only the flags the result justifies reach the caller.
        detail::flag_filter filter;

        const float128 a = fabs(x);
        float128 res;
        // Repeated multiplication is exact for as long as the powers are, and otherwise loses
        // about half an ulp per unit of the exponent; exp(y * log(a)) loses about half an ulp per
        // unit of y * log(a) instead. The multiplication wins for a short exponent, and for any
        // base outside [1/2, 2) where log(a) is at least about one. It is limited to exponents an
        // int32_t holds, which is more than any finite result needs at that size.
        const int32_t a_expo = ilogb(a);
        if (y.is_int() && fabs(y) <= float128(INT32_MAX) && (fabs(y) <= float128(64) || a_expo >= 1 || a_expo <= -2))
            res = pown_magnitude(a, static_cast<int32_t>(y));
        else
            res = exp(y * log(a));

        return filter(negative_odd ? -res : res);
    }
    /**
     * @brief The arguments every logarithm treats alike (IEEE 754-2008 9.2.1).
     *
     * A NaN propagates, a zero is a pole - the division by zero exception, and -inf - any other
     * negative argument is the invalid operation, and +inf is its own logarithm. log2() used to
     * return -inf for a negative argument, a finite value for +inf, and a finite value for a NaN.
     *
     * @param x Argument
     * @param result Receives the answer when there is one
     * @return True when x was one of these and @p result holds the answer.
     */
    [[nodiscard]] FP128_INLINE static constexpr bool log_special(const float128& x, float128& result) noexcept
    {
        if (x.is_nan()) {
            result = propagate_nan(x);
            return true;
        }
        if (x.is_zero()) {
            detail::raise_flags(FE_DIVBYZERO);
            result = -inf();
            return true;
        }
        if (x.is_negative()) {
            result = invalid_operation();
            return true;
        }
        if (x.is_inf()) {
            result = x;
            return true;
        }
        return false;
    }

    /// @brief Largest |t| the log1p_small() series is used for. Chosen so twelve terms suffice.
    [[nodiscard]] FP128_FORCE_INLINE static constexpr float128 log1p_small_limit() noexcept { return float128(0, 0, EXP_BIAS - 4, 0); }
    /// @brief Terms of the log1p_small() series. |s| stays below 2^-4.9, so s^24 is past the last bit.
    static constexpr int32_t LOG1P_TERMS = 12;

    /**
     * @brief Natural logarithm of 1+t for |t| <= 2^-4, accurate relative to the result.
     *
     * log(1+t) = 2 * atanh(t / (2+t)). The substitution replaces the Mercator series in t, which
     * would need a hundred terms at this precision, with one in s^2 where |s| < 2^-4.9, so twelve
     * terms reach the last bit of the mantissa.
     *
     * What matters more than the term count is that the result stays accurate relative to its own
     * size. Near one the logarithm is proportional to t, and t reaches this function exactly - the
     * subtraction that produced it cancelled two values within a factor of two of each other, which
     * is exact. Computing the same thing as a difference of two numbers near one loses a bit for
     * every power of two the argument sits closer to one: at t = 2^-60 barely half the mantissa
     * survives.
     *
     * @param t Offset from one, must satisfy |t| <= 2^-4
     * @return log(1+t)
     */
    [[nodiscard]] FP128_INLINE static float128 log1p_small(const float128& t) noexcept
    {
        if (t.is_zero())
            return t;

        const float128 s = t / (float128(2) + t);
        const float128 w = sqr(s);

        // 1/(2k+1) for the Horner pass below. Built once: the divisions are far more expensive
        // than the series itself, and the values do not depend on the argument.
        static const auto odd_reciprocal = [] {
            std::array<float128, LOG1P_TERMS + 1> table {};
            for (int32_t k = 0; k <= LOG1P_TERMS; ++k)
                table[static_cast<size_t>(k)] = float128::one() / float128(2 * k + 1);
            return table;
        }();

        float128 acc = odd_reciprocal[LOG1P_TERMS];
        for (int32_t k = LOG1P_TERMS - 1; k >= 0; --k)
            acc = odd_reciprocal[static_cast<size_t>(k)] + w * acc;

        return (s * acc) << 1;
    }
    /**
     * @brief Calculates the natural Log (base e) of x: log(x)
     * @param x The number to perform log on.
     * @return log(x)
     */
    [[nodiscard]] friend FP128_INLINE float128 log(float128 x) noexcept
    {
        float128 special;
        if (log_special(x, special))
            return special;
        // From here on only the flags the result justifies reach the caller.
        detail::flag_filter filter;

        // Sterbenz guarantees the subtraction is exact over the range the branch covers, so the
        // series below sees the offset from one with every bit it has.
        const float128 t = x - float128::one();
        if (fabs(t) <= log1p_small_limit())
            return filter(log1p_small(t));

        return filter(log2(x) * float128::ln2());
    }
    /**
     * @brief Calculates the Log base 2 of x: y = log2(x)
     *
     * Accurate relative to its own result, which matters most for an argument close to one: the
     * answer is then proportional to x - 1, and that difference is carried through exactly rather
     * than being rounded onto a grid it would barely register on. That holds on both sides of one;
     * below it the exponent is -1, and the answer is formed without the -1 that would cancel. The
     * earlier implementation built the answer as a fixed point fraction of 112 bits before
     * returning it as a float, which left log2(1 + 2^-112) with no correct significant bits at all.
     *
     * @param x The number to perform log2 on.
     * @return log2(x)
     */
    [[nodiscard]] friend FP128_INLINE float128 log2(float128 x) noexcept
    {
        float128 special;
        if (log_special(x, special))
            return special;
        // From here on only the flags the result justifies reach the caller.
        detail::flag_filter filter;

        // Calculate the log in 2 steps:
        // - The integer part is simple and fast via the exponent.
        // - The fraction part is log2 of the mantissa, by argument reduction and a short series.
        // The result is the sum of the two. Based on the identity:
        // log(x + y) = log(x) + log(y)
        //
        // The exponent comes from get_components(), which normalizes a subnormal: the exponent
        // field of every subnormal is the same and says nothing about its magnitude.
        uint64_t frac_low = 0, frac_high = 0;
        int32_t expo = 0;
        uint32_t mantissa_sign = 0;
        x.get_components(frac_low, frac_high, expo, mantissa_sign);
        frac_high &= UPPER_FRAC_MASK;  // drop the unity bit, leaving f in [0,1) scaled by 2^112

        // x is a power of 2, so the mantissa contributes nothing.
        if ((frac_low | frac_high) == 0) {
            return filter(float128(expo));
        }

        // The leading log2_reduction_bits fraction bits choose the reciprocal to reduce with.
        const size_t j = static_cast<size_t>(frac_high >> (FRAC_BITS - 64 - log2_reduction_bits));

        // x in [1 - 2^-6, 1): the exponent is -1 and the mantissa is within 2^-5 of two, which puts
        // its leading fraction bits in one of the last two table rows.
        constexpr size_t entries = size_t{1} << log2_reduction_bits;
        const bool just_below_one = (expo == -1) && (j >= entries - 2);

        // z, the reduced argument, as a 128 bit fraction. The series below needs |z| <= 2^-6.
        uint64_t z_low = 0, z_high = 0;
        uint32_t z_sign = 0;
        if (j == 0) {
            // Already inside the series' range, so no reduction is applied - and none is wanted.
            // f is exact here, and multiplying it by a rounded reciprocal would cost it every
            // significant bit it has when it is small. That is what makes log2 near 1 accurate:
            // the answer is then proportional to f, and f survives intact.
            z_low = frac_low << 16;
            z_high = (frac_high << 16) | (frac_low >> 48);
        } else if (just_below_one) {
            // The mirror image of the case above. Reducing from the exponent -1 would form the
            // answer as -1 plus the logarithm of a mantissa just under two, and that sum cancels: it
            // loses a bit for every power of two x sits closer to one, half the mantissa at 2^-60
            // below one and nearly all of it at 2^-110. Instead the series is taken at one itself,
            // with z = x - 1. That difference is exact - x is mantissa/2, which as a 128 bit fraction
            // has room for every bit of it - and the answer is then proportional to z, as above.
            const uint64_t m_low = frac_low << 15;
            const uint64_t m_high = (1ull << 63) | (frac_high << 15) | (frac_low >> 49);
            const uint8_t borrow = subborrow_u64(0, 0, m_low, &z_low);
            subborrow_u64(borrow, 0, m_high, &z_high);
            z_sign = 1;
        } else {
            // mantissa/2 as a 128 bit fraction, which is in [0.5,1)
            const uint64_t m_low = frac_low << 15;
            const uint64_t m_high = (1ull << 63) | (frac_high << 15) | (frac_low >> 49);
            uint64_t p_low = 0, p_high = 0;
            mul128_high(m_low, m_high, log2_recip_table[j][1], log2_recip_table[j][0], p_low, p_high);

            // z = 2 * (mantissa/2 * recip) - 1. Normally non negative; a reciprocal that rounded
            // down can take it just below zero.
            constexpr uint64_t one_half = 1ull << 63;
            uint64_t d_low = 0, d_high = 0;
            if (p_high >= one_half) {
                d_low = p_low;
                d_high = p_high - one_half;
            } else {
                d_low = 0ull - p_low;
                d_high = one_half - p_high - ((p_low != 0) ? 1ull : 0ull);
                z_sign = 1;
            }
            z_low = d_low << 1;
            z_high = (d_high << 1) | (d_low >> 63);
        }

        // Horner over 1/(n*ln2), from the last term down, giving A = log2(1+z)/z.
        //
        // Every value here is halved relative to the mathematics - log2_inv_n_table already holds
        // 1/(2n*ln2) - so that all of them stay below one and can be held as plain 128 bit
        // fractions. A is about 1.44 and would not fit otherwise. The accumulator stays inside
        // [0.019, 0.734] throughout, so neither the subtraction nor the addition can leave range.
        uint64_t acc_low = log2_inv_n_table[LOG2_TERMS - 1][1];
        uint64_t acc_high = log2_inv_n_table[LOG2_TERMS - 1][0];
        for (int32_t n = LOG2_TERMS - 1; n >= 1; --n) {
            uint64_t term_low = 0, term_high = 0;
            mul128_high(z_low, z_high, acc_low, acc_high, term_low, term_high);
            const uint64_t q_low = log2_inv_n_table[n - 1][1], q_high = log2_inv_n_table[n - 1][0];
            if (z_sign) {
                const uint8_t carry = addcarryx_u64(0, q_low, term_low, &acc_low);
                addcarryx_u64(carry, q_high, term_high, &acc_high);
            } else {
                const uint8_t borrow = subborrow_u64(0, q_low, term_low, &acc_low);
                subborrow_u64(borrow, q_high, term_high, &acc_high);
            }
        }

        // The rest of (0.5, 1) that the table reduces. The answer is -1 + T + S, where T is the table
        // value below and S the series, and it is as small as 0.0227, so the -1 still cancels up to
        // five and a half bits. Rounding T to a float128 first leaves an absolute error of 2^-114
        // for those bits to expose - measured, up to 9 ulp at 1 - 2^-4 and 42 at 1 - 2^-6 - so the
        // sum is formed on the exact 128 bit fractions instead, where every term is good to a few
        // units of 2^-128, and rounded once. That keeps this range within 0.75 ulp.
        if (expo == -1 && j != 0 && !just_below_one) {
            uint64_t s_low = 0, s_high = 0;
            mul128_high(z_low, z_high, acc_low, acc_high, s_low, s_high);
            s_high = (s_high << 1) | (s_low >> 63);  // undo the halving the table carries
            s_low <<= 1;

            // r = 1 - T - S = -log2(x). T is in (0,1), so 1 - T fits a fraction exactly.
            uint64_t r_low = 0, r_high = 0;
            const uint8_t borrow = subborrow_u64(0, 0, log2_value_table[j][1], &r_low);
            subborrow_u64(borrow, 0, log2_value_table[j][0], &r_high);
            if (z_sign) {
                const uint8_t carry = addcarryx_u64(0, r_low, s_low, &r_low);
                addcarryx_u64(carry, r_high, s_high, &r_high);
            } else {
                const uint8_t borrow_s = subborrow_u64(0, r_low, s_low, &r_low);
                subborrow_u64(borrow_s, r_high, s_high, &r_high);
            }

            float128 res = from_fraction128(r_low, r_high);
            res.set_sign(1);
            return filter(res);
        }

        // The product is formed in float128 rather than in the fraction arithmetic above. A fraction
        // grid is absolute, and log2(1+z) is proportional to z, so holding the product on that grid
        // would leave a small result with only as many significant bits as it has room above 2^-128.
        // Multiplying an exact z by an A that is accurate in its own right keeps the answer accurate
        // relative to its own size, which is what a floating point type is expected to do.
        float128 series = from_fraction128(z_low, z_high) * from_fraction128(acc_low, acc_high);
        series <<= 1;  // undo the halving the table carries
        series.set_sign(z_sign);

        // Taken at one, so the series is the whole answer: neither the exponent nor a table value
        // applies.
        if (just_below_one) {
            return filter(series);
        }

        // -log2(recip), the part of the answer the reduction removed. Zero when nothing was reduced.
        const float128 table_value = (j != 0) ? from_fraction128(log2_value_table[j][1], log2_value_table[j][0]) : float128();

        // The two fraction parts are summed before the exponent is added in. Both are below one, so
        // that first addition rounds against a small value; adding the exponent first would round
        // twice against a number as large as 16000 and cost the answer a bit for nothing.
        return filter(float128(expo) + (table_value + series));
    }
    /**
     * @brief Calculates Log base 10 of x: log10(x)
     * @param x The number to perform log on.
     * @return log10(x)
     */
    [[nodiscard]] friend FP128_INLINE float128 log10(float128 x) noexcept
    {
        float128 special;
        if (log_special(x, special))
            return special;
        // From here on only the flags the result justifies reach the caller.
        detail::flag_filter filter;

        const float128 t = x - float128::one();
        if (fabs(t) <= log1p_small_limit())
            return filter(log1p_small(t) * float128::log10_e());

        return filter(log2(x) * float128::log10_2());
    }
    /**
     * @brief Calculates Log base 2 of x as an integer ignoring the sign of x.
     * Similar to: floor(log2(fabs(x)))
     * @param x The number to perform log on.
     * @return logb(x)
     */
    [[nodiscard]] friend FP128_INLINE float128 logb(float128 x) noexcept
    {
        if (x.is_nan())
            return propagate_nan(x);
        if (x.is_inf())
            return fabs(x);
        if (x.is_zero()) {
            detail::raise_flags(FE_DIVBYZERO);
            return -inf();
        }

        // get_components() normalizes a subnormal, whose stored exponent field is the same for
        // every one of them and therefore says nothing about the value.
        uint64_t l = 0, h = 0;
        int32_t expo = 0;
        uint32_t sign = 0;
        x.get_components(l, h, expo, sign);
        return float128(expo);
    }
    /**
     * @brief Calculates the natural Log (base e) of 1 + x: log(1 + x)
     * @param x The number to perform log on.
     * @return log1p(x)
     */
    [[nodiscard]] friend FP128_INLINE float128 log1p(float128 x) noexcept
    {
        if (x.is_nan())
            return propagate_nan(x);
        if (x.is_zero() || (x.is_inf() && x.is_positive()))
            return x;

        const float128 one_value = float128::one();
        if (x < -one_value)
            return invalid_operation();
        if (x == -one_value) {
            detail::raise_flags(FE_DIVBYZERO);
            return -inf();
        }
        // From here on only the flags the result justifies reach the caller.
        detail::flag_filter filter;

        if (fabs(x) <= log1p_small_limit())
            return filter(log1p_small(x));

        // Above the series' range but still below one, 1+x rounds and the logarithm of a value
        // that close to one magnifies the error. Both subtractions below are exact - the first by
        // Sterbenz, the second because the operands agree to within a rounding - so the correction
        // puts back exactly what the rounding of 1+x took away.
        if (fabs(x) < one_value) {
            const float128 sum = one_value + x;
            const float128 rounding = (sum - one_value) - x;
            return filter(log(sum) - rounding / sum);
        }

        return filter(log(one_value + x));
    }

    //
    // Trigonometric functions
    //

    /**
     * @brief Calculate the sine function over a limited range [-0.5pi, 0.5pi]
     * Using the Maclaurin series expansion, the formula is:
     *              x^3   x^5   x^7
     * sin(x) = x - --- + --- - --- + ...
     *               3!    5!    7!
     * @param x value in Radians in the range [-0.5pi, 0.5pi]
     * @return Sine of x
     */
    [[nodiscard]] friend FP128_INLINE float128 sin1(float128 x) noexcept
    {
        assert(fabs(x) <= float128::half_pi());

        // first part of the series is just 'x'
        const float128 xx = x * x;
        float128 elem_denom, elem_nom = x;

        // Compute the rest of the series, starting with: -(x^3 / 3!). The reciprocal factorials stop
        // at 1/50!, and the argument is reduced long before the terms could need more; the bound
        // only makes sure an unreduced argument cannot keep the loop alive. One once could: past
        // 2^321 the powers overflow, inf * 0 is a NaN, and a NaN never compares equal to zero.
        for (int i = 3, sign = 1; i <= 51; i += 2, sign = 1 - sign) {
            elem_nom *= xx;
            fact_reciprocal(i, elem_denom);
            float128 elem = elem_nom * elem_denom;  // next element in the series
            // precision limit has been hit
            if (!elem)
                break;
            x += (sign) ? -elem : elem;
        }

        return x;
    }
    /**
     * @brief Calculate the cosine function over a limited range [-0.5pi, 0.5pi]
     * Since the sin1 function converges faster, call it with the modified angle.
     * @param x value in Radians in the range [-0.5pi, 0.5pi]
     * @return Cosine of x
     */
    [[nodiscard]] friend FP128_INLINE float128 cos1(const float128& x) noexcept
    {
        constexpr float128 half_pi = float128::half_pi();
        assert(fabs(x) <= half_pi);
        return (x.is_positive()) ? sin1(half_pi - x) : -sin1(-half_pi - x);
    }
    /**
     * @brief Calculate the Sine function
     * Ultimately uses sin() with a reduced range of [-pi/4, pi/4]
     * @param x value in Radians
     * @return Sine of x
     */
    /**
     * @name Unevaluated sums
     *
     * A value held as a pair of float128 whose sum is the quantity meant, with the second one
     * below the last place of the first. That is 226 bits of significand out of two 113 bit
     * values, which is what the trigonometric argument reduction needs: the residual it produces
     * has to be known to the absolute precision of the *result*, and near a zero of the sine the
     * result is many orders of magnitude below the argument it came from.
     * @{
     */

    /**
     * @brief Exact sum of two values, as a rounded sum and the error it made.
     *
     * Dekker's algorithm. Nothing here is an approximation: a + b is exactly sum + err for any two
     * finite operands.
     *
     * @param a First operand
     * @param b Second operand
     * @param sum Receives the rounded sum
     * @param err Receives a + b - sum, exactly
     */
    FP128_INLINE static void two_sum(const float128& a, const float128& b, float128& sum, float128& err) noexcept
    {
        sum = a + b;
        const float128 b_part = sum - a;
        err = (a - (sum - b_part)) + (b - b_part);
    }

    /**
     * @brief Exact sum of two values whose magnitudes are ordered.
     * @param a First operand, must not be smaller in magnitude than b
     * @param b Second operand
     * @param sum Receives the rounded sum
     * @param err Receives a + b - sum, exactly
     */
    FP128_INLINE static void quick_two_sum(const float128& a, const float128& b, float128& sum, float128& err) noexcept
    {
        sum = a + b;
        err = b - (sum - a);
    }
    /// @}

    /**
     * @brief Nearest integer to x / (pi/2), the quadrant the argument falls in.
     *
     * Beyond 2^62 the quadrant no longer fits in the integer the reduction counts with, and an
     * argument that large has fewer bits below the binary point than a quadrant is wide - its
     * sine is not determined by the value in any useful sense. Everything up to there reduces
     * exactly, see reduce_half_pi().
     *
     * @param x Argument, below 2^60 in magnitude; reduce_large() takes the larger ones
     * @return Quadrant index, zero for an argument too large to reduce.
     */
    [[nodiscard]] FP128_INLINE static int64_t quadrant_of(const float128& x) noexcept
    {
        constexpr float128 two_over_pi(0x2A53F84EAFA3EA6A, 0x45F306DC9C88, 0x3FFE, 0);
        const float128 scaled = x * two_over_pi;
        if (fabs(scaled) >= ldexp(float128::one(), 62))
            return 0;
        return llround(scaled);
    }

    /**
     * @brief Reduces an argument of any size: x = n * pi/2 + r, with |r| <= pi/4 (Payne and Hanek).
     *
     * Subtracting multiples of a stored pi/2, as reduce_half_pi() does, needs the multiple n as an
     * integer and pi/2 to as many bits as n has plus the 226 the residual keeps - out of reach once
     * n runs to thousands of bits. Here x * 2/pi is formed instead, in integer arithmetic, from the
     * 113 bit mantissa and a 384 bit window of the binary expansion of 2/pi chosen by the exponent.
     * The bits of 2/pi above the window only add multiples of four to the product, which leave the
     * quadrant where it is; the bits below it only reach past the fraction the window keeps.
     *
     * That fraction has 380 bits or more, which leaves 226 significant ones even when x is closer
     * to a multiple of pi/2 than any binary128 argument can be: for a random argument the fraction
     * leads with about as many zero bits as the logarithm of the count of binary128 values, 127.
     * Rounded to the nearest quadrant and multiplied by pi/2, it becomes the residual.
     *
     * Before this existed sin(), cos() and tan() reduced nothing above 2^62, and returned values
     * such as 2^3240 for sin(1.5e21). Above 2^321 the series they evaluated never terminated.
     *
     * @param x Argument, finite and non zero
     * @param hi Receives the residual
     * @param lo Receives what did not fit in hi
     * @return n mod 4, the quadrant.
     */
    FP128_INLINE static int32_t reduce_large(const float128& x, float128& hi, float128& lo) noexcept
    {
        constexpr float128 half_pi_hi(0x8469898CC51701B8, 0x921FB54442D1, 0x3FFF, 0);
        constexpr float128 half_pi_mid(0xA67CC74020BBEA64, 0xCD129024E088, 0x3F8C, 0);
        constexpr int32_t window_words = 6;
        constexpr int32_t window_bits = window_words * 64;
        constexpr int32_t product_words = window_words + 2;

        uint64_t ml = 0, mh = 0;
        int32_t e = 0;
        uint32_t sign = 0;
        x.get_components(ml, mh, e, sign);

        // |x| = m * 2^(e-112) and bit i of 2/pi is worth 2^-i, so bit i times m lands on 2^(e-112-i).
        // The window starts at the bit that lands on 2^1; the ones above it are multiples of four.
        const int32_t first = (e - 113 > 1) ? (e - 113) : 1;
        const int32_t frac_bits = (first + window_bits - 1) - (e - FRAC_BITS);

        // The window, least significant word first.
        uint64_t window[window_words] {};
        const int32_t word = (first - 1) >> 6;
        const int32_t bit = (first - 1) & 63;
        for (int32_t k = 0; k < window_words; ++k) {
            const uint64_t w0 = detail::two_over_pi_bits[word + k];
            const uint64_t w1 = detail::two_over_pi_bits[word + k + 1];
            window[window_words - 1 - k] = (bit == 0) ? w0 : ((w0 << bit) | (w1 >> (64 - bit)));
        }

        // The product of the mantissa and the window, scaled by 2^frac_bits.
        uint64_t p[product_words] {};
        const uint64_t m[2] = {ml, mh};
        for (int32_t i = 0; i < 2; ++i) {
            uint64_t carry = 0;
            for (int32_t j = 0; j < window_words; ++j) {
                uint64_t high_word = 0;
                uint64_t low_word = mulx_u64(m[i], window[j], &high_word);
                uint8_t c = addcarryx_u64(0, low_word, carry, &low_word);
                high_word += c;
                c = addcarryx_u64(0, p[i + j], low_word, &p[i + j]);
                carry = high_word + c;
            }
            p[i + window_words] = carry;
        }

        // frac_bits is at most 436 for the exponents that reach here, inside the 512 bits of p; the
        // bound keeps that true for any argument.
        const auto bit_at = [&p](int32_t index) {
            return ((index >> 6) < product_words) ? static_cast<uint32_t>((p[index >> 6] >> (index & 63)) & 1) : 0u;
        };
        uint32_t quadrant = bit_at(frac_bits) | (bit_at(frac_bits + 1) << 1);

        // Keep the fraction only.
        for (int32_t i = 0; i < product_words; ++i) {
            const int32_t lowest = i * 64;
            if (lowest >= frac_bits)
                p[i] = 0;
            else if (lowest + 64 > frac_bits)
                p[i] &= FP128_MAX_VALUE_64(frac_bits - lowest);
        }

        // A fraction of one half or more belongs to the next quadrant, as a negative residual.
        const bool negative = bit_at(frac_bits - 1) != 0;
        if (negative) {
            ++quadrant;
            uint8_t borrow = 0;
            for (int32_t i = 0; i < product_words; ++i)
                borrow = subborrow_u64(borrow, 0, p[i], &p[i]);
            for (int32_t i = 0; i < product_words; ++i) {
                const int32_t lowest = i * 64;
                if (lowest >= frac_bits)
                    p[i] = 0;
                else if (lowest + 64 > frac_bits)
                    p[i] &= FP128_MAX_VALUE_64(frac_bits - lowest);
            }
        }

        int32_t msb = -1;
        for (int32_t i = product_words - 1; i >= 0 && msb < 0; --i) {
            if (p[i] != 0)
                msb = i * 64 + 63 - static_cast<int32_t>(lzcnt64(p[i]));
        }

        hi = float128();
        lo = float128();
        if (msb >= 0) {
            // 128 bits of the fraction starting at bit pos, zeros below bit 0.
            const auto extract = [&p](int32_t pos, uint64_t& l, uint64_t& h) {
                const int32_t from = (pos > 0) ? pos : 0;
                const int32_t w = from >> 6;
                const int32_t b = from & 63;
                const uint64_t p0 = (w < product_words) ? p[w] : 0;
                const uint64_t p1 = (w + 1 < product_words) ? p[w + 1] : 0;
                const uint64_t p2 = (w + 2 < product_words) ? p[w + 2] : 0;
                l = (b == 0) ? p0 : ((p0 >> b) | (p1 << (64 - b)));
                h = (b == 0) ? p1 : ((p1 >> b) | (p2 << (64 - b)));
                if (pos < 0)
                    shift_left128_inplace_safe(l, h, -pos);
                h &= MAX_SIG_HIGH;  // 113 bits
            };

            // The top 113 bits of the fraction and the 113 below them, both exact.
            uint64_t l = 0, h = 0;
            extract(msb - FRAC_BITS, l, h);
            float128 f_hi;
            f_hi.set_components(l, h, msb - frac_bits, 0);
            extract(msb - 2 * FRAC_BITS - 1, l, h);
            const float128 f_lo = norm_round_pack(0, msb - FRAC_BITS - 1 - frac_bits, l, h, false);

            // r = f * pi/2, with the error of the leading product recovered exactly.
            const float128 product = f_hi * half_pi_hi;
            const float128 product_err = fma(f_hi, half_pi_hi, -product);
            quick_two_sum(product, product_err + (f_hi * half_pi_mid + f_lo * half_pi_hi), hi, lo);
            if (negative) {
                hi = -hi;
                lo = -lo;
            }
        }

        // sin and cos of -x follow from the reduction of x mirrored.
        if (sign != 0) {
            hi = -hi;
            lo = -lo;
            quadrant = 0u - quadrant;
        }
        return static_cast<int32_t>(quadrant & 3);
    }

    /**
     * @brief Reduces x to [-pi/4, pi/4], as an unevaluated sum, and returns the quadrant.
     * @param x Argument, finite
     * @param hi Receives the residual
     * @param lo Receives what did not fit in hi
     * @return n mod 4, where x = n * pi/2 + hi + lo.
     */
    FP128_INLINE static int32_t reduce_trig(const float128& x, float128& hi, float128& lo) noexcept
    {
        // Up to 2^60 the multiple fits the integer reduce_half_pi() counts with, and its 342 bits
        // of pi/2 leave the residual 226 good bits. Beyond that the multiple is formed in integer
        // arithmetic instead.
        if (x.get_exponent() < 60) {
            const int64_t n = quadrant_of(x);
            reduce_half_pi(x, n, hi, lo);
            return static_cast<int32_t>(n & 3);
        }
        return reduce_large(x, hi, lo);
    }

    /**
     * @brief Subtracts n*(pi/2) from x, leaving the residual as an unevaluated sum.
     *
     * pi/2 is irrational, so no single float128 multiple of it can be subtracted without leaving
     * an error of n ulp behind - and the sine of the residual is only as accurate as the residual
     * itself. Where the answer is small, which is exactly where an argument lands near a multiple
     * of pi/2, that error is the whole answer: sin() used to return a value with 55 correct bits
     * for an argument a few ulp away from pi.
     *
     * Two things fix it. pi/2 is carried in three pieces covering 342 bits rather than 113, and
     * each product n*piece is split into its rounded value and its exact error with fma(), so no
     * multiplication is lost either. The residual accumulates as a pair, giving it 226 bits.
     *
     * @param x Argument
     * @param n Multiple of pi/2 to remove
     * @param hi Receives the residual
     * @param lo Receives what did not fit in hi
     */
    FP128_INLINE static void reduce_half_pi(const float128& x, int64_t n, float128& hi, float128& lo) noexcept
    {
        constexpr float128 half_pi_parts[3] = {
            float128(0x8469898CC51701B8, 0x921FB54442D1, 0x3FFF, 0),
            float128(0xA67CC74020BBEA64, 0xCD129024E088, 0x3F8C, 0),
            float128(0xE19C72FEC8841ABA, 0x3B19376BAD7D, 0x3F1A, 1),
        };

        hi = x;
        lo = float128();
        if (n == 0)
            return;

        const float128 nf(n);
        for (const float128& part : half_pi_parts) {
            // n * part is exactly product + product_err
            const float128 product = nf * part;
            const float128 product_err = fma(nf, part, -product);

            float128 sum, err;
            two_sum(hi, -product, sum, err);
            // Everything that did not fit is collected and folded back in, which keeps the pair
            // normalized for the next piece.
            quick_two_sum(sum, (err + lo) - product_err, hi, lo);
        }
    }

    [[nodiscard]] friend float128 sin(float128 x) noexcept
    {
        if (x.is_nan())
            return propagate_nan(x);
        if (x.is_zero())
            return x;
        if (x.is_inf())
            return invalid_operation();
        // From here on only the flags the result justifies reach the caller.
        detail::flag_filter filter;

        float128 hi, lo;
        const int32_t n = reduce_trig(x, hi, lo);

        switch (n) {
        case 0:  // [-45-45) degrees
            return filter(sin1(hi) + cos1(hi) * lo);
        case 1:  // [45-135) degrees
            return filter(cos1(hi) - sin1(hi) * lo);
        case 2:  // [135-225) degrees
            return filter(-(sin1(hi) + cos1(hi) * lo));
        case 3:  // [225-315) degrees
        default:
            return filter(-(cos1(hi) - sin1(hi) * lo));
        }
    }
    /**
     * @brief Calculate the inverse sine function
     *
     * asin(x) = atan(x / sqrt((1-x)(1+x))). The Newton iteration on sin() this replaced divided by
     * cos(), the derivative, which vanishes at the ends of the domain: it converged slowly there,
     * to about 58 correct bits at asin(1). Near the ends 1-x is exact, so the argument of atan()
     * keeps its accuracy however close to one x gets.
     *
     * @param x value in the range [-1,1]
     * @return Inverse sine of x. Outside the domain the result is a NaN, the invalid operation.
     */
    [[nodiscard]] friend float128 asin(float128 x) noexcept
    {
        if (x.is_nan())
            return propagate_nan(x);
        const float128 one_value = float128::one();
        const float128 a = fabs(x);
        if (a > one_value)
            return invalid_operation();
        if (x.is_zero())
            return x;
        // From here on only the flags the result justifies reach the caller.
        detail::flag_filter filter;

        float128 res;
        if (a == one_value)
            res = float128::half_pi();
        else if (a < ldexp(one_value, -57)) {
            res = a;  // the cubic term is below half an ulp
            detail::raise_flags(FE_INEXACT);
        }
        else
            res = atan(a / sqrt((one_value - a) * (one_value + a)));

        res.set_sign(x.get_sign());
        return filter(res);
    }
    /**
     * @brief Calculate the cosine function
     * Ultimately uses sin1() with a reduced range of [-pi/4, pi/4]
     * Sine's Maclaurin series converges faster than Cosine's.
     * @param x value in Radians
     * @return Cosine of x
     */
    [[nodiscard]] friend float128 cos(float128 x) noexcept
    {
        if (x.is_nan())
            return propagate_nan(x);
        if (x.is_inf())
            return invalid_operation();
        // From here on only the flags the result justifies reach the caller.
        detail::flag_filter filter;

        float128 hi, lo;
        const int32_t n = reduce_trig(x, hi, lo);

        switch (n) {
        case 0:  // [-45-45) degrees
            return filter(cos1(hi) - sin1(hi) * lo);
        case 1:  // [45-135) degrees
            return filter(-(sin1(hi) + cos1(hi) * lo));
        case 2:  // [135-225) degrees
            return filter(-(cos1(hi) - sin1(hi) * lo));
        case 3:  // [225-315) degrees
        default:
            return filter(sin1(hi) + cos1(hi) * lo);
        }
    }
    /**
     * @brief Calculate the inverse cosine function
     *
     * acos(x) = 2 * atan(sqrt((1-x)/(1+x))). The Newton iteration this replaced divided by sin(),
     * which is zero at acos(1): that returned a NaN, and so did acos(1 - 2^-100). Near x = 1 the
     * difference 1-x is exact, and near x = -1 the sum 1+x is, so the quotient keeps its accuracy
     * at both ends.
     *
     * @param x value in the range [-1,1]
     * @return Inverse cosine of x. Outside the domain the result is a NaN, the invalid operation.
     */
    [[nodiscard]] friend float128 acos(float128 x) noexcept
    {
        if (x.is_nan())
            return propagate_nan(x);
        const float128 one_value = float128::one();
        if (fabs(x) > one_value)
            return invalid_operation();
        if (x == one_value)
            return float128();
        if (x == -one_value)
            return float128::pi();
        // From here on only the flags the result justifies reach the caller.
        detail::flag_filter filter;

        return filter(atan(sqrt((one_value - x) / (one_value + x))) << 1);
    }
    /**
     * @brief Calculate the tangent function
     * tan(x) = sin(x)/cos(x)
     * @param x value
     * @return Tangent of x
     */
    [[nodiscard]] friend FP128_FORCE_INLINE float128 tan(float128 x) noexcept { return sin(x) / cos(x); }
    /**
     * @brief Calculate the inverse tangent function
     * @param x value
     * @return Arctangent of x
     */
    [[nodiscard]] friend float128 atan(float128 x) noexcept
    {
        if (x.is_nan())
            return propagate_nan(x);
        if (x.is_zero())
            return x;
        // From here on only the flags the result justifies reach the caller.
        detail::flag_filter filter;

        // constants for segmentation
        constexpr float128 half_pi = float128::half_pi();  // pi / 2
        bool comp = false;
        constexpr int max_iterations = 6;

        // make argument positive, save the sign
        auto sign = x.get_sign();
        x.set_sign(0);

        // limit argument to 0..1. x is greater than one here, so the division cannot be by zero.
        if (x > 1) {
            comp = true;
            x = one() / x;
        }

        // initial step uses the CRT function.
        float128 res = ::atan(static_cast<double>(x));
        const float128 eps = fabs(res >> 110);

        //
        // Xn+1 =  Xn - cos(Xn) * ( sin(Xn) - a * cos(Xn))
        //
        // where 'a' is the argument, each iteration will converge on the result if the initial
        //  estimate is close enough.
        for (int i = 0; i < max_iterations; ++i) {
            float128 cos_xn = cos(res);
            float128 sin_xn = sin(res);
            float128 e = cos_xn * (sin_xn - x * cos_xn);  // this is the iteration estimated error
            res -= e;
            if (fabs(e) <= eps)
                break;
        }

        // restore complement if needed
        if (comp)
            res = half_pi - res;
        // restore sign if needed
        res.set_sign(sign);
        return filter(res);
    }
    /**
     * @brief Calculate the inverse tangent function of the ratio y / x
     * @param y value
     * @param x value
     * @return Arctangent of y / x in the range [-pi, pi]. The zeros and infinities follow
     *         IEEE 754-2008 9.2.1: the sign of a zero y is kept, a negative x (-0 included) puts
     *         the result at +-pi, and two infinities give +-pi/4 or +-3pi/4.
     */
    [[nodiscard]] friend float128 atan2(float128 y, float128 x) noexcept
    {
        // constants for segmentation
        constexpr float128 pi = float128::pi();
        constexpr float128 half_pi = float128::half_pi();  // pi / 2
        constexpr float128 three_quarter_pi(0x234F272993D1414A, 0x2D97C7F3321D, 0x4000, 0);

        if (x.is_nan() || y.is_nan())
            return propagate_nan(y, x);

        // The zeros and the infinities. Each result takes the sign of y, a zero y included; the
        // previous code returned pi for atan2(0, 1) and lost the sign of every zero.
        float128 res;
        if (y.is_zero() || x.is_zero() || y.is_inf() || x.is_inf()) {
            if (y.is_zero())
                res = x.is_negative() ? pi : float128();  // on the axis, the left half including -0
            else if (y.is_inf())
                res = x.is_inf() ? (x.is_negative() ? three_quarter_pi : float128::quarter_pi()) : half_pi;
            else if (x.is_zero())
                res = half_pi;
            else
                res = x.is_negative() ? pi : float128();  // y finite, x infinite
            res.set_sign(y.get_sign());
            return res;
        }
        // From here on only the flags the result justifies reach the caller.
        detail::flag_filter filter;

        // save the signs of x, y
        bool comp = fabs(y) > fabs(x);
        float128 ratio;

        // calculate the ratio keeping it below 1.0
        ratio = (comp) ? x / y : y / x;
        res = atan(ratio);

        if (comp)
            res = (res.is_negative()) ? -half_pi - res : half_pi - res;

        if (x > 0)
            return filter(res);

        // x < 0
        return filter((y < 0) ? res - pi : res + pi);
    }
    /**
     * @brief Calculate the hyperbolic sine function
     * Use the exponent function which produces more accurate results than the power series.
     *           e^x - e^(-x)
     * sinh(x) = ------------
     *                2
     * @param x value
     * @return Sine of x
     */
    [[nodiscard]] friend FP128_INLINE float128 sinh(const float128 x) noexcept
    {
        if (x.is_nan())
            return propagate_nan(x);
        if (x.is_inf() || x.is_zero())
            return x;
        // From here on only the flags the result justifies reach the caller.
        detail::flag_filter filter;

        const float128 a = fabs(x);
        float128 res;
        // Below one, exp(a) and exp(-a) agree to within a factor of e and their difference cancels
        // away most of the mantissa - at a = 2^-60 all of it. expm1() holds the same quantity with
        // no cancellation at all, and t/(t+1) is exp(-a)-1 expressed through it.
        if (a < float128::one()) {
            const float128 t = expm1(a);
            res = (t + t / (t + float128::one())) >> 1;
        } else if (a < float128(11350)) {
            res = (exp(a) - exp(-a)) >> 1;
        } else {
            // exp(a) overflows before sinh(a), which is half of it, does. Squaring exp(a/2) reaches
            // the same value without the intermediate overflow; exp(-a) is far below its last place.
            const float128 half_power = exp(a >> 1);
            res = half_power * (half_power >> 1);
        }

        return filter((x.is_negative()) ? -res : res);
    }
    /**
     * @brief Calculates the inverse hyperbolic sine
     * For positive x:
     * asinh(x) = log(x + sqrt(x^2 + 1))
     * For negative x, the function returns the result with the sign inverted
     * @param x value
     * @return Inverse hyperbolic sine of x
     */
    [[nodiscard]] friend FP128_INLINE float128 asinh(const float128 x) noexcept
    {
        if (x.is_nan())
            return propagate_nan(x);
        if (x.is_inf() || x.is_zero())
            return x;
        // From here on only the flags the result justifies reach the caller.
        detail::flag_filter filter;

        const float128 one_value = float128::one();
        const float128 a = fabs(x);
        float128 res;

        // Three regimes, because the closed form is only usable in the middle one.
        //
        // Below 2^-57 the answer is the argument to within half an ulp, and the closed form would
        // square it to zero and hand back log(1) = 0 instead.
        //
        // Above 2^57 the square overflows long before the argument does, and asinh(a) is log(2a)
        // to well past the last bit.
        //
        // In between, adding a to sqrt(a*a+1) cancels for a negative argument, which is why the
        // sign is taken out first, and the a*a/(1+sqrt(1+a*a)) form is used rather than
        // sqrt(a*a+1)-a: near zero the latter subtracts two values that agree to every bit the
        // result needs.
        if (a < ldexp(one_value, -57)) {
            res = a;
            detail::raise_flags(FE_INEXACT);  // asinh(a) is a hair below a, not equal to it
        } else if (a > ldexp(one_value, 57)) {
            res = log(a) + float128::ln2();
        } else {
            const float128 a2 = sqr(a);
            res = log1p(a + a2 / (one_value + sqrt(one_value + a2)));
        }

        return filter((x.is_positive()) ? res : -res);
    }
    /**
     * @brief Calculate the hyperbolic cosine function over a limited range [-0.5pi, 0.5pi]
     *           e^x + e^(-x)
     * cosh(x) = ------------
     *                2
     * @param x value in Radians in the range [-0.5pi, 0.5pi]
     * @return Sine of x
     */
    [[nodiscard]] friend FP128_INLINE float128 cosh(const float128 x) noexcept
    {
        if (x.is_nan())
            return propagate_nan(x);
        const float128 a = fabs(x);
        if (a.is_inf())
            return a;
        // From here on only the flags the result justifies reach the caller.
        detail::flag_filter filter;
        // As in sinh(): exp(a) overflows a little before cosh(a) = exp(a)/2 does.
        if (a > float128(11350)) {
            const float128 half_power = exp(a >> 1);
            return filter(half_power * (half_power >> 1));
        }
        return filter((exp(a) + exp(-a)) >> 1);
    }
    /**
     * @brief Calculates the inverse hyperbolic cosine
     * For x >= 1:
     * acosh(x) = log(x + sqrt(x^2 - 1))
     * For x < 1, the function return zero
     * @param x value in the range [1, inf]
     * @return Inverse hyperbolic cosine of x
     */
    [[nodiscard]] friend FP128_INLINE float128 acosh(const float128 x) noexcept
    {
        if (x.is_nan())
            return propagate_nan(x);
        // Outside the domain the result is undefined rather than zero, which is what the previous
        // version returned for every argument below one.
        const float128 one_value = float128::one();
        if (x < one_value)
            return invalid_operation();
        if (x.is_inf())
            return x;
        // From here on only the flags the result justifies reach the caller.
        detail::flag_filter filter;

        // Near one, x*x-1 cancels: at x = 1+2^-100 the square rounds back to one and the root
        // comes out zero. The argument of log1p below is built from t = x-1, which is exact over
        // that range, so the answer keeps its significant bits.
        if (x < float128(2)) {
            const float128 t = x - one_value;
            return filter(log1p(t + sqrt((t << 1) + sqr(t))));
        }
        // Beyond 2^57 the square overflows before the argument does and acosh(x) is log(2x).
        if (x > ldexp(one_value, 57))
            return filter(log(x) + float128::ln2());

        return filter(log(x + sqrt(sqr(x) - one_value)));
    }
    /**
     * @brief Calculates the hyperbolic tangent
     *           e^x - e^(-x)
     * tanh(x) = ------------
     *           e^x + e^(-x)
     * @param x value
     * @return hyperbolic tangent of x
     */
    [[nodiscard]] friend FP128_INLINE float128 tanh(const float128 x) noexcept
    {
        if (x.is_nan())
            return propagate_nan(x);
        if (x.is_zero())
            return x;
        // From here on only the flags the result justifies reach the caller.
        detail::flag_filter filter;

        const float128 a = fabs(x);
        float128 res;
        // Past 40 the two exponentials differ by more than the mantissa can express and the
        // quotient is one to the last bit.
        if (x.is_inf() || a > float128(40)) {
            res = float128::one();
        } else {
            //           e^x - e^(-x)         e^(2x) - 1
            // tanh(x) = ------------  =  ----------------
            //           e^x + e^(-x)      (e^(2x) - 1) + 2
            //
            // The right hand form is written in terms of expm1, which holds e^(2x)-1 without the
            // cancellation that subtracting the two exponentials suffers for a small argument.
            const float128 t = expm1(a << 1);
            res = t / (t + float128(2));
        }

        return filter((x.is_negative()) ? -res : res);
    }
    /**
     * @brief Calculates the inverse hyperbolic tangent
     *                       1 + x
     * atanh(x) = 0.5 * log( -----)
     *                       1 - x
     * @param x value in the range (-1, 1)
     * @return Inverse hyperbolic tangent of x
     */
    [[nodiscard]] friend FP128_INLINE float128 atanh(const float128 x) noexcept
    {
        if (x.is_nan())
            return propagate_nan(x);
        if (x.is_zero())
            return x;

        constexpr auto one_value = float128::one();
        const float128 a = fabs(x);
        // The endpoints are poles and outside them the function is undefined. Returning zero for
        // all three, as this used to, gave atanh(1) a finite value and atanh(2) a defined one.
        if (a > one_value)
            return invalid_operation();
        if (a == one_value) {
            detail::raise_flags(FE_DIVBYZERO);
            return x.is_negative() ? -inf() : inf();
        }
        // From here on only the flags the result justifies reach the caller.
        detail::flag_filter filter;

        // 2*atanh(a) = log1p(2a/(1-a)). Splitting the argument as below keeps the numerator from
        // rounding away for a small a, where the answer is proportional to it: the previous form,
        // log((1+x)/(1-x)), formed a quotient that rounds to exactly one and then took its
        // logarithm, which is zero.
        float128 res;
        if (a < float128::half())
            res = log1p((a << 1) + (sqr(a) << 1) / (one_value - a)) >> 1;
        else
            res = log1p((a << 1) / (one_value - a)) >> 1;

        return filter((x.is_negative()) ? -res : res);
    }
    /**
     * @brief Calculates the Maclaurin series constants for the erf function.
     * The array will hold  1 / (n! * (2n + 1))
     * @param a pointer to array that receives the results. The array must be preallocated.
     * @param count Element count in the array
     */
    friend void erf_constants(float128* a, int32_t count)
    {
        if (a == nullptr)
            return;

        a[0] = 1;
        float128 f = 1;  // value of 0!
        for (int32_t i = 1; i < count; ++i) {
            f *= i;
            a[i] = float128::one() / (f * (2 * i + 1));
        }
    }
    /**
     * @brief Computes the error function of a value.
     * The erf function return a value in the range -1.0 to 1.0.
     * There's no error return.
     * @param x A floating-point value.
     * @return The erf functions return the Gauss error function of x.
     */
    [[nodiscard]] friend FP128_INLINE float128 erf(float128 x) noexcept
    {
        if (x.is_nan())
            return propagate_nan(x);
        if (x.is_zero())
            return x;
        // From here on only the flags the result justifies reach the caller.
        detail::flag_filter filter;

        const uint32_t sign = x.get_sign();
        const float128 a = fabs(x);
        float128 res;

        // erfc(9) is 4.1e-37, which is below half the last place of one, so erf is one from here
        // on and the series need never be run outside the range it converges quickly over.
        if (x.is_inf() || a >= float128(9)) {
            res = float128::one();
        } else {
            //                        2                   n
            //           2x        -x     inf         (2x^2)
            // erf(x) = ------ * e     *  sum   ------------------
            //          sqrt(pi)          n=0    1*3*5*...*(2n+1)
            //
            // Every term of this series is positive, unlike the Maclaurin series in x^(2n+1) that
            // this replaced: that one alternates, and its terms peak around n = x^2 at a value
            // 2^24 times the sum by the time x reaches 4, so the answer was formed by cancelling
            // away 24 bits of it. Past x = 2 the old code did not use a series at all - it called
            // the CRT's double erf and returned 53 correct bits.
            const float128 xx = sqr(a);
            // The square rounds, and exp() magnifies that by x^2, which reaches 81 here. Keeping
            // the exact error of the multiply and applying it as a factor puts those bits back.
            const float128 xx_err = fma(a, a, -xx);
            const float128 two_xx = xx << 1;

            // The sum runs to eighty terms near the top of the range and every one of them rounds
            // against a running total, which alone would cost the answer seven bits. two_sum()
            // hands back what each addition lost, so it can be carried and added in at the end.
            float128 term = float128::one();
            float128 sum = term;
            float128 lost;
            for (int32_t n = 1; n < 1024; ++n) {
                term = term * two_xx / float128(2 * n + 1);
                float128 rounded, err;
                two_sum(sum, term, rounded, err);
                sum = rounded;
                lost += err;
                if (term.get_exponent() + FRAC_BITS + 8 < sum.get_exponent())
                    break;
            }
            sum += lost;

            const float128 gaussian = exp(-xx) * (float128::one() - xx_err);
            res = ((a * sum) << 1) * float128::inv_sqrt_pi() * gaussian;
        }

        res.set_sign(sign);
        return filter(res);
    }
    /**
     * @brief Computes the complementary error function, 1 - erf(x).
     *
     * Above one the answer is a small difference of two values close to one, so it cannot be
     * reached through erf() - at x = 6 the subtraction would leave 55 of the 113 bits, and past
     * 13 nothing at all. It is computed directly instead, from Laplace's continued fraction
     *
     *                       2
     *                     -x
     *                    e                 1
     *      erfc(x) = ---------- * ------------------------
     *                 sqrt(pi)     x + (1/2)/(x + 1/(x + ...))
     *
     * evaluated with the modified Lentz algorithm. The fraction converges for every x above one,
     * and fastest where the series form is weakest.
     *
     * @param x Input value
     * @return 1 - erf(x)
     */
    [[nodiscard]] friend FP128_INLINE float128 erfc(float128 x) noexcept
    {
        if (x.is_nan())
            return propagate_nan(x);
        if (x.is_inf())
            return x.is_negative() ? float128(2) : float128();
        if (x.is_zero())
            return float128::one();
        // From here on only the flags the result justifies reach the caller.
        detail::flag_filter filter;

        const float128 one_value = float128::one();
        // erfc is symmetric about one: the negative half is where the answer approaches two and
        // keeps every bit, so it is the reflection that is well conditioned there.
        if (x.is_negative())
            return filter(float128(2) - erfc(-x));
        if (x < one_value)
            return filter(one_value - erf(x));
        // exp(-x*x) reaches the smallest subnormal at 106.9
        if (x > float128(107))
            return filter.inexact(float128());

        constexpr int32_t max_iterations = 4096;
        const float128 tiny = ldexp(one_value, -16000);
        const float128 epsilon = ldexp(one_value, -FRAC_BITS - 8);

        float128 f = tiny;
        float128 c = f;
        float128 d;
        for (int32_t n = 1; n < max_iterations; ++n) {
            // a_1 is one, every a after it is (n-1)/2
            const float128 a = (n == 1) ? one_value : (float128(n - 1) >> 1);

            d = x + a * d;
            if (d.is_zero())
                d = tiny;
            c = x + a / c;
            if (c.is_zero())
                c = tiny;
            d = one_value / d;

            const float128 delta = c * d;
            f *= delta;
            if (fabs(delta - one_value) < epsilon)
                break;
        }

        const float128 xx = sqr(x);
        const float128 xx_err = fma(x, x, -xx);
        const float128 gaussian = exp(-xx) * (one_value - xx_err);
        return filter(gaussian * float128::inv_sqrt_pi() * f);
    }
    /**
     * @brief Tests whether x is finite (not infinite and not NaN).
     * @param x Input value
     * @return Non-zero if finite, zero otherwise.
     */
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr bool isfinite(const float128& x) noexcept { return x.is_finite(); }
    /**
     * @brief Computes (x * y) + z with a single rounding at the end.
     *
     * The whole point of a fused multiply add is that the product is not rounded before the
     * addend reaches it, which is what makes the result correct to the last bit and what lets a
     * caller recover the error of a multiplication as fma(x, y, -x*y). Writing it as x * y + z
     * rounds twice and loses everything the addend cancels against.
     *
     * The exact product of the two 113 bit mantissas is 226 bits, and the addend is another 113.
     * Both are placed in a wide accumulator with their leading bits aligned, added or subtracted
     * there, and the single rounding is applied when the result is read back out. An operand too
     * far below the other one to reach the accumulator cannot change any bit of the result except
     * through the tie break, so it is folded into the sticky bit.
     *
     * @param x The first value to multiply.
     * @param y The second value to multiply.
     * @param z The value to add to the product.
     * @return (x * y) + z, correctly rounded.
     */
    [[nodiscard]] friend FP128_INLINE constexpr float128 fma(float128 x, float128 y, float128 z) noexcept
    {
        const uint32_t product_sign = x.get_sign() ^ y.get_sign();

        // The product of a zero and an infinity is the invalid operation even when the addend is a
        // quiet NaN (IEEE 754-2008 7.2 leaves that case to the implementation; signaling it is what
        // x86 does). Otherwise a NaN operand propagates.
        if ((x.is_inf() && y.is_zero()) || (y.is_inf() && x.is_zero())) {
            if (z.is_nan())
                detail::raise_flags(FE_INVALID);
            return z.is_nan() ? quiet_nan(z) : invalid_operation();
        }
        if (x.is_nan() || y.is_nan() || z.is_nan())
            return propagate_nan(x, y, z);

        // An infinite product decides the result on its own, unless the addend is the opposite
        // infinity and the sum is undefined.
        if (x.is_inf() || y.is_inf()) {
            if (z.is_inf() && z.get_sign() != product_sign)
                return invalid_operation();
            float128 res = inf();
            res.set_sign(product_sign);
            return res;
        }
        if (z.is_inf())
            return z;

        // A zero product leaves the addend. Two zeros give a zero that is negative only when both
        // of them are, or, rounding downward, when either is: the sign rule for an exact zero sum.
        if (x.is_zero() || y.is_zero()) {
            if (!z.is_zero())
                return z;
            float128 res;
            res.set_sign((z.get_sign() == product_sign) ? product_sign : exact_zero_sign());
            return res;
        }
        // A zero addend leaves the product, which the multiply already rounds correctly.
        if (z.is_zero())
            return x * y;

        uint64_t lx = 0, hx = 0, ly = 0, hy = 0, lz = 0, hz = 0;
        int32_t ex = 0, ey = 0, ez = 0;
        uint32_t sx = 0, sy = 0, sz = 0;
        x.get_components(lx, hx, ex, sx);
        y.get_components(ly, hy, ey, sy);
        z.get_components(lz, hz, ez, sz);

        // The exact 226 bit product of the two 113 bit mantissas. Both were normalized to
        // [2^112, 2^113), so the product is in [2^224, 2^226) and its value is
        // product * 2^(ex + ey - 224).
        uint64_t product[4] {};
        mul128to256(lx, hx, ly, hy, product);
        const int32_t product_bits = wide_msb(product, 4);

        // Weight of the leading bit of each operand, which is what they are aligned by.
        const int32_t product_msb = ex + ey - 2 * FRAC_BITS + product_bits;
        const int32_t top = (product_msb > ez) ? product_msb : ez;

        // Bit WIDE_TOP of an accumulator carries the weight 2^top.
        uint64_t product_acc[WIDE_WORDS] {};
        uint64_t addend_acc[WIDE_WORDS] {};
        const uint64_t addend[2] = {lz, hz};
        bool sticky = wide_place(product_acc, product, 4, WIDE_TOP - (top - product_msb) - product_bits);
        sticky = wide_place(addend_acc, addend, 2, WIDE_TOP - (top - ez) - FRAC_BITS) || sticky;

        // Only the operand with the smaller leading bit can have been placed below the
        // accumulator, so a sticky bit always belongs to the smaller of the two.
        uint64_t* acc = product_acc;
        uint32_t sign = product_sign;
        if (product_sign == sz) {
            wide_add(product_acc, addend_acc);
        } else {
            const int32_t cmp = wide_cmp(product_acc, addend_acc);
            // Equal accumulators mean equal values: a dropped tail puts its operand more than 157
            // bits below the other one, which no accumulator of the larger one can match.
            if (cmp == 0)
                return float128(0, static_cast<uint64_t>(exact_zero_sign()) << 63);

            if (cmp > 0) {
                wide_sub(product_acc, addend_acc);
            } else {
                wide_sub(addend_acc, product_acc);
                acc = addend_acc;
                sign = sz;
            }
            // What the subtrahend lost below the accumulator makes the difference smaller by less
            // than one of its last places. Borrowing that one place and marking the remainder
            // sticky says exactly that.
            if (sticky)
                wide_dec(acc);
        }

        const int32_t msb = wide_msb(acc, WIDE_WORDS);
        int32_t expo = top - WIDE_TOP + msb;

        // Read the 113 bit mantissa out of the accumulator, with the bit below it and whether
        // anything below that was set.
        uint64_t l = 0, h = 0;
        uint64_t guard = 0;
        const int32_t lsb = msb - FRAC_BITS;
        if (lsb <= 0) {
            // Fewer bits survived than the mantissa holds, so the result is exact. Cancellation
            // that deep cannot coexist with a dropped tail, which needs the operands far apart.
            l = acc[0];
            h = acc[1];
            shift_left128_inplace_safe(l, h, -lsb);
        } else {
            const int32_t word = lsb >> 6;
            const int32_t bit = lsb & 63;
            l = acc[word] >> bit;
            h = (word + 1 < WIDE_WORDS) ? (acc[word + 1] >> bit) : 0;
            if (bit != 0) {
                if (word + 1 < WIDE_WORDS)
                    l |= acc[word + 1] << (64 - bit);
                if (word + 2 < WIDE_WORDS)
                    h |= acc[word + 2] << (64 - bit);
            }
            h &= FRAC_UNITY | UPPER_FRAC_MASK;

            guard = FP128_GET_BIT(acc[(lsb - 1) >> 6], (lsb - 1) & 63);
            for (int32_t i = 0; i < ((lsb - 1) >> 6); ++i)
                sticky = sticky || acc[i] != 0;
            const int32_t below = (lsb - 1) & 63;
            if (below != 0)
                sticky = sticky || (acc[(lsb - 1) >> 6] & FP128_MAX_VALUE_64(below)) != 0;
        }

        // The single rounding, to the final width: a subnormal result keeps fewer than 113 bits, and
        // rounding to 113 first would round it twice.
        return round_pack(sign, expo, l, h, (guard << 63) | (sticky ? 1 : 0));
    }
    /**
     * @brief Gets the mantissa and exponent of a floating-point number.
     * @param x Floating-point value.
     * @param expptr Floating-point value.
     * @return Mantissa in the [0.5,1) range.
     */
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr float128 frexp(float128 x, int* expptr) noexcept
    {
        if (x.is_special() || x.is_zero()) {
            *expptr = 0;
            return x;
        }

        uint64_t l, h;
        int32_t e;
        uint32_t s;
        x.get_components(l, h, e, s);
        *expptr = e + 1;
        float128 res(l, h, EXP_BIAS - 1, s);
        return res;
    }
    /**
     * @brief Multiplies a floating-point number by an integral power of two.
     * @param x Floating-point value.
     * @param exp Integer exponent.
     * @return The ldexp functions return the value of x * 2^exp if successful. On overflow, and depending on the sign of x, ldexp returns +/- inf
     */
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr float128 ldexp(float128 x, int exp) noexcept
    {
        if (x.is_zero())
            return x;
        if (x.is_special())
            return x.is_nan() ? propagate_nan(x) : x;

        // Clamped before it is negated, which INT_MIN cannot be.
        if (exp > 0)
            return x << ((exp < SCALE_LIMIT) ? exp : SCALE_LIMIT);
        return x >> ((exp > -SCALE_LIMIT) ? -exp : SCALE_LIMIT);
    }

    /**
     * @name Classification and comparison
     *
     * The `<cmath>` predicates. They are functions rather than the macros the C header defines, in
     * the same way the standard library's own C++ overloads are, so they take part in overload
     * resolution and can be found by argument dependent lookup.
     * @{
     */

    /// @brief True when the sign bit is set, including for -0 and for a negative NaN.
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr bool signbit(const float128& x) noexcept { return x.get_sign() != 0; }
    /// @brief True when x is neither zero, subnormal, infinite nor NaN.
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr bool isnormal(const float128& x) noexcept { return x.is_normal() && !x.is_zero(); }
    /// @brief One of FP_NAN, FP_INFINITE, FP_ZERO, FP_SUBNORMAL or FP_NORMAL.
    [[nodiscard]] friend FP128_INLINE constexpr int fpclassify(const float128& x) noexcept
    {
        if (x.is_nan())
            return FP_NAN;
        if (x.is_inf())
            return FP_INFINITE;
        if (x.is_zero())
            return FP_ZERO;
        return x.is_subnormal() ? FP_SUBNORMAL : FP_NORMAL;
    }
    /// @brief True when either operand is a NaN, so that the two do not compare in any order.
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr bool isunordered(const float128& x, const float128& y) noexcept
    {
        return x.is_nan() || y.is_nan();
    }
    /// @brief x > y, false rather than unordered when either is a NaN.
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr bool isgreater(const float128& x, const float128& y) noexcept
    {
        return !isunordered(x, y) && x > y;
    }
    /// @brief x >= y, false when either is a NaN.
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr bool isgreaterequal(const float128& x, const float128& y) noexcept
    {
        return !isunordered(x, y) && x >= y;
    }
    /// @brief x < y, false when either is a NaN.
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr bool isless(const float128& x, const float128& y) noexcept
    {
        return !isunordered(x, y) && x < y;
    }
    /// @brief x <= y, false when either is a NaN.
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr bool islessequal(const float128& x, const float128& y) noexcept
    {
        return !isunordered(x, y) && x <= y;
    }
    /// @brief x < y or x > y, false when either is a NaN.
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr bool islessgreater(const float128& x, const float128& y) noexcept
    {
        return !isunordered(x, y) && (x < y || x > y);
    }
    /// @}

    /**
     * @name IEEE 754-2008 operations beyond <cmath>
     *
     * Operations clause 5 requires that the C library of the time did not name. They are spelled
     * the way ISO/IEC TS 18661-1, and after it C23, binds them to C.
     * @{
     */

    /**
     * @brief The encoding as a key whose unsigned order is the total order of IEEE 754-2008 5.10.
     *
     * Positive encodings already ascend with their value, NaNs above infinity and quiet NaNs above
     * signaling ones; negative encodings descend. Flipping every bit of a negative encoding and
     * only the sign bit of a positive one puts the whole range in one ascending sequence.
     */
    FP128_FORCE_INLINE constexpr void total_order_key(uint64_t& l, uint64_t& h) const noexcept
    {
        const uint64_t flip = (get_sign() != 0) ? UINT64_MAX : 0;
        l = low ^ flip;
        h = high ^ (flip | SIGN_MASK);
    }

    /**
     * @brief totalOrder(x, y): true when x precedes y, or equals it, in the IEEE 754 total order.
     *
     * The order every value takes part in, NaNs included: -NaN < -inf < negative values < -0 <
     * +0 < positive values < +inf < +NaN, a signaling NaN before a quiet one of the same sign, and
     * NaNs of the same kind ordered by payload. It never signals, even for a signaling NaN.
     *
     * @param x First value
     * @param y Second value
     * @return True when x is ordered before or the same as y.
     */
    [[nodiscard]] friend FP128_INLINE constexpr bool totalorder(const float128& x, const float128& y) noexcept
    {
        uint64_t xl = 0, xh = 0, yl = 0, yh = 0;
        x.total_order_key(xl, xh);
        y.total_order_key(yl, yh);
        return (xh != yh) ? (xh < yh) : (xl <= yl);
    }
    /// @brief totalOrderMag(x, y): totalorder() of the absolute values.
    [[nodiscard]] friend FP128_INLINE constexpr bool totalordermag(const float128& x, const float128& y) noexcept
    {
        return totalorder(fabs(x), fabs(y));
    }

    /**
     * @brief maxNumMag: whichever of x and y has the larger magnitude, or fmax() when they tie.
     * @param x First value
     * @param y Second value
     * @return The larger magnitude. A quiet NaN is treated as missing; a signaling one is the
     *         invalid operation and gives a quiet NaN.
     */
    [[nodiscard]] friend FP128_INLINE constexpr float128 fmaxmag(const float128& x, const float128& y) noexcept
    {
        if (!x.is_nan() && !y.is_nan()) {
            const float128 ax = fabs(x), ay = fabs(y);
            if (ax > ay)
                return x;
            if (ay > ax)
                return y;
        }
        return fmax(x, y);
    }
    /// @brief minNumMag: whichever of x and y has the smaller magnitude, or fmin() when they tie.
    /// @copydetails fmaxmag
    [[nodiscard]] friend FP128_INLINE constexpr float128 fminmag(const float128& x, const float128& y) noexcept
    {
        if (!x.is_nan() && !y.is_nan()) {
            const float128 ax = fabs(x), ay = fabs(y);
            if (ax < ay)
                return x;
            if (ay < ax)
                return y;
        }
        return fmin(x, y);
    }

    /// @brief isCanonical: always true, every binary128 encoding is canonical.
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr bool iscanonical(const float128&) noexcept { return true; }
    /// @brief isSignaling: true for a signaling NaN.
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr bool issignaling(const float128& x) noexcept { return x.is_signaling(); }
    /// @brief isSubnormal: true for a subnormal value, which zero is not.
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr bool issubnormal(const float128& x) noexcept { return x.is_subnormal() && !x.is_zero(); }
    /// @brief isZero: true for either zero.
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr bool iszero(const float128& x) noexcept { return x.is_zero(); }
    /// @}

    /**
     * @brief A quiet NaN carrying the payload spelled out in the argument.
     *
     * Mirrors the C library's nan(): the string is read as an unsigned integer, decimal by
     * default, hexadecimal after an 0x prefix, and its low bits become the NaN's fraction. An
     * empty or unparsable string gives the default NaN.
     *
     * @param payload Digits of the payload, may be empty
     * @return A quiet NaN.
     */
    [[nodiscard]] friend FP128_INLINE float128 nan(const char* payload) noexcept
    {
        uint64_t bits = 1;
        if (payload != nullptr && *payload != '\0') {
            const uint64_t parsed = strtoull(payload, nullptr, 0);
            // A zero fraction is an infinity rather than a NaN, so the default payload stands in.
            if (parsed != 0)
                bits = parsed;
        }
        // The leading fraction bit marks the NaN quiet; the payload fills the bits below it.
        float128 res = float128(bits, 0, INF_EXP_BIASED, 0);
        res.high |= QUIET_NAN_BIT;
        return res;
    }

    /// @brief True when x is an integer whose last place is one, which is the parity a tie needs.
    [[nodiscard]] FP128_INLINE static constexpr bool is_odd_int(const float128& x) noexcept
    {
        uint64_t l = 0, h = 0;
        int32_t expo = 0;
        uint32_t sign = 0;
        x.get_components(l, h, expo, sign);
        // Below one there is no units bit; above 2^112 the last place is two or more, so every
        // value that far out is even.
        if (expo < 0 || expo > FRAC_BITS)
            return false;
        const int32_t bit = FRAC_BITS - expo;
        return (bit < 64) ? ((l >> bit) & 1) != 0 : ((h >> (bit - 64)) & 1) != 0;
    }

    /**
     * @name Rounding to the current mode
     *
     * Both round in the current rounding direction, which is to nearest with ties to even unless
     * FP128_IEEE_ENV is defined and fesetround() chose another. They differ in one respect:
     * rint is IEEE 754's roundToIntegralExact and raises the inexact exception when the result
     * differs from its argument, nearbyint never does. round() differs from both: it breaks a tie
     * away from zero whatever the rounding direction.
     * @{
     */
    /// @brief x rounded to an integral value in the current direction, raising inexact when that changes it.
    [[nodiscard]] friend FP128_INLINE constexpr float128 rint(const float128& x) noexcept
    {
        return round_integral(x, detail::current_rounding(), false, true);
    }
    /// @brief x rounded to an integral value in the current direction.
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr float128 nearbyint(const float128& x) noexcept
    {
        return round_integral(x, detail::current_rounding(), false, false);
    }
    /// @}

    /**
     * @name Exponent manipulation
     * @{
     */
    /// @brief x * 2^n, the same operation ldexp performs.
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr float128 scalbn(const float128& x, int n) noexcept { return ldexp(x, n); }
    /// @brief x * 2^n for a long exponent, clamped to what the format can reach.
    [[nodiscard]] friend FP128_INLINE constexpr float128 scalbln(const float128& x, long n) noexcept
    {
        // Anything past the format's range saturates the same way the shift itself would, and
        // clamping keeps the conversion to int well defined.
        constexpr long limit = SCALE_LIMIT;
        const long clamped = (n > limit) ? limit : ((n < -limit) ? -limit : n);
        return ldexp(x, static_cast<int>(clamped));
    }
    /// @}

    /**
     * @brief The representable value next to x in the direction of y.
     * @param x Starting value
     * @param y Value giving the direction
     * @return The neighbour of x towards y, or y itself when the two are equal.
     */
    [[nodiscard]] friend FP128_INLINE constexpr float128 nextafter(const float128& x, const float128& y) noexcept
    {
        if (x.is_nan() || y.is_nan())
            return propagate_nan(x, y);
        if (x == y)
            return y;  // the sign of y is what the standard hands back here
        return (y > x) ? nextUp(x) : nextDown(x);
    }
    /// @brief The neighbour of x towards y. The same as nextafter: there is no wider type to convert from.
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr float128 nexttoward(const float128& x, const float128& y) noexcept
    {
        return nextafter(x, y);
    }

    /**
     * @brief |x| mod |y|, along with the low bits of the quotient it implies.
     *
     * The shared engine behind fmod, remainder and remquo. All three are exact operations on the
     * mantissas, and all three need the same long division; only what they do with the quotient
     * differs.
     *
     * @param x Dividend, must be finite and non zero
     * @param y Divisor, must be finite and non zero
     * @param quotient Receives the low bits of the integer quotient
     * @return The remainder, with a positive sign.
     */
    [[nodiscard]] FP128_INLINE static float128 fmod_magnitude(const float128& x, const float128& y, uint64_t& quotient) noexcept
    {
        quotient = 0;

        uint64_t lx = 0, hx = 0, ly = 0, hy = 0;
        int32_t ex = 0, ey = 0;
        uint32_t sx = 0, sy = 0;
        x.get_components(lx, hx, ex, sx);
        y.get_components(ly, hy, ey, sy);

        if (ex < ey || (ex == ey && uint128_t(lx, hx) < uint128_t(ly, hy)))
            return fabs(x);

        uint128_t rem(lx, hx);
        const uint128_t mod(ly, hy);
        uint128_t block_quotient = rem / mod;
        rem -= block_quotient * mod;
        uint64_t low = 0, high = 0;
        block_quotient.get_components(low, high);
        quotient = low;

        int32_t shift = ex - ey;
        constexpr int32_t block = 15;
        while (shift > 0) {
            const int32_t step = (shift < block) ? shift : block;
            rem <<= step;
            shift -= step;
            block_quotient = rem / mod;
            rem -= block_quotient * mod;
            block_quotient.get_components(low, high);
            // Only the low bits of the quotient are ever wanted, so the overflow is dropped.
            quotient = (quotient << step) + low;
        }

        if (rem.is_zero())
            return float128();

        uint64_t rl = 0, rh = 0;
        rem.get_components(rl, rh);
        const int32_t msb = static_cast<int32_t>(log2(rl, rh));
        shift_left128_inplace_safe(rl, rh, FRAC_BITS - msb);

        float128 res;
        res.set_components(rl, rh, ey - FRAC_BITS + msb, 0);
        return res;
    }

    /**
     * @brief IEEE 754 remainder: x - n*y with n the quotient rounded to the nearest even integer.
     *
     * Unlike fmod, whose result keeps the sign of x, this one lands in [-|y|/2, |y|/2].
     *
     * @param x Dividend
     * @param y Divisor
     * @return The remainder.
     */
    [[nodiscard]] friend FP128_INLINE float128 remainder(const float128& x, const float128& y) noexcept
    {
        int quotient = 0;
        return remquo(x, y, &quotient);
    }

    /**
     * @brief The IEEE remainder together with the low bits of the quotient.
     * @param x Dividend
     * @param y Divisor
     * @param quo Receives at least the low three bits of the quotient, with the sign of x/y
     * @return The remainder.
     */
    [[nodiscard]] friend FP128_INLINE float128 remquo(const float128& x, const float128& y, int* quo) noexcept
    {
        if (quo != nullptr)
            *quo = 0;
        if (x.is_nan() || y.is_nan())
            return propagate_nan(x, y);
        if (x.is_inf() || y.is_zero())
            return invalid_operation();
        if (y.is_inf() || x.is_zero())
            return x;

        uint64_t bits = 0;
        const float128 magnitude = fabs(y);
        float128 rem = fmod_magnitude(x, y, bits);

        // fmod leaves the remainder in [0, |y|); the IEEE one is the nearer of that and
        // |y| - it, with a tie going to whichever leaves an even quotient.
        const float128 twice = rem << 1;
        const bool round_up = (twice > magnitude) || (twice == magnitude && (bits & 1) != 0);
        if (round_up) {
            rem -= magnitude;
            ++bits;
        }

        if (quo != nullptr) {
            // The standard asks for the low bits of the magnitude of the quotient, signed by the
            // signs of the operands. Seven bits is the customary amount and fits any int.
            const int low_bits = static_cast<int>(bits & 0x7F);
            *quo = (x.get_sign() != y.get_sign()) ? -low_bits : low_bits;
        }

        rem.set_sign(rem.is_zero() ? x.get_sign() : (x.get_sign() ^ rem.get_sign()));
        return rem;
    }

    /**
     * @brief Length of the diagonal of a box, computed without an intermediate overflow.
     * @param x First side
     * @param y Second side
     * @param z Third side
     * @return sqrt(x*x + y*y + z*z)
     */
    [[nodiscard]] friend FP128_INLINE float128 hypot(const float128& x, const float128& y, const float128& z) noexcept
    {
        if (x.is_signaling() || y.is_signaling() || z.is_signaling())
            return propagate_nan(x, y, z);
        if (x.is_inf() || y.is_inf() || z.is_inf())
            return inf();
        if (x.is_nan() || y.is_nan() || z.is_nan())
            return propagate_nan(x, y, z);

        // Scaling by the largest term keeps every square inside the format's range, which squaring
        // the values as they came would not: a side above 2^8192 overflows on its own.
        const float128 largest = fmax(fabs(x), fmax(fabs(y), fabs(z)));
        if (largest.is_zero())
            return largest;
        // From here on only the flags the result justifies reach the caller.
        detail::flag_filter filter;

        const int32_t expo = ilogb(largest);
        const float128 a = ldexp(x, -expo);
        const float128 b = ldexp(y, -expo);
        const float128 c = ldexp(z, -expo);
        return filter(ldexp(sqrt(sqr(a) + sqr(b) + sqr(c)), expo));
    }

    /// @brief Terms of the Stirling series lgamma() runs. Sixteen reach 2^-119 for an argument of 40.
    static constexpr int32_t STIRLING_TERMS = 16;

    /**
     * @brief Natural logarithm of the absolute value of the gamma function.
     *
     * Stirling's asymptotic series converges usefully only for a large argument, so a small one is
     * walked up to 40 with the recurrence gamma(x+1) = x*gamma(x) and the logarithms of the
     * factors are subtracted off afterwards. A negative argument goes through the reflection
     * formula, which is why the result is the logarithm of the absolute value: the gamma function
     * alternates sign between the poles at the non positive integers.
     *
     * @param x Input value
     * @return log(|gamma(x)|)
     */
    [[nodiscard]] friend FP128_INLINE float128 lgamma(float128 x) noexcept
    {
        if (x.is_nan())
            return propagate_nan(x);
        if (x.is_inf())
            return inf();
        // The non positive integers are the poles of the gamma function
        if (x.is_zero() || (x.is_negative() && x.is_int())) {
            detail::raise_flags(FE_DIVBYZERO);
            return inf();
        }
        // gamma(1) and gamma(2) are both one. The recurrence and the series below would reach the
        // logarithm of 31! and subtract it from itself, which cancels to a small non zero value
        // rather than to the exact answer.
        if (x == float128::one() || x == float128(2))
            return float128();
        // From here on only the flags the result justifies reach the caller.
        detail::flag_filter filter;

        constexpr float128 coefficients[STIRLING_TERMS] = {
            float128(0x5555555555555555, 0x555555555555, 0x3FFB, 0),  // B2 / (2*1)
            float128(0x6C16C16C16C16C17, 0x6C16C16C16C1, 0x3FF6, 1),  // B4 / (4*3)
            float128(0xA01A01A01A01A01A, 0xA01A01A01A01, 0x3FF4, 0),  // B6 / (6*5)
            float128(0x3813813813813814, 0x381381381381, 0x3FF4, 1),  // B8 / (8*7)
            float128(0x3570EA73806E5479, 0xB951E2B18FF2, 0x3FF4, 0),  // B10 / (10*9)
            float128(0xC81F6AB0D9993C7D, 0xF6AB0D9993C7, 0x3FF5, 1),  // B12 / (12*11)
            float128(0xA41A41A41A41A41A, 0xA41A41A41A41, 0x3FF7, 0),  // B14 / (14*13)
            float128(0x7DC2064A8ED3175C, 0xE4286CB0F539, 0x3FF9, 1),  // B16 / (16*15)
            float128(0xFFA1876FE96381E0, 0x6FE96381E067, 0x3FFC, 0),  // B18 / (18*17)
            float128(0x9EDBDB9CE625987D, 0x6476701181F3, 0x3FFF, 1),  // B20 / (20*19)
            float128(0x5A74F53910C8B380, 0xACE44322CE00, 0x4002, 0),  // B22 / (22*21)
            float128(0xAAB67EE25D73C0F9, 0x39B2525CCCC1, 0x4006, 1),  // B24 / (24*23)
            float128(0x1B4E81B4E81B4E82, 0x12234E81B4E8, 0x400A, 0),  // B26 / (26*25)
            float128(0x7EB3FEDDD8496920, 0x1A198AE1C4AB, 0x400E, 1),  // B28 / (28*27)
            float128(0xA38433DC9FB888D4, 0x51A2089A6E11, 0x4012, 0),  // B30 / (30*29)
            float128(0x77880C2D3577880C, 0xD1089B142D35, 0x4016, 1),  // B32 / (32*31)
        };
        constexpr float128 half_log_two_pi(0x4A69297920028832, 0xD67F1C864BEB, 0x3FFE, 0);
        constexpr float128 stirling_limit(0, 0, EXP_BIAS + 5, 0);  // 32, where sixteen terms suffice

        // gamma(x) * gamma(1-x) = pi / sin(pi*x) reflects the negative half onto the positive one,
        // where the series lives.
        if (x.is_negative()) {
            const float128 reflected = sin(float128::pi() * x);
            return filter(log(float128::pi() / fabs(reflected)) - lgamma(float128::one() - x));
        }

        // Walk up to where the asymptotic series is accurate, remembering what to divide out.
        float128 scale_log;
        while (x < stirling_limit) {
            scale_log += log(x);
            x += float128::one();
        }

        const float128 inv_xx = float128::one() / sqr(x);
        float128 series;
        for (int32_t i = STIRLING_TERMS - 1; i >= 1; --i)
            series = (series + coefficients[i]) * inv_xx;
        series = (series + coefficients[0]) / x;

        return filter((x - float128::half()) * log(x) - x + half_log_two_pi + series - scale_log);
    }

    /**
     * @brief The gamma function.
     *
     * Exact for the small integer arguments, where the answer is a factorial the format holds
     * without rounding, and exp(lgamma) elsewhere. The exponential is what limits the accuracy: a
     * logarithm as large as 7000 carries its own last bit into the seventh place of the result.
     *
     * @param x Input value
     * @return gamma(x)
     */
    [[nodiscard]] friend FP128_INLINE float128 tgamma(float128 x) noexcept
    {
        if (x.is_nan())
            return propagate_nan(x);
        if (x.is_inf())
            return x.is_negative() ? invalid_operation() : x;
        // The poles, and the two zeros of 1/gamma that a signed zero argument picks out
        if (x.is_zero()) {
            detail::raise_flags(FE_DIVBYZERO);
            return x.is_negative() ? -inf() : inf();
        }
        if (x.is_negative() && x.is_int())
            return invalid_operation();
        // From here on only the flags the result justifies reach the caller.
        detail::flag_filter filter;
        // gamma(1756) overflows
        if (x > float128(1756))
            return filter.inexact(inf());

        // 34! is the largest factorial a binary128 holds exactly, so every integer argument up to
        // 35 comes back with no rounding at all.
        if (x.is_int() && x <= float128(35)) {
            float128 result = float128::one();
            for (float128 i = float128(2); i < x; i += float128::one())
                result *= i;
            return filter(result);
        }

        const float128 magnitude = exp(lgamma(x));
        if (x.is_positive())
            return filter(magnitude);

        // Between two poles the gamma function keeps one sign, alternating with every step: it is
        // negative on (-1, 0), positive on (-2, -1), and so on, which is the parity of floor(x).
        const float128 floor_x = floor(x);
        return filter(is_odd_int(floor_x) ? -magnitude : magnitude);
    }

    /// @brief Absolute value, the name `<cmath>` gives the floating point overload alongside fabs.
    [[nodiscard]] friend FP128_FORCE_INLINE constexpr float128 abs(const float128& x) noexcept { return fabs(x); }

    /// @brief User-defined literal implementation for constructing float128 from a string.
    friend float128 operator""_f128(const char* literal) { return float128(literal); }
};

static_assert(sizeof(float128) == sizeof(uint64_t) * 2);

/***********************************************************************************
 *                        Text conversion and standard library integration
 ************************************************************************************/

namespace detail
{
/// @brief Default precision the e, f and g styles use when none is given, as in printf.
inline constexpr int32_t DEFAULT_PRECISION = 6;

/// @brief Appends the decimal exponent of a scientific form, with a sign and at least two digits.
inline void append_exponent(std::string& out, int32_t exponent, char marker)
{
    out += marker;
    out += (exponent < 0) ? '-' : '+';
    // The magnitude is taken in the unsigned domain, where the negation wraps. Negating the signed
    // value is undefined for the most negative one and produces the same bit pattern for the rest.
    uint32_t magnitude = (exponent < 0) ? (0u - static_cast<uint32_t>(exponent)) : static_cast<uint32_t>(exponent);

    // Ten digits is the widest a uint32_t can be. No exponent this function is handed comes close -
    // the largest is 16494, from the hexadecimal form of the smallest subnormal - but the buffer is
    // sized for the type rather than for the callers, so the loops below cannot run past it.
    //
    // The digits are produced least significant first and so are written from the end of the
    // buffer backwards, which leaves them already in order and lets the whole run be appended at
    // once. Filling forwards and reversing afterwards needs three loops over one shared counter,
    // and MSVC's code analysis cannot see that the counter stays in range across all of them.
    constexpr int32_t max_digits = 10;
    char digits[max_digits];
    int32_t index = max_digits;
    while (index > 0) {
        digits[--index] = static_cast<char>('0' + (magnitude % 10));
        magnitude /= 10;
        if (magnitude == 0)
            break;
    }

    // The exponent always shows at least two digits, as it does for a double.
    while (index > max_digits - 2)
        digits[--index] = '0';

    out.append(digits + index, static_cast<size_t>(max_digits - index));
}

/**
 * @brief Fewest significant digits that read back as the same value.
 *
 * The search is a bisection rather than a scan because the property is monotone: a decimal with
 * more digits is at least as close to the value as one with fewer, so if some length reads back
 * correctly then every longer one does too.
 *
 * @param value Finite, non zero value
 * @return Digit count in [1, 36].
 */
[[nodiscard]] inline int32_t shortest_digit_count(const float128& value)
{
    uint64_t low = 0, high = 0;
    int32_t exponent = 0;
    uint32_t sign = 0;
    value.get_components(low, high, exponent, sign);
    uint64_t want_low = 0, want_high = 0;
    value.get_bits(want_low, want_high);
    want_high &= ~(1ull << 63);

    // The trial conversions raise exceptions that belong to none of the caller's operations.
    const flag_quiet quiet;
    char digits[40];
    int32_t lower = 1;
    int32_t upper = 36;  // max_digits10, which always reads back

    while (lower < upper) {
        const int32_t middle = (lower + upper) / 2;
        const int32_t exponent10 = to_decimal_digits(low, high, exponent, middle, digits);

        // Read back the way from_chars() reads, rounded once to the final width, which for a
        // subnormal is fewer than 113 bits: comparing 113 bit mantissas instead demanded more
        // digits than a subnormal needs.
        uint64_t back_low = 0, back_high = 0, extra = 0;
        int32_t back_exponent = 0;
        bool same = false;
        if (from_decimal_digits(digits, middle, exponent10 - middle, back_low, back_high, back_exponent, &extra)) {
            const float128 back = float128::round_pack(0, back_exponent, back_low, back_high, extra, rounding::nearest_even);
            uint64_t got_low = 0, got_high = 0;
            back.get_bits(got_low, got_high);
            same = got_low == want_low && got_high == want_high;
        }
        if (same)
            upper = middle;
        else
            lower = middle + 1;
    }
    return lower;
}

/**
 * @brief Renders a float128 as text.
 *
 * The one place the character form of a value is produced. std::format, the stream inserter,
 * to_chars and the conversion to std::string all come through here so that they cannot disagree
 * with each other.
 *
 * @param value Value to render
 * @param style One of 'a', 'e', 'f' or 'g', or '\0' for the shortest form that reads back exactly
 * @param precision Digits after the point, or significant digits for 'g'. Negative selects the
 *        default for the style
 * @param uppercase Emit the digits, the exponent marker and the special values in upper case
 * @param alternate Keep the decimal point and, for 'g', the trailing zeros
 * @param sign_char Character to emit for a non negative value: '-' means emit nothing
 * @return The rendered text.
 */
[[nodiscard]] inline std::string render(const float128& value, char style, int32_t precision, bool uppercase, bool alternate, char sign_char)
{
    std::string out;
    if (value.is_negative())
        out += '-';
    else if (sign_char != '-')
        out += sign_char;

    if (value.is_nan()) {
        out += uppercase ? "NAN" : "nan";
        return out;
    }
    if (value.is_inf()) {
        out += uppercase ? "INF" : "inf";
        return out;
    }

    const char* const hex_digits = uppercase ? "0123456789ABCDEF" : "0123456789abcdef";

    // The hexadecimal form is a direct reading of the encoding, so it needs none of the decimal
    // machinery and is exact by construction at any precision.
    if (style == 'a') {
        out += uppercase ? "0X" : "0x";
        if (value.is_zero()) {
            out += '0';
            if (alternate || precision > 0) {
                out += '.';
                for (int32_t i = 0; i < precision; ++i)
                    out += '0';
            }
            append_exponent(out, 0, uppercase ? 'P' : 'p');
            return out;
        }

        uint64_t low = 0, high = 0;
        int32_t exponent = 0;
        uint32_t sign = 0;
        value.get_components(low, high, exponent, sign);

        // 28 hex digits hold the whole fraction, and the leading one sits above the point.
        char fraction[28];
        for (int32_t i = 0; i < 28; ++i) {
            const int32_t shift = 108 - 4 * i;
            const uint64_t nibble = (shift >= 64) ? ((high >> (shift - 64)) & 0xF) : ((low >> shift) & 0xF);
            fraction[i] = hex_digits[nibble];
        }

        const auto hex_value = [](char c) {
            if (c >= '0' && c <= '9')
                return c - '0';
            return ((c >= 'a') ? (c - 'a') : (c - 'A')) + 10;
        };

        int32_t kept = (precision >= 0) ? precision : 28;
        if (kept > 28)
            kept = 28;
        // Rounding a shortened fraction, in the current direction: to nearest the first dropped
        // digit decides and a tie goes to the even digit (it used to go up). A carry can run
        // through the fraction and into the leading one, which is the next power of two.
        bool leading_two = false;
        if (precision >= 0 && precision < 28) {
            const int32_t first_dropped = hex_value(fraction[precision]);
            bool rest = false;
            for (int32_t i = precision + 1; i < 28; ++i)
                rest = rest || fraction[i] != '0';
            const int32_t last_kept = (precision > 0) ? hex_value(fraction[precision - 1]) : 1;
            bool up = false;
            switch (digit_rounding_for(sign)) {
            case digit_rounding::nearest_even: up = first_dropped > 8 || (first_dropped == 8 && (rest || (last_kept & 1) != 0)); break;
            case digit_rounding::away_from_zero: up = first_dropped != 0 || rest; break;
            default: break;
            }
            if (first_dropped != 0 || rest)
                raise_flags(FE_INEXACT);
            if (up) {
                int32_t index = precision - 1;
                while (index >= 0) {
                    const int32_t digit = hex_value(fraction[index]);
                    if (digit != 15) {
                        fraction[index] = hex_digits[digit + 1];
                        break;
                    }
                    fraction[index] = '0';
                    --index;
                }
                if (index < 0)
                    leading_two = true;
            }
        }
        if (precision < 0) {
            // Trailing zeros carry no information in the shortest form.
            while (kept > 0 && fraction[kept - 1] == '0')
                --kept;
        }

        out += leading_two ? '2' : '1';
        if (kept > 0 || alternate) {
            out += '.';
            out.append(fraction, static_cast<size_t>(kept));
        }
        append_exponent(out, exponent, uppercase ? 'P' : 'p');
        return out;
    }

    if (value.is_zero()) {
        out += '0';
        // The fixed and scientific layouts always show the digits asked for; the general one shows
        // none unless the alternate form wanted them.
        int32_t zeros = 0;
        if (style == 'f' || style == 'e')
            zeros = (precision >= 0) ? precision : DEFAULT_PRECISION;
        else if (alternate)
            zeros = (precision > 0) ? precision - 1 : DEFAULT_PRECISION - 1;

        if (zeros > 0 || alternate) {
            out += '.';
            out.append(static_cast<size_t>(zeros), '0');
        }
        if (style == 'e')
            append_exponent(out, 0, uppercase ? 'E' : 'e');
        return out;
    }

    uint64_t low = 0, high = 0;
    int32_t exponent = 0;
    uint32_t sign = 0;
    value.get_components(low, high, exponent, sign);

    char digits[MAX_OUTPUT_DIGITS];
    int32_t significant = 0;
    int32_t exponent10 = 0;
    // The digits are rounded in the current direction, except for the shortest form, which is
    // defined by reading back to nearest. Either way a conversion that drops something is inexact.
    const digit_rounding direction = digit_rounding_for(sign);
    bool exact = false;

    if (style == '\0') {
        significant = shortest_digit_count(value);
        exponent10 = to_decimal_digits(low, high, exponent, significant, digits, &exact);
    } else if (style == 'e') {
        significant = ((precision >= 0) ? precision : DEFAULT_PRECISION) + 1;
        if (significant > MAX_OUTPUT_DIGITS)
            significant = MAX_OUTPUT_DIGITS;
        exponent10 = to_decimal_digits(low, high, exponent, significant, digits, &exact, direction);
    } else if (style == 'g') {
        significant = (precision > 0) ? precision : ((precision == 0) ? 1 : DEFAULT_PRECISION);
        if (significant > MAX_OUTPUT_DIGITS)
            significant = MAX_OUTPUT_DIGITS;
        exponent10 = to_decimal_digits(low, high, exponent, significant, digits, &exact, direction);
    } else {
        // Fixed notation asks for a count of digits after the point rather than significant ones,
        // and how many that is depends on where the value sits. One digit is generated first to
        // find that out, then the real request is made.
        const int32_t after_point = (precision >= 0) ? precision : DEFAULT_PRECISION;
        int32_t probe_exponent10 = to_decimal_digits(low, high, exponent, 1, digits, &exact);
        significant = probe_exponent10 + after_point;

        // A single digit probe rounds, and a value just below a power of ten rounds up into the
        // next one: 99 comes back as 1e2 and claims one more digit than it has. Generating the
        // digits reveals the true exponent, and one retry with it is always enough because the
        // rounding can only move the exponent by one.
        if (significant > 0 && significant <= MAX_OUTPUT_DIGITS) {
            const int32_t actual_exponent10 = to_decimal_digits(low, high, exponent, significant, digits, nullptr, direction);
            if (actual_exponent10 != probe_exponent10) {
                probe_exponent10 = actual_exponent10;
                significant = actual_exponent10 + after_point;
            }
        }

        if (significant <= 0) {
            // Every digit asked for is zero, and the value lies below the last of them. To nearest,
            // it rounds up only from above half a place - a tie goes to the even digit, the zero
            // already there; away from zero it always does.
            raise_flags(FE_INEXACT);
            bool round_up = false;
            if (direction == digit_rounding::away_from_zero)
                round_up = true;
            else if (direction == digit_rounding::nearest_even)
                round_up = (significant == 0) && (digits[0] > '5' || (digits[0] == '5' && !exact));
            out += '0';
            if (after_point > 0 || alternate) {
                out += '.';
                for (int32_t i = 0; i < after_point; ++i)
                    out += (round_up && i == after_point - 1) ? '1' : '0';
            }
            return out;
        }

        if (significant > MAX_OUTPUT_DIGITS) {
            // Past MAX_OUTPUT_DIGITS the exact expansion is no longer produced and the rest of
            // the request is filled with zeros. The expansion of a subnormal runs to 16494 places
            // after the point; the cap is set at what writing the widest finite value out in full
            // needs, which is 4933.
            exponent10 = to_decimal_digits(low, high, exponent, MAX_OUTPUT_DIGITS, digits, &exact, direction);
            if (!exact)
                raise_flags(FE_INEXACT);
            const int32_t padding = significant - MAX_OUTPUT_DIGITS;
            significant = MAX_OUTPUT_DIGITS;
            std::string body(digits, static_cast<size_t>(significant));
            body.append(static_cast<size_t>(padding), '0');
            // Lay the point out from the exponent the generator returned.
            if (exponent10 <= 0) {
                out += "0.";
                out.append(static_cast<size_t>(-exponent10), '0');
                out += body;
            } else {
                out.append(body, 0, static_cast<size_t>(exponent10));
                out += '.';
                out.append(body, static_cast<size_t>(exponent10), std::string::npos);
            }
            return out;
        }

        exponent10 = to_decimal_digits(low, high, exponent, significant, digits, &exact, direction);
    }
    if (!exact)
        raise_flags(FE_INEXACT);

    const int32_t scientific_exponent = exponent10 - 1;

    int32_t kept = significant;
    if ((style == 'g' || style == '\0') && !alternate) {
        while (kept > 1 && digits[kept - 1] == '0')
            --kept;
    }

    bool scientific = (style == 'e');
    if (style == 'g') {
        // printf's rule for %g, which std::format keeps for the g type.
        scientific = (scientific_exponent < -4) || (scientific_exponent >= significant);
    } else if (style == '\0') {
        // The default is whichever layout is shorter, with a tie going to the fixed one, which is
        // what std::to_chars produces and therefore what std::format gives a double. The rule for
        // %g would print 100 as 1e+02: it compares the exponent against the digit count, and the
        // shortest form of a round number has very few digits.
        int32_t fixed_length = 0;
        if (exponent10 <= 0)
            fixed_length = 2 - exponent10 + kept;
        else if (exponent10 >= kept)
            fixed_length = exponent10;
        else
            fixed_length = kept + 1;

        int32_t exponent_digits = 2;
        for (int32_t magnitude = (scientific_exponent < 0) ? -scientific_exponent : scientific_exponent; magnitude >= 100; magnitude /= 10)
            ++exponent_digits;
        const int32_t scientific_length = kept + ((kept > 1) ? 1 : 0) + 2 + exponent_digits;

        scientific = scientific_length < fixed_length;
    }

    if (scientific) {
        out += digits[0];
        if (kept > 1 || alternate) {
            out += '.';
            out.append(digits + 1, static_cast<size_t>(kept - 1));
        }
        append_exponent(out, scientific_exponent, uppercase ? 'E' : 'e');
        return out;
    }

    if (exponent10 <= 0) {
        out += "0.";
        out.append(static_cast<size_t>(-exponent10), '0');
        out.append(digits, static_cast<size_t>(kept));
    } else if (exponent10 >= kept) {
        out.append(digits, static_cast<size_t>(kept));
        out.append(static_cast<size_t>(exponent10 - kept), '0');
        if (alternate)
            out += '.';
    } else {
        out.append(digits, static_cast<size_t>(exponent10));
        out += '.';
        out.append(digits + exponent10, static_cast<size_t>(kept - exponent10));
    }

    // Fixed notation pads out to the requested number of digits after the point.
    if (style == 'f') {
        const int32_t after_point = (precision >= 0) ? precision : DEFAULT_PRECISION;
        const size_t point = out.find('.');
        const size_t written = (point == std::string::npos) ? 0 : (out.size() - point - 1);
        if (after_point > 0) {
            if (point == std::string::npos)
                out += '.';
            for (size_t i = written; i < static_cast<size_t>(after_point); ++i)
                out += '0';
        }
    }
    return out;
}
}  // namespace detail

/**
 * @brief The shortest decimal string that reads back as this exact value.
 * @param value Value to convert
 * @return Decimal text, in scientific notation when the exponent is far from zero.
 */
[[nodiscard]] inline std::string to_string(const float128& value)
{
    return detail::render(value, '\0', -1, false, false, '-');
}

namespace detail
{
/**
 * @brief Reads one of the special values: inf, infinity, nan, nan(n-char-sequence) or snan.
 *
 * Case is ignored, as IEEE 754-2008 5.12.1 asks. A payload in parentheses is read the way
 * fp128::nan() reads one, and the parentheses are only consumed when they are closed. "snan" reads
 * as a signaling NaN, which the standard recommends and strtod() does not offer.
 *
 * @param cursor Start of the text, past any sign
 * @param last End of the text
 * @param negative The sign that came before it
 * @param value Receives the value
 * @return Past what was read, or nullptr when the text is none of these.
 */
inline const char* parse_special_value(const char* cursor, const char* last, bool negative, float128& value)
{
    const auto lower = [](char c) { return static_cast<char>((c >= 'A' && c <= 'Z') ? (c - 'A' + 'a') : c); };
    const auto matches = [&](const char* word, size_t length) {
        if (static_cast<size_t>(last - cursor) < length)
            return false;
        for (size_t i = 0; i < length; ++i) {
            if (lower(cursor[i]) != word[i])
                return false;
        }
        return true;
    };

    if (matches("inf", 3)) {
        cursor += 3;
        if (matches("inity", 5))
            cursor += 5;
        value = negative ? -float128::inf() : float128::inf();
        return cursor;
    }
    if (matches("snan", 4)) {
        value = float128::signaling_nan();
        value.set_sign(negative ? 1 : 0);
        return cursor + 4;
    }
    if (matches("nan", 3)) {
        cursor += 3;
        value = float128::nan();
        if (cursor < last && *cursor == '(') {
            const char* close = cursor + 1;
            while (close < last && (isalnum(static_cast<unsigned char>(*close)) || *close == '_'))
                ++close;
            if (close < last && *close == ')') {
                value = fp128::nan(std::string(cursor + 1, close).c_str());
                cursor = close + 1;
            }
        }
        value.set_sign(negative ? 1 : 0);
        return cursor;
    }
    return nullptr;
}

/**
 * @brief Reads a hexadecimal significand with an optional binary exponent, such as 1.8p1.
 *
 * The form printf's %a writes, without its 0x prefix, as std::from_chars reads it for
 * chars_format::hex. IEEE 754-2008 5.12.3 requires it, correctly rounded. Every digit is four bits
 * of the significand, so the first 32 significant digits are collected exactly and any further
 * non zero digit only counts towards the sticky bit; the value is then rounded once.
 *
 * @param first Start of the text
 * @param last One past its end
 * @param value Receives the parsed value, untouched when nothing was parsed
 * @return As from_chars().
 */
inline std::from_chars_result from_hex_chars(const char* first, const char* last, float128& value)
{
    std::from_chars_result result {first, std::errc {}};
    const char* cursor = first;
    if (cursor == last) {
        result.ec = std::errc::invalid_argument;
        return result;
    }

    bool negative = false;
    if (*cursor == '-' || *cursor == '+') {
        negative = (*cursor == '-');
        ++cursor;
    }
    if (const char* end = parse_special_value(cursor, last, negative, value)) {
        result.ptr = end;
        return result;
    }

    const auto hex_digit = [](char c) -> int32_t {
        if (c >= '0' && c <= '9')
            return c - '0';
        if (c >= 'a' && c <= 'f')
            return c - 'a' + 10;
        if (c >= 'A' && c <= 'F')
            return c - 'A' + 10;
        return -1;
    };

    // The value is (h:l) * 2^binary_exponent, plus something below that when sticky is set.
    uint64_t l = 0, h = 0;
    int32_t bits = 0;
    int64_t binary_exponent = 0;
    bool sticky = false, any_digit = false, nonzero = false;
    const auto take = [&](int32_t digit, bool after_point) {
        any_digit = true;
        if (!nonzero && digit == 0) {
            // a leading zero: nothing to keep, but after the point it still moves the exponent
            if (after_point)
                binary_exponent -= 4;
            return;
        }
        nonzero = true;
        if (bits <= 124) {
            h = (h << 4) | (l >> 60);
            l = (l << 4) | static_cast<uint64_t>(digit);
            bits += 4;
            if (after_point)
                binary_exponent -= 4;
        } else {
            sticky = sticky || digit != 0;
            if (!after_point)
                binary_exponent += 4;
        }
    };

    while (cursor < last && hex_digit(*cursor) >= 0)
        take(hex_digit(*cursor++), false);
    if (cursor < last && *cursor == '.') {
        ++cursor;
        while (cursor < last && hex_digit(*cursor) >= 0)
            take(hex_digit(*cursor++), true);
    }
    if (!any_digit) {
        result.ec = std::errc::invalid_argument;
        return result;
    }
    result.ptr = cursor;

    // The binary exponent is only consumed when it is well formed, so that "1p" reads as 1.
    if (cursor < last && (*cursor == 'p' || *cursor == 'P')) {
        const char* probe = cursor + 1;
        bool exponent_negative = false;
        if (probe < last && (*probe == '-' || *probe == '+')) {
            exponent_negative = (*probe == '-');
            ++probe;
        }
        if (probe < last && *probe >= '0' && *probe <= '9') {
            int64_t magnitude = 0;
            while (probe < last && *probe >= '0' && *probe <= '9') {
                if (magnitude < 1000000)
                    magnitude = magnitude * 10 + (*probe - '0');
                ++probe;
            }
            binary_exponent += exponent_negative ? -magnitude : magnitude;
            result.ptr = probe;
        }
    }

    if (!nonzero) {
        value = float128();
        value.set_sign(negative ? 1 : 0);
        return result;
    }

    // Far outside the format's range only an overflow or a zero can come of it, whatever the exact
    // exponent; clamping keeps the arithmetic below in range without changing which.
    const int64_t clamped = (binary_exponent < -100000) ? -100000 : ((binary_exponent > 100000) ? 100000 : binary_exponent);
    value = float128::norm_round_pack(negative ? 1u : 0u, static_cast<int32_t>(clamped) + 112, l, h, sticky);
    if (value.is_inf() || value.is_zero())
        result.ec = std::errc::result_out_of_range;
    return result;
}
}  // namespace detail

/**
 * @brief Converts text to a float128, in the shape of std::from_chars.
 *
 * Reads the longest prefix of [first, last) that forms a number: an optional sign, decimal digits
 * with an optional point and an optional exponent, or one of inf, infinity, nan, nan(payload) and
 * snan. The result is the representable value nearest the one the text names, rounded once in the
 * current rounding direction - including when it is subnormal, where fewer than 113 bits are kept.
 *
 * @param first Start of the text
 * @param last One past the end of the text
 * @param value Receives the parsed value, untouched when nothing was parsed
 * @return ptr points past what was consumed; ec is invalid_argument when no number was found and
 *         result_out_of_range when the value is beyond the format's range, in which case value is
 *         set to the infinity or zero it rounded to.
 */
inline std::from_chars_result from_chars(const char* first, const char* last, float128& value)
{
    std::from_chars_result result {first, std::errc {}};
    const char* cursor = first;
    if (cursor == last) {
        result.ec = std::errc::invalid_argument;
        return result;
    }

    bool negative = false;
    if (*cursor == '-' || *cursor == '+') {
        negative = (*cursor == '-');
        ++cursor;
    }
    if (const char* end = detail::parse_special_value(cursor, last, negative, value)) {
        result.ptr = end;
        return result;
    }

    // The digits are collected without the point, which only decides the exponent.
    char digits[detail::MAX_SIGNIFICANT_DIGITS];
    int32_t count = 0;
    int32_t exponent10 = 0;
    bool any_digit = false;
    bool seen_significant = false;

    // What is being built is `digits` read as an integer, times 10^exponent10. A digit before the
    // point that does not fit raises the exponent, since its place is still there. One after the
    // point lowers it when it is kept, and a leading zero after the point lowers it as well even
    // though nothing is stored for it.
    const auto take_integer = [&](char digit) {
        any_digit = true;
        seen_significant = seen_significant || (digit != '0');
        if (!seen_significant)
            return;
        if (count < detail::MAX_SIGNIFICANT_DIGITS)
            digits[count++] = digit;
        else
            ++exponent10;
    };
    const auto take_fraction = [&](char digit) {
        any_digit = true;
        seen_significant = seen_significant || (digit != '0');
        if (!seen_significant) {
            --exponent10;
            return;
        }
        if (count < detail::MAX_SIGNIFICANT_DIGITS) {
            digits[count++] = digit;
            --exponent10;
        }
    };

    while (cursor < last && *cursor >= '0' && *cursor <= '9')
        take_integer(*cursor++);

    if (cursor < last && *cursor == '.') {
        ++cursor;
        while (cursor < last && *cursor >= '0' && *cursor <= '9')
            take_fraction(*cursor++);
    }

    if (!any_digit) {
        result.ec = std::errc::invalid_argument;
        return result;
    }
    result.ptr = cursor;

    // An exponent is only consumed when it is well formed, so that "1e" reads as 1.
    if (cursor < last && (*cursor == 'e' || *cursor == 'E')) {
        const char* probe = cursor + 1;
        bool exponent_negative = false;
        if (probe < last && (*probe == '-' || *probe == '+')) {
            exponent_negative = (*probe == '-');
            ++probe;
        }
        if (probe < last && *probe >= '0' && *probe <= '9') {
            int64_t magnitude = 0;
            while (probe < last && *probe >= '0' && *probe <= '9') {
                if (magnitude < 1000000)
                    magnitude = magnitude * 10 + (*probe - '0');
                ++probe;
            }
            exponent10 += static_cast<int32_t>(exponent_negative ? -magnitude : magnitude);
            result.ptr = probe;
        }
    }

    const uint32_t sign = negative ? 1u : 0u;
    uint64_t low = 0, high = 0, extra = 0;
    int32_t exponent = 0;
    if (!detail::from_decimal_digits(digits, count, exponent10, low, high, exponent, &extra)) {
        if (exponent == 0) {
            // the digits were all zeros
            value = float128();
            value.set_sign(sign);
            return result;
        }
        // Far beyond the range in one direction or the other. A value of the right size and a
        // sticky bit let the rounding decide the result, which depends on the direction: an
        // overflow can stop at the largest finite value, an underflow at the smallest subnormal.
        value = float128::round_pack(sign, (exponent > 0) ? 20000 : -20000, 0, 1ull << 48, 1);
        result.ec = std::errc::result_out_of_range;
        return result;
    }

    // The mantissa comes back unrounded, so a subnormal result is rounded once, to its own width.
    value = float128::round_pack(sign, exponent, low, high, extra);
    if (value.is_inf() || value.is_zero())
        result.ec = std::errc::result_out_of_range;
    return result;
}

/**
 * @brief Converts text to a float128 in a given format, in the shape of std::from_chars.
 *
 * chars_format::hex reads a hexadecimal significand and a binary exponent without the 0x prefix,
 * as std::from_chars does; every other format reads the decimal forms from_chars() above accepts.
 *
 * @param first Start of the text
 * @param last One past the end of the text
 * @param value Receives the parsed value, untouched when nothing was parsed
 * @param fmt The format to read
 * @return As from_chars() above.
 */
inline std::from_chars_result from_chars(const char* first, const char* last, float128& value, std::chars_format fmt)
{
    if (fmt == std::chars_format::hex)
        return detail::from_hex_chars(first, last, value);
    return from_chars(first, last, value);
}

/**
 * @brief Converts a float128 to text, in the shape of std::to_chars.
 * @param first Start of the output range
 * @param last One past the end of the output range
 * @param value Value to convert
 * @param fmt Layout to use
 * @param precision Digits after the point, or significant digits for `general`
 * @return ptr points past the last character written; ec is value_too_large when the range was too
 *         small, in which case nothing was written.
 */
inline std::to_chars_result to_chars(char* first, char* last, const float128& value, std::chars_format fmt, int precision)
{
    char style = '\0';
    switch (fmt) {
    case std::chars_format::scientific: style = 'e'; break;
    case std::chars_format::fixed:      style = 'f'; break;
    case std::chars_format::hex:        style = 'a'; break;
    default:                            style = 'g'; break;
    }

    const std::string text = detail::render(value, style, precision, false, false, '-');
    if (text.size() > static_cast<size_t>(last - first))
        return {last, std::errc::value_too_large};

    for (size_t i = 0; i < text.size(); ++i)
        first[i] = text[i];
    return {first + text.size(), std::errc {}};
}

/// @brief Converts a float128 to text with the default precision for the layout.
inline std::to_chars_result to_chars(char* first, char* last, const float128& value, std::chars_format fmt)
{
    return to_chars(first, last, value, fmt, -1);
}

/// @brief Converts a float128 to the shortest text that reads back as the same value.
inline std::to_chars_result to_chars(char* first, char* last, const float128& value)
{
    const std::string text = detail::render(value, '\0', -1, false, false, '-');
    if (text.size() > static_cast<size_t>(last - first))
        return {last, std::errc::value_too_large};

    for (size_t i = 0; i < text.size(); ++i)
        first[i] = text[i];
    return {first + text.size(), std::errc {}};
}

/**
 * @brief Writes a float128 to a stream, honouring its formatting state.
 *
 * The stream's precision, its fixed, scientific and hexfloat flags, showpos, showpoint, uppercase,
 * width, fill and adjustfield are all applied, so a float128 behaves in a stream the way a double
 * does.
 */
inline std::ostream& operator<<(std::ostream& os, const float128& value)
{
    const std::ios_base::fmtflags flags = os.flags();
    const std::ios_base::fmtflags style_flags = flags & std::ios_base::floatfield;

    char style = '\0';
    int32_t precision = static_cast<int32_t>(os.precision());
    if (style_flags == std::ios_base::fixed) {
        style = 'f';
    } else if (style_flags == std::ios_base::scientific) {
        style = 'e';
    } else if (style_flags == (std::ios_base::fixed | std::ios_base::scientific)) {
        style = 'a';
        precision = -1;  // hexfloat ignores the precision
    } else {
        // The default floatfield is the general layout, where the precision counts significant
        // digits and zero means one.
        style = 'g';
    }

    const char sign_char = (flags & std::ios_base::showpos) ? '+' : '-';
    const bool uppercase = (flags & std::ios_base::uppercase) != 0;
    const bool alternate = (flags & std::ios_base::showpoint) != 0;
    std::string text = detail::render(value, style, precision, uppercase, alternate, sign_char);

    const std::streamsize width = os.width();
    os.width(0);
    if (static_cast<std::streamsize>(text.size()) < width) {
        const size_t padding = static_cast<size_t>(width) - text.size();
        const std::ios_base::fmtflags adjust = flags & std::ios_base::adjustfield;
        if (adjust == std::ios_base::left) {
            text.append(padding, os.fill());
        } else if (adjust == std::ios_base::internal && !text.empty() && (text[0] == '-' || text[0] == '+')) {
            // internal puts the fill between the sign and the digits
            text.insert(1, padding, os.fill());
        } else {
            text.insert(0, padding, os.fill());
        }
    }
    return os << text;
}

/**
 * @brief Reads a float128 from a stream.
 *
 * Accepts what the constructor from a string does. Failure leaves the value untouched and sets
 * failbit, the way the builtin extractors do.
 */
inline std::istream& operator>>(std::istream& is, float128& value)
{
    std::string token;
    if (!(is >> token))
        return is;

    float128 parsed;
    const std::from_chars_result result = from_chars(token.data(), token.data() + token.size(), parsed);
    if (result.ec == std::errc::invalid_argument || result.ptr != token.data() + token.size()) {
        is.setstate(std::ios_base::failbit);
        return is;
    }

    value = parsed;
    return is;
}

}  // namespace fp128

namespace std
{
/**
 * @brief Numeric properties of fp128::float128, the binary128 interchange format.
 *
 * Specialized so that a function template written against a builtin floating point type - anything
 * reaching for numeric_limits<T>::epsilon() to size a tolerance, or for max() to seed a minimum -
 * compiles and behaves correctly when instantiated with float128.
 *
 * The values are given as encodings rather than computed, which keeps every one of them usable in
 * a constant expression and independent of the string parser.
 *
 * The const, volatile and const volatile forms need no specialization of their own: \<limits\>
 * already defines numeric_limits<cv T> to have the members of numeric_limits<T>.
 */
template <> class numeric_limits<fp128::float128>
{
public:
    static constexpr bool is_specialized = true;
    static constexpr bool is_signed = true;
    static constexpr bool is_integer = false;
    static constexpr bool is_exact = false;
    static constexpr bool has_infinity = true;
    static constexpr bool has_quiet_NaN = true;
    static constexpr bool has_signaling_NaN = true;
    static constexpr bool is_bounded = true;
    static constexpr bool is_modulo = false;
    /// @brief True only with FP128_IEEE_ENV, which adds the rounding directions and the exception
    ///        flags IEEE 754 requires; the format and the default arithmetic conform either way.
    ///        See float128::is754version2008().
    static constexpr bool is_iec559 = fp128::float128::is754version2008();
    /// @brief Exceptions raise flags (with FP128_IEEE_ENV) but never trap.
    static constexpr bool traps = false;
    /// @brief Tininess is detected after rounding, as x86 does for double.
    static constexpr bool tinyness_before = false;
    static constexpr float_round_style round_style = round_to_nearest;

    /// @brief Mantissa bits including the implicit leading one.
    static constexpr int digits = 113;
    /// @brief Decimal digits that survive a round trip through the type.
    static constexpr int digits10 = 33;
    /// @brief Decimal digits needed to distinguish every value of the type.
    static constexpr int max_digits10 = 36;
    static constexpr int radix = 2;
    static constexpr int min_exponent = -16381;
    static constexpr int max_exponent = 16384;
    static constexpr int min_exponent10 = -4931;
    static constexpr int max_exponent10 = 4932;

    // has_denorm and has_denorm_loss are deprecated since C++23 but remain part of the interface a
    // generic caller may read, so they are provided.
    FP128_SUPPRESS_DEPRECATED_BEGIN
    static constexpr float_denorm_style has_denorm = denorm_present;
    static constexpr bool has_denorm_loss = false;
    FP128_SUPPRESS_DEPRECATED_END

    /// @brief Smallest positive normal value, 2^-16382.
    [[nodiscard]] static constexpr fp128::float128 min() noexcept { return fp128::float128(0, 0x0001000000000000ull); }
    /// @brief Largest finite value, (2 - 2^-112) * 2^16383.
    [[nodiscard]] static constexpr fp128::float128 max() noexcept { return fp128::float128(UINT64_MAX, 0x7FFEFFFFFFFFFFFFull); }
    /// @brief Most negative finite value.
    [[nodiscard]] static constexpr fp128::float128 lowest() noexcept { return fp128::float128(UINT64_MAX, 0xFFFEFFFFFFFFFFFFull); }
    /// @brief Difference between one and the next larger value, 2^-112.
    [[nodiscard]] static constexpr fp128::float128 epsilon() noexcept { return fp128::float128(0, 0x3F8F000000000000ull); }
    /// @brief Largest rounding error in units in the last place, one half.
    [[nodiscard]] static constexpr fp128::float128 round_error() noexcept { return fp128::float128(0, 0x3FFE000000000000ull); }
    /// @brief Smallest positive subnormal, 2^-16494.
    [[nodiscard]] static constexpr fp128::float128 denorm_min() noexcept { return fp128::float128(1, 0); }
    [[nodiscard]] static constexpr fp128::float128 infinity() noexcept { return fp128::float128::inf(); }
    [[nodiscard]] static constexpr fp128::float128 quiet_NaN() noexcept { return fp128::float128::nan(); }
    [[nodiscard]] static constexpr fp128::float128 signaling_NaN() noexcept { return fp128::float128::signaling_nan(); }
};

/**
 * @brief std::format support for fp128::float128.
 *
 * Accepts the whole floating point format specification: fill and alignment, a sign, the alternate
 * form, zero padding, a width, a precision, and any of the a, A, e, E, f, F, g and G types. An
 * empty specification gives the shortest text that reads back as the same value, which is what
 * std::format produces for a double.
 *
 * A width or precision given as a nested replacement field is not supported; both have to be
 * written out as digits.
 */
template <> struct formatter<fp128::float128, char>
{
    constexpr auto parse(basic_format_parse_context<char>& context)
    {
        auto it = context.begin();
        const auto end = context.end();
        if (it == end || *it == '}')
            return it;

        // [[fill]align]
        const auto is_align = [](char c) { return c == '<' || c == '>' || c == '^'; };
        if (it + 1 != end && is_align(*(it + 1))) {
            fill = *it;
            align = *(it + 1);
            it += 2;
        } else if (is_align(*it)) {
            align = *it++;
        }

        // [sign]
        if (it != end && (*it == '+' || *it == '-' || *it == ' '))
            sign = *it++;

        // [#]
        if (it != end && *it == '#') {
            alternate = true;
            ++it;
        }

        // [0], which is an alignment rather than a fill in its own right
        if (it != end && *it == '0') {
            zero_pad = true;
            ++it;
        }

        // [width]
        while (it != end && *it >= '0' && *it <= '9')
            width = width * 10 + (*it++ - '0');

        // [.precision]
        if (it != end && *it == '.') {
            ++it;
            precision = 0;
            while (it != end && *it >= '0' && *it <= '9')
                precision = precision * 10 + (*it++ - '0');
        }

        // [type]
        if (it != end && *it != '}') {
            switch (*it) {
            case 'a':
            case 'A':
            case 'e':
            case 'E':
            case 'f':
            case 'F':
            case 'g':
            case 'G':
                type = *it++;
                break;
            default:
                throw format_error("invalid type in the format specification for float128");
            }
        }

        if (it != end && *it != '}')
            throw format_error("unmatched brace in the format specification for float128");
        return it;
    }

    template <typename FormatContext> auto format(const fp128::float128& value, FormatContext& context) const
    {
        const string text = render(value);
        return std::copy(text.begin(), text.end(), context.out());
    }

    /// @brief Builds the text, including the padding the specification asks for.
    [[nodiscard]] string render(const fp128::float128& value) const
    {
        const bool uppercase = (type >= 'A' && type <= 'Z');
        const char style = static_cast<char>(uppercase ? (type - 'A' + 'a') : type);
        string text = fp128::detail::render(value, style, precision, uppercase, alternate, sign);

        if (static_cast<int>(text.size()) >= width)
            return text;

        const size_t padding = static_cast<size_t>(width) - text.size();
        // Zero padding goes between the sign and the digits, and only when no alignment was given.
        // A special value is never zero padded, the way it is not for a double either.
        if (zero_pad && align == 0 && !value.is_nan() && !value.is_inf()) {
            const size_t offset = (!text.empty() && (text[0] == '-' || text[0] == '+' || text[0] == ' ')) ? 1u : 0u;
            text.insert(offset, padding, '0');
            return text;
        }

        switch (align) {
        case '<':
            text.append(padding, fill);
            break;
        case '^':
            text.insert(0, padding / 2, fill);
            text.append(padding - padding / 2, fill);
            break;
        case '>':
        default:
            // A number aligns right by default, as every arithmetic type does.
            text.insert(0, padding, fill);
            break;
        }
        return text;
    }

    char fill = ' ';       ///< Character the padding is made of.
    char align = 0;        ///< One of < > ^, or zero when none was given.
    char sign = '-';       ///< One of + - space.
    char type = 0;         ///< One of a A e E f F g G, or zero for the shortest form.
    bool alternate = false;///< The # flag: keep the point and the trailing zeros.
    bool zero_pad = false; ///< The 0 flag: pad with zeros after the sign.
    int width = 0;         ///< Minimum field width.
    int precision = -1;    ///< Digits after the point, or -1 when none was given.
};

/**
 * @brief Hash support, so a float128 can be a key in an unordered container.
 *
 * The two zeros compare equal and therefore have to hash equal, which the raw encoding would not
 * do: they differ in the sign bit.
 */
template <> struct hash<fp128::float128>
{
    [[nodiscard]] size_t operator()(const fp128::float128& value) const noexcept
    {
        uint64_t low = 0, high = 0;
        value.get_bits(low, high);
        if (value.is_zero())
            low = high = 0;

        // splitmix64's finalizer, applied to the two halves in turn
        uint64_t state = low + 0x9E3779B97F4A7C15ull;
        state = (state ^ (state >> 30)) * 0xBF58476D1CE4E5B9ull;
        state = (state ^ (state >> 27)) * 0x94D049BB133111EBull;
        state ^= high + 0x9E3779B97F4A7C15ull + (state << 6) + (state >> 2);
        state = (state ^ (state >> 30)) * 0xBF58476D1CE4E5B9ull;
        return static_cast<size_t>(state ^ (state >> 31));
    }
};

/// @name Common type
///
/// float128 is wider than every builtin arithmetic type, so a mixed expression produces a
/// float128. Without these the default common_type would try to form the ternary operator over the
/// two, which is ambiguous: float128 converts to double and double converts to float128.
/// @{
template <typename T>
    requires is_arithmetic_v<T>
struct common_type<fp128::float128, T>
{
    using type = fp128::float128;
};
template <typename T>
    requires is_arithmetic_v<T>
struct common_type<T, fp128::float128>
{
    using type = fp128::float128;
};
template <> struct common_type<fp128::float128, fp128::float128>
{
    using type = fp128::float128;
};
/// @}
}  // namespace std

#endif  // FP128_FLOAT128_H
