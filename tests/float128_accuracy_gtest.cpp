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
#include "float128_ref_check.h"

/**********************************************************************
 * float128 accuracy tests
 *
 * Every case here compares a full 113 bit result against the correctly rounded one, which the
 * rest of the suite cannot do: it checks against a double and so leaves the bottom 60 mantissa
 * bits unverified. The bound on each test is the error the implementation actually produces over
 * its reference table, rounded up a little, so it is both a regression guard and the documented
 * accuracy of that function.
 *
 * Run with FP128_PRINT_ULP set in the environment to see the measured error of each one.
 *
 * Three of the bounds are far above the rest and say something about the implementation rather
 * than about rounding. erf sums up to eighty terms, each rounding against a running total, and
 * erfc runs a continued fraction of comparable length. pow is exp(y*log(x)): log is accurate to
 * its own last bit, but multiplying by a y that carries the product up to 11356 turns that
 * relative error into an absolute one, and the exponential magnifies it. Fixing either properly
 * needs the intermediate carried in more than 113 bits.
 ***********************************************************************/

using namespace f128_ref;

TEST(float128_accuracy, sqrt)
{
    CheckUnary("sqrt", sqrt_ref, [](const float128& x) { return sqrt(x); }, 1);
}
TEST(float128_accuracy, cbrt)
{
    CheckUnary("cbrt", cbrt_ref, [](const float128& x) { return cbrt(x); }, 2);
}
TEST(float128_accuracy, exp)
{
    CheckUnary("exp", exp_ref, [](const float128& x) { return exp(x); }, 2);
}
TEST(float128_accuracy, exp2)
{
    CheckUnary("exp2", exp2_ref, [](const float128& x) { return exp2(x); }, 2);
}
TEST(float128_accuracy, expm1)
{
    CheckUnary("expm1", expm1_ref, [](const float128& x) { return expm1(x); }, 2);
}
TEST(float128_accuracy, log)
{
    CheckUnary("log", log_ref, [](const float128& x) { return log(x); }, 2);
}
TEST(float128_accuracy, log2)
{
    CheckUnary("log2", log2_ref, [](const float128& x) { return log2(x); }, 4);
}
TEST(float128_accuracy, log10)
{
    CheckUnary("log10", log10_ref, [](const float128& x) { return log10(x); }, 2);
}
TEST(float128_accuracy, log1p)
{
    CheckUnary("log1p", log1p_ref, [](const float128& x) { return log1p(x); }, 2);
}
TEST(float128_accuracy, sin)
{
    CheckUnary("sin", sin_ref, [](const float128& x) { return sin(x); }, 6);
}
TEST(float128_accuracy, cos)
{
    CheckUnary("cos", cos_ref, [](const float128& x) { return cos(x); }, 8);
}
TEST(float128_accuracy, tan)
{
    CheckUnary("tan", tan_ref, [](const float128& x) { return tan(x); }, 8);
}
TEST(float128_accuracy, asin)
{
    CheckUnary("asin", asin_ref, [](const float128& x) { return asin(x); }, 12);
}
TEST(float128_accuracy, acos)
{
    CheckUnary("acos", acos_ref, [](const float128& x) { return acos(x); }, 32);
}
TEST(float128_accuracy, atan)
{
    CheckUnary("atan", atan_ref, [](const float128& x) { return atan(x); }, 6);
}
TEST(float128_accuracy, sinh)
{
    CheckUnary("sinh", sinh_ref, [](const float128& x) { return sinh(x); }, 2);
}
TEST(float128_accuracy, cosh)
{
    CheckUnary("cosh", cosh_ref, [](const float128& x) { return cosh(x); }, 2);
}
TEST(float128_accuracy, tanh)
{
    CheckUnary("tanh", tanh_ref, [](const float128& x) { return tanh(x); }, 2);
}
TEST(float128_accuracy, asinh)
{
    CheckUnary("asinh", asinh_ref, [](const float128& x) { return asinh(x); }, 2);
}
TEST(float128_accuracy, acosh)
{
    CheckUnary("acosh", acosh_ref, [](const float128& x) { return acosh(x); }, 2);
}
TEST(float128_accuracy, atanh)
{
    CheckUnary("atanh", atanh_ref, [](const float128& x) { return atanh(x); }, 4);
}
TEST(float128_accuracy, erf)
{
    CheckUnary("erf", erf_ref, [](const float128& x) { return erf(x); }, 128);
}
TEST(float128_accuracy, erfc)
{
    CheckUnary("erfc", erfc_ref, [](const float128& x) { return erfc(x); }, 32);
}
TEST(float128_accuracy, atan2)
{
    CheckBinary("atan2", atan2_ref, [](const float128& y, const float128& x) { return atan2(y, x); }, 4);
}
TEST(float128_accuracy, hypot)
{
    CheckBinary("hypot", hypot_ref, [](const float128& x, const float128& y) { return hypot(x, y); }, 2);
}
TEST(float128_accuracy, pow)
{
    CheckBinary("pow", pow_ref, [](const float128& x, const float128& y) { return pow(x, y); }, 256);
}
TEST(float128_accuracy, fmod)
{
    // fmod is exact: the result is a difference of representable values, so anything above zero
    // is a defect rather than a rounding effect.
    CheckBinary("fmod", fmod_ref, [](const float128& x, const float128& y) { return fmod(x, y); }, 0);
}
TEST(float128_accuracy, fma)
{
    // fma is required to round the whole x*y+z once, which makes it exact to the last bit.
    CheckTernary("fma", fma_ref, [](const float128& x, const float128& y, const float128& z) { return fma(x, y, z); }, 0);
}
TEST(float128_accuracy, remainder)
{
    // Like fmod, an exact operation: the result is representable and no rounding takes place.
    CheckBinary("remainder", remainder_ref, [](const float128& x, const float128& y) { return remainder(x, y); }, 0);
}
TEST(float128_accuracy, tgamma)
{
    // exp(lgamma) is what limits this: at the top of the range the logarithm reaches 7000, so its
    // own last bit lands in the seventh place of the result.
    CheckUnary("tgamma", tgamma_ref, [](const float128& x) { return tgamma(x); }, 1u << 15);
}
TEST(float128_accuracy, lgamma)
{
    CheckUnary("lgamma", lgamma_ref, [](const float128& x) { return lgamma(x); }, 1u << 12);
}
