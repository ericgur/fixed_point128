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

// log2_ulp_dump.cpp : dumps log2 results as exact bit patterns, for fixed_point128 and float128.
//
// Half of the accuracy harness; log2_ulp_check.py is the other half and computes the reference.
// See tools/README.md.
//
// Nothing is printed in decimal: these types carry up to 127 fraction bits and a decimal rendering
// would lose exactly the bits being measured. Every value is emitted as
//
//     sign mantissa_high mantissa_low exponent
//
// standing for (-1)^sign * (mantissa_high * 2^64 + mantissa_low) * 2^exponent, which is exact for
// both types. One unit in the last place of a result is then 2^exponent of that result, so the
// same reader handles the fixed grid of fixed_point128 and the floating grid of float128.

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include "fixed_point128.h"
#include "float128.h"

using namespace fp128;

/** @brief Input class the samples currently being emitted belong to; tags each output line. */
static const char* g_group = "random";

/**
 * @brief Fraction bits of binary128, which is the layout float128 implements.
 *
 * Spelled out here rather than taken from float128, whose copy of it is private. These are fixed
 * by IEEE 754 and cannot drift.
 */
constexpr int32_t FLOAT128_FRAC_BITS = 112;

/** @brief Exponent bias of binary128. */
constexpr int32_t FLOAT128_EXP_BIAS = 0x3FFF;

/**
 * @brief Prints one value in the common exact form.
 * @param sign 1 for a negative value.
 * @param mant_high High QWORD of the mantissa.
 * @param mant_low Low QWORD of the mantissa.
 * @param exponent Power of two the mantissa is scaled by.
 */
static void EmitValue(uint32_t sign, uint64_t mant_high, uint64_t mant_low, int32_t exponent)
{
    printf(" %u %016llX %016llX %d", sign, static_cast<unsigned long long>(mant_high), static_cast<unsigned long long>(mant_low), exponent);
}

/**
 * @brief Reassembles the raw QWORDs of a fixed_point128 through the public get_bit() accessor.
 *
 * The members are private and this tool is not a friend, so the bits are read one at a time. It
 * runs once per sample and the cost is irrelevant next to the mpmath side of the comparison.
 *
 * @tparam I Number of integer bits.
 * @param v Value to take apart.
 * @param low Receives the low QWORD of the magnitude.
 * @param high Receives the high QWORD of the magnitude.
 */
template <int32_t I> static void RawBits(const fixed_point128<I>& v, uint64_t& low, uint64_t& high)
{
    low = 0;
    high = 0;
    for (uint32_t b = 0; b < 64; ++b) {
        low |= static_cast<uint64_t>(static_cast<uint32_t>(v.get_bit(b))) << b;
        high |= static_cast<uint64_t>(static_cast<uint32_t>(v.get_bit(b + 64))) << b;
    }
}

/**
 * @brief Emits one fixed_point128 sample: the argument and log2 of it.
 *
 * A fixed_point128<I> holds its value on a fixed grid of 2^-F, so both the mantissa scale and the
 * ulp are 2^-F regardless of how large the value is.
 *
 * @tparam I Number of integer bits.
 * @param x Argument to log2.
 */
template <int32_t I> static void Emit(const fixed_point128<I>& x)
{
    constexpr int32_t F = 128 - I;
    uint64_t x_low = 0, x_high = 0, r_low = 0, r_high = 0;
    RawBits(x, x_low, x_high);
    const fixed_point128<I> result = log2(x);
    RawBits(result, r_low, r_high);

    printf("fixed_point128<%d>/%s", I, g_group);
    EmitValue(static_cast<uint32_t>(x.is_negative()), x_high, x_low, -F);
    EmitValue(static_cast<uint32_t>(result.is_negative()), r_high, r_low, -F);
    printf("\n");
}

/**
 * @brief Emits one float128 sample: the argument and log2 of it.
 *
 * get_components() hands back the 113 bit mantissa with its implicit unity bit in place and the
 * unbiased exponent, having already normalized a subnormal argument, so the pair it returns is the
 * exact value with no decoding left to do here.
 *
 * @param x Argument to log2.
 */
static void EmitFloat(const float128& x)
{
    uint64_t x_low = 0, x_high = 0, r_low = 0, r_high = 0;
    int32_t x_expo = 0, r_expo = 0;
    uint32_t x_sign = 0, r_sign = 0;
    x.get_components(x_low, x_high, x_expo, x_sign);
    const float128 result = log2(x);
    if (result.is_zero() || !isfinite(result)) {
        return;
    }
    result.get_components(r_low, r_high, r_expo, r_sign);

    printf("float128/%s", g_group);
    EmitValue(x_sign, x_high, x_low, x_expo - FLOAT128_FRAC_BITS);
    EmitValue(r_sign, r_high, r_low, r_expo - FLOAT128_FRAC_BITS);
    printf("\n");
}

/**
 * @brief xorshift64*, so a run is reproducible and does not depend on the CRT's generator.
 * @param state Generator state, updated in place.
 * @return The next value.
 */
[[nodiscard]] static uint64_t Next(uint64_t& state)
{
    state ^= state >> 12;
    state ^= state << 25;
    state ^= state >> 27;

    return state * 0x2545F4914F6CDD1Dull;
}

/**
 * @brief Emits samples for one fixed_point128 instantiation, over the input classes that matter.
 *
 * Four groups, because they stress different parts of the implementation:
 * - Values in [1,2), where the series does all its work.
 * - Values over the whole representable range, which exercise the exponent path. Skipped when the
 *   instantiation has too few integer bits to hold the answer: log2 of the smallest representable
 *   value is -F, so 2^(I-1) has to exceed F or the result overflows and the comparison would be
 *   measuring that instead.
 * - Values just above one, the worst case for cancellation, where log2 is near zero.
 * - The reduction boundaries themselves and their immediate neighbours, where the table index
 *   changes and |z| is at its largest.
 *
 * @tparam I Number of integer bits.
 * @param count Samples per group.
 */
template <int32_t I> static void Sweep(uint64_t count)
{
    using fp = fixed_point128<I>;
    constexpr int32_t F = 128 - I;
    uint64_t seed = 0x123456789ABCDEF0ull;

    for (uint64_t n = 0; n < count; ++n) {
        const uint64_t high = Next(seed), low = Next(seed);
        // one() has bit F set; masking the fraction below it keeps the value inside [1,2)
        const fp v = fp::one() + fp(low, high & ((F >= 64) ? ((1ull << (F - 64)) - 1) : 0ull), 0);
        Emit<I>(v);
    }

    if constexpr (I >= 8 && (I >= 40 || (1ll << (I - 1)) > F)) {
        for (uint64_t n = 0; n < count; ++n) {
            const uint64_t high = Next(seed), low = Next(seed);
            const fp v = fp(low, high, 0);
            if (v.is_zero()) {
                continue;
            }
            Emit<I>(v);
        }
    }

    for (uint64_t n = 0; n < count; ++n) {
        const fp v = fp::one() + fp(Next(seed), 0, 0);
        Emit<I>(v);
    }

    constexpr int32_t entries = 1 << log2_reduction_bits;
    for (uint32_t j = 0; j < entries; ++j) {
        const fp step = fp::one() >> log2_reduction_bits;
        const fp boundary = fp::one() + step * j;
        Emit<I>(boundary);
        for (uint64_t d = 1; d <= 3; ++d) {
            Emit<I>(boundary + fp(d, 0, 0));
            if (j != 0) {
                Emit<I>(boundary - fp(d, 0, 0));
            }
        }
    }
}

/**
 * @brief Emits samples for float128, over the same input classes plus the exponent range.
 * @param count Samples per group.
 */
static void SweepFloat(uint64_t count)
{
    uint64_t seed = 0x0FEDCBA987654321ull;
    constexpr int32_t frac_bits = FLOAT128_FRAC_BITS;

    // Mantissas spread over [1,2) at exponent zero, where the series does all its work.
    g_group = "unit-range";
    for (uint64_t n = 0; n < count; ++n) {
        const uint64_t high = Next(seed) & FP128_MAX_VALUE_64(frac_bits - 64);
        EmitFloat(float128(Next(seed), high, FLOAT128_EXP_BIAS, 0));
    }

    // The same, but spread over the exponent range as well.
    g_group = "wide-exponent";
    for (uint64_t n = 0; n < count; ++n) {
        const uint64_t high = Next(seed) & FP128_MAX_VALUE_64(frac_bits - 64);
        const uint32_t expo = static_cast<uint32_t>(Next(seed) % 32000) + 200;
        EmitFloat(float128(Next(seed), high, expo, 0));
    }

    // Mantissas just above one, the worst case for cancellation.
    g_group = "near-one";
    for (uint64_t n = 0; n < count; ++n) {
        EmitFloat(float128(Next(seed), 0, FLOAT128_EXP_BIAS, 0));
    }

    // The reduction boundaries and their immediate neighbours.
    g_group = "boundaries";
    constexpr int32_t entries = 1 << log2_reduction_bits;
    for (uint64_t j = 0; j < entries; ++j) {
        // mantissa = 1 + j/64, i.e. the top log2_reduction_bits fraction bits set to j
        const uint64_t high = j << (frac_bits - 64 - log2_reduction_bits);
        for (uint64_t d = 0; d <= 3; ++d) {
            EmitFloat(float128(d, high, FLOAT128_EXP_BIAS, 0));
            if (d != 0) {
                EmitFloat(float128(0ull - d, high - 1, FLOAT128_EXP_BIAS, 0));
            }
        }
    }
}

int main(int argc, char* argv[])
{
    const uint64_t count = (argc > 1) ? strtoull(argv[1], nullptr, 10) : 1000;

    // A spread of integer bit counts: the two extremes, the value the benchmark uses, and enough
    // in between to catch a mistake that only shows up at a particular fraction width.
    Sweep<1>(count);
    Sweep<2>(count);
    Sweep<3>(count);
    Sweep<4>(count);
    Sweep<6>(count);
    Sweep<8>(count);
    Sweep<10>(count);
    Sweep<20>(count);
    Sweep<32>(count);
    Sweep<64>(count);
    SweepFloat(count);

    return 0;
}
