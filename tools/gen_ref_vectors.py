#!/usr/bin/env python3
"""Generate binary128 reference vectors for the float128 test suite.

The test suite used to check float128 math results against the double returned by the CRT, which
only pins down 53 of the 113 mantissa bits. This script produces the missing 60: every reference
value is computed by mpmath at 240 bits and then rounded to binary128 with round-half-even, so the
committed vector is the correctly rounded result and a test can assert a ulp bound against it.

Inputs are drawn as binary128 bit patterns and decoded exactly, so the value the reference was
computed from is representable and the comparison has no second rounding hiding in it.

Usage:
    python tools/gen_ref_vectors.py [-o tests/float128_ref_data.h]

The output is committed. Re-run it only when the vector set itself changes; the seed is fixed, so
an unchanged configuration reproduces an identical file.
"""

import argparse
import random
import sys
from fractions import Fraction

try:
    import mpmath
    from mpmath import mp, mpf
except ImportError:  # pragma: no cover - a helpful message beats a traceback
    sys.exit("This script needs mpmath: python -m pip install mpmath")

# 240 bits leaves better than 120 bits of headroom over the 113 a binary128 keeps, so the
# double rounding through mpmath cannot move the final result off the correctly rounded one for
# any input that is not astronomically close to a rounding boundary.
mp.prec = 240

SEED = 0x12345678

# binary128 format parameters
P = 113                     # significant bits, including the implicit one
FRAC_BITS = 112
EXP_BIAS = 16383
EXP_MAX = 0x7FFF            # biased exponent of inf/NaN
MIN_SUB_EXP = -16494        # exponent of the least significant bit of the smallest subnormal
MAX_EXP = 16383             # unbiased exponent of the largest finite value
MIN_NORM_EXP = -16382       # unbiased exponent of the smallest normal

ZERO = (0, 0)
NEG_ZERO = (0, 1 << 63)
INF = (0, EXP_MAX << 48)
NEG_INF = (0, (1 << 63) | (EXP_MAX << 48))
NAN = (1, EXP_MAX << 48)


# ---------------------------------------------------------------------------
# binary128 encoding and decoding
# ---------------------------------------------------------------------------

def encode(x):
    """Round an mpmath value to binary128 and return it as a (low, high) QWORD pair."""
    if mpmath.isnan(x):
        return NAN
    if mpmath.isinf(x):
        return NEG_INF if x < 0 else INF

    sign, man, exp, bc = mpf(x)._mpf_
    if man == 0:
        return NEG_ZERO if sign else ZERO

    # The value is man * 2^exp with man holding bc bits, so its leading one has weight 2^msb.
    msb = exp + bc - 1
    # Quantize to the last place of the result: one ulp below the leading bit for a normal value,
    # or the fixed smallest subnormal step when that would fall below the format's floor.
    unit = max(msb - FRAC_BITS, MIN_SUB_EXP)

    shift = exp - unit
    if shift >= 0:
        q = man << shift
    else:
        dropped = -shift
        q = man >> dropped
        rest = man & ((1 << dropped) - 1)
        half = 1 << (dropped - 1)
        # round half to even
        if rest > half or (rest == half and (q & 1)):
            q += 1

    if q == 0:
        return NEG_ZERO if sign else ZERO

    # Rounding up can carry into an extra bit, which is exactly a step to the next power of two.
    if q.bit_length() > P:
        assert q.bit_length() == P + 1 and (q & 1) == 0
        q >>= 1
        unit += 1

    if q.bit_length() < P:
        # Subnormal: the leading one sits below 2^-16382, so there is no implicit bit and the
        # biased exponent field stays zero.
        assert unit == MIN_SUB_EXP
        high = q >> 64
    else:
        biased = unit + FRAC_BITS + EXP_BIAS
        if biased >= EXP_MAX:
            return NEG_INF if sign else INF
        high = ((q >> 64) & ((1 << 48) - 1)) | (biased << 48)

    if sign:
        high |= 1 << 63
    return (q & ((1 << 64) - 1), high)


def decode_fraction(bits):
    """Decode a (low, high) QWORD pair into an exact Fraction.

    The rational is what the two remainder functions need: their results are defined by an exact
    integer quotient, which for arguments a few hundred binary exponents apart runs to more digits
    than any working precision would hold.
    """
    low, high = bits
    sign = high >> 63
    biased = (high >> 48) & EXP_MAX
    man = ((high & ((1 << 48) - 1)) << 64) | low
    if biased == 0:
        value = Fraction(man, 1 << -MIN_SUB_EXP)
    else:
        man |= 1 << FRAC_BITS
        value = Fraction(man, 1) * Fraction(2) ** (biased - EXP_BIAS - FRAC_BITS)
    return -value if sign else value


def decode(bits):
    """Decode a (low, high) QWORD pair into an exact mpmath value.

    Only finite values are decoded; the generators never feed a special value to a reference
    computation, they are covered by dedicated tests on the C++ side instead.
    """
    value = decode_fraction(bits)
    return mpf(value.numerator) / mpf(value.denominator)


def to_fraction(value):
    """Convert an int, float, Fraction or mpf to an exact Fraction.

    mpf has no as_integer_ratio, so Fraction() cannot take one directly; its internal
    (sign, mantissa, exponent) triple is exact and converts without loss.
    """
    if isinstance(value, Fraction):
        return value
    if isinstance(value, (int, float)):
        return Fraction(value)
    sign, man, exp, _ = mpf(value)._mpf_
    result = Fraction(man) * Fraction(2) ** exp
    return -result if sign else result


def from_real(value):
    """Encode a Python int, float, Fraction or mpf as binary128."""
    if isinstance(value, Fraction):
        return encode(mpf(value.numerator) / mpf(value.denominator))
    return encode(mpf(value))


# ---------------------------------------------------------------------------
# Input sampling
# ---------------------------------------------------------------------------

class Sampler:
    """Draws binary128 inputs. One instance per run keeps the whole file reproducible."""

    def __init__(self, seed):
        self.rng = random.Random(seed)

    def bits(self, exp_lo, exp_hi, signed=True):
        """A value with a random 112-bit fraction and an exponent drawn uniformly in the range.

        Sampling the exponent rather than the value covers the format the way it is actually laid
        out: half of all binary128 values are below one, and a uniform draw over a linear range
        would never produce them.
        """
        e = self.rng.randint(exp_lo, exp_hi)
        frac = self.rng.getrandbits(FRAC_BITS)
        sign = self.rng.randint(0, 1) if signed else 0
        biased = e + EXP_BIAS
        return (frac & ((1 << 64) - 1),
                (sign << 63) | (biased << 48) | (frac >> 64))

    def linear(self, lo, hi):
        """A value drawn uniformly from a linear range, as a binary128."""
        # 128 random bits of resolution, so the low mantissa bits are not all zero
        t = Fraction(self.rng.getrandbits(128), 1 << 128)
        lo, hi = to_fraction(lo), to_fraction(hi)
        return from_real(lo + t * (hi - lo))

    def near(self, center, exp_lo, exp_hi, signed=True):
        """A value at a random small offset from a center point.

        This is what exercises the argument reduction of the trigonometric functions and the
        cancellation-prone region of log1p and expm1. Pass signed=False to stay on one side of the
        center, which is what a function whose domain ends there needs.
        """
        delta = decode(self.bits(exp_lo, exp_hi, signed))
        return encode(mpf(center) + delta)


# ---------------------------------------------------------------------------
# Vector set definitions
# ---------------------------------------------------------------------------

UNARY_COUNT = 80
BINARY_COUNT = 64
TERNARY_COUNT = 64


def unary_sets(s):
    """Return {name: (mpmath_reference, [input_bits, ...])} for the single argument functions."""
    def draw(n, fn):
        return [fn() for _ in range(n)]

    # Exponent windows are chosen per function so the samples land where the function is defined
    # and where its implementation actually has to work, rather than saturating to 0 or inf.
    sets = {}

    sets['sqrt'] = (mpmath.sqrt,
                    draw(UNARY_COUNT, lambda: s.bits(-16300, 16300, signed=False)))
    # mpmath.cbrt takes the principal root, which is complex for a negative argument. The C
    # library's cbrt is the real root instead, and it is odd.
    sets['cbrt'] = (lambda x: -mpmath.cbrt(-x) if x < 0 else mpmath.cbrt(x),
                    draw(UNARY_COUNT, lambda: s.bits(-16300, 16300)))
    # exp overflows above ln(max) = 11356.5 and underflows to zero below ln(min_denorm) = -11433
    sets['exp'] = (mpmath.exp,
                   draw(UNARY_COUNT // 2, lambda: s.linear(-11350, 11350)) +
                   draw(UNARY_COUNT // 2, lambda: s.bits(-120, 6)))
    sets['exp2'] = (lambda x: mpmath.power(2, x),
                    draw(UNARY_COUNT // 2, lambda: s.linear(-16380, 16380)) +
                    draw(UNARY_COUNT // 2, lambda: s.bits(-120, 6)))
    # expm1 exists for the region where exp(x)-1 loses everything to cancellation
    sets['expm1'] = (mpmath.expm1,
                     draw(UNARY_COUNT // 2, lambda: s.bits(-140, 0)) +
                     draw(UNARY_COUNT // 2, lambda: s.linear(-40, 40)))
    sets['log'] = (mpmath.log,
                   draw(UNARY_COUNT // 2, lambda: s.bits(-16300, 16300, signed=False)) +
                   draw(UNARY_COUNT // 2, lambda: s.near(1, -60, -1)))
    sets['log2'] = (lambda x: mpmath.log(x, 2),
                    draw(UNARY_COUNT // 2, lambda: s.bits(-16300, 16300, signed=False)) +
                    draw(UNARY_COUNT // 2, lambda: s.near(1, -60, -1)))
    sets['log10'] = (mpmath.log10,
                     draw(UNARY_COUNT // 2, lambda: s.bits(-16300, 16300, signed=False)) +
                     draw(UNARY_COUNT // 2, lambda: s.near(1, -60, -1)))
    sets['log1p'] = (mpmath.log1p,
                     draw(UNARY_COUNT // 2, lambda: s.bits(-140, 0)) +
                     draw(UNARY_COUNT // 2, lambda: s.bits(1, 200, signed=False)))

    sets['sin'] = (mpmath.sin,
                   draw(UNARY_COUNT // 2, lambda: s.linear(-100, 100)) +
                   draw(UNARY_COUNT // 4, lambda: s.bits(-120, 0)) +
                   draw(UNARY_COUNT // 4, lambda: s.near(mpmath.pi, -60, -20)))
    sets['cos'] = (mpmath.cos,
                   draw(UNARY_COUNT // 2, lambda: s.linear(-100, 100)) +
                   draw(UNARY_COUNT // 4, lambda: s.bits(-120, 0)) +
                   draw(UNARY_COUNT // 4, lambda: s.near(mpmath.pi / 2, -60, -20)))
    sets['tan'] = (mpmath.tan,
                   draw(UNARY_COUNT // 2, lambda: s.linear(-100, 100)) +
                   draw(UNARY_COUNT // 2, lambda: s.bits(-120, 0)))
    sets['asin'] = (mpmath.asin,
                    draw(UNARY_COUNT // 2, lambda: s.linear(-1, 1)) +
                    draw(UNARY_COUNT // 2, lambda: s.bits(-120, -1)))
    sets['acos'] = (mpmath.acos,
                    draw(UNARY_COUNT // 2, lambda: s.linear(-1, 1)) +
                    draw(UNARY_COUNT // 2, lambda: s.bits(-120, -1)))
    sets['atan'] = (mpmath.atan,
                    draw(UNARY_COUNT // 2, lambda: s.bits(-16300, 16300)) +
                    draw(UNARY_COUNT // 2, lambda: s.linear(-10, 10)))

    sets['sinh'] = (mpmath.sinh,
                    draw(UNARY_COUNT // 2, lambda: s.linear(-11350, 11350)) +
                    draw(UNARY_COUNT // 2, lambda: s.bits(-140, 2)))
    sets['cosh'] = (mpmath.cosh,
                    draw(UNARY_COUNT // 2, lambda: s.linear(-11350, 11350)) +
                    draw(UNARY_COUNT // 2, lambda: s.bits(-140, 2)))
    sets['tanh'] = (mpmath.tanh,
                    draw(UNARY_COUNT // 2, lambda: s.linear(-30, 30)) +
                    draw(UNARY_COUNT // 2, lambda: s.bits(-140, 2)))
    sets['asinh'] = (mpmath.asinh,
                     draw(UNARY_COUNT // 2, lambda: s.bits(-16300, 8000)) +
                     draw(UNARY_COUNT // 2, lambda: s.bits(-140, 2)))
    sets['acosh'] = (mpmath.acosh,
                     draw(UNARY_COUNT // 2, lambda: s.bits(1, 8000, signed=False)) +
                     draw(UNARY_COUNT // 2, lambda: s.near(1, -60, -1, signed=False)))
    sets['atanh'] = (mpmath.atanh,
                     draw(UNARY_COUNT // 2, lambda: s.bits(-140, -1)) +
                     draw(UNARY_COUNT // 2, lambda: s.linear(mpf(-1) + mpf(2) ** -60,
                                                            mpf(1) - mpf(2) ** -60)))

    sets['erf'] = (mpmath.erf,
                   draw(UNARY_COUNT // 2, lambda: s.linear(-10, 10)) +
                   draw(UNARY_COUNT // 2, lambda: s.bits(-140, 1)))
    # erfc underflows to zero a little above 106
    sets['erfc'] = (mpmath.erfc,
                    draw(UNARY_COUNT // 2, lambda: s.linear(-6, 106)) +
                    draw(UNARY_COUNT // 2, lambda: s.bits(-140, 1)))
    # tgamma overflows above 1755.5
    sets['tgamma'] = (mpmath.gamma,
                      draw(UNARY_COUNT // 2, lambda: s.linear(mpf(2) ** -40, 1755)) +
                      draw(UNARY_COUNT // 4, lambda: s.linear(-20, -1 - mpf(2) ** -20)) +
                      draw(UNARY_COUNT // 4, lambda: s.bits(-140, -1, signed=False)))
    # lgamma is log|gamma|, defined everywhere except the non positive integers
    sets['lgamma'] = (lambda x: mpmath.log(abs(mpmath.gamma(x))),
                      draw(UNARY_COUNT // 2, lambda: s.bits(-60, 60, signed=False)) +
                      draw(UNARY_COUNT // 4, lambda: s.linear(mpf(2) ** -40, 1000)) +
                      draw(UNARY_COUNT // 4, lambda: s.linear(-20, -1 - mpf(2) ** -20)))
    return sets


def binary_sets(s):
    """Return {name: (mpmath_reference, [(x_bits, y_bits), ...])} for the two argument functions."""
    def draw(n, fn):
        return [fn() for _ in range(n)]

    sets = {}
    sets['atan2'] = (mpmath.atan2,
                     draw(BINARY_COUNT // 2, lambda: (s.bits(-200, 200), s.bits(-200, 200))) +
                     draw(BINARY_COUNT // 2, lambda: (s.linear(-10, 10), s.linear(-10, 10))))
    # A large exponent on either side would overflow or flush the result, so both stay moderate.
    sets['hypot'] = (mpmath.hypot,
                     draw(BINARY_COUNT, lambda: (s.bits(-8000, 8000), s.bits(-8000, 8000))))
    sets['pow'] = (lambda x, y: mpmath.power(x, y),
                   draw(BINARY_COUNT // 2, lambda: (s.linear(mpf(2) ** -20, 100),
                                                    s.linear(-30, 30))) +
                   draw(BINARY_COUNT // 2, lambda: (s.bits(-60, 60, signed=False),
                                                    s.linear(-8, 8))))
    # Both remainder functions are exact, and both are emitted from the bit patterns rather than
    # from decoded values: the quotient they are defined through is an integer of up to four
    # hundred bits for arguments this far apart, and mpmath's fmod follows Python's floored
    # convention rather than the truncated one C specifies.
    sets['fmod'] = (None, draw(BINARY_COUNT, lambda: (s.bits(-200, 200), s.bits(-200, 200))))
    sets['remainder'] = (None, draw(BINARY_COUNT, lambda: (s.bits(-200, 200), s.bits(-200, 200))))
    return sets


def exact_fmod(x_bits, y_bits):
    """C's fmod: x - y*n with n the quotient truncated toward zero."""
    x = decode_fraction(x_bits)
    y = decode_fraction(y_bits)
    return from_real(x - y * int(x / y))


def exact_remainder(x_bits, y_bits):
    """C's remainder: x - y*n with n the quotient rounded half to even."""
    x = decode_fraction(x_bits)
    y = decode_fraction(y_bits)
    quotient = x / y
    floor = quotient.numerator // quotient.denominator
    diff = quotient - floor
    if diff > Fraction(1, 2):
        n = floor + 1
    elif diff < Fraction(1, 2):
        n = floor
    else:
        n = floor if floor % 2 == 0 else floor + 1
    return from_real(x - y * n)


EXACT_BINARY = {'fmod': exact_fmod, 'remainder': exact_remainder}


def ternary_sets(s):
    """Return {name: (mpmath_reference, [(x, y, z), ...])} for the three argument functions."""
    def draw(n, fn):
        return [fn() for _ in range(n)]

    # fma has to be exact in the product before the add, so the interesting inputs are the ones
    # where x*y and z very nearly cancel and every bit of the discarded half of the product
    # decides the result.
    def cancelling():
        x = s.bits(-60, 60)
        y = s.bits(-60, 60)
        product = decode(x) * decode(y)
        # perturb the negated product by a few ulp so the sum lands near, but not on, zero
        z = encode(-product * (1 + decode(s.bits(-110, -100))))
        return (x, y, z)

    return {
        'fma': (lambda x, y, z: x * y + z,
                draw(TERNARY_COUNT // 2, cancelling) +
                draw(TERNARY_COUNT // 2, lambda: (s.bits(-200, 200), s.bits(-200, 200),
                                                  s.bits(-200, 200)))),
    }


# ---------------------------------------------------------------------------
# Emission
# ---------------------------------------------------------------------------

def qword(x):
    return f"0x{x:016X}"


def emit_unary(out, name, ref, inputs):
    out.append(f"inline constexpr ref_unary {name}_ref[] = {{")
    for bits in inputs:
        result = encode(ref(decode(bits)))
        out.append(f"    {{{qword(bits[0])}, {qword(bits[1])}, "
                   f"{qword(result[0])}, {qword(result[1])}}},")
    out.append("};")
    out.append("")


def emit_binary(out, name, ref, inputs):
    exact = EXACT_BINARY.get(name)
    out.append(f"inline constexpr ref_binary {name}_ref[] = {{")
    for x, y in inputs:
        result = exact(x, y) if exact else encode(ref(decode(x), decode(y)))
        out.append(f"    {{{qword(x[0])}, {qword(x[1])}, {qword(y[0])}, {qword(y[1])}, "
                   f"{qword(result[0])}, {qword(result[1])}}},")
    out.append("};")
    out.append("")


def emit_ternary(out, name, ref, inputs):
    out.append(f"inline constexpr ref_ternary {name}_ref[] = {{")
    for x, y, z in inputs:
        result = encode(ref(decode(x), decode(y), decode(z)))
        out.append(f"    {{{qword(x[0])}, {qword(x[1])}, {qword(y[0])}, {qword(y[1])}, "
                   f"{qword(z[0])}, {qword(z[1])}, {qword(result[0])}, {qword(result[1])}}},")
    out.append("};")
    out.append("")


HEADER = '''// Generated by tools/gen_ref_vectors.py - do not edit by hand.
//
// Correctly rounded binary128 reference values for the float128 math functions, computed by
// mpmath at 240 bits and rounded to the format with round half to even. Each input is itself a
// binary128 bit pattern, so the only rounding anywhere in a comparison against these tables is
// the one the function under test performs.

#ifndef FP128_FLOAT128_REF_DATA_H
#define FP128_FLOAT128_REF_DATA_H

#include <cstdint>

namespace f128_ref
{
/// @brief One reference case for a single argument function: f(x) == r.
struct ref_unary {
    uint64_t xl, xh;  ///< Argument, low and high QWORD of the binary128 encoding.
    uint64_t rl, rh;  ///< Correctly rounded result.
};

/// @brief One reference case for a two argument function: f(x, y) == r.
struct ref_binary {
    uint64_t xl, xh;  ///< First argument.
    uint64_t yl, yh;  ///< Second argument.
    uint64_t rl, rh;  ///< Correctly rounded result.
};

/// @brief One reference case for a three argument function: f(x, y, z) == r.
struct ref_ternary {
    uint64_t xl, xh;  ///< First argument.
    uint64_t yl, yh;  ///< Second argument.
    uint64_t zl, zh;  ///< Third argument.
    uint64_t rl, rh;  ///< Correctly rounded result.
};

'''

FOOTER = '''}  // namespace f128_ref

#endif  // FP128_FLOAT128_REF_DATA_H
'''


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-o', '--output', default='tests/float128_ref_data.h',
                        help='header file to write (default: tests/float128_ref_data.h)')
    args = parser.parse_args()

    sampler = Sampler(SEED)
    out = []

    for name, (ref, inputs) in unary_sets(sampler).items():
        print(f"  {name} ({len(inputs)} cases)", file=sys.stderr)
        emit_unary(out, name, ref, inputs)
    for name, (ref, inputs) in binary_sets(sampler).items():
        print(f"  {name} ({len(inputs)} cases)", file=sys.stderr)
        emit_binary(out, name, ref, inputs)
    for name, (ref, inputs) in ternary_sets(sampler).items():
        print(f"  {name} ({len(inputs)} cases)", file=sys.stderr)
        emit_ternary(out, name, ref, inputs)

    with open(args.output, 'w', encoding='utf-8', newline='\n') as f:
        f.write(HEADER)
        f.write('\n'.join(out))
        f.write(FOOTER)
    print(f"wrote {args.output}", file=sys.stderr)


if __name__ == '__main__':
    main()
