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
from math import isqrt

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


# ---------------------------------------------------------------------------
# Correctly rounded operations (IEEE 754-2008 clause 5)
# ---------------------------------------------------------------------------
#
# The operations IEEE 754 requires to be correctly rounded are checked bit for bit, sign of zero
# included, rather than within an ulp bound. Their references are computed on exact rationals:
# encode() rounds through a 240 bit mpmath value first, which is plenty for a transcendental but
# not for a sum like 1 + 2^-113 + 2^-500, where the 240 bit intermediate drops the bit that breaks
# the tie. A random draw also rarely lands where rounding is hard, so the inputs are aimed at
# those places: ties and values a sticky bit away from them, deep cancellation, results in the
# subnormal range and at the edge of overflow.
#
# They come from a sampler of their own, so the tables above keep their contents.

EXACT_SEED = 0x7542008
EXACT_COUNT = 160


def encode_exact(q, zero_sign=0):
    """Round an exact Fraction to binary128, ties to even, with no intermediate rounding.

    zero_sign is the sign to give an exactly zero result, which the operation decides.
    """
    if q == 0:
        return NEG_ZERO if zero_sign else ZERO
    sign = 1 if q < 0 else 0
    a = -q if sign else q
    msb = a.numerator.bit_length() - a.denominator.bit_length()
    if Fraction(2) ** msb > a:
        msb -= 1
    unit = max(msb - FRAC_BITS, MIN_SUB_EXP)
    scaled = a / Fraction(2) ** unit
    n = scaled.numerator // scaled.denominator
    rest = scaled - n
    if rest > Fraction(1, 2) or (rest == Fraction(1, 2) and n & 1):
        n += 1
    if n == 0:
        return NEG_ZERO if sign else ZERO
    if n.bit_length() > P:
        n >>= 1
        unit += 1
    if n.bit_length() < P:
        high = n >> 64
    else:
        biased = unit + FRAC_BITS + EXP_BIAS
        if biased >= EXP_MAX:
            return NEG_INF if sign else INF
        high = ((n >> 64) & ((1 << 48) - 1)) | (biased << 48)
    if sign:
        high |= 1 << 63
    return (n & ((1 << 64) - 1), high)


def round_to_binary(q, p, emin, emax):
    """Round an exact Fraction to a binary format of p bits and return it as a Fraction.

    An overflow comes back as a Fraction too large for binary128, which encode_exact() turns into
    an infinity of the right sign.
    """
    if q == 0:
        return q
    sign = -1 if q < 0 else 1
    a = abs(q)
    msb = a.numerator.bit_length() - a.denominator.bit_length()
    if Fraction(2) ** msb > a:
        msb -= 1
    unit = max(msb - (p - 1), emin - (p - 1))
    scaled = a / Fraction(2) ** unit
    n = scaled.numerator // scaled.denominator
    rest = scaled - n
    if rest > Fraction(1, 2) or (rest == Fraction(1, 2) and n & 1):
        n += 1
    if n.bit_length() > p:
        n >>= 1
        unit += 1
    if unit + n.bit_length() - 1 > emax:
        return sign * Fraction(2) ** (MAX_EXP + 1)
    return sign * Fraction(n) * Fraction(2) ** unit


def exact_sqrt(q):
    """sqrt of a positive Fraction whose denominator is a power of two, as something that rounds the same."""
    k = 130 + (q.denominator.bit_length() + 1) // 2
    n = q * Fraction(4) ** k
    r = isqrt(n.numerator)
    if r * r == n.numerator:
        return Fraction(r) / Fraction(2) ** k
    # strictly between r and r + 1, which no rounding boundary of 113 bits separates
    return Fraction(2 * r + 1, 2) / Fraction(2) ** k


class ExactSampler(Sampler):
    """Inputs aimed at where correct rounding is hard."""

    def mantissa(self):
        """A 112 bit fraction: random, sparse (a tie or a sticky bit away from one), or all ones."""
        style = self.rng.random()
        if style < 0.5:
            return self.rng.getrandbits(FRAC_BITS)
        if style < 0.8:
            m = 0
            for _ in range(self.rng.randint(1, 3)):
                m |= 1 << self.rng.randint(0, FRAC_BITS - 1)
            return m
        return ((1 << FRAC_BITS) - 1) ^ self.rng.getrandbits(3)

    def value(self, e, sign=None):
        """A normal value with the exponent e, or a subnormal one below the normal range."""
        s = self.rng.randint(0, 1) if sign is None else sign
        if e < MIN_NORM_EXP:
            frac = max(self.rng.getrandbits(FRAC_BITS) >> self.rng.randint(0, FRAC_BITS - 1), 1)
            return (frac & ((1 << 64) - 1), (s << 63) | (frac >> 64))
        frac = self.mantissa()
        return (frac & ((1 << 64) - 1), (s << 63) | ((e + EXP_BIAS) << 48) | (frac >> 64))

    def exponent_of(self, bits):
        biased = (bits[1] >> 48) & EXP_MAX
        return biased - EXP_BIAS if biased else MIN_NORM_EXP


def exact_sets(s):
    """Return {name: (arity, reference, [inputs])} for the correctly rounded operations."""
    def draw(n, fn):
        return [fn() for _ in range(n)]

    q = EXACT_COUNT // 4

    def aligned():
        x = s.value(s.rng.randint(-100, 100))
        return (x, s.value(s.exponent_of(x) - s.rng.randint(0, 135)))

    def cancelling():
        x = s.value(s.rng.randint(-100, 100))
        y = decode_fraction(x) * (1 + Fraction(s.rng.getrandbits(20) + 1, 1 << s.rng.randint(60, 140)))
        return (x, encode_exact(-y))

    def tiny_pair():
        return (s.value(s.rng.randint(-16500, -16370)), s.value(s.rng.randint(-16500, -16370)))

    def big_pair():
        return (s.value(s.rng.randint(16370, 16383)), s.value(s.rng.randint(16370, 16383)))

    def subnormal_product():
        e = s.rng.randint(-9000, -7000)
        return (s.value(e), s.value(s.rng.randint(-16500, -16370) - e))

    def overflowing_product():
        return (s.value(s.rng.randint(8180, 8200)), s.value(s.rng.randint(8180, 8200)))

    def subnormal_quotient():
        e = s.rng.randint(100, 3000)
        return (s.value(s.rng.randint(-16500, -16370) + e), s.value(e))

    def subnormal_fma():
        e = s.rng.randint(-9000, -7000)
        x, y = s.value(e), s.value(s.rng.randint(-16500, -16380) - e)
        return (x, y, s.value(s.rng.randint(-16494, -16383)))

    def cancelling_fma():
        x, y = s.value(s.rng.randint(-60, 60)), s.value(s.rng.randint(-60, 60))
        p = decode_fraction(x) * decode_fraction(y)
        return (x, y, encode_exact(-p * (1 + Fraction(s.rng.getrandbits(16) + 1, 1 << s.rng.randint(100, 200)))))

    def add(x, y):
        a, b = decode_fraction(x), decode_fraction(y)
        if a + b == 0:
            return encode_exact(Fraction(0), (x[1] >> 63) & (y[1] >> 63))
        return encode_exact(a + b)

    def sub(x, y):
        return add(x, (y[0], y[1] ^ (1 << 63)))

    def mul(x, y):
        return encode_exact(decode_fraction(x) * decode_fraction(y), (x[1] ^ y[1]) >> 63)

    def div(x, y):
        return encode_exact(decode_fraction(x) / decode_fraction(y), (x[1] ^ y[1]) >> 63)

    def fma(x, y, z):
        sum_ = decode_fraction(x) * decode_fraction(y) + decode_fraction(z)
        return encode_exact(sum_)

    def to_double(x):
        return encode_exact(round_to_binary(decode_fraction(x), 53, -1022, 1023), x[1] >> 63)

    def to_float(x):
        return encode_exact(round_to_binary(decode_fraction(x), 24, -126, 127), x[1] >> 63)

    def integral(rounder):
        def ref(x):
            v = decode_fraction(x)
            r = rounder(v)
            return encode_exact(Fraction(r), x[1] >> 63)
        return ref

    def floor_(v):
        return v.numerator // v.denominator

    def ceil_(v):
        return -((-v.numerator) // v.denominator)

    def trunc_(v):
        return floor_(v) if v >= 0 else ceil_(v)

    def round_(v):
        a = abs(v) + Fraction(1, 2)
        n = a.numerator // a.denominator
        return n if v >= 0 else -n

    def rint_(v):
        a = abs(v)
        n = a.numerator // a.denominator
        rest = a - n
        if rest > Fraction(1, 2) or (rest == Fraction(1, 2) and n & 1):
            n += 1
        return n if v >= 0 else -n

    def integral_input():
        return s.value(s.rng.randint(-3, 114))

    def small_negative():
        return s.value(s.rng.randint(-4, -1), sign=1)

    sets = {
        'add_exact': (2, add, draw(EXACT_COUNT - 3 * q, aligned) + draw(q, cancelling) + draw(q, tiny_pair) + draw(q, big_pair)),
        'sub_exact': (2, sub, draw(EXACT_COUNT - 3 * q, aligned) + draw(q, cancelling) + draw(q, tiny_pair) + draw(q, big_pair)),
        'mul_exact': (2, mul, draw(EXACT_COUNT - 2 * q, lambda: (s.value(s.rng.randint(-200, 200)), s.value(s.rng.randint(-200, 200))))
                      + draw(q, subnormal_product) + draw(q, overflowing_product)),
        'div_exact': (2, div, draw(EXACT_COUNT - q, lambda: (s.value(s.rng.randint(-200, 200)), s.value(s.rng.randint(-200, 200))))
                      + draw(q, subnormal_quotient)),
        'sqrt_exact': (1, lambda x: encode_exact(exact_sqrt(decode_fraction(x))),
                       draw(EXACT_COUNT - q, lambda: s.value(s.rng.randint(-16494, 16383), sign=0))
                       + draw(q, lambda: encode_exact(Fraction(s.rng.getrandbits(56) + 1) ** 2 * Fraction(2) ** (2 * s.rng.randint(-200, 200))))),
        'fma_exact': (3, fma, draw(EXACT_COUNT // 2, subnormal_fma) + draw(EXACT_COUNT // 2, cancelling_fma)),
        'to_double_exact': (1, to_double, draw(EXACT_COUNT, lambda: s.value(s.rng.randint(-1080, 1030)))),
        'to_float_exact': (1, to_float, draw(EXACT_COUNT, lambda: s.value(s.rng.randint(-155, 130)))),
        'floor_exact': (1, integral(floor_), draw(EXACT_COUNT - q, integral_input) + draw(q, small_negative)),
        'ceil_exact': (1, integral(ceil_), draw(EXACT_COUNT - q, integral_input) + draw(q, small_negative)),
        'trunc_exact': (1, integral(trunc_), draw(EXACT_COUNT - q, integral_input) + draw(q, small_negative)),
        'round_exact': (1, integral(round_), draw(EXACT_COUNT - q, integral_input) + draw(q, small_negative)),
        'rint_exact': (1, integral(rint_), draw(EXACT_COUNT - q, integral_input) + draw(q, small_negative)),
    }
    return sets


def large_trig_sets(s):
    """Arguments far above 2^62, where the reduction has to read 2/pi to thousands of bits."""
    def draw(n, fn):
        return [fn() for _ in range(n)]

    def large():
        return s.bits(60, 16383)

    return {
        'sin_large': (mpmath.sin, draw(UNARY_COUNT, large)),
        'cos_large': (mpmath.cos, draw(UNARY_COUNT, large)),
        'tan_large': (mpmath.tan, draw(UNARY_COUNT, large)),
    }


def emit_exact(out, name, arity, ref, inputs):
    kind = {1: 'ref_unary', 2: 'ref_binary', 3: 'ref_ternary'}[arity]
    out.append(f"inline constexpr {kind} {name}_ref[] = {{")
    for args in inputs:
        args = (args,) if arity == 1 else args
        result = ref(*args)
        fields = ', '.join(f"{qword(a[0])}, {qword(a[1])}" for a in args)
        out.append(f"    {{{fields}, {qword(result[0])}, {qword(result[1])}}},")
    out.append("};")
    out.append("")

HEADER = '''// Generated by tools/gen_ref_vectors.py - do not edit by hand.
//
// Correctly rounded binary128 reference values for the float128 math functions, computed by
// mpmath at 240 bits and rounded to the format with round half to even. Each input is itself a
// binary128 bit pattern, so the only rounding anywhere in a comparison against these tables is
// the one the function under test performs.
//
// The tables named *_exact are the operations IEEE 754 requires to be correctly rounded. Their
// references are computed on exact rationals and are meant to be matched bit for bit.

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

    # The correctly rounded operations and the large trigonometric arguments draw from a sampler
    # of their own, so adding them left every table above unchanged.
    exact_sampler = ExactSampler(EXACT_SEED)
    for name, (arity, ref, inputs) in exact_sets(exact_sampler).items():
        print(f"  {name} ({len(inputs)} cases)", file=sys.stderr)
        emit_exact(out, name, arity, ref, inputs)
    with mpmath.workprec(17000):
        # sin of an argument near 2^16383 needs the argument reduced against pi to as many bits
        for name, (ref, inputs) in large_trig_sets(exact_sampler).items():
            print(f"  {name} ({len(inputs)} cases)", file=sys.stderr)
            emit_unary(out, name, ref, inputs)

    with open(args.output, 'w', encoding='utf-8', newline='\n') as f:
        f.write(HEADER)
        f.write('\n'.join(out))
        f.write(FOOTER)
    print(f"wrote {args.output}", file=sys.stderr)


if __name__ == '__main__':
    main()
