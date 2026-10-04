#!/usr/bin/env python3
"""IEEE 754-2008 differential check of fp128::float128 against exact arithmetic.

Draws operands aimed at where correct rounding is hard - ties, values a sticky bit away from one,
deep cancellation, results in the subnormal range and at the edge of overflow - has ieee_probe run
each operation, and compares every bit of the result with a reference computed on exact rationals.
The sign of a zero counts, and so does the quiet bit of a NaN.

With --env the probe has to be built with FP128_IEEE_ENV. Every operation is then run in all four
rounding directions, and the exception flags each one raised are compared as well.

The committed reference tables in tests/float128_ref_data.h hold a few hundred cases per operation,
which is what the test suite can afford to run on every build. This runs tens of thousands, drawn
afresh from the seed, and is meant to be run by hand when the arithmetic is touched.

Usage:
    python tools/ieee_check.py --probe out/build/msvc/bin/Release/ieee_probe.exe
    python tools/ieee_check.py --probe out/build/msvc/bin/Release/ieee_probe_env.exe --env
"""

import argparse
import os
import random
import struct
import subprocess
import sys
import tempfile
from fractions import Fraction
from math import isqrt

try:
    import mpmath
    from mpmath import mpf
except ImportError:  # pragma: no cover - a helpful message beats a traceback
    sys.exit("This script needs mpmath: python -m pip install mpmath")

P, EMIN, EMAX, BIAS = 113, -16382, 16383, 16383
M48 = (1 << 48) - 1
M64 = (1 << 64) - 1
QUIET = 1 << 47  # the quiet bit, inside the high QWORD

NEAREST, TOWARD_ZERO, UPWARD, DOWNWARD = range(4)
MODE_NAMES = ['nearest', 'toward_zero', 'upward', 'downward']
INVALID, DIVBYZERO, OVERFLOW, UNDERFLOW, INEXACT = 1, 2, 4, 8, 16
ANY_NAN = 'nan'

INF = (0, 0x7FFF << 48)
NEG_INF = (0, (1 << 63) | (0x7FFF << 48))
QNAN = (0, (0x7FFF << 48) | QUIET)
SNAN = (1, 0x7FFF << 48)
ZERO = (0, 0)
NEG_ZERO = (0, 1 << 63)


# ---------------------------------------------------------------------------
# binary formats
# ---------------------------------------------------------------------------

def decode(bits):
    """(kind, sign, value) of a binary128 encoding; kind is zero, finite, inf, qnan or snan."""
    low, high = bits
    sign = high >> 63
    biased = (high >> 48) & 0x7FFF
    frac = ((high & M48) << 64) | low
    if biased == 0x7FFF:
        if frac == 0:
            return ('inf', sign, None)
        return ('qnan' if frac & (1 << 111) else 'snan', sign, None)
    if biased == 0:
        if frac == 0:
            return ('zero', sign, Fraction(0))
        value = Fraction(frac) * Fraction(2) ** (EMIN - 112)
    else:
        value = Fraction((1 << 112) | frac) * Fraction(2) ** (biased - BIAS - 112)
    return ('finite', sign, -value if sign else value)


def round_binary(q, p, emin, emax, mode):
    """Round a non zero Fraction to a binary format of p bits.

    Returns (kind, sign, n, lsb, flags): kind is finite, zero, inf or max (the largest finite value,
    which the directed modes overflow to), n * 2^lsb the rounded magnitude.
    """
    sign = 1 if q < 0 else 0
    a = -q if sign else q

    def rounded(lsb):
        scaled = a / Fraction(2) ** lsb
        n = scaled.numerator // scaled.denominator
        rest = scaled - n
        if rest != 0:
            if mode == NEAREST:
                n += 1 if (rest > Fraction(1, 2) or (rest == Fraction(1, 2) and n & 1)) else 0
            elif mode == UPWARD:
                n += 0 if sign else 1
            elif mode == DOWNWARD:
                n += 1 if sign else 0
        return n, rest != 0

    e = a.numerator.bit_length() - a.denominator.bit_length()
    if Fraction(2) ** e > a:
        e -= 1
    # Tininess is detected after rounding: rounded with an unbounded exponent, is it below 2^emin?
    unbounded, _ = rounded(e - (p - 1))
    tiny = e + (1 if unbounded == 1 << p else 0) < emin

    lsb = max(e - (p - 1), emin - (p - 1))
    n, inexact = rounded(lsb)
    if n == 1 << p:
        n >>= 1
        lsb += 1
    if n and lsb + n.bit_length() - 1 > emax:
        to_inf = mode == NEAREST or (mode == UPWARD and not sign) or (mode == DOWNWARD and sign)
        return ('inf' if to_inf else 'max', sign, 0, 0, OVERFLOW | INEXACT)
    flags = (INEXACT | (UNDERFLOW if tiny else 0)) if inexact else 0
    return ('finite' if n else 'zero', sign, n, lsb, flags)


def encode128(q, mode, zero_sign=0):
    """The binary128 encoding of an exact Fraction rounded in a direction, and the flags raised."""
    if q == 0:
        return ((0, zero_sign << 63), 0)
    kind, sign, n, lsb, flags = round_binary(q, P, EMIN, EMAX, mode)
    if kind == 'zero':
        return ((0, sign << 63), flags)
    if kind == 'inf':
        return ((0, (sign << 63) | (0x7FFF << 48)), flags)
    if kind == 'max':
        return ((M64, (sign << 63) | (0x7FFE << 48) | M48), flags)
    if n < (1 << 112):
        biased, frac = 0, n
    else:
        biased, frac = lsb + 112 + BIAS, n - (1 << 112)
    return ((frac & M64, (sign << 63) | (biased << 48) | (frac >> 64)), flags)


def encode_narrow(q, p, emin, emax, bias, ebits, mode):
    """The encoding of an exact Fraction rounded to a narrower interchange format."""
    total = 1 + ebits + p - 1
    if q == 0:
        return (0, 0)
    kind, sign, n, lsb, flags = round_binary(q, p, emin, emax, mode)
    top = sign << (total - 1)
    if kind == 'zero':
        return (top, flags)
    if kind == 'inf':
        return (top | (((1 << ebits) - 1) << (p - 1)), flags)
    if kind == 'max':
        return (top | (((1 << ebits) - 2) << (p - 1)) | ((1 << (p - 1)) - 1), flags)
    if n < (1 << (p - 1)):
        return (top | n, flags)
    return (top | ((lsb + (p - 1) + bias) << (p - 1)) | (n - (1 << (p - 1))), flags)


def make(sign, e, frac):
    return (frac & M64, (sign << 63) | ((e + BIAS) << 48) | (frac >> 64))


def neg(bits):
    return (bits[0], bits[1] ^ (1 << 63))


def is_nan(bits):
    return decode(bits)[0] in ('qnan', 'snan')


# ---------------------------------------------------------------------------
# references
# ---------------------------------------------------------------------------

def nan_result(*operands):
    """A NaN operand propagates: a quiet NaN out, invalid when any operand is signaling."""
    return (ANY_NAN, INVALID if any(decode(x)[0] == 'snan' for x in operands) else 0)


def arithmetic(op, a, b, c, mode):
    da, db, dc = decode(a), decode(b), decode(c)
    arity = {'add': 2, 'sub': 2, 'mul': 2, 'div': 2, 'fma': 3, 'sqrt': 1, 'sqr': 1}[op]
    operands = (a, b, c)[:arity]
    if any(is_nan(x) for x in operands):
        return nan_result(*operands)
    exact_zero_sign = 1 if mode == DOWNWARD else 0
    if op in ('add', 'sub'):
        sb = db[1] ^ (1 if op == 'sub' else 0)
        if da[0] == 'inf' or db[0] == 'inf':
            if da[0] == 'inf' and db[0] == 'inf' and da[1] != sb:
                return (ANY_NAN, INVALID)
            return ((0, ((da[1] if da[0] == 'inf' else sb) << 63) | (0x7FFF << 48)), 0)
        total = da[2] + (-db[2] if op == 'sub' else db[2])
        if total == 0:
            same = da[0] == 'zero' and db[0] == 'zero' and da[1] == sb
            return ((0, (da[1] if same else exact_zero_sign) << 63), 0)
        return encode128(total, mode)
    if op in ('mul', 'sqr'):
        if op == 'sqr':
            db = da
        sign = da[1] ^ db[1]
        if da[0] == 'inf' or db[0] == 'inf':
            if da[0] == 'zero' or db[0] == 'zero':
                return (ANY_NAN, INVALID)
            return ((0, (sign << 63) | (0x7FFF << 48)), 0)
        return encode128(da[2] * db[2], mode, sign)
    if op == 'div':
        sign = da[1] ^ db[1]
        if (da[0] == 'inf' and db[0] == 'inf') or (da[0] == 'zero' and db[0] == 'zero'):
            return (ANY_NAN, INVALID)
        if da[0] == 'inf' or db[0] == 'zero':
            return ((0, (sign << 63) | (0x7FFF << 48)), DIVBYZERO if db[0] == 'zero' and da[0] != 'inf' else 0)
        if da[0] == 'zero' or db[0] == 'inf':
            return ((0, sign << 63), 0)
        return encode128(da[2] / db[2], mode, sign)
    if op == 'sqrt':
        if da[0] == 'zero':
            return (a, 0)
        if da[1]:
            return (ANY_NAN, INVALID)
        if da[0] == 'inf':
            return (a, 0)
        k = 130 + (da[2].denominator.bit_length() + 1) // 2
        n = da[2] * Fraction(4) ** k
        r = isqrt(n.numerator)
        # strictly between r and r + 1 when inexact, where no rounding boundary separates them
        root = Fraction(r) if r * r == n.numerator else Fraction(2 * r + 1, 2)
        return encode128(root / Fraction(2) ** k, mode)
    if op == 'fma':
        sign = da[1] ^ db[1]
        if (da[0] == 'inf' and db[0] == 'zero') or (da[0] == 'zero' and db[0] == 'inf'):
            return (ANY_NAN, INVALID)
        if da[0] == 'inf' or db[0] == 'inf':
            if dc[0] == 'inf' and dc[1] != sign:
                return (ANY_NAN, INVALID)
            return ((0, (sign << 63) | (0x7FFF << 48)), 0)
        if dc[0] == 'inf':
            return (c, 0)
        product = da[2] * db[2]
        total = product + dc[2]
        if total == 0:
            # Zeros of the same sign keep it; anything else is an exact zero sum of opposite signs.
            if product == 0 and dc[0] == 'zero' and sign == dc[1]:
                return ((0, sign << 63), 0)
            return ((0, exact_zero_sign << 63), 0)
        return encode128(total, mode)
    raise ValueError(op)


def integral(op, a, mode):
    d = decode(a)
    if d[0] in ('qnan', 'snan'):
        return nan_result(a)
    if d[0] in ('inf', 'zero'):
        return (a, 0)
    v = abs(d[2])
    n = v.numerator // v.denominator
    rest = v - n
    away = (mode == UPWARD and not d[1]) or (mode == DOWNWARD and d[1])
    if rest:
        if op == 'round':
            n += 1 if rest >= Fraction(1, 2) else 0
        elif op == 'floor':
            n += 1 if d[1] else 0
        elif op == 'ceil':
            n += 0 if d[1] else 1
        elif op == 'trunc':
            pass
        elif mode == NEAREST:
            n += 1 if (rest > Fraction(1, 2) or (rest == Fraction(1, 2) and n & 1)) else 0
        elif away:
            n += 1
    bits, _ = encode128(Fraction(-n if d[1] else n), NEAREST, d[1])
    return (bits, INEXACT if (rest and op == 'rint') else 0)


def remainders(op, a, b):
    da, db = decode(a), decode(b)
    if is_nan(a) or is_nan(b):
        return nan_result(a, b)
    if da[0] == 'inf' or db[0] == 'zero':
        return (ANY_NAN, INVALID)
    if db[0] == 'inf' or da[0] == 'zero':
        return (a, 0)
    x, y = da[2], db[2]
    q = x / y
    if op == 'fmod':
        n = int(q)
    else:
        n = abs(q).numerator // abs(q).denominator
        rest = abs(q) - n
        if rest > Fraction(1, 2) or (rest == Fraction(1, 2) and n & 1):
            n += 1
        n = n if q >= 0 else -n
    r = x - n * y
    bits, _ = encode128(r, NEAREST, da[1])
    return (bits, 0)


def reference(op, a, b, c, mode):
    """(expected bits, expected flags) of an operation; expected bits ANY_NAN for a quiet NaN."""
    if op in ('add', 'sub', 'mul', 'div', 'sqrt', 'fma', 'sqr'):
        return arithmetic(op, a, b, c, mode)
    if op in ('rint', 'nearbyint', 'floor', 'ceil', 'trunc', 'round'):
        return integral(op, a, mode)
    if op in ('remainder', 'fmod'):
        return remainders(op, a, b)
    if op == 'ldexp':
        n = b[0] - (1 << 64) if b[0] >> 63 else b[0]
        return encode128(decode(a)[2] * Fraction(2) ** n, mode)
    if op == 'todouble':
        bits, flags = encode_narrow(decode(a)[2], 53, -1022, 1023, 1023, 11, mode)
        return ((bits, 0), flags)
    if op == 'tofloat':
        bits, flags = encode_narrow(decode(a)[2], 24, -126, 127, 127, 8, mode)
        return ((bits, 0), flags)
    if op in ('toint64', 'toint32'):
        width = 64 if op == 'toint64' else 32
        return (((int(decode(a)[2]) & ((1 << width) - 1)), 0), 0)
    raise ValueError(op)


# ---------------------------------------------------------------------------
# inputs
# ---------------------------------------------------------------------------

class Sampler:
    def __init__(self, seed):
        self.rng = random.Random(seed)

    def mantissa(self):
        """Random, sparse (a tie or a sticky bit away from one), or all ones."""
        style = self.rng.random()
        if style < 0.5:
            return self.rng.getrandbits(112)
        if style < 0.8:
            m = 0
            for _ in range(self.rng.randint(1, 3)):
                m |= 1 << self.rng.randint(0, 111)
            return m
        return ((1 << 112) - 1) ^ self.rng.getrandbits(3)

    def value(self, lo, hi, sign=None):
        s = self.rng.getrandbits(1) if sign is None else sign
        e = self.rng.randint(lo, hi)
        if e < EMIN:
            frac = max(self.rng.getrandbits(112) >> self.rng.randint(0, 111), 1)
            return (frac & M64, (s << 63) | (frac >> 64))
        return make(s, e, self.mantissa())

    @staticmethod
    def exponent(bits):
        biased = (bits[1] >> 48) & 0x7FFF
        return biased - BIAS if biased else EMIN


def operations(s, count):
    """(op, a, b, c) for every operation, hard cases first."""
    v = s.value
    cases = []
    for _ in range(count):
        a = v(-200, 200)
        b = make(s.rng.getrandbits(1), s.exponent(a) - s.rng.randint(0, 135), s.mantissa())
        cases += [('add', a, b, ZERO), ('sub', a, b, ZERO)]
        cases += [('mul', v(-300, 300), v(-300, 300), ZERO), ('div', v(-300, 300), v(-300, 300), ZERO)]
        cases.append(('sqrt', v(-16494, 16383, sign=0), ZERO, ZERO))
        x, y = v(-100, 100), v(-100, 100)
        z = make(s.rng.getrandbits(1), s.exponent(x) + s.exponent(y) + s.rng.randint(-120, 120), s.mantissa())
        cases.append(('fma', x, y, z))
        cases.append(('todouble', v(-1100, 1030), ZERO, ZERO))
        cases.append(('tofloat', v(-160, 130), ZERO, ZERO))
        r = v(-4, 116)
        cases += [(op, r, ZERO, ZERO) for op in ('rint', 'nearbyint', 'floor', 'ceil', 'trunc', 'round')]
        cases += [('remainder', v(-200, 200), v(-200, 200), ZERO), ('fmod', v(-200, 200), v(-200, 200), ZERO)]
    for _ in range(count // 4):
        # subnormal results, overflow, sums of tiny values
        e = s.rng.randint(-9000, -7000)
        cases.append(('mul', v(e, e), v(-16500 - e, -16370 - e), ZERO))
        cases.append(('mul', v(8180, 8200), v(8180, 8200), ZERO))
        e = s.rng.randint(100, 3000)
        cases.append(('div', v(-16500 + e, -16370 + e), v(e, e), ZERO))
        cases.append(('add', v(-16500, -16370), v(-16500, -16370), ZERO))
        cases.append(('sqr', v(-8250, -8190), ZERO, ZERO))
        e = s.rng.randint(-9000, -7000)
        cases.append(('fma', v(e, e), v(-16500 - e, -16380 - e), v(-16494, -16383)))
        cases.append(('ldexp', v(-16382, -16000), (s.rng.randint(-150, -1) & M64, 0), ZERO))
        cases.append(('toint64', make(s.rng.getrandbits(1), s.rng.randint(-3, 61), s.mantissa()), ZERO, ZERO))
        cases.append(('toint32', make(s.rng.getrandbits(1), s.rng.randint(-3, 29), s.mantissa()), ZERO, ZERO))
    for op in ('add', 'sub', 'mul', 'div', 'fma'):
        cases += [(op, SNAN, make(0, 0, 0), make(0, 0, 0)), (op, make(0, 0, 0), QNAN, make(0, 0, 0))]
    cases += [('add', INF, NEG_INF, ZERO), ('mul', INF, ZERO, ZERO), ('div', ZERO, ZERO, ZERO), ('div', make(0, 0, 0), NEG_ZERO, ZERO),
              ('sqrt', make(1, 0, 0), ZERO, ZERO), ('sqrt', SNAN, ZERO, ZERO), ('sub', make(0, 0, 0), make(0, 0, 0), ZERO),
              ('fma', make(0, 0, 0), make(0, 0, 0), make(1, 0, 0), )]
    return cases


def special_values():
    """IEEE 754-2008 9.2.1 and friends: (op, a, b, expected). Checked to nearest only."""
    one, two, half, mone = make(0, 0, 0), make(0, 1, 0), make(0, -1, 0), make(1, 0, 0)
    pi = 'pi'
    return [
        ('log', mone, ZERO, ANY_NAN), ('log', NEG_INF, ZERO, ANY_NAN), ('log', INF, ZERO, INF), ('log', NEG_ZERO, ZERO, NEG_INF),
        ('log', one, ZERO, ZERO), ('log2', INF, ZERO, INF), ('log2', QNAN, ZERO, ANY_NAN), ('log10', mone, ZERO, ANY_NAN),
        ('log1p', NEG_ZERO, ZERO, NEG_ZERO), ('exp', NEG_INF, ZERO, ZERO), ('exp', SNAN, ZERO, ANY_NAN), ('expm1', NEG_INF, ZERO, mone),
        ('sin', NEG_ZERO, ZERO, NEG_ZERO), ('sin', INF, ZERO, ANY_NAN), ('cos', INF, ZERO, ANY_NAN), ('tan', NEG_ZERO, ZERO, NEG_ZERO),
        ('asin', two, ZERO, ANY_NAN), ('acos', one, ZERO, ZERO), ('acos', two, ZERO, ANY_NAN), ('atan', NEG_ZERO, ZERO, NEG_ZERO),
        ('atanh', one, ZERO, INF), ('atanh', two, ZERO, ANY_NAN), ('cosh', NEG_INF, ZERO, INF), ('tanh', NEG_INF, ZERO, mone),
        ('logb', ZERO, ZERO, NEG_INF), ('nextup', SNAN, ZERO, ANY_NAN), ('fmin', SNAN, one, ANY_NAN), ('fmin', QNAN, one, one),
        ('atan2', ZERO, ZERO, ZERO), ('atan2', NEG_ZERO, ZERO, NEG_ZERO), ('atan2', ZERO, NEG_ZERO, pi), ('atan2', ZERO, one, ZERO),
        ('atan2', ZERO, mone, pi), ('atan2', INF, INF, 'pi/4'),
        ('pow', QNAN, ZERO, one), ('pow', one, QNAN, one), ('pow', mone, INF, one), ('pow', half, INF, ZERO), ('pow', two, NEG_INF, ZERO),
        ('pow', two, QNAN, ANY_NAN), ('pow', NEG_ZERO, half, ZERO), ('pow', NEG_ZERO, make(1, 1, 1 << 111), NEG_INF),
        ('pow', NEG_INF, half, INF), ('pow', half, make(0, 14, 0x3880 << 96), ZERO),
        ('hypot', INF, QNAN, INF), ('hypot', QNAN, NEG_INF, INF),
    ]


def parse_cases():
    """(text, expected) for the conversion from character sequences."""
    # 40 significant digits just above the halfway point between 1 and its neighbour: exactly what
    # the decimal conversion has to round correctly (its limit is 40 digits; IEEE 754 asks for 39).
    with mpmath.workdps(400):
        above = mpmath.nstr(mpf(1) + mpf(2) ** -113 + mpf(2) ** -300, 40, strip_zeros=False)
    return [
        ('0x1.8p1', encode128(Fraction(3), NEAREST)[0]),
        ('0x1fffffffffffffffffffffffffffffff', encode128(Fraction(2 ** 125), NEAREST)[0]),
        ('0x1.00000000000000000000000000008p0', make(0, 0, 0)),
        ('snan', 'snan'),
        ('nan(0x7b)', (0x7B, (0x7FFF << 48) | QUIET)),
        ('3.3e-4966', (1, 0)),
        ('3.2e-4966', ZERO),
        ('-0', NEG_ZERO),
        ('1e5000', INF),
        ('-1e-5000', NEG_ZERO),
        (above, encode128(Fraction(above), NEAREST)[0]),
    ]


# ---------------------------------------------------------------------------
# running
# ---------------------------------------------------------------------------

def run_probe(probe, lines):
    with tempfile.TemporaryDirectory() as tmp:
        cases_path = os.path.join(tmp, 'cases.txt')
        out_path = os.path.join(tmp, 'out.txt')
        with open(cases_path, 'w') as f:
            f.write('\n'.join(lines) + '\n')
        subprocess.run([probe, cases_path, out_path], check=True)
        with open(out_path) as f:
            return [line.split() for line in f]


def case_line(op, a, b, c):
    return f'{op} {a[0]:x} {a[1]:x} {b[0]:x} {b[1]:x} {c[0]:x} {c[1]:x}'


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--probe', required=True, help='path to the ieee_probe executable')
    parser.add_argument('--env', action='store_true', help='the probe has FP128_IEEE_ENV: check every direction and the flags')
    parser.add_argument('--count', type=int, default=2000, help='draws per operation and direction (default 2000)')
    parser.add_argument('--seed', type=int, default=20261004)
    args = parser.parse_args()

    sampler = Sampler(args.seed)
    modes = range(4) if args.env else (NEAREST,)
    plan = []
    for mode in modes:
        plan.append(('setmode', mode, None))
        for case in operations(sampler, args.count):
            plan.append(('op', mode, case))
    plan.append(('setmode', NEAREST, None))
    specials = special_values()
    parses = parse_cases()

    lines = []
    for kind, mode, case in plan:
        lines.append(case_line('setmode', (mode, 0), ZERO, ZERO) if kind == 'setmode' else case_line(*case))
    lines += [case_line(op, a, b, ZERO) for op, a, b, _ in specials]
    lines += [f'parse {text}' for text, _ in parses]
    rows = run_probe(os.path.abspath(args.probe), lines)
    assert len(rows) == len(lines), (len(rows), len(lines))

    stats, examples = {}, {}

    def record(key, ok, message):
        total_bad = stats.setdefault(key, [0, 0])
        total_bad[0] += 1
        if not ok:
            total_bad[1] += 1
            if len(examples.setdefault(key, [])) < 3:
                examples[key].append(message)

    for (kind, mode, case), row in zip(plan, rows):
        if kind == 'setmode':
            continue
        op, a, b, c = case
        got = (int(row[0], 16), int(row[1], 16))
        got_flags = int(row[2], 16)
        want, want_flags = reference(op, a, b, c, mode)
        if want == ANY_NAN:
            ok = decode(got)[0] == 'qnan'
        else:
            ok = got == want
        if args.env:
            ok = ok and got_flags == want_flags
        record(f'{op}:{MODE_NAMES[mode]}', ok,
               f'in=({a[1]:016x}{a[0]:016x}, {b[1]:016x}{b[0]:016x}, {c[1]:016x}{c[0]:016x}) '
               f'got={got[1]:016x}{got[0]:016x}/{got_flags} want={want if want == ANY_NAN else f"{want[1]:016x}{want[0]:016x}"}/{want_flags}')

    offset = len(plan)
    for (op, a, b, want), row in zip(specials, rows[offset:]):
        got = (int(row[0], 16), int(row[1], 16))
        if want == ANY_NAN:
            ok = decode(got)[0] == 'qnan'
        elif isinstance(want, str):
            with mpmath.workprec(400):
                value = mpmath.pi if want == 'pi' else mpmath.pi / 4
                man, exp = mpf(value).man_exp
            ok = got == encode128(Fraction(int(man)) * Fraction(2) ** int(exp), NEAREST)[0]
        else:
            ok = got == want
        record(f'special:{op}', ok, f'{op}({a[1]:016x}{a[0]:016x}, {b[1]:016x}{b[0]:016x}) got={got[1]:016x}{got[0]:016x}')

    offset += len(specials)
    for (text, want), row in zip(parses, rows[offset:]):
        got = (int(row[0], 16), int(row[1], 16))
        ok = decode(got)[0] == 'snan' if want == 'snan' else got == want
        record('parse', ok, f'"{text[:50]}" got={got[1]:016x}{got[0]:016x}')

    failed = 0
    for key in sorted(stats):
        total, bad = stats[key]
        failed += bad
        print(f'{"FAIL" if bad else "ok  "} {key:28s} {bad:6d}/{total}')
        for message in examples.get(key, []):
            print('        ' + message)
    return 1 if failed else 0


if __name__ == '__main__':
    sys.exit(main())
