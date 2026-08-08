#!/usr/bin/env python3
"""Generates - or verifies - the constant tables fixed_point128::log2() reduces with.

The tables live in include/fp128_shared.h. They are not hand written and must not be
hand edited: every entry depends on the others, and log2_value_table in particular is derived from
the *rounded* reciprocal stored next to it rather than from the exact value that reciprocal
approximates. That is what makes the argument reduction exact, and editing either table on its own
would silently break it.

Usage:
    python tools/gen_log2_tables.py            # print the tables, to paste into the header
    python tools/gen_log2_tables.py --check    # verify the header matches what this would emit

--check is the useful one day to day: it exits non zero if the tables and this generator have
drifted apart, which is the only way that mistake would otherwise be noticed.

Requires mpmath (pip install mpmath). Everything is computed at 600 bits, far beyond the 128 the
tables hold, so the rounding of each entry is the only error in them.

Three scalings are emitted, each chosen so the consumer needs no conversion at runtime:

  log2_recip_table  raw 128 bit fraction, round(v * 2^128) for v in [0,1). The reduction multiply
                    treats both operands as pure fractions, so this form is independent of the
                    fixed_point128 template parameter.
  log2_value_table  raw 128 bit fraction as well. One shift per call converts it to the caller's
                    scaling.
  log2_inv_n_table  already in fixed_point128<1> form, round(v * 2^127). The series reads one entry
                    per iteration, so a per iteration shift would cost more than the multiply the
                    entry is used for.
"""
import argparse
import pathlib
import re
import sys

try:
    from mpmath import mp, mpf, log
except ImportError:  # pragma: no cover - depends on the environment
    sys.exit("this script needs mpmath: pip install mpmath")

mp.prec = 600

REDUCTION_BITS = 6                      # must match log2_reduction_bits in the header
ENTRIES = 1 << REDUCTION_BITS
MAX_TERMS = 24                          # enough for F = 127 plus guard bits
TWO128 = mpf(2) ** 128
TWO127 = mpf(2) ** 127

HEADER = pathlib.Path(__file__).resolve().parent.parent / "include" / "fp128_shared.h"


def split(scaled):
    """Splits a 128 bit integer into (high, low) QWORDs."""
    if not 0 <= scaled < (1 << 128):
        raise ValueError(f"value does not fit 128 bits: {scaled}")
    return (scaled >> 64) & 0xFFFFFFFFFFFFFFFF, scaled & 0xFFFFFFFFFFFFFFFF


def as_fraction128(value):
    """v in [0,1) -> round(v * 2^128), the raw 128 bit fraction form."""
    return split(int(mp.nint(value * TWO128)))


def as_fp1(value):
    """v in [0,2) -> round(v * 2^127), the raw form of fixed_point128<1>."""
    return split(int(mp.nint(value * TWO127)))


def build():
    """Returns {table name: [(high, low), ...]} plus the widest |z| the reduction leaves."""
    ln2 = log(mpf(2))
    recip, value = [], []
    max_z = mpf(0)

    for j in range(ENTRIES):
        c = 1 + mpf(j) / ENTRIES
        if j == 0:
            # 1.0 is not a 128 bit fraction. The largest value below it makes the reduction a near
            # no-op instead of an exact one, and the matching table entry absorbs the difference.
            high, low = 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF
        else:
            high, low = as_fraction128(1 / c)
        recip.append((high, low))

        stored = mpf((high << 64) | low) / TWO128
        # Derived from the stored reciprocal, so that log2(m) = log2(m * stored) - log2(stored)
        # holds exactly whatever the reciprocal rounded to.
        t_value = -log(stored) / ln2
        if not 0 <= t_value < 1:
            raise ValueError(f"log2_value_table[{j}] outside [0,1): {t_value}")
        value.append(as_fraction128(t_value))

        for m in (c, c + mpf(1) / ENTRIES):
            max_z = max(max_z, abs(m * stored - 1))

    # 1/(n*ln2): the division by ln(2) is folded into the series rather than applied at the end,
    # which removes both a multiply and a rounding step from every call. Entry 0 is 1/ln2 = 1.4427,
    # which still fits fixed_point128<1>, and the accumulator peaks around 1.47.
    inv_n = [as_fp1(1 / (n * ln2)) for n in range(1, MAX_TERMS + 1)]

    return {"log2_recip_table": recip, "log2_value_table": value, "log2_inv_n_table": inv_n}, max_z


def render(name, pairs):
    lines = [f"inline constexpr uint64_t {name}[][2] = {{"]
    lines += [f"    {{0x{high:016X}ull, 0x{low:016X}ull}}," for high, low in pairs]
    lines.append("};")
    return "\n".join(lines)


def parse_header(text, name):
    """Extracts the (high, low) pairs of one table from the header source."""
    match = re.search(rf"inline constexpr uint64_t {name}\[\]\[2\] = \{{(.*?)\n\}};", text, re.S)
    if not match:
        return None
    return [(int(h, 16), int(l, 16))
            for h, l in re.findall(r"\{0x([0-9A-Fa-f]+)ull, 0x([0-9A-Fa-f]+)ull\}", match.group(1))]


def check(tables):
    text = HEADER.read_text(encoding="utf-8")

    bits = re.search(r"log2_reduction_bits = (\d+)", text)
    if not bits or int(bits.group(1)) != REDUCTION_BITS:
        print(f"FAIL log2_reduction_bits: header has {bits.group(1) if bits else '?'}, "
              f"this generator assumes {REDUCTION_BITS}")
        return 1

    failures = 0
    for name, expected in tables.items():
        found = parse_header(text, name)
        if found is None:
            print(f"FAIL {name}: not found in {HEADER.name}")
            failures += 1
        elif found != expected:
            print(f"FAIL {name}: header has {len(found)} entries, generator produced "
                  f"{len(expected)}" if len(found) != len(expected) else
                  f"FAIL {name}: {sum(a != b for a, b in zip(found, expected))} entries differ")
            failures += 1
        else:
            print(f"ok   {name}: {len(found)} entries match")

    return 1 if failures else 0


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--check", action="store_true",
                        help="verify the header matches this generator instead of printing")
    args = parser.parse_args()

    tables, max_z = build()
    if args.check:
        return check(tables)

    ln2 = log(mpf(2))
    print(f"// generated by tools/gen_log2_tables.py - {REDUCTION_BITS} reduction bits, {ENTRIES} entries")
    print(f"// widest |z| the reduction leaves: 2^{mp.nstr(log(max_z) / ln2, 6)}")
    print()
    for name, pairs in tables.items():
        print(render(name, pairs))
        print()

    return 0


if __name__ == "__main__":
    sys.exit(main())
