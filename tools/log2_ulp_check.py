#!/usr/bin/env python3
"""Reports the error of log2 in ulps, against an mpmath reference.

The other half of the harness; log2_ulp_dump reproduces the library's answers as exact bit
patterns and this computes what they should have been. See tools/README.md.

Usage:
    log2_ulp_dump 4000 > dump.txt
    python tools/log2_ulp_check.py dump.txt
    python tools/log2_ulp_check.py dump.txt --baseline old.txt   # compare two implementations

Each dump line is a label followed by two values, the argument and log2 of it:

    <label> <sign> <mant_high> <mant_low> <exp>  <sign> <mant_high> <mant_low> <exp>

standing for (-1)^sign * (mant_high * 2^64 + mant_low) * 2^exp. One ulp of a result is 2^exp of
that result, which is the fixed grid spacing for fixed_point128 and the local grid spacing for
float128, so the figures are comparable across both.

Requires mpmath (pip install mpmath). The reference is computed at 400 bits, far beyond the 128
under test, so the reference itself contributes nothing to the numbers reported.
"""
import argparse
import sys
from collections import defaultdict

try:
    from mpmath import mp, mpf, log
except ImportError:  # pragma: no cover - depends on the environment
    sys.exit("this script needs mpmath: pip install mpmath")

mp.prec = 400


def value(sign, mant_high, mant_low, exponent):
    """Rebuilds the exact value a dumped quadruple stands for."""
    magnitude = mpf((mant_high << 64) | mant_low) * mpf(2) ** exponent
    return -magnitude if sign else magnitude


def measure(path):
    """Returns {label: (samples, mean ulp, max ulp, worst argument)}, preserving first-seen order."""
    worst = {}
    total = defaultdict(float)
    count = defaultdict(int)

    with open(path, encoding="utf-8-sig") as source:
        for line in source:
            parts = line.split()
            if len(parts) != 9:
                continue
            label = parts[0]
            x = value(int(parts[1]), int(parts[2], 16), int(parts[3], 16), int(parts[4]))
            got = value(int(parts[5]), int(parts[6], 16), int(parts[7], 16), int(parts[8]))
            if x <= 0:
                continue

            ulp = mpf(2) ** int(parts[8])
            ulps = float(abs(got - log(x) / log(mpf(2))) / ulp)
            total[label] += ulps
            count[label] += 1
            if label not in worst or ulps > worst[label][0]:
                worst[label] = (ulps, f"{parts[2]}:{parts[3]}")

    return {k: (count[k], total[k] / count[k], worst[k][0], worst[k][1]) for k in count}


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("dump", help="output of log2_ulp_dump")
    parser.add_argument("--baseline", help="a second dump to compare against, e.g. from before a change")
    args = parser.parse_args()

    results = measure(args.dump)
    baseline = measure(args.baseline) if args.baseline else {}

    width = max(len(k) for k in results) if results else 20
    if baseline:
        print(f"{'type':<{width}} {'samples':>9} {'mean ulp':>11} {'was':>11} {'max ulp':>11} {'was':>11}")
        for label, (samples, mean, worst, _) in results.items():
            old = baseline.get(label)
            old_mean = f"{old[1]:.4f}" if old else "-"
            old_max = f"{old[2]:.4f}" if old else "-"
            print(f"{label:<{width}} {samples:>9} {mean:>11.4f} {old_mean:>11} {worst:>11.4f} {old_max:>11}")
    else:
        print(f"{'type':<{width}} {'samples':>9} {'mean ulp':>11} {'max ulp':>11}   worst argument (high:low)")
        for label, (samples, mean, worst, arg) in results.items():
            print(f"{label:<{width}} {samples:>9} {mean:>11.4f} {worst:>11.4f}   {arg}")

    overall = max(worst for _, _, worst, _ in results.values())
    print(f"\noverall max error: {overall:.4f} ulp")
    if baseline:
        was = max(worst for _, _, worst, _ in baseline.values())
        print(f"baseline max error: {was:.4f} ulp")

    return 0


if __name__ == "__main__":
    sys.exit(main())
