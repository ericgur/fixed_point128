#!/usr/bin/env python3
"""Generate - or verify - the 2/pi table float128's trigonometric argument reduction uses.

sin, cos and tan of a large argument need x * 2/pi modulo 4 to well over a hundred bits below the
binary point, and for an argument near the top of the binary128 range that takes about sixteen
thousand bits of 2/pi: the bits above the window a given argument reads only add multiples of four,
the bits below it are past the precision the reduction keeps. The table holds the binary expansion
of 2/pi after the binary point, most significant word first.

Usage:
    python tools/gen_two_over_pi.py           # print the table to paste into include/float128.h
    python tools/gen_two_over_pi.py --check   # verify the table in the header

The check exits non zero on a mismatch.
"""

import argparse
import re
import sys

try:
    from mpmath import mp, mpf
except ImportError:  # pragma: no cover - a helpful message beats a traceback
    sys.exit("This script needs mpmath: python -m pip install mpmath")

# The reduction reads WINDOW_BITS bits starting at most at bit 16383 - 113 of the expansion, so the
# table has to reach past bit 16270 + 384. The extra words leave room for the window to change.
WORDS = 264
HEADER = 'include/float128.h'
ARRAY = 'two_over_pi_bits'


def words():
    mp.prec = WORDS * 64 + 128
    scaled = int(mp.floor((mpf(2) / mp.pi) * (mpf(2) ** (WORDS * 64))))
    return [(scaled >> (64 * (WORDS - 1 - i))) & ((1 << 64) - 1) for i in range(WORDS)]


def render(values):
    lines = []
    for i in range(0, len(values), 4):
        lines.append('    ' + ', '.join(f'0x{v:016X}' for v in values[i:i + 4]) + ',')
    return '\n'.join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--check', action='store_true', help=f'verify the table in {HEADER}')
    args = parser.parse_args()

    values = words()
    if not args.check:
        print(render(values))
        return 0

    text = open(HEADER, encoding='utf-8').read()
    match = re.search(ARRAY + r'\[\]\s*=\s*\{(.*?)\};', text, re.S)
    if match is None:
        print(f'FAIL {ARRAY}: not found in {HEADER}')
        return 1
    found = [int(v, 16) for v in re.findall(r'0x([0-9A-Fa-f]{16})', match.group(1))]
    if found != values:
        print(f'FAIL {ARRAY}: {len(found)} entries, expected {len(values)} matching ones')
        return 1
    print(f'ok   {ARRAY}: {len(values)} entries match')
    return 0


if __name__ == '__main__':
    sys.exit(main())
