# Maintenance tools

Two things live here, both concerned with `log2()` - the one in `fixed_point128` and the one in
`float128`, which share their reduction tables. Neither is part of the library or of a normal build.

| | what it is |
| --- | --- |
| `gen_log2_tables.py` | generates - or verifies - the constant tables `log2()` reduces with |
| `log2_ulp_dump.cpp` + `log2_ulp_check.py` | measures the error of `log2()` in ulps against a high precision reference |

Both Python scripts need [mpmath](https://mpmath.org/):

```sh
pip install mpmath
```

## The constant tables

`log2()` reduces its argument using three tables in `include/fp128_shared.h`:
`log2_recip_table`, `log2_value_table` and `log2_inv_n_table`. They are generated, not hand
written, and **must not be hand edited**. The entries depend on each other: `log2_value_table[j]`
holds the logarithm of the *rounded* reciprocal stored in `log2_recip_table[j]`, not of the round
number that reciprocal approximates. That is what makes the argument reduction exact rather than
merely close, and changing either table on its own silently breaks it.

To check that the header and the generator still agree:

```sh
python tools/gen_log2_tables.py --check
```

```
ok   log2_recip_table: 64 entries match
ok   log2_value_table: 64 entries match
ok   log2_inv_n_table: 24 entries match
```

It exits non zero on a mismatch, which is the only way that mistake would otherwise come to light.
To regenerate - after changing `REDUCTION_BITS`, say - run it without `--check` and paste the
output over the three arrays in the header. `log2_reduction_bits` in the header has to be changed
to match; `--check` verifies that too.

## Measuring accuracy

`log2_ulp_dump` prints the library's answers as exact bit patterns, and `log2_ulp_check.py`
computes what they should have been at 400 bits of precision and reports the difference in ulps.
Nothing is printed in decimal anywhere in between: a `fixed_point128` carries up to 127 fraction
bits, and rendering one as decimal loses exactly the bits being measured.

Build the dump tool and run the pair:

```sh
cmake --build --preset msvc-release --target log2_ulp_dump
./out/build/msvc/bin/Release/log2_ulp_dump 4000 > dump.txt
python tools/log2_ulp_check.py dump.txt
```

```
type                        samples    mean ulp     max ulp   worst argument (high:low)
fixed_point128<1>/unit-range   2045      0.8946      1.9421   8000000000000000:F21351A3DCFBBC83
...
float128/near-one               800      0.3806      1.2091   0001000000000000:00B13AAA419A85B5
float128/boundaries             447      0.3880      1.2235   0001000000000000:0000000000000001

overall max error: 1.9421 ulp
```

The argument is the number of samples per input class per type; a few thousand takes a minute or
two, most of it in mpmath.

Rows are split by input class, because they break different things: values in `[1,2)` where the
series does its work, values across the exponent range, values just above one where cancellation is
worst, and every argument reduction boundary with its immediate neighbours. Keeping them apart
matters - the `near-one` class is the one that catches a `log2` which is accurate in absolute terms
but not relative to its own result, and averaging it in with the rest hides that completely.

One ulp always means one unit in the last place *of that result*: the fixed grid 2^-F for
`fixed_point128`, and the local grid for `float128`. So the columns are comparable across every row
even though the two types measure precision differently.

### Comparing two implementations

This is the useful mode when changing `log2()`. Dump before and after, then:

```sh
python tools/log2_ulp_check.py after.txt --baseline before.txt
```

```
   I     F   samples    mean ulp         was     max ulp         was
   1   127      8259      0.9798           -      1.9425           -
   2   126      8259      0.5779      0.5937      1.5956      1.5351
  ...
```

A change to `log2()` should not make any row worse. Note that the mean matters as much as the
maximum: the maximum is an extreme value over a finite sample and moves around by a few percent
between runs, while the mean is stable.

If you build the two dumps from different revisions, build each one fresh - and remember that on
Windows the first run of a newly linked binary is scanned by the OS, so its timings (though not its
output, which is what matters here) are worthless.

### Note on the tools' own build

`log2_ulp_dump` is excluded from the default build. Either name the target explicitly, as above, or
configure with `-DFP128_BUILD_TOOLS=ON`. It is only wired into the CMake build - the MSVC solution
in `msvc/` does not carry a project for it, since these are run occasionally and from a shell.
