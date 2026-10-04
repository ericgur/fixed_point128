# Maintenance tools

None of these is part of the library or of a normal build. They are run by hand, when the code they
check is touched.

| | what it is |
| --- | --- |
| `gen_log2_tables.py` | generates - or verifies - the constant tables `log2()` reduces with |
| `log2_ulp_dump.cpp` + `log2_ulp_check.py` | measures the error of `log2()` in ulps against a high precision reference |
| `ieee_probe.cpp` + `ieee_check.py` | checks the correctly rounded `float128` operations bit for bit, in every rounding direction |
| `gen_two_over_pi.py` | generates - or verifies - the 2/pi table the trigonometric functions reduce large arguments with |
| `gen_ref_vectors.py` | generates the reference tables the test suite checks against, `tests/float128_ref_data.h` |

The first two are concerned with `log2()` - the one in `fixed_point128` and the one in `float128`,
which share their reduction tables. The rest concern `float128` alone.

The Python scripts need [mpmath](https://mpmath.org/):

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
Nothing is printed in decimal anywhere in between: a `fixed_point128` carries up to 126 fraction
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
series does its work, values across the exponent range, values just above one and just below it
where cancellation is worst, and every argument reduction boundary with its immediate neighbours.
Keeping them apart matters - the `near-one` and `below-one` classes are the ones that catch a `log2`
which is accurate in absolute terms but not relative to its own result, and averaging them in with
the rest hides that completely.

The two sides of one need a class each because `float128` takes a different path on each. Above one
the exponent is zero and nothing cancels; below it the exponent is -1, and adding that to the
logarithm of a mantissa just under two cancels. Before 0.10.0.0 `log2` did exactly that and lost a
bit of its result for every power of two the argument sat closer to one - half the mantissa at
2^-60, nearly all of it at 2^-110 - and with only `near-one` in the harness nothing showed it.

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

## The IEEE 754 differential check

IEEE 754 requires addition, subtraction, multiplication, division, square root, fused multiply-add
and the conversions to be correctly rounded: the result has to be the exact one, rounded once, in
the current rounding direction. A comparison within an ulp cannot tell that apart from an almost
right implementation, and a random operand rarely lands where the difference shows - next to a tie,
a sticky bit away from one, in the subnormal range. `ieee_check.py` aims its operands at exactly
those places, has `ieee_probe` run each operation, and compares every bit of the result with a
reference computed on exact rationals. The sign of a zero counts, and so does the quiet bit of a NaN.

The test suite carries a few hundred such cases per operation (`tests/float128_ieee_gtest.cpp`). This
draws tens of thousands afresh from a seed, and is the thing to run after touching the arithmetic.

```sh
cmake --build out/build/msvc/tools --config Release --target ieee_probe ieee_probe_env
python tools/ieee_check.py --probe out/build/msvc/bin/Release/ieee_probe.exe
python tools/ieee_check.py --probe out/build/msvc/bin/Release/ieee_probe_env.exe --env
```

```
ok   add:nearest                   0/2302
ok   div:nearest                   0/2502
...
ok   sqrt:downward                 0/300
```

`ieee_probe_env` is the same program built with `FP128_IEEE_ENV`; with `--env` every operation runs
in all four rounding directions, and the exception flags it raised are compared too. `--count` sets
the draws per operation and direction, `--seed` the seed. It exits non zero on any mismatch.

(The Visual Studio generator cannot build these through `cmake --build --preset ... --target`, which
looks for the project at the top of the build tree; building from the `tools` directory of the build
tree, as above, works for every generator.)

## The 2/pi table

`sin`, `cos` and `tan` reduce an argument above 2^60 by multiplying it with a 384 bit window of the
binary expansion of 2/pi, read as deep as the argument's exponent requires - about 16650 bits for one
near the top of the range. The expansion is `two_over_pi_bits` in `include/float128.h`:

```sh
python tools/gen_two_over_pi.py --check   # verify the table in the header
python tools/gen_two_over_pi.py           # print it, to paste over the array
```
