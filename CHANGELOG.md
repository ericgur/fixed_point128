# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to a four-part `MAJOR.MINOR.PATCH.BUILD` version scheme
exposed through the `FP128_VERSION*` macros and the `fp128::version*` constants in
[`include/fp128_shared.h`](include/fp128_shared.h).

## [0.11.0.0] - unreleased

### Added

- **Minimum integer bits are checked at compile time.** A constant or function that cannot work in
  a narrow instantiation now rejects it with a `static_assert` naming itself and the bound, instead
  of overflowing quietly or failing to link:

  | Needs | Constants and functions |
  |---|---|
  | `I >= 2` | `pi()`, `e()`, `exp()`, `exp2()`, `expm1()`, `pow()`, `tanh()` |
  | `I >= 3` | `pi2()` |
  | `I >= 4` | `sin1()`, `cos1()`, `sin()`, `cos()`, `tan()`, `asin()`, `acos()`, `atan()`, `atan2()`, `sinh()`, `cosh()` |

  Each bound is the smallest `I` that works rather than a safe margin, and `tests/fixed_point128_gtest.cpp`
  pins that down by evaluating each one at exactly its minimum.

  The seven functions that already carried `requires (I >= 4)` keep the same bound but now diagnose
  properly. The constraint never reached the caller: the namespace scope forward declarations that
  shadow the CRT are unconstrained, so a violation bound to one of those and surfaced as
  `error LNK2019: unresolved external symbol "fp128::sin<3>"` at link time. Constraining those
  declarations as well would have been worse - `fixed_point128` converts implicitly to `double`, so
  the call would have silently reached `::sin(double)`. A `static_assert` in the body leaves the
  friend the best match, so it is selected, instantiated, and rejected at the call site.

- `asin()`, `acos()`, `atan()` and `atan2()` had no bound at all despite evaluating `sin()` and
  `cos()`; below `I` of 4 they were the same link error.

- **`div_128bit()`, a 64 bit limb long division**, and with it a benchmark row that finally covers
  the integer types' 128 bit divisor path — `Division by 128-bit value above 2^64`. Every row that
  existed before divided by a value that fits in one QWORD, so nothing measured the long division
  the integer types actually run.

  Every 128 bit type reaches long division once the divisor needs more than 64 bits: `int128_t` and
  `uint128_t` divide 128 by 128, `fixed_point128` and `float128` divide 256 by 128 because they
  scale the numerator by 2^128 first. `div_32bit()` answers all of those in 32 bit limbs, which
  makes them m=4,n=4 and m=8,n=4 — up to five main loop passes with a four iteration
  multiply-subtract in each. In 64 bit limbs the same divisions are one and three passes of two.

  Measured on the 256/128 shape, pinned to one core, across three divisor magnitudes:

  | | MSVC 19.51 | clang-cl 22.1.3 |
  |---|---|---|
  | `div_32bit`, 32 bit limbs | 217-227 cycles | 353-375 cycles |
  | `div_128bit`, 64 bit limbs | 68-73 cycles | 44-46 cycles |
  | | **3.1x** | **7.9x** |

  Clang used to pay a 1.6x penalty against MSVC on `div_32bit`, which four measurements had failed
  to explain — it is not the DIV instruction (16.0 against 16.3 cycles), not the `alloca`, not
  auto-vectorization, and not branch prediction. It was the 32 bit limb code itself: on
  `div_128bit` clang is the faster of the two toolchains.

  End to end over the whole benchmark suite, against the previous commit, best of three interleaved
  runs pinned to one core — **+8.2% geometric mean on MSVC and +10.9% on clang-cl over all 92 rows**,
  with no row outside the noise moving backwards:

  | Benchmark | MSVC | clang-cl |
  |---|---|---|
  | `float128` division by a 128 bit value | +137% | +307% |
  | `uint128_t` / `int128_t` division by a value above 2^64 | +97% / +99% | +195% / +175% |
  | `fixed_point128` division by a 128 bit fractional value | +124% | +110% |
  | `float128` / `fixed_point128` `sqrt` | +56% / +46% | +95% / +33% |
  | `asinh`, `acosh`, `atanh` | +9% to +27% | +8% to +54% |

  `div_32bit()` is left in place but no longer has a caller inside the library.

### Changed

- **`FP128_USE_RECIPROCAL_FOR_DIVISION` now defaults to `0`**, so `fixed_point128` divides by long
  division rather than by multiplying by a reciprocal. The default was `1` because the reciprocal
  was 1.4x-1.8x faster, and it cost up to 1.7 ulp to get that. `div_128bit()` inverts the trade:
  the long division now runs at **2.0x** (MSVC) and **2.7x** (clang-cl) the rate of the reciprocal
  over a table of 256 random divisors, and it was already the more accurate of the two — the
  residual `|a - q * b|` puts it closer to the exact quotient every single time the two disagree.
  Set the macro to `1` to get the old behaviour back.

- **Five call sites that multiplied by a reciprocal now divide.** `fixed_point128::sqrt()` ran a
  Newton loop whose every iteration called `reciprocal()`, itself a Newton loop —
  `float128::sqrt()` had always divided instead. `atan()` and `pow()` in both types inverted their
  argument or result the same way. All five are faster and none are less accurate, since the
  division is the more accurate operation. The two `pow()` sites keep their old answer for a zero:
  `fixed_point128` has no infinity, so `reciprocal()` answers zero there where `operator/=` throws,
  and that case is now guarded rather than divided.

- **`uint128_t` and `int128_t` are no longer trivially copyable**, and are measurably faster for
  it. Their copy and move members were `= default`; they now assign the two QWORDs individually.
  MSVC implements a *trivial* 16-byte copy as a single `vmovups` once `/arch:AVX2` is on, and
  every operator in the class produces its result as two QWORD stores — so that one wide load
  overlaps both stores, fails store-to-load forwarding, and stalls for roughly 15 cycles. It fires
  on any expression whose 128-bit result reaches memory, `*out = *a * *b` included.

  Measured on MSVC with plain array code (`out[i] = in[i] op x`, 64K elements), interleaving the
  two builds so drift and throttling hit both equally: `sqrt` 1.7x, division by a 128-bit value
  1.4x, `std::vector` growth 1.35x, addition 1.05x, and no measurable change to multiplication,
  `log10`, `pow` or division by a 64-bit integer. Nothing measured slower. The `vector` result is
  the surprising one — losing `memmove` relocation costs less than the member-wise copy saves, at
  this size.

  What this breaks: `std::is_trivially_copyable_v` is now `false` for both types, so `std::bit_cast`
  no longer accepts them and a container of them relocates with a copy loop instead of a
  `memmove`. Both remain standard layout, and `tests/{u,}int128_t_gtest.cpp` now assert the
  absence of trivial copyability so the members cannot drift back to `= default`.
  `fixed_point128` and `float128` have always been written this way.

- **`fixed_point128` is now 128 bits wide, not 129.** The separate `uint32_t sign` member is
  gone; the object is the QWORD pair `high:low` read as a two's complement integer, with the
  sign in the MSB of `high`, and the value it stands for is that integer divided by
  2<sup>F</sup>:

  ```
  value = (int128)(high:low) / 2^F,  where F = 127 - I
  ```

  The template parameter **I** keeps its meaning - the integer bits, not counting the sign -
  and its range narrows from `[1, 64]` to `[1, 63]`. The fraction gives up the bit the sign
  now occupies, so `F` is one smaller for every instantiation and results carry one bit less
  precision. The value range becomes `[-2^I, 2^I - 2^-F]`, asymmetric like every two's
  complement type: the most negative value has no positive counterpart, so negating it - or
  taking `fabs` of it - wraps back to itself, exactly as `-INT64_MIN` does.

  Consequences worth knowing:
  - **A signed zero no longer exists.** It used to be representable and compared unequal to
    plain zero; every operation now produces the one zero there is.
  - **`operator>>` shifts arithmetically**, replicating the sign bit, and keeps rounding to
    nearest as before.
  - **The bitwise operators cover all 128 bits**, the sign included. `~x` is now `-x - epsilon()`.
  - **`operator uint64_t` and `operator uint32_t` wrap** a negative value modulo 2^64 or 2^32,
    the way a builtin signed to unsigned conversion does, rather than returning the magnitude.
    The signed conversions still truncate towards zero.
  - Multiplication and division work on magnitudes internally, so `(-a) * b` stays bit
    identical to `-(a * b)`, as it was before.

- **Breaking API change** — the raw component accessors drop their sign parameter:
  `get_components(low, high)` and the constructor `fixed_point128(low, high)` take the two
  QWORDs of the value, sign included. The three argument forms are gone rather than adapted,
  so existing calls fail to compile instead of silently meaning something else.

- **`std::numeric_limits`** — `digits` is 127 rather than 128, `min_exponent` is `I - 127`,
  and `lowest()` is `-2^I`, one step below `-max()`.

- **Performance** — measured against the previous representation with `bench -t fixed_point128`
  on MSVC: the comparisons run at 1.6x, addition and subtraction at 1.5x, multiplication by a
  32-bit integer at 1.5x, division by a power of two at 1.6x, and `log`/`log2`/`log10` at
  1.1x, none of which need to look at a sign field any more. Multiplication by a 128-bit value
  is unchanged and the Mandelbrot composite is at 0.95x; what pays for the rest is `exp` and
  its family at 0.82x-0.87x and the trigonometric functions at 0.85x-0.91x, which multiply
  values of alternating sign and so take the magnitude round trip on every term.

- The debugger visualizer in `fixed_point128.natvis` decodes the new layout.

### Fixed

- **The benchmark understated `uint128_t` and `int128_t` by up to 6x on MSVC.** `BenchBinary()`
  timed `result = lhs op rhs` into an escaped result. MSVC builds the operator's return value in a
  stack temporary with two QWORD stores and then copies it out with a single 16-byte load; that
  load overlaps both stores, cannot be store-forwarded, and waits for them to reach L1 — roughly
  15 cycles on every iteration. Only the two integer types were charged for it, because the copy
  MSVC is free to widen is the compiler-generated one, and `fixed_point128` and `float128` assign
  their two QWORDs individually. The result was a table in which `fixed_point128<10>` multiplied
  two 128-bit values faster than `uint128_t` did — 491M/s against 258M/s, for four multiplies plus
  a shift and a rounding step against three multiplies.

  `BenchBinary()` now applies the operation in place, as `BenchAccumulate()` already did for the
  same reason. The 128-bit multiplication reads 1.56G/s for `uint128_t` and 490M/s for
  `fixed_point128<10>`, and clang-cl — which never emitted the wide copy, and so was measuring the
  multiply all along — agrees with both. The division measurements gain between 1.15x and 1.7x,
  the stall being a smaller share of an operation that is dominated by a `div`. `BenchUnary()` was
  checked and is not affected, since every unary function the integer types have returns a
  `uint64_t`.

- **`udiv128()` was a `__udivti3` library call on Clang and GCC, not a divide instruction.** The
  portable body divides a `__uint128_t` by a `uint64_t`; neither compiler can prove the quotient
  fits in 64 bits, so both lower it to the compiler-rt/libgcc helper — a full software 128÷128
  division — where MSVC reaches the `DIV` instruction through its `_udiv128` intrinsic. This is
  what made MSVC look ~2.8x faster than clang-cl at 128-bit division; the gap was Clang's, not an
  MSVC measurement artifact.

  x86-64 under a GCC-style frontend now issues the instruction directly as inline assembly, behind
  the new `FP128_X64` detection macro and the `std::is_constant_evaluated()` guard the other
  assembly helpers already use. AArch64 keeps the portable path — it has no 128/64 divide
  instruction. Interleaved measurements with clang-cl over 64K-element arrays: division by a
  64-bit integer 2.7x (92 → 250 M/s), by a 128-bit value 2.5x (70 → 178 M/s), `sqrt` 1.9x
  (16.5 → 31 M/s). clang-cl now edges out MSVC on division by a 64-bit integer.

  The precondition is `hi_dividend < divisor` — `DIV` raises #DE rather than truncating when the
  quotient overflows — which is the contract `_udiv128` already carried, and which both call sites
  satisfy: `div_32bit()` passes a zero high half, `div_64bit()` passes the remainder of a previous
  division by the same divisor. `tests/uint128_t_gtest.cpp` now checks the division identity over
  edge cases and 65536 random ones, on whichever spelling the toolchain selected.

  Side effect: an x86-64 clang-cl build no longer references `__udivti3` at all, so it no longer
  needs `clang_rt.builtins` on the link line. `cmake/FP128ClangRuntime.cmake` stays, since the
  portable path is still what ARM64 uses.

- **`uint128_t::operator/=(T)` copied its numerator and quotient through stack arrays**, and paid
  for both. `div_64bit()` was handed a `{low, high}` copy as the numerator and wrote the quotient
  into the members, or — while this was being tracked down — a local quotient copied back over
  them. MSVC vectorizes either copy into a 16-byte access that overlaps the QWORD stores on the
  other side of it and cannot forward from them, putting a stall immediately in front of the two
  divisions and another behind them.

  Numerator and quotient are now the members themselves, which are adjacent QWORDs in a standard
  layout class. `div_64bit()` reads `u[j]` into a local before writing `q[j]` and walks `j`
  downwards, so the two may alias; every early return that would leave the quotient unwritten is
  already excluded by the trivial-case checks in the caller.

  This is what keeps the copy-member change above from costing anything here. With the numerator
  copy still in place that change made `out[i] = in[i] / 5` run 3.4x *slower*, because the stall
  simply moved from the assignment into the copy; removing both copies restores the operation to
  its previous rate.

## [0.10.0.0] - 2026-08-10

### Added

- **`float128` completeness** — the 25 remaining `<cmath>` functions:
  `rint`/`nearbyint`, `scalbn`/`scalbln`, `nextafter`/`nexttoward`,
  `remainder`/`remquo`, three-argument `hypot`, `lgamma`/`tgamma`, `cbrt`, `abs`,
  `nan(const char*)` and the classification predicates.
- **Exact decimal conversion** — new [`include/fp128_decimal.h`](include/fp128_decimal.h)
  implements exact binary128 ↔ decimal, so a value round trips through its text form
  and any requested precision is correctly rounded. This backs the string constructor,
  `to_chars`/`from_chars` and `std::formatter`.
- **Accuracy test harness** — `tests/float128_ref_check.h`, with generated reference
  vectors in `tests/float128_ref_data.h` produced by `tools/gen_ref_vectors.py`,
  measures every function against correctly rounded references and states its error
  as a ulp bound.
- New test suites `tests/float128_accuracy_gtest.cpp` and
  `tests/float128_format_gtest.cpp`.

### Changed

- Precision improved in several existing `float128` functions.
- **Documentation** — GitHub Pages now publishes the full Doxygen API reference rather
  than a Jekyll-rendered `README.md`, built and deployed by the new
  `.github/workflows/docs.yml`.
- Version bumped to `0.10.0.0` (`FP128_VERSION_MINOR` 9 → 10).

### Fixed

- The GitHub Pages build and deployment failure introduced by the GoogleTest submodule,
  whose own Jekyll site broke the root-level build.

## [0.9.0.0] - 2026-08-08

The first release to carry a version number; this entry summarizes the work that led
up to it.

### Added

- **Library versioning** — the `FP128_VERSION*` macros and `fp128::version*` constants,
  which the benchmark prints and records in its JSON report.
- The missing scalar-on-the-left arithmetic operators for all four types, constrained
  on `std::is_arithmetic_v`.
- `sqr()`, which squares faster than `x * x`.
- A CMake build (`CMakeLists.txt`, `CMakePresets.json`) alongside the Visual Studio
  solution, which was converted to the `.slnx` format.
- New `tools/` directory: `gen_log2_tables.py`, `log2_ulp_check.py` and
  `log2_ulp_dump.cpp` generate and ulp-check the three log2 reduction tables.
- `tests/diffcalls.py`, which diffs inlining between two builds.

### Changed

- **constexpr conversions** — `Double`/`Float` now hold a raw bit pattern converted with
  `std::bit_cast` behind `f()`/`e()`/`s()`/`val()` accessors plus a `make(s, e, f)`
  factory, replacing the bit-field overlay that constant evaluation forbids; constexpr
  coverage was then pushed across all four 128-bit classes.
- **Accuracy** — `log2` and multiplication reworked for `float128` and `fixed_point128`,
  backed by the three generated reduction tables.
- **Performance** — add/sub optimization, two force-inline passes, and a 6.7x fix on
  128-bit multiply: binary operators return the named parameter (`return lhs OP= rhs`)
  to avoid an MSVC store-forwarding stall. `FP128_USE_RECIPROCAL_FOR_DIVISION` became a
  build-overridable macro. Also earlier `fixed_point128` addition and multiplication
  improvements, and ARM64 assembly improvements.
- **Benchmark harness** — largely rewritten: `Escape()`/`Barrier()` instead of an opaque
  call (clang had been vectorizing the addition chain to a bogus 10G/s), batch sizing
  with fastest-batch reporting, dead-code fixes, coverage for all four types, and JSON
  reports.
- Renamed `fixed_point128_shared.h` → `fp128_shared.h` across headers, VS projects,
  tests, tools and docs, and reorganized the library into `include/`, `tests/`, `bench/`,
  `tools/` and `msvc/`.
- `int128_t` and `uint128_t` restructured onto a single shared base.
- Improved conversion to string, and an improved `.natvis` file.
- Hardened all GoogleTest suites, with GoogleTest itself pulled in as a submodule.

### Fixed

- `fixed_point128` took magnitudes by negating a signed value — undefined behavior at
  `INT_MIN`/`INT64_MIN` — in two constructors and in `operator*=`/`operator/=`; the
  negation is now done in the unsigned domain.
- `asin`/`acos` now reject a Newton step that fails to improve their residual.
- An incorrect oracle in the `int128_t`/`uint128_t` sqrt tests.
- Division and addition with special numbers, and multiplication with NaN.
- `fixed_point128` when the integer part is 64 bits, `atan2`, and the round-to-even
  implementation.
- clang compatibility for `fixed_point128` and for the benchmark's add/sub/mul/div
  conditions.

[0.10.0.0]: https://github.com/ericgur/fixed_point128/releases/tag/v0.10.0.0
[0.9.0.0]: https://github.com/ericgur/fixed_point128/releases/tag/v0.9.0.0
