# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to a four-part `MAJOR.MINOR.PATCH.BUILD` version scheme
exposed through the `FP128_VERSION*` macros and the `fp128::version*` constants in
[`include/fp128_shared.h`](include/fp128_shared.h).

## [0.11.0.0] - unreleased

### Changed

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
