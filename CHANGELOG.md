# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to a four-part `MAJOR.MINOR.PATCH.BUILD` version scheme
exposed through the `FP128_VERSION*` macros and the `fp128::version*` constants in
[`include/fp128_shared.h`](include/fp128_shared.h).

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
