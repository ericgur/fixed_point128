# fixed_point128 library
 A 128 bit fixed-point class template for fast, high precision calculations.
 This code is used in my [Mandelbrot just-for-fun project](https://github.com/ericgur/Mandelbrot). With double precision floats I could zoom the image to 2<sup>44</sup>, with the fixed_point128, 2<sup>113</sup> is possible.
 
 The `inc` directory contains the header files for the **fixed_point128** library, a header-only C++20 library providing 128-bit integer, fixed-point, and floating-point arithmetic types. 
 
 All types reside in the `fp128` namespace.

## fixed_point128 Template class Highlights
 - Most operations are very fast. 1-10x slower than double precision. ~10x faster than MPIR at similar precision.
 - Up to 38 fraction digits (decimal) are supported.
 - Has a superset of integer and floating point functions including all standard C/C++ operators.
 - The single template parameter **\<I\>** allows the user to specify 1-64 bits for the integer part, the rest are allocated to the fraction.
 - An object can be created from all int/float types as well as from strings representing a float.
 - Supports conversions from one template instance to another (2 instances with different **\<I\>** parameter).

## float128 class Highlights
 - Based on the IEEE 754 binary128 format.
 - Implements most of the standard library.
 
 ## Dependencies and Prerequisites
 - C++20
 - Standard C++ library.
 - 64 bit builds only
 - A supported compiler:
   - **MSVC** - Visual Studio 2019+ (Windows).
   - **Clang 17+** - Windows (Clang toolset), Linux or macOS.
   - **GCC** - Linux or macOS. GCC takes the same code path as Clang (`__uint128_t` and `__builtin_*`).

## Building

The library itself is header-only: add the `inc` directory to your include path and include the header you need. There is nothing to compile or link. The sections below cover the benchmark and test suite that ship with the repository.

### Benchmark

**GCC / Clang** - via CMake (3.20 or newer). Builds `src/Bench.cpp` into `bin/`:

```sh
cmake -S . -B build
cmake --build build
./bin/bench
```

**MSVC** - open `fixed_point128.slnx` in Visual Studio 2019+ and build.

### Tests

The GoogleTest suite lives in `gtest`. CMake must be installed and on the `PATH`. Run the scripts below from the `gtest` directory.

**MSVC** - run `setup.bat` once to generate `build\fixed_point128_gtest.sln`, then after each code change:
- `build.bat` - build the test app.
- `test.bat` - run all tests.
- `test_failed.bat` - re-run only the tests that failed.

**Clang** - uses the Clang toolset that ships with Visual Studio (install the "C++ Clang tools for Windows" component). It builds into `build_clang\`, so both toolchains can be kept side by side. Run `setup_clang.bat` once, then `build_clang.bat` and `test_clang.bat`.

 ## Dependency Graph

```
fixed_point128_shared.h
    |
    +--- int128_shared.h
    |        |
    |        +--- int128_t.h
    |        |
    |        +--- uint128_t.h
    |                 |
    |                 +--- float128.h
    |
    +--- fixed_point128.h
```

All headers depend on `fixed_point128_shared.h`. The `float128` class additionally depends on `uint128_t.h`. Do not include `fixed_point128_shared.h` or `int128_shared.h` directly; they are pulled in automatically by the other headers.

---

## File Details

### fixed_point128.h

Template class `fixed_point128<I>` where **I** is the number of integer bits (range `[1, 64]`). The remaining `128 - I` bits store the fractional part. This gives compile-time control over the trade-off between range and precision.

**Data layout:** `uint64_t low` + `uint64_t high` + `uint32_t sign` (separate sign bit).

**Features:**
- Construction from integer types, `double`, C strings (accurate to 37 decimal digits), and raw components.
- Full arithmetic operators including optimized 64-bit multiply path.
- Cross-template assignment and conversion between different `I` values.
- Comprehensive math library:
  - **Basic:** `fabs`, `floor`, `ceil`, `trunc`, `round`, `copysign`, `fmod`, `modf`, `fdim`, `fmin`, `fmax`.
  - **Power / Root:** `sqr`, `sqrt`, `pow`, `hypot`, `reciprocal`.
  - **Exponential / Logarithmic:** `exp`, `exp2`, `expm1`, `log`, `log2`, `log10`, `log1p`, `logb`.
  - **Trigonometric:** `sin`, `cos`, `tan`, `asin`, `acos`, `atan`, `atan2`.
  - **Hyperbolic:** `sinh`, `cosh`, `tanh`, `asinh`, `acosh`, `atanh`.
- Built-in constants: `pi()`, `pi2()`, `half_pi()`, `e()`, `sqrt_2()`, `golden_ratio()`, `one()`, `half()`, `epsilon()`.

### float128.h

IEEE 754-2008 binary128 (quadruple-precision) floating-point type, aligned to 16 bytes.

**Bit layout:** 1-bit sign, 15-bit exponent (bias 16383), 112-bit fraction.

**Data layout:** `uint64_t low` + `uint64_t high`.

**Features:**
- Full IEEE 754 special-value handling: NaN propagation, infinity arithmetic, subnormals.
- Classification queries: `is_zero`, `is_finite`, `is_normal`, `is_subnormal`, `is_nan`, `is_signaling`, `is_inf`, `is_special`, `is_int`, `is_negative`, `is_positive`, `is_exponent_of_2`.
- Construction from `float`, `double`, integer types, and C strings (including scientific notation and special values).
- Arithmetic operators: `+`, `-`, `*`, `/`, `<<`, `>>`.
- Comprehensive math library (50+ functions):
  - **Basic:** `fabs`, `floor`, `ceil`, `trunc`, `round`, `copysign`, `fmod`, `modf`, `fdim`, `fmin`, `fmax`.
  - **Power / Root:** `sqr`, `sqrt`, `cbrt`, `pow`, `hypot`.
  - **Exponential / Logarithmic:** `exp`, `exp2`, `expm1`, `log`, `log2`, `log10`, `log1p`, `logb`.
  - **Trigonometric:** `sin`, `cos`, `tan`, `asin`, `acos`, `atan`, `atan2`.
  - **Hyperbolic:** `sinh`, `cosh`, `tanh`, `asinh`, `acosh`, `atanh`.
  - **Error functions:** `erf`, `erfc`.
  - **Rounding:** `llrint`, `llround`, `lrint`, `lround`.
  - **Other:** `frexp`, `ldexp`, `ilogb`, `reciprocal`, `double_factorial`.
- Built-in constants: `pi()`, `half_pi()`, `e()`, `sqrt_2()`, `tenth()`.
- User-defined literal: `_f128` (e.g. `3.14_f128`).

### int128_t.h

Signed 128-bit integer stored in two's complement representation, aligned to 16 bytes.

**Data layout:** `uint64_t low` (lower 64 bits) + `uint64_t high` (upper 64 bits; MSB is the sign bit).

**Features:**
- Construction from integer types, `double`, C strings, `std::string`, and raw `low`/`high` components.
- Full arithmetic operators: `+`, `-`, `*`, `/`, `%`, `<<`, `>>`, `&`, `|`, `^`, and their assignment variants.
- Comparison operators: `==`, `!=`, `<`, `<=`, `>`, `>=`.
- Math functions: `abs`, `sqrt`, `log`, `log2`, `log10`, `pow`.
- Conversions to `int64_t`, `uint64_t`, `float`, `double`, `long double`, and strings.
- User-defined literal: `_int128` (e.g. `12345_int128`).

### uint128_t.h

Unsigned 128-bit integer, aligned to 16 bytes. Mirrors the API surface of `int128_t` with unsigned semantics.

**Data layout:** `uint64_t low` + `uint64_t high`.

**Features:**
- Same constructor set and operator suite as `int128_t`, adapted for unsigned arithmetic.
- Math functions: `sqrt`, `log`, `log2`, `log10`, `pow`.
- User-defined literal: `_uint128` (e.g. `99999_uint128`).

### int128_shared.h

Implementation shared by both 128-bit integer types. `int128_t` and `uint128_t` are aliases of a single class template defined here, `int128_base<bool IsSigned>`, so the constructors, operators and math functions listed above physically live in this header; `int128_t.h` and `uint128_t.h` only declare the aliases and their user-defined literals.

**Key contents:**
- **`int128_base<IsSigned>`** - 16-byte aligned class template holding `uint64_t low` + `uint64_t high`. `int128_t` is `int128_base<true>`, `uint128_t` is `int128_base<false>`. Values are stored in two's complement, so the bit patterns of the shared operations match exactly between the signed and unsigned types.
- **Constructors** - default, copy, move, `double`, any builtin integral type, and C strings (decimal or hexadecimal with an optional sign).
- **Operator suite** - arithmetic, bitwise, shift and comparison operators, each with a template overload accepting any type convertible to `int128_base`.
- **Free functions** - `abs`, `sqr`, `sqrt`, `log`, `log2`, `log10`, `pow`, `lzcnt128`.
- **Conversions** - `operator double`, `operator float`, `operator std::string`.

Because `int128_t` and `uint128_t` are aliases rather than distinct classes, neither can be forward declared; include `int128_t.h` or `uint128_t.h` instead.

### fixed_point128_shared.h

Foundation header providing platform-specific intrinsic wrappers and common helper functions used by the rest of the library.

**Key contents:**
- **Build configuration macros** - Compiler detection (`FP128_MSVC`, `FP128_CLANG`), inline control (`FP128_INLINE`, `FP128_FORCE_INLINE`), and feature flags (`FP128_CPP_STYLE_MODULO`, `FP128_USE_RECIPROCAL_FOR_DIVISION`).
- **Intrinsic wrappers** - Portable wrappers for `lzcnt`, `popcnt`, `mulx`, `addcarryx`, and `udiv128` covering both MSVC and GCC/Clang.
- **128-bit shift functions** - `shift_right128`, `shift_left128`, and rounding variants.
- **Multi-word division** - `div_32bit` and `div_64bit`, derived from *Hacker's Delight* by Henry S. Warren Jr.
- **Bit manipulation** - `lzcnt128`, `popcnt128`, `log2`, and `twos_complement128`.
- **IEEE 754 unions** - `Double` and `Float` structs for accessing bit fields of native floating-point values.

## Platform Support

The library targets 64-bit platforms and provides optimized intrinsic paths for:
- **MSVC** (Windows) - Uses compiler intrinsics such as `_umul128`, `_udiv128`, `__lzcnt64`, `_addcarry_u64`.
- **GCC / Clang** (Linux, macOS, Windows with Clang toolset) - Uses `__uint128_t` and `__builtin_*` functions.

 
 ## Visual Studio Debugger Support
 Add the file **fixed_point128.natvis** to your Visual Studio project.
 Enables pretty print of the objects value (limited to double precision).

 ## Code Examples

 ### Mandelbrot main loop
 Modified code snippet from my Mandelbrot C++/QT6 project: https://github.com/ericgur/Mandelbrot
 Allows plotting the Mandelbrot (or Julia set) with a zoom of 2<sup>113</sup>.

```cpp
constexpr int MAX_ITERATION = 1000;

fixed_point128<8> xc, yc;  // coordinate arguments.
fixed_point128<8> radius = 2, radius_sq = sqr(radius);
fixed_point128<8> usq = 0, vsq = 0, u = 0, v = 0, tmp, modulus_sq = 0;
int iter = 0;

// Find how many iterations are needed to have the (xc,yc) coordinates diverge (absolute value > 2)
while (modulus_sq < radius_sq && ++iter < MAX_ITERATION) {
    // real
    tmp = usq - vsq + xc;

    // imaginary
    // v = 2.0 * (u * v) + yc;
    v = ((u * v) << 1) + yc;
    u = tmp;
    usq = sqr(u);
    vsq = sqr(v);
    // compare the squared magnitude against radius^2, avoiding a square root
    modulus_sq = usq + vsq;
}
```

## Acknowledgements
- `div_32bit` (multi-precision integer division) is derived from the book *"Hacker's Delight"* 2nd Edition by Henry S. Warren Jr. 
It was converted to 32 bit operations and modified a bit. The algorithm is an implementation of Knuth's "Algorithm D" from the book *"The Art of Computer Programming"*.
- Logarithm functions are derived from [Dan Moulding's log2fix](https://github.com/dmoulding/log2fix).
- Square root uses Newton-Raphson iteration based on *Math Toolkit for Real Time Programming* by Jack W. Crenshaw.

## License

Released under the MIT License, Copyright (c) 2022 Eric Gur. See [LICENSE](LICENSE) for the full text.
