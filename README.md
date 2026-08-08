# fixed_point128 library
 A 128 bit fixed-point class template for fast, high precision calculations.
 This code is used in my [Mandelbrot just-for-fun project](https://github.com/ericgur/Mandelbrot). With double precision floats I could zoom the image to 2<sup>44</sup>, with the fixed_point128, 2<sup>113</sup> is possible.
 
 The `include` directory contains the header files for the **fixed_point128** library, a header-only C++20 library providing 128-bit integer, fixed-point, and floating-point arithmetic types. 
 
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

 Using the library needs nothing beyond the above. Building the benchmark or the tests additionally
 needs CMake 3.21+ (or Visual Studio), and the tests need the GoogleTest submodule described below.

## Repository Layout

```
include/     the library - header only, this is all a consumer needs
bench/       benchmark program
tests/       GoogleTest suite
tools/       maintenance tools for log2: table generator and accuracy harness
external/    third party dependencies (GoogleTest, as a git submodule)
cmake/       CMake helper modules
msvc/        Visual Studio project files (the solution lives at the repository root)
DoxyGen/     documentation configuration
```

## Getting the Sources

GoogleTest is tracked as a git submodule, so clone recursively:

```sh
git clone --recurse-submodules https://github.com/ericgur/fixed_point128.git
```

If the repository was cloned without `--recurse-submodules`, fetch the submodule afterwards:

```sh
git submodule update --init --recursive
```

Only the tests need the submodule. Consumers of the library, and the benchmark, do not.

## Building

The library itself is header-only: add the `include` directory to your include path and include the
header you need. There is nothing to compile or link. The sections below cover the benchmark and the
test suite that ship with the repository.

Two build systems are maintained in parallel and both build the same sources:

- **CMake** - Windows, macOS and Linux, with MSVC, clang-cl, Clang or GCC. This is the command line path.
- **Visual Studio solution** (`fixed_point128.slnx`) - Windows only, with both the MSVC and Clang toolsets.

### CMake

Requires CMake 3.21 or newer (3.25 for the presets below), and Ninja for the `clang` and `gcc` presets.
Configure once per toolchain, then build and test:

```sh
cmake --preset msvc                   # configure; or: clang-cl, clang, gcc
cmake --build --preset msvc-release   # build; or: msvc-debug
ctest --preset msvc-release -j        # test;  or: msvc-debug
```

Every configure preset uses a multi-config generator, so the configure step is configuration agnostic:
there is no `msvc-debug` *configure* preset and `cmake --preset msvc-debug` is an error. Debug versus
Release is selected at build and test time, which is why `-debug` and `-release` appear only on
`cmake --build --preset` and `ctest --preset`. Building the directory directly rather than through a
build preset needs an explicit `--config`, otherwise CMake defaults to Debug:

```sh
cmake --build out/build/msvc --config Release
```

| Preset | Toolchain | Platform |
| --- | --- | --- |
| `msvc` | Visual Studio generator, MSVC toolset | Windows |
| `clang-cl` | Visual Studio generator, ClangCL toolset | Windows |
| `clang` | Ninja + `clang++` | macOS, Linux, Windows |
| `gcc` | Ninja + `g++` | Linux |

The `msvc` and `clang-cl` presets do not pin a generator, so CMake selects the newest Visual Studio
installed on the machine. The `clang-cl` preset needs the "C++ Clang tools for Windows" Visual Studio
component.

Everything lands under `out/build/<preset>/`. Useful extras:

```sh
ctest --preset msvc-debug --rerun-failed          # only what failed last time
ctest --preset msvc-debug -R float128             # only matching tests
cmake --preset msvc -DFP128_BUILD_TESTS=OFF       # skip the tests (no submodule needed)
./out/build/clang/bin/Release/bench               # run the benchmark
```

### Visual Studio

Open `fixed_point128.slnx` in Visual Studio 2019 or newer. The solution holds three projects:

| Project | Output |
| --- | --- |
| `bench` | the benchmark executable |
| `fixed_point128_tests` | the GoogleTest suite |
| `googletest` | GoogleTest built as a static library from the submodule |

Four configurations are available - `Debug`, `Release`, and `Debug_clang` / `Release_clang` for the
Clang toolset - all for `x64`. Binaries are written to `bin/`.

#### Running and debugging the tests

Set **fixed_point128_tests** as the startup project and press F5. The project reference to
`googletest` builds and links the static library automatically; there is nothing to configure.

- Breakpoints work in the test sources and throughout the header-only library, which is compiled into
  the test binary with full debug information.
- `fixed_point128.natvis` is compiled into the test PDB, so `fixed_point128<I>`, `float128`,
  `int128_t` and `uint128_t` show their values in Locals and Watch instead of raw QWORDs.
- **Test Explorer** lists the individual cases through the "Test Adapter for Google Test" component of
  the C++ workload. Right-click a single test and choose *Debug* to break inside just that case.
- To narrow down an F5 run, add arguments under Project Properties -> Debugging -> Command Arguments,
  for example `--gtest_filter=float128.Add* --gtest_break_on_failure`. The latter drops into the
  debugger at the first failing assertion. These are stored in the untracked `.vcxproj.user` file.

Switching the configuration to *Debug_clang* gives the same experience against the Clang-built binaries.

## Updating GoogleTest

The submodule is pinned to a release tag rather than tracking a branch, so that a given commit of this
repository always builds against a known GoogleTest version. To move the pin:

```sh
cd external/googletest
git fetch --tags
git checkout v1.18.0     # the desired release tag
cd ../..
git add external/googletest
git commit -m "Update GoogleTest to v1.18.0"
```

Everyone else picks the new version up with `git submodule update --init --recursive` after pulling.

`git submodule update --remote external/googletest` would instead advance to the tip of GoogleTest's
default branch. That is deliberately not used here, because it pins an arbitrary untagged commit.

 ## Dependency Graph

```
fp128_shared.h
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

All headers depend on `fp128_shared.h`. The `float128` class additionally depends on `uint128_t.h`. Do not include `fp128_shared.h` or `int128_shared.h` directly; they are pulled in automatically by the other headers.

---

## File Details

### fixed_point128.h

Template class `fixed_point128<I>` where **I** is the number of integer bits (range `[1, 64]`). The remaining `128 - I` bits store the fractional part. This gives compile-time control over the trade-off between range and precision.

**Data layout:** `uint64_t low` + `uint64_t high` + `uint32_t sign` (separate sign bit).

**Features:**
- Construction from integer types, `double`, C strings (accurate to 37 decimal digits), and raw components.
- Full arithmetic operators including optimized 64-bit multiply path.
- A builtin scalar may appear on either side of any binary operator. The scalar is widened rather than the object narrowed, so the result is `fixed_point128<I>` and `1 - f` keeps the fraction of `f`.
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
- Arithmetic operators: `+`, `-`, `*`, `/`, `<<`, `>>`. A builtin scalar may appear on either side, and the result is a `float128` either way.
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
- A builtin scalar may appear on either side of any of the above. The scalar is widened, so the result is 128 bit and the comparisons stay exact for values no builtin type can hold. This also makes `1 << n` produce the expected power of two where a builtin shift that wide would be undefined behavior.
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
- **Operator suite** - arithmetic, bitwise, shift and comparison operators, each with a template overload accepting any type convertible to `int128_base`, plus a mirror overload constrained to `std::is_arithmetic_v` for a builtin scalar on the left. The constraint is what keeps the two sets from being an equally good match for each other.
- **Free functions** - `abs`, `sqr`, `sqrt`, `log`, `log2`, `log10`, `pow`, `lzcnt128`.
- **Conversions** - `operator double`, `operator float`, `operator std::string`.

Because `int128_t` and `uint128_t` are aliases rather than distinct classes, neither can be forward declared; include `int128_t.h` or `uint128_t.h` instead.

### fp128_shared.h

Foundation header providing platform-specific intrinsic wrappers and common helper functions used by the rest of the library.

**Key contents:**
- **Build configuration macros** - Compiler detection (`FP128_MSVC`, `FP128_CLANG`), inline control (`FP128_INLINE`, `FP128_FORCE_INLINE`), and feature flags (`FP128_CPP_STYLE_MODULO`, `FP128_USE_RECIPROCAL_FOR_DIVISION`).
- **Intrinsic wrappers** - Portable wrappers for `lzcnt`, `popcnt`, `mulx`, `addcarryx`, and `udiv128` covering both MSVC and GCC/Clang.
- **128-bit shift functions** - `shift_right128`, `shift_left128`, and rounding variants.
- **Multi-word division** - `div_32bit` and `div_64bit`, derived from *Hacker's Delight* by Henry S. Warren Jr.
- **Bit manipulation** - `lzcnt128`, `popcnt128`, `log2`, and `twos_complement128`.
- **IEEE 754 unions** - `Double` and `Float` structs for accessing bit fields of native floating-point values.

**Overridable macros:** two of the build configuration macros are meant to be set from the build, either on the command line or by defining them before the first include of a library header.

| Macro | Default | Effect |
| --- | --- | --- |
| `FP128_DISABLE_INLINE` | `0` | Set it to `1` to turn every `FP128_INLINE` and `FP128_FORCE_INLINE` into `noinline`, so that a profile attributes time to the function it was actually spent in. |
| `FP128_USE_RECIPROCAL_FOR_DIVISION` | `1` | Selects how `fixed_point128` divides by a value that is neither a power of two nor an integer: `a * reciprocal(b)` when non-zero, long division when zero. The reciprocal is 1.4x-1.8x faster and up to 1.7 ulp less accurate. Only `fixed_point128` reads it - for `float128` the reciprocal measures both slower and less accurate, and the integer types have no reciprocal to multiply by. The comment on the macro carries the measurements. |

```sh
cl  /DFP128_USE_RECIPROCAL_FOR_DIVISION=0 ...   # MSVC, clang-cl
c++ -DFP128_USE_RECIPROCAL_FOR_DIVISION=0 ...   # Clang, GCC
```

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
