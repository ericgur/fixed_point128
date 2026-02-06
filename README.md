# fixed_point128 library
 A 128 bit fixed-point class template for fast, high precision calculations.
 This code is used in my [Mandelbrot just-for-fun project](https://github.com/ericgur/Mandelbrot). With double precision floats I could zoom the image to 2^44^, with the fixed_point128, 2^113^ is possible.
 
 This `inc` directory contains the header files for the **fixed_point128** library, a header-only C++20 library providing 128-bit integer, fixed-point, and floating-point arithmetic types. 
 
 All types reside in the `fp128` namespace.

## fixed_point128 Template class Highlights
 - Most operations are very fast. 1-10x slower than double precision. ~10x faster than MPIR at similar precision.
 - Up to 38 fraction digits (decimal) are supported.
 - Has a superset of integer and floating point functions including all standard C/C++ operators.
 - The single template paramter **\<I\>** allows the user to specify 1-64 bits for the integer part, the rest are allocated to the fraction.
 - An object can be created from all int/float types as well as from strings representing a float.
 - Supports converions from one template instance to another (2 instances with different **\<I\>** parameter).

## float128 class Highlights
 - Based on the IEEE 754 binary128 format.
 - Implementes most of the standard library.
 
 ## Dependencies and Perquisites
 - Visual Studio 2019+ (MSFT Compiler)
 - C++20
 - Standard C++ library.
 - 64 bit builds only
 
 ## Dependency Graph

```
fixed_point128_shared.h
    |
    +--- int128_t.h
    |
    +--- uint128_t.h
    |        |
    |        +--- float128.h
    |
    +--- fixed_point128.h
```

All headers depend on `fixed_point128_shared.h`. The `float128` class additionally depends on `uint128_t.h`. Do not include `fixed_point128_shared.h` directly; it is pulled in automatically by the other headers.

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
  - **Power / Root:** `sqrt`, `pow`, `hypot`, `reciprocal`.
  - **Exponential / Logarithmic:** `exp`, `exp2`, `expm1`, `log`, `log2`, `log10`, `log1p`, `logb`.
  - **Trigonometric:** `sin`, `cos`, `tan`, `asin`, `acos`, `atan`, `atan2`.
  - **Hyperbolic:** `sinh`, `cosh`, `tanh`, `asinh`, `acosh`, `atanh`.
- Built-in constants: `pi`, `e`, `sqrt2`, `golden_ratio`, `one()`, `half()`, `epsilon()`.

### float128.h

IEEE 754-2008 binary128 (quadruple-precision) floating-point type, aligned to 16 bytes.

**Bit layout:** 1-bit sign, 15-bit exponent (bias 16383), 112-bit fraction.

**Data layout:** `uint64_t low` + `uint64_t high`.

**Features:**
- Full IEEE 754 special-value handling: NaN propagation, infinity arithmetic, subnormals.
- Classification queries: `is_zero`, `is_finite`, `is_normal`, `is_subnormal`, `is_nan`, `is_signaling_nan`, `is_infinite`, `fpclassify`.
- Construction from `float`, `double`, integer types, and C strings (including scientific notation and special values).
- Arithmetic operators: `+`, `-`, `*`, `/`, `<<`, `>>`.
- Comprehensive math library (50+ functions):
  - **Basic:** `fabs`, `floor`, `ceil`, `trunc`, `round`, `copysign`, `fmod`, `modf`, `fdim`, `fmin`, `fmax`.
  - **Power / Root:** `sqrt`, `cbrt`, `pow`, `hypot`.
  - **Exponential / Logarithmic:** `exp`, `exp2`, `expm1`, `log`, `log2`, `log10`, `log1p`, `logb`.
  - **Trigonometric:** `sin`, `cos`, `tan`, `asin`, `acos`, `atan`, `atan2`.
  - **Hyperbolic:** `sinh`, `cosh`, `tanh`, `asinh`, `acosh`, `atanh`.
  - **Error functions:** `erf`, `erfc`.
  - **Rounding:** `llrint`, `llround`, `lrint`, `lround`.
  - **Other:** `frexp`, `ldexp`, `ilogb`, `nextafter`, `reciprocal`, `factorial`.
- Built-in constants: `pi`, `e`, `sqrt2`.
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
 Allows plotting the Mandelbrot (or Julia set) with a zoom of 2^113^.

    fixed_point128<8> x, y; // coordinate arguments.
    fixed_point128<8> radius = 2, radius_sq = radius * radius;
    fixed_point128<8> usq = 0, vsq = 0, u = 0, v = 0, tmp, uv, modulus = 0;
    
    // Find how many iterations are needed to have the (x,y) coordinates diverge (absolute value > 2)
    while (modulus < radius_sq && ++iter < MAX_ITERATION) {
        // real
        tmp = usq - vsq + xc;

        // imaginary
        //v = 2.0 * (u * v) + y;
        v = ((u * v) << 1) + yc;
        u = tmp;
        usq = u * u;
        vsq = v * v;
        // check uv vector amplitude is smaller than 2
        modulus = usq + vsq;
    }

## Acknologements
- `div_32bit` (multi-precision integer division) is derived from the book *"Hacker's Delight"* 2nd Edition by Henry S. Warren Jr. 
It was converted to 32 bit operations and mdified a bit. The algorithm is an implementation of Knuth's "Algorithm D" from the book *"The Art of Computer Pogramming"*.
- Logarithm functions are derived from [Dan Moulding's log2fix](https://github.com/dmoulding/log2fix).
- Square root uses Newton-Raphson iteration based on *Math Toolkit for Real Time Programming* by Jack W. Crenshaw.
**Acknowledgements:**  
