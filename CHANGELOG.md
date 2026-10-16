# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to a four-part `MAJOR.MINOR.PATCH.BUILD` version scheme
exposed through the `FP128_VERSION*` macros and the `fp128::version*` constants in
[`include/fp128_shared.h`](include/fp128_shared.h).

## [1.0.1.0] - 2026-10-16

### Fixed

- **Qt projects could not include the library.** The decimal formatter in `fp128_decimal.h` named a
  local lambda `emit`, which Qt defines as an empty macro, so every header that includes the
  formatter failed to compile in a Qt project. The lambda is renamed. `tests/qt_keywords_gtest.cpp`
  now builds every header with each unprefixed Qt macro (`emit`, `signals`, `slots`, `foreach`,
  `forever`, `SIGNAL`, `qDebug`, `M_PI` and the rest) defined as a compile error, so a library
  identifier that matches one of them fails the build.

## [1.0.0.0] - 2026-10-15

### Fixed

`float128` now conforms to IEEE 754-2008. An audit against exact references found the following,
each now checked by `tests/float128_ieee_gtest.cpp` and `tools/ieee_check.py`:

- **Division truncated** whenever the dividend's mantissa was the smaller of the two: 29% of random
  quotients came out one ulp low. The quotient's normalization assumed a fixed position for its
  leading bit and so skipped the rounding.
- **Addition and subtraction** decided the rounding from a three bit window, without a sticky bit
  for the rest: about 5% of results were an ulp off.
- **Every result in the subnormal range was rounded twice** - multiplication, division, `fma`,
  `ldexp`, the shifts and decimal input - and about one in ten came out an ulp off. Every operation
  now ends in a single rounding step that knows the destination's width.
- **`sqrt`** is correctly rounded; a quarter of its results used to be an ulp off.
- **Conversions:**
  - to `double` and `float`: correctly rounded (`float` used to round twice, through `double`);
  - overflow and underflow keep their sign;
  - 2^-1074 no longer converts to zero;
  - NaNs keep their sign and payload;
  - `float128(-0.0)` is negative zero;
  - a quiet NaN from a `double` stays quiet, where it used to become a signaling one;
  - `long double` converts exactly where it is wider than `double`.
- **Rounding to integral:** `ceil`, `trunc`, `round`, `rint`, `nearbyint` and `modf` keep the sign
  of a zero result. `round(0.5 - 2^-114)` is 0 rather than 1, and `round(2^112 + 1)` is itself.
- **NaNs:** every operation quiets a signaling NaN and delivers the payload of its first NaN
  operand. `get_class()` reports `signalingNaN`, and `fmin`/`fmax` treat a signaling NaN as invalid.
- **Math functions:**
  - `log`, `log2` and `log10` of a negative argument are NaN rather than -inf, of +inf are +inf,
    and of a NaN are NaN (they were finite);
  - `log2` of a subnormal uses its real exponent;
  - `pow` follows every special case of IEEE 754-2008 9.2.1 and no longer overflows for an integer
    exponent above 11355 (`pow(1 + 2^-40, 12000)` was infinite);
  - `atan2(0, 1)` is 0 rather than pi, and every zero and infinity follows the standard;
  - `asin`/`acos` outside [-1, 1] are NaN rather than 0, `acos(1)` is 0 rather than NaN, and both
    are accurate near the ends of the domain (they were off by up to 2^54 ulp);
  - `sin`, `cos` and `tan` reduce arguments of any size (they returned garbage above 2^62, and
    never returned at all above 2^321);
  - `hypot` scales instead of overflowing, and an infinity wins over a quiet NaN;
  - `cosh` and `sinh` no longer overflow just below their own overflow point.
- **Hexadecimal text:** `0x1.8p1` is read correctly (it used to read as 1), with correct rounding;
  `from_chars` accepts `chars_format::hex`. Shortened hexadecimal output rounds ties to even.
- **`inf / 0`** no longer counts as a division by zero.

A review of the test suite, with coverage and mutation testing, found the following. Each is now
checked by the test file of its type.

`float128`:

- **`remainder` and `remquo` near the top of the range** doubled the remainder to decide the rounding,
  which overflows above half of `max`: the exact `remainder(1.5 * 2^16383, max)` raised overflow and
  inexact, and rounding towards zero or downward returned `x` itself with a quotient of 0.
- `isless()`, `isgreater()` and the rest of the quiet comparisons, and `isunordered()`, raise the
  invalid exception for a signaling NaN, as IEEE 754 5.11 requires.
- `ilogb()` of 0, an infinity or a NaN raises the invalid exception; `float128(float)` of a signaling
  NaN raises it too, as the `double` conversion did; `frexp()` of a signaling NaN quiets it.
- `is_subnormal()` is false for a zero and `is_int()` is true for one, as documented.
- `double_factorial()` of a negative argument read before its table; (-1)!! is now 1 and other
  negative arguments are NaN. The table is built at compile time, which also removes a data race on
  its first use.
- `exp10()` of a negative exponent was up to 1809 ulp off (3453 in the directed rounding modes). It
  is now correctly rounded except for `exp10(-4931)` rounded upward, which is 1 ulp off.
- **`cos(0)` was 1 - 2^-113** rather than 1, as was the cosine of every tiny argument; `sin(pi/2)`,
  `cos(pi)` and `atan(1)` now round to their constants as well.
- **`lgamma` and `tgamma`** lost every bit near the zeros of `lgamma` at 1 and 2 and next to the
  poles of `tgamma` (`lgamma(2 - 2^-100)` was 6e32 ulp off, `tgamma(-1 + 2^-40)` 4.7e11 ulp). Both are
  rewritten; `tgamma` is exact for every integer up to 38, the largest factorial a float128 holds
  exactly being 37!, where it used to stop at 35. They are also about three times faster.
- `pow(x, 0.5)` is `sqrt(x)`, so an exact root is exact; `pow(denorm_min, 0.5)` was 4078 ulp off.
  `log1p(+-denorm_min)` returned zero.
- **Decimal text input** dropped every digit past the 40th, so a long input could round the wrong
  way. It is now correctly rounded at any length.
- **Fixed notation of a value below the last requested digit** was wrong: `{:.0f}` of 0.7 was "0"
  and `{:.1f}` of 0.046 was "0.1", through `std::format`, the streams and `to_chars` alike.
- `to_chars` in `chars_format::hex` wrote "0x1.8p+01", which `from_chars` could not read back; it now
  writes "1.8p+1", as `std::to_chars` does. `{:.Na}` pads with zeros past 28 digits, and
  `float128("-0x1p100000")` is negative infinity rather than positive.
- `operator>>` reads the number and stops where it ends, as it does for a `double` ("1.5,2.5" reads
  1.5), accepts hexadecimal and the special values, and sets `failbit` on bad input.
- The documented accuracy of `pow` (4 + 2.5 |y log2 x| ulp), `erfc` (56 ulp), `lgamma` and `tgamma`
  is now what they achieve; the old bounds were exceeded outside the regions the tests sampled.

`fixed_point128`:

- **`pow(1, y)` was 0** for every `y`.
- **`x % y` and `fmod()`** were computed from the quotient rounded to nearest, so a quotient that
  rounded up to an integer gave 0: `(9 - epsilon) % 3` was 0 rather than `3 - epsilon`. They now
  return the exact remainder with the sign of `x`. `lowest() % -1` crashed the process with an
  integer overflow trap at `I = 63`.
- **`hypot()`** overflowed its intermediate sum of squares: `hypot(300, 400)` was 0 at `I = 16`.
- **`sqrt()`** lost about I/2 bits, up to 2.9e9 ulp at `I = 63`. It is now correctly rounded, as are
  `hypot()` and `reciprocal()`.
- **`sin()` and `cos()`**: the series wrapped around at `I = 16`, and `cos(0)` was not 1 at `I = 63`.
  `tan`, `atan`, `asin` and `acos`, which are built on them, improve with them; `asin` and `acos`
  were also off by up to 8e16 ulp near +-1.
- **Overflow inside a function**: `exp` of a negative argument and `tanh` returned garbage once
  `e^|x|` overflowed, and `asinh`, `acosh` and `atanh` overflowed and terminated the program.
  `atan2` lost its sign when `x / y` underflowed.
- **`log`, `log10` and `exp2`** used ln 2 and log10(2) constants with wrong digits, up to 7015 ulp
  off at `I = 1`.
- **`0 / 0` and `0 % 0`** returned 0 instead of throwing, as every other division by zero does.
- `fixed_point128<1>` to `<3>` converted to text with digits like `"0.:6:::"`.
- **The text form was not one-to-one**: it printed a digit too few, so `epsilon()` printed as "0",
  `max()` as 2^I, which is out of range, and three quarters of the values at `I = 32` read back as
  another value. Every value now survives a round trip through text, and the string constructor is
  correctly rounded at every width (it was up to 22.8 ulp off at `I = 1`).
- `operator>>` reports bad and out of range input through `failbit`, and
  `numeric_limits<fixed_point128>::is_modulo` is `true`, since the arithmetic wraps.
- `reciprocal()` returned garbage rather than the documented 0 when the result is out of range;
  `logb()` of a negative argument threw, and `ilogb(0)` returned `I - 128` instead of `FP_ILOGB0`.
- The accuracy of every math function is measured and stated in its documentation.

The integer types:

- **`uint128_t` and `int128_t` compared with a floating point value** truncated it first, so
  `int128_t(2) == 2.5` held and `int128_t(-2) > -2.5` did not, and saturated it, so the largest
  value equalled NaN, infinity and `1e39`. The two are now compared as real numbers: a fraction
  orders correctly, a NaN is unordered (only `!=` holds), and a value beyond the range orders beyond
  every value of the type. Integral operands are compared as before.
- The binary `&`, `|` and `^` are `noexcept`, like their compound forms.
- `shift_right128()` and `shift_left128()` return the right value for a shift of 0. No caller in the
  library passed one.
- `operator>>` behaves like it does for a builtin integer: input without digits sets `failbit` and
  stores 0, a value out of range sets `failbit` and stores the end of the range (it used to wrap
  silently), reading stops where the number ends, and `std::hex` and `std::oct` are honoured.

### Added

- **`FP128_IEEE_ENV`**: opt in IEEE 754 rounding directions and exception flags, through
  `fp128::fegetround`, `fesetround`, `fetestexcept`, `feclearexcept`, `feraiseexcept`,
  `fegetexceptflag` and `fesetexceptflag`. Without the macro nothing changes in speed.
- `totalorder`, `totalordermag`, `fmaxmag`, `fminmag`, `iscanonical`, `issignaling`,
  `issubnormal` and `iszero`, as named in ISO/IEC TS 18661-1.
- Text input reads `snan` and `nan(payload)`.
- `tools/ieee_check.py` with `tools/ieee_probe.cpp`, a differential check of the correctly rounded
  operations in every rounding direction, and `tools/gen_two_over_pi.py` for the trigonometric
  reduction table.
- A test program of its own for the environment, `fixed_point128_env_tests`.

### Changed

- **Conversion from `float128` to the integer types truncates towards zero** at every magnitude, as
  the builtin conversions do. It used to round to nearest from one up and truncate below one, so
  `(int)float128(1.75)` was 2 and `(int)float128(0.75)` was 0. A value out of range saturates and a
  NaN converts to zero (a NaN used to saturate). `llrint`/`lrint` return 0 when out of range, as
  their documentation always said.
- **`std::numeric_limits<float128>::is_iec559` and `float128::is754version2008()`** are `true` only
  with `FP128_IEEE_ENV`, since the standard requires the environment the macro adds.
- The `iterations` argument of `sqrt()` no longer has an effect, for `float128` and `fixed_point128`
  alike: the root is correctly rounded either way.
- The relational operators raise the invalid exception for a NaN operand, as IEEE 754's signaling
  comparisons do; `==`, `!=` and `isless()` and friends raise it only for a signaling NaN.
- **float128 is faster.** Rounding once, branch free, and the reformulated functions measure ahead of
  0.12.0.0 on 31 of the 32 float128 benchmarks with MSVC and 30 with Clang (P-core; median
  +31% and +18%): addition +22% / +28%, multiplication +29% / +3%, division +18% / +14%,
  `sqrt` +15% / +21%, `exp` +41% / +23%, `sin` +56% / +36%, Mandelbrot +7% / +4% (MSVC / Clang).
  The exception is Clang on an E-core, where addition and subtraction measure 23-24% slower.
- **Multiplying a `uint128_t` or `int128_t` by the other 128 bit integer type, or by `float128`, no
  longer compiles.** It used to compile and silently use only the low 64 bits of the right operand:
  `int128_t(3) * uint128_t(2^64 + 1)` was 3. The other mixed operators were already rejected.
- **`fixed_point128` rounds to nearest, ties away from zero, in every operation.** Multiplication
  and division by a `fixed_point128` already did. Right shifts, narrowing conversions between
  widths and division by a power of two now do as well; they decided the rounding from a single
  dropped bit, an error of up to 0.75 ulp, and broke ties to even. Division by an integer used to
  truncate. One consequence: `epsilon() >> 1` is `epsilon()`, so a loop that halves a value until it
  reaches zero does not end, as was already the case for `x *= 0.5`.
- `float128::reciprocal()` is `1 / x`: correctly rounded, and about 3.2 times faster than the Newton
  iteration it replaces. `reciprocal(0)` is infinity, as its documentation now says.
- **Long division is specialized for the width of its numerator.** `div_128bit()` takes the word
  count as a template argument, so each type gets a division unrolled for its own shape whether or
  not the compiler inlines it; with the count an ordinary argument, clang++ compiled one generic copy
  for every type. Measured against the same code without it (MSVC and Clang, P- and E-core): division of
  the 128 bit integers by a divisor above 2^64 +22% to +63%, `float128` division +7% to +53% and
  `reciprocal` +10% to +69%, `fixed_point128` division by a fraction +11% to +44% (-6% with MSVC on
  a P-core). `div_128bit(q, r, u, v, m)` still takes the count at run time.
- **The reworked math functions are faster as well**, measured against the code before the test
  review (MSVC and Clang, P- and E-core): `float128` `sin`, `cos`, `tan`, `asin`, `acos` and `atan`
  +17% to +45%; `fixed_point128` `tanh` +63% to +84%, `sin` +34% to +44%, `tan` +19% to +30% and
  `acos` +12% to +23%. `fixed_point128::atanh` is 24% to 36% slower, the price of no longer
  overflowing and of its accuracy near +-1.

## [0.12.0.0] - 2026-10-03

### Fixed

- **`float128::log2()` near one** — the error just below one (up to 8.7 ulp) and just above it
  (up to 3.8 ulp) is now at most 1.5 ulp on both sides, and arguments within 2^-4 of one run 3-6x
  faster. `log`, `log10` and `log1p` gain 4-13% from the same change.
- **`tools/log2_ulp_dump`** builds again against the two's complement `fixed_point128`, emits
  negative results correctly, labels its input classes, and gains a `below-one` class.
- **`std::numeric_limits<fixed_point128<I>>`** conforms to the standard and to the arithmetic:
  - every member function is `constexpr`, as the standard requires of a specialization;
  - `min_exponent` is `I - 126`, so that 2^(min_exponent-1) is the smallest positive value, the
    way `max_exponent = I` already read; it was one too low;
  - `max_digits10` is 40 rather than 39 for the fifteen `I` where 39 digits print some
    neighbouring values identically (4, 7, 10, 14, 17, 20, 24, 27, 30, 37, 40, 47, 50, 57, 60);
  - `round_style` is `round_to_nearest`, which is what multiplication and division do;
  - `traps` is `true`: division by zero throws, the same as for the integer types;
  - `denorm_min()` returns the smallest positive value instead of zero.
- **`std::numeric_limits<int128_t>::is_modulo`** is `true`. The signed type wraps on overflow
  exactly as the unsigned one does.
- **No deprecation warnings under C++23 and C++26.** Clang reported two per specialization for the
  `has_denorm` members, which are now declared with the warning suppressed locally.

### Removed

- The `numeric_limits` specializations for `const`, `volatile` and `const volatile` qualified
  types. `<limits>` already provides them, forwarding to the unqualified specialization.

### Changed

- The MSVC warning suppressions shared by every first-party target live in one place:
  `fp128_set_warnings()` for CMake and `msvc/fixed_point128.warnings.props` for MSBuild, which
  also silences code-analysis warning C6262 for the benchmark's 16 KB argument arrays.

## [0.11.0.0] - 2026-08-22

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

[1.0.1.0]: https://github.com/ericgur/fixed_point128/releases/tag/v1.0.1.0
[1.0.0.0]: https://github.com/ericgur/fixed_point128/releases/tag/v1.0.0.0
[0.12.0.0]: https://github.com/ericgur/fixed_point128/releases/tag/v0.12.0.0
[0.11.0.0]: https://github.com/ericgur/fixed_point128/releases/tag/v0.11.0.0
[0.10.0.0]: https://github.com/ericgur/fixed_point128/releases/tag/v0.10.0.0
[0.9.0.0]: https://github.com/ericgur/fixed_point128/releases/tag/v0.9.0.0
