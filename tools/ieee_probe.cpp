/**
 * @file ieee_probe.cpp
 * @brief The C++ half of the IEEE 754 differential check, see ieee_check.py and README.md.
 *
 * Reads one operation per line and writes its result as bits, so that the Python half can compare
 * every bit - the sign of a zero and the quiet bit of a NaN included - against a reference computed
 * on exact rationals. Nothing passes through decimal on the way, which would round away exactly
 * the bits being checked.
 *
 * Input lines are "op a_low a_high b_low b_high c_low c_high" in hex, or "parse text". Each output
 * line is "low high flags", where flags is the set of float128 exceptions the operation raised
 * (1 invalid, 2 division by zero, 4 overflow, 8 underflow, 16 inexact), always zero unless the
 * probe was built with FP128_IEEE_ENV. "setmode m" selects the rounding direction: 0 to nearest,
 * 1 towards zero, 2 upward, 3 downward.
 *
 * Usage: ieee_probe <input file> <output file>
 */

#include <cstdio>
#include <cstring>
#include <string>
#include "float128.h"

using namespace fp128;

namespace
{
/// @brief The exceptions raised since the last call, as the bits the Python half expects, cleared.
int TakeFlags()
{
    const int raised = fp128::fetestexcept(FE_ALL_EXCEPT);
    int bits = 0;
    if (raised & FE_INVALID)
        bits |= 1;
    if (raised & FE_DIVBYZERO)
        bits |= 2;
    if (raised & FE_OVERFLOW)
        bits |= 4;
    if (raised & FE_UNDERFLOW)
        bits |= 8;
    if (raised & FE_INEXACT)
        bits |= 16;
    fp128::feclearexcept(FE_ALL_EXCEPT);
    return bits;
}

/// @brief Writes a float128 result.
void Put(FILE* out, const float128& r)
{
    uint64_t low = 0, high = 0;
    r.get_bits(low, high);
    fprintf(out, "%016llx %016llx %x\n", static_cast<unsigned long long>(low), static_cast<unsigned long long>(high), TakeFlags());
}

/// @brief Writes a result of up to 64 bits: a double, a float or an integer.
void PutBits(FILE* out, uint64_t bits)
{
    fprintf(out, "%016llx 0 %x\n", static_cast<unsigned long long>(bits), TakeFlags());
}
}  // namespace

int main(int argc, char** argv)
{
    if (argc != 3) {
        fprintf(stderr, "usage: ieee_probe <input file> <output file>\n");
        return 2;
    }
    FILE* in = fopen(argv[1], "r");
    FILE* out = fopen(argv[2], "w");
    if (in == nullptr || out == nullptr) {
        fprintf(stderr, "cannot open the files\n");
        return 2;
    }

    char op[64];
    unsigned long long a0 = 0, a1 = 0, b0 = 0, b1 = 0, c0 = 0, c1 = 0;
    while (fscanf(in, "%63s", op) == 1) {
        if (strcmp(op, "parse") == 0) {
            static char text[1024];
            if (fscanf(in, "%1023s", text) != 1)
                break;
            Put(out, float128(static_cast<const char*>(text)));
            continue;
        }
        if (fscanf(in, "%llx %llx %llx %llx %llx %llx", &a0, &a1, &b0, &b1, &c0, &c1) != 6)
            break;

        const float128 a(a0, a1), b(b0, b1), c(c0, c1);
        const std::string o(op);
        if (o == "setmode") {
            static const int modes[4] = {FE_TONEAREST, FE_TOWARDZERO, FE_UPWARD, FE_DOWNWARD};
            fp128::fesetround(modes[a0 & 3]);
            TakeFlags();
            fprintf(out, "0 0 0\n");
        } else if (o == "add") {
            Put(out, a + b);
        } else if (o == "sub") {
            Put(out, a - b);
        } else if (o == "mul") {
            Put(out, a * b);
        } else if (o == "div") {
            Put(out, a / b);
        } else if (o == "sqrt") {
            Put(out, sqrt(a));
        } else if (o == "fma") {
            Put(out, fma(a, b, c));
        } else if (o == "sqr") {
            Put(out, sqr(a));
        } else if (o == "rint") {
            Put(out, rint(a));
        } else if (o == "nearbyint") {
            Put(out, nearbyint(a));
        } else if (o == "floor") {
            Put(out, floor(a));
        } else if (o == "ceil") {
            Put(out, ceil(a));
        } else if (o == "trunc") {
            Put(out, trunc(a));
        } else if (o == "round") {
            Put(out, round(a));
        } else if (o == "remainder") {
            Put(out, remainder(a, b));
        } else if (o == "fmod") {
            Put(out, fmod(a, b));
        } else if (o == "ldexp") {
            Put(out, ldexp(a, static_cast<int>(static_cast<int64_t>(b0))));
        } else if (o == "fromdouble") {
            double d = 0;
            memcpy(&d, &a0, sizeof(d));
            Put(out, float128(d));
        } else if (o == "todouble") {
            const double d = static_cast<double>(a);
            uint64_t bits = 0;
            memcpy(&bits, &d, sizeof(d));
            PutBits(out, bits);
        } else if (o == "tofloat") {
            const float f = static_cast<float>(a);
            uint32_t bits = 0;
            memcpy(&bits, &f, sizeof(f));
            PutBits(out, bits);
        } else if (o == "toint64") {
            PutBits(out, static_cast<uint64_t>(static_cast<int64_t>(a)));
        } else if (o == "toint32") {
            PutBits(out, static_cast<uint32_t>(static_cast<int32_t>(a)));
        } else if (o == "touint64") {
            PutBits(out, static_cast<uint64_t>(a));
        } else if (o == "less") {
            PutBits(out, a < b);
        } else if (o == "equal") {
            PutBits(out, a == b);
        } else if (o == "isless") {
            PutBits(out, isless(a, b));
        } else if (o == "log") {
            Put(out, log(a));
        } else if (o == "log2") {
            Put(out, log2(a));
        } else if (o == "log10") {
            Put(out, log10(a));
        } else if (o == "log1p") {
            Put(out, log1p(a));
        } else if (o == "exp") {
            Put(out, exp(a));
        } else if (o == "exp2") {
            Put(out, exp2(a));
        } else if (o == "expm1") {
            Put(out, expm1(a));
        } else if (o == "sin") {
            Put(out, sin(a));
        } else if (o == "cos") {
            Put(out, cos(a));
        } else if (o == "tan") {
            Put(out, tan(a));
        } else if (o == "asin") {
            Put(out, asin(a));
        } else if (o == "acos") {
            Put(out, acos(a));
        } else if (o == "atan") {
            Put(out, atan(a));
        } else if (o == "sinh") {
            Put(out, sinh(a));
        } else if (o == "cosh") {
            Put(out, cosh(a));
        } else if (o == "tanh") {
            Put(out, tanh(a));
        } else if (o == "atanh") {
            Put(out, atanh(a));
        } else if (o == "logb") {
            Put(out, logb(a));
        } else if (o == "nextup") {
            Put(out, nextafter(a, float128::inf()));
        } else if (o == "atan2") {
            Put(out, atan2(a, b));
        } else if (o == "pow") {
            Put(out, pow(a, b));
        } else if (o == "hypot") {
            Put(out, hypot(a, b));
        } else if (o == "fmin") {
            Put(out, fmin(a, b));
        } else if (o == "fmax") {
            Put(out, fmax(a, b));
        } else {
            fprintf(out, "unknown\n");
        }
    }
    fclose(out);
    fclose(in);
    return 0;
}
