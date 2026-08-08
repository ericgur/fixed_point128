/***********************************************************************************
    MIT License

    Copyright (c) 2025 Eric Gur (ericgur@iname.com)

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.
************************************************************************************/

/**
 * @file fp128_decimal.h
 * @brief Exact conversion between binary128 and decimal.
 *
 * A binary128 value is a 113 bit integer times a power of two, so its decimal expansion is finite:
 * at most 4933 digits before the point and 16494 after it. Everything here works with that
 * expansion in full rather than with an approximation of it, which is what lets the digits be
 * correctly rounded at any requested precision and lets a decimal string be read back to the
 * nearest representable value.
 *
 * The alternative - scaling by a power of ten held in the type itself - cannot do either. A
 * negative power of ten is not representable in binary, so every scaling rounds, and the error
 * compounds across the digits. That is what the previous conversions did, and it left the string
 * form of a value about an ulp away from the value and unable to survive a round trip.
 *
 * These are implementation details of float128, not part of the interface.
 *
 * @see float128.h
 */

#ifndef FP128_DECIMAL_H
#define FP128_DECIMAL_H

#include <cstdint>
#include "fp128_shared.h"

namespace fp128
{
namespace detail
{

/**
 * @brief Fixed capacity unsigned big integer, 32 bit limbs, least significant first.
 *
 * Sized for the widest intermediate any conversion below needs: reading the smallest subnormal
 * back from its decimal form scales a 40 digit mantissa by 2^16716, which is 527 limbs.
 *
 * The operations are only the ones the conversions use - multiply and divide by a value that fits
 * in a limb, shift, compare - so none of them needs a general multi word division.
 */
class big_uint
{
public:
    static constexpr int32_t CAPACITY = 560;         ///< Limbs, 17920 bits.
    static constexpr int32_t LIMB_BITS = 32;         ///< Bits per limb.
    static constexpr uint64_t LIMB_BASE = 1ull << LIMB_BITS;  ///< One past the largest limb value.

    /// @brief Constructs zero.
    big_uint() noexcept : used(0), limb {} {}

    /// @brief Sets the value to a 128 bit integer.
    void set(uint64_t low, uint64_t high) noexcept
    {
        clear();
        limb[0] = static_cast<uint32_t>(low);
        limb[1] = static_cast<uint32_t>(low >> 32);
        limb[2] = static_cast<uint32_t>(high);
        limb[3] = static_cast<uint32_t>(high >> 32);
        used = 4;
        trim();
    }

    /// @brief Sets the value to zero.
    void clear() noexcept
    {
        for (int32_t i = 0; i < used; ++i)
            limb[i] = 0;
        used = 0;
    }

    [[nodiscard]] bool is_zero() const noexcept { return used == 0; }
    [[nodiscard]] int32_t limbs() const noexcept { return used; }
    [[nodiscard]] uint32_t at(int32_t index) const noexcept { return (index < used) ? limb[index] : 0u; }

    /// @brief Index of the highest set bit, or -1 when the value is zero.
    [[nodiscard]] int32_t bit_length() const noexcept
    {
        if (used == 0)
            return 0;
        const uint32_t top = limb[used - 1];
        int32_t bit = 0;
        for (uint32_t probe = top; probe != 0; probe >>= 1)
            ++bit;
        return (used - 1) * LIMB_BITS + bit;
    }

    /// @brief Multiplies by a value that fits in a limb.
    /// @param multiplier Value to multiply by
    /// @return False when the product did not fit in the capacity.
    bool mul_small(uint32_t multiplier) noexcept
    {
        uint64_t carry = 0;
        for (int32_t i = 0; i < used; ++i) {
            const uint64_t product = static_cast<uint64_t>(limb[i]) * multiplier + carry;
            limb[i] = static_cast<uint32_t>(product);
            carry = product >> LIMB_BITS;
        }
        while (carry != 0) {
            if (used >= CAPACITY)
                return false;
            limb[used++] = static_cast<uint32_t>(carry);
            carry >>= LIMB_BITS;
        }
        return true;
    }

    /// @brief Adds a value that fits in a limb.
    void add_small(uint32_t addend) noexcept
    {
        uint64_t carry = addend;
        for (int32_t i = 0; carry != 0; ++i) {
            if (i >= used) {
                if (i >= CAPACITY)
                    return;
                limb[i] = 0;
                used = i + 1;
            }
            const uint64_t sum = static_cast<uint64_t>(limb[i]) + carry;
            limb[i] = static_cast<uint32_t>(sum);
            carry = sum >> LIMB_BITS;
        }
    }

    /// @brief Divides by a value that fits in a limb.
    /// @param divisor Value to divide by, must not be zero
    /// @return The remainder.
    uint32_t divmod_small(uint32_t divisor) noexcept
    {
        uint64_t remainder = 0;
        for (int32_t i = used - 1; i >= 0; --i) {
            const uint64_t dividend = (remainder << LIMB_BITS) | limb[i];
            limb[i] = static_cast<uint32_t>(dividend / divisor);
            remainder = dividend % divisor;
        }
        trim();
        return static_cast<uint32_t>(remainder);
    }

    /// @brief Multiplies in place by 10^count, in blocks that fit in a limb.
    /// @return False when the product did not fit in the capacity.
    bool mul_pow10(int32_t count) noexcept
    {
        // 10^9 is the largest power of ten below 2^32
        while (count >= 9) {
            if (!mul_small(1000000000u))
                return false;
            count -= 9;
        }
        return (count <= 0) || mul_small(POWERS_OF_10[count]);
    }

    /// @brief Divides in place by 10^count, reporting whether anything was discarded.
    /// @param count Power of ten to divide by
    /// @param sticky Set to true when a division left a non zero remainder
    void div_pow10(int32_t count, bool& sticky) noexcept
    {
        while (count >= 9) {
            sticky = (divmod_small(1000000000u) != 0) || sticky;
            count -= 9;
        }
        if (count > 0)
            sticky = (divmod_small(POWERS_OF_10[count]) != 0) || sticky;
    }

    /// @brief Shifts left, growing the value.
    /// @return False when the result did not fit in the capacity.
    bool shift_left(int32_t bits) noexcept
    {
        if (bits <= 0 || used == 0)
            return true;

        const int32_t limb_shift = bits / LIMB_BITS;
        const int32_t bit_shift = bits % LIMB_BITS;
        const int32_t extra = (bit_shift != 0) ? 1 : 0;
        if (used + limb_shift + extra > CAPACITY)
            return false;

        // Copied from the top down so a source limb is read before it is overwritten.
        for (int32_t i = used - 1 + extra; i >= 0; --i) {
            uint32_t value = 0;
            if (bit_shift == 0) {
                value = (i < used) ? limb[i] : 0u;
            } else {
                const uint32_t high = (i < used) ? limb[i] : 0u;
                const uint32_t low = (i > 0) ? limb[i - 1] : 0u;
                value = static_cast<uint32_t>((high << bit_shift) | (low >> (LIMB_BITS - bit_shift)));
            }
            limb[i + limb_shift] = value;
        }
        for (int32_t i = 0; i < limb_shift; ++i)
            limb[i] = 0;

        used += limb_shift + extra;
        trim();
        return true;
    }

    /// @brief Shifts right, discarding the bits that fall off the bottom.
    void shift_right(int32_t bits) noexcept
    {
        if (bits <= 0 || used == 0)
            return;
        const int32_t limb_shift = bits / LIMB_BITS;
        const int32_t bit_shift = bits % LIMB_BITS;
        if (limb_shift >= used) {
            clear();
            return;
        }

        for (int32_t i = 0; i + limb_shift < used; ++i) {
            const uint32_t low = limb[i + limb_shift];
            if (bit_shift == 0) {
                limb[i] = low;
            } else {
                const uint32_t high = (i + limb_shift + 1 < used) ? limb[i + limb_shift + 1] : 0u;
                limb[i] = static_cast<uint32_t>((low >> bit_shift) | (high << (LIMB_BITS - bit_shift)));
            }
        }
        const int32_t remaining = used - limb_shift;
        for (int32_t i = remaining; i < used; ++i)
            limb[i] = 0;
        used = remaining;
        trim();
    }

    /// @brief Clears every bit at or above the given position.
    void keep_low_bits(int32_t bits) noexcept
    {
        if (bits <= 0) {
            clear();
            return;
        }
        const int32_t limb_index = bits / LIMB_BITS;
        const int32_t bit_index = bits % LIMB_BITS;
        if (limb_index >= used)
            return;

        for (int32_t i = limb_index + ((bit_index != 0) ? 1 : 0); i < used; ++i)
            limb[i] = 0;
        if (bit_index != 0)
            limb[limb_index] &= (1u << bit_index) - 1u;
        trim();
    }

    /**
     * @brief Treats the value as a fraction, multiplies it by a small value and peels off the
     *        integer part that crosses the point.
     *
     * The buffer holds a fraction scaled by 2^(width*32), so multiplying it produces digits above
     * the point that no longer belong in it. Those are what the caller wants.
     *
     * @param multiplier Value to multiply by, small enough that the overflow fits in a limb
     * @param width Limb count the fraction occupies
     * @return The integer part the multiply produced.
     */
    uint32_t fraction_mul_extract(uint32_t multiplier, int32_t width) noexcept
    {
        uint64_t carry = 0;
        for (int32_t i = 0; i < width; ++i) {
            const uint64_t product = static_cast<uint64_t>(limb[i]) * multiplier + carry;
            limb[i] = static_cast<uint32_t>(product);
            carry = product >> LIMB_BITS;
        }
        used = width;
        trim();
        return static_cast<uint32_t>(carry);
    }

    /// @brief Compares the value against one half of 2^(width*32), the fraction's midpoint.
    /// @return 1 above the midpoint, -1 below it, 0 exactly on it.
    [[nodiscard]] int32_t compare_half(int32_t width) const noexcept
    {
        const uint32_t top = at(width - 1);
        constexpr uint32_t half_limb = 1u << (LIMB_BITS - 1);
        if (top > half_limb)
            return 1;
        if (top < half_limb)
            return -1;
        for (int32_t i = 0; i < width - 1; ++i) {
            if (at(i) != 0)
                return 1;
        }
        return 0;
    }

    /// @brief Powers of ten that fit in a limb, indexed by the exponent.
    static constexpr uint32_t POWERS_OF_10[10] = {1u,      10u,      100u,      1000u,      10000u,
                                                  100000u, 1000000u, 10000000u, 100000000u, 1000000000u};

private:
    /// @brief Drops the leading zero limbs so that used is the significant limb count.
    void trim() noexcept
    {
        while (used > 0 && limb[used - 1] == 0)
            --used;
    }

    int32_t used;              ///< Significant limb count.
    uint32_t limb[CAPACITY];   ///< Limbs, least significant first.
};

/// @brief Largest number of significant digits any conversion here produces.
///
/// A binary128 needs 36 to survive a round trip, and one more is kept so that a caller asking for
/// the maximum can still see the digit that decides the rounding.
inline constexpr int32_t MAX_SIGNIFICANT_DIGITS = 40;

/// @brief Largest number of digits to_decimal_digits() produces.
///
/// Enough to write any finite binary128 out in full in fixed notation, the widest of which has
/// 4933 digits before the point. The exact expansion of a subnormal runs to 16494 places after the
/// point instead, and a request past this cap is filled with zeros rather than continuing it.
inline constexpr int32_t MAX_OUTPUT_DIGITS = 5000;

/// @brief Limbs of the fraction buffer. 16640 bits reaches below the smallest subnormal's last bit.
inline constexpr int32_t FRACTION_LIMBS = 520;

/**
 * @brief Correctly rounded decimal digits of a finite, non zero binary128 value.
 *
 * The value is decomposed into an exact integer part and an exact fractional part, and the digits
 * of each are produced by arithmetic that cannot round: repeated division by a power of ten for
 * the integer part, repeated multiplication for the fraction. Only the last digit kept is rounded,
 * once, from the exact remainder.
 *
 * @param mantissa_low Low QWORD of the 113 bit mantissa, its leading one at bit 112
 * @param mantissa_high High QWORD of the mantissa
 * @param exponent Unbiased exponent, so the value is mantissa * 2^(exponent-112)
 * @param requested Significant digits wanted, clamped to MAX_SIGNIFICANT_DIGITS
 * @param digits Output buffer, receives `requested` digits as characters with no terminator
 * @param exact Optional, set to true when the digits are the whole value and nothing was discarded
 * @return The decimal exponent: the value is 0.<digits> * 10^result.
 */
inline int32_t to_decimal_digits(uint64_t mantissa_low, uint64_t mantissa_high, int32_t exponent, int32_t requested, char* digits,
                                 bool* exact = nullptr) noexcept
{
    if (requested > MAX_OUTPUT_DIGITS)
        requested = MAX_OUTPUT_DIGITS;
    if (requested < 1)
        requested = 1;

    // The value is mantissa * 2^binary_exponent.
    const int32_t binary_exponent = exponent - 112;

    big_uint integer_part;
    big_uint fraction;
    integer_part.set(mantissa_low, mantissa_high);

    if (binary_exponent >= 0) {
        integer_part.shift_left(binary_exponent);
    } else {
        const int32_t below_point = -binary_exponent;
        // The fraction buffer holds the fractional part scaled by 2^(FRACTION_LIMBS*32). That
        // scale reaches 2^-16640, below the 2^-16494 of the smallest subnormal's last bit, so the
        // shift is always to the left and nothing is lost.
        fraction.set(mantissa_low, mantissa_high);
        fraction.keep_low_bits(below_point);
        fraction.shift_left(FRACTION_LIMBS * big_uint::LIMB_BITS - below_point);
        integer_part.shift_right(below_point);
    }

    // Integer digits come out least significant first, nine at a time.
    char reversed[5000];
    int32_t integer_digits = 0;
    while (!integer_part.is_zero() && integer_digits + 9 <= static_cast<int32_t>(sizeof(reversed))) {
        uint32_t block = integer_part.divmod_small(1000000000u);
        for (int32_t i = 0; i < 9; ++i) {
            reversed[integer_digits++] = static_cast<char>('0' + (block % 10));
            block /= 10;
        }
    }
    // The last block is zero padded to nine digits, and those zeros are not part of the value.
    while (integer_digits > 0 && reversed[integer_digits - 1] == '0')
        --integer_digits;

    // The value is 0.<digits> * 10^decimal_exponent. A value below one has no integer digits, and
    // the zeros the fraction starts with are not significant: each one moves the exponent instead.
    int32_t decimal_exponent = integer_digits;
    int32_t emitted = 0;
    char rounding_digit = '0';
    bool have_rounding_digit = false;
    bool rest_nonzero = false;
    bool seen_significant = integer_digits > 0;

    const auto emit = [&](char digit) {
        if (emitted < requested)
            digits[emitted++] = digit;
        else if (!have_rounding_digit) {
            rounding_digit = digit;
            have_rounding_digit = true;
        } else {
            rest_nonzero = rest_nonzero || (digit != '0');
        }
    };

    for (int32_t i = integer_digits - 1; i >= 0; --i)
        emit(reversed[i]);

    // Fraction digits, nine at a time, until the request and its rounding digit are both covered.
    while (!fraction.is_zero() && (emitted < requested || !have_rounding_digit)) {
        uint32_t block = fraction.fraction_mul_extract(1000000000u, FRACTION_LIMBS);
        char block_digits[9];
        for (int32_t i = 8; i >= 0; --i) {
            block_digits[i] = static_cast<char>('0' + (block % 10));
            block /= 10;
        }
        for (int32_t i = 0; i < 9; ++i) {
            if (!seen_significant) {
                if (block_digits[i] == '0') {
                    --decimal_exponent;
                    continue;
                }
                seen_significant = true;
            }
            emit(block_digits[i]);
        }
    }
    // Anything still left in the buffer sits below every digit produced.
    rest_nonzero = rest_nonzero || !fraction.is_zero();

    // A value whose expansion ended early is padded, which is exact.
    while (emitted < requested)
        digits[emitted++] = '0';

    if (exact != nullptr)
        *exact = !rest_nonzero && (!have_rounding_digit || rounding_digit == '0');

    // Round half to even against the exact remainder.
    bool round_up = false;
    if (rounding_digit > '5')
        round_up = true;
    else if (rounding_digit == '5')
        round_up = rest_nonzero || ((digits[requested - 1] - '0') % 2 != 0);

    if (round_up) {
        int32_t index = requested - 1;
        while (index >= 0 && digits[index] == '9')
            digits[index--] = '0';
        if (index >= 0) {
            ++digits[index];
        } else {
            // Every digit was a nine, so the carry created a new leading one.
            digits[0] = '1';
            ++decimal_exponent;
        }
    }

    return decimal_exponent;
}

/**
 * @brief Reads a decimal mantissa and exponent back to the nearest binary128.
 *
 * The value is digits * 10^exponent10. Both directions are computed exactly - a positive power of
 * ten by multiplication, a negative one by division with the discarded remainder kept as a sticky
 * bit - so the 113 bit result is the correctly rounded one and never two roundings deep.
 *
 * @param digits Decimal digits, most significant first, no sign and no point
 * @param count Number of digits
 * @param exponent10 Power of ten the digits are scaled by
 * @param mantissa_low Receives the low QWORD of the 113 bit mantissa
 * @param mantissa_high Receives the high QWORD
 * @param exponent Receives the unbiased exponent, so the value is mantissa * 2^(exponent-112)
 * @return False when the value is zero, or too large or too small for the format to hold, in which
 *         case the exponent says which: below -100000 for an underflow, above 100000 for an
 *         overflow.
 */
inline bool from_decimal_digits(const char* digits, int32_t count, int32_t exponent10, uint64_t& mantissa_low, uint64_t& mantissa_high,
                                int32_t& exponent) noexcept
{
    mantissa_low = mantissa_high = 0;
    exponent = 0;

    // Leading zeros carry no information, trailing ones move into the exponent.
    int32_t first = 0;
    while (first < count && digits[first] == '0')
        ++first;
    if (first == count)
        return false;  // the value is zero
    while (count > first && digits[count - 1] == '0') {
        --count;
        ++exponent10;
    }

    // More digits than this cannot change which binary128 is nearest: the format separates
    // neighbours by one part in 2^112, and 40 digits resolve one part in 10^40.
    int32_t kept = count - first;
    bool sticky = false;
    if (kept > MAX_SIGNIFICANT_DIGITS) {
        for (int32_t i = first + MAX_SIGNIFICANT_DIGITS; i < count; ++i)
            sticky = sticky || (digits[i] != '0');
        exponent10 += kept - MAX_SIGNIFICANT_DIGITS;
        kept = MAX_SIGNIFICANT_DIGITS;
    }

    big_uint value;
    for (int32_t i = 0; i < kept; ++i) {
        value.mul_small(10u);
        value.add_small(static_cast<uint32_t>(digits[first + i] - '0'));
    }

    // A rough decimal magnitude is enough to reject what cannot be represented before building an
    // intermediate for it. log10 of the largest finite value is 4932, of the smallest subnormal
    // -4966, and one extra digit of slack on each side keeps every borderline case on the slow
    // path where it is decided exactly.
    const int32_t magnitude = exponent10 + kept;
    if (magnitude > 4935) {
        exponent = 200000;
        return false;
    }
    if (magnitude < -4970) {
        exponent = -200000;
        return false;
    }

    int32_t binary_exponent = 0;
    if (exponent10 >= 0) {
        if (!value.mul_pow10(exponent10)) {
            exponent = 200000;
            return false;
        }
    } else {
        // Scaling up before the division leaves room for the 113 bits wanted plus a guard bit.
        // A decimal place is log2(10) = 3.3219 bits and 10/3 is a hair more, so the shift always
        // covers the division with the fixed 160 bits of slack left over. Rounding that factor up
        // to four instead overshoots by 3000 limbs at the bottom of the exponent range, which is
        // past the capacity: the conversion then failed rather than being merely wasteful.
        const int32_t headroom = (-exponent10) * 10 / 3 + 160;
        if (!value.shift_left(headroom)) {
            exponent = -200000;
            return false;
        }
        binary_exponent -= headroom;
        value.div_pow10(-exponent10, sticky);
    }

    if (value.is_zero()) {
        exponent = -200000;
        return false;
    }

    // Take the top 113 bits, with the bit below them and everything under that.
    const int32_t bits = value.bit_length();
    const int32_t drop = bits - 113;
    if (drop > 0) {
        // The guard bit and the sticky bits have to be read before the shift discards them.
        const int32_t guard_index = drop - 1;
        const uint32_t guard = (value.at(guard_index / big_uint::LIMB_BITS) >> (guard_index % big_uint::LIMB_BITS)) & 1u;
        for (int32_t i = 0; i < guard_index && !sticky; ++i)
            sticky = ((value.at(i / big_uint::LIMB_BITS) >> (i % big_uint::LIMB_BITS)) & 1u) != 0;

        for (int32_t i = 0; i < drop; ++i)
            value.divmod_small(2);
        binary_exponent += drop;

        // round half to even
        if (guard != 0 && (sticky || (value.at(0) & 1u) != 0)) {
            value.add_small(1u);
            // A carry out of 113 bits lands on the next power of two.
            if (value.bit_length() > 113) {
                value.divmod_small(2);
                ++binary_exponent;
            }
        }
    } else if (drop < 0) {
        value.shift_left(-drop);
        binary_exponent += drop;
    }

    mantissa_low = static_cast<uint64_t>(value.at(0)) | (static_cast<uint64_t>(value.at(1)) << 32);
    mantissa_high = static_cast<uint64_t>(value.at(2)) | (static_cast<uint64_t>(value.at(3)) << 32);
    // The caller's convention puts the leading one at bit 112, so the exponent counts from there.
    exponent = binary_exponent + 112;
    return true;
}

}  // namespace detail
}  // namespace fp128

#endif  // FP128_DECIMAL_H
