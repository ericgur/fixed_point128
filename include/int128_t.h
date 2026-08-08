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

/***********************************************************************************
                                Acknologements

************************************************************************************/

/**
 * @file int128_t.h
 * @brief 128-bit signed integer type.
 *
 * Provides @ref fp128::int128_t, a software-emulated 128-bit signed integer using
 * two's complement representation. Supports full arithmetic, bitwise, comparison,
 * and conversion operators, as well as math functions (abs, sqrt, log, log2,
 * log10, pow).
 * All methods are inline for maximum performance.
 *
 * @see int128_shared.h for the class template both 128 bit integer types alias.
 * @see fp128_shared.h for supporting intrinsics and utilities.
 */

#ifndef FP128_INT128_T_H
#define FP128_INT128_T_H

#include "int128_shared.h"

namespace fp128
{

/***********************************************************************************
 *                                  Main Code
 ************************************************************************************/

/**
 * @typedef int128_t
 * @brief 128 bit signed integer type.
 *
 * An alias of the signed instantiation of @ref fp128::int128_base, which holds the entire
 * implementation. This type implements the standard operators a signed integer type has.<BR>
 * All of int128_t's methods are inline for maximum performance.
 *
 * <B>Implementation notes:</B>
 * <UL>
 * <LI>Values are stored in two's complement, so the bit patterns of the shared operations match
 *     the unsigned type exactly.</LI>
 * <LI>Overflow is handled silently, similar to builtin integer operations.</LI>
 * <LI>An int128_t object is not thread safe. Accessing a const object from multiple threads is safe.</LI>
 * <LI>int128_t is <B>conditionally safe</B>, 2 different non const objects can be accessed concurrently.</LI>
 * <LI>Only 64 bit builds are supported.</LI>
 * <LI>Being an alias rather than a distinct class, int128_t cannot be forward declared. Include
 *     this header instead.</LI>
 * </UL>
 */
using int128_t = int128_base<true>;

/**
 * @brief User-defined literal for constructing int128_t from a string.
 *
 * Defined at namespace scope rather than as a friend of the class template: its signature does
 * not depend on the template parameter, so a friend definition inside the template would be
 * emitted once per instantiation and collide. Literal operators are found by ordinary unqualified
 * lookup only, never by ADL, so namespace scope is where they belong anyway.
 *
 * @param literal Decimal or hexadecimal digits with an optional sign, see the
 *                int128_t(const char*) constructor.
 * @return The parsed value.
 */
[[nodiscard]] FP128_INLINE int128_t operator""_int128(const char* literal)
{
    return int128_t(literal);
}

}  // namespace fp128

#endif  // FP128_INT128_T_H
