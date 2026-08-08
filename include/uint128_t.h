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
 * @file uint128_t.h
 * @brief 128-bit unsigned integer type.
 *
 * Provides @ref fp128::uint128_t, a software-emulated 128-bit unsigned integer
 * with full arithmetic, bitwise, comparison, and conversion operators, as well as
 * math functions (sqrt, log, log2, log10, pow).
 * All methods are inline for maximum performance.
 *
 * @see int128_shared.h for the class template both 128 bit integer types alias.
 * @see fp128_shared.h for supporting intrinsics and utilities.
 */

#ifndef FP128_UINT128_T_H
#define FP128_UINT128_T_H

#include "int128_shared.h"

namespace fp128
{

/***********************************************************************************
 *                                  Main Code
 ************************************************************************************/

/**
 * @typedef uint128_t
 * @brief 128 bit unsigned integer type.
 *
 * An alias of the unsigned instantiation of @ref fp128::int128_base, which holds the entire
 * implementation. This type implements the standard operators an unsigned integer type has.<BR>
 * All of uint128_t's methods are inline for maximum performance.
 *
 * <B>Implementation notes:</B>
 * <UL>
 * <LI>Overflow is handled silently, similar to builtin integer operations.</LI>
 * <LI>A uint128_t object is not thread safe. Accessing a const object from multiple threads is safe.</LI>
 * <LI>uint128_t is <B>conditionally safe</B>, 2 different non const objects can be accessed concurrently.</LI>
 * <LI>Only 64 bit builds are supported.</LI>
 * <LI>Being an alias rather than a distinct class, uint128_t cannot be forward declared. Include
 *     this header instead.</LI>
 * </UL>
 */
using uint128_t = int128_base<false>;

/**
 * @brief User-defined literal for constructing uint128_t from a string.
 *
 * Defined at namespace scope rather than as a friend of the class template: its signature does
 * not depend on the template parameter, so a friend definition inside the template would be
 * emitted once per instantiation and collide. Literal operators are found by ordinary unqualified
 * lookup only, never by ADL, so namespace scope is where they belong anyway.
 *
 * @param literal Decimal or hexadecimal digits, see the uint128_t(const char*) constructor.
 * @return The parsed value.
 */
[[nodiscard]] FP128_INLINE uint128_t operator""_uint128(const char* literal)
{
    return uint128_t(literal);
}

}  // namespace fp128

#endif  // FP128_UINT128_T_H
