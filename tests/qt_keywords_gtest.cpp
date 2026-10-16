/**
 * @file qt_keywords_gtest.cpp
 * @brief Builds every library header the way a Qt project sees it: with Qt's macros defined.
 *
 * Qt defines its keywords (`emit`, `signals`, `slots`, `foreach`, `forever`) and a set of other
 * unprefixed names as macros, and a Qt project includes those definitions before it includes ours.
 * A library identifier that matches one of them is rewritten by the preprocessor, which breaks the
 * build or, worse, quietly changes the code: Qt's `emit` expands to nothing, so a call `emit(x)`
 * becomes the expression `(x)`.
 *
 * This translation unit defines each of those names as a stray `)` and then includes the library
 * headers. No declaration, statement or expression can absorb a stray `)`, so any use of a Qt name
 * as an identifier fails the build: Clang and GCC name the macro in a note under the error, MSVC
 * reports a syntax error on the line that used it. The check runs at compile time; the test case
 * exists so that it shows up in the test list.
 *
 * The names are every unprefixed macro the public headers of Qt 6.12 define for user code, apart
 * from two that only exist on obsolete toolchains (`for` and `main`). Names starting with `Q_`,
 * `QT_` or `qt_` cannot be listed exhaustively; CLAUDE.md rules them out for library identifiers.
 *
 * Only the code each compiler actually builds is checked: an identifier inside a branch of
 * fp128_shared.h for another compiler or architecture is checked when that compiler builds this
 * file.
 */

// remove warnings from gtest itself
#if defined(_MSC_VER)
#pragma warning(push)
#pragma warning(disable : 26439)
#pragma warning(disable : 26495)
#endif
#include <gtest/gtest.h>
#if defined(_MSC_VER)
#pragma warning(pop)
#endif

// The standard headers the library includes come in before the Qt names are defined. They are not
// ours to check, and <cmath> defines the M_ constants itself on some platforms.
#include <algorithm>
#include <array>
#include <bit>
#include <cassert>
#include <cctype>
#include <cfenv>
#include <cfloat>
#include <charconv>
#include <climits>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <format>
#include <functional>
#include <istream>
#include <limits>
#include <memory>
#include <new>
#include <ostream>
#include <stdexcept>
#include <string>
#include <type_traits>

/// @brief What every Qt name expands to in this translation unit: a token nothing can absorb.
#define FP128_QT_NAME_USED_AS_IDENTIFIER )

// Keywords: qtmetamacros.h and qforeach.h.
#define signals FP128_QT_NAME_USED_AS_IDENTIFIER
#define slots FP128_QT_NAME_USED_AS_IDENTIFIER
#define emit FP128_QT_NAME_USED_AS_IDENTIFIER
#define foreach FP128_QT_NAME_USED_AS_IDENTIFIER
#define forever FP128_QT_NAME_USED_AS_IDENTIFIER

// Meta-object connection strings: qobjectdefs.h. Qt defines these as function-like macros; they are
// object-like here so that the bare name is caught too.
#define SIGNAL FP128_QT_NAME_USED_AS_IDENTIFIER
#define SLOT FP128_QT_NAME_USED_AS_IDENTIFIER
#define METHOD FP128_QT_NAME_USED_AS_IDENTIFIER

// Logging: qlogging.h and qloggingcategory.h.
#define qDebug FP128_QT_NAME_USED_AS_IDENTIFIER
#define qInfo FP128_QT_NAME_USED_AS_IDENTIFIER
#define qWarning FP128_QT_NAME_USED_AS_IDENTIFIER
#define qCritical FP128_QT_NAME_USED_AS_IDENTIFIER
#define qFatal FP128_QT_NAME_USED_AS_IDENTIFIER
#define qErrnoWarning FP128_QT_NAME_USED_AS_IDENTIFIER
#define qCDebug FP128_QT_NAME_USED_AS_IDENTIFIER
#define qCInfo FP128_QT_NAME_USED_AS_IDENTIFIER
#define qCWarning FP128_QT_NAME_USED_AS_IDENTIFIER
#define qCCritical FP128_QT_NAME_USED_AS_IDENTIFIER
#define qCFatal FP128_QT_NAME_USED_AS_IDENTIFIER

// Strings, moves and the application instance: qstring.h, qcompilerdetection.h,
// qcoreapplication.h and qguiapplication.h.
#define qPrintable FP128_QT_NAME_USED_AS_IDENTIFIER
#define qUtf8Printable FP128_QT_NAME_USED_AS_IDENTIFIER
#define qUtf16Printable FP128_QT_NAME_USED_AS_IDENTIFIER
#define qMove FP128_QT_NAME_USED_AS_IDENTIFIER
#define qApp FP128_QT_NAME_USED_AS_IDENTIFIER
#define qGuiApp FP128_QT_NAME_USED_AS_IDENTIFIER

// Math constants: qmath.h defines each one that <cmath> has not, so <cmath>'s are replaced here.
#undef M_E
#undef M_LOG2E
#undef M_LOG10E
#undef M_LN2
#undef M_LN10
#undef M_PI
#undef M_PI_2
#undef M_PI_4
#undef M_1_PI
#undef M_2_PI
#undef M_2_SQRTPI
#undef M_SQRT2
#undef M_SQRT1_2
#define M_E FP128_QT_NAME_USED_AS_IDENTIFIER
#define M_LOG2E FP128_QT_NAME_USED_AS_IDENTIFIER
#define M_LOG10E FP128_QT_NAME_USED_AS_IDENTIFIER
#define M_LN2 FP128_QT_NAME_USED_AS_IDENTIFIER
#define M_LN10 FP128_QT_NAME_USED_AS_IDENTIFIER
#define M_PI FP128_QT_NAME_USED_AS_IDENTIFIER
#define M_PI_2 FP128_QT_NAME_USED_AS_IDENTIFIER
#define M_PI_4 FP128_QT_NAME_USED_AS_IDENTIFIER
#define M_1_PI FP128_QT_NAME_USED_AS_IDENTIFIER
#define M_2_PI FP128_QT_NAME_USED_AS_IDENTIFIER
#define M_2_SQRTPI FP128_QT_NAME_USED_AS_IDENTIFIER
#define M_SQRT2 FP128_QT_NAME_USED_AS_IDENTIFIER
#define M_SQRT1_2 FP128_QT_NAME_USED_AS_IDENTIFIER

// Every public header, which between them include every other library header.
#include "fixed_point128.h"
#include "float128.h"
#include "int128_t.h"
#include "uint128_t.h"

/// @brief The headers built with every Qt name defined, and the decimal formatter, whose digit sink
///        was once named `emit`, still produces the digits.
TEST(QtKeywords, HeadersBuildWithQtMacrosDefined)
{
    char buffer[64];
    const std::to_chars_result result =
        fp128::to_chars(buffer, buffer + sizeof(buffer), fp128::float128(0.375), std::chars_format::fixed, 3);
    ASSERT_EQ(result.ec, std::errc {});
    EXPECT_EQ(std::string(buffer, result.ptr), "0.375");
}
