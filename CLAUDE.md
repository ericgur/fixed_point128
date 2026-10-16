# Documentation
- Document code using doxygen style comments.
- Clear and precise explanations.

# Build Configuration

- C++20 standard, C17 for C files
- OpenMP enabled for parallel processing
- AVX2 instruction set on Windows (MSVC)
- Code must compile without warnings on MSVC and Clang.

# Code Style Rules
- C++20 standard, C17 for C files
- Write clear and concise comments for each function.
- Document code using Doxygen style
- Write verbose documentation for class headers.
- Use PascalCase for global functions, private and protected class/struct methods.
- Use camelCase for global functions that start with the letter 'q'. i.e. QT style.
- Never name anything after a Qt keyword or macro: no variable, parameter, lambda, function, type,
  member, namespace, enumerator or macro. Qt projects define these names as macros before they include
  our headers, so any library identifier that matches one breaks their build. A local lambda named
  `emit` did exactly that. The names, taken from the Qt 6 public headers:
  - Keywords: `signals`, `slots`, `emit`, `foreach`, `forever`.
  - Meta-object: `SIGNAL`, `SLOT`, `METHOD`.
  - Logging: `qDebug`, `qInfo`, `qWarning`, `qCritical`, `qFatal`, `qErrnoWarning`, `qCDebug`,
    `qCInfo`, `qCWarning`, `qCCritical`, `qCFatal`.
  - Other: `qPrintable`, `qUtf8Printable`, `qUtf16Printable`, `qMove`, `qApp`, `qGuiApp`.
  - Math constants from qmath.h: `M_E`, `M_LOG2E`, `M_LOG10E`, `M_LN2`, `M_LN10`, `M_PI`, `M_PI_2`,
    `M_PI_4`, `M_1_PI`, `M_2_PI`, `M_2_SQRTPI`, `M_SQRT2`, `M_SQRT1_2`. Avoid every `M_` name.
  - Prefixes: no identifier may start with `Q_`, `QT_` or `qt_`. Qt names its own macros that way.

  `tests/qt_keywords_gtest.cpp` defines each listed name as a compile error. When this list
  changes, change that test to match.
- Use camelCase for public methods, struct methods/data members and local variables.
- Private and protected class/struct data members use camelCase. Example: "dataMember".
- No newlines before the opening curly brace of any code block, just a single space  
  (such as after `if`, `for`, `while`, `foreach`, `using`, `try`, etc.).
- Code blocks containing a macro or function call must use curly braces.
  (such as after `if`, `for`, `while`, `foreach`, `using`, `try`, etc.).
- Ensure that the final `return` statement of a method is on its own line.
- All functions that return a value that is not an error code or boolean success value must be set as [[nodiscard]].
- Indentation is 4 spaces. No tabs.
- All source files must be UTF-8 encoded. Convert if they don't meet this requirement.
- Line ending must be LF (Unix style). Convert if they don't meet this requirement.
