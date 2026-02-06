# Documentation
- Document code using doxygen style comments.
- Clear and precide explanations.

# Build Configuration

- C++20 standard, C17 for C files
- OpenMP enabled for parallel processing
- AVX2 instruction set on Windows (MSVC)
/clee- Qt6 with widgets, core, and gui modules

# Code Style Rules
- C++20 standard, C17 for C files
- Write clear and concise comments for each function.
- Document code using Doxygen style
- Write verbose documentation for class headers.
- Use PascalCase for global functions, private and protected class/struct methods.
- Use camelCase for global functions that start with the letter 'q'. i.e. QT style.
- Use camelCase for public methods, struct methods/data members and local variables.
- Private and protected class/struct data members should start with an underscore. The rest of the identifier is camelCase. Example: "_dataMember".
- No newlines before the opening curly brace of any code block, just a single space  
  (such as after `if`, `for`, `while`, `foreach`, `using`, `try`, etc.).
- Code blocks containing a macro or function call must use curly braces.
  (such as after `if`, `for`, `while`, `foreach`, `using`, `try`, etc.).
- Ensure that the final `return` statement of a method is on its own line.
- All functions that return a value that is not an error code or boolean success value must be set as [[nodiscard]].
- Indentation is 4 spaces. No tabs.