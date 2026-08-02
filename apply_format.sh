#!/bin/bash
find include bench tests -type f -name "*.cpp" -o -name "*.hpp" -o -name "*.h" -exec clang-format -i {}
