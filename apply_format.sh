#!/bin/bash
find src -type f -name "*.cpp" -o -name "*.hpp" -o -name "*.h" -exec clang-format -i {}
