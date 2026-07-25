Install CMake and add it to the path (install option).

MSVC build
One time:
1) Setup the build system by running setup.bat. This creates a VS solution in build\fixed_point128_gtest.sln

After code changes:
1) Build the test app via build.bat
2a) Run all the tests via test.bat
2b) Run only the failed tests with run_failed.bat

Clang build
Uses the Clang toolset that ships with Visual Studio (install the "C++ Clang tools for Windows"
component if -T ClangCL fails to configure). It builds into build_clang\, so both toolchains can
be kept side by side.
One time:
1) Setup the build system by running setup_clang.bat

After code changes:
1) Build the test app via build_clang.bat
2) Run all the tests via test_clang.bat

The headers need no instruction set flags under Clang: fixed_point128_shared.h only uses the x86
intrinsics on MSVC and provides portable implementations for Clang. The one extra requirement is
the Clang builtins library, which supplies __udivti3 for the __uint128_t division in udiv128().
CMakeLists.txt locates it next to the compiler and adds it to the link line automatically.
