##
# @file CMakeLists.txt
# @brief Top level build script for the fixed_point128 library, benchmark and unit tests.
#
# fixed_point128 is header only, so the library itself is exposed as an INTERFACE target that simply
# carries the include directory and the C++20 requirement. The benchmark and the GoogleTest suite are
# ordinary executables that link against it.
#
# Supported toolchains: MSVC, clang-cl, Clang and GCC, on Windows, macOS and Linux. Use the presets in
# CMakePresets.json rather than invoking cmake by hand:
#
#     cmake --preset msvc                  # configure; or clang-cl / clang / gcc
#     cmake --build --preset msvc-release  # build; or msvc-debug
#     ctest --preset msvc-release          # test;  or msvc-debug
#
# Every configure preset uses a multi-config generator (Visual Studio, or Ninja Multi-Config), so the
# configure step is configuration agnostic - there is no `msvc-debug` *configure* preset, and asking for
# one is an error. Debug versus Release is chosen at build and test time, which is why `-debug` and
# `-release` only ever appear on `cmake --build --preset` and `ctest --preset`.
#
# For the same reason, building the directory directly instead of through a build preset needs an
# explicit configuration: bare `cmake --build out/build/msvc` silently builds Debug, the first entry in
# CMAKE_CONFIGURATION_TYPES. Pass `--config Release` if you go that route.
#
# The unit tests need the GoogleTest submodule:
#
#     git submodule update --init --recursive
#
# The MSVC solution in msvc/ builds the same sources through MSBuild and is maintained in parallel.
##

cmake_minimum_required(VERSION 3.21)

project(fixed_point128
        VERSION 1.0.0
        DESCRIPTION "128 bit integer, fixed-point and floating-point arithmetic types"
        HOMEPAGE_URL "https://github.com/ericgur/fixed_point128"
        LANGUAGES CXX)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(FP128CompilerOptions)
include(FP128ClangRuntime)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(FP128_BUILD_BENCH "Build the fixed_point128 benchmark" ${PROJECT_IS_TOP_LEVEL})
option(FP128_BUILD_TESTS "Build the fixed_point128 unit tests" ${PROJECT_IS_TOP_LEVEL})
option(FP128_BUILD_TOOLS "Build the fixed_point128 maintenance tools (see tools/README.md)" OFF)

##
# The library target.
#
# Header only, so INTERFACE: consumers inherit the include directory and the C++20 requirement and
# nothing is compiled here. This is also what lets the tests and the benchmark use portable
# `#include <fixed_point128.h>` instead of relative paths into the source tree.
##
add_library(fixed_point128 INTERFACE)
add_library(fixed_point128::fixed_point128 ALIAS fixed_point128)

target_include_directories(fixed_point128 INTERFACE
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
    "$<INSTALL_INTERFACE:include>")

target_compile_features(fixed_point128 INTERFACE cxx_std_20)

if(FP128_BUILD_BENCH)
    add_subdirectory(bench)
endif()

if(FP128_BUILD_TESTS)
    enable_testing()
    add_subdirectory(external)
    add_subdirectory(tests)
endif()

##
# The tools are added unconditionally so the target exists and can be built on demand, but they are
# excluded from the default build unless FP128_BUILD_TOOLS is on. They are run by hand when log2()
# changes, and there is no reason to pay for them on an ordinary build.
##
add_subdirectory(tools EXCLUDE_FROM_ALL)
if(FP128_BUILD_TOOLS)
    set_target_properties(log2_ulp_dump PROPERTIES EXCLUDE_FROM_ALL OFF)
endif()
