##
# @file CMakeLists.txt
# @brief Build script for the fixed_point128 benchmark executable (GCC/Clang).
#
# MSVC builds are covered by fixed_point128.slnx.
#
#   cmake -S . -B build
#   cmake --build build
#   ./bin/bench
##

cmake_minimum_required(VERSION 3.20)

project(fixed_point128_bench LANGUAGES CXX)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type" FORCE)
endif()

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/bin")

add_executable(bench src/Bench.cpp)
target_include_directories(bench PRIVATE inc)
target_compile_options(bench PRIVATE -Wall -Wextra -O3)

# Bench.cpp uses __declspec(noinline), which needs explicit opt-in outside MSVC.
include(CheckCXXCompilerFlag)
check_cxx_compiler_flag("-fdeclspec" HAVE_FDECLSPEC)
if(HAVE_FDECLSPEC)
    target_compile_options(bench PRIVATE -fdeclspec)
endif()

if(CMAKE_SYSTEM_PROCESSOR MATCHES "^(x86|X86|AMD64|amd64|x86_64|i[3-6]86)$")
    target_compile_options(bench PRIVATE -mavx2 -mbmi -mbmi2 -mlzcnt)
endif()
