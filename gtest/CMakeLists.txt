##
# @file CMakeLists.txt
# @brief Build script for the fixed_point128 unit tests.
#
# Supports three toolchains on Windows:
#   MSVC      cmake -S . -B build        -G "Visual Studio 18 2026"
#   clang-cl  cmake -S . -B build-clang  -G "Visual Studio 18 2026" -T ClangCL
#   clang++   cmake -S . -B build-clang  -G Ninja -DCMAKE_CXX_COMPILER=clang++
#
# See setup.bat (MSVC) and setup_clang.bat (clang-cl).
##

cmake_minimum_required(VERSION 3.14)
project(fixed_point128_gtest CXX)

# GoogleTest requires at least C++14, this library requires C++20
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
# use static CRT
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

set(BUILD_GMOCK ON)
add_subdirectory(external/googletest)

enable_testing()

add_executable(
  fixed_point128_gtest
  fixed_point128_gtest.cpp
  uint128_t_gtest.cpp
  int128_t_gtest.cpp
  float128_gtest.cpp
)

# The debugger visualizer is consumed by the Visual Studio generators. Other generators treat an
# unknown extension as a header and silently ignore it, so only add it where it does something.
if(CMAKE_GENERATOR MATCHES "Visual Studio")
  target_sources(fixed_point128_gtest PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../fixed_point128.natvis)
endif()

##
# Clang support.
#
# No instruction set flags are needed: fixed_point128_shared.h only reaches for the x86 intrinsics
# under MSVC and provides portable implementations of all of them for Clang, so the headers compile
# with the baseline target.
#
# What Clang does need is the compiler runtime. udiv128() divides a __uint128_t, which lowers to a
# call to __udivti3 and lives in compiler-rt rather than in the MSVC runtime, so linking fails with
# an undefined symbol unless the builtins library is on the link line.
##
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  get_filename_component(FP128_CLANG_BIN "${CMAKE_CXX_COMPILER}" DIRECTORY)
  file(GLOB FP128_CLANG_RT_CANDIDATES
       "${FP128_CLANG_BIN}/../lib/clang/*/lib/windows/clang_rt.builtins-*.lib"
       "${FP128_CLANG_BIN}/../lib/clang/*/lib/*/libclang_rt.builtins.a"
       "${FP128_CLANG_BIN}/../lib/clang/*/lib/linux/libclang_rt.builtins-*.a")

  if(FP128_CLANG_RT_CANDIDATES)
    list(GET FP128_CLANG_RT_CANDIDATES 0 FP128_CLANG_RT)
    message(STATUS "fixed_point128: linking Clang builtins for __udivti3: ${FP128_CLANG_RT}")
    target_link_libraries(fixed_point128_gtest PRIVATE "${FP128_CLANG_RT}")
  elseif(NOT WIN32)
    # on a GNU style toolchain the driver links its own runtime, nothing to do
    message(STATUS "fixed_point128: relying on the default compiler runtime for __udivti3")
  else()
    message(WARNING
      "fixed_point128: could not locate clang_rt.builtins next to ${CMAKE_CXX_COMPILER}. "
      "The link may fail with an undefined __udivti3.")
  endif()
endif()

target_link_libraries(
  fixed_point128_gtest
  PRIVATE
  GTest::gtest_main
)

include(GoogleTest)

gtest_discover_tests(fixed_point128_gtest)
