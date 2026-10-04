##
# @file CMakeLists.txt
# @brief Build script for the fixed_point128 GoogleTest suite.
#
# The tests are deliberately built for the baseline instruction set rather than with the benchmark's
# AVX2 flags: fp128_shared.h only reaches for x86 intrinsics under MSVC and supplies portable
# implementations of all of them elsewhere, so the correctness suite should exercise the code the way
# the widest range of consumers will compile it.
#
# The equivalent MSBuild project is msvc/fixed_point128_tests.vcxproj.
##

add_executable(fixed_point128_tests
    fixed_point128_gtest.cpp
    float128_gtest.cpp
    float128_accuracy_gtest.cpp
    float128_format_gtest.cpp
    float128_ieee_gtest.cpp
    int128_t_gtest.cpp
    uint128_t_gtest.cpp
    gtest_shared.h
    float128_ref_check.h
    float128_ref_data.h)

target_link_libraries(fixed_point128_tests
    PRIVATE
        fixed_point128::fixed_point128
        GTest::gtest_main)

fp128_set_warnings(fixed_point128_tests)
fp128_link_clang_builtins(fixed_point128_tests)

# gtest_discover_tests below runs this executable at build time, so a runtime the loader cannot find
# would break the build itself rather than merely the test run.
fp128_link_static_runtime(fixed_point128_tests)

# The string round-trip tests copy between fixed size buffers they own with strcpy/strlen, which the
# Windows CRT marks deprecated in favour of its non-portable _s variants. Using those would make the
# tests Windows-only for no benefit.
fp128_uses_msvc_frontend(FP128_TESTS_MSVC_FRONTEND)
if(FP128_TESTS_MSVC_FRONTEND)
    target_compile_definitions(fixed_point128_tests PRIVATE _CRT_SECURE_NO_WARNINGS)
endif()

# The debugger visualizer is consumed by the Visual Studio generators. Other generators treat an
# unknown extension as a header and silently ignore it, so only add it where it does something.
if(CMAKE_GENERATOR MATCHES "Visual Studio")
    target_sources(fixed_point128_tests PRIVATE "${CMAKE_SOURCE_DIR}/fixed_point128.natvis")
endif()

set_target_properties(fixed_point128_tests PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

include(GoogleTest)
gtest_discover_tests(fixed_point128_tests)

##
# The floating point environment of float128 is opt in, and the macro that enables it changes the
# definitions of inline functions. A translation unit built with it cannot share a program with one
# built without it, so its tests are a program of their own.
#
# The equivalent MSBuild project is msvc/fixed_point128_env_tests.vcxproj.
##

find_package(Threads REQUIRED)

add_executable(fixed_point128_env_tests
    float128_env_gtest.cpp
    gtest_shared.h
    float128_ref_check.h
    float128_ref_data.h)

target_compile_definitions(fixed_point128_env_tests PRIVATE FP128_IEEE_ENV=1)

target_link_libraries(fixed_point128_env_tests
    PRIVATE
        fixed_point128::fixed_point128
        GTest::gtest_main
        Threads::Threads)

fp128_set_warnings(fixed_point128_env_tests)
fp128_link_clang_builtins(fixed_point128_env_tests)
fp128_link_static_runtime(fixed_point128_env_tests)

if(FP128_TESTS_MSVC_FRONTEND)
    target_compile_definitions(fixed_point128_env_tests PRIVATE _CRT_SECURE_NO_WARNINGS)
endif()

set_target_properties(fixed_point128_env_tests PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

gtest_discover_tests(fixed_point128_env_tests)
