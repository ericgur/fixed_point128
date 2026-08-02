##
# @file CMakeLists.txt
# @brief Build script for the fixed_point128 benchmark executable.
#
# Unlike the unit tests, the benchmark opts into the AVX2 instruction set family: it exists to measure
# the library at its fastest, and the resulting binary is not expected to run on pre-Haswell hardware.
##

add_executable(bench Bench.cpp)

target_link_libraries(bench PRIVATE fixed_point128::fixed_point128)

fp128_set_warnings(bench)
fp128_set_arch_flags(bench)
fp128_enable_declspec(bench)
fp128_link_clang_builtins(bench)
fp128_link_static_runtime(bench)

set_target_properties(bench PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
