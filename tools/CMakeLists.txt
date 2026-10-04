##
# @file CMakeLists.txt
# @brief Build script for the fixed_point128 maintenance tools.
#
# Off by default: these exist to be run by hand when log2() or the float128 arithmetic is touched,
# not on every build. Enable with -DFP128_BUILD_TOOLS=ON, or build the target directly:
#
#     cmake --build --preset msvc-release --target log2_ulp_dump
#     cmake --build --preset msvc-release --target ieee_probe ieee_probe_env
#
# The Python halves of the harness need no build step at all. See README.md in this directory.
##

add_executable(log2_ulp_dump log2_ulp_dump.cpp)

target_link_libraries(log2_ulp_dump PRIVATE fixed_point128::fixed_point128)

fp128_set_warnings(log2_ulp_dump)
fp128_enable_declspec(log2_ulp_dump)
fp128_link_clang_builtins(log2_ulp_dump)
fp128_link_static_runtime(log2_ulp_dump)

set_target_properties(log2_ulp_dump PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

# The C++ half of the IEEE 754 differential check, twice: as the library is built by default, and
# with the opt in floating point environment, whose rounding directions and flags ieee_check.py --env
# exercises. The macro changes inline functions, so the two have to be separate programs.
foreach(probe ieee_probe ieee_probe_env)
    add_executable(${probe} ieee_probe.cpp)
    target_link_libraries(${probe} PRIVATE fixed_point128::fixed_point128)
    fp128_set_warnings(${probe})
    fp128_link_clang_builtins(${probe})
    fp128_link_static_runtime(${probe})
    set_target_properties(${probe} PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
endforeach()
target_compile_definitions(ieee_probe_env PRIVATE FP128_IEEE_ENV=1)
