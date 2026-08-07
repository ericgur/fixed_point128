##
# @file CMakeLists.txt
# @brief Build script for the fixed_point128 maintenance tools.
#
# Off by default: these exist to be run by hand when log2() is touched, not on every build. Enable
# with -DFP128_BUILD_TOOLS=ON, or build the target directly:
#
#     cmake --build --preset msvc-release --target log2_ulp_dump
#
# The Python halves of the harness need no build step at all. See README.md in this directory.
##

add_executable(log2_ulp_dump log2_ulp_dump.cpp)

target_link_libraries(log2_ulp_dump PRIVATE fixed_point128::fixed_point128)

fp128_set_warnings(log2_ulp_dump)
fp128_enable_declspec(log2_ulp_dump)
fp128_link_clang_builtins(log2_ulp_dump)
fp128_link_static_runtime(log2_ulp_dump)

##
# The same two MSVC /W4 diagnostics the benchmark and the tests suppress, reached for the same
# reasons - see bench/CMakeLists.txt for what they are and why they are deliberate.
##
fp128_uses_msvc_frontend(FP128_TOOLS_MSVC_FRONTEND)
if(FP128_TOOLS_MSVC_FRONTEND)
    target_compile_options(log2_ulp_dump PRIVATE /wd4201 /wd4702)
endif()

set_target_properties(log2_ulp_dump PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
