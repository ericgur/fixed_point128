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

##
# Two MSVC /W4 diagnostics fire on deliberate constructs in the library headers, both reached because
# the benchmark measures every 128 bit type. They are suppressed here rather than on the
# fixed_point128 INTERFACE target, so that consumers of the library keep full control of their own
# warning settings. The unit tests suppress the same pair, for the same reasons.
#
#   C4201  float128 stores its 128 bits in an anonymous struct wrapping an anonymous union. The layout
#          is the point of the type, and the extension is universally supported.
#   C4702  MSVC reports the code after an `if constexpr` chain as unreachable for the instantiations
#          where an earlier branch returns, even though it is live for all the others. A compiler
#          limitation rather than dead code; only the debug build reaches it.
##
fp128_uses_msvc_frontend(FP128_BENCH_MSVC_FRONTEND)
if(FP128_BENCH_MSVC_FRONTEND)
    target_compile_options(bench PRIVATE /wd4201 /wd4702)
endif()

set_target_properties(bench PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
