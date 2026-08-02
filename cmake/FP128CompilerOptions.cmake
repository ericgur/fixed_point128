##
# @file FP128CompilerOptions.cmake
# @brief Toolchain aware compiler option helpers shared by the benchmark and the unit tests.
#
# The four supported toolchains do not agree on flag syntax, and two of them disagree with their own
# compiler id: clang-cl reports `CMAKE_CXX_COMPILER_ID` as "Clang" but accepts MSVC style flags. Every
# helper here therefore branches on `CMAKE_CXX_COMPILER_FRONTEND_VARIANT` (the flag *syntax*) rather
# than on the compiler id (the code *generator*), which is the only distinction that matters when
# deciding between `/W4` and `-Wall`.
#
# Provided functions:
#   fp128_set_warnings(<target>)     enable the project warning level
#   fp128_set_arch_flags(<target>)   opt into AVX2/BMI/LZCNT on x86 hosts
#   fp128_enable_declspec(<target>)  allow __declspec on non-MSVC frontends
##

include_guard(GLOBAL)
include(CheckCXXCompilerFlag)

##
# @brief Returns TRUE when the compiler accepts MSVC style command line flags.
#
# True for both cl.exe and clang-cl. `CMAKE_CXX_COMPILER_FRONTEND_VARIANT` is empty for compilers that
# only ever had one frontend, so fall back to the MSVC variable in that case.
#
# @param outVar Name of the variable to receive the result.
##
function(fp128_uses_msvc_frontend outVar)
    if(CMAKE_CXX_COMPILER_FRONTEND_VARIANT)
        if(CMAKE_CXX_COMPILER_FRONTEND_VARIANT STREQUAL "MSVC")
            set(${outVar} TRUE PARENT_SCOPE)
        else()
            set(${outVar} FALSE PARENT_SCOPE)
        endif()
    elseif(MSVC)
        set(${outVar} TRUE PARENT_SCOPE)
    else()
        set(${outVar} FALSE PARENT_SCOPE)
    endif()
endfunction()

##
# @brief Enables the project warning level on a target.
#
# The codebase is required to compile warning free on both MSVC and Clang, so this is deliberately
# strict. `/utf-8` is passed on the MSVC frontend because all sources in this repository are UTF-8 and
# MSVC otherwise decodes them using the system ANSI code page.
#
# @param target Target to apply the options to.
##
function(fp128_set_warnings target)
    fp128_uses_msvc_frontend(msvcFrontend)
    if(msvcFrontend)
        target_compile_options(${target} PRIVATE /W4 /utf-8)
    else()
        target_compile_options(${target} PRIVATE -Wall -Wextra)
    endif()
endfunction()

##
# @brief Opts a target into the AVX2 instruction set family on x86 hosts.
#
# The headers never require these instructions - fixed_point128_shared.h provides a portable fallback
# for every intrinsic it uses - so this is a performance opt-in only, and is applied to the benchmark
# rather than to the unit tests. Non-x86 targets (notably Apple Silicon) are skipped: the flags do not
# exist there.
#
# MSVC's `/arch:AVX2` implies BMI, BMI2, LZCNT and FMA. clang-cl maps `/arch:AVX2` to `-mavx2` alone,
# so the GNU style flags are used for every Clang frontend to get the same instruction set.
#
# @param target Target to apply the options to.
##
function(fp128_set_arch_flags target)
    if(NOT CMAKE_SYSTEM_PROCESSOR MATCHES "^(x86|X86|AMD64|amd64|x86_64|i[3-6]86)$")
        return()
    endif()

    if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
        target_compile_options(${target} PRIVATE /arch:AVX2)
    else()
        target_compile_options(${target} PRIVATE -mavx2 -mbmi -mbmi2 -mlzcnt)
    endif()
endfunction()

##
# @brief Enables __declspec support on compilers where it is not native.
#
# Bench.cpp annotates its measurement helpers with __declspec(noinline) to keep the optimizer from
# folding away the code being timed. Clang understands the attribute but only accepts it behind
# -fdeclspec; the flag is probed rather than assumed because GCC has no equivalent.
#
# @param target Target to apply the option to.
##
function(fp128_enable_declspec target)
    check_cxx_compiler_flag("-fdeclspec" FP128_HAVE_FDECLSPEC)
    if(FP128_HAVE_FDECLSPEC)
        target_compile_options(${target} PRIVATE -fdeclspec)
    endif()
endfunction()
