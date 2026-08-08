##
# @file FP128ClangRuntime.cmake
# @brief Runtime library handling for the Clang toolchains.
#
# The two frontends need opposite things, and both of them only at link time:
#
#   fp128_link_clang_builtins()  clang-cl, to resolve __udivti3
#   fp128_link_static_runtime()  MinGW targeting drivers, so the binaries run outside the toolchain
#
# udiv128() divides a __uint128_t. On the Clang code path that lowers to a call to __udivti3, which
# lives in compiler-rt rather than in the C runtime.
#
# Whether anything has to be done about that depends on the driver, not on the platform. A GNU style
# driver links its own compiler runtime and resolves the call unaided - that covers Clang on Linux and
# macOS, and equally the MinGW targeting Clang distributions on Windows. clang-cl targeting the MSVC
# ABI links no runtime of its own, so there the link fails with an undefined __udivti3 unless the
# builtins library is named on the link line explicitly.
#
# That same GNU style driver on Windows creates the opposite problem. It resolves the call, but it does
# so against DLLs - libc++.dll and libunwind.dll for llvm-mingw, libstdc++-6.dll and libgcc_s_seh-1.dll
# for MinGW GCC - that ship inside the toolchain's own bin directory and are not on PATH. The link
# succeeds and the executable then dies at startup with STATUS_DLL_NOT_FOUND (0xC0000135) for anyone
# who has not put the compiler on their PATH, ctest included. MSVC and clang-cl are unaffected: they
# import only the UCRT, which is part of the operating system.
##

include_guard(GLOBAL)

include(CheckLinkerFlag)
include("${CMAKE_CURRENT_LIST_DIR}/FP128CompilerOptions.cmake")

##
# @brief Resolves the full path of the compiler-rt builtins library for the current toolchain.
#
# The driver is asked for the answer with `-print-libgcc-file-name`, which reports the library for the
# target the driver actually defaults to and so never confuses a 32 bit build with a 64 bit one. Only
# if that fails - an old driver that does not understand the option, or a distribution that ships the
# compiler without compiler-rt - is the compiler-rt directory beside the compiler searched by hand.
#
# The lookup is cached, because the function runs once per target it is applied to and the answer
# cannot change unless the compiler does.
#
# @param outVar Name of the variable to receive the library path, or a false value when none was found.
##
function(fp128_find_clang_builtins outVar)
    if(DEFINED FP128_CLANG_BUILTINS)
        set(${outVar} "${FP128_CLANG_BUILTINS}" PARENT_SCOPE)
        return()
    endif()

    set(builtins "")

    execute_process(
        COMMAND "${CMAKE_CXX_COMPILER}" -print-libgcc-file-name --rtlib=compiler-rt
        RESULT_VARIABLE printResult
        OUTPUT_VARIABLE printedPath
        ERROR_QUIET
        OUTPUT_STRIP_TRAILING_WHITESPACE)

    if(printResult EQUAL 0 AND EXISTS "${printedPath}")
        set(builtins "${printedPath}")
    else()
        get_filename_component(clangBin "${CMAKE_CXX_COMPILER}" DIRECTORY)
        file(GLOB candidates
             "${clangBin}/../lib/clang/*/lib/windows/clang_rt.builtins-*.lib"
             "${clangBin}/../lib/clang/*/lib/*/clang_rt.builtins.lib")

        # Several distributions ship the 32 and 64 bit builtins side by side, and an alphabetical pick
        # would take i386 for a 64 bit build, so discard the candidates of the wrong width first.
        if(CMAKE_SIZEOF_VOID_P EQUAL 8)
            list(FILTER candidates EXCLUDE REGEX "-(i[3-6]86|arm)\\.lib$")
        else()
            list(FILTER candidates EXCLUDE REGEX "-(x86_64|aarch64)\\.lib$")
        endif()

        if(candidates)
            list(GET candidates 0 builtins)
        endif()
    endif()

    set(FP128_CLANG_BUILTINS "${builtins}" CACHE INTERNAL
        "compiler-rt builtins library supplying __udivti3, empty when none was located")
    set(${outVar} "${builtins}" PARENT_SCOPE)
endfunction()

##
# @brief Links the Clang compiler-rt builtins library into a target when the toolchain requires it.
#
# No-op for MSVC, and for the GNU style Clang drivers that already link their own runtime.
#
# @param target Target to add the library to.
##
function(fp128_link_clang_builtins target)
    if(NOT CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        return()
    endif()

    fp128_uses_msvc_frontend(msvcFrontend)
    if(NOT msvcFrontend)
        # A GNU style driver links its own runtime, on Windows just as much as anywhere else.
        return()
    endif()

    fp128_find_clang_builtins(clangRt)
    if(clangRt)
        target_link_libraries(${target} PRIVATE "${clangRt}")
    else()
        message(WARNING
            "fixed_point128: could not locate clang_rt.builtins for ${CMAKE_CXX_COMPILER}. "
            "The link may fail with an undefined __udivti3.")
    endif()
endfunction()

##
# @brief Links the C++ runtime statically for the MinGW targeting drivers on Windows.
#
# Applied only where the runtime would otherwise be a set of DLLs private to the toolchain: a GNU style
# driver producing a Windows binary. The condition is deliberately on the frontend rather than on the
# compiler id, so that it covers llvm-mingw Clang and MinGW GCC alike while leaving clang-cl and MSVC -
# which have no such problem - untouched.
#
# `-static` rather than `-static-libstdc++`, because the latter leaves the unwinder behind as a DLL and
# the executable still fails to start. Elsewhere this is skipped entirely: on Linux and macOS the
# runtime is a system library, and a fully static link there is either discouraged or unsupported.
#
# The flag is probed instead of assumed, since a distribution that ships only the import libraries
# cannot honour it, and a link error at that point would be far less obvious than skipping it here.
#
# @param target Target to apply the link option to.
##
function(fp128_link_static_runtime target)
    if(NOT WIN32)
        return()
    endif()

    fp128_uses_msvc_frontend(msvcFrontend)
    if(msvcFrontend)
        return()
    endif()

    check_linker_flag(CXX "-static" FP128_HAVE_STATIC_RUNTIME)
    if(FP128_HAVE_STATIC_RUNTIME)
        target_link_options(${target} PRIVATE -static)
    else()
        message(WARNING
            "fixed_point128: ${CMAKE_CXX_COMPILER} rejected -static. ${target} will import the "
            "toolchain's runtime DLLs and will only run with the compiler's bin directory on PATH.")
    endif()
endfunction()
