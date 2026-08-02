##
# @file FP128ClangRuntime.cmake
# @brief Locates the Clang builtins library needed by the 128 bit division path.
#
# udiv128() divides a __uint128_t. On the Clang code path that lowers to a call to __udivti3, which
# lives in compiler-rt rather than in the C runtime. On a GNU style toolchain the driver links its own
# runtime and there is nothing to do, but clang-cl targeting the MSVC ABI does not, so the link fails
# with an undefined __udivti3 unless the builtins library is placed on the link line explicitly.
##

include_guard(GLOBAL)

##
# @brief Links the Clang compiler-rt builtins library into a target when the toolchain requires it.
#
# No-op for MSVC, and for GNU style Clang/GCC drivers that already link their own runtime.
#
# @param target Target to add the library to.
##
function(fp128_link_clang_builtins target)
    if(NOT CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        return()
    endif()

    get_filename_component(clangBin "${CMAKE_CXX_COMPILER}" DIRECTORY)
    file(GLOB clangRtCandidates
         "${clangBin}/../lib/clang/*/lib/windows/clang_rt.builtins-*.lib"
         "${clangBin}/../lib/clang/*/lib/*/libclang_rt.builtins.a"
         "${clangBin}/../lib/clang/*/lib/linux/libclang_rt.builtins-*.a")

    if(clangRtCandidates)
        list(GET clangRtCandidates 0 clangRt)
        message(STATUS "fixed_point128: linking Clang builtins for __udivti3: ${clangRt}")
        target_link_libraries(${target} PRIVATE "${clangRt}")
    elseif(NOT WIN32)
        # On a GNU style toolchain the driver links its own runtime, nothing to do.
        message(STATUS "fixed_point128: relying on the default compiler runtime for __udivti3")
    else()
        message(WARNING
            "fixed_point128: could not locate clang_rt.builtins next to ${CMAKE_CXX_COMPILER}. "
            "The link may fail with an undefined __udivti3.")
    endif()
endfunction()
