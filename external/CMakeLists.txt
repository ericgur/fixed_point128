##
# @file CMakeLists.txt
# @brief Brings the vendored third party dependencies into the build.
#
# GoogleTest is tracked as a git submodule pinned to a release tag rather than fetched at configure
# time, so that builds are reproducible and work without network access. See the README for how to
# update the pin.
##

if(NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/googletest/CMakeLists.txt")
    message(FATAL_ERROR
        "The GoogleTest submodule is not initialized, so the unit tests cannot be built.\n"
        "Run:\n"
        "    git submodule update --init --recursive\n"
        "Or configure with -DFP128_BUILD_TESTS=OFF to skip the tests entirely.")
endif()

# Use the shared (DLL) CRT, matching the /MD default of the MSBuild projects in msvc/. Mixing CRTs
# between the GoogleTest library and the test executable causes link errors and heap corruption.
set(gtest_force_shared_crt ON  CACHE BOOL "" FORCE)

# No test includes a gmock header, so building it is pure cost.
set(BUILD_GMOCK           OFF CACHE BOOL "" FORCE)

# This is a leaf project; never install GoogleTest as a side effect of installing fixed_point128.
set(INSTALL_GTEST         OFF CACHE BOOL "" FORCE)

add_subdirectory(googletest EXCLUDE_FROM_ALL)

##
# Silence GoogleTest's own diagnostics.
#
# internal_utils.cmake compiles GoogleTest with -W4 -WX on the MSVC frontend, so any new compiler
# diagnostic becomes a hard error in third-party code we do not maintain. clang-cl 21 does exactly
# that: -Wcharacter-conversion fires inside gtest-printers.h and breaks the build against a pinned,
# unmodified GoogleTest release. Since we cannot act on warnings in vendored code, they are noise at
# best, so the warning level is dropped to zero for these targets only.
#
# The flag is appended through COMPILE_OPTIONS, which CMake places after the COMPILE_FLAGS property
# that GoogleTest sets, so it wins.
##
fp128_uses_msvc_frontend(FP128_EXTERNAL_MSVC_FRONTEND)
foreach(gtestTarget IN ITEMS gtest gtest_main)
    if(TARGET ${gtestTarget})
        if(FP128_EXTERNAL_MSVC_FRONTEND)
            target_compile_options(${gtestTarget} PRIVATE /w /WX-)
        else()
            target_compile_options(${gtestTarget} PRIVATE -w -Wno-error)
        endif()
    endif()
endforeach()
