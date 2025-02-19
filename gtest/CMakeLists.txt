cmake_minimum_required(VERSION 3.14)
project(fixed_point128_gtest)

# GoogleTest requires at least C++14
set(CMAKE_CXX_STANDARD 20)
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

target_link_libraries(
  fixed_point128_gtest
  GTest::gtest_main
)

include(GoogleTest)

gtest_discover_tests(fixed_point128_gtest)
