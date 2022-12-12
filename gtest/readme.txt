Install CMake and add it to the path (install option).
One time:
1) Setup the build system by running setup.bat
2) Open the generated solution build\fixed_point128_gtest.sln 
3) Modify the runtime library for fixed_point128_gtest to use the static variant (not DLL)

After code changes:
1) Build the test app via build.bat
2) Run the tests via test.bat

