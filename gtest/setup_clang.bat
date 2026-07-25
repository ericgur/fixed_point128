@echo off
rem Configures the unit tests against the Clang toolset that ships with Visual Studio.
rem Uses a separate build directory so the MSVC build (setup.bat) stays intact.
cd /d "%~dp0"
rmdir /S /Q build_clang
cmake -S . -B build_clang -Wno-dev -G "Visual Studio 18 2026" -T ClangCL
