@echo off
rmdir /S /Q build
cmake -S . -B build -Wno-dev -G "Visual Studio 18 2026"
