@echo off
rmdir /S /Q build
cmake -S . -B build -Wno-dev
