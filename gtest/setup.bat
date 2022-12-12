@echo off
rmdir /S /Q build
cmake -S . -B build -Wno-dev
:: copy the natvis files for easy debugging of the test app
copy ../*.natvis build
