@echo off
cls
set cur_dir=%CD%
cd build

ctest --output-on-failure --parallel %NUMBER_OF_PROCESSORS%

cd %cur_dir%
