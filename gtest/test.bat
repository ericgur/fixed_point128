@echo off
cls
set cur_dir=%CD%
cd build

ctest --output-on-failure

cd %cur_dir%
