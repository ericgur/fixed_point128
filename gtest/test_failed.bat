@echo off
cls
set cur_dir=%CD%
cd build

ctest --output-on-failure --rerun-failed

cd %cur_dir%
