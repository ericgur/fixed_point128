@echo off
set cur_dir=%CD%
cd build

ctest

cd %cur_dir%
