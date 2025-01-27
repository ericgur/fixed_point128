@echo off
cls
pushd build
ctest --output-on-failure --parallel %NUMBER_OF_PROCESSORS%
popd
