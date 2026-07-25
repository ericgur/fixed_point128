@echo off
cls
cd /d "%~dp0"
pushd build_clang
ctest --output-on-failure --parallel %NUMBER_OF_PROCESSORS% %*
popd
