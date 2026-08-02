Get-ChildItem -Path include, bench, tests -Include *.cpp, *.hpp, *.h -File -Recurse | ForEach-Object { clang-format -i $_.FullName }
