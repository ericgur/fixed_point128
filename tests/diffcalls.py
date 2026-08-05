#!/usr/bin/env python3
"""Compare two builds of the same binary and report differences in what got inlined.

Every function in this library is a header only `inline`, which is a *hint*: the compiler is free to
ignore it and emit a real call instead. It does exactly that under MSVC's whole program optimization
(/GL), where the link time code generator re-decides inlining and can turn a small hot operator into
an out-of-line call. The resulting slowdown does not show up as a source change, only as a benchmark
regression, which is easy to mistake for code layout noise.

This tool tells the two apart. It disassembles both binaries, and for every benchmark function it
compares the set of fp128 symbols still reached through a `call`. A symbol that is called in one
build but not the other is an inlining difference, and is usually fixable by promoting the callee
from FP128_INLINE to FP128_FORCE_INLINE. If the call sets are identical, any timing difference is
code placement, and no source change will remove it.

Typical use, comparing whole program optimization against the same build without it:

    msbuild msvc/bench.vcxproj /p:Configuration=Release /p:Platform=x64 /t:Rebuild ^
        /p:WholeProgramOptimization=false /p:OutDir=%TEMP%/noltcg/ /p:IntDir=%TEMP%/noltcg_obj/
    python tests/diffcalls.py %TEMP%/noltcg/bench.exe bin/x64_Release/bench.exe

Requires dumpbin.exe and undname.exe from a Visual Studio installation. They are located through
vswhere, the PATH, or the FP128_DUMPBIN environment variable.
"""
import os
import pathlib
import re
import shutil
import subprocess
import sys

# Symbols that say nothing about the library: the harness and the standard library.
NOISE = ("print_ips", "DoNotOptimize", "chrono", "logic_error", "domain_error",
         "Query_perf", "basic_string", "operator new", "operator delete")

# Only functions whose name contains this are compared. The benchmark names every measurement
# function bench_<operation>, so this keeps setup and reporting code out of the report.
SCOPE = "bench_"


def find_tool(name):
    """Locates a Visual Studio build tool, preferring an explicit override."""
    override = os.environ.get("FP128_DUMPBIN")
    if override:
        candidate = pathlib.Path(override).parent / name
        if candidate.is_file():
            return str(candidate)

    found = shutil.which(name)
    if found:
        return found

    vswhere = pathlib.Path(os.environ.get("ProgramFiles(x86)", r"C:\Program Files (x86)")) \
        / "Microsoft Visual Studio" / "Installer" / "vswhere.exe"
    if vswhere.is_file():
        result = subprocess.run(
            [str(vswhere), "-latest", "-products", "*", "-property", "installationPath"],
            capture_output=True, text=True).stdout.strip()
        if result:
            tools = pathlib.Path(result) / "VC" / "Tools" / "MSVC"
            if tools.is_dir():
                for version in sorted(tools.iterdir(), reverse=True):
                    candidate = version / "bin" / "Hostx64" / "x64" / name
                    if candidate.is_file():
                        return str(candidate)

    sys.exit(f"Could not find {name}. Set FP128_DUMPBIN to its full path, or run from a "
             f"Developer Command Prompt.")


DUMPBIN = find_tool("dumpbin.exe")
UNDNAME = find_tool("undname.exe")


def undecorate(symbol):
    """Turns a decorated MSVC symbol into something readable, or returns it unchanged."""
    output = subprocess.run([UNDNAME, symbol], capture_output=True, text=True,
                            errors="replace").stdout
    match = re.search(r'is :- "(.*)"', output, re.S)
    return match.group(1).strip() if match else symbol


def call_sets(binary):
    """Returns {function: {symbols it calls}} for every function matching SCOPE."""
    path = pathlib.Path(binary)
    if not path.is_file():
        sys.exit(f"No such file: {binary}")

    listing = subprocess.run([DUMPBIN, "/disasm:nobytes", str(path)],
                             capture_output=True, text=True, errors="replace").stdout
    sets, current = {}, None
    for line in listing.splitlines():
        # A function starts at a non indented line ending in a colon.
        if line and not line.startswith(" ") and line.rstrip().endswith(":"):
            current = line.rstrip()[:-1]
            if SCOPE in current:
                sets[current] = set()
            else:
                current = None
        elif current:
            match = re.search(r"\bcall\s+(\?\S+)", line)
            if match and not any(n in match.group(1) for n in NOISE):
                sets[current].add(match.group(1))
    return sets


def report(title, sets_a, sets_b):
    """Prints symbols called in sets_b but not in sets_a, most widespread first."""
    extra = {}
    for function in sorted(set(sets_a) & set(sets_b)):
        for symbol in sets_b[function] - sets_a[function]:
            extra.setdefault(symbol, []).append(function)

    print(f"\n{title}: {len(extra)} symbol(s)")
    if not extra:
        print("  none - the two builds inline the same set of functions")
        return
    for symbol, functions in sorted(extra.items(), key=lambda kv: -len(kv[1])):
        print(f"  {undecorate(symbol)}")
        print(f"      in {len(functions)} benchmark(s), e.g. {undecorate(functions[0])}")


def main():
    if len(sys.argv) != 3:
        sys.exit(f"Usage: {pathlib.Path(sys.argv[0]).name} <baseline.exe> <candidate.exe>")

    baseline, candidate = sys.argv[1], sys.argv[2]
    base_sets, cand_sets = call_sets(baseline), call_sets(candidate)
    shared = set(base_sets) & set(cand_sets)
    print(f"baseline : {baseline}  ({len(base_sets)} '{SCOPE}' functions)")
    print(f"candidate: {candidate}  ({len(cand_sets)} '{SCOPE}' functions)")
    print(f"comparing {len(shared)} function(s) present in both")

    report("Inlined in the baseline but called in the candidate (candidate regressions)",
           base_sets, cand_sets)
    report("Inlined in the candidate but called in the baseline (candidate improvements)",
           cand_sets, base_sets)

    return 0


if __name__ == "__main__":
    sys.exit(main())
