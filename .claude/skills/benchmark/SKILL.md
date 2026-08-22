---
name: benchmark
description: Build and run the fixed_point128 release benchmarks for MSVC and Clang, compare them against the stored baseline, and generate an HTML report of the differences. Use when asked to benchmark, measure performance, check for a performance regression, compare against the baseline, or update the benchmark baseline.
---

# fixed_point128 benchmark

Measures the two release configurations the project tracks on Windows and diffs them against a
baseline held outside version control.

| configuration | preset          | compiler                                       |
| ------------- | --------------- | ---------------------------------------------- |
| `msvc`        | `msvc-release`  | cl.exe, Visual Studio generator                |
| `clang`       | `clang-release` | command line clang++, Ninja Multi-Config       |

clang-cl is not measured. It has never diverged from clang on this benchmark by more than the run to
run noise, so measuring it costs a third of the wall clock time and adds nothing.

## Running it

```powershell
.claude\skills\benchmark\scripts\Invoke-Fp128Benchmark.ps1
```

Builds both configurations and measures each of them on **every core type the CPU has** - on a hybrid
part that means one run pinned to a P-core and one to an E-core, reported separately. Writes
`bench/results/report.html`. Roughly 4 minutes per configuration per core type at the default 5 runs.

A P-core and an E-core are different machines for this purpose (2.6x apart on this benchmark), so
their numbers are never mixed. On a uniform CPU there is a single set of runs and a single section.

| flag                      | effect                                                              |
| ------------------------- | ------------------------------------------------------------------- |
| `-SaveBaseline`           | adopt this run as the new baseline                                   |
| `-Toolchain msvc`         | measure one configuration only                                       |
| `-CoreType P` / `E`       | measure one core type only; default `Auto` measures every type        |
| `-Rounds <n>`             | runs per configuration, default 5 (minimum 3)                        |
| `-Cpu <n>`                | pin to one exact logical CPU, overriding `-CoreType`                 |
| `-NoBuild`                | measure what is already built, skipping the build step               |
| `-ReportOnly`             | re-render the report from the JSON already in `bench/results/`       |

The script finds clang++ and ninja itself; neither needs to be on PATH.

## Where things go

Everything lives under `bench/results/`, which `.gitignore` excludes:

```
bench/results/baseline/bench_{msvc,clang}_{pcore,ecore}_release.json   the baseline
bench/results/current/bench_{msvc,clang}_{pcore,ecore}_release.json    the run just measured
bench/results/runs/{msvc,clang}_{pcore,ecore}_round<n>.json            all 5 runs behind it
bench/results/report.html                                              the report
```

The report carries, for each core type, one baseline-vs-current table per toolchain and one
**MSVC vs Clang** table for the current build, ordered by the size of the disagreement so the
benchmarks where the two code generators differ most come first.

Each report records the branch, commit and dirty state it was measured at, and the report shows both
sides, so a comparison can never be read without knowing which two states of the tree produced it.

## Establishing a baseline for a change

The baseline is only meaningful if it was measured on this machine, under this protocol. To measure
what a working tree change did:

1. `git stash push -- <the changed files>`
2. run with `-SaveBaseline`
3. `git stash pop`
4. run again

For a stronger comparison, build both variants first, copy the two `bench.exe` files aside, and run
them round-robin in one session. Interleaving spreads any machine drift evenly across both sides
instead of letting it accumulate on the second one measured.

## Reading the result

- Changes under 2% are reported as unchanged. That is this machine's run to run noise, not a result.
- The `noise` column is the spread across the three runs that survive the trim. A change that is not
  several times larger than it is unproven, whatever its sign.
- A benchmark that appears or disappears is shown as `new` or `removed`, never silently dropped.
- **Check that the change could have reached the benchmark before believing it.** Any edit shifts
  every function's address, and a hot loop landing on a different 32-byte phase can move 20% with
  byte-identical instructions. Measured on this repo: a one-line change inside
  `fixed_point128::operator*=` moved the `uint128_t` comparison benchmark by 22%, in a function whose
  machine code was unchanged - it had simply been pushed 96 bytes down the image. Confirm with an
  asm diff (`cl /FAsc` on both variants, compare PROC bodies) before attributing a number to a change.

## How a score is produced

Each configuration is run **5 times with no warm-up**. Per benchmark the fastest and the slowest
result are discarded and the score is the **mean of the remaining 3**; the spread across those same 3
is the reported noise.

The trim replaces the warm-up rather than merely tolerating its absence. Windows malware-scans a
freshly linked binary on its first execution and the scan lands inside the measurement - historically
worth over 100% on a single result - so run 1 is the slowest and is exactly what dropping the worst
removes. If you time something by hand instead, either discard the first run or trim the same way.

Nothing sleeps between runs: the benchmark is single threaded and will not push a desktop part into
thermal throttling, so a cool down would buy nothing but wall clock time.

## Rules that keep the numbers honest

These are not style preferences; each is a way this machine has produced a wrong answer before.
- **Never run unpinned.** The dev box is a hybrid i9-12900K (CPUs 0-7 P-cores, 8-15 E-cores). An
  unpinned thread migrates mid-run and measures the same binary 2.6x apart. Pinning takes the run to
  run spread from a 27% median to 0.4%.
- **Do not compare against a baseline from another machine or another protocol.** Rates are absolute,
  not normalized.
- When a result is surprising, check whether both configurations moved. A large change on one
  compiler only usually localizes a codegen difference; a change on both is usually real.
