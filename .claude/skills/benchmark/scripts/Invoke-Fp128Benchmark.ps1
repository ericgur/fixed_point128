<#
.SYNOPSIS
    Builds and runs the fixed_point128 release benchmarks, then compares them against a stored baseline.

.DESCRIPTION
    Drives the two release configurations this project measures on Windows - MSVC (Visual Studio
    generator) and Clang (Ninja Multi-Config, the command line clang++). clang-cl is deliberately not
    measured: it has never diverged from clang on this benchmark by more than the run to run noise.

    On a hybrid CPU every configuration is measured twice, once pinned to a performance core and once
    to an efficiency core, and the two are reported separately. They are different machines for this
    purpose - a P-core runs this benchmark about 2.6x faster than an E-core - and averaging them, or
    letting Windows choose between them, produces a number that describes neither.

    The rest of the protocol is not incidental either. Each part of it exists because of a specific way
    this machine has previously produced wrong numbers:

    - Every configuration is run five times and each benchmark's fastest and slowest result is thrown
      away; the score is the mean of the three that remain. Discarding the slowest is also what makes
      a separate warm-up unnecessary: the first execution of a freshly linked binary is malware
      scanned on Windows, and that run is the slow one the trim removes.
    - Every run is pinned. An unpinned thread migrates between core types mid-run and measures the
      same binary 2.6x apart; pinning takes the run to run spread from a 27% median to 0.4%.
    - Nothing waits between runs. The benchmark is single threaded, which will not drive a desktop
      part into thermal throttling, so a cool down would only add wall clock time.

.PARAMETER Rounds
    Runs per configuration. Default 5: the fastest and slowest are discarded and the score is the mean
    of the rest, so this must be at least 3 and is only meaningful from 5 upwards.

.PARAMETER CoreType
    Which core types to measure on. Auto (the default) measures every type the CPU has: both P and E
    on a hybrid part, a single set on a uniform one.

.PARAMETER Cpu
    Pin to this exact logical CPU and measure nothing else, overriding -CoreType. For investigating one
    particular core; the automatic choice is better for ordinary use.

.PARAMETER Toolchain
    Which configurations to measure. Default both.

.PARAMETER SaveBaseline
    Promote the run just measured to the baseline, replacing whatever was there.

.PARAMETER NoBuild
    Measure the binaries already in the build tree instead of rebuilding first.

.PARAMETER ReportOnly
    Skip building and measuring entirely and re-render the report from the JSON already in the results
    directory. Useful after editing the report, and for reporting on results collected some other way.

.PARAMETER Label
    Free text recorded with the run and shown in the report's build stamp. Use it when the two sides of
    a comparison come from the same commit - measuring a variant binary with -NoBuild, say - because
    the commit alone then cannot tell the reader which side is which.

.PARAMETER ResultsDir
    Where results and the report live. Default <repo>/bench/results, which is git ignored.

.EXAMPLE
    .\Invoke-Fp128Benchmark.ps1
    Build both configurations, measure them on every core type, and report against the baseline.

.EXAMPLE
    .\Invoke-Fp128Benchmark.ps1 -SaveBaseline
    Same, and then adopt this run as the new baseline.

.EXAMPLE
    .\Invoke-Fp128Benchmark.ps1 -CoreType P -Toolchain msvc -Rounds 3
    A quick single configuration check on a performance core, at the lowest useful run count.
#>
[CmdletBinding()]
param(
    [ValidateRange(3, 20)]
    [int]$Rounds = 5,

    [ValidateSet('Auto', 'P', 'E')]
    [string[]]$CoreType = @('Auto'),

    [ValidateRange(-1, 63)]
    [int]$Cpu = -1,

    [ValidateSet('msvc', 'clang')]
    [string[]]$Toolchain = @('msvc', 'clang'),

    [switch]$SaveBaseline,

    [switch]$NoBuild,

    [switch]$ReportOnly,

    [string]$Label,

    [string]$ResultsDir
)

$ErrorActionPreference = 'Stop'

# ---------------------------------------------------------------------------------------------
# Layout
# ---------------------------------------------------------------------------------------------

# scripts -> benchmark -> skills -> .claude -> repo root
$repoRoot = (Resolve-Path (Join-Path $PSScriptRoot '..\..\..\..')).Path

if (-not $ResultsDir) {
    $ResultsDir = Join-Path $repoRoot 'bench\results'
}
$baselineDir = Join-Path $ResultsDir 'baseline'
$currentDir  = Join-Path $ResultsDir 'current'
$runsDir     = Join-Path $ResultsDir 'runs'
$reportPath  = Join-Path $ResultsDir 'report.html'

foreach ($dir in @($ResultsDir, $baselineDir, $currentDir, $runsDir)) {
    if (-not (Test-Path $dir)) {
        New-Item -ItemType Directory -Path $dir -Force | Out-Null
    }
}

##
# @brief The measured configurations.
#
# ConfigurePreset and BuildPreset name entries in CMakePresets.json. ExeDir is where that preset's
# generator drops the release benchmark, which differs between the Visual Studio and the Ninja
# Multi-Config generator.
##
$configs = [ordered]@{
    msvc = @{
        Name             = 'msvc'
        DisplayName      = 'MSVC'
        ConfigurePreset  = 'msvc'
        BuildPreset      = 'msvc-release'
        ExeDir           = Join-Path $repoRoot 'out\build\msvc\bin\Release'
        NeedsClangOnPath = $false
    }
    clang = @{
        Name             = 'clang'
        DisplayName      = 'Clang'
        ConfigurePreset  = 'clang'
        BuildPreset      = 'clang-release'
        ExeDir           = Join-Path $repoRoot 'out\build\clang\bin\Release'
        NeedsClangOnPath = $true
    }
}

# Rates that move by less than this are reported as unchanged. The pinned, best-of-N protocol has a
# run to run spread of roughly 0.4% median and 1.4% at p90, so 2% sits comfortably outside it.
#
# Run to run noise is not the only reason a rate moves without the code changing. Any edit shifts every
# function's address, and a hot loop that lands on a different alignment can move by 20% or more with
# byte-identical instructions. The report cannot detect that; when a benchmark moves in code the change
# could not have reached, suspect layout before believing the number.
$noiseThresholdPercent = 2.0

# ---------------------------------------------------------------------------------------------
# CPU topology
# ---------------------------------------------------------------------------------------------

##
# @brief Maps each logical CPU to its efficiency class, via GetLogicalProcessorInformationEx.
#
# Windows records a per core "efficiency class" - higher is faster - which is how a hybrid part
# distinguishes performance cores from efficiency cores. Neither WMI nor the .NET framework exposes it,
# so this reads the native structure directly.
#
# The layout being walked, on x64:
#   SYSTEM_LOGICAL_PROCESSOR_INFORMATION_EX { DWORD Relationship; DWORD Size; union { ... } }
#   PROCESSOR_RELATIONSHIP { BYTE Flags; BYTE EfficiencyClass; BYTE Reserved[20]; WORD GroupCount;
#                            GROUP_AFFINITY GroupMask[]; }
# which puts EfficiencyClass at +9, GroupCount at +30, and the first 16 byte GROUP_AFFINITY at +32.
#
# @return Hashtable of logical CPU number to efficiency class, or $null if it could not be read.
##
function Get-CpuEfficiencyClasses {
    $source = @'
using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Runtime.InteropServices;

public static class Fp128Topology
{
    [DllImport("kernel32.dll", SetLastError = true)]
    private static extern bool GetLogicalProcessorInformationEx(int relationship, byte[] buffer, ref int length);

    private const int RelationProcessorCore = 0;

    public static Dictionary<int, int> EfficiencyClasses()
    {
        int length = 0;
        GetLogicalProcessorInformationEx(RelationProcessorCore, null, ref length);
        byte[] buffer = new byte[length];
        if (!GetLogicalProcessorInformationEx(RelationProcessorCore, buffer, ref length))
            throw new Win32Exception(Marshal.GetLastWin32Error());

        var classes = new Dictionary<int, int>();
        int offset = 0;
        while (offset + 8 <= length)
        {
            int size = BitConverter.ToInt32(buffer, offset + 4);
            if (size <= 0) break;
            int efficiencyClass = buffer[offset + 9];
            int groupCount = BitConverter.ToUInt16(buffer, offset + 30);
            for (int g = 0; g < groupCount; g++)
            {
                int groupOffset = offset + 32 + g * 16;
                ulong mask = BitConverter.ToUInt64(buffer, groupOffset);
                int group = BitConverter.ToUInt16(buffer, groupOffset + 8);
                for (int bit = 0; bit < 64; bit++)
                    if ((mask & (1UL << bit)) != 0)
                        classes[group * 64 + bit] = efficiencyClass;
            }
            offset += size;
        }
        return classes;
    }
}
'@

    try {
        if (-not ('Fp128Topology' -as [type])) {
            Add-Type -TypeDefinition $source -ErrorAction Stop
        }
        $classes = [Fp128Topology]::EfficiencyClasses()
        if ($classes.Count -eq 0) {
            return $null
        }

        $result = @{}
        foreach ($key in $classes.Keys) {
            $result[[int]$key] = [int]$classes[$key]
        }
        return $result
    } catch {
        Write-Warning ('Could not read the CPU topology ({0}); treating the CPU as uniform.' -f $_.Exception.Message)
        return $null
    }
}

##
# @brief Chooses which logical CPUs to measure on.
#
# On a hybrid CPU the highest efficiency class is the performance cores and the lowest is the efficiency
# cores. CPU 0 is avoided wherever there is an alternative: it services interrupts, and on this class of
# machine that alone is worth a few percent.
#
# @return Target descriptors, each with Label (used in file names), Description, Cpu and Class.
##
function Get-CoreTargets {
    if ($Cpu -ge 0) {
        return @(@{ Label = "cpu$Cpu"; Description = "CPU $Cpu"; Cpu = $Cpu; Class = 'explicit' })
    }

    $classes = Get-CpuEfficiencyClasses
    if (-not $classes) {
        return @(@{ Label = 'cpu'; Description = 'default CPU'; Cpu = 0; Class = 'uniform' })
    }

    $distinct = @($classes.Values | Sort-Object -Unique)
    $cpus = @($classes.Keys | Sort-Object)

    if ($distinct.Count -eq 1) {
        $chosen = @($cpus | Where-Object { $_ -ne 0 })
        if ($chosen.Count -eq 0) {
            $chosen = $cpus
        }
        return @(@{ Label = 'cpu'; Description = "CPU $($chosen[0]) (uniform CPU)"; Cpu = $chosen[0]; Class = 'uniform' })
    }

    $performanceCpus = @($cpus | Where-Object { $classes[$_] -eq $distinct[-1] -and $_ -ne 0 })
    $efficiencyCpus  = @($cpus | Where-Object { $classes[$_] -eq $distinct[0]  -and $_ -ne 0 })

    $wanted = $CoreType
    if ($wanted -contains 'Auto') {
        $wanted = @('P', 'E')
    }

    $targets = @()
    if ($wanted -contains 'P' -and $performanceCpus.Count -gt 0) {
        # CPU 2 wherever it qualifies, so results stay comparable with every baseline taken before the
        # script picked the core itself.
        $pick = $performanceCpus[0]
        if ($performanceCpus -contains 2) {
            $pick = 2
        }
        $targets += @{ Label = 'pcore'; Description = "P-core (CPU $pick)"; Cpu = $pick; Class = 'P' }
    }
    if ($wanted -contains 'E' -and $efficiencyCpus.Count -gt 0) {
        $targets += @{ Label = 'ecore'; Description = "E-core (CPU $($efficiencyCpus[0]))"; Cpu = $efficiencyCpus[0]; Class = 'E' }
    }

    if ($targets.Count -eq 0) {
        throw "No logical CPU matches -CoreType $($CoreType -join ',')."
    }
    return $targets
}

# ---------------------------------------------------------------------------------------------
# Toolchain discovery
# ---------------------------------------------------------------------------------------------

##
# @brief Returns the directory holding an executable, searching PATH and then known install roots.
#
# The Clang preset drives clang++ and Ninja directly and neither is on PATH on a stock Visual Studio
# machine - both ship inside other products. Wildcards in @p Candidates are expanded, highest match
# first, so a versioned directory such as llvm-mingw1706_64 does not have to be named exactly.
#
# @param Name       Executable to find, without the .exe suffix.
# @param Candidates Directories to search when the executable is not already on PATH.
# @return The containing directory, or $null when the executable could not be found.
##
function Find-ToolDirectory {
    param(
        [Parameter(Mandatory)][string]$Name,
        [Parameter(Mandatory)][string[]]$Candidates
    )

    $onPath = Get-Command $Name -ErrorAction SilentlyContinue
    if ($onPath) {
        return (Split-Path -Parent $onPath.Source)
    }

    foreach ($candidate in $Candidates) {
        $found = @(Get-Item (Join-Path $candidate "$Name.exe") -ErrorAction SilentlyContinue |
                   Sort-Object FullName -Descending)
        if ($found.Count -gt 0) {
            return (Split-Path -Parent $found[0].FullName)
        }
    }

    return $null
}

##
# @brief Prepends clang++ and ninja to PATH for this process, so the Clang preset can be configured.
#
# The Visual Studio bundled LLVM is searched last on purpose: its clang++ targets the MSVC ABI, which
# would turn the "clang" configuration into a second clang-cl rather than the independent code
# generator it is measured for.
#
# @return $true when both tools were found.
##
function Enable-ClangToolchain {
    $clangDir = Find-ToolDirectory -Name 'clang++' -Candidates @(
        'C:\Qt\Tools\llvm-mingw*\bin',
        'C:\Program Files\LLVM\bin',
        'C:\Program Files\Microsoft Visual Studio\*\*\VC\Tools\Llvm\x64\bin')

    $ninjaDir = Find-ToolDirectory -Name 'ninja' -Candidates @(
        'C:\Qt\Tools\Ninja',
        'C:\Program Files\Microsoft Visual Studio\*\*\Common7\IDE\CommonExtensions\Microsoft\CMake\Ninja')

    if (-not $clangDir) {
        Write-Warning 'clang++.exe not found; skipping the clang configuration.'
        return $false
    }
    if (-not $ninjaDir) {
        Write-Warning 'ninja.exe not found; skipping the clang configuration.'
        return $false
    }

    foreach ($dir in @($clangDir, $ninjaDir)) {
        if (($env:PATH -split ';') -notcontains $dir) {
            $env:PATH = "$dir;$env:PATH"
        }
    }

    Write-Host ('  clang++ : {0}' -f (Join-Path $clangDir 'clang++.exe'))
    Write-Host ('  ninja   : {0}' -f (Join-Path $ninjaDir 'ninja.exe'))
    return $true
}

# ---------------------------------------------------------------------------------------------
# Build
# ---------------------------------------------------------------------------------------------

##
# @brief Configures a preset if its build tree does not exist yet, then builds the bench target.
#
# @param Config One entry of $configs.
##
function Build-BenchTarget {
    param([Parameter(Mandatory)][hashtable]$Config)

    $buildTree = Join-Path $repoRoot ('out\build\{0}' -f $Config.ConfigurePreset)
    if (-not (Test-Path (Join-Path $buildTree 'CMakeCache.txt'))) {
        Write-Host ("Configuring preset '{0}'..." -f $Config.ConfigurePreset)
        & cmake --preset $Config.ConfigurePreset | Out-Null
        if ($LASTEXITCODE -ne 0) {
            throw "cmake --preset $($Config.ConfigurePreset) failed"
        }
    }

    Write-Host ("Building preset '{0}'..." -f $Config.BuildPreset)
    # Merging stderr into the output stream would otherwise raise a terminating NativeCommandError
    # under ErrorActionPreference = Stop, and MSBuild writes to stderr on builds that succeed.
    $previousPreference = $ErrorActionPreference
    $ErrorActionPreference = 'Continue'
    $output = & cmake --build --preset $Config.BuildPreset --target bench 2>&1
    $ErrorActionPreference = $previousPreference
    if ($LASTEXITCODE -ne 0) {
        $output | ForEach-Object { Write-Host $_ }
        throw "cmake --build --preset $($Config.BuildPreset) failed"
    }
}

# ---------------------------------------------------------------------------------------------
# Measurement
# ---------------------------------------------------------------------------------------------

##
# @brief Runs bench.exe once, pinned to a single logical CPU at High priority.
#
# The affinity mask and the priority are applied to the live process, which is the only way to set them
# from PowerShell. Both land long before the process finishes sizing its first batch, so no measured
# work escapes them.
#
# Reading the process Handle is what keeps ExitCode readable: without it the object loses its native
# handle when the process ends and reports a null exit code, which would silently pass off a crashed run
# as a good one.
#
# @param ExeDir    Directory holding bench.exe; also the working directory, so any JSON report the run
#                  writes lands there.
# @param BenchArgs Command line for bench.exe.
# @param TargetCpu Logical CPU to pin to.
##
function Invoke-PinnedBench {
    param(
        [Parameter(Mandatory)][string]$ExeDir,
        [Parameter(Mandatory)][string[]]$BenchArgs,
        [Parameter(Mandatory)][int]$TargetCpu
    )

    $exe = Join-Path $ExeDir 'bench.exe'
    if (-not (Test-Path $exe)) {
        throw "Benchmark binary not found: $exe"
    }

    $stdout = Join-Path $ExeDir 'bench_stdout.txt'
    $process = Start-Process -FilePath $exe -ArgumentList $BenchArgs -WorkingDirectory $ExeDir `
                             -PassThru -NoNewWindow -RedirectStandardOutput $stdout
    $null = $process.Handle
    $process.ProcessorAffinity = [IntPtr]([int64]1 -shl $TargetCpu)
    $process.PriorityClass = 'High'
    $process.WaitForExit()

    if ($process.ExitCode -ne 0) {
        throw "bench.exe exited with $($process.ExitCode); see $stdout"
    }
}

##
# @brief Returns the commit the working tree is at, and whether it has uncommitted changes.
#
# Recorded in every report so that a comparison can never be read without knowing which two states of
# the tree produced it.
##
function Get-GitDescription {
    $commit = & git -C $repoRoot rev-parse --short HEAD 2>$null
    $branch = & git -C $repoRoot rev-parse --abbrev-ref HEAD 2>$null
    $status = & git -C $repoRoot status --porcelain 2>$null

    $description = [ordered]@{ commit = 'unknown'; branch = 'unknown'; dirty = [bool]$status }
    if ($commit) { $description.commit = "$commit".Trim() }
    if ($branch) { $description.branch = "$branch".Trim() }
    return $description
}

##
# @brief Measures one configuration on one core and returns the fastest result per benchmark.
#
# Runs the binary $Rounds times with nothing discarded up front, then scores each benchmark as the
# mean of its results with the fastest and the slowest removed. The trim is doing two jobs: it drops
# the one run that interference made slow, and it drops the malware scan that lands inside the first
# execution of a freshly linked binary, which is why there is no separate warm-up.
#
# The spread across the runs that survive the trim is kept alongside the score. That is the run to run
# variation of the number being reported, and it is what tells a reader of the report whether a
# difference is larger than this machine's noise.
#
# @param Config One entry of $configs.
# @param Target One entry from Get-CoreTargets.
# @return A report object ready to be serialized.
##
function Measure-Configuration {
    param(
        [Parameter(Mandatory)][hashtable]$Config,
        [Parameter(Mandatory)][hashtable]$Target
    )

    Write-Host ('Measuring {0} on {1}, {2} runs...' -f $Config.DisplayName, $Target.Description, $Rounds)

    $runs = @()
    for ($round = 1; $round -le $Rounds; $round++) {
        $stopwatch = [Diagnostics.Stopwatch]::StartNew()
        Invoke-PinnedBench -ExeDir $Config.ExeDir -BenchArgs @('-j') -TargetCpu $Target.Cpu

        $produced = @(Get-ChildItem (Join-Path $Config.ExeDir 'bench_*_release.json'))
        if ($produced.Count -eq 0) {
            throw "The run produced no JSON report in $($Config.ExeDir)"
        }

        $saved = Join-Path $runsDir ('{0}_{1}_round{2}.json' -f $Config.Name, $Target.Label, $round)
        Move-Item $produced[0].FullName $saved -Force
        $runs += (Get-Content $saved -Raw | ConvertFrom-Json)

        Write-Host ('  round {0}/{1}  {2,5:N1}s' -f $round, $Rounds, $stopwatch.Elapsed.TotalSeconds)
    }

    # Key on type/group/name rather than on position: a benchmark added or removed between two builds
    # would otherwise silently shift every row after it.
    $byKey = [ordered]@{}
    foreach ($run in $runs) {
        foreach ($result in $run.results) {
            $key = '{0}|{1}|{2}' -f $result.type, $result.group, $result.name
            if (-not $byKey.Contains($key)) {
                $byKey[$key] = [ordered]@{ type = $result.type; group = $result.group; name = $result.name; rates = @() }
            }
            $byKey[$key].rates += [double]$result.iterationsPerSecond
        }
    }

    $results = @()
    foreach ($entry in $byKey.Values) {
        # Sort, drop one from each end, and average what is left. With five runs that keeps three:
        # the slowest is whatever the machine interfered with, the fastest is the one measurement most
        # likely to have caught an unrepresentatively quiet moment.
        $sorted = @($entry.rates | Sort-Object)
        $kept = @($sorted[1..($sorted.Count - 2)])

        $score = ($kept | Measure-Object -Average).Average
        $spread = 0.0
        if ($score -gt 0) {
            $high = ($kept | Measure-Object -Maximum).Maximum
            $low = ($kept | Measure-Object -Minimum).Minimum
            $spread = 100.0 * ($high - $low) / $score
        }

        $results += [ordered]@{
            type                = $entry.type
            group               = $entry.group
            name                = $entry.name
            iterationsPerSecond = $score
            rates               = $entry.rates
            keptRates           = $kept
            spreadPercent       = $spread
        }
    }

    $first = $runs[0]
    $libraryVersion = 'unknown'
    if ($first.PSObject.Properties.Name -contains 'libraryVersion') {
        $libraryVersion = $first.libraryVersion
    }

    return [ordered]@{
        toolchain       = $Config.Name
        toolchainName   = $Config.DisplayName
        coreLabel       = $Target.Label
        coreDescription = $Target.Description
        coreClass       = $Target.Class
        compiler        = $first.compiler
        compilerTag     = $first.compilerTag
        libraryVersion  = $libraryVersion
        build           = $first.build
        timestamp       = $first.timestamp
        rounds          = $Rounds
        cpu             = $Target.Cpu
        sampling        = "mean of $($Rounds - 2) of $Rounds pinned runs, fastest and slowest discarded"
        label           = $Label
        git             = Get-GitDescription
        results         = $results
    }
}

# ---------------------------------------------------------------------------------------------
# Comparison
# ---------------------------------------------------------------------------------------------

##
# @brief Reads the spread field from a result, tolerating reports written before it existed.
##
function Get-SpreadPercent {
    param([Parameter(Mandatory)]$Result)

    if ($Result.PSObject.Properties.Name -contains 'spreadPercent') {
        return [double]$Result.spreadPercent
    }
    return $null
}

##
# @brief Indexes a report's results by type, group and name.
##
function Get-ResultIndex {
    param([Parameter(Mandatory)]$Report)

    $index = @{}
    foreach ($result in $Report.results) {
        $index['{0}|{1}|{2}' -f $result.type, $result.group, $result.name] = $result
    }
    return $index
}

##
# @brief Joins a baseline report and a current report into per benchmark comparison rows.
#
# Benchmarks present on only one side are kept, with a null rate on the missing side, so that adding or
# removing a benchmark shows up in the report instead of vanishing from it.
#
# @return Comparison rows in baseline order, with current-only rows appended.
##
function Compare-Reports {
    param(
        [Parameter(Mandatory)]$Baseline,
        [Parameter(Mandatory)]$Current
    )

    $currentByKey = Get-ResultIndex -Report $Current
    $rows = @()
    $seen = @{}

    foreach ($baseResult in $Baseline.results) {
        $key = '{0}|{1}|{2}' -f $baseResult.type, $baseResult.group, $baseResult.name
        $seen[$key] = $true

        $baseRate = [double]$baseResult.iterationsPerSecond
        $currentRate = $null
        $currentNoise = $null
        $delta = $null
        if ($currentByKey.ContainsKey($key)) {
            $currentResult = $currentByKey[$key]
            $currentRate = [double]$currentResult.iterationsPerSecond
            $currentNoise = Get-SpreadPercent -Result $currentResult
            if ($baseRate -gt 0) {
                $delta = 100.0 * ($currentRate - $baseRate) / $baseRate
            }
        }

        $rows += [pscustomobject]@{
            Type          = $baseResult.type
            Group         = $baseResult.group
            Name          = $baseResult.name
            BaselineRate  = $baseRate
            CurrentRate   = $currentRate
            BaselineNoise = Get-SpreadPercent -Result $baseResult
            CurrentNoise  = $currentNoise
            DeltaPercent  = $delta
        }
    }

    foreach ($result in $Current.results) {
        $key = '{0}|{1}|{2}' -f $result.type, $result.group, $result.name
        if ($seen.ContainsKey($key)) {
            continue
        }
        $rows += [pscustomobject]@{
            Type          = $result.type
            Group         = $result.group
            Name          = $result.name
            BaselineRate  = $null
            CurrentRate   = [double]$result.iterationsPerSecond
            BaselineNoise = $null
            CurrentNoise  = Get-SpreadPercent -Result $result
            DeltaPercent  = $null
        }
    }

    return $rows
}

##
# @brief Compares the current build under two toolchains, benchmark by benchmark.
#
# The sign convention is "@p Right relative to @p Left": positive means the right hand toolchain
# produced the faster code. Rows come back ordered by the size of the disagreement regardless of its
# sign, because where the two code generators differ most is the only reason to read this table.
#
# @return Comparison rows, largest absolute difference first.
##
function Compare-Toolchains {
    param(
        [Parameter(Mandatory)]$Left,
        [Parameter(Mandatory)]$Right
    )

    $rightByKey = Get-ResultIndex -Report $Right
    $rows = @()

    foreach ($leftResult in $Left.results) {
        $key = '{0}|{1}|{2}' -f $leftResult.type, $leftResult.group, $leftResult.name
        if (-not $rightByKey.ContainsKey($key)) {
            continue
        }

        $leftRate = [double]$leftResult.iterationsPerSecond
        $rightRate = [double]$rightByKey[$key].iterationsPerSecond
        if ($leftRate -le 0 -or $rightRate -le 0) {
            continue
        }

        $rows += [pscustomobject]@{
            Type         = $leftResult.type
            Group        = $leftResult.group
            Name         = $leftResult.name
            LeftRate     = $leftRate
            RightRate    = $rightRate
            DeltaPercent = 100.0 * ($rightRate - $leftRate) / $leftRate
        }
    }

    return @($rows | Sort-Object { - [Math]::Abs($_.DeltaPercent) })
}

# ---------------------------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------------------------

##
# @brief Escapes text for inclusion in HTML.
#
# Benchmark names carry both angle brackets and ampersands - "fixed_point128<10>" and
# "Operators >, >=, <, <=" - so this is load bearing, not defensive.
##
function ConvertTo-HtmlText {
    param([Parameter(Mandatory)][AllowEmptyString()][string]$Text)

    return $Text.Replace('&', '&amp;').Replace('<', '&lt;').Replace('>', '&gt;').Replace('"', '&quot;')
}

##
# @brief Formats a number with the invariant culture.
#
# Used for every number that ends up inside CSS or in a data attribute the page's script parses. A
# machine on a comma-decimal locale would otherwise emit "width:1,5%", which no browser accepts.
##
function Format-Invariant {
    param(
        [Parameter(Mandatory)][double]$Value,
        [Parameter(Mandatory)][string]$Format
    )

    return $Value.ToString($Format, [System.Globalization.CultureInfo]::InvariantCulture)
}

##
# @brief Formats a rate in iterations per second as millions per second.
##
function Format-Rate {
    param([Parameter(Mandatory)][double]$Rate)

    return '{0:N1}' -f ($Rate / 1e6)
}

##
# @brief Renders a signed percentage with an explicit plus sign.
##
function Format-Delta {
    param([Parameter(Mandatory)][double]$Value)

    $sign = ''
    if ($Value -ge 0) {
        $sign = '+'
    }
    return '{0}{1:N1}%' -f $sign, $Value
}

##
# @brief Renders the "branch @ commit - version" stamp identifying one side of a comparison.
##
function Format-BuildStamp {
    param([Parameter(Mandatory)]$Report)

    $stamp = 'unknown commit'
    if ($Report.PSObject.Properties.Name -contains 'git') {
        $stamp = '{0} @ {1}' -f $Report.git.branch, $Report.git.commit
        if ($Report.git.dirty) {
            $stamp += ' + uncommitted changes'
        }
    }
    if ($Report.PSObject.Properties.Name -contains 'label' -and $Report.label) {
        $stamp = '{0} - {1}' -f $Report.label, $stamp
    }
    if ($Report.PSObject.Properties.Name -contains 'libraryVersion') {
        $stamp += ' - ' + $Report.libraryVersion
    }
    return $stamp
}

##
# @brief Trims a compiler's self-reported name down to something a heading can hold.
#
# Clang reports itself as "Clang 17.0.6 (https://github.com/llvm/llvm-project.git 6009708b...)"; the
# provenance belongs in the JSON, not across the top of a table.
##
function Format-CompilerName {
    param([Parameter(Mandatory)][string]$Name)

    $cut = $Name.IndexOf(' (')
    if ($cut -gt 0) {
        return $Name.Substring(0, $cut)
    }
    return $Name
}

##
# @brief Renders the two sided bar that visualises a signed percentage.
#
# The bar occupies one half of the track, so a 100% change fills it; anything larger clamps.
##
function Format-Bar {
    param(
        [Parameter(Mandatory)][double]$Value,
        [Parameter(Mandatory)][string]$Class
    )

    $direction = 'right'
    if ($Value -lt 0) {
        $direction = 'left'
    }
    $width = [Math]::Min(100.0, [Math]::Abs($Value)) / 2.0
    return '<span class="track {0}"><span class="bar {1}" style="width:{2}%"></span></span>' -f `
           $direction, $Class, (Format-Invariant -Value $width -Format '0.0')
}

##
# @brief Renders one summary card.
##
function New-SummaryCard {
    param(
        [Parameter(Mandatory)][string]$Label,
        [Parameter(Mandatory)]$Value,
        [Parameter(Mandatory)][string]$Class,
        [string]$Note = ''
    )

    $noteHtml = ''
    if ($Note) {
        $noteHtml = '<span class="note">{0}</span>' -f (ConvertTo-HtmlText $Note)
    }
    return '<div class="card"><span class="value {0}">{1}</span><span class="label">{2}</span>{3}</div>' -f `
           $Class, (ConvertTo-HtmlText ([string]$Value)), (ConvertTo-HtmlText $Label), $noteHtml
}

##
# @brief Renders one baseline versus current row.
##
function New-ResultRow {
    param([Parameter(Mandatory)]$Row)

    $delta = $Row.DeltaPercent
    if ($null -eq $delta) {
        $class = 'missing'
        $deltaText = 'new'
        if ($null -eq $Row.CurrentRate) {
            $deltaText = 'removed'
        }
        $deltaAttribute = '0'
        $bar = '<span class="track right"></span>'
    } else {
        $class = 'flat'
        if ($delta -gt $noiseThresholdPercent) {
            $class = 'good'
        } elseif ($delta -lt -$noiseThresholdPercent) {
            $class = 'bad'
        }
        $deltaText = Format-Delta -Value $delta
        $deltaAttribute = Format-Invariant -Value $delta -Format '0.0000'
        $bar = Format-Bar -Value $delta -Class $class
    }

    $noise = ''
    if ($null -ne $Row.CurrentNoise) {
        $noise = '{0:N1}%' -f $Row.CurrentNoise
    }

    $baselineText = '-'
    if ($null -ne $Row.BaselineRate) {
        $baselineText = Format-Rate -Rate $Row.BaselineRate
    }
    $currentText = '-'
    if ($null -ne $Row.CurrentRate) {
        $currentText = Format-Rate -Rate $Row.CurrentRate
    }

    return ('<tr class="{0}" data-delta="{1}"><td class="type">{2}</td><td>{3}</td>' +
            '<td class="num">{4}</td><td class="num">{5}</td><td class="num noise">{6}</td>' +
            '<td class="num delta {0}">{7}</td><td class="barcol">{8}</td></tr>') -f `
           $class, $deltaAttribute, (ConvertTo-HtmlText $Row.Type), (ConvertTo-HtmlText $Row.Name),
           $baselineText, $currentText, $noise, $deltaText, $bar
}

##
# @brief Renders one toolchain versus toolchain row.
##
function New-ToolchainRow {
    param(
        [Parameter(Mandatory)]$Row,
        [Parameter(Mandatory)][string]$LeftName,
        [Parameter(Mandatory)][string]$RightName
    )

    $delta = $Row.DeltaPercent
    $class = 'right-ahead'
    $winner = $RightName
    if ($delta -lt 0) {
        $class = 'left-ahead'
        $winner = $LeftName
    }
    if ([Math]::Abs($delta) -le $noiseThresholdPercent) {
        $class = 'flat'
        $winner = 'tie'
    }

    return ('<tr class="{0}" data-delta="{1}"><td class="type">{2}</td><td>{3}</td>' +
            '<td class="num">{4}</td><td class="num">{5}</td><td class="winner {0}">{6}</td>' +
            '<td class="num delta {0}">{7}</td><td class="barcol">{8}</td></tr>') -f `
           $class,
           (Format-Invariant -Value $delta -Format '0.0000'),
           (ConvertTo-HtmlText $Row.Type),
           (ConvertTo-HtmlText $Row.Name),
           (Format-Rate -Rate $Row.LeftRate),
           (Format-Rate -Rate $Row.RightRate),
           (ConvertTo-HtmlText $winner),
           (Format-Delta -Value $delta),
           (Format-Bar -Value $delta -Class $class)
}

##
# @brief Renders one configuration on one core: build metadata, a summary, and the benchmark table.
##
function New-ComparisonSection {
    param([Parameter(Mandatory)][hashtable]$Comparison)

    $baseline = $Comparison.Baseline
    $current  = $Comparison.Current
    $rows     = @($Comparison.Rows)

    $compared = @($rows | Where-Object { $null -ne $_.DeltaPercent })
    $faster   = @($compared | Where-Object { $_.DeltaPercent -gt $noiseThresholdPercent })
    $slower   = @($compared | Where-Object { $_.DeltaPercent -lt -$noiseThresholdPercent })

    $section = New-Object System.Text.StringBuilder
    [void]$section.AppendLine(('<h3>{0}</h3>' -f (ConvertTo-HtmlText (Format-CompilerName -Name $current.compiler))))

    [void]$section.AppendLine('<table class="meta"><tbody>')
    [void]$section.AppendLine(('<tr><th>baseline</th><td>{0}</td><td>{1}</td></tr>' -f `
                               (ConvertTo-HtmlText (Format-BuildStamp -Report $baseline)),
                               (ConvertTo-HtmlText $baseline.timestamp)))
    [void]$section.AppendLine(('<tr><th>current</th><td>{0}</td><td>{1}</td></tr>' -f `
                               (ConvertTo-HtmlText (Format-BuildStamp -Report $current)),
                               (ConvertTo-HtmlText $current.timestamp)))
    [void]$section.AppendLine('</tbody></table>')

    [void]$section.AppendLine('<div class="cards">')
    [void]$section.AppendLine((New-SummaryCard -Label 'faster' -Value $faster.Count -Class 'good'))
    [void]$section.AppendLine((New-SummaryCard -Label 'slower' -Value $slower.Count -Class 'bad'))
    [void]$section.AppendLine((New-SummaryCard -Label 'unchanged' -Value ($compared.Count - $faster.Count - $slower.Count) -Class 'flat'))
    if ($faster.Count -gt 0) {
        $best = @($faster | Sort-Object DeltaPercent -Descending)[0]
        [void]$section.AppendLine((New-SummaryCard -Label 'largest gain' -Value ('{0:N0}%' -f $best.DeltaPercent) -Class 'good' -Note $best.Name))
    }
    if ($slower.Count -gt 0) {
        $worst = @($slower | Sort-Object DeltaPercent)[0]
        [void]$section.AppendLine((New-SummaryCard -Label 'largest loss' -Value ('{0:N0}%' -f $worst.DeltaPercent) -Class 'bad' -Note $worst.Name))
    }
    [void]$section.AppendLine('</div>')

    [void]$section.AppendLine('<div class="tablewrap"><table class="results"><thead><tr>')
    [void]$section.AppendLine('<th>type</th><th>benchmark</th><th class="num">baseline M/s</th><th class="num">current M/s</th>')
    [void]$section.AppendLine('<th class="num">noise</th><th class="num">change</th><th class="barcol"></th>')
    [void]$section.AppendLine('</tr></thead><tbody>')
    foreach ($row in $rows) {
        [void]$section.AppendLine((New-ResultRow -Row $row))
    }
    [void]$section.AppendLine('</tbody></table></div>')

    return $section.ToString()
}

##
# @brief Renders the toolchain against toolchain table for one core.
##
function New-ToolchainDiffSection {
    param([Parameter(Mandatory)][hashtable]$Diff)

    $rows = @($Diff.Rows)
    $leftName = $Diff.LeftName
    $rightName = $Diff.RightName

    $rightAhead = @($rows | Where-Object { $_.DeltaPercent -gt $noiseThresholdPercent })
    $leftAhead  = @($rows | Where-Object { $_.DeltaPercent -lt -$noiseThresholdPercent })

    $section = New-Object System.Text.StringBuilder
    [void]$section.AppendLine(('<h3>{0} vs {1}</h3>' -f (ConvertTo-HtmlText $leftName), (ConvertTo-HtmlText $rightName)))
    [void]$section.AppendLine(('<p class="sub">Current build under both toolchains, largest disagreement first. ' +
                               'A positive difference means {0} produced the faster code.</p>') -f `
                              (ConvertTo-HtmlText $rightName))

    [void]$section.AppendLine('<div class="cards">')
    [void]$section.AppendLine((New-SummaryCard -Label ('{0} ahead' -f $rightName) -Value $rightAhead.Count -Class 'right-ahead'))
    [void]$section.AppendLine((New-SummaryCard -Label ('{0} ahead' -f $leftName) -Value $leftAhead.Count -Class 'left-ahead'))
    [void]$section.AppendLine((New-SummaryCard -Label 'within noise' -Value ($rows.Count - $rightAhead.Count - $leftAhead.Count) -Class 'flat'))
    if ($rows.Count -gt 0) {
        $widest = $rows[0]
        $widestClass = 'right-ahead'
        if ($widest.DeltaPercent -lt 0) {
            $widestClass = 'left-ahead'
        }
        [void]$section.AppendLine((New-SummaryCard -Label 'widest gap' -Value (Format-Delta -Value $widest.DeltaPercent) `
                                                  -Class $widestClass -Note ('{0} - {1}' -f $widest.Type, $widest.Name)))
    }
    [void]$section.AppendLine('</div>')

    [void]$section.AppendLine('<div class="tablewrap"><table class="results toolchains"><thead><tr>')
    [void]$section.AppendLine(('<th>type</th><th>benchmark</th><th class="num">{0} M/s</th><th class="num">{1} M/s</th>' -f `
                               (ConvertTo-HtmlText $leftName), (ConvertTo-HtmlText $rightName)))
    [void]$section.AppendLine('<th>faster</th><th class="num">difference</th><th class="barcol"></th>')
    [void]$section.AppendLine('</tr></thead><tbody>')
    foreach ($row in $rows) {
        [void]$section.AppendLine((New-ToolchainRow -Row $row -LeftName $leftName -RightName $rightName))
    }
    [void]$section.AppendLine('</tbody></table></div>')

    return $section.ToString()
}

##
# @brief The report stylesheet. Self contained, and follows the browser's light/dark preference.
##
function Get-ReportStyle {
    return @'
<style>
:root {
    --bg: #ffffff; --fg: #1a1a1a; --muted: #6b7280; --line: #e5e7eb; --panel: #f9fafb;
    --good: #15803d; --bad: #b91c1c; --flat: #6b7280;
    --good-bar: #86efac; --bad-bar: #fca5a5;
    --left: #b45309; --right: #1d4ed8; --left-bar: #fcd34d; --right-bar: #93c5fd;
}
@media (prefers-color-scheme: dark) {
    :root {
        --bg: #16181d; --fg: #e8e8ea; --muted: #9ba1ac; --line: #2c3038; --panel: #1d2027;
        --good: #4ade80; --bad: #f87171; --flat: #9ba1ac;
        --good-bar: #166534; --bad-bar: #7f1d1d;
        --left: #fbbf24; --right: #60a5fa; --left-bar: #78350f; --right-bar: #1e3a8a;
    }
}
* { box-sizing: border-box; }
body { margin: 0; background: var(--bg); color: var(--fg);
       font: 15px/1.5 -apple-system, "Segoe UI", Roboto, Helvetica, Arial, sans-serif; }
main { max-width: 1120px; margin: 0 auto; padding: 32px 20px 64px; }
h1 { font-size: 26px; margin: 0 0 4px; letter-spacing: -0.01em; }
h2 { font-size: 21px; margin: 46px 0 4px; padding-top: 26px; border-top: 2px solid var(--line); }
h3 { font-size: 16px; margin: 32px 0 10px; }
section:first-of-type h2 { border-top: none; padding-top: 0; margin-top: 28px; }
.sub { color: var(--muted); font-size: 13px; margin: 0 0 8px; }
.legend { margin-top: 32px; }
table { border-collapse: collapse; width: 100%; }
table.meta { width: auto; margin: 0 0 16px; font-size: 13px; color: var(--muted); }
table.meta th { text-align: left; font-weight: 600; padding: 2px 16px 2px 0; color: var(--fg); }
table.meta td { padding: 2px 16px 2px 0; }
.cards { display: flex; flex-wrap: wrap; gap: 10px; margin: 0 0 16px; }
.card { background: var(--panel); border: 1px solid var(--line); border-radius: 8px;
        padding: 10px 14px; min-width: 108px; display: flex; flex-direction: column; }
.card .value { font-size: 21px; font-weight: 650; font-variant-numeric: tabular-nums; }
.card .label { font-size: 11px; text-transform: uppercase; letter-spacing: 0.05em; color: var(--muted); }
.card .note { font-size: 11px; color: var(--muted); margin-top: 4px; max-width: 230px; }
.tablewrap { overflow-x: auto; }
table.results { font-size: 13px; min-width: 760px; }
table.results th { text-align: left; font-weight: 600; font-size: 11px; text-transform: uppercase;
                   letter-spacing: 0.04em; color: var(--muted); padding: 6px 10px;
                   border-bottom: 1px solid var(--line); }
table.results td { padding: 5px 10px; border-bottom: 1px solid var(--line); }
table.results tr:hover td { background: var(--panel); }
.num { text-align: right; font-variant-numeric: tabular-nums; white-space: nowrap; }
.type { color: var(--muted); white-space: nowrap; }
.noise { color: var(--muted); font-size: 12px; }
.delta, .winner { font-weight: 650; }
.winner { white-space: nowrap; }
.good { color: var(--good); }
.bad { color: var(--bad); }
.flat { color: var(--flat); font-weight: 400; }
.left-ahead { color: var(--left); }
.right-ahead { color: var(--right); }
tr.missing td { color: var(--muted); font-style: italic; }
.barcol { width: 180px; padding: 0 10px; }
.track { display: block; width: 100%; height: 9px; position: relative; }
.track::before { content: ""; position: absolute; left: 50%; top: -1px; width: 1px; height: 11px; background: var(--line); }
.bar { position: absolute; top: 0; height: 9px; border-radius: 2px; }
.track.right .bar { left: 50%; }
.track.left .bar { right: 50%; }
.bar.good { background: var(--good-bar); }
.bar.bad { background: var(--bad-bar); }
.bar.left-ahead { background: var(--left-bar); }
.bar.right-ahead { background: var(--right-bar); }
.bar.flat { background: var(--line); }
button { font: inherit; font-size: 12px; background: var(--panel); color: var(--fg); border: 1px solid var(--line);
         border-radius: 6px; padding: 4px 10px; margin-bottom: 8px; cursor: pointer; }
button:hover { border-color: var(--muted); }
</style>
'@
}

##
# @brief A sort control per table.
#
# The rendered order is the meaningful default in both cases - benchmark order for a baseline table,
# largest disagreement first for a toolchain table - so the toggle offers signed order against it and
# restores the original by replaying the row list captured at load.
##
function Get-ReportScript {
    return @'
<script>
document.querySelectorAll("table.results").forEach(function (table) {
    var body = table.tBodies[0];
    var original = Array.prototype.slice.call(body.rows);
    var sorted = false;
    var button = document.createElement("button");
    button.textContent = "sort by change";
    button.addEventListener("click", function () {
        sorted = !sorted;
        var rows = original.slice();
        if (sorted) {
            rows.sort(function (a, b) {
                return parseFloat(a.dataset.delta) - parseFloat(b.dataset.delta);
            });
        }
        rows.forEach(function (row) { body.appendChild(row); });
        button.textContent = sorted ? "original order" : "sort by change";
    });
    table.parentNode.parentNode.insertBefore(button, table.parentNode);
});
</script>
'@
}

##
# @brief Renders every measured core into one standalone HTML page.
#
# @param Cores Ordered map of core label to a hashtable of Description, Sampling, Comparisons and Diff.
# @param Path  File to write.
##
function New-ComparisonReport {
    param(
        [Parameter(Mandatory)]$Cores,
        [Parameter(Mandatory)][string]$Path
    )

    $html = New-Object System.Text.StringBuilder
    [void]$html.AppendLine('<!DOCTYPE html>')
    [void]$html.AppendLine('<html lang="en"><head><meta charset="utf-8">')
    [void]$html.AppendLine('<meta name="viewport" content="width=device-width, initial-scale=1">')
    [void]$html.AppendLine('<title>fixed_point128 benchmark comparison</title>')
    [void]$html.AppendLine((Get-ReportStyle))
    [void]$html.AppendLine('</head><body><main>')
    [void]$html.AppendLine('<h1>fixed_point128 benchmark comparison</h1>')
    [void]$html.AppendLine(('<p class="sub">Release builds. Generated {0}.</p>' -f `
                            (ConvertTo-HtmlText (Get-Date -Format 'yyyy-MM-dd HH:mm'))))

    foreach ($label in $Cores.Keys) {
        $core = $Cores[$label]
        [void]$html.AppendLine('<section>')
        [void]$html.AppendLine(('<h2>{0}</h2>' -f (ConvertTo-HtmlText $core.Description)))
        [void]$html.AppendLine(('<p class="sub">{0}</p>' -f (ConvertTo-HtmlText $core.Sampling)))

        foreach ($toolchain in $core.Comparisons.Keys) {
            [void]$html.AppendLine((New-ComparisonSection -Comparison $core.Comparisons[$toolchain]))
        }
        if ($core.Diff) {
            [void]$html.AppendLine((New-ToolchainDiffSection -Diff $core.Diff))
        }
        [void]$html.AppendLine('</section>')
    }

    [void]$html.AppendLine(('<p class="sub legend">Rates are millions of iterations per second; higher is better. ' +
                            'Changes smaller than {0}% are shown as unchanged - that is this machine''s run to run ' +
                            'noise, not a measured result. The noise column is the run to run variation of the ' +
                            'current figure, measured across the runs that survive the trim; treat any change that ' +
                            'is not several times larger than it as unproven. ' +
                            'Noise is not the only trap: any edit moves every function''s address, and a hot loop ' +
                            'landing on a different alignment can shift 20% with byte-identical instructions. A ' +
                            'benchmark that moves in code the change could not reach is layout, not a result.</p>') -f `
                           (Format-Invariant -Value $noiseThresholdPercent -Format '0.#'))
    [void]$html.AppendLine('</main>')
    [void]$html.AppendLine((Get-ReportScript))
    [void]$html.AppendLine('</body></html>')

    [System.IO.File]::WriteAllText($Path, $html.ToString(), (New-Object System.Text.UTF8Encoding($false)))
}

##
# @brief The result file name for one configuration on one core.
##
function Get-ResultFileName {
    param(
        [Parameter(Mandatory)][string]$ToolchainName,
        [Parameter(Mandatory)][string]$CoreLabel
    )

    return 'bench_{0}_{1}_release.json' -f $ToolchainName, $CoreLabel
}

# ---------------------------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------------------------

$coreTargets = @(Get-CoreTargets)
Write-Host ('Measuring on: {0}' -f (($coreTargets | ForEach-Object { $_.Description }) -join ', '))

$selected = @()
foreach ($name in $Toolchain) {
    $config = $configs[$name]
    if ($config.NeedsClangOnPath -and -not $ReportOnly) {
        Write-Host 'Locating the clang toolchain...'
        if (-not (Enable-ClangToolchain)) {
            continue
        }
    }
    $selected += $config
}

if ($selected.Count -eq 0) {
    throw 'No usable configuration; nothing to measure.'
}

if (-not $ReportOnly) {
    foreach ($config in $selected) {
        if (-not $NoBuild) {
            Build-BenchTarget -Config $config
        }

        foreach ($target in $coreTargets) {
            $report = Measure-Configuration -Config $config -Target $target
            $path = Join-Path $currentDir (Get-ResultFileName -ToolchainName $config.Name -CoreLabel $target.Label)
            [System.IO.File]::WriteAllText($path, ($report | ConvertTo-Json -Depth 6), (New-Object System.Text.UTF8Encoding($false)))
            Write-Host ('Wrote {0}' -f $path)
        }
    }
}

if ($SaveBaseline) {
    foreach ($config in $selected) {
        foreach ($target in $coreTargets) {
            $name = Get-ResultFileName -ToolchainName $config.Name -CoreLabel $target.Label
            $source = Join-Path $currentDir $name
            if (Test-Path $source) {
                Copy-Item $source (Join-Path $baselineDir $name) -Force
            }
        }
    }
    Write-Host ('Baseline updated from this run: {0}' -f $baselineDir)
}

$cores = [ordered]@{}
foreach ($target in $coreTargets) {
    $comparisons = [ordered]@{}
    $currentReports = [ordered]@{}

    foreach ($config in $selected) {
        $name = Get-ResultFileName -ToolchainName $config.Name -CoreLabel $target.Label
        $baselinePath = Join-Path $baselineDir $name
        $currentPath = Join-Path $currentDir $name

        if (-not (Test-Path $currentPath)) {
            Write-Warning ('No current result for {0} on {1}; nothing to compare.' -f $config.DisplayName, $target.Description)
            continue
        }
        $current = Get-Content $currentPath -Raw | ConvertFrom-Json
        $currentReports[$config.Name] = $current

        if (-not (Test-Path $baselinePath)) {
            Write-Warning ('No baseline for {0} on {1}; rerun with -SaveBaseline to create one.' -f $config.DisplayName, $target.Description)
            continue
        }
        $baseline = Get-Content $baselinePath -Raw | ConvertFrom-Json
        $comparisons[$config.Name] = @{
            Baseline = $baseline
            Current  = $current
            Rows     = @(Compare-Reports -Baseline $baseline -Current $current)
        }
    }

    if ($currentReports.Count -eq 0) {
        continue
    }

    # The toolchain table needs two current runs on the same core; with only one measured there is
    # nothing to diff, and the section is dropped rather than rendered empty.
    $diff = $null
    $names = @($currentReports.Keys)
    if ($names.Count -ge 2) {
        $left = $currentReports[$names[0]]
        $right = $currentReports[$names[1]]
        $diff = @{
            LeftName  = $left.toolchainName
            RightName = $right.toolchainName
            Rows      = @(Compare-Toolchains -Left $left -Right $right)
        }
    }

    $first = $currentReports[$names[0]]
    $description = $target.Description
    if ($first.PSObject.Properties.Name -contains 'coreDescription' -and $first.coreDescription) {
        $description = $first.coreDescription
    }

    $sampling = 'Pinned to CPU {0}.' -f $first.cpu
    if ($first.PSObject.Properties.Name -contains 'sampling' -and $first.sampling) {
        $sampling = '{0}, on CPU {1}.' -f $first.sampling, $first.cpu
    }

    $cores[$target.Label] = @{
        Description = $description
        Sampling    = $sampling
        Comparisons = $comparisons
        Diff        = $diff
    }
}

if ($cores.Count -eq 0) {
    Write-Host 'Nothing to report.'
    return
}

New-ComparisonReport -Cores $cores -Path $reportPath

Write-Host ''
foreach ($label in $cores.Keys) {
    Write-Host ('{0}:' -f $cores[$label].Description)
    foreach ($name in $cores[$label].Comparisons.Keys) {
        $rows   = @($cores[$label].Comparisons[$name].Rows | Where-Object { $null -ne $_.DeltaPercent })
        $faster = @($rows | Where-Object { $_.DeltaPercent -gt $noiseThresholdPercent }).Count
        $slower = @($rows | Where-Object { $_.DeltaPercent -lt -$noiseThresholdPercent }).Count
        Write-Host ('  {0,-6} {1} faster, {2} slower, {3} unchanged (of {4})' -f `
                    $name, $faster, $slower, ($rows.Count - $faster - $slower), $rows.Count)
    }
    if ($cores[$label].Diff) {
        $diffRows = @($cores[$label].Diff.Rows)
        if ($diffRows.Count -gt 0) {
            Write-Host ('  {0} vs {1}: widest gap {2} on {3} {4}' -f `
                        $cores[$label].Diff.LeftName, $cores[$label].Diff.RightName,
                        (Format-Delta -Value $diffRows[0].DeltaPercent), $diffRows[0].Type, $diffRows[0].Name)
        }
    }
}
Write-Host ''
Write-Host ('Report: {0}' -f $reportPath)
