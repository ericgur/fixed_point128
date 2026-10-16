<#
.SYNOPSIS
    The pre-release gate for fixed_point128: every test passes, every build is free of warnings, and
    the documentation generates without a warning under the same Doxygen the CI uses.

.DESCRIPTION
    Builds the two test programs - fixed_point128_tests and fixed_point128_env_tests - from a clean
    state with each configuration the project releases on, runs them, and counts the compiler
    diagnostics. A clean build is deliberate: an incremental build skips the translation units that
    did not change, and a warning in one of those would not be printed at all.

    Then builds the documentation with DoxyGen/Doxyfile. The Doxygen version is read from
    .github/workflows/docs.yml so this check can never drift from the CI: an older Doxygen does not
    understand WARN_AS_ERROR = FAIL_ON_WARNINGS_PRINT and misses warnings the CI fails on. When the
    Doxygen on PATH (or in Program Files) is not that version, the official Windows build of it is
    downloaded once into out\doxygen-<version>, which is git-ignored.

    Always runs from the repository root, wherever it is started from: cmake --preset reads
    CMakePresets.json from the current directory.

.PARAMETER SkipTests
    Skip the builds and the test runs.

.PARAMETER SkipDocs
    Skip the documentation build.

.PARAMETER Incremental
    Build incrementally instead of from clean. Faster, but a warning in an unchanged file is not
    seen, so a release check should not use it.

.EXAMPLE
    .\Invoke-Fp128ReleaseCheck.ps1
    The full gate. Exit code 0 only when every step passed.
#>
param(
    [switch]$SkipTests,
    [switch]$SkipDocs,
    [switch]$Incremental
)

$ErrorActionPreference = 'Stop'

# scripts -> release -> skills -> .claude -> repo root
$repoRoot = (Resolve-Path (Join-Path $PSScriptRoot '..\..\..\..')).Path
$logDir = Join-Path $repoRoot 'out\release-check'
New-Item -ItemType Directory -Force $logDir | Out-Null

# The configurations a release is checked on: the two Windows toolchains, and MSVC in Debug as well
# because the Debug build is where the FP128_ASSERT checks and the unoptimized code paths run.
$configs = @(
    @{ Preset = 'msvc';     Config = 'Release' },
    @{ Preset = 'msvc';     Config = 'Debug' },
    @{ Preset = 'clang-cl'; Config = 'Release' }
)
$testPrograms = @('fixed_point128_tests', 'fixed_point128_env_tests')
$results = New-Object System.Collections.Generic.List[object]

##
# @brief Runs a native command with stderr merged into the output, without letting MSBuild's habit of
#        writing to stderr on a successful build raise a terminating error.
##
function Invoke-Native {
    param([Parameter(Mandatory)][scriptblock]$Command)
    $previous = $ErrorActionPreference
    $ErrorActionPreference = 'Continue'
    try {
        $output = & $Command 2>&1 | ForEach-Object { "$_" }
    } finally {
        $ErrorActionPreference = $previous
    }
    return , @($output)
}

##
# @brief Returns the path of a doxygen.exe of the given version, downloading it when necessary.
##
function Get-Doxygen {
    param([Parameter(Mandatory)][string]$Version)

    $candidates = @()
    $onPath = Get-Command doxygen -ErrorAction SilentlyContinue
    if ($onPath) { $candidates += $onPath.Source }
    $candidates += Join-Path $env:ProgramFiles 'doxygen\bin\doxygen.exe'
    $cached = Join-Path $repoRoot "out\doxygen-$Version\doxygen.exe"
    $candidates += $cached

    foreach ($exe in $candidates) {
        if (Test-Path $exe) {
            $found = (& $exe --version 2>$null | Select-Object -First 1) -replace '\s.*$', ''
            if ($found -eq $Version) {
                return $exe
            }
        }
    }

    $url = "https://www.doxygen.nl/files/doxygen-$Version.windows.x64.bin.zip"
    $zip = Join-Path $logDir "doxygen-$Version.zip"
    Write-Host "Downloading Doxygen $Version (the CI's version) from $url"
    Invoke-WebRequest -Uri $url -OutFile $zip -UseBasicParsing
    Expand-Archive -Path $zip -DestinationPath (Split-Path $cached) -Force
    Remove-Item $zip -Force
    if (-not (Test-Path $cached)) {
        throw "doxygen.exe not found in the downloaded archive at $cached"
    }
    return $cached
}

Push-Location $repoRoot
try {
    if (-not $SkipTests) {
        foreach ($c in $configs) {
            $tag = '{0}-{1}' -f $c.Preset, $c.Config
            $buildTree = Join-Path $repoRoot ('out\build\{0}' -f $c.Preset)
            if (-not (Test-Path (Join-Path $buildTree 'CMakeCache.txt'))) {
                Write-Host "Configuring preset '$($c.Preset)'..."
                $out = Invoke-Native { cmake --preset $c.Preset }
                if ($LASTEXITCODE -ne 0) {
                    $out | Out-File -Encoding utf8 (Join-Path $logDir "$tag-configure.txt")
                    throw "cmake --preset $($c.Preset) failed"
                }
            }

            Write-Host "Building $tag$(if (-not $Incremental) { ' from clean' })..."
            $buildArgs = @('--build', $buildTree, '--config', $c.Config, '--target') + $testPrograms
            if (-not $Incremental) { $buildArgs += '--clean-first' }
            $buildArgs += @('--', '/m', '/nologo', '/v:minimal')
            $out = Invoke-Native { cmake @buildArgs }
            $buildRc = $LASTEXITCODE
            $out | Out-File -Encoding utf8 (Join-Path $logDir "$tag-build.txt")
            $diagnostics = @($out | Where-Object { $_ -match 'warning C\d+|error C\d+|: warning:|: error:' } | Sort-Object -Unique)
            $results.Add([pscustomobject]@{ Step = "build $tag"; Ok = ($buildRc -eq 0 -and $diagnostics.Count -eq 0); Detail = "rc=$buildRc, $($diagnostics.Count) compiler diagnostics" })
            $diagnostics | Select-Object -First 10 | ForEach-Object { Write-Host "  $_" }
            if ($buildRc -ne 0) {
                continue
            }

            foreach ($program in $testPrograms) {
                $exe = Join-Path $buildTree ('bin\{0}\{1}.exe' -f $c.Config, $program)
                $out = Invoke-Native { & $exe --gtest_brief=1 }
                $testRc = $LASTEXITCODE
                $out | Out-File -Encoding utf8 (Join-Path $logDir "$tag-$program.txt")
                $summary = ($out | Where-Object { $_ -match '^\[  (PASSED|FAILED)  \]|tests? from .* ran' }) -join ' | '
                $results.Add([pscustomobject]@{ Step = "test $tag $program"; Ok = ($testRc -eq 0); Detail = $summary })
            }
        }
    }

    if (-not $SkipDocs) {
        $workflow = Get-Content (Join-Path $repoRoot '.github\workflows\docs.yml') -Raw
        if ($workflow -notmatch 'DOXYGEN_VERSION:\s*([0-9.]+)') {
            throw 'DOXYGEN_VERSION not found in .github/workflows/docs.yml'
        }
        $doxygenVersion = $Matches[1]
        $doxygen = Get-Doxygen -Version $doxygenVersion

        # PROJECT_NUMBER comes from this variable, exactly as the workflow derives it
        $shared = Get-Content (Join-Path $repoRoot 'include\fp128_shared.h') -Raw
        $parts = foreach ($field in 'MAJOR', 'MINOR', 'PATCH', 'BUILD') {
            if ($shared -match "#define\s+FP128_VERSION_$field\s+(\d+)") { $Matches[1] }
        }
        $env:FP128_DOC_VERSION = $parts -join '.'

        Write-Host "Building the documentation with Doxygen $doxygenVersion (fixed_point128 $env:FP128_DOC_VERSION)..."
        $out = Invoke-Native { & $doxygen 'DoxyGen/Doxyfile' }
        $docRc = $LASTEXITCODE
        $out | Out-File -Encoding utf8 (Join-Path $logDir 'doxygen.txt')
        $warnings = @($out | Where-Object { $_ -match 'warning:|error:' })
        $warnings | Select-Object -First 10 | ForEach-Object { Write-Host "  $_" }
        $results.Add([pscustomobject]@{ Step = "docs (Doxygen $doxygenVersion)"; Ok = ($docRc -eq 0 -and $warnings.Count -eq 0); Detail = "rc=$docRc, $($warnings.Count) warnings" })
    }
} finally {
    Pop-Location
}

Write-Host ''
foreach ($r in $results) {
    Write-Host ('{0,-5} {1,-50} {2}' -f $(if ($r.Ok) { 'ok' } else { 'FAIL' }), $r.Step, $r.Detail)
}
Write-Host "Logs: $logDir"

if ($results.Count -eq 0 -or @($results | Where-Object { -not $_.Ok }).Count -gt 0) {
    Write-Host 'RELEASE CHECK FAILED'
    exit 1
}
Write-Host 'RELEASE CHECK PASSED'
exit 0
