<#
.SYNOPSIS
    Writes the GitHub release notes of a fixed_point128 version from its CHANGELOG.md section.

.DESCRIPTION
    Takes the "## [<Version>] - <date>" section of CHANGELOG.md as it is at a git ref, and formats it
    the way the GitHub releases are written: the "###" headings become "##", the lines of a wrapped
    paragraph or bullet are joined into one (GitHub renders the release body with hard line breaks),
    nested bullets keep their own lines, and a footer links the full CHANGELOG at the release tag.

.PARAMETER Version
    Four part version, e.g. 1.1.0.0. The section header must already carry its date.

.PARAMETER Ref
    Git ref to read CHANGELOG.md from. Defaults to HEAD; use origin/main once the release is merged.

.PARAMETER OutFile
    Where to write the notes. Defaults to out\release-notes-<Version>.md under the repository root.

.EXAMPLE
    .\New-Fp128ReleaseNotes.ps1 -Version 1.1.0.0 -Ref origin/main
#>
param(
    [Parameter(Mandatory)][ValidatePattern('^\d+\.\d+\.\d+\.\d+$')][string]$Version,
    [string]$Ref = 'HEAD',
    [string]$OutFile
)

$ErrorActionPreference = 'Stop'

# scripts -> release -> skills -> .claude -> repo root
$repoRoot = (Resolve-Path (Join-Path $PSScriptRoot '..\..\..\..')).Path
if (-not $OutFile) {
    $OutFile = Join-Path $repoRoot "out\release-notes-$Version.md"
}

$changelog = (& git -C $repoRoot show "${Ref}:CHANGELOG.md" | Out-String) -replace "`r", ''
if ($LASTEXITCODE -ne 0) {
    throw "CHANGELOG.md not found at $Ref"
}

$escaped = [regex]::Escape($Version)
$match = [regex]::Match($changelog, "(?ms)^## \[$escaped\] - \d{4}-\d{2}-\d{2}\n(.*?)(?=^## \[|\z)")
if (-not $match.Success) {
    throw "No dated '## [$Version] - <date>' section in CHANGELOG.md at $Ref"
}
$section = $match.Groups[1].Value -replace '(?m)^\[[^\]]+\]: .*$', ''  # the link references below the last entry

$lines = New-Object System.Collections.Generic.List[string]
foreach ($line in ($section -split "`n")) {
    if ($line.Trim() -eq '') {
        $lines.Add('')
        continue
    }
    $startsNew = $line.StartsWith('#') -or $line.StartsWith('- ') -or ($line -match '^\s+- ')
    $previous = if ($lines.Count -gt 0) { $lines[$lines.Count - 1] } else { '' }
    if (-not $startsNew -and $previous -ne '' -and -not $previous.StartsWith('#')) {
        $lines[$lines.Count - 1] = $previous + ' ' + $line.Trim()
    } elseif ($line.StartsWith('###')) {
        $lines.Add($line.Substring(1))
    } else {
        $lines.Add($line)
    }
}

$body = (($lines -join "`n") -replace "\n{3,}", "`n`n").Trim()
$body += "`n`nThe full entry is in [CHANGELOG.md](https://github.com/ericgur/fixed_point128/blob/v$Version/CHANGELOG.md).`n"

New-Item -ItemType Directory -Force (Split-Path $OutFile) | Out-Null
[System.IO.File]::WriteAllText($OutFile, $body, (New-Object System.Text.UTF8Encoding($false)))
Write-Host "Wrote $OutFile ($(($body -split "`n").Count) lines)"
