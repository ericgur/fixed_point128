# Notification hook (auto-installed by the Claude Code VS extension). The CLI fires this when Claude
# needs the user's attention (a permission prompt, or it went idle waiting for input); we forward the
# message to the VS bridge's /notify endpoint so the extension can raise an in-IDE notification
# (InfoBar + taskbar flash). Observe-only: always exits 0 so the CLI is never blocked.
$ErrorActionPreference = 'Stop'
try {
    # Read stdin as UTF-8 (default console input encoding garbles non-ASCII).
    $stdin = New-Object System.IO.StreamReader([Console]::OpenStandardInput(), [System.Text.Encoding]::UTF8)
    $p = $stdin.ReadToEnd() | ConvertFrom-Json
    $msg = [string]$p.message
    if (-not $msg) { exit 0 }

    # Find the Visual Studio bridge: the MOST-SPECIFIC workspace match (longest workspaceFolders prefix
    # of this cwd) whose port is actually listening. Avoids two failure modes: a parent-folder instance
    # (e.g. the repo root) shadowing a subfolder, and a stale "zombie" lockfile (dead instance, recycled
    # PID) whose port no longer answers.
    function Test-Port([int]$pt) {
        try {
            $c = New-Object System.Net.Sockets.TcpClient
            $live = $c.BeginConnect('127.0.0.1', $pt, $null, $null).AsyncWaitHandle.WaitOne(300) -and $c.Connected
            $c.Close(); return $live
        } catch { return $false }
    }
    $ideDir = Join-Path $env:USERPROFILE '.claude\ide'
    $cands = @()
    foreach ($f in Get-ChildItem $ideDir -Filter *.lock -ErrorAction SilentlyContinue) {
        try {
            $j = Get-Content -Raw $f.FullName | ConvertFrom-Json
            if ($j.ideName -ne 'Visual Studio') { continue }
            $ws = if ($j.workspaceFolders) { [string]$j.workspaceFolders[0] } else { '' }
            $match = [bool]($ws -and $p.cwd -and ($p.cwd -like ($ws + '*')))
            $cands += [pscustomobject]@{ Port = [int]$f.BaseName; Token = $j.authToken; Score = (([int]$match) * 1000000 + $ws.Length) }
        } catch { }
    }
    $port = $null; $token = $null
    foreach ($cand in ($cands | Sort-Object Score -Descending)) {
        if (Test-Port $cand.Port) { $port = $cand.Port; $token = $cand.Token; break }
    }
    if (-not $port) { exit 0 }

    $body = @{ message = $msg } | ConvertTo-Json -Compress
    $bytes = [System.Text.Encoding]::UTF8.GetBytes($body)
    Invoke-RestMethod -Uri "http://127.0.0.1:$port/notify" -Method Post `
        -ContentType 'application/json; charset=utf-8' `
        -Headers @{ 'x-claude-code-ide-authorization' = $token } `
        -Body $bytes -TimeoutSec 5 | Out-Null
} catch { }
exit 0
