param(
    [string]$Revision = 'working',
    [string]$Label = 'working',
    [int]$Side = 100,
    [int]$Repetitions = 3,
    [int[]]$Threads = @(1, 2, 4, 8),
    [string]$OutputDirectory = '.local/native-benchmark',
    [string]$Compiler = '',
    [switch]$BuildOnly,
    [switch]$RunOnly
)

# This helper snapshots headers, instruments only that temporary copy, and
# tracks requested live C++ allocation payloads. It does not measure RSS, R
# memory, malloc/OpenMP runtime allocations, allocator overhead, or setup peaks.
# Example: ./tools/benchmark-native.ps1 -Revision 6cb9d30 -Label baseline
#          ./tools/benchmark-native.ps1 -Revision working -Label optimized
# Use -BuildOnly / -RunOnly to keep measurement periods free of build activity.
$ErrorActionPreference = 'Stop'
$repository = [IO.Path]::GetFullPath((Join-Path $PSScriptRoot '..'))
$destination = [IO.Path]::GetFullPath((Join-Path $repository $OutputDirectory))
if ($Label -notmatch '^[A-Za-z0-9_-]+$') { throw 'Label may contain letters, digits, underscore, and hyphen only' }
if ($BuildOnly -and $RunOnly) { throw 'BuildOnly and RunOnly cannot both be set' }
$snapshot = Join-Path $destination $Label
$null = New-Item -ItemType Directory -Force -Path $snapshot
$utf8 = New-Object System.Text.UTF8Encoding($false)
$executable = Join-Path $snapshot 'benchmark-native.exe'

if (-not $RunOnly) {
    if (-not $Compiler) {
        if (Test-Path 'C:/rtools45/x86_64-w64-mingw32.static.posix/bin/g++.exe') {
            $Compiler = 'C:/rtools45/x86_64-w64-mingw32.static.posix/bin/g++.exe'
        } else {
            $candidate = Get-Command g++ -ErrorAction SilentlyContinue
            if ($candidate) { $Compiler = $candidate.Source }
            else { throw 'Supply -Compiler with an OpenMP-capable C++ compiler' }
        }
    }
    if ($Revision -eq 'working') {
        $include = Join-Path $snapshot 'inst/include'
        $null = New-Item -ItemType Directory -Force -Path $include
        Copy-Item -LiteralPath (Join-Path $repository 'inst/include/RcppAMGCL.h') -Destination $include -Force
        Copy-Item -LiteralPath (Join-Path $repository 'inst/include/amgcl') -Destination $include -Recurse -Force
    } else {
        $archive = Join-Path $snapshot 'headers.zip'
        & git -C $repository archive --format=zip "--output=$archive" $Revision inst/include
        if ($LASTEXITCODE -ne 0) { throw 'Could not snapshot requested Git revision' }
        Expand-Archive -LiteralPath $archive -DestinationPath $snapshot -Force
        $include = Join-Path $snapshot 'inst/include'
    }
    $header = Join-Path $include 'RcppAMGCL.h'
    $originalHeader = Join-Path $snapshot 'RcppAMGCL.original.h'
    Copy-Item -LiteralPath $header -Destination $originalHeader -Force
    $headerHash = (Get-FileHash -LiteralPath $originalHeader -Algorithm SHA256).Hash.ToLowerInvariant()
    $benchmarkHash = (Get-FileHash -LiteralPath (Join-Path $repository 'tools/benchmark-native.cpp') -Algorithm SHA256).Hash.ToLowerInvariant()
    $sourceRevision = if ($Revision -eq 'working') { 'HEAD' } else { $Revision }
    $snapshotKind = if ($Revision -eq 'working') { 'working-tree snapshot; commit is base history, header hash identifies exact source' } else { 'git revision' }
    $sourceCommit = (& git -C $repository rev-parse $sourceRevision | Out-String).Trim()
    if ($LASTEXITCODE -ne 0) { throw 'Could not resolve source revision' }
    $source = [IO.File]::ReadAllText($header)
    $needle = 'result.residual = residual(rhs, result.x);'
    if (($source.Split(@($needle), [StringSplitOptions]::None).Length - 1) -ne 1) {
        throw 'Expected exactly one residual call for temporary instrumentation'
    }
    $replacement = @'
const auto benchmark_residual_begin = std::chrono::steady_clock::now();
        result.residual = residual(rhs, result.x);
        rcppamgcl_benchmark::residual_seconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - benchmark_residual_begin).count();
'@
    [IO.File]::WriteAllText($header, $source.Replace($needle, $replacement), $utf8)
    $arguments = @('-std=c++14', '-O2', '-Wall', '-Wextra', '-fopenmp', "-I$include",
                   (Join-Path $repository 'tools/benchmark-native.cpp'), '-o', $executable)
    & $Compiler @arguments
    if ($LASTEXITCODE -ne 0) { throw 'Native benchmark compilation failed' }
    $metadata = @(
        "revision=$Revision", "label=$Label", "compiler=$Compiler",
        "source_commit=$sourceCommit", "snapshot_kind=$snapshotKind",
        "original_header_sha256=$headerHash",
        "benchmark_cpp_sha256=$benchmarkHash",
        "flags=-std=c++14 -O2 -Wall -Wextra -fopenmp",
        "built_utc=$([DateTime]::UtcNow.ToString('o'))",
        'allocation_metric=requested live C++ new/new[] payload bytes after setup, input CSR moved',
        'excluded=R memory; malloc/OpenMP runtime allocations; allocator overhead; peak setup memory',
        'timing=residual call instrumented in temporary header only; one solve warmup before repetitions',
        (& $Compiler --version | Out-String).Trim()
    )
    [IO.File]::WriteAllLines((Join-Path $snapshot 'metadata.txt'), $metadata, $utf8)
}

if (-not $BuildOnly) {
    if (-not (Test-Path -LiteralPath $executable)) { throw 'Benchmark executable is missing; build it first' }
    $runMetadata = @(
        "started_utc=$([DateTime]::UtcNow.ToString('o'))",
        "side=$Side", "repetitions=$Repetitions", "threads=$($Threads -join ',')",
        "OMP_THREAD_LIMIT=$env:OMP_THREAD_LIMIT", "OMP_NUM_THREADS=$env:OMP_NUM_THREADS",
        "OMP_DYNAMIC=$env:OMP_DYNAMIC", "OMP_PROC_BIND=$env:OMP_PROC_BIND",
        "OMP_PLACES=$env:OMP_PLACES", "PROCESSOR_IDENTIFIER=$env:PROCESSOR_IDENTIFIER",
        "NUMBER_OF_PROCESSORS=$env:NUMBER_OF_PROCESSORS",
        "os=$([Environment]::OSVersion.VersionString)"
    )
    $records = New-Object 'System.Collections.Generic.List[string]'
    foreach ($threadCount in $Threads) {
        Write-Host "Measuring $Label, side=$Side, threads=$threadCount, repetitions=$Repetitions"
        $lines = @(& $executable $Side $threadCount $Repetitions)
        if ($LASTEXITCODE -ne 0) { throw "Benchmark failed with threads=$threadCount" }
        [IO.File]::WriteAllLines((Join-Path $snapshot "threads-$threadCount.csv"), $lines, $utf8)
        if ($records.Count -eq 0) { $records.Add($lines[0]) }
        for ($index = 1; $index -lt $lines.Length; ++$index) { $records.Add($lines[$index]) }
    }
    $output = Join-Path $snapshot 'results.csv'
    [IO.File]::WriteAllLines($output, $records, $utf8)
    $runMetadata += "finished_utc=$([DateTime]::UtcNow.ToString('o'))"
    [IO.File]::WriteAllLines((Join-Path $snapshot 'run-metadata.txt'), $runMetadata, $utf8)
    Write-Host "Wrote $output"
}
