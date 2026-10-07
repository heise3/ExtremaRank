param(
    [Parameter(Mandatory=$true)][string]$CurrentLibrary,
    [Parameter(Mandatory=$true)][string]$LegacyLibrary,
    [string]$OutputDirectory="results/windows-local"
)
$ErrorActionPreference = "Stop"
$rscript = (Get-Command Rscript -ErrorAction Stop).Source
if (!(Test-Path "benchmarks/mac_pseudobulk.R")) { throw "Run from the repository root" }
$CurrentLibrary = (Resolve-Path $CurrentLibrary).Path
$LegacyLibrary = (Resolve-Path $LegacyLibrary).Path
New-Item -ItemType Directory -Force $OutputDirectory | Out-Null
$outputPath = (Resolve-Path $OutputDirectory).Path
$previousRlibs = $env:R_LIBS
try {
    $env:R_LIBS = $CurrentLibrary
    & $rscript r/extremarank/tests/backends.R
    if ($LASTEXITCODE -ne 0) { throw "Native count backend correctness tests failed" }
} finally { $env:R_LIBS = $previousRlibs }
# Start-Process uses a Windows command line. Reject embedded quotes before quoting paths.
function Quoted([string]$value) {
    if ($value.Contains('"')) { throw "Paths cannot contain a double quote" }
    return '"' + $value + '"'
}
foreach ($caseName in @("sparse-100k", "sparse-250k")) {
    $inputPath = Join-Path $outputPath "$caseName.rds"
    & $rscript benchmarks/mac_pseudobulk.R prepare $caseName $inputPath
    if ($LASTEXITCODE -ne 0) { throw "Input preparation failed" }
    foreach ($method in @("released-0.5", "matrix-operator", "native")) {
        $library = if ($method -eq "released-0.5") { $LegacyLibrary } else { $CurrentLibrary }
        $reportPath = Join-Path $outputPath "$caseName-$method.json"
        $arguments = @("benchmarks/mac_pseudobulk.R", "run", (Quoted $inputPath), $method, (Quoted $reportPath), (Quoted $library))
        $process = Start-Process $rscript -ArgumentList $arguments -PassThru -NoNewWindow `
            -RedirectStandardOutput (Join-Path $outputPath "$caseName-$method.log") `
            -RedirectStandardError (Join-Path $outputPath "$caseName-$method.stderr")
        $null = $process.Handle
        [long]$peak = 0
        while (!$process.HasExited) {
            $process.Refresh()
            try { $peak = [Math]::Max($peak, $process.PeakWorkingSet64) } catch {}
            Start-Sleep -Milliseconds 20
        }
        $process.WaitForExit()
        if ($process.ExitCode -ne 0) { throw "$caseName $method failed; inspect stderr" }
        # Get the final OS high-water mark when still accessible; retain the largest sampled mark.
        try { $process.Refresh(); $peak = [Math]::Max($peak, $process.PeakWorkingSet64) } catch {}
        if ($peak -le 0) { throw "Peak process working set unavailable" }
        $report = Get-Content -Raw $reportPath | ConvertFrom-Json
        $report | Add-Member -NotePropertyName peak_working_set_bytes -NotePropertyValue $peak
        $report | Add-Member -NotePropertyName memory_measurement -NotePropertyValue "Windows OS PeakWorkingSet64 high-water mark sampled every 20 ms; final process lifetime tail may be missed"
        $report | ConvertTo-Json -Depth 20 | Set-Content -Encoding utf8 $reportPath
        Write-Host "$caseName $method exact hash passed; peak working set $peak bytes"
    }
}
