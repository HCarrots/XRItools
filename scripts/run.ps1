param(
    [string]$File = '',
    [string]$QtRoot = "$env:USERPROFILE\miniconda3\envs\xrsabre\Library",
    [string]$BuildDirectory
)
$ErrorActionPreference = 'Stop'
$projectRoot = Split-Path $PSScriptRoot -Parent
if ($File) {
    if (![IO.Path]::IsPathRooted($File)) { $File = Join-Path $projectRoot $File }
    if (!(Test-Path -LiteralPath $File -PathType Leaf)) { throw "Input file not found: $File" }
} else {
    $sample = '15469_grid_scan_20261014231918_1.nxs'
    foreach ($directory in @('workspace\Data', 'Data')) {
        $candidate = Join-Path (Join-Path $projectRoot $directory) $sample
        if (Test-Path -LiteralPath $candidate -PathType Leaf) {
            $File = $candidate
            break
        }
    }
}
$portable = Join-Path $projectRoot 'dist\XRI\XRI.exe'
if (!$BuildDirectory -and (Test-Path -LiteralPath $portable -PathType Leaf)) {
    $executable = $portable
} else {
    if (!$BuildDirectory) { $BuildDirectory = 'build' }
    if (![IO.Path]::IsPathRooted($BuildDirectory)) { $BuildDirectory = Join-Path $projectRoot $BuildDirectory }
    $executable = Join-Path $BuildDirectory 'XRI.exe'
    $env:PATH = "$QtRoot\bin;$QtRoot\lib\qt6\bin;$env:PATH"
    $env:QT_PLUGIN_PATH = "$QtRoot\plugins;$QtRoot\lib\qt6\plugins"
}
if (!(Test-Path -LiteralPath $executable -PathType Leaf)) {
    throw "XRI executable not found at $executable. Run scripts\build.ps1 first."
}
$launchOptions = @{ FilePath = $executable; WindowStyle = 'Normal' }
if ($File) { $launchOptions.ArgumentList = '"' + $File + '"' }
Start-Process @launchOptions
