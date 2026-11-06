param(
    [string]$File = '',
    [string]$QtRoot = "$env:USERPROFILE\miniconda3\envs\xrsabre\Library"
)
$ErrorActionPreference = 'Stop'
$projectRoot = Split-Path $PSScriptRoot -Parent
if (!$File) { $File = Join-Path $projectRoot 'Data\15469_grid_scan_20261014231918_1.nxs' }
if (![IO.Path]::IsPathRooted($File)) { $File = Join-Path $projectRoot $File }
$portable = Join-Path $projectRoot 'dist\XRI\XRI.exe'
if (Test-Path -LiteralPath $portable) {
    Start-Process -FilePath $portable -ArgumentList ('"' + $File + '"') -WindowStyle Normal
} else {
    $env:PATH = "$QtRoot\bin;$QtRoot\lib\qt6\bin;$env:PATH"
    $env:QT_PLUGIN_PATH = "$QtRoot\plugins;$QtRoot\lib\qt6\plugins"
    Start-Process -FilePath (Join-Path $projectRoot 'build\XRI.exe') -ArgumentList ('"' + $File + '"') -WindowStyle Normal
}
