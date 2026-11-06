param(
    [string]$QtRoot = "$env:USERPROFILE\miniconda3\envs\xrsabre\Library",
    [string]$Hdf5Root = $QtRoot,
    [string]$BuildDirectory = "build",
    [switch]$Clean,
    [switch]$Test,
    [switch]$Run
)
$ErrorActionPreference = 'Stop'
$projectRoot = Split-Path $PSScriptRoot -Parent
if (!(Test-Path -LiteralPath "$QtRoot\lib\cmake\Qt6\Qt6Config.cmake")) {
    throw "Qt development package missing. Supply -QtRoot <Qt SDK prefix> (Qt Widgets, Concurrent)."
}
$vswhere = "${env:ProgramFiles(x86)}\Microsoft Visual Studio\Installer\vswhere.exe"
$vsRoot = if (Test-Path -LiteralPath $vswhere) { & $vswhere -latest -products '*' -requires Microsoft.VisualStudio.Component.VC.Tools.x86.x64 -property installationPath }
if (!$vsRoot) {
    $candidate = Get-ChildItem -LiteralPath "$env:ProgramFiles\Microsoft Visual Studio" -Filter VsDevCmd.bat -Recurse -ErrorAction SilentlyContinue | Select-Object -First 1
    if (!$candidate) { throw 'MSVC C++ development tools were not found.' }
    $vsRoot = $candidate.Directory.Parent.Parent.FullName
}
$devCommand = Join-Path $vsRoot 'Common7\Tools\VsDevCmd.bat'
# Import MSVC's environment into this process; no global PATH changes.
$startInfo = [Diagnostics.ProcessStartInfo]::new()
$startInfo.FileName = $env:ComSpec
$startInfo.Arguments = '/d /s /c ""' + $devCommand + '" -no_logo -arch=x64 -host_arch=x64 >nul && set"'
$startInfo.UseShellExecute = $false
$startInfo.RedirectStandardOutput = $true
$startInfo.CreateNoWindow = $true
$process = [Diagnostics.Process]::Start($startInfo)
$environmentLines = $process.StandardOutput.ReadToEnd() -split "`r?`n"
$process.WaitForExit()
if ($process.ExitCode) { throw 'Could not initialize the MSVC environment.' }
foreach ($line in $environmentLines) {
    if ($line -match '^([^=]+)=(.*)$') { [Environment]::SetEnvironmentVariable($matches[1], $matches[2], 'Process') }
}
$cmake = Join-Path $vsRoot 'Common7\IDE\CommonExtensions\Microsoft\CMake\CMake\bin\cmake.exe'
$ctest = Join-Path (Split-Path $cmake) 'ctest.exe'
$ninja = Join-Path $vsRoot 'Common7\IDE\CommonExtensions\Microsoft\CMake\Ninja'
$env:PATH = "$QtRoot\bin;$QtRoot\lib\qt6\bin;$Hdf5Root\bin;$ninja;$env:PATH"
$env:QT_PLUGIN_PATH = "$QtRoot\plugins;$QtRoot\lib\qt6\plugins"
$env:VSLANG = '1033'
if (![IO.Path]::IsPathRooted($BuildDirectory)) { $BuildDirectory = Join-Path $projectRoot $BuildDirectory }
& $cmake -S $projectRoot -B $BuildDirectory -G Ninja '-DCMAKE_BUILD_TYPE=Release' "-DCMAKE_PREFIX_PATH=$QtRoot" "-DHDF5_ROOT=$Hdf5Root" '-DHDF5_USE_STATIC_LIBRARIES=OFF'
if ($LASTEXITCODE) { throw 'CMake configuration failed.' }
if ($Clean) {
    & $cmake --build $BuildDirectory --target clean
    if ($LASTEXITCODE) { throw 'Build cleanup failed.' }
}
& $cmake --build $BuildDirectory --parallel
if ($LASTEXITCODE) { throw 'Build failed.' }
if ($Test) {
    & $ctest --test-dir $BuildDirectory --output-on-failure
    if ($LASTEXITCODE) { throw 'Regression tests failed.' }
}
if ($Run) {
    Start-Process -FilePath (Join-Path $BuildDirectory 'XRI.exe') -ArgumentList ('"' + (Join-Path $projectRoot 'Data\15469_grid_scan_20261014231918_1.nxs') + '"') -WindowStyle Normal
}
