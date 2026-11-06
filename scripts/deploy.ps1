param(
    [string]$QtRoot = "$env:USERPROFILE\miniconda3\envs\xrsabre\Library",
    [string]$Hdf5Root = $QtRoot,
    [string]$BuildDirectory = 'build',
    [string]$OutputDirectory = 'dist\XRI'
)
$ErrorActionPreference = 'Stop'
$projectRoot = Split-Path $PSScriptRoot -Parent
& "$PSScriptRoot\build.ps1" -QtRoot $QtRoot -Hdf5Root $Hdf5Root -BuildDirectory $BuildDirectory -Test
if (![IO.Path]::IsPathRooted($BuildDirectory)) { $BuildDirectory = Join-Path $projectRoot $BuildDirectory }
if (![IO.Path]::IsPathRooted($OutputDirectory)) { $OutputDirectory = Join-Path $projectRoot $OutputDirectory }
New-Item -ItemType Directory -Path $OutputDirectory -Force | Out-Null
Copy-Item -LiteralPath "$BuildDirectory\XRI.exe" -Destination $OutputDirectory -Force
$deployTool = @("$QtRoot\bin\windeployqt.exe", "$QtRoot\lib\qt6\bin\windeployqt.exe") | Where-Object { Test-Path -LiteralPath $_ } | Select-Object -First 1
if (!$deployTool) { throw 'windeployqt.exe was not found in the Qt SDK.' }
$qtpaths = @("$QtRoot\bin\qtpaths6.exe", "$QtRoot\bin\qtpaths.exe", "$QtRoot\lib\qt6\bin\qtpaths6.exe") | Where-Object { Test-Path -LiteralPath $_ } | Select-Object -First 1
& $deployTool --qtpaths $qtpaths --release --no-translations --no-opengl-sw --compiler-runtime "$OutputDirectory\XRI.exe"
if ($LASTEXITCODE) { throw 'Qt deployment failed.' }
$cmake = (Get-Command cmake.exe -ErrorAction SilentlyContinue).Source
if (!$cmake) { $cmake = (Get-ChildItem -LiteralPath "$env:ProgramFiles\Microsoft Visual Studio" -Filter cmake.exe -Recurse | Select-Object -First 1).FullName }
& $cmake "-DAPP=$OutputDirectory/XRI.exe" "-DOUTPUT=$OutputDirectory" "-DSEARCH_PATHS=$QtRoot/bin;$Hdf5Root/bin;$QtRoot/lib/qt6/bin" -P "$PSScriptRoot\copy-runtime.cmake"
if ($LASTEXITCODE) { throw 'Runtime dependency deployment failed.' }
Copy-Item -LiteralPath "$projectRoot\README.md", "$projectRoot\THIRD_PARTY.md" -Destination $OutputDirectory -Force
# Conda packages retain their upstream notices in the package cache.
$binaryNames = [Collections.Generic.HashSet[string]]::new([StringComparer]::OrdinalIgnoreCase)
Get-ChildItem -LiteralPath $OutputDirectory -Recurse -Filter '*.dll' | ForEach-Object { [void]$binaryNames.Add($_.Name) }
foreach ($prefix in @((Split-Path $QtRoot -Parent), (Split-Path $Hdf5Root -Parent)) | Select-Object -Unique) {
    $metadataDirectory = Join-Path $prefix 'conda-meta'
    if (!(Test-Path -LiteralPath $metadataDirectory)) { continue }
    foreach ($record in Get-ChildItem -LiteralPath $metadataDirectory -Filter '*.json') {
        $metadata = Get-Content -LiteralPath $record.FullName -Raw | ConvertFrom-Json
        $matchesBinary = $false
        foreach ($entry in $metadata.files) {
            if ($binaryNames.Contains([IO.Path]::GetFileName($entry))) { $matchesBinary = $true; break }
        }
        if (!$matchesBinary) { continue }
        $cachePath = $metadata.extracted_package_dir
        if (!$cachePath) { $cachePath = $metadata.link.source }
        if (!$cachePath) { continue }
        $notices = Join-Path $cachePath 'info\licenses'
        if (Test-Path -LiteralPath $notices) {
            $destination = Join-Path $OutputDirectory ('licenses\' + $metadata.name)
            New-Item -ItemType Directory -Path $destination -Force | Out-Null
            Get-ChildItem -LiteralPath $notices | Copy-Item -Destination $destination -Recurse -Force
        }
    }
}
$eigenNotices = Join-Path $OutputDirectory 'licenses\eigen3'
New-Item -ItemType Directory -Path $eigenNotices -Force | Out-Null
Copy-Item -LiteralPath "$projectRoot\third_party\eigen3\COPYING.MPL2", "$projectRoot\third_party\eigen3\COPYING.README", "$projectRoot\third_party\eigen3\COPYING.BSD" -Destination $eigenNotices -Force
$imagingDocs = Join-Path $OutputDirectory 'docs'
New-Item -ItemType Directory -Path $imagingDocs -Force | Out-Null
Copy-Item -LiteralPath "$projectRoot\docs\imaging-migration.md", "$projectRoot\docs\imaging-comparison.json" -Destination $imagingDocs -Force
Copy-Item -LiteralPath "$projectRoot\examples\DT_Be_ROI.txt" -Destination $OutputDirectory -Force
Copy-Item -LiteralPath "$projectRoot\examples\imaging.15469.json" -Destination $OutputDirectory -Force
Write-Host "Ready: $OutputDirectory\XRI.exe"
