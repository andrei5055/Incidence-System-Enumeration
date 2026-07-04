# K18A2 Bug Tracking & Tracing Documentation

This document records the exact backup version that introduced the K18A2 classification bug (which caused $L=10$ to generate only 5 isomorphism classes instead of 6) and the script used to trace the regression.

## 1. Tracing Summary

* **Working Reference ZIP**: `ZIP/TripleSys_260619_184453.tar` (June 19, 2026, 18:44:53) -> Generates all **6** isomorphism classes.
* **First ZIP with Bug**: `ZIP/TripleSys_260620_000844.tar` (June 20, 2026, 00:08:44) -> Generates only **5** isomorphism classes.

The bug was introduced in the changes committed between **June 19, 2026, 18:44:53** and **June 20, 2026, 00:08:44**.

---

## 2. Test Automation Script (`test_zip.ps1`)

We used the following PowerShell script to automate extracting the source files, setting `min_cycle_length = 10` in `k18a2.h`, rebuilding the project, and counting the generated P-matrix results:

```powershell
param(
    [string]$tarName
)

$tarPath = "c:\TripleSys\Tt4\ZIP\$tarName"
if (-not (Test-Path $tarPath)) {
    Write-Error "Tar file not found: $tarPath"
    exit 1
}

# 1. Clean previous results
$folders = @(
    "c:\TripleSys\Logs_18-17-A2_Tmp",
    "c:\TripleSys\Logs_18-17-A2_Tmp1",
    "c:\TripleSys\Logs_18-17-A2_Result"
)
foreach ($f in $folders) {
    if (Test-Path $f) {
        Remove-Item -Recurse -Force $f | Out-Null
    }
}

# 2. Extract files
cd c:\TripleSys\Tt4\scratch
if (Test-Path c:\TripleSys\Tt4\scratch\TripleSys) {
    Remove-Item -Recurse -Force c:\TripleSys\Tt4\scratch\TripleSys | Out-Null
}
tar -xf $tarPath ./TripleSys/Include/k18a2.h ./TripleSys/Source/k18a2support.cpp

# 3. Copy to workspace
Copy-Item -Force c:\TripleSys\Tt4\scratch\TripleSys\Include\k18a2.h c:\TripleSys\Tt4\TripleSys\Include\k18a2.h
Copy-Item -Force c:\TripleSys\Tt4\scratch\TripleSys\Source\k18a2support.cpp c:\TripleSys\Tt4\TripleSys\Source\k18a2support.cpp

# Touch files to force MSBuild compile
(Get-Item c:\TripleSys\Tt4\TripleSys\Include\k18a2.h).LastWriteTime = Get-Date
(Get-Item c:\TripleSys\Tt4\TripleSys\Source\k18a2support.cpp).LastWriteTime = Get-Date

# 4. Modify min_cycle_length to 10 in k18a2.h
$h_content = Get-Content c:\TripleSys\Tt4\TripleSys\Include\k18a2.h
$h_content = $h_content -replace 'min_cycle_length\s*=\s*\d+', 'min_cycle_length = 10'
$h_content | Set-Content c:\TripleSys\Tt4\TripleSys\Include\k18a2.h

# 5. Build
cd c:\TripleSys\Tt4
cmd.exe /c "build.bat nopause" | Out-Null

# 6. Run test
cd c:\TripleSys\Tt4\NewTests
$out = cmd.exe /c "run18a2.bat nopause"

# 7. Count files
$files17 = @()
$files16 = @()
$p17 = "c:\TripleSys\Logs_18-17-A2_Result\Complete_graphs\18\18x17x2_18"
$p16 = "c:\TripleSys\Logs_18-17-A2_Result\Complete_graphs\18\18x16x2_18"

if (Test-Path $p17) {
    $files17 = Get-ChildItem -Path $p17 -Filter "P00*.txt"
}
if (Test-Path $p16) {
    $files16 = Get-ChildItem -Path $p16 -Filter "P00*.txt"
}

$total_files = $files17.Count + $files16.Count
Write-Host "TAR: $tarName | Files in 18x17: $($files17.Count) | Files in 18x16: $($files16.Count) | Total Results: $total_files"
```
