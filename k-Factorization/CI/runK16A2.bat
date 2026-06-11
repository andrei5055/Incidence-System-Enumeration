@echo off
cd /d "%~dp0"
setlocal enabledelayedexpansion

SET nopause=0
:parse
IF "%~1" == "" GOTO end_parse
IF /I "%~1" == "nopause" (
    SET nopause=1
)
SHIFT
GOTO parse
:end_parse

SET "ROOT=..\x64\Release"
copy "%ROOT%\k-Sys.exe" . >nul
copy "%ROOT%\*.dll" . >nul

k-Sys.exe paramK16A2.txt

:: Check if execution failed
if %ERRORLEVEL% NEQ 0 (
    echo k-Sys execution failed.
    if %nopause% == 0 pause
    exit /b %ERRORLEVEL%
)

SET "FILE1=.\ExpectedResults\Complete_graphs\16\16x15x2_16\P_OrderedMatrices.txt"
SET "FILE2=..\Logs_CK_A2\Complete_graphs\16\16x15x2_16\P_OrderedMatrices.txt"

if not exist "%FILE2%" (
    echo "Error: Output file %FILE2% was not generated."
    if %nopause% == 0 pause
    exit /b 1
)

fc "%FILE1%" "%FILE2%" >nul
if %ERRORLEVEL% EQU 0 (
    echo "Success: Files are identical"
    if %nopause% == 0 pause
    exit /b 0
) else (
    echo "Differences found between K16A2 output and expected results."
    fc "%FILE1%" "%FILE2%"
    if %nopause% == 0 pause
    exit /b 1
)
