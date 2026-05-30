@echo off
cd /d "%~dp0"
setlocal enabledelayedexpansion

SET "ROOT=..\x64\Release"
copy "%ROOT%\k-Sys.exe" .
copy "%ROOT%\*.dll" .

k-Sys.exe paramK1111p1f.txt

:: Check if copy failed
if %ERRORLEVEL% NEQ 0 (
    echo Copy failed.
    pause
    exit /b
)

SET "WINMERGE=C:\Users\%username%\AppData\Local\Programs\WinMerge\WinMergeU.exe"
SET "FILE1=.\ExpectedResults\2-Partite_graphs\11\22x11x2_22\PC0000000001.txt"
SET "FILE2=..\LOGS_CK\2-Partite_graphs\11\22x11x2_22\PC0000000001.txt"

:: Run WinMerge directly
"%WINMERGE%" /noninteractive /minimize /enableexitcode "%FILE1%" "%FILE2%"

if %ERRORLEVEL% EQU 0 (
    echo "Success: Files are identical"
    pause
    exit /b
)

:: If it gets here, ERRORLEVEL was not 0
echo Differences found (Code %ERRORLEVEL%). Opening WinMerge for review.
"%WINMERGE%" "%FILE1%" "%FILE2%"

pause
