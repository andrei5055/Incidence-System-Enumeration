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

k-Sys.exe param16A2.txt

    if %nopause% == 0 pause
    exit /b %ERRORLEVEL%

