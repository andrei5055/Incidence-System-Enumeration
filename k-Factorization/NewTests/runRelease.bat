@echo off
SETLOCAL

SET ROOT=..\x64\Release
copy "%ROOT%\k-Sys.exe" . >nul
copy "%ROOT%\*.dll" . >nul

cd /d "%~dp0"

SET paramRel=.\paramU1F.txt
SET nopause=0

:parse
IF "%~1" == "" GOTO end_parse
IF /I "%~1" == "nopause" (
    SET nopause=1
) ELSE (
    SET paramRel=%~1
)
SHIFT
GOTO parse
:end_parse

echo Using parameter file: "%paramRel%"
k-Sys.exe "%paramRel%"

IF %nopause% == 0 pause
