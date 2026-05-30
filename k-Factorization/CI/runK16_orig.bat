@echo off
cd /d "%~dp0"
setlocal enabledelayedexpansion

SET "ROOT=..\x64\Release"
copy "%ROOT%\k-Sys.exe" .
copy "%ROOT%\*.dll" .

k-Sys.exe paramK16p1f.txt
pause
