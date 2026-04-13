SETLOCAL

SET ROOT=..\x64\Release
copy %ROOT%\k-Sys.exe .
copy %ROOT%\*.dll .

@echo off
cd /d "%~dp0"

SET paramRel=kparam_rel.txt

IF NOT "%1" == "" (
    SET paramRel=%1
)

k-Sys.exe "%paramRel%"
if "%1" == "" pause
