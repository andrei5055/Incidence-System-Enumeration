SETLOCAL
SET ROOT=..\..\CDTools
if not exist %ROOT%\ (
    SET ROOT=..
)

SET ROOT=%ROOT%\x64\Release
copy %ROOT%\k-Sys.exe .
copy %ROOT%\*.dll .

@echo off

SET paramRel=..\..\param.txt

IF NOT "%1" == "" (
    SET paramRel=%1
)
cmd /k k-Sys.exe %paramRel%

