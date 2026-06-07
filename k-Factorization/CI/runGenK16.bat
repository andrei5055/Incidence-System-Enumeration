SETLOCAL

SET ROOT=..\x64\Release
copy %ROOT%\k-Sys.exe .
copy %ROOT%\*.dll .

@echo off

SET paramRel=.\paramGenK16.txt

IF NOT "%1" == "" (
    SET paramRel=%1
)

 k-Sys.exe "%paramRel%"
pause

