rmdir "../Logs_CI/" /s /q
copy ..\..\CDTools\x64\Release\k-Sys.exe .
copy ..\..\CDTools\x64\Release\*.dll .
@echo off
SET _free=1000
SET _used=2000
k-Sys.exe paramU1F.txt
if %ERRORLEVEL% EQU 0 (
cmd /k k-Sys.exe param.txt
) else (
cmd /k echo "'param.txt' test set skipped because of error in the 'paramU1F.txt' test set" 
)
