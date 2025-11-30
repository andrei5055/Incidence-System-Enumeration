copy ..\x64\Release\k-Sys.exe .
copy ..\x64\Release\*.dll .

@echo off
cmd /k k-Sys.exe ../../param.txt

