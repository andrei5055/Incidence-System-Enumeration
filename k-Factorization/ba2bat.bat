@echo off
for /r %%f in (*.ba) do (
    move "%%f" "%%~dpnf.bat" >nul
)