@echo off
cd /d "%~dp0"

:: Copy all *.bat into *.ba in this folder
for %%f in (*.bat) do (
    copy "%%f" "%%~dpnf.ba" >nul
)

:: Get timestamp using PowerShell (Format: yyMMdd_HHmmss)
for /f %%i in ('powershell -NoProfile -Command "Get-Date -Format yyMMdd_HHmmss"') do set TIMESTAMP=%%i

:: Archive name
set ARCHIVE_NAME=GenKn_%TIMESTAMP%.zip

echo Creating archive %ARCHIVE_NAME%...

:: Use PowerShell to zip the source and script files
powershell -NoProfile -Command "Compress-Archive -Path main.cpp, generateKn.cpp, generateKn.h, run.ba, zip.ba -DestinationPath %ARCHIVE_NAME% -Force"

if %ERRORLEVEL% neq 0 (
    echo [ERROR] Failed to create archive.
    del *.ba >nul 2>&1
    exit /b 1
)

:: Cleanup: delete all *.ba files created in this folder
del *.ba >nul 2>&1

echo Archive %ARCHIVE_NAME% created successfully.
