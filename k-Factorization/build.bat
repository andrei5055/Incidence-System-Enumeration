@echo off
taskkill /F /IM k-Sys.exe /T >nul 2>nul
ping -n 3 127.0.0.1 >nul

SET nopause=0
:parse
IF "%~1" == "" GOTO end_parse
IF /I "%~1" == "nopause" (
    SET nopause=1
)
SHIFT
GOTO parse
:end_parse

set VSWHERE="%ProgramFiles(x86)%\Microsoft Visual Studio\Installer\vswhere.exe"
if not exist %VSWHERE% (
    echo "Error: vswhere not found. Please install Visual Studio."
    exit /b 1
)

for /f "usebackq tokens=*" %%i in (`%VSWHERE% -latest -products * -requires Microsoft.Component.MSBuild -find MSBuild\**\Bin\MSBuild.exe`) do (
  set MSBUILD="%%i"
)

if not defined MSBUILD (
    echo "Error: MSBuild not found."
    exit /b 1
)

echo Building TripleSys.sln...
%MSBUILD% TripleSys.sln /p:Configuration=Release /p:Platform=x64
if %ERRORLEVEL% neq 0 (
    echo "Build failed."
    if %nopause% == 0 pause
    exit /b 1
)

echo "Build successful."
    if %nopause% == 0 pause
    exit /b %ERRORLEVEL%

