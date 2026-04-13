@echo off
taskkill /F /IM k-Sys.exe /T >nul 2>nul
del /q x64\Release\k-Sys.exe >nul 2>nul
del /q x64\Release\TripleSys.lib >nul 2>nul
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
    exit /b 1
)

echo "Build successful."
