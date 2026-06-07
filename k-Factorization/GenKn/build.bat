@echo off
setlocal enabledelayedexpansion

:: Locate vswhere
set VSWHERE="%ProgramFiles(x86)%\Microsoft Visual Studio\Installer\vswhere.exe"
if not exist %VSWHERE% (
    set VSWHERE="%ProgramFiles%\Microsoft Visual Studio\Installer\vswhere.exe"
)

if not exist %VSWHERE% (
    echo [ERROR] vswhere.exe not found. Visual Studio is required to build this project.
    exit /b 1
)

:: Get VS installation path
for /f "usebackq tokens=*" %%i in (`%VSWHERE% -latest -property installationPath`) do (
    set VS_PATH=%%i
)

if not defined VS_PATH (
    echo [ERROR] Could not find Visual Studio installation path.
    exit /b 1
)

:: Check for vcvars64.bat or vcvarsall.bat
set VCVARS=
if exist "!VS_PATH!\VC\Auxiliary\Build\vcvars64.bat" (
    set VCVARS="!VS_PATH!\VC\Auxiliary\Build\vcvars64.bat"
) else if exist "!VS_PATH!\VC\Auxiliary\Build\vcvarsall.bat" (
    set VCVARS="!VS_PATH!\VC\Auxiliary\Build\vcvarsall.bat" x64
)

if "!VCVARS!" == "" (
    echo [ERROR] Could not find vcvars64.bat or vcvarsall.bat.
    exit /b 1
)

echo [INFO] Setting up MSVC environment...
call !VCVARS! >nul
echo [INFO] Compiling GenKn.exe with AVX2 and O2 optimizations...
cl /O2 /arch:AVX2 /std:c++20 /EHsc "%~dp0main.cpp" "%~dp0generateKn.cpp" /Fe:"%~dp0GenKn.exe"
if !ERRORLEVEL! neq 0 (
    echo [ERROR] Compilation failed.
    exit /b 1
)
echo [INFO] Build successful.
