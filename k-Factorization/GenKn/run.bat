@echo off
setlocal enabledelayedexpansion

:: Default values
set "PAUSE_AT_END=1"
set "REBUILD=0"
set "N=16"

:: Parse command line arguments
set "ARGS="

:arg_loop
if "%1" == "" goto arg_done
if "%1" == "--build" (
    set REBUILD=1
    shift
    goto arg_loop
)
if "%1" == "nopause" (
    set PAUSE_AT_END=0
    shift
    goto arg_loop
)

:: Append argument to ARGS
if not defined ARGS (
    set "ARGS=%1"
) else (
    set "ARGS=!ARGS! %1"
)
shift
goto arg_loop

:arg_done
if not defined ARGS (
    set "ARGS=16"
)

:: If binary does not exist, force a build
if not exist "%~dp0GenKn.exe" set REBUILD=1

:: Locate vswhere
set VSWHERE="%ProgramFiles(x86)%\Microsoft Visual Studio\Installer\vswhere.exe"
if not exist %VSWHERE% (
    set VSWHERE="%ProgramFiles%\Microsoft Visual Studio\Installer\vswhere.exe"
)

if not exist %VSWHERE% (
    if exist "%~dp0GenKn.exe" (
        goto run
    )
    echo [ERROR] vswhere.exe not found. Visual Studio is required to build this project.
    goto end
)

:: Get VS installation path
for /f "usebackq tokens=*" %%i in (`%VSWHERE% -latest -property installationPath`) do (
    set VS_PATH=%%i
)

if not defined VS_PATH (
    if exist "%~dp0GenKn.exe" (
        goto run
    )
    echo [ERROR] Could not find Visual Studio installation path.
    goto end
)

:: Check for vcvars64.bat or vcvarsall.bat
set VCVARS=
if exist "!VS_PATH!\VC\Auxiliary\Build\vcvars64.bat" (
    set VCVARS="!VS_PATH!\VC\Auxiliary\Build\vcvars64.bat"
) else if exist "!VS_PATH!\VC\Auxiliary\Build\vcvarsall.bat" (
    set VCVARS="!VS_PATH!\VC\Auxiliary\Build\vcvarsall.bat" x64
)

if "!VCVARS!" == "" (
    if exist "%~dp0GenKn.exe" (
        goto run
    )
    echo [ERROR] Could not find vcvars64.bat or vcvarsall.bat.
    goto end
)

if %REBUILD% == 1 (
    echo [INFO] Setting up MSVC environment...
    call !VCVARS! >nul
    echo [INFO] Compiling GenKn.exe with AVX2 and O2 optimizations...
    cl /O2 /arch:AVX2 /std:c++20 /EHsc "%~dp0main.cpp" "%~dp0generateKn.cpp" /Fe:"%~dp0GenKn.exe"
    if !ERRORLEVEL! neq 0 (
        echo [ERROR] Compilation failed.
        goto end
    )
    echo [INFO] Build successful.
)

:run
echo [INFO] Running GenKn.exe with %ARGS%
"%~dp0GenKn.exe" %ARGS%

:end
if "%PAUSE_AT_END%" == "1" (
    echo.
    echo Press any key to exit...
    pause >nul
)
