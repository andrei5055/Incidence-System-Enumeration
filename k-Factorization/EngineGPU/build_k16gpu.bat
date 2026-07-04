@echo off
REM ============================================================================
REM Build K16GPU.dll (CUDA engine for the K16 P1F full enumeration, k16gpu.cu).
REM
REM No GPU is needed to COMPILE. Default output embeds PTX (compute_90) which
REM the driver JIT-compiles on the target card (RTX 5090 = sm_120 Blackwell).
REM On the 5090 box with CUDA 12.8+/13.x, ALSO emit native SASS for max speed:
REM   set GENCODE=-gencode arch=compute_120,code=sm_120 -gencode arch=compute_90,code=compute_90
REM   build_k16gpu.bat
REM
REM NOTE: compile and link are two separate nvcc steps ON PURPOSE -- the
REM single-step "--shared" build crashes cudafe++ with CUDA 12.5 + MSVC 19.44.
REM -allow-unsupported-compiler covers minor CUDA/MSVC version skew.
REM ============================================================================
cd /d "%~dp0"
if "%GENCODE%"=="" set GENCODE=-gencode arch=compute_90,code=compute_90

REM Locate MSVC host compiler for nvcc (vswhere), fall back to cl on PATH.
set "CCBIN="
for /f "usebackq delims=" %%i in (`"%ProgramFiles(x86)%\Microsoft Visual Studio\Installer\vswhere.exe" -latest -find **\Hostx64\x64\cl.exe 2^>nul`) do set "CCBIN=%%~dpi"
if defined CCBIN (set "CCARG=-ccbin "%CCBIN%"") else set "CCARG="

nvcc -O3 -c k16gpu.cu -o k16gpu.obj %GENCODE% -allow-unsupported-compiler %CCARG% %*
if %ERRORLEVEL% NEQ 0 ( echo Compile FAILED. & exit /b 1 )
nvcc --shared k16gpu.obj -o ..\x64\Release\K16GPU.dll -allow-unsupported-compiler %CCARG%
if %ERRORLEVEL% NEQ 0 ( echo Link FAILED. & exit /b 1 )
echo Built ..\x64\Release\K16GPU.dll
