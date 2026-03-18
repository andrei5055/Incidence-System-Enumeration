@echo off
cd /d "%~dp0"

REM Archive name
set ARX=TripleSys

REM Use PowerShell for robust date (yyyyMMdd)
for /f %%i in ('powershell -NoProfile -Command "Get-Date -Format yyMMdd"') do set TODAY=%%i

REM Copy all *.bat into *.ba
for /r %%f in (*.bat) do (
   copy "%%f" "%%~dpnf.ba" >nul
)

REM Build list of files for tar
set TARFILES=./k-Sys/*.vcxproj* ./k-Sys/*.cpp ./TripleSys/Source ./TripleSys/Include ./TripleSys/Manuals
set TARFILES=%TARFILES% *.sln ./TripleSys/TripleSys.* ./EngineGPU/EngineGPU.* ./EngineGPU/*.cu
set TARFILES=%TARFILES% ./Utils/Include ./Utils/Source ./Utils/Utils.* *.ba
set TARFILES=%TARFILES% ./OneApp/*.cpp ./OneApp/*.h ./OneApp/*.vcxproj ./OneApp/sycl_target_flags.props

REM Conditionally add CI files if they exist
set CI_FILES=
if exist "./CI/param*.txt" set CI_FILES=%CI_FILES% ./CI/param*.txt
if exist "./CI/*.ba" set CI_FILES=%CI_FILES% ./CI/*.ba

REM Append CI files to archive list
if defined CI_FILES set TARFILES=%TARFILES% %CI_FILES%

REM Create archive, ignore missing files
tar -c -J -f "%ARX%_%TODAY%.tar" %TARFILES%

REM Cleanup: delete all *.ba files
del /s /q /f *.ba >nul 2>&1
