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
set TARFILES=./*.md
set TARFILES=./k-Sys/*.vcxproj* ./k-Sys/*.cpp ./TripleSys/Source ./TripleSys/Include ./TripleSys/Manuals
set TARFILES=%TARFILES% *.sln ./TripleSys/TripleSys.* ./EngineGPU/EngineGPU.* ./EngineGPU/*.cu
set TARFILES=%TARFILES% ./Utils/Include ./Utils/Source ./Utils/Utils.* *.ba
set TARFILES=%TARFILES% ./OneApp/*.cpp ./OneApp/*.h ./OneApp/*.vcxproj ./OneApp/sycl_target_flags.props
set TARFILES=%TARFILES% ./OldTests/ExpectedResults ./OldTests/*.ba ./OldTests/*.txt
set TARFILES=%TARFILES% ./NewTests/*.ba ./NewTests/*.txt
set TARFILES=%TARFILES% ./CI/*.ba ./CI/*.txt

REM Create archive, ignore missing files
tar -c -J -f "%ARX%_%TODAY%.tar" %TARFILES%

REM Cleanup: delete all *.ba files
del /s /q /f *.ba >nul 2>&1
