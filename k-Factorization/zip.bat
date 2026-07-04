@echo off
rem Define the true ESCAPE character
for /F %%a in ('echo prompt $E ^| cmd') do set "ESC=%%a"

cd /d "%~dp0"

REM Archive name
set ARX=TripleSys

REM Use PowerShell for robust date (yyyyMMdd)
for /f %%i in ('powershell -NoProfile -Command "Get-Date -Format yyMMdd"') do set TODAY=%%i

SET "clean_time=%TIME: =0%"
SET "hhmmss=%clean_time:~0,2%%clean_time:~3,2%%clean_time:~6,2%"

REM Copy all *.bat into *.ba
for /r %%f in (*.bat) do (
   copy "%%f" "%%~dpnf.ba" >nul
)

REM Build list of files for tar
set TARFILES=./k-Sys/*.vcxproj* ./k-Sys/*.cpp ./TripleSys/Source ./TripleSys/Include ./TripleSys/Manuals ./MD/*.tex
REM NOTE: use ".vcxproj*" not ".*" here. On Windows "TripleSys.*" ALSO matches an
REM extensionless name, so a stray build dir "TripleSys\TripleSys\" gets matched and
REM tar recurses its whole x64\Release tree (bloated the 260627 zip to 24MB). The
REM ".vcxproj*" form grabs the .vcxproj/.filters/.user files and never a directory.
set TARFILES=%TARFILES% *.sln ./TripleSys/TripleSys.vcxproj* ./EngineGPU/EngineGPU.vcxproj* ./EngineGPU/*.cu ./EngineGPU/*.ba
set TARFILES=%TARFILES% ./Utils/Include ./Utils/Source ./Utils/Utils.* *.ba
set TARFILES=%TARFILES% ./OneApp/*.cpp ./OneApp/*.h ./OneApp/*.vcxproj ./OneApp/sycl_target_flags.props
set TARFILES=%TARFILES% ./OldTests/ExpectedResults ./OldTests/*.ba ./OldTests/*.txt
set TARFILES=%TARFILES% ./NewTests/*.ba ./NewTests/*.txt
set TARFILES=%TARFILES% ./CI/*.ba ./CI/*.txt

REM Create archive, ignore missing files
if not exist ".\ZIP" mkdir ".\ZIP"
tar -c -J -f ".\ZIP\%ARX%_%TODAY%_%hhmmss%.tar" %TARFILES%
rem Check the exit status and print colored results
if %errorlevel% equ 0 (
    echo %ESC%[32m[ OK ] File created: ".\ZIP\%ARX%_%TODAY%_%hhmmss%.tar" %ESC%[0m
) else (
    echo %ESC%[31m[ FAULT ] Tar command failed.%ESC%[0m
)

rem Cleanup: delete all *.ba files
del /s /q /f *.ba >nul 2>&1

pause