start /wait C:\Users\%username%\AppData\Local\Programs\WinMerge\WinMergeU.exe -noninteractive -minimize /f *.*  /enableexitcode /r ../../LOGS_CI_ST ../LOGS_CI
if %ERRORLEVEL% EQU 0 (
   echo Success
) else (
   C:\Users\%username%\AppData\Local\Programs\WinMerge\WinMergeU.exe /r /f *.* ../LOGS_CI ../../LOGS_CI_ST

   rem exit /b %errorlevel%
)