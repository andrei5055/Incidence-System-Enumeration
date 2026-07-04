@echo off

set D=".\k-Factorization"
set ARX=CD_Tool

set YEAR=%DATE:~12,2%
set MONTH=%DATE:~4,2%
set DAY=%DATE:~7,2%

rem cd "%D%"
rem zip.bat
rem cd ..

copy /Y zip.bat zip.ba

tar -c -J -f "%ARX%_%YEAR%%MONTH%%DAY%.tar" "./BIBD_list/BIBD_list*" "./CDTools/CDTools.sln" "./CDTools/CDTools.vcxproj*" "./CDTools/CMake*" "./CDTools_GPU/*" "./CombiLogic/CombiLogic*" "./Source" "%D%\*.tar" ".\zip.ba" 

del zip.ba

