@echo off

set ARX=TripleSys
set YEAR=%DATE:~12,2%
set MONTH=%DATE:~4,2%
set DAY=%DATE:~7,2%
copy /Y zip.bat zip.ba
copy /Y .\CI\run.bat .\CI\run.ba
copy /Y .\CI\run_st.bat .\CI\run_st.ba

tar -c -J -f "%ARX%_%YEAR%%MONTH%%DAY%.tar" "./k-Sys/*.vcxproj*" "./k-Sys/*.cpp"  "./TripleSys/Source" "./TripleSys/Include" "./TripleSys/Manuals" "*.sln" ^
"./TripleSys/TripleSys.*" "./EngineGPU/EngineGPU.*" "./EngineGPU/*.cu" "./Utils/Include" "./Utils/Source" "./Utils/Utils.*" "zip.ba" "./CI/param*.txt" "./CI/*.ba" "./OneApp/*.cpp" "./OneApp/*.h" "./OneApp/*.vcxproj" ./OneApp/sycl_target_flags.props

del zip.ba
if exist .\CI\run.bat del .\CI\run.ba
if exist .\CI\run_st.bat del .\CI\run_st.ba