@echo off
REM ============================================================================
REM  run14a2.bat - K14 VALIDATION ORACLE via the representative method.
REM
REM  K14 (2n=14 == 2 mod 4, like K18) is the structurally-correct validator for
REM  the K18/K20 engine: it has the fixed-point-free even cycle types (2^7,
REM  4^3.2^1) that K16 (== 0 mod 4) lacks, and being small it COMPLETES order-2
REM  (the case intractable at K16/K18/K20). Ground truth (this project's own
REM  solver, Logs_CI/.../14x13x2): 23 P1Fs; the rep method reproduces the 21 with
REM  |Aut|>1.
REM
REM  ORDER SET = {2,3,7,13}: the primes dividing any K14 |Aut| (2,3,7,13; note 5
REM  never divides a K14 |Aut|). Their union = the full |Aut|>1 classification.
REM  EXPECTED (verified by the pipeline's independent Aut-checker):
REM    order-2  -> 16  {2:3,4:1,6:5,12:5,84:1,156:1}   (PROVEN; fpf 2^7 terminates)
REM    order-3  -> 17  {3:5,6:5,12:5,84:1,156:1}
REM    order-7  -> 1   {84:1}
REM    order-13 -> 1   {156:1}
REM    COMBINED -> 21  {2:3,3:5,4:1,6:5,12:5,84:1,156:1}
REM  (order-4 -> 5 {12:4,156:1} counts order-4 ELEMENTS, not |Aut| div by 4;
REM   composite orders are not clean checks -- use primes + the |Aut|>1 union.)
REM
REM  Override the orders with the REP_ORDERS / REP_ORDER env var before running
REM  (e.g.  set REP_ORDER=2  for just the order-2 completion check).
REM  Set worker count via KThreads in param14a2.txt.
REM ============================================================================
cd /d "%~dp0"
setlocal enabledelayedexpansion

SET nopause=0
:parse
IF "%~1" == "" GOTO end_parse
IF /I "%~1" == "nopause" (
    SET nopause=1
)
SHIFT
GOTO parse
:end_parse

REM The rep method is unseeded but the dispatch needs ONE 14-player 3-row starter
REM to trigger the run-once solve. Stage it from the CI output if not present.
SET "STARTDIR=..\..\Logs_14-3-A1\Complete_graphs\14\14x3x2_14"
SET "STARTSRC=..\Logs_CI\Complete_graphs\14\14x3x2\P0000000001.txt"
IF NOT EXIST "%STARTDIR%\P0000000001.txt" (
    IF EXIST "%STARTSRC%" (
        if not exist "%STARTDIR%" mkdir "%STARTDIR%"
        copy "%STARTSRC%" "%STARTDIR%\" >nul
        echo Staged K14 3-row starter from CI output.
    ) ELSE (
        echo WARNING: no K14 starter at "%STARTDIR%" and no CI source at "%STARTSRC%".
        echo          Run CI\Run.bat first to generate the 14x3x2 starters.
    )
)

SET "ROOT=..\x64\Release"
copy "%ROOT%\k-Sys.exe" . >nul
copy "%ROOT%\*.dll" . >nul

REM Default validation sweep (overridden if REP_ORDERS / REP_ORDER is already set).
IF "%REP_ORDERS%%REP_ORDER%" == "" SET REP_ORDERS=2,3,7,13

k-Sys.exe param14a2.txt

    if %nopause% == 0 pause
    exit /b %ERRORLEVEL%
