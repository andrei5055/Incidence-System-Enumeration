@echo off
REM ============================================================================
REM  run16a2full.bat - FULL K16 |Aut|>2 classification via the representative
REM  method (the K16 analog of run18a2full.bat). Runs automorphism orders 3,5,7
REM  in one pass; all results flow into the same pipeline, which canonicalizes
REM  the union -> one combined, deduplicated classification in the result folder.
REM
REM  Expected (literature / validation oracle): 30 distinct P1Fs with |Aut|>2,
REM    by |Aut|: {3:19, 5:5, 7:4, 14:1, 15:1} = 30.
REM  Only orders 3,5,7 are needed: K16 has no |Aut| divisible by 4, and the
REM  |Aut|=14 / |Aut|=15 classes are caught via order 7 and orders 3/5
REM  respectively (their lower-order powers). Running 4,6,8,... would only add
REM  zero-result passes.
REM
REM  NOTE: order 2 (|Aut| = a power of 2, e.g. the 59 |Aut|=2 classes) is
REM  omitted -- type 2^k is over the centralizer cap and awaits the
REM  generator-based per-level dedup (the same piece K18 L=2 needs). |Aut|=2 is
REM  < 3 anyway, so it is outside the |Aut|>2 target.
REM
REM  Set worker count via KThreads in param16a2.txt.
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

REM Automorphism orders to classify (comma-separated, no spaces).
SET REP_ORDERS=2,3,5,7

SET "ROOT=..\x64\Release"
copy "%ROOT%\k-Sys.exe" . >nul
copy "%ROOT%\*.dll" . >nul

k-Sys.exe param16a2.txt

    if %nopause% == 0 pause
    exit /b %ERRORLEVEL%
