@echo off
REM ============================================================================
REM  run18a2full.bat - FULL K18 |Aut|>1 classification via the representative
REM  method. All orders below run in one pass; results flow into the same
REM  pipeline, which canonicalizes the union -> one combined, deduplicated
REM  classification in the result folder.
REM
REM  ORDER SET = {3,4,5,7,11,13,17}: the prime orders {3,5,7,11,13,17} (Cauchy:
REM  every |Aut|>1 has an element of prime order) PLUS order-4 as the tractable
REM  PROXY for order-2 (order-2 itself is intractable -- see below -- so order-4
REM  is how the order-4-having 2-power classes get caught). Composite orders
REM  6,8,9,10,12,14,15,16 are intentionally DROPPED: every class they could find
REM  already has a prime-order or order-4 element in this set, so they only add
REM  redundant zero/duplicate passes. This set is COMPLETE for |Aut|>1 except
REM  the order-2-only part documented next.
REM
REM  ORDER-2 IS OMITTED -- IT IS INTRACTABLE FOR THIS METHOD. order-2 has TWO
REM  dense over-cap involution types that do not terminate: 2^8.1^2 and 2^9
REM  (verified 2026-06-27: 2^8.1^2 alone expands a 223k-task frontier with no
REM  plateau, like 2^9). The only |Aut|>1 classes this misses are those whose
REM  automorphism group is an elementary-abelian 2-group (|Aut|=2, or Klein-type
REM  with all involutions). To ATTEMPT order-2 anyway, run a SEPARATE pass with
REM  REP_ORDER=2 (the 7 small types 2^1.1^16..2^7.1^4 finish in seconds at +0,
REM  then it hangs on 2^8.1^2 / 2^9 -- kill it manually). Do NOT add 2 here: it
REM  would block the pipeline from ever producing the combined result.
REM
REM  Set worker count via KThreads in param18a2.txt (e.g. 32 on the big PC).
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
REM Complete tractable |Aut|>1 set: primes {3,5,7,11,13,17} + order-4 proxy.
REM order-2 deliberately excluded (intractable: 2^8.1^2 and 2^9). See header.
SET REP_ORDERS=3,4,5,7,11,13,17

SET "ROOT=..\x64\Release"
copy "%ROOT%\k-Sys.exe" . >nul
copy "%ROOT%\*.dll" . >nul

k-Sys.exe param18a2.txt

    if %nopause% == 0 pause
    exit /b %ERRORLEVEL%
