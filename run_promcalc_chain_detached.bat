@echo off
REM Detached launcher for the stage-B -> validate -> stage-C chain.
REM
REM The chain itself needs nothing but stdlib, so it runs under the plain interpreter; stage C is
REM invoked inside the chain through `uv run` because it needs the project env.
REM
REM Sentinel discipline, learned the hard way: the redirect MUST precede the echo. Written as
REM `echo EXIT=%RC%>>file`, cmd reads the digit before `>>` as a file handle and the line vanishes.

setlocal
cd /d "%~dp0"

set "SENTINEL=data\promcalc\chain_sentinel.txt"

python scripts\promcalc_chain_stage_c.py
set RC=%ERRORLEVEL%

>>"%SENTINEL%" echo CHAIN_EXIT=%RC% at %DATE% %TIME%
endlocal
