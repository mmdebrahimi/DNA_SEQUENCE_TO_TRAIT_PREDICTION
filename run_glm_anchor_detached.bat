@echo off
REM Detached runner for the GLM sweep anchor and its temperature grid points.
REM
REM WHY DETACHED: each point takes ~22 min of generation (11.2 s/seq x 120 on CPU) and
REM TWO consecutive attempts died at session end -- one at 95/120, writing nothing. The
REM job was never defective; the session boundary was. Same reason
REM run_mechanism_audit_detached.bat exists for the ~1 hr mechanism audit.
REM
REM Usage:
REM   run_glm_anchor_detached.bat <out-json> <sentinel> <log> [temperature]
REM
REM With no temperature: the DEFAULT-sampling anchor (temperature stays 1.0, no knob
REM verification needed -- there is no non-default knob to verify).
REM With a temperature: adds --verify-knob, so the point REFUSES unless the knob is
REM proven to move the output. That is the 2026-10-20 dead-knob lesson as a gate.
REM
REM --compare-to is always passed: it checks the frozen geometry against the recorded
REM baseline BEFORE the 22 minutes, which is the whole point of having built it. The
REM bins=6 here is the recorded value; leaving --bins at its default 10 is exactly the
REM mistake that voided a completed run.
REM
REM Writes EXIT=<code> as the last line of <sentinel>, so completion is distinguishable
REM from death. Do NOT poll it in an until-loop -- a waiter whose job was killed outlives
REM the session (measured: two bash processes alive ~10 h with 39.5 s CPU).

setlocal
cd /d "%~dp0"

set "OUT=%~1"
set "SENTINEL=%~2"
set "LOG=%~3"
set "TEMP_ARG="
if not "%~4"=="" set "TEMP_ARG=--temperature %~4 --verify-knob"

set UV_CACHE_DIR=E:/uv_cache
set HF_HOME=E:/hf_cache
set HF_HUB_OFFLINE=1
set TRANSFORMERS_OFFLINE=1

if exist "%SENTINEL%" del "%SENTINEL%"

uv run python scripts/glm_generator_falsifier.py ^
  --model auto --cache-dir E:/hf_cache --device cpu ^
  --conditional --context 300 --natural-set promoter ^
  --n 120 --length 150 --limit-natural 1200 ^
  --modes kmer,both --markov-order 3 --bins 6 ^
  --compare-to wiki/glm_generator_falsifier_REAL_2026-10-07.json ^
  %TEMP_ARG% ^
  --out "%OUT%" > "%LOG%" 2>&1

echo EXIT=%ERRORLEVEL%>> "%SENTINEL%"
endlocal
