@echo off
REM Detached runner for the GLM sweep anchor and its temperature grid points.
REM
REM WHY DETACHED: each point takes ~22 min of generation (11.2 s/seq x 120 on CPU) and TWO
REM consecutive attempts died at session end -- one at 95/120, writing nothing. The job was
REM never defective; the session boundary was. Same reason run_mechanism_audit_detached.bat
REM exists for the ~1 hr mechanism audit.
REM
REM Usage:
REM   run_glm_anchor_detached.bat <out-json> <sentinel> <log> [temperature]
REM
REM With no temperature: the DEFAULT-sampling anchor (temperature stays 1.0; there is no
REM non-default knob to verify, so --verify-knob would be a no-op).
REM With a temperature: adds --verify-knob, so the point REFUSES unless the knob is proven to
REM move the output. That is the 2026-10-20 dead-knob lesson as a gate.
REM
REM --compare-to is always passed: it checks the frozen geometry against the recorded baseline
REM BEFORE the 22 minutes. bins=6 here is the RECORDED value; leaving --bins at its default 10
REM is exactly the mistake that voided a completed run.

setlocal
cd /d "%~dp0"

set "OUT=%~1"
set "SENTINEL=%~2"
set "LOG=%~3"
set "TEMP_ARG="
if not "%~4"=="" set "TEMP_ARG=--temperature %~4 --verify-knob"

set "LOCKFILE=%~dp0.glm_run.lock"

set UV_CACHE_DIR=E:/uv_cache
set HF_HOME=E:/hf_cache
set HF_HUB_OFFLINE=1
set TRANSFORMERS_OFFLINE=1

if exist "%SENTINEL%" del "%SENTINEL%"

REM --- One job at a time. Two concurrent torch processes is this host's documented
REM --- access-violation hazard, and nothing used to stop a second launch. Liveness is by PID,
REM --- NOT by the sentinel: the sentinel is written only at exit, so it says nothing at all
REM --- while a job is running. A lock whose PID is gone is reclaimed, so a killed job cannot
REM --- wedge the next launch.
REM ---
REM --- `run` wraps the job rather than `acquire` ... `release` around it, and that distinction
REM --- is load-bearing: a launcher's `acquire` process EXITS immediately, so the pid it
REM --- recorded is already dead when a second launch looks, the lock is correctly reclaimed as
REM --- stale, and the second job starts anyway. Found by the two-launch manual check; the
REM --- in-process unit tests could not see it. `run` holds the lock for the CHILD's lifetime.
REM --- Exit 9 means refused-to-start, which is distinguishable from any child exit code.
REM --- Through `uv run python`, NOT bare `python`. The `--`-separator defect this depends on
REM --- is INTERPRETER-DEPENDENT: it bites under the venv (3.11.17, what uv and the test suite
REM --- use) and not under the 3.14.7 on PATH. With the fix either works, but pinning the
REM --- interpreter means the launcher no longer depends on whichever python is first on PATH.
uv run python "%~dp0scripts\glm_run_lock.py" run "%LOCKFILE%" --label "glm anchor %~1" -- ^
  uv run python scripts/glm_generator_falsifier.py ^
  --model auto --cache-dir E:/hf_cache --device cpu ^
  --conditional --context 300 --natural-set promoter ^
  --n 120 --length 150 --limit-natural 1200 ^
  --modes kmer,both --markov-order 3 --bins 6 ^
  --compare-to wiki/glm_generator_falsifier_REAL_2026-10-07.json ^
  %TEMP_ARG% ^
  --out "%OUT%" > "%LOG%" 2>&1

set "RC=%ERRORLEVEL%"

REM --- The redirect is written BEFORE the echo on purpose. `echo EXIT=%RC%>> file` expands to
REM --- `echo EXIT=0>> file`, and cmd reads a digit immediately preceding `>>` as a FILE HANDLE
REM --- number -- so it redirected handle 0 (stdin) and left the sentinel EMPTY while `EXIT=`
REM --- went to the console. That is why the first detached run's completion was undetectable
REM --- and had to be noticed from `git status` instead. Measured, not theorised.
>>"%SENTINEL%" echo EXIT=%RC%

endlocal
exit /b %RC%
