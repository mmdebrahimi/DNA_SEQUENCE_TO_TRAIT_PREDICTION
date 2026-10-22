@echo off
REM ============================================================================
REM Stage B of the promcalc chain: per-site mechanistic decomposition of all
REM 44,106 GSE144621 peak tiles. Detached, resume-safe, ~5.2 h.
REM
REM WHY DETACHED. The previous full run of this exact stage was killed by
REM SESSION TEARDOWN at 88.4% (38,974 of 44,106 tiles). Three GLM sweep runs
REM died the same way in the same week. A session-parented shell is the known
REM cause, so this launches via:
REM
REM     Start-Process -FilePath <repo>\run_promcalc_decompose_detached.bat ^
REM       -WorkingDirectory <repo> -WindowStyle Hidden
REM
REM COST, MEASURED HERE 2026-10-22 -- NOT the recorded figure.
REM   The wiki records 0.26 s/tile single-core and 785 tiles/min on 7 workers
REM   (~56 min). Measured on promotercalculator 1.2.4 on this host, the solve
REM   is ~3.0 s/tile, about 11x that. The constructor is FREE (0.000 s), so the
REM   cost is in Promoter_Calculator.run() and is NOT amortizable by reusing the
REM   object. 44,106 x 3.0 s / 7 workers ~= 18,900 s ~= 5.2 h.
REM   Whether the recorded 0.26 s/tile came from a faster package version or was
REM   measured differently is UNRESOLVED; 3.0 s/tile is what this host does.
REM
REM DO NOT POLL PROGRESS BY DECOMPRESSING THE OUTPUT. Recorded trap: doing that
REM while 7 workers competed for 8 cores gave 138 tiles/min and a 5.0 h
REM projection against a real 785 tiles/min -- observing the sweep slowed the
REM sweep, and the slowed rate was then quoted as the sweep's rate. Read the
REM LOG's tail or the sentinel instead; never re-read the growing gzip.
REM
REM TWO VERSION BREAKS WERE FIXED BEFORE THIS LAUNCH, both caught by the
REM bounded --smoke gate in 24 s rather than after 5 h:
REM   1. `from promoter_calculator import Promoter_Calculator` -- 1.2.4 no longer
REM      re-exports the CLASS at the top level (only the lowercase FUNCTION).
REM   2. each per-TSS value is now a PromoCalcResults object with __getitem__ but
REM      NO .get(), which the worker called.
REM The CLASS ITSELF IS UNCHANGED (1.2.4's own wrapper.py drives it with the same
REM ctor and reads the same output key), so this is an import-PATH fix and stage
REM C's frozen reconcile anchors remain the right acceptance test.
REM
REM ON DEATH: re-launch with --resume appended. It skips tile_ids present across
REM the main file AND every shard, drops a truncated shard's LAST tile (its row
REM set may be partial), and writes a new decomposition.part-NN.csv.gz, because
REM a gzip stream cannot be safely appended.
REM
REM SENTINEL: EXIT=0 clean. EXIT=-1 ran and crashed in native code. ABSENT
REM sentinel means killed before reaching the exit line (the teardown signature).
REM ============================================================================

set "OUTDIR=%~dp0data\promcalc"
set "LOCKFILE=%OUTDIR%\decompose.lock"
set "LOG=%OUTDIR%\decompose_run.log"
set "SENTINEL=%OUTDIR%\decompose_run.sentinel"
set "PY=C:\Users\Laptop Portable\tools\promcalc_env\Scripts\python.exe"

if not exist "%OUTDIR%" mkdir "%OUTDIR%"
del "%SENTINEL%" 2>nul

REM --- ONE BLAS THREAD PER WORKER. Measured 2026-10-22 and worth ~10 h on this job.
REM --- The solve is BLAS-parallel: one tile costs 8.7 core-seconds pinned to a single thread but
REM --- only ~3.0 s wall unpinned, so an unpinned process quietly spreads across ~3 cores. With 7
REM --- worker PROCESSES each ALSO spawning its own thread pool, 8 cores get oversubscribed several
REM --- times over and most of the machine goes into context switching: the first launch of this
REM --- exact job measured 28 tiles/min (~26 h projected). Pinning each worker to one thread makes
REM --- the process-level parallelism the only parallelism, which is the level that scales here.
REM --- promotercalculator pulls numpy + scipy + scikit-learn, so all four vendor knobs are set.
REM --- Env vars are inherited by the spawned Pool children, which is what makes this reach them.
set OMP_NUM_THREADS=1
set OPENBLAS_NUM_THREADS=1
set MKL_NUM_THREADS=1
set NUMEXPR_NUM_THREADS=1

REM --- The ISOLATED interpreter, not `uv run python`. Stage B must not see the
REM --- project env: promotercalculator pulls sklearn 1.9.1 while the project is
REM --- pinned at 1.8.0, and `uv.lock` is a frozen reproducibility unit so the
REM --- dependency cannot be added to it. The separate venv is the only route
REM --- that satisfies both.
>>"%LOG%" echo ==== start %DATE% %TIME% (args: %*) ====
uv run python "%~dp0scripts\glm_run_lock.py" run "%LOCKFILE%" --label "promcalc decompose" -- ^
  "%PY%" "%~dp0scripts\promcalc_decompose.py" ^
    --in-dir "%OUTDIR%" ^
    --workers 7 ^
    %* >>"%LOG%" 2>&1

set RC=%ERRORLEVEL%
>>"%LOG%" echo ==== end %DATE% %TIME% rc=%RC% ====

REM --- Redirect BEFORE the echo. `echo EXIT=%RC%>> file` expands to
REM --- `echo EXIT=0>> file` and cmd reads the digit immediately preceding `>>`
REM --- as a FILE HANDLE number, so it redirects handle 0 (stdin) and leaves the
REM --- sentinel holding a bare `EXIT=`.
>>"%SENTINEL%" echo EXIT=%RC%
exit /b %RC%
