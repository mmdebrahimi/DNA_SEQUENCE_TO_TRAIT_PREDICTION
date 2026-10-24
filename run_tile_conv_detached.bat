@echo off
REM Long-running builder for the `conv` column of the per-tile prediction tables.
REM
REM WHY THIS IS THE RIGHT LONG RUN FOR AN OVERNIGHT WINDOW. The 3+-channel
REM orthogonal-stack arm (`orthogonal-modality-stack-promcalc-plus-N-channels`) is the one
REM `planned` arm on the ceiling criterion that is RUNNABLE: its substrate was regenerated
REM 2026-10-23 and nothing is missing except this compute. `scripts/glm_tile_hybrid_arm.py`
REM already generalises over channel count (`COMBOS`, `combine(preds, names)`,
REM `orthogonality(a, b)`, a shuffled control) and refuses with
REM INDETERMINATE_MISSING_SCORER until `conv` exists in the tables. This produces it.
REM
REM TWO-PHASE BY DESIGN, not by accident. glm_tile_predictions.py's own docstring says to
REM "run the cheap scorers first, verify the anchors, then commit the expensive conv run",
REM and `merge_payload` exists to union scorer columns into an existing seed table.
REM PHASE A ALREADY RAN (2026-10-23, attended) with --scorers gc,max_tx,decomposition and
REM the reconcile gate PASSED, reproducing the published arm exactly:
REM   seed 0: gc +0.3000  max_tx +0.3973  decomposition +0.4021
REM   ...all 10 seeds identical to wiki/promcalc_decomposition_arm_2026-10-23.json per_seed.
REM So the substrate is verified BEFORE spending the hours. This phase adds `conv` only.
REM
REM MERGE SAFETY: merge_payload unions columns for the SAME seed and REFUSES on a split
REM mismatch rather than aligning one scorer's predictions against another's tiles. Do NOT
REM pass --no-merge here; that would rewrite the seed files from scratch and destroy the
REM Phase-A columns. The flag's own help text records that footgun.
REM
REM THE 6-DEATH HISTORY DOES NOT APPLY THE SAME WAY, and that is the reason this is
REM launchable unattended at all. Deaths #1-#4 were STATUS_IN_PAGE_ERROR (0xc0000006) with
REM torch mmapping the 4.33 GB GENERator weights off the USB volume; deaths #5/#6 produced
REM NO event on any of 118 logs and remain unexplained. This job touches NO external volume:
REM inputs are data/promcalc (159 MB) + data/mpra (96 MB) on C:, and the conv is ~13k
REM parameters. The one long job that ever SURVIVED on this host (promcalc stage B, ~12 h to
REM rc=0) also had its data on C:. That is a favourable-but-not-proven difference, so this
REM still writes a sentinel and a resumable log rather than assuming success.
REM
REM Usage (launch DETACHED from PowerShell so it outlives the session):
REM   Start-Process -FilePath <repo>\run_tile_conv_detached.bat -WorkingDirectory <repo> -WindowStyle Hidden
REM
REM Progress is per-seed: each completed seed is merged into its table immediately, so a
REM mid-run death costs only the seed in flight. Read data/promcalc/tile_conv.log for the
REM per-seed lines and tile_conv.sentinel for the exit code.

setlocal
cd /d "%~dp0"

set "LOG=%~dp0data\promcalc\tile_conv.log"
set "SENTINEL=%~dp0data\promcalc\tile_conv.sentinel"

if exist "%SENTINEL%" del "%SENTINEL%"

REM --- No HF_HOME / no external cache: this job loads no pretrained weights. The conv is
REM --- trained from scratch per seed, which is exactly why it is CPU-feasible here.
set PYTHONUNBUFFERED=1

echo ==== start %DATE% %TIME% ==== > "%LOG%"

uv run python scripts/glm_tile_predictions.py ^
  --seeds 10 --scorers conv ^
  --in-dir data/promcalc --mpra-cache data/mpra ^
  --out-dir data/promcalc/tile_predictions ^
  --device cpu >> "%LOG%" 2>&1

set "RC=%ERRORLEVEL%"

echo ==== end %DATE% %TIME% rc=%RC% ==== >> "%LOG%"

REM --- The redirect is written BEFORE the echo on purpose. `echo EXIT=%RC%>> file` expands
REM --- to `echo EXIT=0>> file`, and cmd reads a digit immediately preceding `>>` as a FILE
REM --- HANDLE number -- it would redirect handle 0 (stdin) and leave the sentinel EMPTY.
REM --- Measured on run_glm_anchor_detached.bat, not theorised.
>>"%SENTINEL%" echo EXIT=%RC%

endlocal
exit /b %RC%
