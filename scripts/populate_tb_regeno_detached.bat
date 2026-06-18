@echo off
REM Detached launcher for the ~1.6 TB RIF regeno cohort fetch (the multi-DAY v1b callability job).
REM Routes to D: (4.4 TB free). Restartable / skip-existing — safe to re-run after an interruption.
REM Run from the repo root:  scripts\populate_tb_regeno_detached.bat
REM Monitor:  type %TEMP%\tb_regeno_pop.log   (or check D:\dna_decode_cache\cryptic\regeno file count)
cd /d %~dp0\..
echo Starting RIF regeno populate (~1.6 TB -> D:). This runs for DAYS. >> %TEMP%\tb_regeno_pop.log
uv run python scripts/populate_tb_cohort.py --drug rifampicin --kind regeno --cache D:/dna_decode_cache/cryptic >> %TEMP%\tb_regeno_pop.log 2>&1
echo DONE rc=%ERRORLEVEL% >> %TEMP%\tb_regeno_pop.log
