# Watchdog for the detached stage-B promcalc regeneration.
#
# WHY THIS EXISTS. Four long runs died unattended on this host in one week. Three were session
# teardown and one was a USB paging fault; the previous full run of THIS EXACT stage was killed by
# teardown at 88.4% (38,974 of 44,106 tiles). Stage B is resume-safe, but only if something actually
# re-launches it -- and the whole point of a detached ~11 h job is that nobody is watching.
#
# WHAT IT DOES. Polls for liveness every $PollSeconds. When the job is NOT alive it decides between
# three states and acts on exactly one:
#   complete        -> coverage reached $ExpectTiles. Stop. Never re-launch.
#   died, retries   -> re-launch the SAME launcher with --resume, which skips every tile already
#                      present across the main file and all shards.
#   died, exhausted -> stop and leave the evidence. Do NOT loop forever.
#
# THE RETRY CAP IS THE SAFETY PROPERTY. An unbounded re-launcher on a job that fails instantly is a
# hot loop that burns the machine all night and buries the real error under thousands of log lines.
# $MaxRetries bounds it, and the coverage check means a COMPLETED run can never be mistaken for a
# dead one and restarted over its own finished output.
#
# COVERAGE IS CHECKED ONLY WHILE THE JOB IS IDLE, and that is deliberate. Reading coverage means
# decompressing a growing multi-hundred-MB gzip, and doing that DURING the run is the recorded trap
# that turned a ~1 h job into a 5.0 h projection -- observing the sweep slowed the sweep, and the
# slowed rate was then quoted as the sweep's rate. Idle, it costs nothing and competes with nothing.
#
# NO EMBEDDED PYTHON, on purpose. The first version inlined the coverage probe as a here-string whose
# closing '@ sat indented inside a function -- a parse error, so the watchdog silently never started
# and wrote no log. The probe lives in scripts/promcalc_coverage.py and is called as a file.
#
# Launch (itself detached, so it outlives the session that started it):
#   Start-Process powershell -WindowStyle Hidden -ArgumentList `
#     '-NoProfile','-ExecutionPolicy','Bypass','-File','<repo>\scripts\promcalc_watchdog.ps1'

param(
  [string]$Repo        = "C:\Users\Laptop Portable\pythonprojects\dna_decode",
  [int]   $PollSeconds = 300,
  [int]   $MaxRetries  = 8,
  [int]   $ExpectTiles = 44106
)

$ErrorActionPreference = "Continue"

$OutDir   = Join-Path $Repo "data\promcalc"
$Launcher = Join-Path $Repo "run_promcalc_decompose_detached.bat"
$Prober   = Join-Path $Repo "scripts\promcalc_coverage.py"
$WdLog    = Join-Path $OutDir "watchdog.log"
$PyIso    = "C:\Users\Laptop Portable\tools\promcalc_env\Scripts\python.exe"

function Write-Wd([string]$msg) {
  $line = "{0}  {1}" -f (Get-Date -Format "yyyy-MM-dd HH:mm:ss"), $msg
  Add-Content -Path $WdLog -Value $line -Encoding utf8
}

function Test-JobAlive {
  # The Pool children are spawn-based, so their command line names multiprocessing rather than the
  # script. Counting ONLY the named script would read a live job as dead whenever the short-lived
  # parent is between states, so both shapes count.
  $procs = @(Get-CimInstance Win32_Process -Filter "Name='python.exe'" -ErrorAction SilentlyContinue |
             Where-Object { $_.CommandLine -like '*promcalc_decompose*' -or
                            $_.CommandLine -like '*multiprocessing*spawn*' })
  return $procs.Count
}

function Get-Coverage {
  # Idle-only. Returns -1 when coverage cannot be determined, which is treated as "not complete"
  # rather than as completion -- failing toward another resume is recoverable, declaring a partial
  # run finished is not.
  try {
    $out = & $PyIso -I $Prober $OutDir 2>$null
    if ($LASTEXITCODE -eq 0 -and "$out".Trim() -match '^\d+$') { return [int]"$out".Trim() }
    return -1
  } catch { return -1 }
}

Write-Wd "watchdog start (poll ${PollSeconds}s, max ${MaxRetries} retries, expect ${ExpectTiles} tiles)"

$retries = 0
while ($true) {
  Start-Sleep -Seconds $PollSeconds

  if ((Test-JobAlive) -gt 0) { continue }   # healthy: say nothing, keep the log readable

  # --- not alive: decide which of the three states this is -----------------------------------
  $sentinel = Join-Path $OutDir "decompose_run.sentinel"
  if (Test-Path $sentinel) { $exitTxt = (Get-Content $sentinel -Raw).Trim() } else { $exitTxt = "ABSENT" }
  $cov = Get-Coverage
  Write-Wd "job not alive. sentinel=[$exitTxt] coverage=$cov/$ExpectTiles"

  if ($cov -ge $ExpectTiles) {
    Write-Wd "COMPLETE at $cov tiles. Stopping; no further launches."
    break
  }
  if ($cov -lt 0) {
    Write-Wd "coverage UNDETERMINED -- treating as incomplete and resuming (safer than declaring done)"
  }
  if ($retries -ge $MaxRetries) {
    Write-Wd "retry cap ${MaxRetries} reached at coverage $cov. Stopping and leaving the evidence."
    break
  }

  $retries++
  Write-Wd "relaunch ${retries}/${MaxRetries} with --resume from coverage $cov"
  # ABSENT sentinel means killed before the exit line (the teardown signature); EXIT=-1 means it ran
  # and crashed in native code. Both resume identically, so the distinction is recorded, not branched.
  Start-Process -FilePath $Launcher -ArgumentList "--resume" -WorkingDirectory $Repo -WindowStyle Hidden
  Start-Sleep -Seconds 90
  if ((Test-JobAlive) -eq 0) { Write-Wd "WARNING: relaunch ${retries} did not come up within 90s" }
}

Write-Wd "watchdog exit after $retries relaunch(es)"
