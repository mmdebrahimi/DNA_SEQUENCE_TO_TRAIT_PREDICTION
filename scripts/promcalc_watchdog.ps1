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
  [string]$Repo         = "C:\Users\Laptop Portable\pythonprojects\dna_decode",
  [int]   $PollSeconds  = 300,
  [int]   $MaxRetries   = 8,
  [int]   $ExpectTiles  = 44106,
  [int]   $StallMinutes = 20
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
  #
  # **THIS MATCH IS DELIBERATELY OVER-BROAD, and Test-Progressing is what makes that safe.**
  # Tightening it to the run-lock's recorded PID looks stricter and is WORSE: on 2026-10-22 killing
  # the parent left six worker children running at 100% CPU with the lock PID already dead, so a
  # PID-only check would have read "dead", relaunched, and put a second job alongside the orphans.
  # An unrelated multiprocessing job can therefore read as "promcalc alive" here -- fail-safe toward
  # alive, because a spurious NO-OP is cheap and a double-launch is not. The cost of that choice is
  # that liveness alone cannot distinguish healthy from hung, which is exactly the gap the stall
  # detector below closes.
  $procs = @(Get-CimInstance Win32_Process -Filter "Name='python.exe'" -ErrorAction SilentlyContinue |
             Where-Object { $_.CommandLine -like '*promcalc_decompose*' -or
                            $_.CommandLine -like '*multiprocessing*spawn*' })
  return $procs.Count
}

function Get-OutputBytes {
  # Total bytes across the main file and every resume shard. Cheap: a size read, NEVER a decompress.
  # Decompressing the growing gzip during a live run is the recorded trap that turned a ~1 h job into
  # a 5.0 h projection, so progress is measured by SIZE while the job runs and by actual coverage only
  # when it is idle.
  $files = @(Get-ChildItem -Path (Join-Path $OutDir "decomposition.csv.gz"),
                                 (Join-Path $OutDir "decomposition.part-*.csv.gz") `
                           -ErrorAction SilentlyContinue)
  if (-not $files) { return -1 }
  return ($files | Measure-Object -Property Length -Sum).Sum
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

$retries    = 0
$lastBytes  = -1
$lastGrowth = Get-Date

while ($true) {
  Start-Sleep -Seconds $PollSeconds

  $alive   = Test-JobAlive
  $stalled = $false

  if ($alive -gt 0) {
    # --- PROGRESS CHECK. Liveness alone cannot tell healthy from hung, and the whole point of an
    # --- unattended 12 h run is that nobody is watching: without this, "7 processes alive, zero
    # --- tiles completing" waits forever and the night is lost.
    $bytes = Get-OutputBytes
    if ($bytes -gt $lastBytes) {
      $lastBytes  = $bytes
      $lastGrowth = Get-Date
      continue                      # alive AND progressing: stay quiet, keep the log readable
    }
    $stallMin = ((Get-Date) - $lastGrowth).TotalMinutes
    if ($stallMin -lt $StallMinutes) { continue }   # gzip writes are lumpy; not a stall yet
    Write-Wd ("STALLED: {0} process(es) alive but output has not grown for {1:N0} min (threshold {2}). Treating as dead." -f $alive, $stallMin, $StallMinutes)
    $stalled = $true
    # --- DELIBERATELY NOT KILLING ANYTHING HERE. The liveness match above is over-broad on purpose,
    # --- so an automated kill could take down an unrelated multiprocessing job. The relaunch below is
    # --- guarded by the run-lock's own PID-liveness check, which refuses a second job while a live
    # --- holder exists -- so the safe outcome is reached without this script choosing what to kill.
  }

  # --- not alive (or stalled): decide which of the three states this is ----------------------
  $sentinel = Join-Path $OutDir "decompose_run.sentinel"
  if (Test-Path $sentinel) { $exitTxt = (Get-Content $sentinel -Raw).Trim() } else { $exitTxt = "ABSENT" }
  $cov = Get-Coverage
  $why = if ($stalled) { "STALLED" } else { "not alive" }
  Write-Wd "job ${why}. sentinel=[$exitTxt] coverage=$cov/$ExpectTiles"

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

  # Reset the stall tracker so the NEW job is judged on its own progress. Without this, the pre-relaunch
  # stall timestamp carries over and the next poll would immediately re-declare a stall, burning every
  # retry in a few minutes -- the retry cap would be spent without any of them having been given time.
  $lastBytes  = Get-OutputBytes
  $lastGrowth = Get-Date
}

Write-Wd "watchdog exit after $retries relaunch(es)"
