# Launch a command FULLY DETACHED from the caller's process tree.
#
# WHY THIS EXISTS. `Start-Process -WindowStyle Hidden` is what run_glm_anchor_detached.bat's own
# comment recommends, and its comment records that form as "VERIFIED 2026-10-21: the job outlived a
# real session boundary". On 2026-10-23 it did NOT: the T=0.6 point died ~5 min in at cond-gen
# 25/120 with the run-lock STILL HELD, no sentinel, D: still attached, and ZERO disk id-51 paging
# errors and ZERO python.exe crash events in the event log -- so the child did not crash and the
# drive did not drop. The whole tree was terminated externally, which is the launcher's documented
# cause (2), session teardown, surviving a Start-Process launch.
#
# Win32_Process.Create spawns through the WMI provider host, so the new process's parent is NOT the
# caller and it is not in the caller's job object. That is the difference that matters here.
#
# Written as a FILE, not an inline -Command: passing \\.\root\cimv2 through bash -> powershell eats
# the leading backslashes and the cast fails with 'Invalid parameter'. A file has no shell in between.

param(
  [Parameter(Mandatory = $true)][string]$CommandLine,
  [Parameter(Mandatory = $true)][string]$WorkingDirectory
)

$ErrorActionPreference = 'Stop'

$cls = [wmiclass]'root\cimv2:Win32_Process'
$res = $cls.Create($CommandLine, $WorkingDirectory, $null)

# ReturnValue 0 == success. Anything else and the PID is meaningless, so say so rather than
# printing a null PID that reads like a launch.
if ($res.ReturnValue -ne 0) {
  Write-Output ("LAUNCH_FAILED ReturnValue={0}" -f $res.ReturnValue)
  exit 1
}

$pidNew = $res.ProcessId
Write-Output ("LAUNCHED pid={0}" -f $pidNew)

# Report the parent so a caller can VERIFY the detachment rather than trust it.
Start-Sleep -Seconds 2
$p = Get-CimInstance Win32_Process -Filter "ProcessId=$pidNew" -ErrorAction SilentlyContinue
if ($p) {
  Write-Output ("parent={0} name={1}" -f $p.ParentProcessId, $p.Name)
} else {
  Write-Output "already exited or not visible (a short-lived cmd wrapper is normal)"
}
