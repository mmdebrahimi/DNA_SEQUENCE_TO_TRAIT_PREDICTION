@echo off
REM Long-running runner for the GLM sweep anchor and its temperature grid points.
REM
REM =============================================================================
REM AMENDMENT 2026-10-23 -- READ THIS BEFORE THE BLOCK BELOW. Four of its claims are
REM measured WRONG and the scoreboard is stale. The old text is left in place on purpose
REM (supersede, do not delete) so the record of what was believed when stays readable.
REM
REM SCOREBOARD IS NOW 6 DEATHS / 1 COMPLETED POINT, and the deaths split into TWO
REM MECHANISTICALLY DIFFERENT CLASSES rather than the one this file describes.
REM
REM (A) DEATHS #1-#4 ARE CRASHES, exactly as block (1) below says. Unchanged, still right.
REM
REM (B) DEATHS #5 AND #6 ARE NOT CRASHES, and this is MEASURED with a positive control
REM     rather than inferred from an absence. A python crash DOES write Application
REM     id-1000 on this host -- proven by #3 and #4 on this very machine, at 10-21
REM     22:23:24 and 10-22 09:57:48, both 'Faulting application name: python.exe'.
REM     Deaths #5 (~11:29) and #6 (~14:42:19) wrote NO python id-1000 at all. A full
REM     sweep of EVERY event log carrying records (118 of them) across the #6 window
REM     returned only Microsoft-Windows-Store licensing noise about Dolby Access.
REM     So #5/#6 are EXTERNAL TERMINATIONS, not faults. Consequence: a PYTHONFAULTHANDLER
REM     or -u change buys nothing here -- there is no fault to dump and the log already
REM     ends on a clean flushed newline mid-progress.
REM
REM (C) WRONG: 'THE USB DRIVE, and this is the one that actually bites.' It bit 4 of 6.
REM     Deaths #5/#6 logged ZERO disk id-51/id-153 events in their windows (the disk
REM     burst on 10-23 is at 10:35:54, over an hour before #5), and D: was attached.
REM
REM (D) WRONG, and it was MY correction that was wrong first: cause (2) below reads an
REM     ABSENT sentinel as session teardown. Teardown is REFUTED twice over. Death #6 was
REM     launched through Win32_Process.Create (scripts/launch_detached_wmi.ps1), verified
REM     by parentage to sit under WmiPrvSE and outside this session's tree, and it died
REM     the same way. And promcalc stage B ran ~12 h to 'rc=0' plus its chain ~6 h ACROSS
REM     a session boundary, both via the same Start-Process form. So detached jobs DO
REM     survive here and the launch METHOD is not the variable. An absent sentinel means
REM     KILLED -- it does not tell you by what. THE CAUSE OF #5/#6 IS UNNAMED. Treat it as
REM     an open hypothesis, not a diagnosis. The one un-eliminated candidate is that torch
REM     mmaps 4.33 GB off the USB volume while promcalc's data sits on C:, and that is a
REM     CANDIDATE: it does not explain the total absence of disk events.
REM
REM (E) WRONG REASON, right verdict: 'Checkpointing ... would NOT have saved this one,
REM     which died before the first sequence.' True of #4 only; #5/#6 died MID-RUN at
REM     25/120 and 45/120, exactly the case checkpointing addresses. The real blocker is
REM     different and harder: dna_decode/glm/generate.py:208 generate_conditional calls
REM     torch.manual_seed(self.seed) ONCE at entry and then loops, so it is a SINGLE RNG
REM     STREAM -- sequence i's draws depend on every draw before it, and a naive resume
REM     hands sequence 45 the draws sequence 0 got. Exact resume needs per-sequence derived
REM     seeds, which CHANGES the generation procedure after two grid points are banked.
REM     That is an authority call, not a technical one.
REM
REM (F) STALE: option (b) says server1's Claude sessions read 'offline ... so nothing can
REM     be dispatched there yet'. They are reachable as of 2026-10-23 and the T=0.6 point
REM     is dispatched. If you run it there, pass --device cpu: the anchor and T=0.8 were
REM     both CPU-generated and CUDA floating-point yields different sampled sequences, so
REM     a CUDA point would be procedurally inconsistent inside a pre-registered grid.
REM =============================================================================
REM
REM TWO SEPARATE THINGS KILL THIS JOB. Both are diagnosed; neither is a defect in the job.
REM
REM (1) THE USB DRIVE, and this is the one that actually bites. MEASURED 2026-10-21, died at
REM     cond-gen 65/120 with sentinel EXIT=-1 and NO Python traceback anywhere in the log.
REM     torch MEMORY-MAPS the 4.4 GB model.safetensors off D:. When that drive hiccups, the
REM     mapped pages fault at the OS level and the process dies in NATIVE code below Python,
REM     so no traceback can ever be printed. Windows event-log signature, 30 s apart:
REM       disk id 153  'IO operation at LBA ... for Disk 1 was retried'
REM       disk id  51  'error detected on \Device\Harddisk1\DR1 during a PAGING operation'
REM       app  id 1000 'Faulting application python.exe, faulting module VCRUNTIME'
REM     RAM was NOT the cause (15.9 GB total, 8.7 GB free). The drive had already detached
REM     once earlier the same day, mid-test-suite. If a run dies with EXIT=-1 and a clean log,
REM     CHECK THE DISK EVENT LOG before suspecting the model or the code.
REM     Mitigation, NOT applied here because it collides with the standing 'heavy files on D:,
REM     never C:' rule and is the user's call: stage the weights on local storage and point
REM     HF_HOME at that copy, so nothing is memory-mapped off USB.
REM
REM     CONFIRMED AND SHARPENED 2026-10-22, death #4, and the mechanism is now named by the OS
REM     rather than inferred. T=0.6 attempt 2 died ~2 min in, during the weight load itself (the
REM     log stops at 'Loading weights: 100%'), sentinel EXIT=-1:
REM       app  id 1000 'python.exe, faulting module torch_cpu.dll, Exception code: 0xc0000006'
REM       disk id  51  x4, same second, '\Device\Harddisk1\DR1 during a paging operation'
REM     0xc0000006 IS STATUS_IN_PAGE_ERROR -- the literal 'a memory-mapped page could not be read
REM     from its backing store' code. The earlier run blamed VCRUNTIME; this one names torch_cpu
REM     and the exception code, so the mmap-off-USB account is no longer an inference.
REM
REM     AND THE DRIVE DID NOT COME BACK. Measured minutes later: `Get-Disk` lists only the NVMe
REM     and the empty card reader, `Get-PSDrive D` errors 'A drive with the name D does not
REM     exist', and the morning's Ntfs health events name D: under TWO different device paths --
REM     HarddiskVolume9 at 08:46, HarddiskVolume6 at 09:38. The volume re-enumerated between the
REM     two sweep attempts and then left. So the 2026-10-21 letter pin WORKED and is not the
REM     issue: the letter is stable, the physical USB link is not.
REM
REM     WALL CLASS: external-wall, not code-closable. No launcher change makes a 4.4 GB mmap survive a
REM     disk that detaches, and the weights themselves live on the missing volume, so the model
REM     cannot be loaded at all right now. Checkpointing the generation loop -- the right fix for
REM     the 65/120 and 35/120 mid-run deaths -- would NOT have saved this one, which died before
REM     the first sequence. Scoreboard: 4 deaths / 1 completed point (T=0.8).
REM     The three real options, two of which are the user's to pick:
REM       (a) replug / different port / different cable, then relaunch     <- user, cheapest
REM       (b) run it on server1 (GPU + local storage)                      <- user; its Claude
REM           sessions read `offline` in ListAgents as of 2026-10-22 09:xx, so nothing can be
REM           dispatched there yet
REM       (c) stage the 4.33 GiB weights on local storage (C: had 58.8 GB free) <- AUTHORITY call:
REM           it contradicts the standing 'never C:' rule, so it is not taken unilaterally
REM
REM (2) THIS FILE DOES NOT DETACH ANYTHING -- detaching is the CALLER's job. The name says
REM     'detached' and that describes the MOTIVE; there is no `start`, no `/b`, nothing that
REM     reparents. Launched from an agent's background shell it stays parented to a shell that
REM     is torn down at session end, which is how the 95/120 and 15/120 deaths happened.
REM     LAUNCH IT LIKE THIS (PowerShell) so it outlives the session:
REM       Start-Process -FilePath <repo>\run_glm_anchor_detached.bat ^
REM         -ArgumentList '<out-json>','<sentinel>','<log>','<temperature>' ^
REM         -WorkingDirectory <repo> -WindowStyle Hidden
REM     VERIFIED 2026-10-21: launched that way, the job outlived a real session boundary that
REM     had killed its predecessor -- then died to cause (1) sixteen minutes later.
REM
REM Cost per point MEASURED 2026-10-21: ~25-30 min, not 22. 12.4 s/seq x 120 plus ~3 min of
REM model load and a 173 s cold-start sequence before the knob check.
REM The out-JSON is written ONLY at the end, so any mid-run death costs elapsed time and
REM corrupts nothing. Read the sentinel to tell the two cases apart: EXIT=-1 means it RAN and
REM crashed (cause 1); an ABSENT sentinel means it was KILLED before the .bat reached its exit
REM line (cause 2). Same considerations apply to run_mechanism_audit_detached.bat (~1 hr).
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

REM --- D:, not E:, and that is now CORRECT rather than aspirational. 2026-10-21 the
REM --- External1TB volume was pinned to D: permanently (its MountedDevices entry is a
REM --- 24-byte disk-signature binding, so it survives reboots and replugs); the card
REM --- reader that had been squatting on D: was evicted with `mountvol D: /D`. These
REM --- lines said E: because the letter had drifted, which is the problem that got fixed
REM --- rather than worked around. `D:/hf_cache` holds the 4.33 GiB GENERator-1.2B weights.
set UV_CACHE_DIR=D:/uv_cache
set HF_HOME=D:/hf_cache
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
  --model auto --cache-dir D:/hf_cache --device cpu ^
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
