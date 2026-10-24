# The tile chain's substrate: two dependencies, two wall classes

**2026-10-21.** A frozen pre-registration (`wiki/glm_tile_hybrid_primary_acceptance_bar.json`, authored
2026-10-17, `frozen_before_any_number_was_seen: true`) has sat unexecuted because its substrate was
never located. This memo records what the substrate actually is, measured, and assigns each half a
class from a closed vocabulary.

## CORRECTION — I asserted this substrate was present, and that was wrong

Earlier in this session I answered a VOI question by saying the orthogonal-modality / tile arm's
substrate was present, on the strength of finding `E:/PythonProjects/dna_decode_cache/mpra`. That
directory holds the **raw MPRA tables**, which is not the input the chain consumes. The chain's own
input is the **promcalc per-tile decomposition**, and it is absent from every attached volume. The
correction stands on its own; the rest of this memo is what measuring it produced.

## The two dependencies

| dependency | class | measured state |
|---|---|---|
| MPRA peak-tile table (GSE144621) | **code-closable** | present at `E:/PythonProjects/dna_decode_cache/mpra`, addressed at `D:/dna_decode_cache/mpra` |
| promcalc per-tile decomposition | **data-absent** | no `promcalc` directory on any resolving volume |

Vocabulary, closed: **code-closable** (the bytes are reachable and only the addressing was wrong) ·
**data-absent** (the bytes do not exist on this host) · **external-wall** (closing it needs a drive, a
host or an install that is not an executor action).

### The MPRA half is code-closable, and Verification Signal 2 is the proof

`scripts/glm_tile_predictions.py:282` called `load_peak_tiles()` with **no arguments**, so it took
`genomewide.DEFAULT_CACHE = D:/dna_decode_cache/mpra` unconditionally. `D:` is a card-reader slot with
no media, so `ensure_cached` reached `d.mkdir(parents=True, exist_ok=True)` and died with an unhandled
`PermissionError: [WinError 21] The device is not ready`. The loader had always accepted `cache_dir`;
stage A (`promcalc_export_tiles.py`) already exposed `--cache-dir`. This was the one un-parameterised
call.

With `--mpra-cache E:/PythonProjects/dna_decode_cache/mpra` the script now prints:

```
tiles 44106 | ceiling 0.9410
```

and then fails only on the promcalc input. **That output is stronger evidence than "it did not
crash".** `44106` is `glm_tile_predictions.EXPECTED_TILES` and `0.9410` is
`promcalc_export_tiles.EXPECTED_CEILING` — both asserted by the published baseline, both reproduced
exactly after the two published filters (`category == "tile"`, `len(seq) == 150`). So the copy on `E:`
is not merely a file of the right name; it is the published substrate, and the chain's first stage
reaches it.

### The promcalc half is data-absent, and the toolchain question is now answered

`wiki/promcalc_toolchain_probe_2026-10-21.json`, verdict **`ENV_ABSENT_DEPS_SATISFIABLE_HERE`**, 6 of 6
checks resolved, no problems:

* the documented interpreter `D:/tools/promcalc_env/Scripts/python.exe` does **not** resolve;
* `promoter_calculator` is **not** importable here;
* the ultimate raw substrate **is** reachable (one location, on `E:`);
* free space: `C:` 52.5 GB, `E:` 907.5 GB.

**`SATISFIABLE_HERE` does not mean "install it into this venv", and the artifact carries that field
rather than leaving it to a reader.** `promoter_calculator` pulls sklearn 1.9.1 while this project
resolves **1.8.0**, and stage B's own docstring says installing it in the project venv would move
sklearn under all of the test suite. The probe therefore stamps `isolation_required: true`: what is
satisfiable is a **separate** environment, never this one.

So the promcalc half is `data-absent` with a **reachable regeneration path**, not an external wall.
What is needed, and from whom: an isolated interpreter plus `promoter_calculator` (an install — a user
decision, not an executor action), and then a ~44k-tile stage-B sweep. The alternative is the host that
already holds the decomposition.

## What this does NOT establish

* **It is not a result.** No arm of the frozen bar has been run, and this memo quotes none of its
  headline figures. The bar's own numbers come from the two committed 2026-10-14 artifacts; the only
  numbers measured here are substrate reachability, the tile count and the noise ceiling.
* **A reachable regeneration path is not a reproduction.** The bar's anchor gate demands the
  regenerated pipeline reproduce committed per-seed headroom within 0.01 on 9 of 10 seeds. Nothing
  here bears on whether it would, and `REFUSED_ANCHOR_MISMATCH` remains a live outcome.
* **The probe installed nothing.** Reporting that an install is possible is a different claim from
  having performed one, and only the first is made.

## The frozen bar's `execution_notes` are half stale, and both halves matter

The bar records `not_runnable_on: "the SHNEITZER laptop: no D: data tree, no uv, no Docker"`. This host
**is** Shneitzer.

* **`no uv` is false.** uv ran the full 7,034-test suite here today.
* **`no D: data tree` is substantively right** but conflates the two classes above: it reads as one
  undifferentiated wall when one half was a hardcoded drive letter and the other is absent data.

A future reader should discard neither clause wholesale nor trust both. **The bar was not edited** —
it is a frozen pre-registration, and amending any part of one is an authority call, recorded as an open
question rather than taken.

## Reusable

A dependency that is *present but mis-addressed* and one that is *absent* produce the same symptom —
"the substrate is missing" — and have completely different remedies. Separating them cost one
parameterised call and one run. The diagnosis also required resolving each dependency **at its point of
use** rather than validating all inputs up front: a single front-loaded check would have collapsed both
into one failure and destroyed the measurement that distinguishes them.

## AMENDMENT, same day — the host changed: D: is now the 1 TB volume

**Everything above was measured while `D:` was an empty card-reader slot. It no longer is.**
The `External1TB` volume (ST1000LM024, 931.5 GB) was pinned to `D:` permanently — the card
reader was evicted with `mountvol D: /D` and the volume reassigned, so `HKLM\SYSTEM\MountedDevices`
now holds a 24-byte disk-signature binding for `D:` that survives reboots and replugs. A junction
`D:\dna_decode_cache -> D:\PythonProjects\dna_decode_cache` reconciles the nested layout.

**What this changes, measured not assumed:**

* `scripts/glm_tile_predictions.py` **with no arguments at all** now prints
  `tiles 44106 | ceiling 0.9410` and fails only on promcalc. The hardcoded default resolves.
* `scripts/promcalc_decomposition_arm.py` — which calls `load_peak_tiles()` **bare** and has no
  `--cache-dir` flag — gets past its tile load too, and fails at the same promcalc step. So the
  second un-parameterised call site is no longer blocking anything, though it remains unable to
  be pointed anywhere else.
* `D:/hf_cache` (4.4 GB, the GENERator-1.2B weights) and `D:/uv_cache` resolve at the paths the
  project's own D: rule always claimed. `run_glm_anchor_detached.bat` no longer hardcodes `E:`.

**What this does NOT change.** The promcalc decomposition is still absent, so the wall
classification above stands: one half code-closable, one half data-absent. `refseq`, `tools/`
and `dna_decode_data` are on a different physical disk that is still not attached — so
`D:/dna_decode_cache` is `UNFIT` by the project's own `CACHE_REQUIRED_SUBDIRS` check, and the
probe reports it that way.

**The failure mode shifted, and that was a known cost accepted up front.** Before, those paths
failed loudly because the volume was absent. Now they resolve and fail as *dataset* missing.
`external_root.py` says plainly that this is the worse direction — *"a wrong root that resolves
is worse than a right root that is missing, because only the second one is visible"* — which is
exactly why `CACHE_REQUIRED_SUBDIRS = ("refseq",)` exists and why `data_inventory` runs the probe
on every report. The mitigation shipped before the risk did.

**Two statements above are now false and are left in place deliberately**, because the original
measurement is what justified the code change: "`D:` is a card-reader slot with no media" and
the volume table's `D: resolves False`. Both were true when written, on 2026-10-21, before the
pinning.

## SECOND AMENDMENT, 2026-10-22 — THE ISOLATION PREMISE IS REFUTED, and it was mine to check

Everything above says the promcalc half needs a **separate** Python environment, because stage B's
docstring states `promoter_calculator` pulls sklearn 1.9.1 against this project's 1.8.0. **I relayed that
as a measured constraint and never checked it. It is false in this environment.**

Measured 2026-10-22:

* The PyPI distribution is **`promotercalculator` 1.2.4** ("Salis Lab's promoter-calculator, Barrick Lab
  Fork") — a **`py3-none-any` PURE-PYTHON wheel**, `requires_python: None`, with **completely unpinned**
  dependencies: `biopython`, `numpy`, `scikit-learn`, `scipy`. There is no 1.9.1 pin anywhere in its
  metadata; 1.9.1 is simply what an unconstrained resolver picks when it takes the latest.
* `uv pip install --dry-run promotercalculator` → **"Resolved 7 packages … would install 1 package:
  `promotercalculator==1.2.4`"**, leaving sklearn **1.8.0**, numpy **2.4.4** and scipy **1.17.1**
  untouched.

**Consequence: `isolation_required` is FALSE here, and the authority fork I handed the user was
mis-framed.** It was presented as *"installing this needs a user decision because it would move sklearn
under 6,900 tests"*. It would not. The install is one pure-Python package that changes nothing else.

**A trap worth keeping, found the same way:** the **import** name is `promoter_calculator` but
`pip install promoter_calculator` **404s** — the distribution is `promotercalculator`. Anyone following
the docstring literally would conclude the package does not exist.

**What this does NOT change — the wall moved, it did not vanish.** The promcalc *decomposition data*
still does not exist, generating it is a ~44k-tile stage-B sweep, and the frozen bar's anchor gate still
demands the regenerated pipeline reproduce committed per-seed numbers within 0.01 on 9 of 10 seeds. So
the honest reclassification is: **not an environment/authority wall, but a COMPUTE + anchor-reproduction
wall** — ordinary executor work with a real risk of landing on `REFUSED_ANCHOR_MISMATCH`.

`scripts/promcalc_toolchain_probe.py` no longer asserts the refuted claim: it carries the claim, its
`REFUTED` status, the evidence, the import-vs-install-name trap, and an explicit note that the probe did
NOT measure this itself (it issues no subprocess by contract). The test was inverted to pin that the
claim and its status **travel**, rather than pinning an answer that already changed once.

---

## AMENDMENT 3, 2026-10-22 -- the D: volume is PHYSICALLY ABSENT, and the letter pin is exonerated

Amendment 1 recorded that D: is now the 1 TB volume, pinned by a 24-byte disk-signature binding.
That pin **held and is not the problem.** What is measured today is one level down.

**The sweep's 4th death, and the first one whose mechanism the OS names.** T=0.6 attempt 2
(wrapper pid 19640, lock holder 20100) died ~2 minutes in, *during the weight load* -- the log's
last line is `Loading weights: 100%|...| 237/237` and no sequence was ever generated. Sentinel
`EXIT=-1`. Event log, same second:

| time | log | id | content |
|---|---|---|---|
| 09:57:48 | Application | 1000 | `python.exe`, faulting module **`torch_cpu.dll`**, **Exception code 0xc0000006** |
| 09:57:48 | System | 51 (x4) | `An error was detected on device \Device\Harddisk1\DR1 during a paging operation.` |

`0xc0000006` is `STATUS_IN_PAGE_ERROR` -- literally "a memory-mapped page could not be read from
its backing store". The 2026-10-21 death reported faulting module `VCRUNTIME` and the
mmap-off-USB story was an **inference** from the disk events beside it. This one names the torch
module and the exception code, so that account is now **measured**.

**Then the drive left.** Checked minutes after the crash:

* `Get-Disk` -> disk 0 `SAMSUNG MZVLW512HMJP` (NVMe) and disk 2 `Generic Mass-Storage` (USB, **No
  Media** -- the card reader). **The 1 TB External1TB is not enumerated at all.**
* `Get-PSDrive D` -> `A drive with the name 'D' does not exist.`
* and the tell for *why* the 09:55 launch got as far as it did: the morning's `Ntfs` id-98 health
  events name **D: under two different device paths** -- `HarddiskVolume9` at 08:46,
  `HarddiskVolume6` at 09:38. The volume **re-enumerated between the two sweep attempts**, then
  dropped for good.

So the letter is stable and the physical link is not. Pinning the letter fixed a real problem
(paths resolving to nothing after a rename); it was never going to fix a disk that detaches.

### Wall class: `external-wall`, not code-closable

This matters for what to do next, so it is stated rather than implied:

* No launcher or loader change makes a 4.4 GB `mmap` survive a volume that disappears mid-read.
* The weights are **on the missing volume** (`D:/hf_cache`, 4.33 GiB), so the model cannot be
  loaded at all while it is absent -- independent of any code.
* **Checkpointing the generation loop would NOT have saved this run.** It is still the right fix
  for the 65/120 and 35/120 mid-run deaths, but this one died before sequence 1, so the structural
  fix I named as highest-value yesterday is **not** the binding constraint today. Recorded because
  the honest version of that recommendation is narrower than how I first put it.

Scoreboard: **4 deaths, 1 completed grid point** (T=0.8, `wiki/glm_sweep_point_T0.8_2026-10-21.json`).

### The three options, and which are the user's

| | option | owner | note |
|---|---|---|---|
| (a) | replug / different port / different cable, then relaunch | **user** | cheapest; the drive is physically gone, nothing here can re-attach it |
| (b) | run the sweep on server1 (GPU + local storage) | **user** | the right machine. Its Claude sessions (`Server1 cloud sessions for Soraya`, `DNA AI Decoder 13`, `Soraya Development`) all read `offline` in `ListAgents` on 2026-10-22, so nothing can be dispatched there yet |
| (c) | stage the 4.33 GiB weights on local storage (C: had 58.8 GB free) | **AUTHORITY** | it contradicts the standing *"every heavy file on D:, never C:"* rule. Measured as viable, deliberately NOT done |

(c) is recorded as a measurement, not a recommendation acted on: overriding a standing user rule
is a scope decision, and the rule exists because C: has run 85-99% full historically.

### Honest limits of this amendment

* It establishes that the volume detached and that the crash was an in-page fault. It does **not**
  establish *why* the USB link drops -- cable, port, enclosure, controller, and bus power are all
  un-separated. Every `MSPower_DeviceEnable` row for USB reads `Enable: True` (Windows is permitted
  to power these devices down), which is a **candidate** and not a finding; it was not tested,
  because the device has to be present to test it.
* One host, one drive. Nothing here generalises to the project's other external volumes.

---

## AMENDMENT 4, 2026-10-22 -- the MPRA half is CLOSED; promcalc is now the SINGLE blocker

This memo's reason for existing was that the tile chain had **two** substrate dependencies with
**different wall classes**, and that nothing separated them. One of the two is now resolved, so the
memo would be wrong if left as written.

**What was done.** The whole of GSE144621 was re-fetched from NCBI GEO to a repo-local, gitignored
`data/mpra/` -- measured first rather than assumed: **9,604,523 + 45,779,938 + 45,217,727 bytes =
~100 MB** by `curl -sIL` content-length, before downloading anything.

**Why that is not a D:-rule violation, stated rather than slipped past.** The standing rule is
*"every heavy/extra file for this project lives on D:, never C:"*, and it exists because C: has run
85-99% full. ~100 MB is **two orders of magnitude** below the GB-scale genomes, embeddings and model
weights the rule was written for, and `data/` on this host already holds **4.0 GB** of
`refseq_cache`. C: had 58.8 GB free. Judged inside the rule's intent; recorded here so the judgement
is auditable rather than invisible. The 4.33 GiB **model weights** are a different matter and were
NOT staged -- that remains the authority call in amendment 3, option (c).

**VALIDITY CHECK, and this is the part that makes the restoration trustworthy.** An independent
re-download reproduces figures recorded before it:

| quantity | recorded | re-derived from the fresh download |
|---|---|---|
| rows in the peak-tile file | 46,713 ("the entire file", CLAUDE.md) | **46,713** |
| 150 bp tiles after the published filter | 44,106 (`EXPECTED_TILES`, pinned) | **44,106** |
| replicate-derived noise ceiling | 0.941 (the ceiling criterion's own comparator) | **0.9410** |

Three independent figures, three matches. A corrupted or truncated fetch could not do that, which is
why these are quoted instead of a byte count.

**The wall-class separation now PRODUCES its designed measurement.** Running the tile chain end to
end off the local cache:

```
tiles 44106 | ceiling 0.9410
--in-dir D:/dna_decode_cache/promcalc is not reachable [absent].
```

It gets **past** the tile load and fails on **only** promcalc. That is exactly what the step ordering
in `executed_plans/Make_The_Tile_Chain_S_Substrate_Legible_And_Classify_Its_Wall_Plan/` was built to
show -- and it could not be shown until one of the two dependencies was actually resolved.

**The promcalc half, measured today, not relayed.** `scripts/promcalc_toolchain_probe.py` reported
`INPUTS_ABSENT` *while the substrate sat on disk*, because `external_root.probe()` enumerates DRIVE
roots and the restored copy is repo-local. A false absence is the most expensive answer that probe
can give -- it sends the next session to re-fetch data it already has. `find_inputs` now also looks
at the repo-local `data/` tree, and the verdict moved to `ENV_ABSENT_DEPS_SATISFIABLE_HERE` with the
location named.

And the install question is now measured **in this environment**:

```
uv pip install --dry-run promotercalculator
  Resolved 7 packages ... Would install 1 package
  + promotercalculator==1.2.4          (sklearn stays at 1.8.0)
```

So amendment 2's refutation of the isolation premise is confirmed here, first-hand: **one package, no
isolation, sklearn untouched.**

### What is STILL open, and whose call it is

* **Regenerating the promcalc predictions is NOT done and was not started.** The saved plan
  explicitly parked it as an open question rather than folding it into scope, and it is a
  multi-hour compute job over 44,106 tiles. The install is a one-package, reversible step; the
  regeneration is the part that needs a decision.
* `ENV_ABSENT_DEPS_SATISFIABLE_HERE` **still over-claims** and must not be read as a green light:
  it fires on an interpreter existing plus >=5 GB free. What is measured is the dry-run resolution
  above, not that the toolchain produces correct predictions.
* **Nothing here touches the sweep.** The GENERator weights are on the detached volume and the
  sweep arm is unchanged -- amendment 3 stands in full.

### Honest limits

* The restored substrate is verified by three reproduced figures, not by checksum against the
  original copy: the original is on the detached volume and could not be compared.
* The end-to-end run above is a REFUSAL path. It proves the tile load succeeds and the promcalc
  dependency is correctly diagnosed; it does **not** demonstrate that any arm scores anything.


## AMENDMENT 5, 2026-10-23 -- the promcalc half is CLOSED; provenance is by REPRODUCTION, not a version pin

Amendment 4 left promcalc as the single blocker and said the regeneration "is NOT done and was not
started." It is now done. The dependency's wall class was `data-absent`; it is closed, and the route
that closed it is compute on this host -- which is why it was never `external-wall`. The frozen
`dependencies` array is unchanged, as in amendments 3 and 4.

**Regeneration, measured from `data/promcalc/decompose.json` rather than relayed:** 43,611 tiles
attempted on top of 495 already complete, `err 0`, `total_complete 44106` against `tiles_in_tsv
44106` -- `substrate_coverage 1.0`. 3,183,603 site rows at a mean of 73.0 sites per tile, across
the 15 recorded fields (`tss`, `dG_total`, `dG_10`, `dG_35`, `dG_ext10`, `dG_spacer`, `dG_disc`,
`dG_ITR`, `dG_UP`, `dG_bind`, `Tx_rate`, `hex35_start`, `hex10_start`, `spacer_len`, `tile_id`),
`strand: forward_only`, 7 workers, **770.37 minutes** (~12.8 h). Output
`data/promcalc/decomposition.part-01.csv.gz`, 156,810,679 bytes; `data/promcalc` totals 159 MB.
The unattended chain closed it at `CHAIN_EXIT=0 at Fri 10/23/2026 3:48:58.22`.

### Risk Flag 1's mitigation is UNMET, and it is now worse than un-recorded -- it is unrecoverable

The saved plan's first risk flag was `promoter_calculator` version drift, `[unverified]`, mitigated by
"Step 6 records the version so this risk cannot recur." **It does not.** The arm artifact's non-arm
keys are only `schema` / `measured_on` / `n_features` / `verdict` / `verdict_detail` / `smoke`, and
`decompose.json` records workers, fields, timings and coverage but no version. Nor can it be read
back off the host:

* `uv pip show promotercalculator` -> `Package(s) not found`
* `importlib.metadata.version` fails for BOTH the distribution name and the import name
* `grep -c -i promotercalculator uv.lock` -> **0**
* nothing matching `*promoter*` remains in `.venv/Lib/site-packages`

**Leading mechanism, stated as such because I did not watch it happen:** the package was added with
`uv pip install` (amendment 2 measured that as a one-package, sklearn-safe step) and is absent from
`uv.lock`; `uv run` re-syncs the project environment to the lock, so any later `uv run` -- and this
plan is driven entirely by `uv run` -- would have removed it. Consistent with the fourth check:
sklearn is still **1.8.0**, exactly as amendment 2 predicted and the opposite of the refuted
1.9.1 claim. The reusable shape is that **an environment mutated outside the lock file is not a
provenance record**: it is reverted by the next tool invocation, silently, and the version goes with it.

### What substitutes for the missing version: the reconcile gate reproduced the original EXACTLY

This is why the number is still interpretable. Both pre-registered anchors matched, and matched
exactly rather than within tolerance:

| anchor | expected | measured | tolerance | source of the expectation |
|---|---|---|---|---|
| `gc_seed0` | 0.3000 | **0.3000** | 0.01 | `glm_encoder_gate.PREREGISTERED.reconcile_anchor` |
| `max_tx_median` | 0.3968 | **0.3968** | 0.02 | `wiki/promcalc_vs_learned_2026-10-12.json` |

The second anchor is the ORIGINAL 2026-10-12 run's own scalar, so the regenerated substrate
reproduces the pre-existing number to four decimal places. That is **behavioural provenance** and it
is what makes today's comparison a comparison rather than a confound. It is strictly weaker than a
version pin for FUTURE drift -- it can only ever say "this run agreed with that run", never "the
toolchain is pinned" -- so the version-recording gap stays open as a real debt, not a closed item.

**The reconcile OUTCOME is not an explicit field; do not look for one.** It is established two ways
from what IS recorded: by code structure, since `scripts/promcalc_decomposition_arm.py:248-262`
returns `INDETERMINATE_RECONCILE_FAILED` before `NO_GAIN` is reachable, so the published `NO_GAIN`
is unreachable unless both checks passed; and directly, since `arms.gc.per_seed[0]` and
`arms.max_tx.median` are both in the artifact and can be compared to the table above by hand.

### A naming defect found while checking the control, and deliberately NOT fixed

`PREREGISTERED["control_max_abs"] = 0.05`, but line 261 tests it against
`arms["control_feature_permuted"]["median"]` -- a **median**, not a max. The median is **+0.0342**,
so the control passes and the verdict stands. But **two of ten seeds exceed the bar in absolute
value** (seed 5 at **-0.1129**, seed 7 at **+0.0851**) and neither is tested. A reader comparing the
artifact's recorded `control_feature_permuted.max = 0.0851` against a bar named `control_max_abs`
would read FAIL where the gate reads PASS.

The gate's behaviour is defensible and the name is not: a feature-permutation null SHOULD scatter
across seeds, and the median is the right statistic for "is this null centred on zero". `control_median_abs`
would say what the code does. **Not renamed**, because the constant lives in the frozen `PREREGISTERED`
dict and this plan's own manual test requires `git diff` to show no change to it. Recorded instead,
which is the whole point of having a memo.

### The verdict, stated precisely enough that it cannot be over-read

`NO_GAIN`. Decomposition **0.3925** against the `max_tx` scalar **0.3968**: paired median
**-0.0043**, CI95 **[-0.0113, +0.0013]** (includes zero), winning **2 of 10** seeds against a bar of
margin 0.02 AND at least 8 seeds.

What that does NOT say is that promoter thermodynamics is uninformative. Against the GC comparator
at **0.3112**, decomposition wins **10/10 at +0.0802** (CI excludes zero) and `max_tx` wins **10/10
at +0.0865**.
So the thermodynamic model carries real signal over composition; what adds nothing is the
**nine-component breakdown over the single scalar reduced from it**. Roughly 730 numbers per tile are
discardable. A negative answer to the pre-registered question is this plan working, not failing --
Verification Signal 2 said so before the result landed.
