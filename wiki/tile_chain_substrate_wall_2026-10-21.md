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
