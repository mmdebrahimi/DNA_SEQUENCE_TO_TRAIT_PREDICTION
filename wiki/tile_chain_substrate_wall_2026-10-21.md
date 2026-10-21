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
