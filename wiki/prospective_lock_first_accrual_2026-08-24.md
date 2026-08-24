# Prospective-lock: FIRST ACCRUAL SCORED (2026-08-24)

**The reproducibility freeze's forward-path #2 produced its first numbers.** The prospective lock
(`wiki/prospective_lock_manifest_2026-06-22.json`, `lock_date=2026-06-13`) was built so that any isolate
becoming public *strictly after* the lock is a leakage-free test case **by construction** — the frozen
decoder cannot have been tuned to data that did not yet exist. Until now every accrual sweep returned
`eligible=0` (a genuine ACCRUING zero, last checked 2026-07-10). This is the first non-zero one.

`verify_lock` passed on every run: the decoder is **byte-identical to the 2026-06-13 lock** (5 pinned
files), so these scores are of the frozen surface, unchanged.

## The cohort

`scripts/fetch_prospective_cohort.py --source ncbi_pd`, run twice (deterministic — the funnel reproduced
exactly):

| stage | n |
|---|---|
| NCBI-PD metadata rows | 918,436 |
| with measured AST + downloadable assembly | 13,859 |
| SRA-dated pre-lock (skipped) | 9,740 |
| **prospective-eligible (strictly post-lock)** | **63 isolates / 215 rows** |

`overall_status=OK` — every scope queried cleanly, so this is a real accrual, not an outage reported as
data. All 63 isolates resolved to `Escherichia_coli_Shigella`.

## Results

| cell | N | R/S | acc | sens | spec | TP/FP/TN/FN |
|---|---|---|---|---|---|---|
| E. coli x **ciprofloxacin** | 61 | 24/37 | **0.967** | 0.917 | **1.000** | 22/0/37/2 |
| E. coli x **gentamicin** | 62 | 49/13 | **0.532** | **0.429** | 0.923 | 21/1/12/28 |

Zero abstentions in both. Ceftriaxone (49R/**1S**) and tetracycline (42R/**0S**) accrued too but are
single-class-starved and were **not** scored — reporting a spec on 1 S isolate would be theatre.

Honest tier: a **temporal** prospective stress test of the frozen decoder. Leakage-free by construction,
but NOT lineage-independent and NOT clinical validation. N is small and will grow as the cohort accrues.

## The gentamicin result is a located catalog gap, NOT decoder decay

sens 0.429 with spec 0.923 means the rule **under-calls** — it misses determinants rather than
mislabelling. Diagnosed on the cached AMRFinder outputs:

* **24 of the 28 false negatives carry a 16S rRNA methyltransferase** (`rmtE1`, `rmtE`, `armA`-family).
* **0 of the 28 carry `aac(3)`** — the aminoglycoside-modifying family the rule does count.

The mechanism is exact. The frozen rule is `subclass_any=['GENTAMICIN'], threshold=1`, i.e. a literal
match on AMRFinder's `Subclass`. In these runs:

| determinant | AMRFinder Subclass | counted? |
|---|---|---|
| `aac(3)-VIa`, `aac(3)-IIe`, `aac(3)-IId` | `GENTAMICIN` | yes |
| `ant(2'')-Ia` | `GENTAMICIN/KANAMYCIN/TOBRAMYCIN` | yes |
| **`rmtE1`** | **`AMINOGLYCOSIDE`** (generic) | **NO** |
| `aph(3')-Ia` | `KANAMYCIN` | no (correct) |
| `aph(3'')-Ib`, `aph(6)-Id`, `aadA*` | `STREPTOMYCIN` | no (correct) |

16S rRNA methyltransferases confer high-level resistance to the 4,6-disubstituted deoxystreptamines —
gentamicin, tobramycin, amikacin. AMRFinder files them under the generic `AMINOGLYCOSIDE` subclass with
no drug named, so a `GENTAMICIN`-substring rule cannot see them.

**This is not drift.** The rule was always blind to `rmt`; it simply never mattered, because the frozen
in-distribution cohort (N=128, acc 0.945 / sens 0.893) evidently contained no `rmt` carriers. The
post-lock accrual is enriched for them. The rule's exclusion of `aph`/`aadA` is *correct and deliberate*
(streptomycin/kanamycin genes that do not confer gentamicin resistance) — the gap is specifically the
generic-subclass methyltransferases.

**Why the labels are probably not the problem:** a label artefact would scatter determinant profiles
randomly across the false negatives. Instead they are coherent — nearly every FN carries an `rmt` and
none carries an `aac(3)`. That is the signature of a missing rule, not a noisy label. (Contrast the
project's own high-sens/low-spec heuristic, where suspecting the label IS the right first move; this is
the mirror case.)

## What this does and does not license

* **NOT changed:** the frozen surface. `amr_rules.py` and `calibrated_amr_rules.json` are byte-identical
  to the lock and must stay so — editing them would invalidate the lock and the reproducibility freeze.
* **Recommendation for a future UNFROZEN revision:** count 16S rRNA methyltransferases toward
  aminoglycoside cells (they are pan-aminoglycoside), either by an explicit gene-family clause or by
  admitting the generic `AMINOGLYCOSIDE` subclass for gentamicin. Any such change needs its own
  validation and would begin a NEW lock — it cannot be retrofitted into this one.
* **The cipro cell holds prospectively** (0.967 / sens 0.917 / spec 1.000), which is the reassuring half
  of the same experiment and makes the gentamicin contrast interpretable rather than a global failure.

## Two pipeline defects fixed to get here

Both survived because every prior accrual returned zero rows — **a pipeline whose output has only ever
been empty has never had its output validated.**

1. **The cohort recorded no `organism`.** The sweep spans Campylobacter + E. coli/Shigella + Klebsiella
   and `call_resistance` applies a per-organism rule, so a mixed cohort would be silently mis-scored.
   *Honest scope:* all 63 rows turned out to be E. coli, so for THIS cohort the defect was latent — the
   default organism was right by luck. What was actually broken is that the artifact could not prove its
   own provenance.
2. **The scorer did not filter by drug.** The first scoring run reported `n_scored=215, R=164, S=51` —
   exactly the *all-drugs* label sum — because it compared a ciprofloxacin prediction against every row's
   own (ceftriaxone / gentamicin / tetracycline) label. That artifact was a meaningless number; it was
   deleted before being committed, and `_scope_to_drug` now refuses rather than scoring a mismatched set.

## Reproduce

```bash
uv run python -m scripts.fetch_prospective_cohort --source ncbi_pd        # ~50 min, network
uv run python -m scripts.prospective_lock_validate \
  --cohort-tsv <out>/prospective_cohort.tsv \
  --drug ciprofloxacin --organism Escherichia_coli_Shigella --amrfinder-organism Escherichia
```
Needs Docker (AMRFinder) for the first pass; subsequent drugs reuse the cached runs.
Artifacts: `wiki/prospective_lock_validation_Escherichia_coli_Shigella_{ciprofloxacin,gentamicin}_2026-08-24.json`.
