# CABBAGE scored cell — *Campylobacter jejuni* × ciprofloxacin, N=4,265

**2026-10-10.** First cell scored through the EMBL-EBI AMR Portal's **pre-computed** AMRFinderPlus,
so the genotype side cost ~1 KB/genome instead of ~95 s/strain of local Docker. Artifacts:
`wiki/cabbage_scored_cell_campylobacter_ciprofloxacin_2026-10-10.json`,
`wiki/cabbage_fn_anatomy_campylobacter_ciprofloxacin_2026-10-10.json`.
Reproduce: `uv run python scripts/cabbage_score_cell.py --max-isolates 5000` then
`uv run python scripts/cabbage_fn_anatomy.py`.

## Result

Frozen `eval/amr_rules.call_resistance` (byte-unchanged), labels from portal MIC through this repo's
own `external_mic_labels.tier_for_isolate` on **STRICT tiers only**, operator-aware.

| | |
|---|---|
| strict-labelled | 4,712 |
| excluded for leakage | 60 |
| scoreable | 4,652 |
| **scored** | **4,265** |
| confusion | tp 997 / fn 21 / tn 3,236 / fp 11 |
| **accuracy** | **0.9925** |
| **sensitivity** | **0.9794** [0.969, 0.986] |
| **specificity** | **0.9966** [0.994, 0.998] |
| fetch failures | 387 (8.3%), excluded and counted, never scored wrong |

Funnel and failure count reproduced identically across two independent runs. An earlier
`--max-isolates 1500` pass scored 1,444 at sens 0.9971 (1 FN) — the **full set is harder**, so the
subset was optimistic and running the whole thing mattered.

## Both error classes diagnosed

**All 21 false negatives carry ZERO counted determinants.** That is the rule's own declared blind spot,
not decay: `call_resistance` emits `undetectable_mechanisms: [efflux, porin_loss, regulatory]` and
states that an S call cannot rule them out. Same shape as the single FN in the 2026-09-28 prospective
Campylobacter cell — abstention-by-design leaking into a hard S call.

**All 11 false positives carry a determinant, and the cause is located.** Per-determinant PPV over the
whole scored set:

| determinant | carriers | truth R | truth S | PPV | Wilson CI |
|---|---|---|---|---|---|
| `gyrA_T86I` | 1,000 | 995 | 5 | 0.995 | [0.988, 0.998] |
| **`gyrA_P104S`** | **6** | **0** | **6** | **0.000** | **[0.0, 0.39]** |
| `gyrA_D90N` | 1 | 1 | 0 | 1.000 | [0.207, 1.0] |
| `gyrA_T86V` | 1 | 1 | 0 | 1.000 | [0.207, 1.0] |
| `gyrA_T86K` | 1 | 1 | 0 | 1.000 | [0.207, 1.0] |
| `gyrA_D90Y` | 1 | 1 | 0 | 1.000 | [0.207, 1.0] |

`gyrA_P104S` appears nowhere among the 997 true positives and accounts for 6 of the 11 false positives;
the other 5 FP are `T86I` carriers that are phenotypically S, consistent with that determinant's 0.995
PPV. **This is the MIRROR of the gentamicin `rmt` case** — there a determinant the rule could not *see*
produced false negatives; here one the rule *counts* produces false positives.

**Attribution, verified not inferred:** on `GCA_004836905.1`, `gyrA_P104S` is the **sole** counted
determinant (`n_determinants=1`, `counter='qrdr_point'`, `threshold=1`) and **AMRFinderPlus itself**
files it as `Class=QUINOLONE` / `Subclass=QUINOLONE`, named "Campylobacter quinolone resistant GyrA".
The frozen rule counts QUINOLONE-class gyrA point mutations, so it is behaving exactly as specified —
the resistance-conferring claim for P104S originates **upstream, in AMRFinder's catalog**. A fix
belongs there, or in an L2 doubt flag, not in the counting logic. Changing what the rule counts would
edit the sha256-pinned frozen surface and retire the active prospective lock, so it is a **user
authority call and is deliberately not made here**.

## The probe's own first answer was wrong, and the guard that caught it is now permanent

The first run of the anatomy probe read `res["hits"]`, a key this record does not have, so **every
bucket came back determinant-free** — and it printed a confident
`ALL_FN_CARRY_NO_DETERMINANT_DECLARED_BLIND_SPOT`, which is a *plausible* verdict that happens to match
the prospective cell's story. The tell was the TP bucket: 997 true positives with no counted
determinant is impossible, since a TP has one by definition. The probe now **refuses a verdict (exit 3)
when the TP bucket is determinant-free**, that being a plumbing signature rather than a result. The FN
conclusion survived the fix; the FP finding only became visible after it.

## Honest limits

- **GEOGRAPHICALLY CONCENTRATED — do not quote this as source-diverse.** Against this repo's own bar
  (`scripts/source_diverse_validate.py`: ≥5 sources AND largest ≤0.60) the cohort **fails on country:
  96.4% United States of America** (5 countries). It passes on `database` label (7 distinct, largest
  NARMS 0.452), but that column is an overlapping **aggregator** label, not a BioProject. Verdict in the
  artifact: `source_concentrated`. A cohort that fails this bar cannot reveal a determinant-family blind
  spot it does not contain — the mechanism that made the gentamicin `rmt` gap structurally invisible.
- **NOT clonality-corrected (G8 open).** These are isolate-level counts and *Campylobacter* surveillance
  collections are clonally structured. The repo's canonical lineage route is Mash via Docker, and
  **Docker's daemon is down on this host with no native `mash` binary**; the portal metadata carries no
  MLST/ST field, so the Mash-free lineage unit used for TB (barcode) and E. coli (ST) is unavailable
  here. Computing MLST locally would mean fetching ~4,652 assemblies (~9 GB) and running BLAST per
  isolate. Recorded as blocked rather than approximated with an unsound proxy.
- **Mechanistically expected, so not evidence of broad catalog coverage.** Cohort R is overwhelmingly a
  single well-catalogued substitution (`gyrA_T86I`, 1,000 of 1,006 counted determinants), so
  near-perfect numbers on this cell test one point mutation, not the catalog.
- **Per-determinant PPV is MARGINAL, not causal** — a carrier may hold several determinants, so credit
  is not attributable to one symbol from these counts alone. And `0 of 6` is suggestive with an upper CI
  of 0.39; it does **not** demonstrate P104S is benign.
- **Leakage exclusion is accession-string based** (version-insensitive) via the fail-closed
  `cohort_manifest` (1,509 prior accessions, manifest complete). It does not resolve BioSample identity,
  so an isolate deposited under two assemblies could survive.
- **387 fetch failures are uncensused by organism.** Whether the portal's `genomes/` tree covers all
  CABBAGE assemblies remains the registered constraint `cabbage-genomes-tree-coverage-uncensused`; this
  run measures 8.3% absent for *this* species only, which is a datapoint toward it, not a census.
- Collection years span 2007–2021 (15 distinct), so the cohort is **not** temporally locked and these
  isolates are **not** prospective relative to the active 2026-08-31 lock.
