# The L2 doubt layer — measured, wired, and registered (F-A, 2026-08-31)

**One sentence:** across the entire deployed AMR surface — 1,818 genomes, 1,279 determinant families
the rules cannot represent, six drugs — exactly **one** family survives a family-wise correction, and
it is the independently confirmed gap.

Steps 2–4 of F-A (`plans/Hybrid_Decoder_Architecture_Plan.md`). Step 1 (the screen) shipped earlier
the same day. Frozen AMR surface byte-unchanged; prospective lock re-verified.

---

## What the full-index run changed

The step-1 run was capped at 220 genomes/drug. At full index (8.3×) the raw `rmt_like` signature
(>=3 R carriers, 0 S) fires on **five** families, not one:

| drug | family | R/S | verdict |
|---|---|---|---|
| gentamicin | `rmtE1` | 36 / 0 | the known gap |
| ciprofloxacin | `qnrA1` | 4 / 0 | noise |
| ciprofloxacin | `oqxA10` | 3 / 0 | noise |
| ceftriaxone | `blaTEM` | 9 / 0 | noise |
| ceftriaxone | `blaOXA` | 6 / 0 | noise |

At n=4 against a cohort base susceptible-rate of 0.583, "all four resistant" happens by chance with
p = 0.030 — against ~125 families screened for that drug. **Purity alone cannot tell a gap from a
coincidence**, and shipping those five side by side would have made the layer noise.

## Two nulls, and the obvious one is wrong

The first correction attempt scored each family by a lower-tail binomial on its observed susceptible
count. It **called `aph(6)-Id` (62R/28S) STRONG at p ≈ 5e-5** — and `aph(6)-Id` is a *correct*
exclusion, a streptomycin determinant that travels with gentamicin resistance by linkage.

That is not a tuning problem. **Every co-occurring determinant is R-enriched**, so an enrichment null
floods the layer with correct exclusions — the same way raw volume ranking buried the answer in step 1.
The `rmt` signature is **purity**: zero susceptible carriers among a well-powered set. One susceptible
carrier is positive evidence the exclusion is deliberate, so it *ends* the signal rather than weakening
it. Pinned by `test_an_enrichment_null_WOULD_have_called_a_correct_exclusion_strong`.

## The result

`P(zero susceptible carriers | n, cohort base rate)`, Bonferroni-corrected over the families screened
for that drug:

| drug | uncounted families | raw signature | **STRONG** | weak | base S-rate |
|---|---:|---:|---:|---:|---:|
| ceftriaxone | 216 | 2 | **0** | 0 | 0.202 |
| ciprofloxacin | 125 | 2 | **0** | 1 | 0.583 |
| **gentamicin** | 131 | 1 | **1** — `rmtE1`, p = 4.11e-12 | 0 | 0.517 |
| meropenem | 317 | 0 | **0** | 1 | 0.870 |
| oxacillin | 401 | — | no labels for this drug | | |
| tetracycline | 89 | 0 | **0** | 0 | 0.281 |

The correction drops 4 of the 5 raw hits and keeps the true one.

**Honest limits.** The AMR arm has exactly **one** independently confirmed gap, so "recovered 1 of 1"
is a single case and **not a rate** — it bounds nothing about gaps never confirmed. Ranking depends on
labels, and only 200 of 1,818 cached genomes carry NCBI-PD calls. Oxacillin is uninformative here.
The two arms below are **not pooled**: one is per-determinant-family, the other per-isolate.

## The position-novelty arm (reproduced, not recomputed)

The incumbent, from its committed artifact: **median sensitivity 0.604** on the EFV catalog-negative
blind spot, lift 4.69. Per drug: doravirine 0.714 · etravirine 0.688 · efavirenz 0.604 ·
rilpivirine 0.571 · nevirapine 0.515.

## Wired into the record

All six target-site paths funnel through one seam (`_target_site_record`), so the block lands once.
Verified on the real CLI:

```
dna-amr --drug efavirenz --observed RT:K103N   ->  R  | doubt: none
dna-amr --drug efavirenz --observed RT:K103R   ->  S  | doubt: WEAK
dna-amr --drug lamivudine --observed RT:M184V  ->  R  | doubt: none, applicable=False
```

The middle case is the point: the prediction stays **S** (the rule is untouched) while the record now
says that S is least trustworthy. That is the 53-isolate HIV blind spot made visible at call time.

**Three states, never collapsed.** `not-applicable` (a position-based catalog — every substitution at
a catalogued position is *already* called, so the flag could never fire) · `not-assessable` (the genome
path does not surface observed substitutions) · `assessed`. Reporting "no doubt" for either of the
first two would be a false clean bill — the exact failure this layer exists to surface.

**The constraint is enforced, not documented.** `DoubtBlock.as_dict()` runs `assert_no_call` on its own
output and raises rather than emit anything call-shaped — recursively, keys and values. A doubt signal
may qualify a call and explain itself; it may never overrule L1 or emit one. That is what keeps L2 out
of the learned-predictor regime that is 0-for-5 de-confounded.

## Registered augment-only

`trust_surface.doubt_layer_for(drug)` attaches under its own `doubt_layer` key. "Does this cell have a
known completeness gap?" and "how well is this cell validated?" are different questions; merging them
is the shared-key silent-overwrite trap. The guard compares each badge computed with the layer on
versus off and requires every pre-existing field to be byte-identical — with a non-vacuity test, since
a layer that attached nothing would make that guard prove nothing.

## Artifacts

`scripts/doubt_layer_per_cell.py` → `wiki/doubt_layer_per_cell_2026-08-31.json` ·
`dna_decode/eval/doubt.py` · 30 tests across `test_doubt_layer.py`,
`test_doubt_record_guard.py`, `test_doubt_trust_surface_augment_only.py`. Suite: 4029 passed, 0 failed.
