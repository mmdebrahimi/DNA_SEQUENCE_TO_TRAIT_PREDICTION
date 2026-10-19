# F3 de-confounder calibration on 1001 Genomes — `VALID_DIFFERENTIAL_RESPONSE`

**2026-10-10.** Closes `blocked:positive_control_not_fold_compatible`, open since 2026-10-09.
Artifacts: `wiki/f3_deconfound_calibration_1001g_2026-10-10.json`,
`wiki/f3_ibs_crosscheck_1001g_2026-10-10.json`, `wiki/f3_rung_census_1001g_2026-10-10.json`.
Reproduce: `uv run python scripts/f3_calibrate_1001g.py --n-markers 1500 --n-perm 100 --max-traits 8 --coarse-percentile 10`.

## What was blocked, and what unblocked it

`cv_ridge_gp`/`cv_model_gp` gained a keyword-only `groups=` on 2026-10-09 (bit-for-bit no-op when
`None`), but a grouped number from it was **not yet evidence of de-confounding** — nothing had shown the
instrument responds to structure when structure is present and stays quiet when it is absent. The
blocker was the positive control: `dgrp` had no rung both fold-compatible and not mostly singletons.

The substrate that works is the **1001 Genomes** Arabidopsis panel (1,135 imputed natural accessions,
already on `D:`) with **AraPheno** phenotypes fetched and ranked by overlap
(`scripts/fetch_arapheno_traits.py` → `wiki/arapheno_trait_overlap_2026-10-10.json`; 536 phenotypes,
0 fetch failures, 163 usable at ≥250 overlap).

## Result

Controls are selected **by construction** from each panel's `design` field, never by a structure
statistic — `pc1` was measured and rejected for this because a bi-parental cross *maximises* LD and
would rank as the most structured panel, which is precisely the negative control.

| arm | design | rung | clusters | median drop | n dropped | sign-test p |
|---|---|---|---|---|---|---|
| 1001G (positive) | natural population | fine p=5 | 157 | **+0.0441** | 7/8 | 0.0703 |
| 1001G (positive) | natural population | coarse p=10 | 90 | **+0.0334** | 7/8 | 0.0703 |
| bloom (negative) | constructed bi-parental | fine p=2 | 170 | **0.0000** | 4/8 | 1.000 |
| bloom (negative) | constructed bi-parental | coarse p=5 | 74 | **−0.0052** | 1/8 | 0.0703 |

**Symmetric differential (positive − negative, each arm at the same granularity): +0.0441 (fine),
+0.0386 (coarse).** Every arm's grouped fit still beats its own grouped permutation null, which is the
frozen criterion's other half.

The negative control is the part worth looking at: median drop **0.0000** with four traits up and four
down. Grouped CV has no systematic effect on a cross with no population structure — which is what makes
the positive arm's +0.0441 attributable to structure rather than to the act of grouping.

Per-trait, the positive arm's drop is largest exactly where the apparent signal was weakest. FT10 and
DTF1 (n≈1,000, r≈0.70) drop 0.024/0.034 and still beat their grouped nulls — real trait signal survives
de-confounding. Raffinose and inter-specific pollination drop 0.128/0.102 and **fail** their grouped
nulls — on those traits the random-CV signal was essentially all structure.

## Two things that cut against the convenient reading

**1. The sign test does not reach 0.05, and the verdict does not rest on it.** 7 of 8 is two-sided
p = 0.0703. The frozen criterion is *differential* — the negative control keeps beating its own grouped
null AND the positive's drop exceeds the negative's — and that is met. **Do not quote this arm as
significant at 0.05.** With 8 traits the floor for a two-sided sign test is 0.0078 (8/8), so 7/8 cannot
clear 0.05 at this n; more traits, not a different test, is what would change it.

**2. The bracket answered the opposite of what it was built to check — and refuted a reading I formed
mid-run.** `select_threshold` takes the *smallest* passing rung, which under-merges, so the worry was
that 157 clusters over 1,135 individuals makes grouped CV nearly random CV and *understates* the
de-confounding. Coarsening to 90 clusters **reduces** the median drop (+0.0441 → +0.0334; paired median
delta −0.0156, only 3 of 8 traits dropping more). So the under-merging bias is **not** masking a larger
effect on this panel. While the run was in flight I had seen only the first three cells — FT10, DTF1 and
rSurvival, all of which do drop more at the coarse rung — and wrote that coarsening increased the drop
"consistently so far". Those three are exactly the 3 of 8 exceptions; the other five reverse it.

A third oddity, reported rather than explained: at the coarse rung the negative control's drop is
mildly **negative** (−0.0052, 7 of 8 traits going that way), i.e. grouped CV slightly *outperforms*
random CV on a structureless cross. That is why the criterion is differential rather than "the negative
must be zero".

## Three defects found and fixed while getting this number

**The coarse bracket endpoint was mis-specified, and the refusal that now catches it did not exist.**
`sweep_thresholds` evaluates the one-fold share on the **full 1,135-accession panel**, but grouped CV
folds each trait's **non-missing subset**. The coarsest rung the sweep passes (p=20, largest 0.1947 on
the full panel) reaches 0.2104 on FT10's 1,003 rows and is **inadmissible for 6 of 8 traits**. The first
coarse run therefore scored only the 2 traits that happened to clear it — and those were the
*lowest-n, largest-drop* ones, so an unpaired median would have credited the coarse rung for
**attrition**. `wiki/f3_rung_census_1001g_2026-10-10.json` measures the admissible set per trait: only
**p=1, 2, 5, 10** are admissible for all 8, making p=10 the correct coarse end. The run now refuses an
inadmissible endpoint and names the traits; the refusal is verified non-vacuous (default refuses and
lists all 6, `--coarse-percentile 10` passes 8/8). Fine-vs-coarse is additionally compared **paired** on
traits scoreable at both rungs.

**Marker budget was asymmetric by an order of magnitude.** The committed `run_arm` loads bloom at
stride 1 (11,623 markers) while the positive arm ran on 1,500. Striding bloom to a comparable budget
moved the negative control's drop from 0.0164 to ≈0 — a validity fix, not a tuning one.

**Percentile matching would have compared different operations.** The same ladder rung buys very
different fold geometry on two differently-structured panels: on 1001G p=10 gives 90 clusters over
1,135, while on bloom p=20 gives 22 over 1,008 — a 12.6-member cluster against a 46-member one. The
arms are matched on **clusters per individual**, a rule fixed before any bloom drop was read
(`matched_rung`), which selects bloom p=2 (mean size 5.93 vs the positive's 7.23) and p=5 (13.62 vs
12.61). Noted in passing: the negative control's admissible range is **wide** (p=1..20) and the
positive's **narrow** (p=5..10), because a cross forms no dominant cluster — the width of the admissible
range is itself a structure signal, observed here and deliberately not used as a test.

## Honest limits

- **Bounds the INSTRUMENT, not any downstream claim.** A calibrated de-confounder does not validate any
  organism-level prediction, and this says nothing about whether a learned representation helps.
- **One positive panel, one negative panel.** Two panels do not establish that the calibration transfers
  to another substrate; `eval/regime.py::DECONFOUNDER_CALIBRATION` should cite this, not generalise it.
- **The two arms have different traits**, so the between-arm comparison is two-sample; only the
  within-arm sign tests are paired. Cross-arm pairing is impossible by construction (Arabidopsis vs
  yeast).
- **Trait independence is the best this substrate allows, not good.** Traits are selected **one per
  study** because 163 usable AraPheno traits span only 8 studies and 121 of those come from a single
  rainfall factorial — a top-N-by-overlap selection would have handed the sign test six correlated
  replicates of one experiment. Even so, **FT10 and DTF1 are both flowering time** and are not
  independent observations.
- **Genotypes are IMPUTED**, which can flatten or sharpen the very relatedness structure being measured.
- **Markers are a MAF≥0.05 subsample** read from chunk-aligned blocks (9,442 passing of 60,000 read),
  not an LD-pruned genome-wide set; bloom is 1,661 markers to the positive's 1,500, so "comparable
  budget" is comparable, not identical.
- `n_perm=100` bounds the resolution of every `beats_null` decision; `MIN_TRAITS_FOR_SIGN_TEST=5` is
  asserted, not derived.
- Clustering is greedy-representative (chaining-resistant) on this repo's own IBS. The project's own
  kinship (IBS, MAC≥5, all markers) gives a different partition — same-pair Jaccard **0.2532** — and the
  calibration deliberately uses the repo IBS, so the verdict does not depend on the converted kinship.
- `r_random` is reused across rungs for the same (panel, trait) because the ungrouped fit cannot see the
  grouping. That is verified rather than assumed — FT10 0.692, DTF1 0.7046, rSurvival 0.2412 and
  Raffinose 0.1817 are byte-identical at p=5 and p=10, and a reuse reproduces its source exactly.
