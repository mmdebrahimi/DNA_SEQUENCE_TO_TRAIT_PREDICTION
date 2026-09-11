# Acquired determinants ARE less clade-coupled than chromosomal ones — measured directly

**Verdict: `SUPPORTED`** on all four pre-registered bands, and on a ceiling-normalized secondary band
the frozen bar did not anticipate. Bar frozen before any score existed:
`wiki/determinant_clade_coupling_acceptance_bar.json`. Run:
`uv run python scripts/determinant_clade_coupling.py` (read-only, offline, no Docker, no model).

This is the direct test of `COUPLING-not-popdesign`, the last surviving mechanism claim from the
2026-09-10 `/innovate` sweep. It is the **first survivor of that sweep to survive a *stronger* test
rather than fall** — the other two that were acted on (`POOLING-dilution`, `DOUBT-bacterial-gap`) were
both refuted once their claim was tested rather than their proxy.

## The claim, and why the kill-test that let it survive was too weak

Horizontal transfer is nature's own randomization: a mobile determinant is decoupled from clade the way
meiosis decouples a locus in a cross. So **acquired** determinants should be less clade-coupled than
**chromosomal** point mutations.

The kill-test *inferred* coupling from a proxy — a sensitivity drop under lineage collapse — across 10
cells with a single drug in the chromosomal arm. This measures the named mechanism itself: how
concentrated each determinant family's carriers are across lineages, with **both arms drawn from the
same 621 genomes and the same lineage partition**.

## Result

621 E. coli genomes, 104 lineages (MLST, single `ecoli_achtman_4` scheme), 92 acquired families and 400
chromosomal families at ≥5 carriers.

| band | point | acquired | gap | perm max (1000) | verdict |
|---|---|---|---|---|---|
| called point, prevalence-matched | 0.4708 (n=24) | 0.0788 (n=92) | 0.3920 | 0.1762 | **SUPPORTED** |
| called point, unrestricted | 0.4383 (n=26) | 0.0788 | 0.3595 | 0.1627 | **SUPPORTED** |
| screened point, prevalence-matched | 0.4019 (n=369) | 0.0788 | 0.3230 | 0.1192 | **SUPPORTED** |
| screened point, unrestricted | 0.3719 (n=400) | 0.0788 | 0.2931 | 0.1296 | **SUPPORTED** |

Every gap exceeds the **maximum** of 1000 arm-label permutations, not merely the 95th percentile. The
headline is the **called** arm; the screened bands keep every position AMRFinder merely screened and are
reported so the pollution is visible, never as the headline.

## The correction that nearly did not happen

Raw coupling falls sharply with prevalence **in both arms**:

| prevalence band | POINT | ACQUIRED |
|---|---|---|
| rare (<0.10) | 0.5366 (n=20) | 0.0937 (n=66) |
| mid (0.10–0.33) | 0.1690 (n=3) | 0.0524 (n=18) |
| common (≥0.33) | 0.0521 (n=3) | 0.0161 (n=8) |

The reading I had written, and was one step from publishing: *the arm-level claim holds, but the
determinants that actually drive most resistance — `gyrA_S83L` (0.041), `parC_S80I` (0.056),
`gyrA_D87N` (0.061) — are individually near-chance, so the mechanism fails where it matters most.*

**That reading is an artifact of the metric, and it is wrong.** The null fixes the *expectation* at a
given carrier count; it does not fix the *maximum*. The largest lineage holds **53 of 621** genomes, so
a determinant with 400 carriers **cannot exceed** 53/400 = 0.133 largest-lineage fraction however
coupled it is — leaving roughly **0.048** of coupling available to it in total. Scored against what was
achievable:

| determinant | n | coupling | ceiling | **normalized** |
|---|---|---|---|---|
| `gyrA_S83L` | 400 | 0.043 | 0.133 | **0.895** |
| `parC_S80I` | 335 | 0.055 | 0.158 | **0.755** |
| `gyrA_D87N` | 307 | 0.060 | 0.173 | **0.698** |
| `parE_I529L` | 48 | 0.873 | 1.000 | 0.977 |
| `ptsI_V25I` | 55 | 0.858 | 0.964 | 1.000 |
| `blaTEM-1` | 282 | −0.014 | 0.188 | **−0.134** |
| `sul2` | 301 | −0.012 | 0.176 | **−0.132** |
| `tet(A)` | 271 | 0.029 | 0.196 | 0.261 |

The dominant chromosomal determinants are **near-saturated, not near-chance**, while acquired
determinants at comparable prevalence sit at or below zero. Arm means on the normalized metric:
**point 0.6415 vs acquired 0.1188, gap 0.5227 against a permutation max of 0.2184.**

**The arm gap is larger under normalization, not smaller.** A corroborating signature: the
*unrestricted* normalized mean (0.6555) is **higher** than the prevalence-matched one (0.6415) —
prevalence matching was only ever compensating for the ceiling, and once the ceiling is handled the
one-sided restriction stops being needed. That is what a correct fix should do.

The frozen bar's verdict rule is defined on the **raw** score and was applied to the raw score
unchanged. Normalization ships as a secondary band, never as a retro-fit of the rule.

## A second defect: the filter discarded real determinants

Reconciling the script's synonymous count (156,471) against the bar's (137,288) exposed an unrelated
bug. The two numbers differ **exactly**: 137,288 true wildtype rows + **19,183 unparseable symbols** =
156,471. The counter conflated "wildtype screen row" with "symbol my regex cannot parse".

The unparseable bucket is not all wildtype. The single-residue regex `^(.*?)_([A-Z])(\d+)([A-Z*])$`
cannot express indels (`ftsI_N337NYRIN`), frameshifts (`cirA_S90YfsTer15`), nonsense (`nfsA_Q67Ter`),
deletions (`23S_G543DEL`) or **promoter positions with negative coordinates** (`ampC_C-42T`). **12 of
the discarded symbols are ones AMRFinder actually CALLS**, six clearing the 5-carrier floor — about a
fifth of the chromosomal arm.

The bar specified this filter to remove *wildtype screen rows*; dropping indels was never its intent, so
the fix is a correction rather than a re-spec. It ships as a **sensitivity band** rather than a silent
reshaping of the frozen substrate
(`wiki/determinant_clade_coupling_filter_sensitivity_2026-09-10.json`):

| filter | point families | point mean (raw) | gap | verdict |
|---|---|---|---|---|
| frozen, single-residue | 24 | 0.4708 | 0.3920 | SUPPORTED |
| corrected, all variant forms | **30** | 0.4794 | 0.4006 | SUPPORTED |

Normalized: 0.6415 → 0.6425. **A real defect that is not load-bearing for the verdict** — which is
precisely what a sensitivity band exists to establish, in either direction.

## What this does and does not establish

**Does:** the mechanism's *premise* — that acquired and chromosomal determinants differ in clade
coupling, by a wide and controlled margin, measured directly rather than inferred from a proxy. The
kill-test assumed this; it is now measured.

**Does not:** that coupling **causes** the de-confounded learned-decoding failures. Population design
and inheritance-coupling remain perfectly correlated across all five failures and three successes, so
the attribution in `dna_decode/eval/regime.py` stays under-determined. Separating them needs the regime
comparison itself, not this measurement.

## Honest limits

- **One organism, one lineage definition.** E. coli, 7-locus MLST. ST is a **coarse** unit; a finer
  phylogeny could move every family in both arms.
- **AMRFinder's own `Subtype` is the mobility proxy.** `POINT` is a chromosomal substitution by
  definition, but this is the tool's classification, not an independent determination of mobility.
- **The normalized metric divides by a small denominator at high prevalence** (~0.048 at n=400), so a
  normalized score there is far noisier than one for a rare family. It corrects a real bias in the raw
  score; it does not make the two ends equally precise.
- **`ceiling_largest_fraction` is a theoretical bound.** It assumes carriers could be packed into the
  biggest lineage arbitrarily — right for a metric ceiling, but no real determinant is free to
  distribute itself.
- **Prevalence matching is one-sided** (removes 31 POINT families and 0 ACQUIRED), which is why the
  unrestricted band ships beside it.
- **Genomes are AMR-cohort leftovers**, enriched for resistance and not a natural population.

## Reusable

**A null that holds N fixed controls the expectation, not the ceiling.** Where the maximum attainable
value of a statistic varies with the same quantity the null is matching on, a deviation-from-chance
score is still confounded — and the confound runs in the direction that makes the highest-prevalence,
most-important cases look inert. Report the achievable range alongside the score, or the metric will
manufacture a finding about exactly the cases you care most about.

**Reconcile a count against its own pre-registration before trusting either.** The ceiling error was
caught by reading a stratification; the filter error was caught only because the bar recorded a number
the script later disagreed with. A frozen bar is a checksum on the substrate, not just on the verdict.
