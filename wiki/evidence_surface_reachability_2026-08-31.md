# The evaluation machinery was card-only — reachability audit (F-C, 2026-08-31)

**One sentence:** four per-cell disclosure layers rendered on the standing report card and only **two**
reached a decoder call; the missing pair included the one caveat that explains this project's own
worst-known metric gap.

F-C actions 1–3 of `plans/Hybrid_Decoder_Architecture_Plan.md`. Frozen AMR surface byte-unchanged.

---

## The audit

Reachability fails independently at three levels, so all three are checked:

| level | question |
|---|---|
| **card** | the layer renders on `wiki/decoder_validation_report_card.json` |
| **record** | it reaches `trust_block` → the JSON a decoder call emits |
| **human** | it reaches the *printed* CLI output |

**Before:**

| layer | on card | in record | human |
|---|---:|---:|---|
| `doubt_layer` | 16 | 16 | ❌ *(JSON only — shipped that morning)* |
| `prospective` | 2 | 2 | only when it **contradicted** |
| `lineage` | 10 | **0** | ❌ |
| `source_concentration` | 10 | **0** | ❌ |

**After:** all four at all three levels.

## Why it mattered, concretely

`escherichia_coli_shigella × gentamicin` reports **sens 0.893** from a cohort that is **95% one
BioProject** and contains **zero `rmt` carriers**; source-diverse measurements of the same cell with
the same frozen rule report **0.523**. A caller deciding on 0.893 could not see the caveat that
explains it. That call now prints three honest layers together:

```
validation: INDEPENDENT_MEASURED -- acc 0.987 (N=15697)  (... || PROSPECTIVE REGRESSION:
            post-lock sens 0.429 on N=62 isolates public after 2026-06-13 ...)
source concentration: SINGLE-SOURCE -- 4 BioProject(s), one dominant (95% of the cohort).
            The metric above describes that source's isolates; it is a narrow estimate,
            in either direction, not necessarily an inflated one
```

Klebsiella × ceftriaxone (36 BioProjects) correctly prints no concentration line — the caveat is
decision-relevant, not boilerplate.

## Four judgements worth carrying

**A JSON-only disclosure is not a disclosure.** The `doubt` block shipped in the record in the morning
and printed nowhere; the human output showed a *static* list saying an S call "can't rule out an
uncatalogued substitution" while never saying that *this* genotype has one. Found by auditing my own
work hours after shipping it.

**Silence is honest in exactly one case.** A doubt line prints unless the flag was assessed and found
nothing. `not-applicable` (position-based catalog — the flag could never fire) and `not-assessable`
(the input path never surfaced the substitutions) both print, because reporting either as silence is a
false clean bill.

**The CI is the result, so a compact renderer must not drop it.** `lineage_one_line` prints the
effective-N disclosure always but withholds the weighted point estimate whenever its interval is
missing — mirroring `_assert_weighted_renderable`, which refuses to render such a point at all.
Effective lineage N is tiny; a bare 0.889 would mislead.

**The error is not directional.** A 97%-single-source ciprofloxacin cell reads *pessimistic*
(spec 0.700 vs 0.988 on 8 BioProjects). The caveat says the estimate is **narrow**, never inflated.

## Honest limits

Reachability, not correctness — this says a caller *can see* a caveat, not that they will act on it.
`prospective` is on 2 cells only. The three non-report-card trust surfaces (HIV, TB, pgx cards) were
**not** audited and are the named follow-on. `lineage` is rendered at the finest available Mash
threshold; the coarser rung collapses harder and is in the record but not the printed line.

## Artifacts

`scripts/evidence_surface_inventory.py` → `wiki/evidence_surface_layer_inventory.json` (derived, never
hand-listed) · `trust_surface.{DISCLOSURE_LAYERS, _cell_layer_for, lineage_one_line,
concentration_one_line}` · 12 tests in `tests/test_evidence_surface_reachable.py`.
