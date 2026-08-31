# The blind-spot framing sweep killed my own recommendation, and the winning framing is curation

`/innovate` run, 2026-08-31. Self-invoke trigger **(b)** — ≥2 approaches failed the same way — so the
**framing sweep** ran rather than a single-framing pass. Four framings, nine candidates, every verdict
produced by an **executed** kill-test through the gated runner.

## Result

| framing | survived | killed |
|---|---:|---:|
| **F4 — curation, not computation** ← **WINNER** | **3** | 0 |
| F3 — ship self-awareness, not prediction | 2 | 0 |
| F2 — the comparator is the free flag, not the catalog | 1 | 1 |
| **F1 — INCUMBENT: find a better scorer** | **0** | **2** |

**The incumbent framing produced zero survivors.** Both of its candidates died to executed tests.

## The two kills, and the first one is mine

**`F1-ddg-cheapest` — KILLED.** Two turns ago I recommended a ΔΔG_bind pilot as *"the cheapest decisive
experiment."* The disproof: `wiki/hiv_blindspot_position_novelty_2026-07-11.json` records a **zero-tool
deterministic position-novelty flag already recovering 60.4% of the EFV blind spot** (lift 4.69, no model,
no structures, no GPU). Something cheaper already exists **and already works**. The recommendation was
made without reading the artifact that answers it.

**`F1-not-pocket` — KILLED.** The plausible theory that ESM2 fails because the blind spot is *not*
pocket-mediated (and so structurally invisible to any binding scorer) is disproved by
`wiki/hiv_blindspot_pocket_localization_2026-07-09.md`, which returns **VERDICT: GO** — blind-spot
resistant isolates are **3.05× burden-adjusted enriched** for functional NNRTI-pocket mutations.

Note what those two together mean: the ΔΔG premise is **sound** (it is pocket-mediated), but ΔΔG is **not
the cheap next move** (a free flag already covers most of it). Those are different questions and I had
merged them.

## The winning framing — F4: the blind spot is a CURATION gap

Three survivors, all executed-verified:

**`F4-drivers-absent`** — the named blind-spot drivers are **absent from the deployed catalog**. The
pocket artifact lists them with counts: **V179D ×12, A98G ×10, H221Y ×7, F227C ×5, V108I ×4, V179E ×3** —
all described there as known functional NNRTI mutations. The shipped `NNRTI_RT_MAJOR_DRMS` covers only the
eight Stanford major positions {100, 101, 103, 106, 181, 188, 190, 230}. Positions 179, 98, 221, 227 and
108 are **not in it**. The blind spot is, in substantial part, a catalog that was scoped to majors
by design.

**`F4-hiv-not-frozen`** — `hiv_amr.py` is **not** pinned by `prospective_lock_manifest_2026-06-22.json`
(which pins `amr_rules.py`, `calibrated_amr_rules.json`, `mic_tiers.py`, `shipped_decoder_surface.py`,
`cohort_manifest.py`). So **curating the HIV catalog does not invalidate the prospective lock or the
reproducibility freeze** — unlike the gentamicin `rmt` fix, which does. This is the structural reason the
HIV gap is cheap to close and the bacterial one is not.

**`F4-same-class-as-rmt`** — no cross-cell catalog-completeness screen exists. The HIV blind spot and the
gentamicin `rmt` gap are the **same failure class**: a curated catalog missing a known determinant family,
found only when an independent label set exposed it. One screen would serve both.

## F3's two survivors are a second, independent product move

**`F3-flag-not-shipped`** — `position_novel` appears in `dna_decode/eval/position_novelty.py` but **not**
in `dna_decode/data/hiv_amr.py` or `dna_decode/cli.py`. A working self-awareness signal sits in the
evaluation layer and never reaches a user.

**`F3-abstention-tier`** — `AbstentionVocab` already exists in `cell_registry.py`, so surfacing the flag
needs **no new evidence-tier concept**, only wiring.

The artifact itself frames the flag correctly: *"a 'catalog call may be incomplete' self-awareness flag,
NOT a resistance predictor."* That is exactly the shape this project's honesty rails want.

## F2's survivor sets the bar for anything future

**`F2-nothing-beats-flag`** — no method recorded in this repo has been shown to beat the flag's 0.604 on
the blind spot. **So the comparator for any new method is 0.604 at zero cost, not the catalog's 0.962.**
(`F2-ddg-marginal` was killed: an artifact *does* already frame ΔΔG against the flag.)

## Method note — two FALSE SURVIVORS were caught by verifying, not by trusting

`F1-ddg-cheapest` "survived" **twice** before it was killed. Both times its kill-test failed for the wrong
reason — first an `AttributeError` from over-clever JSON parsing, then a key-name guess that matched
nothing (the real key is `median_flag_sens_on_blindspot`). The engine correctly read each non-zero exit as
*"disproof not found"*; the claim was simply never tested.

**A survivor whose test failed for the wrong reason is indistinguishable from a real survivor in the
output.** What caught it was checking a verdict that looked wrong against a fact I already knew. The
third version carries a guard assert that fails loudly when nothing parses, so the test can no longer
return a silent false survivor.

This is the same defect class the H1 discrimination controls exist for — and it arrived through the one
predicate kind (`test-exit-0`) where controls are deferred in v0.

## Handoff — not executed

The top survivor is **F4**. It is a **curation** move on an unfrozen catalog, not a build:

1. verify each of the six drivers against a citable source (Stanford HIVDB / literature) — **no fabrication**;
2. decide whether they enter as majors or as a separate accessory tier with its own evidence label;
3. re-score the blind spot to measure what curation alone recovers;
4. only then ask what ΔΔG adds **over** curation + the flag.

Whether to edit a shipped catalog is a **scope decision**, so it stops here.

Reproduce: kill-tests and the framings ledger are in this session's scratchpad; the sweep is
`python ~/.claude/skills/soraya/scripts/framing_sweep.py --sweep <ledger.json> --cwd .`
