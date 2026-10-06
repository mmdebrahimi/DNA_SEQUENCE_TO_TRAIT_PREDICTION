# Open threads — what a fresh session should pick up

Short by design. **Transient** state only: what is in flight, what is waiting on the user, what was just
learned that isn't durable yet. Durable findings belong in `CLAUDE.md` / `wiki/`; scope facts belong in
`scripts/project_status.py` (derived, never written).

Prune aggressively. A stale entry here is worse than an empty file.

_Last updated: 2026-10-05._

> **Pruned 2026-10-05.** This file had sat at 2026-09-28 / 428 lines, and its own "known-stale" section
> had itself gone stale in the opposite direction (it warned that 46 traits was really 44; the live count
> is **46** again after `dna-identify` and `dna-tb` shipped). Everything removed was already durable in
> `CLAUDE.md` (gentamicin v2 lock, HCMV contracting, the declined NNRTI curation, PEAR's reclassification,
> the doubt-layer firing rate) or in `project_state/` (the five project families, F-C/F-D). Checked before
> deleting, not assumed.

---

## IN FLIGHT — deep-research campaign F1 / F3 / F4

Plan `plans/Deep_Research_Missing_Part_2026-09-28_Plan.md`. User ratified full breadth (F1–F4).
Screening spine is built and green (`scripts/deep_research_screen.py`, 13 tests; runs both
`rejection_gates` G1–G10 and `regime.screen_proposal`, and refuses to call anything promising until they
have run).

| family | state |
|---|---|
| **F2** condition coverage | **ANSWERED NEGATIVELY.** "Public K-12" adds samples, not condition breadth — glucose share 60% → **82.7%** of annotated rows, **42%** carry no carbon annotation, only **3** carbon sources reach 10 samples. Computed from the repo's own metadata; the "contact us" wall covers only `log_tpm`. |
| **F1** label substrates | one shallow search. DMS axis is a clean negative. **Not searched:** measured-AST collections, AllTheBacteria, EBI AMR Portal, Zenodo/Figshare thesis deposits. |
| **F3** de-confounded population designs | **NEVER SEARCHED. On the critical path** — population DESIGN, not organism complexity, is the measured barrier to the north star's step 4. |
| **F4** cross-assay dosage transfer | **NEVER SEARCHED.** |

**Relaunch SEQUENTIALLY (or 2 at a time) and forbid sub-agent spawning.** The first attempt spawned 12
concurrent agents; all died on a session rate limit with **zero** final reports. First-party search worked
fine immediately after — the fan-out was the constraint, not the session.

## IN FLIGHT — the cross-organism transfer ladder (drafted 2026-10-05, NOT saved)

A 9-step plan to turn the project's **single** step-4 data point (E. coli-tuned essentiality decoder applied
UNCHANGED to human, AUROC 0.5805) into a 5-rung ladder, to separate **(A)** decay-with-distance from
**(B)** a eukaryote-boundary cliff. Drafted + adversarially reviewed in conversation; **never written to
`plans/`** — deliberately, because the review mutated it (the save-once rule).

**Three findings that change it, all measured, none in the plan as drafted:**

1. **AUROC is the WRONG primary metric — it is tie-dominated.** 83.1% of human essentials score EXACTLY 0
   (E. coli 62.1%), and 11 of 13 patterns have P(hit|non-essential)=0.0000 — the decoder is **silent**, not
   wrong. Coverage and precision-where-it-fires move in OPPOSITE directions across the two rungs
   (coverage 0.3789→0.1689; precision 0.3684→**0.9583**) and AUROC blends them. **Use coverage.**
2. **The miss is ~2.5× host-specific vs phrasing, and the ladder survives** —
   `wiki/essentiality_missed_vocabulary_2026-10-05.md`. Mechanism (iii) floor 141/566 (24.9%) vs
   mechanism (ii) floor 57/566 (10.1%). So a cheap pre-test does **not** pre-empt the ladder, but ~10% is a
   catalogue PHRASING gap unrelated to phylogeny and must be separated out before any rung is scored.
3. **The FBA walls do not apply.** `fba/essentiality_labels.py` marks S. aureus `LABEL_WALLED` and
   P. aeruginosa `MODEL_WALLED`, but both are **model-gene-join** walls; this decoder joins to a genome
   annotation, so two "blocked" rungs are reachable. Reachability is still a HYPOTHESIS until fetched.

**Open AUTHORITY forks:** re-run `/technical-plan` around coverage, or `/save-plan` as-is? (The review also
wants a W0 schema probe before any parser, and the functional-class step demoted to a provenance-gated
audit — it is otherwise a second adjudication system.)

## The system design — drafted, awaiting ratification

`plans/Hybrid_Decoder_Architecture_Plan.md`. The hybrid is **not** "catalog + ML predictor" (that framing
scored 0 survivors). Measured shape: **CALL / DOUBT / EVIDENCE**.

| layer | status |
|---|---|
| **L1 CALL** — deterministic curated rules | shipped, 128 cells |
| **L2 DOUBT** — "this call may be incomplete, and why"; never a competing call | **shipped and firing** (`dna_decode/eval/doubt.py`); nobody in the field ships this |
| **L3 EVIDENCE** — de-confounding, nulls, denominators, leakage, provenance | built; **four** cohort-evidence layers now render (lineage, source-concentration, prospective, species-composition) |
| **L4 LEARNED** — forward/inverse, orthogonal modalities; molecular + constructed ONLY | shipped, bounded |

## Waiting on the user — authority calls, not executor tasks

1. **Re-score `Klebsiella × meropenem` on *K. pneumoniae* only?** (new 2026-10-05)
   Measured: the cohort is **38% *K. aerogenes*** and **all 16 false negatives are those off-species
   isolates** — on *K. pneumoniae* there are zero misses (sens 12/12). The pneumoniae-only matrix is
   **sens 1.000 / spec 0.880 / acc 0.919, N=37** against the published 0.467 / 0.900 / 0.683.
   Publishing it edits a **frozen unit** of the reproducibility freeze and changes a SCORED number on the
   trust surface. Current answer is **disclose** (`species_composition` layer, augment-only).
   `wiki/meropenem_fn16_diagnosis_2026-10-05.md`.

2. **Drop lone-porin counting from the meropenem rule?** *(evidence added 2026-10-05, pointing the
   other way — see the closed lever below)*
   All 3 false positives on the AR Bank cohort are porin-only, and the rule's `threshold=1` on a lone
   `ompK35`/`ompK36` truncation over-calls. But it edits sha256-pinned `amr_rules.py` → retires the
   active 2026-08-31 v2 lock and restarts the prospective clock for every cell, motivated by **3
   isolates** on a cohort whose own powering gate HARD-FAILED (n=5 susceptibles).

3. **Whether a single-source cell warrants more than disclosure.** 3 of 10 SCORED AMR cells rest on one
   BioProject. Current answer is *disclose*, namespace-separate. Demoting them is a scope decision.

4. **Whether to compress the other long CLAUDE.md bullets.** The file loads every session; several long
   bullets cite a resolvable memo, so their derivations could become pointers. Measure first:
   `uv run python scripts/claude_md_weight.py`. Two long bullets have no external store and must stay
   whole — the tool already protects them.

## Cheap untried levers (executor work, no authority needed)

- ~~**`dna-identify` cross-check on the 16 meropenem FN.**~~ **DONE 2026-10-05** — ran it, and it
  CORROBORATES: 16/16 FN abstain (`above_max_distance`, 0.1129–0.1339, inside the 0.10849–0.13386 band
  `thresholds.py` records for K. aerogenes), 0 called `klebsiella_pneumoniae`, controls 12/14 called at
  0.0011–0.0090. **Partition agreement with the GenBank labels: 30/30, zero disagreements.** Scope: a
  CLOSED-SET router establishes NOT-pneumoniae, not *aerogenes* specifically.
  `wiki/meropenem_fn_identify_crosscheck_2026-10-05.json`.
- ~~**Re-run AMRFinder under the correct `-O` for the 23 off-species isolates.**~~ **ILL-POSED — CLOSED
  2026-10-05 without running it.** Measured from the pinned image (`amrfinder -l`): the organism list
  offers `Enterobacter_asburiae` / `Enterobacter_cloacae` / `Klebsiella_oxytoca` /
  `Klebsiella_pneumoniae` and **no `Klebsiella_aerogenes`**. So there is no "correct flag" to re-run
  under, and the species contamination in that cohort is **not remediable by re-running** — a structural
  limit of the tool's supported set, not a to-do. (The nearest available options sit under *Enterobacter*,
  which is where this organism's former name sits; choosing one would be a judgement about an unsupported
  species, not a correction.) This is also why the 23 isolates cannot simply be re-scored "properly".
- ~~**The 2 *K. aerogenes* true positives are unexplained.**~~ **DONE 2026-10-05, and they are NOT one
  class.** `GCA_003951185.1` is called R by `ompK35_G41TfsTer32` — a **lone porin truncation**, its only
  counted carbapenem determinant. `GCA_003951605.1` is called R by `blaNMC-A`, a real class-A
  carbapenemase. **This bears on authority fork 2 and points the other way:** a lone-porin call that is
  RIGHT, so dropping lone-porin counting converts a true positive into a false negative. n=1 against the
  3 AR Bank FPs — not settling, but the fix now has a measured cost on the same arm.
- **The 3 Klebsiella cohorts with incomplete cached AMRFinder runs** (gentamicin 3/60, tetracycline 33/60,
  ceftriaxone 54/60) have their species cross-tab withheld for that reason. Completing them needs Docker.
- **The two *K. aerogenes* true positives are unexplained** — they were called R, so they carry some
  CARBAPENEM-subclass determinant. Offline, cheap, not chased.

## The FBA switch cell, as it now stands

Direction is fine, **silence** is the problem: 61–76% of genes emit one identical ratio for every
condition. Do **not** build the relative ranking rule (oracle ceiling +4 genes, and a deployable version
must infer k = the original problem restated). **Do not quote a raw exact-set rate** — a conditionally
essential gene is two-sided by construction, so a constant call can never match one. The measured
bottleneck is fixed in both datasets: PRECISE-1K ∩ Keio carbon = **11 of 28**, 621 of 1,035 samples
glucose.

## Known-stale / do not trust without re-deriving

- Any **cell count or trait count written in prose**, including in this file.
  `uv run python scripts/project_status.py` is the authority — it is derived, never written.
- `wiki/project_distillation_2026-08-29.md` says 46 CLI traits. That was wrong when written (it was 44)
  and is now right again by coincidence, which is the cleanest possible illustration of why a written
  count is not evidence.
