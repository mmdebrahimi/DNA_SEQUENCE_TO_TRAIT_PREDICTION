# AMR-Portal Tier 3/4 curation — big-picture plan (2026-07-01)

How to tackle the ~44 unscored AMR-Portal cells (Tier 3 + Tier 4 of `wiki/amr_portal_unscored_triage_2026-06-28.md`)
without falling into the intrinsic-gene over-call trap. This is a **methodology + phased backlog**, NOT a sweep —
each cell is a deliberate curated-rule effort, and an honest ABSTAIN is a valid, common outcome.

## The load-bearing idea: mechanism-DECODABILITY, not organism, decides feasibility

A determinant caller can only decode resistance that is **encoded in gene/allele presence or a target-site
point mutation**. The 44 cells sort into 4 mechanism classes — this is the real triage axis:

| Class | Decodable? | Rule shape (existing) | Examples in the 44 |
|---|---|---|---|
| **A. Target-site point mutation** in a known gene | YES (clean) | `qrdr_point` (gyrA/parC), `rpoB` point | Neisseria gonorrhoeae cipro (gyrA), Staph cipro (gyrA/grlA) + rif (rpoB) |
| **B. Acquired gene** (horizontally-transferred) | YES *if* intrinsic look-alikes are excluded | `subclass_any` acquired / `gene_prefixes` (+ acquired-only curation) | Enterobacter/Proteus/Serratia cipro-qnr/gent-aac/tet-tetA/TMP-SMX; Enterococcus van/aac; Acinetobacter acquired-only |
| **C. Intrinsic / chromosomal-derepression** | NO → ABSTAIN | `EXPRESSION_FLOOR` / ABSTAIN | intrinsic AmpC cef (Enterobacter/Serratia/Proteus), OXA-51 carbapenem (Acinetobacter), intrinsic-gent (Enterococcus low-level) |
| **D. PBP / porin / efflux-expression** | NO (not gene-presence) → ABSTAIN | ABSTAIN | Strep pneumoniae β-lactams (PBP mutation), Pseudomonas porin/efflux, Neisseria penA-mosaic cef (borderline) |

**The trap (why curated, not swept):** applying the generic E. coli `DRUG_RULE` to a class-C/D organism×drug
over-calls (spec → 0) because the organism carries the class gene intrinsically. Every cell must pass an
intrinsic-exclusion review + a specificity guard before it can be endorsed.

## The reusable per-cell curation RECIPE (7 steps — the same for every cell)

1. **Determinant sourcing** — the organism×drug mechanism(s) from CARD + AMRFinder (DB already installed) +
   literature (`/research`; WHO/CLSI where a catalogue exists). Output: the acquired genes / point mutations
   that CAUSE resistance in *that* organism.
2. **Mechanism-class triage** — tag each mechanism A/B/C/D (table above). If the dominant mechanism is C/D →
   **stop, ABSTAIN** (document as `ABSTAINS_BY_DESIGN`, no score).
3. **Intrinsic-exclusion** — list the organism's intrinsic/look-alike genes that would over-call (AmpC,
   OXA-51, OqxAB, intrinsic tet/aac) → **exclude** them; rule counts acquired/strength-specific determinants
   only.
4. **Rule draft** — reuse an existing shape with an organism-specific determinant list: `qrdr_point` (class A),
   acquired `subclass_any`/`gene_prefixes` (class B). New organisms get an `organism_rules/<org>_amr.py`
   module (the `pneumo_amr.py` pattern already exists) or a `calibrated_amr_rules.json` registry entry
   (OPT-IN, not frozen-default).
5. **Independent powering check** — `amr_portal_score_independent.py` on the provenance-disjoint AMR-Portal
   isolates (≥10 R AND ≥10 S). Under-powered → report, don't endorse.
6. **Validate + specificity guard** — score sens/spec/acc; run the `organism_drug_validate.py` pattern that
   catches the spec-collapse over-call. **Endorse only if spec clears a pre-set bar** (draft: spec ≥ 0.85 on
   the powered independent set) AND the mechanism strata reproduce.
7. **Gate + ship** — namespace-separate artifact; add to the AMR-Portal report card as `SCORED` (endorsed) or
   `ABSTAINS_BY_DESIGN` / `INDETERMINATE`. Frozen surface stays byte-unchanged (registry entries are OPT-IN).

## Phased backlog (decompose — ranked VOI ÷ effort; ~44 cells)

- **Phase A — point-mutation quick wins (class A; highest VOI, lowest effort, reuses `qrdr_point`).**
  Neisseria gonorrhoeae **cipro** (gyrA S91/D95) + tet; **Staphylococcus cipro** (gyrA/grlA) + **rif** (rpoB).
  ~4-5 cells. These are the closest analogue to the already-validated Campylobacter/Salmonella cipro cells.
- **Phase B — acquired-gene Enterobacterales (class B; broadest cell count).** Enterobacter (10), Proteus (3),
  Serratia (3) for **cipro / gent / tet / TMP-SMX** — reuse the E. coli acquired determinant lists + the
  intrinsic-**AmpC exclusion for cef** (cef → likely ABSTAIN per class C). ~14-16 candidate cells; expect a
  mix of SCORED (cipro/gent/tet/TMP-SMX) + ABSTAIN (cef).
- **Phase C — Gram-positive determinant (class A/B; organism-specific).** Staph gent/tet/TMP-SMX; Enterococcus
  cipro/gent(high-level markers)/tet. ~6-8 cells; the oxacillin/mecA + cefoxitin-surrogate lesson governs
  Staph β-lactams.
- **Phase D — ABSTAIN documentation (class C/D; no scoring, just honest non-cells).** Strep pneumoniae
  β-lactams (cef/mero — PBP), intrinsic AmpC cef across Enterobacterales, Pseudomonas non-mero efflux,
  C. difficile / anaerobe. ~10-12 cells → documented `ABSTAINS_BY_DESIGN` on the card.

**Tier 3 (Campylobacter tet/gent):** tet = `tetO`/`tetM` ribosomal-protection (class B, curatable, different
from E. coli efflux `tetA`); gent = aac acquired (class B). Both are Phase-B-shaped — a curated Campylobacter
rule, then powering + validate.

## Tooling to reuse (no new infra needed)
- `scripts/amr_portal_score_independent.py` (+ `RULE_ORGANISM` map) — the powered provenance-disjoint scorer.
- `scripts/organism_drug_validate.py` — the specificity-collapse guard (the intrinsic over-call detector).
- `dna_decode/organism_rules/` (pneumo_amr already started) + `calibrated_amr_rules.json` (OPT-IN registry).
- `dna_decode/data/mic_tiers.py` loci catalogs + `scripts/drug_mechanism_audit.py` (per-drug mechanism audit).
- `dna_decode/data/experimental_drug_rules.py` — the pattern for a non-frozen scorer-local overlay.
- Data: AMR-Portal parquet (on D:, present) + AMRFinder DB (installed) — determinant sourcing + labels in hand.

## Guardrails + falsifiers (per cell, pre-committed)
- **Intrinsic-exclusion review is mandatory** before any endorsement (the trap).
- **Specificity falsifier:** spec < 0.85 on the powered independent set ⇒ the rule over-calls ⇒ do NOT endorse
  (ABSTAIN or re-curate). Sens-only success is not enough (the high-sens/low-spec → suspect-the-label lesson).
- **Provenance-disjoint + powering** (≥10R/≥10S) enforced by the scorer; under-powered ⇒ report-not-endorse.
- **Namespace separation** — new cells go to the AMR-Portal card (opt-in registry), NEVER the frozen
  deployed surface (`amr_rules.py`/`calibrated_amr_rules.json` byte-unchanged unless a deliberate ratify-first
  freeze amendment).
- **Honest ABSTAIN is a WIN** — a documented `ABSTAINS_BY_DESIGN` (class C/D) is a correct outcome, not a gap.

## Honest expectation (set before starting)
Of ~44 cells, a realistic split: **~15-20 SCORED** (Phase A + B + C decodable), **~20-25 ABSTAIN** (Phase D +
intrinsic cef). The value is a broader *honestly-tiered* AMR-Portal card, not "44 scored." Each cell is ~a
half-to-one `--until-mvp` session; the whole backlog is many sessions — pursue Phase A first, one cell at a
time, and re-rank after each (the Campylobacter/Neisseria cipro cells will calibrate the effort estimate).

> **DONE 2026-07-01 (Phase A, cell 1):** Neisseria gonorrhoeae ciprofloxacin SCORED -- gyrA-QRDR (Ser91/Asp95) rule, acc 0.968 / sens 0.943 / spec 0.99 on 5618R/6406S provenance-disjoint AMR-Portal measured-AST (spec>>0.85 floor; strata gyrA-present 0.988 vs absent 0.048). organism_rules/neisseria_amr.py + scripts/neisseria_cipro_amr_portal_validate.py + 6 tests; wired into the AMR-Portal card overlay section. Frozen surface byte-unchanged. The 7-step recipe transfers cleanly; per-cell cost ~1 attended run. Next: Staphylococcus cipro (gyrA/grlA) or the acquired-gene Enterobacterales (Phase B).

## Recommended first executable move (next session, not now)
**Phase A, Neisseria gonorrhoeae ciprofloxacin** — gyrA point mutation (S91F/D95), the cleanest class-A cell,
directly reuses the validated `qrdr_point` machinery. One curated rule + powering + validate → the first new
Tier-3/4 SCORED cell, and it calibrates the per-cell recipe cost for everything after.
