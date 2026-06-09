# Wider-AMR transferability synthesis — 2026-06-08

> Where does the deterministic AMR decoder's rule transfer, and where does it break? Synthesis across all
> organism×drug validations run via `scripts/organism_drug_validate.py` (deployed DRUG_RULE applied
> UNCHANGED; only AMRFinder `-O` varies). Each: NCBI Pathogen-Detection cohort, ≥15R/15S where available,
> NCBI labels. North-star framing: this map IS the decoder's honest scope-of-trust statement.

## The map

| organism (phylum) | drug | acc | sens | spec | verdict | boundary |
|---|---|---:|---:|---:|---|---|
| E. coli (Pseudomonadota) | cipro/cef/tet/gent | 0.83–0.93 | — | — | **transfers** | — |
| Klebsiella pneumoniae | ciprofloxacin | **1.000** | 1.000 | 1.000 | **VALIDATED** | — |
| Klebsiella pneumoniae | ceftriaxone / tet / gent / meropenem | 0.93–1.0 | — | — | **transfers** | — |
| Acinetobacter baumannii | meropenem | 0.500 | 1.000 | 0.000 | FAILS | **CONTENT** |
| Pseudomonas aeruginosa | meropenem | 0.500 | 1.000 | 0.000 | FAILS | **CONTENT (+EXPRESSION)** |
| Campylobacter (Campylobacterota) | ciprofloxacin | 0.500 | 0.000 | 1.000 | FAILS | **TUNING** |
| Enterobacter cloacae | ceftriaxone | 0.455 | 0.375 | 0.667 | FAILS* | **EXPRESSION** |
| Salmonella | ciprofloxacin | 0.567 | 0.133 | 1.000 | FAILS | **TUNING + CONTENT** |

\* Enterobacter N=11 (8R/3S), label-limited — directional only.

## Three boundary flavors (the load-bearing taxonomy)

The deterministic rule is `R iff (#curated drug-relevant determinants) ≥ threshold`, with per-drug
subclass/gene refinements. It fails in exactly three distinguishable ways:

1. **CONTENT** — counts the *wrong* genes. The organism carries determinants that match the drug's class
   but don't confer the phenotype. Acinetobacter (intrinsic blaOXA-51 *gene presence*) and Pseudomonas
   (intrinsic nalC/oprD *point-mutation variants*) both → spec 0 on meropenem. **Fix = organism-specific
   content curation** (strength tiers / exclude intrinsics). Family-level auto-intrinsic-flag SURVIVES as
   an automated fix (see `self_calibration_falsifier_2026-06-08`).

2. **TUNING** — right genes, wrong integer. The mechanism caller transfers perfectly but the count
   *threshold* is organism-specific. Campylobacter cipro: single gyrA T86I confers clinical R, but the
   E. coli-tuned threshold=2 misses it → sens 0. **Fix = per-organism threshold.** Auto-threshold LOO
   selection SURVIVES as an automated fix (same falsifier doc).

3. **EXPRESSION** — right genes, can't see their *regulation*. Resistance is via over-expression /
   derepression of an intrinsic gene, invisible to gene-presence. Enterobacter cef (derepressed AmpC) and
   the residual Acinetobacter/Pseudomonas FN (ISAba1-driven OXA-51, efflux up-regulation). **No
   presence-based fix exists** — needs promoter/regulatory-region inference from the assembly (the open
   `/hypothesise` IS-element-upstream strand). This is the hard floor.

## What this tells the tool

- **Trust zone is rule×organism-specific, NOT just "Enterobacterales".** cipro/cef/tet/gent transfer
  across E. coli + Klebsiella. But Salmonella cipro FAILS under the deployed rule (qnr + single-gyrA, needs
  a broad counter @ threshold 1) and Enterobacter cef FAILS (derepressed AmpC) — both *are* Enterobacterales.
  Lesson: even within a family, the right counter+threshold is organism-specific. The Klebsiella-vs-Salmonella
  cipro contrast is the proof — the SAME QRDR-point-only design choice makes Klebsiella perfect and Salmonella
  fail. There is no family-wide rule; calibrate per organism.
- **Carbapenems on non-fermenters (Acinetobacter, Pseudomonas) are out of trust zone** — intrinsic content
  + expression-driven R defeat presence-based calling. Honest output here = abstain / flag, not predict.
- **CONTENT and TUNING are cheaply auto-fixable** from a ≥15R/15S cohort (both falsifiers survived). The
  `calibrate_organism(cohort)` build (tracked hypothesis H3) would convert each new organism from a
  hand-curation task into an automated calibration — modulo the EXPRESSION floor, which it cannot cross.
- **EXPRESSION is the frontier.** Every residual FN across organisms traces to expression-level mechanisms.
  The single highest-value capability extension is reading regulatory context (IS-elements, promoter/
  derepression mutations) from the assembly the decoder already has.

## Method honesty (applies to every row)

1 organism × 1 drug × N≈30 (some smaller), NCBI labels (different source/curation, NOT a different-lab
study). Verdicts are directional cross-organism signals, not graded clinical-accuracy claims. No refinement
has been wired into the deployed rule — all fixes documented as candidates pending independent cohorts.
