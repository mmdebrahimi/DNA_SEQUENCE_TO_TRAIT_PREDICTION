# dna-amr — Multi-Organism AMR Decoder — Capstone (v0.4.0) — 2026-06-07

> Consolidates the AMR arc into one deliverable. The deterministic mechanism-feature AMR decoder is
> validated across **6 drugs × 4 organisms × 4 mechanism classes, spanning the gram divide**, deployed as
> `dna-amr` / `dna-decode`. This is the milestone; further breadth is diminishing-returns.

## What the tool is

`dna-amr` (in the unified `dna-decode` CLI): a deterministic, interpretable antibiotic-resistance R/S
caller. It reads AMRFinderPlus's CURATED determinants from a genome and applies a per-drug rule
(`amr_rules.py::DRUG_RULE`). NOT embeddings, NOT a black box, NOT a clinical tool — a transparent,
auditable decoder that names the determinants driving every call + its own blind spots.

## The validated matrix

| Drug | Mechanism class | Rule | Organisms validated |
|---|---|---|---|
| ciprofloxacin | QRDR target point-mutations | `qrdr_point` ≥2 (gyrA/parC/parE) | E. coli 0.925 · Klebsiella 1.0 · Pseudomonas 0.867 |
| ceftriaxone | acquired ESBL/AmpC | ≥1, CEPHALOSPORIN/CARBAPENEM subclass | E. coli 0.933 · Klebsiella 0.80 |
| gentamicin | acquired AME | ≥1, GENTAMICIN subclass | E. coli 0.945 · Klebsiella 0.867 |
| meropenem | acquired carbapenemase | ≥1, CARBAPENEM subclass | Klebsiella 0.867 |
| tetracycline | acquired tet efflux/RPP | ≥1, `tet*` gene-prefix | E. coli 0.917 · Klebsiella 0.80 (efflux-limited) |
| oxacillin | acquired mecA (MRSA) | ≥1, METHICILLIN subclass | S. aureus: genotype sens 1.0 (label-limited) |

Organisms: **E. coli, Klebsiella pneumoniae, Pseudomonas aeruginosa, Staphylococcus aureus** (gram-neg ×3 +
gram-pos ×1). Validation tiers: in-cohort + held-out (cipro 0.862) + cross-source (NCBI, independent) +
cross-organism. Every per-drug rule beats naive AMRFinder ("any drug-class determinant → R") on independent
data.

## The one engineering principle (held across every case)

**Count the drug's SPECIFIC resistance determinants, not the broad drug-class bag.** The broad count
over-calls because organisms carry intrinsic/co-class genes that don't confer the specific resistance:
- cipro/tet: K. pneumoniae intrinsic OqxAB efflux (tagged QUINOLONE + TETRACYCLINE) → use QRDR point-muts /
  acquired `tet*` only.
- cef: intrinsic narrow β-lactamases (blaTEM-1) → use extended-spectrum subclass only.
- gent: streptomycin/kanamycin AMEs (aph/aadA) → use GENTAMICIN subclass only.
The refinement is mechanistically canonical AND cross-organism-robust (it also improved E. coli: cipro
cross-source 0.955→1.0, tet 0.833→0.917).

## The honest limits (the recurring binding constraint)

Every frontier confirmed the same thing: **the method generalizes; the de-confounded, reliably-labeled
substrate is the binding constraint.**
- **Blind spots (named in every S call's `undetectable_mechanisms`):** efflux-overexpression, porin-loss,
  regulatory — expression phenotypes invisible to curated determinants. Cause the FN at tet/Klebsiella
  (efflux) and the mero/cef FN tail.
- **Label quality:** S. aureus oxacillin spec 0.333 was oxacillin-AST label noise (mecA genotype was right;
  cefoxitin is the proper surrogate, substrate-sparse here).
- **Embedding frontier closed 0-for-3** (pathotype/AMR/carbon-util) — deterministic mechanism features are
  the product; embeddings need a YES/YES/YES substrate that doesn't exist at solo scale.

## Why this is the milestone (not "one more organism")

Adding more organisms/drugs now re-confirms the two findings above without new information. The decoder is
complete, validated, deployed, honest about its blind spots. The genuinely high-value next leaps are
different in kind and need resources only the user can supply:
- cross-LAB / clinical-partner validation (independent data access),
- a non-AMR sampling-independent labeled substrate at depth (acquisition/funding),
- multimodal / eukaryotic (roadmap Phase 5/6 — own `/idea-anchor` cycle + compute).

## Provenance

Per-drug + per-organism artifacts: `wiki/{dna_amr_multidrug,dna_amr_xsource,klebsiella_drug_matrix,
klebsiella_*_validate,pseudomonas_*_validate,staphylococcus_*_validate}_2026-06-*.md`. Rules:
`dna_decode/eval/amr_rules.py::DRUG_RULE`. Generalized validator: `scripts/organism_drug_validate.py`
(any NCBI organism × drug = one command). 108 tests green. Tags v0.2.0 → v0.4.0.
