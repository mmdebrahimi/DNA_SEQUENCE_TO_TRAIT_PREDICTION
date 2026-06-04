# v0 cross-axis ExPEC rule — external validity on independent labels (2026-06-04)

**Closes** the `dna_decode/pathotype/markers.py` scope-limit "K=1 cross-axis tested in-sample on N=24."
**Method:** `scripts/expec_rule_external_validity.py` applies the SHIPPED rule (`meets_cross_axis_support`:
≥1 iron-acquisition gene AND ≥1 capsule/serum gene) to the recorded virulence-gene calls
(`Other_Vir_genes`) of the Horesh 2021 isolation-source-labelled strains. ExPEC = blood/urine =
independent of the resolver's marker rules.

## Result

| group | N | gene-called | cross-axis fires | rate | strong≥2 | rate |
|---|---|---|---|---|---|---|
| **ExPEC (independent positive)** | 1,574 | 1,209 | 640 | 0.407 (all) / **0.529 (gene-called)** | 291 | 0.185 |
| EPEC (intestinal; DEC-gated neg) | 269 | 94 | 2 | **0.007** | 0 | 0.0 |
| ETEC (intestinal; DEC-gated neg) | 183 | 11 | 0 | **0.000** | 0 | 0.0 |
| Not determined (commensal-ish neg) | 2,774 | 999 | 348 | 0.125 | 84 | 0.03 |

## Verdict — specificity holds, recall is modest (and optimistic in-sample)

- **Specificity STRONG:** the cross-axis rule fires on only 0.7% of EPEC and 0% of ETEC. It does NOT
  spuriously call intestinal/diarrheagenic strains ExPEC. (Even the 0.7% EPEC is harmless live — the
  DEC-module gate resolves LEE/eae ABOVE the ExPEC branch.)
- **Recall MODEST out-of-cohort: ~0.53 (lower bound)** vs the in-sample 0.917 on N=24. The in-sample
  number was optimistic. Two non-exclusive causes:
  1. **Detection sparsity (measurement floor):** Horesh's `Other_Vir_genes` lists notable VFs, not the
     resolver's full VF-DB from-FASTA scan → 0.53 is a LOWER bound; the live from-FASTA resolver likely
     recalls higher. 365/1,574 ExPEC have no gene calls at all (excluded from the recall denominator).
  2. **Genuine ExPEC heterogeneity:** many blood/urine isolates are "accidental"/opportunistic and lack
     the canonical iron+capsule ExPEC signature the cross-axis rule requires.

## Action taken — DOCUMENT, do not re-tune

The rule is NOT loosened. Re-tuning K or the cross-axis bar to chase Horesh recall would (a) repeat the
overfit trap the user already refused (the 0.833-over-0.917 decision), (b) optimise against
detection-limited gene calls, and (c) trade away the strong specificity. The honest move is to update the
shipped scope-limit from "in-sample N=24, recall 0.917" to the out-of-cohort numbers below.

**Updated scope-limit (canonical):** cross-axis ExPEC rule — out-of-cohort on N=1,209 independent
(isolation-source) ExPEC: **recall ≈0.53 (lower bound, detection-limited); specificity ≈0.99 vs
intestinal pathotypes (EPEC/ETEC).** The rule is conservative (high precision/specificity, modest
recall) — appropriate for a compatibility resolver with abstention.

## Provenance
`data/horesh_2021/F1_genome_metadata.csv` (Figshare 13270073). Rule reused verbatim from
`dna_decode.pathotype.expec_score.meets_cross_axis_support`. Rule-logic conformance on an
independent-LABEL population holding gene-detection constant — NOT a from-FASTA re-run.
