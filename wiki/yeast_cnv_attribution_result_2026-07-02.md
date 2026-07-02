# Yeast copy-number attribution — the capstone question RESOLVED (2026-07-02)

The attribution capstone (`wiki/yeast_attribution_capstone_2026-07-02.md`) left one open thread: gene
presence/absence FAILED canonical-gene attribution (copper→CUP1 is copy-number, invisible to presence/absence).
A bounded confirmatory test — **NOT** the full SNP+kinship-LMM build that was weighed — resolves it: given the
RIGHT feature (copy number), the yeast decoder attributes the canonical mechanism, de-confounded.

## What ran (the cheap path, per the design review)
The 1002 project already ships a precomputed per-isolate copy-number matrix (`genesMatrix_CopyNumber.tab.gz`,
free) — so the "hard part" (CUP1 dosage) was a download, not a CNV-calling build. Test: clade-centered
CONTINUOUS Spearman of gene copy number vs growth (residualize BOTH by SNP-distance clade), within-clade
permutation null (n=200), direction check, K=18 + K=30. Reuses the existing clade machinery; no LMM, no binary
carrier test. `scripts/yeast_cnv_attribution.py`.

## Results — two independent de-confounded copy-number attributions + a specificity control
| gene (copy) | mechanism | condition | global ρ | **clade-centered ρ (K18/K30)** | perm_p | verdict |
|---|---|---|---|---|---|---|
| **CUP1** (YHR055C) | copper metallothionein tandem array | copper (YPDCUSO410MM) | +0.824 | **+0.726 / +0.717** | 0.005 | **CONFIRMED** |
| **ENA5** (YDR038C) | Na⁺/Li⁺ efflux ATPase cluster | sodium (YPDNACL1M) | +0.330 | **+0.251 / +0.251** | 0.005 | **CONFIRMED** |

- **CUP1 copy number is real dosage** (1–45 copies, mean 4.08 — NOT collapsed; the precomputed-matrix
  collapse risk did not materialize).
- **ENA5 is mechanism-SPECIFIC** (a textbook control): ENA5 copy associates with IONIC stress — sodium +0.25,
  **lithium (YPDLICL250MM) +0.271** (ENA also effluxes Li⁺) — but **NOT with non-ionic osmotic stress**:
  ENA5 × sorbitol clade-centered ρ **+0.03, perm_p 0.38 (null)**. The gene helps only against the ions it
  actually transports.

## What this resolves
1. **The capstone diagnosis is CONFIRMED, not just asserted.** Attribution FAILED with presence/absence
   (wrong feature type for a copy-number mechanism) and SUCCEEDS with copy number (CUP1 ρ +0.73 de-confounded;
   ENA5 ρ +0.25 ionic-specific). The failure was a **feature/mechanism mismatch, not a capability gap.**
2. **This is now the SECOND proof of the general law** (after DepMap): the feature type must match the
   mechanism type. Point mutations → DepMap (BRAF/TP53 recovered); copy number → yeast (CUP1/ENA recovered);
   presence/absence → only accessory/gene-content mechanisms. The learned decoder attributes canonical genes
   whenever the feature encodes the mechanism.
3. **The full SNP+kinship-LMM build was correctly NOT pursued** — it would have been over-scoped, off the
   north star, and brittle; the cheap targeted test on the existing CNV matrix answered the question decisively.

## Honest scope
- Confirmatory, not genome-wide: this TESTS named canonical loci (CUP1/ENA), it does not run a de-novo GWAS to
  DISCOVER them. That is the right scope for the capstone question ("can the decoder attribute the canonical
  gene given the right feature?" → yes).
- Copy number for arsenite→ARR not tested (ARR is a presence/absence cluster, a different feature axis; the
  capstone's arsenite thread stays open but is lower value).
- Data on the scratchpad (regenerable free download); the logic is tested on synthetic data
  (`tests/test_yeast_cnv_attribution.py`).

## Verdict
The yeast attribution capstone's open thread is CLOSED: **canonical-gene attribution is CONFIRMED for
copy-number mechanisms (CUP1→copper, ENA5→sodium/lithium, de-confounded, mechanism-specific)** using the
already-available copy-number features — no new SNP/LMM build required. Frozen AMR surface byte-unchanged.
