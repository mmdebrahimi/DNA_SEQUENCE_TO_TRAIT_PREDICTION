# Data hunt round 2 — widen census · gnomAD benign · at-scale AM holdout (2026-07-23)

**Status:** ✅ three sequential data-hunt moves (user directive). All in the R2 molecular regime (clears the 8
rejection gates by construction). Frozen AMR surface byte-unchanged (READ-only). Per-move JSON:
`wiki/clinical_gene_landscape_census_2026-07-23.json` (move 1), `wiki/clinical_gnomad_benign_census_2026-07-23.json`
(move 2), `wiki/mavedb_am_holdout_2026-07-23.json` (move 3).

## Move 1 — widen the census (ClinVar-only): confirms the benign-sparsity wall at scale

Ran the offset-auto census over 20 more clinically-relevant genes (CYP2C19, HMBS, ASPA, CPOX, RASA1, MPL, PRKN,
TYK2, kinases, SNCA, RHO). **Result: 0 new AUROC-viable genes.** Every widened gene was SINGLE_CLASS (benign-
sparse: HMBS 48P/3B, ASPA 61P/6B, PRKN 17P/14B) or NO_JOIN (offset artifacts on domain-numbered kinase assays:
RASA1 off=280, SRC off=269, LYN off=63 — the assayed domain carries few ClinVar variants). This is the honest
confirmation that **ClinVar-benign sparsity — not the join — is the binding constraint**, exactly what move 2
targets. (`scripts/clinical_gene_landscape_census.py --genes …`.)

## Move 2 — gnomAD frequency-benign (the lever): unlocks 5 new genes for the DMS-ceiling, non-circularly

Supplied the missing BENIGN class from gnomAD r4 common variants (AF ≥ 1e-4 = the ACMG frequency-benign
principle; `scripts/gnomad_benign.py`), positives = ClinVar-pathogenic. **7 genes SCORED (5 genuinely new)** on
the fair, non-circular decoders:

| Gene | **DMS-ceiling AUROC** (fitness-align) | BLOSUM floor | disease | new? |
|---|---|---|---|---|
| **KCNH2** | **0.954** | 0.397 | long-QT (cardiac) | **NEW** |
| MSH2 | 0.945 | 0.881 | Lynch | (had ClinVar) |
| **TSC2** | **0.861** | 0.680 | tuberous sclerosis | **NEW** |
| MTHFR | 0.732 | 0.633 | homocystinuria | **NEW** |
| PRKN | 0.731 | 0.692 | Parkinson | **NEW** |
| LDLR | 0.727 | 0.778 | hypercholesterolemia | (had ClinVar) |
| G6PD | 0.673 | 0.702 | G6PD deficiency | **NEW** |

**THE LOAD-BEARING CIRCULARITY (surfaced, not hidden):** AlphaMissense was TRAINED using population-common
variants as its benign weak-label, so **AM's AUROC on a gnomAD-benign set is CIRCULAR** — it is reported under
`am_auroc_CIRCULAR` (0.87–0.99, trivially high) and NEVER headlined. The fair non-circular decoders on
gnomAD-benign are the **DMS-itself ceiling** (wet-lab; never saw gnomAD) and **BLOSUM** (both above). ESM2 is
also fair (self-supervised) but GPU-slow at census scale → deferred.

**What this buys + honest limits:**
- gnomAD unlocks the **fitness-alignment CEILING** (does the molecular assay separate clinical-pathogenic from
  population-common?) on 5 constrained genes ClinVar alone could not score — confirming R2 applies broadly
  (KCNH2 0.954, TSC2 0.861 strong; G6PD 0.673 weak — likely a property mismatch, echoing the CYP2C9
  abundance-vs-activity finding).
- It does **NOT** validate the deployable AM decoder (circular on this benign source) — the deployable AM
  clinical number stays the ClinVar-only 4-gene result.
- The **most-constrained genes stay STILL_UNDERPOWERED even with gnomAD** (KRAS 0 benign, GCK 4, PTEN 6, TP53
  10) — highly-constrained oncogenes/suppressors genuinely have too few common variants. A real biological
  wall, not a code gap. Lowering to AF ≥ 1e-5 would clear them numerically but weakens the benign proxy
  (a 1-in-100k variant can be pathogenic) — NOT taken; rigor over volume.
- ClinVar-benign and gnomAD-benign are COMPLEMENTARY: TP53/F9 have curated ClinVar benign but few gnomAD-common
  (constrained); KCNH2/TSC2/G6PD/MTHFR/PRKN are the reverse.

## Move 3 — AlphaMissense on the held-out MaveDB DMS set: deployable decoder generalizes at ESM2 level

Placed AlphaMissense (deployable, no-GPU) on the same held-out task as the ESM2 prospective holdout — for each
held-out human MaveDB assay (gene NOT in ProteinGym), |Spearman| of AM pathogenicity vs the wet-lab DMS score,
offset-aligned. **AM median |Spearman| = 0.502 over 57 held-out human assays, vs ESM2's 0.492 human**
(`wiki/mavedb_full_esm2_2026-07-22`). The deployable AM decoder **matches/edges ESM2** on held-out molecular-
fitness ranking — it is a competitive deployable choice for fitness ranking, not only clinical pathogenicity.

**Leakage framing (honest):** AM is NOT sequence-held-out (trained proteome-wide, saw these sequences), but the
**DMS-fitness labels are independent of AM's training** (AM trained on pathogenicity weak-labels, not DMS
scores), so AM-vs-DMS-fitness is a fair LABEL-independent generalization test comparable to ESM2's zero-shot
number. The full **leakage-free ESM2+ProSST hybrid at 2383-assay scale remains the Kaggle follow-up** (GPU-bound).

## Net of the round

- **Move 1**: the benign-sparsity wall is real at scale (ClinVar-only widen → 0 new viable).
- **Move 2**: gnomAD breaks it for the CEILING on +5 constrained genes (KCNH2 0.954 …) but is CIRCULAR for AM;
  the most-constrained genes stay walled (a biological, not code, limit).
- **Move 3**: the deployable AM decoder generalizes to held-out MaveDB DMS at ESM2 level (0.502 ≈ 0.492).
- Binding constraint remains **labels** (benign-class availability + AM-circularity), not code or model —
  consistent with the project north star.

Reproduce: `scripts/clinical_gene_landscape_census.py` · `scripts/gnomad_benign.py` +
`scripts/clinical_gnomad_benign_census.py` · `scripts/mavedb_am_holdout.py`. Tests:
`tests/test_gnomad_benign.py` (6) + `tests/test_clinical_gene_landscape_census.py` (4). Memory:
`feedback_g2p_decoder_regime_boundary`.
