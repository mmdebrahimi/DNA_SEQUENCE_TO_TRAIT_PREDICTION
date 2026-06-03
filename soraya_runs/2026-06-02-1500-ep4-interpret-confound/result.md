<!-- result.md for 2026-06-02-1500-ep4-interpret-confound -->
# Result — EP-4 model interpretation + study-confound probe

Artifact: `research_outputs/pathotype_model_interpret_confound_2026-06-02.json`.

## Finding 1 — the confound is STRUCTURAL on this substrate
`class × study = {(ExPEC, Salipante): 12, (EPEC, Hazen): 12}` → **study == class exactly**. Biology and study/assembler batch are the SAME partition here. No feature can separate them. The 0.729 (and every AUROC on this 24-genome subset) **cannot be attributed to pathotype biology** — it is pathotype-OR-batch, indistinguishable.

## Finding 2 — the signal is SPARSE (concentration math)
Full presence/absence LOSO = 0.736 (matches the diagnostic's 0.729). Only **43 / 65,536 k-mers carry any importance**; 9 reach 50%, 16 reach 80%. Within-fold top-K LOSO is monotone-decreasing — fewer k-mers do BETTER:

| top-K (within-fold selected) | LOSO AUROC |
|---|---|
| 10 | **0.806** |
| 25 | 0.771 |
| 50 | 0.757 |
| 100 | 0.750 |
| 250 | 0.729 |
| full (~70k) | 0.736 |

So ~10 k-mers carry the whole signal. The script's auto-verdict on concentration alone: `SPARSE → gene-presence-like → hopeful`.

## Finding 3 — ANALYST OVERRIDE: sparse, but the features look like BATCH, not biology
The top-10 discriminative k-mers, with per-class presence (/12):

```
TCTAGGCA imp .106  ExPEC 12/12  EPEC 6/12      CTTCTAGG imp .052  11/12  6/12
CCTAGGCA imp .070   2/12        6/12           CCTAGACG imp .047   9/12  6/12
GGGTCTAG imp .065  11/12        7/12           CTTCTAGA imp .047  12/12  7/12
TCTAGGAC imp .054  12/12        6/12           TACCTAGA imp .046   8/12  11/12
GATCTAGA imp .052  10/12        6/12           CCTAGGAC imp .046   2/12  5/12
```

**Every one contains a CTAG core** (CTAGG / TCTAGA / CCTAGG = XbaI / AvrII restriction-site motifs). **CTAG is the most under-represented tetranucleotide in E. coli** → these are RARE 8-mers, and rare-k-mer presence/absence is dominated by assembly completeness / coverage / assembler version, NOT gene content. Presence is **partial** (the "absent" class still has 5–7/12), unlike the clean ~12/0 of a discrete virulence locus. With study==class and Salipante/Hazen using different assembly pipelines, the **most parsimonious explanation is assembler/sequencing batch artifact, not pathotype virulence biology.**

→ Override the auto-"hopeful": **the learned k-mer signal (incl. the 0.729) is NOT trustworthy as pathotype biology on this substrate.** This is on the WORRYING side, and it reinforces — does not weaken — the case for abandoning the learned track in favor of the v0 known-marker resolver.

## Net for the project
The whole learned-representation bake-off (k-mer 0.514/0.604/0.729 + NT 0.38) was run on a **study==class** substrate, so none of those numbers can be read as biology. The honest, north-star-aligned path is the ledger-LOCKED **v0 deterministic virulence-gene-cluster resolver** (eae/bfp/stx/LT-ST/AAF-aggR/papC-afa…): it searches for SPECIFIC known biology, so it is immune to "whatever separates the two assembly batches." See `recommendation.md`.
