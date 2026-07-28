# Essentiality gold-standard LABEL wall — exhaustive negative (2026-07-28)

**Finding: free, machine-readable, GOLD-STANDARD gene-essentiality LABELS are not programmatically
fetchable from this host, across every route tried. This gates v0.1 (AUROC), E3 (learned complement),
and E4 (cross-organism transfer) — but NOT v0 (the conserved-core decoder, which is label-independent
and shipped, row 582).**

## Routes tried (all walled) — do NOT re-hunt without a new lead
| source | organism(s) | result |
|---|---|---|
| OGEE (v3.ogee.info / ogee.medgenius.info / ogeedb.embl.de) | multi | ALL hosts DOWN (000) |
| DEG (tubic.org/deg) | multi-bacterial | bulk-download dir form-gated (403); org file → HTML |
| DepMap portal API | human CRISPR | bot-verification HTML challenge (not JSON) |
| DeeplyEssential (GitHub) | 30 bacteria | repo has only the species INDEX; labels in a Google-Drive snapshot |
| UniProt REST disruption-phenotype | E. coli | reachable but 2.5% coverage + curation-biased (64% ess vs true ~7%) + noisy → NOT gold-standard |
| Figshare search | — | noisy, no clean deposit surfaced |
| GitHub raw candidates | — | 404 |

## Why E3 is the HARDEST-gated (R2 pre-bar check outcome)
E3 (the learned complement) is a **SUPERVISED** task: training a model to predict essentiality needs
essentiality labels to TRAIN on — not merely to validate. There is no validated zero-shot ESM2→essentiality
(unlike the DMS-validated zero-shot ESM2→variant-effect cell). So a Kaggle GPU run producing ESM2 embeddings
without training labels would be MOTION (unusable embeddings). E3 is blocked at its core, not just its eval.

## The unblock (a single external input — user authority/data fork)
ONE gold-standard essential-gene label file unblocks v0.1 + E3 + E4 simultaneously. Any of:
- **Keio essential-gene list** (303 E. coli genes; Baba 2006 Mol Syst Biol supplementary) — the canonical set.
- **The 248-consensus** (Goodall 2018 mBio, Keio∩PEC∩TraDIS) — highest-confidence.
- **DepMap common-essentials CSV** (exported from the DepMap browser, which bypasses the API bot-wall) — human.
- **DEG access** (fill the download form once) — multi-bacterial.

Format needed: a plain list of essential gene symbols/locus-tags (+ ideally the full gene universe or a
non-essential list). The v0 decoder already has the E. coli gene universe (4318 genes) on D:, so a bare
essential-gene list is enough to compute the real AUROC + train E3.

## Status
v0 conserved-core decoder SHIPPED + validated by composition (row 582). Everything beyond is label-gated.
This is a genuine data fork, NOT a code wall — surfaced for the user, not fabricated around.
