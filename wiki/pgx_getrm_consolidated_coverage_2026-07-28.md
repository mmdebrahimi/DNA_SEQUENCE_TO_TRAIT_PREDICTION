# Full GeT-RM ingest — coverage gain + "are the other cells useful?" (2026-07-28)

**Yes — emphatically.** The full GeT-RM Consolidated table (363 samples) + 137-PGx-sample table (137 + ENA
sequencing URLs) ingested to `tests/data/pgx_getrm/getrm_consolidated_truth.tsv` (**4,193 truth rows / 333
samples / 37 genes/loci**, public-domain). Parser: `scripts/getrm_consolidated_ingest.py`. This is a much
bigger unlock than "more samples for genes we already do" — it splits into four tiers of value.

## Tier A — SCORE cells we already BUILT but never validated (highest VOI, no new decoder needed)

We have working decoder cells for these but **zero concordance number** — GeT-RM now gives free truth:

| gene | GeT-RM samples | today |
|---|---|---|
| **UGT1A1** | **197** | cell exists, unscored → irinotecan dosing |
| **SLCO1B1** | **137** | cell exists, unscored → statin myopathy risk |
| **CYP4F2** | **96** | cell exists, unscored → warfarin dose |

Scoring these = **3 new validated cells from code we already shipped.** Pure win.

## Tier B — FIX an underpowered cell (CYP3A5: n=8 → 137)

CYP3A5 is currently flagged **UNDERPOWERED** (only ~8 GeT-RM samples overlapped our source). GeT-RM has
**137**. This can take CYP3A5 (tacrolimus dosing) from "populated but unvalidatable" to properly scored.

## Tier C — EXPAND N on already-scored genes (more validation confidence)

| gene | GeT-RM samples | we currently score on |
|---|---|---|
| CYP2D6 | 272 | ~54 SNP-decodable |
| CYP2C9 / CYP2C19 | 248 each | ~87 |
| TPMT | 147 | ~85–98 |
| CYP2B6 / CYP2C8 | 137 each | ~87–88 |

More samples per gene → tighter, more credible numbers. (Caveat below.)

## Tier D — ROADMAP: 15 new pharmacogenes with FREE truth (no decoder yet)

Each is a candidate for a NEW cell that would be **free to validate the moment it's built** —
**CYP3A4 (164)**, CYP1A2, CYP2A6, CYP2E1, NAT1/NAT2, the GST family (GSTM1/P1/T1), the SLC transporters
(SLC15A2/SLC22A2/SLCO2B1), the UGT2B family (UGT2B7/2B15/2B17). Plus **11 HLA loci** (HLA-A/B/C/DRB1/DQ…)
that could validate/expand our HLA cell. This is the prioritized build list for future PGx cells.

## The load-bearing caveat (honest — R2)

GeT-RM truth is a **star-allele diplotype** per sample; to SCORE a cell we also need that sample's **genotype
from our source**. The 137-table ships **ENA sequencing URLs** (fastq/CRAM) per sample, so the genotypes ARE
fetchable — via the read-level CRAM tool (`scripts/pgx_cram_genotype.py`, shipped today) or the 1000G VCFs.
So "more truth samples" only converts to "more scored samples" for the subset whose genotype we can fetch +
whose sites the panel covers. The upside is real but each cell needs a fetch+score pass (the next `--until-mvp`),
not an instant number.

## Recommended next moves (ranked)

1. **Score UGT1A1 + SLCO1B1 + CYP4F2** (Tier A) — 3 shipped-but-unscored cells → validated, using existing code.
2. **Re-score CYP3A5** with the 137-sample truth (Tier B) — clear the UNDERPOWERED flag.
3. Expand N on the scored genes (Tier C) — confidence.
4. Pick 1–2 roadmap genes to BUILD (Tier D) — CYP3A4 is the highest-coverage new candidate.

Files: `scripts/getrm_consolidated_ingest.py` · `tests/data/pgx_getrm/getrm_consolidated_truth.tsv` ·
`wiki/pgx_getrm_consolidated_coverage_2026-07-28.json`. Frozen AMR/forward surfaces byte-unchanged (data +
parser only). Source xlsx are user browser-downloads (CDC public-domain), not committed.
