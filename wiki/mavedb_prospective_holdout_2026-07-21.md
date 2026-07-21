# MaveDB prospective-holdout harness — R2 molecular cell (2026-07-21)

**Status:** ✅ harness BUILT + proven end-to-end on the real MaveDB API + real ProteinGym catalog + real CPU
scoring. The R2 analog of the AMR prospective-lock: a leakage-free held-out validation set for the frozen
`forward` hybrid. `scripts/mavedb_prospective_holdout.py` · frozen surface byte-unchanged.

## What it is
MaveDB deep-mutational-scanning score sets whose target gene is NOT in the ProteinGym v1.1 benchmark (the set
the frozen hybrid was tuned + validated on) are HELD OUT by construction. R2 has NO population-structure/
clonality confound (a designed mutant library has no lineage axis), so benchmark overlap is the only leakage
vector — closed by a conservative gene-symbol dedup vs the committed `proteingym_v1.1_substitutions_catalog.tsv`.

## Real run (in-session, real network + CPU)
- ProteinGym benchmark gene symbols: **178**.
- MaveDB fetched: 100 (the search endpoint's page cap — full manifest needs pagination, a follow-up).
- **HELD OUT (protein_coding, not in ProteinGym): 86** — 58 *Homo sapiens*, 37 published ≥2024.
- BLOSUM62 scoring proof (CPU; `|Spearman|` — see caveat) on 5 held-out human assays:

| URN | gene | n missense | \|ρ\| (BLOSUM62) |
|---|---|---|---|
| 00000062-b-1 | CYP2C19 | 121 | 0.358 |
| 00000062-a-1 | CYP2C9 | 106 | 0.290 |
| 00000107-a-1 | PSAT1 | 197 | 0.296 |
| 00000107-b-1 | PSAT1 | 1914 | 0.332 |
| 00000113-a-3 | APP | 751 | 0.192 |

Sane BLOSUM-vs-DMS magnitudes (conservativeness weakly tracks functional effect, as the project's own
BLOSUM-baseline findings predict). The ESM2+ProSST hybrid on this same held-out set is the **Kaggle GPU
follow-up** (the forward cell's established pattern; expected to substantially beat BLOSUM ~0.3, given the
hybrid's ProteinGym median ~0.5+). CYP2C19/CYP2C9 are pharmacogenomics genes — this holdout also seeds the R1
human catalog lane.

## Honest caveats
- **Score-direction:** MaveDB does NOT standardize functional-score direction per assay (higher can mean more
  OR less fit — the curation ProteinGym adds). So a raw per-assay Spearman SIGN is uninterpretable without the
  assay's direction metadata; `|Spearman|` is a HARNESS PROOF only. A real headline requires per-assay
  direction normalization (from the MaveDB record's methodText / a sign heuristic vs a reference).
- **Page cap:** 100-record sample; the full held-out manifest needs search pagination.
- Scoring proof is the CPU BLOSUM62 baseline; the deployed hybrid is GPU (Kaggle).

## Next step
Kaggle run: score the 86-assay held-out manifest with the frozen ESM2+ProSST hybrid + per-assay
direction-normalized Spearman → the first prospective (leakage-free-by-benchmark-exclusion) number for the R2
molecular cell. Scoping context: `wiki/mavedb_forward_cell_scoping_2026-07-21.md`.
