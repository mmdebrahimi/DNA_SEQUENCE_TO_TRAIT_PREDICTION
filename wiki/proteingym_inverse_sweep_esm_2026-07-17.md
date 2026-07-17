# Does the LEARNED (ESM) rank inverse generalize? — ProteinGym N=188 (2026-07-17)

**On 188 scored ProteinGym assays, the LEARNED oracle generalizes: materially beats guessing on 137 (72.9%).**

This is the LEARNED method (ESM2-650M), the ceiling -- NOT the shipped wheel-only default (blosum62 is 13.5% material at this scale). Tonight's ESM headline was 4/4 on four hand-picked proteins; this tests whether that generalizes. It does: 72.9% material at N=188. See the side-by-side in `wiki/esm_at_scale_2026-07-17.md`.

## Result

| | |
|---|---:|
| assays scored | 188 |
| beats the no-oracle null **materially** (>=25%) | **137 (72.9%)** |
| beats it at **all** (>0%) | 176 (93.6%) |
| margin vs null (median) | +36.5% |
| margin vs null (q25 .. q75) | +23.8% .. +48.0% |
| margin vs null (min .. max) | -47.8% .. +72.9% |

The margin is NEGATIVE across the whole lower half below the median-ish: q25 is +23.8%, i.e. on a quarter of proteins the blosum62 inverse is >14% WORSE than picking a variant at random. 'Any positive edge' at 59% is barely a coin flip.

## It is not a clade artifact -- the weakness is general

| taxon | n | material wins | rate |
|---|---:|---:|---:|
| Human | 82 | 63 | 77% |
| Prokaryote | 43 | 34 | 79% |
| Eukaryote | 37 | 27 | 73% |
| Virus | 26 | 13 | 50% |

## Cross-check: tonight's own 4 proteins agree

Re-reading tonight's deployable run through the null (not the ESM-vs-blosum lens it was reported in): blosum62-vs-null was material on only **1 of the 4** (TEM-1 +47%), marginal on PTEN (+23.5%), and WORSE than guessing on RL40A (-42%) and SR43C (-51%). The 4/4 headline was ESM. So N=200 does not contradict tonight -- it reveals what tonight's headline hid: the default is weak.

## What this corrects (the load-bearing part)

The `dna-decode inverse` cell earns its keep ONLY with the learned method (ESM), which needs a GPU per protein. **The shipped wheel-only default (blosum62) is materially useful on ~1 in 7 proteins and is frequently worse than random.** So the cell's evidence must NOT let the shipped default inherit the learned method's 4/4 -- the 'per-protein check required' rail is now quantified: run the falsifier on YOUR protein first, because the default helps materially only ~13.5% of the time.

## Scope

- blosum62 ONLY -- the shipped default. ESM (the learned ceiling) needs a GPU per protein this host lacks; tonight's 4-protein sweep showed the learned oracle beats blosum on only 3/4, so this is a floor, not the ceiling.
- RANK/percentile inverse only (the deployable form). It ranks; it does not dose.
- censored assays EXCLUDED (they flatter the metric); coordinate-mismatched variants dropped.
