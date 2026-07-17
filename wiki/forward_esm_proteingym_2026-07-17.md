# Is the shipped forward default (blosum62) 'modest' at scale? — ProteinGym N=194 (2026-07-17)

**Median |Spearman| = 0.493 (mean 0.439) across 194 assays.** The forward cell quotes blosum62 on two proteins (TEM-1 0.3465, PTEN 0.182) and calls it 'REAL but modest'. At scale 'modest' is accurate for the MEDIAN -- but the cell LEADS with TEM-1's 0.35, which is **top-13% (153/194 proteins reach 0.30), not typical**. PTEN's 0.18 is near the median.

| |Spearman| | value |
|---|---:|
| min | 0.012 |
| q25 | 0.328 |
| **median** | **0.493** |
| q75 | 0.588 |
| max | 0.738 |
| below 0.15 (near useless) | 21/194 |
| at/above 0.30 | 153/194 |

So the shipped wheel-only forward default is a weak-to-modest ranker on the typical protein; ESM (the learned method, ~0.49 median on ProteinGym, GPU) is what makes the forward cell strong, which the cell already states. This is not an overclaim correction -- the cell says 'modest' -- but it recontextualises the headline 0.35 as top-decile rather than representative.

blosum62 is a signed severity score; the cell reports |rho| per protein. This measures the shipped default's rank quality at scale -- ESM (the learned method) is separately ~0.49 median on ProteinGym and beats blosum everywhere (the cell already says so).
