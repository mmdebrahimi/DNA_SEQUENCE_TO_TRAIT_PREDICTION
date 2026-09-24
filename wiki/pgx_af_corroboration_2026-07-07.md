# Multi-gene PGx AF-corroboration (2026-09-24)

_KNOWLEDGE_BASELINE AF-corroboration (like DPYD/VKORC1/SLCO1B1) — confirms each new cell's actionable variant is real + correctly-positioned at CPIC-expected population frequency; NOT an independent per-sample concordance number. NOT clinical._

**Verdict: AF_CORROBORATED** — 6/6 actionable variants in band across ABCG2, CYP4F2, NUDT15, UGT1A1.

- **GeT-RM concordance:** PARTLY SUPERSEDED 2026-09-22/24: the 'EXTERNAL_WALL for all' claim this field used to make is FALSE for UGT1A1 and CYP4F2 -- both are now scored against the GeT-RM CDC CONSOLIDATED multi-lab consensus (n=65 and n=64; wiki/ugt1a1_getrm_concordance_2026-09-22, wiki/cyp4f2_getrm_concordance_2026-09-24). The wall was the CYP-only ursaPGx benchmark, not GeT-RM as such. Still EXTERNAL_WALL for NUDT15 (n=12 truth overlap, underpowered) and ABCG2 (absent from the consolidated table)

| gene | allele | rsid | ALT | pop | ALT AF | expected band | verdict |
|---|---|---|---|---|---|---|---|
| NUDT15 | *3 | rs116855232 | T | EAS | 0.0952 | 0.05-0.15 | IN_BAND |
| UGT1A1 | *80 | rs887829 | T | EUR | 0.2982 | 0.22-0.42 | IN_BAND |
| UGT1A1 | *6 | rs4148323 | A | EAS | 0.1379 | 0.08-0.2 | IN_BAND |
| CYP4F2 | *3 | rs2108622 | T | EUR | 0.2903 | 0.2-0.4 | IN_BAND |
| CYP4F2 | *3 | rs2108622 | T | EAS | 0.2128 | 0.15-0.3 | IN_BAND |
| ABCG2 | 141K | rs2231142 | T | EAS | 0.2907 | 0.2-0.4 | IN_BAND |

_AFs: Ensembl 1000G phase3 (ALT = functional allele, per population), 2026-07-07. Expected bands: CPIC allele-functionality + gnomAD. NOT a clinical tool._
