# TB coordinate-alignment probe — ALIGNED (2026-06-16)

Verifies the TB-decoder plan's load-bearing assumption: the WHO catalogue v2 (pinned `0bb3914348c5`) and the cached CRyPTIC VCFs share reference **NC_000962.3** and align at the variant level. Data-only (no network). Sibling of the CRyPTIC feasibility probe.

- WHO catalogue: **30699** distinct variants, **438** grade-1/2 (Assoc w R).
- Cached CRyPTIC VCFs scanned: **30**.
- Verdict: **ALIGNED**.

## Sentinels

| sentinel | grade-1/2 | catalogue coords | VCF matches (PASS+GT≥1, exact ref/alt) | aligned |
|---|---|---|---|---|
| RIF rpoB S450L (rpoB_p.Ser450Leu) | True | 761155 C>T | 10 | **True** |
| INH katG S315T (katG_p.Ser315Thr) | True | 2155168 C>G | 10 | **True** |

## Honesty

Verifies coordinate alignment on 2 canonical sentinels against the cached VCF subset. Proves the catalogue<->VCF reference frame matches; does NOT establish cohort sens/spec or biological independence (WHO catalogue built partly from CRyPTIC -> any CRyPTIC score is a KNOWLEDGE_BASELINE).
