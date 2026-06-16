# CRyPTIC feasibility probe — GREEN: free, all-gates-passing TB measured-MIC substrate (2026-06-16)

First executor step on the **acquisition** strategic fork. The 2026-06-16 `/research` shortlist
(`research_outputs/acquirable-label-sources-2026-06-16.md`) ranked CRyPTIC #1; this probe tests it end to
end **without building the full TB decoder or any real-world acquisition.** Verdict: **GREEN — proceed.**

## What was tested (scripts/cryptic_feasibility_probe.py)
Downloaded the CRyPTIC reuse table (`data/raw/cryptic/CRyPTIC_reuse_table_20240917.csv`, 4.7 MB, 12,287
isolates, free EBI FTP, gitignored). Four feasibility axes:

### 1. POWERING — ALL 13 drugs clear ≥20/class (R/S at HIGH phenotype quality)
| drug | R_high | S_high | drug | R_high | S_high |
|---|---|---|---|---|---|
| isoniazid | **4467** | 5051 | kanamycin | 759 | 8573 |
| rifampicin | **3448** | 5507 | amikacin | 597 | 8375 |
| rifabutin | 3195 | 6846 | clofazimine | 106 | 7656 |
| ethambutol | 1460 | 5154 | delamanid | 85 | 8009 |
| ethionamide | 1189 | 7062 | linezolid | 76 | 7064 |
| levofloxacin | 1178 | 6595 | bedaquiline | 71 | 8464 |
| moxifloxacin | 991 | 5793 | | | |

This **dwarfs the E. coli substrate** — the project's deepest free-public source, where most drugs were
underpowered for provenance-disjoint cohorts. CRyPTIC gives thousands of R per first-line drug.

### 2. GENOTYPE ACCESS — clean, no assembly needed
Per-isolate VCF (variant calls vs H37Rv NC_000962.3) is fetchable from the FTP + parses. The genotype is
the **pre-computed variant call** — so the TB determinant rule reads VCF positions/codons directly, with
**no assembly / Docker / AMRFinder step** (cleaner plumbing than the E. coli path).

### 3. DETERMINANT-RULE PROOF-OF-CONCEPT (rifampicin) — the rule transfers
PoC rule: predict R iff the VCF carries a variant in the rpoB RRDR (NC_000962.3:761055-761165, the
canonical RIF region). On a balanced 15R/15S HIGH-quality sample: **sens 1.0 / spec 1.0** (TP15/FP0/TN15/FN0,
0 VCF-fetch misses). The deterministic POINT-mutation rule — the decoder's strongest regime (cf. cipro QRDR
0.925, Klebsiella 1.0) — works on TB out of the box.

### 4. CLONALITY proxy — 11 collection sites
11 distinct collection sites (top 2643 / 1828 / 1631 …). Coarse but a real provenance spread.

## Honest caveats (what the probe does NOT yet establish)
- **The PoC number is not a validated metric.** N=30, first-in-file sample, and a WINDOW proxy (any RRDR
  variant). The real rule needs the **WHO TB mutation catalogue** (codon-level; excludes synonymous/benign
  RRDR variants that would be FPs, and includes out-of-window R determinants like rpoB I491F that the window
  misses → FNs). Expect the honest full-cohort number to be below 1.0/1.0.
- **TB is globally clonal (lineages 1–4).** Raw N is huge, but the project's Mash/SNP clonality correction
  is ESSENTIAL here — effective lineages will be far fewer than 12,287, and the lineage-disclosure layer must
  run before any cell is called validated. The 11-site spread is only a coarse proxy.
- **New organism = a modeling build, not an acquisition.** Needs a TB determinant catalog (WHO catalogue /
  tb-profiler / AMRFinder TB) + an organism route in the decoder. This is the genome→phenotype kingdom-jump
  the project already proved on fungal C. auris — engineering, not procurement.
- **This expands the DETERMINISTIC decoder, not the learned-embedding niche.** TB AMR has a curated catalog,
  so (per the embedding-niche three-part test) embeddings would lose to the determinant rule here too. If the
  goal is reopening the *learned* decoder, CRyPTIC is not it (that still needs a catalog-free phenotype).

## Verdict + recommended next step
**GREEN.** CRyPTIC clears all 4 gates, is free + immediately fetchable (no MTA), genotype-accessible, and
the determinant rule transfers. The acquisition fork's best target needs **no acquisition action** — just a
decode build. Recommended next step (a real plan, not a probe): **build the TB determinant rule from the WHO
catalogue + score rifampicin (then isoniazid) on a clonality-corrected CRyPTIC cohort** — the first TB cell,
extending the decoder to a new organism in its strongest mechanism class.

## Artifacts
- Probe: `scripts/cryptic_feasibility_probe.py` (+ `tests/test_cryptic_feasibility_probe.py`).
- Machine-readable: `wiki/cryptic_feasibility_probe_2026-06-16.json`.
- Raw (gitignored): `data/raw/cryptic/CRyPTIC_reuse_table_20240917.csv`.
- Upstream research: `research_outputs/acquirable-label-sources-2026-06-16.md`.
