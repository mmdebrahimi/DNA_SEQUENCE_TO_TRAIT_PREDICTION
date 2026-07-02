# FinnGen summary-stats ingest — free host-genetics locus priors (2026-07-02)

The free, no-auth forward-path from the catalog plan (Track C, item 1): FinnGen publishes GWAS summary
statistics for 2,470 disease endpoints openly (no application, no fee). Fetched + a real ingest built.
**Honest scope up front: this is HOST (human) genetics — a locus-priors / feature layer banked for a possible
future human-variant arm. It does NOT feed the pathogen/organism decoder (that is pathogen genetics, a
different axis). Not load-bearing; not overclaimed.**

## What was done
- Fetched the **R12 manifest** (2,470 endpoints; `phenocode / phenotype / category / num_cases / num_controls
  / path_https`) from the public GCS bucket `finngen-public-data-r12` — free, no auth.
- Built `scripts/finngen_ingest.py`: streams a per-endpoint `.gz` (schema `#chrom pos ref alt rsids
  nearest_genes pval mlogp beta sebeta af_*`), filters to genome-wide-significant variants (p < 5e-8), clumps
  to the **lead variant per nearest gene**, emits a compact priors table (gene → lead variant, rsid, p, beta,
  risk/protective direction) + honest metadata (cases/controls, n_variants, n_gw_sig).
- Ran it on the most project-relevant endpoint: **AB1_TUBERCULOSIS** (host susceptibility to TB;
  3,063 cases / 497,285 controls).

## Result — and it recovers textbook biology
`wiki/finngen_priors_AB1_TUBERCULOSIS.json`: 21,327,062 variants → **233 genome-wide-significant** → clumped
to **3 gene loci, ALL in the HLA class-II region on chr6 (~32.3–32.5 Mb)**:

| gene | lead variant | rsid | mlogp | direction |
|---|---|---|---|---|
| TSBP1-AS1 | 6:32373621:A:G | rs9391858 | 8.66 | protective |
| HLA-DRB9 | 6:32459827:G:C | rs13207140 | 8.16 | protective |
| HLA-DRA | 6:32447985:A:G | rs12529078 | 8.12 | protective |

This is the **canonical host-TB-susceptibility signal** — the HLA-DR (immune antigen-presentation) locus — so
the ingest is validated: the free FinnGen path works and extracts real, biologically-correct loci. (233
variants collapsing to one HLA peak within ~90 kb is exactly the expected LD structure.)

## Honest scope + notes
- **HOST genetics, not the pathogen decoder.** The TB host-susceptibility HLA signal is about which *humans*
  get TB, NOT about the *M. tuberculosis* genome (which the deterministic decoder scores). Different axis.
  This is a prior/feature layer for a future human-variant arm — banked, not wired into the decoder.
- The 773 MB per-endpoint `.gz` is gitignored-class (lives on `D:/dna_decode_cache/finngen/`); only the tiny
  priors JSON + the ingest script + tests are committed. The manifest URL is documented (not committed).
- The ingest generalizes to any of the 2,470 endpoints:
  `uv run python -m scripts.finngen_ingest --gz <endpoint>.gz --endpoint <PHENOCODE> --manifest <R12_manifest>`.
- TB at 3,063 cases is modestly powered (233 GW-sig, all one HLA peak); larger endpoints (e.g.
  AB1_BACTINF_NOS, 25,506 cases) would yield richer loci if a human arm is opened.

## Verdict
Track C item 1 COMPLETE: the free FinnGen summary-stats path is real, ingested, and validated (recovers the
HLA-DR host-TB signal). Banked as an honestly-scoped host-genetics priors layer. Frozen AMR surface
byte-unchanged; 3 offline tests (`tests/test_finngen_ingest.py`).
