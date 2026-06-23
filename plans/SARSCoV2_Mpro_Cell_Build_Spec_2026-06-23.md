# SARS-CoV-2 antiviral cell — build spec (the next free-independent-label cell, 2026-06-23)

Direction-1 build. The label gate is CLEARED (`wiki/next_independent_label_cell_feasibility_2026-06-23.md`):
Stanford **CoV-RDB** is free/open/downloadable, lab-measured, variant-level — the HIVDB of coronaviruses.
This spec scopes the v0 cell, mirroring the HIV / fungal target-site engine pattern.

## v0 scope (decided from the data, not assumed)
- **Primary: nirmatrelvir × Mpro (3CLpro / nsp5).** Best-powered (115 CoV-RDB selection records) AND cleanest
  reference (Mpro = a single clean ORF1ab segment, no frameshift). This is the v0 cell.
- **Secondary (same gene, sparse): ensitrelvir / lufotrelvir / GC376 × Mpro** — fold into the Mpro catalog;
  per-drug differential is v0.1.
- **Deferred to v0.1: remdesivir × RdRp (nsp12).** Catalog is thin (~12 records: L480, V557…) AND nsp12
  spans a −1 ribosomal frameshift (NC_045512.2 13442–13468 → 13468–16236) → the reference CDS is non-trivial.
  Do it after Mpro proves the pattern.

## Data substrate (acquired / located)
- **Catalog source (IN HAND):** `data/raw/sarscov2/invitro_selection_antiviral.csv` — assembled from CoV-RDB
  `tables/invitro_selection_results/*.csv`, filtered to Mpro/RdRp/PLpro, 172 records / 126 unique
  `(gene,pos,aa,drug)`. Schema `(ref_name, rx_name, gene, position, amino_acid)`. Nirmatrelvir/Mpro = 115.
- **Independent validation labels (THE make-or-break, NOT yet pulled):** CoV-RDB `tables/rx_fold/*.csv` =
  measured fold-change, but RELATIONAL — keyed by `iso_name` + `rx_name`; the isolate→Mpro-mutation join is
  in `isolates.d/` + `isolate_mutations.d/`. Two ways to get the join: (a) the CoV-RDB **SQLite dump**
  (covdb.stanford.edu download / Stanford CDN — GitHub releases API returned empty, find the real URL) and
  one SQL query `rx_fold ⋈ isolate_mutations WHERE gene=_3CLpro AND rx=nirmatrelvir`; or (b) the public
  **CoV-RDB API**. **Powering reality to measure FIRST:** how many isolates carry an Mpro mutation AND have a
  measured nirmatrelvir fold-change. Clinical nirmatrelvir resistance is rare → expect a SMALL validation N
  (honest-underpowered, like the CAB / 2nd-gen-INSTI cells). That is a powering limit, NOT circularity.

## Engine reuse (mirror HIV / fungal — minimal new code)
- **Catalog:** `dna_decode/data/sarscov2_amr.py` — curated Mpro resistance substitutions from the CoV-RDB
  selection set (+ a `MproTargetClass`), mirroring `hiv_amr.py` (`call_target_class`, position-based v0 OR
  mutant-level). `--observed Mpro:E166V,L50F` wheel-only path.
- **Reference + caller:** commit `data/sarscov2_ref/SARSCoV2_Mpro_NC045512_cds.fna` (NC_045512.2 nsp5
  10055–10972, 306 codons) + `scripts/sarscov2_caller.py` reusing
  `fungal_erg11_caller.observed_substitutions` codon-mapper (the SAME shared mapper HIV uses). **Reference
  integrity gate** (the HIV pattern): a test self-checks the committed reference translation == catalog WT at
  every catalog position before any call.
- **CLI:** route `--drug nirmatrelvir` (+ ensitrelvir/lufotrelvir) → `_sarscov2_main`; `--organism` relabel
  `Escherichia`→`SARS-CoV-2`; genome-FASTA mode via the committed Mpro ref; same `amr-mechanism-call-v1`
  record. Offline (no BLAST) → INDETERMINATE.
- **Validation:** `scripts/sarscov2_mpro_validate.py` vs CoV-RDB measured fold-change (independent), honest
  powering gate (ACCRUING / UNDERPOWERED / SCORED). Its OWN trust surface — NAMESPACE-SEPARATE from the
  bacterial NCBI-PD card and the HIV card (in-distribution ≠ provenance-disjoint; do not conflate).

## Acceptance bar (draft — ratify)
1. `sarscov2_amr` Mpro catalog sourced from the committed CoV-RDB selection data + `--observed` path.
2. Committed Wuhan-Hu-1 Mpro reference + reference-integrity self-check test (translate==WT at catalog posns).
3. CLI `--drug nirmatrelvir` wheel-only + genome-FASTA modes; `amr-mechanism-call-v1` record.
4. Validation vs CoV-RDB measured fold-change with an HONEST powering verdict (likely UNDERPOWERED v0).
5. Tests (catalog / caller / CLI / validate); FROZEN bacterial AMR surface byte-unchanged.

## Honesty rails (carry from HIV/TB)
- CoV-RDB is in-distribution for the catalog (the catalog is partly built from it) UNLESS we validate on a
  held-out / later fold-change set — label the cell accordingly (KNOWLEDGE_BASELINE vs INDEPENDENT).
- Report powering honestly; a rare-resistance drug yields a small N — say UNDERPOWERED, never inflate.
- Per-drug differential (nirmatrelvir vs ensitrelvir) + RdRp/remdesivir + the frameshift = v0.1.

## Provenance
Feasibility `wiki/next_independent_label_cell_feasibility_2026-06-23.md`; catalog substrate
`data/raw/sarscov2/invitro_selection_antiviral.csv` (CoV-RDB `invitro_selection_results`, MIT). Ledger row 165.
