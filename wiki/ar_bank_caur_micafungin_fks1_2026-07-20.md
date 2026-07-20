# AR Bank C. auris micafungin (FKS1) — echinocandin arm: reference verified + S-side clean, R-side SRA-accruing

**Date:** 2026-07-20 · **Status:** ✅ **POWERED — SCORED_ENDORSED** (5R/8S, sens 0.60, spec 1.00, acc 0.846)
· **Cell:** `dna_decode/data/fungal_amr` FKS1 echinocandin catalog (NON-FROZEN) via `fungal_erg11_caller`
(gene-agnostic) · **Frozen surface:** byte-unchanged (verify_lock OK).

## POWERED RESULT (5R/8S) + the uncatalogued-variant findings

All 6 R were pursued via SRA read-mapping; 5 landed (the 6th, SRR26666804 @ 3.4 GB, deferred on host-disk
I/O — not needed for powering). The per-R FKS1 target-site calls:

| R isolate | FKS1 variant | catalog call | interpretation |
|---|---|---|---|
| SAMN38094192 | **S639P** | R ✓ | canonical HS1 — caught |
| SAMN38094212 | **S639Y** | R ✓ | canonical HS1 — caught |
| SAMN38094223 | **S639F** | R ✓ | canonical HS1 — caught |
| SAMN38094197 | **F635C** | S (miss) | **substitution at a CATALOGUED HOTSPOT position** (F635); catalog lists only `F635del` → defensible catalog-completion candidate |
| SAMN38094230 | **W691L** | S (miss) | non-hotspot position; only 1 R carries it → isolate-specific, NOT a systematic blind spot |

**Determinant-scan TRANSFERS to the echinocandin target: sens 0.60, spec 1.00 (8/8 S correct), POWERED,
SCORED_ENDORSED.** 3/5 R caught cleanly via canonical S639 substitutions (S639F/P/Y) — the catalog works on
the dominant FKS1 HS1 mechanism.

**Two misses, honestly distinct** (both DISCLOSED, neither auto-added — the ERG11 clade-IV over-call rail):
- **F635C** (SAMN38094197) is the stronger finding: F635 IS a documented HS1 hotspot, but the catalog only
  encodes the *deletion* `F635del`, so the *substitution* F635C slips through. Treating any F635 substitution
  as resistance (not just the del) would lift sens 0.60 → 0.80 — a defensible catalog completion, flagged for
  independent review, NOT auto-applied on n=1.
- **W691L** (SAMN38094230) is a non-hotspot variant carried by a single R → the "candidate blind spot" from
  the n=1 partial resolves as **isolate-specific**, not a mechanism the catalog systematically misses.

The distinction is the whole point of the uncatalogued-variant disclosure: it surfaced BOTH, and let the
hotspot-vs-non-hotspot position separate a real catalog gap (F635C) from an isolate quirk (W691L) — instead
of silently scoring both R as S.

## Why this arm

The fungal cell's azole/ERG11 arm is externally validated (fluconazole POWERED, `ar_bank_caur_powered_result_2026-07-20.md`).
**Echinocandin/FKS1 is a genuinely different drug class AND target gene** — validating it tests whether the
determinant-scan method generalizes to a second fungal mechanism, not just a second azole. AR Bank C. auris
has micafungin MIC labels: **6R / 26S** (CDC tentative breakpoint micafungin ≥4 → R).

## The load-bearing derisk — FKS1 numbering VERIFIED (R2 derive-don't-assert)

The catalog encodes C. auris FKS1 hotspots **S639F/P/Y (HS1), F635del, R1354H (HS2)**. A wrong reference
frame/coordinate would silently mis-number every echinocandin call. Verified against the real B8441
reference (GCF_002759435.1), where **FKS1 is annotated `GSC1`** (the S. cerevisiae name; locus B9J08_02922,
1888 aa — a plain "FKS1" grep misses it):

| position | reference residue | catalog expects |
|---|---|---|
| 635 | **F** | F635del ✓ |
| 639 | **S** | S639F/P/Y ✓ |
| 1354 | **R** | R1354H ✓ |

All three WT residues match exactly. Committed reference `data/fungal_ref/Cauris_FKS1_cds.fna` (5667 nt,
in-frame); self-consistency confirmed (BLAST of the FKS1 ref vs its own B8441 genome → zero substitutions =
WT). Pinned by `tests/test_fungal_fks1_reference.py` (3 tests: exists+in-frame, WT residues, catalog
coverage) — a future frame error fails loudly.

## Result so far

| view | n | R/S | sens | spec | acc | verdict |
|---|---|---|---|---|---|---|
| **S-side (assembled)** | 8 | 0R/8S | — | **1.00** | **1.00** | clean; R-side pending |

All 8 micafungin-S isolates that have a downloadable assembly (overlap the azole-assembled set) carry **WT
FKS1** (zero substitutions) → called S. Clean specificity, 0 errors — the FKS1 caller works end-to-end on
the real genomes.

## Why R-side is SRA-accruing (data-availability ceiling, not a cell gap)

All **6 micafungin-R** isolates are **SRA-reads-only** (0 downloadable assemblies) — the same AR-Bank
Gram-positive/fungal assembly ceiling seen on the azole arm. R-side is being powered by targeted FKS1
read-mapping (`assemble_sra_cohort --method map --erg11-ref Cauris_FKS1_cds.fna`, 6 R resolved to SRR
runs SRR26666757/776/783/804/795/808). When the FKS1 consensuses land, re-run the scorer to fold them in;
expect FKS1 HS1 mutations (S639F/P/Y) in the R isolates → R, giving a powered two-sided result.

## Update (R-side accruing) — uncatalogued-variant disclosure + a candidate blind spot

`scripts/ar_bank_caur_micafungin_finalize.py` combines the assembled S + SRA-mapped R and — applying the
ERG11 clade-IV lineage-marker lesson from earlier today to the echinocandin target — records the **raw**
FKS1 substitutions (catalogued + uncatalogued) for BOTH R and S, then classifies each uncatalogued variant:

- **R-only** uncatalogued variant → `CANDIDATE_BLIND_SPOT` (possible determinant OR R-lineage marker; needs
  independent evidence — **NOT auto-added to the catalog**, exactly the over-call the ERG11 lesson warns against).
- present in any **S** isolate → `BENIGN_POLYMORPHISM` / `LINEAGE_MARKER` (non-discriminative, not causal).

Current partial (1R/8S — 1 of 6 R landed, rest SRA-accruing): spec **1.00** (8/8 S correct), sens 0.0 (the
one R scored S). The lone R (**SAMN38094230**, MIC>8) carries FKS1 **W691L** — a clean, full-length,
non-canonical variant (S639/F635/R1354 all WT; 0 consensus gaps) that the catalog does not list → flagged
`CANDIDATE_BLIND_SPOT`. Two S-only variants (L1572I, I1817V) are correctly flagged `BENIGN_POLYMORPHISM`.
**Verdict deferred to the full R set:** if the other 5 R carry canonical S639F/P/Y the catalog works and
W691L is an isolate-specific quirk; if they share W691L / other non-hotspot variants, that is a genuine FKS1
catalog blind spot (the echinocandin analogue of the distributed-mechanism tet failure). Not concluded on
n=1.

## Scope

- CURATED_NONFROZEN (fungal FKS1 echinocandin catalog); genotype = BLAST FKS1 (annotated GSC1) target-site;
  phenotype = CDC micafungin MIC; provenance-disjoint. The caller is gene-agnostic
  (`observed_substitutions(genome, cds_ref, gene="FKS1")`); the scorer `ar_bank_caur_validate.py` is now
  generalized (`--drug/--gene/--cds-ref`, defaults fluconazole/ERG11).
- Artifact: `wiki/ar_bank_caur_validation_micafungin_*.json`. Frozen surface `verify_lock` OK.
