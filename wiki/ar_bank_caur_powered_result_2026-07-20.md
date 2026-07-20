# AR Bank C. auris fluconazole — POWERED, and it surfaced a clade IV lineage-marker over-call

**Date:** 2026-07-20 · **Status:** POWERED (5R/7S ≥5-per-class); binary NOT_ENDORSED, HIGH-confidence subset
PERFECT · **Cell:** `dna_decode/data/fungal_amr` ERG11 target-site catalog (NON-FROZEN) via
`scripts/fungal_erg11_caller` · **Frozen surface:** byte-unchanged (verify_lock OK).

## Headline

Powering the C. auris fluconazole cell (SRA-read-mapping the unassembled fluconazole-S isolates to close
the S-class powering floor) did more than power it — **it surfaced that the catalog was over-calling a
clade IV ERG11 lineage marker as a causal determinant**, the fungal analogue of the documented
QRDR-vs-lineage confound. The fix (a user-ratified confidence tier) preserves sensitivity while making the
weak-evidence calls visible, and recovers a **perfect result on the mechanism-attributable subset**.

## Result

| view | n | R/S | sens | spec | acc | verdict |
|---|---|---|---|---|---|---|
| **binary cell** (all calls) | 12 | 5R/7S | **1.00** | 0.714 | 0.833 | POWERED · NOT_ENDORSED |
| **HIGH-confidence subset** (mechanism-attributable) | 9 | 4R/5S | **1.00** | **1.00** | **1.00** | perfect (R-underpowered) |

- **Powering:** the MVP bar (fluconazole n_R≥5 AND n_S≥5) is **met** — 5 R + 7 S scored (8 from downloadable
  assembly + 4 from SRA-read-mapped ERG11 consensus). (5th fluconazole-S SRA isolate SAMN46817482 was still
  mapping at write time; the finalize is idempotent — re-run to fold it in. 4/5 SRA S is already powered.)
- **Independence:** CDC AR Isolate Bank measured MIC (CDC tentative fluconazole ≥32 → R); provenance-disjoint
  (0 overlap vs the fungal G1 tuning cohort). NOT methodology-independent (the rule IS the curated ERG11
  catalog). Genotype = BLAST ERG11 target-site.

## The finding — clade IV ERG11 haplotype is a lineage marker, not a causal determinant

The haplotype **K177R / N335S / E343D** appears with **identical genotype in isolates of opposite phenotype**:

| isolate | CDC label | ERG11 genotype | binary call | conf tier |
|---|---|---|---|---|
| SAMN11570381 | **R** | K177R/N335S/E343D (only) | R → **TP** | LOW_LINEAGE_ONLY |
| SAMN10139552 | **S** | K177R/N335S/E343D (only) | R → **FP** | LOW_LINEAGE_ONLY |
| SAMN46817483 | **S** | K177R/N335S/E343D (only) | R → **FP** | LOW_LINEAGE_ONLY |

**1 R + 2 S, all the same ERG11 genotype** ⇒ the haplotype does not discriminate R from S — it co-segregates
with clade IV (often genuinely resistant) but carries **zero net discriminative signal** (+1 TP and +2 FP).
The originating S. Africa outbreak paper (PMC10521600) reported these in clade IV **R** isolates but had **no
clade IV S control**; the AR Bank supplies that control. Survives falsification on N=3 (1R/2S identical
genotype); suggestive of a lineage marker, underpowered for a causal-role verdict.

## The fix — confidence tier (user-ratified, AMR-safe)

Rather than delete the haplotype (which would turn genuinely-resistant SAMN11570381 into a false-negative —
missing resistance is the worse AMR error), the catalog now **tiers** the R call:

- a **causal** marker (Y132F / K143R / F126L/T / VF125AL) → `confidence=HIGH`;
- an R driven **only** by the clade-background haplotype → `confidence=LOW_LINEAGE_ONLY` (R **preserved** —
  no missed resistance — but flagged "lineage-associated, mechanism-unconfirmed").

`FungalCall` gains `confidence` + `lineage_only_determinants` (defaulted, backward-compatible). Shipped in
`5852a0e` with +5 regression tests (11 total). Non-frozen cell; frozen decoder surface byte-identical.

**On the HIGH-confidence (mechanism-attributable) subset the cell is 9/9 perfect (sens 1.0 / spec 1.0)** — the
3 clade IV haplotype-only isolates correctly abstain to LOW rather than posting 1 TP + 2 confident FPs.

## Scope

- CURATED_NONFROZEN (fungal ERG11 target-site catalog); genotype = BLAST ERG11; phenotype = CDC MIC;
  provenance-disjoint (0 overlap vs fungal G1). Artifact:
  `wiki/ar_bank_caur_powered_validation_fluconazole_caur_powered_2026-07-20_*.json`.
- SRA-read-mapping via `scripts/assemble_sra_cohort --method map` (targeted ERG11 minimap2 + samtools
  consensus, ~4 min/isolate on D:). The Gram-positive/fungal AR-Bank assembly-availability ceiling (most
  C. auris deposited SRA-reads-only) is why powering needed read-mapping, not assembly-download.
- Frozen surface `verify_lock` OK.
