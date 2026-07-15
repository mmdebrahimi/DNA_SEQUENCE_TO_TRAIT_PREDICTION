# Pre-resistance "one-nt-from-resistance" base-rate census (2026-07-15)

The `/idea-validation-council` verdict on the **escape / pre-resistance forecaster** idea was **Pursue test**:
four of six lenses converged on one crux — catalogued point-mutations are single-nt-accessible from wildtype
BY CONSTRUCTION, so the tool's "self-validation" is circular and it only adds information beyond the catalogue
if REAL genomes carry a non-WT, non-resistant **intermediate** codon that sits **strictly closer** to a
catalogued resistance codon than wildtype does. This census runs that test. Deterministic, offline, no GPU;
frozen decoder surface byte-unchanged. `scripts/pre_resistance_base_rate_census.py`.

**Substrate:** HIV — the cleanest available (mutant-level WT→R catalogues NNRTI-RT + CAI-capsid; committed
HXB2 reference CDS per gene; position-columnar per-isolate amino-acid panels in `data/raw/hiv/*_DataSet.txt`).

## Result

| Part | Metric | Value | Council expectation |
|---|---|---:|---|
| **A** accessibility | catalogued DRMs single-nt from resistance | **19/28 = 67.9%** | "≥90% ⇒ self-validation vacuous" |
| **A** | multi-step DRMs (WT ≥2 nt away) | **9/28** | assumed ~0 |
| **B** intermediate-carriers | isolates carrying a *strictly-closer, non-resistant* residue | **231** | "~0 ⇒ NO-GO" |
| **B** | broad (any non-WT, non-resistant residue at a DRM position) | 441 | — |
| — | **verdict at the ≥5 bar** | **GO** | leaned NO-GO |

Multi-step DRMs: `G190Q, G190S, K101P, V106M, Y181A, Y181I, Y188L, K70H, K70S`. Biggest carrier positions:
RT101 (141 isolates one nt from K101P), RT106 (85), RT181 (15), RT190 (3), RT188 (2), CA70 (0).

## The finding (refines the council, does not simply confirm it)

1. **Self-validation is only PARTIALLY vacuous.** On this curated *major*-DRM set, only 68% of DRMs are
   single-nt from resistance — a full **third are genuinely 2-nt (multi-step)**. The council's "≥90%
   single-nt" was a general intuition; the real major-NNRTI/CAI catalogue is enriched for multi-step
   mutations, which is exactly where a "closer intermediate" has room to exist. That 2-nt gap *is* the
   published **genetic-barrier** signal.
2. **The substrate is NOT empty.** Hundreds of real HIVDB isolates carry a non-WT, non-resistant residue that
   sits strictly closer to a catalogued resistance allele than wildtype (231 after excluding every
   already-resistant residue at the position). The `~0 → NO-GO` branch did **not** fire.

## Honest scope — what this census does and does NOT establish

- **Does establish:** there is a non-trivial substrate to flag — the forecaster's core claim ("some genomes
  are one nt from a catalogued resistance allele where wildtype is two") is real on multi-step DRMs.
- **Upper bound (load-bearing caveat):** the panels are amino-acid-level, so Part B counts a carrier if ANY
  codon of the observed residue is one nt from resistance — the isolate's *exact* codon is unknown. 231 is an
  UPPER BOUND; codon-exact confirmation needs nucleotide data (TB VCFs on `D:` / cipro N=147 QRDR) — the
  named follow-up. An upper bound of 231 (≫5) still clears the `≥5 vs ~0` decision decisively on the GO side.
- **Does NOT touch the council's other objections** (untouched by this census, still live):
  - **Novelty (Competitor lens):** "genetic barrier to resistance" is a 15-year published metric (Kliemann
    2016, Maïga 2009); "pre-resistance in *M. tuberculosis*" is a named 2021 *Nat Commun* paper (Torres Ortiz).
    A non-empty substrate does not make the concept novel.
  - **Actionability (Contrarian/Outsider):** no surveillance workflow is keyed to "proximity"; the actionable
    unit is confirmed resistance + transmission.
  - **Phenotype refinement (not run):** a cleaner "pre-resistance in a *susceptible* genome" count would
    restrict carriers to low-fold-change isolates (the panels carry fold-change columns) — deferred.

## Verdict

**Census: GO** on its narrow question (the substrate exists; the `~0` NO-GO did not fire). This does **not**
auto-authorize building the forecaster — that remains a strategic go/no-go the council's novelty +
actionability objections bear on. The census converts "is there anything to flag?" from *assumed-empty* to
*measured-non-empty (231 upper-bound carriers, on the 32% multi-step DRM fraction)*.

Artifact: `wiki/pre_resistance_base_rate_census_2026-07-15.json`. Run:
`uv run python scripts/pre_resistance_base_rate_census.py`.
