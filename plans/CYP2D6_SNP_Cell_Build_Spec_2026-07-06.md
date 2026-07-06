# CYP2D6 SNP-defined star-allele cell — build spec (the last major pharmacogene, 2026-07-06)

**Status: PLAN — awaiting user ratification before execution (Planning-STOP).** Decompose + `--plan` produced;
no build code written yet.

## R2 pre-bar framing check — the load-bearing decision (derive, don't assert)

The user asked for a "CYP2D6 structural caller." **That bar is INFEASIBLE from a VCF and is the wrong bar.**
Measured from the real data this session:
- CYP2D6 GeT-RM truth (local `star-allele-comparison_common.tsv`, `CYP2D6_getrm_cons`): **87 samples.**
- **56/87 are SNP-only-truth (VCF-decodable).** **31/87 are STRUCTURAL** — gene deletions (`*5`), duplications
  (`*36x2`, `xN`), and CYP2D6-CYP2D7 hybrids (`*13`, `*36`, `*68`). Structural detection needs **read-depth /
  breakpoint analysis on a BAM/CRAM** (the Cyrius/Aldy/StellarPGx approach) — it is NOT recoverable from the
  1000G/PGP-UK phased VCFs (no coverage track).
- The 56 SNP-only samples are **well-powered** (comparable to CYP2C8 82 / TPMT 85 / CYP3A5 8) and dominated
  by clinically-central alleles: **\*4 (19), \*2 (19), \*17 (6), \*41 (5), \*29 (5), \*35 (4), \*3 (3)**, +
  \*6, \*9, \*10, \*21, \*40, \*14, \*15, \*45, \*46.

**Reframed bar (feasible + validatable + honest):** a **CYP2D6 SNP-defined star-allele cell**, validated on the
**56-sample SNP-only GeT-RM subset**, with the **31 structural samples as a documented, WITHHELD (not
mis-called) BAM-required blind spot** — exactly the "core-comparable + non-core residual" pattern the CYP2C
cells already use, plus a structural-withhold sentinel. This is the reachable form of "the last major
pharmacogene," not a Cyrius reimplementation.

## Decompose — families + flow-down + critical path

Terminal: **CYP2D6 is a registered, GeT-RM-validated PGx cell (SNP surface) with an honest structural blind
spot**, on the same trust surface as the other 8 PGx genes.

- **F1 — Catalog** (`dna_decode/pgx/cyp2d6_catalog.py`): verified GRCh38 (chr22, ~42.12–42.13 Mb) defining
  variants for the SNP-defined core {*2,*3,*4,*6,*9,*10,*17,*29,*35,*41}; CPIC activity-value → phenotype
  (CYP2D6 uses the activity-score model like CYP2C9). Multi-SNP alleles (e.g. *2 = 2851C>T + 4181G>C; *17,
  *41) → reuse `compound_caller`. **DEP: none.**
- **F2 — Structural withhold sentinel**: a sentinel layer that WITHHOLDS the phenotype when the VCF evidence
  is consistent with a structural allele the SNP proxy can't resolve (mirrors the CYP2C19 *4/*35 sentinel).
  v0-honest fallback: score core-comparable only on SNP-only truth; structural truth = non-core residual.
  **DEP: F1.**
- **F3 — Wiring**: `runner.call_cyp2d6` + `PGX_GENES` + `dna-pgx --gene cyp2d6` dispatch + registry
  `CellContract` + `cli_routable_manifest`. **DEP: F1.**
- **F4 — Validation**: fetch chr22 CYP2D6 region (existing `scripts/fetch_1000g_region.py`) + a `cyp2d6`
  config in `scripts/pgx_getrm_concordance.py` → real core-comparable concordance on the 56-sample subset.
  **DEP: F1, F3, data.**
- **F5 — Report card + docs**: add CYP2D6 to `pgx_report_card` + regenerate the certification capstone (→ 74
  cells) + README note. **DEP: F4.**

**Critical path:** F1 → F3 → F4 → F5 (F2 parallel to F3; F2's full sentinel can be a v0.1 follow-on).

## Acceptance bar (DRAFT — ratify before execution)

MVP `--until-mvp` criteria (all checkable):
1. `dna_decode/pgx/cyp2d6_catalog.py` — every defining coord VERIFIED via Ensembl REST + empirically
   AF-confirmed on 1000G (the *4/*10 orientation-flip guard); unit tests green (`test-exit-0 pytest`).
2. CYP2D6 scored vs `CYP2D6_getrm_cons` on the SNP-only subset → a committed
   `wiki/pgx_getrm_concordance_cyp2d6_*.{md,json}` (`file-exists`).
3. Full pgx + registry suite green incl. a CYP2D6 coverage-in-`PGX_GENES` test (`test-exit-0 pytest`).
4. Frozen bacterial/viral/fungal/TB AMR surface byte-unchanged (`test-exit-0` leak guard).

## Verdict-time pre-commitments (before results land — per project pattern)

Core-comparable concordance on the 56-sample SNP-only subset:
- **≥0.90 → CLEAN PASS** (ships as NEAR_INDEPENDENT, like CYP2C8/3A5/TPMT).
- **0.75–0.90 → NOISY PASS** — inspect the mis-calls; likely a missing multi-SNP allele definition (e.g. *2
  needs both SNPs) → add the component, re-score. Ships with the residual documented.
- **<0.75 → FAIL / re-scope** — probable structural-truth leakage into the "SNP-only" filter (a `+`/`x`/`*5`
  the regex missed) OR a strand/orientation error. Diagnose before shipping; do NOT ship a low number.

## Engine reuse (minimal new code — mirror CYP2C8/TPMT)

Catalog dataclasses + `call_diplotype`/`assemble_compound_diplotype` + `pgx_getrm_concordance.py` GENES-config
+ `fetch_1000g_region.py` + the registry `_PGX_CONTRACTS` pattern are ALL reused. New code ≈ one catalog +
one runner fn + one concordance config + tests. No new machinery.

## Honesty rails (carry from the PGx cells)

- **Structural blind spot is load-bearing**: *5/*13/*36/*68/xN are NOT VCF-decodable → WITHHELD, never
  mis-called; the concordance is explicitly "SNP-surface only, 56/87 GeT-RM samples; structural = BAM-required
  (Cyrius-class), out of scope." Never imply full CYP2D6 typing.
- CALLING validatable vs GeT-RM consensus; PHENOTYPE faithful-to-CPIC (activity-score). NOT a clinical tool.
- Frozen AMR surface untouched (non-frozen pgx package).

## Provenance
GeT-RM subset counts (56 SNP-only / 31 structural) + the SNP-only star distribution measured 2026-07-06 from
`tests/data/pgx_getrm/star-allele-comparison_common.tsv`. Coords to be Ensembl-verified at F1 (rs3892097 *4,
rs1065852 *10, rs35742686 *3, rs5030655 *6, rs28371706 *17, rs28371725 *41, rs16947 *2, …). Est. ~18–25
executable steps (a CYP2C8-class build).

## Alternative considered (the "different organism/trait" fork)
A new organism/trait (e.g. a 2nd Gram-positive AMR drug, or a new human trait) was the other user-offered
option. CYP2D6-SNP is recommended over it: highest-value remaining pharmacogene, GeT-RM-validatable NOW on a
well-powered subset, ~100% machinery reuse (near-zero new infra risk). A new organism/trait carries dataset +
substrate unknowns (the label wall) that CYP2D6-SNP does not.
