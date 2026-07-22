# R1 × R2 pharmacogene synthesis — corroboration is PROPERTY-SPECIFIC (2026-07-22)

**Status:** ✅ an honest cross-regime finding (not a win, not a wall). For CYP2C9 — the one gene with BOTH a
mature R1 catalog cell (star-allele → CPIC phenotype, GeT-RM 1.0) AND today's R2 molecular substrate (MaveDB
DMS) — the two regimes **do NOT corroborate for the clinically-important alleles**, and the *reason why* is the
value: **R1×R2 corroboration is property-specific.** Frozen surface byte-unchanged.

## The question
The project decodes CYP2C9 in two independent regimes:
- **R1 (deterministic catalog):** star allele → CPIC metabolizer phenotype. `dna_decode/pgx/cyp2c9_catalog.py`,
  validated **73/73 vs GeT-RM** (the CDC free consensus panel) — the R1 cell is real + independent.
- **R2 (learned molecular):** ESM2 / measured DMS of CYP2C9 variants (today's MaveDB work; ESM2 pgx 0.547).

Do they agree? For each star-allele-defining **missense** variant with a known CPIC function, does its measured
DMS score put it in the damaging tail?

## Result (real MaveDB fetch; `scripts/pgx_r1r2_synthesis.py`)
| star | variant | CPIC function (R1) | abundance-DMS percentile (R2) | in damaging tail? |
|---|---|---|---|---|
| **\*2** | R144C | decreased (AS 0.5) | 0.769 | **no** (stable) |
| **\*3** | I359L | **no function (AS 0.0)** | 0.525 | **no** (stable) |
| \*8 | R150H | decreased | 0.482 | no |
| \*9 | H251R | decreased | 0.856 | no |
| **\*11** | R335W | decreased/no | 0.150 | **YES** |
| \*5 | D360E | no function | (not in this assay) | — |

**Only 1/5 CPIC-reduced-function alleles land in the R2 damaging tail.** The clinically-dominant \*2 and \*3 —
the alleles CPIC's warfarin guidance turns on — read as *stable*.

## Why (the load-bearing mechanism — data-grounded, not a guess)
**Every free CYP2C9 MaveDB assay is an ABUNDANCE assay** (the score-set titles literally say "CYP2C9 abundance
scores as measured by VAMP-seq"; verified across `00000095-a-1`, `00000095-b-1`, `00000062-a-1`). But **CYP2C9
\*2 (R144C) and \*3 (I359L) are CATALYTIC-ACTIVITY defects — the protein is expressed and stable, it just
metabolizes substrate more slowly** (established CYP2C9 pharmacology; \*3 reduces intrinsic clearance with
retained protein). An abundance readout **cannot see an activity defect**. The one concordant allele, \*11
(R335W), is a *destabilizing* variant — the property the assay actually measures.

## The lesson (reusable)
- **R1×R2 corroboration is PROPERTY-SPECIFIC.** The R2 molecular assay validates the R1 clinical call ONLY when
  it measures the SAME molecular property the phenotype depends on (activity vs abundance vs binding). The free
  CYP2C9 DMS measures the *wrong property* for the metabolizer phenotype, so it cannot validate \*2/\*3.
- This is the **g2p regime boundary at finer grain**: earlier the boundary was "learned wins only when
  fitness-aligned"; here it sharpens to "aligned to the SAME molecular property." It echoes the HIV-DRM finding
  (resistance reached via conservative substitutions at conserved sites → abundance/likelihood scorers call
  them benign) — the same failure shape, a different gene family.
- **Consequence for the decoder:** the R1 catalog (GeT-RM-validated) remains the trustworthy CYP2C9 metabolizer
  cell; the free R2 substrate does NOT independently validate it. To cross-validate CYP2C9 clinical function
  with a molecular assay you'd need a CYP2C9 **activity** DMS (4-MU / luciferin turnover), which is not in the
  free MaveDB set. Named, not built.

## Honest scope
- CYP2C9 only. **CYP2C19 is excluded** — its core alleles are non-missense (\*2 splice / \*3 stop / \*17
  promoter), invisible to any missense DMS or ESM.
- Abundance DMS only (the free CYP2C9 substrate); an activity DMS would be the property-matched test.
- The percentile/damaging-tail orientation is set from the nonsense-vs-all median (lower = damaging, verified).
- Kept as a negative — NOT tuned to manufacture concordance.

Reproduce: `uv run python scripts/pgx_r1r2_synthesis.py` → `wiki/pgx_r1r2_synthesis_2026-07-22.json`.
Related: `wiki/pgx_report_card.md` (R1), `wiki/mavedb_pgx_esm2_2026-07-21.md` (R2), memory
`feedback_g2p_decoder_regime_boundary`.
