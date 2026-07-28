# Handoff → DNA-11: propagate the PGx sentinel-withhold layer to the 7 leaking cells

**From:** Soraya (new-phenotype-cells session) · **To:** DNA-11 (owns the PGx / decoder-science lane)
**Date:** 2026-07-28 · **Type:** verified handoff (R4) — grounded in the real caller code, NOT a guessed runbook.
**Status:** DNA-11's lane to execute; this session did NOT modify any pgx cell (non-colliding).

## Why this is the highest-VOI human-precision lever

Precision = *do not emit wrong calls*. Today, across most PGx cells, a rare **non-core** star-allele is
**silently mis-called `*1`** (a confident wrong answer — worse than an abstention). The fix is the
**sentinel-withhold** pattern: a defining SNP that proves a non-core allele the core proxy cannot resolve
→ **withhold the phenotype** (`phenotype_status=phenotype_withheld`) instead of a confident mis-call. This
raises precision *by construction* on the exact samples that currently fail, and it generalizes across the
whole panel at once.

## The grounded finding — the mechanism EXISTS; only 2 of 9 star-allele cells use it

The sentinel-withhold logic is **already generic** in `dna_decode/pgx/caller.py` (verified 2026-07-28):

- `call_diplotype(..., sentinels=SENTINELS, sentinel_counts=...)` — gene-parameterized (`caller.py:163-174`).
- `caller.py:240-253` — a proven non-core sentinel hit sets `phenotype_status="phenotype_withheld"`,
  `phenotype=None`, and appends the `non_core_allele_sentinel` flag. `_scan_sentinel_counts` (`caller.py:272`)
  drives it from the VCF; an `alt="*"` wildcard = any non-ref ALT at that site signals the non-core allele.

So the machinery is done. The leak is that **only 2 catalogs populate a non-empty `SENTINELS` list**:

| gene | `SENTINELS` | guarded? | note |
|---|---|---|---|
| **CYP2C19** | 2 (rs28399504→*4, rs12769205→*35) | ✅ reference impl | the `*4b`-under-`*17` + `*35`-vs-`*2` traps |
| **CYP2C9**  | 4 (*5/*8/*9/*11) | ✅ reference impl | clean `alt="*"` wildcard-site pattern |
| **CYP2C8**  | `[]` | ❌ **LEAKS** | demotion_rule: "rare non-core mis-called *1 (no sentinel v0)" |
| **CYP2B6**  | `[]` | ❌ **LEAKS** | single-SNP *6-proxy; can't split *6/*9 — sentinel especially valuable |
| **CYP3A5**  | `[]` | ❌ **LEAKS** | non-core non-expressor alleles → silently `*1` (expressor) = wrong direction |
| **TPMT**    | `[]` | ❌ **LEAKS** | demotion_rule names *2/*8/*16 mis-called *1 — clinically high-stakes (thiopurine) |
| **NUDT15**  | `[]` | ❌ **LEAKS** | rarer non-core → *1; clinically high-stakes (thiopurine) |
| **UGT1A1**  | `[]` | ❌ **LEAKS** | non-core (e.g. *27/*37) → *1 (irinotecan) |
| **DPYD**    | `[]` | ⚠️ **weak case** | CPIC only doses the 4 actionable haplotypes → non-actionable-by-design; LOWEST priority |
| CYP2D6 | — (no `SENTINELS`) | separate | SNP-surface; non-core SNP alleles (*14/*15/*21/*40/*46) — a distinct, harder case (structural already handled by the CRAM stack) |
| ABCG2 / SLCO1B1 / CYP4F2 / VKORC1 | — | **N/A** | single-defining-variant cells — no star-allele proxy, so no "non-core proxy mis-call" failure mode |

## The task (per leaking star-allele cell)

For each of the 7 leaking cells, populate its catalog `SENTINELS: list[SentinelVariant]` with the
defining site(s) of the non-core alleles the core proxy cannot resolve — mirroring CYP2C9's
`SentinelVariant(rsid, chrom, pos, ref, alt="*", implies="*N", note=...)` wildcard-site pattern.

**Which alleles need a guard** (grounded in each catalog's own `ALLELE_FUNCTION`/`ACTIVITY_VALUE` +
its demotion_rule; ordered by clinical stakes):

1. **TPMT** — `*2`, `*8`, `*16` (the demotion_rule already names these). No-function → mis-called `*1` (normal) is the dangerous direction (thiopurine over-dose risk).
2. **NUDT15** — the reduced/no-function alleles beyond the modeled `*3` (`*2` already collapses to `*3`; guard the rarer no-function alleles).
3. **UGT1A1** — non-core reduced-function alleles beyond `*80`/`*6` (e.g. `*27`, `*37`).
4. **CYP3A5** — non-expressor alleles beyond `*3/*6/*7` (guard against a non-core non-expressor reading as the `*1` expressor).
5. **CYP2C8** — non-core beyond `*2/*3/*4`.
6. **CYP2B6** — the `*6/*9`-split site (its documented blind spot) + other non-core.
7. **DPYD** — LOWEST priority / possibly skip: CPIC's own posture is non-actionable outside the 4 haplotypes, so a `*1` call there is CPIC-faithful, not a precision defect. Decide explicitly rather than reflexively.

## ⚠️ Anti-fabrication rail (load-bearing)

**Do NOT invent rsIDs / coordinates / ALT alleles.** Each `SentinelVariant`'s `(rsid, chrom, pos, ref, alt)`
must be sourced **verbatim from PharmVar** (the allele-definition authority) for that gene, exactly as the
CYP2C19/CYP2C9 entries were. This handoff names the *alleles to guard* (the WHAT) + the *source* (PharmVar);
it deliberately does **not** supply coordinates, to avoid fabricated catalog data.

## Acceptance criteria (checkable)

1. **No regression:** GeT-RM core-diplotype concordance UNCHANGED for each edited cell (`scripts/pgx_getrm_concordance.py` — the existing falsifier). A sentinel must never fire on a core call.
2. **Real withhold:** a held-out sample carrying a guarded non-core allele → `phenotype_status=phenotype_withheld`, `phenotype=None` (not a confident mis-call). Reuse the CYP2C19/CYP2C9 sentinel test as the template.
3. **Registry truth-up:** update each edited cell's `cell_registry` `demotion_rule` from "mis-called *1 (no sentinel v0)" → the sentinel-guarded wording (so the trust surface reflects the new precision).
4. Frozen AMR + forward surfaces byte-unchanged (pgx is outside them; assert anyway).

## Effort / ordering

~6 star-allele cells × (source PharmVar sites + populate `SENTINELS` + 1 test + registry wording). TPMT /
NUDT15 / UGT1A1 first (highest clinical stakes + clearest non-core sets); CYP3A5 / CYP2C8 / CYP2B6 next;
DPYD is a decide-explicitly (likely skip). The mechanism + 2 reference implementations already exist, so
each cell is a small, well-patterned addition — not new architecture.

## Provenance

Grounded against `dna_decode/pgx/caller.py` (sentinel mechanism) + the 13 per-gene catalogs +
`dna_decode/data/cell_registry.py` demotion_rules, all read 2026-07-28. No pgx code modified by this
session (handoff only). See the session-board `handoffs-out` row (Soraya new-phenotype-cells).
