# CYP2D6 hybrid IDENTITY (*13/*36/*68) — Phase B technical plan (drafted 2026-07-07)

**Status: EXECUTED 2026-07-06 (via `--advance`) — Phase A GO + Phase B classifier SHIPPED + full-N validated.**
Phase A (`scripts/cyp2d6_psv_evidence.py`) proved the signal; Phase B (`scripts/cyp2d6_psv_phaseb_falsifier.py`
+ `dna_decode/pgx/cyp2d6_hybrid_identity.py`) validated at FULL N (spec 1.0, *68 4/4, *36 6/8; *13 n=1
unpowered) and shipped the three-level abstaining classifier. `wiki/cyp2d6_hybrid_identity_2026-07-06.md`.
Remaining residuals (subtle *36 conversion, *13 power) are long-read / more-truth, NOT code. This plan folded
in the pre-exec `/brainstorm` (2 rounds) refinements.

## Where we are
- SHIPPED CYP2D6 surface: SNP diplotype (46/47) + copy-number (*5/*xN, 26/26) + hybrid-PRESENCE (CYP2D7
  depth, sens 0.62/spec 1.0) + **hybrid-identity Phase A evidence table** (read-level PSV D6-fraction profile;
  the falsifier proved the signal reproduces the Cyrius signatures).
- The remaining gap: assign WHICH hybrid (*13 vs *36 vs *68) — Phase B.

## Corrected facts (from the brainstorm — do not regress)
1. **PSV source = the real Cyrius `data/cyp2d6_psv/CYP2D6_SNP_38.txt`** (117 PSVs, GRCh38-native,
   paired-coord `chr pos_CYP2D6 base_CYP2D6 pos_CYP2D7 base_CYP2D7 annotation`). NOT `CYP2D6.json`; NO liftover.
2. **Count BOTH pos_CYP2D6 AND pos_CYP2D7 per PSV** — D6-coordinate-only measures aligner placement, not
   biology. (Phase A already does this.)
3. **The validation bar is per-allele with honest denominators, NOT ">=90% on 3 anchors"** — the committed
   hybrid set is 13 rows (*68=4, *36=8, **\*13=1 [NA19785 only]**); anchors NA24631/HG01161 are NOT in the
   data. *13 is single-sample-UNPOWERED and must be labelled so.

## Phase B — the abstaining identity classifier (the build)
- **B1. Full-set evidence extraction.** Run the Phase A extractor (mpileup D6+D7 regions) over the full 13
  committed hybrid rows + a matched non-hybrid set (normals + pure-dups + deletions) + the fetchable Cyrius
  anchors (NA24631=*36, HG01161=*13 — fetch to grow *13/*36 power). Emit the per-sample regional profile.
- **B2. Copy-normalized evidence (brainstorm R2).** Per-site fractions alone over-call on multi-copy
  arrangements (*36x2, mixed). Normalize the D6-fraction by the sample's total CYP2D6+CYP2D7 copy number
  (already available from the shipped CN caller) so a duplication's flat-but-elevated profile is not mistaken
  for a shift.
- **B3. Three-level classifier (brainstorm R2).** Output exactly one of: **resolved identity** (*13/*36/*68)
  | **hybrid_present_identity_unresolved** | **evidence_not_callable**. NEVER force a hybrid-positive into a
  specific allele. Decision features (from Phase A, validated): directional 5'-3' shift (+ -> *68-like; - ->
  *13-like) + exon9-tip dip (*36-like) + the breakpoint region (contiguous run of shifted PSVs).
- **B4. Callability gate (brainstorm R2 open-Q).** Require a minimum callable-PSV count + regional coverage
  (derive the threshold on B1 data) before allowing a resolved call; else `evidence_not_callable`.
- **B5. Per-allele honest validation.** Report resolved-identity concordance PER ALLELE with honest
  denominators; **flag *13 as single-sample-unpowered** unless the anchors add power. Target: beat "presence-
  only" on the well-powered alleles (*68, *36); do not claim a *13 number on n=1-2.

## Open tradeoffs (ratify at build time)
- **mpileup flag contract.** Phase A used permissive `-B -q0 -Q0` (matches the depth path) and the falsifier
  passed. Decision: keep permissive, or add MAPQ/baseQ filtering for cleaner paralog fractions — derive
  empirically on the B1 set (the flags that best separate), don't assume.
- **Aligner-portability.** PSV fractions may be aligner/reference-build dependent -> scope the cell to
  GRCh38-aligned CRAMs with the benchmark mapper, or invest in cross-aligner robustness (larger).

## Documented residuals (do NOT try to close in v0.4)
- The **\*36 embedded-exon-9 gene-conversion** blind spot: a pure conversion yields no CN change and only a
  tiny exon-9-tip signal -> may stay `hybrid_present_identity_unresolved`. Honest; matches Cyrius's own caveat.
- **\*36x2 / multi-copy** explicit copy-number modelling -> defer until read-level linkage evidence is strong.

## Falsifier for Phase B (pre-committed)
Resolved-identity concordance on the well-powered alleles (*68 n>=4, *36 n>=8) must beat presence-only AND the
three-level classifier must never emit a WRONG specific allele on a non-hybrid (spec-first, like the shipped
presence detector). If B1's copy-normalized profiles do NOT cleanly separate *68 from *36 from non-hybrids at
the full N -> STOP; presence-only + Phase A evidence table stay the honest ceiling and identity is documented
as long-read-required.

## Reuse (near-zero new infra beyond Phase A)
The Phase A extractor + the shipped CN caller + the remote-CRAM tooling + the registry/report-card pattern are
all reused. New code ~= one classifier fn (features -> three-level call) + a validation script + tests.
