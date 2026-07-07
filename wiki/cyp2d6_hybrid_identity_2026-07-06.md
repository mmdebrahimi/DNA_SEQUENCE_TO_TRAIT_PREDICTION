# CYP2D6 hybrid IDENTITY (Phase B) — read-level PSV classifier, full-N validated (2026-07-06)

**The flagship pharmacogene's last gap closed (for the well-powered alleles).** Phase A proved the read-level
CYP2D6-CYP2D7 PSV D6-fraction signal reproduces the Cyrius hybrid signatures (n=1/type). Phase B validated it
at FULL N on the structural panel and shipped a three-level abstaining identity classifier.

## Full-N falsifier (`wiki/cyp2d6_psv_phaseb_falsifier.json`)
39 samples (13 hybrid / 26 non-hybrid), 117 Cyrius PSVs, mpileup base-counts at BOTH pos_CYP2D6 + pos_CYP2D7.

| metric | value |
|---|---|
| **verdict** | **GO_BUILD_CLASSIFIER** |
| **non-hybrid specificity** (the confound test) | **1.00 (26/26)** — no normal/pure-dup/deletion mis-called a hybrid |
| **\*68** identity signal | **4/4 (1.0)** |
| **\*36** identity signal | **6/8 (0.75)** — the 2 misses are subtle exon-9 gene-conversions (NA18617 *36/*36, NA18572) |
| **\*13** identity signal | 1/1 — **UNPOWERED (n=1, NA19785)** |
| overall hybrid sensitivity | 0.85 (11/13) |

## The classifier (`dna_decode.pgx.cyp2d6_hybrid_identity.classify_hybrid_identity`)
Three levels (never forces a call): **resolved** `*68` (5' CYP2D6 / 3' CYP2D7, intron-1) / `*13` (opposite,
5' CYP2D7) / `*36` (exon-9 conversion) · **`hybrid_present_identity_unresolved`** (weak/ambiguous — e.g. a
subtle *36) · **`evidence_not_callable`** (< 50 callable PSVs). Features: the directional 5'-3' D6-fraction
shift + the exon-9-tip dip. HIGH-SPECIFICITY (a resolved call is trustworthy). The D6-fraction is inherently
copy-normalized (a fraction, not absolute depth), so a duplication's flat-elevated profile does not trip it.

## Honest scope
- **Resolves *68 cleanly + *36 partially; abstains on subtle *36 conversions** (documented residual — matches
  Cyrius's own gene-conversion caveat) and on ambiguous profiles.
- **\*13 is single-sample-validated (UNPOWERED)** — HG01161's CRAM fetch failed, so *13 stays n=1; the call
  fires but carries no statistical weight.
- Read-depth/pileup, GRCh38-aligned CRAMs; NOT long-read; NOT a clinical tool. Reference tool: Cyrius (WGS).

## CYP2D6 is now decoder-complete for a VCF+CRAM
SNP diplotype 46/47 · copy-number *5/*xN 26/26 · hybrid-presence sens 0.62/spec 1.0 · **hybrid-identity
*68 4/4 / *36 6/8 / spec 1.0** — everything short-read WGS can resolve; the remaining residuals (subtle *36
conversion, *13 power) are long-read / more-truth-data, not code.
