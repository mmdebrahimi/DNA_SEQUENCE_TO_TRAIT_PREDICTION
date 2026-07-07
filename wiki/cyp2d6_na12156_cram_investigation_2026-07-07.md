# Real-human CRAM acquisition + the CYP2D6 46/47 miss (NA12156) — 2026-07-07

**Directive:** acquire real-human CRAMs to attack the one CYP2D6 GeT-RM SNP miss (46/47), which the
concordance flagged as a *structural confound* (`cnv_hybrid_unassessed`). Real-surface, end-to-end.

## What was acquired (the reusable capability)

`scripts/resolve_1000g_cram.py` — resolves any 1000G sample → its 30x high-coverage `*.final.cram` URL by
parsing the canonical `1000G_2504_high_coverage` sequence indices (2504 panel + the 698-related index),
`ftp://`→`http://` normalized for the remote-slice path. This turns the prior ad-hoc curl into a committed,
tested (6 offline tests), live-verified capability that feeds `cyp2d6_structural_probe.py --compute`'s
`sample/truth/cram_url` cohort TSVs. Live-verified: `NA12156 → ERR3239310`, `NA12878 → ERR3239334`.

## The NA12156 finding (structural surface, real CRAM)

Sliced NA12156's CRAM remotely (Docker samtools depth, ENA reference auto-fetch, no full-reference download)
over the CYP2D6 body + single-copy control:

| sample | GeT-RM truth | SNP-surface call | CYP2D6/control depth | ratio | structural CN |
|---|---|---|---|---|---|
| NA12156 | `*1/*4` | `*4/*4` (the miss) | 45.23 / 37.89 | **1.194** | **CN = 2 (NORMAL)** |

**This contradicted the CNV hypothesis.** NA12156 is copy-number-**normal** (ratio 1.194 ≈ the NORMAL
baseline 1.26; no `*5` deletion, no `*xN` duplication). So the SNP `*4/*4`-vs-truth-`*1/*4` miss is **NOT a
copy-number event.**

## Narrowed diagnosis

With CN ruled out, the remaining explanations for a `*1/*4` individual reading homozygous-`*4` at rs3892097
(at CN-normal) are:
1. **A CYP2D6–CYP2D7 hybrid on the `*1` background** (e.g. `*1/(*68+*4)` — `*68` is a common tandem partner
   of `*4`). The CN surface is blind to hybrids (they can be CN-normal); the **read-level PSV hybrid-identity
   surface** (`scripts/cyp2d6_psv_evidence.py` + `dna_decode.pgx.cyp2d6_hybrid_identity`) is the tool that
   resolves this. **← the concrete next step.**
2. A 1000G-panel genotyping artifact at rs3892097 (called `1/1` where truly `0/1`).

## Concrete next step (named, not done here)

Run the read-level PSV hybrid surface on NA12156's now-acquirable CRAM: generate the 117-PSV D6+D7 mpileup
(the generation loop is currently ad-hoc — committing it as a helper is part of this next step) → feed
`cyp2d6_psv_evidence` → check for a `*68`-like directional 5′-high/3′-low signal. If `*68+*4` is present, the
structural+hybrid surfaces TOGETHER explain the miss (and the honest resolution is `*1/(*68+*4)`, not a SNP
error); if flat, the miss is a panel artifact. Either outcome closes the 46/47 open question honestly.

## Honesty

- Real CRAM, real depth, real boundary (R3). CN=2 is a measured fact, not an assertion.
- This did NOT (yet) move 46/47 → 47/47 — it **removed the wrong hypothesis** (CNV) and pointed at the right
  tool (hybrid surface). The committed 40-row structural ratios set is unchanged (NA12156 measured separately).
- Frozen AMR surface byte-unchanged.
