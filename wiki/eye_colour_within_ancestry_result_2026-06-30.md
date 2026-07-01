# Eye-colour M2-full — within-ancestry re-score DISENTANGLES the confound (2026-06-30)

M2's structural first cut (`wiki/eye_colour_ancestry_confound_2026-06-30.md`) showed rs12913832 is strongly
ancestry-informative (EUR blue-allele 0.636 vs EAS 0.002, 318×) — raising the worry that the v0 0.993 could
be the SNP *tagging* European ancestry rather than mechanism. M2-full resolves it with the within-ancestry
re-score now that D: is reconnected.

## Method

Stratify OpenSNP users by SELF-REPORTED ancestry (keyword classifier over the `ethnicity` + `Ancestry` free-
text columns → EUROPEAN / NON_EUROPEAN / MIXED_UNKNOWN; conservative — mixed/ambiguous/Ashkenazi excluded).
Re-score v0 (rs12913832 binary calls) within each stratum. **The logic:** within a ~homogeneous-ancestry
group there is little ancestry variance for the SNP to tag, so if accuracy HOLDS there, the signal is
mechanistic.

## Result

| Stratum | N (binary) | accuracy | brown-sens | blue-spec |
|---|---|---|---|---|
| Full cohort | 298 | 0.993 | 0.98 | 1.00 |
| **Self-reported EUROPEAN** | **45** | **1.00** | 1.00 | 1.00 |
| Self-reported NON_EUROPEAN | 9 | 1.00 | 1.00 | — (no blue) |
| MIXED / unknown | 244 | 0.992 | — | — |

**Verdict: `MECHANISTIC_HOLDS_WITHIN_EUROPEAN`.** Among self-reported Europeans, rs12913832 separates
blue/brown 45/45 — accuracy does NOT fall relative to the full cohort. The confound does not drive the
result: rs12913832 predicts eye colour by mechanism (it is the causal HERC2/OCA2 variant), not merely by
tagging ancestry. The non-European stratum (n=9) is all-brown/all-correct — consistent with the blue allele
being near-absent outside Europe (they carry AA/brown and are called brown), though it can't test blue-spec.

## Honest limits

- **Self-reported ancestry via a keyword proxy** (not genetic-ancestry inference from AIMs). A stronger
  version would infer ancestry from genome-wide markers; deferred.
- **n=45 European is modest** (most OpenSNP users left `ethnicity`/`Ancestry` blank → MIXED_UNKNOWN). The
  verdict clears the ≥30 powering floor but is not large.
- **NON_EUROPEAN n=9** is too small to test within-non-European specificity (expected — OpenSNP is
  Euro-dominated; that skew is the whole reason the confound worry existed).
- A verify-in-batch bug was caught + fixed here: "asian" was substring-matching inside "cauc**asian**",
  mis-bucketing ~13 European users into MIXED_UNKNOWN → fixed with word-boundary matching (regression-pinned).

## Disposition

M2 is complete: the confound was (a) quantified structurally (318× EUR/EAS) AND (b) disentangled empirically
(within-European accuracy holds at 1.0). The eye-colour cell's ancestry caveat is now resolved to "confound
present in structure but NOT driving the accuracy — mechanistic within Europeans (self-report proxy, n=45)."

## Reproduce
`uv run python scripts/eye_colour_within_ancestry_validate.py` (reads the D: OpenSNP zip). Classifier:
`dna_decode/data/eye_colour_ancestry.py::classify_self_reported_ancestry`. Tests: `tests/test_eye_colour_ancestry.py`.
