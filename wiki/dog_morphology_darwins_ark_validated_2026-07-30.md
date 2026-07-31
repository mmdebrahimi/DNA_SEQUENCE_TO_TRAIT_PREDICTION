# Dog morphology relative-signal validation — Darwin's Ark (2026-07-30)

**Verdict: HEIGHT + EAR are strong validated known-SNP morph traits; the 4 covariate-adjusted "rerun" traits
the directive named are NOT single-SNP-mappable on the classic loci.** Step 2 of the dog morphology cell:
the relative-signal validation (does dosage track the owner-reported morphology ordinal), exactly like the
coat concordance. Reproducer: `scripts/dog_morphology_darwins_ark_validate.py` (→ this `.md` + a `.json`).

## 1. Height (Q121) — the 4-locus body-size polygenic score

Reconfirmed with the committed catalog (`dna_decode.pigment.dog_body_size.SIZE_LOCI`): the sum of the 4
big-allele dosages (0–8) vs the Q121 height z-score, **N=3276**:

- **polygenic r = +0.619, R² = 0.383** — ~38% of cross-breed height variance from 4 SNPs. Clean monotonic
  gradient (score 0 → −1.15, score 8 → +0.77). (Pinning + per-locus detail: `dog_body_size_darwins_ark_pinned_2026-07-30.md`.)

## 2. Ear type (Q125) — MSRB3, a NEW validated known-SNP morph trait

There is **no codebook** for the Q-numbers in this dataset, so trait identity was derived FUNCTIONALLY (which
morphology question a known locus's dosage tracks). Scanning the classic single-SNP morphology loci against
all 9 morphology questions surfaced one strong, unambiguous hit:

- **Q125 ↔ MSRB3 ear locus `chr10:8612500:A:G`, r = +0.543, N=2834.** This SNP is the **exact published
  canFam4 ear lead variant** (chr10:8,612,500; Sci Rep 2025). Monotonic dose-response: dose 0 → −1.13,
  dose 1 → −0.32, dose 2 → +0.43.
- **Cleanly resolved from body size:** the ear lead sits 91 kb from the HMGA2 body-size SNP (chr10:8,703,415),
  and Q125's correlation with the *size* SNP is only −0.127 (opposite sign). This is the exact MSRB3-vs-HMGA2
  confound Morrill 2022 had to untangle in the diverse cohort — here it separates cleanly (|r_ear| > 2·|r_size|).

Pinned as `dog_body_size.MORPH_LOCI['EAR']` (single-SNP morphology, separate from the polygenic size score).

## 3. The 4 "rerun" traits (Q124/127/128/245) — honest negative

The four covariate-adjusted rerun morphology traits the directive named do **NOT** map to a classic single
SNP. Best |r| across FGF5 (coat length), KRT71 (curl), and MSRB3 (ear):

| trait | FGF5 | KRT71 | MSRB3/ear | best |
|---|---|---|---|---|
| Q124_rerun | −0.05 | +0.04 | −0.17 | 0.17 |
| Q127_rerun | −0.09 | −0.08 | +0.12 | 0.12 |
| Q128_rerun | +0.04 | −0.04 | +0.05 | 0.05 |
| Q245_recoded | −0.09 | +0.02 | −0.12 | 0.12 |

**Max |r| = 0.213** — no strong single-SNP signal. Two honest reads (can't disambiguate without the codebook):
(a) these traits are **SV/indel-caused** (leg-length FGF4 retrogene, tail T-box, furnishings RSPO2 insertion)
→ substrate-limited exactly like the coat indels; or (b) they are different traits than coat-length/curl/ear.
The signal is **not** a z-score artifact — the raw 34Q values are the same rank-preserving normalized ordinals.

Also: **FGF5 (coat length) + KRT71 (curl) show NO strong signal on any of the 9 morphology questions**
(max |r| ~0.14) — those coat-texture traits are simply not among the measured Darwin's Ark morph Qs.

## Honest conclusion

Step 2 validates **two** strong relative-signal morph traits on this free imputed-SNV substrate — **height**
(polygenic, r=0.619) and **ear type** (single-SNP MSRB3, r=0.543, exact literature lead, resolved from size) —
both the SNP-based regime where the substrate works. It also **corrects the directive's premise**: the 4
rerun traits are not the known-SNP ones (Q125 ear, among the other 5 morph Qs, is). This is a relative-signal
validation (dosage tracks the owner-reported ordinal), NOT a calibrated absolute predictor.

Shipped: `scripts/dog_morphology_darwins_ark_validate.py` (reproducer) + `tests/test_dog_morphology_darwins_ark.py`
(offline) + `dog_body_size.MORPH_LOCI['EAR']` (pinned ear locus). Frozen AMR/forward surfaces byte-unchanged.
Step 3 (ship `typing:dog:morphology` — polygenic height + ear, both validated) is the natural next increment;
the 4 rerun traits ABSTAIN (no known-SNP mapping) rather than a fabricated call.
