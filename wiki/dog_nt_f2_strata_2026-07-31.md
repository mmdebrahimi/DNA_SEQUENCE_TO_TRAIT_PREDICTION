# Dog masked reconstruction F2 — structured-signal vs calibration (2026-07-31)

The decisive test that closes the reconstruction track. The scale finding (bigger NT reconstructs
better; `wiki/dog_nt_scale_2026-07-31.md`) raised one question the brainstorm flagged as decisive: is
the NT-vs-Markov win **STRUCTURED** (concentrated in coding/functional regions = a real signal) or
**FLAT** (mere calibration)? Answered on Kaggle T4 (free), NT-v2-500M, per-base marginal NLL, tokens
bucketed by canFam4 CDS annotation (NCBI efetch; UCSC unreachable from this host).

## Result

Substrate: 34 windows across canFam4 chr1 20–21 Mb (636 coding + 6130 intergenic 6-mer tokens);
disjoint 200 kb Markov-train region (chr1 22.0–22.2 Mb); Markov given its best k (=3, hardest baseline).

| stratum | n bases | NT-500M accuracy | NT NLL | Markov NLL (k=3) | NLL delta (Markov − NT) |
|---|---|---|---|---|---|
| **coding** | 3,816 | 0.291 | 1.5865 | 1.3588 | **−0.228 (NT LOSES)** |
| **intergenic** | 36,780 | 0.408 | 1.3708 | 1.3184 | **−0.052 (~tie/slight loss)** |

**structured-signal gap (coding delta − intergenic delta) = −0.175.**

## Interpretation — NOT a structured signal (calibration, confirmed)

1. **NT is WORSE, not better, in coding [controlled — same windows, same Markov].** The "captures
   functional structure" hypothesis predicted NT would be MORE advantaged in coding; the opposite holds.
   NT-500M's reconstruction skill is concentrated in **low-information, repetitive INTERGENIC** sequence
   (accuracy 0.408) and is near-chance in high-information CODING sequence (0.291, vs the 0.25 null).
   A low-order Markov captures coding codon-composition better than NT does. This is the brainstorm's
   "global reconstruction win = calibration, not biology" prediction, confirmed and sharpened: the skill
   is anti-correlated with functional constraint.

2. **The scale-test win is WINDOW-dependent, NOT baseline-fragile [RESOLVED 2026-07-31 by a controlled
   test].** My initial hypothesis (the scale win was inflated by a weak 4 kb Markov) was FALSIFIED. Holding
   the scale window fixed and varying ONLY the Markov training set: NT-500M NLL 1.0668; Markov-4kb NLL
   **1.3561** (delta +0.289 — reproduces the scale test exactly); Markov-**200kb** NLL **1.3384** (delta
   +0.272). The pure baseline-strength effect is only **−0.018 nats** — a 200 kb Markov barely beats a
   4 kb one on this window, and NT still wins by +0.27. So the scale-test's arbitrary first-1200 bp window
   was simply a FAVORABLE (easy, low-complexity) window; the broader F2 annotation-selected windows are
   genuinely harder (NT's own NLL 1.07 → 1.37). The scale/F2 disagreement is a **window/sampling effect**,
   not baseline strength. What is SOLID: the within-run scale ORDERING (bigger NT > smaller NT) stands;
   what does NOT generalize is "NT-500M beats Markov" — it beats it on favorable windows and loses on a
   representative sample. Controlled result: `wiki/dog_nt_scale_window_control_2026-07-31.json`.

   **Reusable lesson (corrected):** a single-window reconstruction number does NOT generalize — sample
   windows representatively before claiming a win. (Baseline strength was HYPOTHESIZED, tested, and was
   NOT the driver here; the general "check baseline strength" hygiene still holds, but window-sampling was
   the real confound.)

## Verdict — the reconstruction track is CLOSED (honest negative for decoding)

The world model's masked-reconstruction skill (a) does not robustly beat a fairly-trained Markov
baseline, and (b) is concentrated in repetitive low-information regions, NOT functional/coding ones →
**it is not a useful decoding signal.** This is fully consistent with — and does not reopen — the closed
embedding-vs-phenotype negatives: native-objective reconstruction skill ≠ phenotype-decoding skill, and
here even the reconstruction skill is calibration rather than functional structure.

**Caveat:** some coding tokens sit near exon boundaries where flanking context is intergenic (less
coding context → harder for NT); this is a minor fraction and cannot explain a −0.175 gap, but a
CDS-interior-only re-run would tighten it (deferred).

## Scope + reproducibility
NT-v2-500M only (the scale-test winner); one 1 Mb region; Kaggle T4 (free, no money). Notebook
`scripts/kaggle_dog_nt_f2_strata.py` (embeds the strata windows + train region); strata built by
`scripts/build_f2_strata.py` (NCBI efetch gbwithparts → CDS intervals → coding/intergenic windows).
Result `wiki/dog_nt_f2_strata_2026-07-31.json`. Frozen AMR/forward surfaces byte-unchanged.
