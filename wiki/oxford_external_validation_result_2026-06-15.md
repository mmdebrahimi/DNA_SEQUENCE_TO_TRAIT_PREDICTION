# Oxford external validation — the frozen v0.5.0 AMR decoder holds on independent measured MIC (2026-06-15)

First fully-INDEPENDENT external validation of the shipped deterministic decoder
(`call_resistance`). Substrate: the Oxford E. coli cohort (Lipworth et al., Lancet Microbe
2025; repo `github.com/samlipworth/ecoli_mic_arg`) — **2897 isolates** with linked
**measured MIC** (broth microdilution) + the cohort's **own AMRFinder genotype**, both keyed
by study `guuid`. Scored by `scripts/oxford_score.py` (pure join → frozen `call_resistance` →
`independent_cohort_validate._conf`; no genome download, no Docker).

## Result (binary R/S, the right scorer for clean measured MIC)

| drug | n scored | accuracy | sensitivity (R) | specificity (S) | R / S |
|---|---|---|---|---|---|
| **ciprofloxacin** | 2841 | **0.960** | **0.935** | **0.963** | 370 / 2471 |
| **gentamicin** | 2873 | **0.990** | **0.922** | **0.995** | 192 / 2681 |
| ceftriaxone | 2868 | 0.729 | 0.945 | 0.709 | 238 / 2630 |

Strict-tier (HIGH_R/HIGH_S, 4× margin) where computable: cipro acc 0.97 / sens 0.935 / spec 0.991 (n=1005).

## Verdict
- **cipro + gent GENERALIZE STRONGLY** to an independent, large, UK clinical measured-MIC cohort — the strongest trust signal the decoder has received. cipro spec actually *improves* vs the frozen NCBI-PD cell (0.991 vs 0.70); gent near-perfect both sides.
- **cef spec on the Oxford genotype layer is 0.709 — but this is an AMRFinder-VERSION CONFOUND, NOT a deployed-decoder gap (CORRECTED 2026-06-15 after /brainstorm).** On the Oxford run, of 765 cef false-positives, 649 (85%) are counted blaEC-5 (intrinsic chromosomal AmpC, near-universal in E. coli, ceftriaxone-S). The catch: **the Oxford run scored the cohort's OWN OLDER AMRFinder output, where blaEC is `Subclass=CEPHALOSPORIN`** and thus passes the cef extended-spectrum refinement. **The DEPLOYED decoder runs AMRFinder v4.2.7 (`ncbi/amr:4.2.7`), where `fam.tsv` classes blaEC as `BETA-LACTAM/BETA-LACTAM`** → the cef refinement already EXCLUDES it (verified: zero counted blaEC in cached v4.2.7 runs). So the deployed decoder does NOT have the blaEC over-call; Oxford cef 0.709 reflects the paper's genotype-caller version, not the shipped rule.
- **Deployed cef is already strong; no frozen-rule change warranted.** On the v4.2.7 NCBI-PD cef cells the deployed rule scores **E. coli 0.967/0.967/0.967** and **Klebsiella 0.944 acc (sens 1.0 / spec 0.9)**. The one residual deployed-path candidate — blaOXA-1 (a narrow class-D penicillinase mis-subclassed `CEPHALOSPORIN` on v4.2.7) — was tested by what-if on the cached v4.2.7 cef cohorts: excluding it changes **nothing** (both cells identical). → **cef thread CLOSED: leave `amr_rules.py` untouched.** Do NOT exclude blaEC in the frozen rule (it would fix an old-AMRFinder adapter artifact, not the shipped decoder, and retrofit the reproducibility-frozen rule for zero deployed gain).
- **Honest external cef number (deferred, optional):** to state how the *deployed* decoder performs on Oxford cef, the Oxford reads (deposited reads-only) would need assembling + a v4.2.7 AMRFinder re-score (blaEC would drop out → expected ~0.95). Not required to justify leaving the rule alone; cipro+gent already stand as the independent-validation headline.

## Independence + leakage (why this counts)
- **Independence:** Oxford's own measured MIC AND its own AMRFinder genotype — MORE independent than the NCBI-PD provenance-disjoint cells (independent genotype caller too).
- **Caveat (honest):** the paper's AMRFinder version/DB differs from the frozen Docker pipeline; the adapter header-normalizes `Gene symbol → Element symbol` so the frozen parser reads it. What is validated is the decoder's **`call_resistance` RULE**, on fully independent genotype+phenotype.
- **Leakage = DISJOINT (project-level, by construction):** PRJNA604975 + PRJNA1007570 deposit **0 NCBI assemblies** (reads-only); the decoder's 1039 tuning/validation accessions are **100% GCA_/GCF_** → an Oxford isolate (no GCA) cannot be in the tuning set.

## Process note (root-cause caught, not shipped as a finding)
The first run showed cipro **sens=0.0** — a plumbing artifact (the frozen parser reads `Element symbol`; the paper's older AMRFinder names it `Gene symbol`), NOT a decoder failure. Caught by the "suspect plumbing before the finding" discipline (370/370 join coverage + a cipro-R isolate with 4 QRDR mutations returning S). Fixed adapter-side; frozen `amr_rules.py` byte-unchanged. The cef-spec gap, by contrast, survived the same scrutiny (cipro+gent clean, cef sens correct) → it is a real finding.

Artifact: `wiki/external_validation_oxford_2026-06-15.json`. Scorer: `scripts/oxford_score.py`. Data (gitignored): `data/raw/oxford/`.
