# The over-call risk is no longer untested — and it is real, in Klebsiella

**The independent archive was found, the counter-examples exist, and this time they survive the control.**
BV-BRC yields **67 gentamicin-susceptible `rmt` carriers with measured MICs**, 162 of 169 carriers new
relative to the NCBI-PD sweep. Unlike PD's 60, these are **not** a label artifact.

**The result is organism-stratified, and that split is the whole finding:**

| organism | R | S | PPV(`rmt` → R) |
|---|---|---|---|
| **Klebsiella** | 58 | **64** | **0.475** |
| Salmonella | 15 | 3 | 0.833 |
| **E. coli / Shigella** — *the deployed rule's validated scope* | 12 | **0** | **1.000** |
| Pseudomonas / Acinetobacter / Enterobacter / Citrobacter | 11 | 0 | 1.000 |

**In Klebsiella the deployed rescue would over-call on more than half of carriers. In E. coli — the only
organism it was validated on — not one susceptible carrier has been found, here or anywhere.**

---

## Why BV-BRC, after PD was exhausted

The [archive search](rmt_independent_archive_search_2026-09-03.md) eliminated the obvious candidates on
structure, not effort: `rmt` surveillance studies **ascertain on high-level aminoglycoside resistance**
(so they hold zero susceptible carriers by construction), NARMS **feeds NCBI-PD weekly** (not
independent), and ATLAS has ~634k MICs but **no paired sequence**.

BV-BRC is different on **both** axes, which is what makes it a real second opinion:

- **Phenotype** — BioSample antibiograms *plus* ~300 hand-curated publications, each row carrying `pmid`,
  `laboratory_typing_method` and `testing_standard`. The publication-curated part is content PD's
  `AST_phenotypes` field does not hold. Filtered to `evidence = Laboratory Method` throughout: BV-BRC also
  ships **ML-predicted** phenotypes, and scoring a deterministic rule against a model's output is gate G1
  firing.
- **Genotype** — `sp_gene` is **CARD/BLAT**, a different caller with different thresholds from both our
  AMRFinder and PD's. The carrier call is genuinely re-determined, not re-served.

**Overlap measured, not assumed:** 7 of 48 resolvable carriers are shared with the PD sweep → **162 new**.
An "independent archive" returning the same isolates would be independent in name only.

## Two controls, both of which the counter-examples survive

**1. The carrier call is not the artifact.** MIC ≤1 for an `rmtB` carrier is biologically extraordinary —
G1405 methylation shifts aminoglycoside MICs 256–512×. The obvious explanation is a bad gene call. It
isn't:

| | CARD identity | CARD coverage | partial hits |
|---|---|---|---|
| **susceptible carriers** (n=67) | median 100, min **99** | median 100, min **100** | **0** |
| resistant carriers (n=129) | median 100, min 93 | median 100, min 82 | 17 |

The susceptible carriers have **better** gene calls than the resistant ones. These are full-length,
near-identical `rmtB` genes in isolates measured at gentamicin MIC ≤1.

**2. The labels are not the artifact.** 94% of the susceptible carriers come from one study
(pmid 36801013) — the same concentration signature that condemned PD's set. So the same `aac(3)` yardstick
decides it, with the decision rule written into the code before the numbers were read:

| within pmid 36801013 | R | S | % R |
|---|---|---|---|
| carries `aac(3)`, no `rmt` | 154 | 2 | **99%** |
| **carries `rmt`** | **18** | **63** | **22%** |
| no known gentamicin determinant | 4 | 110 | **4%** |
| *the same `aac(3)` stratum outside this study* | *147* | *31* | *83%* |

**Verdict: `SPECIFIC_TO_RMT`.** This study's gentamicin column tracks genotype better than the archive
average — 99% R for the undisputed determinant, 4% R when no determinant is present. It is not calling
everything susceptible; it contributes 18 **resistant** `rmt` carriers of its own. Only `rmt` is anomalous.

Contrast PD's PRJNA1322038, which failed the identical test at 2% vs 97% and 0R/60S. **Same control, same
threshold, opposite verdict** — which is exactly what a control is for.

## What this changes

- **The standing honest limit is retired.** "The over-call risk is UNTESTED" has been true since the v2
  lock. It is now **measured**: in Klebsiella, PPV 0.475.
- **The deployed rule's own scope survives.** E. coli is 12/12 here and 146/146 on PD. The v2 validation
  (N=131 E. coli, sens 0.523 → 0.892 at unchanged specificity) is not contradicted by any of this.
- **The rule must not be extended to Klebsiella**, and CLAUDE.md's note that non-E. coli carriers are
  "outside the validated organism scope" turns out to have been load-bearing rather than cautious
  boilerplate.
- **`gene context: core | acquired` is the shape of the explanation.** This is exactly the pattern the
  [AMRrules schema](prior_art_decoder_landscape_2026-09-03.md) encodes as a first-class field and that our
  own recorded lesson names — a determinant's phenotypic consequence is organism-dependent. Whether
  Klebsiella `rmtB` is silent, poorly expressed, or plasmid-context-dependent is not settled here.

## What is NOT settled

- **The mechanism.** A full-length `rmtB` at MIC ≤1 is unexplained. Silencing, promoter loss, low copy
  number and expression context are all candidates; none was tested. This memo reports a phenotype
  discordance, not its cause.
- **Whether the E. coli scope is genuinely safe or merely under-sampled.** Twelve carriers is not many. The
  honest statement is "no susceptible E. coli carrier has been found in two independent archives", not
  "none exists".
- **CARD vs AMRFinder.** Carriers here are CARD calls. A carrier under CARD is not guaranteed to be one
  under AMRFinder, which is what the deployed rule actually consumes. Re-calling these 67 genomes with
  AMRFinder is the obvious next check and was not done.
- **One dominant study.** 94% concentration survived the control, but it remains one lab's isolates; a
  second Klebsiella source would strengthen it considerably.

## Reproduce

```bash
uv run python scripts/gentamicin_rmt_bvbrc_hunt.py      # -> wiki/gentamicin_rmt_bvbrc_hunt.json
uv run python scripts/gentamicin_rmt_bvbrc_control.py   # -> wiki/gentamicin_rmt_bvbrc_control.json
```

Network-only; no frozen file is touched. BV-BRC 403s the default `python-urllib` User-Agent and answers
an outage with **HTTP 200 wrapping a 503 envelope** — both are handled explicitly, because either one
would otherwise look like an honest empty result.
