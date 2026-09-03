# The `rmt` counter-examples exist, all 60 come from one submission, and they cannot test the rule

**The deployed gentamicin v2 rule's one untested risk, hunted properly for the first time.** The rule
widens the frozen `subclass_any={GENTAMICIN}` with `symbol_rescue=^(rmt[A-H]\d*|npmA\d*)$`. Its
sensitivity gain is measured (0.523 → 0.892). Its **specificity** claim never was: no S-labelled `rmt`
carrier had ever been found, so "specificity unchanged" was an *absence of counter-examples*, not a bound.

**Result: 60 S-labelled carriers found — the first ever — and every one fails gate G2. The specificity
risk remains UNTESTED, now for a measured reason rather than an absence. Separately, the evidence FOR the
rule is much stronger than it was: 146/146, up from 62/63.**

---

## Why the previous search could not have found them

`gentamicin_rmt_label_hunt.py` keeps only PD rows whose accession is **already in the local AMRFinder
cache**. It therefore answers *"of my 109 carriers, which are labelled?"* — never *"of all labelled
isolates, which carry `rmt`?"*. A counter-example outside that cache is structurally invisible to it, and
outside the cache is exactly where one would live. That is a limit of the question, not of effort.

**The inversion:** PD metadata ships `AMR_genotypes` — NCBI's own AMRFinder calls — in the same row as the
measured `AST_phenotypes`. So the question can be asked the other way round over every labelled isolate,
with no local AMRFinder run and no cache bound. Population: **20,816 gentamicin-labelled isolates across
11 organism groups**, vs 63 before.

The rescue regex is **imported from the deployed rule**, never retyped — a hunt scored against a re-typed
pattern would not be testing what ships.

## What came back

| | R | S | PPV(rmt→R) |
|---|---|---|---|
| All labelled `rmt` carriers | 146 | 60 | 0.709 |
| **Inside BioProject PRJNA1322038** | **0** | **60** | **0.000** |
| **Outside it (23 other BioProjects)** | **146** | **0** | **1.000** |

Taken at face value the first row retracts the rule. The second and third rows are why it does not: a
determinant cannot be 100% R across 23 independent submissions and 0% R in one. **Perfect separation by
submission is gate G2 (study == class) at its maximum**, and the largest-source share among the S carriers
is **100%** (University of Queensland).

## The control — because "one project" is not by itself a verdict

Two hypotheses survive concentration alone, with opposite consequences: the project's labels are odd
(the rule stays untested), or the project genuinely holds susceptible `rmt` carriers (the rule is
defective). The discriminator was **stated in the code before the numbers were seen**: use `aac(3)` — the
classic gentamicin-modifying enzyme, the determinant the *frozen* rule already counted — as an internal
yardstick that nobody disputes.

| within PRJNA1322038 | R | S | % R |
|---|---|---|---|
| carries `aac(3)`, no `rmt` | 3 | 130 | **2%** |
| carries `rmt` | 0 | 60 | 0% |
| **no known gentamicin gene** | 219 | 37 | **86%** |
| *the same `aac(3)` comparison outside* | *1,452* | *50* | ***97%*** |

**That project calls its `aac(3)` carriers susceptible 98% of the time, against 97% resistant everywhere
else.** And its isolates carrying *no* gentamicin determinant are called resistant 86% of the time. The
gentamicin column is **anti-correlated with genotype** — consistent with an inverted or mis-mapped R/S
encoding in that one submission, and not comparable to the rest of PD.

**Verdict: `LABEL_ARTIFACT`.** These 60 isolates cannot serve as the specificity test.

## What this does and does not settle

- **It does NOT vindicate the rule.** A label artifact means these counter-examples cannot test it, not
  that none exists. The over-call risk stays **UNTESTED** — the honest limit in `CLAUDE.md` stands,
  sharpened rather than retired.
- **It does strengthen the evidence for the rule considerably.** Previously 62/63 on a cache-bounded set
  of 63. Now **146/146 = 1.000** on a cache-independent population spanning 23 BioProjects and 11 organism
  groups. Every `rmt` carrier with a label from an ordinary submission is resistant.
- **The deployed validation is uncontaminated.** Zero overlap between PRJNA1322038 and the 24 accessions
  in `wiki/gentamicin_v2_validation_2026-08-31.json` — checked, not assumed.
- **Nothing was changed in the frozen surface**, and the v2 lock is untouched. This is measurement.

## Update 2026-09-03 — SOLO grading strengthens this, and the pooled number was the wrong denominator

A prior-art scan (`wiki/prior_art_decoder_landscape_2026-09-03.md`) surfaced the WHO TB catalogue's
**SOLO** method: grade a determinant only on isolates where it appears WITHOUT another known determinant
for the same drug. By that standard the pooled 146/206 above credits `rmt` for **47 resistant carriers
that also carry `aac(3)`** — the classic gentamicin-modifying enzyme — which cannot tell you what `rmt`
did.

| stratum | R | S | PPV | 95% CI (Wilson) |
|---|---|---|---|---|
| pooled (as first reported) | 146 | 60 | 0.709 | [0.643, 0.767] |
| **SOLO — no `aac(3)` co-carriage** | 99 | 21 | 0.825 | [0.747, 0.883] |
| co-carriage with `aac(3)` | 47 | 39 | 0.547 | [0.442, 0.647] |
| **SOLO, artifact project excluded** | **99** | **0** | **1.000** | **[0.963, 1.000]** |

The co-carriage stratum at 0.547 is where the pooled number's dilution came from. The properly-controlled
evidence for the rule is **99 solo carriers, all resistant**, clearing WHO's grade-1 bar (≥5 solo
occurrences, PPV CI lower ≥ 0.25) by a wide margin — stronger and better controlled than what was first
published here. It does **not** touch the specificity question: zero susceptible solo carriers outside the
artifact project remains an absence of counter-examples, not a bound. Reproduce with
`uv run python scripts/solo_ppv.py`.

*Solo is the conservative estimator — it discards data to buy an uncontaminated denominator; a Nature
Communications 2025 analysis found penalised multivariable regression beats SOLO on sensitivity for TB.
And "solo" here means no co-carried `aac(3)` only, because that is the single co-carriage flag the hunt
recorded; a wider co-determinant screen would lower the solo count further.*

## Honest limits

- The **carrier** call is NCBI's AMRFinder (PD's `AMR_genotypes`), a different version/DB from ours —
  independent of our pipeline, but still a tool-derived feature. Only the phenotype is measured.
- Same archive as the earlier hunt (NCBI-PD), different **population**. It is a new search, not a new
  data source; a genuinely independent archive (e.g. a clinical MIC collection) would be stronger.
- The control tests one project's **internal consistency**. It cannot prove those labels wrong — only
  that they do not behave like any other submission. The inverted-encoding reading is an inference.
- `I` (intermediate) never appeared among carriers (0 of 206), so nothing here bears on borderline calls.
- Klebsiella and *P. aeruginosa* carriers are outside the rule's E. coli-validated organism scope; they
  are reported for completeness, and the 146/146 figure spans all groups.

## Reproduce

```bash
uv run python scripts/gentamicin_rmt_specificity_hunt.py     # -> wiki/gentamicin_rmt_specificity_hunt.json
uv run python scripts/gentamicin_rmt_project_control.py      # -> wiki/gentamicin_rmt_project_control.json
```

Both are network-only and touch no frozen file. A run with fetch errors or `--row-cap` records
`complete=false` and cannot be quoted as a bound.
