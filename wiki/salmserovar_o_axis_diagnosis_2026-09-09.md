# The O-axis gap is reference-DB content — and my own threshold change is exonerated twice

The SeqSero2 comparison (2026-09-08) localized the serovar caller's entire −0.175 gap to one axis: H1
agreed 200/200, H2 146/146, **O only 159/200**, with 85% of the O errors in two strings. That was
recorded as a **lead, not a diagnosis**. This is the diagnosis.

## The hypothesis I started with was mine, and it was wrong

On 2026-09-04 I lowered this caller's O-antigen **coverage threshold 80 → 40**. The obvious suspicion is
that partial alignments now clear the bar and produce the over-calls. **Tested twice, falsified twice.**

- **`9,46` over-calls:** the winning allele is `O__9,46__239` at **99.70 identity / 100.00 coverage** — a
  full-length, near-perfect hit. No threshold change produced that.
- **`1,3,19` over-calls:** `O__1,3,19__126` hits **121/200 at every coverage cut from 40 through 90**,
  because where it hits its coverage is **99.2–100%**. The cut is irrelevant to it.

My change is not the cause of either. Stating that plainly matters more than it would have if the answer
had gone the other way.

## Two reference-DB defects, both structural

### 1. There is no plain O9 allele

`grep -c ">O__9__"` on the shipped DB returns **0**. The DB carries `9,46` and `9,46,27` but **not `9`**.

O group 9 (D1 — Enteritidis, Typhi, Dublin) is among the most prevalent in *Salmonella*, and the
reference tool calls it on **17 of 187** cohort genomes with a resolved O. Our caller **cannot emit it at
all**, so those genomes are forced onto the nearest available allele — which is exactly the `9,46 → 9`
pattern, 15 of the 33 O disagreements.

**No threshold, selection rule, or code change can fix a missing reference sequence.**

### 2. `O__1,3,19__126` is a 130 bp near-universal fragment

It is **130 bp against an O-allele median of 1155** — the shortest in the DB — and it aligns at 99.2–100%
coverage on **121/200 = 60%** of a diverse cohort, while the reference calls O=1,3,19 on only **2%**.
That is **28× its antigen's true prevalence**. A 130 bp sequence present in most *Salmonella* is not
typing an antigen. Our caller emits `1,3,19` on 14 genomes where the reference says 4.

## The check that was wrong before the allele was

The first version of the audit used a flat "hits >25% of the cohort" bar and flagged **`O__4__231`** too.
That was **the check being wrong, not the allele**: O=4 genuinely *is* 32% of this cohort, so a 32%
hit-rate is correct behaviour. A flat bar cannot separate promiscuity from prevalence.

The audit now measures each allele's hit-rate **against its own antigen's reference prevalence**. `O__4__231`
clears at 1×; `O__1,3,19__126` stands alone at 28×. The false flag was removed before publishing, not after.

## What is deliberately NOT done

**Neither defect is fixed here.** Both fixes need *sourced* reference sequences — a genuine O9 `wzx/wzy`,
a full-length `1,3,19` — and writing biological reference data from memory is the fabrication hazard this
project guards against elsewhere.

**Deleting the bad allele is also not proposed.** `O__1,3,19__126` is the DB's **only** `1,3,19` entry, so
removing it would trade ~10 false positives for the loss of all E4 calling, which the reference says is
correct on 4 cohort genomes. That is a curation decision on shipped data, not a cleanup.

## Sized impact

| defect | bound |
|---|---|
| missing O9 | **17/187 (9.1%)** of genomes structurally uncallable on O |
| `1,3,19` over-fire | emitted on **14**, correct on **4** → ~10 false |
| together | ~27 of the 33 O disagreements, consistent with the 85%-two-strings finding |

Fixing both is the plausible path from O 0.795 toward the reference tool's level. That is an estimate
from disagreement counts, **not a measured post-fix result**.

## Honest limits

- The missing-antigen check is **data-driven, not curated**: it lists only antigens the reference tool
  actually emitted on our own 200-genome cohort. An antigen neither the cohort nor the tool exercised
  would not be flagged, so **absence of a flag is not proof of completeness**.
- The 300 bp / 5× / 10% bars are **hygiene tripwires** separating clear outliers from the bulk, not
  derived constants. They flag candidates for inspection; they never adjudicate an allele.
- Promiscuity is measured against SeqSero2's calls as the prevalence baseline — a **tool**, so this
  inherits that tool's own error rate as a baseline uncertainty.
- One cohort, one DB build. The cohort caps each serovar at 12, so its antigen prevalence is deliberately
  flattened and is not a population distribution.

## Reproduce

```bash
uv run python scripts/salmserovar_db_audit.py
```

Offline from the committed DB plus the cached blastn sweep. Frozen AMR surface byte-unchanged.
See [`salmserovar_db_audit_2026-09-09.json`](salmserovar_db_audit_2026-09-09.json).
