# A shipped bug the sibling cell had already fixed — found by testing a second typing cell

Extending the wet-lab validation from Salmonella to E. coli O:H serotype found a **live defect in the
shipped caller**, diagnosed it from a fix the sibling cell already carried, and closed it:

| axis | before (coverage-only) | after (identity-primary) | change |
|---|---|---|---|
| **H accuracy** | 0.770 — 34 misses | **0.926** — 11 misses | **+0.155**, misses −68% |
| **O accuracy** | 0.931 — 9 misses | **0.962** — 5 misses | +0.031 |
| resolution (both axes) | unchanged | unchanged | the fix picks *which* allele, not *whether* |

Same 150 isolates, same labels, same DB. Only the selection rule changed.

---

## The hypothesis that led here, and how it survived contact

The Salmonella serovar cell scored 0.702 against a wet-lab label, and its failure was **not uniform** —
the H antigen was frequently right where O was unresolved or mis-grouped. That suggested something
narrower than "the caller is weak": **the O-antigen axis might be the systematic weakness of the
serotyping cells.** E. coli O:H is the natural second test — independent cell, different allele DB, and
a label that decomposes into exactly the two axes under test.

**The hypothesis was half right, and the half that was wrong is the more useful half.**

- **O reproduced as a *resolution* weakness** — 87.3% resolved vs H's 98.7%, the same shape as Salmonella.
- **On *accuracy* the asymmetry INVERTED.** H resolved almost always but was only 77% right, while O was
  93% right when it committed.

So the two axes were failing in **different ways**: O **abstains**, H is **confidently wrong**. That
distinction matters more than the original hypothesis, because an over-confident wrong call is far more
dangerous downstream than an abstention — and it is what pointed at the real defect.

## The defect

H misses were **concentrated, not diffuse**: `H21→H8` alone was 9 of 34 (26%), and the top four pairs
were 65%. That is the signature of systematic allele confusion.

The salmserovar report card records a bug **already fixed there**: *"flagellin alleles cross-hybridize at
full coverage, so `_best_per_axis` selected by coverage-only picked the WRONG H antigen."* Inspection of
`dna_decode/serotype/runner.py` confirmed the E. coli caller was **still selecting by coverage only**, in
two places — and E. coli H antigens are fliC flagellins, the same biology.

**The fix was applied to salmserovar and never propagated.** It also predicts the observed asymmetry: O
antigens use wzx/wzy, which cross-hybridize far less, so O errors stayed scattered (max 3 for one pair).

Diagnosis confirmed **by code inspection before the fix**, not inferred from the numbers alone.

## Same pattern elsewhere — surfaced, deliberately NOT fixed

| caller | selection | status |
|---|---|---|
| `ktype` | identity-primary | already correct |
| `pneumoserotype` | `(coverage, identity)` — coverage-**primary** | **pattern present, not fixed** |
| `plasmid` | `cov > cur[...]` — coverage-**only** | **pattern present, not fixed** |
| `mlst` | exact allele match | not applicable |

Neither carries a rationale comment, so the choice looks unexamined rather than deliberate. **They were
not changed**, for two reasons: no wet-lab validation cohort exists for either, and the biology differs —
pneumococcal cps is whole-locus reference matching and plasmid Inc typing is replicon-family matching,
neither of which is per-antigen flagellin. Identity-primary is **not automatically correct** there.
Changing them blind would be an unvalidated behaviour change of exactly the kind this work exists to
catch. **Building a Quellung-labelled pneumo cohort is the named next test.**

## A false verdict this run produced, and the guard added

The first pilot printed a confident `NO_AXIS_ASYMMETRY_HYPOTHESIS_NOT_SUPPORTED` from **0 resolved calls
on both axes** — because the scorer read `O_type`/`H_type` keys the caller does not emit (it returns one
combined `"O?:H7"` string). A broken extraction wearing the costume of a clean negative.

The scorer now **refuses to emit a verdict** (exit 3) when both axes resolve zero, since that is a
plumbing signature rather than a result.

## Honest limits

- **No in-silico incumbent.** NCBI-PD does not populate `computed_types` for E. coli, so unlike the
  Salmonella run there is no comparator and the **absolute** accuracy is much less interpretable. What
  this run supports is the internally-controlled **before/after** and **O-vs-H** contrasts — same
  isolates, same labels, same caller.
- **The before/after is not blinded.** The fix was chosen *after* seeing the failure pattern. It is a
  strong result (+0.155 on a pre-existing cohort, mechanism confirmed in code, matching a documented
  sibling fix) but it is not a pre-registered test.
- **Label provenance is unprovable per isolate.** E. coli O:H serotyping is traditionally agglutination,
  but here there is no incumbent score to use as a circularity probe the way the Salmonella run had.
- O-only labels score the O axis alone, so **the two axes have different denominators** — reported per
  axis, never pooled.
- `H-` (non-motile) is a **meaningful** label, not a missing value, and is scored as one.
- Per-O-group and per-BioProject caps flatten prevalence away from O157: 55 O-groups, 33 BioProjects,
  largest source share 0.167 — clears the project's own 0.60 diversity bar.

## Reproduce

```bash
uv run python scripts/serotype_oh_validate.py --target 150
```

Needs blastn + network. The frozen AMR surface is byte-unchanged — this is a typing cell.
