# The third gate screen found a defect in the gates — Oxford E. coli bacteraemia

**Verdict: `REJECTED`, gate G7.** Run: `uv run python scripts/screen_candidate_gates.py --candidate
oxford` (offline, seconds). Artifact: `wiki/oxford_gate_screen_2026-09-11.json`.

This closes candidate row 1 of `project_state/evidence-surface-2026-08-31.md` — *"screen a THIRD
candidate through the gate screen (n=2 worked examples is the schema's stated limit)"*. The schema's own
docstring named n=2 as its limit. The third example found a real conflation in G2, and produced a
failure mode neither prior example had.

## Why Oxford, specifically

The two existing screens fill one diagonal each: PEAR is **L4 → CLEARS**, HBV is **L1 → REJECTED on
G1**. Nothing had ever exercised the L1 gates *past* G1 — HBV dies at the first gate because no measured
HBV phenotype exists at all. So the informative third candidate is an **L1 candidate with an
unimpeachable label**, to find out whether the L1 path can reach a verdict on the other nine gates.

Oxford is that, and it is also the sharper choice for a second reason: **this project already trusted
it.** It was scored for external re-validation on 2026-06-15 (gentamicin acc 0.990) and re-read for the
`rmt` prevalence probe on 2026-09-03. Pointing our own gates at our own accepted cohort is the same
move as `source_diverse_validate` refusing six of our own SCORED cells.

Every number below is **measured from the committed deposit** (`data/raw/oxford/`), not transcribed from
a memo.

## The defect: G2 tested a different quantity than the memo defines

`wiki/negative_results_map_2026-06-13.md` defines G2 as:

> **Study == class** — the label is confounded with the source study / submitter (**one BioProject
> supplies most of one class**) — *contingency table of class × BioProject/submitter; a dominant cell =
> trip*.

The implementation tested `largest_source_share` **of the whole cohort** against a 60% bar. That is a
different quantity, and the two diverge hardest exactly where it matters most:

**Oxford is a single-source cohort.** One study, supplying *both* classes. Cohort share = 1.00, so the
implemented rule tripped. But with one source there is **no source variation at all**, so source cannot
explain any label variance — the confound G2 names is **structurally impossible**. The screen rejected a
cohort whose wet-lab MICs we had already validated against, for a confound it cannot have.

Worse, the rule made G2 **unsatisfiable by any single-source cohort whatsoever**, regardless of label
quality.

**Fixed** (`dna_decode/eval/rejection_gates.py`): G2 now returns `not_applicable` when `n_sources <= 1`,
prefers the memo's actual contingency (`max_class_share_from_one_source`) when supplied, and falls back
to the cohort-wide share only as an explicitly-labelled weaker proxy.

**Single-source is not harmless, and the fix does not pretend otherwise.** It is a *generalizability*
limit, not a *validity* one — and it already has its own machinery: the `source_concentration` disclosure
layer and `source_diverse_validate.py`'s 60% refusal. Oxford is itself the proof the concern is real:
**zero `rmt` carriers among 4,979 isolates**, so it structurally cannot test a rule keyed on them. The
concern is routed to the layer that owns it; the caveat rides in G2's reason string.

**Both committed screens still reproduce** (`--verify`: PEAR `CLEARS`, HBV `REJECTED`). PEAR is
unaffected (constructed variation → G2 already `not_applicable`); HBV dies at G1 regardless.

## The result: excellent labels, unusable metadata

| gate | verdict | measured |
|---|---|---|
| G1 circular label | pass | Oxford hospital-lab broth-microdilution MIC; no genomic tool produced it |
| G2 study == class | **n/a** | single source — confound structurally impossible (was a false trip) |
| G3 sampling-defined | pass | an MIC reading; ascertained on **bacteraemia**, independent of resistance |
| G4 surveillance domination | pass | smallest class n=**192** (gentamicin 192 R / 2,681 S), bar 20 |
| G5 assembly attrition | pass | 2,897 records with MIC; 4,979 carry an AMRFinder scan |
| G6 phenotype censoring | pass | **0** of 2,873 gentamicin MICs fail to resolve R vs S |
| G7 provenance not separable | **TRIP** | submitter/centre/collection populated on **0%** of records |
| G8 dedup collapses balance | insufficient | no MLST or Mash distance ships in the deposit |
| G9 / G10 decoder scoreability | pass | the deployed rules already score these drugs on this cohort |

**This is a third distinct failure mode.** HBV fails because *no measured phenotype exists*. Oxford's
labels are as good as this project has: wet-lab MIC, ascertained independently of the phenotype, both
classes well-powered, essentially zero breakpoint censoring. It fails because the **deposit carries only
`guuid`** — no submitter, centre or collection field on any record — so a leakage-clean
provenance-disjoint split cannot be built from it as shipped.

## Why G7 is reported and G2 was fixed — the distinction that matters

Both gates tripped on Oxford. Only one was a bug, and conflating them would have turned a real finding
into special pleading:

- **G2 was a schema defect.** It tripped on a confound that is *logically impossible* given the cohort's
  design. No measurement could change that; the rule was testing the wrong quantity.
- **G7 is a true finding.** It tripped on a capability that is *genuinely absent* — the fields really
  are not there. The honest response is to report it, not to declare it inapplicable.

Having just found one over-firing gate, the temptation is to read every inconvenient gate the same way.
That is the motivated-reasoning failure this repo has recorded before
(`[[feedback_apply_your_own_standard_to_your_own_cohort]]`). The test applied here: *is the concern
structurally impossible, or merely unmeasured?*

## A unit trap that would have silently corrupted the screen

The MIC columns are `_lower` / `_upper` pairs holding **log2 dilution indices**, not mg/L — values are
negative (`-3.0 → 0.0`). Applying CLSI breakpoints to them directly returned **`R = 0` for all three
drugs**, which reads as a clean, powered, all-susceptible cohort rather than as an error.

The authoritative convention is recorded in `scripts/oxford_score.py:14` — **MIC = 2\*\*upper** (the
paper encodes MIC=8 as `lower=log2(4), upper=log2(8)`). Under it, what the naive reading called
"straddling the breakpoint" are precisely the resistant isolates: gentamicin 192 R, ciprofloxacin 370 R,
ceftriaxone 238 R. **The interpretation was read out of the existing consumer rather than guessed** —
inventing a plausible one would have produced a fully self-consistent, entirely wrong screen.

## Honest limits

- **This screens the DEPOSIT, not the cohort in principle.** The submitter/centre fields may well be
  fetchable from ENA for `PRJNA604975`. That has **not** been done, so G7's trip is a statement about
  the files on disk. If those fields were fetched, G7 could flip and the verdict would turn on G8.
- **G8 is unmeasured.** No MLST or Mash distance ships in the deposit, so effective lineage count is
  unknown. G7 is decisive on its own, so nothing here hinges on it — but the screen is not complete.
- **A REJECTED verdict bounds the LABEL question only.** It does not say Oxford is unusable: the project
  used it correctly as a *scoped external check against an already-frozen decoder*, where disjointness
  is established at the accession/BioSample level by `cohort_manifest`, not by an internal provenance
  split. What the gates say is that it cannot serve as a **general** label source.
- **n=3 worked examples.** Better than the n=2 the schema itself flagged, still small. G2's fix is
  verified against all three; the fallback path is now exercised by none of them.

## Reusable

**A gate derived from prose can silently test a different quantity than the prose defines, and the
divergence shows up only at an extreme.** G2's cohort-share proxy agrees with the memo's class × source
contingency across most multi-source cohorts and inverts completely at one source. Screening candidates
that resemble the ones the schema was built from will never surface it — the third example was chosen
precisely because it sat where neither prior one did.

**When a screen rejects something you already trust, that is the highest-information outcome available**
— either the screen has a bug or the trust was misplaced. Here it was one of each: a real bug in G2, and
a real limitation of the deposit in G7.
