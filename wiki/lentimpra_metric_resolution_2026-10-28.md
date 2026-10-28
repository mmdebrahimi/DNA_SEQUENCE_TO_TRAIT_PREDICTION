# lentiMPRA published-metric resolution: ATTEMPTED and REFUSED

Date: 2026-10-28 · Run: `2026-10-27-2342-lentimpra-metric-ceiling-unblock`

## CORRECTIONS to this memo, made 2026-10-28

Two claims in the version committed as `b28fa8dc` were wrong. Both are corrected below; the
superseded wording is recorded here verbatim so the change is legible rather than a silent
overwrite.

1. **The decisive next check could not have worked.** It read: *"read this workbook's own
   header/units row for the 'Benchmarking results' sheet -- blocked today because the substrate
   volume (D:) is detached"*. The parser proves the header is `test fold | <method names>` and
   carries no coefficient, so reading it could never settle the metric. The rows that can settle it
   are the ones the parser **discards**, and they are now reachable via `unparsed_rows()`.
2. **The arm's blocker was overstated.** It read: *"the arm
   `substrate-change-constructed-eukaryotic-or-human-panel-UNRUN` stays blocked on provenance, not
   on compute"*. That conflates cannot-COMPARE with cannot-RUN. Nothing gates running the arm; the
   metric gates comparing its result against the published incumbents. The arms ledger never made
   this error -- its own rationale reads *"the acquisition wall is DISCHARGED; the arm is still
   PLANNED because it has not RUN"*, with the metric clause scoped to comparison: *"no number may
   be compared against it until that is resolved"*.

A third claim is **softened, not retracted**: see the corroboration section's wording on how weakly
a two-decimal agreement discriminates between coefficients.

## Verdict

`UNRESOLVED_TABLE_IDENTITY_UNVERIFIED`. `published_incumbents()` still returns
`metric="UNSPECIFIED_IN_SOURCE"`. **No number may be compared against the published incumbents.**

**Table identity is ONE of three findings, not the load-bearing one.** The refusal does not depend
on it: `paper_reports_both` -- the publication reporting Pearson *and* Spearman for model
predictions, so there is no single global benchmarking metric to copy -- is independent and stands
whatever the table numbering turns out to be. Were table identity settled tomorrow, the metric would
still be unspecified for this sheet.

This is the outcome the plan's Step 7 was written to allow. The failure mode it guards against is a
GUESSED metric, not an unresolved one: a plausible-sounding choice here yields a fully
self-consistent comparison that is entirely wrong, with nothing downstream able to catch it -- the
same shape as the refused ECOFF values.

## A correction to my own VOI estimate, made the same day

Earlier on 2026-10-27 I ranked this **first** on cost-adjusted value and described it as "one read of
the paper's methods", "costs minutes". **That was wrong**, and the error was in the premise, not the
estimate: I assumed the metric was stated somewhere and merely unread. Three findings say otherwise.

## What was found, from the source

Read the open-access publication text directly (36 pages extracted; 30 occurrences of "Pearson",
6 of "Spearman") rather than a search-tool paraphrase.

1. **The paper reports BOTH coefficients for model predictions.** Its MPRALegNet panels carry
   *"corresponding Pearson (r) and Spearman (rho) correlation values"*. So there is no single global
   benchmarking metric in the paper to copy onto our sheet.
2. **The table identities do not evidently match, and this is the load-bearing finding.** The
   paper's IN-TEXT *Supplementary Table 6* is a biochemical-feature / bigWig manifest -- it is cited
   for *"1,506 bigWig files for HepG2 cells; 1,206 files for K562 cells; and 277 files for
   WTC11/H1-ESCs"* and for *"a compilation of thousands of biochemical features"*. Our `S6_FILE` is
   `Table S6 - folds and performance.xlsx`, carrying a `Benchmarking results` sheet **the paper never
   names**. So "the paper's Table S6 states the metric" is a question about a table the paper does
   not describe, and citing the paper's Table S6 as the provenance of our incumbent numbers is
   UNSUPPORTED. `wiki/lentimpra_acquisition_2026-10-27.md` refers to "Table S6" in a way that invites
   exactly that conflation.
3. **The workbook itself is unreachable.** `D:/dna_decode_cache/mpra_human` and the `E:` equivalent
   are both absent (measured), so the rows that could state the metric for THIS sheet could not be
   inspected. **CORRECTED 2026-10-28:** this finding originally named the sheet's own "header/units
   row" as *the one place that could state the metric*. That is false twice over -- the header is
   `test fold | <method names>` and carries no coefficient, and it is not the only candidate row. The
   candidates are the rows the parser DISCARDS (a caption, a units note), now reachable via
   `unparsed_rows()`. What the detached volume blocks is reading those.

## Sourced corroboration, recorded as NOT the answer

The paper states MPRALegNet reached *"Pearson r ~ 0.83"* combining folds; our HepG2 MPRALegNet
fold-mean is **0.8333**. **A two-decimal agreement discriminates poorly between the two
coefficients**, which are often close to one another on standardized, near-linear data -- so this is
much weaker evidence for Pearson than it first reads as. (Stated qualitatively on purpose: putting
an unsourced numeric interval inside a memo whose whole argument is that unsourced values must be
refused would undercut the memo.) It is **not** the
answer: it compares a COMBINED-fold figure to a PER-CELL fold-mean, and says nothing about the
sheet. It is recorded in code so a future session does not re-derive it, and labelled so a future
session does not promote it.

## Honest limits

- Scoped to the open-access version fetched today; the publisher's final supplementary manifest was
  not consulted and may number its tables differently from the paper's in-text citations. **The
  numbering discrepancy may be a renumbering between versions rather than a different table** -- that
  is unresolved either way, which is the point.
- A failed read is evidence about a ROUTE, not a property of the source (the six-routes lesson).
  **Three routes have now been tried, and only one of the three results is evidence about the
  source at all:**

  | route | result | evidence about the source? |
  |---|---|---|
  | search index | never retrieved the supplementary listing | no -- uninformative |
  | `nature.com/articles/s41586-024-08430-9` | `303` redirect to `idp.nature.com` | **yes** -- the article HTML is behind authentication |
  | Europe PMC REST (DOI, then title) | `503 "Search service is temporarily unavailable. Please retry later."` | no -- service down, explicitly retryable |

  **The Europe PMC result was nearly mis-filed as evidence.** Queried by DOI and then by title it
  returned what parsed as `hits: 0`, which reads exactly like "no such record". A control query for
  a paper that certainly is indexed returned `0` as well, and the raw body turned out to be the
  service's own 503 error object being parsed as an empty result set. A 503 wearing the costume of
  an absent record. **Retry this route before concluding anything from it** -- and the lesson
  generalises: when a lookup returns nothing, run one control lookup that must succeed before
  treating the nothing as a finding.

## The decisive next check -- now runnable, not a prose instruction

Call **`unparsed_rows()`** and read the rows the parser DISCARDS. The header row cannot settle it:
the header is `test fold | <method names>` and carries no coefficient. A caption or units row is
where a metric statement would sit, and those are exactly the rows `published_incumbents()` skips
without a record -- which is why no amount of reading its output could ever have found one.

`unparsed_rows(data_dir, sheet="Benchmarking results")` returns every skipped row with its 1-based
index and non-empty cells, and raises `LentiMpraDataError` when the workbook is absent rather than
returning `[]` (an empty list would read as "there are no unparsed rows" when the truth is "could
not look" -- the same conflation this memo is about). Attaching the substrate volume is a
PREREQUISITE for running it, not the reason the previous check was wrong.

**What this blocks, stated at its true scope.** The metric blocks **comparison against the published
incumbents**. It does **not** block running the arm: the acquisition wall is discharged, and
`substrate-change-constructed-eukaryotic-or-human-panel-UNRUN` is `planned` because it has not RUN,
on compute and substrate-attachment like its siblings.
