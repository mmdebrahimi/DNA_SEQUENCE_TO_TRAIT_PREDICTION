# lentiMPRA published-metric resolution: ATTEMPTED and REFUSED

Date: 2026-10-28 · Run: `2026-10-27-2342-lentimpra-metric-ceiling-unblock`

## Verdict

`UNRESOLVED_TABLE_IDENTITY_UNVERIFIED`. `published_incumbents()` still returns
`metric="UNSPECIFIED_IN_SOURCE"`. **No number may be compared against the published incumbents.**

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
   are both absent (measured), so the sheet's own header/units row -- the one place that could state
   the metric for THIS sheet -- could not be inspected.

## Sourced corroboration, recorded as NOT the answer

The paper states MPRALegNet reached *"Pearson r ~ 0.83"* combining folds; our HepG2 MPRALegNet
fold-mean is **0.8333**. Agreement at two decimals is suggestive of Pearson. It is **not** the
answer: it compares a COMBINED-fold figure to a PER-CELL fold-mean, and says nothing about the
sheet. It is recorded in code so a future session does not re-derive it, and labelled so a future
session does not promote it.

## Honest limits

- Scoped to the open-access version fetched today; the publisher's final supplementary manifest was
  not consulted and may number its tables differently from the paper's in-text citations. **The
  numbering discrepancy may be a renumbering between versions rather than a different table** -- that
  is unresolved either way, which is the point.
- A failed read is evidence about a ROUTE, not a property of the source (the six-routes lesson).
  Two routes were tried: the search index, then the full text. Neither is exhaustive.

## The decisive next check, and it is cheap once unblocked

Attach the substrate volume and read the `Benchmarking results` sheet's own header row. If it names
the coefficient, this resolves in one command. **Until then the arm
`substrate-change-constructed-eukaryotic-or-human-panel-UNRUN` stays blocked on provenance, not on
compute** -- which is a narrower and more actionable blocker than it had this morning.
