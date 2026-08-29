# The rescue is now MEASURED: +0.369 sensitivity on 131 leakage-gated isolates that had never been scored

Two runs ago this was "well-motivated, in-distribution-safe, rescue unproven — the prospective cohort's
per-isolate determinant calls are not on disk". They are not. But a different set of isolates was, and
nobody had noticed.

## Where the data came from — free, and already paid for

AMRFinder is the expensive step (~95 s/genome, Docker). This repo had **1,818 genomes with cached
AMRFinder output**. Against the fail-closed accession manifest, **311 of them are genuinely disjoint** from
every cohort that ever tuned a rule. NCBI Pathogen Detection carries an `AST_phenotypes` call for **200** of
those. Three (organism, drug) cells clear 20R/20S. **Zero marginal compute.**

## The result, on E. coli × gentamicin (n=131, 65R/66S, all disjoint)

| | accuracy | sensitivity | specificity |
|---|---:|---:|---:|
| frozen rule | 0.756 | **0.523** | 0.985 |
| candidate (`rmt*`/`npmA`) | **0.939** | **0.892** | **0.985 — unchanged** |

**+0.369 sensitivity at zero specificity cost.** 24 of 31 false negatives rescued.

And the mechanism reproduces independently: **24 of the 31 false negatives (77%) carry an `rmt`-family
gene**, against the prospective accrual's "24 of 28". Two disjoint isolate sets, same signature.

The frozen rule's 0.523 here also independently corroborates the prospective accrual's 0.429 — both far
below the frozen cohort's 0.893. Three datasets now say the same thing.

## The other two cells (frozen rule, previously unscored, disjoint)

> **CORRECTION (2026-08-29).** The Klebsiella × ceftriaxone row below should not have been reported. The
> `source_diverse_validate.py` arm applies a source-diversity bar to its OWN cohorts, and that cohort
> fails it: **5 BioProjects but 68% from one**, over the 60% bar. Reporting a number from a cohort as
> concentrated as the ones this work criticises is exactly the error being exposed. It is now emitted as
> `status: source_concentrated` with **no metrics**. The E. coli cipro row stands (8 sources, 31% share).

| cell | n | acc | sens | spec |
|---|---:|---:|---:|---:|
| E. coli × ciprofloxacin | 131 | 0.977 | 0.958 | 0.988 |
| ~~Klebsiella × ceftriaxone~~ | ~~66~~ | — | — | — (retracted: source-concentrated) |

Cipro holds up strongly on data it has never seen. These are **not** added to the report card — they are
scored here, not censused through the provdisjoint arm, and conflating the two is the shared-key trap.

## What is STILL not established, and it is the same thing

**Zero S-labelled `rmt` carriers.** Not in the 150 local labelled isolates, not among 63 publicly-labelled
carriers, and not among these 66 disjoint S isolates. So specificity "unchanged at 0.985" remains
**arithmetic, not evidence** — the candidate cannot produce a false positive where no S isolate carries the
gene.

Three independent datasets failing to produce a single S-labelled `rmt` carrier is consistent with the
mechanism (16S methyltransferases confer *high-level* aminoglycoside resistance, so carriage without
resistance should be rare) — but an absence is not a bound, and I am not going to promote it into one.

## The leakage gate earned its keep, and my cheap version did not

The first pass filtered on "never appeared in a `selected.tsv`" and found 956 candidates / 294 labelled /
**5** powered cells. Run against the real gate (`cohort_manifest.prior_accessions`, which also scans the
parquet cohorts) the pool lost **645 accessions — two thirds** — leaving 311 / 200 / **3** cells. E. coli
ceftriaxone and tetracycline dropped out entirely once leakage was removed.

A hand-rolled exclusion check beside the data under-covers. The fail-closed manifest exists so nobody has
to re-derive it, and this is the third time that lesson has been paid for in this repo.

**A near-miss worth recording:** my edit applying the manifest gate silently no-op'd (an unasserted
string replacement), and the re-run printed byte-identical numbers. I nearly reported the ungated 5-cell
result as gated. What caught it was that identical output after a real change is itself a signal.

## Status

- The `rmt` gap is now **quantified on independent data**: it costs the frozen rule 0.369 of sensitivity on
  E. coli gentamicin.
- The candidate recovers essentially all of it (0.892 vs the frozen cohort's 0.893) with no measurable
  specificity cost.
- **Still not deployed.** Changing the frozen surface invalidates the prospective lock and the
  reproducibility freeze. A v2 lock is a user authority call — and the evidence for making that call is now
  substantially stronger than it was two runs ago.

Reproduce: `uv run python scripts/unscored_genome_label_census.py` (network, writes the labels sidecar),
then score offline against `dna_decode/eval/amr_rules.call_resistance`.
