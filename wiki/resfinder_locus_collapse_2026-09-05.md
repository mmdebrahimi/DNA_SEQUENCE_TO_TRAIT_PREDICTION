# `dna-resfinder` reported 165 β-lactamase genes per genome. The real number is about two.

The fourth `FAITHFUL_TO_TOOL` cell to be pointed at real data, and the second live shipped defect the
exercise has found. This one is clinically meaningful rather than cosmetic.

## What it did

One cached E. coli genome, run through the shipped caller:

| | β-lactam | aminoglycoside |
|---|---|---|
| AMRFinder | **3** (`blaOXA-1`, `blaTEMp_G162T`, `blaTEM-1`) | **3** |
| `dna-resfinder` | **190** | **19** |

AMRFinder's calls were a *subset*. The genome carries `blaTEM-1B` at 100.0% identity / 100.0% coverage —
a narrow-spectrum penicillinase. The shipped caller **also** reported `blaTEM-52B`, `blaTEM-52C`,
`blaTEM-12`, `blaTEM-10` and `blaTEM-24`, at 99.0–99.8% identity.

Those are ESBLs. A reader of that output would conclude the genome encodes extended-spectrum
β-lactamase activity. It does not. **That is a wrong clinical reading of a genome, not a counting
problem.**

## Mechanism — exact, and the third instance of one failure

TEM β-lactamases differ by one to three point mutations; the ESBL phenotype comes from a couple of
amino-acid substitutions. At the nucleotide level `blaTEM-52` is >99% identical to `blaTEM-1`. The
caller's bar is 90% identity / 60% coverage, so **every TEM in the catalog clears it against a single
TEM locus** — and the output was keyed on the *allele name*, so each cleared independently and all were
reported present. Measured on the reference genome: **191 called alleles resolving to 2 actual loci**,
183 of them stacked on one.

This is the same underlying failure as the two typing bugs found this week — closely-related alleles
cross-hybridize and the selection rule fails to pick a winner — in its third form:

| cell | what the rule did |
|---|---|
| `salmserovar` | picked the **wrong** winner (coverage-only tiebreak) — fixed previously |
| `serotype` | picked the **wrong** winner (same defect, never propagated) — fixed 2026-09-04 |
| **`resfinder`** | **picked no winner at all — reported every candidate as present** |

## The fix: collapse by position, not by name

**Not by gene name.** `blaOXA-1` and `blaOXA-48` share a prefix, are functionally unrelated
(narrow-spectrum vs carbapenemase), and a genome can genuinely carry both. Name-based grouping would
merge two real genes.

**By locus.** Alleles hitting the same genomic position are one gene; the best-matching allele there is
the call. Selection is **identity-primary** — within a locus every variant sits at ~100% coverage, so a
coverage-first tiebreak is settled by dict iteration order rather than by sequence, which is precisely
the serotype defect.

Clustering is **greedy-representative, not single-linkage**. Single-linkage chains, and this project has
already seen a real 7-copy tandem `blaTEM` array; each hit is compared against its cluster's
highest-identity representative, so a tandem copy that does not overlap the representative starts its
own locus. `n_alleles_at_locus` ships in the output, so the collapse stays auditable.

On the reference genome the fix returns exactly what AMRFinder returned: `blaOXA-1` + `blaTEM-1B`
(191 alleles → 2 loci) and `aph(3'')-Ib` + `aadA1` + `aph(6)-Id` (33 → 3).

## Measured on 648 genomes

Every cached assembly that also has a committed AMRFinder run. Both rules scored from **one blastn pass
per genome**, so only the grouping differs.

| | β-lactam | aminoglycoside |
|---|---|---|
| mean genes/genome **before** | **165.52** | 9.16 |
| mean genes/genome **after** | **1.59** | 2.44 |
| AMRFinder | 1.94 | 2.46 |
| Jaccard vs AMRFinder, before → after | **0.0125 → 0.7754** | **0.4126 → 0.7786** |
| genomes where one locus was multi-reported | **522 / 648** | 400 / 648 |

**Verdict `LOCUS_COLLAPSE_IMPROVES_AGREEMENT_WITH_AN_INDEPENDENT_CALLER`** — agreement rises in every
class measured, on identical inputs. The defect was not an edge case: it fired on **80% of genomes** in
the β-lactam class.

**Reported, not buried:** after the fix the caller sits slightly *below* AMRFinder (1.59 vs 1.94
β-lactam genes/genome), so some under-calling remains. And the β-lactam **exact** Jaccard is 0.4269
against the normalized 0.7754 — the gap is naming (`blaTEM-1B` vs AMRFinder's `blaTEM-1`). The
normalized figure is the **lenient** reading and both are in the artifact.

## The sibling cells were measured, not assumed

`disinfinder` shares the defective pattern *literally* — it imports `resfinder.gene_of` and keys on it
the same way. It was still measured before being touched, and it is **inert**: the disinfectant DB holds
**16 alleles** with no dense variant families, and across 40 genomes old equals new (34 loci) with **0**
multi-reported loci. **Left unchanged**, matching how `plasmid` and `pneumoserotype` were handled — but
re-measure if that DB ever grows. `plasmid` already groups by `replicon_family` (a real family function)
and `pointfinder` has a different architecture entirely.

One real bug out of four suspects, again. **The pattern-match is a lead, not a diagnosis.**

## Honest limits

- **The comparator is a TOOL, not a wet-lab label.** This measures agreement with an independent curated
  implementation, not correctness — both callers could be wrong together. **The cell stays
  `FAITHFUL_TO_TOOL`**; what changed is which tool it is faithful to, and how closely.
- Only two ResFinder class DBs are committed (aminoglycoside, β-lactam). This says nothing about the
  other classes the full ResFinder DB covers.
- Genomes are whatever this project had cached, drawn from AMR cohorts — **enriched for resistance**, not
  a random sample.
- AMRFinder `POINT` rows are excluded from the comparison: a ResFinder allele DB cannot represent a point
  mutation, so scoring against them would charge the caller with missing something it cannot express.
- The 50% reciprocal-overlap bar means a tandem array whose copies overlap more than that would still
  collapse into one call.

## Reproduce

```bash
uv run python scripts/resfinder_locus_collapse_validate.py
```

Needs blastn + cached assemblies + committed AMRFinder runs. Frozen AMR surface byte-unchanged and the
2026-08-31 prospective lock re-verified — this is a finder cell, not the frozen decoder.
See [`resfinder_locus_collapse_2026-09-05.json`](resfinder_locus_collapse_2026-09-05.json).
