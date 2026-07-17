# The flowering cell's FLC route — the distinctive claim, tested (2026-07-17)

**Verdict: `FLC_ROUTE_VALIDATED_ADDS_OVER_FRI_ONLY`** — n=106.

The Table S3 scoring run could only test the **FRI route**, because S3 carries no FLC — so the cell
collapsed to exactly the naive FRI-only rule its Da(1)-12 anchor exists to catch, and its
distinctive two-locus claim went **untested**. This tests it, by joining FLC *expression* from an
independent free source (AraPheno phenotype 29, Atwell 2010) to Table S3's FRI status + flowering time.

## The rule's four cells — all called correctly

*late iff FT16 > 61.0d ; strong FLC iff expression > 0.857 (cohort medians)*

| FRI | FLC | n | % late observed | the cell calls | |
|---|---|---:|---:|---|---|
| functional | strong | 47 | 85% | `late` ✅ |  |
| functional | weak | 23 | 39% | `early` ✅ | **← Da(1)-12 class** |
| lof | strong | 6 | 17% | `early` ✅ | **← Lz-0 class** |
| lof | weak | 30 | 10% | `early` ✅ |  |

The **Da(1)-12 class** (functional FRI + weak FLC) is the cell's signature prediction: a naive
FRI-only rule calls these LATE, and only **39%** of them are. The **Lz-0 class**
(FRI-LoF yet late, via FRI-independent FLC upregulation) is real but rare here — 1 of 6 — which is
why the cell caps that branch at MEDIUM confidence rather than calling it wrong.

## Does the second locus earn its place?

Measured where a FRI-only rule commits: the **70 functional-FRI accessions, all of which it calls LATE**.

- strong FLC → **85% late** · weak FLC → **39% late**
- calls **rescued** by FLC: **14** · **broken**: **9** · net **+5**
- accuracy: two-locus **0.811** vs FRI-only 0.764 (+0.047); null 0.500

It is a **net** gain, not a clean one: FLC fixes 14 calls and breaks 9.

## The honest headline: within ancestry

FRI genotype tracks ancestry (the S3 run measured a +23pp pooled advantage collapsing to +3.4pp
within-ancestry), so the within-ancestry figure is the one that means anything:

| rule | within-ancestry accuracy | vs its own null |
|---|---:|---:|
| FRI-only | 0.767 | +0.016 |
| **two-locus (FRI + FLC)** | **0.803** | **+0.052** |
| (null = guess each group's majority) | 0.751 | — |

So across 5 ancestry groups the FLC route roughly **triples** the within-ancestry
advantage the FRI-only rule manages. That is the strongest honest statement here.

## The caveat that must ship with the number: it rides on the threshold

| FLC cut (quantile) | q20 | q30 | q40 | q50 | q60 | q70 |
|---|---|---|---|---|---|---|
| FLC's benefit | +0.028 | +0.066 | +0.047 | +0.047 | +0.000 | -0.085 |

**The benefit dies at q60 and reverses at q70.** The median cut is not biological: Werner 2005
reports weak/null FLC alleles are *rare*, which a median split cannot represent — so the plausible
range is the low quantiles, where the benefit holds (+0.028 to +0.066). But a reader must know the
gain is not threshold-robust: over-call weak FLC and the second locus becomes actively harmful.

## Scope

- IN-DISTRIBUTION: the cell's catalogue and both label sources trace to the same literature.
- FLC expression is a PROXY for allele status (see honest_mapping).
- n=105 -- the join of S3 (1,017) with AraPheno FLC (167).
- HABIT/direction only, not days-to-flower.

## Provenance

- **fri_status_and_flowering_time**: Zhang L & Jimenez-Gomez JM (2020) Plant Journal 103:154-165, doi:10.1111/tpj.14716, Table S3 (CC-BY 4.0)
- **flc_expression**: Atwell S et al. (2010) Nature 465:627-631 -- FLC expression, via AraPheno phenotype 29 (https://arapheno.1001genomes.org/rest/phenotype/29/values.csv)
- mapping: the cell wants an ALLELE status; AraPheno gives measured FLC EXPRESSION. Expression is a downstream PROXY for allele strength (a weak/null allele yields low steady-state mRNA -- how Michaels 2003 defined the weak alleles) but is NOT the same measurement
