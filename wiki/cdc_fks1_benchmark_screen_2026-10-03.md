# CDC FKS1 *C. auris* benchmark: reachability INCOMPLETE, but the screen found a catalog gap offline (2026-10-03)

**Verdict: `G5_INCOMPLETE` on reachability — and that is the honest status, not a pass.** I could not
retrieve the per-isolate accession list or MIC table with the tools available here. But the screen was
not wasted: it produced a **concrete, offline-verified gap in our own deployed FKS1 catalog** that is
actionable regardless of whether the cohort is ever fetched.

## The dataset

*A benchmark dataset for validating FKS1 mutations in Candida auris* — Misas et al., **Microbiol
Spectr** 2025, [doi:10.1128/spectrum.03147-24](https://journals.asm.org/doi/10.1128/spectrum.03147-24),
open access at **PMC12323592**. CDC Mycotic Diseases Branch.

**100 WGS *C. auris* isolates: 53 echinocandin-susceptible / 47 resistant**, four clades (I–IV),
analysed with MycoSNP-nf. For comparison, our existing AR Bank micafungin cohort is **32 (6R/26S)** —
so this is roughly **8× the resistant class**.

## Reachability — what actually blocked, named precisely

| route | result |
|---|---|
| publisher full text (`journals.asm.org`) | **HTTP 403** |
| PMC (`pmc.ncbi.nlm.nih.gov/articles/PMC12323592/`) | cookie-consent wall; no article body returned |
| NCBI BioProject search (E-utilities, 3 query forms) | **count 0** — not discoverable by title |

The likely explanation for the BioProject miss is that these isolates sit under an existing CDC
surveillance umbrella project rather than a new one, with the accession list in the supplementary
table. **"Three routes failed" is evidence about three routes, not a property of the source** — this
repo has that lesson on record from the ECOFF episode, where six failed routes were wrongly reported as
a human-only wall and the real cause was a wrong query parameter. So this is INCOMPLETE, not blocked.

**A reachability risk worth flagging before anyone spends labour:** the paper's own figure descriptions
put per-isolate MICs in a **heatmap (Fig. 1)** and **boxplots (Fig. 2)**. If per-isolate MIC values
exist only in figures, that is the PEAR problem again — numbers locked in plots — and the extraction
cost is real. Unverified either way.

**Unblock:** a human opening PMC12323592's supplementary files and reporting (a) the BioProject/SRA
accession and (b) whether a per-isolate MIC + FKS1-genotype table exists as a file. Everything
downstream already exists: `scripts/build_fungal_cohort.py`, `scripts/assemble_sra_cohort.py --method
map` (targeted ERG11/FKS1 read-mapping, ~4 min/isolate), and the committed `data/fungal_ref/` CDS
references.

## What the screen DID establish, offline and verified

**Our deployed FKS1 catalog carries three mutations, and all three are at one position.**

```
FUNGAL_RESISTANCE_MUTATIONS -> FKS1 entries: S639F, S639P, S639Y
```

Searched the module for `W691`, `689`, `691`, `696`, `HS3`, `third` — **all absent**.

The benchmark paper reports resistant isolates carrying non-synonymous FKS1 hotspot mutations in
**HS1 (F635–P643)**, **HS2 (D1350–L1357)**, and a **potential third hotspot at W691 (L689–N696)**, with
44 of 47 resistant isolates carrying one, and it names **D642Y** among the commonest.

So our catalog is **HS1-S639-only** and would structurally miss:

- other HS1 positions in F635–P643, including the named **D642Y**
- **all of HS2** (D1350–L1357)
- the **third hotspot** (W691)

> **Scope care:** our catalog's contents are **verified offline** (read from the module). The hotspot
> coordinates and mutation names are taken from a **secondary summary** of the paper, not from the paper
> itself, which I could not fetch. Treat the region boundaries as paper-reported-via-search until the
> full text is read. The gap in *our* catalog does not depend on them — S639-only is S639-only.

This is the fungal arm's analogue of the completeness question the L2 doubt layer exists for, and the
record currently lists the fungal cells as UNMEASURED because no free phenotype source was known. A
catalog carrying one position out of three hotspot regions is a **predicted blind spot** that can be
stated now and tested the moment a cohort lands.

## Gate screen, as far as it can be taken

| gate | status |
|---|---|
| G1 circular label | **PASS** — AFST MIC is a wet-lab reading, not a genotype-tool call |
| G2 study == class | **unscreened** (needs the accession/source table) |
| G3 sampling-defined | **PASS** — the label is an MIC, not a collection context |
| G4 surveillance domination | **RISK** — it is a CDC surveillance-derived panel by construction |
| G5 assembly attrition | **INCOMPLETE** — the blocker above |
| G6 censoring | **unscreened** — echinocandin MICs commonly report as `>8`/`<=0.06`; operator-aware handling exists (`external_mic_labels.MicValue`) but the censored fraction is unmeasured |
| G8 clonal collapse | **RISK, bounded** — 4 clades named, but *C. auris* is clade-structured and clonal within clades |

**A PASS on the label question would not make this a build recommendation.** Per the standing rule,
artifact reachability and regime fit are separate questions — and reachability is precisely what is
open.

## Recommendation

Do **not** spend labour on this cohort until the two reachability facts are in hand. **Do** record the
S639-only catalog gap now: it is free, verified, and predicts where the fungal cell fails whether or not
this particular benchmark is ever used.
