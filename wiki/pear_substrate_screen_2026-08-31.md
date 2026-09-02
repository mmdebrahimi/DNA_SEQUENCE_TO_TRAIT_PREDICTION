# PEAR screened — it is not what F-E assumed, and the blocker is not a gate

**Verdict: PEAR is a constructed-variation DMS for the L4 forward cell, NOT a label source for the
deterministic AMR decoder — and NOT an acquisition. It is public and free. **It clears every rejection
gate that applies to it, G6 now included** — see the update note below. The blocker is ARTIFACT FORMAT.**

> **Update 2026-09-01 — G6 is CLOSED, and the artifact blocker is GONE.** This headline was corrected
> *downward* earlier today (to "incomplete") because G6 was unscreened and the mechanical screen said so.
> It is now corrected *upward*, on measurement rather than on assertion. R was installed on D: (4.1 TB
> free, which is what made the "C: at 99%" blocker evaporate), the authors' `.RData` was read, and the
> per-variant fitness values were extracted. **G6 passes:** mode-share 0.0019 over 2,106 distinct levels
> on 2,114 cefotaxime measurements — nowhere near the >25% / <20-level degeneracy bars. The ten-gate
> screen now returns **CLEARS**. Detail, including a trap that nearly produced the opposite verdict:
> `wiki/pear_g6_screen_2026-09-01.md`.

Screened 2026-08-31 per F-E's recorded discipline: *resolve accessions before recommending, never after.*
That discipline paid — but not in the way it was written for.

---

## What PEAR actually is (verified, not from memory)

Rhee-style deep mutational scan in all but name — [Zhang et al. 2022, *Mol Biol Evol* 39(5):msac086](https://academic.oup.com/mbe/article/39/5/msac086/6575838):

- ~200,000 constructed *E. coli* MG1655 strains, each carrying one variant of the **single-copy**
  `blaCTX-M-14` gene at a **fixed chromosomal site** (`att`HK022)
- **relative growth** measured for **~23,000** strains under cefotaxime and ceftazidime selection
- covers **>90% of single mutations** in the gene
- readout is barcode frequency by Illumina — a continuous fitness score, not an R/S call

## The reclassification, and why it matters

F-E's ledger scoped PEAR as a candidate **label source to clear the AMR label wall**. That is wrong on
two counts, and both change what should be done with it:

**1. It is not an L1 substrate.** These are constructed lab strains with a continuous growth readout —
not natural isolates with clinical AST. `dna_decode/eval/regime.py` classifies
`constructed + molecular + supervised` as **`WORKS`**, evidenced by the TEM-1 genome-edit path at
Spearman **0.761**. PEAR sits in that regime, so it **cannot** retire the binding uncertainty the
reproducibility freeze names (labels for the deterministic decoder).

**2. It is not an acquisition.** It is public: BioProject in SRA, code and figure data on GitHub, no DUA,
no money, no partner. The user-authority gate F-E attached to it was a consequence of the
misclassification, not of anything about the dataset.

**What it IS:** the natural external replication of the one working learned regime. Our forward cell was
validated on **TEM-1 + ampicillin**; PEAR is the same shape on a **different β-lactamase (CTX-M-14) with
different drugs (cefotaxime, ceftazidime)**. Same family, independent lab, independent measurement.

## Rejection-gate screen

Scored against `wiki/negative_results_map_2026-06-13.md` **as an L4 constructed-molecular substrate** —
the intended layer decides which gates even apply, so it is stated first.

| gate | verdict | why |
|---|---|---|
| G1 circular label | **pass** | wet-lab growth measurement; no genomic tool produced the label |
| G2 study == class | **n/a** | no class-vs-study confound — variation is constructed, not sampled |
| G3 sampling-defined label | **pass** | relative growth is an assay reading, not a collection context |
| G4 surveillance domination | **n/a** | not a surveillance corpus |
| G5 assembly attrition | **n/a** | no assemblies needed; variants defined by construction |
| G6 phenotype censoring | **OPEN — the one real risk** | see below |
| G7 provenance not separable | **n/a** | no provenance split needed; split by POSITION, as the forward cell already does |
| G8 dedup collapses balance | **n/a** | no clonality — ancestry randomised by construction, which is the point |
| G9 causal variant unrecorded | **pass** | every variant is known by construction |
| G10 variant class off-panel | **pass** | single substitutions in one gene, directly representable |

**PEAR clears the gates.** That is a genuinely rare result in this project — and it is exactly what the
regime map predicts, since the gates were written for natural-population label sources and constructed
variation sidesteps most of them by construction.

### G6 is the one to check first, and this repo has been burned by it

A selection-based growth assay has a **floor**: every non-functional variant dies at the same rate, so its
score piles up at one value. That is the censored-assay trap already recorded here — CcdB was **79.3%
tied at its ceiling** and posted the *best* number in the whole forward/inverse sweep because quantile
targets collapsed onto the tie. `assay_degeneracy()` in `scripts/forward_inverse_roundtrip.py` exists to
gate exactly this and **must be run on PEAR before any score is believed** (mode-share > 25% or fewer than
20 distinct levels → exclude).

## The actual blocker: artifact format, not availability

This is where the probe earned its keep. Verified end to end:

| artifact | status |
|---|---|
| BioProject `PRJNA687219` | **resolves** — *E. coli* K-12 MG1655 (taxid 511145), correct title, 45 SRA experiments, Sun Yat-sen University |
| SRA payload | **478 Gbases / ~0.2 TB of RAW READS** — barcode sequencing, not a fitness table |
| GitHub `woson2020/CTXM-14` | **exists**; entire repo is 2 `.RData` files + R/notebook scripts + the barcode-calling pipeline |
| plain-text data (CSV/TSV/parquet) | **none anywhere in the repo** |
| `Data_for_Figure2.RData` (3.3 MB) | **a serialized ggplot2 PLOT OBJECT**, not a table |

That last line is the finding. The file is uncompressed `RDX3`; `pyreadr` fails
(*"unsupported features"*) and `rdata` fails on `RObjectType.WEAKREF`. A string scan explains both: the
object is named **`Figure.2A`** and the file is full of ggplot2 internals (`legend.position`,
`position_dodge`, `strip.position`, `non_position_scales`, "Suppressing axis rendering when
strip.position…") around an axis label `"Relative growth"`. Plot objects carry environments and closures,
which is precisely what both parsers choke on.

The numbers are very likely embedded in the plot's data slot — but extracting them needs **R**, which is
**not installed on this host**, and `C:` is at **99% (3.0 GB free)**.

## Wall classification: CODE-CLOSABLE, not external

Nothing here needs a partner, a purchase, or a permission. Three routes, cheapest first:

1. **Supplementary tables** — UNRESOLVED. The data-availability statement could not be read (PMC is
   cookie-gated, the Oxford page truncates, Europe PMC is JS-rendered). Papers of this type usually
   deposit the variant-fitness table as supplementary Excel/TSV. **Check this before installing anything.**
2. **Install R** (~150–200 MB) and `saveRDS(ggplot_build(Figure.2A)$data, ...)` → CSV. Feasible but a Care
   concern at 3.0 GB free.
3. **Re-derive from raw reads** with their own pipeline — **infeasible here** (0.2 TB against 3.0 GB free).

## Honest limits

The reclassification is verified; the *usability* is not. I have not read one PEAR fitness value. Whether
the ggplot object carries all ~23,000 variants or only the subset plotted in Figure 2A is **unknown** —
and a figure-scoped subset would be a much weaker substrate than the paper's headline number suggests.
G6 is unscreened. Licence/reuse terms were not checked.

## What this changes

- **F-E's PEAR entry is reclassified**: not an acquisition, not an L1 label source, no authority gate.
- **It does not retire the label wall.** The binding constraint on the deterministic AMR track is
  unchanged, and PEAR was never going to change it.
- **It is a live, free, gate-clearing candidate for extending the forward cell** — pending the extraction
  route and a G6 degeneracy screen.
