# The Moderna/Merck neoantigen readout, screened against our own 8 rejection gates

**Date:** 2026-08-23 · **Type:** external-landscape + candidate-substrate screen
**Trigger:** user-supplied video — *"How Big of a Deal is this Cancer Vaccine?"*, Hank Green
(`youtube.com/watch?v=DgeRVUSuGvQ`), full caption track retrieved and read.
**Screens against:** `wiki/negative_results_map_2026-06-13.md` → the 8 rejection gates.

This is the first external-landscape doc in the repo. It exists because the Merck/Moderna result looks
like a counterexample to this project's central finding, and it is not — it is a **confirmation of the
regime boundary**, from the opposite direction.

## 1 · What was actually announced (and two things circulating that are wrong)

**INTerpath-001**, Phase 3, announced 2026-08-19: **intismeran autogene** (= V940 = mRNA-4157) +
**pembrolizumab** in completely-resected stage IIB–IV melanoma. Met primary endpoint (recurrence-free
survival) and a key secondary (distant metastasis-free survival). First positive Phase 3 for an
individualized neoantigen therapy and for any mRNA cancer therapy.

**Correction 1 — the hazard ratios everyone is quoting are Phase 2b, not Phase 3.** The Phase 3 topline
does **not** disclose an HR. `HR 0.51` (RFS) and `HR 0.411` (DMFS) are the 5-year KEYNOTE-942 Phase 2b
numbers from ASCO 2026.

**Correction 2 — "EchoNeo" is not verifiably Moderna's.** It is AACR 2026 Abstract 4376 (Liu, Mao, Wu,
Qiu); no Moderna affiliation appears in the primary source. A third-party blog made the attribution and it
has propagated. The abstract's *content* is still worth having: it states that current neoantigen
strategies "rely primarily on MHC binding affinity, leading to limited accuracy in immunogenicity and a
high false-positive rate."

**The pipeline**, as far as public sources describe it: tumour + blood NGS (**DNA and RNA**, so expression
informs design; blood = germline comparator) → ML ranks candidate neoantigens by predicted MHC
presentation and immunogenicity **against that patient's HLA type** → up to **34** epitopes concatenated
onto a single mRNA → LNP. Moderna's surrounding infrastructure: **mRNA Design Studio** (AWS-hosted
sequence designer feeding automated production) and **"Lucy"** (July 2026, closed-loop ML + lab
automation). Exact model architecture, features and training data are proprietary.

**Not AI-invented.** The video's framing is the accurate one: a therapy that *leverages* ML at one step,
not one an AI thought up.

## 2 · The screen — the gates that closed our tracks are the ones this substrate passes

The ML step is a genotype→phenotype problem: mutation → peptide → will MHC-I present it. Its label source
is **mass-spec immunopeptidomics**. Screened against our 8 gates:

| gate | verdict | evidence |
|---|---|---|
| **G1 circular label** | **PASS** | the label is a physical measurement — peptides eluted off MHC and read by LC-MS/MS — not a gene-call from a tool the model would compete with. *Caveat: peptide ID uses a database search engine, so FDR/search choices shape which peptides are called.* |
| **G2 study == class** | **PARTIAL TRIP** | not the classic form, but the shape is inverted and worse: positives are **observed** peptides, negatives are **assumed-unobserved** decoys. Detection limits, not biology, define the negative class. |
| **G3 sampling-defined label** | **PASS (cleanly)** | the label is an assay reading, not a description of where the sample came from. This is the gate that killed pathotype, and it does not fire here. |
| **G4 surveillance domination** | **PASS, with an allele caveat** | multi-lab aggregation (SysteMHC v2.0: 2,447 samples, 303 allotypes). But the field's own known limit is *insufficient MS data for many HLA alleles* — per-allele, not per-corpus, is where it thins. |
| **G5 assembly attrition** | **N/A — passes trivially** | no genome assembly is needed. The 96%-drop failure mode that dominated our BV-BRC work simply does not exist here. |
| **G6 phenotype censoring** | **PARTIAL TRIP** | no MIC-style interval censoring, but the negative side *is* the censored side — non-detection ≠ non-presentation. Same defect as G2, seen from the other end. |
| **G7 provenance not separable** | **LIKELY PASS (unverified)** | HLA Ligand Atlas carries tissue + donor annotations; SysteMHC carries sample/study. Donor-disjoint and allele-disjoint splits look constructible. Not verified against the actual files. |
| **G8 dedup collapses balance** | **NEEDS CHECKING** | the analogue of lineage collapse is peptide/source-protein redundancy — nested and overlapping peptides from the same protein region. The field handles this with source-protein-level splits; we have not measured it. |

**Free and reachable** (all HTTP 200, checked): IEDB · HLA Ligand Atlas (90,428 HLA-I ligands, 51 HLA-I
alleles, 30 benign tissues) · SysteMHC v2.0 (**~1.0M unique class-I peptides, 303 allotypes**; 2,028,964
unique peptides total at 1% FDR) · MHCflurry · NetMHCpan.

## 3 · What this settles

> **This is the "PDB for genomics" that our own analysis said does not exist — except it exists for one
> narrow molecular question, and that is exactly where the learned model wins.**

Set against the project's regime boundary:

| precondition | our AMR / essentiality work | MHC-I presentation |
|---|---|---|
| large independent wet-lab label corpus | **absent** — 11 of 27 decoder cells have *no free phenotype source* | **~1M peptides across 303 allotypes** |
| label from an assay, not a sampling context | often violated (G3 closed pathotype) | **satisfied** |
| population-structure / clonality confound | **everywhere** — TB raw 0.92 → lineage-collapsed 0.44 | **absent at the peptide level** |

So the Merck/Moderna result is not evidence that our conclusions were wrong. It is evidence that the
binding constraint we identified — **labels, not models** — is the correct one, because where the labels
exist, the learned model works.

**Two further parallels worth recording:**

1. **Their pipeline fails where ours does, for the same reason.** *Presentation* has a large MS corpus and
   is well predicted. *Immunogenicity* — will a presented peptide actually provoke T cells — has small,
   biased labels, and that is where the false-positive rate lives (the AACR abstract says so outright).
   BioNTech's pancreatic trial is that gap in clinical form: **8 of 16 patients mounted no immune response
   at all**, and only the 8 who did had better outcomes.
2. **The video's "one last wrench" has since been challenged, and the challenge is our own lesson.** The
   2025 Nature observational result (NSCLC 3-year OS **55.7% vs 30.8%** for mRNA-COVID-vaccinated patients
   on checkpoint inhibitors) was reanalysed by MSK (Jee et al., **8,368** ICI patients): the benefit was
   **not specific to immunotherapy** (it appeared with chemotherapy too), was strongest early-pandemic and
   faded after 2021, and **shrank or vanished under landmark analysis** — read as selection bias, healthier
   patients being likelier to get vaccinated. A dramatic effect that survives 39 covariate adjustments and
   dies under a better-designed split is structurally identical to the clonality inflation this project has
   spent a week correcting for. A randomized Phase 3 is being designed.

## 4 · Recommendation: do NOT build here — and the reason is not the gates

The screen comes back mostly clean, which makes the honest recommendation less obvious and worth stating
plainly: **this is not a good direction for this project.**

- **It is not the north star.** The stated goal is an *AI DNA decoder*. MHC-I presentation is
  peptide→presentation — protein-level, downstream of the DNA question, and adjacent at best to the
  existing `forward/` variant-effect cell.
- **The field is mature and we would bring no data advantage.** NetMHCpan, MHCflurry, SHERPA, neoMS,
  MHCnuggets and others are established on the *same* public corpora. Entering with the same data and less
  domain infrastructure is not a differentiator.
- **The genuinely open sub-problem is immunogenicity, and it fails our gates** — small, biased,
  assay-inconsistent labels. That is the same wall, one layer up.

**What it is good for:** it is a **positive control for the negative-results map.** Every prior use of
those 8 gates was to reject something. This is the first worked example of what a substrate looks like when
it *passes* — which makes the map a usable screen in both directions rather than only a rejection list.

## 5 · Honest limits

1. **The video is a secondary source.** Everything load-bearing above was checked against Merck/Moderna
   primary releases or the underlying papers; the video contributed framing and the pancreatic/BioNTech
   pointer, not facts.
2. **Moderna's actual model is proprietary.** The pipeline description is assembled from press material,
   an AWS case study and trial documentation. No architecture, feature set or training corpus is public.
3. **G7 and G8 are unverified** — asserted from dataset documentation, not measured against downloaded
   files. If this were ever pursued, those two are the first things to check, not the last.
4. **The MSK reanalysis is itself one preprint**, not a settled overturn. Both the original and the
   reanalysis are observational.
5. **No claim is made about clinical efficacy.** The Phase 3 HR is not public.

## Sources

Merck press release · Moderna newsroom · Merck INTerpath backgrounder (PDF) · AACR 2026 Abstract 4376 ·
Moderna mRNA platform + AI blog · AWS Moderna case study · SysteMHC Atlas v2.0 (NAR 2024) · HLA Ligand
Atlas (PMC8054196) · SHERPA / composite MHC presentation modeling (MCP 2023) · MD Anderson ESMO 2025 ·
MSK reanalysis (medRxiv 2026-01) · YouTube `DgeRVUSuGvQ` caption track.
