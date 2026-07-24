# Free label-acquisition — deep-research round 2 (2026-07-23)

**Status:** DECISION AID (research-only; the build is a user go/no-go + a lane-coordination call). **Method:**
4 parallel web-research agents, each screening candidates against the 8 rejection gates
(`wiki/negative_results_map_2026-06-13.md`) in territory the prior round (`wiki/data_acquisition_voi_memo_2026-07-18.md`)
did NOT deeply cover. **Grounding:** the winning free-independent paradigm = a curated **target-site / determinant
catalog** scored against a **FREE, INDEPENDENT, WET-LAB-MEASURED** phenotype (HIV HIVDB fold-change, SARS-CoV-2
CoV-RDB, CRyPTIC TB MIC) — NOT a rules-interpreter's R/S call.

## Headline: the "free sources are exhausted" verdict is REFUTED (scoped)

The 2026-07-18 memo was **correct for its scope — bacterial-AMR-MIC surveillance** (BV-BRC / NCBI-PD / KlebNET
all gate-trip; and this round CONFIRMS the surveillance-MIC platforms Vivli/ATLAS/SIDERO are dead too — they carry
**no linked genomes**, G5). But three adjacent territories it did not sweep each yield genuinely-new, gate-clearing,
free label sources:
1. **New pathogen families in the winning viral-catalog paradigm** — HCMV (herpesvirus) is a near-exact analogue of the shipped HIV/SARS-CoV-2 cells.
2. **Independent (CRyPTIC-disjoint) TB DST** — TB Portals gives the out-of-distribution TB number the project explicitly wanted.
3. **A whole new NON-AMR measured trait axis** — phage host-range + BacDive biochemical phenotypes, both measured + genome-linked, structurally exempt from the AMR label wall.

## Ranked candidates (all clear the free + measured + genome/genotype-linked bar unless noted)

| # | Candidate | Phenotype (measured) | Free? | Fit to paradigm | Verdict | Top caveat |
|---|---|---|---|---|---|---|
| **1** | **HCMV antiviral** (UL97/UL54/UL56 → GCV/CDV/FOS/letermovir) | recombinant marker-transfer **fold-change EC50** (Chou compilations, open-access PMC/ASM) | **YES** (open PMC tables) | **near-exact** — target-site catalog + measured FC + built-in phenotyped-benign sensitive class | **GO (top)** | catalog is dozens–hundreds of mutations (curate from tables), not a single download |
| **2** | **TB Portals (NIAID)** | MGIT/LJ **measured DST** + WGS in SRA (no DUA), joinable by `condition_id` | **YES** (open) | reuses shipped TB cell → **independent** out-of-distribution number | **GO** | per-isolate CRyPTIC/WHO-v2 disjointness must be BioSample-checked (tooling exists) |
| **3** | **Phage host-range** (Picard/Guelin, Nat Microbiol 2024) | **38,688 wet-lab lysis interactions** (403 sequenced *E. coli* × 96 phages) | **YES** (open + Zenodo/BioProject) | NEW non-AMR axis; matches AMR track's within-*E. coli* lineage depth (163 STs, 8 phylogroups) | **GO (new axis)** | receptor/adsorption catalog is the baseline (their genome model AUROC 0.86) → classify regime first |
| **4** | **Antimalarial lab-evolution** (Winzeler 2024 `adk9893`; Cowell 2018 `aan4472`) | **measured EC50 fold-shift** (448 + 262 clones vs matched parents), WGS in SRA | **YES** (open, no DUA) | partial — mechanism-discovery shape | **GO w/ caveats** | mostly EXPERIMENTAL compounds (not clinical kelch13/pfcrt drugs); limited lab-founder lineages (G8) |
| 5 | **BacDive (DSMZ)** | measured biochemical/physiological +/− (enzyme, carbon-util, growth-temp) on 50,588 linked genomes | **YES** (open REST API, no DUA — cleanest access) | NEW non-AMR axis | **MARGINAL-GO** | type-strain-dominated → verify within-single-organism phenotype∩genome depth |
| 6 | **HCV replicon fold-change** (NS3/NS5A/NS5B) | subgenomic-replicon **FC** (drug labels + in-vitro papers) | YES (scattered) | good | **MARGINAL-GO** | no clean download — manual extraction; subtype-stratified |
| 7 | Pooled post-2023 TB (India PZA `PRJNA1155695`, Peru, Singapore) | measured pDST + free SRA | YES | independent top-up | MARGINAL→GO in aggregate | individually small / drug-patchy |

## Confirmed CLOSED (do not pursue as a label; reasons verified this round)

- **Surveillance MIC — Vivli AMR Register / Pfizer ATLAS / SIDERO-WT / SENTRY:** measured MIC but **NO linked genomes** (β-lactamase gene flags only) → **G5**. Worse than BV-BRC.
- **AllTheBacteria (2024):** AMR labels are **AMRFinderPlus calls** → **G1 circular**. BUT keep as a free 2.44M-genome **POOL** (infrastructure to fetch genomes for any measured label obtained elsewhere).
- **WWARN in-vitro IC50:** content clears the gates but is **DUA-gated** (DAC review + signed agreement) → not free-open.
- **MalariaGEN Pf7/Pf8:** genotype-only; "resistance" is marker-INFERRED → **G1 circular**.
- **WHO TB catalogue v2 dataset:** it IS the project's TB rule (CRyPTIC + Seq&Treat) → **circular / not independent**.
- **ReSeqTB:** platform discontinued (~2021) — dead.
- **geno2pheno[hbv]/[hcv] interpretations, HSV, HBV:** interpreter output = **G1 circular**; HSV resistance is loss-of-function (TK-null) → shape-mismatch for a substitution catalog; HBV measured data too scattered.

## Recommendation

**Build #1 — the HCMV antiviral cell — as the highest-VOI, lowest-risk next capability.** It reuses the EXACT
architecture already shipped + validated for HIV-1 and SARS-CoV-2 (curated target-site catalog → free open-access
wet-lab-measured recombinant fold-change → R/S with a genuine phenotyped-benign sensitive class), extends the decoder
to a NEW pathogen family (herpesvirus), and clears every gate cleanly. Lowest build risk because the cell shape is
proven.

**Then #2 — TB Portals** — to produce a genuinely INDEPENDENT (CRyPTIC-disjoint) number for the SHIPPED TB decoder,
reusing the existing TB cell + `biosample_resolver`/`cohort_manifest` disjointness tooling. High scientific value
(external validation of a deployed claim), moderate build.

**Optional new-axis bet — phage host-range** — the one genuinely NEW measured phenotype with AMR-track-grade
within-*E. coli* depth. Worth a cheap regime-classification probe first (their genome-only model already hits 0.86
AUROC → likely a "curated-catalog-wins-again" regime, not the embedding-friendly one — but that itself is a clean,
publishable determinant-decoder result).

## Honesty rails

- **Unverified numbers** (per agents): exact HCMV catalog size; TB-Portals WGS∩DST isolate count + its per-isolate
  CRyPTIC disjointness; the Zenodo DOI for the phage data; antimalarial per-clinical-drug R/S counts; BacDive
  within-species depth. Each is a cheap pre-build check, not a blocker.
- **Lane coordination (R4):** the shipped viral cells (HIV/SARS-CoV-2) live in this repo; a new HCMV cell + the TB
  Portals independent run touch the decoder surface — coordinate with DNA-11 (decoder cells) / the acquisition lane
  before building. This memo is a decision aid; the build is a **user go/no-go**.
- This does NOT reopen the closed negatives (embeddings 0-for-4; MIC-continuous; organism-polygenic). It adds
  NEW label substrates in the DETERMINISTIC-catalog paradigm that already wins.

## Sources (per-candidate URLs in the agent transcripts)

HCMV: PMC3262590 / PMC3773841 / AAC 10.1128/aac.49.7.2710 / AAC 10.1128/aac.00922-18 · TB Portals:
tbportals.niaid.nih.gov + datasharing.tbportals.niaid.nih.gov · Phage: nature.com/articles/s41564-024-01832-5 ·
Antimalarial: PMC11809290 (Winzeler 2024) / PMC5925756 (Cowell 2018) · BacDive: bacdive.dsmz.de + NAR D1/D748 ·
Vivli/ATLAS (closed, G5): amr.vivli.org · AllTheBacteria (pool): allthebacteria.org.
