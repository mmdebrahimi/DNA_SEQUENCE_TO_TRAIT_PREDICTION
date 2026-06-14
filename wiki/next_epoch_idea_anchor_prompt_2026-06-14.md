# Next-epoch `/idea-anchor` prompt — first decode move on a CONFIRMED non-circular label source (2026-06-14)

**Supersedes** `wiki/next_epoch_idea_anchor_prompt_2026-06-13.md`. That one framed the question as "does a
non-circular label source exist?" — **the `/research` run (2026-06-14) + `--advance` answered it: YES.**
The Oxford E. coli cohort is a confirmed free, open-deposit, sampling-independent, provenance-disjoint
label source — the labels wall is broken. The question is now **which decode move to make first**, not
whether one is possible.

**Why this exists:** the AMR decoder is banked/frozen at v0.5.0; every public-label expansion was closed as
circular/confounded/saturated. The 2026-06-14 research found the first clean external label source since
AMR MIC itself. This anchor frames the project's reopening around it.

---

## What changed (the facts the anchor is built on)

- **TOP candidate — Oxford E. coli cohort (re-validation substrate).** Lipworth et al., Lancet Microbe
  2025 / medRxiv 2024.05.15.24307162: **2,875 isolates** with linked WGS + **measured clinical MIC** across
  **8 drugs incl. ciprofloxacin / ceftriaxone / gentamicin** (direct overlap with the frozen decoder's
  cells). UK Oxford University Hospitals provenance — **outside** the US NARMS/CDC/GenomeTrakr ecosystem the
  decoder was tuned on. Blood-culture WGS (~2,410 of the 2,875) is **openly deposited at ENA/NCBI BioProject
  `PRJNA604975`** (verbatim data-availability statement, Genome Medicine 2021 / PMC8414751) → **free
  download, NOT MTA**. Clears all 4 criteria; criterion-3 (downloadable assemblies) CONFIRMED.
- **RUNNER-UP — ecoref / eLife panel (new-trait substrate).** Galardini et al. 2017: **894 strains**, growth
  measured across **214 conditions**, **696 public genomes** (ecoref repo). A genuinely NEW, non-AMR,
  sampling-independent decodable trait (growth/fitness) at depth, free download.
- **DISQUALIFIED (recorded so they aren't re-investigated):** Pfizer ATLAS (6.5M measured MIC but NO genome
  assemblies → fails criterion 3); NARMS raw MIC (IS the excluded surveillance ecosystem → fails criterion 4).
- Full scorecard: `research_outputs/acquirable-ecoli-phenotype-label-sources-2026-06-14.md`.

---

## Paste-ready command

```
/idea-anchor First decode-epoch move on the Oxford E. coli cohort — external re-validation of the frozen v0.5.0 AMR decoder on independent, clinically-measured MIC.

PROBLEM: The deterministic E. coli AMR decoder is validated + frozen at v0.5.0 on PROVENANCE-DISJOINT but US-surveillance-adjacent NCBI-PD cohorts (10 SCORED cells, lineage-disclosed). It has never been tested on a fully INDEPENDENT, non-US, clinically-measured-MIC cohort. The 2026-06-14 research found exactly such a source — the first clean non-circular label source since AMR MIC itself.

THE SUBSTRATE (confirmed this session): the Oxford E. coli cohort (Lipworth et al., Lancet Microbe 2025) — 2,875 isolates with linked WGS + MEASURED clinical MIC across 8 drugs incl. cipro/cef/gent; UK OUH provenance (outside NARMS/CDC); blood-culture genomes (~2,410) OPENLY deposited at ENA/NCBI BioProject PRJNA604975 (free download, not MTA). It clears all 4 non-negotiable label criteria (sampling-independent / non-circular / ≥100-150 same-organism with downloadable assemblies / provenance-disjoint-feasible).

GOAL: scope the SMALLEST credible first move — run the SHIPPED deterministic call_resistance() rules (cipro/cef/gent — the cells that overlap) against the Oxford cohort's measured MIC, on the openly-deposited genomes, and report acc/sens/spec on a genuinely independent external cohort with the clonality/lineage discipline already built. A PASS = the strongest possible trust upgrade for the frozen decoder (it generalizes outside its tuning provenance). A FAIL = a real, honest finding about where it breaks. This is RE-VALIDATION of the existing tool on independent labels — NOT a new classifier, NOT embeddings, NOT a new trait.

OPEN FORKS to resolve in framing:
  A. PRIMARY = re-validate the existing decoder on Oxford (executor-doable, free, highest trust-value). Is this the target?
  B. SECONDARY = a NEW non-AMR trait on the ecoref panel (894 strains / 214 growth conditions / 696 public genomes) — new science, but growth phenotype may be lineage-confounded (the embedding-0-for-4 risk). Defer behind A, or pursue in parallel?

VERIFICATION STILL NEEDED (pre-build, cheap — none blocks the genome download):
  - Is the MIC value table in a downloadable supplement, or must it be extracted from the paper?
  - Exact E. coli count with BOTH genome AND MIC after the join (~2,410 blood expected).
  - The 465 urine isolates (2020 subset) accession — may be a separate/later deposit than PRJNA604975.
  - Does the Oxford MIC method (clinical instrument vs broth microdilution) match the decoder's CLSI/EUCAST breakpoint assumptions in mic_tiers.py?

SUCCESS CRITERION the downstream chain should pin: a per-drug (cipro/cef/gent) sens/spec on a leakage-checked, clonality-corrected Oxford slice, with an honest PASS/FAIL verdict + the provenance caveat (different country/lab/method than the tuning data).

CONSTRAINTS: solo hobby; two machines (laptop + Precision 7780); money needs approval (this path is FREE — open ENA download + local AMRFinder via Docker, same toolchain as the existing 10 cells); north star is a DECODER TOOL not papers — re-validation matters because it tells the user whether to TRUST the shipped tool outside its tuning provenance. Do NOT reintroduce circular/sampling-defined labels; do NOT restart embeddings/MIC-continuous (closed with recorded negatives).
```

---

## What to expect / how it chains

- `/idea-anchor` stops after framing. Likely ≤3 foundational questions: (i) confirm PRIMARY = Oxford
  re-validation vs new-trait; (ii) re-validate the cipro/cef/gent overlap only, or all 8 Oxford drugs (4 are
  new cells); (iii) acceptance bar for "generalizes" (e.g. sens/spec within X pp of the frozen cells).
- Recommended next step it will likely emit: **`/probe`** (repo-grounded — this names the decoder's own
  `call_resistance`, `mic_tiers.py`, `cohort_manifest.py`, Docker AMRFinder path) → then `/technical-plan`
  the Oxford-slice pilot.
- **Honest flag:** unlike the 2026-06-13 anchor, this fork's bottleneck is NO LONGER acquisition — the
  genomes are free + open. The bottleneck is now ordinary executor work (download a slice from PRJNA604975,
  join to the paper's MIC, run AMRFinder + call_resistance, score). That makes this the first genuinely
  executor-doable forward move since v0.5.0 — the project is back in motion.
