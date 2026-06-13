# Next-epoch `/idea-anchor` prompt — new non-circular label source (2026-06-13)

**Why this exists:** As of 2026-06-13 every executor-reachable thread is terminal — AMR provdisjoint grid
SATURATED (10 SCORED cells), embedding bet CLOSED 0-for-4, eukaryotic cycle COMPLETE, pathotype
label-blocked, MIC-continuous closed (circular + censored + underpowered). The single recurring
meta-constraint behind all of it is **labels, not models**. Forward motion requires acquiring a *new*
phenotype label source that breaks the constraint — an acquisition decision only the user can make. This
file holds the drafted `/idea-anchor` to open that epoch (user-only skill → paste-ready).

---

## Paste-ready command

```
/idea-anchor New non-circular phenotype label source to reopen the decoder beyond the saturated AMR grid.

PROBLEM: The E. coli genotype→phenotype decoder is at a validated milestone but forward motion is blocked
by a single recurring meta-constraint — LABELS, not models. Every track died on the same wall: pathotype
labels are circular (EnteroBase=BlastFrost gene-calls, NCBI-PD=AMRFinder-derived) or study==class
confounded; foundation-model embeddings went 0-for-4 (learn lineage/population-structure, not mechanism,
even on the best-designed eukaryotic test); MIC-continuous is dead (91% of BV-BRC MIC is computational/
XGBoost-from-genome = circular, ~70% of the rest censored at the breakpoint); and the provenance-disjoint
AMR validation grid is SATURATED (10 SCORED cells; Salmonella + other deep organisms are
surveillance-ecosystem-dominated so the disjoint R-pool collapses). The ONE clean label the project ever
had — AMR MIC R/S from broth microdilution — is sampling-independent but exhausted at free-NCBI-PD scale.

GOAL: Identify and ACQUIRE a new phenotype label source that breaks this constraint — reopening either
(a) validation of the existing deterministic decoder on genuinely independent labels, or (b) a new
decodable trait — by clearing ALL FOUR criteria the project learned are non-negotiable:
  1. SAMPLING-INDEPENDENT — the label is a lab measurement / assay, NOT the sampling context itself
     (clinical isolation-site IS the confound; an MIC/growth/biochemical assay is not).
  2. NON-CIRCULAR — NOT derived by a genomic tool the decoder would compete against (no AMRFinder/
     BlastFrost/XGBoost-from-genome provenance); ideally wet-lab or clinically measured.
  3. ORGANISM-DEPTH — ≥~100-150 same-organism isolates carrying the label AND downloadable assemblies
     (the historical ~96% assembly-availability drop is the real bottleneck, not raw record count).
  4. PROVENANCE-DISJOINT-FEASIBLE — enough label-bearing isolates OUTSIDE the NARMS/CDC/FDA/GenomeTrakr/
     PulseNet/USDA surveillance ecosystem to build a leakage-clean ≥20/class cohort, AND not clonally
     dominated (must survive the Mash-lineage clonality correction already built).

CANDIDATE SOURCES TO WEIGH (rank acquirability vs the 4 criteria; not exhaustive):
  - Raw clinical/reference MIC sets with MEASURED (not computational) values: NARMS raw MIC, EUCAST MIC
    distributions, hospital/clinical-lab AST exports, CARD-adjacent measured panels.
  - Curated literature cohorts already flagged label-independent: Whittam DECA pathotype, von Mentzer 2021
    ETEC reference (both need direct-contact/MTA — an acquisition, not a download).
  - Quantitative-assay phenotypes with public genome links: growth/fitness, biochemical/metabolic,
    challenge phenotypes — IF they reach single-organism depth.

DELIVERABLE the downstream chain should produce (NOT this idea-anchor — just framing): a RANKED SHORTLIST
of acquirable label sources, each scored pass/fail on the 4 criteria + an acquisition path (free download
/ public API / MTA / direct-contact / paid) + a rough N-after-all-filters estimate.

CONSTRAINTS: solo hobby project; money spend requires explicit approval (free sources first); two-machine
setup (laptop + Precision 7780); north star is a DECODER TOOL, not papers — a new label only matters if it
makes the tool decode something it currently can't, honestly. Do NOT propose anything that reintroduces
circular or sampling-defined labels; that's the exact trap being escaped.
```

---

## What to expect / how it chains

- `/idea-anchor` stops after framing (planning-pipeline `Each step STOPS`). Likely ≤3 foundational
  questions: (i) validate-existing-decoder vs new-trait — which is the real target; (ii) is paid data
  access on the table; (iii) is direct-contact/MTA acquisition acceptable or download-only.
- Recommended next step it will likely emit: **`/research`** (survey + score candidate sources), NOT
  `/probe` — this is a data-acquisition question, not a repo-code question.
- **Honest flag:** this fork's bottleneck is ACQUISITION, not analysis. The highest-leverage candidates
  (Whittam DECA, von Mentzer) gate on the user making contact / accepting an MTA — a real-world action no
  skill can execute. Expect the chain to land a decision back on the user, not a thing the executor drives
  to completion inline.
