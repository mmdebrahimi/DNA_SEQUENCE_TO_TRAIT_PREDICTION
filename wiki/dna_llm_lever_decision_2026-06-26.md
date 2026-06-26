# DNA-LLM forward-lever decision (2026-06-26)

Acting on the two non-foreclosed levers from the functional-alphabet closeout (cohort expansion; a
distributed-mechanism drug). Read-only feasibility census of the raw BV-BRC tables on disk.

## Headline: BOTH levers are data-FEASIBLE. The wall is download + AMRFinder (Docker, ~hours), NOT GPU.

| Drug | downloadable R/S (assembly + MLST) | shared R+S lineages (ceiling) | current probe |
|---|---|---|---|
| ciprofloxacin | 6124 / 15378 | **317** | 6 lineages / 43 pairs (N=147, p=0.0565) |
| tetracycline | 10770 / 10595 | **685** | — |

There is enormous headroom: 317 / 685 shared R+S lineages vs the current 6. A powered within-lineage probe
needs only ~20 lineages.

## Reconciliation with the strict-MIC "infeasible" census (2026-05-18)

This is NOT a contradiction — it's a label-bar difference, and it matters:
- The **strict-MIC census** used a 4x-safety-margin MIC-VALUE label (paper-grade). That threw away >95% of
  strains → cipro 17R/4S, tet 1R/0S "infeasible". Correct, for a publication-grade cohort.
- This census uses the **relaxed binary phenotype** (BV-BRC's own Resistant/Susceptible call). That is the
  noisier label the project deliberately avoided for PAPER claims — but it is exactly the label the
  **audit framework** (mechanism x MIC merge + SUSPEND gate) was built to handle, and the v0 decoder ethos
  already endorses ("train on categorical labels, propagate the audit verdict"). For a within-lineage
  DIAGNOSTIC (not a paper claim), the relaxed label + audit gate is the right bar.

## Concrete next artifact (built this run)

`wiki/dna_llm_shared_lineage_manifest_tetracycline_2026-06-26.json` — a balanced **tet shared-lineage
cohort: 20 MLST lineages x 120 strains (60R / 60S)**, real downloadable accessions. This is the cohort a
powered tet within-lineage probe would fetch. (`--manifest ciprofloxacin` produces the cipro analog.)

## Decision + recommendation

- **Lever 1 (cipro expansion): FEASIBLE.** A re-selected ~20-lineage cohort easily powers p<0.05.
- **Lever 2 (tet, distributed-mechanism): FEASIBLE + higher-value.** tet is the regime where curated
  determinants are incomplete, so the functional-alphabet probe on tet is the real diagnostic: if the
  curated-determinant alphabet FAILS within-lineage on tet (expected), that is the first concrete evidence
  that a learned representation has *headroom* over the deterministic decoder — the only honest case for
  the DNA-LLM #3 build.

**Recommended:** run the **tet** cohort first. The remaining work is a single ~3-4 hr Docker batch
(fetch 120 genomes via NCBI Datasets + AMRFinder each at ~95s + re-run `functional_alphabet_probe.py` on
the new cohort). It is **code-closable, attended, no money, no GPU** — but it is a real ~4 hr resource
commit + it adopts the relaxed-label cohort bar, so it is a USER greenlight decision, not an
auto-executed one.

## Wall classification

**Code-closable / external-data-acquisition.** Needs: network (NCBI Datasets fetch) + Docker (AMRFinder) —
both present on this host. No paid resource, no GPU. The only judgment call is the relaxed-label bar (above).
