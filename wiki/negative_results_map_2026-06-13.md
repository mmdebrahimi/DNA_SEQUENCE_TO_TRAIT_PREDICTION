# Operational negative-results map — why each G2P expansion was rejected (2026-06-13)

**Purpose:** This is NOT a prose graveyard of failed attempts. It is an *operational* map: each rejected
expansion is encoded as a reusable GATE so future work (future-you, or anyone forking this tool) can
screen a candidate dataset/trait against the failure modes BEFORE spending labor. The recurring lesson —
verified across every closed track this project ran — is that **honest public bacterial genotype→phenotype
decoding is bounded by LABELS, not models.** A candidate that trips any gate below cannot support an honest
learned decoder on available public data.

This map is the scientific contribution that accompanies the shipped deterministic AMR decoder: it states
the boundary of what public G2P can honestly support, so the tool is never overextended into a dishonest
regime.

## The 8 rejection gates (screen every candidate against these)

| # | Gate | Trips when | How to check (cheap) |
|---|---|---|---|
| G1 | **Circular label** | the phenotype label is itself produced by a genomic tool the decoder would compete against (AMRFinder / BlastFrost / any genome→label model) | inspect label provenance field / dataset methods — is the label wet-lab/clinical or gene-call-derived? |
| G2 | **Study == class** | the label is confounded with the source study / submitter (one BioProject supplies most of one class) | contingency table of class × BioProject/submitter; a dominant cell = trip |
| G3 | **Sampling-defined label** | the label IS the sampling context, not a measurement (e.g. blood-vs-feces = pathotype by definition) | ask: is the label an assay reading, or a description of where/why the isolate was collected? |
| G4 | **Surveillance domination** | excluding the surveillance ecosystem (NARMS/CDC/FDA/GenomeTrakr/PulseNet/USDA) collapses the resistant/positive pool below a usable cohort | census the OTHER (non-ecosystem) class counts; <20/class = trip |
| G5 | **Assembly attrition** | label-bearing records exist but lack downloadable assemblies (the historical ~96% drop) | count records WITH `assembly_accession` that NCBI Datasets can fetch, not raw record count |
| G6 | **Phenotype censoring** | the quantitative label is interval-censored exactly where it matters (MIC `>X` / `<=X` at the breakpoint) | tally exact vs censored values among the in-class subset; majority-censored = trip |
| G7 | **Provenance not separable** | metadata is too thin to build a leakage-clean provenance-disjoint split | are submitter / center / collection fields populated per-record? absent = trip |
| G8 | **Dedup collapses balance** | after Mash-lineage clonality correction, one class drops below usable effective-N (clonally dominated) | greedy-representative Mash cluster per class; <~3 effective lineages = trip |

## The verified failure record (each row = a closed track + the gates it tripped)

| Track | Verdict | Gates tripped | Evidence |
|---|---|---|---|
| Pathotype (EnteroBase / NCBI-PD labels) | label-blocked | G1 (BlastFrost/AMRFinder-derived), G3 (isolation-site), G2 (study==class on the 24-genome ExPEC/EPEC) | `research_outputs/horesh-f1-label-provenance-audit-2026-06-04.md` (H1 falsified; curated-independent fraction 20.5%) |
| Foundation-model embeddings (0-for-4) | embedding learns lineage not mechanism | (not a label gate — a MODEL ceiling) | `wiki/embedding_niche_cross_domain_synthesis_2026-06-12.md` (cipro within-lineage=chance; Arabidopsis within-group r2 −0.13) |
| MIC-continuous (graded resistance) | not-feasible | G1 (91% of BV-BRC MIC is XGBoost-from-genome), G6 (~70% breakpoint-censored), G8 (within-R exact N=6) | ledger action 91 (`project_state/dna-decode-2026-05-11.md`); probe 2026-06-13 |
| AMR grid — Salmonella tet/gent | underpowered, infeasible | G4 (ecosystem-dominated: tet 4871R ecosystem → 4R disjoint) | ledger action 79; `wiki/provdisjoint_census_results.json` |
| AMR grid — Acinetobacter/Pseudomonas/Klebsiella-class broad expansion | intrinsic-gene degeneracy | (mechanism ceiling — intrinsic class-genes over-call, spec→0) | `~/.claude/...memory/feedback_intrinsic_genes_break_broad_amr_class_rules.md` |

## What SURVIVED (the shipped product, and why it cleared the gates)

The deterministic AMR decoder on **AMR MIC R/S (broth-microdilution)** — the single label that cleared all
8 gates at free-public scale: sampling-independent lab measurement (not G1/G3), provenance-separable
(not G7), and powered + lineage-robust on the acquired-gene mechanisms (β-lactamase, tet efflux,
aminoglycoside-modifying) after clonality disclosure (survives G4/G8 on 10 organism×drug cells). The cipro
chromosomal-QRDR cells are clonality-inflated at the isolate level (disclosed, not hidden — see the report
card's lineage table). See `wiki/decoder_validation_report_card.md`.

## How to use this map

Before proposing any new trait/organism/label source, screen it against G1–G8. A candidate that trips any
gate is not a viable honest-decoder substrate on available public data — do NOT spend acquisition or
modeling labor on it. The only ways forward that the map does NOT foreclose:
1. A **non-public** label source (clinical/biobank/collaborator wet-lab measurements) that clears G1/G3/G7
   by construction — an ACQUISITION decision, gated on a concrete named source in hand.
2. **Prospective-lock validation** of the existing decoder against later-arriving independent data (needs
   no new label today; see `wiki/reproducibility_freeze_2026-06-13.md`).
