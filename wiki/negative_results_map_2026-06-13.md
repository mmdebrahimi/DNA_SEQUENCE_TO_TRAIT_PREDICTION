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

## The 10 rejection gates (screen every candidate against these)

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
| G9 | **Causal variant unrecorded** (added 2026-08-26) | the decoder's rule names allele SYMBOLS and dominance order but never the causal VARIANT, so no genotype file can be scored against it — unvalidatable as written, and no cohort fixes it | for each locus, does the catalog record a concrete variant (`c.` / `p.` / a named indel)? Count unrecorded loci; a majority = trip. Measured: **40 of 65 colour-cell loci (62%) record none**, and 7 of 19 cells record none for ANY locus (`wiki/colour_cell_substrate_screen_2026-08-26.md`) |
| G10 | **Variant class off-panel** (added 2026-08-26) | the causal variants ARE recorded but are indels / structural / CNV, which a SNP array or imputed biallelic-SNV panel cannot represent — so the commonest free genotype cohorts cannot score the rule however good it is | classify each recorded causal variant as SNV vs indel/structural; count the off-panel share. Measured: **14 of 25 recorded colour-cell variants (56%) are indel/structural**, and this is what sank the dog cell empirically — black 0.994 but every other base colour unscorable (`wiki/dog_coat_darwins_ark_measured_2026-07-30.md`) |

## The verified failure record (each row = a closed track + the gates it tripped)

| Track | Verdict | Gates tripped | Evidence |
|---|---|---|---|
| Pathotype (EnteroBase / NCBI-PD labels) | label-blocked | G1 (BlastFrost/AMRFinder-derived), G3 (isolation-site), G2 (study==class on the 24-genome ExPEC/EPEC) | `research_outputs/horesh-f1-label-provenance-audit-2026-06-04.md` (H1 falsified; curated-independent fraction 20.5%) |
| Foundation-model embeddings (0-for-4) | embedding learns lineage not mechanism | (not a label gate — a MODEL ceiling) | `wiki/embedding_niche_cross_domain_synthesis_2026-06-12.md` (cipro within-lineage=chance; Arabidopsis within-group r2 −0.13) |
| MIC-continuous (graded resistance) | not-feasible | G1 (91% of BV-BRC MIC is XGBoost-from-genome), G6 (~70% breakpoint-censored), G8 (within-R exact N=6) | ledger action 91 (`project_state/dna-decode-2026-05-11.md`); probe 2026-06-13 |
| AMR grid — Salmonella tet/gent | underpowered, infeasible | G4 (ecosystem-dominated: tet 4871R ecosystem → 4R disjoint) | ledger action 79; `wiki/provdisjoint_census_results.json` |
| AMR grid — Acinetobacter/Pseudomonas/Klebsiella-class broad expansion | intrinsic-gene degeneracy | (mechanism ceiling — intrinsic class-genes over-call, spec→0) | `~/.claude/...memory/feedback_intrinsic_genes_break_broad_amr_class_rules.md` |
| Cross-organism catalog transfer (g8r2 — "zero target labels" bypass) | label-wall NOT bypassed + already-built | (not a new label gate — the naive unchanged transfer FAILS per-organism on threshold/content, so target labels are still needed to calibrate) | `wiki/wider_amr_transferability_synthesis_2026-06-08.md` (Campylobacter/Salmonella cipro FAIL: TUNING/CONTENT); `wiki/innovate_round2_frontier_closeout_2026-07-13.md` |
| Self-supervised catalog (g5r2 — train an effect-predictor on the catalog itself) | no negative class / circular where one exists | (not a new label gate — 3/4 cells are all-positive DRM sets; WHO-TB grades ARE distilled phenotype → self-prediction is circular) | `wiki/innovate_round2_frontier_closeout_2026-07-13.md` (HIV/SARS/fungal all-positive; WHO 457 R-assoc / 550 benign / 47,139 uncertain) |

## What SURVIVED (the shipped product, and why it cleared the gates)

The deterministic AMR decoder on **AMR MIC R/S (broth-microdilution)** — the single label that cleared all
8 gates at free-public scale: sampling-independent lab measurement (not G1/G3), provenance-separable
(not G7), and powered + lineage-robust on the acquired-gene mechanisms (β-lactamase, tet efflux,
aminoglycoside-modifying) after clonality disclosure (survives G4/G8 on 10 organism×drug cells). The cipro
chromosomal-QRDR cells are clonality-inflated at the isolate level (disclosed, not hidden — see the report
card's lineage table). See `wiki/decoder_validation_report_card.md`.

## How to use this map

Before proposing any new trait/organism/label source, screen it against G1–G10. A candidate that trips any
gate is not a viable honest-decoder substrate on available public data — do NOT spend acquisition or
modeling labor on it. The only ways forward that the map does NOT foreclose:
1. A **non-public** label source (clinical/biobank/collaborator wet-lab measurements) that clears G1/G3/G7
   by construction — an ACQUISITION decision, gated on a concrete named source in hand.
2. **Prospective-lock validation** of the existing decoder against later-arriving independent data (needs
   no new label today; see `wiki/reproducibility_freeze_2026-06-13.md`).

### G1–G8 gate the LABEL; G9–G10 gate the DECODER (added 2026-08-26)

G9/G10 (added 2026-08-26) are the DECODER-SIDE gates: G1–G8 all ask whether a usable LABEL exists, while G9/G10 ask whether the decoder's own rule is scoreable against a genotype at all. Screen a new curated-catalog cell against them BEFORE building it — the animal-colour family reached 19 cells before anyone checked, and 7 of them cannot be validated on any substrate as written.

### The map runs in BOTH directions (added 2026-08-23)

Every use of G1–G8 above was a rejection. The first worked example of a substrate that **passes** is
recorded at `wiki/moderna_int_landscape_and_gate_screen_2026-08-23.md` — MHC-I peptide presentation, the
label layer under the Merck/Moderna neoantigen therapy. It passes **G1 / G3 / G5** cleanly, which are
exactly the gates that closed pathotype and the BV-BRC cohorts, and its residual risk sits on the
*negative* class (G2/G6: MS non-detection ≠ non-presentation).

Two things that makes usable:

- **A passing screen is not a build recommendation.** That substrate was screened and then *declined* —
  wrong layer for the north star (peptide→presentation, not DNA→trait), a mature field with strong
  incumbents on the identical public corpora, and no data advantage. Record the screen; decide separately.
- **It calibrates what "clears the gates" actually looks like**, so a future candidate can be compared
  against a real positive rather than only against the failure record.
