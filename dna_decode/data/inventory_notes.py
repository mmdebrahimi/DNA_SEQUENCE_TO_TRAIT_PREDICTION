"""What each dataset IS -- the part of the data inventory that cannot be derived.

`scripts/data_inventory.py` derives everything mechanical: what exists, how big, which scripts reference
it, which artifacts report on it. What it CANNOT derive is meaning -- what a dataset is, what it was
already used for, why it is blocked, what it could still support. That lives here.

THIS FILE IS DELIBERATELY SMALL AND DELIBERATELY INCOMPLETE. Curating all 248 datasets would recreate
the failure the inventory exists to prevent: a large hand-written surface that goes stale silently. The
inventory REPORTS its own coverage (`documented` vs `undocumented`), so an absent note is visible rather
than hidden, and the derived fields still answer "has anything used this?" for every dataset regardless.

GUARDED IN BOTH DIRECTIONS:
  * a note whose dataset no longer exists FAILS `--self-check` (the staleness guard);
  * a dataset with no note is listed as UNDOCUMENTED, never silently dropped (the coverage guard).

Add a note when you learn something about a dataset that the next reader would otherwise rediscover.
Write what you VERIFIED, not what you remember -- the whole point is that this is the layer a reader
trusts, so an unverified claim here is worse than no claim.

Fields:
  what          -- one line: what the data physically is
  used_for      -- what has already been done with it (past tense, verified)
  status        -- live / closed-negative / blocked / superseded / unused
  could_support -- what it could still be used for (judgement, explicitly speculative)
  gotcha        -- the trap a first-time user hits (optional)
"""
from __future__ import annotations

NOTES: dict[str, dict] = {
    # ---------------------------------------------------------------- reference / registries
    "external": {
        "what": "External-cohort validation inputs + crosswalks (Oxford, PROBAC, AR Bank, Sci234).",
        "used_for": "The external-cohort re-validation arm; TMP-SMX experimental overlay scoring.",
        "status": "live",
        "could_support": "Any new independent-cohort scoring that must stay namespace-separate from "
                         "the frozen provenance-disjoint arm.",
        "gotcha": "Results here MUST land in wiki/external_validation_* -- writing the provdisjoint "
                  "schema would silently overwrite a frozen cell (the shared-key trap).",
    },
    "amrfinder_runs": {
        "what": "Cached AMRFinderPlus output (main.tsv + mutations.tsv) per assembly.",
        "used_for": "Every deterministic AMR call; the resfinder/pointfinder concordance runs.",
        "status": "live",
        "could_support": "Free re-scoring of any cached genome at zero marginal AMRFinder cost -- "
                         "1,818 genomes are cached and 311 are leakage-disjoint.",
        "gotcha": "The POINT screen is NOT in main.tsv (zero POINT rows across all 1,818 runs); it is "
                  "in a separate mutations.tsv in --mutation_all form including [WILDTYPE] rows.",
    },
    "amrfinder_db": {"what": "AMRFinderPlus reference database.", "used_for": "All AMRFinder runs.",
                     "status": "live", "could_support": "n/a -- infrastructure."},
    "clinvar": {"what": "ClinVar variant-classification release.",
                "used_for": "The clinvar / germline-pathogenicity cell.", "status": "live",
                "could_support": "Any Mendelian-classification work."},
    "hiv_ref": {"what": "Committed HXB2 CDS references (RT / PR / IN / CA).",
                "used_for": "HIV genome-mode calling for all 5 drug classes.", "status": "live",
                "could_support": "Any new HIV gene target.",
                "gotcha": "Reference integrity is self-checked against catalog WT at every DRM position; "
                          "a frame/coordinate error fails loudly rather than silently mis-calling."},
    "virulencefinder_db": {"what": "VirulenceFinder allele DB (E. coli-scoped).",
                           "used_for": "The genome-map virulence overlay + pathotype resolver.",
                           "status": "live",
                           "could_support": "Non-E. coli VF DBs are the named unbuilt extension."},
    "arabidopsis": {"what": "Arabidopsis flowering-time (FT10) phenotype + accession data.",
                    "used_for": "The Path B embedding test at Gate G2.",
                    "status": "closed-negative",
                    "could_support": "Nothing on the embedding axis -- H2 was FALSIFIED under "
                                     "de-confounding (embeddings learned population structure).",
                    "gotcha": "Do NOT reopen on a bigger GPU. A negative de-confounded metric is a "
                              "signal-vs-structure problem, not a compute problem."},
    "horesh_2021": {"what": "Horesh H4 E. coli pathotype cohort (24 genomes).",
                    "used_for": "The ExPEC recall gate (10/12 = 0.833, precision 1.0).",
                    "status": "blocked",
                    "could_support": "Nothing without a new label source -- ExPEC comes from isolation "
                                     "SITE (gate G3, sampling-defined) and cannot be de-confounded by "
                                     "a bigger cohort."},

    # ---------------------------------------------------------------- cohorts / raw
    "cryptic": {"what": "CRyPTIC M. tuberculosis compendium (12,287 isolates, measured BMD MIC).",
                "used_for": "The TB RIF/INH in-distribution baseline via the VARIANTS.parquet adapter.",
                "status": "live",
                "could_support": "Other TB drugs in the WHO catalogue.",
                "gotcha": "A CRyPTIC-scored number is IN-DISTRIBUTION (the WHO catalogue was built "
                          "partly from CRyPTIC) -- never call it independent validation."},
    "tb_goldset": {"what": "EBI AMR-Portal provenance-disjoint TB cohort (39,193 isolates w/ DST).",
                   "used_for": "The first genuinely independent TB number (RIF 0.920/0.955 raw).",
                   "status": "live",
                   "could_support": "Additional TB drugs with DST coverage in the same portal."},
    "tb_indep": {"what": "Per-isolate TB independent-run VCFs + results checkpoint (2,845 isolates).",
                 "used_for": "The lineage-collapsed independent TB headline (RIF 0.444/0.979).",
                 "status": "live",
                 "could_support": "Re-scoring under a revised catalogue with no refetch."},
    "revalidation_2026-07-11": {
        "what": "A 5.3GB revalidation cohort staged 2026-07-11.",
        "used_for": "UNVERIFIED -- referenced by one artifact and no code.",
        "status": "unused",
        "could_support": "Unknown until inspected. Flagged because 5.3GB of staged cohort that nothing "
                         "consumes is exactly the shape of a forgotten asset.",
    },

    # ---------------------------------------------------------------- processed
    "models": {"what": "Trained classifier pickles + their provenance blocks.",
               "used_for": "pipeline predict; the v0/v0.1 decoder path.", "status": "superseded",
               "could_support": "Little -- the embedding-classifier track is a closed negative; the "
                                "shipped product is the deterministic catalog."},
    "embeddings": {"what": "HDF5 foundation-model embedding caches (NT, mostly cipro cohorts).",
                   "used_for": "Every NT-embedding experiment (gate B, N=40, N=147).",
                   "status": "closed-negative",
                   "could_support": "Nothing on the AMR axis -- 0-for-5 de-confounded. Retained as the "
                                    "reproducibility record, not as live substrate.",
                   "gotcha": "populate(skip_existing=True) skips at GENE level, so a crash leaves a "
                             "partial strain that mean-pools into a valid-LOOKING embedding. Every "
                             "consumer must call verify_complete() first."},

    # ---------------------------------------------------------------- caches (heavy, regenerable)
    "refseq": {"what": "NCBI RefSeq genome cache (FASTA + GFF3 + GBK per accession).",
               "used_for": "Nearly every genome-mode call in the repo (57 consumers).",
               "status": "live", "could_support": "Any new organism cell without refetching."},
    "proteingym": {"what": "ProteinGym DMS substitution benchmark + per-assay reference scores.",
                   "used_for": "The forward variant-effect validation; the ESM2 scale finding.",
                   "status": "live",
                   "could_support": "Any new variant-effect scorer -- this is the benchmark that "
                                    "settled 650M > 3B > 15B.",
                   "gotcha": "The published Average_Spearman is a mean of 5 CATEGORY averages, not a "
                             "per-assay mean; never compare a per-assay median against it."},
    "prospective": {"what": "Prospective-lock cohort sweeps + accrual state.",
                    "used_for": "The 2026-08-24 first accrual (E. coli cipro 0.967 / gent 0.532).",
                    "status": "blocked",
                    "could_support": "The project's STRONGEST tier -- leakage-free by construction. "
                                     "Blocked only on post-lock isolates accruing.",
                    "gotcha": "Both E. coli cells are superseded_by_surface_change since the "
                              "2026-08-31 gentamicin v2 lock; their prior numbers describe a RETIRED "
                              "rule and the clock restarted."},
    "essentiality": {"what": "Keio / Fitness Browser conditional-essentiality data.",
                     "used_for": "The FBA conditional-essentiality line (~25 wiki artifacts).",
                     "status": "live",
                     "could_support": "The condition-SWITCH cell is the one genuinely open regime; the "
                                      "bottleneck is MEASURED (11 of 28 conditions overlap PRECISE-1K).",
                     "gotcha": "An `infeasible` deletion solve is GENUINE essentiality, not a solver "
                               "bug -- do not 'fix' the NaN-to-essential coding."},
    "ss2_checkpoint": {"what": "Per-isolate SeqSero2 + salmserovar policy-variant checkpoints (200 isolates).",
                       "used_for": "Every salmserovar scoring run; the 2026-09-09 resolution stratification.",
                       "status": "live",
                       "could_support": "Any re-scoring of the salmserovar cohort with NO blastn and NO "
                                        "Docker -- the variants differ only in the name lookup."},
    "ecoli_sero_asm": {"what": "E. coli serotype cohort assemblies + per-isolate results.jsonl.",
                       "used_for": "The O:H validation and the identity-primary selection fix.",
                       "status": "live",
                       "could_support": "Re-scoring offline. NOT sufficient for a resolution "
                                        "stratification: the rows carry only final O_call/H_call, no "
                                        "allele candidates and no second rule variant."},
    "ktype_kaptive": {"what": "Per-isolate ktype-vs-Kaptive concordance rows (307 genomes).",
                      "used_for": "The 2026-09-08 near-ceiling verdict (0.8788 strict).",
                      "status": "live",
                      "could_support": "Little that is new -- the wzi_50 multi-KL composition is "
                                       "ALREADY recorded in the artifact; re-deriving it is packaging."},
    "caur_sra": {"what": "13.2GB of C. auris SRA reads (reads/ + asm/ + assemblies/).",
                 "used_for": "UNVERIFIED -- nothing in the tree references this directory by name.",
                 "status": "unused",
                 "could_support": "Probably the raw inputs behind the fungal G1 cohort, which was "
                                  "built by targeted ERG11 read-mapping. Verify before deleting AND "
                                  "before re-downloading -- it is the single largest orphan on disk.",
                 "gotcha": "13.2GB on a disk-tight host. Its status is a genuine open question, which "
                           "is exactly why it is recorded rather than assumed dead."},
    "nt-gate-b": {"what": "3.8GB NT embedding cache for the gate-B cohort.",
                  "used_for": "UNVERIFIED by name -- scripts reference nt_gate_b_cohort_67.h5 "
                              "(underscores) rather than this directory (hyphens).",
                  "status": "closed-negative",
                  "could_support": "Nothing live -- the NT-embedding track is a closed negative. "
                                   "Strong reclaim candidate if disk is needed."},
    "data files donwload": {
        "what": "3.4GB directory, name misspelled (`donwload`), contents unclassified.",
        "used_for": "Nothing references it.",
        "status": "unused",
        "could_support": "Unknown. Recorded because an unexamined multi-GB directory is precisely "
                         "what gets re-downloaded later.",
    },
}
