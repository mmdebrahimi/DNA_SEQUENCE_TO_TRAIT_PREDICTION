"""Genome Map v1 — the "Bakta honesty report".

A single-genome, evidence-tiered function/QC map: it RE-TIERS Bakta's own
db-light annotation for honesty (4 tiers) + overlays the existing AMRFinder
curated determinant cells behind a HARD join-quality gate, and reports a
DB-labelled unknown rate. This is the achievable, label-free form of the
dna_decode north star ("DNA -> what its parts do") — NOT a learned
genotype->phenotype predictor (that arm is a closed negative on free data;
see wiki/north_star_distance_brainstorm_2026-06-17.md).

v1 is a SPIKE: a 3-bacterial-genome prototype that ends in a GO/NO-GO verdict
on whether to invest further. Bacterial-AMR-only (TB's VCF-vs-GFF contract is
out of the v1 spike).

Package layout:
  ingest.py            — the shared ##FASTA-stripping GFF loader (Bakta + offline both go through it)
  annotate.py          — thin Bakta (db-light) Docker runner
  amrfinder.py         — thin AMRFinderPlus Docker runner (explicit organism)
  tiers.py / tier_vocab.py  — the 4-tier classifier (vocab seeded from the Step-2 manifest)
  phenotype_overlay.py — DeterminantHit + the hard join-quality gate
  build_map.py         — the map assembler (raw-field schema + phenotype wall + unknown rate)
  gate.py              — the G1/G2 prevent-wrong-inference gate + GO/NO-GO verdict

The phenotype tier READS the frozen AMR surface (dna_decode/eval/amr_rules.py +
dna_decode/data/calibrated_amr_rules.json); it NEVER modifies it.
"""
from __future__ import annotations

# The four v1 evidence tiers, in precedence order (highest -> lowest). Only the
# first may carry a phenotype claim (the phenotype wall). `pathway-module` is
# DEFERRED (no KEGG in v1).
TIER_DETERMINANT_PHENOTYPE = "determinant-phenotype"
TIER_CURATED_FUNCTION = "curated-molecular-function"
TIER_HOMOLOGY_HYPOTHESIS = "homology-only-hypothesis"
TIER_UNKNOWN = "unknown"

TIER_PRECEDENCE = (
    TIER_DETERMINANT_PHENOTYPE,
    TIER_CURATED_FUNCTION,
    TIER_HOMOLOGY_HYPOTHESIS,
    TIER_UNKNOWN,
)

# The single tier that may emit a phenotype/property claim (R1, the phenotype wall).
PHENOTYPE_TIER = TIER_DETERMINANT_PHENOTYPE

__all__ = [
    "TIER_DETERMINANT_PHENOTYPE",
    "TIER_CURATED_FUNCTION",
    "TIER_HOMOLOGY_HYPOTHESIS",
    "TIER_UNKNOWN",
    "TIER_PRECEDENCE",
    "PHENOTYPE_TIER",
]
