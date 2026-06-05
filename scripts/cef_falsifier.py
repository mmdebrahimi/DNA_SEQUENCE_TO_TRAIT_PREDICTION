"""Cef falsifier — thin shim over the drug-agnostic scripts/amr_falsifier.py.

Superseded by amr_falsifier.py (2026-06-04) which adds the cohort de-confound gate + CI-aware
verdict. Kept as a convenience entry with cef defaults; the de-confound gate now BLOCKS the cef
gate_b_cohort (CONFOUNDED: R≈USA, S≈Africa/India, 1 shared lineage) — see
plans/cef_falsifier_brainstorm.md. Pass extra args through to override.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from scripts.amr_falsifier import main

CEF_DEFAULTS = [
    "--drug", "ceftriaxone",
    "--cohort", "data/processed/gate_b_cohort.parquet",
    "--nt-cache", "data/processed/embeddings/nt_gate_b_cohort_67.h5",
]

if __name__ == "__main__":
    argv = sys.argv[1:] if len(sys.argv) > 1 else CEF_DEFAULTS
    sys.exit(main(argv))
