"""Step 1 of the expression_context plan — build the INDEPENDENT Acinetobacter meropenem cohort.

The expression_context signal (ISAba1-upstream-of-blaOXA-51) was falsified IN-SAMPLE on the cached N=30
cohort (data/raw/acinetobacter_meropenem/). Before the ABSTAIN->R override can be promoted (even to opt-in),
it must hold on a DISJOINT cohort. This script fetches up to 15R/15S Acinetobacter meropenem strains from
NCBI Pathogen Detection AST, EXCLUDING the cohort-1 accessions, downloads assemblies, and runs AMRFinder
(cached/restartable) — reusing the existing validator primitives (organism_drug_validate +
independent_cohort_validate). No money; Docker + network only.

HARD KILL GATE (the anchor's binding risk): if fewer than 15S OR 10R have downloadable assemblies free,
prints COHORT_INFEASIBLE and exits 1 — documenting "floor not independently validatable on free data"
rather than faking a cohort. Exit 0 = cohort built (>=10R AND >=15S, disjoint from cohort-1).

Usage: .venv/Scripts/python.exe scripts/build_acinetobacter_indep_cohort.py
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from scripts.independent_cohort_validate import _exclusion_set, select_independent_cohort
from scripts.organism_drug_validate import _run_dir, ensure_run

GROUP = "Acinetobacter"
AMR_ORG = "Acinetobacter_baumannii"
DRUG = "meropenem"
FIRST_SLUG = "acinetobacter_meropenem"          # cohort-1 dir (data/raw/<FIRST_SLUG>/selected.tsv)
MIN_R, MIN_S = 10, 15                            # promotion gate needs n_S>=15; R floor 10


def main() -> int:
    base = Path("data/raw/acinetobacter_meropenem_indep")
    own_runs = base / "amrfinder_runs"
    gcache = base / "refseq"
    reuse_glob = "data/raw/acinetobacter_*/amrfinder_runs"
    exclude = _exclusion_set(FIRST_SLUG)
    print(f"Excluding {len(exclude)} cohort-1 accessions from {FIRST_SLUG}/selected.tsv")

    sel = select_independent_cohort(GROUP, DRUG, exclude, base / "selected.tsv")
    overlap = set(sel) & exclude
    assert not overlap, f"independent cohort overlaps cohort-1: {overlap}"
    nR = sum(sel.values())
    nS = len(sel) - nR
    print(f"Selected (label-level, pre-download): {len(sel)} ({nR}R/{nS}S)")

    # Download + AMRFinder per strain (cached/restartable). Count strains that actually got a run.
    got = {}
    for i, (acc, y) in enumerate(sel.items(), 1):
        if _run_dir(acc, own_runs, reuse_glob):
            got[acc] = y
            continue
        print(f"  [{i}/{len(sel)}] {acc} (label {'R' if y else 'S'}) ...", flush=True)
        mt = ensure_run(acc, own_runs, gcache, AMR_ORG, reuse_glob)
        if mt is not None:
            got[acc] = y

    gR = sum(got.values())
    gS = len(got) - gR
    print(f"\nWith downloadable assembly + AMRFinder run: {len(got)} ({gR}R/{gS}S)")

    if gR < MIN_R or gS < MIN_S:
        print(f"COHORT_INFEASIBLE: need >={MIN_R}R and >={MIN_S}S with free downloadable assemblies; "
              f"got {gR}R/{gS}S. The EXPRESSION floor is NOT independently validatable on free data for "
              f"Acinetobacter x meropenem at this bar.")
        return 1
    print(f"COHORT_OK: {gR}R/{gS}S independent + disjoint from cohort-1. Ready for Step 6 eval.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
