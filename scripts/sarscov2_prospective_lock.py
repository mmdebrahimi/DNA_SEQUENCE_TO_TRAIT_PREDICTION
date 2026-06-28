"""SARS-CoV-2 Mpro prospective-lock — the free, leakage-free-by-construction path to an INDEPENDENT
nirmatrelvir number (the only viable SARS v0.1 lever; acquisition is blocked, see
wiki/sarscov2_mpro_v0.1_independence_infeasible_2026-06-27.md).

Mirrors the AMR prospective-lock (`dna_decode/eval/prospective_lock.py`) for the SARS Mpro cell, REUSING its
generic primitives (`_sha256_file`, `is_prospective_eligible`) — touches NO AMR-surface file and does NOT
modify `prospective_lock.py`. Pins the frozen SARS Mpro decoder surface NOW + stamps a lock DATE; any future
CoV-RDB (or external) Mpro-inhibitor fold record whose study/publication date is STRICTLY AFTER the lock is a
provably-leakage-free independent test case (the frozen catalog could not have been tuned to data that didn't
exist). The census (2026-06-27) showed 0 held-out studies TODAY, so this ACCRUES — it produces the mechanism,
not a number yet. `verify_lock` re-hashes the live surface so a drifted catalog is tamper-evident.
"""
from __future__ import annotations

import json
import sys
from datetime import date as _date
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from dna_decode.eval.prospective_lock import _sha256_file, is_prospective_eligible  # noqa: E402

# The frozen SARS Mpro decoder surface (catalog + genome-mode caller). NOT the AMR FROZEN_SURFACE_FILES.
SARS_MPRO_SURFACE_FILES: tuple[str, ...] = (
    "dna_decode/data/sarscov2_amr.py",     # MPRO_MAJOR_DRMS catalog + call (the decision surface)
    "scripts/sarscov2_caller.py",          # genome-mode BLAST->codon-map caller
)
SCHEMA = "sarscov2-mpro-prospective-lock-v1"


def surface_hashes(repo: Path = REPO) -> dict[str, str]:
    return {rel: _sha256_file(rel, repo) for rel in SARS_MPRO_SURFACE_FILES}


def compute_manifest(lock_date: str, repo: Path = REPO) -> dict:
    return {
        "schema": SCHEMA,
        "cell": "sarscov2_mpro_nirmatrelvir",
        "lock_date": lock_date,
        "manifest_created": lock_date,
        "eligibility_rule": ("a Mpro-inhibitor fold record is an eligible INDEPENDENT test iff its study's "
                             "earliest public/publication date is STRICTLY AFTER lock_date (leakage-free by "
                             "construction — the frozen catalog could not have been tuned to it)"),
        "cutoff_justification": ("the SARS Mpro catalog (sarscov2_amr.MPRO_MAJOR_DRMS) shipped 2026-06-23 "
                                 "verbatim from CoV-RDB invitro_selection_results; lock_date is set on/after "
                                 "that freeze so the pinned surface is provably unchanged since"),
        "surface_sha256": surface_hashes(repo),
        "accrual_status": ("ACCRUING — the 2026-06-27 census found 0 held-out CoV-RDB studies; all nirmatrelvir/"
                           "ensitrelvir fold comes from the 2 catalog-building studies (FDA23-NTV, Krismer23). "
                           "This manifest makes any LATER study a provable independent test; it is the mechanism, "
                           "not a number"),
        "scoring_path": ("re-run scripts/sarscov2_mpro_validate.py, filter rx_fold to refs whose study date > "
                         "lock_date (is_prospective_eligible), score the frozen catalog on that held-out set"),
        "what_this_is_not": ("NOT an independent number today (cohort empty); NOT a frozen-surface change "
                             "(read-only pin); NOT the AMR prospective-lock (separate cell, separate manifest)"),
    }


def verify_lock(manifest: dict, repo: Path = REPO) -> tuple[bool, list[str]]:
    """Re-hash the live SARS Mpro surface vs the manifest -> (ok, drifted_files). Tamper-evident."""
    pinned = manifest.get("surface_sha256", {})
    drifted = [rel for rel, want in pinned.items() if _sha256_file(rel, repo) != want]
    return (not drifted, drifted)


def main(argv=None) -> int:
    import argparse
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--lock-date", default=_date.today().isoformat(), help="eligibility cutoff YYYY-MM-DD")
    ap.add_argument("--out", type=Path, default=None)
    ap.add_argument("--verify", type=Path, default=None, help="verify an existing manifest instead of writing")
    a = ap.parse_args(argv)
    if a.verify:
        man = json.loads(a.verify.read_text(encoding="utf-8"))
        ok, drifted = verify_lock(man)
        print("LOCK INTACT" if ok else f"LOCK DRIFTED: {drifted}")
        return 0 if ok else 2
    man = compute_manifest(a.lock_date)
    out = a.out or (REPO / "wiki" / f"sarscov2_prospective_lock_manifest_{_date.today().isoformat()}.json")
    out.write_text(json.dumps(man, indent=2), encoding="utf-8")
    print(f"SARS Mpro prospective-lock pinned (lock_date={man['lock_date']}); surface:")
    for rel, h in man["surface_sha256"].items():
        print(f"  {rel}  {h[:16]}…")
    print(f"-> {out}  (ACCRUING; 0 held-out studies today)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
