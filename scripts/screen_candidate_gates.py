"""Run the ten rejection gates over a candidate, and CHECK the screen against the hand-made verdicts.

Two candidates have been screened by hand and committed as memos. They are encoded here as evidence
packets plus the per-gate verdict each memo actually recorded, so `--verify` re-derives them mechanically.
If the screen and a memo disagree, that is a finding to diagnose -- one of them is wrong -- not a number
to report as passing.

    uv run python scripts/screen_candidate_gates.py --verify        # reproduce both hand verdicts
    uv run python scripts/screen_candidate_gates.py --candidate pear
    uv run python scripts/screen_candidate_gates.py --packet my_candidate.json

Offline, read-only, seconds. Writes nothing unless --out is given.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from dna_decode.eval.rejection_gates import (  # noqa: E402
    INSUFFICIENT_DATA, L1_AMR_RS, L4_FORWARD_CONTINUOUS, NOT_APPLICABLE, PASS, TRIP,
    screen_candidate,
)

# --- the two hand-worked candidates ---------------------------------------------------------------
# Evidence is transcribed from the committed memos. Where a memo did not measure something, the field is
# ABSENT -- never guessed to make the screen come out clean.

PEAR = {
    "candidate": "PEAR (Zhang 2022 blaCTX-M-14 DMS)",
    "memo": "wiki/pear_substrate_screen_2026-08-31.md",
    "intended_layer": L4_FORWARD_CONTINUOUS,
    "evidence": {
        "label_provenance_evidence": "relative growth of ~23,000 constructed strains under cefotaxime/"
                                     "ceftazidime selection, read out by barcode sequencing. No genomic "
                                     "tool produced the label.",
        "label_is_measured": True,
        "label_semantics_evidence": "relative growth is an assay reading, not a collection context",
        "label_is_assay_reading": True,
        "variation_is_constructed": True,
        "genotype_defined_by_construction": True,
        "loci_without_recorded_variant_fraction": 0.0,   # every variant known by construction
        "off_panel_variant_fraction": 0.0,               # single substitutions in one gene
        # mode_share / n_distinct_values DELIBERATELY ABSENT: the memo's own honest limit is
        # "I have not read one PEAR fitness value. G6 is unscreened."
    },
    # what the memo's gate table actually recorded
    "expected": {"G1": PASS, "G2": NOT_APPLICABLE, "G3": PASS, "G4": NOT_APPLICABLE,
                 "G5": NOT_APPLICABLE, "G6": INSUFFICIENT_DATA, "G7": NOT_APPLICABLE,
                 "G8": NOT_APPLICABLE, "G9": PASS, "G10": PASS},
    "expected_verdict": "INCOMPLETE",
}

HBV = {
    "candidate": "HBV RT (nucleos(t)ide-analogue resistance)",
    "memo": "wiki/hbv_cell_gate_screen_2026-09-01.md",
    "intended_layer": L1_AMR_RS,
    "evidence": {
        "label_provenance_evidence": "every free HBV resource is an interpretation SYSTEM "
                                     "(geno2pheno[HBV], HIV-GRADE) or a mutation-PREVALENCE annotation "
                                     "(Stanford HBVseq/HBVrtDB, dormant since 2012). No free compiled "
                                     "measured-phenotype source found. Cause is stated in the field: HBV "
                                     "has no simple cell-culture system, so phenotypic testing is scarce.",
        "label_is_measured": False,
        # Nothing else was measured, and nothing else needed to be: G1 is decisive on its own.
    },
    "expected": {"G1": TRIP},
    "expected_verdict": "REJECTED",
}

CANDIDATES = {"pear": PEAR, "hbv": HBV}


def run_one(spec: dict) -> dict:
    res = screen_candidate(spec["candidate"], spec["intended_layer"], spec["evidence"])
    out = res.as_dict()
    out["memo"] = spec.get("memo")
    return out


def check_one(spec: dict) -> tuple[bool, list[str]]:
    """Does the screen reproduce what the memo recorded? Returns (ok, discrepancies)."""
    res = screen_candidate(spec["candidate"], spec["intended_layer"], spec["evidence"])
    got = {g.gate: g.verdict for g in res.gates}
    bad = []
    for gate, want in spec["expected"].items():
        if got.get(gate) != want:
            bad.append(f"{gate}: memo says {want!r}, screen says {got.get(gate)!r}")
    if res.verdict != spec["expected_verdict"]:
        bad.append(f"overall: memo implies {spec['expected_verdict']!r}, screen says {res.verdict!r}")
    return (not bad), bad


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--verify", action="store_true",
                    help="re-derive both committed hand verdicts and report any disagreement")
    ap.add_argument("--candidate", choices=sorted(CANDIDATES), help="screen one committed candidate")
    ap.add_argument("--packet", type=Path,
                    help="screen a candidate from a JSON packet {candidate, intended_layer, evidence}")
    ap.add_argument("--out", type=Path, help="write the screen result as JSON")
    a = ap.parse_args()

    if a.packet:
        spec = json.loads(a.packet.read_text(encoding="utf-8"))
        res = run_one(spec)
    elif a.candidate:
        res = run_one(CANDIDATES[a.candidate])
    elif a.verify:
        res = None
    else:
        ap.error("pass --verify, --candidate, or --packet")

    if a.verify:
        print("Reproducing the committed hand verdicts\n" + "=" * 60)
        failed = 0
        for key, spec in CANDIDATES.items():
            ok, bad = check_one(spec)
            print(f"\n{key}: {'REPRODUCED' if ok else 'DISAGREES'}  ({spec['memo']})")
            for line in bad:
                print(f"    ! {line}")
            failed += 0 if ok else 1
        print("\n" + "=" * 60)
        if failed:
            print(f"{failed} candidate(s) disagree with their memo. One of the two is wrong -- "
                  "diagnose before trusting either.")
            return 1
        print("Both hand verdicts re-derived mechanically.")
        return 0

    print(json.dumps(res, indent=2))
    if a.out:
        a.out.write_text(json.dumps(res, indent=2) + "\n", encoding="utf-8")
        print(f"\nwrote {a.out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
