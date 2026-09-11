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
        # G6 MEASURED 2026-09-01 (was absent while no PEAR fitness value had been read). Values from
        # wiki/pear_g6_screen_2026-09-01.json, Figure3.A cefotaxime, 2,114 per-variant effect sizes
        # extracted from the authors' .RData via R. Wild-type baseline rows excluded -- left in, they
        # alone put mode-share at 0.2678 and trip the bar.
        "mode_share": 0.0019,
        "n_distinct_values": 2106.0,
    },
    # what the memo's gate table records, with G6 now closed by measurement
    "expected": {"G1": PASS, "G2": NOT_APPLICABLE, "G3": PASS, "G4": NOT_APPLICABLE,
                 "G5": NOT_APPLICABLE, "G6": PASS, "G7": NOT_APPLICABLE,
                 "G8": NOT_APPLICABLE, "G9": PASS, "G10": PASS},
    "expected_verdict": "CLEARS",
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

# THE THIRD WORKED EXAMPLE, and the one that filled the schema's empty diagonal: an L1 candidate that
# might CLEAR. PEAR is L4/CLEARS, HBV is L1/REJECTED-on-G1, so nothing had ever exercised the L1 gates
# past G1. Every number below is MEASURED from the committed deposit (data/raw/oxford/), not transcribed
# from a memo -- MIC = 2**upper per scripts/oxford_score.py, the paper's log2-dilution encoding.
OXFORD = {
    "candidate": "Oxford E. coli bacteraemia cohort (PRJNA604975 deposit)",
    "memo": "wiki/oxford_gate_screen_2026-09-11.md",
    "intended_layer": L1_AMR_RS,
    "evidence": {
        "label_provenance_evidence": "Oxford clinical-microbiology broth-microdilution MIC, measured in "
                                     "the hospital lab and deposited as per-drug log2 dilution intervals. "
                                     "No genomic tool produced the label; the cohort's AMRFinder run is a "
                                     "separate genotype file that never feeds the phenotype.",
        "label_is_measured": True,
        "label_semantics_evidence": "an MIC dilution reading, not a description of where or why the "
                                    "isolate was collected. Ascertainment is on BACTERAEMIA (a clinical "
                                    "syndrome), which is independent of the resistance phenotype scored.",
        "label_is_assay_reading": True,
        # G4: measured per drug; gentamicin is the most constrained of the three (192 R / 2,681 S).
        "non_ecosystem_min_class_n": 192.0,
        # G5: 4,979 isolates carry an AMRFinder scan in the deposit; 2,897 carry a MIC row.
        "n_fetchable_assemblies": 2897.0,
        # G6 (L1 form): breakpoint censoring. Of 2,873 gentamicin MICs, ZERO fail to resolve R vs S
        # against the CLSI breakpoints; ceftriaxone leaves 6 of 2,874 unresolved.
        "censored_fraction": 0.0021,
        # G2: ONE study. This is the field that exposes the conflation -- see the memo.
        "largest_source_share": 1.0,
        "n_sources": 1.0,
        # G9/G10: the deployed rules already score these drugs on this cohort (gent acc 0.990,
        # 2026-06-15), so the rule is scoreable against the genotype by demonstration.
        "loci_without_recorded_variant_fraction": 0.0,
        "off_panel_variant_fraction": 0.0,
        # G7 MEASURED, and it TRIPS. The deposit carries `guuid` and nothing else -- no submitter,
        # centre or collection field on any record -- so a leakage-clean provenance-disjoint split
        # cannot be built FROM THE DEPOSIT AS SHIPPED. Reported rather than explained away: unlike G2,
        # this is a capability genuinely ABSENT, not a confound structurally impossible. (The fields may
        # be fetchable from ENA for PRJNA604975; that has NOT been done, so this scores the deposit.)
        "provenance_field_populated_fraction": 0.0,
        # G8 deliberately ABSENT: no MLST or Mash distance ships in the deposit, so effective lineage
        # count is not measurable from it. Guessing it to make the screen resolve is what this module
        # refuses to do -- and G7 is decisive on its own, so nothing hinges on it.
    },
    "expected": {"G1": PASS, "G2": NOT_APPLICABLE, "G3": PASS, "G4": PASS, "G5": PASS,
                 "G6": PASS, "G7": TRIP, "G9": PASS, "G10": PASS},
    "expected_verdict": "REJECTED",
}

CANDIDATES = {"pear": PEAR, "hbv": HBV, "oxford": OXFORD}


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
