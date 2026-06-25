"""Validate the S. pneumoniae gene-presence AMR rule vs WET-LAB measured AST (GPS Poland cohort).

Scope = the clean gene-presence drugs (macrolide erm/mef, tetracycline tet). β-lactams are deferred (PBP
engine; see the go/no-go memo). The measured label is the GPS pipeline-paper Supplementary Data 1 disc/agar
AST (CLSI S. pneumoniae disc-diffusion zone breakpoints); the determinant genotype is the GPS pipeline's own
determinant calls (Supplementary Data 2).

HONESTY TIER: the WET-LAB measured AST is independent of any caller (clears the circularity gate). The
determinant CALLS here are GPS's (gene-presence BLAST) — for erm/mef/tet that is ~AMRFinder-equivalent, so
this is a NEAR-INDEPENDENT measured-label number; a truly independent run swaps in OUR AMRFinder (organism
Streptococcus_pneumoniae) on the assemblies (Docker; deferred). Emitted with `genotype_source` so the tier
is never overstated.

    uv run python scripts/pneumo_amr_validate.py   # reads the GPS CSVs from D:/dna_decode_cache/pneumo_gps
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))
from dna_decode.organism_rules.pneumo_amr import call_drug  # noqa: E402

GPS = Path("D:/dna_decode_cache/pneumo_gps")

# CLSI S. pneumoniae disc-diffusion zone-diameter breakpoints (mm): higher zone = more susceptible.
DISC_BREAKPOINTS = {
    "erythromycin": {"col": "Erythromycin", "S": 21, "R": 15},   # 15 µg disc: S>=21, R<=15
    "tetracycline": {"col": "Tetracycline", "S": 23, "R": 18},   # 30 µg disc: S>=23, R<=18
}
# determinant column in Supplementary Data 2 (the GPS pipeline genotype)
DET_COL = {"erythromycin": "ERY_Determinant", "tetracycline": "TET_Determinant"}


def _zone(v):
    try:
        return float((v or "").strip())
    except ValueError:
        return None


def _measured_rs(drug: str, zone: float | None) -> str | None:
    if zone is None:
        return None
    bp = DISC_BREAKPOINTS[drug]
    return "S" if zone >= bp["S"] else "R" if zone <= bp["R"] else "I"


def main(argv=None) -> int:
    m3 = {r["ERR"].strip(): r for r in csv.DictReader(open(GPS / "sd_3.csv", encoding="utf-8")) if r.get("ERR")}
    m4 = {r["Sample_ID"].strip(): r for r in csv.DictReader(open(GPS / "sd_4.csv", encoding="utf-8")) if r.get("Sample_ID")}
    both = set(m3) & set(m4)

    out = {"schema": "pneumo-amr-genepresence-validation-v1",
           "label": "WET-LAB measured AST (GPS Poland; CLSI disc-diffusion zone breakpoints)",
           "genotype_source": "GPS pipeline determinant calls (gene-presence BLAST; ~AMRFinder-equivalent for erm/mef/tet)",
           "honesty_tier": "NEAR_INDEPENDENT_MEASURED (label independent; determinant-calling = GPS, our-AMRFinder swap deferred)",
           "scope": "gene-presence drugs only (macrolide, tetracycline); beta-lactam PBP cell deferred",
           "drugs": {}}
    for drug in ("erythromycin", "tetracycline"):
        bp = DISC_BREAKPOINTS[drug]
        tp = tn = fp = fn = skip = 0
        for k in both:
            meas = _measured_rs(drug, _zone(m3[k].get(bp["col"])))
            if meas in (None, "I"):
                skip += 1
                continue
            call = call_drug(drug, [m4[k].get(DET_COL[drug], "")])
            pred = call.prediction if call else "S"
            if meas == "R" and pred == "R":
                tp += 1
            elif meas == "S" and pred == "S":
                tn += 1
            elif meas == "S" and pred == "R":
                fp += 1
            else:
                fn += 1
        n = tp + tn + fp + fn
        out["drugs"][drug] = {
            "n": n, "skipped_I_or_missing": skip,
            "accuracy": round((tp + tn) / n, 3) if n else None,
            "sensitivity": round(tp / (tp + fn), 3) if (tp + fn) else None,
            "specificity": round(tn / (tn + fp), 3) if (tn + fp) else None,
            "tp": tp, "fp": fp, "tn": tn, "fn": fn,
        }
        print(f"{drug}: n={n} acc={out['drugs'][drug]['accuracy']} "
              f"sens={out['drugs'][drug]['sensitivity']} spec={out['drugs'][drug]['specificity']} "
              f"TP{tp} FP{fp} TN{tn} FN{fn}")
    (REPO / "wiki" / "pneumo_amr_genepresence_validation.json").write_text(json.dumps(out, indent=2), encoding="utf-8")
    print("-> wiki/pneumo_amr_genepresence_validation.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
