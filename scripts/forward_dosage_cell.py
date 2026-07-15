"""Forward-cell dosage head — real validation on PTEN: calibrate the AlphaMissense score to the measured
DMS magnitude + split-conformal prediction intervals, report held-out coverage + interval narrowing.

The forward cell's methods rank variants; this turns the score into a CALIBRATED MAGNITUDE ("this edit's
effect is X +/- q on the DMS scale") with a pre-registered coverage target. Uses AlphaMissense on PTEN
(instant, no GPU) as the score; the calibration is method-agnostic (any forward-cell score works).

PRE-REGISTERED BAR: CALIBRATED_DOSAGE iff |mean held-out coverage - TARGET| <= COVER_TOL (0.05) AND the
score meaningfully narrows the interval vs the predict-the-mean marginal (interval_narrowing > NARROW_MIN).
The narrowing gate is the load-bearing honesty check (conformal coverage holds even for a useless model).
"""
from __future__ import annotations

import csv
import json
import sys
from datetime import date as _date
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from dna_decode.forward import am_table_for_mutants, load_am_for_uniprot  # noqa: E402
from dna_decode.forward.dosage import conformal_q, evaluate_dosage  # noqa: E402

PG = Path("D:/dna_decode_cache/proteingym")
ASSAY = "PTEN_HUMAN_Mighell_2018"
UNIPROT = "P60484"
TARGET = 0.80
COVER_TOL = 0.05
NARROW_MIN = 0.02
REPEATS = 20
SEED = 7


def load_xy():
    """x = AlphaMissense benign score (1 - pathogenicity), y = measured DMS_score, over shared variants."""
    dms = {}
    mutants = []
    for r in csv.DictReader(open(PG / "pg_dms" / "DMS_ProteinGym_substitutions" / f"{ASSAY}.csv",
                                 encoding="utf-8")):
        m = r["mutant"].strip()
        if ":" in m:
            continue
        mutants.append(m)
        try:
            dms[m] = float(r["DMS_score"])
        except (TypeError, ValueError):
            pass
    am = am_table_for_mutants(load_am_for_uniprot(PG / "am_filtered.tsv", UNIPROT), 0, mutants)
    xs, ys = [], []
    for m, a in am.items():
        if m in dms:
            xs.append(1.0 - a)     # higher = benign (aligns polarity with DMS function score)
            ys.append(dms[m])
    return xs, ys


def main() -> int:
    import numpy as np
    # --- non-duplication check: my conformal_q == J2's canonical _conformal_q ---
    import hiv_quantitative_calibration as qc  # noqa: E402
    probe = np.array([0.1, 0.5, 0.3, 0.9, 0.2, 0.7])
    assert abs(conformal_q(probe, 0.2) - qc._conformal_q(probe, 0.2)) < 1e-12, "conformal_q drift vs J2 helper"

    xs, ys = load_xy()
    x = np.asarray(xs, float); y = np.asarray(ys, float)
    n = len(x)
    rng = np.random.RandomState(SEED)
    covs, widths, narrows, psp = [], [], [], []
    for _ in range(REPEATS):
        idx = rng.permutation(n)
        a, b = n // 2, (3 * n) // 4       # 50% fit / 25% calib / 25% test
        fi, ci, ti = idx[:a], idx[a:b], idx[b:]
        res = evaluate_dosage(x[fi], y[fi], x[ci], y[ci], x[ti], y[ti], coverage=TARGET)
        covs.append(res.coverage); widths.append(res.halfwidth)
        narrows.append(res.interval_narrowing); psp.append(res.point_spearman)

    mean_cov = float(np.mean(covs)); mean_w = float(np.mean(widths))
    mean_narrow = float(np.mean(narrows)); mean_psp = float(np.mean(psp))
    calibrated = abs(mean_cov - TARGET) <= COVER_TOL
    informative = mean_narrow > NARROW_MIN
    verdict = "CALIBRATED_DOSAGE" if (calibrated and informative) else (
        "CALIBRATED_BUT_UNINFORMATIVE" if calibrated else "MISCALIBRATED")

    res = {
        "cell": "forward_dosage_head", "assay": ASSAY, "score_method": "alphamissense",
        "n_variants": n, "target_coverage": TARGET, "repeats": REPEATS,
        "mean_held_out_coverage": round(mean_cov, 4),
        "mean_interval_halfwidth": round(mean_w, 4),
        "mean_interval_narrowing_vs_marginal": round(mean_narrow, 4),
        "mean_point_spearman": round(mean_psp, 4),
        "calibrated": calibrated, "informative": informative, "verdict": verdict,
        "honesty": ("Split-conformal coverage is guaranteed even for a useless model (interval widens to the "
                    "marginal); the narrowing>0 gate is what shows the score actually pins the DMS magnitude. "
                    "conformal_q verified byte-equal to J2's hiv_quantitative_calibration._conformal_q."),
        "interval_units": "DMS score units; interval = calibrated_point +/- halfwidth",
    }
    out = REPO / "wiki" / f"forward_dosage_cell_pten_{_date.today().isoformat()}.json"
    out.write_text(json.dumps(res, indent=2), encoding="utf-8")
    print(f"[forward-dosage] {ASSAY} via AlphaMissense | n={n} target={TARGET}")
    print(f"  mean held-out coverage = {mean_cov:.4f} (|.-target|={abs(mean_cov-TARGET):.4f}, tol {COVER_TOL})")
    print(f"  mean interval halfwidth = {mean_w:.4f} DMS units | narrowing vs marginal = {mean_narrow:.4f}")
    print(f"  mean point Spearman = {mean_psp:.4f} | VERDICT = {verdict}")
    print(f"  artifact -> {out}")
    return 0 if verdict == "CALIBRATED_DOSAGE" else 1


if __name__ == "__main__":
    raise SystemExit(main())
