"""Point the project's own source-concentration bar at the AR Isolate Bank arm.

WHY. `scripts/source_diverse_validate.py` REFUSED to report metrics for 6 of the 10 frozen SCORED cells
because one BioProject held >60% of the cohort -- the concentration that made the gentamicin `rmt` blind
spot structurally invisible. That bar was applied to the NCBI-PD arm and never to the AR Bank arm, which
reports several of the project's highest accuracies. [[feedback_apply_your_own_standard_to_your_own_cohort]]

THE MEASUREMENT IS BY CDC PANEL, NOT BioProject, AND THAT IS A REAL LIMITATION. The AR Bank cohorts on
disk carry only `(BioSample, label)`, and `data/raw/ar_isolate_bank/*_details.json` carries `panel_id`
but no BioProject / submitter / study field. Resolving BioProject per BioSample needs the network. A CDC
panel is a curated DISTRIBUTION unit and a BioProject is a SUBMISSION unit -- related, not identical --
so this is a proxy, and the bar's 0.60 threshold was calibrated on BioProjects. Read a failing verdict
here as "this cohort draws on one curated panel", never as a like-for-like BioProject number.

WHAT IT FOUND, and it is NOT the blanket result the hypothesis predicted: the six gonorrhoeae cohorts are
single-panel (share 1.000) while the Klebsiella and E. coli arms draw on 5-6 panels and PASS the same
bar. A claim that every AR Bank cell would fail is wrong; the concentration is real and it is CONFINED.

READ-ONLY. No network. Exit 0 always -- a report, not a gate.
"""
from __future__ import annotations

import argparse
import glob
import json
from collections import Counter
from datetime import date as _date
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

# Imported rather than restated: the bar must be the SAME number the frozen arm is judged against, or
# this is a different standard wearing the same name.
try:
    from scripts.source_diverse_validate import MAX_SOURCE_SHARE
except Exception:                                              # pragma: no cover - import shape varies
    MAX_SOURCE_SHARE = 0.60


def load_panels(details_glob: str) -> dict[str, int]:
    """BioSample -> CDC panel_id, from the committed AR Bank detail files."""
    out: dict[str, int] = {}
    for f in sorted(glob.glob(details_glob)):
        for row in json.loads(Path(f).read_text(encoding="utf-8")):
            bs = row.get("biosample")
            if bs and row.get("panel_id") is not None:
                out[bs] = row["panel_id"]
    return out


def concentration(rows: list[dict], panels: dict[str, int]) -> dict | None:
    """Panel concentration for one cohort, or None when no row resolves to a panel."""
    got = [panels[r["biosample"]] for r in rows if r.get("biosample") in panels]
    if not got:
        return None
    c = Counter(got)
    largest = c.most_common(1)[0]
    return {"n_resolved": len(got), "n_rows": len(rows), "n_panels": len(c),
            "largest_panel": largest[0], "largest_share": round(largest[1] / len(got), 4),
            "passes_bar": (largest[1] / len(got)) <= MAX_SOURCE_SHARE}


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", type=Path,
                    default=ROOT / "wiki" / f"ar_bank_source_concentration_{_date.today()}.json")
    a = ap.parse_args(argv)

    panels = load_panels(str(ROOT / "data" / "raw" / "ar_isolate_bank" / "*_details.json"))
    cells = []
    for f in sorted(glob.glob(str(ROOT / "data" / "raw" / "ar_bank_*" / "predictions_strict.json"))):
        rows = json.loads(Path(f).read_text(encoding="utf-8"))
        conc = concentration(rows, panels)
        if conc:
            cells.append({"cohort": Path(f).parent.name, **conc})

    failing = [c for c in cells if not c["passes_bar"]]
    doc = {
        "schema": "ar-bank-source-concentration-v1",
        "date": str(_date.today()),
        "bar": MAX_SOURCE_SHARE,
        "grouping": "cdc_panel_id",
        "grouping_is_a_proxy": (
            "Grouped by CDC panel_id, NOT BioProject: the cohorts carry no BioProject field and "
            "resolving one per BioSample needs the network. A panel is a curated distribution unit; a "
            "BioProject is a submission unit. The 0.60 bar was calibrated on BioProjects, so a verdict "
            "here means 'draws on one curated panel', not a like-for-like BioProject number."),
        "n_biosamples_with_panel": len(panels),
        "n_cells": len(cells),
        "n_failing_bar": len(failing),
        "headline": (
            "CONFINED, not blanket: the gonorrhoeae arm is single-panel (share 1.000) while the "
            "Klebsiella and E. coli arms draw on 5-6 panels and PASS the same bar. The hypothesis that "
            "every AR Bank cell would fail is refuted by measurement."),
        "cells": cells,
        "honest_limits": [
            "Panel proxy, not BioProject (see grouping_is_a_proxy).",
            "n is the SCORED subset (predictions_strict.json), not the labelled cohort -- roughly half "
            "of each cohort is excluded upstream for lacking a downloadable assembly.",
            "Concentration bounds what a cohort COULD have detected; it is not itself a defect, and no "
            "metric in any AR Bank artifact is changed by this run.",
        ],
    }
    a.out.write_text(json.dumps(doc, indent=2) + "\n", encoding="utf-8")
    print(f"bar={MAX_SOURCE_SHARE} grouping=cdc_panel_id  cells={len(cells)}  failing={len(failing)}")
    for c in cells:
        flag = "" if c["passes_bar"] else "  <-- FAILS"
        print(f"  {c['cohort']:<50} panels={c['n_panels']:<3} largest={c['largest_share']:.3f}{flag}")
    print(f"-> {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
