"""Paired naive-AMRFinder baseline vs the FROZEN call_resistance rule on the Oxford cohort.

The validate-wrapper-vs-underlying-tool rail (committed lesson): a curated POLICY LAYER over
an external curated-DB tool must BEAT *naive use of that tool* on INDEPENDENT data — else the
in-cohort number only proves the tool works, not that the layer adds value. The 2026-06-15
Oxford run reported the frozen decoder number alone; this closes that gap by scoring a paired
naive baseline on the SAME measured-MIC labels and reporting the delta.

- naive = "ANY AMRFinder determinant whose Class matches the drug's AMR class -> R" (the
  non-expert use of the raw tool; the broad `mic_tiers.amrfinder_classes_for(drug)` match,
  WITHOUT the curated subclass/point/threshold refinement). No abstain — naive calls everything.
- frozen = the shipped `call_resistance` DRUG_RULE (subclass_any / qrdr_point / threshold curation).

Reuses data/raw/oxford/ + scripts/oxford_score.py helpers. No frozen-surface change; no download.
"""
from __future__ import annotations

import json
import sys
import tempfile
from datetime import date as _date
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from dna_decode.data.mic_tiers import amrfinder_classes_for
from dna_decode.eval.amr_rules import call_resistance
from scripts.independent_cohort_validate import _conf
from scripts.oxford_score import (
    PILOT_DRUGS,
    REGISTRY_ORGANISM,
    group_amrfinder,
    load_mic_labels,
)


def _class_col_index(header: str) -> int:
    cols = header.split("\t")
    return cols.index("Class")


def naive_predict(guuid: str, drug: str, amr: dict[str, list[str]], class_idx: int) -> str:
    """Non-expert use of AMRFinder: R iff ANY determinant row's Class is in the drug's
    AMR-class set. No subclass refinement, no point-mutation specificity, no abstain."""
    wanted = amrfinder_classes_for(drug)
    for ln in amr.get(guuid, []):
        parts = ln.split("\t")
        if class_idx < len(parts) and parts[class_idx].strip().upper() in wanted:
            return "R"
    return "S"


def main() -> int:
    base = Path("data/raw/oxford")
    _, binary_labels, n_total = load_mic_labels(base / "main_data.csv")
    header, amr = group_amrfinder(base / "amrfinder.tsv")
    class_idx = _class_col_index(header)
    print(f"loaded {n_total} MIC rows; {len(amr)} guuids with AMRFinder hits; Class col @ {class_idx}")

    tmpdir = Path(tempfile.mkdtemp(prefix="oxford_naive_"))
    frozen_cache: dict[tuple[str, str], str] = {}

    def frozen_predict(guuid: str, drug: str) -> str:
        key = (guuid, drug)
        if key in frozen_cache:
            return frozen_cache[key]
        p = tmpdir / f"{guuid}.tsv"
        if not p.exists():
            p.write_text(header + "\n" + "\n".join(amr.get(guuid, [])) + "\n", encoding="utf-8")
        pred = call_resistance(p, drug, organism=REGISTRY_ORGANISM)["prediction"]
        frozen_cache[key] = pred
        return pred

    drugs_out = {}
    for drug in PILOT_DRUGS:
        labelmap = binary_labels[drug]
        frozen_pairs = [(frozen_predict(g, drug), 1 if rs == "R" else 0) for g, rs in labelmap.items()]
        naive_pairs = [(naive_predict(g, drug, amr, class_idx), 1 if rs == "R" else 0)
                       for g, rs in labelmap.items()]
        fz, nv = _conf(frozen_pairs), _conf(naive_pairs)
        fz_balacc = round((fz["sens"] + fz["spec"]) / 2, 4)
        nv_balacc = round((nv["sens"] + nv["spec"]) / 2, 4)
        delta = {
            "acc": round(fz["acc"] - nv["acc"], 4),
            "sens": round(fz["sens"] - nv["sens"], 4),
            "spec": round(fz["spec"] - nv["spec"], 4),
            "balacc": round(fz_balacc - nv_balacc, 4),
        }
        # value-add verdict on BALANCED accuracy (the honest net metric): a naive baseline can
        # game one axis (call everything R -> sens~1, spec~0), so a per-axis sens guard is the
        # wrong test. balacc = (sens+spec)/2 nets out the over-call. The curated layer ADDS
        # VALUE iff it strictly beats naive on balanced accuracy by a clear margin.
        if delta["balacc"] >= 0.03:
            verdict = "CURATED_LAYER_ADDS_VALUE"
        elif delta["balacc"] <= -0.03:
            verdict = "NAIVE_BEATS_CURATED"
        else:
            verdict = "NAIVE_TIES_CURATED"
        drugs_out[drug] = {"frozen": {**fz, "balacc": fz_balacc}, "naive": {**nv, "balacc": nv_balacc},
                           "delta_frozen_minus_naive": delta, "value_add_verdict": verdict}
        print(f"\n{drug}: n={fz['n_scored']}/{nv['n_scored']} (frozen/naive)")
        print(f"  frozen  acc={fz['acc']} sens={fz['sens']} spec={fz['spec']} balacc={fz_balacc}")
        print(f"  naive   acc={nv['acc']} sens={nv['sens']} spec={nv['spec']} balacc={nv_balacc}")
        print(f"  delta   balacc={delta['balacc']:+} spec={delta['spec']:+} sens={delta['sens']:+}  -> {verdict}")

    out = {
        "_schema": "external-validation-naive-comparator-v1",
        "date": _date.today().isoformat(),
        "cohort": "oxford_lipworth_ecoli_mic_arg",
        "organism": REGISTRY_ORGANISM,
        "baseline_definition": ("naive = R iff ANY AMRFinder determinant Class in "
                                "mic_tiers.amrfinder_classes_for(drug); no subclass/point/threshold "
                                "refinement, no abstain. frozen = shipped call_resistance DRUG_RULE."),
        "rail": "validate-wrapper-vs-underlying-tool: the curated layer must beat naive tool use on independent data.",
        "drugs": drugs_out,
        "source": "github.com/samlipworth/ecoli_mic_arg (PRJNA604975 + PRJNA1007570)",
        "frozen_surface_changed": False,
    }
    outp = Path(f"wiki/external_validation_oxford_naive_comparator_{_date.today().isoformat()}.json")
    outp.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(f"\nartifact -> {outp}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
