"""Build the shippable HIV supervised blind-spot COMPLEMENT models — NNRTI + PI + INSTI (2026-07-12).

Recommendation 2 (integrate), extended across HIV genes. The catalog fold-in was REJECTED (hard rules trade
sens for spec); the value is the WEIGHTED continuous score, so we ship THAT: train the full-sequence logistic
on the free Stanford fold-change label for a strong well-powered drug per class, then SERIALIZE the learned
per-mutation weights (`{"<pos><aa>": coef}` + intercept) to a small committable JSON. The complement scorer
(`dna_decode/data/hiv_supervised_complement.py`) loads a class's JSON and scores a genotype OFFLINE — no
training data, no sklearn at inference.

Classes: NNRTI (RT, train EFV), PI (protease, train LPV), INSTI (integrase, train RAL). Each is a RANKING
complement, NOT a hard R/S rule and NOT a catalog replacement. Frozen decoder surface + hiv_amr catalog
untouched.

DEPLOYABILITY (corrected 2026-09-30). This docstring used to call all three "deployment-validated" and
quote leave-study-out 0.81 / 0.89 / 0.89. Only the first was measured — `hiv_supervised_deployability.py`
hardcodes DRUG="EFV" and the NNRTI dataset, so it covered RT alone; the PI and INSTI figures were literals
in the CLASSES dict with no artifact behind them, stamped into each shipped model and rendered to users by
`dna_decode/eval/doubt.py`. `scripts/hiv_deployability_per_gene.py` measures every gene and
`deployability_block()` reads its result, so:

  NNRTI  measured 0.8102 on 1883 blind-spot isolates / 201 R over 84 studies — the literal was the real
         value rounded, and it is also the harness ANCHOR (it reproduces the committed 2026-07-12 artifact
         exactly, n and R identical).
  INSTI  measured 0.8923 on 844 / 59 R over 69 studies, PASSES the same pre-registered bar. The unsourced
         0.89 turns out ACCURATE — this is a confirmation, not a catch. It is also the SECOND gene for the
         supervised blind-spot result, whose evidence had been RT-only.
  PI     NO NUMBER. Its own position-based catalog calls 74.2% of isolates resistant, which empties its
         blind spot to 3 R of 910 — below the 10-per-class floor an AUROC needs. So the 0.89 was not only
         unsourced, it is not measurable on this substrate; the field now carries None plus the reason.
"""
from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import sys
from datetime import date as _date
from pathlib import Path

import numpy as np
from sklearn.linear_model import LogisticRegression

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO)); sys.path.insert(0, str(REPO / "scripts"))

spec = importlib.util.spec_from_file_location(
    "hiv_supervised_deployability", REPO / "scripts" / "hiv_supervised_deployability.py")
D = importlib.util.module_from_spec(spec); spec.loader.exec_module(D)
S, ev = D.S, D.ev

RAW = REPO / "data" / "raw" / "hiv"
REF = REPO / "data" / "hiv_ref"
CUTOFF, C = 3.0, 1.0

# class -> (dataset, training drug, protein source, deployability blind-spot AUROC, output JSON)
#
# THE `leave_study_out` VALUES WERE HARDCODED LITERALS UNTIL 2026-09-30, and the docstring above called
# all three "deployment-validated" when only NNRTI was: `hiv_supervised_deployability.py` hardcodes
# DRUG="EFV" + the NNRTI dataset, so it measured RT only. `scripts/hiv_deployability_per_gene.py` now
# measures every gene (its RT anchor reproduces the committed 0.8102 exactly), and `deployability_block()`
# below reads the MEASURED value out of that artifact rather than restating a number here -- so a rebuild
# cannot silently reintroduce an unsourced figure.
DEPLOYABILITY_ARTIFACT = REPO / "wiki" / "hiv_deployability_per_gene_2026-09-30.json"

CLASSES = {
    "NNRTI": {"dataset": "NNRTI_DataSet.Full.txt", "drug": "EFV", "protein": "rt",
              "out": "hiv_nnrti_supervised_complement.json"},
    "PI": {"dataset": "PI_DataSet.Full.txt", "drug": "LPV", "protein": "HIV1_PR_HXB2_cds.fna",
           "out": "hiv_pi_supervised_complement.json"},
    "INSTI": {"dataset": "INI_DataSet.Full.txt", "drug": "RAL", "protein": "HIV1_IN_HXB2_cds.fna",
              "out": "hiv_insti_supervised_complement.json"},
}


def deployability_block(cls: str) -> dict:
    """The MEASURED deployability block for `cls`, read from the per-gene artifact.

    A gene whose catalog-negative blind spot cannot support an AUROC (the PI case: its own catalog
    over-calls, leaving 3 R of 910) carries `leave_study_out_blindspot_auroc: None` plus the reason --
    never a number. RAISES if the artifact is missing or does not cover the class, because emitting a
    model with no deployability provenance is how the unsourced literals shipped in the first place.
    """
    if not DEPLOYABILITY_ARTIFACT.exists():
        raise FileNotFoundError(
            f"deployability artifact absent: {DEPLOYABILITY_ARTIFACT}. Run "
            "`uv run python scripts/hiv_deployability_per_gene.py` first -- a complement model must not "
            "ship an unsourced deployability number.")
    art = json.loads(DEPLOYABILITY_ARTIFACT.read_text(encoding="utf-8"))
    if not (art.get("anchor") or {}).get("reproduces"):
        raise ValueError(f"{DEPLOYABILITY_ARTIFACT.name}: the RT anchor did not bind; refusing to stamp "
                         "deployability from a harness that cannot reproduce a known result")
    cell = (art.get("cells") or {}).get(cls)
    if cell is None:
        raise KeyError(f"{DEPLOYABILITY_ARTIFACT.name} does not cover class {cls}")
    rel = str(DEPLOYABILITY_ARTIFACT.relative_to(REPO)).replace("\\", "/")
    auroc = cell.get("leave_study_out_blindspot_auroc_MEASURED")
    blk = {"leave_study_out_blindspot_auroc": auroc, "measured": auroc is not None, "artifact": rel,
           "blind_spot_n": cell.get("blind_spot_n"), "blind_spot_R": cell.get("blind_spot_R"),
           "n_studies": cell.get("n_studies")}
    if auroc is None:
        blk["blind_spot_status"] = "UNSCOREABLE"
        blk["blind_spot_S"] = cell.get("blind_spot_S")
        blk["catalog_called_positive_fraction"] = cell.get("catalog_called_positive_fraction")
        blk["reason"] = cell.get("refusal")
    else:
        blk["pass_rule"] = art.get("pass_rule")
        blk["passes"] = (cell.get("blind_spot") or {}).get("pass")

    # The class that anchors the harness cites the independent prior artifact it reproduces.
    anchor = art.get("anchor") or {}
    if cell.get("class") == "NNRTI" and anchor.get("artifact"):
        blk["corroborating_artifact"] = anchor["artifact"]

    # Provenance history, DERIVED from the artifact rather than hand-written: what this field used to
    # claim, and whether measuring it confirmed or corrected that claim.
    shipped = (art.get("shipped_literals") or {}).get(cls)
    if shipped is not None:
        if auroc is None:
            blk["superseded_literal"] = (
                f"this field shipped an unsourced literal {shipped} until {art.get('date')}; it is NOT "
                "measurable on this substrate (see `reason`), so no number replaces it")
        else:
            delta = art.get("measured_minus_shipped", {}).get(cls)
            verb = ("CONFIRMS it" if delta is not None and abs(delta) <= 0.01
                    else f"CORRECTS it (delta {delta:+})")
            blk["superseded_literal"] = (
                f"this field shipped the literal {shipped} until {art.get('date')}; measurement {verb}")
    return blk


def _protein(src):
    if src == "rt":
        return ev.rt_protein()
    seq = "".join(l.strip() for l in open(REF / src) if not l.startswith(">"))
    return "".join(ev.CODON.get(seq[i:i + 3], "X") for i in range(0, len(seq) - 2, 3))


def build(cls: str) -> Path:
    cfg = CLASSES[cls]
    prot = _protein(cfg["protein"])
    rows = list(csv.DictReader(open(RAW / cfg["dataset"], encoding="utf-8"), delimiter="\t"))
    pcols = [c for c in rows[0] if c.startswith("P") and c[1:].isdigit()]
    drifted = set()
    for c in pcols:
        p = int(c[1:])
        if p <= len(prot) and sum(1 for r in rows if (r[c] or "").strip() == prot[p - 1]) > 0.01 * len(rows):
            drifted.add(p)
    drug = cfg["drug"]
    have = [r for r in rows if r.get(drug) not in ("NA", "", "-", None)]
    X, feats = S.build_onehot(have, pcols, drifted, prot)
    y = np.array([1 if float(r[drug]) >= CUTOFF else 0 for r in have])
    clf = LogisticRegression(max_iter=3000, C=C, solver="liblinear").fit(X, y)
    weights = {f"{p}{aa}": round(float(w), 5) for (p, aa), w in zip(feats, clf.coef_[0])}
    model = {
        "schema": "hiv-supervised-complement-v1", "drug_class": cls, "date": str(_date.today()),
        "gene": {"NNRTI": "RT", "PI": "PR", "INSTI": "IN"}[cls],
        "drug_trained": drug, "cutoff_fold": CUTOFF, "solver": "liblinear-logistic-L2", "C": C,
        "label_source": "Stanford HIVDB PhenoSense fold-change (free, independent)",
        "n_train": len(have), "n_features": len(weights), "intercept": round(float(clf.intercept_[0]), 5),
        "feature_key": "<pos><aa> of a NON-WT residue (e.g. '103N'); risk = sigmoid(intercept + sum coef)",
        "deployability": deployability_block(cls),
        "honest_scope": (f"Blind-spot RANKING complement for {cls} (trained on {drug}). NOT a hard R/S rule; "
                         "complements — does not replace — the deployed catalog. In-distribution to Stanford; "
                         "supervised (needs the free label to have been trained)."),
        "weights": weights,
    }
    out = REF / cfg["out"]
    out.write_text(json.dumps(model, indent=1), encoding="utf-8")
    top = sorted(weights.items(), key=lambda kv: -kv[1])[:6]
    print(f"[{cls}] trained {drug} n={len(have)}; {len(weights)} weights; intercept {model['intercept']}; "
          f"top {top}; -> {out} ({out.stat().st_size//1024} KB)")
    return out


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--class", dest="cls", choices=list(CLASSES) + ["all"], default="all")
    a = ap.parse_args(argv)
    for cls in (list(CLASSES) if a.cls == "all" else [a.cls]):
        build(cls)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
