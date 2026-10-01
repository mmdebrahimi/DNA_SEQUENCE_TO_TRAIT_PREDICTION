"""Leave-one-STUDY-out blind-spot deployability PER GENE — RT, protease, integrase (2026-09-30).

WHY THIS EXISTS. `scripts/build_hiv_complement_model.py` ships three complement models and stamps each with
a `deployability.leave_study_out_blindspot_auroc`, which `dna_decode/eval/doubt.py` renders to a user. Its
docstring calls all three "deployment-validated". Only ONE is: `hiv_supervised_deployability.py` hardcodes
`DRUG = "EFV"` + `NNRTI_DataSet.Full.txt`, so it measured RT only, and the committed artifact
`wiki/hiv_supervised_deployability_2026-07-12.json` reports 0.8102 — the shipped 0.81 for NNRTI is REAL.
The **PI and INSTI 0.89 are hardcoded literals in a config dict with no artifact behind them.**

This measures them, which ALSO closes the largest honest limit on the 2026-09-30 GLM headroom result: every
drug there shared ONE gene (RT), so it generalized across DRUGS, not across GENES. Protease and integrase
are different genes with their own catalogs and their own Stanford datasets.

WHAT IS REUSED, NOT RE-DERIVED — this is what lets the anchor bind. The split (GroupKFold by `RefID`), the
blind-spot metric, the pre-registered pass rule (AUROC >= 0.65 AND > mutation-burden AND shuffled null <
0.55), the one-hot featurizer and the 1%-consensus drift detection are all IMPORTED from
`hiv_supervised_deployability.py` / `hiv_supervised_vs_catalog.py`. The comparator is the DEPLOYED
`hiv_amr.call_hiv_observed(drug, {gene: subs})` for every gene.

THE ANCHOR. RT/EFV is re-measured here and required to reproduce the committed 0.8102. A missed anchor
REFUSES the PI/INSTI numbers (exit 3) — a number produced by a harness that cannot reproduce a known result
is not a measurement. Because the original scored the NNRTI majors set directly while this uses the deployed
caller, their per-isolate agreement is also reported, so the anchor is interpretable rather than a bare pass.

THE PI EXPECTATION IS PRE-STATED, so a null there cannot be dressed up afterwards: `PI` already carries
`THRESHOLD_NEVER_FIRES` in `wiki/supervised_complement_nonvacuity_2026-09-30.json` — the position-based PI
catalog over-calls, which EMPTIES its own blind spot of resistant isolates. A blind spot with too few R
cannot be scored, so this REFUSES an AUROC below `MIN_PER_CLASS` rather than reporting an unpowered one.

Run: uv run python scripts/hiv_deployability_per_gene.py [--self-check]
Exits 2 when the gitignored Stanford data is absent, 3 when the RT anchor does not bind.
"""
from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import sys
from datetime import date
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
for p in (str(REPO), str(REPO / "scripts")):
    if p not in sys.path:
        sys.path.insert(0, p)

RAW = REPO / "data" / "raw" / "hiv"
REF = REPO / "data" / "hiv_ref"
COMMITTED_ANCHOR = REPO / "wiki" / "hiv_supervised_deployability_2026-07-12.json"
ANCHOR_TOL = 0.02
MIN_PER_CLASS = 10          # a blind spot with fewer R than this is not scoreable
CUTOFF = 3.0

# class -> dataset, training drug (the DEPLOYED caller's own drug name), gene key, protein reference
GENES = {
    "NNRTI": {"dataset": "NNRTI_DataSet.Full.txt", "drug": "efavirenz", "col": "EFV",
              "gene": "RT", "protein": "rt", "role": "anchor"},
    "PI": {"dataset": "PI_DataSet.Full.txt", "drug": "lopinavir", "col": "LPV",
           "gene": "PR", "protein": "HIV1_PR_HXB2_cds.fna", "role": "measure"},
    "INSTI": {"dataset": "INI_DataSet.Full.txt", "drug": "raltegravir", "col": "RAL",
              "gene": "IN", "protein": "HIV1_IN_HXB2_cds.fna", "role": "measure"},
}


class PerGeneRefused(RuntimeError):
    """Raised instead of returning a number that cannot be trusted."""


def _harness():
    """The NNRTI deployability harness -- imported for its split + metric + pass rule, not re-written."""
    spec = importlib.util.spec_from_file_location(
        "hiv_supervised_deployability", REPO / "scripts" / "hiv_supervised_deployability.py")
    D = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(D)
    return D


def _protein(src: str, ev) -> str:
    if src == "rt":
        return ev.rt_protein()
    path = REF / src
    if not path.exists():
        raise PerGeneRefused(f"reference absent: {path}")
    seq = "".join(l.strip() for l in open(path, encoding="utf-8") if not l.startswith(">"))
    return "".join(ev.CODON.get(seq[i:i + 3], "X") for i in range(0, len(seq) - 2, 3))


def measure_gene(cls: str, cfg: dict, D) -> dict:
    """One gene: leave-one-STUDY-out CV, then the blind-spot metric on catalog-NEGATIVE isolates."""
    import numpy as np

    import dna_decode.data.hiv_amr as H
    ev, S = D.ev, D.S

    path = RAW / cfg["dataset"]
    if not path.exists():
        raise PerGeneRefused(f"dataset absent: {path}")
    prot = _protein(cfg["protein"], ev)
    rows = list(csv.DictReader(open(path, encoding="utf-8"), delimiter="\t"))
    pcols = [c for c in rows[0] if c.startswith("P") and c[1:].isdigit()]
    col = cfg["col"]
    if col not in rows[0]:
        raise PerGeneRefused(f"fold column {col} absent in {path.name}")
    for meta in ("PtID", "RefID", "Subtype"):
        if meta not in rows[0]:
            raise PerGeneRefused(f"{path.name} lacks {meta} -- leave-study-out needs it")

    have = [r for r in rows if r.get(col) not in ("NA", "", "-", None)]
    if len(have) < 50:
        raise PerGeneRefused(f"{cls}: only {len(have)} rows carry a {col} fold")

    # 1%-consensus drift detection, generalized from the NNRTI harness
    drifted = set()
    for c in pcols:
        p = int(c[1:])
        if p <= len(prot) and sum(1 for r in have if (r[c] or "").strip() == prot[p - 1]) > 0.01 * len(have):
            drifted.add(p)

    gene = cfg["gene"]

    def subs_of(r):
        return {f"{prot[p-1]}{p}{aa}" for p, aa in ev.isolate_muts(r, pcols)
                if p <= len(prot) and p not in drifted and prot[p - 1] != aa}

    # THE DEPLOYED comparator, for every gene
    cat = np.array([1 if H.call_hiv_observed(cfg["drug"], {gene: subs_of(r)}).prediction == "R" else 0
                    for r in have])
    neg = (cat == 0)
    y = np.array([1 if float(r[col]) >= CUTOFF else 0 for r in have])
    burden = [len(ev.isolate_muts(r, pcols)) for r in have]

    X, feats = S.build_onehot(have, pcols, drifted, prot)
    refid = np.array([r["RefID"] for r in have])
    n_studies = len(set(refid.tolist()))

    out = {
        "class": cls, "gene": gene, "drug": cfg["drug"], "fold_column": col,
        "dataset": cfg["dataset"], "protein_len": len(prot), "n_positions_screened": len(pcols),
        "n_rows_with_fold": len(have), "n_R": int(y.sum()), "n_features": len(feats),
        "n_studies": n_studies, "n_patients": len(set(r["PtID"] for r in have)),
        "n_drift_positions": len(drifted),
        "catalog_called_positive_fraction": round(float(cat.mean()), 4),
        "blind_spot_n": int(neg.sum()), "blind_spot_R": int(y[neg].sum()),
        "blind_spot_S": int((neg.sum() - y[neg].sum())),
    }
    if n_studies < 2:
        out["status"] = "NO_STUDY_STRUCTURE"
        return out

    # THE REFUSAL: an over-calling catalog empties its own blind spot of resistant isolates.
    n_r, n_s = out["blind_spot_R"], out["blind_spot_S"]
    if n_r < MIN_PER_CLASS or n_s < MIN_PER_CLASS:
        out["status"] = "BLIND_SPOT_UNSCOREABLE"
        out["refusal"] = (
            f"the catalog-negative subset holds {n_r} R and {n_s} S; fewer than {MIN_PER_CLASS} of a class "
            "cannot support an AUROC, so NO number is reported (an over-calling catalog empties its own "
            "blind spot -- this is a property of the incumbent, not a failure of the complement)")
        return out

    oof = D._fit_predict_grouped(X, y, refid)
    out["overall_auroc"] = round(ev.auroc(y.tolist(), oof.tolist()), 4)
    out["blind_spot"] = D._blindspot_metrics(y.tolist(), oof.tolist(), neg.tolist(), burden)
    out["leave_study_out_blindspot_auroc_MEASURED"] = out["blind_spot"]["auroc"]
    out["status"] = "SCORED"

    if cls == "NNRTI":
        # The original harness scored the NNRTI majors set directly; report agreement so the anchor is
        # interpretable rather than a bare pass/fail.
        majors = H.NNRTI_RT_MAJOR_DRMS
        old = np.array([1 if any(f"{H._RT_WT.get(p,'?')}{p}{aa}" in majors
                                 for p, aa in ev.isolate_muts(r, pcols)) else 0 for r in have])
        out["deployed_vs_majors_agreement"] = f"{int((old == cat).sum())}/{len(have)}"
    return out


def check_anchor(rt: dict) -> dict:
    """The RT cell must reproduce the committed deployability artifact, or PI/INSTI are refused."""
    if not COMMITTED_ANCHOR.exists():
        return {"reproduces": False, "reason": f"committed anchor absent: {COMMITTED_ANCHOR}"}
    d = json.loads(COMMITTED_ANCHOR.read_text(encoding="utf-8"))
    b = (d.get("results", d).get("B_leave_study_out") or {}).get("blind_spot") or {}
    want, got = b.get("auroc"), rt.get("leave_study_out_blindspot_auroc_MEASURED")
    ok = (want is not None and got is not None and abs(got - want) <= ANCHOR_TOL
          and b.get("n") == rt.get("blind_spot_n") and b.get("R") == rt.get("blind_spot_R"))
    return {"reproduces": bool(ok), "committed_auroc": want, "measured_auroc": got,
            "committed_n": b.get("n"), "measured_n": rt.get("blind_spot_n"),
            "committed_R": b.get("R"), "measured_R": rt.get("blind_spot_R"),
            "tolerance": ANCHOR_TOL, "artifact": str(COMMITTED_ANCHOR.relative_to(REPO)).replace("\\", "/")}


# The values `build_hiv_complement_model.py` stamped into every shipped complement model until this
# measurement replaced them. They are HISTORY, not configuration -- `deployability_block()` there now reads
# the measured artifact, so these live here solely to report what is being superseded (and by how much).
# Recorded from the pre-correction builder at commit b67ee2f.
SUPERSEDED_LITERALS = {"NNRTI": 0.81, "PI": 0.89, "INSTI": 0.89}


def shipped_claims() -> dict:
    """What each model's deployability field claimed BEFORE this measurement.

    Also asserts the builder no longer carries those literals, so this cannot quietly become a second live
    source of the same number -- which is the defect being fixed.
    """
    spec = importlib.util.spec_from_file_location(
        "build_hiv_complement_model", REPO / "scripts" / "build_hiv_complement_model.py")
    B = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(B)
    live = {cls: cfg for cls, cfg in B.CLASSES.items() if "leave_study_out" in cfg}
    if live:
        raise PerGeneRefused(
            f"the builder still carries deployability literals for {sorted(live)}; it must DERIVE them "
            "from the measured artifact, or a rebuild reintroduces an unsourced number")
    return dict(SUPERSEDED_LITERALS)


def self_check() -> None:
    D = _harness()
    assert hasattr(D, "_fit_predict_grouped") and hasattr(D, "_blindspot_metrics"), \
        "the split + metric must be imported from the NNRTI harness, never re-written"
    assert D.CUTOFF == CUTOFF, "the fold cutoff must match the harness"
    ship = shipped_claims()
    assert ship.get("PI") == 0.89 and ship.get("INSTI") == 0.89, \
        f"expected the two unsourced literals; got {ship}"
    assert ship.get("NNRTI") == 0.81, f"NNRTI literal changed: {ship}"
    assert MIN_PER_CLASS == 10
    print("self-check OK -- harness imported; shipped literals read:", ship)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", default=str(REPO / "wiki" / f"hiv_deployability_per_gene_{date.today()}.json"))
    ap.add_argument("--self-check", action="store_true")
    a = ap.parse_args(argv)
    if a.self_check:
        self_check()
        return 0

    D = _harness()
    cells, refused = {}, {}
    for cls, cfg in GENES.items():
        try:
            cells[cls] = measure_gene(cls, cfg, D)
        except PerGeneRefused as exc:
            refused[cls] = str(exc)
            print(f"{cls}: REFUSED -- {exc}")

    if "NNRTI" not in cells:
        print("the RT anchor could not be measured; refusing to report any gene")
        return 2

    anchor = check_anchor(cells["NNRTI"])
    print(f"\nANCHOR (RT/EFV vs {anchor.get('artifact')}): "
          f"committed {anchor.get('committed_auroc')} n={anchor.get('committed_n')}/"
          f"R={anchor.get('committed_R')} vs measured {anchor.get('measured_auroc')} "
          f"n={anchor.get('measured_n')}/R={anchor.get('measured_R')} -> "
          f"{'REPRODUCES' if anchor['reproduces'] else 'MISMATCH'}")
    if cells["NNRTI"].get("deployed_vs_majors_agreement"):
        print(f"  deployed call_hiv_observed vs the original majors set: "
              f"{cells['NNRTI']['deployed_vs_majors_agreement']}")
    if not anchor["reproduces"]:
        print("REFUSING the PI/INSTI numbers: a harness that cannot reproduce a known result is not a "
              "measurement.")
        return 3

    ship = shipped_claims()
    print(f"\n{'class':7s} {'gene':5s} {'n':>5s} {'catPos':>7s} {'bsN':>5s} {'bsR':>5s} {'bsS':>5s} "
          f"{'MEASURED':>9s} {'SHIPPED':>8s}  status")
    for cls in GENES:
        c = cells.get(cls)
        if c is None:
            print(f"{cls:7s} {'':5s} {'':>5s} {'':>7s} {'':>5s} {'':>5s} {'':>5s} {'':>9s} "
                  f"{ship.get(cls, 0):8.2f}  REFUSED")
            continue
        m = c.get("leave_study_out_blindspot_auroc_MEASURED")
        print(f"{cls:7s} {c['gene']:5s} {c['n_rows_with_fold']:5d} "
              f"{c['catalog_called_positive_fraction']:7.3f} {c['blind_spot_n']:5d} "
              f"{c['blind_spot_R']:5d} {c['blind_spot_S']:5d} "
              f"{(f'{m:.4f}' if m is not None else 'withheld'):>9s} "
              f"{ship.get(cls, float('nan')):8.2f}  {c['status']}")

    deltas = {}
    for cls, c in cells.items():
        m = c.get("leave_study_out_blindspot_auroc_MEASURED")
        if m is not None and cls in ship:
            deltas[cls] = round(m - ship[cls], 4)
    if deltas:
        print("\nMEASURED minus SHIPPED:", deltas)
    for cls, c in cells.items():
        if c.get("refusal"):
            print(f"\n{cls}: NO NUMBER -- {c['refusal']}")

    art = {
        "schema": "hiv-deployability-per-gene-v1", "date": str(date.today()),
        "question": ("is the shipped leave-one-STUDY-out blind-spot deployability number MEASURED for each "
                     "HIV gene, and does the supervised blind-spot result hold on a SECOND and THIRD gene?"),
        "provenance_of_the_shipped_numbers": (
            "scripts/build_hiv_complement_model.py CLASSES[...]['leave_study_out'] are hardcoded literals "
            "stamped into each model's deployability field and rendered by dna_decode/eval/doubt.py. Only "
            "NNRTI had an artifact behind it (hiv_supervised_deployability.py hardcodes DRUG='EFV' and the "
            "NNRTI dataset, so it measured RT only); PI and INSTI were unsourced."),
        "shipped_literals": ship, "anchor": anchor,
        "measured_minus_shipped": deltas,
        "pass_rule": ("imported verbatim from hiv_supervised_deployability._blindspot_metrics: AUROC >= 0.65 "
                      "AND > mutation-burden AND shuffled null < 0.55"),
        "min_per_class": MIN_PER_CLASS,
        "honest_scope": (
            "IN-DISTRIBUTION to the Stanford knowledge base and ~96% subtype B. Leave-one-STUDY-out is the "
            "out-of-distribution split, NOT an independent cohort. A linear model over one-hot substitution "
            "tokens is the WEAKEST supervised family member, so a pass is a FLOOR on that family. A gene "
            "whose catalog-negative subset holds too few R of either class is REFUSED, not scored."),
        "cells": cells, "refused": refused,
    }
    Path(a.out).write_text(json.dumps(art, indent=2) + "\n", encoding="utf-8")
    print(f"\nwrote {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
