"""Per-class non-vacuity of the HIV supervised blind-spot complement (2026-09-30).

THE QUESTION. `dna_decode/data/hiv_supervised_complement.py` is wired into the L2 doubt layer for NNRTI /
PI / INSTI. The wiring measured NNRTI only. Does the shipped 0.5 threshold actually separate anything on
the OTHER two classes' blind spots?

METHOD. Partition isolates by the DEPLOYED `call_hiv_observed`; the catalog-SUSCEPTIBLE ones ARE the blind
spot. Ask whether the complement's risk separates the truly-resistant ones (measured Stanford fold >= the
builder's own cutoff) from the rest, and report enrichment over the blind spot's base rate.

Genotype extraction is lifted from the complement's OWN builder (`build_hiv_complement_model.build`:
drifted-position detection -> `ev.isolate_muts` -> `prot[p-1] != aa`) so the tokens scored here are built
the same way as the features the model was trained on.

THE PARTITION GUARD is the point of the script's shape, not decoration. Two earlier hand-written versions
of this scan each produced a confident WRONG design conclusion from a format mismatch, and both failures
had the same signature: every isolate landed in ONE bucket. `--self-check` asserts the guard fires.

TWO DATASETS, deliberately, because they are different questions:
  training   `*_DataSet.Full.txt` -- what the model was FIT on (NNRTI n_train 4222 matches exactly)
  validation `*_DataSet.txt`      -- what the rest of this repo's HIV validation uses
BOTH are in-distribution. The deployable claim stays the shipped leave-one-STUDY-out AUROC; the in-sample
AUROC this script prints is NOT performance and is labelled so in the output.

Run: uv run python scripts/supervised_complement_nonvacuity.py [--dataset validation|training|both]
Data (`data/raw/hiv/*_DataSet*.txt`) is gitignored; exits 2 when absent rather than reporting an empty scan.
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
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))
    sys.path.insert(0, str(REPO / "scripts"))

RAW = REPO / "data" / "raw" / "hiv"
REF = REPO / "data" / "hiv_ref"

# class -> (fold column + trained drug, full drug name, protein source). Mirrors the builder's CLASSES.
CASES: dict[str, dict[str, str]] = {
    "NNRTI": {"stem": "NNRTI_DataSet", "col": "EFV", "drug": "efavirenz", "protein": "rt"},
    "PI": {"stem": "PI_DataSet", "col": "LPV", "drug": "lopinavir", "protein": "HIV1_PR_HXB2_cds.fna"},
    "INSTI": {"stem": "INI_DataSet", "col": "RAL", "drug": "raltegravir", "protein": "HIV1_IN_HXB2_cds.fna"},
}
DATASETS = {"training": ".Full.txt", "validation": ".txt"}


class ScanRefused(RuntimeError):
    """The scan cannot be trusted -- raised instead of returning a number."""


def partition_is_non_vacuous(n_r: int, n_s: int, lo: float = 0.05, hi: float = 0.95) -> bool:
    """A filter that puts (almost) everything in ONE bucket is broken, not a finding.

    THE GUARD BOTH EARLIER VERSIONS OF THIS SCAN LACKED. One swallowed an AttributeError on every row and
    reported '0 catalog-susceptible isolates'; the other used `103N` tokens where the catalog needs `K103N`
    so EVERY isolate became blind-spot and a 0.996 AUROC was really the full cohort. Same signature both
    times. Pure + cheap on purpose: call it before reading any number off the partition.
    """
    if n_r <= 0 or n_s <= 0:
        return False
    return lo < n_r / (n_r + n_s) < hi


def enrichment(precision: float, base_rate: float) -> float:
    """Precision over the blind spot's own base rate -- how much better than flagging at random."""
    if base_rate <= 0:
        raise ScanRefused("base rate 0: no resistant isolate in the blind spot, enrichment is undefined")
    return precision / base_rate


def _load_ev_and_builder():
    spec = importlib.util.spec_from_file_location(
        "hiv_supervised_deployability", REPO / "scripts" / "hiv_supervised_deployability.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.ev


def _protein(src: str, ev):
    if src == "rt":
        return ev.rt_protein()
    seq = "".join(l.strip() for l in open(REF / src) if not l.startswith(">"))
    return "".join(ev.CODON.get(seq[i:i + 3], "X") for i in range(0, len(seq) - 2, 3))


def scan(cls: str, dataset: str) -> dict:
    """One (class, dataset) cell. Raises ScanRefused rather than returning an untrustworthy number."""
    import numpy as np

    from dna_decode.data import hiv_supervised_complement as C
    from dna_decode.data.hiv_amr import call_hiv_observed, gene_for_hiv_drug

    ev = _load_ev_and_builder()
    cfg = CASES[cls]
    path = RAW / (cfg["stem"] + DATASETS[dataset])
    if not path.exists():
        raise ScanRefused(f"dataset absent: {path}")

    cutoff = float(C.model_info(cls).get("cutoff_fold") or 3.0)
    prot = _protein(cfg["protein"], ev)
    rows = list(csv.DictReader(open(path, encoding="utf-8"), delimiter="\t"))
    pcols = [c for c in rows[0] if c.startswith("P") and c[1:].isdigit()]
    drifted = {int(c[1:]) for c in pcols
               if int(c[1:]) <= len(prot)
               and sum(1 for r in rows if (r[c] or "").strip() == prot[int(c[1:]) - 1]) > 0.01 * len(rows)}
    keep = {p for p in range(1, min(max(int(c[1:]) for c in pcols), len(prot)) + 1) if p not in drifted}
    gene = gene_for_hiv_drug(cfg["drug"])

    have = [r for r in rows if r.get(cfg["col"]) not in ("NA", "", "-", None)]
    sus, n_r, n_err, n_wt = [], 0, 0, 0
    for r in have:
        muts = {(p, aa) for p, aa in ev.isolate_muts(r, pcols)
                if p in keep and p <= len(prot) and prot[p - 1] != aa}
        if not muts:
            n_wt += 1
            continue
        try:
            pred = call_hiv_observed(cfg["drug"], {gene: {f"{prot[p - 1]}{p}{aa}" for p, aa in muts}}).prediction
        except Exception:
            n_err += 1
            continue
        if pred == "R":
            n_r += 1
            continue
        sus.append((C.blind_spot_risk({f"{p}{aa}" for p, aa in muts}, drug_class=cls), float(r[cfg["col"]])))

    if n_err:
        raise ScanRefused(f"{n_err} rows errored -- fix the call before reading any number")
    if not partition_is_non_vacuous(n_r, len(sus)):
        raise ScanRefused(
            f"the catalog partition is degenerate (R={n_r}, S={len(sus)}) -- that is the signature of a "
            "format mismatch, which produced two confident wrong conclusions before. Not a finding.")

    risk = np.array([s[0] for s in sus]); fold = np.array([s[1] for s in sus])
    truly = fold >= cutoff
    base = float(truly.mean())
    out = {
        "class": cls, "dataset": dataset, "trained_drug": cfg["col"], "gene": gene,
        "cutoff_fold": cutoff, "n_rows_with_fold": len(have), "n_wildtype_only_skipped": n_wt,
        "n_catalog_r": n_r, "catalog_r_fraction": round(n_r / (n_r + len(sus)), 4),
        "n_blind_spot": len(sus), "n_truly_resistant": int(truly.sum()), "base_rate": round(base, 4),
        "threshold": C.DEFAULT_THRESHOLD, "risk_max": round(float(risk.max()), 4),
        "leave_study_out_blindspot_auroc":
            (C.model_info(cls).get("deployability") or {}).get("leave_study_out_blindspot_auroc"),
    }
    sel = risk >= C.DEFAULT_THRESHOLD
    out["n_flagged"] = int(sel.sum())
    if not sel.sum():
        out["verdict"] = "THRESHOLD_NEVER_FIRES"
    elif truly.sum() == 0:
        out["verdict"] = "NO_RESISTANT_ISOLATE_IN_BLIND_SPOT"
    else:
        prec = float(truly[sel].mean())
        out.update(precision=round(prec, 4),
                   recall=round(float(truly[sel].sum() / truly.sum()), 4),
                   enrichment=round(enrichment(prec, base), 2),
                   verdict="NON_VACUOUS" if prec > base else "NO_BETTER_THAN_BASE_RATE")
    if truly.sum():
        try:
            from sklearn.metrics import roc_auc_score
            # IN-SAMPLE. Named so in the key -- it is NOT the deployable claim.
            out["blindspot_auroc_IN_SAMPLE_not_performance"] = round(float(roc_auc_score(truly, risk)), 4)
        except Exception:
            pass
    return out


def self_check() -> None:
    assert partition_is_non_vacuous(500, 500)
    assert not partition_is_non_vacuous(0, 900), "an empty R bucket must refuse"
    assert not partition_is_non_vacuous(900, 0), "an empty S bucket must refuse"
    assert not partition_is_non_vacuous(999, 1), "a 99.9% bucket must refuse"
    assert abs(enrichment(0.729, 0.047) - 15.51) < 0.02
    try:
        enrichment(0.5, 0.0)
    except ScanRefused:
        pass
    else:                                                     # pragma: no cover
        raise AssertionError("a zero base rate must refuse, not divide")
    print("self-check OK")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dataset", choices=["validation", "training", "both"], default="both")
    ap.add_argument("--out", default=str(REPO / "wiki" /
                                        f"supervised_complement_nonvacuity_{date.today()}.json"))
    ap.add_argument("--self-check", action="store_true")
    a = ap.parse_args(argv)
    if a.self_check:
        self_check()
        return 0

    wanted = ["validation", "training"] if a.dataset == "both" else [a.dataset]
    cells, refused = [], []
    for ds in wanted:
        for cls in CASES:
            try:
                cells.append(scan(cls, ds))
            except ScanRefused as exc:
                refused.append({"class": cls, "dataset": ds, "refused": str(exc)})
                print(f"REFUSED {cls}/{ds}: {exc}")
    if not cells:
        print("no cell produced a trustworthy number (gitignored data absent?)")
        return 2

    for c in cells:
        line = (f"{c['class']:6s} {c['dataset']:11s} blind spot n={c['n_blind_spot']:5d}  "
                f"R={c['n_truly_resistant']:4d}  base={c['base_rate']:.3f}  "
                f"flagged={c['n_flagged']:4d}  ")
        line += (f"prec={c['precision']:.3f} recall={c['recall']:.3f} "
                 f"ENRICH={c['enrichment']:.2f}x" if "enrichment" in c else f"-- {c['verdict']}")
        print(line)
    art = {"schema": "supervised-complement-nonvacuity-v1", "date": str(date.today()),
           "note": ("Per-class non-vacuity of the HIV supervised blind-spot complement on the catalog's own "
                    "blind spot. EVERY number is IN-DISTRIBUTION to the Stanford knowledge base (the "
                    "'training' rows literally so); the deployable claim remains the shipped "
                    "leave-one-STUDY-out AUROC, NOT the in-sample AUROC reported here."),
           "cells": cells, "refused": refused}
    Path(a.out).write_text(json.dumps(art, indent=2) + "\n", encoding="utf-8")
    print(f"\nwrote {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
