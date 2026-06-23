"""Validate the SARS-CoV-2 Mpro catalog against Stanford CoV-RDB MEASURED fold-change (deliverable: the
independent-label test for the coronavirus cell).

Joins the frozen Mpro catalog (`dna_decode/data/sarscov2_amr.MPRO_MAJOR_DRMS`) against CoV-RDB's measured
nirmatrelvir fold-change (`rx_fold` ⋈ `isolate_mutations`, from the committed payload), scores R/S at a stated
fold threshold, and reports sens/spec + a per-mutation measured-fold table.

HONESTY RAILS (load-bearing, mirror the TB CRyPTIC baseline):
  1. IN-DISTRIBUTION, NOT independent. The catalog is built from CoV-RDB selection records and the fold-change
     is ALSO CoV-RDB -> status `COV_RDB_IN_DISTRIBUTION_KNOWLEDGE_BASELINE`. A truly-independent number needs
     fold-change from a source NOT used to build the catalog (held-out study / clinical isolates) = v0.1.
  2. UNDERPOWERED. Clinical nirmatrelvir resistance is rare; the set is ~37 Mpro-mutant isolates -> a powering
     verdict, never a confident headline. Wilson CI reported.
  3. OPERATOR-AWARE fold censoring (the MIC-censoring lesson): a '>' bound is never falsely called S, a '<'
     bound never falsely R.
  4. The catalog is SELECTION-DERIVED -> expect over-call (weak passengers like T21I/A173V have fold ~1 but
     are catalogued); the per-mutation table surfaces the v0.1 pruning candidates.

Pure logic (`binarize_fold` / `catalog_predict` / `score`) is unit-tested; the CoV-RDB load is the live part.
"""
from __future__ import annotations

import csv
import io
import sys
import zipfile
from collections import Counter
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from dna_decode.data.sarscov2_amr import MPRO_MAJOR_DRMS  # noqa: E402

DEFAULT_PAYLOAD = Path("D:/dna_decode_cache/data files donwload/covid-drdb-payload/payload.zip")
MPRO_GENES = ("_3CLpro", "3CLpro", "Mpro")
# catalog as a (position, mutant_aa) set — match isolate_mutations (gene,pos,aa) directly (numbering-robust)
_CATALOG_POSMUT = {(int("".join(c for c in s[1:-1] if c.isdigit())), s[-1]) for s in MPRO_MAJOR_DRMS}


def binarize_fold(fold: str, cmp: str, threshold: float) -> str | None:
    """Measured fold-change -> 'R' / 'S' / None, OPERATOR-AWARE (the MIC-censoring lesson).

    '=' v: R iff v>=threshold. '>' v (lower bound): R if v>=threshold, else AMBIGUOUS (None) — a '>2' with
    threshold 2.5 can't be called S. '<' v (upper bound): S if v<threshold, else None."""
    try:
        v = float(str(fold).strip())
    except (TypeError, ValueError):
        return None
    op = (cmp or "=").strip()
    if op in ("=", "~", ""):
        return "R" if v >= threshold else "S"
    if op in (">", ">="):
        return "R" if v >= threshold else None          # lower bound below threshold -> can't call S
    if op in ("<", "<="):
        return "S" if v < threshold else None            # upper bound above threshold -> can't call R
    return None


def catalog_predict(mpro_posmut: set[tuple[int, str]]) -> str:
    """R iff the isolate carries >=1 catalogued Mpro major substitution (by position+mutant)."""
    return "R" if (mpro_posmut & _CATALOG_POSMUT) else "S"


def score(records: list[dict], threshold: float) -> dict:
    """records: [{iso, fold, cmp, mpro_posmut}]. Returns confusion + sens/spec + per-mutation fold table."""
    tp = fp = tn = fn = 0
    scored = 0
    permut: dict[tuple, list] = {}
    for r in records:
        pheno = binarize_fold(r["fold"], r.get("cmp", "="), threshold)
        for pm in r["mpro_posmut"]:
            permut.setdefault(pm, []).append(r["fold"])
        if pheno is None:
            continue
        scored += 1
        pred = catalog_predict(r["mpro_posmut"])
        if pred == "R" and pheno == "R": tp += 1
        elif pred == "R" and pheno == "S": fp += 1
        elif pred == "S" and pheno == "S": tn += 1
        else: fn += 1
    sens = tp / (tp + fn) if (tp + fn) else None
    spec = tn / (tn + fp) if (tn + fp) else None
    return {"threshold": threshold, "n_records": len(records), "n_scored": scored,
            "tp": tp, "fp": fp, "tn": tn, "fn": fn, "sens": sens, "spec": spec,
            "n_R": tp + fn, "n_S": tn + fp,
            "per_mutation_fold": {f"{p}{a}": sorted(set(v)) for (p, a), v in sorted(permut.items())}}


def _read_csv(zf: zipfile.ZipFile, name: str):
    return csv.DictReader(io.StringIO(zf.read(name).decode("utf-8-sig")))


def load_records(payload_zip: Path, drug: str = "Nirmatrelvir") -> list[dict]:
    """CoV-RDB payload -> [{iso, fold, cmp, mpro_posmut}] for `drug` fold records on Mpro-mutant isolates."""
    zf = zipfile.ZipFile(str(payload_zip))
    names = zf.namelist()
    isomut: dict[str, set] = {}
    for n in names:
        if "/isolate_mutations.d/" in n and n.endswith(".csv"):
            for r in _read_csv(zf, n):
                if r.get("gene") in MPRO_GENES:
                    try:
                        isomut.setdefault(r["iso_name"], set()).add((int(r["position"]), r["amino_acid"].strip()))
                    except (ValueError, KeyError):
                        pass
    out = []
    for n in names:
        if "/rx_fold/" in n and n.endswith(".csv"):
            for r in _read_csv(zf, n):
                if r.get("rx_name") == drug and r.get("iso_name") in isomut:
                    out.append({"iso": r["iso_name"], "fold": r.get("fold"),
                                "cmp": r.get("fold_cmp", "="), "mpro_posmut": isomut[r["iso_name"]]})
    return out


def main(argv=None) -> int:
    import argparse
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--payload", type=Path, default=DEFAULT_PAYLOAD, help="CoV-RDB payload zip")
    ap.add_argument("--drug", default="Nirmatrelvir")
    ap.add_argument("--threshold", type=float, default=2.5, help="fold-change R cutoff (reduced suscept.)")
    ap.add_argument("--min-per-class", type=int, default=10)
    a = ap.parse_args(argv)
    if not a.payload.exists():
        print(f"ERROR: CoV-RDB payload not found at {a.payload} (download "
              "codeload.github.com/hivdb/covid-drdb-payload/zip/refs/heads/master)", file=sys.stderr)
        return 2
    recs = load_records(a.payload, a.drug)
    res = score(recs, a.threshold)
    powered = (res["n_R"] >= a.min_per_class and res["n_S"] >= a.min_per_class)
    res["status"] = "COV_RDB_IN_DISTRIBUTION_KNOWLEDGE_BASELINE"
    res["powering"] = "POWERED" if powered else "UNDERPOWERED"
    print(f"[sarscov2-mpro-validate] {a.drug} @ fold>={a.threshold} | status={res['status']} "
          f"powering={res['powering']}")
    print(f"  records={res['n_records']} scored={res['n_scored']}  R={res['n_R']} S={res['n_S']}")
    print(f"  confusion: TP={res['tp']} FP={res['fp']} TN={res['tn']} FN={res['fn']}")
    print(f"  sens={res['sens']} spec={res['spec']}")
    print("  per-mutation measured fold (selection-derived catalog over-call surfaces here):")
    for m, folds in res["per_mutation_fold"].items():
        print(f"    {m}: {folds}")
    import json
    out = REPO / "wiki" / "sarscov2_mpro_cov_rdb_validation.json"
    out.write_text(json.dumps({k: v for k, v in res.items() if k != "per_mutation_fold"}
                              | {"per_mutation_fold": res["per_mutation_fold"]}, indent=2, default=str),
                   encoding="utf-8")
    print(f"artifact -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
