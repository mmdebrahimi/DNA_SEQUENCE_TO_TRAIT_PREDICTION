"""Validate the rs12913832 eye-colour decoder on OpenSNP (free per-individual genotype + self-reported colour).

The flagship off-pathogen cell: the deterministic gene->trait pattern (curated rule + free MEASURED-ish label)
applied to a human visible trait. OpenSNP (opensnp.org) provides, per consenting user, a raw DTC genotype file
(23andMe / AncestryDNA format) + self-reported phenotypes incl. eye colour. This scorer:
  1. parses each raw genotype file -> the rs12913832 call (the SHARED `eye_colour.call_eye_colour`),
  2. bins the free-text eye-colour phenotype -> blue / brown / other (intermediate/green excluded from the
     binary score, reported separately),
  3. scores blue-vs-brown accuracy + a confusion table, with the honest caveats baked in.

HONESTY RAILS (mirroring the AMR independence discipline):
  - SELF-REPORTED label (a `near-independent` tier, NOT a lab assay; still independent of the genome ->
    non-circular). The label noise is real + disclosed.
  - ANCESTRY-CONFOUNDED: rs12913832 is European-calibrated; OpenSNP is ancestry-mixed. The headline is the
    WITHIN-blue/brown binary on self-reported European-typical colours; a within-ancestry split is the v0.1
    de-confound (needs ancestry calls, deferred — flagged, not faked).
  - Intermediate/green eyes are a real third class the single locus only weakly resolves -> reported, not
    forced into the binary.

DATA (gitignored): data/raw/opensnp/ — phenotype CSV + per-user genotype files. ABSENT -> status NOT_FETCHED
(this scorer + the caller are the durable deliverable; the real-data score is the bonus when data lands).
No fabrication: if no data, it reports NOT_FETCHED, never a synthetic number.
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

from dna_decode.data.eye_colour import EYE_COLOUR_LOCUS, call_eye_colour  # noqa: E402

DATA = REPO / "data" / "raw" / "opensnp"

# free-text eye-colour -> binary class. Deliberately CONSERVATIVE: ambiguous/green/hazel -> "other"
# (excluded from the binary score, counted separately). Lowercased substring match.
_BLUE_TERMS = ("blue", "light blue", "blue-grey", "blue-gray", "blue/grey", "grey", "gray")
_BROWN_TERMS = ("brown", "dark brown", "light brown", "dark")
_OTHER_TERMS = ("green", "hazel", "amber", "mixed", "blue-green", "green-blue")


def bin_eye_colour(text: str) -> str | None:
    """free-text self-reported eye colour -> 'blue' / 'brown' / 'other' / None (blank)."""
    t = (text or "").strip().lower()
    if not t:
        return None
    if any(o in t for o in _OTHER_TERMS):           # green/hazel first (a 'blue-green' is 'other', not 'blue')
        return "other"
    if any(b in t for b in _BLUE_TERMS):
        return "blue"
    if any(b in t for b in _BROWN_TERMS):
        return "brown"
    return "other"


def genotype_from_dtc_file(path: Path, rsid: str = EYE_COLOUR_LOCUS) -> str | None:
    """Extract one rsID's genotype from a 23andMe / AncestryDNA raw file. Returns the 2-allele string, or None.

    23andMe: tab-sep `rsid  chrom  pos  genotype`. AncestryDNA: tab-sep `rsid chrom pos allele1 allele2`.
    Both are commented with leading '#'. Robust to either shape.
    """
    try:
        for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
            if not line or line.startswith("#") or line.startswith("rsid"):
                continue
            parts = line.replace(",", "\t").split("\t")
            if parts and parts[0].strip() == rsid:
                if len(parts) >= 5:                  # AncestryDNA: allele1 + allele2 in cols 4,5
                    return (parts[3].strip() + parts[4].strip()).upper()
                if len(parts) == 4:                  # 23andMe: genotype in col 4
                    return parts[3].strip().upper()
    except OSError:
        return None
    return None


def score(data_dir: Path = DATA) -> dict:
    """Score blue-vs-brown on the fetched OpenSNP cohort. NOT_FETCHED if the data dir is absent/empty."""
    pheno_csv = data_dir / "eye_colour_phenotypes.csv"   # columns: user_id, genotype_filename, eye_colour
    if not pheno_csv.exists():
        return {"status": "NOT_FETCHED", "data_dir": str(data_dir),
                "note": ("OpenSNP eye-colour data not present. Fetch per the runbook (a SMALL slice: the "
                         "eye-colour phenotype table + the matching per-user genotype files) into "
                         f"{data_dir}/ then re-run. The caller + this harness are the durable deliverable.")}
    rows = list(csv.DictReader(pheno_csv.open(encoding="utf-8")))
    conf = {"TP": 0, "FP": 0, "TN": 0, "FN": 0}      # R=brown(positive, melanin-present) convention is arbitrary;
    strata = {"blue": {"blue": 0, "brown": 0, "intermediate": 0, "indet": 0},
              "brown": {"blue": 0, "brown": 0, "intermediate": 0, "indet": 0}}
    n_other = n_nogeno = n_nolabel = 0
    for r in rows:
        label = bin_eye_colour(r.get("eye_colour", ""))
        if label is None:
            n_nolabel += 1; continue
        if label == "other":
            n_other += 1; continue
        gfile = data_dir / (r.get("genotype_filename") or "")
        gt = genotype_from_dtc_file(gfile) if gfile.name else None
        if not gt:
            n_nogeno += 1; continue
        pred = call_eye_colour(gt)["prediction"]
        bucket = pred if pred in ("blue", "brown") else ("intermediate" if pred == "intermediate" else "indet")
        strata[label][bucket] += 1
    # binary blue/brown confusion (intermediate + indeterminate excluded, reported separately)
    tp = strata["brown"]["brown"]; fn = strata["brown"]["blue"]
    tn = strata["blue"]["blue"]; fp = strata["blue"]["brown"]
    n = tp + fn + tn + fp
    return {
        "status": "SCORED" if n else "NO_BINARY_PAIRS",
        "rule": f"{EYE_COLOUR_LOCUS} single-locus v0 (strand-agnostic)",
        "label_tier": "self-reported (near-independent; non-circular; noisier than a lab assay)",
        "n_users": len(rows), "n_binary_scored": n,
        "n_other_excluded": n_other, "n_no_genotype": n_nogeno, "n_no_label": n_nolabel,
        "confusion_brown_positive": {"TP": tp, "FP": fp, "TN": tn, "FN": fn},
        "accuracy": round((tp + tn) / n, 3) if n else None,
        "brown_sens": round(tp / (tp + fn), 3) if (tp + fn) else None,
        "blue_spec": round(tn / (tn + fp), 3) if (tn + fp) else None,
        "strata_pred_by_label": strata,
        "caveats": ["self-reported label (not a lab assay)",
                    "ancestry-confounded: rs12913832 is European-calibrated; within-ancestry split is v0.1 (deferred)",
                    "intermediate/green/hazel excluded from the binary (reported in strata)",
                    "single-locus v0; IrisPlex 6-SNP = v0.1 (needs sourced Walsh-2011 coefficients)"],
    }


def main(argv=None) -> int:
    import argparse
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data-dir", type=Path, default=DATA)
    a = ap.parse_args(argv)
    res = score(a.data_dir)
    out = REPO / "wiki" / f"eye_colour_opensnp_validation_{_date.today().isoformat()}.json"
    out.write_text(json.dumps(res, indent=2), encoding="utf-8")
    print(json.dumps(res, indent=2))
    print(f"\n-> {out}")
    return 0 if res.get("status") in ("SCORED", "NOT_FETCHED") else 1


if __name__ == "__main__":
    raise SystemExit(main())
