"""M2-full: within-European re-score of the eye-colour decoder (disentangle the ancestry confound).

The confound worry (quantified in wiki/eye_colour_ancestry_confound_2026-06-30.md): rs12913832 is strongly
ancestry-informative (EUR blue-allele 0.636 vs EAS 0.002), so in a Europe-majority cohort the SNP could be
predicting eye colour by TAGGING European ancestry rather than by mechanism. The disentangler: restrict to
users who SELF-REPORT European ancestry (a ~homogeneous-ancestry stratum -> little ancestry variance to tag)
and re-score. If accuracy HOLDS within Europeans, the signal is MECHANISTIC, not an ancestry artifact.

Reuses the zip-native OpenSNP helpers + the v0 (rs12913832) rule + the keyword ancestry classifier. HONEST:
self-reported ancestry (a keyword proxy, not genetic inference); MIXED/ambiguous excluded; small non-European
N is a known limit (the whole point is that OpenSNP is Euro-dominated). Emits
wiki/eye_colour_within_ancestry_validation_<date>.json.
"""
from __future__ import annotations

import csv
import io
import json
import sys
import zipfile
from datetime import date as _date
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from dna_decode.data.eye_colour import call_eye_colour  # noqa: E402
from dna_decode.data.eye_colour_ancestry import classify_self_reported_ancestry  # noqa: E402
from scripts.eye_colour_opensnp_ingest import (  # noqa: E402
    DEFAULT_ZIP, _find_phenotype_member, _genotype_members_by_user, _pick_member, _rs_from_member,
)
from scripts.eye_colour_opensnp_validate import bin_eye_colour  # noqa: E402


def _accuracy(pairs: list[tuple[str, str]]) -> dict:
    tp = sum(1 for l, p in pairs if l == "brown" and p == "brown")
    fn = sum(1 for l, p in pairs if l == "brown" and p == "blue")
    tn = sum(1 for l, p in pairs if l == "blue" and p == "blue")
    fp = sum(1 for l, p in pairs if l == "blue" and p == "brown")
    n = tp + fn + tn + fp
    return {"n_binary": n, "TP": tp, "FP": fp, "TN": tn, "FN": fn,
            "accuracy": round((tp + tn) / n, 3) if n else None,
            "brown_sens": round(tp / (tp + fn), 3) if (tp + fn) else None,
            "blue_spec": round(tn / (tn + fp), 3) if (tn + fp) else None}


def run(zip_path: Path) -> dict:
    if not zip_path.exists():
        return {"status": "ZIP_NOT_PRESENT", "zip": str(zip_path)}
    zf = zipfile.ZipFile(str(zip_path))
    pheno = _find_phenotype_member(zf)
    if not pheno:
        return {"status": "NO_PHENOTYPE_CSV"}
    raw = zf.read(pheno).decode("utf-8", errors="replace")
    rdr = csv.reader(io.StringIO(raw), delimiter=";")
    hdr = next(rdr, [])
    hl = [h.strip().lower() for h in hdr]
    i_uid = next((i for i, h in enumerate(hl) if h in ("user_id", "id", "user")), 0)
    i_eye = next((i for i, h in enumerate(hl) if h.strip() in ("eye color", "eye colour")), None)
    i_eth = hl.index("ethnicity") if "ethnicity" in hl else None
    i_anc = hl.index("ancestry") if "ancestry" in hl else None
    if i_eye is None or (i_eth is None and i_anc is None):
        return {"status": "MISSING_COLUMNS", "have_eye": i_eye is not None}

    # user -> (label, ancestry_class)
    users: dict[str, tuple[str, str]] = {}
    mx = max(x for x in (i_uid, i_eye, i_eth, i_anc) if x is not None)
    for row in rdr:
        if len(row) <= mx:
            continue
        uid = row[i_uid].strip()
        label = bin_eye_colour(row[i_eye].strip())
        if not uid or label is None or label == "other":
            continue
        eth = row[i_eth].strip() if i_eth is not None else ""
        anc = row[i_anc].strip() if i_anc is not None else ""
        users[uid] = (label, classify_self_reported_ancestry(eth, anc))

    geno_by_uid = _genotype_members_by_user(zf)
    strata: dict[str, list[tuple[str, str]]] = {"EUROPEAN": [], "NON_EUROPEAN": [], "MIXED_UNKNOWN": []}
    all_pairs: list[tuple[str, str]] = []
    n_nogeno = n_nors = 0
    for uid, (label, anc_class) in users.items():
        member = _pick_member(geno_by_uid.get(uid, []))
        if not member:
            n_nogeno += 1
            continue
        gt = _rs_from_member(zf, member)  # rs12913832 (v0)
        if not gt:
            n_nors += 1
            continue
        pred = call_eye_colour(gt)["prediction"]
        if pred not in ("blue", "brown"):
            continue  # v0 abstains on heterozygotes; within-ancestry test is on the v0 binary calls
        strata[anc_class].append((label, pred))
        all_pairs.append((label, pred))

    eur = _accuracy(strata["EUROPEAN"])
    allc = _accuracy(all_pairs)
    disentangled = (
        eur["n_binary"] >= 30 and eur["accuracy"] is not None and allc["accuracy"] is not None
        and eur["accuracy"] >= allc["accuracy"] - 0.03
    )
    return {
        "status": "SCORED" if all_pairs else "NO_USERS_SCORED",
        "schema": "eye-colour-within-ancestry-v1", "date": _date.today().isoformat(),
        "rule": "v0 rs12913832 (binary calls only), stratified by SELF-REPORTED ancestry",
        "n_no_genotype_file": n_nogeno, "n_rs12913832_missing": n_nors,
        "v0_all_cohort": allc,
        "v0_within_european": eur,
        "v0_non_european": _accuracy(strata["NON_EUROPEAN"]),
        "v0_mixed_unknown": _accuracy(strata["MIXED_UNKNOWN"]),
        "confound_verdict": "CONFOUND_REDUCED_NOT_RESOLVED" if disentangled else
                            ("UNDERPOWERED_WITHIN_EUROPEAN" if eur["n_binary"] < 30 else "ACCURACY_DROPS_WITHIN_EUROPEAN"),
        "interpretation": (
            "within-European accuracy ~= full-cohort accuracy -> CONSISTENT WITH a mechanistic signal within "
            "self-reported Europeans; ancestry confounding is REDUCED, NOT RESOLVED (n small, self-report "
            "keyword proxy not genetic ancestry, homozygote binary calls only). A genome-wide/AIM ancestry "
            "inference + CIs would be needed to resolve it."
            if disentangled else
            "within-European stratum underpowered or accuracy drops -- see n_binary; cannot disentangle."
        ),
        "caveats": [
            "verdict is REDUCED-NOT-RESOLVED: consistent-with-mechanism, not proof the confound is absent",
            "self-reported ancestry via keyword proxy (not genetic ancestry inference)",
            "MIXED/ambiguous/Ashkenazi excluded from EUROPEAN/NON_EUROPEAN strata",
            "OpenSNP is Euro-dominated -> NON_EUROPEAN stratum is small (expected; that IS the point)",
            "v0 binary calls only (homozygotes); heterozygotes abstained as in v0",
        ],
    }


def main(argv=None) -> int:
    import argparse
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--zip", type=Path, default=DEFAULT_ZIP)
    a = ap.parse_args(argv)
    res = run(a.zip)
    if res.get("status") == "SCORED":
        out = REPO / "wiki" / f"eye_colour_within_ancestry_validation_{_date.today().isoformat()}.json"
        out.write_text(json.dumps(res, indent=2), encoding="utf-8")
        print(f"[wrote {out}]")
    print(json.dumps(res, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
