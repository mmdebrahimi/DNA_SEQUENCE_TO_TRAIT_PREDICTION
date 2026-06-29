"""Ingest the OpenSNP archive dump (Internet Archive mirror) → score the rs12913832 eye-colour decoder.

The flagship off-pathogen validation, user-ratified 2026-06-28 (live OpenSNP was deleted 2025-04; this uses
the archived dump at archive.org/details/opensnp_data_dumps). Reads DIRECTLY from the 21 GB zip via
`zipfile` (selective member streaming — the phenotype CSV + only each eye-colour user's genotype member;
NEVER extracts 21 GB). No further network (one download already on D:).

Robust-by-design (the real zip's column/filename conventions are verified at run, not assumed): fuzzy
eye-colour column match + a genotype-filename → user-id regex; `--inspect` prints the real structure first.

Reuses the SHARED rule (`eye_colour.call_eye_colour`, strand-agnostic) + the binner/scorer math
(`eye_colour_opensnp_validate`). HONESTY: self-reported label (near-independent, non-circular, noisy);
ancestry-confounded (rs12913832 European-calibrated) → within-ancestry split is v0.1 (deferred, flagged).
"""
from __future__ import annotations

import csv
import io
import json
import re
import sys
import zipfile
from datetime import date as _date
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from dna_decode.data.eye_colour import EYE_COLOUR_LOCUS, call_eye_colour  # noqa: E402
from scripts.eye_colour_opensnp_validate import bin_eye_colour  # noqa: E402

DEFAULT_ZIP = Path("D:/dna_decode_cache/opensnp/opensnp_datadump.2017-12-08.zip")
_USER_RE = re.compile(r"user(\d+)_", re.IGNORECASE)
# prefer 23andme files (clean rsid\tchrom\tpos\tgenotype); ancestry/ftdna handled too.
_GENO_EXT_PRIORITY = (".23andme.txt", ".23andme-exome-vcf.txt", ".ancestry.txt", ".ftdna-illumina.txt", ".txt")


def _find_phenotype_member(zf: zipfile.ZipFile) -> str | None:
    cands = [n for n in zf.namelist() if "phenotype" in n.lower() and n.lower().endswith(".csv")]
    return sorted(cands, key=len)[0] if cands else None


def _pick_eye_column(hl: list[str]) -> int | None:
    """Pick the free-text eye-COLOUR column, not a decoy.

    The real OpenSNP 2017 phenotype CSV has SIX 'eye'-containing headers — e.g. 'Eye pigmentation'
    (appears BEFORE 'Eye color'), 'Hair and eye color', "Mother's eye color". A naive first-'eye'
    match grabs 'Eye pigmentation' (wrong). Prefer an exact 'eye colo(u)r', then a clean
    colour-bearing 'eye' header (excluding hair/mother/relatives), then any 'eye' header.
    """
    exact = {"eye color", "eye colour", "eye colors", "eye colours"}
    for i, h in enumerate(hl):
        if h.strip() in exact:
            return i
    for i, h in enumerate(hl):
        if "eye" in h and ("color" in h or "colour" in h) and "hair" not in h and "mother" not in h:
            return i
    return next((i for i, h in enumerate(hl) if "eye" in h), None)


def _eye_colour_by_user(zf: zipfile.ZipFile, pheno_member: str) -> dict[str, str]:
    """Parse the OpenSNP phenotype CSV → {user_id: raw eye-colour string} for users with a value.
    OpenSNP CSV is ';'-separated; one row per user; a column whose header contains 'eye'."""
    raw = zf.read(pheno_member).decode("utf-8", errors="replace")
    # sniff delimiter (OpenSNP uses ';'); fall back to ','
    delim = ";" if raw.split("\n", 1)[0].count(";") >= raw.split("\n", 1)[0].count(",") else ","
    rdr = csv.reader(io.StringIO(raw), delimiter=delim)
    header = next(rdr, [])
    hl = [h.strip().lower() for h in header]
    uid_i = next((i for i, h in enumerate(hl) if h in ("user_id", "id", "user")), 0)
    eye_i = _pick_eye_column(hl)
    if eye_i is None:
        return {}
    out = {}
    for row in rdr:
        if len(row) <= max(uid_i, eye_i):
            continue
        uid, val = row[uid_i].strip(), row[eye_i].strip()
        if uid and val and val.lower() not in ("-", "rather not say", "unknown", ""):
            out[uid] = val
    return out


def _genotype_members_by_user(zf: zipfile.ZipFile) -> dict[str, list[str]]:
    by_uid: dict[str, list[str]] = {}
    for n in zf.namelist():
        m = _USER_RE.search(Path(n).name)
        if m and n.lower().endswith(".txt") and "phenotype" not in n.lower():
            by_uid.setdefault(m.group(1), []).append(n)
    return by_uid


def _pick_member(members: list[str]) -> str | None:
    for ext in _GENO_EXT_PRIORITY:
        for n in members:
            if n.lower().endswith(ext):
                return n
    return members[0] if members else None


def _rs_from_member(zf: zipfile.ZipFile, member: str, rsid: str = EYE_COLOUR_LOCUS) -> str | None:
    """Stream a genotype member, return the 2-allele genotype for rsid (handles 23andMe + AncestryDNA shapes)."""
    try:
        with zf.open(member) as fh:
            for bline in fh:
                line = bline.decode("utf-8", errors="replace")
                if not line or line[0] == "#" or line.startswith("rsid") or line.startswith('"rsid'):
                    continue
                if rsid in line:
                    parts = line.rstrip("\n").replace('"', "").replace(",", "\t").split("\t")
                    if parts and parts[0].strip() == rsid:
                        if len(parts) >= 5:
                            return (parts[3].strip() + parts[4].strip()).upper()
                        if len(parts) >= 4:
                            return parts[3].strip().upper()
    except (KeyError, zipfile.BadZipFile, OSError):
        return None
    return None


def run(zip_path: Path, limit: int | None = None, inspect: bool = False) -> dict:
    if not zip_path.exists():
        return {"status": "ZIP_NOT_PRESENT", "zip": str(zip_path),
                "note": "download not complete; re-run when the D: dump finishes."}
    zf = zipfile.ZipFile(str(zip_path))
    pheno_member = _find_phenotype_member(zf)
    geno_by_uid = _genotype_members_by_user(zf)
    if inspect:
        names = zf.namelist()
        return {"status": "INSPECT", "n_members": len(names), "phenotype_member": pheno_member,
                "n_users_with_genotype": len(geno_by_uid),
                "sample_members": names[:8],
                "phenotype_header": (zf.read(pheno_member)[:400].decode("utf-8", "replace") if pheno_member else None)}
    if not pheno_member:
        return {"status": "NO_PHENOTYPE_CSV", "members_sample": zf.namelist()[:10]}
    eye = _eye_colour_by_user(zf, pheno_member)
    strata = {"blue": {"blue": 0, "brown": 0, "intermediate": 0, "indet": 0},
              "brown": {"blue": 0, "brown": 0, "intermediate": 0, "indet": 0}}
    n_other = n_nogeno = n_nors = 0
    scored_users = 0
    for uid, raw_colour in eye.items():
        label = bin_eye_colour(raw_colour)
        if label is None or label == "other":
            n_other += 1; continue
        member = _pick_member(geno_by_uid.get(uid, []))
        if not member:
            n_nogeno += 1; continue
        gt = _rs_from_member(zf, member)
        if not gt:
            n_nors += 1; continue
        pred = call_eye_colour(gt)["prediction"]
        bucket = pred if pred in ("blue", "brown", "intermediate") else "indet"
        strata[label][bucket] += 1
        scored_users += 1
        if limit and scored_users >= limit:
            break
    tp = strata["brown"]["brown"]; fn = strata["brown"]["blue"]
    tn = strata["blue"]["blue"]; fp = strata["blue"]["brown"]
    n = tp + fn + tn + fp
    return {
        "status": "SCORED" if n else "NO_BINARY_PAIRS",
        "schema": "eye-colour-opensnp-archive-v1", "date": _date.today().isoformat(),
        "source": "OpenSNP archive dump (archive.org/opensnp_data_dumps, 2017-12-08); live site deleted 2025-04",
        "rule": f"{EYE_COLOUR_LOCUS} single-locus v0 (strand-agnostic)",
        "label_tier": "self-reported (near-independent, non-circular, noisy)",
        "n_users_eye_labelled": len(eye), "n_binary_scored": n,
        "n_other_excluded": n_other, "n_no_genotype_file": n_nogeno, "n_rs12913832_missing": n_nors,
        "confusion_brown_positive": {"TP": tp, "FP": fp, "TN": tn, "FN": fn},
        "accuracy": round((tp + tn) / n, 3) if n else None,
        "brown_sens": round(tp / (tp + fn), 3) if (tp + fn) else None,
        "blue_spec": round(tn / (tn + fp), 3) if (tn + fp) else None,
        "strata_pred_by_label": strata,
        "caveats": ["self-reported label (not a lab assay)",
                    "ancestry-confounded: rs12913832 European-calibrated; within-ancestry split = v0.1 (deferred)",
                    "intermediate/green/hazel excluded from the binary (reported in strata)",
                    "2017 dump (the last archived snapshot); single-locus v0 (IrisPlex 6-SNP = v0.1)"],
    }


def main(argv=None) -> int:
    import argparse
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--zip", type=Path, default=DEFAULT_ZIP)
    ap.add_argument("--limit", type=int, default=None, help="cap scored users (for a fast first pass)")
    ap.add_argument("--inspect", action="store_true", help="print the zip's real structure + exit")
    a = ap.parse_args(argv)
    res = run(a.zip, limit=a.limit, inspect=a.inspect)
    if not a.inspect and res.get("status") == "SCORED":
        out = REPO / "wiki" / f"eye_colour_opensnp_validation_{_date.today().isoformat()}.json"
        out.write_text(json.dumps(res, indent=2), encoding="utf-8")
        print(f"[wrote {out}]")
    print(json.dumps(res, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
