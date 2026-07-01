"""M3: cross-cohort REPLICATION of the eye-colour decoder on Personal Genome Project (PGP-Harvard) data.

PGP is the ethically-clean OpenSNP successor (open-consent, CC0 public-domain, participants explicitly
consented to public release). This runs the SAME v0 (rs12913832) + v0.1 (IrisPlex 6-SNP) callers on a
SECOND, independent cohort -> does the eye-colour rule replicate? Uses PGP's "Basic Phenotypes 2015" survey
(google_surveys/19) for the eye-colour label + each participant's donated 23andMe/AncestryDNA file (linked
from their profile). BOUNDED PILOT by design (--limit): replication of a near-deterministic rule doesn't need
the full 1127-participant cohort, and a bounded run is polite to PGP's servers (rate-limited, cached to D:).

HONEST: self-reported eye colour (as OpenSNP); complete-case for v0.1 (all 6 SNPs); a pilot subset (N bounded)
-> a replication SIGNAL, not a full-cohort number. Emits wiki/eye_colour_pgp_validation_<date>.json.
"""
from __future__ import annotations

import csv
import io
import json
import re
import sys
import time
import urllib.request
import zipfile
from datetime import date as _date
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from dna_decode.data.eye_colour import EYE_COLOUR_LOCUS, call_eye_colour  # noqa: E402
from dna_decode.data.eye_colour_irisplex import _SNP_ORDER, predict_irisplex  # noqa: E402
from scripts.eye_colour_opensnp_validate import bin_eye_colour  # noqa: E402

# PRIVACY / RETENTION POSTURE (public-but-sensitive human genotype data):
# - Downloaded DTC files are cached to D:/dna_decode_cache/pgp/ -- OUTSIDE the repo (never git-tracked;
#   .gitignore's /data/* + the external drive both keep raw genomes out of version control).
# - Source is PGP-Harvard CC0 open-consent data (participants explicitly consented to public release).
# - Local-only research cache; not synced cross-machine. Purge with:  rm -rf D:/dna_decode_cache/pgp/
#   (safe -- files are re-fetchable from PGP). Do NOT commit the cache or copy it off the local host.
BASE = "https://my.pgp-hms.org"
SURVEY_URL = f"{BASE}/google_surveys/19/download"
CACHE = Path("D:/dna_decode_cache/pgp")
_UA = {"User-Agent": "dna-decode-eye-colour-replication/1.0 (research; polite; rate-limited)"}
_DL_RE = re.compile(r'data-summarize-as="name">([^<]+)</td>\s*<td data-summarize-as="size">\s*'
                    r'<a href="/user_file/download/(\d+)"', re.S)


def _get(url: str, timeout: int = 90) -> bytes:
    req = urllib.request.Request(url, headers=_UA)
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read()


def survey_eye_labels(csv_bytes: bytes) -> dict[str, str]:
    """{huID: binned label (blue/brown/other)} from PGP Basic Phenotypes 2015 (left eye, else right)."""
    rows = list(csv.reader(io.StringIO(csv_bytes.decode("utf-8", errors="replace"))))
    hdr = rows[0]
    li = next(i for i, h in enumerate(hdr) if "Left Eye Color - Text" in h)
    ri = next(i for i, h in enumerate(hdr) if "Right Eye Color - Text" in h)
    out: dict[str, str] = {}
    for r in rows[1:]:
        if len(r) <= max(li, ri):
            continue
        hu = r[0].strip()
        raw = (r[li].strip() or r[ri].strip())
        if hu and raw:
            lab = bin_eye_colour(raw)
            if lab:
                out[hu] = lab
    return out


def profile_genotype_file(hu: str) -> tuple[int, str] | None:
    """Return (download_id, kind) for a participant's best DTC genotype file, prefer 23andMe then Ancestry."""
    try:
        html = _get(f"{BASE}/profile/{hu}").decode("utf-8", errors="replace")
    except Exception:
        return None
    files = [(name.lower(), int(fid)) for name, fid in _DL_RE.findall(html)]
    for name, fid in files:
        if "23andme" in name:
            return fid, "23andme"
    for name, fid in files:
        if "ancestry" in name or (name.endswith(".zip") and "dna" in name):
            return fid, "ancestry"
    return None


def _snps_from_lines(lines, rsids: set[str]) -> dict[str, str]:
    out: dict[str, str] = {}
    want = set(rsids)
    for line in lines:
        if not want:
            break
        if not line or line[0] == "#" or line.startswith("rsid") or line.startswith('"rsid'):
            continue
        parts = line.rstrip("\n").replace('"', "").replace(",", "\t").split("\t")
        if not parts:
            continue
        rid = parts[0].strip()
        if rid in want:
            if len(parts) >= 5:
                out[rid] = (parts[3].strip() + parts[4].strip()).upper()
            elif len(parts) >= 4:
                out[rid] = parts[3].strip().upper()
            want.discard(rid)
    return out


def fetch_snps(fid: int, kind: str, rsids: set[str]) -> dict[str, str]:
    """Download (cached) a participant genotype file + extract the requested rsids. Handles zipped files."""
    CACHE.mkdir(parents=True, exist_ok=True)
    cached = CACHE / f"{fid}.bin"
    if cached.exists():
        data = cached.read_bytes()
    else:
        data = _get(f"{BASE}/user_file/download/{fid}")
        cached.write_bytes(data)
        time.sleep(1.0)  # polite: only sleep on an actual network fetch
    if data[:2] == b"PK":  # zip (AncestryDNA export, sometimes 23andMe)
        try:
            zf = zipfile.ZipFile(io.BytesIO(data))
            member = next((n for n in zf.namelist() if n.lower().endswith(".txt")), None)
            if not member:
                return {}
            text = zf.read(member).decode("utf-8", errors="replace")
        except zipfile.BadZipFile:
            return {}
    else:
        text = data.decode("utf-8", errors="replace")
    return _snps_from_lines(text.splitlines(), rsids)


def _acc(pairs):
    tp = sum(1 for l, p in pairs if l == "brown" and p == "brown")
    fn = sum(1 for l, p in pairs if l == "brown" and p == "blue")
    tn = sum(1 for l, p in pairs if l == "blue" and p == "blue")
    fp = sum(1 for l, p in pairs if l == "blue" and p == "brown")
    n = tp + fn + tn + fp
    return {"n_binary": n, "TP": tp, "FP": fp, "TN": tn, "FN": fn,
            "accuracy": round((tp + tn) / n, 3) if n else None}


def run(limit: int = 60, survey_bytes: bytes | None = None) -> dict:
    try:
        sb = survey_bytes if survey_bytes is not None else _get(SURVEY_URL)
    except Exception as e:
        return {"status": "SURVEY_FETCH_FAILED", "error": str(e)}
    labels = survey_eye_labels(sb)
    targets = [(hu, lab) for hu, lab in labels.items() if lab in ("blue", "brown")]

    want = set(_SNP_ORDER)
    v0_pairs, v01_pairs = [], []
    n_nofile = n_nors = n_incomplete = attempted = 0
    for hu, label in targets:
        if attempted >= limit:
            break
        attempted += 1
        link = profile_genotype_file(hu)
        time.sleep(0.5)  # polite between profile fetches
        if not link:
            n_nofile += 1
            continue
        fid, kind = link
        try:
            gts = fetch_snps(fid, kind, want)
        except Exception:
            n_nofile += 1
            continue
        v0_gt = gts.get(EYE_COLOUR_LOCUS)
        if not v0_gt:
            n_nors += 1
            continue
        v0_call = call_eye_colour(v0_gt)["prediction"]
        if v0_call in ("blue", "brown"):
            v0_pairs.append((label, v0_call))
        v01 = predict_irisplex(gts)
        if v01["status"] == "PREDICTED":
            if v01["prediction"] in ("blue", "brown"):
                v01_pairs.append((label, v01["prediction"]))
        else:
            n_incomplete += 1

    return {
        "status": "SCORED" if (v0_pairs or v01_pairs) else "NO_USERS_SCORED",
        "schema": "eye-colour-pgp-replication-v1", "date": _date.today().isoformat(),
        "cohort": "Personal Genome Project (PGP-Harvard); Basic Phenotypes 2015 survey; CC0 open-consent",
        "label_tier": "self-reported (independent 2nd cohort)",
        "is_pilot": True, "limit": limit, "n_profiles_attempted": attempted,
        "n_eye_labelled_total": len(labels), "n_blue_brown_available": len(targets),
        "n_no_genotype_file": n_nofile, "n_rs12913832_missing": n_nors, "n_incomplete_6snp": n_incomplete,
        "v0_rs12913832": _acc(v0_pairs),
        "v01_irisplex": _acc(v01_pairs),
        "replication_verdict": "REPLICATES" if (_acc(v0_pairs)["n_binary"] >= 20
                               and (_acc(v0_pairs)["accuracy"] or 0) >= 0.90) else "PILOT_UNDERPOWERED",
        "caveats": [
            "self-reported label (2nd independent cohort; non-circular)",
            "BOUNDED PILOT (--limit) -> replication signal, not a full-cohort number; full run = raise --limit",
            "complete-case for v0.1 (all 6 SNPs); v0 needs only rs12913832",
            "PGP is also Euro-majority -> same ancestry caveat as OpenSNP",
        ],
    }


def main(argv=None) -> int:
    import argparse
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--limit", type=int, default=60)
    a = ap.parse_args(argv)
    res = run(limit=a.limit)
    if res.get("status") == "SCORED":
        out = REPO / "wiki" / f"eye_colour_pgp_validation_{_date.today().isoformat()}.json"
        out.write_text(json.dumps(res, indent=2), encoding="utf-8")
        print(f"[wrote {out}]")
    print(json.dumps(res, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
