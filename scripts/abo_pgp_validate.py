"""M4: ABO blood-type (O vs non-O) decoder validated on PGP — the serological Mendelian cell.

Reuses the M3 PGP harness (survey 19 + profile->genotype-file linking, files cached on D:) but for the
BLOOD TYPE label (survey col '1.1 - Blood Type') + the rs8176719 O-status rule. Deterministic curated rule
(c.261delG deletion -> O) x free independent self-report label. Same architecture as eye colour, a different
trait class (serological). BOUNDED PILOT (--limit); rate-limited; genotype files reuse the M3 D: cache.

Emits wiki/abo_pgp_validation_<date>.json.
"""
from __future__ import annotations

import csv
import io
import json
import sys
import time
from datetime import date as _date
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from dna_decode.data.abo_blood import bin_blood_type, call_abo_o_status  # noqa: E402
from scripts.eye_colour_pgp_validate import (  # noqa: E402
    SURVEY_URL, _get, fetch_snps, profile_genotype_file,
)

RS = "rs8176719"


def blood_labels(csv_bytes: bytes) -> dict[str, str]:
    rows = list(csv.reader(io.StringIO(csv_bytes.decode("utf-8", errors="replace"))))
    hdr = rows[0]
    bi = next(i for i, h in enumerate(hdr) if "Blood Type" in h)
    out: dict[str, str] = {}
    for r in rows[1:]:
        if len(r) <= bi:
            continue
        hu = r[0].strip()
        lab = bin_blood_type(r[bi].strip())
        if hu and lab:
            out[hu] = lab
    return out


def run(limit: int = 200, survey_bytes: bytes | None = None) -> dict:
    try:
        sb = survey_bytes if survey_bytes is not None else _get(SURVEY_URL)
    except Exception as e:
        return {"status": "SURVEY_FETCH_FAILED", "error": str(e)}
    labels = blood_labels(sb)
    targets = list(labels.items())

    pairs = []  # (label, pred) over O/non-O
    n_nofile = n_nors = attempted = 0
    for hu, label in targets:
        if attempted >= limit:
            break
        attempted += 1
        link = profile_genotype_file(hu)
        time.sleep(0.5)
        if not link:
            n_nofile += 1
            continue
        fid, kind = link
        try:
            gts = fetch_snps(fid, kind, {RS})
        except Exception:
            n_nofile += 1
            continue
        gt = gts.get(RS)
        pred = call_abo_o_status(gt) if gt else "INDETERMINATE"
        if pred in ("O", "non-O"):
            pairs.append((label, pred))
        else:
            n_nors += 1

    # O = positive class
    tp = sum(1 for l, p in pairs if l == "O" and p == "O")
    fn = sum(1 for l, p in pairs if l == "O" and p == "non-O")
    tn = sum(1 for l, p in pairs if l == "non-O" and p == "non-O")
    fp = sum(1 for l, p in pairs if l == "non-O" and p == "O")
    n = tp + fn + tn + fp
    return {
        "status": "SCORED" if n else "NO_USERS_SCORED",
        "schema": "abo-pgp-o-status-v1", "date": _date.today().isoformat(),
        "cohort": "PGP-Harvard Basic Phenotypes 2015; CC0 open-consent",
        "rule": "rs8176719 c.261delG homozygous deletion (DD) -> O, else non-O (SOURCED; O/non-O only)",
        "label_tier": "self-reported blood type (independent)",
        "is_pilot": True, "limit": limit, "n_profiles_attempted": attempted,
        "n_blood_labelled_total": len(labels),
        "n_no_genotype_file": n_nofile, "n_rs8176719_missing_or_nocall": n_nors,
        "confusion_O_positive": {"TP": tp, "FP": fp, "TN": tn, "FN": fn, "n": n},
        "accuracy": round((tp + tn) / n, 3) if n else None,
        "O_sens": round(tp / (tp + fn), 3) if (tp + fn) else None,
        "nonO_spec": round(tn / (tn + fp), 3) if (tn + fp) else None,
        "caveats": [
            "O vs non-O only; A/B distinction deferred (tag-SNP allele coding not fabricated)",
            "self-reported blood type (independent, noisy) -- some self-reports are wrong",
            "BOUNDED PILOT (--limit); low PGP survey->genotype linkage yield (~5%, as M3)",
            "rs8176719 I/D coding; '--' no-call -> INDETERMINATE (never guessed O)",
        ],
    }


def main(argv=None) -> int:
    import argparse
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--limit", type=int, default=200)
    a = ap.parse_args(argv)
    res = run(limit=a.limit)
    if res.get("status") == "SCORED":
        out = REPO / "wiki" / f"abo_pgp_validation_{_date.today().isoformat()}.json"
        out.write_text(json.dumps(res, indent=2), encoding="utf-8")
        print(f"[wrote {out}]")
    print(json.dumps(res, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
