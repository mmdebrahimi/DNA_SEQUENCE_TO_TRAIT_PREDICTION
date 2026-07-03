"""Generalized single-SNP openSNP validator — scores any dna_decode.data.single_snp_traits cell.

Extends the eye-colour openSNP pattern to earwax (ABCC11) / lactase (LCT) / cilantro (OR6A2), REUSING the
zip-streaming machinery (`eye_colour_opensnp_ingest` helpers: phenotype-member finder, per-user genotype
members, member picker, one-pass rsid extractor). No new download — the archived 21 GB dump on D:.

Contract (mirrors the eye-colour cell): self-reported label (non-circular, noisy) -> PILOT/DEMO tier. Per
trait: bin the free-text phenotype, apply the deterministic strand-agnostic rule, build a confusion matrix
on the trait's positive class. INDETERMINATE calls + unbinnable labels are reported, never guessed.
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

from dna_decode.data.single_snp_traits import TRAITS, SingleSnpTrait  # noqa: E402
from scripts.eye_colour_opensnp_ingest import (  # noqa: E402
    _find_phenotype_member, _genotype_members_by_user, _pick_member, _rs_from_member,
)

DEFAULT_ZIP = Path("D:/dna_decode_cache/opensnp/opensnp_datadump.2017-12-08.zip")


def _pheno_rows(zf: zipfile.ZipFile, member: str):
    raw = zf.read(member).decode("utf-8", errors="replace")
    delim = ";" if raw.split("\n", 1)[0].count(";") >= raw.split("\n", 1)[0].count(",") else ","
    rdr = csv.reader(io.StringIO(raw), delimiter=delim)
    header = [h.strip() for h in next(rdr, [])]
    return header, list(rdr)


def _label_by_user(header, rows, trait: SingleSnpTrait) -> dict[str, str]:
    """Map user_id -> binned label for the RICHEST matching phenotype column (most binnable values)."""
    hl = [h.lower() for h in header]
    uid_i = next((i for i, h in enumerate(hl) if h in ("user_id", "id", "user")), 0)
    cand_cols = [i for i, h in enumerate(hl) if any(k in h for k in trait.phenotype_keywords)]
    best: dict[str, str] = {}
    for ci in cand_cols:
        got: dict[str, str] = {}
        for r in rows:
            if len(r) <= max(uid_i, ci):
                continue
            uid = r[uid_i].strip()
            lab = trait.binner(r[ci].strip())
            if uid and lab is not None:
                got[uid] = lab
        if len(got) > len(best):
            best = got
    return best


def run(zip_path: Path, trait_key: str, limit: int | None = None) -> dict:
    trait = TRAITS[trait_key]
    if not zip_path.exists():
        return {"status": "ZIP_NOT_PRESENT", "trait": trait_key, "zip": str(zip_path)}
    zf = zipfile.ZipFile(str(zip_path))
    pheno_member = _find_phenotype_member(zf)
    if not pheno_member:
        return {"status": "NO_PHENOTYPE_CSV", "trait": trait_key}
    header, rows = _pheno_rows(zf, pheno_member)
    labels = _label_by_user(header, rows, trait)
    geno_by_uid = _genotype_members_by_user(zf)

    pos, neg = trait.positive_label, trait.negative_label
    tp = fp = tn = fn = 0
    n_indet = n_nogeno = n_nors = scored = 0
    for uid, label in labels.items():
        member = _pick_member(geno_by_uid.get(uid, []))
        if not member:
            n_nogeno += 1
            continue
        gt = _rs_from_member(zf, member, rsid=trait.rsid)
        if not gt:
            n_nors += 1
            continue
        pred = trait.call(gt)
        if pred not in (pos, neg):
            n_indet += 1
            continue
        if label == pos and pred == pos:
            tp += 1
        elif label == neg and pred == pos:
            fp += 1
        elif label == neg and pred == neg:
            tn += 1
        else:
            fn += 1
        scored += 1
        if limit and scored >= limit:
            break
    n = tp + fp + tn + fn
    return {
        "status": "SCORED" if n else "NO_PAIRS",
        "schema": "single-snp-opensnp-v1", "date": _date.today().isoformat(),
        "trait": trait_key, "rsid": trait.rsid, "gene": trait.gene, "tier": trait.tier,
        "rule_source": trait.source,
        "source": "OpenSNP archive dump (archive.org/opensnp_data_dumps, 2017-12-08)",
        "label_tier": "self-reported (near-independent, non-circular, noisy) — PILOT/DEMO",
        "positive_class": pos, "negative_class": neg,
        "n_users_labelled": len(labels), "n_scored": n,
        "n_indeterminate_call": n_indet, "n_no_genotype_file": n_nogeno, "n_rsid_missing": n_nors,
        "confusion_positive_%s" % pos: {"TP": tp, "FP": fp, "TN": tn, "FN": fn},
        "accuracy": round((tp + tn) / n, 3) if n else None,
        "pos_sens": round(tp / (tp + fn), 3) if (tp + fn) else None,
        "neg_spec": round(tn / (tn + fp), 3) if (tn + fp) else None,
        "caveats": [
            "self-reported label (not a lab assay) — PILOT/DEMO tier, like the eye-colour openSNP cell",
            "deployed textbook rule -> a WIN validates the rule, it is not new biological signal",
            ("WEAK-EFFECT association by design (not Mendelian) -> near-chance is the CORRECT, honest outcome"
             if trait.tier == "WEAK_ASSOCIATION_CONTRAST"
             else "ancestry-confounded where the derived allele is population-specific (e.g. lactase European "
                  "-13910*T); within-ancestry split deferred"),
        ],
    }


def main(argv=None) -> int:
    import argparse
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--trait", choices=[*TRAITS, "all"], default="all")
    ap.add_argument("--zip", type=Path, default=DEFAULT_ZIP)
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--out-dir", type=Path, default=REPO / "wiki")
    a = ap.parse_args(argv)
    keys = list(TRAITS) if a.trait == "all" else [a.trait]
    rc = 0
    for k in keys:
        res = run(a.zip, k, limit=a.limit)
        print(json.dumps(res, indent=2))
        if res.get("status") == "SCORED":
            out = a.out_dir / f"single_snp_{k}_opensnp_{_date.today().isoformat()}.json"
            out.write_text(json.dumps(res, indent=2), encoding="utf-8")
            print(f"[wrote {out}]")
        else:
            rc = 1
        print("=" * 70)
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
