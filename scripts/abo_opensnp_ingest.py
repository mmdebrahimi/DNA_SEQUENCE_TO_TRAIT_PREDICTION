"""Stage the J3-ABO substrate from the OpenSNP 2017-12-08 dump — per-sample ABO genotype + self-report label.

This is the CPU-only, no-fabrication INPUT that BOTH J3 embedding arms (frozen_fm + learned_rep) consume,
plus the de-confoundable label for the J3 downstream falsifier (`j3-abo-embedding-v1`). It does NOT compute
embeddings (that is the workhorse-GPU step) — it produces the verified substrate + a manifest so the
workhorse emits embeddings aligned to a clean, deterministic sample index.

Reuses the eye-colour OpenSNP ingest machinery (selective zip streaming — NEVER extracts the 21 GB dump):
`_find_phenotype_member` / `_genotype_members_by_user` / `_pick_member` / `rsids_from_member`.

ABO locus (chr9q34.2, ABO gene):
  * rs8176719  c.261delG  — the O-allele frameshift DELETION. 23andMe reports D (deletion) / I (insertion).
      DD = homozygous deletion → blood type O (deterministic).  DI / II → at least one functional allele → non-O.
  * rs8176746  (A/B-distinguishing, p.Leu266Met)  — A vs B tag.
  * rs8176747  (A/B-distinguishing, p.Gly268Ala)  — A vs B tag.

HONESTY: self-reported ABO label (near-independent, non-circular, NOISY — self-report ~15% erroneous +
non-deletional O alleles rs8176719 misses). The deterministic O-vs-non-O call is the KNOWN-mechanism
baseline the learned rep must BEAT on a de-confounded slice; the label is the falsifier target, not truth.
"""
from __future__ import annotations

import argparse
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

from scripts.eye_colour_opensnp_ingest import (  # noqa: E402
    _find_phenotype_member, _genotype_members_by_user, _pick_member, rsids_from_member,
)

DEFAULT_ZIP = Path("D:/dna_decode_cache/opensnp/opensnp_datadump.2017-12-08.zip")
ABO_SNPS = {"rs8176719", "rs8176746", "rs8176747"}

# EXACT ABO-relevant phenotype headers (normalized lower/stripped) — NOT fuzzy substring (avoids
# metABOlite / asparagus / Kell-Kidd-Diego decoys). Verified against the real 2017 dump header.
_ABO_HEADERS = {
    "abh blood group (antigens)", "abo rh", "blood type", "blood type a1, rh positive",
}


def _norm(h: str) -> str:
    return " ".join(h.strip().lower().split())


def parse_abo_group(raw: str) -> str | None:
    """Self-reported blood-type string -> ABO group in {O, A, B, AB}, or None if no ABO info.

    Strips Rh (+/-, positive/negative/rhesus), maps written '0'->O, resolves genotype-style (AO->A, BO->B).
    Pure-Rh statements ('Rh negative') and junk -> None."""
    if not raw:
        return None
    s = raw.strip().lower()
    for junk in ("rather not say", "unknown", "n/a", "na", "none", "-", "can't", "cant "):
        if s == junk or s.startswith(junk):
            return None
    # drop Rh words + signs
    for tok in ("rhesus", "rh", "positive", "negative", "pos ", "neg ", "+", "-", "positiv", "negativ", "/"):
        s = s.replace(tok, " ")
    s = s.replace("0", "o")          # people write zero for O
    s = " ".join(s.split())
    letters = "".join(ch for ch in s.upper() if ch in "OAB")
    if not letters:
        return None
    st = set(letters)
    if "A" in st and "B" in st:      # AB (or genotype AB) -> AB
        return "AB"
    if "A" in st:                    # A, AA, AO -> A
        return "A"
    if "B" in st:                    # B, BB, BO -> B
        return "B"
    if "O" in st:                    # O, OO -> O
        return "O"
    return None


def deterministic_o_call(rs8176719_gt: str | None) -> str | None:
    """rs8176719 (D=deletion=O allele): DD -> O; DI/II -> non_O; missing/other -> None (abstain)."""
    if not rs8176719_gt:
        return None
    g = rs8176719_gt.upper().replace("-", "").replace(" ", "")
    # 23andMe reports I/D; some report the actual alleles. Map to D-count.
    if g in ("DD",):
        return "O"
    if g in ("DI", "ID", "II"):
        return "non_O"
    # allele-form fallback: a homozygous single-base (deletion often shown as '--' or 'DD') — abstain otherwise
    return None


def _abo_label_by_user(zf: zipfile.ZipFile, pheno_member: str) -> dict[str, str]:
    raw = zf.read(pheno_member).decode("utf-8", errors="replace")
    rows = list(csv.reader(io.StringIO(raw), delimiter=";"))
    if not rows:
        return {}
    hdr = rows[0]
    abo_cols = [i for i, h in enumerate(hdr) if _norm(h) in _ABO_HEADERS]
    out: dict[str, str] = {}
    for row in rows[1:]:
        uid = row[0].strip() if row else ""
        if not uid or uid in out:
            continue
        for ci in abo_cols:
            if len(row) > ci and parse_abo_group(row[ci]) is not None:
                out[uid] = row[ci].strip()
                break
    return out


def run(zip_path: Path, limit: int | None = None) -> dict:
    if not zip_path.exists():
        return {"status": "ZIP_NOT_PRESENT", "zip": str(zip_path)}
    zf = zipfile.ZipFile(str(zip_path))
    pheno_member = _find_phenotype_member(zf)
    if not pheno_member:
        return {"status": "NO_PHENOTYPE_CSV"}
    labels = _abo_label_by_user(zf, pheno_member)
    geno = _genotype_members_by_user(zf)

    samples = []
    n_no_geno = n_no_rs = n_concordant = n_discordant = 0
    for uid, raw_label in labels.items():
        member = _pick_member(geno.get(uid, []))
        if not member:
            n_no_geno += 1
            continue
        gts = rsids_from_member(zf, member, ABO_SNPS)
        if "rs8176719" not in gts:
            n_no_rs += 1
            continue
        abo = parse_abo_group(raw_label)
        det = deterministic_o_call(gts.get("rs8176719"))
        label_ovn = "O" if abo == "O" else ("non_O" if abo in ("A", "B", "AB") else None)
        concordant = (det is not None and label_ovn is not None and det == label_ovn)
        if det is not None and label_ovn is not None:
            n_concordant += int(concordant)
            n_discordant += int(not concordant)
        samples.append({
            "user_id": uid, "genotype_member": member,
            "raw_label": raw_label, "abo_group": abo, "label_o_vs_nonO": label_ovn,
            "rs8176719": gts.get("rs8176719"), "rs8176746": gts.get("rs8176746"),
            "rs8176747": gts.get("rs8176747"),
            "deterministic_o_call": det, "det_vs_label_concordant": concordant,
        })
        if limit and len(samples) >= limit:
            break
    n = len(samples)
    scored = n_concordant + n_discordant
    return {
        "status": "STAGED" if n else "NO_SAMPLES",
        "schema": "j3-abo-substrate-v1", "date": _date.today().isoformat(),
        "source": "OpenSNP archive dump 2017-12-08 (D:/dna_decode_cache/opensnp)",
        "locus": "ABO chr9q34.2; rs8176719 (O-deletion) + rs8176746/rs8176747 (A/B tags)",
        "label_tier": "self-reported ABO (near-independent, non-circular, ~15% noisy; non-deletional O missed)",
        "n_users_abo_labelled": len(labels), "n_substrate_samples": n,
        "n_no_genotype_file": n_no_geno, "n_rs8176719_missing": n_no_rs,
        "deterministic_baseline": {
            "n_scored_o_vs_nonO": scored, "n_concordant": n_concordant, "n_discordant": n_discordant,
            "concordance": round(n_concordant / scored, 3) if scored else None,
            "note": "rs8176719 DD->O vs self-report O-vs-nonO; this is the KNOWN-mechanism baseline the "
                    "learned_rep must BEAT on a de-confounded slice (else no embedding niche).",
        },
        "samples": samples,
    }


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--zip", type=Path, default=DEFAULT_ZIP)
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--out", type=Path, default=REPO / "data" / "j3_abo" / "j3_abo_substrate.json")
    ap.add_argument("--tsv", type=Path, default=REPO / "data" / "j3_abo" / "j3_abo_substrate.tsv")
    a = ap.parse_args(argv)
    res = run(a.zip, limit=a.limit)
    if res.get("status") == "STAGED":
        a.out.parent.mkdir(parents=True, exist_ok=True)
        a.out.write_text(json.dumps(res, indent=2), encoding="utf-8")
        cols = ["user_id", "raw_label", "abo_group", "label_o_vs_nonO", "rs8176719", "rs8176746",
                "rs8176747", "deterministic_o_call", "det_vs_label_concordant", "genotype_member"]
        with a.tsv.open("w", encoding="utf-8", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(cols)
            for s in res["samples"]:
                w.writerow([s[c] for c in cols])
        print(f"[wrote {a.out} + {a.tsv}]")
    summary = {k: v for k, v in res.items() if k != "samples"}
    print(json.dumps(summary, indent=2))
    return 0 if res.get("status") == "STAGED" else 1


if __name__ == "__main__":
    raise SystemExit(main())
