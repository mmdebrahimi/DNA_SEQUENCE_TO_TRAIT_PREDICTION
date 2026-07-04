"""Capture 1000 Genomes ancestry-stratified pairwise LD (Ensembl REST) — validates the imputation layer's
LD basis on an INDEPENDENT multi-population panel + quantifies its ancestry limit.

1000G is NOT a trait substrate (phenotype = ancestry/sex/family only) — its value here is as the LD /
ancestry REFERENCE the imputation pre-processor needs. The openSNP-derived imputation maps are
European-dominated; this queries 1000G r²/D' per super-population for a determinant↔tag pair so the map can
carry an INDEPENDENT, ancestry-stratified LD annotation (valid where r² >= threshold, abstain-worthy where
low). Lightweight: Ensembl REST, no bulk VCF download (C: is disk-tight). Free, no auth.
"""
from __future__ import annotations

import json
import sys
import urllib.request
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
_ENSEMBL = "https://rest.ensembl.org"
SUPERPOPS = ["EUR", "AFR", "EAS", "SAS", "AMR"]
VALID_R2 = 0.90                                       # r2 >= this => the tag is a reliable proxy in that pop


def _rest(path: str) -> list | dict | None:
    req = urllib.request.Request(_ENSEMBL + path, headers={"Content-Type": "application/json"})
    try:
        with urllib.request.urlopen(req, timeout=30) as r:
            return json.loads(r.read().decode("utf-8"))
    except Exception:
        return None


def pairwise_ld(tag: str, target: str) -> dict:
    out = {}
    for sp in SUPERPOPS:
        d = _rest(f"/ld/human/pairwise/{tag}/{target}?population_name=1000GENOMES:phase_3:{sp}")
        if isinstance(d, list) and d:
            out[sp] = {"r2": round(float(d[0]["r2"]), 3), "d_prime": round(float(d[0].get("d_prime", 0)), 3)}
        else:
            out[sp] = {"r2": None, "d_prime": None, "note": "absent/monomorphic in this population"}
    valid = [sp for sp, v in out.items() if isinstance(v.get("r2"), float) and v["r2"] >= VALID_R2]
    return {"schema": "1000g-ld-annotation-v1", "tag": tag, "target": target,
            "source": "Ensembl REST /ld/human/pairwise (1000 Genomes phase 3 super-populations)",
            "valid_r2_threshold": VALID_R2, "by_superpopulation": out,
            "valid_populations": valid,
            "ancestry_caveat": (f"tag is a reliable proxy (r2>={VALID_R2}) in {valid or 'NO'} population(s); "
                                "imputation from this tag is UNRELIABLE elsewhere — a map derived in one "
                                "ancestry must NOT be applied across ancestries without this check")}


def main(argv=None) -> int:
    import argparse
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--tag", default="rs657152")
    ap.add_argument("--target", default="rs8176719")
    ap.add_argument("--out", type=Path, default=None)
    a = ap.parse_args(argv)
    res = pairwise_ld(a.tag, a.target)
    if res["by_superpopulation"].get("EUR", {}).get("r2") is None and not res["valid_populations"]:
        print("WARN: no 1000G LD returned (REST unreachable?) — nothing written", file=sys.stderr)
    out = a.out or (REPO / "data" / "imputation" / f"{a.target}_from_{a.tag}_1000g_ld.json")
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(res, indent=2), encoding="utf-8")
    print(json.dumps(res, indent=2))
    print(f"\n[wrote {out}]")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
