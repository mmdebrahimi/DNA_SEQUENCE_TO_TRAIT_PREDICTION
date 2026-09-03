"""Does the deployed `rmt` rescue miss real carriers of the wider 16S-RMTase family?

THE MECHANISM, sourced not recalled. Acquired 16S rRNA methyltransferases fall in two families:
  - m7G1405: ArmA, RmtA..RmtI -- confer resistance to 4,6-disubstituted deoxystreptamines, a class that
    INCLUDES gentamicin.
  - m1A1408: NpmA, NpmB, NpmC, WarA -- broader spectrum, also covering gentamicin.
So EVERY member of both families confers gentamicin resistance. The deployed v2 rescue is
`^(rmt[A-H]\\d*|npmA\\d*)$`, whose character range stops one letter short of `rmtI` and whose npm branch
covers only `npmA`. `armA` needs no rescue -- AMRFinder files it under Subclass GENTAMICIN already.

So the rule under-covers by construction: `rmtI`, `npmB`, `npmC`, `warA`. That is a SAFE-DIRECTION gap --
it can only cost sensitivity, never specificity -- and it does NOT bear on the untested over-call risk.
This probe measures whether the gap is REAL in the data or merely nomenclatural.

Counts carriers of the wider family in NCBI-PD's own `AMR_genotypes` field, split by whether the deployed
rescue would have caught them, with the measured gentamicin phenotype where one exists.

Network-only; writes wiki/rmt_family_coverage_probe.json.
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from dna_decode.data.pd_ast import ast_label_for  # noqa: E402
from gentamicin_rmt_specificity_hunt import (  # noqa: E402
    RESCUE_RE, latest_metadata_url, parse_amr_genotypes,
)
import urllib.request  # noqa: E402

# The full acquired 16S-RMTase family. `armA` is deliberately EXCLUDED: it is already counted by the
# frozen Subclass=GENTAMICIN rule, so including it here would double-count and misattribute coverage.
FAMILY_RE = re.compile(r"^(rmt[A-Z]\d*|npm[A-Z]\d*|warA\d*)$")

GROUPS = ("Klebsiella", "Escherichia_coli_Shigella", "Acinetobacter",
          "Enterobacter_hormaechei", "Pseudomonas_aeruginosa", "Salmonella")


def scan(group: str, drug: str) -> dict:
    r = urllib.request.urlopen(latest_metadata_url(group), timeout=600)
    cols = r.readline().decode("utf8", "replace").rstrip("\n").split("\t")
    idx = {c: i for i, c in enumerate(cols)}
    gi, ai = idx["AMR_genotypes"], idx["AST_phenotypes"]

    seen: dict[str, dict] = {}
    n_rows = 0
    for line in r:
        n_rows += 1
        f = line.decode("utf8", "replace").rstrip("\n").split("\t")
        if len(f) <= max(gi, ai):
            continue
        for sym in parse_amr_genotypes(f[gi]):
            if not FAMILY_RE.match(sym):
                continue
            rec = seen.setdefault(sym, {"symbol": sym, "rescued_by_deployed_rule": bool(RESCUE_RE.match(sym)),
                                        "n_genomes": 0, "labelled": {"R": 0, "S": 0, "I": 0}})
            rec["n_genomes"] += 1
            lab = ast_label_for(f[ai], drug)
            if lab in ("R", "S", "I"):
                rec["labelled"][lab] += 1
    return {"group": group, "n_rows": n_rows, "symbols": sorted(seen.values(), key=lambda x: -x["n_genomes"])}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--groups", default=",".join(GROUPS))
    ap.add_argument("--drug", default="gentamicin")
    ap.add_argument("--out", type=Path, default=ROOT / "wiki" / "rmt_family_coverage_probe.json")
    a = ap.parse_args()

    per_group, errors = [], {}
    for g in [x.strip() for x in a.groups.split(",") if x.strip()]:
        try:
            res = scan(g, a.drug)
        except Exception as e:
            errors[g] = f"{type(e).__name__}: {e}"
            print(f"{g:26} ERROR {errors[g][:60]}", flush=True)
            continue
        per_group.append(res)
        miss = [s for s in res["symbols"] if not s["rescued_by_deployed_rule"]]
        print(f"{g:26} family symbols={len(res['symbols']):2d}  "
              f"NOT rescued: {[s['symbol'] for s in miss] or 'none'}", flush=True)

    agg: dict[str, dict] = {}
    for gres in per_group:
        for s in gres["symbols"]:
            a_ = agg.setdefault(s["symbol"], {"symbol": s["symbol"],
                                              "rescued_by_deployed_rule": s["rescued_by_deployed_rule"],
                                              "n_genomes": 0, "labelled": {"R": 0, "S": 0, "I": 0}})
            a_["n_genomes"] += s["n_genomes"]
            for k in ("R", "S", "I"):
                a_["labelled"][k] += s["labelled"][k]

    missed = {k: v for k, v in agg.items() if not v["rescued_by_deployed_rule"]}
    n_missed_genomes = sum(v["n_genomes"] for v in missed.values())
    n_missed_labelled = sum(sum(v["labelled"].values()) for v in missed.values())

    print(f"\nfamily symbols seen: {sorted(agg)}")
    print(f"NOT rescued by the deployed rule: {sorted(missed) or 'NONE'}")
    print(f"  genomes carrying an unrescued family member : {n_missed_genomes}")
    print(f"  of those, with a gentamicin label           : {n_missed_labelled}")
    for v in sorted(missed.values(), key=lambda x: -x["n_genomes"]):
        print(f"    {v['symbol']:8} genomes={v['n_genomes']:5d}  labelled R/S/I="
              f"{v['labelled']['R']}/{v['labelled']['S']}/{v['labelled']['I']}")

    out = {"schema": "rmt-family-coverage-probe-v1", "drug": a.drug,
           "deployed_rescue": RESCUE_RE.pattern, "family_pattern": FAMILY_RE.pattern,
           "note": ("armA is deliberately excluded from the family pattern -- AMRFinder files it under "
                    "Subclass GENTAMICIN, so the frozen rule already counts it and including it would "
                    "misattribute coverage."),
           "per_group": per_group, "aggregate": sorted(agg.values(), key=lambda x: -x["n_genomes"]),
           "unrescued": sorted(missed.values(), key=lambda x: -x["n_genomes"]),
           "n_unrescued_genomes": n_missed_genomes, "n_unrescued_labelled": n_missed_labelled,
           "errors": errors, "complete": not errors,
           "honest_limits": [
               "This measures COVERAGE (does the rule see the determinant), not correctness. Every "
               "member of both RMTase families confers gentamicin resistance by mechanism, so a missed "
               "carrier is a SENSITIVITY loss -- it can never cause an over-call.",
               "It therefore does NOT bear on the deployed rule's untested over-call/specificity risk.",
               "Carrier calls are NCBI's own AMRFinder (PD AMR_genotypes), a tool-derived feature.",
               "A symbol absent here may still exist in nature; absence in PD is not absence in the world.",
           ]}
    a.out.write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")
    print(f"\ncomplete={not errors}  wrote {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
