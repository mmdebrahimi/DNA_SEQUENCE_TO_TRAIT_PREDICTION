"""S. aureus rifampicin curated rule scored on the EBI AMR Portal (Tier-3/4 cell 3, provenance-disjoint).

Applies the non-frozen rpoB-RRDR rule (`organism_rules/staphylococcus_amr.call_sa_rifampicin`) to the AMR
Portal's provenance-disjoint measured-AST S. aureus rifampicin isolates. NAMESPACE-SEPARATE artifact
(`wiki/amr_portal_staphylococcus_rif_<date>.json`), CURATED_NONFROZEN. Endorsement falsifier: spec >= 0.85 on
the powered (>=10 R/S) provenance-disjoint set AND rpoB-present R-rate > rpoB-absent R-rate.
"""
from __future__ import annotations

import json
import sys
from datetime import date as _date
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from dna_decode.organism_rules.staphylococcus_amr import DRUG_RIF, ORGANISM, call_sa_rifampicin  # noqa: E402
from scripts.amr_portal_score_independent import (  # noqa: E402
    DEFAULT_CRYPTIC, DEFAULT_PARQUET, GENO_PARQUET,
    _load_genotype, _load_leak_set, _load_phenotype_cells, wilson_ci,
)
from scripts.independent_cohort_validate import _conf  # noqa: E402

SPEC_FLOOR = 0.85
MIN_PER_CLASS = 10


def _symbols(dets):
    return [(d.get("amr_element_symbol") or d.get("gene_symbol") or "") for d in dets]


def score(isolates, geno) -> dict:
    strata = {"rpoB_present": {"R": 0, "S": 0}, "rpoB_absent": {"R": 0, "S": 0}}
    pairs = []
    for bs, leaked, sir in isolates:
        if leaked:
            continue
        call = call_sa_rifampicin(_symbols(geno.get(bs, [])))
        pairs.append((call["prediction"], 1 if sir == "R" else 0))
        k = "rpoB_present" if call["matched_rpoB"] else "rpoB_absent"
        strata[k]["R" if sir == "R" else "S"] += 1
    conf = _conf(pairs)
    tbl = {}
    for k, v in strata.items():
        n = v["R"] + v["S"]
        tbl[k] = {"n": n, "R": v["R"], "S": v["S"], "r_rate": round(v["R"] / n, 3) if n else None}
    n_R, n_S = conf["tp"] + conf["fn"], conf["tn"] + conf["fp"]
    spec = conf["spec"]
    pres, absnt = tbl["rpoB_present"]["r_rate"], tbl["rpoB_absent"]["r_rate"]
    strata_ok = bool(pres is not None and absnt is not None and pres > absnt)
    powered = n_R >= MIN_PER_CLASS and n_S >= MIN_PER_CLASS
    endorsed = bool(powered and spec is not None and spec >= SPEC_FLOOR and strata_ok)
    return {"binary": conf, "n_R": n_R, "n_S": n_S, "powered": powered,
            "sens_wilson95": wilson_ci(conf["tp"], n_R), "spec_wilson95": wilson_ci(conf["tn"], n_S),
            "strata": tbl, "strata_reproduced": strata_ok, "spec_floor": SPEC_FLOOR,
            "headline": "SCORED" if endorsed else ("UNDERPOWERED" if not powered else "INDETERMINATE")}


def main(argv=None) -> int:
    import argparse
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--pheno", type=Path, default=DEFAULT_PARQUET)
    ap.add_argument("--geno", type=Path, default=GENO_PARQUET)
    ap.add_argument("--cryptic", type=Path, default=DEFAULT_CRYPTIC)
    a = ap.parse_args(argv)
    for p in (a.pheno, a.geno):
        if not p.exists():
            print(f"ERROR: AMR Portal parquet not found: {p}", file=sys.stderr)
            return 2
    leak = _load_leak_set(a.cryptic)
    geno = _load_genotype(a.geno, {ORGANISM})
    pheno = _load_phenotype_cells(a.pheno, [(ORGANISM, DRUG_RIF)], leak)
    r = score(pheno.get((ORGANISM, DRUG_RIF), []), geno)
    out = {
        "_schema": "amr-portal-staphylococcus-rif-v1", "date": _date.today().isoformat(),
        "organism": ORGANISM, "drug": DRUG_RIF,
        "rule_status": "CURATED_NONFROZEN", "rule_scope": "scorer_local", "not_in_shipped_surface": True,
        "rule_text": call_sa_rifampicin([])["rule"],
        "independence_tier": ("EBI AMR Portal provenance-disjoint; genotype = the Portal's AMRFinderPlus POINT "
                              "calls; phenotype = wet-lab measured AST (binary R/S, non-circular)"),
        "endorsement_gate": f"spec >= {SPEC_FLOOR} on powered (>=10 R/S) provenance-disjoint + strata reproduce",
        "note_rare_R": "rifampicin resistance is rare in S. aureus -> R class small; powered but S-skewed",
        "frozen_surface_changed": False, **r,
    }
    outp = REPO / "wiki" / f"amr_portal_staphylococcus_rif_{_date.today().isoformat()}.json"
    outp.write_text(json.dumps(out, indent=2, default=str), encoding="utf-8")
    b = r["binary"]
    print(f"SA rif: nR={r['n_R']} nS={r['n_S']} acc={b['acc']} sens={b['sens']} spec={b['spec']} "
          f"strata_repro={r['strata_reproduced']} -> {r['headline']}")
    print(f"  strata: {r['strata']}")
    print(f"  -> {outp}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
