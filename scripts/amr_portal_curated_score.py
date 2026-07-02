"""Generic curated-cell scorer for the AMR-Portal Tier-3/4 backlog (the reusable workhorse).

Replaces the per-cell scorer proliferation (neisseria_cipro / staph_cipro / staph_rif were bespoke). A
curated cell = (organism, drug, rule_fn) where `rule_fn(determinant_dicts) -> {"prediction": "R"/"S",
"determinant_present": bool}`. This module loads the AMR-Portal provenance-disjoint measured-AST isolates,
applies the rule, computes binary sens/spec + Wilson CIs + the determinant-present/absent strata + the
spec>=0.85 endorsement falsifier, and writes a NAMESPACE-SEPARATE `wiki/amr_portal_<slug>_<date>.json`
(CURATED_NONFROZEN, not the frozen deployed surface). Cells register in `CURATED_CELLS`.
"""
from __future__ import annotations

import json
import sys
from datetime import date as _date
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from scripts.amr_portal_score_independent import (  # noqa: E402
    DEFAULT_CRYPTIC, DEFAULT_PARQUET, GENO_PARQUET,
    _load_genotype, _load_leak_set, _load_phenotype_cells, wilson_ci,
)
from scripts.independent_cohort_validate import _conf  # noqa: E402

SPEC_FLOOR = 0.85
MIN_PER_CLASS = 10


def score_cell(isolates, geno, rule_fn) -> dict:
    """isolates=[(biosample, leaked, SIR)]; geno={bs:[determinant dicts]}; rule_fn(dets)->{prediction, determinant_present}."""
    strata = {"determinant_present": {"R": 0, "S": 0}, "determinant_absent": {"R": 0, "S": 0}}
    pairs = []
    for bs, leaked, sir in isolates:
        if leaked:
            continue
        r = rule_fn(geno.get(bs, []))
        pairs.append((r["prediction"], 1 if sir == "R" else 0))          # _conf wants int label
        k = "determinant_present" if r["determinant_present"] else "determinant_absent"
        strata[k]["R" if sir == "R" else "S"] += 1
    conf = _conf(pairs)
    tbl = {}
    for k, v in strata.items():
        n = v["R"] + v["S"]
        tbl[k] = {"n": n, "R": v["R"], "S": v["S"], "r_rate": round(v["R"] / n, 3) if n else None}
    n_R, n_S = conf["tp"] + conf["fn"], conf["tn"] + conf["fp"]
    spec = conf["spec"]
    pres, absnt = tbl["determinant_present"]["r_rate"], tbl["determinant_absent"]["r_rate"]
    strata_ok = bool(pres is not None and absnt is not None and pres > absnt)
    powered = n_R >= MIN_PER_CLASS and n_S >= MIN_PER_CLASS
    endorsed = bool(powered and spec is not None and spec >= SPEC_FLOOR and strata_ok)
    return {"binary": conf, "n_R": n_R, "n_S": n_S, "powered": powered,
            "sens_wilson95": wilson_ci(conf["tp"], n_R), "spec_wilson95": wilson_ci(conf["tn"], n_S),
            "strata": tbl, "strata_reproduced": strata_ok, "spec_floor": SPEC_FLOOR,
            "headline": "SCORED" if endorsed else ("UNDERPOWERED" if not powered else "INDETERMINATE")}


def run_cell(organism, drug, rule_fn, slug, rule_text, *,
             pheno=DEFAULT_PARQUET, geno=GENO_PARQUET, cryptic=DEFAULT_CRYPTIC, write=True) -> dict:
    leak = _load_leak_set(cryptic)
    g = _load_genotype(geno, {organism})
    p = _load_phenotype_cells(pheno, [(organism, drug)], leak)
    r = score_cell(p.get((organism, drug), []), g, rule_fn)
    out = {"_schema": "amr-portal-curated-cell-v1", "date": _date.today().isoformat(),
           "organism": organism, "drug": drug, "rule_status": "CURATED_NONFROZEN",
           "rule_scope": "scorer_local", "not_in_shipped_surface": True, "rule_text": rule_text,
           "independence_tier": "EBI AMR Portal provenance-disjoint; AMRFinderPlus genotype x measured AST",
           "endorsement_gate": f"spec >= {SPEC_FLOOR} on powered (>=10 R/S) + strata reproduce",
           "frozen_surface_changed": False, **r}
    if write:
        (REPO / "wiki" / f"amr_portal_{slug}_{_date.today().isoformat()}.json").write_text(
            json.dumps(out, indent=2, default=str), encoding="utf-8")
    b = r["binary"]
    print(f"{organism} {drug}: nR={r['n_R']} nS={r['n_S']} acc={b['acc']} sens={b['sens']} spec={b['spec']} "
          f"repro={r['strata_reproduced']} -> {r['headline']}")
    return out


# ---- curated cell registry ---------------------------------------------------------------------------------
def _syms(dets):
    return [(d.get("amr_element_symbol") or d.get("gene_symbol") or "") for d in dets]


def _ng_tet(dets):
    from dna_decode.organism_rules.neisseria_amr import call_ng_tetracycline
    c = call_ng_tetracycline(_syms(dets))
    return {"prediction": c["prediction"], "determinant_present": bool(c["matched_tetM"])}


CURATED_CELLS = {
    "neisseria_tet": ("Neisseria gonorrhoeae", "tetracycline", _ng_tet, "neisseria_tet",
                      "tet(M) -> R (rpsJ V57M accessory-only)"),
}


def main(argv=None) -> int:
    import argparse
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("cell", choices=sorted(CURATED_CELLS))
    a = ap.parse_args(argv)
    org, drug, fn, slug, txt = CURATED_CELLS[a.cell]
    run_cell(org, drug, fn, slug, txt)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
