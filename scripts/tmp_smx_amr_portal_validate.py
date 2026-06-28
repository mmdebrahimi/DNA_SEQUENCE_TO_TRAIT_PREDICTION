"""TMP-SMX EXPERIMENTAL overlay scored on the EBI AMR Portal (provenance-disjoint, measured AST).

Tier-1 of the AMR-Portal unscored-cell triage (`wiki/amr_portal_unscored_triage_2026-06-28.md`): the
5 powered TMP-SMX cells (E. coli / Klebsiella / Salmonella / Shigella sonnei / Shigella flexneri) have a
REAL rule — the experimental `(>=1 acquired sul) AND (>=1 acquired dfr)` overlay (NON-frozen,
`dna_decode/data/experimental_drug_rules.py`), already validated on Oxford + Sci234. This scores it on the
AMR Portal's free measured-AST isolates.

NAMESPACE-SEPARATE by design (the shared-key silent-overwrite trap): this writes its OWN
`wiki/amr_portal_tmpsmx_experimental_<date>.json` artifact, branded EXPERIMENTAL_SCORED / scorer_local /
not_in_shipped_surface — it does NOT touch the frozen `amr_portal_independent_*` card (which is the
deployed-surface validation). The FROZEN AMR surface is untouched (this rule is not in it).

The AMR Portal gives BINARY measured AST (R/S), not MIC -> only the binary metric + the 4-genotype-strata
table + the strata-reproduction gate apply (no strict/relaxed MIC tiers). Reuses the portal scorer's
genotype/phenotype/leakage loaders verbatim; reuses the experimental rule's tmp_smx_call.
"""
from __future__ import annotations

import json
import sys
from datetime import date as _date
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from dna_decode.data.experimental_drug_rules import (  # noqa: E402
    DRUG, RULE_SCOPE, RULE_STATUS, tmp_smx_call,
)
from scripts.amr_portal_score_independent import (  # noqa: E402
    DEFAULT_CRYPTIC, DEFAULT_PARQUET, GENO_PARQUET,
    _load_genotype, _load_leak_set, _load_phenotype_cells, wilson_ci,
)
from scripts.independent_cohort_validate import _conf  # noqa: E402

# the 5 TMP-SMX-powered AMR-Portal organisms (organism-column names) -> all route to the E. coli-default
# acquired-gene rule shape (the overlay is gene-family presence, organism-agnostic for Enterobacterales +
# Shigella). Campylobacter etc. NOT here (no TMP-SMX powering / different flora).
TMP_SMX_ORGS = ("Escherichia coli", "Klebsiella pneumoniae", "Salmonella enterica",
                "Shigella sonnei", "Shigella flexneri")
CELLS = [(org, DRUG) for org in TMP_SMX_ORGS]


def _symbols(dets: list[dict]) -> list[str]:
    return [(d.get("amr_element_symbol") or d.get("gene_symbol") or "") for d in dets]


def score_cell_sir(isolates, geno) -> dict:
    """isolates = [(biosample, leaked, SIR)]; geno = {biosample: [determinant dicts]}.
    Bins disjoint isolates into the 4 genotype strata + computes binary confusion + the
    strata-reproduction gate (sul+dfr is the highest-R stratum AND sul-only R-rate < 0.5)."""
    strata = {k: {"R": 0, "S": 0} for k in ("sul+dfr", "sul-only", "dfr-only", "neither")}
    pairs = []
    for bs, leaked, sir in isolates:
        if leaked:
            continue                                   # provenance-disjoint only
        call = tmp_smx_call(_symbols(geno.get(bs, [])))
        has_sul, has_dfr = bool(call["matched_sul"]), bool(call["matched_dfr"])
        stratum = ("sul+dfr" if (has_sul and has_dfr) else
                   "dfr-only" if has_dfr else "sul-only" if has_sul else "neither")
        y = 1 if sir == "R" else 0
        pairs.append((call["prediction"], y))
        strata[stratum]["R" if sir == "R" else "S"] += 1
    tbl = {}
    for k, v in strata.items():
        n = v["R"] + v["S"]
        tbl[k] = {"n": n, "R": v["R"], "S": v["S"], "r_rate": round(v["R"] / n, 3) if n else None}
    rates = {k: t["r_rate"] for k, t in tbl.items() if t["r_rate"] is not None}
    sd, so = tbl["sul+dfr"]["r_rate"], tbl["sul-only"]["r_rate"]
    reproduced = bool(rates and sd is not None and sd == max(rates.values()) and (so is None or so < 0.5))
    conf = _conf(pairs)
    return {"binary": conf, "strata": tbl, "strata_reproduced": reproduced,
            "n_R": conf["tp"] + conf["fn"], "n_S": conf["tn"] + conf["fp"],
            "sens_wilson95": wilson_ci(conf["tp"], conf["tp"] + conf["fn"]),
            "spec_wilson95": wilson_ci(conf["tn"], conf["tn"] + conf["fp"]),
            "headline": "SCORED" if reproduced else "INDETERMINATE"}


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
    geno = _load_genotype(a.geno, {o for o, _ in CELLS})
    pheno = _load_phenotype_cells(a.pheno, CELLS, leak)
    cells = {}
    for (org, drug), isolates in pheno.items():
        r = score_cell_sir(isolates, geno)
        cells[f"{org}|{drug}"] = {"organism": org, "drug": drug, **r}
    n_scored = sum(1 for c in cells.values() if c["headline"] == "SCORED")
    out = {
        "_schema": "amr-portal-tmpsmx-experimental-v1", "date": _date.today().isoformat(),
        "drug": DRUG, "rule_status": RULE_STATUS, "rule_scope": RULE_SCOPE,
        "not_in_shipped_surface": True,
        "rule_text": tmp_smx_call([])["rule_text"],
        "independence_tier": ("EBI AMR Portal provenance-disjoint (BioSample/ERS/GCA disjoint vs CRyPTIC + "
                              "our tuning cohorts); genotype = the Portal's own AMRFinderPlus run; phenotype "
                              "= wet-lab measured AST (binary R/S, non-circular)"),
        "metric_note": ("binary measured AST only (Portal has no MIC -> no strict/relaxed tier); the "
                        "4-genotype-strata table + strata-reproduction gate is the experimental-honesty check "
                        "(sul+dfr highest-R AND sul-only R-rate < 0.5, reproducing the Sci234/Oxford pattern)"),
        "n_cells": len(cells), "n_scored": n_scored, "cells": cells,
        "frozen_surface_changed": False,
    }
    outp = REPO / "wiki" / f"amr_portal_tmpsmx_experimental_{_date.today().isoformat()}.json"
    outp.write_text(json.dumps(out, indent=2, default=str), encoding="utf-8")
    print(f"{'cell':<42} {'nR':>6} {'nS':>6} {'acc':>6} {'sens':>6} {'spec':>6} strata-repro headline")
    for k, c in cells.items():
        b = c["binary"]
        s = lambda x: f"{x:.3f}" if isinstance(x, float) else "  -  "
        print(f"{k:<42} {c['n_R']:>6} {c['n_S']:>6} {s(b['acc']):>6} {s(b['sens']):>6} {s(b['spec']):>6} "
              f"{str(c['strata_reproduced']):>11} {c['headline']}")
    print(f"\n{n_scored}/{len(cells)} cells SCORED (strata-reproduced) -> {outp}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
