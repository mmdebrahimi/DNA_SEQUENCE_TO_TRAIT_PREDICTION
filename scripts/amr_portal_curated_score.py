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

import tempfile  # noqa: E402

from scripts.amr_portal_score_independent import (  # noqa: E402
    DEFAULT_CRYPTIC, DEFAULT_PARQUET, GENO_PARQUET,
    _load_genotype, _load_leak_set, _load_phenotype_cells, genotype_to_main_tsv, wilson_ci,
)
from dna_decode.eval.amr_rules import call_resistance  # noqa: E402  (FROZEN rule, called unchanged)
from scripts.independent_cohort_validate import _conf  # noqa: E402

_TMP = Path(tempfile.gettempdir()) / "_amrportal_curated_main.tsv"


def default_rule_fn(drug: str, organism_route: str | None = None):
    """rule_fn that applies the FROZEN deployed `call_resistance` default rule (organism_route=None -> the
    E. coli DRUG_RULE default) to an isolate's determinants. For Enterobacterales acquired-gene cells whose
    flora shares E. coli's acquired determinants; the spec>=0.85 guard catches the intrinsic over-call
    (e.g. AmpC cef). ABSTAIN (expression-floor drugs like meropenem) -> excluded from the binary."""
    if drug == "trimethoprim-sulfamethoxazole":
        # TMP-SMX is not in the frozen DRUG_RULE (call_resistance) -> use the EXPERIMENTAL sul-AND-dfr overlay.
        from dna_decode.data.experimental_drug_rules import tmp_smx_call

        def _fn(dets):
            syms = [(d.get("amr_element_symbol") or d.get("gene_symbol") or "") for d in dets]
            c = tmp_smx_call(syms)
            return {"prediction": c["prediction"], "determinant_present": c["prediction"] == "R"}
        return _fn

    def _fn(dets):
        _TMP.write_text(genotype_to_main_tsv(dets), encoding="utf-8")
        pred = call_resistance(_TMP, drug, organism=organism_route)["prediction"]
        return {"prediction": pred if pred in ("R", "S") else "ABSTAIN",
                "determinant_present": pred == "R"}
    return _fn

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


def _efm(fn_name, key):
    def _fn(dets):
        from dna_decode.organism_rules import enterococcus_amr as e
        c = getattr(e, fn_name)(_syms(dets))
        return {"prediction": c["prediction"], "determinant_present": bool(c[key])}
    return _fn


def _cj(fn_name, key):
    def _fn(dets):
        from dna_decode.organism_rules import campylobacter_amr as c
        r = getattr(c, fn_name)(_syms(dets))
        return {"prediction": r["prediction"], "determinant_present": bool(r[key])}
    return _fn


CURATED_CELLS = {
    "neisseria_tet": ("Neisseria gonorrhoeae", "tetracycline", _ng_tet, "neisseria_tet",
                      "tet(M) -> R (rpsJ V57M accessory-only)"),
    "efm_cipro": ("Enterococcus faecium", "ciprofloxacin", _efm("call_efm_ciprofloxacin", "matched_qrdr"),
                  "enterococcus_faecium_cipro", "gyrA/parC QRDR -> R"),
    "efm_tet": ("Enterococcus faecium", "tetracycline", _efm("call_efm_tetracycline", "matched_tet"),
                "enterococcus_faecium_tet", "acquired tet gene -> R"),
    "efm_gent": ("Enterococcus faecium", "gentamicin", _efm("call_efm_gentamicin", "matched_aph2"),
                 "enterococcus_faecium_gent", "aph(2'') high-level -> R (intrinsic aac(6')-Ii excluded)"),
    "cj_tet": ("Campylobacter jejuni", "tetracycline", _cj("call_cj_tetracycline", "matched_tetO"),
               "campylobacter_jejuni_tet", "tet(O)-family ribosomal protection -> R"),
    "cj_gent": ("Campylobacter jejuni", "gentamicin", _cj("call_cj_gentamicin", "matched_gent"),
                "campylobacter_jejuni_gent", "true gent enzyme aph(2'')/aac(3) -> R (aad9/spw non-gent excluded)"),
    "cc_tet": ("Campylobacter coli", "tetracycline", _cj("call_cj_tetracycline", "matched_tetO"),
               "campylobacter_coli_tet", "tet(O)-family ribosomal protection -> R"),
    "cc_gent": ("Campylobacter coli", "gentamicin", _cj("call_cj_gentamicin", "matched_gent"),
                "campylobacter_coli_gent", "true gent enzyme aph(2'')/aac(3) -> R (aad9/spw non-gent excluded)"),
}


_DRUG_ABBR = {"ciprofloxacin": "cipro", "gentamicin": "gent", "tetracycline": "tet",
              "trimethoprim-sulfamethoxazole": "tmpsmx", "ceftriaxone": "cef", "meropenem": "mero"}


def _slug(org, drug):
    return org.lower().replace(" ", "_") + "_" + _DRUG_ABBR.get(drug, drug.replace("-", "")[:6])


# Phase B — Enterobacterales acquired-gene cells scored via the FROZEN default rule + spec-guard.
PHASE_B = ([("Enterobacter cloacae", d) for d in
            ("ciprofloxacin", "gentamicin", "tetracycline", "trimethoprim-sulfamethoxazole",
             "ceftriaxone", "meropenem")]
           + [("Proteus mirabilis", d) for d in ("gentamicin", "trimethoprim-sulfamethoxazole")]
           + [("Serratia marcescens", d) for d in ("ceftriaxone", "ciprofloxacin", "meropenem")])


def run_batch(cells) -> None:
    for org, drug in cells:
        run_cell(org, drug, default_rule_fn(drug), _slug(org, drug),
                 f"FROZEN deployed DRUG_RULE default (organism-agnostic) on {org} — spec-guard endorses")


def main(argv=None) -> int:
    import argparse
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("cell", nargs="?", choices=sorted(CURATED_CELLS))
    ap.add_argument("--batch", choices=["phaseB"])
    a = ap.parse_args(argv)
    if a.batch == "phaseB":
        run_batch(PHASE_B)
        return 0
    org, drug, fn, slug, txt = CURATED_CELLS[a.cell]
    run_cell(org, drug, fn, slug, txt)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
