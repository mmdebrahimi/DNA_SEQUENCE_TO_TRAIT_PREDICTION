"""GENERIC no-compute external validation of NON-FROZEN organism cells on NCBI Pathogen Detection.

The gonococcus run (2026-07-21) discovered that NCBI-PD publishes its OWN AMRFinderPlus calls (the
`AMR_genotypes` field, point mutations included) for every organism it covers, and the isolates have
downloadable assemblies -> external validation of any determinant-based cell is a pure metadata JOIN + score,
OFFLINE, no Docker/Kaggle/assembly. This is the REUSABLE substrate: point it at any organism's cohort +
determinants (built from `<PDG>.amr.metadata.tsv`) + its NON-FROZEN cell.

Each organism config: the cell module's per-drug call functions. Score vs measured AST; DEGENERATE guard
(a cell predicting all-one-class is never endorsed even at spec/sens 1.0).

  uv run python scripts/score_ncbipd_extval.py --organism campylobacter
  uv run python scripts/score_ncbipd_extval.py --organism gono

HONEST CAVEATS (same as any cell): provenance-disjoint (frozen accessions excluded at cohort-build) but NOT
methodology-independent (same AMRFinderPlus + same cell); RAW sens/spec is CLONALITY-INFLATED (lineage
collapse via Mash on the assemblies is the follow-up).
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from datetime import date as _date
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import re  # noqa: E402

from dna_decode.eval.clonality import cluster_weighted_confusion  # noqa: E402
from dna_decode.organism_rules import campylobacter_amr as CJ  # noqa: E402
from dna_decode.organism_rules import neisseria_amr as NG  # noqa: E402
from dna_decode.organism_rules import pneumo_amr as PN  # noqa: E402
from dna_decode.organism_rules import staphylococcus_amr as SA  # noqa: E402
from scripts.amr_portal_score_independent import wilson_ci  # noqa: E402
from scripts.independent_cohort_validate import _conf  # noqa: E402


def _pneumo(drug):
    """Adapter: pneumo_amr.call_drug takes normalized gene tokens (its rule keys are 'ermb'/'mefa'/'tetm'
    with no parens), but NCBI-PD emits `erm(B)`/`tet(M)` -> strip non-alphanumerics so the substring match
    fires. Returns the scorer's {'prediction','rule'} dict; out-of-scope drug -> INDETERMINATE."""
    def call(symbols):
        norm = [re.sub(r"[^a-z0-9]", "", (s or "").lower()) for s in symbols]
        c = PN.call_drug(drug, norm)
        if c is None:
            return {"prediction": "INDETERMINATE", "rule": f"pneumo {drug} out of v0 scope"}
        return {"prediction": c.prediction, "rule": f"pneumo {drug}: {'|'.join(c.rule_tokens)} gene-presence -> R"}
    return call

SPEC_FLOOR = 0.85
MIN_PER_CLASS = 5

# organism -> (display name, cohort dir, {drug: call function})
ORGANISMS = {
    "gono": ("Neisseria gonorrhoeae", "data/raw/gono_ncbipd_extval", {
        "ciprofloxacin": NG.call_ng_ciprofloxacin, "cefixime": NG.call_ng_cefixime,
        "ceftriaxone": NG.call_ng_ceftriaxone, "azithromycin": NG.call_ng_azithromycin,
        "penicillin": NG.call_ng_penicillin, "tetracycline": NG.call_ng_tetracycline,
    }),
    "campylobacter": ("Campylobacter", "data/raw/campy_ncbipd_extval", {
        "tetracycline": CJ.call_cj_tetracycline, "gentamicin": CJ.call_cj_gentamicin,
    }),
    "staph": ("Staphylococcus aureus", "data/raw/staph_ncbipd_extval", {
        "ciprofloxacin": SA.call_sa_ciprofloxacin, "rifampin": SA.call_sa_rifampicin,
    }),
    "pneumo": ("Streptococcus pneumoniae", "data/raw/pneumo_ncbipd_extval", {
        "erythromycin": _pneumo("erythromycin"), "tetracycline": _pneumo("tetracycline"),
    }),
}


def score_cell(name: str, cohort_dir: str, drugs: dict) -> dict:
    labels = {r["biosample"]: r for r in csv.DictReader(open(f"{cohort_dir}/cohort.tsv"), delimiter="\t")}
    dets = {r["biosample"]: [s for s in (r["determinants"] or "").split(";") if s]
            for r in csv.DictReader(open(f"{cohort_dir}/determinants.tsv"), delimiter="\t")}
    # NCBI-PD SNP clusters (biosample -> PDS); NULL -> a unique singleton lineage. No Mash/Docker needed.
    clusters, _sing = {}, [10 ** 7]
    cpath = Path(f"{cohort_dir}/clusters.tsv")
    if cpath.exists():
        for r in csv.DictReader(open(cpath), delimiter="\t"):
            pds = r.get("PDS_acc") or "NULL"
            if pds and pds != "NULL":
                clusters[r["biosample"]] = hash(pds) % (10 ** 6)
            else:
                _sing[0] += 1
                clusters[r["biosample"]] = _sing[0]
    print(f"{name}: {len(labels)} NCBI-PD isolates (provenance-disjoint; frozen accessions excluded)"
          + (f"; {len(clusters)} SNP-clustered -> lineage-collapse ON" if clusters else ""))
    results = {}
    for drug, fn in drugs.items():
        scored, preds, labs = [], {}, {}
        for bs, row in labels.items():
            rs = row.get(drug, "")
            if rs not in ("R", "S"):
                continue
            pred = fn(dets.get(bs, []))["prediction"]
            if str(pred).upper() in ("R", "S"):
                scored.append((pred, 1 if rs == "R" else 0))
                preds[bs], labs[bs] = pred, rs
        conf = _conf(scored)
        n_R, n_S = conf["tp"] + conf["fn"], conf["tn"] + conf["fp"]
        powered = n_R >= MIN_PER_CLASS and n_S >= MIN_PER_CLASS
        degenerate = ((n_R >= MIN_PER_CLASS and conf["sens"] == 0.0)
                      or (n_S >= MIN_PER_CLASS and conf["spec"] == 0.0))
        endorsed = bool(powered and not degenerate and conf["spec"] is not None
                        and conf["spec"] >= SPEC_FLOOR and (conf["sens"] is None or conf["sens"] >= 0.5))
        headline = ("DEGENERATE_NOT_ENDORSED" if degenerate else "SCORED_ENDORSED" if endorsed
                    else "UNDERPOWERED" if not powered else "SCORED_NOT_ENDORSED")
        # lineage-collapsed (clonality-corrected) — one vote per SNP cluster; mixed clusters = DISCORDANT
        lineage = None
        if clusters and preds:
            clus = {bs: clusters[bs] for bs in preds if bs in clusters}
            lineage = cluster_weighted_confusion({bs: preds[bs] for bs in clus},
                                                 {bs: labs[bs] for bs in clus}, clus)
        results[drug] = {"drug": drug, "rule": fn(["_probe_"])["rule"], "binary": conf,
                         "n_R": n_R, "n_S": n_S, "sens_wilson95": wilson_ci(conf["tp"], n_R),
                         "spec_wilson95": wilson_ci(conf["tn"], n_S), "powered": powered,
                         "headline": headline, "lineage_collapsed": lineage}
        lstr = (f" | LINEAGE sens={lineage.get('sens')} spec={lineage.get('spec')} "
                f"discordant={lineage.get('n_discordant')}") if lineage else ""
        print(f"  {drug:14s}: n={conf['n_scored']:3d} ({n_R}R/{n_S}S) acc={conf['acc']} "
              f"sens={conf['sens']} spec={conf['spec']} -> {headline}{lstr}")
    return results


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--organism", choices=list(ORGANISMS), required=True)
    a = ap.parse_args()
    name, cohort_dir, drugs = ORGANISMS[a.organism]
    results = score_cell(name, cohort_dir, drugs)
    artifact = {
        "_schema": "ncbipd-extval-v1", "date": _date.today().isoformat(), "organism": name,
        "cell": f"{a.organism}_amr (NON-FROZEN)", "cohort_dir": cohort_dir,
        "independence": "provenance-disjoint (frozen accessions excluded) but NOT methodology-independent",
        "clonality_caveat": "RAW is clonality-inflated; lineage-collapse (Mash) is the follow-up",
        "results": results, "frozen_surface_changed": False,
    }
    out = Path(f"wiki/{a.organism}_ncbipd_extval_{_date.today().isoformat()}.json")
    out.write_text(json.dumps(artifact, indent=2), encoding="utf-8")
    print(f"artifact: {out}")
    return 0 if any(r["powered"] for r in results.values()) else 1


if __name__ == "__main__":
    raise SystemExit(main())
