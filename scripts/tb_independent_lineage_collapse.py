"""Lineage-collapsed (publication-grade) independent TB number from the AMR-Portal cohort VCFs.

The raw sens/spec from `run_tb_independent_amr_portal.py` is CLONALITY-INFLATED — TB R classes are clonally
dominated, so one over-sampled clone carries the metric. The honest headline is the LINEAGE-COLLAPSED number
(each same-lineage same-prediction set counts once), exactly as the in-distribution CRyPTIC baseline reports.

This is a thin post-processor that reuses the FROZEN lineage-collapse path UNCHANGED:
  per-isolate masked VCF (already produced by the runner; FILTER .->PASS) + measured label + the WHO
  determinants + the pinned Napier barcode  ->  `score_tb_cryptic.run_v1b`  ->  raw + lineage-collapsed
  sens/spec + Wilson CI + effective-lineage-N + n_discordant, via `clonality.cluster_weighted_confusion`.

NO Docker, NO re-fetch: the barcode lineage is read from each isolate's VCF calls (the same masked VCFs on
disk). Runs on whatever VCFs exist (the smoke 60 now; the full 2,845 when the background align finishes).
Honest rails identical to the raw memo: accession-level provenance-disjoint; measured phenotype; WHO rule
UNCHANGED; no regeno -> callability unassessed.
"""
from __future__ import annotations

import csv
import json
import os
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

COHORT = REPO / "data" / "raw" / "tb_goldset" / "amr_portal_tb_disjoint_cohort.tsv"
WORK = Path(os.environ.get("TB_INDEP_WORK", str(REPO / "data" / "raw" / "tb_indep")))


def _pass_mask(vcf_text: str) -> str:
    return "\n".join(
        ln if ln.startswith("#") else
        "\t".join(p if i != 6 or p != "." else "PASS" for i, p in enumerate(ln.split("\t")))
        for ln in vcf_text.splitlines())


def load_cohort_labels() -> dict[str, dict]:
    out = {}
    with open(COHORT, encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["leaked"] == "0" and r["assembly"].startswith("GCA_"):
                out[r["strain_id"]] = {"rif": r["rif_label"], "inh": r["inh_label"]}
    return out


def load_masked(work: Path, strain_ids) -> dict[str, str]:
    """{strain_id: PASS-masked VCF text} for the isolates whose VCF exists on disk."""
    vcfdir = work / "vcf"
    masked = {}
    for sid in strain_ids:
        vp = vcfdir / f"{sid}.vcf"
        if vp.exists() and vp.stat().st_size > 0:
            masked[sid] = _pass_mask(vp.read_text(encoding="utf-8"))
    return masked


def collapse(work: Path = WORK) -> dict:
    from dna_decode.data import tb_lineage_barcode, tb_who_catalogue
    from scripts.score_tb_cryptic import run_v1b
    tb_who_catalogue.verify_pins()
    barcode = tb_lineage_barcode.load_barcode()
    dets = {"rifampicin": tb_who_catalogue.load_determinants("rifampicin"),
            "isoniazid": tb_who_catalogue.load_determinants("isoniazid")}
    labels = load_cohort_labels()
    masked = load_masked(work, labels)
    result = {"n_isolates_with_vcf": len(masked), "drugs": {}}
    for drug, code in (("rifampicin", "rif"), ("isoniazid", "inh")):
        strain_label = {sid: labels[sid][code] for sid in masked if labels[sid][code] in ("R", "S")}
        strain_masked = {sid: masked[sid] for sid in strain_label}
        res = run_v1b(strain_masked, strain_label, dets[drug], barcode,
                      drug=drug, cohort_complete=False)
        result["drugs"][drug] = res
    return result


def main(argv=None) -> int:
    import argparse
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--work", type=Path, default=WORK)
    a = ap.parse_args(argv)
    res = collapse(a.work)
    out = REPO / "wiki" / "tb_independent_amr_portal_lineage_collapsed.json"
    out.write_text(json.dumps(res, indent=2, default=str), encoding="utf-8")
    print(f"[tb-indep-lineage] {res['n_isolates_with_vcf']} isolates with VCF")
    for drug, r in res["drugs"].items():
        st = r.get("status")
        raw = r.get("raw", {})
        lc = r.get("lineage_collapsed", {})
        print(f"  {drug}: status={st}")
        print(f"    raw      : sens={raw.get('sens')} spec={raw.get('spec')} (nR={raw.get('n_R')} nS={raw.get('n_S')})")
        if lc:
            print(f"    lineage  : sens={lc.get('sens')} spec={lc.get('spec')} "
                  f"(R-lineages={lc.get('n_clusters_R')} S-lineages={lc.get('n_clusters_S')} "
                  f"discordant={lc.get('n_discordant')})")
    print(f"artifact -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
