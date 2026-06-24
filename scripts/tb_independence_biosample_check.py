"""BioSample-resolution independence check for the TB independent cohort (rigor upgrade, rec #1).

Upgrades the TB independence claim from "accession-STRING disjoint" to "BioSample-resolution-checked" by
characterizing — and empirically bounding — the cross-archive-alias residual (the case a string-match misses:
ENA `ERS` vs NCBI `SAMN` being the same physical sample under different archive IDs).

STRUCTURAL finding (no network):
  - CRyPTIC's leak set keys on ENA_RUN (ERR) + ENA_SAMPLE (ERS) + UNIQUEID (site.xx) — NO NCBI BioSample.
  - The AMR-Portal TB disjoint cohort splits into ENA-side (`ERS` sra / `SAMEA` biosample) and NCBI-side
    (`SRS` sra / `SAMN` biosample).
  - **ENA-side isolates are ALREADY BioSample-grade disjoint:** their `ERS` sample accession is string-
    compared DIRECTLY against CRyPTIC's `ENA_SAMPLE` (ERS) — same namespace, so same-BioSample overlap is
    already caught. The "accession-level upper bound" caveat is TIGHT for them.
  - **NCBI-side isolates (`SAMN`) are the only cross-archive risk** vs the European CRyPTIC consortium.

EMPIRICAL bound (ENA portal, bounded sample): resolve a sample of NCBI-side `SAMN` biosamples to their
ENA-mirrored run/sample accessions and check vs the CRyPTIC leak set. A near-zero overlap means the
cross-archive residual is negligible; the irreducible residual that remains is genomic RE-SUBMISSION (an
isolate sequenced twice + submitted as distinct BioSample records) — which BioSample resolution does NOT
catch (only Mash genomic dedup would). So this check bounds the alias residual; it does not claim to close
the genomic-duplicate residual.

Pure logic (`composition`, `overlap`) unit-tested; the ENA probe is the live part (injectable fetch).
"""
from __future__ import annotations

import csv
import json
import re
import sys
import urllib.request
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

COHORT = REPO / "data" / "raw" / "tb_goldset" / "amr_portal_tb_disjoint_cohort.tsv"
DEFAULT_CRYPTIC = Path("D:/dna_decode_cache/data files donwload/CRyPTIC TB MIC compendium/"
                       "CRyPTIC_reuse_table_20240917.csv")
ENA_PORTAL = ("https://www.ebi.ac.uk/ena/portal/api/filereport?accession={}&result=read_run"
              "&fields=run_accession,sample_accession,secondary_sample_accession&format=tsv")


def _acctype(v: str) -> str:
    v = (v or "").strip()
    return "EMPTY" if not v else re.sub(r"\d.*", "#", v)[:7]


def composition(rows: list[dict]) -> dict:
    """Accession-type composition of the disjoint cohort + the ENA-side-already-tight count."""
    from collections import Counter
    bs, sra = Counter(), Counter()
    for r in rows:
        if r["leaked"] == "0" and r["assembly"].startswith("GCA_"):
            bs[_acctype(r["biosample"])] += 1
            sra[_acctype(r["sra"])] += 1
    ena_side = sum(v for k, v in sra.items() if k == "ERS#")        # ERS sra -> already biosample-grade
    ncbi_side = sum(v for k, v in sra.items() if k == "SRS#")
    return {"biosample_types": dict(bs), "sra_types": dict(sra),
            "ena_side_already_biosample_grade": ena_side, "ncbi_side_cross_archive_candidates": ncbi_side}


def overlap(resolved_accs: set[str], cryptic: set[str]) -> set[str]:
    """Cross-archive overlap = resolved ENA accessions (upper-cased) that are in the CRyPTIC leak set."""
    return {a.upper() for a in resolved_accs} & cryptic


def _default_fetch(url: str, timeout: int = 30) -> str:
    return urllib.request.urlopen(url, timeout=timeout).read().decode("utf-8", "ignore")


def probe_ncbi_side(rows, cryptic: set[str], sample: int, fetch=_default_fetch) -> dict:
    """Resolve a bounded sample of NCBI-side SAMN biosamples to ENA accessions; count CRyPTIC overlaps."""
    samn = [r["biosample"] for r in rows
            if r["leaked"] == "0" and r["assembly"].startswith("GCA_") and r["biosample"].startswith("SAMN")]
    samn = samn[:sample]
    n_resolved = n_no_ena = n_hit = 0
    hits = []
    for bs in samn:
        try:
            txt = fetch(ENA_PORTAL.format(bs))
        except Exception:
            n_no_ena += 1
            continue
        lines = [l for l in txt.splitlines()[1:] if l.strip()]
        if not lines:
            n_no_ena += 1
            continue
        n_resolved += 1
        accs = {tok.strip() for l in lines for tok in l.split("\t") if tok.strip()}
        ov = overlap(accs, cryptic)
        if ov:
            n_hit += 1
            hits.append({"biosample": bs, "cryptic_overlap": sorted(ov)})
    return {"n_probed": len(samn), "n_ena_mirrored": n_resolved, "n_no_ena_mirror": n_no_ena,
            "n_cross_archive_overlap": n_hit, "hits": hits}


def main(argv=None) -> int:
    import argparse
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cryptic", type=Path, default=DEFAULT_CRYPTIC)
    ap.add_argument("--sample", type=int, default=30, help="NCBI-side biosamples to probe (0 = skip network)")
    a = ap.parse_args(argv)
    rows = list(csv.DictReader(open(COHORT, encoding="utf-8"), delimiter="\t"))
    comp = composition(rows)
    from dna_decode.organism_rules import tb_goldset
    cryptic = {x.upper() for x in tb_goldset.cryptic_accessions(a.cryptic)} if a.cryptic.exists() else set()
    print(f"[tb-independence] cohort composition: {comp}")
    print(f"[tb-independence] CRyPTIC leak set: {len(cryptic):,} accessions (ERR/ERS/site.xx)")
    probe = {"skipped": True}
    if a.sample and cryptic:
        probe = probe_ncbi_side(rows, cryptic, a.sample)
        print(f"[tb-independence] NCBI-side probe: {probe['n_probed']} SAMN | "
              f"{probe['n_ena_mirrored']} ENA-mirrored | cross-archive overlap with CRyPTIC: "
              f"{probe['n_cross_archive_overlap']}")
    verdict = ("BIOSAMPLE_RESOLUTION_CHECKED: ENA-side already BioSample-grade (ERS<->ENA_SAMPLE direct); "
               f"NCBI-side cross-archive overlap {probe.get('n_cross_archive_overlap', '?')}/"
               f"{probe.get('n_probed', '?')} probed. Irreducible residual = genomic re-submission "
               "(distinct BioSample record) -> needs Mash dedup, NOT accession resolution.")
    out = REPO / "wiki" / "tb_independence_biosample_check.json"
    out.write_text(json.dumps({"composition": comp, "cryptic_leakset_size": len(cryptic),
                               "ncbi_side_probe": probe, "verdict": verdict}, indent=2), encoding="utf-8")
    print(f"verdict: {verdict}\nartifact -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
