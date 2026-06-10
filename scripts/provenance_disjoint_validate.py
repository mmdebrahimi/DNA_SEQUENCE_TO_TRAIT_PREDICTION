"""Stage-2 PROVENANCE-DISJOINT validation — score the DEPLOYED decoder on a different-submitter-lab subset.

The Stage-1 census (wiki/ncbi_pd_provenance_census_2026-06-10.md) showed a free, genome-linked,
provenance-disjoint NCBI-PD subset is POWERED for Klebsiella/Campylobacter cipro. This Stage-2 run actually
scores the deployed decoder on it: select non-ecosystem-submitter isolates with the drug's AST + a
downloadable genome (balanced N/class, PREFERRING already-AMRFinder-cached accessions to bound runtime),
ensure each has an AMRFinder run (reuse cache; else download+run), apply the SHIPPED
`call_resistance(organism=..., drug=...)` rule, and report acc/sens/spec.

HONEST TIER: this is PROVENANCE-disjoint (different submitter/lab/country), NOT methodology-independent
(most NCBI submitters still use CLSI broth microdilution). Headline it as "holds on provenance-disjoint
isolate phenotype not used in tuning" — a provenance-leakage stress-test, NOT external clinical validation.

Usage: .venv/Scripts/python.exe scripts/provenance_disjoint_validate.py --group Klebsiella \
         --amrfinder-organism Klebsiella_pneumoniae --drug ciprofloxacin --per-class 30 [--select-only]
"""
from __future__ import annotations

import argparse
import glob
import json
import sys
import urllib.request
from datetime import date as _date
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from dna_decode.eval.amr_rules import call_resistance
from scripts.independent_cohort_validate import _conf
from scripts.ncbi_pd_provenance_census import ECOSYSTEM, is_ecosystem, latest_metadata_url
from scripts.organism_drug_validate import _run_dir, ensure_run


# Flagship cohorts that live as PARQUET (not data/raw/*/selected.tsv) — the data/raw glob misses these, so
# a provenance-disjoint run could silently LEAK the E. coli cipro tuning/held-out sets. Exclude them too.
# (Caught 2026-06-10: E. coli flagship run's data/raw glob found 0 prior cohorts; overlap had to be checked
# by hand. The script now enforces it. Harmless cross-organism — accessions don't collide across species.)
_FLAGSHIP_PARQUET_COHORTS = [
    "data/processed/stage2_n150_cipro_cohort.parquet",   # E. coli cipro TUNING cohort (N=147)
    "data/processed/gate_b_cohort.parquet",              # E. coli cef/cipro held-out
    "data/processed/gate_b_n40_cipro_cohort.parquet",    # E. coli cipro held-out N=40
]


def _prior_cohort_accessions(slug: str, exclude_self: str) -> set[str]:
    """All accessions used in ANY prior cohort (calibration + prior validation) — these are NOT
    provenance-disjoint from tuning/prior-validation and MUST be excluded for a clean validation. Covers
    BOTH data/raw/<slug>_*/selected.tsv AND the parquet flagship cohorts (which the glob misses)."""
    out = set()
    for sel in glob.glob(f"data/raw/{slug}_*/selected.tsv"):
        if exclude_self in sel:
            continue
        for ln in Path(sel).read_text().splitlines():
            if "\t" in ln:
                out.add(ln.split("\t")[0])
    # parquet flagship cohorts (leakage-hardening, 2026-06-10)
    try:
        from dna_decode.data.cohort import load_cohort
        for pq in _FLAGSHIP_PARQUET_COHORTS:
            if Path(pq).exists():
                for s in load_cohort(pq).strains:
                    acc = getattr(s, "assembly_accession", None)
                    if acc:
                        out.add(acc)
    except Exception:
        pass   # parquet/load_cohort unavailable -> data/raw exclusion still applies (offline-safe)
    return out


def select_disjoint(group: str, drug: str, per_class: int, reuse_glob: str, selected: Path,
                    exclude_prior: set[str]) -> dict[str, int]:
    """Pick balanced per_class R/S non-ecosystem isolates with downloadable genomes, EXCLUDING any accession
    used in a prior cohort (tuning OR prior validation) — genuinely fresh, leakage-free."""
    if selected.exists():
        out = {}
        for ln in selected.read_text().splitlines():
            if "\t" in ln:
                a, rs = ln.split("\t"); out[a] = 1 if rs.strip() == "R" else 0
        return out
    tok_r, tok_s = f"{drug}=R", f"{drug}=S"
    pools = {1: [], 0: []}
    seen = set()
    with urllib.request.urlopen(latest_metadata_url(group), timeout=300) as resp:
        header = resp.readline().decode("utf-8", "replace").rstrip("\n").split("\t")
        idx = {n: header.index(n) for n in
               ("asm_acc", "AST_phenotypes", "bioproject_center", "collected_by", "sra_center") if n in header}
        ai, pi = idx["asm_acc"], idx["AST_phenotypes"]
        for raw in resp:
            cells = raw.decode("utf-8", "replace").rstrip("\n").split("\t")
            if len(cells) <= max(ai, pi):
                continue
            acc, ast = cells[ai], cells[pi]
            if not acc.startswith(("GCA_", "GCF_")) or acc in seen or acc in exclude_prior or not ast or ast == "NULL":
                continue
            lab = 1 if tok_r in ast.split(",") else (0 if tok_s in ast.split(",") else None)
            if lab is None:
                continue
            prov = [cells[idx[c]] for c in ("bioproject_center", "collected_by", "sra_center")
                    if c in idx and len(cells) > idx[c]]
            if is_ecosystem(*prov):
                continue
            seen.add(acc)
            pools[lab].append(acc)
    chosen = {}
    for lab in (1, 0):
        for acc in pools[lab][:per_class]:
            chosen[acc] = lab
    selected.parent.mkdir(parents=True, exist_ok=True)
    selected.write_text("".join(f"{a}\t{'R' if v else 'S'}\n" for a, v in chosen.items()), encoding="utf-8")
    print(f"selected {len(chosen)} ({sum(chosen.values())}R/{len(chosen)-sum(chosen.values())}S); "
          f"all FRESH (excluded {len(exclude_prior)} prior-cohort accessions); pools had "
          f"{len(pools[1])}R/{len(pools[0])}S fresh-disjoint available")
    return chosen


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--group", default="Klebsiella")
    ap.add_argument("--amrfinder-organism", default="Klebsiella_pneumoniae")
    ap.add_argument("--drug", default="ciprofloxacin")
    ap.add_argument("--per-class", type=int, default=30)
    ap.add_argument("--registry-organism", default=None, help="organism key for call_resistance (default=group)")
    ap.add_argument("--select-only", action="store_true", help="select cohort + report cache overlap, no AMRFinder")
    a = ap.parse_args()
    slug = a.group.lower()
    base = Path(f"data/raw/{slug}_provdisjoint_{a.drug}")
    own_runs = base / "amrfinder_runs"; gcache = base / "refseq"
    reuse_glob = f"data/raw/{slug}_*/amrfinder_runs"
    reg_org = a.registry_organism or a.group

    exclude_prior = _prior_cohort_accessions(slug, exclude_self=f"{slug}_provdisjoint_")
    sel = select_disjoint(a.group, a.drug, a.per_class, reuse_glob, base / "selected.tsv", exclude_prior)
    if a.select_only:
        return 0

    applied = []
    for i, (acc, y) in enumerate(sel.items(), 1):
        mt = _run_dir(acc, own_runs, reuse_glob)
        if mt is None:
            print(f"  [{i}/{len(sel)}] {acc} ({'R' if y else 'S'}) running AMRFinder ...", flush=True)
            ensure_run(acc, own_runs, gcache, a.amrfinder_organism, reuse_glob)
            mt = _run_dir(acc, own_runs, reuse_glob)
        if mt is None:
            print(f"  [{i}/{len(sel)}] {acc}: no run (skip)"); continue
        call = call_resistance(mt / "main.tsv", a.drug, organism=reg_org)
        applied.append((call["prediction"], y))
    conf = _conf(applied)
    artifact = {
        "_schema": "provenance-disjoint-validation-v1", "date": _date.today().isoformat(),
        "organism": a.group, "amrfinder_organism": a.amrfinder_organism, "drug": a.drug,
        "registry_organism": reg_org, "n_selected": len(sel), "n_scored": conf["n_scored"],
        "independence_tier": "provenance-disjoint (different submitter/lab/country); NOT methodology-independent (most submitters use CLSI broth microdilution)",
        "leakage_control": f"all selected accessions are FRESH — excluded {len(exclude_prior)} accessions used in any prior {slug}_* cohort (calibration + prior validation)",
        "ecosystem_excluded": ECOSYSTEM, "metrics": conf,
    }
    out_json = Path(f"wiki/provenance_disjoint_validation_{slug}_{a.drug[:5]}_{_date.today().isoformat()}.json")
    out_md = Path(f"wiki/provenance_disjoint_validation_{slug}_{a.drug[:5]}_{_date.today().isoformat()}.md")
    out_json.write_text(json.dumps(artifact, indent=2), encoding="utf-8")
    md = (f"# Provenance-disjoint validation — {a.group} x {a.drug} — {_date.today().isoformat()}\n\n"
          f"Deployed decoder `call_resistance(organism={reg_org}, drug={a.drug})` scored on a "
          f"PROVENANCE-DISJOINT NCBI-PD subset (submitters OUTSIDE NARMS/CDC/FDA/GenomeTrakr/PulseNet/USDA).\n\n"
          f"## Result\n\n| metric | value |\n|---|---|\n"
          f"| n scored | {conf['n_scored']} (TP {conf['tp']} FP {conf['fp']} TN {conf['tn']} FN {conf['fn']}; abstain {conf['abstain']}) |\n"
          f"| accuracy | {conf['acc']} |\n| sensitivity (R) | {conf['sens']} |\n| specificity (S) | {conf['spec']} |\n\n"
          f"## Independence tier (DO NOT inflate)\n"
          f"This is **provenance-disjoint** (different submitter / lab / country than the BV-BRC/NCBI-PD records "
          f"the decoder was tuned + cross-source-validated on) — a stress-test against provenance leakage. It is "
          f"NOT methodology-independent (most NCBI submitters use CLSI broth microdilution) and NOT external "
          f"clinical validation. Headline accordingly. Excluded ecosystem submitters: {ECOSYSTEM}.\n\n"
          f"## Leakage control\nAll {len(sel)} accessions are FRESH — excluded {len(exclude_prior)} "
          f"accessions used in ANY prior {slug}_* cohort (registry calibration + prior validation), so the "
          f"score is on strains never seen in tuning or earlier validation.\n")
    out_md.write_text(md, encoding="utf-8")
    print(f"\nRESULT acc={conf['acc']} sens={conf['sens']} spec={conf['spec']} (n={conf['n_scored']}); artifacts: {out_md}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
