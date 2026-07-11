"""Fresh-cohort RE-VALIDATION of the 10 frozen SCORED cells — the shipped decoder on UNSEEN genomes.

WHY (offered a long compute window 2026-07-11): every move that produces a genuinely-NEW validated number
is closed or acquisition-gated, so the highest honest value from a long run is VERIFICATION. This re-scores
each deployed SCORED cell on a SECOND, freshly-selected provenance-disjoint cohort — genomes EXCLUDED from
every prior cohort (tuning + the frozen validation), i.e. NOT the same isolates the committed number came
from. If the cell's sens/spec holds on brand-new genomes, that is the strongest reproducibility-of-the-CLAIM
evidence short of external clinical validation; if it drifts, that is a real finding to surface before any
acquisition/publication.

SAFETY (load-bearing — this runs UNATTENDED):
  * Everything is routed into ISOLATED dirs: cohorts → `data/raw/revalidation_<date>/<cell>/`, summary →
    `wiki/revalidation_<date>/`. It NEVER writes `data/raw/*_provdisjoint_*/selected.tsv`, the committed
    `wiki/provenance_disjoint_validation_*.json`, or the report card — so a re-run cannot mutate the frozen
    reproducibility surface, and the report-card top-level glob (`wiki/provenance_disjoint_validation_*.json`)
    never sees these (they live in a subdir).
  * The isolated cohorts sit at `data/raw/revalidation_<date>/<cell>/selected.tsv` (TWO levels deep), so the
    leakage manifest (`data/raw/*/selected.tsv`, one level) never picks them up as a prior cohort either.
  * Checkpointed per cell (restartable): a wedged Docker mid-run leaves a partial isolated dir + a summary
    that says which cells completed; nothing frozen is touched.
  * Composes the FROZEN `provenance_disjoint_validate` functions unchanged (select_disjoint / ensure_run /
    _run_dir / call_resistance / _conf); zero edits to that script.

Run:  uv run python scripts/revalidate_scored_cells.py            # all 10 cells (hours; Docker + network)
      uv run python scripts/revalidate_scored_cells.py --smoke    # per_class=2 sanity (a few genomes)
"""
from __future__ import annotations

import argparse
import glob
import json
import sys
from datetime import date as _date
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from dna_decode.eval.amr_rules import call_resistance  # noqa: E402
from dna_decode.eval.cohort_manifest import build_manifest, prior_accessions  # noqa: E402
from scripts.independent_cohort_validate import _conf  # noqa: E402
from scripts.organism_drug_validate import _run_dir, ensure_run  # noqa: E402
from scripts.provenance_disjoint_validate import select_disjoint  # noqa: E402

# The 10 deployed SCORED cells: (group, amrfinder_organism, drug, registry_organism).
CELLS = [
    ("Campylobacter", "Campylobacter", "ciprofloxacin", "Campylobacter"),
    ("Escherichia_coli_Shigella", "Escherichia", "ceftriaxone", "Escherichia_coli_Shigella"),
    ("Escherichia_coli_Shigella", "Escherichia", "ciprofloxacin", "Escherichia_coli_Shigella"),
    ("Escherichia_coli_Shigella", "Escherichia", "gentamicin", "Escherichia_coli_Shigella"),
    ("Escherichia_coli_Shigella", "Escherichia", "tetracycline", "Escherichia_coli_Shigella"),
    ("Klebsiella", "Klebsiella_pneumoniae", "ceftriaxone", "Klebsiella"),
    ("Klebsiella", "Klebsiella_pneumoniae", "ciprofloxacin", "Klebsiella"),
    ("Klebsiella", "Klebsiella_pneumoniae", "gentamicin", "Klebsiella"),
    ("Klebsiella", "Klebsiella_pneumoniae", "meropenem", "Klebsiella"),
    ("Klebsiella", "Klebsiella_pneumoniae", "tetracycline", "Klebsiella"),
]

DATE = _date.today().isoformat()
ISO_ROOT = REPO / "data" / "raw" / f"revalidation_{DATE}"
OUT_DIR = REPO / "wiki" / f"revalidation_{DATE}"


def frozen_number(group: str, drug: str) -> dict | None:
    """The committed frozen provenance-disjoint sens/spec for this cell (latest dated JSON at wiki/ TOP level)."""
    slug = group.lower()
    cands = sorted(glob.glob(str(REPO / "wiki" / f"provenance_disjoint_validation_{slug}_{drug[:5]}_*.json")))
    if not cands:
        return None
    d = json.loads(Path(cands[-1]).read_text(encoding="utf-8"))
    m = d.get("metrics", {})
    return {"sens": m.get("sens"), "spec": m.get("spec"), "n_scored": m.get("n_scored"),
            "date": d.get("date"), "file": Path(cands[-1]).name}


def revalidate_cell(group: str, amrfinder_org: str, drug: str, reg_org: str, per_class: int) -> dict:
    slug = group.lower()
    base = ISO_ROOT / f"{slug}_{drug}"
    base.mkdir(parents=True, exist_ok=True)
    own_runs = base / "amrfinder_runs"
    gcache = base / "refseq"
    selected = base / "selected.tsv"
    # reuse cached AMRFinder runs from ANY frozen cohort (read-only) to save time on the rare shared accession;
    # the fresh selection excludes frozen accessions, so most will be genuinely new runs.
    reuse_glob = "data/raw/*_provdisjoint_*/amrfinder_runs"

    # Fresh cohort: exclude EVERY prior accession (frozen cohorts + tuning). exclude_cohort=<my isolated
    # base name>, which is NOT in the manifest, so nothing is un-excluded -> genuinely fresh genomes.
    manifest = build_manifest()
    exclude_prior = prior_accessions(manifest, exclude_cohort=base.name)
    sel = select_disjoint(group, drug, per_class, reuse_glob, selected, exclude_prior)

    applied = []
    for i, (acc, y) in enumerate(sel.items(), 1):
        mt = _run_dir(acc, own_runs, reuse_glob)
        if mt is None:
            print(f"    [{group} {drug}] [{i}/{len(sel)}] {acc} ({'R' if y else 'S'}) AMRFinder ...", flush=True)
            ensure_run(acc, own_runs, gcache, amrfinder_org, reuse_glob)
            mt = _run_dir(acc, own_runs, reuse_glob)
        if mt is None:
            continue
        call = call_resistance(mt / "main.tsv", drug, organism=reg_org)
        applied.append((call["prediction"], y))
    conf = _conf(applied)

    fro = frozen_number(group, drug)
    fresh = {"sens": conf.get("sens"), "spec": conf.get("spec"), "n_scored": conf.get("n_scored")}
    drift = None
    if fro and fro.get("sens") is not None and fresh["sens"] is not None:
        drift = {"d_sens": round(fresh["sens"] - fro["sens"], 3),
                 "d_spec": round((fresh["spec"] or 0) - (fro["spec"] or 0), 3)}
    return {"group": group, "drug": drug, "registry_organism": reg_org,
            "n_selected": len(sel), "fresh": fresh, "frozen": fro, "drift": drift,
            "cohort_disjoint_from_frozen": True, "excluded_prior_accessions": len(exclude_prior)}


def load_checkpoint(path: Path) -> dict[str, dict]:
    if not path.exists():
        return {}
    return {f"{r['group']}|{r['drug']}": r for r in
            (json.loads(ln) for ln in path.read_text(encoding="utf-8").splitlines() if ln.strip())}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--per-class", type=int, default=30)
    ap.add_argument("--smoke", action="store_true", help="per_class=2 sanity run (a few genomes/cell)")
    a = ap.parse_args(argv)
    per_class = 2 if a.smoke else a.per_class

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    ckpt = OUT_DIR / "checkpoint.jsonl"
    done = load_checkpoint(ckpt)
    print(f"[revalidate] {len(CELLS)} SCORED cells, per_class={per_class}, isolated -> {ISO_ROOT}\n"
          f"[revalidate] {len(done)} already checkpointed", flush=True)

    results = list(done.values())
    with open(ckpt, "a", encoding="utf-8") as fh:
        for group, amrf, drug, reg in CELLS:
            key = f"{group}|{drug}"
            if key in done:
                continue
            print(f"[revalidate] === {group} x {drug} ===", flush=True)
            try:
                rec = revalidate_cell(group, amrf, drug, reg, per_class)
            except Exception as e:  # noqa: BLE001 — one cell failure must not abort the sweep
                print(f"    FAILED {group} {drug}: {type(e).__name__}: {e}", flush=True)
                rec = {"group": group, "drug": drug, "error": f"{type(e).__name__}: {e}"}
            fh.write(json.dumps(rec) + "\n")
            fh.flush()
            results.append(rec)
            fr = rec.get("fresh", {}); fz = (rec.get("frozen") or {})
            print(f"    fresh sens={fr.get('sens')} spec={fr.get('spec')} (n={fr.get('n_scored')}) | "
                  f"frozen sens={fz.get('sens')} spec={fz.get('spec')} | drift={rec.get('drift')}", flush=True)

    # summary (SUBDIR only — never the report-card top-level glob)
    (OUT_DIR / "revalidation_summary.json").write_text(json.dumps(
        {"_schema": "scored-cell-revalidation-v1", "date": DATE, "per_class": per_class,
         "note": "SECOND provenance-disjoint cohort per cell (genomes disjoint from the frozen cohort). "
                 "Isolated from all frozen artifacts; NOT a report-card input.",
         "cells": results}, indent=2), encoding="utf-8")

    lines = [f"# Fresh-cohort re-validation of the 10 SCORED cells ({DATE})", "",
             "Each deployed cell re-scored on a SECOND provenance-disjoint cohort whose genomes are DISJOINT "
             "from the frozen validation cohort (excluded via the accession manifest). Frozen artifacts + the "
             "report card are UNTOUCHED — this is an independent confidence check, not a re-baseline.", "",
             "| cell | fresh sens/spec (n) | frozen sens/spec | Δsens | Δspec |",
             "|---|---|---|---:|---:|"]
    for r in results:
        if r.get("error"):
            lines.append(f"| {r['group']} × {r['drug']} | ERROR: {r['error']} | — | — | — |")
            continue
        fr, fz, dr = r["fresh"], (r.get("frozen") or {}), (r.get("drift") or {})
        lines.append(f"| {r['group']} × {r['drug']} | {fr.get('sens')}/{fr.get('spec')} ({fr.get('n_scored')}) | "
                     f"{fz.get('sens')}/{fz.get('spec')} | {dr.get('d_sens', '—')} | {dr.get('d_spec', '—')} |")
    lines += ["", "Frozen numbers from the committed `wiki/provenance_disjoint_validation_*.json`. A large Δ "
              "means the cell does NOT reproduce on fresh genomes — investigate before trusting the headline.",
              "", "Generated by `scripts/revalidate_scored_cells.py`."]
    (OUT_DIR / "revalidation_summary.md").write_text("\n".join(lines), encoding="utf-8")
    print(f"[revalidate] summary -> {OUT_DIR / 'revalidation_summary.md'}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
