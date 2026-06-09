"""Out-of-sample validation of the calibrated AMR registry on INDEPENDENT (disjoint) cohorts.

The registry (dna_decode/data/calibrated_amr_rules.json) was calibrated IN-SAMPLE (N~30 each). Before any
calibrated config can be promoted from opt-in to a default, it must hold on a DIFFERENT set of strains.
This script, per CALIBRATED organism×drug:

  1. builds a SECOND cohort DISJOINT from the first (excludes the cohort-1 accessions in
     data/raw/<first_slug>/selected.tsv), up to 15R/15S, written to data/raw/<slug>_indep_<drug>/,
  2. downloads + runs AMRFinder (C:-local; reuses the validator's tested ensure_run),
  3. applies the IN-SAMPLE registry rule to the independent cohort via call_resistance(organism=...)
     -> out-of-sample acc/sens/spec (the generalization test),
  4. RE-CALIBRATES on the independent cohort -> does it recover the same (counter, threshold)?
     (config-stability test),
  5. writes wiki/calibrated_registry_independent_validation_<date>.{md,json}.

Robust + restartable: per-organism try/except (one failure doesn't kill the others); AMRFinder runs are
cached so re-invocation resumes. Designed for unattended overnight operation. No money, no D: dependency.

Usage: .venv/Scripts/python.exe scripts/independent_cohort_validate.py
"""
from __future__ import annotations

import json
import re
import traceback
import urllib.request
from datetime import date as _date
from pathlib import Path

from dna_decode.eval.amr_rules import call_resistance, load_calibrated_registry
from dna_decode.eval.calibrate_organism import calibrate, features_from_main_tsv
from scripts.organism_drug_validate import _run_dir, ensure_run, latest_metadata_url

TARGET_PER_CLASS = 15

# (ncbi_group, amrfinder_org, drug, first_cohort_slug, registry_key)
TARGETS = [
    ("Campylobacter", "Campylobacter", "ciprofloxacin", "campylobacter_ciprofloxacin", "Campylobacter|ciprofloxacin"),
    ("Klebsiella", "Klebsiella_pneumoniae", "ciprofloxacin", "klebsiella_cipro", "Klebsiella|ciprofloxacin"),
    ("Salmonella", "Salmonella", "ciprofloxacin", "salmonella_ciprofloxacin", "Salmonella|ciprofloxacin"),
]


def _exclusion_set(first_slug: str) -> set[str]:
    sel = Path(f"data/raw/{first_slug}/selected.tsv")
    out = set()
    if sel.exists():
        for ln in sel.read_text().splitlines():
            if "\t" in ln:
                out.add(ln.split("\t")[0])
    return out


def select_independent_cohort(group: str, drug: str, exclude: set[str], selected: Path) -> dict[str, int]:
    """Pick up to 15R/15S strains for `drug` EXCLUDING `exclude` (cohort-1 accessions). Idempotent."""
    if selected.exists():
        out = {}
        for ln in selected.read_text().splitlines():
            if "\t" in ln:
                a, rs = ln.split("\t"); out[a] = 1 if rs == "R" else 0
        return out
    tok_r, tok_s = f"{drug}=R", f"{drug}=S"
    chosen: dict[str, int] = {}
    with urllib.request.urlopen(latest_metadata_url(group), timeout=400) as resp:
        header = resp.readline().decode().rstrip("\n").split("\t")
        ai = header.index("asm_acc"); pi = header.index("AST_phenotypes")
        seen = set()
        for raw in resp:
            cells = raw.decode("utf-8", "replace").rstrip("\n").split("\t")
            if len(cells) <= max(ai, pi):
                continue
            acc, ast = cells[ai], cells[pi]
            if not acc.startswith("GCA_") or acc in seen or acc in exclude or not ast or ast == "NULL":
                continue
            for kv in ast.split(","):
                if kv == tok_r or kv == tok_s:
                    lab = 1 if kv == tok_r else 0
                    have = sum(1 for v in chosen.values() if v == lab)
                    if have < TARGET_PER_CLASS:
                        chosen[acc] = lab
                    seen.add(acc); break
            if (sum(1 for v in chosen.values() if v == 1) >= TARGET_PER_CLASS and
                    sum(1 for v in chosen.values() if v == 0) >= TARGET_PER_CLASS):
                break
    selected.parent.mkdir(parents=True, exist_ok=True)
    selected.write_text("".join(f"{a}\t{'R' if v else 'S'}\n" for a, v in chosen.items()), encoding="utf-8")
    return chosen


def _conf(pairs):
    tp = sum(1 for p, y in pairs if p == "R" and y == 1); fp = sum(1 for p, y in pairs if p == "R" and y == 0)
    tn = sum(1 for p, y in pairs if p == "S" and y == 0); fn = sum(1 for p, y in pairs if p == "S" and y == 1)
    ab = sum(1 for p, y in pairs if p == "ABSTAIN")
    n = tp + fp + tn + fn
    return {"n_scored": n, "tp": tp, "fp": fp, "tn": tn, "fn": fn, "abstain": ab,
            "acc": round((tp + tn) / n, 3) if n else None,
            "sens": round(tp / (tp + fn), 3) if (tp + fn) else None,
            "spec": round(tn / (tn + fp), 3) if (tn + fp) else None}


def validate_one(group, amr_org, drug, first_slug, reg_key, registry) -> dict:
    slug = group.lower().replace(" ", "_")
    base = Path(f"data/raw/{slug}_indep_{drug}")
    own_runs = base / "amrfinder_runs"; gcache = base / "refseq"
    reuse_glob = f"data/raw/{slug}_*/amrfinder_runs"
    exclude = _exclusion_set(first_slug)
    sel = select_independent_cohort(group, drug, exclude, base / "selected.tsv")
    nR = sum(sel.values()); nS = len(sel) - nR
    print(f"\n=== {group} {drug} INDEPENDENT cohort: {len(sel)} ({nR}R/{nS}S); excluded {len(exclude)} cohort-1 accs ===")
    # disjointness assertion
    overlap = set(sel) & exclude
    assert not overlap, f"independent cohort overlaps cohort-1: {overlap}"
    for i, acc in enumerate(sel, 1):
        if _run_dir(acc, own_runs, reuse_glob):
            continue
        print(f"  [{i}/{len(sel)}] {acc} ...")
        ensure_run(acc, own_runs, gcache, amr_org, reuse_glob)
    # build features + labels for strains with runs
    strains, labels, applied = [], [], []
    reg_cfg = registry.get("rules", {}).get(reg_key, {})
    for acc, y in sel.items():
        rd = _run_dir(acc, own_runs, reuse_glob)
        if not rd:
            continue
        mt = rd / "main.tsv"
        strains.append(features_from_main_tsv(mt, drug)); labels.append("R" if y else "S")
        # apply the IN-SAMPLE registry rule out-of-sample
        c = call_resistance(mt, drug, organism=group)
        applied.append((c["prediction"], y))
    oos = _conf(applied)
    # re-calibrate on the independent cohort -> config stability
    recal = calibrate(strains, labels, drug) if strains else None
    same_config = (recal is not None and reg_cfg
                   and recal.counter == reg_cfg.get("counter")
                   and recal.threshold == reg_cfg.get("threshold"))
    result = {
        "group": group, "drug": drug, "n": len(sel), "nR": nR, "nS": nS,
        "n_with_runs": len(strains), "excluded_cohort1": len(exclude),
        "in_sample_config": {k: reg_cfg.get(k) for k in ("counter", "threshold", "intrinsic_families_excluded", "loo_balanced_accuracy")},
        "out_of_sample_applying_in_sample_rule": oos,
        "recalibrated_on_independent": (None if recal is None else {
            "counter": recal.counter, "threshold": recal.threshold,
            "intrinsic_families_excluded": recal.intrinsic_families_excluded,
            "loo_balanced_accuracy": recal.loo_balanced_accuracy, "verdict": recal.verdict}),
        "config_recovered": bool(same_config),
    }
    print(f"  OOS (in-sample rule on indep): acc={oos['acc']} sens={oos['sens']} spec={oos['spec']} (N={oos['n_scored']})")
    if recal:
        print(f"  re-calibrated: {recal.counter}@{recal.threshold} (in-sample was {reg_cfg.get('counter')}@{reg_cfg.get('threshold')}) -> recovered={same_config}")
    return result


def main() -> int:
    registry = load_calibrated_registry()
    results = []
    for group, amr_org, drug, first_slug, reg_key in TARGETS:
        try:
            results.append(validate_one(group, amr_org, drug, first_slug, reg_key, registry))
        except Exception as e:
            print(f"  !! {group} {drug} FAILED: {type(e).__name__}: {e}")
            traceback.print_exc()
            results.append({"group": group, "drug": drug, "error": f"{type(e).__name__}: {e}"})
    d = _date.today().isoformat()
    # markdown
    md = [f"# Calibrated AMR registry — INDEPENDENT-cohort out-of-sample validation — {d}", "",
          "> The registry was calibrated IN-SAMPLE (N~30). This applies each in-sample rule to a DISJOINT",
          "> second cohort (cohort-1 accessions excluded) + re-calibrates on it. Promotion gate: a config",
          "> that holds out-of-sample (acc & sens >= 0.80) AND is recovered on the independent cohort is",
          "> eligible to become a default; otherwise it stays opt-in. NCBI labels; AMRFinder pinned image.",
          "", "| organism | drug | indep N | in-sample cfg | OOS acc | OOS sens | OOS spec | re-cal cfg | recovered |",
          "|---|---|---:|---|---:|---:|---:|---|---|"]
    for r in results:
        if "error" in r:
            md.append(f"| {r['group']} | {r['drug']} | — | — | ERROR | {r['error']} | | | |"); continue
        ic = r["in_sample_config"]; oos = r["out_of_sample_applying_in_sample_rule"]; rc = r["recalibrated_on_independent"]
        rcs = f"{rc['counter']}@{rc['threshold']}" if rc else "—"
        md.append(f"| {r['group']} | {r['drug']} | {r['n_with_runs']}/{r['n']} | "
                  f"{ic['counter']}@{ic['threshold']} | {oos['acc']} | {oos['sens']} | {oos['spec']} | "
                  f"{rcs} | {'YES' if r['config_recovered'] else 'no'} |")
    md += ["", "## Reading", "- **OOS acc/sens** = the in-sample registry rule applied to strains it was NOT",
           "  calibrated on. High => the calibrated config generalizes (promotion-eligible).",
           "- **recovered** = re-calibrating from scratch on the independent cohort picks the SAME",
           "  (counter, threshold). YES => the config choice is stable, not a cohort-1 artifact.",
           "", "## Honest scope", "Second cohort is disjoint by accession but same label source (NCBI AST).",
           "A held-out NCBI cohort is a stronger test than in-sample but still not a different-lab study.",
           "Promotion of any config from opt-in to default remains a deliberate decision on this evidence."]
    Path("wiki").mkdir(exist_ok=True)
    out_md = Path("wiki") / f"calibrated_registry_independent_validation_{d}.md"
    out_json = Path("wiki") / f"calibrated_registry_independent_validation_{d}.json"
    out_md.write_text("\n".join(md), encoding="utf-8")
    out_json.write_text(json.dumps({"date": d, "results": results}, indent=2), encoding="utf-8")
    print(f"\nWrote {out_md} + {out_json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
