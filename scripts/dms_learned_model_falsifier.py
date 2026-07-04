#!/usr/bin/env python
"""dms_learned_model_falsifier.py - does a LEARNED protein model capture the causal molecular signal?

THE QUESTION (the JEPA/CLIP thesis, at the layer where it can WORK): given a protein, can a learned
variant-effect model "discover which variants drive the measured phenotype"? We test it end-to-end on
dna_decode's OWN cached data: correlate a learned model's per-variant scores (AlphaMissense) against
WET-LAB deep-mutational-scan (DMS) measured effects (ProteinGym), across every joinable assay, with a
SHUFFLED negative control + a pre-registered bar. Then contextualise against the full ProteinGym
leaderboard (217 assays x ~100 learned models).

HONEST SCOPE (integrity rail): this validates the MOLECULAR-phenotype variant-effect direction (protein
function) - exactly where learned representations succeed. It does NOT rescue the COMPLEX-ORGANISMAL
phenotype direction (dna_decode's 0-for-5 de-confounded negative). AlphaMissense is a supervised
structure+sequence learned model (a proxy for the learned-representation family JEPA belongs to), not
JEPA itself; the point is that a LEARNED REPRESENTATION captures causal molecular signal - the green
light + substrate for building a JEPA/CLIP variant on protein DMS.

Pre-registered falsifier:
  PASS iff median |Spearman(AlphaMissense, DMS)| over joinable assays >= 0.30 AND the shuffled control
  median is ~0 (< 0.05). (Literature: AlphaMissense vs DMS ~0.45; strong PLMs ~0.4-0.5.)

Run:  python scripts/dms_learned_model_falsifier.py [--max-assays N]
Data: D:/dna_decode_cache/proteingym (cached; no download, no GPU, no money).
"""
from __future__ import annotations

import csv
import json
import os
import sys

import numpy as np

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

D = "D:/dna_decode_cache/proteingym"
DMS_DIR = f"{D}/pg_dms/DMS_ProteinGym_substitutions"
PASS_BAR = 0.30
SHUFFLE_MAX = 0.05


def spearman(x, y):
    """Spearman rho via rank + Pearson (no scipy dependency)."""
    x = np.asarray(x, float); y = np.asarray(y, float)
    if len(x) < 5:
        return float("nan")
    rx = np.argsort(np.argsort(x)).astype(float)
    ry = np.argsort(np.argsort(y)).astype(float)
    rx -= rx.mean(); ry -= ry.mean()
    denom = np.sqrt((rx * rx).sum() * (ry * ry).sum())
    return float((rx * ry).sum() / denom) if denom > 0 else float("nan")


def load_alphamissense(path=f"{D}/am_pg.tsv"):
    """{accession: {mutant: am_score}} keyed by wt+pos+mut (e.g. 'M1A')."""
    am = {}
    with open(path) as f:
        for ln in f:
            p = ln.rstrip("\n").split("\t")
            if len(p) < 3:
                continue
            acc, mut, score = p[0], p[1], p[2]
            try:
                am.setdefault(acc, {})[mut] = float(score)
            except ValueError:
                continue
    return am


def joinable_assays(am_acc):
    """List of (dms_id, dms_filename, accession) whose protein has AlphaMissense scores."""
    e2a = json.load(open(f"{D}/pg_entryname_to_accession.json"))
    out = []
    with open(f"{D}/pg_reference.csv", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            acc = e2a.get(row["UniProt_ID"]) or e2a.get(row["DMS_id"])
            if acc in am_acc:
                out.append((row["DMS_id"], row["DMS_filename"], acc))
    return out


def load_dms(path):
    """{mutant: DMS_score} for SINGLE mutants only (skip the big mutated_sequence column)."""
    out = {}
    with open(path, encoding="utf-8") as f:
        r = csv.DictReader(f)
        for row in r:
            m = row.get("mutant", "")
            if not m or ":" in m or ";" in m:   # single mutants only
                continue
            try:
                out[m] = float(row["DMS_score"])
            except (ValueError, KeyError, TypeError):
                continue
    return out


def run_assay(am_prot, dms_path):
    dms = load_dms(dms_path)
    xs, ys = [], []
    for mut, ds in dms.items():
        if mut in am_prot:
            xs.append(am_prot[mut]); ys.append(ds)
    if len(xs) < 20:
        return None
    rho = spearman(xs, ys)
    # shuffled negative control (deterministic permutation)
    rng = np.random.default_rng(0)
    yss = np.array(ys); rng.shuffle(yss)
    rho_shuf = spearman(xs, yss)
    return {"n": len(xs), "rho": rho, "rho_shuf": rho_shuf}


def leaderboard_context(models=("ESM2 (650M)", "EVE (ensemble)", "GEMME", "TranceptEVE L", "ESM-1v (ensemble)")):
    """Median Spearman of strong learned models across all 217 ProteinGym DMS assays (field context)."""
    out = {}
    try:
        with open(f"{D}/pg_spearman_dms.csv", encoding="utf-8") as f:
            r = csv.DictReader(f)
            cols = {m: [] for m in models if m in r.fieldnames}
            for row in r:
                for m in cols:
                    try:
                        cols[m].append(float(row[m]))
                    except (ValueError, TypeError):
                        pass
        for m, v in cols.items():
            out[m] = (float(np.median(v)), len(v))
    except Exception as e:  # noqa: BLE001
        out["_error"] = str(e)
    return out


def main(argv=None):
    argv = argv if argv is not None else sys.argv[1:]
    max_assays = 40
    if "--max-assays" in argv:
        max_assays = int(argv[argv.index("--max-assays") + 1])

    am = load_alphamissense()
    assays = joinable_assays(set(am))
    print(f"=== DMS learned-model falsifier (AlphaMissense vs wet-lab DMS) ===")
    print(f"AlphaMissense proteins: {len(am)} | joinable DMS assays: {len(assays)} "
          f"| evaluating up to {max_assays}\n")

    results = []
    for dms_id, fn, acc in assays[:max_assays]:
        path = os.path.join(DMS_DIR, fn)
        if not os.path.exists(path):
            continue
        r = run_assay(am[acc], path)
        if r:
            r["dms_id"] = dms_id
            results.append(r)

    if not results:
        print("FAIL: no assays produced a usable join"); return 1

    abs_rho = [abs(r["rho"]) for r in results]
    abs_shuf = [abs(r["rho_shuf"]) for r in results]
    med, med_shuf = float(np.median(abs_rho)), float(np.median(abs_shuf))

    # show the strongest few, signed (AlphaMissense pathogenic HIGHER -> expect NEGATIVE vs DMS fitness)
    top = sorted(results, key=lambda r: -abs(r["rho"]))[:6]
    print("strongest joins (signed Spearman; negative = pathogenic variants reduce measured fitness):")
    for r in top:
        print(f"  {r['dms_id'][:44]:44s} n={r['n']:5d}  rho={r['rho']:+.3f}  shuffled={r['rho_shuf']:+.3f}")

    print(f"\nAlphaMissense vs DMS over {len(results)} assays: median |Spearman| = {med:.3f}  "
          f"(shuffled control median = {med_shuf:.3f})")

    print("\nfield context - median Spearman across all 217 ProteinGym DMS assays (cached leaderboard):")
    for m, (val, n) in leaderboard_context().items():
        print(f"  {m:22s} {val:.3f}  (n={n} assays)")

    ok = med >= PASS_BAR and med_shuf < SHUFFLE_MAX
    print()
    if ok:
        print(f"PASS: a LEARNED model captures the causal molecular signal (median |rho| {med:.3f} >= "
              f"{PASS_BAR}, shuffled ~0).")
        print("  => The learned-REPRESENTATION thesis WORKS at the molecular-phenotype layer -> green light")
        print("     + substrate for a JEPA/CLIP protein model. Does NOT rescue complex-phenotype (0-for-5).")
    else:
        print(f"FAIL: median |rho| {med:.3f} < {PASS_BAR} or shuffled {med_shuf:.3f} >= {SHUFFLE_MAX}.")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
