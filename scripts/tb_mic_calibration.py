"""CRyPTIC TB MIC calibration — Family B (MIC) extension of the genome world-model plan (2026-07-11).

Family B calibrated HIV fold-change. The user asked to extend it to MIC. This builds coverage-valid
prediction intervals for MEASURED CRyPTIC BMD-MIC (RIF + INH) on the doubling-dilution scale — the
quantitative TB decoder ("MIC within +/- N dilutions with verified 90% coverage").

WHY this is a valid MIC substrate (NOT the closed BV-BRC negative): CRyPTIC MIC is a WET-LAB broth-
microdilution measurement (not the G1-circular XGBoost-from-genome BV-BRC MIC that was rejected). It IS
end-censored on the dilution ladder (RIF 58% / INH 36% `<=`/`>` values) — handled honestly below.

METHOD:
  * Features: per-isolate WHO grade-1/2 determinant PRESENCE, from streaming `VARIANTS.parquet` once and
    reusing the FROZEN matcher `tb_amr._matched_determinants` (genomic SNV/MNV/indel match). Determinants
    present in >= MIN_SUPPORT isolates are kept. Cached to JSON so re-runs skip the 2.9 GB stream.
  * Target: y = log2(MIC) = the doubling-dilution index. Censored `<=X`/`>X` parsed to the bound + a flag.
  * Model: additive ElasticNet (reuse `hiv_epistasis.nested_oof`) -> held-out OOF log2-MIC.
  * CENSORING-AWARE split-conformal (reuse `hiv_quantitative_calibration._conformal_q`): the conformal
    quantile q is set on RESOLVED (uncensored) isolates only (where y is exact); the interval is oof +/- q.
    - coverage_resolved = held-out fraction of RESOLVED isolates whose exact y is in the interval (the honest
      coverage number; averaged over REPEATS shuffles).
    - consistency_censored = fraction of censored isolates whose interval is CONSISTENT with the censoring
      bound (`<=X`: interval lower <= log2(X); `>X`: interval upper >= log2(X)) — uses the censored data honestly.
    - width = 2q in log2 = the interval is MIC x/div 2^q; a WIDE interval honestly says determinant presence
      does not pin the exact MIC rung (conformal stays VALID regardless of point-model quality).

PRE-REGISTERED BAR: CALIBRATED_MIC_INTERVALS iff |coverage_resolved - 0.90| <= COVER_TOL (=0.05) on BOTH
RIF and INH; else MIC_CENSORING_OR_MODEL_DOMINATES (report which). Honest scope: split-conformal gives
MARGINAL coverage; the resolved subset is enriched for mid-ladder MICs (the censored ends are dropped from the
coverage number but scored for consistency). Frozen AMR + TB surfaces READ-only (uses tb_amr as-is).
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from collections import Counter
from datetime import date as _date
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "scripts"))

import hiv_epistasis as he  # noqa: E402  (nested_oof + _r2)
import hiv_quantitative_calibration as qc  # noqa: E402  (_conformal_q)
from dna_decode.data import tb_who_catalogue  # noqa: E402
from dna_decode.organism_rules import tb_amr  # noqa: E402
from scripts.score_tb_cryptic_parquet import (  # noqa: E402
    DEFAULT_DUMP, DEFAULT_REUSE, DRUG_CODE, indel_determinant_targets, load_calls_by_strain,
)

MIN_SUPPORT = 10
COVER_TOL = 0.05
REPEATS = 20
TARGET = 0.90
SEED = 0
_DEFAULT_CACHE = REPO / "data" / "processed" / "tb_mic_features_cache.json"


def parse_mic(s: str):
    """CRyPTIC MIC cell -> (log2_index, censored_kind in {'exact','left','right'}) or None on NA."""
    s = (s or "").strip()
    if not s or s.upper() == "NA":
        return None
    kind = "exact"
    if s.startswith("<="): kind, s = "left", s[2:]
    elif s.startswith("<"): kind, s = "left", s[1:]
    elif s.startswith(">="): kind, s = "right", s[2:]
    elif s.startswith(">"): kind, s = "right", s[1:]
    try:
        v = float(s)
    except ValueError:
        return None
    if v <= 0:
        return None
    return math.log2(v), kind


def build_features(dump: Path, reuse: Path, cache: Path, force=False):
    """{drug: {uid: [det_symbol,...]}} + {drug: {uid: mic_string}} — cached (streams VARIANTS.parquet once)."""
    if cache.exists() and not force:
        d = json.loads(cache.read_text(encoding="utf-8"))
        return d["features"], d["mics"]
    import csv
    tb_who_catalogue.verify_pins()
    dets = {drug: tb_who_catalogue.load_determinants(drug) for drug in DRUG_CODE}
    wanted = {d.pos for drug in dets for d in dets[drug]}
    targets = {}
    for drug in dets:
        targets.update(indel_determinant_targets(dets[drug]))
    print(f"[tb-mic] streaming VARIANTS.parquet for {len(wanted)} determinant positions ...", flush=True)
    calls, indel_hits = load_calls_by_strain(dump / "VARIANTS.parquet", wanted, targets)
    # per-isolate matched determinant SYMBOLS per drug (reuse the frozen matcher)
    features = {drug: {} for drug in DRUG_CODE}
    for drug, ds in dets.items():
        for uid, c in calls.items():
            matched = tb_amr._matched_determinants(c, ds)
            syms = sorted({m.variant for m in matched})   # Determinant.variant = the catalogue key
            if syms:
                features[drug][uid] = syms
    # MIC strings per drug from the reuse table
    mics = {drug: {} for drug in DRUG_CODE}
    with open(reuse, encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            uid = (r.get("UNIQUEID") or "").strip()
            if not uid:
                continue
            for drug, code in DRUG_CODE.items():
                v = (r.get(f"{code}_MIC") or "").strip()
                if v and v.upper() != "NA":
                    mics[drug][uid] = v
    cache.parent.mkdir(parents=True, exist_ok=True)
    cache.write_text(json.dumps({"features": features, "mics": mics}), encoding="utf-8")
    print(f"[tb-mic] cached features for {len(mics['rifampicin'])} RIF / {len(mics['isoniazid'])} INH MIC isolates",
          flush=True)
    return features, mics


def calibrate_drug(feat_syms: dict, mic_strings: dict, seed=SEED):
    """Build determinant-presence X, log2-MIC y (+censoring), OOF ElasticNet, censoring-aware conformal."""
    uids = [u for u in mic_strings if parse_mic(mic_strings[u]) is not None]
    parsed = {u: parse_mic(mic_strings[u]) for u in uids}
    # determinant vocabulary with support (present in >= MIN_SUPPORT isolates that HAVE a MIC)
    counts = Counter()
    for u in uids:
        counts.update(feat_syms.get(u, []))
    vocab = sorted(d for d, n in counts.items() if n >= MIN_SUPPORT)
    vidx = {d: j for j, d in enumerate(vocab)}
    X = np.zeros((len(uids), len(vocab)), float)
    for i, u in enumerate(uids):
        for d in feat_syms.get(u, []):
            j = vidx.get(d)
            if j is not None:
                X[i, j] = 1.0
    y = np.array([parsed[u][0] for u in uids], float)
    kinds = np.array([parsed[u][1] for u in uids])
    if len(uids) < he.N_MIN or float(np.var(y)) < 1e-9 or X.shape[1] == 0:
        return {"n": len(uids), "powered": False, "note": "too few / no features / degenerate MIC"}
    oof = he.nested_oof(X, y, seed)
    resolved = kinds == "exact"
    n_res = int(resolved.sum())
    if n_res < he.N_MIN:
        return {"n": len(uids), "n_resolved": n_res, "powered": False, "note": "too few resolved MICs"}
    res_idx = np.flatnonzero(resolved)
    abs_res = np.abs(y[res_idx] - oof[res_idx])
    rng = np.random.default_rng(seed)
    covs, qs = [], []
    for _ in range(REPEATS):
        perm = rng.permutation(len(res_idx))
        half = len(perm) // 2
        calib, test = res_idx[perm[:half]], res_idx[perm[half:]]
        q = qc._conformal_q(np.abs(y[calib] - oof[calib]), alpha=1 - TARGET)
        covs.append(float(np.mean(np.abs(y[test] - oof[test]) <= q)))
        qs.append(q)
    cover_resolved = float(np.mean(covs)); q = float(np.mean(qs))
    # consistency on censored isolates
    cens_idx = np.flatnonzero(~resolved)
    consistent = 0
    for i in cens_idx:
        lo, hi = oof[i] - q, oof[i] + q
        if kinds[i] == "left" and lo <= y[i]:
            consistent += 1
        elif kinds[i] == "right" and hi >= y[i]:
            consistent += 1
    return {
        "n": len(uids), "n_resolved": n_res, "n_censored": int(len(cens_idx)),
        "n_determinant_features": len(vocab), "powered": True,
        "r2_oof_resolved": round(he._r2(y[res_idx], oof[res_idx]), 4),
        "cover_resolved_90": round(cover_resolved, 4),
        "halfwidth_log2": round(q, 3), "interval_fold_factor": round(2 ** q, 2),
        "consistency_censored": round(consistent / len(cens_idx), 4) if len(cens_idx) else None,
        "calibrated": bool(abs(cover_resolved - TARGET) <= COVER_TOL),
    }


def run(dump: Path, reuse: Path, cache: Path, force=False, seed=SEED):
    features, mics = build_features(dump, reuse, cache, force=force)
    per_drug = {drug: calibrate_drug(features.get(drug, {}), mics.get(drug, {}), seed) for drug in DRUG_CODE}
    powered = [m for m in per_drug.values() if m.get("powered")]
    calibrated_all = powered and all(m["calibrated"] for m in powered)
    verdict = ("CALIBRATED_MIC_INTERVALS" if calibrated_all
               else ("MIC_CENSORING_OR_MODEL_DOMINATES" if powered else "NO_POWERED_DRUGS"))
    return {
        "artifact": "tb_mic_calibration", "schema": "tb-mic-calibration-v1",
        "question": "Do split-conformal intervals on CRyPTIC BMD-MIC (RIF/INH) achieve nominal held-out "
                    "coverage on resolved MICs (a calibrated quantitative TB decoder)?",
        "substrate": "CRyPTIC reuse-table measured BMD-MIC (wet-lab, NOT G1-circular BV-BRC) + WHO grade-1/2 "
                     "determinant presence from VARIANTS.parquet (frozen tb_amr matcher)",
        "prereg": {"MIN_SUPPORT": MIN_SUPPORT, "COVER_TOL": COVER_TOL, "REPEATS": REPEATS, "target": TARGET,
                   "N_MIN": he.N_MIN, "seed": seed,
                   "censoring": "conformal quantile set on RESOLVED MICs only; censored scored for consistency"},
        "verdict": verdict,
        "honest_caveats": [
            "CRyPTIC MIC is end-censored (RIF 58% / INH 36% <=/> values) — coverage is on the RESOLVED subset "
            "(mid-ladder-enriched), censored isolates scored only for consistency.",
            "Features = WHO grade-1/2 determinant presence (the catalogue that defines R/S) -> the point model "
            "explains R/S better than the exact MIC rung; the conformal WIDTH honestly reflects that.",
            "Split-conformal gives MARGINAL (not per-genotype) coverage; measured wet-lab MIC, in-distribution "
            "vs the WHO catalogue (not independent).",
        ],
        "citation": "CRyPTIC Consortium 2022; WHO TB mutation catalogue v2 (2023)",
        "per_drug": per_drug,
    }


def render_md(res, generated):
    L = [f"# CRyPTIC TB MIC calibration — are the MIC prediction intervals honest? ({generated})", "",
         f"**Verdict: {res['verdict']}**", "",
         f"{res['question']} Substrate: {res['substrate']}.", "",
         "`cover_resolved_90` = held-out coverage on RESOLVED (uncensored) MICs (target 0.90). "
         "`interval_fold_factor` = MIC within ×/÷ this factor (2^halfwidth). `consistency_censored` = fraction "
         "of censored isolates whose interval respects the censoring bound.", "",
         "| drug | n | resolved | censored | dets | R2(res) | **cover_90** | MIC ± | consistency | calibrated |",
         "|---|---|---|---|---|---|---|---|---|---|"]
    for drug, m in res["per_drug"].items():
        if not m.get("powered"):
            L.append(f"| {drug} | {m.get('n')} | — | — | — | — | {m.get('note','')} | — | — | — |")
            continue
        L.append(f"| {drug} | {m['n']} | {m['n_resolved']} | {m['n_censored']} | {m['n_determinant_features']} | "
                 f"{m['r2_oof_resolved']} | **{m['cover_resolved_90']}** | ×/÷{m['interval_fold_factor']} | "
                 f"{m['consistency_censored']} | {'YES' if m['calibrated'] else 'no'} |")
    L += ["", "## Honest caveats"] + [f"- {c}" for c in res["honest_caveats"]]
    L += ["", f"Citation: {res['citation']}."]
    return "\n".join(L)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dump-dir", type=Path, default=DEFAULT_DUMP)
    ap.add_argument("--reuse-csv", type=Path, default=DEFAULT_REUSE)
    ap.add_argument("--cache", type=Path, default=_DEFAULT_CACHE)
    ap.add_argument("--force-stream", action="store_true")
    ap.add_argument("--out", type=Path, default=None)
    a = ap.parse_args(argv)
    if not a.cache.exists() and not (a.dump_dir / "VARIANTS.parquet").exists():
        print(f"ERROR: no feature cache and VARIANTS.parquet absent under {a.dump_dir}", file=sys.stderr)
        return 2
    today = _date.today().isoformat()
    res = run(a.dump_dir, a.reuse_csv, a.cache, force=a.force_stream)
    out = a.out or (REPO / "wiki" / f"tb_mic_calibration_{today}.json")
    out.write_text(json.dumps(res, indent=2), encoding="utf-8")
    out.with_suffix(".md").write_text(render_md(res, today), encoding="utf-8")
    print(render_md(res, today))
    print(f"\n[wrote {out} + .md]  verdict={res['verdict']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
