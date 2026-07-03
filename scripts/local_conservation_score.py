"""FROZEN-SPEC local independent-sites conservation score — the project computes its OWN deterministic
conservation, rather than citing ProteinGym's Site-Independent number.

Purist path (a pre-build review's Alternative 2): reproduce the independent-sites baseline with a PINNED spec,
using ProteinGym's OWN per-protein MSA + its OWN redundancy-based sequence weights (so we don't double-correct
redundancy). Target: reproduce ~0.43 median Spearman on FUNCTION assays and correlate per-assay with ProteinGym's
published Site-Independent (a correct-implementation check).

FROZEN SPEC (declared before looking at outcomes; hashed at runtime):
  * MSA: ProteinGym a2m (`MSA_filename`). Match columns = per-sequence chars after stripping lowercase + '.'
    (a2m insert states); the remaining uppercase/'-' are the match-state residues (length = #match columns).
  * Weights: ProteinGym's precomputed per-sequence weights `<UniProt_ID>_theta_<MSA_theta>.npy` (NO extra filter).
  * Frequency: f_j(a) = (sum_seq w * [res_j==a] + LAMBDA) / (sum_seq w * [res_j is a std AA] + 20*LAMBDA),
    over the 20 standard AAs only (gaps excluded from numerator + denominator). LAMBDA = 0.5.
  * Score(wt->mut at protein position p): col = p - MSA_start (0-based); score = log2 f_col(mut) - log2 f_col(wt)
    (delta-log-odds -- cancels family composition, the review's flagged primary). Positions outside
    [MSA_start, MSA_end] ABSTAIN. Higher score = more tolerated -> positive Spearman vs DMS_score.
Data (ProteinGym MSAs + weights + DMS/zero-shot) on D:.
"""
from __future__ import annotations

import hashlib
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

AAS = "ACDEFGHIKLMNPQRSTVWY"
AA_IDX = {a: i for i, a in enumerate(AAS)}
LAMBDA = 0.5
_STRIP = {ord(c): None for c in "abcdefghijklmnopqrstuvwxyz."}   # a2m insert states
_SINGLE = re.compile(r"^([A-Z])(\d+)([A-Z])$")


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()[:16]


def focus_pos2col(a2m_path: Path, msa_start: int) -> dict:
    """Map protein position -> match-column index by walking the FOCUS (first) sequence:
    uppercase = match residue (consumes a protein position AND a match column), lowercase = insert residue
    (position only), '-' = match gap (column only), '.' = insert gap (neither). This is the correct a2m
    coordinate map -- match columns do NOT equal protein positions when the focus carries inserts/trims."""
    focus = []
    started = False
    with open(a2m_path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                if started:
                    break
                started = True
            elif started:
                focus.append(line.strip())
    seq = "".join(focus)
    pos2col = {}
    pp, col = int(msa_start), 0
    for c in seq:
        if c.isupper():
            pos2col[pp] = col; pp += 1; col += 1
        elif c == "-":
            col += 1
        elif c.islower():
            pp += 1
        # '.' -> neither
    return pos2col


def weighted_freqs(a2m_path: Path, weights: np.ndarray | None) -> np.ndarray:
    """Stream the a2m; accumulate weighted per-match-column AA counts -> frequency matrix [ncol x 20].
    weights=None -> uniform (unweighted) frequencies (the redundancy-uncorrected fallback)."""
    counts = None
    seq_i = -1
    cur = []
    with open(a2m_path, encoding="utf-8", errors="replace") as fh:
        def flush(seqchars):
            nonlocal counts, seq_i
            if not seqchars:
                return
            seq_i += 1
            w = 1.0 if weights is None else (float(weights[seq_i]) if seq_i < len(weights) else 1.0)
            match = "".join(seqchars).translate(_STRIP).upper()          # match-state residues
            if counts is None:
                counts = np.zeros((len(match), 20), dtype="float64")
            arr = np.frombuffer(match.encode("ascii", "replace"), dtype=np.uint8)
            n = min(len(arr), counts.shape[0])
            for j in range(n):
                idx = AA_IDX.get(chr(arr[j]))
                if idx is not None:
                    counts[j, idx] += w
        for line in fh:
            if line.startswith(">"):
                flush(cur); cur = []
            else:
                cur.append(line.strip())
        flush(cur)
    denom = counts.sum(axis=1, keepdims=True) + 20 * LAMBDA
    return (counts + LAMBDA) / denom


def score_assay(dms_csv: Path, a2m_path: Path, weights_path: Path, msa_start: int, msa_end: int) -> dict | None:
    if not a2m_path.exists():
        return None
    weighted = weights_path.exists()
    weights = np.load(weights_path) if weighted else None                # None -> unweighted fallback
    freqs = weighted_freqs(a2m_path, weights)
    logf = np.log2(freqs)
    pos2col = focus_pos2col(a2m_path, int(msa_start))                     # correct a2m coordinate map
    dms = pd.read_csv(dms_csv)
    dms = dms[~dms["mutant"].astype(str).str.contains(":")]
    dms["DMS_score"] = pd.to_numeric(dms["DMS_score"], errors="coerce")
    p = dms["mutant"].astype(str).str.extract(_SINGLE)
    dms = dms.assign(wt=p[0], pos=pd.to_numeric(p[1], errors="coerce"), to=p[2]).dropna(subset=["wt", "to", "pos", "DMS_score"])
    scores, keep = [], []
    for wt, pos, to in zip(dms["wt"], dms["pos"].astype(int), dms["to"]):
        col = pos2col.get(pos, -1)
        if 0 <= col < freqs.shape[0] and wt in AA_IDX and to in AA_IDX:
            scores.append(float(logf[col, AA_IDX[to]] - logf[col, AA_IDX[wt]]))
            keep.append(True)
        else:
            keep.append(False)
    d = dms[keep]
    if len(d) < 20:
        return None
    rho = float(spearmanr(scores, d["DMS_score"])[0])
    # ProteinGym's own Site_Independent Spearman on the same rows (reproduction check), if present in the DMS csv
    si_repro = None
    if "Site_Independent" in dms.columns:
        dd = pd.read_csv(dms_csv)
        dd = dd[dd["mutant"].isin(d["mutant"])]
        v = pd.to_numeric(dd["Site_Independent"], errors="coerce")
        y = pd.to_numeric(dd["DMS_score"], errors="coerce")
        m = v.notna() & y.notna()
        if m.sum() >= 20:
            si_repro = round(float(spearmanr(v[m], y[m])[0]), 4)
    return {"n_scored": len(d), "local_conservation_spearman": round(rho, 4),
            "proteingym_site_independent_spearman": si_repro, "weighted": weighted,
            "msa_sha": _sha256(a2m_path), "weights_sha": _sha256(weights_path) if weighted else None,
            "ncol": int(freqs.shape[0]), "n_seqs": int(len(weights)) if weighted else None}


def run(reference_csv: Path, dms_dir: Path, msa_dir: Path, weights_dir: Path, selections=None) -> dict:
    ref = pd.read_csv(reference_csv)
    ref = ref[ref["taxon"].astype(str).str.contains("Human", case=False, na=False)]
    if selections:
        ref = ref[ref["coarse_selection_type"].isin(selections)]
    rows = []
    for _, r in ref.iterrows():
        did = r["DMS_id"]
        dcsv = Path(dms_dir) / f"{did}.csv"
        a2m = Path(msa_dir) / str(r["MSA_filename"])
        wpath = Path(weights_dir) / f"{r['UniProt_ID']}_theta_{r['MSA_theta']}.npy"
        if not dcsv.exists() or not a2m.exists():                       # weights optional (unweighted fallback)
            continue
        res = score_assay(dcsv, a2m, wpath, int(r["MSA_start"]), int(r["MSA_end"]))
        if res:
            res.update({"DMS_id": did, "selection": r["coarse_selection_type"]})
            rows.append(res)
    df = pd.DataFrame(rows)
    out = {"spec": {"lambda": LAMBDA, "formula": "weighted delta-log2-odds; gaps excluded; ProteinGym weights",
                    "n_assays": len(df)}, "assays": rows}
    if len(df):
        for st, g in df.groupby("selection"):
            out.setdefault("by_selection", {})[st] = {
                "n": len(g), "local_conservation": round(float(g["local_conservation_spearman"].median()), 4),
                "proteingym_site_independent": (round(float(g["proteingym_site_independent_spearman"].median()), 4)
                                                if g["proteingym_site_independent_spearman"].notna().any() else None)}
        out["overall_median_local"] = round(float(df["local_conservation_spearman"].median()), 4)
        wdf = df[df["weighted"]]                                         # ProteinGym-weights subset (rigorous)
        out["weighted_subset"] = {
            "n": len(wdf),
            "local_median": round(float(wdf["local_conservation_spearman"].median()), 4) if len(wdf) else None,
            "proteingym_si_median": (round(float(wdf["proteingym_site_independent_spearman"].median()), 4)
                                     if wdf["proteingym_site_independent_spearman"].notna().any() else None)}
        rep = wdf.dropna(subset=["proteingym_site_independent_spearman"])
        if len(rep):
            out["weighted_subset"]["reproduction_delta_median"] = round(float(
                (rep["local_conservation_spearman"] - rep["proteingym_site_independent_spearman"]).abs().median()), 4)
        udf = df[~df["weighted"]]
        out["unweighted_fallback"] = {
            "n": len(udf),
            "local_median": round(float(udf["local_conservation_spearman"].median()), 4) if len(udf) else None}
    return out


def main(argv=None) -> int:
    import argparse
    DG = Path("D:/dna_decode_cache/proteingym")
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--reference", type=Path, default=DG / "pg_reference.csv")
    ap.add_argument("--dms-dir", type=Path, default=DG / "pg_zeroshot")
    ap.add_argument("--msa-dir", type=Path, default=DG / "pg_msa" / "DMS_msa_files")
    ap.add_argument("--weights-dir", type=Path, default=DG / "pg_weights" / "DMS_msa_weights")
    ap.add_argument("--selections", nargs="*", default=["Activity"])
    ap.add_argument("--out", type=Path, default=REPO / "wiki" / "local_conservation_scores.json")
    a = ap.parse_args(argv)
    res = run(a.reference, a.dms_dir, a.msa_dir, a.weights_dir, a.selections)
    a.out.write_text(json.dumps(res, indent=2), encoding="utf-8")
    print(f"assays scored={res['spec']['n_assays']} | overall median local={res.get('overall_median_local')}")
    print(f"reproduction |local - ProteinGym-SI| median={res.get('reproduction_delta_median')}")
    for st, v in res.get("by_selection", {}).items():
        print(f"  {st:14} n={v['n']:>2} local={v['local_conservation']} proteingym_SI={v['proteingym_site_independent']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
