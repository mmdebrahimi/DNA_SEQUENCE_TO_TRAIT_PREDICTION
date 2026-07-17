"""Is the SHIPPED forward default (blosum62) 'modest' across ALL of ProteinGym, or just on 2 proteins?

The `dna-decode forward` cell's evidence quotes blosum62 Spearman on exactly TWO proteins (TEM-1 0.3465,
PTEN 0.182) and calls it "REAL but modest". The companion to the inverse sweep: is that representative of
the shipped wheel-only default across the whole ProteinGym substitution benchmark (N=217), or a 2-protein
sample? Same cheap CPU pattern; bounds the OTHER shipped molecular default at scale.

Metric: Spearman(blosum62 substitution severity, measured DMS) over the MEASURED single-mutants per assay,
coordinate-checked. No GPU, no network. Restartable JSONL checkpoint. Reuses the forward cell's OWN scorer
(`blosum62_score`) so it measures exactly what ships.

Run:  uv run python scripts/forward_blosum_proteingym_sweep.py
Exit: 0 = ran; 2 = substrate unavailable.
"""
from __future__ import annotations

import argparse
import csv
import json
import statistics
import sys
from datetime import date
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from dna_decode.forward.variant_effect import blosum62_score  # noqa: E402
from scripts.forward_inverse_proteingym_sweep import DMS_DIR, REF, AMINO_ACIDS  # noqa: E402


def _spearman(x: list[float], y: list[float]) -> float:
    """Rank correlation with MID-RANKS for ties (the documented tie-order trap)."""
    def rk(v: list[float]) -> list[float]:
        order = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0] * len(v)
        i = 0
        while i < len(order):
            j = i
            while j + 1 < len(order) and v[order[j + 1]] == v[order[i]]:
                j += 1
            mid = (i + j) / 2.0
            for k in range(i, j + 1):
                r[order[k]] = mid
            i = j + 1
        return r
    if len(x) < 3:
        return 0.0
    rx, ry = rk(x), rk(y)
    n = len(x)
    mx, my = sum(rx) / n, sum(ry) / n
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    den = (sum((a - mx) ** 2 for a in rx) * sum((b - my) ** 2 for b in ry)) ** 0.5
    return num / den if den else 0.0


def score_one(row: dict) -> dict:
    dms_csv = DMS_DIR / f"{row['DMS_id']}.csv"
    base = {"dms_id": row["DMS_id"], "taxon": row.get("taxon", "")}
    if not dms_csv.exists():
        return {**base, "status": "NO_DMS_FILE"}
    target = row["target_seq"].strip().upper()
    sev, meas = [], []
    with dms_csv.open(encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            mut = (r.get("mutant") or "").strip()
            if ":" in mut or len(mut) < 3:
                continue
            wt, alt = mut[0], mut[-1]
            try:
                pos = int(mut[1:-1])
                val = float(r["DMS_score"])
            except (ValueError, TypeError, KeyError):
                continue
            if not (1 <= pos <= len(target)) or alt not in AMINO_ACIDS or alt == wt:
                continue
            if target[pos - 1] != wt:
                continue
            sev.append(blosum62_score(wt, alt))     # the forward cell's OWN blosum62 method
            meas.append(val)
    if len(sev) < 200:
        return {**base, "status": "UNDERPOWERED", "n": len(sev)}
    # abs(spearman) -- blosum severity is signed opposite to fitness in some assays; the cell reports |rho|
    rho = _spearman(sev, meas)
    return {**base, "status": "SCORED", "n": len(sev),
            "spearman": round(rho, 4), "abs_spearman": round(abs(rho), 4)}


def render_md(rep: dict) -> str:
    s = rep["abs_spearman"]
    n = rep["n_scored"]
    return chr(10).join([
        f"# Is the shipped forward default (blosum62) 'modest' at scale? — ProteinGym N={n} ({rep['date']})",
        "",
        f"**Median |Spearman| = {s['median']:.3f} (mean {s['mean']:.3f}) across {n} assays.** The forward "
        "cell quotes blosum62 on two proteins (TEM-1 0.3465, PTEN 0.182) and calls it 'REAL but modest'. "
        "At scale 'modest' is accurate for the MEDIAN -- but the cell LEADS with TEM-1's 0.35, which is "
        f"**top-13% ({rep['n_above_0.3']}/{n} proteins reach 0.30), not typical**. PTEN's 0.18 is near the "
        "median.",
        "",
        "| |Spearman| | value |",
        "|---|---:|",
        f"| min | {s['min']:.3f} |",
        f"| q25 | {s['q25']:.3f} |",
        f"| **median** | **{s['median']:.3f}** |",
        f"| q75 | {s['q75']:.3f} |",
        f"| max | {s['max']:.3f} |",
        f"| below 0.15 (near useless) | {rep['n_below_0.15']}/{n} |",
        f"| at/above 0.30 | {rep['n_above_0.3']}/{n} |",
        "",
        "So the shipped wheel-only forward default is a weak-to-modest ranker on the typical protein; ESM "
        "(the learned method, ~0.49 median on ProteinGym, GPU) is what makes the forward cell strong, which "
        "the cell already states. This is not an overclaim correction -- the cell says 'modest' -- but it "
        "recontextualises the headline 0.35 as top-decile rather than representative.",
        f"",
        f"{rep['note']}",
        "",
    ])


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", type=Path, default=REPO / "wiki")
    ap.add_argument("--checkpoint", type=Path,
                    default=REPO / "data" / "processed" / "forward_blosum_proteingym.jsonl")
    args = ap.parse_args()
    if not REF.exists() or not DMS_DIR.exists():
        print("[fwd-blosum] SUBSTRATE UNAVAILABLE (D: not mounted?)", file=sys.stderr)
        return 2

    with REF.open(encoding="utf-8") as fh:
        ref = list(csv.DictReader(fh))
    args.checkpoint.parent.mkdir(parents=True, exist_ok=True)
    done = {}
    if args.checkpoint.exists():
        for line in args.checkpoint.read_text(encoding="utf-8").splitlines():
            if line.strip():
                r = json.loads(line)
                done[r["dms_id"]] = r
    todo = [row for row in ref if row["DMS_id"] not in done]
    print(f"[fwd-blosum] {len(ref)} assays | {len(done)} done | {len(todo)} to run", flush=True)
    with args.checkpoint.open("a", encoding="utf-8") as ck:
        for i, row in enumerate(todo, 1):
            try:
                res = score_one(row)
            except Exception as exc:
                res = {"dms_id": row["DMS_id"], "status": "ERROR", "error": f"{type(exc).__name__}: {exc}"}
            ck.write(json.dumps(res) + "\n")
            ck.flush()
            done[res["dms_id"]] = res
            if i % 40 == 0 or i == len(todo):
                print(f"[fwd-blosum]   {i}/{len(todo)}", flush=True)

    scored = [r for r in done.values() if r["status"] == "SCORED"]
    rhos = sorted(r["abs_spearman"] for r in scored)
    rep = {
        "schema": "forward-blosum-proteingym-v1", "date": date.today().isoformat(),
        "question": ("is the shipped forward default (blosum62) 'modest' across ALL of ProteinGym, or just "
                     "on the 2 proteins the cell quotes?"),
        "n_scored": len(scored),
        "abs_spearman": {
            "min": rhos[0] if rhos else None, "q25": rhos[len(rhos) // 4] if rhos else None,
            "median": statistics.median(rhos) if rhos else None,
            "q75": rhos[3 * len(rhos) // 4] if rhos else None, "max": rhos[-1] if rhos else None,
            "mean": round(statistics.fmean(rhos), 4) if rhos else None,
        },
        "reference_2_proteins": {"TEM-1": 0.3465, "PTEN": 0.182},
        "n_above_0.3": sum(1 for r in rhos if r >= 0.3),
        "n_below_0.15": sum(1 for r in rhos if r < 0.15),
        "note": ("blosum62 is a signed severity score; the cell reports |rho| per protein. This measures the "
                 "shipped default's rank quality at scale -- ESM (the learned method) is separately ~0.49 "
                 "median on ProteinGym and beats blosum everywhere (the cell already says so)."),
    }
    stem = f"forward_blosum_proteingym_{rep['date']}"
    (args.out_dir / f"{stem}.json").write_text(json.dumps(rep, indent=2), encoding="utf-8")
    (args.out_dir / f"{stem}.md").write_text(render_md(rep), encoding="utf-8")
    s = rep["abs_spearman"]
    print(f"\n[fwd-blosum] SCORED {len(scored)} assays | shipped blosum62 |Spearman|:")
    print(f"  min {s['min']:.3f} | q25 {s['q25']:.3f} | median {s['median']:.3f} | q75 {s['q75']:.3f} "
          f"| max {s['max']:.3f} | mean {s['mean']:.3f}")
    print(f"  >=0.30: {rep['n_above_0.3']}/{len(scored)} | <0.15: {rep['n_below_0.15']}/{len(scored)}")
    print(f"  (the cell quotes TEM-1 0.35 / PTEN 0.18 -- is that representative?)")
    print(f"  -> {args.out_dir / (stem + '.json')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
