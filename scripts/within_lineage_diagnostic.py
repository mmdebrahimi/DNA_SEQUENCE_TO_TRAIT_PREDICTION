"""Within-lineage conditional ranking diagnostic — does the embedding edge survive CONDITIONING on lineage?

The de-confound gate + leakage-safe CV say the falsifier's AUROC isn't a batch artifact, but a held-out
strain can still have same-MLST relatives in train, so a raw +Xpp NT-vs-k-mer gap could partly mean "NT
encodes lineage context richer than bag-of-k-mers" rather than "NT learned the resistance MECHANISM".

This asks the sharper question directly: WITHIN each shared lineage (an MLST carrying BOTH R and S),
across all R/S strain pairs, how often does the model score the R strain above the S strain? That is a
within-lineage concordance (≈ a lineage-conditioned AUC). If NT's advantage persists inside shared
lineages, the signal is mechanism, not lineage. Cheap (reads the persisted scores sidecar; no model rerun).

Input: the `<drug>_falsifier_<date>.scores.json` written by scripts/amr_falsifier.py
  {strain_ids, y_true, mlst, scores: {variant: [per-strain score|null]}, auroc}

Run:  uv run python scripts/within_lineage_diagnostic.py --scores wiki/ciprofloxacin_falsifier_2026-06-05.scores.json
"""
from __future__ import annotations

import argparse
import json
import sys
from collections import defaultdict
from datetime import date
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent


def within_lineage_concordance(scores, y, mlst):
    """Pooled within-MLST R-vs-S concordance: over every (R,S) pair that SHARES an MLST,
    fraction with score(R) > score(S) (ties = 0.5). Returns (concordance, n_pairs, n_shared_lineages)."""
    by_mlst = defaultdict(list)
    for i, m in enumerate(mlst):
        if scores[i] is not None and m not in (None, "None", ""):
            by_mlst[m].append(i)
    wins = ties = total = 0
    shared = 0
    for m, idxs in by_mlst.items():
        r = [i for i in idxs if y[i] == 1]
        s = [i for i in idxs if y[i] == 0]
        if not r or not s:
            continue
        shared += 1
        for ri in r:
            for si in s:
                total += 1
                if scores[ri] > scores[si]:
                    wins += 1
                elif scores[ri] == scores[si]:
                    ties += 1
    if total == 0:
        return float("nan"), 0, shared
    return (wins + 0.5 * ties) / total, total, shared


def permutation_ci(scores, y, mlst, *, n_perm=2000, seed=42):
    """Permute labels WITHIN each shared MLST; null distribution of within-lineage concordance.
    Returns (observed, null_mean, p_value_one_sided_gt)."""
    rng = np.random.default_rng(seed)
    obs, _, _ = within_lineage_concordance(scores, y, mlst)
    if not np.isfinite(obs):
        return obs, float("nan"), float("nan")
    by_mlst = defaultdict(list)
    for i, m in enumerate(mlst):
        if scores[i] is not None and m not in (None, "None", ""):
            by_mlst[m].append(i)
    shared = [idxs for idxs in by_mlst.values()
              if any(y[i] == 1 for i in idxs) and any(y[i] == 0 for i in idxs)]
    null = []
    for _ in range(n_perm):
        yp = list(y)
        for idxs in shared:
            lab = [y[i] for i in idxs]
            rng.shuffle(lab)
            for i, l in zip(idxs, lab):
                yp[i] = l
        c, _, _ = within_lineage_concordance(scores, yp, mlst)
        null.append(c)
    null = np.array(null)
    p = float((null >= obs).mean())
    return obs, float(null.mean()), p


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Within-lineage conditional ranking diagnostic")
    ap.add_argument("--scores", type=Path, required=True, help="<drug>_falsifier_<date>.scores.json")
    ap.add_argument("--output", type=Path, default=None)
    ap.add_argument("--n-perm", type=int, default=2000)
    args = ap.parse_args(argv)

    if not args.scores.exists():
        print(f"PARKED: scores sidecar not found at {args.scores} "
              f"(run scripts/amr_falsifier.py first).", file=sys.stderr)
        return 2
    d = json.loads(args.scores.read_text(encoding="utf-8"))
    y = [int(v) for v in d["y_true"]]
    mlst = d["mlst"]
    drug = d.get("drug", "?")

    rows = []
    for name, sc in d["scores"].items():
        conc, npairs, nshared = within_lineage_concordance(sc, y, mlst)
        obs, null_mean, p = permutation_ci(sc, y, mlst, n_perm=args.n_perm)
        rows.append((name, d["auroc"].get(name, float("nan")), conc, npairs, nshared, null_mean, p))

    nt = max((r for r in rows if r[0].startswith("NT")), key=lambda r: r[1])
    km = next((r for r in rows if r[0] == "k-mer-XGB"), None)
    delta = (nt[2] - km[2]) if km else float("nan")

    today = date.today().isoformat()
    lines = [
        f"# {drug} within-lineage conditional-ranking diagnostic ({today})", "",
        "> Conditions on lineage: within each MLST carrying BOTH R and S, fraction of R/S pairs scored",
        "> R>S (≈ lineage-conditioned AUC). If NT's edge persists here, the signal is mechanism, not lineage.",
        f"> Source scores: `{args.scores.name}`", "",
        f"**Shared lineages used:** {nt[4]} · **within-lineage R/S pairs:** {nt[3]}", "",
        "| Variant | overall AUROC | within-lineage concordance | perm-null mean | p(≥obs) |",
        "|---|---:|---:|---:|---:|",
    ]
    for name, auroc, conc, npairs, nshared, null_mean, p in rows:
        lines.append(f"| {name} | {auroc:.3f} | {conc:.3f} | {null_mean:.3f} | {p:.3f} |")
    lines += [
        "",
        f"**NT-best within-lineage concordance {nt[2]:.3f} vs k-mer {km[2]:.3f} → "
        f"Δ {delta:+.3f}**" if km else "",
        "",
        "## Reading",
        "- concordance ~0.5 = no within-lineage discrimination (signal was lineage, not mechanism).",
        "- concordance >>0.5 with low p = the model ranks R above S EVEN within the same lineage = mechanism.",
        f"- NT edge persists within lineage if NT concordance > k-mer AND NT p is small. "
        f"(Δ {delta:+.3f}, NT p {nt[6]:.3f}).",
        "- Small n_pairs ⇒ low power; treat as directional diagnostic, not a gate.",
    ]
    out = args.output or (ROOT / f"wiki/{drug}_within_lineage_diagnostic_{today}.md")
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("\n".join(lines), encoding="utf-8")
    print(f"[within-lineage] {drug}: NT-best conc {nt[2]:.3f} (p {nt[6]:.3f}) vs k-mer {km[2]:.3f}; "
          f"Δ {delta:+.3f} over {nt[3]} pairs in {nt[4]} lineages -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
