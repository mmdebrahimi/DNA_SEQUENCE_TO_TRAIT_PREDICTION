"""Where is the k-mer ceiling, and is the committed order-3 chain a STRONG null or a weak one?

WHY THIS EXISTS. The committed scorer calls its order-3 chain "the STRONG 4-mer null" and every
beyond-k-mer verdict in this project is stated against it. Published practice disagrees: the one
DNA-language-model evaluation that reports a Markov baseline at all (BEND 2.0, arXiv 2611.03806)
uses a TUNED chain of **order 10 with smoothing strength 64**, selected on a held-out split. If the
k-mer ceiling actually sits well below order 3's cross-entropy, then "the model loses to a 4-mer
chain" understates the problem -- the model is losing to a weak baseline and the real k-mer ceiling
is further away still.

THE BAR IS FROZEN HERE, BEFORE ANY NUMBER EXISTS, and it is deliberately two-sided:

* WEAK_NULL  -- the best tuned order beats order 3 by >= 0.0100 nats held-out. Then order 3 is NOT
  the k-mer ceiling, and every existing "beats / loses to order-3" verdict must be restated against
  the tuned chain before it means anything.
* AT_PLATEAU -- the best tuned order is within 0.0100 nats of order 3. Then order 3 sits at the
  k-mer plateau, the project's framing stands, and the published order-10 convention simply does
  not bind at this window length.

0.0100 nats is one third of the ~0.03-nat total margin this task has above chance, so it is a
difference large enough to change a reading rather than a difference merely detectable.

HONEST LIMITS, stated here rather than discovered later.
* This sweeps a HOMOGENEOUS chain only. The companion probe (glm_frame_headroom_probe.py) measures
  what a positional variable adds; the two are separate axes and are not combined here.
* Smoothing is tuned per order on a disjoint TUNE split, never on the scored split. Without that,
  high orders win in-sample by construction and the sweep would be meaningless.
* Promoter windows are rebuilt here as the 150 bp upstream of each annotated gene start,
  strand-aware. That is the standard definition and matches the committed corpus geometry
  (window_length_bp 150), but it is REBUILT, not the committed corpus's own window list -- so the
  absolute cross-entropies are not expected to match the committed artifacts to the last digit. The
  ORDER COMPARISON is internal to this run and unaffected.
* An unseen context backs off to UNIFORM, never to zero probability.
"""
from __future__ import annotations

import argparse
import json
import math
import random
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

MG1655 = ROOT / "data" / "cache" / "refseq" / "GCF_000005845.2"

BASES = "ACGT"
_IDX = {b: i for i, b in enumerate(BASES)}
LN4 = math.log(4)
WINDOW = 150

#: Frozen before any number existed. See the module docstring.
PREREGISTERED = {
    "margin_nats": 0.0100,
    "verdicts": {
        "WEAK_NULL": "best tuned order beats order 3 by >= 0.0100 nats held-out",
        "AT_PLATEAU": "best tuned order is within 0.0100 nats of order 3",
    },
    "consequence_if_weak_null": ("every existing beats/loses-to-order-3 verdict must be restated "
                                 "against the tuned chain before it means anything"),
    "reference": "BEND 2.0 arXiv 2611.03806 uses a tuned order-10 chain with smoothing 64",
}

ALPHAS = (0.01, 0.1, 0.5, 1.0, 4.0, 16.0, 64.0)


class SweepError(RuntimeError):
    pass


_COMP = str.maketrans("ACGTN", "TGCAN")


def revcomp(s: str) -> str:
    return s.translate(_COMP)[::-1]


def read_fasta(path: Path) -> dict[str, str]:
    seqs: dict[str, str] = {}
    name = None
    buf: list[str] = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith(">"):
            if name is not None:
                seqs[name] = "".join(buf).upper()
            name = line[1:].split()[0]
            buf = []
        else:
            buf.append(line.strip())
    if name is not None:
        seqs[name] = "".join(buf).upper()
    return seqs


def promoter_windows(gff: Path, fna: Path, *, width: int = WINDOW) -> list[str]:
    """The `width` bp immediately upstream of each annotated gene start, strand-aware."""
    contigs = read_fasta(fna)
    if not contigs:
        raise SweepError(f"no sequence parsed from {fna}")
    out: list[str] = []
    for line in gff.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line or line.startswith("#"):
            continue
        f = line.split("\t")
        if len(f) < 9 or f[2] != "gene":
            continue
        contig, start, end, strand = f[0], int(f[3]), int(f[4]), f[6]
        seq = contigs.get(contig)
        if seq is None:
            continue
        if strand == "+":
            lo = start - 1 - width
            if lo < 0:
                continue
            win = seq[lo:start - 1]
        else:
            win = seq[end:end + width]
            if len(win) < width:
                continue
            win = revcomp(win)
        if len(win) != width or any(c not in _IDX for c in win):
            continue
        out.append(win)
    if not out:
        raise SweepError("no promoter windows built")
    return out


def fit(seqs: list[str], order: int) -> dict[str, list[int]]:
    table: dict[str, list[int]] = {}
    for s in seqs:
        for i in range(order, len(s)):
            ctx = s[i - order:i]
            row = table.get(ctx)
            if row is None:
                row = [0, 0, 0, 0]
                table[ctx] = row
            row[_IDX[s[i]]] += 1
    return table


def cross_entropy(table: dict[str, list[int]], seqs: list[str], order: int, alpha: float) -> tuple[float, int]:
    tot = 0.0
    n = 0
    for s in seqs:
        for i in range(order, len(s)):
            row = table.get(s[i - order:i])
            if row is None:
                tot += LN4
            else:
                tot += -math.log((row[_IDX[s[i]]] + alpha) / (sum(row) + 4.0 * alpha))
            n += 1
    if n == 0:
        raise SweepError("no scoreable positions")
    return tot / n, n


def three_way_split(seqs: list[str], seed: int) -> tuple[list[str], list[str], list[str]]:
    """FIT / TUNE / SCORE. Smoothing is chosen on TUNE and never on SCORE."""
    idx = list(range(len(seqs)))
    random.Random(seed).shuffle(idx)
    a = len(idx) * 50 // 100
    b = len(idx) * 75 // 100
    return ([seqs[i] for i in idx[:a]],
            [seqs[i] for i in idx[a:b]],
            [seqs[i] for i in idx[b:]])


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--max-order", type=int, default=8)
    ap.add_argument("--max-windows", type=int, default=None)
    ap.add_argument("--out", default=None)
    a = ap.parse_args()

    gff, fna = MG1655 / "annotations.gff3", MG1655 / "genome.fna"
    for p in (gff, fna):
        if not p.exists():
            print(f"REFUSING: missing input {p}", file=sys.stderr)
            return 3

    wins = promoter_windows(gff, fna)
    if a.max_windows:
        wins = wins[:a.max_windows]
    fit_s, tune_s, score_s = three_way_split(wins, a.seed)
    print(f"MG1655 promoter windows: {len(wins)} total "
          f"(fit {len(fit_s)} / tune {len(tune_s)} / score {len(score_s)})", flush=True)

    rows: list[dict] = []
    for order in range(1, a.max_order + 1):
        table = fit(fit_s, order)
        best_alpha, best_tune = None, None
        for alpha in ALPHAS:
            ce, _ = cross_entropy(table, tune_s, order, alpha)
            if best_tune is None or ce < best_tune:
                best_alpha, best_tune = alpha, ce
        score_ce, n = cross_entropy(table, score_s, order, best_alpha)
        rows.append({"order": order, "tuned_alpha": best_alpha,
                     "tune_ce_nats": round(best_tune, 6),
                     "score_ce_nats": round(score_ce, 6),
                     "margin_below_uniform_nats": round(LN4 - score_ce, 6),
                     "n_contexts_fitted": len(table), "n_scored_positions": n})
        print(f"  order {order}: alpha {best_alpha:<5} score CE {score_ce:.6f} "
              f"(margin {LN4 - score_ce:+.6f}, {len(table)} contexts)", flush=True)

    by_order = {r["order"]: r for r in rows}
    order3 = by_order.get(3)
    if order3 is None:
        raise SweepError("order 3 not swept -- the comparison is undefined")
    best = min(rows, key=lambda r: r["score_ce_nats"])
    delta = order3["score_ce_nats"] - best["score_ce_nats"]
    verdict = "WEAK_NULL" if delta >= PREREGISTERED["margin_nats"] else "AT_PLATEAU"

    res = {
        "record": "glm-markov-order-sweep-v1",
        "question": "where is the k-mer ceiling for MG1655 promoter windows, and is order 3 at it?",
        "PREREGISTERED": PREREGISTERED,
        "geometry": {"window_bp": WINDOW, "n_windows": len(wins),
                     "split": "50 fit / 25 tune / 25 score, sequence-level, seeded"},
        "uniform_ce_nats": round(LN4, 6),
        "alphas_searched": list(ALPHAS),
        "sweep": rows,
        "order3": order3,
        "best": best,
        "delta_order3_minus_best_nats": round(delta, 6),
        "VERDICT": verdict,
        "reading": (f"order 3 scores {order3['score_ce_nats']:.6f}; the best tuned order is "
                    f"{best['order']} at {best['score_ce_nats']:.6f}; difference "
                    f"{delta:+.6f} nats against a frozen bar of {PREREGISTERED['margin_nats']}"),
        "honest_limits": [
            "Homogeneous chains only -- a positional variable is a separate axis, measured by "
            "scripts/glm_frame_headroom_probe.py.",
            "Promoter windows are REBUILT as 150 bp upstream of annotated gene starts, not the "
            "committed corpus's own window list, so absolute CE is not expected to match the "
            "committed artifacts digit-for-digit. The order comparison is internal to this run.",
            "Smoothing is tuned on a disjoint TUNE split; an unseen context backs off to uniform.",
            "This bounds what K-MER STATISTICS can reach. It says nothing about any model, and it "
            "is not evidence the diffusion generator can or cannot beat it.",
        ],
    }

    out = Path(a.out) if a.out else ROOT / "wiki" / "glm_markov_order_sweep_2026-11-25.json"
    out.write_text(json.dumps(res, indent=1) + "\n", encoding="utf-8")
    print(f"\nVERDICT {verdict} | order3 {order3['score_ce_nats']:.6f} vs best "
          f"(order {best['order']}) {best['score_ce_nats']:.6f} | delta {delta:+.6f}")
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
