"""How much structure can a HOMOGENEOUS order-3 Markov chain NOT represent?

THE QUESTION THIS SETTLES. Two shipped hypotheses rest on premises nobody has measured:

* The codon-frame positive control assumes protein-coding DNA carries frame-periodic structure a
  homogeneous order-3 chain cannot encode (it conditions on the previous 3 bases regardless of
  reading frame). If that gap is near zero, the control is dead on its own terms -- a frame-aware
  model would have nothing to win -- and the whole "is it the model or the substrate" diagnostic
  evaporates before the domain-shift objection is even reached.
* The position-concentrated-scoring hypothesis assumes promoter windows carry POSITION-specific
  structure that uniform position sampling dilutes. Same shape of question, different conditioning
  variable.

Both reduce to one quantity: the conditional mutual information between the base and a
POSITIONAL variable, GIVEN the 3-base left context. That is exactly the information a homogeneous
order-3 chain is structurally blind to, because its sufficient statistic is the 4-mer conditional
and carries no position index.

    gap = H(base | prev3) - H(base | prev3, positional_var)

WHY HELD-OUT ESTIMATION IS NOT OPTIONAL. Conditioning on more variables NEVER increases training
entropy, so an in-sample gap is positive by construction and measures nothing. Every number here is
a held-out cross-entropy under add-alpha smoothing fitted on a disjoint train half. A gap that
survives held-out estimation is real structure; one that does not is overfitting. This is the same
discipline the project applies elsewhere and the reason a naive version of this probe would have
produced a confident, meaningless positive.

WHAT IT DOES NOT SHOW. This measures what a CHAIN cannot represent. It does NOT show the diffusion
model captures it, and it is not a claim about any model -- it bounds the PRIZE available to a
positionally-aware model, nothing more. A large gap means the diagnostic is worth running; it does
not mean it will pass.

THE COMPARISON IS PAIRED BY CONSTRUCTION. Both substrates are scored by the same estimator, the
same smoothing, the same split rule and the same held-out protocol, so the CDS-vs-promoter contrast
cannot be an artifact of differing method.
"""
from __future__ import annotations

import argparse
import json
import math
import random
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

MG1655 = ROOT / "data" / "cache" / "refseq" / "GCF_000005845.2"
PROMOTER_CORPUS = ROOT / "wiki" / "glm_training_corpus_representative_2026-11-22.txt"

ORDER = 3           # pinned to the committed scorer's chain order -- the whole point is this chain
ALPHA = 0.5         # add-alpha (Jeffreys) smoothing; reported, not tuned
BASES = "ACGT"
_IDX = {b: i for i, b in enumerate(BASES)}
LN4 = math.log(4)


class ProbeError(RuntimeError):
    pass


# ---------------------------------------------------------------------------- sequence loading

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


_COMP = str.maketrans("ACGTN", "TGCAN")


def revcomp(s: str) -> str:
    return s.translate(_COMP)[::-1]


def load_cds(gff: Path, fna: Path, *, min_len: int = 300) -> list[str]:
    """In-frame CDS sequences, strand-corrected so index 0 IS codon position 0.

    Only phase-0, length-divisible-by-3, single-segment CDS rows are kept. A JOINED CDS (the
    programmed-frameshift shape this repo has documented) would have two rows sharing an id; those
    are dropped rather than concatenated, because a wrong join would silently corrupt the frame --
    which is the one variable this probe measures.
    """
    contigs = read_fasta(fna)
    if not contigs:
        raise ProbeError(f"no sequence parsed from {fna}")
    seen: dict[str, int] = defaultdict(int)
    rows: list[tuple[str, int, int, str, str]] = []
    for line in gff.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line or line.startswith("#"):
            continue
        f = line.split("\t")
        if len(f) < 9 or f[2] != "CDS":
            continue
        attrs = f[8]
        key = ""
        for part in attrs.split(";"):
            if part.startswith("ID="):
                key = part[3:]
                break
        seen[key] += 1
        rows.append((f[0], int(f[3]), int(f[4]), f[6], key))

    out: list[str] = []
    for contig, start, end, strand, key in rows:
        if key and seen[key] > 1:
            continue                     # joined / multi-segment CDS -> frame not reconstructable here
        seq = contigs.get(contig)
        if seq is None:
            continue
        sub = seq[start - 1:end]
        if len(sub) < min_len or len(sub) % 3 != 0:
            continue
        if strand == "-":
            sub = revcomp(sub)
        if any(c not in _IDX for c in sub):
            continue
        out.append(sub)
    if not out:
        raise ProbeError("no usable in-frame CDS extracted")
    return out


def load_promoters(path: Path, *, limit: int | None = None) -> list[str]:
    seqs = [ln.strip().upper() for ln in path.read_text(encoding="utf-8").splitlines() if ln.strip()]
    seqs = [s for s in seqs if all(c in _IDX for c in s)]
    if not seqs:
        raise ProbeError(f"no usable promoter windows in {path}")
    if limit is not None:
        seqs = seqs[:limit]
    return seqs


# ---------------------------------------------------------------------------- the estimator

def _fit(seqs: list[str], posfn) -> dict:
    """Count table keyed (context, positional_value). posfn(i, seq) -> hashable, or None to pool."""
    table: dict[tuple, list[int]] = {}
    for s in seqs:
        for i in range(ORDER, len(s)):
            ctx = s[i - ORDER:i]
            pv = posfn(i, s)
            key = (ctx, pv)
            row = table.get(key)
            if row is None:
                row = [0, 0, 0, 0]
                table[key] = row
            row[_IDX[s[i]]] += 1
    return table


def _cross_entropy(table: dict, seqs: list[str], posfn, *, alpha: float = ALPHA) -> tuple[float, int]:
    """Held-out mean negative log-likelihood in nats. An unseen cell backs off to UNIFORM, never zero."""
    tot = 0.0
    n = 0
    for s in seqs:
        for i in range(ORDER, len(s)):
            ctx = s[i - ORDER:i]
            key = (ctx, posfn(i, s))
            row = table.get(key)
            if row is None:
                tot += LN4
            else:
                denom = sum(row) + 4.0 * alpha
                tot += -math.log((row[_IDX[s[i]]] + alpha) / denom)
            n += 1
    if n == 0:
        raise ProbeError("no scoreable positions")
    return tot / n, n


def split(seqs: list[str], seed: int) -> tuple[list[str], list[str]]:
    idx = list(range(len(seqs)))
    random.Random(seed).shuffle(idx)
    cut = len(idx) // 2
    return [seqs[i] for i in idx[:cut]], [seqs[i] for i in idx[cut:]]


def measure(seqs: list[str], posfn, label: str, seed: int) -> dict:
    tr, te = split(seqs, seed)
    tbl = _fit(tr, posfn)
    ce, n = _cross_entropy(tbl, te, posfn)
    return {"label": label, "held_out_ce_nats": round(ce, 6), "n_predictions": n,
            "n_train_seqs": len(tr), "n_test_seqs": len(te), "n_cells_fitted": len(tbl),
            "margin_below_uniform_nats": round(LN4 - ce, 6)}


# ---------------------------------------------------------------------------- positional variables

def pos_pooled(i, s):          # the homogeneous chain: no positional variable at all
    return None


def pos_frame(i, s):           # codon position -- defined ONLY for in-frame CDS
    return i % 3


def pos_window(i, s):          # absolute index in a fixed-length promoter window
    return i


def pos_window_bucket10(i, s):  # coarse position -- fewer cells, far less overfitting pressure
    return (i * 10) // max(len(s), 1)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--max-cds", type=int, default=2000)
    ap.add_argument("--max-promoters", type=int, default=20000)
    ap.add_argument("--out", default=None)
    a = ap.parse_args()

    gff = MG1655 / "annotations.gff3"
    fna = MG1655 / "genome.fna"
    for p in (gff, fna, PROMOTER_CORPUS):
        if not p.exists():
            print(f"REFUSING: missing input {p}", file=sys.stderr)
            return 3

    cds = load_cds(gff, fna)
    random.Random(a.seed).shuffle(cds)
    cds = cds[:a.max_cds]
    prom = load_promoters(PROMOTER_CORPUS, limit=a.max_promoters)
    print(f"CDS: {len(cds)} in-frame sequences | promoters: {len(prom)} windows", flush=True)

    res: dict = {
        "record": "glm-frame-headroom-probe-v1",
        "question": ("how much structure can a HOMOGENEOUS order-3 Markov chain NOT represent, "
                     "on coding sequence (frame) vs promoter windows (window position)?"),
        "order": ORDER, "alpha": ALPHA, "seed": a.seed,
        "uniform_ce_nats": round(LN4, 6),
        "estimator": ("held-out cross-entropy, 50/50 sequence-level split, add-alpha smoothing, "
                      "unseen cell backs off to uniform; identical protocol for both substrates"),
    }

    # --- CDS: does frame-awareness buy anything a homogeneous order-3 chain cannot get?
    cds_h = measure(cds, pos_pooled, "cds_homogeneous_order3", a.seed)
    cds_f = measure(cds, pos_frame, "cds_frame_aware_order3", a.seed)
    res["cds"] = {"homogeneous": cds_h, "frame_aware": cds_f,
                  "gap_nats": round(cds_h["held_out_ce_nats"] - cds_f["held_out_ce_nats"], 6)}
    print(f"  CDS  homogeneous {cds_h['held_out_ce_nats']:.4f} | frame-aware "
          f"{cds_f['held_out_ce_nats']:.4f} | gap {res['cds']['gap_nats']:+.4f}", flush=True)

    # --- promoters: does position-awareness buy anything? (two granularities)
    pr_h = measure(prom, pos_pooled, "promoter_homogeneous_order3", a.seed)
    pr_p = measure(prom, pos_window, "promoter_position_aware_order3", a.seed)
    pr_b = measure(prom, pos_window_bucket10, "promoter_bucket10_order3", a.seed)
    res["promoter"] = {"homogeneous": pr_h, "position_aware": pr_p, "bucket10": pr_b,
                       "gap_nats_full_position": round(pr_h["held_out_ce_nats"] - pr_p["held_out_ce_nats"], 6),
                       "gap_nats_bucket10": round(pr_h["held_out_ce_nats"] - pr_b["held_out_ce_nats"], 6)}
    print(f"  PROM homogeneous {pr_h['held_out_ce_nats']:.4f} | position-aware "
          f"{pr_p['held_out_ce_nats']:.4f} | gap {res['promoter']['gap_nats_full_position']:+.4f}", flush=True)
    print(f"  PROM bucket10 {pr_b['held_out_ce_nats']:.4f} | gap "
          f"{res['promoter']['gap_nats_bucket10']:+.4f}", flush=True)

    # --- the headroom contrast that actually decides the ranking
    res["headroom"] = {
        "cds_total_margin_below_uniform_nats": cds_h["margin_below_uniform_nats"],
        "promoter_total_margin_below_uniform_nats": pr_h["margin_below_uniform_nats"],
        "reading": ("the margin is how much an order-3 chain itself extracts; the gap is what a "
                    "positionally-aware model could add ON TOP. A substrate with a large margin AND "
                    "a large gap is where a positional model has room to win."),
    }

    out = Path(a.out) if a.out else ROOT / "wiki" / "glm_frame_headroom_probe_2026-11-25.json"
    out.write_text(json.dumps(res, indent=1) + "\n", encoding="utf-8")
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
