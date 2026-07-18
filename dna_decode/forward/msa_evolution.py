"""MSA -> per-variant evolutionary-score table: the run-time evolution component of the modality hybrid.

The modality-hybrid finding (`wiki/forward_modality_hybrid_2026-07-17.md`) is that combining ESM2 (learned
sequence) with an ORTHOGONAL evolution signal beats ESM2-650M alone. This module is the deployable scaffold
for the evolution component on a NOVEL protein: parse an MSA, reweight sequences, and emit a per-variant
score table oriented `higher = preserved` -- ready to drop into `variant_effect.rank_average_hybrid`.

**Honest tier (measured, `scripts/forward_modality_hybrid_sweep.py` R2 scan 2026-07-17):** the
SITE-INDEPENDENT profile model this ships is the evolution FLOOR -- on its own it does NOT lift ESM2 in the
hybrid (Site_Independent+ESM2 delta -0.003, win 47%, p=0.68). The lift needs GEMME-grade coevolution
(+0.022) or MSA-Transformer (+0.013); EVmutation is a marginal +0.005. So this module's VALUE is the
reusable pipeline (MSA parse -> reweight -> per-variant score -> hybrid adapter) with a PLUGGABLE evolution
model: `evolution_table_from_scores` accepts ANY {mutation: score} table (a precomputed GEMME / MSA-T run),
and `site_independent_table` is the built-in floor. Swap the model, keep the pipeline.

a2m format: the FIRST record is the focus/query (full protein); UPPERCASE = match columns, lowercase = query
inserts (no aligned column -> unscorable), '.' = insert gaps in other seqs. Match columns are extracted by
keeping uppercase + '-'. A mutation at a query position that maps to an insert gets NO evolution score
(exactly as ProteinGym's Site_Independent leaves it NaN).
"""
from __future__ import annotations

import math
from pathlib import Path

AA = "ACDEFGHIKLMNPQRSTVWY"
AA_SET = frozenset(AA)


def parse_a2m(path: str | Path) -> tuple[str, str, list[str]]:
    """Return (focus_name, focus_raw, match_column_sequences). Each match-column sequence keeps only the
    match columns (uppercase + '-' in the alignment), so all are the same length as the focus's match set."""
    name, buf, records = None, "", []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n").rstrip("\r")
            if line.startswith(">"):
                if name is not None:
                    records.append((name, buf))
                name, buf = line[1:], ""
            else:
                buf += line.strip()
    if name is not None:
        records.append((name, buf))
    if not records:
        raise ValueError(f"empty MSA: {path}")
    focus_name, focus_raw = records[0]

    def match_only(s: str) -> str:
        return "".join(c for c in s if c.isupper() or c == "-")

    match_cols = [match_only(s) for s in (r[1] for r in records)]
    ncol = len(match_cols[0])
    if any(len(s) != ncol for s in match_cols):
        raise ValueError("MSA match-column extraction produced ragged rows (non-standard a2m)")
    return focus_name, focus_raw, match_cols


def query_pos_to_col(focus_raw: str) -> dict[int, int]:
    """1-based query residue position -> match-column index. Insert (lowercase) query residues map to NO
    column (absent from the dict); focus gaps '-' consume a column without a query residue."""
    qpos, col, out = 0, 0, {}
    for c in focus_raw:
        if c == "-":                 # match column, no query residue
            col += 1
        elif c.isupper():            # match column + query residue
            qpos += 1
            out[qpos] = col
            col += 1
        elif c.islower():            # query insert residue, no column
            qpos += 1
    return out


def sequence_weights(match_cols: list[str], theta: float = 0.2, max_seqs: int | None = None) -> list[float]:
    """Reweight each sequence by 1/(cluster size at (1-theta) identity) over MATCH columns. O(N^2) -- for a
    large MSA pass ProteinGym's precomputed .npy instead (load_weights_npy). `max_seqs` caps N for a novel
    protein's on-the-fly compute (a documented approximation, not used when a canonical .npy is supplied)."""
    cols = match_cols if max_seqs is None else match_cols[:max_seqs]
    n = len(cols)
    L = len(cols[0]) if cols else 0
    thresh = (1.0 - theta) * L
    counts = [1.0] * n
    for i in range(n):
        si = cols[i]
        for j in range(i + 1, n):
            sj = cols[j]
            ident = sum(1 for a, b in zip(si, sj) if a == b)
            if ident >= thresh:
                counts[i] += 1.0
                counts[j] += 1.0
    return [1.0 / c for c in counts]


def load_weights_npy(path: str | Path) -> list[float]:
    """Load ProteinGym's precomputed per-sequence weights (.npy) -- the canonical reweighting."""
    import numpy as np
    return np.load(path).astype(float).tolist()


def site_independent_table(msa_path: str | Path, *, weights: list[float] | None = None,
                           weights_npy: str | Path | None = None, theta: float = 0.2,
                           pseudocount: float = 1.0, max_seqs_for_weights: int | None = 20000
                           ) -> dict[str, float]:
    """Weighted site-independent log-odds table: {mutation: score}, score = log f(mut,col) - log f(wt,col),
    higher = more preserved (the BLOSUM/ESM2 sign). Only MATCH-column positions get scored (inserts have no
    column). Every non-wt standard AA at every scorable position is emitted (e.g. "M69L")."""
    focus_name, focus_raw, match_cols = parse_a2m(msa_path)
    n = len(match_cols)
    if weights is None:
        weights = (load_weights_npy(weights_npy) if weights_npy is not None
                   else sequence_weights(match_cols, theta=theta, max_seqs=max_seqs_for_weights))
    if len(weights) != n:
        raise ValueError(f"weights length {len(weights)} != MSA depth {n}")

    L = len(match_cols[0])
    # weighted amino-acid counts per column
    wcount = [{a: 0.0 for a in AA} for _ in range(L)]
    for k in range(n):
        seq, wk = match_cols[k], weights[k]
        for j, ch in enumerate(seq):
            if ch in AA_SET:
                wcount[j][ch] += wk
    # weighted frequencies with a symmetric pseudocount over the 20-aa alphabet
    freq = []
    for j in range(L):
        tot = sum(wcount[j].values()) + pseudocount * len(AA)
        freq.append({a: (wcount[j][a] + pseudocount) / tot for a in AA})

    q2c = query_pos_to_col(focus_raw)
    # the query residue at each scorable position (mirrors query_pos_to_col's positional walk)
    qpos, qres = 0, {}
    for c in focus_raw:
        if c.isupper():
            qpos += 1
            qres[qpos] = c.upper()
        elif c.islower():
            qpos += 1
            qres[qpos] = c.upper()
        # '-' consumes no query residue

    out: dict[str, float] = {}
    for pos, col in q2c.items():
        wt = qres.get(pos)
        if wt not in AA_SET:
            continue
        fwt = freq[col][wt]
        base = math.log(fwt)
        for alt in AA:
            if alt == wt:
                continue
            out[f"{wt}{pos}{alt}"] = math.log(freq[col][alt]) - base
    return out


def evolution_table_from_scores(scores: dict[str, float]) -> dict[str, float]:
    """Pass-through adapter for a PRECOMPUTED evolution model (GEMME / MSA-Transformer / EVmutation): a
    {mutation: score} table already oriented higher=preserved, ready for `rank_average_hybrid`. The pluggable
    upgrade slot -- swap the floor site-independent model for a coevolution model without touching the pipe."""
    if not scores:
        raise ValueError("empty evolution score table")
    return dict(scores)
