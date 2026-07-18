"""Regime-B forward variant-effect predictor: (protein sequence, point mutation) -> predicted effect on
molecular fitness, with an honest confidence + regime tag + abstention.

The load-bearing validated quantity is the RANK CORRELATION between the deterministic substitution score
and the wet-lab DMS-measured effect (see scripts/tem1_forward_cell.py). BLOSUM62 is the deterministic,
no-GPU, no-network baseline the project found to be strong (ProteinGym context); ESM2 zero-shot
(scripts/esm_zeroshot_dms.py, median Spearman ~0.49) is a drop-in `method` for a later upgrade.

HONEST RAILS (hard-won):
  - This predicts MOLECULAR fitness (enzyme activity / stability / fold), the fitness-aligned regime where
    substitution-severity + learned models WORK. It is NOT an organism-level polygenic predictor (Regime C,
    a closed negative) — callers outside the molecular regime get `regime="C_organismal"` + `abstain=True`.
  - For the ANTAGONISTICALLY-SELECTED resistance direction (does an edit confer clinical resistance),
    raw likelihood/exchangeability scorers FAIL (the resistance-conservativeness finding) — use the
    Regime-A determinant catalogue there, NOT this module.
  - The per-variant tier is a coarse read of a continuous score; the VALIDATED claim is the rank
    correlation on a given protein, reported alongside every batch.
"""
from __future__ import annotations

from dataclasses import dataclass, field

# 20 standard amino acids (one-letter). '*' = stop/nonsense.
_AA = set("ACDEFGHIKLMNPQRSTVWY")

_BLOSUM62 = None


def _blosum62():
    """Lazy-load the authoritative BLOSUM62 matrix from Biopython (same source as dms_variant_effect_benchmark)."""
    global _BLOSUM62
    if _BLOSUM62 is None:
        from Bio.Align import substitution_matrices  # lazy: keeps import light for non-scoring callers
        _BLOSUM62 = substitution_matrices.load("BLOSUM62")
    return _BLOSUM62


def parse_mutation(mut: str) -> tuple[str, int, str]:
    """'M69L' / ' A42G ' -> ('M', 69, 'L'). Raises ValueError on a malformed token.

    Accepts a trailing '*' (nonsense) or 'X' as the alt. Position is 1-based (ProteinGym / biology convention).
    """
    s = (mut or "").strip()
    if len(s) < 3:
        raise ValueError(f"malformed mutation token: {mut!r}")
    wt, alt = s[0].upper(), s[-1].upper()
    pos_str = s[1:-1]
    if not pos_str.isdigit():
        raise ValueError(f"malformed mutation token (no integer position): {mut!r}")
    pos = int(pos_str)
    if wt not in _AA:
        raise ValueError(f"WT residue {wt!r} is not a standard amino acid in {mut!r}")
    if alt not in _AA and alt not in ("*", "X"):
        raise ValueError(f"ALT residue {alt!r} is not a standard amino acid / stop in {mut!r}")
    if pos < 1:
        raise ValueError(f"position must be 1-based positive in {mut!r}")
    return wt, pos, alt


def blosum62_score(wt: str, alt: str) -> float:
    """BLOSUM62 substitution score wt->alt. HIGHER = more conservative substitution = more likely to PRESERVE
    function. A nonsense/stop ('*'/'X') returns the most-damaging floor (below any real substitution)."""
    wt, alt = wt.upper(), alt.upper()
    if alt in ("*", "X"):
        return -10.0  # nonsense: maximally damaging, below BLOSUM62's minimum real entry (-4)
    m = _blosum62()
    try:
        return float(m[wt, alt])
    except (KeyError, IndexError):
        return float(m[alt, wt])  # matrix is symmetric; guard order


# BLOSUM62 diagonal (identity) ranges ~4..11; off-diagonal ~ -4..3. A synonymous/near-neutral swap scores
# high; a radical swap scores negative. These coarse tiers read the continuous score into a label.
_PRESERVED_AT = 0      # score >= 0  -> conservative -> predicted preserved
_DAMAGING_BELOW = -2   # score <= -2 -> radical      -> predicted damaging (else uncertain)


@dataclass
class ForwardPrediction:
    mutation: str
    wt: str
    pos: int
    alt: str
    protein: str
    regime: str                       # "B_molecular" | "C_organismal"
    method: str                       # "blosum62" (ESM2 drop-in later)
    raw_score: float                  # continuous; higher = more preserved function (the VALIDATED quantity)
    predicted_effect: str             # "preserved" | "damaging" | "uncertain" | "abstain"
    confidence: str                   # "high" | "medium" | "low"
    phenotype_axis: str               # what the DMS fitness measures for this protein
    abstain: bool
    notes: list[str] = field(default_factory=list)

    def as_dict(self) -> dict:
        return {
            "mutation": self.mutation, "wt": self.wt, "pos": self.pos, "alt": self.alt,
            "protein": self.protein, "regime": self.regime, "method": self.method,
            "raw_score": self.raw_score, "predicted_effect": self.predicted_effect,
            "confidence": self.confidence, "phenotype_axis": self.phenotype_axis,
            "abstain": self.abstain, "notes": self.notes,
        }


# ESM2 zero-shot delta (logP(alt)-logP(wt)) tiers — its own scale (benign ~ >= -1; damaging << 0)
_ESM_PRESERVED = -1.0
_ESM_DAMAGING = -5.0

# Hybrid rank thresholds — normalized combined rank in [0,1], higher = more preserved
_HYB_PRESERVED = 0.66
_HYB_DAMAGING = 0.33


def _midranks(vals: list[float]) -> list[float]:
    """Mid-ranks (ties share the average rank) — the documented tie-order trap."""
    order = sorted(range(len(vals)), key=lambda i: vals[i])
    r = [0.0] * len(vals)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and vals[order[j + 1]] == vals[order[i]]:
            j += 1
        mid = (i + j) / 2.0
        for k in range(i, j + 1):
            r[order[k]] = mid
        i = j + 1
    return r


def rank_average_hybrid(tables: list[dict[str, float]]) -> dict[str, float]:
    """Combine >=2 per-protein variant-effect score tables into ONE label-free rank score per variant.

    Each input table maps mutation -> a raw score already oriented so **higher = more preserved** (the
    BLOSUM/ESM2 sign). The combine is a **rank-average** over the intersection of variants: within each table
    the variants are mid-ranked, ranks are averaged across tables, and the result is normalized to [0,1]
    (higher = more preserved). This is the deployable form of the modality-hybrid finding
    (`wiki/forward_modality_hybrid_2026-07-17.md`): a naive rank-average of orthogonal modalities
    (sequence ESM2 (+) evolution GEMME (+) structure ProSST) beats ESM2-650M on 84-90% of ProteinGym proteins.

    It RANKS — it does not dose. No label, no calibrator, no fitting: the same deployability class as the
    inverse cell. Raises on <2 tables or an empty variant intersection.
    """
    if len(tables) < 2:
        raise ValueError("rank_average_hybrid needs >=2 score tables (one per modality)")
    shared = set(tables[0])
    for t in tables[1:]:
        shared &= set(t)
    if not shared:
        raise ValueError("no variant is present in all score tables — nothing to combine")
    keys = sorted(shared)
    n = len(keys)
    summed = [0.0] * n
    for t in tables:
        ranks = _midranks([t[k] for k in keys])
        for i, rk in enumerate(ranks):
            summed[i] += rk
    denom = (n - 1) * len(tables) if n > 1 else 1.0
    return {k: summed[i] / denom for i, k in enumerate(keys)}


def predict_effect(protein_seq: str, mutation: str, *, protein: str = "protein",
                   phenotype_axis: str = "molecular fitness (DMS-measured)",
                   method: str = "blosum62", regime: str = "B_molecular",
                   esm_table: dict | None = None, am_table: dict | None = None,
                   esm_if_table: dict | None = None,
                   hybrid_tables: list[dict] | None = None) -> ForwardPrediction:
    """Forward edit -> predicted phenotype effect for ONE point mutation on ONE protein.

    - method="blosum62" (default): deterministic substitution severity — no model, no network.
    - method="esm2": learned ESM2 zero-shot delta from a precomputed `esm_table`
      (dna_decode/forward/esm_scorer.esm2_logp_table) — pass the table so the model runs ONCE per protein.
    - Verifies the WT residue matches `protein_seq[pos-1]` (1-based) — a coordinate/frame error fails LOUDLY
      (the same discipline as the HIV/TB reference-integrity gates), never a silent wrong call.
    - Regime C (organismal-polygenic) -> abstain by construction (closed negative).
    - Returns a continuous `raw_score` (the validated rank quantity) + a coarse tier + confidence.
    """
    wt, pos, alt = parse_mutation(mutation)
    notes: list[str] = []

    if regime == "C_organismal":
        return ForwardPrediction(mutation, wt, pos, alt, protein, "C_organismal", method,
                                 raw_score=float("nan"), predicted_effect="abstain", confidence="low",
                                 phenotype_axis=phenotype_axis, abstain=True,
                                 notes=["organism-level polygenic trait — closed-negative regime; abstaining"])

    if protein_seq:
        if pos > len(protein_seq):
            raise ValueError(f"position {pos} beyond protein length {len(protein_seq)} for {mutation!r}")
        ref = protein_seq[pos - 1].upper()
        if ref != wt:
            raise ValueError(
                f"WT mismatch for {mutation!r}: protein has {ref!r} at position {pos}, mutation asserts {wt!r} "
                f"(coordinate/frame error — refusing to score)")
    else:
        notes.append("no protein sequence supplied — WT residue not verified against a reference")

    if method == "blosum62":
        score = blosum62_score(wt, alt)
        if alt in ("*", "X"):
            effect, conf = "damaging", "high"
            notes.append("nonsense/stop substitution — truncation, maximally damaging")
        elif score >= _PRESERVED_AT:
            effect, conf = "preserved", "medium"
        elif score <= _DAMAGING_BELOW:
            effect, conf = "damaging", "medium"
        else:
            effect, conf = "uncertain", "low"
    elif method == "esm2":
        from .esm_scorer import esm2_delta
        if esm_table is None:
            raise ValueError("method='esm2' requires esm_table= (build once via esm_scorer.esm2_logp_table)")
        score = esm2_delta(esm_table, wt, pos, alt)
        if alt in ("*", "X"):
            effect, conf = "damaging", "high"
            notes.append("nonsense/stop substitution — truncation, maximally damaging")
        elif score >= _ESM_PRESERVED:
            effect, conf = "preserved", "medium"
        elif score <= _ESM_DAMAGING:
            effect, conf = "damaging", "medium"
        else:
            effect, conf = "uncertain", "low"
    elif method == "alphamissense":
        from .am_scorer import am_tier
        if am_table is None or mutation not in am_table:
            raise ValueError(f"method='alphamissense' needs am_table with {mutation!r} "
                             f"(human-proteome only; variant not AlphaMissense-covered)")
        am = am_table[mutation]                 # AM pathogenicity in [0,1]; higher = damaging
        score = 1.0 - am                        # flip so higher = benign/preserved (BLOSUM/ESM sign)
        effect = am_tier(am)
        conf = "medium" if effect != "uncertain" else "low"
    elif method == "esm_if":
        from .structure_scorer import esm_if_tier
        if esm_if_table is None or mutation not in esm_if_table:
            raise ValueError(f"method='esm_if' needs esm_if_table with {mutation!r} "
                             f"(structure-based; build via structure_scorer.esm_if_variant_table)")
        score = esm_if_table[mutation]          # ESM-IF conditional-LL delta; higher = structure-compatible
        effect = esm_if_tier(score)
        conf = "medium" if effect != "uncertain" else "low"
    elif method == "hybrid":
        if not hybrid_tables or len(hybrid_tables) < 2:
            raise ValueError("method='hybrid' needs hybrid_tables= (>=2 per-protein score tables, one per "
                             "modality, each oriented higher=preserved)")
        combined = rank_average_hybrid(hybrid_tables)
        if mutation not in combined:
            raise ValueError(f"method='hybrid': {mutation!r} not present in all modality tables "
                             f"(rank-hybrid scores only the shared candidate set)")
        score = combined[mutation]              # normalized combined rank in [0,1]; higher = preserved
        if score >= _HYB_PRESERVED:
            effect, conf = "preserved", "medium"
        elif score <= _HYB_DAMAGING:
            effect, conf = "damaging", "medium"
        else:
            effect, conf = "uncertain", "low"
        notes.append("hybrid = rank-average of orthogonal modalities; RANKS, does not dose "
                     "(wiki/forward_modality_hybrid_2026-07-17.md)")
    else:
        raise NotImplementedError(
            f"method {method!r} not supported; use 'blosum62' / 'esm2' / 'alphamissense' / 'esm_if' / 'hybrid'")

    return ForwardPrediction(mutation, wt, pos, alt, protein, "B_molecular", method,
                             raw_score=score, predicted_effect=effect, confidence=conf,
                             phenotype_axis=phenotype_axis, abstain=False, notes=notes)
