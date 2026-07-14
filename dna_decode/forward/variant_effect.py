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


def predict_effect(protein_seq: str, mutation: str, *, protein: str = "protein",
                   phenotype_axis: str = "molecular fitness (DMS-measured)",
                   method: str = "blosum62", regime: str = "B_molecular",
                   esm_table: dict | None = None, am_table: dict | None = None) -> ForwardPrediction:
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
    else:
        raise NotImplementedError(f"method {method!r} not supported; use 'blosum62' / 'esm2' / 'alphamissense'")

    return ForwardPrediction(mutation, wt, pos, alt, protein, "B_molecular", method,
                             raw_score=score, predicted_effect=effect, confidence=conf,
                             phenotype_axis=phenotype_axis, abstain=False, notes=notes)
