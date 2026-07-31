"""Masked-genome reconstruction harness -- the dog "world-model" probe (F1/F1' engine).

Masks k-mer tokens of a genomic sequence, runs a DNA masked-LM to reconstruct them, and scores the
result AGAINST an order-k Markov baseline over the IDENTICAL masked base positions.

F1' clean-comparison discipline (after the 2026-07-31 adversarial review found the F1 smoke's
comparison biased against NT on two axes):
- PRIMARY endpoint = per-base NEGATIVE LOG-LIKELIHOOD delta (Markov NLL - NT NLL; positive = NT
  better). NLL is the native surface for a probability model; accuracy is secondary. NT's per-base
  distribution is the MARGINAL over the vocab softmax (P(base@j) summed over all tokens with that base
  at offset j), NOT the single argmax 6-mer (which discards the distribution and biased against NT).
- The Markov baseline is either fit on a DISJOINT region (no leakage) or LEAVE-ONE-OUT scored (each
  target's own count excluded) -- never the F1 smoke's same-slice-fit-and-score transductive leakage.
- Both sides scored on the EXACT same masked base indices.

Self-supervised: the genome is its own label -- no phenotype label anywhere.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

from dna_decode.interp.markov_baseline import MarkovModel
from dna_decode.models.foundation import FoundationModel, TokenPrediction

_BASES = ("A", "C", "G", "T")
_BI = {b: i for i, b in enumerate(_BASES)}


def _masked_base_indices(preds: list[TokenPrediction]) -> list[int]:
    idxs: list[int] = []
    for p in preds:
        idxs.extend(range(p.base_start, p.base_start + len(p.true_kmer)))
    return idxs


def _nt_per_base(preds: list[TokenPrediction]) -> tuple[float, float, int]:
    """(mean per-base NLL, per-base marginal argmax accuracy, n_bases) from NT base_marginals."""
    total_nll = 0.0
    correct = 0
    n = 0
    for p in preds:
        if not p.base_marginals:
            raise ValueError("model did not provide per-base marginals; cannot score per-base NLL")
        for j, base in enumerate(p.true_kmer):
            if base not in _BI:
                continue
            marg = p.base_marginals[j]
            n += 1
            total_nll += -math.log(max(marg[_BI[base]], 1e-12))
            if max(range(4), key=lambda i: marg[i]) == _BI[base]:
                correct += 1
    return (total_nll / n if n else float("inf"), correct / n if n else 0.0, n)


@dataclass
class ReconstructionScore:
    """One clean masked-reconstruction result. PRIMARY = nll_delta (Markov NLL - NT NLL; >0 = NT
    beats the baseline). Accuracy metrics are secondary; the argmax-token metric is retained only for
    continuity with the deprecated F1 smoke."""

    region_label: str
    n_tokens_masked: int
    n_bases_scored: int
    # primary
    nt_per_base_nll: float
    markov_per_base_nll: float
    nll_delta: float                   # markov_nll - nt_nll ; POSITIVE = NT better  <-- HEADLINE
    # secondary
    nt_per_base_marginal_accuracy: float
    markov_per_base_accuracy: float
    accuracy_delta: float              # nt_marginal_acc - markov_acc
    markov_k: int
    markov_mode: str                   # "disjoint-fit" | "leave-one-out" | "same-slice(LEAKY)"
    model_name: str

    def as_dict(self) -> dict:
        return {
            "region_label": self.region_label,
            "PRIMARY_nll_delta_markov_minus_nt": round(self.nll_delta, 4),
            "nt_per_base_nll": round(self.nt_per_base_nll, 4),
            "markov_per_base_nll": round(self.markov_per_base_nll, 4),
            "nt_per_base_marginal_accuracy": round(self.nt_per_base_marginal_accuracy, 4),
            "markov_per_base_accuracy": round(self.markov_per_base_accuracy, 4),
            "accuracy_delta": round(self.accuracy_delta, 4),
            "markov_k": self.markov_k,
            "markov_mode": self.markov_mode,
            "n_tokens_masked": self.n_tokens_masked,
            "n_bases_scored": self.n_bases_scored,
            "model_name": self.model_name,
            "note": "PRIMARY = per-base NLL delta (Markov - NT), positive = NT beats the baseline; "
                    "NT uses per-base MARGINALS not argmax-token; Markov is disjoint-fit or LOO "
                    "(no same-slice leakage); both scored on identical masked bases",
        }


def score_reconstruction(
    model: FoundationModel,
    sequence: str,
    *,
    markov_k: int = 5,
    positions=None,
    region_label: str = "region",
    markov_fit_sequences=None,
    leave_one_out: bool | None = None,
    markov_alpha: float = 1.0,
    strict: bool = False,
) -> ReconstructionScore:
    """Clean masked-reconstruction score (per-base NLL primary).

    Markov leakage control:
    - `markov_fit_sequences` given AND disjoint from `sequence` -> disjoint-fit (leave_one_out
      defaults False).
    - `markov_fit_sequences` None -> Markov fit on `sequence` itself; leave_one_out defaults TRUE
      (exclude each target's own count) so the same-slice transductive leakage is removed.
    Pass `leave_one_out` explicitly to override the default.
    """
    if not getattr(model, "supports_mlm", False):
        raise ValueError(f"model {model.name!r} does not support masked reconstruction")
    seq = sequence.upper()
    preds = model.masked_token_predictions(seq, positions=positions, strict=strict)
    if not preds:
        raise ValueError("no maskable full-length tokens in sequence")
    return score_from_predictions(
        preds, seq, model_name=model.name, markov_k=markov_k, region_label=region_label,
        markov_fit_sequences=markov_fit_sequences, leave_one_out=leave_one_out,
        markov_alpha=markov_alpha,
    )


def score_from_predictions(
    preds: list[TokenPrediction],
    sequence: str,
    *,
    model_name: str,
    markov_k: int = 5,
    region_label: str = "region",
    markov_fit_sequences=None,
    leave_one_out: bool | None = None,
    markov_alpha: float = 1.0,
) -> ReconstructionScore:
    """Score already-computed NT predictions against a Markov baseline (lets a caller compute the
    expensive NT forward ONCE and sweep Markov k cheaply). Same leakage-control semantics as
    `score_reconstruction`."""
    seq = sequence.upper()
    nt_nll, nt_acc, n_bases = _nt_per_base(preds)
    masked_idx = _masked_base_indices(preds)

    disjoint = markov_fit_sequences is not None
    loo = (not disjoint) if leave_one_out is None else leave_one_out
    markov = MarkovModel.fit(markov_fit_sequences or [seq], k=markov_k)
    markov_nll = markov.nll_on_masked(seq, masked_idx, leave_one_out=loo, alpha=markov_alpha)
    markov_acc = markov.accuracy_on_masked(seq, masked_idx, leave_one_out=loo)
    mode = "disjoint-fit" if disjoint else ("leave-one-out" if loo else "same-slice(LEAKY)")

    return ReconstructionScore(
        region_label=region_label,
        n_tokens_masked=len(preds),
        n_bases_scored=n_bases,
        nt_per_base_nll=nt_nll,
        markov_per_base_nll=markov_nll,
        nll_delta=markov_nll - nt_nll,
        nt_per_base_marginal_accuracy=nt_acc,
        markov_per_base_accuracy=markov_acc,
        accuracy_delta=nt_acc - markov_acc,
        markov_k=markov_k,
        markov_mode=mode,
        model_name=model_name,
    )
