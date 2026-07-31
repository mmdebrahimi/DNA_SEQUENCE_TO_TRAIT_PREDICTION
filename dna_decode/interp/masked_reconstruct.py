"""Masked-genome reconstruction harness -- the dog "world-model" probe (F1 engine).

Masks k-mer tokens of a genomic sequence, runs a DNA masked-LM to reconstruct them, and scores the
result AGAINST a cheap order-k Markov baseline over the IDENTICAL masked base positions. The honest
headline is the DELTA (LM per-base accuracy - Markov per-base accuracy), never raw accuracy: a
repetitive genome is reconstructed well by a Markov model with no "understanding", so raw accuracy
would flatter the world model. A null or negative delta is a valid, bankable finding.

Self-supervised: the genome is its own label -- no phenotype label anywhere (this is what lets the
probe dodge the label wall that closed the AMR/embedding tracks).

Granularity note (load-bearing honesty): NT tokenizes into non-overlapping 6-mers, so masking ONE
token hides a 6-base span. The LM predicts those 6 bases bidirectionally from flanking 6-mers; the
Markov baseline is teacher-forced (handed the TRUE left context) -> a conservative, strong baseline.
The LM's bidirectional context vs the Markov's causal-only context is a real asymmetry, stated not
hidden; a bidirectional baseline is an F2 hardening.
"""
from __future__ import annotations

from dataclasses import dataclass, field

from dna_decode.interp.markov_baseline import MarkovModel
from dna_decode.models.foundation import FoundationModel, TokenPrediction

_BASES = frozenset("ACGT")


def _lm_base_counts(preds: list[TokenPrediction]) -> tuple[int, int]:
    """(correct_bases, total_bases) comparing pred_kmer vs true_kmer base-by-base."""
    correct = total = 0
    for p in preds:
        n = min(len(p.true_kmer), len(p.pred_kmer))
        for a, b in zip(p.true_kmer[:n], p.pred_kmer[:n]):
            if a not in _BASES:
                continue
            total += 1
            if a == b:
                correct += 1
    return correct, total


def _masked_base_indices(preds: list[TokenPrediction]) -> list[int]:
    """Absolute base indices the LM was asked to reconstruct (union over masked tokens)."""
    idxs: list[int] = []
    for p in preds:
        idxs.extend(range(p.base_start, p.base_start + len(p.true_kmer)))
    return idxs


@dataclass
class ReconstructionScore:
    """Result of one masked-reconstruction run. The DELTA is the headline; raw accuracies are
    reported for transparency but are explicitly NOT the claim."""

    region_label: str
    n_tokens_masked: int
    n_bases_scored: int
    lm_base_accuracy: float
    markov_base_accuracy: float
    delta: float                       # lm_base_accuracy - markov_base_accuracy  <-- HEADLINE
    lm_token_accuracy: float           # exact 6-mer match rate
    lm_mean_true_prob: float           # mean softmax prob the LM put on the true token
    markov_k: int
    model_name: str

    def as_dict(self) -> dict:
        return {
            "region_label": self.region_label,
            "headline_delta_lm_minus_markov": round(self.delta, 4),
            "lm_base_accuracy": round(self.lm_base_accuracy, 4),
            "markov_base_accuracy": round(self.markov_base_accuracy, 4),
            "markov_k": self.markov_k,
            "lm_token_accuracy": round(self.lm_token_accuracy, 4),
            "lm_mean_true_prob": round(self.lm_mean_true_prob, 4),
            "n_tokens_masked": self.n_tokens_masked,
            "n_bases_scored": self.n_bases_scored,
            "model_name": self.model_name,
            "note": "headline is delta vs order-k Markov; raw accuracy is NOT the claim; "
                    "masking granularity is 6-base tokens; Markov is teacher-forced (conservative)",
        }


def score_reconstruction(
    model: FoundationModel,
    sequence: str,
    *,
    markov_k: int = 5,
    positions=None,
    region_label: str = "region",
    markov_fit_sequences=None,
) -> ReconstructionScore:
    """Run masked reconstruction + the Markov baseline over the identical masked base set.

    `markov_fit_sequences` defaults to `[sequence]` (fit on the region). The LM and Markov are
    scored on the EXACT same absolute base indices, so the delta is a fair head-to-head.
    """
    if not getattr(model, "supports_mlm", False):
        raise ValueError(f"model {model.name!r} does not support masked reconstruction")
    seq = sequence.upper()

    preds = model.masked_token_predictions(seq, positions=positions)
    if not preds:
        raise ValueError("no maskable full-length tokens in sequence")

    lm_correct, n_bases = _lm_base_counts(preds)
    lm_base_acc = lm_correct / n_bases if n_bases else 0.0
    token_correct = sum(1 for p in preds if p.true_kmer == p.pred_kmer)
    lm_token_acc = token_correct / len(preds)
    mean_true_prob = sum(p.true_prob for p in preds) / len(preds)

    masked_idx = _masked_base_indices(preds)
    markov = MarkovModel.fit(markov_fit_sequences or [seq], k=markov_k)
    markov_base_acc = markov.accuracy_on_masked(seq, masked_idx)

    return ReconstructionScore(
        region_label=region_label,
        n_tokens_masked=len(preds),
        n_bases_scored=n_bases,
        lm_base_accuracy=lm_base_acc,
        markov_base_accuracy=markov_base_acc,
        delta=lm_base_acc - markov_base_acc,
        lm_token_accuracy=lm_token_acc,
        lm_mean_true_prob=mean_true_prob,
        markov_k=markov_k,
        model_name=model.name,
    )
