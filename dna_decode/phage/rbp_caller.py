"""RBP-level phage receptor caller for the mixed (RBP-variable) clades.

The v0 genome-homology caller ABSTAINS / mis-calls on the RBP-variable clades (T-even Tequatrovirus,
Drexlerviridae) because whole genomes are similar across DIFFERENT receptors there — receptor is
determined by the RECEPTOR-BINDING PROTEIN (tail fiber tip), not the genome backbone. This caller
transfers the receptor from the nearest RBP by protein k-mer similarity (the GenoPHI-validated,
BLAST-free approach), so two genomically-similar T-even phages with different tail fibers get different
receptors.

Pure functions here (k-mers, similarity, nearest-transfer, leave-one-out) are offline-testable; the
.gbk RBP extraction + the LBNL cohort run live in scripts/rbp_receptor_validate.py.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

# committed RBP reference (built by scripts/build_rbp_reference.py from the MIT-licensed LBNL Phage
# Datasheets); makes the --rbp-fasta caller wheel-only (no repo clone at run time).
DEFAULT_RBP_REFERENCE = Path(__file__).resolve().parents[2] / "data" / "phage_ref" / "rbp_reference.faa"


def protein_kmers(seq: str, k: int = 4) -> frozenset[str]:
    """The set of length-k substrings of a protein sequence (upper-cased, non-AA stripped)."""
    s = "".join(c for c in seq.upper() if c.isalpha())
    if len(s) < k:
        return frozenset()
    return frozenset(s[i:i + k] for i in range(len(s) - k + 1))


def kmer_similarity(a: frozenset[str], b: frozenset[str]) -> float:
    """Jaccard similarity of two k-mer sets (0 if either empty)."""
    if not a or not b:
        return 0.0
    inter = len(a & b)
    if inter == 0:
        return 0.0
    return inter / len(a | b)


@dataclass(frozen=True)
class RbpCall:
    status: str                 # "CALLED" | "INDETERMINATE"
    predicted_receptor: str | None
    nearest_phage: str | None
    similarity: float | None
    reason: str = ""


def nearest_rbp_receptor(query_kmers: frozenset[str], ref_kmers: dict[str, frozenset[str]],
                         receptors: dict[str, str], *, min_similarity: float = 0.05,
                         exclude: str | None = None) -> RbpCall:
    """Transfer the receptor of the RBP most k-mer-similar to the query. INDETERMINATE (abstain, never
    a fabricated call) when no reference RBP clears `min_similarity`."""
    best_label, best_sim = None, 0.0
    for label, ks in ref_kmers.items():
        if label == exclude:
            continue
        sim = kmer_similarity(query_kmers, ks)
        if sim > best_sim:
            best_sim, best_label = sim, label
    if best_label is None or best_sim < min_similarity:
        return RbpCall("INDETERMINATE", None, best_label, best_sim,
                       reason=f"no reference RBP cleared min_similarity={min_similarity}")
    return RbpCall("CALLED", receptors.get(best_label), best_label, best_sim)


@dataclass
class RbpLOOResult:
    n_total: int = 0
    n_called: int = 0
    n_correct: int = 0
    per_receptor: dict[str, list[int]] = field(default_factory=dict)  # receptor -> [correct, called]
    predictions: list[dict] = field(default_factory=list)

    @property
    def accuracy(self) -> float | None:
        return (self.n_correct / self.n_called) if self.n_called else None


def load_rbp_reference(path: str | Path | None = None, k: int = 4):
    """Load the committed RBP reference FASTA (`>{phage}|{receptor}`). Returns (kmers{phage->set},
    receptors{phage->class}). Comment lines start with ';'."""
    fp = Path(path) if path else DEFAULT_RBP_REFERENCE
    kmers: dict[str, frozenset[str]] = {}
    receptors: dict[str, str] = {}
    label = receptor = None
    seq_parts: list[str] = []

    def _flush():
        if label and seq_parts:
            kmers[label] = protein_kmers("".join(seq_parts), k=k)
            receptors[label] = receptor
    with open(fp, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith(";"):
                continue
            if line.startswith(">"):
                _flush()
                head = line[1:]
                label, _, receptor = head.partition("|")
                label, receptor = label.strip(), receptor.strip()
                seq_parts = []
            else:
                seq_parts.append(line.strip())
    _flush()
    return kmers, receptors


def call_rbp_from_protein(protein: str, *, ref_path: str | Path | None = None, k: int = 4,
                          min_similarity: float = 0.05) -> RbpCall:
    """Predict a phage's receptor from its tail-fiber RBP protein sequence, by k-mer transfer against the
    committed RBP reference. Wheel-only. INDETERMINATE when no reference RBP clears `min_similarity`."""
    ref_kmers, receptors = load_rbp_reference(ref_path, k=k)
    if not ref_kmers:
        return RbpCall("INDETERMINATE", None, None, None, reason="empty RBP reference")
    return nearest_rbp_receptor(protein_kmers(protein, k=k), ref_kmers, receptors,
                                min_similarity=min_similarity)


def leave_one_out_rbp(rbp_kmers: dict[str, frozenset[str]], receptors: dict[str, str],
                      *, k: int = 4, min_similarity: float = 0.05) -> RbpLOOResult:
    """Leave-one-out RBP-homology receptor transfer over the labelled RBP set."""
    res = RbpLOOResult()
    for label, true_receptor in receptors.items():
        res.n_total += 1
        call = nearest_rbp_receptor(rbp_kmers[label], rbp_kmers, receptors,
                                    min_similarity=min_similarity, exclude=label)
        correct = call.status == "CALLED" and call.predicted_receptor == true_receptor
        if call.status == "CALLED":
            res.n_called += 1
            bucket = res.per_receptor.setdefault(true_receptor, [0, 0])
            bucket[1] += 1
            if correct:
                res.n_correct += 1
                bucket[0] += 1
        res.predictions.append({
            "phage": label, "true": true_receptor, "predicted": call.predicted_receptor,
            "nearest": call.nearest_phage, "similarity": round(call.similarity or 0.0, 3),
            "status": call.status, "correct": bool(correct)})
    return res
