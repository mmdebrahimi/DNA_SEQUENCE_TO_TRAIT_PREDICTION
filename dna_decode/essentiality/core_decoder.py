"""The R1 conserved-core essentiality decoder — deterministic, offline, label-independent.

Premise (established E. coli/universal-core biology): essential genes concentrate in a small set of
housekeeping FUNCTIONS. This decoder scores each gene by how strongly its product matches those core
functions and predicts essential above a threshold. It reads FUNCTION (product text), never a label,
so it is a genuine determinant-catalogue decoder (the AMR paradigm, for essentiality).

HONESTY: this is the conserved-core PRIOR. It predicts the ~universal essential core well; it will MISS
organism-specific / conditionally-essential genes (the tail the R2 learned complement targets, Family E3).
"""
from __future__ import annotations
import re
from dataclasses import dataclass, field

# Core-essential FUNCTION catalogue (curated from the universal-essential-gene literature: translation,
# DNA replication, transcription, cell-envelope/division, and a few core cofactor/lipid pathways).
# Each entry: (weight, compiled product/name pattern). Higher weight = more reliably essential.
_CORE = [
    (3.0, r"\bribosomal (protein|subunit)\b|\b30S\b|\b50S\b"),                 # ribosome
    (3.0, r"--tRNA ligase\b|\btRNA synthetase\b|aminoacyl-tRNA"),              # tRNA synthetases
    (3.0, r"DNA polymerase III|DNA-directed DNA polymerase|replicative"),      # replicative DNA pol
    (3.0, r"DNA-directed RNA polymerase\b|\bRNA polymerase (subunit|core)"),   # RNA polymerase
    (2.5, r"DNA gyrase|topoisomerase|DNA primase|replicative DNA helicase|DnaA|replication initiat"),
    (2.5, r"cell division protein Fts|divisome|septal|\bFtsZ\b|\bFtsA\b"),     # cell division
    (2.5, r"translation initiation factor|elongation factor|release factor|ribosome-binding factor"),
    (2.0, r"outer membrane protein assembly|BAM complex|\bBamA\b|Lpt|lipopolysaccharide (transport|assembly)"),
    (2.0, r"lipid (A|II)\b|UDP-N-acetyl|peptidoglycan|\bMur[A-G]\b|cell wall"),  # envelope biogenesis
    (2.0, r"signal recognition particle|preprotein translocase|Sec(A|Y|E)\b"),  # protein export core
    (1.5, r"rRNA (methyltransferase|processing|maturation)|ribosome (assembly|biogenesis)"),
    (1.5, r"undecaprenyl|isoprenoid|\bIspA?\b|polyprenyl"),                     # essential lipid carrier
    (1.5, r"chaperone (GroEL|GroES|DnaK)|folding catalyst|peptidyl-prolyl"),    # essential chaperones
]
_CORE = [(w, re.compile(p, re.I)) for w, p in _CORE]

# strong NEGATIVE signals (products that are almost never singly-essential in rich medium)
_NON = re.compile(r"transposase|prophage|\bIS[0-9]|pseudogene|hypothetical protein|transporter for|"
                  r"catabolism|degradation|utilization|resistance|toxin-antitoxin|fimbria|pilus|flagell",
                  re.I)


@dataclass
class EssentialityCall:
    gene: str
    product: str
    core_score: float
    matched: list[str] = field(default_factory=list)
    prediction: str = "non_essential"   # essential | non_essential | uncertain
    def as_dict(self):
        return {"gene": self.gene, "product": self.product, "core_score": round(self.core_score, 2),
                "matched": self.matched, "prediction": self.prediction}


def score_gene(gene: str, product: str, threshold: float = 2.0) -> EssentialityCall:
    """Score one gene's essentiality by conserved-core FUNCTION match (deterministic, label-free)."""
    text = f"{gene or ''} {product or ''}"
    score = 0.0; matched = []
    for w, pat in _CORE:
        if pat.search(text):
            score += w; matched.append(pat.pattern[:28])
    if _NON.search(product or ""):
        score -= 2.0
    pred = "essential" if score >= threshold else ("uncertain" if score > 0 else "non_essential")
    return EssentialityCall(gene, product, score, matched, pred)


def decode_genome(genes: list[tuple[str, str]], threshold: float = 2.0) -> list[EssentialityCall]:
    """genes = [(gene_symbol, product), ...] -> per-gene essentiality calls."""
    return [score_gene(g, p, threshold) for g, p in genes]
