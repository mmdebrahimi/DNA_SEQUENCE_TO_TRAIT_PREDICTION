"""Forward variant-effect prediction: given an E. coli protein + a minor edit (point mutation),
predict the change in phenotype. The FORWARD direction of the decoder (edit -> phenotype), validated
per-variant against FREE deep-mutational-scan (DMS) measured fitness (ProteinGym).

Scope (per the project's G2P regime boundary, feedback_g2p_decoder_regime_boundary):
  Regime A (edit hits a curated determinant)  -> deterministic AMR catalogue (amr_rules / tb_amr), already shipped.
  Regime B (edit changes protein molecular fitness/stability) -> THIS module (DMS-validated variant effect).
  Regime C (organism-level polygenic trait)   -> closed negative; abstain.
"""
from .genome_edit import (  # noqa: F401
    GenomeEditPrediction,
    cds_point_edit,
    predict_genome_edit,
    translate_codon,
)
from .variant_effect import (  # noqa: F401
    ForwardPrediction,
    blosum62_score,
    parse_mutation,
    predict_effect,
)

# esm_scorer is imported lazily (torch/transformers) — not re-exported at package import time.
