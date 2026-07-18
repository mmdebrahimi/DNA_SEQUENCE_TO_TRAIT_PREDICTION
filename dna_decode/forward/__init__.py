"""Forward variant-effect prediction: given an E. coli protein + a minor edit (point mutation),
predict the change in phenotype. The FORWARD direction of the decoder (edit -> phenotype), validated
per-variant against FREE deep-mutational-scan (DMS) measured fitness (ProteinGym).

Scope (per the project's G2P regime boundary, feedback_g2p_decoder_regime_boundary):
  Regime A (edit hits a curated determinant)  -> deterministic AMR catalogue (amr_rules / tb_amr), already shipped.
  Regime B (edit changes protein molecular fitness/stability) -> THIS module (DMS-validated variant effect).
  Regime C (organism-level polygenic trait)   -> closed negative; abstain.
"""
from .am_scorer import (  # noqa: F401
    am_table_for_mutants,
    am_tier,
    load_am_for_uniprot,
)
from .dosage import (  # noqa: F401
    DosageResult,
    conformal_q,
    dosage_intervals,
    evaluate_dosage,
)
from .genome_edit import (  # noqa: F401
    GenomeEditPrediction,
    cds_point_edit,
    predict_genome_edit,
    translate_codon,
)
from .structure_scorer import (  # noqa: F401
    StructureMethodUnavailable,
    alphafold_pdb_url,
    esm_if_tier,
    esm_if_variant_table,
)
from .prosst_scorer import (  # noqa: F401
    prosst_tier,
    prosst_variant_table,
    quantize_structure,
)
from .router import (  # noqa: F401
    REGIME_A,
    REGIME_B,
    REGIME_C,
    REGIME_UNKNOWN,
    catalogue_call,
    classify_edit,
    predict_edit,
    variant_key,
)
from .variant_effect import (  # noqa: F401
    ForwardPrediction,
    blosum62_score,
    parse_mutation,
    predict_effect,
)

# esm_scorer is imported lazily (torch/transformers) — not re-exported at package import time.
