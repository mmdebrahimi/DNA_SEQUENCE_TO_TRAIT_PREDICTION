"""One-call novel-protein predictor — sequence + mutation -> the VALIDATED ESM2+ProSST hybrid, end to end.

Ties together the pieces validated separately (all locally, no Kaggle): ESM2 masked-marginal
(`esm_scorer.esm2_logp_table`, proven == the cached GPU table to ~0.01 logP) + ProSST from a 3D structure
(`prosst_scorer.prosst_variant_table` on tokens from `quantize_structure` / `fetch_alphafold_pdb`, self-
quantized == ProteinGym's pre-quantized 217/217) -> `variant_effect.predict_variant_hybrid`
(+0.067 vs ESM2 alone, N=56, `wiki/prosst_lift_2026-07-18.md`).

Two input shapes per modality (compute-or-supply), so you pay the heavy forward passes ONCE per protein:
  - ESM2: pass a precomputed `esm_table` ({pos:{aa:logp}} or {mut:score}), else it is computed from `seq`
    (WARNING: ESM2-650M masked-marginal is L forward passes — minutes on CPU for a few-hundred-aa protein).
  - ProSST: pass precomputed `prosst_table` OR `structure_tokens`, else `pdb_path`, else `uniprot`
    (fetched from AlphaFold + quantized locally).
"""
from __future__ import annotations

from pathlib import Path

from .variant_effect import ForwardPrediction, predict_variant_hybrid


def _all_single_mutants(seq: str) -> list[str]:
    aa = "ACDEFGHIKLMNPQRSTVWY"
    return [f"{seq[i]}{i + 1}{a}" for i in range(len(seq)) for a in aa if a != seq[i].upper()]


def predict_hybrid_from_sequence(protein_seq: str, mutation: str, *, protein: str = "protein",
                                 esm_table: dict | None = None, esm_model: str = "facebook/esm2_t33_650M_UR50D",
                                 prosst_table: dict | None = None, structure_tokens: list[int] | None = None,
                                 pdb_path: str | Path | None = None, uniprot: str | None = None,
                                 vocab: int = 2048, struct_cache: str | Path = "D:/dna_decode_cache/alphafold",
                                 phenotype_axis: str = "molecular fitness (DMS-measured)") -> ForwardPrediction:
    """Score ONE mutation with the validated ESM2+ProSST hybrid, computing whichever component tables aren't
    supplied. Returns the ForwardPrediction (rank + percentile among the protein's variants). Raises with a
    clear message if ProSST's structure inputs / stack are unavailable."""
    muts = _all_single_mutants(protein_seq)

    # --- ESM2 component ---
    if esm_table is None:
        from .esm_scorer import esm2_logp_table
        esm_table = esm2_logp_table(protein_seq, model_name=esm_model)   # {pos:{aa:logp}}; heavy, computed once

    # --- ProSST component ---
    if prosst_table is None:
        from .prosst_scorer import prosst_variant_table, quantize_structure, fetch_alphafold_pdb
        if structure_tokens is None:
            if pdb_path is None:
                if uniprot is None:
                    raise ValueError("need one of prosst_table / structure_tokens / pdb_path / uniprot for ProSST")
                pdb_path = fetch_alphafold_pdb(uniprot, Path(struct_cache))
            structure_tokens = quantize_structure(pdb_path, vocab)
        prosst_table = prosst_variant_table(protein_seq, muts, structure_tokens=structure_tokens, vocab=vocab)

    return predict_variant_hybrid(protein_seq, mutation, esm_table=esm_table, prosst_table=prosst_table,
                                  protein=protein, phenotype_axis=phenotype_axis)
