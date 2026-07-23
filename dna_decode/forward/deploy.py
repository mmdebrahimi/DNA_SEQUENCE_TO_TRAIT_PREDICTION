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

from .variant_effect import ForwardPrediction, predict_effect, predict_variant_hybrid


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


# ---------------------------------------------------------------------------
# Deployable one-call path: resolve the requested method against host capabilities, compute whatever
# tables that method needs, degrade honestly when a dependency is missing, and stamp full provenance.
# This is what `dna-decode forward --method <strong>` calls. Heavy computes go through module-level names
# (_esm2_logp_table / _prosst_table_from_structure / _gemme_table) so tests can inject cheap stand-ins.
# ---------------------------------------------------------------------------

def _esm2_logp_table(seq: str, model: str = "facebook/esm2_t33_650M_UR50D") -> dict:
    from .esm_scorer import esm2_logp_table
    return esm2_logp_table(seq, model_name=model)


def _prosst_table_from_structure(seq: str, muts: list[str], *, uniprot: str | None, pdb_path=None,
                                 structure_tokens: list[int] | None = None, vocab: int = 2048,
                                 struct_cache: str | Path = "D:/dna_decode_cache/alphafold") -> dict:
    from .prosst_scorer import prosst_variant_table, quantize_structure, fetch_alphafold_pdb
    if structure_tokens is None:
        if pdb_path is None:
            if uniprot is None:
                raise ValueError("ProSST needs a structure: pass --uniprot, --pdb, or --structure-tokens")
            pdb_path = fetch_alphafold_pdb(uniprot, Path(struct_cache))
        structure_tokens = quantize_structure(pdb_path, vocab)
    return prosst_variant_table(seq, muts, structure_tokens=structure_tokens, vocab=vocab)


def _gemme_table(seq: str, muts: list[str], msa) -> dict:
    from .gemme_scorer import run_gemme, gemme_table_from_column
    col = run_gemme(msa, seq)
    return gemme_table_from_column(col, muts)


def predict_effect_deployable(protein_seq: str, mutation: str, *, method: str = "auto",
                              protein: str = "protein", uniprot: str | None = None,
                              pdb_path=None, structure_tokens: list[int] | None = None,
                              esm_table: dict | None = None, prosst_table: dict | None = None,
                              gemme_table: dict | None = None, am_table: dict | None = None,
                              msa=None, degrade: bool = True, caps: dict | None = None,
                              phenotype_axis: str = "molecular fitness (DMS-measured)") -> dict:
    """Run the strongest requested-and-runnable forward method end to end, with honest provenance.

    method="auto"  -> pick the strongest method runnable on this host.
    A requested method whose dependency is absent -> degrade to the strongest runnable method (degrade=True)
    or raise RuntimeError (degrade=False). The returned dict records method_requested / method_used /
    degraded / degrade_reason so the caller never mistakes a fallback for the real thing.
    """
    from .capabilities import runnable_methods

    esm_table = _coerce_esm_pos_keys(esm_table)
    rm = runnable_methods(caps)
    has_seq = bool(protein_seq)
    has_structure = bool(uniprot or pdb_path or structure_tokens or prosst_table is not None)
    has_msa = bool(msa is not None or gemme_table is not None)

    # "computable NOW" = installed dep is present AND the call-time inputs this method needs were supplied
    # (a precomputed table also satisfies it). This is input-AWARE, so `auto` never picks a method whose
    # structure/MSA is missing on this call.
    def _computable(m: str) -> tuple[bool, str]:
        if m == "blosum62":
            return True, ""
        if m == "esm2":
            if esm_table is not None:
                return True, ""
            if not rm["esm2"][0]:
                return False, rm["esm2"][1]
            return (has_seq, "" if has_seq else "esm2 needs a protein sequence")
        if m == "prosst":
            if prosst_table is not None:
                return has_seq, "" if has_seq else "prosst needs a protein sequence"
            if not rm["prosst"][0]:
                return False, rm["prosst"][1]
            if not has_structure:
                return False, "prosst needs a structure (--uniprot / --pdb / --structure-tokens)"
            return (has_seq, "" if has_seq else "prosst needs a protein sequence")
        if m == "gemme":
            if gemme_table is not None:
                return has_seq, "" if has_seq else "gemme needs a protein sequence"
            if not rm["gemme"][0]:
                return False, rm["gemme"][1]
            if not has_msa:
                return False, "gemme needs an MSA (--msa)"
            return (has_seq, "" if has_seq else "gemme needs a protein sequence")
        if m == "hybrid":
            n = sum(_computable(x)[0] for x in ("esm2", "prosst", "gemme"))
            return (n >= 2, "" if n >= 2 else f"hybrid needs >=2 computable modalities here (have {n})")
        return False, f"unknown method {m!r}"

    def _strongest_computable() -> str:
        return next((m for m in ("hybrid", "prosst", "esm2", "gemme", "blosum62") if _computable(m)[0]),
                    "blosum62")

    requested = method
    if method == "auto":
        method = _strongest_computable()

    degraded = False
    degrade_reason = None
    ok, why = _computable(method)
    if not ok:
        if not degrade:
            raise RuntimeError(f"method {method!r} not runnable here: {why}")
        fallback = _strongest_computable()
        degraded = True
        degrade_reason = f"requested {method!r} not runnable ({why}); used {fallback!r}"
        method = fallback

    muts = _all_single_mutants(protein_seq) if protein_seq else [mutation]

    # compute whatever tables the resolved method needs (unless supplied)
    if method in ("esm2", "hybrid") and esm_table is None and protein_seq:
        esm_table = _esm2_logp_table(protein_seq)
    if method in ("prosst", "hybrid") and prosst_table is None and protein_seq:
        prosst_table = _prosst_table_from_structure(protein_seq, muts, uniprot=uniprot, pdb_path=pdb_path,
                                                    structure_tokens=structure_tokens)
    if method in ("gemme", "hybrid") and gemme_table is None and msa is not None and protein_seq:
        gemme_table = _gemme_table(protein_seq, muts, msa)

    if method == "hybrid":
        tables = [t for t in (esm_table_to_variant(esm_table, protein_seq), prosst_table, gemme_table)
                  if t is not None]
        pred = predict_effect(protein_seq, mutation, protein=protein, phenotype_axis=phenotype_axis,
                              method="hybrid", hybrid_tables=tables)
    else:
        pred = predict_effect(protein_seq, mutation, protein=protein, phenotype_axis=phenotype_axis,
                              method=method, esm_table=esm_table, prosst_table=prosst_table,
                              gemme_table=gemme_table, am_table=am_table)

    d = pred.as_dict()
    d.update({"method_requested": requested, "method_used": method, "degraded": degraded,
              "degrade_reason": degrade_reason})
    return d


def _coerce_esm_pos_keys(esm_table: dict | None) -> dict | None:
    """A pos-table loaded from JSON has string keys ({'1':{...}}); esm2_delta needs int keys. Coerce
    numeric-string top-level keys to int for a pos-table; leave a {mut:score} table untouched."""
    if not esm_table:
        return esm_table
    sample_v = next(iter(esm_table.values()))
    if isinstance(sample_v, dict):   # pos-table shape
        out = {}
        for k, v in esm_table.items():
            out[int(k) if isinstance(k, str) and k.lstrip("-").isdigit() else k] = v
        return out
    return esm_table


def esm_table_to_variant(esm_table: dict | None, protein_seq: str) -> dict | None:
    """Normalize an ESM2 position-table ({pos:{aa:logp}}) into a {mutation: score} table for the hybrid.
    A {mut: score} table (or None) passes through unchanged."""
    if esm_table is None:
        return None
    # position-table shape -> convert; already-variant shape -> passthrough
    sample = next(iter(esm_table), None)
    if isinstance(esm_table.get(sample), dict):
        from .variant_effect import esm_pos_table_to_variant_table
        return esm_pos_table_to_variant_table(esm_table, protein_seq)
    return esm_table
