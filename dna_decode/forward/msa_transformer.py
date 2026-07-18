"""MSA-Transformer evolution scores: the LIFTING coevolution model for the modality hybrid.

`msa_evolution.site_independent_table` is the evolution FLOOR (a profile model — does not lift ESM2). This
module is the validated LIFTING upgrade: it runs MSA-Transformer (`esm_msa1b_t12_100M_UR50S`) on an MSA and
emits a per-variant score table oriented `higher = preserved`, ready for `variant_effect.rank_average_hybrid`
via `msa_evolution.evolution_table_from_scores`. Measured (`wiki/msa_transformer_lift_2026-07-17.md`): our
own MSA-Transformer scores reproduce ProteinGym's `MSA_Transformer` column (median Spearman ~0.91) AND the
ESM2 (+) MSA-Transformer hybrid lifts vs ESM2-650M paired — the deployable form of the +0.013 finding.

Zero-shot convention: the fast **wt-marginal** (one forward pass; read the query row's per-position
log-probs; score = logP(alt) − logP(wt)). ProteinGym uses the slower masked-marginal, so this is a fast
approximation — validated to reproduce it in rank at ~0.91. Match-column positions only (MSA inserts have no
aligned column), exactly like the site-independent floor.

**Heavy optional deps** (`torch` + `fair-esm`, lazy-imported; the 115M model caches to `TORCH_HOME`). The
whole module mirrors the ESM2 path in `variant_effect`: learned method, precompute-once, table-fed.
"""
from __future__ import annotations

from pathlib import Path

from .msa_evolution import parse_a2m, query_pos_to_col

_BUNDLE = None   # cached (model, alphabet, batch_converter)


def load_msa_transformer():
    """Lazy-load + cache MSA-Transformer (esm_msa1b_t12_100M_UR50S). Needs torch + fair-esm."""
    global _BUNDLE
    if _BUNDLE is None:
        import esm
        model, alphabet = esm.pretrained.esm_msa1b_t12_100M_UR50S()
        model.eval()
        _BUNDLE = (model, alphabet, alphabet.get_batch_converter())
    return _BUNDLE


def subsample_msa(match_cols: list[str], n_seqs: int) -> list[str]:
    """Keep the focus (row 0) + a deterministic evenly-strided subsample of the rest to depth n_seqs."""
    if len(match_cols) <= n_seqs:
        return match_cols
    rest = len(match_cols) - 1
    stride = max(1, rest // (n_seqs - 1))
    idx = [0] + list(range(1, len(match_cols), stride))[: n_seqs - 1]
    return [match_cols[i] for i in idx]


def _focus_residues(focus_raw: str) -> dict[int, str]:
    qpos, out = 0, {}
    for c in focus_raw:
        if c.isupper() or c.islower():
            qpos += 1
            out[qpos] = c.upper()
    return out


def msa_transformer_table(msa_path: str | Path, *, n_seqs: int = 128, model_bundle=None) -> dict[str, float]:
    """Run MSA-Transformer on `msa_path` and return {mutation: score}, higher = preserved, for every standard
    substitution at every MATCH-column position of the query. `model_bundle` reuses a preloaded
    (model, alphabet, batch_converter) across proteins; else loads once and caches."""
    import torch
    AA = "ACDEFGHIKLMNPQRSTVWY"
    model, alphabet, bc = model_bundle or load_msa_transformer()
    # bound memory: read at most 4000 rows (closest homologs first) before the depth-n_seqs subsample --
    # some MSAs have 100k+ rows and loading them all is GB-scale (it OOM-killed the first lift run).
    focus_name, focus_raw, match_cols = parse_a2m(msa_path, max_rows=max(4000, n_seqs * 10))
    sub = subsample_msa(match_cols, n_seqs)
    _, _, toks = bc([(f"s{i}", s) for i, s in enumerate(sub)])
    with torch.no_grad():
        logits = model(toks, repr_layers=[], return_contacts=False)["logits"]
    lp = torch.log_softmax(logits[0, 0].float(), dim=-1)      # query row 0 -> (L+1, vocab); BOS at 0

    q2c = query_pos_to_col(focus_raw)
    qres = _focus_residues(focus_raw)
    idx = alphabet.get_idx
    out: dict[str, float] = {}
    for pos, col in q2c.items():
        wt = qres.get(pos)
        if wt not in AA:
            continue
        base = lp[col + 1, idx(wt)].item()
        for alt in AA:
            if alt != wt:
                out[f"{wt}{pos}{alt}"] = lp[col + 1, idx(alt)].item() - base
    return out
