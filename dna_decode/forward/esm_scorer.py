"""ESM2 zero-shot masked-marginal scoring for the forward variant-effect predictor — the learned upgrade
over the deterministic BLOSUM62 baseline (`variant_effect.py`).

For ONE protein, `esm2_logp_table(seq)` runs the model ONCE per residue position (mask that residue, read
the log-softmax over the 20 amino acids) — a {pos: {aa: log-prob}} table from which ANY point mutation's
zero-shot score is instant: score(wt->alt @ pos) = logP(alt|context) - logP(wt|context). Higher = the model
prefers the mutant = more likely benign = correlates POSITIVELY with fitness (same sign as BLOSUM: higher =
preserved). This is the standard ESM zero-shot variant-effect signal (published ESM2-650M median Spearman
~0.49 on ProteinGym).

Lazy + CPU-safe (weights cached under HF_HOME); the 650M model on a single ~300-aa protein is ~L masked
forward passes, feasible on CPU. Mirrors scripts/esm_zeroshot_dms.py's masked-marginal method.
"""
from __future__ import annotations

_AA = "ACDEFGHIKLMNPQRSTVWY"
_CACHE: dict[str, tuple] = {}


def _load(model_name: str):
    if model_name not in _CACHE:
        import torch  # noqa: F401  (imported so a missing-torch env fails here, clearly)
        from transformers import AutoModelForMaskedLM, AutoTokenizer
        tok = AutoTokenizer.from_pretrained(model_name)
        model = AutoModelForMaskedLM.from_pretrained(model_name).eval()
        _CACHE[model_name] = (tok, model)
    return _CACHE[model_name]


def esm2_logp_table(seq: str, model_name: str = "facebook/esm2_t33_650M_UR50D",
                    positions=None, batch: int = 8) -> dict[int, dict[str, float]]:
    """{pos(1-based): {aa: log-prob}} via ESM2 masked-marginals over `positions` (default all residues).

    A batch of masked copies (one masked residue each) is forwarded together; at each masked token the
    log-softmax over the 20 standard AAs is recorded. Deterministic (eval mode, no sampling).
    """
    import torch
    tok, model = _load(model_name)
    L = len(seq)
    if positions is None:
        positions = list(range(1, L + 1))
    enc = tok(seq, return_tensors="pt")
    ids = enc["input_ids"][0]               # [L+2]: <cls> residues... <eos>; residue i (1-based) at index i
    mask_id = tok.mask_token_id
    aa_ids = tok.convert_tokens_to_ids(list(_AA))
    table: dict[int, dict[str, float]] = {}
    with torch.no_grad():
        for start in range(0, len(positions), batch):
            chunk = positions[start:start + batch]
            stack = ids.repeat(len(chunk), 1).clone()   # [b, L+2]
            for r, p in enumerate(chunk):
                stack[r, p] = mask_id                    # mask residue p (token index == 1-based pos)
            logits = model(stack).logits                 # [b, L+2, V]
            for r, p in enumerate(chunk):
                lp = torch.log_softmax(logits[r, p].float(), dim=-1)
                table[p] = {aa: float(lp[i]) for aa, i in zip(_AA, aa_ids)}
    return table


def esm2_delta(table: dict[int, dict[str, float]], wt: str, pos: int, alt: str) -> float:
    """Zero-shot score for wt->alt at pos: logP(alt) - logP(wt). Higher = benign (preserved). Nonsense/X
    is not in the AA table -> damaging floor. Raises KeyError if `pos` was not scored into the table."""
    row = table[pos]
    if alt in ("*", "X") or alt not in row:
        return -20.0
    if wt not in row:
        return row.get(alt, -20.0)
    return row[alt] - row[wt]
