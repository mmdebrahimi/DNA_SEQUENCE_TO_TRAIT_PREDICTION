"""Rung-2 mutation-effect predictor for a protein point substitution (2026-07-13).

Point at a protein sequence + a single point mutation ("A123V") and get back:
  (1) the DETERMINISTIC amino-acid sequence change (certain — no model),
  (2) a molecular-effect RANK score from ESM2-650M masked marginals, honestly labelled.

The honesty rails (from the design review) are baked in:
  * `damage_llr = logP(wt) - logP(mut)` — FROZEN sign contract (higher = MORE damaging). Two other repo
    scripts use opposite conventions, so this module pins ONE and documents it in the output.
  * The score is a RANK score, NOT a per-mutation probability. It is presented with a `position_percentile`
    (how damaging vs. the other 18 substitutions at that position) and a fixed honest caveat.
  * `direction_hint` is a 3-bucket HINT derived from the percentile, never a resistant/susceptible or
    phenotype CALL.
  * Scope = rung-2 MOLECULAR function (stability/activity), NOT rung-4 cellular phenotype. The model's
    accuracy ceiling is the ProteinGym benchmark (~0.49 Spearman median, ~0.52 stability) — a zero-shot
    FACE-VALIDITY level, and the target E. coli proteins are likely in ESM's pretraining set (so this is
    face-validity, not prospective validation).

This module is separate from the AMR decoder surface (frozen surface untouched). The ESM masked-marginal
loop reuses the exact pattern from scripts/hiv_esm_vs_catalog.py; the pure functions need no torch.
"""
from __future__ import annotations

import json
import re
from pathlib import Path

AA = list("ACDEFGHIKLMNPQRSTVWY")
MODEL = "facebook/esm2_t33_650M_UR50D"
_MUT_RE = re.compile(r"^([A-Z])(\d+)([A-Z])$")

HONEST_CAVEAT = (
    "damage_llr is a zero-shot ESM2-650M RANK score (logP(wt)-logP(mut); higher = more damaging), NOT a "
    "per-mutation probability and NOT a resistant/susceptible or phenotype call. Accuracy ceiling is the "
    "ProteinGym benchmark (~0.49 Spearman median / ~0.52 on stability assays) — a zero-shot FACE-VALIDITY "
    "level; the protein is likely in ESM's pretraining set, so this is face-validity, not prospective "
    "validation. Scope: rung-2 molecular function (stability/activity), NOT rung-4 cellular phenotype."
)


class MutationParseError(ValueError):
    """Raised when a mutation string is malformed or inconsistent with the sequence."""


def parse_mutation(mutation: str):
    """'A123V' -> (wt='A', pos=123 (1-based), mut='V'). Raises on malformed input."""
    m = _MUT_RE.match(mutation.strip().upper())
    if not m:
        raise MutationParseError(f"malformed mutation {mutation!r}; expected <WT><pos><MUT> e.g. 'A123V'")
    wt, pos, mut = m.group(1), int(m.group(2)), m.group(3)
    if wt not in AA or mut not in AA:
        raise MutationParseError(f"non-standard residue in {mutation!r} (expected 20 canonical AAs)")
    if wt == mut:
        raise MutationParseError(f"{mutation!r} is a no-op (wt == mut)")
    return wt, pos, mut


def apply_edit(seq: str, mutation: str):
    """Deterministic: apply the point substitution -> {wt, pos, mut, wt_seq, mut_seq}. Certain, no model."""
    seq = seq.strip().upper()
    wt, pos, mut = parse_mutation(mutation)
    if pos < 1 or pos > len(seq):
        raise MutationParseError(f"position {pos} out of range for sequence length {len(seq)}")
    if seq[pos - 1] != wt:
        raise MutationParseError(f"mutation {mutation!r} says WT={wt} at {pos} but sequence has {seq[pos-1]}")
    mut_seq = seq[: pos - 1] + mut + seq[pos:]
    return {"wt": wt, "pos": pos, "mut": mut, "wt_seq": seq, "mut_seq": mut_seq}


def damage_llr(logp: dict, pos: int, wt: str, mut: str) -> float:
    """FROZEN contract: logP(wt) - logP(mut) at `pos` (1-based). Higher = more damaging."""
    col = logp[pos]
    return float(col[wt] - col[mut])


def position_percentile(logp: dict, pos: int, wt: str, mut: str) -> float:
    """Fraction of the 19 non-WT substitutions at `pos` that are LESS-or-equally damaging than `mut`.
    1.0 => this mutation is the most damaging option at the position; 0.0 => the least."""
    col = logp[pos]
    this = col[wt] - col[mut]
    others = [col[wt] - col[a] for a in AA if a != wt]
    if not others:
        return 0.0
    return round(sum(1 for d in others if d <= this) / len(others), 4)


def direction_hint(percentile: float) -> str:
    """A 3-bucket HINT from the position-percentile — NOT a phenotype call."""
    if percentile >= 0.8:
        return "likely-deleterious"
    if percentile <= 0.2:
        return "likely-tolerated"
    return "uncertain"


def predict(seq: str, mutation: str, logp: dict) -> dict:
    """Assemble the honest output record. `logp` = {pos(1-based): {aa: log-prob}} from ESM masked marginals."""
    edit = apply_edit(seq, mutation)
    pos, wt, mut = edit["pos"], edit["wt"], edit["mut"]
    llr = round(damage_llr(logp, pos, wt, mut), 4)
    pct = position_percentile(logp, pos, wt, mut)
    return {
        "artifact": "protein_mutation_effect", "schema": "protein-mutation-effect-v1",
        "mutation": f"{wt}{pos}{mut}",
        "sequence_change": {"wt_residue": wt, "position": pos, "mut_residue": mut,
                            "length": len(seq), "certain": True},
        "damage_llr": llr, "damage_llr_definition": "logP(wt) - logP(mut); higher = more damaging",
        "position_percentile": pct,
        "direction_hint": direction_hint(pct),
        "honest_caveat": HONEST_CAVEAT,
        "provenance": {"model": MODEL, "method": "zero-shot masked-marginal log-likelihood"},
    }


# ---- ESM masked-marginal computation (lazy torch; reuses the hiv_esm_vs_catalog pattern) ----

def masked_marginals(seq: str, cache_path: str | Path | None = None, batch: int = 8,
                     progress: bool = False) -> dict:
    """{pos(1-based): {aa: log-prob}} for every position of `seq` via ESM2-650M masked marginals (CPU).
    Cached to `cache_path` (JSON) keyed by the sequence; regenerates on sequence mismatch."""
    seq = seq.strip().upper()
    if cache_path is not None:
        cp = Path(cache_path)
        if cp.exists():
            d = json.loads(cp.read_text(encoding="utf-8"))
            if d.get("sequence") == seq:
                return {int(k): v for k, v in d["logp"].items()}
    import torch
    from transformers import AutoModelForMaskedLM, AutoTokenizer
    tok = AutoTokenizer.from_pretrained(MODEL)
    model = AutoModelForMaskedLM.from_pretrained(MODEL).eval()
    tok_ids = {a: tok.convert_tokens_to_ids(a) for a in AA}
    ids = tok(seq, return_tensors="pt")["input_ids"]     # CLS at index 0 -> residue i is token i
    positions = list(range(1, len(seq) + 1))
    logp: dict[int, dict] = {}
    with torch.no_grad():
        for s in range(0, len(positions), batch):
            chunk = positions[s:s + batch]
            b = ids.repeat(len(chunk), 1)
            for r, p in enumerate(chunk):
                b[r, p] = tok.mask_token_id
            out = model(b).logits.float().log_softmax(-1)
            for r, p in enumerate(chunk):
                logp[p] = {a: float(out[r, p, i]) for a, i in tok_ids.items()}
            if progress:
                print(f"  masked {s + len(chunk)}/{len(positions)}", flush=True)
    if cache_path is not None:
        Path(cache_path).parent.mkdir(parents=True, exist_ok=True)
        Path(cache_path).write_text(json.dumps({"sequence": seq, "model": MODEL, "logp": logp}),
                                    encoding="utf-8")
    return logp
