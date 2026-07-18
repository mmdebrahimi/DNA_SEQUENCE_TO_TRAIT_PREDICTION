"""Structure-based method for the forward variant-effect predictor — ProSST (quantized-structure LM).

ProSST (Li et al., NeurIPS 2024) is the STRONG structure model: it quantizes a 3D structure into discrete
tokens (a GVP encoder + k-means), then a sequence-structure disentangled-attention transformer predicts
per-position amino-acid distributions conditioned on BOTH sequence and structure. Its zero-shot variant score
is the log-likelihood ratio `logP(alt | seq, struct) − logP(wt | seq, struct)` at the mutated position
(higher = more preserved), the same sign convention as BLOSUM / ESM2 / (1−AM) / ESM-IF.

WHY ProSST over the existing ESM-IF seam: on our own N=95 sweep `ESM2 + ProSST` gave +0.05 paired vs
ESM2-650M (win 87%) — ~4x the evolution/MSA path and the biggest single modality lever; ProSST is the
current ProteinGym structure-tier leader (zero-shot 0.504) whereas ESM-IF (0.479) lost to sequence-only ESM2.

DEPENDENCY REALITY (this host, 2026-07-18): the `PdbQuantizer` structure encoder needs `torch_geometric`
(the same stack that keeps ESM-IF seam-only on Windows/CPU) + the `prosst` package; the ProSST transformer
needs `transformers` (+ `trust_remote_code`). So this module LAZY-imports and raises
`StructureMethodUnavailable` when the stack is absent. The SEAM (predict_effect method='prosst' + the
leaderboard column + the pluggable rank_average_hybrid table) is complete + mock-tested; the real forward
pass runs on a Linux/GPU host (Kaggle T4), and — because ProteinGym ships PRE-QUANTIZED structures — the
transformer-only path (pass `structure_tokens=`) can skip `torch_geometric` entirely. Structures come from
the AlphaFold DB by UniProt via the shared `fetch_alphafold_pdb`.
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

# Reuse the AlphaFold fetch + the shared unavailable-signal from the ESM-IF structure module.
from .structure_scorer import (  # noqa: F401
    StructureMethodUnavailable,
    alphafold_pdb_url,
    fetch_alphafold_pdb,
)

# Pin the HF model revision — trust_remote_code executes the repo's modeling code (supply-chain).
# Leave None to track main; set a commit SHA to pin. Surfaced in the plan's Risk Flags.
PROSST_REVISION = None

# ProSST log-ratio tiers (log-likelihood-ratio scale, like ESM-IF's masked-marginal delta).
_PROSST_PRESERVED = -1.0
_PROSST_DAMAGING = -3.0

_AA = "ACDEFGHIKLMNPQRSTVWY"

_BUNDLE: dict = {}   # cache (model, tokenizer) per vocab size


def prosst_model_name(vocab: int = 2048) -> str:
    """HuggingFace model id for a ProSST structure-vocabulary size (20/128/512/1024/2048/4096; 2048 = best)."""
    return f"AI4Protein/ProSST-{vocab}"


def _load_prosst(vocab: int = 2048):
    """Lazy-load + cache the ProSST masked-LM + tokenizer. Needs `transformers` (+ trust_remote_code)."""
    if vocab in _BUNDLE:
        return _BUNDLE[vocab]
    try:
        from transformers import AutoModelForMaskedLM, AutoTokenizer
    except Exception as e:  # ModuleNotFoundError: transformers
        raise StructureMethodUnavailable(
            f"ProSST unavailable — `transformers` is not installed ({type(e).__name__}: {e}). "
            f"pip install transformers on a compatible host, then re-run.") from e
    name = prosst_model_name(vocab)
    model = AutoModelForMaskedLM.from_pretrained(name, trust_remote_code=True, revision=PROSST_REVISION)
    tokenizer = AutoTokenizer.from_pretrained(name, trust_remote_code=True, revision=PROSST_REVISION)
    # The ProSST checkpoint OMITS cls.predictions.decoder.weight, expecting it tied to the input embeddings.
    # Newer `transformers` no longer auto-ties it (the config uses a legacy tie key), so it loads RANDOM ->
    # garbage logits (reproduction ~0). Force-tie it (verified: CASP3 reproduction 0.013 -> 1.0). See
    # LESSONS_LEARNED 2026-07-18.
    try:
        model.cls.predictions.decoder.weight = model.get_input_embeddings().weight
    except AttributeError:
        pass   # a future ProSST variant may name the head differently; leave as-loaded
    _BUNDLE[vocab] = (model.eval(), tokenizer)
    return _BUNDLE[vocab]


# The ProSST quantizer repo (cloned; bundles the GVP AE.pt + {vocab}.joblib checkpoints).
PROSST_REPO = os.environ.get("PROSST_REPO", "D:/prosst_repo")


class _SerialPool:
    """Drop-in for pathos Pool/ThreadPool that runs in-process — the repo's Windows multiprocessing `spawn`
    deadlocks/EOFErrors on this host, and a single PDB needs no parallelism."""
    def __init__(self, *a, **k):
        pass
    def map(self, f, *iters):
        return list(map(f, *iters))
    imap = imap_unordered = map
    def close(self):
        pass
    join = terminate = restart = clear = close


def quantize_structure(pdb_path: str | Path, vocab: int = 2048) -> list[int]:
    """PDB backbone -> per-residue ProSST structure-token list (the novel-protein deploy path).

    Runs LOCALLY (validated on this Windows/CPU host, 2026-07-18: self-quantized GRB2 == ProteinGym's
    pre-quantized tokens 217/217) via the cloned ProSST repo (`$PROSST_REPO`) + three shims: a pure-python
    `torch_scatter` package, biotite 1.x's renamed `filter_backbone`, and serial pathos + `num_workers=0`
    (Windows `spawn` is unusable). Raises StructureMethodUnavailable if the repo/`torch_geometric` are absent
    — on ProteinGym you can skip this and pass the PRE-QUANTIZED `structure_tokens=` directly instead.
    """
    repo = Path(PROSST_REPO)
    if not (repo / "prosst" / "structure" / "get_sst_seq.py").exists():
        raise StructureMethodUnavailable(
            f"ProSST quantizer repo not found at {repo} — clone https://github.com/ai4protein/ProSST to "
            f"$PROSST_REPO (bundles the GVP AE.pt + {vocab}.joblib), or pass pre-quantized structure_tokens=.")
    if str(repo) not in sys.path:
        sys.path.insert(0, str(repo))
    try:
        import pathos.multiprocessing as _mp
        import pathos.threading as _th
        _mp.Pool = _mp.ProcessPool = _SerialPool          # patch BEFORE the repo binds them at import
        _th.ThreadPool = _SerialPool
        from prosst.structure.get_sst_seq import SSTPredictor   # needs torch_geometric + the torch_scatter shim
    except Exception as e:
        raise StructureMethodUnavailable(
            f"ProSST quantizer stack unavailable ({type(e).__name__}: {e}) — needs torch_geometric + the "
            f"torch_scatter shim + biotite; or pass pre-quantized structure_tokens=.") from e
    predictor = SSTPredictor(structure_vocab_size=vocab, num_processes=0, num_threads=1)  # 0 => no spawn
    res = predictor.predict_from_pdb(str(pdb_path))
    rec = res[0] if isinstance(res, list) else res
    if not isinstance(rec, dict):
        return list(rec)
    suffix = f"{vocab}_sst_seq"                            # key may be prefixed by the cluster-model path
    key = next((k for k in rec if str(k).endswith(suffix)), None)
    if key is None:
        raise KeyError(f"no '*{suffix}' key in quantizer result (keys: {list(rec)})")
    return list(rec[key])


def prosst_variant_table(wt_seq: str, mutants, *, structure_tokens: list[int] | None = None,
                         pdb_path: str | Path | None = None, vocab: int = 2048,
                         model_bundle=None) -> dict[str, float]:
    """{DMS mutation 'wt{pos}alt' -> ProSST log-ratio (logP(alt) − logP(wt)) at the mutated position}.

    Provide EITHER `structure_tokens` (pre-quantized — transformer-only, no torch_geometric) OR a `pdb_path`
    (quantized here via `quantize_structure`). WT-marginal: one forward on the wild-type sequence + structure,
    read per-position log-probs (ProSST's reference zero-shot protocol). Higher = preserved.

    Raises StructureMethodUnavailable if the ProSST stack is absent (this host). The exact tensor wiring
    (structure-token/ss_input alignment + CLS offset) is finalized against the Step-6 Kaggle run, mirroring
    how the ESM-IF sign was validated on its real run — CI exercises the seam, not this forward.
    """
    if structure_tokens is None and pdb_path is None:
        raise ValueError("prosst_variant_table needs structure_tokens= or pdb_path=")   # guard BEFORE model load
    import torch
    model, tokenizer = model_bundle or _load_prosst(vocab)
    if structure_tokens is None:
        structure_tokens = quantize_structure(pdb_path, vocab)
    if len(structure_tokens) != len(wt_seq):
        raise ValueError(f"structure tokens ({len(structure_tokens)}) != sequence length ({len(wt_seq)}) "
                         f"— quantizer/sequence mismatch")

    enc = tokenizer([wt_seq], return_tensors="pt")
    # ss_input_ids: ProSST's canonical wiring (zero_shot/proteingym_benchmark.py) shifts each structure
    # token +3 (past <pad>/<cls>/<eos>) and wraps with <cls>=1 / <eos>=2, matching the AA tokenization.
    ss = torch.tensor([[1, *[t + 3 for t in structure_tokens], 2]], dtype=torch.long)
    with torch.no_grad():
        logits = model(input_ids=enc["input_ids"], attention_mask=enc["attention_mask"],
                       ss_input_ids=ss).logits
    # strip CLS/EOS so residue at 1-based pos maps to 0-based index pos-1 (canonical slice logits[:, 1:-1])
    lp = torch.log_softmax(logits[0, 1:-1].float(), dim=-1)

    vocab_map = tokenizer.get_vocab()
    table: dict[str, float] = {}
    for m in mutants:
        m = m.strip()
        if ":" in m or len(m) < 3 or not m[1:-1].isdigit():
            continue
        wt, pos, alt = m[0], int(m[1:-1]), m[-1]
        if pos > len(wt_seq) or wt_seq[pos - 1] != wt or wt not in _AA or alt not in _AA:
            continue
        table[m] = float(lp[pos - 1, vocab_map[alt]] - lp[pos - 1, vocab_map[wt]])
    return table


def prosst_tier(delta: float) -> str:
    """ProSST log-ratio -> forward-cell tier (higher = structure-compatible = preserved)."""
    if delta >= _PROSST_PRESERVED:
        return "preserved"
    if delta <= _PROSST_DAMAGING:
        return "damaging"
    return "uncertain"
