"""Step 7 — Foundation-model wrappers (Evo / DNABERT-2 / Nucleotide Transformer / GENA-LM).

Uniform interface across pretrained DNA foundation models. Frozen embedding
mode only in Phase 1; fine-tuning is Phase 2.

Lazy loading: model weights load on first embed call, not at construction.
This keeps tests + CLI startup cheap.

Sliding window: sequences longer than the model's context get split into
overlapping windows (stride = context/2) and embeddings are mean-pooled.
"""
from __future__ import annotations

import hashlib
import os
from abc import ABC, abstractmethod
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import yaml

DEFAULT_CONFIG_PATH = Path("config/datasources.yaml")
DEFAULT_DEVICE_ENV_VAR = "DNA_DECODE_DEVICE"
DEFAULT_MODEL_SOURCE_ENV_VAR = "DNA_DECODE_MODEL_DIR"


class FoundationModelError(Exception):
    """Wrapper for foundation-model load / inference failures."""


@dataclass(frozen=True)
class ModelMetadata:
    """Static metadata for a foundation model (from config/datasources.yaml)."""

    name: str
    huggingface_id: str
    embedding_dim: int
    max_context: int


@dataclass(frozen=True)
class TokenPrediction:
    """One masked-token reconstruction result (for the dog world-model probe).

    A model that tokenizes into k-mers reconstructs at k-mer granularity: masking one token hides
    `len(true_kmer)` bases at once. `pred_kmer` is the argmax token's decoded string; base-level
    accuracy compares `pred_kmer` to `true_kmer` position-by-position (min length).
    """

    token_index: int          # 0-based index into the model's k-mer token stream (specials excluded)
    base_start: int           # 0-based absolute start of this token's bases in the input sequence
    true_kmer: str
    pred_kmer: str
    true_prob: float          # softmax prob the model assigned to the TRUE token at the masked slot
    pred_prob: float          # softmax prob of the argmax (predicted) token
    # Per-base MARGINAL distribution at each offset j in the masked token: P(base@j) summed over all
    # vocab tokens whose j-th base == that base. The FAIR per-base surface (vs the argmax-token proxy):
    # a k-tuple of (pA, pC, pG, pT). Empty when the model doesn't provide it.
    base_marginals: tuple = ()


class FoundationModel(ABC):
    """Abstract base for DNA foundation-model wrappers.

    Subclasses provide load + tokenize + forward primitives; this class
    handles batching + sliding-window aggregation uniformly.
    """

    def __init__(self, metadata: ModelMetadata, device: str | None = None):
        self.metadata = metadata
        self._device = device or os.environ.get(DEFAULT_DEVICE_ENV_VAR, "cuda")
        self._loaded = False

    @property
    def device(self) -> str:
        return self._device

    @property
    def name(self) -> str:
        return self.metadata.name

    @property
    def embedding_dim(self) -> int:
        return self.metadata.embedding_dim

    @property
    def max_context(self) -> int:
        return self.metadata.max_context

    def _model_source(self) -> str:
        """Return actual pretrained source for weight loading."""
        return _resolve_model_source(self.metadata.name, self.metadata.huggingface_id)

    def _ensure_loaded(self) -> None:
        if not self._loaded:
            self._load_weights()
            self._loaded = True

    @abstractmethod
    def _load_weights(self) -> None:
        """Subclass-specific weight loading (HuggingFace transformers / Evo SDK / etc)."""

    @abstractmethod
    def _embed_window(self, sequence: str) -> np.ndarray:
        """Compute the embedding for a sequence that fits inside `max_context`.

        Returns a 1-D array of shape (embedding_dim,).
        """

    def embed(self, sequence: str) -> np.ndarray:
        """Embed a single sequence. Sliding-window aggregation for long sequences.

        Returns shape (n_windows_or_1, embedding_dim).
        """
        if not sequence:
            raise ValueError("Cannot embed empty sequence")

        self._ensure_loaded()
        windows = self._slide_windows(sequence)
        if len(windows) == 1:
            return self._embed_window(windows[0]).reshape(1, -1)

        # Mean-pool overlapping windows into a single representation
        per_window = np.stack([self._embed_window(w) for w in windows])
        return per_window  # caller can mean-pool or per-window as needed

    def _embed_window_batch(self, sequences: list[str]) -> np.ndarray:
        """Embed a batch of sequences that each fit inside max_context.

        Default impl: per-sequence loop over `_embed_window`. Subclasses
        with batchable forward passes (NT, DNABERT-2) override this for
        true GPU batching.

        Returns shape (len(sequences), embedding_dim).
        """
        if not sequences:
            return np.empty((0, self.embedding_dim), dtype=np.float32)
        return np.stack([self._embed_window(s) for s in sequences])

    def embed_batch(self, sequences: list[str]) -> np.ndarray:
        """Embed a batch of sequences. Returns shape (len(sequences), embedding_dim).

        Fast path (all sequences ≤ max_context): single batched forward pass
        via `_embed_window_batch`. Slow path (some sequences need windowing):
        per-sequence loop with mean-pool across windows.
        """
        if not sequences:
            return np.empty((0, self.embedding_dim), dtype=np.float32)

        self._ensure_loaded()
        # Fast path: every sequence fits in one window → batched forward pass.
        if all(len(s) <= self.max_context for s in sequences):
            return self._embed_window_batch(sequences)

        # Slow path: some sequences need windowing → per-sequence loop.
        out = np.empty((len(sequences), self.embedding_dim), dtype=np.float32)
        for i, seq in enumerate(sequences):
            per_window = self.embed(seq)
            out[i] = per_window.mean(axis=0)
        return out

    def _slide_windows(self, sequence: str) -> list[str]:
        """Split a sequence into overlapping windows (stride = max_context / 2).

        Single window if `len(sequence) <= max_context`.
        """
        if len(sequence) <= self.max_context:
            return [sequence]

        stride = max(1, self.max_context // 2)
        windows: list[str] = []
        start = 0
        while start < len(sequence):
            end = min(start + self.max_context, len(sequence))
            windows.append(sequence[start:end])
            if end == len(sequence):
                break
            start += stride
        return windows

    # --- masked-reconstruction (world-model probe) -------------------------------------------
    # Capability flag: only masked-LM models (NT, DNABERT, mock) implement the forward below.
    supports_mlm: bool = False
    kmer_size: int = 1

    def masked_token_predictions(
        self, sequence: str, positions=None, batch_size: int = 16, strict: bool = False
    ) -> list["TokenPrediction"]:
        """Mask each selected k-mer token, run the MLM head, return per-token reconstructions.

        `positions` = 0-based indices into the model's k-mer token stream (specials excluded);
        None = every full-length token. `strict=True` raises if a requested position was dropped.
        Base sequence must fit `max_context`. Only masked-LM subclasses (`supports_mlm`) implement it.
        """
        raise FoundationModelError(
            f"{self.name} does not support masked reconstruction (supports_mlm is False)"
        )


class MockFoundationModel(FoundationModel):
    """Deterministic hash-based mock for tests + smoke pipeline.

    Each window's embedding is a hash-derived seed-determined random vector.
    Reproducible across runs; no GPU, no network, no model weights.
    """

    def __init__(self, metadata: ModelMetadata | None = None, device: str = "cpu"):
        super().__init__(
            metadata
            or ModelMetadata(
                name="mock", huggingface_id="mock://mock", embedding_dim=128, max_context=512
            ),
            device=device,
        )

    def _load_weights(self) -> None:
        pass  # nothing to load

    def _embed_window(self, sequence: str) -> np.ndarray:
        # Hash → seed → deterministic random vector
        h = hashlib.sha256(sequence.encode("utf-8")).digest()
        seed = int.from_bytes(h[:8], "big")
        rng = np.random.default_rng(seed)
        return rng.standard_normal(self.embedding_dim).astype(np.float32)

    # Deterministic 6-mer masked-LM mock: reconstructs EVEN token positions correctly and
    # corrupts ODD ones (first base rotated A->C->G->T->A). Gives a KNOWN base-accuracy so the
    # reconstruction scorer can be pinned offline with no weights.
    supports_mlm: bool = True
    kmer_size: int = 6

    def masked_token_predictions(self, sequence, positions=None, batch_size: int = 16,
                                 strict: bool = False):
        s = sequence.upper()
        k = self.kmer_size
        tokens = [s[i:i + k] for i in range(0, len(s) - k + 1, k)]  # non-overlapping full k-mers
        _bi = {"A": 0, "C": 1, "G": 2, "T": 3}
        requested = list(range(len(tokens))) if positions is None else list(positions)
        idxs = [t for t in requested if 0 <= t < len(tokens) and all(c in _bi for c in tokens[t])]
        if strict and len(idxs) != len(requested):
            from_ = [t for t in requested if t not in set(idxs)]
            raise FoundationModelError(f"mock masked_token_predictions(strict): dropped {from_[:5]}")
        _rot = {"A": "C", "C": "G", "G": "T", "T": "A"}

        def _marg(true_base: str, correct: bool):
            if correct:
                v = [0.1 / 3] * 4
                v[_bi[true_base]] = 0.9
            else:  # mass on the rotated wrong base -> argmax != true, true marginal 0.1
                v = [0.0] * 4
                v[_bi[true_base]] = 0.1
                v[_bi[_rot[true_base]]] = 0.9
            return tuple(v)

        out = []
        for ti in idxs:
            true_kmer = tokens[ti]
            correct_token = ti % 2 == 0
            if correct_token:
                pred_kmer, true_p, pred_p = true_kmer, 0.9, 0.9
                marg = tuple(_marg(b, True) for b in true_kmer)
            else:  # first base wrong (rotated), rest correct
                pred_kmer, true_p, pred_p = _rot.get(true_kmer[0], "N") + true_kmer[1:], 0.1, 0.6
                marg = tuple(_marg(b, j != 0) for j, b in enumerate(true_kmer))
            out.append(
                TokenPrediction(
                    token_index=ti, base_start=ti * k, true_kmer=true_kmer, pred_kmer=pred_kmer,
                    true_prob=true_p, pred_prob=pred_p, base_marginals=marg,
                )
            )
        return out


class EvoModel(FoundationModel):
    """Evo (Together AI / Stanford) wrapper — microbial DNA LM, 7B params, 131K context.

    HuggingFace ID: togethercomputer/evo-1-131k-base. Real load is gated on
    transformers + bitsandbytes for 4-bit quantization (originally planned for
    Phase 1 on RTX 4090; project's actual hardware is GTX 860M CC=5.0 so the
    4-bit Evo path is unreachable here — bitsandbytes requires CC ≥ 7.0).
    Wrapper kept as scaffolding for future compute-upgrade scenarios.
    """

    def _load_weights(self) -> None:
        try:
            from transformers import AutoModel, AutoTokenizer
        except ImportError as e:
            raise FoundationModelError(
                "transformers not installed; run `uv sync` to install Phase 1 deps"
            ) from e
        try:
            self._tokenizer = AutoTokenizer.from_pretrained(
                self._model_source(), trust_remote_code=True
            )
            # Real loading would use load_in_4bit=True via bitsandbytes; deferred to
            # first real-data run to avoid dragging a multi-GB model download into tests.
            self._model = AutoModel.from_pretrained(
                self._model_source(), trust_remote_code=True
            ).to(self._device).eval()
        except Exception as e:
            raise FoundationModelError(f"Failed to load Evo weights: {e}") from e

    def _embed_window(self, sequence: str) -> np.ndarray:
        import torch

        inputs = self._tokenizer(sequence, return_tensors="pt", truncation=True).to(
            self._device
        )
        with torch.no_grad():
            outputs = self._model(**inputs)
        # Mean-pool over the sequence dim for the final hidden state
        hidden = outputs.last_hidden_state.squeeze(0).mean(dim=0)
        return hidden.cpu().float().numpy()


class DNABERT2Model(FoundationModel):
    """DNABERT-2 wrapper — multi-species BPE-tokenized DNA model, 117M params."""

    def _load_weights(self) -> None:
        try:
            from transformers import AutoConfig, AutoModel, AutoTokenizer
        except ImportError as e:
            raise FoundationModelError("transformers not installed") from e
        try:
            model_source = self._model_source()
            self._tokenizer = AutoTokenizer.from_pretrained(
                model_source, trust_remote_code=True
            )
            config = AutoConfig.from_pretrained(model_source, trust_remote_code=True)
            # Force DNABERT2 onto the safe PyTorch attention path locally.
            if getattr(config, "attention_probs_dropout_prob", 0.0) == 0.0:
                config.attention_probs_dropout_prob = 1e-8
            self._model = AutoModel.from_pretrained(
                model_source, trust_remote_code=True, config=config
            ).to(self._device).eval()
        except Exception as e:
            raise FoundationModelError(f"Failed to load DNABERT-2 weights: {e}") from e

    def _embed_window(self, sequence: str) -> np.ndarray:
        import torch

        inputs = self._tokenizer(sequence, return_tensors="pt", truncation=True).to(
            self._device
        )
        with torch.no_grad():
            outputs = self._model(**inputs)
        hidden = outputs[0].squeeze(0).mean(dim=0)
        return hidden.cpu().float().numpy()


class NucleotideTransformerModel(FoundationModel):
    """Nucleotide Transformer v2 multi-species wrapper."""

    def _load_weights(self) -> None:
        try:
            from transformers import AutoModelForMaskedLM, AutoTokenizer
        except ImportError as e:
            raise FoundationModelError("transformers not installed") from e
        try:
            self._tokenizer = AutoTokenizer.from_pretrained(
                self._model_source(), trust_remote_code=True
            )
            self._model = AutoModelForMaskedLM.from_pretrained(
                self._model_source(), trust_remote_code=True
            ).to(self._device).eval()
        except Exception as e:
            raise FoundationModelError(f"Failed to load NT weights: {e}") from e

    def _embed_window(self, sequence: str) -> np.ndarray:
        import torch

        inputs = self._tokenizer(sequence, return_tensors="pt", truncation=True).to(
            self._device
        )
        with torch.no_grad():
            outputs = self._model(**inputs, output_hidden_states=True)
        hidden = outputs.hidden_states[-1].squeeze(0).mean(dim=0)
        return hidden.cpu().float().numpy()

    def _embed_window_batch(self, sequences: list[str]) -> np.ndarray:
        """Batched forward pass for NT with mask-aware mean pooling.

        Pads sequences to longest-in-batch via tokenizer padding=True, then
        mean-pools final hidden states using attention_mask so padding tokens
        do not contribute. Numerically equivalent to per-sequence calls for
        single-sequence input (mask is all-ones).
        """
        import torch

        if not sequences:
            return np.empty((0, self.embedding_dim), dtype=np.float32)
        inputs = self._tokenizer(
            sequences, return_tensors="pt", padding=True, truncation=True
        ).to(self._device)
        with torch.no_grad():
            outputs = self._model(**inputs, output_hidden_states=True)
        hidden = outputs.hidden_states[-1]  # (B, T, D)
        mask = inputs["attention_mask"].unsqueeze(-1).float()  # (B, T, 1)
        pooled = (hidden * mask).sum(dim=1) / mask.sum(dim=1).clamp(min=1)  # (B, D)
        return pooled.cpu().float().numpy()

    # --- masked reconstruction (world-model probe) -------------------------------------------
    supports_mlm: bool = True
    kmer_size: int = 6  # NT-v2 is a non-overlapping 6-mer tokenizer

    def _base_marginal_maps(self):
        """Cache: for each k-mer offset j, a LongTensor mapping vocab-token-id -> base index
        (0=A,1=C,2=G,3=T), or -1 for any token whose j-th char is not ACGT / is too short. Lets the
        per-base marginal P(base@j) be an index_add scatter over the full softmax."""
        import torch

        if getattr(self, "_bmm_cache", None) is not None:
            return self._bmm_cache
        vocab = self._tokenizer.get_vocab()  # token_str -> id
        vsize = max(vocab.values()) + 1
        base_idx = {"A": 0, "C": 1, "G": 2, "T": 3}
        maps = []
        for j in range(self.kmer_size):
            m = torch.full((vsize,), -1, dtype=torch.long)
            for tok, tid in vocab.items():
                if len(tok) == self.kmer_size and tok[j] in base_idx:
                    m[tid] = base_idx[tok[j]]
            maps.append(m)
        self._bmm_cache = maps
        return maps

    def masked_token_predictions(self, sequence, positions=None, batch_size: int = 16,
                                 strict: bool = False):
        """Mask each selected 6-mer token, run the MLM head, return argmax + per-base MARGINALS.

        Reads `.logits` (the MLM head already present via AutoModelForMaskedLM). Only full-length
        ACGT 6-mer tokens are maskable (tail remainder + specials skipped). `strict=True` RAISES if
        any explicitly-requested position was dropped (out of range / truncated / non-ACGT) rather
        than silently scoring a subset -- the F2 full-chromosome-sweep guard. `base_marginals` on
        each result marginalizes the full vocab softmax to per-base P(base@j) (the fair surface).
        """
        import torch

        self._ensure_loaded()
        seq = sequence.upper()
        enc = self._tokenizer(seq, return_tensors="pt", truncation=True)
        input_ids = enc["input_ids"][0]  # includes special tokens (e.g. <cls>)
        tok_strs = self._tokenizer.convert_ids_to_tokens(input_ids.tolist())
        special = set(self._tokenizer.all_special_ids)
        mask_id = self._tokenizer.mask_token_id
        if mask_id is None:
            raise FoundationModelError("NT tokenizer has no mask token")

        kmers = []  # (abs_idx, kmer_str, base_start)
        base_cursor = 0
        for abs_idx, tid in enumerate(input_ids.tolist()):
            if tid in special:
                continue
            s = tok_strs[abs_idx]
            kmers.append((abs_idx, s, base_cursor))
            base_cursor += len(s)

        requested = list(range(len(kmers))) if positions is None else list(positions)
        targets = [t for t in requested if 0 <= t < len(kmers)
                   and len(kmers[t][1]) == self.kmer_size and all(c in "ACGT" for c in kmers[t][1])]
        if strict and len(targets) != len(requested):
            dropped = [t for t in requested if t not in set(targets)]
            raise FoundationModelError(
                f"masked_token_predictions(strict): {len(dropped)} requested position(s) dropped "
                f"(out-of-range/truncated/non-ACGT), e.g. {dropped[:5]}; n_tokens={len(kmers)}"
            )

        bmm = self._base_marginal_maps()
        preds: list[TokenPrediction] = []
        for start in range(0, len(targets), batch_size):
            chunk = targets[start:start + batch_size]
            batch = input_ids.unsqueeze(0).repeat(len(chunk), 1).clone()
            for row, tpos in enumerate(chunk):
                batch[row, kmers[tpos][0]] = mask_id
            with torch.no_grad():
                logits = self._model(input_ids=batch.to(self._device)).logits  # (B, T, V)
            probs = torch.softmax(logits.float(), dim=-1).cpu()
            for row, tpos in enumerate(chunk):
                abs_idx, true_kmer, base_start = kmers[tpos]
                dist = probs[row, abs_idx]  # (V,)
                pred_id = int(torch.argmax(dist))
                true_id = int(input_ids[abs_idx])
                pred_kmer = self._tokenizer.convert_ids_to_tokens([pred_id])[0]
                marginals = []
                for j in range(self.kmer_size):
                    acc = torch.zeros(4)
                    m = bmm[j]
                    keep = m >= 0
                    acc.index_add_(0, m[keep], dist[keep])
                    marginals.append(tuple(float(x) for x in acc))
                preds.append(
                    TokenPrediction(
                        token_index=tpos, base_start=base_start, true_kmer=true_kmer,
                        pred_kmer=pred_kmer, true_prob=float(dist[true_id]),
                        pred_prob=float(dist[pred_id]), base_marginals=tuple(marginals),
                    )
                )
        return preds


class GenaLMModel(FoundationModel):
    """GENA-LM wrapper — multi-species DNA LM (AIRI)."""

    def _load_weights(self) -> None:
        try:
            from transformers import AutoModel, AutoTokenizer
        except ImportError as e:
            raise FoundationModelError("transformers not installed") from e
        try:
            self._tokenizer = AutoTokenizer.from_pretrained(self._model_source())
            self._model = (
                AutoModel.from_pretrained(self._model_source()).to(self._device).eval()
            )
        except Exception as e:
            raise FoundationModelError(f"Failed to load GENA-LM weights: {e}") from e

    def _embed_window(self, sequence: str) -> np.ndarray:
        import torch

        inputs = self._tokenizer(sequence, return_tensors="pt", truncation=True).to(
            self._device
        )
        with torch.no_grad():
            outputs = self._model(**inputs)
        hidden = outputs[0].squeeze(0).mean(dim=0)
        return hidden.cpu().float().numpy()


_MODEL_REGISTRY: dict[str, type[FoundationModel]] = {
    "mock": MockFoundationModel,
    "evo": EvoModel,
    "dnabert2": DNABERT2Model,
    "nucleotide_transformer": NucleotideTransformerModel,
    "gena_lm": GenaLMModel,
}


def _model_source_env_var(name: str) -> str:
    """Per-model override env var for a local weights directory."""
    return f"{DEFAULT_MODEL_SOURCE_ENV_VAR}_{name.upper()}"


def _resolve_model_source(name: str, huggingface_id: str) -> str:
    """Resolve model source from env override or config value."""
    env_name = _model_source_env_var(name)
    override = os.environ.get(env_name, "").strip()
    if override:
        return override
    global_override = os.environ.get(DEFAULT_MODEL_SOURCE_ENV_VAR, "").strip()
    if global_override:
        candidate = Path(global_override) / name
        if candidate.exists():
            return str(candidate)
    return huggingface_id


def model_factory(
    name: str,
    config_path: Path | str = DEFAULT_CONFIG_PATH,
    device: str | None = None,
) -> FoundationModel:
    """Construct a foundation-model wrapper by name. Reads metadata from config."""
    if name not in _MODEL_REGISTRY:
        raise FoundationModelError(
            f"Unknown foundation model: {name!r}. Known: {sorted(_MODEL_REGISTRY.keys())}"
        )

    if name == "mock":
        return MockFoundationModel(device=device or "cpu")

    path = Path(config_path)
    if not path.exists():
        raise FoundationModelError(f"Config not found: {path}")
    with open(path) as f:
        cfg = yaml.safe_load(f)
    models_cfg = cfg.get("foundation_models", {})
    if name not in models_cfg:
        raise FoundationModelError(f"Model {name!r} not in config['foundation_models']")
    meta_dict = models_cfg[name]
    metadata = ModelMetadata(
        name=name,
        huggingface_id=meta_dict["huggingface_id"],
        embedding_dim=meta_dict["embedding_dim"],
        max_context=meta_dict["max_context"],
    )
    return _MODEL_REGISTRY[name](metadata, device=device)
