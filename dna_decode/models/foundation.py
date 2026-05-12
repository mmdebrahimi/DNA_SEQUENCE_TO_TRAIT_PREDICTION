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


class FoundationModelError(Exception):
    """Wrapper for foundation-model load / inference failures."""


@dataclass(frozen=True)
class ModelMetadata:
    """Static metadata for a foundation model (from config/datasources.yaml)."""

    name: str
    huggingface_id: str
    embedding_dim: int
    max_context: int


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

    def embed_batch(self, sequences: list[str]) -> np.ndarray:
        """Embed a batch of sequences. Returns shape (len(sequences), embedding_dim).

        Long sequences are mean-pooled across windows before stacking.
        """
        if not sequences:
            return np.empty((0, self.embedding_dim), dtype=np.float32)

        self._ensure_loaded()
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


class EvoModel(FoundationModel):
    """Evo (Together AI / Stanford) wrapper — microbial DNA LM, 7B params, 131K context.

    HuggingFace ID: togethercomputer/evo-1-131k-base. Real load is gated on
    transformers + bitsandbytes for 4-bit quantization (Phase 1 RTX 4090 default).
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
                self.metadata.huggingface_id, trust_remote_code=True
            )
            # Real loading would use load_in_4bit=True via bitsandbytes; deferred to
            # first real-data run to avoid dragging a multi-GB model download into tests.
            self._model = AutoModel.from_pretrained(
                self.metadata.huggingface_id, trust_remote_code=True
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
            from transformers import AutoModel, AutoTokenizer
        except ImportError as e:
            raise FoundationModelError("transformers not installed") from e
        try:
            self._tokenizer = AutoTokenizer.from_pretrained(
                self.metadata.huggingface_id, trust_remote_code=True
            )
            self._model = AutoModel.from_pretrained(
                self.metadata.huggingface_id, trust_remote_code=True
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
                self.metadata.huggingface_id, trust_remote_code=True
            )
            self._model = AutoModelForMaskedLM.from_pretrained(
                self.metadata.huggingface_id, trust_remote_code=True
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


class GenaLMModel(FoundationModel):
    """GENA-LM wrapper — multi-species DNA LM (AIRI)."""

    def _load_weights(self) -> None:
        try:
            from transformers import AutoModel, AutoTokenizer
        except ImportError as e:
            raise FoundationModelError("transformers not installed") from e
        try:
            self._tokenizer = AutoTokenizer.from_pretrained(self.metadata.huggingface_id)
            self._model = (
                AutoModel.from_pretrained(self.metadata.huggingface_id).to(self._device).eval()
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
