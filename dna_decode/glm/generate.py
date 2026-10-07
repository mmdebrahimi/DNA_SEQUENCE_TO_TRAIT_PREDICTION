"""The generator seam — one interface, so the falsifier never knows which generator it is scoring.

Any object with `.name` and `.generate(n, length) -> list[str]` is a generator here. That is what lets the
whole measurement apparatus be built and validated against pure-Python baselines (`baselines.py`) before a
single weight is downloaded, and it is why swapping in a real model changes one line rather than the harness.

`HFGenerator` wraps a HuggingFace causal DNA LM (the adopted one being
`GenerTeam/GENERator-v2-prokaryote-1.2b-base` — MIT, 4.33 GiB measured, verified against the HF API rather
than taken from a survey). Everything heavy is lazily imported so this module stays importable, and the
test suite runnable, on a host with no torch and no weights.

**The fp16 trap, encoded here rather than left in a document.** `recommended_dtype()` returns fp32 on any
GPU with compute capability < 7.0. The reason is specific and bites hard: **consumer Pascal (GTX 10xx,
CC 6.1) executes fp16 at 1/64 of its fp32 rate** — only the GP100/Tesla P100 got fast fp16 — so the
instinctive "load in half precision to save VRAM" makes a 1070 roughly unusable rather than faster. Pascal
also has no bf16 at all, and bitsandbytes needs CC >= 7.0. This repo has already paid for two sibling traps
(AMD Polaris being dead for modern PyTorch; Kaggle's `enable_gpu` provisioning a P100 whose fp16 kernels do
not exist), so the rule lives in code where it cannot be forgotten.
"""
from __future__ import annotations

from typing import Protocol, runtime_checkable

#: Compute capability at/above which fp16 is genuinely fast (Volta/Turing+ have real fp16 units).
FAST_FP16_MIN_CC = (7, 0)

#: The adopted generator. Verified against the HF API 2026-10-06: exists, license=mit, 4.33 GiB.
DEFAULT_MODEL = "GenerTeam/GENERator-v2-prokaryote-1.2b-base"


@runtime_checkable
class Generator(Protocol):
    """Anything the falsifier can score."""

    name: str

    def generate(self, n: int, length: int) -> list[str]:
        ...


def recommended_dtype(device: str = "auto") -> str:
    """PURE-ish: the dtype to actually use on this host. Returns 'float32' or 'float16'.

    fp16 is recommended ONLY on a GPU with CC >= 7.0. On Pascal it would be a 64x slowdown, and on CPU
    fp16 is emulated and also slower, so both resolve to fp32.
    """
    try:
        import torch
    except ImportError:
        return "float32"
    if device == "cpu" or not torch.cuda.is_available():
        return "float32"
    major, minor = torch.cuda.get_device_capability(0)
    return "float16" if (major, minor) >= FAST_FP16_MIN_CC else "float32"


def device_report() -> dict:
    """What this host can actually do. Read-only; safe to call anywhere, torch or not."""
    out: dict = {"torch": None, "cuda": False, "device_name": None, "compute_capability": None,
                 "vram_gib": None, "recommended_dtype": "float32", "fp16_fast": False, "notes": []}
    try:
        import torch
    except ImportError:
        out["notes"].append("torch not installed; generation unavailable, baselines still work")
        return out
    out["torch"] = torch.__version__
    out["cuda"] = bool(torch.cuda.is_available())
    if out["cuda"]:
        cc = torch.cuda.get_device_capability(0)
        out["device_name"] = torch.cuda.get_device_name(0)
        out["compute_capability"] = f"{cc[0]}.{cc[1]}"
        out["vram_gib"] = round(torch.cuda.get_device_properties(0).total_memory / 2**30, 2)
        out["fp16_fast"] = cc >= FAST_FP16_MIN_CC
        if not out["fp16_fast"]:
            out["notes"].append(
                f"CC {cc[0]}.{cc[1]} < 7.0: fp16 is NOT fast here (consumer Pascal runs it at ~1/64 of "
                "fp32; Maxwell has no fp16 units). Using fp32. bitsandbytes and bf16 are also unavailable."
            )
    else:
        out["notes"].append("no CUDA device; CPU fp32")
    out["recommended_dtype"] = recommended_dtype()
    return out


class HFGenerator:
    """A HuggingFace causal DNA LM as a `Generator`. All heavy imports are lazy."""

    def __init__(self, model_id: str = DEFAULT_MODEL, *, device: str = "auto",
                 dtype: str | None = None, cache_dir: str | None = None, seed: int = 0,
                 prompt: str = "", temperature: float = 1.0, top_k: int = 0):
        self.model_id = model_id
        self.device = device
        self.dtype = dtype or recommended_dtype(device)
        self.cache_dir = cache_dir
        self.seed = seed
        self.prompt = prompt
        self.temperature = temperature
        self.top_k = top_k
        self._tok = None
        self._model = None
        #: Set by `load()` when the GPU was too small and CPU was used instead. Reported, never silent.
        self.fallback_reason: str | None = None

    def _estimated_vram_gib(self, auto_config) -> float | None:
        """Parameter bytes implied by the model config, from its `num_parameters` or dims.

        Returns None when the config does not expose enough to estimate — in which case the guard does
        NOT fire, because refusing a load on a number we could not compute would be worse than trying.
        """
        bytes_per = 2 if self.dtype in ("float16", "bfloat16") else 4
        try:
            cfg = auto_config.from_pretrained(self.model_id, cache_dir=self.cache_dir,
                                              trust_remote_code=True)
        except Exception:  # noqa: BLE001 - an un-inspectable config must not block the load
            return None
        n = getattr(cfg, "num_parameters", None)
        if not n:
            h = getattr(cfg, "hidden_size", None)
            layers = getattr(cfg, "num_hidden_layers", None)
            vocab = getattr(cfg, "vocab_size", None)
            if not (h and layers):
                return None
            # transformer block ~12*h^2 params, plus embeddings
            n = 12 * layers * h * h + (vocab or 0) * h
        return n * bytes_per / 2**30

    @property
    def name(self) -> str:
        return f"hf:{self.model_id.split('/')[-1]}"

    def load(self):
        """Materialize tokenizer + model. Separated from __init__ so construction is free and testable."""
        if self._model is not None:
            return self
        import torch
        from transformers import AutoConfig, AutoModelForCausalLM, AutoTokenizer
        dev = self.device
        if dev == "auto":
            dev = "cuda" if torch.cuda.is_available() else "cpu"
            # GUARD BEFORE THE HEAVY LOAD -- measured, not assumed. `device_report()` already knows this
            # GPU's VRAM, so OOMing is avoidable: a 1.2B model at fp32 needs ~4.8 GB and the local 860M
            # has 4.0 GiB, which failed with "tried to allocate 44.00 MiB ... 0 bytes free" AFTER paying
            # for the full download and load. Falling back to CPU is slow; crashing is useless.
            if dev == "cuda":
                need = self._estimated_vram_gib(AutoConfig)
                have = torch.cuda.get_device_properties(0).total_memory / 2**30
                if need is not None and need > have * 0.85:   # 15% headroom for activations + KV cache
                    self.fallback_reason = (
                        f"model needs ~{need:.1f} GiB at {self.dtype} but this GPU has {have:.1f} GiB; "
                        "fell back to CPU rather than OOM after loading"
                    )
                    dev = "cpu"
        td = getattr(torch, self.dtype)
        self._tok = AutoTokenizer.from_pretrained(self.model_id, cache_dir=self.cache_dir,
                                                  trust_remote_code=True)
        self._model = AutoModelForCausalLM.from_pretrained(
            self.model_id, cache_dir=self.cache_dir, dtype=td, trust_remote_code=True,
        ).to(dev).eval()
        self._dev = dev
        return self

    def generate_conditional(self, prompts: list[str], length: int, *,
                             progress_every: int = 5) -> list[str]:
        """Continue each real genomic prompt for `length` bases. ONE output per prompt.

        **Only the CONTINUATION is returned** — the prompt is stripped — so the scored region contains no
        real sequence and there is no leakage into the discriminator. This is the meaningful task, because
        the only reference set with discriminating power is the ALIGNED promoter set (see
        `corpus.window_context`), and an unconditional sample cannot be fairly scored against it.
        """
        import sys
        import time

        import torch
        self.load()
        torch.manual_seed(self.seed)
        t0 = time.time()
        max_new = max(4, length // 6 + 4)
        out: list[str] = []
        with torch.no_grad():
            for i, pr in enumerate(prompts):
                ids = self._tok(pr, return_tensors="pt").input_ids.to(self._dev)
                n_prompt_tok = ids.shape[1]
                g = self._model.generate(
                    ids, max_new_tokens=max_new, do_sample=True,
                    temperature=self.temperature,
                    top_k=self.top_k if self.top_k else None,
                    pad_token_id=self._tok.pad_token_id or self._tok.eos_token_id or 0,
                )
                # decode ONLY the newly generated tokens -- never the prompt
                new_ids = g[0][n_prompt_tok:]
                txt = self._tok.decode(new_ids, skip_special_tokens=True)
                s = "".join(c for c in txt.upper() if c in "ACGT")
                out.append(s[:length])
                if progress_every and ((i + 1) % progress_every == 0 or i == 0):
                    el = time.time() - t0
                    rate = el / (i + 1)
                    print(f"    cond-gen {i + 1}/{len(prompts)}  {rate:.1f}s/seq  "
                          f"elapsed {el / 60:.1f}m  ETA {(len(prompts) - i - 1) * rate / 60:.1f}m",
                          flush=True)
                    sys.stdout.flush()
        return out

    def generate(self, n: int, length: int, *, progress_every: int = 5) -> list[str]:
        """Sample `n` sequences of about `length` bases.

        GENERator tokenizes as 6-mers, so the emitted length is a multiple of 6 and is TRIMMED to exactly
        `length`. Trimming rather than padding keeps the length distribution exact, which matters because
        the discriminator's features are length-sensitive.

        **Progress is PRINTED, and that was earned.** A silent version of this loop burned 52 CPU-minutes
        on a 1.2B fp32 CPU run with no way to tell sequence 5 from sequence 115, so it had to be killed
        blind. A long job that reports only its verdict forces every failure to be re-created by hand;
        `progress_every` emits a rate and an ETA so a run can be sized, or abandoned, on evidence.
        """
        import sys
        import time

        import torch
        self.load()
        torch.manual_seed(self.seed)
        t0 = time.time()
        # 6-mer tokenizer -> ~length/6 new tokens, plus slack for the prompt and any special tokens.
        max_new = max(4, length // 6 + 4)
        ids = self._tok(self.prompt, return_tensors="pt").input_ids.to(self._dev) if self.prompt \
            else torch.full((1, 1), self._tok.bos_token_id or 0, device=self._dev, dtype=torch.long)
        out: list[str] = []
        with torch.no_grad():
            for i in range(n):
                g = self._model.generate(
                    ids, max_new_tokens=max_new, do_sample=True,
                    temperature=self.temperature,
                    top_k=self.top_k if self.top_k else None,
                    pad_token_id=self._tok.pad_token_id or self._tok.eos_token_id or 0,
                )
                txt = self._tok.decode(g[0], skip_special_tokens=True)
                s = "".join(c for c in txt.upper() if c in "ACGT")
                out.append(s[:length])
                if progress_every and ((i + 1) % progress_every == 0 or i == 0):
                    el = time.time() - t0
                    rate = el / (i + 1)
                    print(f"    gen {i + 1}/{n}  {rate:.1f}s/seq  elapsed {el / 60:.1f}m  "
                          f"ETA {(n - i - 1) * rate / 60:.1f}m", flush=True)
                    sys.stdout.flush()
        return out
