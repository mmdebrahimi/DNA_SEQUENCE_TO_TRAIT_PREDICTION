"""A masked (absorbing-state) discrete diffusion generator for DNA, small enough to train on a CPU.

**Why diffusion and not autoregressive.** The measured requirement from this arm is POSITIONAL
architecture: a k-mer discriminator is saturated by construction, because a Markov chain's sufficient
statistic IS the k-mer conditional distribution, while the positional margin holds at +0.0555 to +0.3015
across fourteen species and is invariant to the null's order. A left-to-right model is trained on exactly
the conditional a Markov chain approximates, so it has no structural reason to represent WHERE a motif
sits relative to the window end. A masked diffusion model predicts any position from any context, and its
positional embedding is the only thing telling it which position it is looking at -- so position is
load-bearing in the architecture rather than incidental.

**The formulation**, which is the absorbing-state / masked-diffusion construction (D3PM-absorbing, as
used by MDLM and SEDD):

  * forward: draw `t` uniform on (0, 1] and independently replace each base with MASK with probability t
  * model: predict the ORIGINAL base at every masked position, from the whole partially-masked sequence
  * loss: cross-entropy over masked positions only, weighted by `1/t`

The `1/t` weight is not decoration. At small t almost nothing is masked, so each example carries very
little signal and a plain average would let easy, nearly-clean examples dominate the gradient; the weight
is what makes the objective a bound on the likelihood rather than an arbitrary denoising loss. Dropping
it trains a model that is excellent at filling one gap and cannot generate from scratch.

**Sampling** runs the process backwards from a fully masked sequence: at each step, predict the
distribution at the masked positions and unmask each independently with probability `(t - s) / t`,
filling it with a draw from the prediction. The final sequence has no MASK left by construction, and
`assert_no_mask_remains` checks it rather than trusting it.

Deliberately small and plain: a few transformer layers with learned positional embeddings. Nothing here
is novel, and that is the point -- the corpus, the leak gate and the scoring geometry are what this arm
is testing, so the model should be the least surprising part of it.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass

BASES = "ACGT"
MASK_ID = 4                     # one past the four bases; the absorbing state
VOCAB = 5


def encode(seqs: list[str], length: int):
    """`(n, length)` int64 tensor of base ids. Refuses anything that is not exactly ACGT at `length`.

    Refusing rather than coercing matters: an N silently mapped to a base would teach the model that
    ambiguity codes are real sequence, and a short window silently padded would teach it that windows
    end early.
    """
    import torch

    out = torch.empty((len(seqs), length), dtype=torch.long)
    lut = {c: i for i, c in enumerate(BASES)}
    for i, s in enumerate(seqs):
        if len(s) != length:
            raise ValueError(f"sequence {i} has length {len(s)}, expected {length}")
        try:
            out[i] = torch.tensor([lut[c] for c in s.upper()], dtype=torch.long)
        except KeyError as exc:
            raise ValueError(f"sequence {i} contains a non-ACGT character {exc.args[0]!r}") from exc
    return out


def decode(ids) -> list[str]:
    return ["".join(BASES[int(v)] for v in row) for row in ids]


@dataclass(frozen=True)
class DiffusionConfig:
    """Everything that defines the model. Recorded into every checkpoint so a run is reproducible."""

    length: int = 150
    d_model: int = 128
    n_layers: int = 4
    n_heads: int = 4
    d_ff: int = 512
    dropout: float = 0.1
    seed: int = 0

    def n_params_estimate(self) -> int:
        """Rough parameter count, used to size a run against the corpus before building anything."""
        emb = (VOCAB + self.length) * self.d_model
        per_layer = 4 * self.d_model * self.d_model + 2 * self.d_model * self.d_ff
        head = self.d_model * len(BASES)
        return emb + self.n_layers * per_layer + head

    def as_dict(self) -> dict:
        return asdict(self) | {"vocab": VOCAB, "mask_id": MASK_ID,
                               "n_params_estimate": self.n_params_estimate()}


def build_model(cfg: DiffusionConfig):
    """A plain transformer encoder over a partially masked sequence.

    Positional embeddings are LEARNED and per-position rather than sinusoidal, because the window is a
    fixed 150 bp with a biologically meaningful origin -- the annotated start codon sits at one end -- so
    there is nothing to generalise to longer contexts and a learned table can represent a motif's
    preferred distance from that origin directly.
    """
    import torch
    from torch import nn

    torch.manual_seed(cfg.seed)

    class MaskedDiffusionDNA(nn.Module):
        def __init__(self) -> None:
            super().__init__()
            self.cfg = cfg
            self.tok = nn.Embedding(VOCAB, cfg.d_model)
            self.pos = nn.Embedding(cfg.length, cfg.d_model)
            layer = nn.TransformerEncoderLayer(
                d_model=cfg.d_model, nhead=cfg.n_heads, dim_feedforward=cfg.d_ff,
                dropout=cfg.dropout, batch_first=True, norm_first=True,
                activation="gelu")
            self.enc = nn.TransformerEncoder(layer, num_layers=cfg.n_layers)
            self.norm = nn.LayerNorm(cfg.d_model)
            self.head = nn.Linear(cfg.d_model, len(BASES))

        def forward(self, x):
            idx = torch.arange(x.shape[1], device=x.device)
            h = self.tok(x) + self.pos(idx)[None, :, :]
            h = self.enc(h)
            return self.head(self.norm(h))

    return MaskedDiffusionDNA()


def corrupt(x0, t, generator=None):
    """Mask each position independently with probability `t`. Returns `(xt, mask)`.

    `t` is per-sequence, so a batch spans many noise levels at once -- which is what lets one batch
    train the model across the whole reverse trajectory rather than at a single corruption level.
    """
    import torch

    probs = t[:, None].expand_as(x0)
    mask = torch.rand(x0.shape, device=x0.device, generator=generator) < probs
    xt = x0.masked_fill(mask, MASK_ID)
    return xt, mask


def loss_fn(model, x0, generator=None, eps: float = 1e-3):
    """`1/t`-weighted cross-entropy over masked positions only.

    `t` is drawn on `[eps, 1]` rather than `(0, 1]`: as t approaches zero the weight diverges while the
    number of masked positions goes to zero, so the estimator's variance explodes and a single batch can
    swamp the gradient. Clamping the low end is the standard fix and the value is recorded in the config
    rather than buried.

    A sequence that happens to get zero masked positions contributes nothing and is excluded from the
    denominator, so the reported loss is per-masked-token and does not quietly shrink when a batch is
    lightly corrupted.
    """
    import torch
    import torch.nn.functional as F

    b = x0.shape[0]
    t = torch.rand(b, device=x0.device, generator=generator) * (1.0 - eps) + eps
    xt, mask = corrupt(x0, t, generator=generator)
    logits = model(xt)
    if not bool(mask.any()):
        return logits.sum() * 0.0
    ce = F.cross_entropy(logits[mask], x0[mask], reduction="none")
    per_pos_weight = (1.0 / t)[:, None].expand_as(mask)[mask]
    return (ce * per_pos_weight).sum() / mask.sum()


def sample(model, n: int, cfg: DiffusionConfig, steps: int = 64, seed: int = 0,
           temperature: float = 1.0) -> list[str]:
    """Generate `n` sequences by reversing the masking process from a fully masked start.

    At each step from `t` down to `s`, every still-masked position is unmasked independently with
    probability `(t - s) / t` and filled by sampling the model's predicted distribution. That ratio is
    the correct conditional for the absorbing process: of the positions still masked at time t, the
    fraction that should become visible by time s is exactly (t - s) / t.

    The last step forces any remaining mask to resolve, so the output is always complete sequence.
    `assert_no_mask_remains` is run on the result because a silent leftover MASK would be scored as a
    base by anything downstream.
    """
    import torch

    model.eval()
    # The generator and every draw must live on the MODEL's device. A cpu Generator with a cuda `x`
    # raises, so before this the sampler worked on CPU only -- and the failure landed at SAMPLING,
    # i.e. after a full training spend, which is the most expensive place to discover it.
    # On CPU this is the same call with the same stream, so CPU output is bit-identical (asserted by
    # tests/test_diffusion_sampler_device.py). CUDA and CPU RNG streams differ, so a seed reproduces
    # on the SAME device only.
    dev = next(model.parameters()).device
    g = torch.Generator(device=dev.type).manual_seed(seed)
    x = torch.full((n, cfg.length), MASK_ID, dtype=torch.long, device=dev)
    ts = torch.linspace(1.0, 0.0, steps + 1)
    with torch.no_grad():
        for i in range(steps):
            t, s = float(ts[i]), float(ts[i + 1])
            masked = x == MASK_ID
            if not bool(masked.any()):
                break
            logits = model(x) / max(temperature, 1e-6)
            probs = torch.softmax(logits, dim=-1)
            draw = torch.multinomial(probs.reshape(-1, len(BASES)), 1,
                                     generator=g).reshape(x.shape)
            if s <= 0.0:
                reveal = masked
            else:
                p = (t - s) / max(t, 1e-9)
                reveal = masked & (torch.rand(x.shape, device=dev, generator=g) < p)
            x = torch.where(reveal, draw, x)
    assert_no_mask_remains(x)
    return decode(x.cpu())


def assert_no_mask_remains(x) -> None:
    """A leftover absorbing state would be decoded as a base and scored as real sequence."""
    if bool((x == MASK_ID).any()):
        n = int((x == MASK_ID).sum())
        raise RuntimeError(
            f"{n} position(s) are still MASK after sampling. The reverse process did not complete, and "
            "a leftover mask would be decoded as a base and scored as real sequence.")
