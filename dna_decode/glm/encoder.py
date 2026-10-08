"""A small learned convolutional sequence encoder — G-A's candidate representation.

**Built against a stop gate that predicted it would fail, on an explicit user decision (2026-10-08).**
`scripts/glm_tile_headroom.py` measured that adding POSITION to composition on real genomic promoter tiles
buys **−0.0068** against a +0.02 bar, and a conv encoder is fundamentally a composition-plus-position
model, so `wiki/glm_tile_headroom_2026-10-08.json` records it as a *predicted*-negative and recommended
ratification before spending the build. The user ratified building it anyway, and that is a defensible
call for a reason the gate itself names: the gate's two honest limits are exactly the two things an
encoder could beat.

1. Positional GC in 10 bins is the **minimal** positional feature. A motif at a precise offset is
   invisible to a bulk shift in one bin.
2. The gate's ridge is **linear in its features**, so a composition×position INTERACTION is not
   representable there at all.

So the encoder is now the test of its own gate: if it also fails, the negative stops being a prediction
from a weak proxy and becomes a measurement from the model class the proxy stood in for. That is a
stronger result than the gate could produce alone, and it is why this module exists.

**What is deliberately small.** 32 channels per branch, three branches. Capacity has already HURT once on
this substrate — 6-mers (4,096 features) scored WORST of every arm (0.2022) — so the risk here is not too
little capacity but too much, and the shuffled-label null is what detects it.

**Three kernel widths, each with a stated reason**, not a hyperparameter sweep:

| width | scale it can see |
|---|---|
| 6 | a core motif box (a −10/−35 element is 6 bp) |
| 20 | motif-plus-spacer, i.e. a two-box promoter with its gap |
| 75 | half the 150 bp tile — whole-promoter architecture |

Widths produce different output lengths (`L − k + 1`), so **pooling before concatenation is structural,
not stylistic** — the branches cannot be concatenated along length.

**torch is lazy-imported inside `fit_encoder`**, mirroring `generate.py`, so every pure helper here and
its tests stay runnable on a default install where torch is absent.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass

#: Fixed alphabet order. One-hot rows are indexed by this, so changing it changes every trained model.
BASES = "ACGT"

#: The three branch widths, with their rationale in the module docstring.
DEFAULT_WIDTHS = (6, 20, 75)


class EncoderDataError(ValueError):
    """Raised on an input shape that could only be handled by silently inventing data."""


@dataclass(frozen=True)
class TrainingProtocol:
    """The FROZEN training protocol. A dataclass, not loose kwargs, so a post-hoc tweak shows in a diff.

    This is the thing Step 7's `registered_protocol_used` flag attests to. If a run needs a different
    protocol, the honest move is a new named instance recorded in the artifact — never an ad-hoc override
    that leaves the artifact claiming the registered one was used.
    """

    optimizer: str = "adam"
    lr: float = 1e-3
    weight_decay: float = 0.0
    batch_size: int = 256
    max_epochs: int = 40
    early_stop_patience: int = 5
    channels: int = 32
    val_fraction: float = 0.2

    def as_dict(self) -> dict:
        return asdict(self)


#: The registered protocol. Step 7 asserts the run used exactly this.
REGISTERED_PROTOCOL = TrainingProtocol()


# ---------------------------------------------------------------------------------------------------
# pure helpers -- no torch, so these stay testable on a default install
# ---------------------------------------------------------------------------------------------------
def conv_output_length(length: int, width: int) -> int:
    """`L - k + 1`. Returns <= 0 for a kernel wider than the sequence, which callers must refuse."""
    return length - width + 1


def encoder_feature_width(widths=DEFAULT_WIDTHS, *, channels: int = 32) -> int:
    """Width of the concatenated pooled vector: two pools (max, mean) per channel per branch."""
    return 2 * channels * len(widths)


def require_uniform_length(seqs: list[str]) -> int:
    """Return the common length, or REFUSE.

    A ragged batch has exactly two silent resolutions — pad or truncate — and both invent sequence that
    was never measured, then return a plausible number. The tile substrate is a fixed 150 bp, so a ragged
    input here means the caller mixed in the frag library, which is a real bug worth stopping on.
    """
    if not seqs:
        raise EncoderDataError("no sequences")
    L = len(seqs[0])
    if L == 0:
        raise EncoderDataError("zero-length sequence")
    bad = [i for i, s in enumerate(seqs) if len(s) != L]
    if bad:
        got = sorted({len(seqs[i]) for i in bad})
        raise EncoderDataError(
            f"ragged input: {len(bad)} of {len(seqs)} sequences differ from the first ({L} bp); "
            f"other lengths {got[:5]}. Refusing rather than padding or truncating.")
    return L


def one_hot(seqs: list[str]):
    """`(N, 4, L)` float32, the layout `Conv1d` requires (N, C, L) — verified, not assumed.

    An unknown base (N, or lowercase soft-mask) becomes an ALL-ZERO column rather than a uniform 0.25 or a
    silent 'A'. Zero says "no evidence of any base here", which is what an N means; 0.25 would assert
    equal evidence for all four, and picking 'A' would fabricate a base.
    """
    import numpy as np
    L = require_uniform_length(seqs)
    idx = {b: i for i, b in enumerate(BASES)}
    out = np.zeros((len(seqs), len(BASES), L), dtype=np.float32)
    for n, s in enumerate(seqs):
        for p, ch in enumerate(s.upper()):
            i = idx.get(ch)
            if i is not None:
                out[n, i, p] = 1.0
    return out


def validate_widths(widths, length: int) -> tuple[int, ...]:
    """Every width must fit inside the sequence, and the set must be non-empty and unique."""
    w = tuple(int(x) for x in widths)
    if not w:
        raise EncoderDataError("no kernel widths given")
    if len(set(w)) != len(w):
        raise EncoderDataError(f"duplicate kernel widths {w}")
    if any(x < 1 for x in w):
        raise EncoderDataError(f"kernel widths must be >= 1, got {w}")
    too_wide = [x for x in w if conv_output_length(length, x) < 1]
    if too_wide:
        raise EncoderDataError(
            f"kernel width(s) {too_wide} exceed the sequence length {length}; a conv branch would have "
            f"no valid output position")
    return w


def peak_blocked_val_split(groups: list[str], *, val_fraction: float, seed: int):
    """Carve an inner validation set out of TRAIN by WHOLE GROUP, for early stopping.

    Group-blocked for the same measured reason every outer split in this project is: on tiles the group is
    the PEAK (tiles inside one peak are offset by tens of bp and overlap heavily), and on the designed grid
    it is the ELEMENT identity. A random inner split would early-stop against near-duplicates of its own
    training rows and systematically over-train. Membership is a pure function of the group label, so
    disjointness is BY CONSTRUCTION.

    The blocking variable is supplied by the caller (`fit_encoder(group_key=...)`) rather than hardcoded,
    so a new substrate must NAME its leakage unit instead of silently inheriting one that does not apply.

    Returns `(train_positions, val_positions)` as index lists into `groups`.
    """
    import random as _r
    uniq = sorted(set(groups))
    if len(uniq) < 2:
        raise EncoderDataError(f"need >= 2 training groups to carve a validation split, got {len(uniq)}")
    ps = list(uniq)
    _r.Random(10_000 + seed).shuffle(ps)
    k = max(1, min(len(ps) - 1, int(round(len(ps) * val_fraction))))
    held = set(ps[:k])
    tr = [i for i, p in enumerate(groups) if p not in held]
    va = [i for i, p in enumerate(groups) if p in held]
    if not tr or not va:
        raise EncoderDataError("degenerate inner split")
    return tr, va


def tile_peak_key(t) -> str:
    """Default grouping: a peak tile's peak id. REFUSES an object with no peak rather than defaulting.

    A silent `""` default would collapse every row into one group, and `peak_blocked_val_split` would then
    refuse with a confusing "need >= 2 groups" — or worse, a future edit adding a fallback would produce a
    RANDOM inner split on a substrate that needs a blocked one. Naming the failure here keeps it legible.
    """
    p = getattr(t, "peak", None)
    if not p:
        raise EncoderDataError(
            "this record has no non-empty `peak`, so the default tile grouping does not apply; pass an "
            "explicit group_key naming this substrate's leakage unit (e.g. the element identity on the "
            "designed grid)")
    return str(p)


def count_parameters(widths=DEFAULT_WIDTHS, *, channels: int = 32) -> int:
    """Parameter count, derived arithmetically so it is checkable without building the model.

    Per branch: `channels * 4 * width` weights + `channels` biases. Head: `feature_width + 1`.
    """
    conv = sum(channels * len(BASES) * w + channels for w in widths)
    head = encoder_feature_width(widths, channels=channels) + 1
    return conv + head


# ---------------------------------------------------------------------------------------------------
# the model + training -- torch lazy-imported
# ---------------------------------------------------------------------------------------------------
def build_encoder(widths=DEFAULT_WIDTHS, *, channels: int = 32):
    """The `SequenceEncoder` module. Imports torch, so it is not importable on a default install."""
    import torch
    from torch import nn

    class SequenceEncoder(nn.Module):
        """One-hot -> parallel Conv1d branches -> ReLU -> global max+mean pool -> concat -> linear.

        Pooling collapses the length axis per branch, which is what makes branches of DIFFERENT kernel
        widths concatenable at all. Max pool answers "is this pattern present anywhere" (a motif detector
        is translation-invariant by nature); mean pool answers "how much of it is there" (a composition
        signal). Both are kept because the gate's finding is precisely that composition carries this
        substrate and position does not — so the model must be able to express composition cleanly, or a
        negative would be about the architecture starving the winning signal rather than about position.
        """

        def __init__(self) -> None:
            super().__init__()
            self.widths = tuple(widths)
            self.branches = nn.ModuleList(
                [nn.Conv1d(len(BASES), channels, kernel_size=w) for w in self.widths])
            self.head = nn.Linear(encoder_feature_width(self.widths, channels=channels), 1)

        def forward(self, x):
            pooled = []
            for conv in self.branches:
                h = torch.relu(conv(x))
                pooled.append(h.amax(dim=2))
                pooled.append(h.mean(dim=2))
            return self.head(torch.cat(pooled, dim=1)).squeeze(1)

    return SequenceEncoder()


class _Determinism:
    """Set CPU determinism for the duration of a fit, then RESTORE it.

    Restoring matters: `use_deterministic_algorithms` is a GLOBAL torch flag, and leaving it on would
    change the behaviour of every later test in the same process — a fit that silently reconfigures the
    interpreter is a cross-test contaminant, not a determinism guarantee.
    """

    def __init__(self, enable: bool) -> None:
        self.enable = enable
        self.prev = None

    def __enter__(self):
        if not self.enable:
            return self
        import torch
        self.prev = torch.are_deterministic_algorithms_enabled()
        torch.use_deterministic_algorithms(True)
        return self

    def __exit__(self, *exc):
        if self.enable and self.prev is not None:
            import torch
            torch.use_deterministic_algorithms(self.prev)
        return False


def resolve_device(device: str) -> str:
    """`cpu` | `cuda` | `auto`. CPU is the determinism REFERENCE; CUDA convolutions may be nondeterministic."""
    if device not in ("cpu", "cuda", "auto"):
        raise EncoderDataError(f"device must be cpu/cuda/auto, got {device!r}")
    if device == "cpu":
        return "cpu"
    try:
        import torch
    except ImportError as e:
        raise EncoderDataError("torch is not installed; device must be 'cpu' only for pure helpers") from e
    if device == "cuda":
        if not torch.cuda.is_available():
            raise EncoderDataError("device='cuda' requested but torch reports no CUDA device")
        return "cuda"
    return "cuda" if torch.cuda.is_available() else "cpu"


def fit_encoder(train, test, *, widths=DEFAULT_WIDTHS, protocol: TrainingProtocol = REGISTERED_PROTOCOL,
                seed: int = 0, device: str = "cpu", shuffle_labels: bool = False,
                group_key=tile_peak_key) -> dict:
    """Train the encoder on `train`, score Spearman on `test`. Mirrors `_fit`'s role in the gate scripts.

    `shuffle_labels` permutes the TRAINING targets only — the validity null. If a model with its labels
    destroyed still scores, the pipeline fabricates signal and nothing may be graded. The test labels are
    never touched, so the null measures what it claims to.

    Early stopping runs against a PEAK-BLOCKED inner split of `train`; `test` is never seen until the
    final score. The returned dict records the protocol actually used, so Step 7's gate can assert it
    rather than trust it.
    """
    import numpy as np
    import torch
    from torch import nn

    from .genomewide import spearman

    if not train:
        raise EncoderDataError("empty training set")
    if not test:
        raise EncoderDataError("empty test set")

    dev = resolve_device(device)
    seqs_tr = [t.seq for t in train]
    L = require_uniform_length(seqs_tr)
    # the test set must share the length -- a different length would mean a different substrate
    if require_uniform_length([t.seq for t in test]) != L:
        raise EncoderDataError("train and test sequence lengths differ; refusing to pad or truncate")
    w = validate_widths(widths, L)

    groups = [group_key(t) for t in train]
    tr_pos, va_pos = peak_blocked_val_split(groups, val_fraction=protocol.val_fraction, seed=seed)

    X_all = one_hot(seqs_tr)
    y_all = np.log1p(np.asarray([t.expression for t in train], dtype=np.float64))
    if shuffle_labels:
        rng = np.random.default_rng(20_000 + seed)
        y_all = rng.permutation(y_all)
    # standardise the target on the INNER TRAIN rows only -- not on all of train, and never on test
    mu = float(y_all[tr_pos].mean())
    sd = float(y_all[tr_pos].std()) or 1.0
    y_all = (y_all - mu) / sd

    Xtr = torch.from_numpy(X_all[tr_pos])
    ytr = torch.from_numpy(y_all[tr_pos].astype(np.float32))
    Xva = torch.from_numpy(X_all[va_pos]).to(dev)
    yva = torch.from_numpy(y_all[va_pos].astype(np.float32)).to(dev)
    Xte = torch.from_numpy(one_hot([t.seq for t in test])).to(dev)
    yte_raw = [t.expression for t in test]

    with _Determinism(dev == "cpu"):
        torch.manual_seed(seed)
        model = build_encoder(w, channels=protocol.channels).to(dev)
        if protocol.optimizer != "adam":
            raise EncoderDataError(f"unsupported optimizer {protocol.optimizer!r}")
        opt = torch.optim.Adam(model.parameters(), lr=protocol.lr,
                              weight_decay=protocol.weight_decay)
        lossf = nn.MSELoss()
        gen = torch.Generator().manual_seed(seed)

        best = float("inf")
        best_state = None
        bad_epochs = 0
        epochs_ran = 0
        n = Xtr.shape[0]
        for epoch in range(protocol.max_epochs):
            model.train()
            order = torch.randperm(n, generator=gen)
            for s in range(0, n, protocol.batch_size):
                b = order[s:s + protocol.batch_size]
                xb, yb = Xtr[b].to(dev), ytr[b].to(dev)
                opt.zero_grad(set_to_none=True)
                lossf(model(xb), yb).backward()
                opt.step()
            model.eval()
            with torch.no_grad():
                vloss = float(lossf(model(Xva), yva))
            epochs_ran = epoch + 1
            if vloss < best - 1e-6:
                best, bad_epochs = vloss, 0
                best_state = {k: v.detach().clone() for k, v in model.state_dict().items()}
            else:
                bad_epochs += 1
                if bad_epochs >= protocol.early_stop_patience:
                    break
        if best_state is not None:
            model.load_state_dict(best_state)
        model.eval()
        with torch.no_grad():
            pred = [float(v) for v in model(Xte).cpu()]

    return {
        "spearman": spearman(pred, yte_raw),
        # the raw test predictions, so a caller can compute a SECOND framing (e.g. AUROC on an
        # active/inactive column) without retraining. Callers must drop this before serialising -- it is
        # one float per test row.
        "predictions": pred,
        "n_train": len(tr_pos),
        "n_inner_val": len(va_pos),
        "n_test": len(test),
        "epochs_ran": epochs_ran,
        "early_stopped": epochs_ran < protocol.max_epochs,
        "best_inner_val_mse": round(best, 6),
        "device": dev,
        "widths": list(w),
        "n_parameters": count_parameters(w, channels=protocol.channels),
        "shuffle_labels": shuffle_labels,
        "protocol": protocol.as_dict(),
        "determinism_reference": dev == "cpu",
        "n_train_groups": len(set(groups)),
    }
