"""The five shipped PromoterAtlas expression heads: how to load them, and how to run them over a
latent that is the WRONG LENGTH for them.

The frozen base encoder reaches median Spearman 0.5076 on the tile-chain rung by masked-mean-pooling
its final latent. The clone ships five trained prediction heads beside that base, and none has ever
been scored on this rung. This module owns every fact about them so the feasibility probe and the
scorer cannot drift apart -- the sibling arms' defect history is why this is one module and not two
copies of the same arithmetic.

**MEASURED, not inherited (2026-11-02, re-asserted live by the probe).** All FIVE heads -- not four,
as an earlier note had it -- are identical in the two ways that matter: every `fc1.weight` is
`(64, 11008)`, and every head's `base_model.*` tensors are byte-identical to `promoteratlas-base.pt`
at 124 of 124 tensors, maximum absolute deviation exactly `0.000e+00`. That identity is what makes
this a FROZEN-ENCODER arm and what lets one forward pass feed all five heads. If a head ever
diverges it stops being a frozen-encoder head and becomes a fine-tuned-encoder one, which is a
different claim and a different arm, so `load_head` REFUSES rather than proceeding.

**The obstruction is arithmetic.** `11008 = 128 * 86`, so a head's `fc1` expects an 86-position
latent while this substrate is 150 positions (which would flatten to `128 * 150 = 19200`). A head
therefore cannot be fed a tile directly; it must be run over 86-wide WINDOWS, and
`150 - 86 + 1 = 65` of them fit. Which windowing is used is an acceptance-bar decision frozen from a
measurement, not a preference, which is why this module offers BOTH routes and never picks one:

* `head_forward_full` runs the convolutional trunk ONCE over the full 150-position latent and only
  then slices 86-wide windows before `fc1`. Cheap -- one trunk evaluation per tile.
* `head_forward_per_window` slices the latent FIRST and recomputes the trunk per window. Expensive
  -- 65 trunk evaluations per tile -- and it is what the head would actually see if handed an 86 bp
  sequence.

**They are not the same function, and that is the point.** `conv1`/`conv2` use `padding='same'`, so
near a window's edge the cheap route sees real neighbouring context where the true route sees zero
padding. The two agree in a window's interior and differ at its edges. A test asserting they are
equal everywhere would be asserting a bug.

**TRAP 1 -- `padding='same'` on an EVEN kernel is ASYMMETRIC, and guessing the wrong way is silent.**
`kernel_size=8` means 7 total padding, and torch splits it `left=3, right=4` (floor on the left).
Measured 2026-11-02 against `F.conv1d(..., padding="same")` with an impulse input: `(3, 4)` matches
and `(4, 3)` does not. A symmetric-looking `(4, 3)` produces a plausible tensor of the right shape
and wrong values, so it would not raise -- it would just move the arm's number.

**TRAP 2 -- the parameter count depends on which question you ask, and the repo's gate asks one of
them.** The constants below are STATE-DICTIONARY ELEMENT counts. The base's comparable figure is
**1,847,492**, which is NOT the `atlas_embed.EXPECTED_PARAMS = 1,847,428` that gate enforces as a
hard refusal. Both are right: the 64-element gap is the eight
`blocks.N.rotary_attention_block.rotary_emb.inv_freq` tensors at 8 elements each, which are
registered BUFFERS, and `EXPECTED_PARAMS` counts `sum(p.numel() for p in model.parameters())`, which
excludes buffers. A probe that sums the state dictionary and compares to `EXPECTED_PARAMS` reports a
checkpoint mismatch that does not exist.

**TRAP 3 -- a broken reduction must not read as a clean null.** `genomewide.spearman` maps NaN to
0.0, so a head output carrying NaN would be scored as "no signal" rather than as an error. Every
entry point here refuses non-finite input up front, exactly as `positional_readout` does.

**TRAP 4 -- in float32 the two routes differ EVERYWHERE, well below the edge effect, and a bar whose
tolerance sits under that floor would be measuring arithmetic noise.** `F.conv1d` selects different
blocking for a 150-position input than for an 86-position one, so the two routes disagree by order
1e-7 even at columns that are mathematically identical; `MaxPool1d` then amplifies it wherever two
channels are nearly tied, because a 1e-7 difference flips the argmax. Measured 2026-11-02: at a
1e-9 threshold ALL 86 columns differ in float32, while in float64 the differing set collapses to
exactly `contaminated_columns(...)` -- zero extra, zero missing, at offsets 0, 30 and 64. So the
edge effect is the deterministic signal, the 1e-7 floor is noise, and the probe must report both so
the bar can place its tolerance above the floor rather than inside it.

**Loading is deliberately routed through `atlas_embed.audit_checkpoint`.** That function
disassembles every pickle with `pickletools.genops` WITHOUT executing it and refuses any global
outside its allowlist; only then is the checkpoint loaded. Measured: all five heads PASS that audit
unchanged, each naming six globals, every one already covered -- so no allowlist edit is needed. The
tempting alternative, `weights_only=True` with `torch.serialization.safe_globals`, provably CANNOT
work on these files: torch keys its allowlist on `f"{__module__}.{__name__}"`, these were pickled
under numpy 1.x as `numpy.core.multiarray.scalar`, and numpy 2.4.4 exposes that object at
`numpy._core.multiarray`, so the name the unpickler asks for and the name the allowlist registers
never meet.

**The clone's own loader is UNUSABLE here and must not be reached for.**
`promoter_atlas.models.expression_predictor.ExpressionPredictor.__init__` calls `torch.load` with no
`weights_only` argument, and torch 2.6 changed that default to `True`, so instantiating it raises
`UnpicklingError` on the embedded Linux `PosixPath` before anything else happens -- independently of
a second defect in the same two lines, which is that its base-weights path is RELATIVE and only
resolves with the current directory set to the clone root. Measured 2026-11-02 under torch 2.7.1.
"""

from __future__ import annotations

import os
import pathlib
from dataclasses import dataclass, field

#: Positions a head's `fc1` expects. `11008 = 128 * 86` is the whole obstruction.
HEAD_SEQ_LENGTH = 86
#: Channels in the base encoder's final latent, i.e. `DNATransformer(dim=128)`.
LATENT_CHANNELS = 128
#: `fc1.weight` must be exactly this shape on every head, or it is not one of these heads.
FC1_SHAPE = (64, LATENT_CHANNELS * HEAD_SEQ_LENGTH)

CONV_KERNEL = 8
#: Measured against `F.conv1d(..., padding="same")`, NOT assumed. See TRAP 1.
SAME_PAD_LEFT, SAME_PAD_RIGHT = 3, 4
POOL_KERNEL, POOL_STRIDE, POOL_PADDING = 7, 1, 3

#: How many pooled columns a ZERO-PADDED window boundary contaminates, in closed form.
#:
#: Derivation, left boundary: `conv1` at column `c` reads input `[c - PAD_L, c + PAD_R]`, so it
#: differs iff `c < PAD_L`; `conv2` reads conv1's output over the same span, extending that to
#: `c < 2 * PAD_L`; the pool reads `[c - POOL_PADDING, c + POOL_PADDING]`, extending it again.
#: The right boundary is the same argument with `PAD_R`, which is why the two are NOT equal --
#: the asymmetric `same` split makes the right edge two columns wider than the left.
#:
#: MEASURED against the real functions in float64: the differing-column set equals this prediction
#: EXACTLY, zero extra and zero missing, at window offsets 0, 30 and 64.
EDGE_COLS_LEFT = 2 * SAME_PAD_LEFT + POOL_PADDING      # 9
EDGE_COLS_RIGHT = 2 * SAME_PAD_RIGHT + POOL_PADDING    # 11

#: STATE-DICTIONARY ELEMENT counts, not parameter counts. See TRAP 2 before comparing either of
#: these to `atlas_embed.EXPECTED_PARAMS`.
EXPECTED_HEAD_ELEMENTS = 2_814_533
EXPECTED_HEAD_ONLY = 967_041
#: The base's state-dictionary element count. `EXPECTED_PARAMS` is 1,847,428; the 64-element gap is
#: eight 8-element rotary `inv_freq` BUFFERS, which `model.parameters()` does not count.
BASE_STATE_DICT_ELEMENTS = 1_847_492
BASE_STATE_DICT_TENSORS = 124

#: Exactly the keys a head carries above the base. Any other extra key means a different head shape.
HEAD_ONLY_KEYS = (
    "conv1.weight", "conv1.bias", "conv2.weight", "conv2.bias",
    "fc1.weight", "fc1.bias", "fc2.weight", "fc2.bias",
)

_BASE_PREFIX = "base_model."

#: Endpoint tokens the clone uses in its filenames. `trspred` predicts TRANSCRIPTION; `tslpred`
#: predicts TRANSLATION, which is a DIFFERENT QUANTITY from promoter activity -- recorded as a
#: stated limit on that head rather than as grounds for dropping it.
TRANSCRIPTION, TRANSLATION = "trspred", "tslpred"


class AtlasHeadsError(RuntimeError):
    """A head cannot be loaded, failed its audit, or was asked for something it refuses to guess."""


@dataclass(frozen=True)
class Head:
    """One loaded head. `endpoint` and `log_transform` are DERIVED, never restated."""

    name: str
    path: str
    endpoint: str
    state_dict: dict
    normalization: dict = field(default_factory=dict)

    @property
    def log_transform(self) -> bool:
        """True only where the checkpoint itself says so (`tslpred-kosuri2013` does; the others
        do not). Spearman is invariant under the monotone denormalisation, so this cannot move a
        rank-based verdict -- but an ABSOLUTE predicted value would be wrong without it."""
        return bool(self.normalization.get("log_transform", False))

    @property
    def predicts_transcription(self) -> bool:
        return self.endpoint == TRANSCRIPTION


# ---------------------------------------------------------------------------------------------------
# endpoint + discovery -- derived from what is on disk, not from a hand-maintained table
# ---------------------------------------------------------------------------------------------------
def head_endpoint(filename: str) -> str:
    """Derive the endpoint from the clone's own filename token. Refuses an unrecognised one."""
    stem = os.path.basename(str(filename))
    for token in (TRANSCRIPTION, TRANSLATION):
        if "-%s-" % token in stem:
            return token
    raise AtlasHeadsError(
        "cannot derive an endpoint from %r: it names neither %r nor %r, so whether it predicts "
        "transcription or translation is unknown and must not be guessed."
        % (stem, TRANSCRIPTION, TRANSLATION)
    )


def head_paths(expression_dir: pathlib.Path | str) -> list[pathlib.Path]:
    """Every head checkpoint in the clone's `trained_weights/expression_prediction` directory.

    Discovered by glob rather than hand-listed, so a clone that ships a sixth head is SEEN instead
    of silently ignored -- the hand-enumerated-list failure mode this repo has hit five times.
    """
    d = pathlib.Path(expression_dir)
    if not d.is_dir():
        raise AtlasHeadsError("expression-head directory absent: %s" % d)
    found = sorted(d.glob("*.pt"))
    if not found:
        raise AtlasHeadsError("no .pt checkpoints under %s" % d)
    return found


# ---------------------------------------------------------------------------------------------------
# loading -- audit statically, THEN load; never the clone's own loader
# ---------------------------------------------------------------------------------------------------
def _alias_posixpath_on_windows() -> None:
    """Make a Linux-saved `PosixPath` instantiable here.

    Idempotent and deliberately NOT conditional on the alias being absent: `atlas_embed.load_atlas`
    performs the same rebinding process-globally and never restores it, so this module must not
    assume `pathlib.PosixPath` is still its real self when it is called.
    """
    if os.name == "nt":
        pathlib.PosixPath = pathlib.WindowsPath  # type: ignore[misc, assignment]


def load_base_state_dict(pa_repo: pathlib.Path | str | None = None) -> dict:
    """The frozen base encoder's state dictionary, audited before it is unpickled."""
    import torch

    from dna_decode.glm.atlas_embed import DEFAULT_PA_REPO, audit_checkpoint

    pa = pathlib.Path(pa_repo if pa_repo is not None else DEFAULT_PA_REPO)
    weights = pa / "trained_weights/base_model/promoteratlas-base.pt"
    if not weights.is_file():
        raise AtlasHeadsError("base checkpoint absent: %s" % weights)
    audit_checkpoint(weights)
    _alias_posixpath_on_windows()
    ck = torch.load(weights, map_location="cpu", weights_only=False)
    return ck["model_state_dict"] if "model_state_dict" in ck else ck


def base_identity(head_state_dict: dict, base_state_dict: dict) -> dict:
    """How far a head's embedded base weights have drifted from the frozen base.

    Pure: takes two dictionaries, touches no disk. Returns `n_checked`, `n_identical`, `max_abs_dev`
    and `missing` so a caller can RECORD the measurement rather than only act on a boolean.
    """
    embedded = {k[len(_BASE_PREFIX):]: v for k, v in head_state_dict.items()
                if k.startswith(_BASE_PREFIX)}
    n_identical, max_dev, missing = 0, 0.0, []
    for k, v in embedded.items():
        if k not in base_state_dict:
            missing.append(k)
            continue
        dev = float((v.float() - base_state_dict[k].float()).abs().max().item())
        max_dev = max(max_dev, dev)
        if dev == 0.0:
            n_identical += 1
    return {"n_checked": len(embedded), "n_identical": n_identical,
            "max_abs_dev": max_dev, "missing": sorted(missing)}


def load_head(path: pathlib.Path | str, *, base_state_dict: dict) -> Head:
    """Load one head, refusing anything that would change what the arm is measuring.

    Three refusals, each because the alternative is a silently wrong number rather than an error:
    a wrong `fc1.weight` shape (not one of these heads at all), an unexpected head-only key set
    (a different architecture), and ANY drift in the embedded base weights (a fine-tuned-encoder
    head masquerading as a frozen-encoder one).
    """
    import torch

    from dna_decode.glm.atlas_embed import audit_checkpoint

    p = pathlib.Path(path)
    if not p.is_file():
        raise AtlasHeadsError("head checkpoint absent: %s" % p)
    endpoint = head_endpoint(p.name)

    audit_checkpoint(p)
    _alias_posixpath_on_windows()
    ck = torch.load(p, map_location="cpu", weights_only=False)
    sd = ck["model_state_dict"] if "model_state_dict" in ck else ck

    fc1 = sd.get("fc1.weight")
    if fc1 is None or tuple(fc1.shape) != FC1_SHAPE:
        raise AtlasHeadsError(
            "%s has fc1.weight %r, the arm requires %r (= 64 x 128 x %d). A head of a different "
            "geometry is not one of these heads."
            % (p.name, None if fc1 is None else tuple(fc1.shape), FC1_SHAPE, HEAD_SEQ_LENGTH)
        )

    extra = tuple(sorted(k for k in sd if not k.startswith(_BASE_PREFIX)))
    if extra != tuple(sorted(HEAD_ONLY_KEYS)):
        raise AtlasHeadsError(
            "%s carries head-only keys %r; the arm expects exactly %r. A different key set is a "
            "different architecture." % (p.name, extra, tuple(sorted(HEAD_ONLY_KEYS)))
        )

    ident = base_identity(sd, base_state_dict)
    if ident["missing"] or ident["n_identical"] != ident["n_checked"]:
        raise AtlasHeadsError(
            "%s: embedded base weights are NOT byte-identical to the frozen base "
            "(%d/%d identical, max |dev| %.3e, missing %r). This head would make the arm a "
            "FINE-TUNED-encoder arm, which is a different claim and a different arm."
            % (p.name, ident["n_identical"], ident["n_checked"], ident["max_abs_dev"],
               ident["missing"])
        )

    norm = dict(ck.get("normalization_params") or {})
    return Head(name=p.stem, path=str(p), endpoint=endpoint, state_dict=sd, normalization=norm)


def load_heads(expression_dir: pathlib.Path | str, *, base_state_dict: dict) -> list[Head]:
    """Every head in the directory, in sorted order. One bad head fails the whole set."""
    return [load_head(p, base_state_dict=base_state_dict) for p in head_paths(expression_dir)]


# ---------------------------------------------------------------------------------------------------
# the windowing arithmetic -- pure, and the one place the 65 comes from
# ---------------------------------------------------------------------------------------------------
def window_starts(length: int, width: int = HEAD_SEQ_LENGTH) -> tuple[int, ...]:
    """Every start offset of a `width`-wide window inside `length` positions, stride 1.

    `window_starts(150, 86)` is 65 windows, 0 through 64. Refuses a window WIDER than the latent
    rather than truncating or padding to fit, because either silently changes what the head is fed.
    """
    if not isinstance(length, int) or not isinstance(width, int):
        raise AtlasHeadsError("length and width must be ints, got %r and %r"
                              % (type(length).__name__, type(width).__name__))
    if width <= 0 or length <= 0:
        raise AtlasHeadsError("length and width must be positive, got %d and %d" % (length, width))
    if width > length:
        raise AtlasHeadsError(
            "a %d-wide window does not fit in %d positions. Refusing rather than padding: a padded "
            "window feeds the head zeros it was never trained on." % (width, length)
        )
    return tuple(range(length - width + 1))


def n_windows(length: int, width: int = HEAD_SEQ_LENGTH) -> int:
    return len(window_starts(length, width))


def contaminated_columns(start: int, length: int,
                         width: int = HEAD_SEQ_LENGTH) -> frozenset[int]:
    """Which columns of the window at `start` the two routes MUST disagree on, in closed form.

    A boundary contaminates columns only when it is NOT shared with the full latent's own boundary.
    At `start == 0` the left edge is shared -- both routes see the same zero padding there -- so only
    the right edge differs; at `start == length - width` only the left edge differs; in between,
    both do. That asymmetry is the whole mechanical content of the cheap-versus-true comparison, and
    it is why the degenerate `length == width` case makes the two routes coincide exactly.

    Verified against the real functions in float64 at offsets 0, 30 and 64: predicted == observed.
    """
    starts = window_starts(length, width)
    if start not in starts:
        raise AtlasHeadsError("start %d is not a valid window offset for length %d width %d"
                              % (start, length, width))
    cols: set[int] = set()
    if start > 0:
        cols |= set(range(0, min(EDGE_COLS_LEFT, width)))
    if start + width < length:
        cols |= set(range(max(0, width - EDGE_COLS_RIGHT), width))
    return frozenset(cols)


def _require_finite(t, name: str):
    """Refuse non-finite input. See TRAP 3: NaN would be scored as a clean null, not as an error."""
    import torch

    if not torch.isfinite(t).all():
        raise AtlasHeadsError(
            "%s carries non-finite values. Refusing: the scorer maps NaN to 0.0, so this would be "
            "reported as 'no signal' rather than as a fault." % name
        )
    return t


def _trunk(latent, sd):
    """`conv1 -> relu -> conv2 -> relu -> pool` over whatever length it is given.

    Length-preserving: `padding='same'` keeps it, and `MaxPool1d(7, stride=1, padding=3)` keeps it
    (`L + 2*3 - 7 + 1 == L`). The pad is applied EXPLICITLY at the measured `(3, 4)` split rather
    than via `padding='same'`, so the asymmetry is visible in the code instead of being a library
    default this module merely hopes is stable.
    """
    import torch.nn.functional as F

    x = latent
    for w, b in (("conv1", "conv1"), ("conv2", "conv2")):
        x = F.conv1d(F.pad(x, (SAME_PAD_LEFT, SAME_PAD_RIGHT)),
                     sd["%s.weight" % w], sd["%s.bias" % b])
        x = F.relu(x)
    return F.max_pool1d(x, kernel_size=POOL_KERNEL, stride=POOL_STRIDE, padding=POOL_PADDING)


def _classify(pooled_window, sd):
    """`flatten -> fc1 -> relu -> fc2` on one 86-wide pooled window -> (B,).

    `flatten(1)` on `(B, 128, 86)` is CHANNEL-MAJOR in C order, which is why the window must be
    sliced on the POSITION axis before flattening and never after.
    """
    import torch.nn.functional as F

    flat = pooled_window.flatten(1)
    if flat.shape[1] != FC1_SHAPE[1]:
        raise AtlasHeadsError(
            "flattened window is %d wide, fc1 expects %d. The window was sliced on the wrong axis "
            "or at the wrong width." % (flat.shape[1], FC1_SHAPE[1])
        )
    return F.linear(F.relu(F.linear(flat, sd["fc1.weight"], sd["fc1.bias"])),
                    sd["fc2.weight"], sd["fc2.bias"]).squeeze(-1)


def head_forward_full(latent, sd, *, width: int = HEAD_SEQ_LENGTH):
    """THE CHEAP ROUTE. One trunk evaluation over the full latent, then slice windows before `fc1`.

    `latent` is `(B, 128, L)`; returns `(B, n_windows)`. Each window therefore sees REAL
    neighbouring context at its edges, because the convolutions ran before the slice.
    """
    import torch

    _require_finite(latent, "latent")
    pooled = _trunk(latent, sd)
    cols = [_classify(pooled[:, :, s:s + width], sd)
            for s in window_starts(int(latent.shape[-1]), width)]
    return _require_finite(torch.stack(cols, dim=-1), "head_forward_full output")


def head_forward_per_window(latent, sd, *, width: int = HEAD_SEQ_LENGTH):
    """THE TRUE ROUTE. Slice the latent first, then recompute the trunk per window.

    `latent` is `(B, 128, L)`; returns `(B, n_windows)`. This is what the head would see if handed
    an 86 bp sequence: each window's edges are ZERO-PADDED, not filled with neighbouring context.
    Costs `n_windows` trunk evaluations per tile instead of one.
    """
    import torch

    _require_finite(latent, "latent")
    cols = [_classify(_trunk(latent[:, :, s:s + width], sd), sd)
            for s in window_starts(int(latent.shape[-1]), width)]
    return _require_finite(torch.stack(cols, dim=-1), "head_forward_per_window output")


def reduce_windows(values, how: str):
    """Reduce `(B, n_windows)` over the window axis -> `(B,)`.

    `mean` is the arm's declared primary; `max` is a pre-declared secondary that the bar marks
    INELIGIBLE as the arm's result. Refuses any other name rather than defaulting to one, because a
    silent default would decide a declared primary by accident.
    """
    import torch

    _require_finite(values, "window values")
    if values.ndim != 2:
        raise AtlasHeadsError("expected (B, n_windows), got shape %r" % (tuple(values.shape),))
    if how == "mean":
        return values.mean(dim=1)
    if how == "max":
        return values.max(dim=1).values
    raise AtlasHeadsError("unknown reduction %r; the bar declares only 'mean' and 'max'" % (how,))


def route_deviation(cheap, true_) -> dict:
    """Summarise how far the cheap route sits from the true one. Pure; no scoring, no verdict.

    Reports max and median absolute deviation plus the correlation of the two, so an acceptance bar
    can CHOOSE a route from a number instead of from a preference.
    """
    import torch

    _require_finite(cheap, "cheap-route values")
    _require_finite(true_, "true-route values")
    if tuple(cheap.shape) != tuple(true_.shape):
        raise AtlasHeadsError("shape mismatch: cheap %r vs true %r"
                              % (tuple(cheap.shape), tuple(true_.shape)))
    d = (cheap.double() - true_.double()).abs()
    a, b = cheap.double().reshape(-1), true_.double().reshape(-1)
    if a.numel() < 2 or float(a.std()) == 0.0 or float(b.std()) == 0.0:
        pearson = None
    else:
        pearson = float(torch.corrcoef(torch.stack([a, b]))[0, 1])
    return {"n": int(d.numel()), "max_abs_dev": float(d.max()),
            "median_abs_dev": float(d.median()), "mean_abs_dev": float(d.mean()),
            "pearson": pearson}
