"""Load the frozen PromoterAtlas base encoder and masked-mean-pool it to 128 dims, LOCALLY.

The arm that produced `arms.native150 = 0.5076` ran this computation inside
`scripts/kaggle/promoteratlas_tile_arm_kernel.py` on a Kaggle T4. The pooling block is lifted here
verbatim in behaviour so the local and Kaggle paths keep ONE definition of the representation; a
second, subtly different pooling would make every local number incomparable with the committed
tables, which is exactly what the reconcile gate exists to catch.

**It runs on CPU and that is sized by measurement, not estimate.** Measured 2026-11-01 on this host:
1,847,428 parameters, 8 threads, batch 64 at 150 bp -> 0.544 s -> 117.6 seq/s, so a full pass over
all 44,106 tiles is about 6.3 minutes. `torch` reports `cuda False` here. No Kaggle push, no GPU, no
network.

**`weights_only=False` is justified by a committed audit, not by convenience.** The checkpoint embeds
a `pathlib.PosixPath` (an unused training-config path), which is off torch's allowlist, so
`weights_only=True` cannot load it. `scripts/promoteratlas_pickle_audit.py` disassembles the pickle
WITHOUT executing it and found exactly four referenced globals, all inside a torch/pathlib/stdlib
allowlist. That audit is re-run here before every load rather than cited.
"""

from __future__ import annotations

import importlib.util
import io
import os
import pathlib
import pickletools
import zipfile
from pathlib import Path

#: Pinned by the arm's acceptance bar. A different count means a different checkpoint.
EXPECTED_PARAMS = 1_847_428

#: Nucleotide -> channel. A masked base in pretraining was all-four-channels-zero, which is why
#: zero-padding is NOT a neutral operation on this encoder.
NUC = {"A": 0, "C": 1, "G": 2, "T": 3}

SAFE_GLOBAL_PREFIXES = (
    "torch", "collections.OrderedDict", "pathlib.PosixPath", "pathlib.WindowsPath",
    "numpy.core.multiarray", "numpy", "builtins.set", "_codecs.encode",
)

#: Where the MIT-licensed clone lives on this host. Overridable by env so a different checkout works.
DEFAULT_PA_REPO = Path(
    os.environ.get("PA_REPO", r"C:/Users/LAPTOP~1/AppData/Local/Temp/pa_probe/PromoterAtlas")
)


class AtlasEmbedError(RuntimeError):
    """The encoder cannot be loaded, or the checkpoint failed its audit."""


def audit_checkpoint(path: Path) -> list[str]:
    """Disassemble the checkpoint's pickles WITHOUT executing them; refuse an off-allowlist global.

    Returns the sorted global list so a caller can record it in an artifact.
    """
    found: dict[str, int] = {}
    with zipfile.ZipFile(path) as z:
        for pk in [n for n in z.namelist() if n.endswith(".pkl")]:
            strings: list[str] = []
            for op, arg, _pos in pickletools.genops(io.BytesIO(z.read(pk))):
                if op.name in ("SHORT_BINUNICODE", "BINUNICODE", "UNICODE"):
                    strings.append(arg)
                elif op.name == "STACK_GLOBAL" and len(strings) >= 2:
                    found["%s.%s" % (strings[-2], strings[-1])] = 1
                elif op.name in ("GLOBAL", "INST"):
                    found[str(arg).replace(" ", ".")] = 1
    bad = [g for g in found if not any(g == p or g.startswith(p) for p in SAFE_GLOBAL_PREFIXES)]
    if bad:
        raise AtlasEmbedError(
            "checkpoint names off-allowlist globals %s; refusing to load" % sorted(bad)
        )
    return sorted(found)


def load_atlas(pa_repo: Path | str = DEFAULT_PA_REPO):
    """Load the frozen base encoder in eval mode. Refuses on a parameter-count mismatch."""
    import torch

    pa = Path(pa_repo)
    src = pa / "src/promoter_atlas/models/dna_transformer.py"
    weights = pa / "trained_weights/base_model/promoteratlas-base.pt"
    for p in (src, weights):
        if not p.is_file():
            raise AtlasEmbedError(
                "PromoterAtlas clone incomplete: %s absent (PA_REPO=%s)" % (p, pa)
            )

    spec = importlib.util.spec_from_file_location("pa_dt_local", src)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)

    audit_checkpoint(weights)
    # The checkpoint embeds a pathlib.PosixPath. On Windows, instantiating one raises
    # NotImplementedError, so alias it HERE and only here. The kernel carries the same guard; the
    # dry run caught a first version that saved the original and never wrote the assignment.
    if os.name == "nt":
        pathlib.PosixPath = pathlib.WindowsPath
    ck = torch.load(weights, map_location="cpu", weights_only=False)
    model = mod.DNATransformer(dim=128, n_blocks=8, n_heads=8)
    model.load_state_dict(ck["model_state_dict"], strict=True)
    n = sum(p.numel() for p in model.parameters())
    if n != EXPECTED_PARAMS:
        raise AtlasEmbedError("%d params, the arm's bar expects %d" % (n, EXPECTED_PARAMS))
    model.eval()
    return model


def masked_mean_pool(model, seqs, *, width=None, device="cpu", batch=64, progress=None):
    """Masked mean-pool of the final latent over the UNPADDED positions -> (N, 128) float32.

    `width is None` means the NATIVE geometry (`width == len(seq)`), which is the arm's primary cell.
    A `width` greater than the sequence length zero-pads, which is how pretraining encoded a MASKED
    base -- measured to perturb the pooled embedding even under masked pooling, because attention
    lets those columns reach the real positions. Available so the declared secondary cell stays
    reproducible; never the default.

    All sequences in one call must share a length: a ragged batch is REFUSED rather than padded,
    because a silent pad would produce a plausible meaningless number.
    """
    import numpy as np
    import torch

    if not seqs:
        raise AtlasEmbedError("masked_mean_pool got zero sequences")
    lengths = {len(s) for s in seqs}
    if len(lengths) != 1:
        raise AtlasEmbedError(
            "ragged input (%d distinct lengths, %d-%d); this function refuses to pad -- split by "
            "length or pass a fixed-geometry cohort" % (len(lengths), min(lengths), max(lengths))
        )
    L = lengths.pop()
    w = L if width is None else int(width)
    if w < L:
        raise AtlasEmbedError(
            "width %d is shorter than the sequences (%d bp); that would truncate" % (w, L)
        )

    out = np.zeros((len(seqs), 128), dtype=np.float32)
    with torch.no_grad():
        for i in range(0, len(seqs), batch):
            chunk = seqs[i:i + batch]
            X = torch.zeros((len(chunk), 4, w))
            for b, s in enumerate(chunk):
                for j, c in enumerate(s):
                    k = NUC.get(c)
                    if k is not None:
                        X[b, k, j] = 1.0
            X = X.to(device)
            _, _, z = model(X, return_latent=True)
            out[i:i + len(chunk)] = z[:, :, :L].mean(dim=2).float().cpu().numpy()
            if progress is not None:
                progress(min(i + batch, len(seqs)), len(seqs))
    return out


def ridge_spearman(Xtr, ytr_raw, Xte, yte_raw, *, alpha=1.0, shuffle_seed=None):
    """Mirrors `dna_decode.glm.expression.fit_predict_ridge` EXACTLY.

    StandardScaler on the train matrix, Ridge(alpha=1.0), fit on `log1p(expression)`, score Spearman
    against RAW expression. Any deviation makes the number incomparable with the committed incumbent
    tables, which is the entire point of the reconcile gate.
    """
    import numpy as np
    from sklearn.linear_model import Ridge
    from sklearn.preprocessing import StandardScaler

    from dna_decode.glm.genomewide import spearman

    ytr = np.log1p(np.asarray(ytr_raw, dtype=float))
    if shuffle_seed is not None:
        rng = np.random.default_rng(90_000 + shuffle_seed)
        ytr = ytr.copy()
        rng.shuffle(ytr)
    sc = StandardScaler().fit(Xtr)
    m = Ridge(alpha=alpha).fit(sc.transform(Xtr), ytr)
    pred = m.predict(sc.transform(Xte))
    return float(spearman(list(pred), list(yte_raw))), [float(v) for v in pred]
