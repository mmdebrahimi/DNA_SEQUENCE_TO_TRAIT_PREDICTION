"""iML1515 loading + FBA gene-KO -> growth/essentiality (the mechanistic edit->trait engine).

cobrapy is a LAZY import (heavy, optional extra) — the pure logic (`call_essential`) imports
with no deps so the honesty rails + unit tests run wheel-only.

Model resolution (offline-first, then network, then cache):
  1. $DNA_DECODE_IML1515  (explicit override)
  2. data/fba/iML1515.xml.gz  (dev checkout)
  3. <cache>/iML1515.xml.gz   (~/.cache/dna_decode/fba or $DNA_DECODE_CACHE)
  4. download from BiGG Models -> cache
"""
from __future__ import annotations

import os
import urllib.request
from pathlib import Path

MODEL_NAME = "iML1515"
_BIGG_URL = "http://bigg.ucsd.edu/static/models/iML1515.xml.gz"

# module-level model cache so repeated CLI calls in one process don't re-parse (~10s)
_MODEL = None


def _cache_dir() -> Path:
    env = os.environ.get("DNA_DECODE_CACHE")
    base = Path(env) if env else (Path.home() / ".cache" / "dna_decode")
    return base / "fba"


def resolve_model_path(download: bool = True) -> Path:
    """Locate iML1515.xml.gz; optionally download from BiGG into the cache. Raises if unresolved."""
    override = os.environ.get("DNA_DECODE_IML1515")
    if override and Path(override).exists():
        return Path(override)
    # dev checkout copy
    dev = Path(__file__).resolve().parent.parent.parent / "data" / "fba" / "iML1515.xml.gz"
    if dev.exists():
        return dev
    cached = _cache_dir() / "iML1515.xml.gz"
    if cached.exists():
        return cached
    if not download:
        raise FileNotFoundError(
            "iML1515 model not found. Set $DNA_DECODE_IML1515, place it at data/fba/iML1515.xml.gz, "
            "or allow the BiGG download."
        )
    cached.parent.mkdir(parents=True, exist_ok=True)
    urllib.request.urlretrieve(_BIGG_URL, cached)  # ~620 KB, free (BiGG Models)
    return cached


def load_model(path: str | Path | None = None, reuse: bool = True):
    """Load (and cache) the iML1515 cobra model. Lazy-imports cobra."""
    global _MODEL
    if reuse and _MODEL is not None and path is None:
        return _MODEL
    try:
        import cobra  # noqa: PLC0415  (lazy heavy import)
    except ImportError as e:  # pragma: no cover - env-dependent
        raise ImportError(
            "cobrapy is required for the FBA cell: `pip install 'dna-decode[fba]'` or `uv pip install cobra`."
        ) from e
    p = Path(path) if path else resolve_model_path()
    model = cobra.io.read_sbml_model(str(p))
    if path is None:
        _MODEL = model
    return model


def call_essential(growth: float | None, wt_growth: float, frac: float = 0.01) -> bool:
    """PURE: is a KO essential? growth < frac x wild-type (or ~0) -> essential.

    frac=0.01 (1% of WT) is the standard growth/no-growth cutoff for in-silico gene essentiality.
    A None / NaN / infeasible growth is essential (the KO broke the model's ability to grow).
    """
    if growth is None:
        return True
    if growth != growth:  # NaN
        return True
    if wt_growth <= 0:
        return growth <= 1e-9
    return growth < max(frac * wt_growth, 1e-6)


def wildtype_growth(model) -> float:
    """Max growth rate (/h) on the model's current medium (default: glucose M9 aerobic)."""
    return float(model.slim_optimize())


def knockout_growth(model, gene_ids) -> float:
    """Growth rate (/h) after knocking out one or more genes (context-managed; model unchanged)."""
    if isinstance(gene_ids, str):
        gene_ids = [gene_ids]
    with model:
        for g in gene_ids:
            model.genes.get_by_id(g).knock_out()
        val = model.slim_optimize()
    return 0.0 if (val is None or val != val) else float(val)


def gene_essentiality(model, frac: float = 0.01) -> dict[str, tuple[float, bool]]:
    """Genome-wide single-gene-deletion essentiality. Returns {gene_id: (growth, is_essential)}."""
    from cobra.flux_analysis import single_gene_deletion  # noqa: PLC0415

    wt = wildtype_growth(model)
    res = single_gene_deletion(model, processes=1).reset_index()
    idcol = "ids" if "ids" in res.columns else res.columns[0]
    out: dict[str, tuple[float, bool]] = {}
    for _, row in res.iterrows():
        gid = row[idcol]
        if isinstance(gid, (frozenset, set)):
            gid = next(iter(gid)) if gid else None
        if gid is None:
            continue
        growth = row["growth"]
        growth = 0.0 if (growth is None or growth != growth) else float(growth)
        out[str(gid)] = (growth, call_essential(growth, wt, frac))
    return out
