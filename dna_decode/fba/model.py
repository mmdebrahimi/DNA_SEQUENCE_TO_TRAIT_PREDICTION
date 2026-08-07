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

# BiGG model id -> the organism that model ACTUALLY reconstructs, verified against the BiGG Models
# API (http://bigg.ucsd.edu/api/v2/models/<id>) on 2026-08-07. This is the PROVENANCE SOURCE OF
# TRUTH: every emitted record stamps its organism from here, never from a hardcoded literal.
#
# History (see wiki/fba_wrong_organism_model_bug_2026-08-07.md): v0.11.0-v0.12.0 mapped
# `saureus` -> iYS1720 (actually a *Salmonella* pan-reactome; 1262/1707 gene ids carry the
# Salmonella Typhimurium `STM` prefix) and `paeruginosa` -> iJN1463 (actually *P. putida*).
# Both are corrected below; the models themselves are kept under their TRUE organism aliases.
MODEL_ORGANISM: dict[str, str] = {
    "iML1515": "Escherichia coli K-12 MG1655",
    "iYS854": "Staphylococcus aureus USA300_TCH1516",
    "iSB619": "Staphylococcus aureus N315",
    "iJN1463": "Pseudomonas putida KT2440",
    "iYS1720": "Salmonella (pan-reactome)",
    "iMM904": "Saccharomyces cerevisiae S288C",
}

# Cross-organism genome-scale models on BiGG (the engine is organism-agnostic; only E. coli is
# Keio-validated -- other organisms are v0 "engine generalizes" with their own essentiality gold
# standard deferred). organism alias -> BiGG model id.
_BIGG_MODELS: dict[str, str] = {
    "escherichia_coli": "iML1515", "ecoli": "iML1515", "e_coli": "iML1515", "escherichia": "iML1515",
    "staphylococcus_aureus": "iYS854", "saureus": "iYS854", "s_aureus": "iYS854",
    "salmonella": "iYS1720", "salmonella_enterica": "iYS1720",
    "pseudomonas_putida": "iJN1463", "pputida": "iJN1463", "p_putida": "iJN1463",
    "saccharomyces_cerevisiae": "iMM904", "yeast": "iMM904", "scerevisiae": "iMM904",
}
_DEFAULT_MODEL_ID = "iML1515"

# Organisms we are ASKED for but have no genome-scale model for in BiGG. Fail LOUDLY here rather
# than silently handing back a different species' model (the v0.11.0 defect). Queried the full
# BiGG model list 2026-08-07: zero hits for "aeruginosa".
_NO_BIGG_MODEL: dict[str, str] = {
    "pseudomonas_aeruginosa": "P. aeruginosa", "paeruginosa": "P. aeruginosa",
    "p_aeruginosa": "P. aeruginosa",
}


def organism_for(model_id: str) -> str:
    """BiGG model id -> the organism it actually reconstructs. Unknown ids are reported as such."""
    return MODEL_ORGANISM.get(model_id, f"unknown organism (BiGG model {model_id})")

# module-level per-model cache so repeated CLI calls in one process don't re-parse (~10s)
_MODELS: dict[str, object] = {}


def _cache_dir() -> Path:
    env = os.environ.get("DNA_DECODE_CACHE")
    base = Path(env) if env else (Path.home() / ".cache" / "dna_decode")
    return base / "fba"


def resolve_model_id(organism: str | None) -> str:
    """organism alias -> BiGG model id (default iML1515). Raises on an unknown organism."""
    if not organism:
        return _DEFAULT_MODEL_ID
    key = organism.strip().lower().replace(" ", "_").replace(".", "").replace("-", "_")
    if key in _BIGG_MODELS:
        return _BIGG_MODELS[key]
    if key in _NO_BIGG_MODEL:
        raise ValueError(
            f"no genome-scale metabolic model is available for {_NO_BIGG_MODEL[key]} in BiGG "
            f"(checked 2026-08-07). Refusing to substitute another organism's model. "
            f"Supply one explicitly with --model / --model-id if you have a reconstruction. "
            f"NOTE: dna-decode v0.11.0-v0.12.0 silently loaded iJN1463 (*Pseudomonas putida*) "
            f"for this alias -- any result from those versions is P. putida, not P. aeruginosa."
        )
    if organism in _BIGG_MODELS.values() or organism in MODEL_ORGANISM:  # a raw model id was passed
        return organism
    raise ValueError(
        f"unknown organism '{organism}'. Known: {sorted(set(_BIGG_MODELS))} "
        f"(or pass a BiGG model id via --model-id)."
    )


def resolve_model_path(model_id: str = _DEFAULT_MODEL_ID, download: bool = True) -> Path:
    """Locate <model_id>.xml.gz; optionally download from BiGG into the cache. Raises if unresolved."""
    if model_id == _DEFAULT_MODEL_ID:
        override = os.environ.get("DNA_DECODE_IML1515")
        if override and Path(override).exists():
            return Path(override)
        dev = Path(__file__).resolve().parent.parent.parent / "data" / "fba" / "iML1515.xml.gz"
        if dev.exists():
            return dev
    cached = _cache_dir() / f"{model_id}.xml.gz"
    if cached.exists():
        return cached
    if not download:
        raise FileNotFoundError(
            f"{model_id} model not found. Place it at data/fba/{model_id}.xml.gz or allow the BiGG download."
        )
    cached.parent.mkdir(parents=True, exist_ok=True)
    urllib.request.urlretrieve(  # free (BiGG Models)
        f"http://bigg.ucsd.edu/static/models/{model_id}.xml.gz", cached
    )
    return cached


def load_model(
    path: str | Path | None = None,
    *,
    organism: str | None = None,
    model_id: str | None = None,
    reuse: bool = True,
):
    """Load (and cache) a genome-scale cobra model. Default = iML1515 (E. coli). Lazy-imports cobra.

    `organism` (alias) or `model_id` (BiGG id) selects a cross-organism model; `path` overrides both.
    """
    mid = model_id or resolve_model_id(organism)
    if reuse and path is None and mid in _MODELS:
        return _MODELS[mid]
    try:
        import cobra  # noqa: PLC0415  (lazy heavy import)
    except ImportError as e:  # pragma: no cover - env-dependent
        raise ImportError(
            "cobrapy is required for the FBA cell: `pip install 'dna-decode[fba]'` or `uv pip install cobra`."
        ) from e
    p = Path(path) if path else resolve_model_path(mid)
    model = cobra.io.read_sbml_model(str(p))
    if path is None:
        _MODELS[mid] = model
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


def synthetic_lethality(model, gene_a: str, gene_b: str, frac: float = 0.01) -> dict:
    """Two-gene edit -> is the PAIR synthetic-lethal? (neither alone lethal, but the double is).

    Synthetic lethality = KO(a) viable AND KO(b) viable AND KO(a,b) lethal. It is how metabolic
    drug-target PAIRS are found (each single is buffered by an isozyme / alternate route; the double
    breaks it). Returns a `fba-synthetic-lethality-v1` record.
    """
    wt = wildtype_growth(model)
    ga = knockout_growth(model, gene_a)
    gb = knockout_growth(model, gene_b)
    gab = knockout_growth(model, [gene_a, gene_b])
    a_ess = call_essential(ga, wt, frac)
    b_ess = call_essential(gb, wt, frac)
    ab_ess = call_essential(gab, wt, frac)
    is_sl = (not a_ess) and (not b_ess) and ab_ess
    return {
        "record": "fba-synthetic-lethality-v1",
        "gene_a": gene_a,
        "gene_b": gene_b,
        "wildtype_growth_per_h": round(wt, 4),
        "ko_a_growth_per_h": round(ga, 4),
        "ko_b_growth_per_h": round(gb, 4),
        "double_ko_growth_per_h": round(gab, 4),
        "single_a_essential": a_ess,
        "single_b_essential": b_ess,
        "double_essential": ab_ess,
        "synthetic_lethal": is_sl,
        "verdict": (
            "SYNTHETIC-LETHAL (each single viable; the PAIR is lethal)" if is_sl
            else "not synthetic-lethal"
            + (" (a single is already essential)" if (a_ess or b_ess) else " (double is viable)")
        ),
    }


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
