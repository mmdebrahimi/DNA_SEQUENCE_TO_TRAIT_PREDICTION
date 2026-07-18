"""Evolution modality for the forward hybrid — GEMME (a deterministic evolutionary-conservation model).

GEMME (Laine et al. 2019) scores a variant by how well it fits the evolutionary history in an MSA — a global
epistatic model with ZERO learned parameters. On ProteinGym it is a top-tier evolution predictor (~0.455)
and, combined with sequence (ESM2) + structure (ProSST), forms the sweep-top 3-way
(`wiki/forward_modality_hybrid_2026-07-17.md`: `ESM2+GEMME+ProSST` 0.547).

This module is the GEMME component of the hybrid, mirroring `prosst_scorer` / `structure_scorer`:
- `gemme_table_from_column(rows, col="GEMME")` — the DEPLOYABLE-NOW path: adapt a precomputed GEMME column
  (e.g. ProteinGym's `pg_zeroshot`) into a `{mutation: score}` table, oriented higher=preserved, ready for
  `variant_effect.rank_average_hybrid`. Because GEMME is DETERMINISTIC (0 params), a precomputed column IS
  canonical GEMME output — using it is the same move as using ProteinGym's pre-quantized ProSST structures.
- `run_gemme(msa_path, query_seq)` — the NOVEL-protein path: run GEMME locally. GEMME needs JET2 (Java) + R +
  an MSA; it is Windows-hostile (the same class as ESM-IF/ProSST's `torch_geometric` quantizer). So it
  LAZY-checks the toolchain and raises `GemmeUnavailable` when absent; run it on a Linux host (feed the MSA
  from `msa_fetch.fetch_msa`). Deferred — the seam is complete, the real run needs the toolchain.

Sign: GEMME's native score is a NEGATIVE conservation delta (more negative = more deleterious). ProteinGym's
column is already oriented higher=fitter, so the adapter passes it through; a raw GEMME run is oriented in
`run_gemme`.
"""
from __future__ import annotations

import shutil
from pathlib import Path


class GemmeUnavailable(RuntimeError):
    """Raised when the GEMME toolchain (JET2/Java + R) is absent (this host)."""


# GEMME score tiers on the ProteinGym-oriented scale (higher = preserved); coarse, like the other modalities.
_GEMME_PRESERVED = -0.5
_GEMME_DAMAGING = -4.0


def gemme_table_from_column(rows, col: str = "GEMME", mutant_key: str = "mutant") -> dict[str, float]:
    """Adapt a precomputed GEMME column into {mutation: score} (higher=preserved). `rows` is an iterable of
    dict rows (e.g. csv.DictReader over ProteinGym's pg_zeroshot). Skips rows with a missing/NA GEMME cell.

    GEMME is deterministic (0 params) -> a precomputed column is canonical GEMME output, so this is the
    deployable evolution component for any protein ProteinGym (or a GEMME run) has already scored."""
    out: dict[str, float] = {}
    for r in rows:
        m = r.get(mutant_key)
        v = r.get(col)
        if not m or v in (None, "", "NA"):
            continue
        try:
            out[m] = float(v)
        except (TypeError, ValueError):
            continue
    if not out:
        raise ValueError(f"no usable {col!r} values in the supplied rows")
    return out


def find_jet2() -> str | None:
    """Resolve the JET2 launcher (GEMME's structure/MSA engine) — env override -> PATH."""
    import os
    env = os.environ.get("JET2_BIN")
    if env and Path(env).exists():
        return env
    return shutil.which("JET2") or shutil.which("jet2")


def run_gemme(msa_path: str | Path, query_seq: str, *, work_dir: str | Path | None = None) -> dict[str, float]:
    """Run GEMME on a real MSA -> {mutation: score} for every single substitution (higher=preserved).

    NOVEL-protein path. Needs the GEMME toolchain (JET2/Java + R) — Windows-hostile; run on a Linux host,
    feeding the MSA from `msa_fetch.fetch_msa`. Raises GemmeUnavailable when the toolchain is absent. The
    subprocess wiring is finalized against the first real Linux run (mirrors ESM-IF/ProSST's deferral)."""
    if find_jet2() is None or shutil.which("R") is None:
        raise GemmeUnavailable(
            "GEMME toolchain not found — needs JET2 (Java) + R (Windows-hostile). Set $JET2_BIN and install "
            "R, or run on a Linux host. For VALIDATION use gemme_table_from_column with a precomputed column.")
    raise GemmeUnavailable(  # pragma: no cover - real subprocess wiring lands on the first Linux run
        "GEMME subprocess wiring is deferred to the first real Linux run; the toolchain resolved but the "
        "runner is not yet finalized. Use gemme_table_from_column for the precomputed/validation path.")


def gemme_tier(score: float) -> str:
    """GEMME (ProteinGym-oriented, higher=preserved) -> forward-cell tier."""
    if score >= _GEMME_PRESERVED:
        return "preserved"
    if score <= _GEMME_DAMAGING:
        return "damaging"
    return "uncertain"
