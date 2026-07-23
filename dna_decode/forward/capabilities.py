"""Forward variant-effect capability preflight - what can this host actually RUN?

The strong forward methods (ESM2 / ProSST / GEMME / the hybrid) are validated but each needs a different
heavy dependency: torch+transformers for ESM2, +the ProSST library +a 3D structure for ProSST, Docker (or a
native JET2+R) +an MSA for GEMME. BLOSUM62 alone is wheel-only. This module is the honest "what's installed
-> which methods are runnable HERE" probe so `dna-decode forward` can (a) tell a user what they can run and
(b) degrade gracefully to the strongest available method instead of erroring.

Pure + cheap: it uses `importlib.util.find_spec` (does NOT import torch) + `shutil.which` (does NOT run
Docker), so calling it costs nothing and never triggers a model download.
"""
from __future__ import annotations

import importlib.util
import shutil

# The ordering IS the strength ranking (measured, wiki/mavedb_holdout_hybrid_2026-07-23.md):
# hybrid ~ prosst > esm2 > alphamissense >> blosum62. `auto` picks the strongest RUNNABLE one.
METHOD_STRENGTH = ["hybrid", "prosst", "esm2", "gemme", "alphamissense", "blosum62"]

# per-method install hint shown when a method is NOT runnable
_INSTALL_HINT = {
    "esm2": "pip install 'dna-decode[forward]'  (torch + transformers>=5)",
    "prosst": ("pip install 'dna-decode[forward]' + torch_geometric + torch_scatter + the AI4Protein/ProSST "
               "repo (structure quantizer; heavy/platform-specific)  (+ a structure: --uniprot / --pdb)"),
    "gemme": "install Docker Desktop (image elodielaine/gemme:gemme)  (+ an MSA: --msa)",
    "hybrid": "needs >=2 of {esm2, prosst, gemme} runnable",
    "alphamissense": "supply --am-table (precomputed AlphaMissense scores)",
}


def _has(mod: str) -> bool:
    try:
        return importlib.util.find_spec(mod) is not None
    except (ImportError, ValueError):
        return False


def probe_capabilities() -> dict:
    """Detect the heavy dependencies present on this host. Cheap: no imports of torch, no Docker calls."""
    return {
        "torch": _has("torch"),
        "transformers": _has("transformers"),
        "prosst": _has("prosst"),
        "docker": shutil.which("docker") is not None,
    }


def runnable_methods(caps: dict | None = None) -> dict:
    """Map each forward method -> (runnable: bool, reason: str), given the host capabilities.

    Structure/MSA are call-time inputs (a method can be *installed* yet need a --uniprot/--msa at run time);
    those are reported as a note, not a hard blocker here - the runner resolves them per call.
    """
    caps = caps if caps is not None else probe_capabilities()
    esm2_ok = caps["torch"] and caps["transformers"]
    prosst_ok = esm2_ok and caps["prosst"]
    gemme_ok = caps["docker"]
    n_hybrid = sum([esm2_ok, prosst_ok, gemme_ok])
    hybrid_ok = n_hybrid >= 2

    def _reason(ok: bool, method: str, extra: str = "") -> str:
        if ok:
            return "runnable" + (f" ({extra})" if extra else "")
        return f"not runnable - {_INSTALL_HINT.get(method, 'unavailable')}"

    return {
        "blosum62": (True, "runnable (wheel-only, deterministic)"),
        "esm2": (esm2_ok, _reason(esm2_ok, "esm2")),
        "prosst": (prosst_ok, _reason(prosst_ok, "prosst", "needs a structure at run time")),
        "gemme": (gemme_ok, _reason(gemme_ok, "gemme", "needs an MSA at run time")),
        "hybrid": (hybrid_ok, _reason(hybrid_ok, "hybrid",
                                      f"{n_hybrid} of esm2/prosst/gemme runnable")),
        "alphamissense": (False, _reason(False, "alphamissense")),
    }


def strongest_runnable(caps: dict | None = None) -> str:
    """The strongest method that is runnable on this host (always at least blosum62)."""
    rm = runnable_methods(caps)
    for m in METHOD_STRENGTH:
        if rm.get(m, (False, ""))[0]:
            return m
    return "blosum62"


def render_capabilities(caps: dict | None = None) -> str:
    """Human-readable preflight report: deps present + which methods are runnable + install hints."""
    caps = caps if caps is not None else probe_capabilities()
    rm = runnable_methods(caps)
    lines = ["forward variant-effect - capability preflight", ""]
    lines.append("  installed dependencies:")
    for dep in ("torch", "transformers", "prosst", "docker"):
        lines.append(f"    {'YES' if caps[dep] else 'no ':>4}  {dep}")
    lines.append("")
    lines.append("  runnable methods (strongest first):")
    for m in METHOD_STRENGTH:
        ok, reason = rm[m]
        mark = " OK " if ok else "  - "
        lines.append(f"    [{mark}] {m:<14} {reason}")
    lines.append("")
    lines.append(f"  --method auto (given the required run-time inputs) would use: {strongest_runnable(caps)}")
    lines.append("  strength ranking (measured, wiki/mavedb_holdout_hybrid_2026-07-23.md): "
                 "hybrid ~ prosst > esm2 > alphamissense >> blosum62")
    return "\n".join(lines)
