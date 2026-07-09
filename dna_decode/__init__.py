"""DNA → trait prediction platform.

Phase 1: E. coli antibiotic resistance prediction with biologically
interpretable attribution. See plans/Ecoli_G2P_Platform_Technical_Plan.md.
"""

# Derive from the installed package metadata so this can never drift from
# pyproject.toml (it was pinned to a stale "0.0.1" while pyproject was 0.6.5).
# Mirrors dna_decode.cli._version(). Falls back only for a source tree with no
# metadata (e.g. running uninstalled).
try:
    from importlib.metadata import PackageNotFoundError, version as _pkg_version

    try:
        __version__ = _pkg_version("dna_decode")
    except PackageNotFoundError:  # not installed (bare source checkout)
        __version__ = "0+unknown"
except ImportError:  # pragma: no cover - importlib.metadata always present on 3.8+
    __version__ = "0+unknown"
