"""Regression guard: no test module may hard-import a heavy optional dependency at module scope.

WHY THIS EXISTS. CI installs only `.[dev]` (the deterministic decoder + dev tools) — NOT torch / the
`[ml]` extra. A bare module-level `import torch` (or xgboost, transformers, …) therefore fails at
COLLECTION with ModuleNotFoundError, which pytest treats as a hard error → `Interrupted: N errors during
collection` → the WHOLE job fails (exit 2), regardless of how many real tests would pass. This has broken
CI — and flooded the maintainer's inbox with failure emails — twice: xgboost in the predict-e2e files, and
torch in `test_genome_jepa_prototype.py`.

The fix each time is the pytest-canonical guard: `pytest.importorskip("<dep>")` BEFORE the import, so the
module SKIPS cleanly on a runner without the dep. This meta-test enforces that guard for every test file so
the failure mode cannot recur silently.

Pure + offline: a line scan of the test files' own source. No heavy import is performed here.
"""
from __future__ import annotations

import re
from pathlib import Path

import pytest

# Optional deps NOT in the CI `.[dev]` install. A module-scope import of any of these must be guarded.
_HEAVY = ("torch", "xgboost", "transformers", "sentence_transformers", "bitsandbytes")
# line-start (module scope, not indented) `import X` / `import X.y` / `from X import …`
_MODULE_IMPORT = re.compile(
    r"^(?:import|from)\s+(" + "|".join(_HEAVY) + r")(?:[.\s]|$)", re.MULTILINE)

_TESTS_DIR = Path(__file__).resolve().parent


def _test_files() -> list[Path]:
    return sorted(p for p in _TESTS_DIR.glob("test_*.py") if p.name != Path(__file__).name)


def test_no_unguarded_heavy_import_at_module_scope():
    offenders: list[str] = []
    for f in _test_files():
        src = f.read_text(encoding="utf-8")
        hits = {m.group(1) for m in _MODULE_IMPORT.finditer(src)}
        if hits and "importorskip" not in src:
            offenders.append(f"{f.name}: module-scope import of {sorted(hits)} without pytest.importorskip")
    assert not offenders, (
        "These test modules hard-import a heavy optional dep at module scope; CI (`.[dev]`, no torch/ml) "
        "will FAIL AT COLLECTION on them. Add `pytest.importorskip(\"<dep>\")` before the import:\n  "
        + "\n  ".join(offenders))


def test_the_guard_regex_actually_matches_a_bare_import():
    """Pin the matcher so a future refactor can't silently neuter the guard."""
    assert _MODULE_IMPORT.search("import torch\n")
    assert _MODULE_IMPORT.search("from torch import nn\n")
    assert _MODULE_IMPORT.search("import xgboost as xgb\n")
    # indented (inside a function/fixture) is fine — not module scope
    assert not _MODULE_IMPORT.search("    import torch\n")
    # a substring in an unrelated name must not false-positive
    assert not _MODULE_IMPORT.search("import torchvision_helper_that_is_not_torch\n".replace("torch", "xtorch"))


def test_known_ml_test_files_carry_the_guard():
    """A positive check: the files that legitimately use torch/xgboost do guard it."""
    for name in ("test_genome_jepa_prototype.py",):
        p = _TESTS_DIR / name
        if p.exists():
            assert "importorskip" in p.read_text(encoding="utf-8"), f"{name} lost its importorskip guard"


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(pytest.main([__file__, "-q"]))
