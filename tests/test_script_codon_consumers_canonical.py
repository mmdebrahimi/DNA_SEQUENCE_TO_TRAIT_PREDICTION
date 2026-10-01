"""Step 7 guard — script-level codon consumers resolve to the canonical table (2026-10-01).

Verified STATICALLY (via `ast`) rather than by importing, because `scripts/forward_inverse_roundtrip.py`
has import-time side effects that make in-process exec fail. That failure is PRE-EXISTING: the HEAD
version fails identically to the post-edit version, so it is a property of the script and not of this
change. A static check is also the stronger guard here — it proves the literal is gone rather than
merely shadowed at runtime.
"""

from __future__ import annotations

import ast
import pathlib

import pytest

from dna_decode.constraints import codon_tables as T

MIGRATED = {
    "scripts/fungal_erg11_caller.py": "_CODON",
    "scripts/forward_inverse_roundtrip.py": "CODON_TABLE",
    "scripts/build_holdout_hybrid_manifest.py": "_CODON",
}

#: Executes standalone on Kaggle with no access to the package, so it CANNOT import and keeps a private
#: table. Step 9's discovery guard pins it to equality with the canonical table.
STANDALONE_EXEMPT = {"scripts/kaggle/pear_prosst_kernel.py"}


@pytest.mark.parametrize("path,alias", sorted(MIGRATED.items()))
def test_each_script_imports_the_canonical_table_under_its_original_alias(path, alias):
    tree = ast.parse(pathlib.Path(path).read_text(encoding="utf-8"))
    found = [n for n in ast.walk(tree) if isinstance(n, ast.ImportFrom)
             and n.module == "dna_decode.constraints.codon_tables"
             and any(a.name == "STANDARD" and a.asname == alias for a in n.names)]
    assert found, f"{path} must import STANDARD as {alias} (alias preserved, call sites unchanged)"


@pytest.mark.parametrize("path", sorted(MIGRATED))
def test_no_migrated_script_still_declares_a_codon_literal(path):
    """Neither a 64-entry dict nor the compressed AAs-string construction may remain."""
    src = pathlib.Path(path).read_text(encoding="utf-8")
    tree = ast.parse(src)
    for node in ast.walk(tree):
        if isinstance(node, ast.Dict) and len(node.keys) >= 60:
            raise AssertionError(f"{path} still declares a {len(node.keys)}-entry dict literal")
    assert "FFLLSSSSYY**CC*WLLLL" not in src, f"{path} still builds a table from a compressed AAs string"


def test_the_loadable_scripts_resolve_to_the_canonical_object():
    """Two of the three import cleanly in-process; assert object identity for those."""
    import importlib.util
    checked = 0
    for path, alias in MIGRATED.items():
        if "forward_inverse_roundtrip" in path:
            continue                                  # pre-existing import-time side effects
        spec = importlib.util.spec_from_file_location(pathlib.Path(path).stem + "_probe", path)
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        assert getattr(mod, alias) is T.STANDARD, path
        checked += 1
    assert checked == 2, "expected two in-process-loadable migrated scripts"


def test_the_high_blast_radius_caller_is_an_IMPORT_ONLY_change():
    """`fungal_erg11_caller`'s codon mapper is reused by the HIV RT/PR/IN and SARS-CoV-2 Mpro callers,
    so its call sites must be untouched — only the table's provenance changed."""
    src = pathlib.Path("scripts/fungal_erg11_caller.py").read_text(encoding="utf-8")
    assert "_CODON.get(" in src, "the original call sites must remain"
    assert "HIGH BLAST RADIUS" in src, "the risk must be documented at the edit site"


def test_the_standalone_kaggle_kernel_is_exempt_and_the_reason_is_recorded():
    """It cannot import the package, so it is allow-listed here and equality-pinned in Step 9."""
    for path in STANDALONE_EXEMPT:
        p = pathlib.Path(path)
        if not p.exists():
            pytest.skip(f"{path} absent")
        src = p.read_text(encoding="utf-8")
        assert "dna_decode.constraints" not in src, (
            f"{path} now imports the package; remove it from STANDALONE_EXEMPT")
