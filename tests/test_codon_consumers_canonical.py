"""Step 6 guard — in-package codon consumers resolve to the ONE canonical table (2026-10-01).

The three in-package tables were verified byte-identical before consolidation, so this refactor must be
behaviourally inert. These tests pin both halves: one object (not three copies), and the same assignments
as before the change.
"""

from __future__ import annotations

from dna_decode.constraints import codon_tables as T
from dna_decode.forward import genome_edit as GE
from dna_decode.forward import inverse as INV
from dna_decode.typing import codon_map as CM

# The PRE-CHANGE assignments, pinned as literals so the refactor is provably inert. These are the
# values all three copies carried (NCBI table 1) before consolidation.
PRE_CHANGE_SAMPLE = {
    "TGA": "*", "TAA": "*", "TAG": "*", "ATG": "M", "TGG": "W",
    "AGA": "R", "AGG": "R", "ATA": "I", "TCG": "S", "TTG": "L", "GGG": "G",
}


def test_all_three_consumers_are_the_SAME_OBJECT_as_canonical():
    """One object, not three copies — that is what makes drift impossible rather than merely unlikely."""
    assert GE._CODON is T.STANDARD
    assert INV.CODON_TABLE is T.STANDARD
    assert CM.CODON is T.STANDARD


def test_the_refactor_is_behaviourally_inert():
    for codon, aa in PRE_CHANGE_SAMPLE.items():
        assert GE._CODON[codon] == aa, codon
    assert len(GE._CODON) == 64
    assert sorted(c for c, a in GE._CODON.items() if a == "*") == ["TAA", "TAG", "TGA"]


def test_no_in_package_module_still_declares_its_own_64_codon_literal():
    """Static check: the literals are gone, not merely shadowed."""
    import ast
    import pathlib
    for path in ("dna_decode/forward/genome_edit.py", "dna_decode/forward/inverse.py",
                 "dna_decode/typing/codon_map.py"):
        tree = ast.parse(pathlib.Path(path).read_text(encoding="utf-8"))
        for node in ast.walk(tree):
            if isinstance(node, ast.Dict) and len(node.keys) >= 60:
                raise AssertionError(f"{path} still declares a {len(node.keys)}-entry dict literal")


def test_each_module_imports_from_the_canonical_source():
    import ast
    import pathlib
    for path, alias in (("dna_decode/forward/genome_edit.py", "_CODON"),
                        ("dna_decode/forward/inverse.py", "CODON_TABLE"),
                        ("dna_decode/typing/codon_map.py", "CODON")):
        tree = ast.parse(pathlib.Path(path).read_text(encoding="utf-8"))
        found = [n for n in ast.walk(tree) if isinstance(n, ast.ImportFrom)
                 and n.module == "dna_decode.constraints.codon_tables"
                 and any(a.asname == alias and a.name == "STANDARD" for a in n.names)]
        assert found, f"{path} does not import STANDARD as {alias}"


def test_pointfinders_import_surface_is_unchanged():
    """`pointfinder/runner.py` imports `typing.codon_map`; Step 6 deliberately did not edit it, so the
    public names it relies on must still resolve."""
    from dna_decode.typing.codon_map import CODON, subject_aa_by_codon, translate
    assert callable(translate) and callable(subject_aa_by_codon)
    assert translate("ATGTGG") == "MW"
    assert CODON is T.STANDARD
    import ast
    import pathlib
    tree = ast.parse(pathlib.Path("dna_decode/pointfinder/runner.py").read_text(encoding="utf-8"))
    assert any(isinstance(n, ast.ImportFrom) and n.module == "dna_decode.typing.codon_map"
               for n in ast.walk(tree)), "pointfinder must still import typing.codon_map"


def test_the_clade_tables_are_reachable_from_the_same_module():
    """The point of consolidation: a consumer can now ask for a non-standard code."""
    assert T.table_for("bacteria.mycoplasma")["TGA"] == "W"
    assert GE._CODON["TGA"] == "*", "the default path must remain the standard code"
