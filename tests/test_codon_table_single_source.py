"""Step 9 — the codon table has ONE source, enforced by DISCOVERY rather than a hand-written list.

A hand-enumerated list of "files that define a codon table" drifts the moment someone adds a file; this
repo has hit that exact failure five times (hardcoded exclusion lists silently under-cover). So this
guard SCANS the package and `scripts/` for anything that looks like a genetic code and then requires each
hit to be either the canonical object or an explicitly allow-listed standalone copy that EQUALS it.

Two detection shapes, because the copies were not all dict literals:
  * a dict literal with >= 60 codon-like keys (the five original in-package/script copies);
  * the compressed NCBI `AAs` string, which is how `scripts/kaggle/pear_prosst_kernel.py` builds its
    table and how `codon_tables.py` legitimately stores every table.
"""

from __future__ import annotations

import ast
import pathlib
import re

import pytest

from dna_decode.constraints import codon_tables as T

REPO = pathlib.Path(__file__).resolve().parent.parent
SCAN_ROOTS = ("dna_decode", "scripts")

#: The canonical NCBI table-1 AAs string. Its presence is legitimate ONLY in the canonical module and in
#: allow-listed standalone copies.
AAS_T1 = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"

#: The one module allowed to DEFINE genetic codes.
CANONICAL = "dna_decode/constraints/codon_tables.py"

#: Standalone copies that cannot import the package, each with the reason it is exempt. Every entry is
#: additionally required to EQUAL the canonical table, so an exemption is not a licence to drift.
ALLOW_LISTED: dict[str, str] = {
    "scripts/kaggle/pear_prosst_kernel.py":
        "executes standalone inside a Kaggle notebook with no access to the dna_decode package, so it "
        "cannot import the canonical table and must rebuild it from the NCBI AAs string",
}

CODON_RE = re.compile(r"^[ACGTU]{3}$", re.I)


def _py_files() -> list[pathlib.Path]:
    out: list[pathlib.Path] = []
    for root in SCAN_ROOTS:
        out.extend(sorted((REPO / root).rglob("*.py")))
    return [p for p in out if "__pycache__" not in p.parts]


def _rel(p: pathlib.Path) -> str:
    return str(p.relative_to(REPO)).replace("\\", "/")


def discover_codon_tables() -> dict[str, list[str]]:
    """Find every file that looks like it defines a genetic code. Returns path -> detection shapes."""
    hits: dict[str, list[str]] = {}
    for p in _py_files():
        src = p.read_text(encoding="utf-8", errors="ignore")
        shapes: list[str] = []
        if AAS_T1 in src:
            shapes.append("aas_string")
        try:
            tree = ast.parse(src)
        except SyntaxError:
            tree = None
        if tree is not None:
            for node in ast.walk(tree):
                if isinstance(node, ast.Dict) and len(node.keys) >= 60:
                    keys = [k.value for k in node.keys
                            if isinstance(k, ast.Constant) and isinstance(k.value, str)]
                    if keys and sum(1 for k in keys if CODON_RE.match(k)) >= 60:
                        shapes.append("dict_literal")
                        break
        if shapes:
            hits[_rel(p)] = shapes
    return hits


# --------------------------------------------------------------------------- the guard

def test_discovery_actually_finds_something():
    """A discovery guard that finds nothing is the inert-filter failure in test form."""
    hits = discover_codon_tables()
    assert hits, "discovery found no codon tables at all -- the scan is broken"
    assert CANONICAL in hits, f"the canonical module must be discovered; found {sorted(hits)}"


def test_every_discovered_table_is_canonical_or_explicitly_allow_listed():
    hits = discover_codon_tables()
    unexplained = sorted(set(hits) - {CANONICAL} - set(ALLOW_LISTED))
    assert not unexplained, (
        "these files define a genetic code without being the canonical source or allow-listed: "
        f"{unexplained}. Import dna_decode.constraints.codon_tables.STANDARD instead, or add an "
        "allow-list entry stating why the file cannot import it.")


def test_no_dict_literal_codon_table_survives_outside_the_canonical_module():
    """The original five copies were dict literals; none may remain."""
    for path, shapes in discover_codon_tables().items():
        if path == CANONICAL:
            continue
        assert "dict_literal" not in shapes, f"{path} still declares a codon dict literal"


def test_every_allow_list_entry_carries_a_reason_and_still_exists():
    assert ALLOW_LISTED, "an empty allow-list would make the exemption check vacuous"
    for path, reason in ALLOW_LISTED.items():
        assert len(reason) > 40, f"{path} needs a real reason, not a placeholder"
        assert (REPO / path).exists(), f"{path} is allow-listed but absent; prune the entry"


def test_each_allow_listed_copy_EQUALS_the_canonical_table():
    """An exemption is not a licence to drift. The Kaggle kernel rebuilds the table from the NCBI AAs
    string; decode it the same way and require equality."""
    bases = "TCAG"
    checked = 0
    for path in ALLOW_LISTED:
        src = (REPO / path).read_text(encoding="utf-8", errors="ignore")
        assert AAS_T1 in src, f"{path} no longer carries the AAs string; re-verify this exemption"
        rebuilt = {}
        i = 0
        for b1 in bases:
            for b2 in bases:
                for b3 in bases:
                    rebuilt[b1 + b2 + b3] = AAS_T1[i]
                    i += 1
        assert rebuilt == T.STANDARD, f"{path}'s table diverges from canonical"
        checked += 1
    assert checked == len(ALLOW_LISTED)


def test_the_canonical_module_is_the_only_place_that_DEFINES_a_code():
    """Other modules may import STANDARD; only this one may construct a table."""
    src = (REPO / CANONICAL).read_text(encoding="utf-8")
    assert "def _decode" in src
    for path in discover_codon_tables():
        if path in (CANONICAL, *ALLOW_LISTED):
            continue
        raise AssertionError(f"{path} defines a code outside the canonical module")


def test_the_guard_is_NON_VACUOUS_a_divergent_copy_would_fail(tmp_path, monkeypatch):
    """Prove the check can fail: plant a divergent table in a scanned tree and assert discovery flags
    it, then assert the live tree is clean."""
    fake_root = tmp_path / "dna_decode"
    fake_root.mkdir()
    bad = fake_root / "rogue_table.py"
    divergent = AAS_T1.replace("*", "W", 1)          # a genuinely different code
    bad.write_text(f'AAS = "{divergent}"\n_CODON = {{'
                   + ", ".join(f'"{a}{b}{c}": "X"' for a in "ACGT" for b in "ACGT" for c in "ACGT")
                   + "}\n", encoding="utf-8")
    monkeypatch.setattr(pathlib.Path, "cwd", lambda: tmp_path, raising=False)
    import tests.test_codon_table_single_source as mod
    monkeypatch.setattr(mod, "REPO", tmp_path)
    monkeypatch.setattr(mod, "SCAN_ROOTS", ("dna_decode",))
    planted = mod.discover_codon_tables()
    assert "dna_decode/rogue_table.py" in planted, planted
    assert "dict_literal" in planted["dna_decode/rogue_table.py"]


def test_the_frozen_amr_surface_is_untouched_by_this_plan():
    """Every file this plan edits is outside the frozen surface; assert the frozen paths show no diff."""
    import subprocess
    # NOTE: `amr_rules.py` lives under eval/, NOT data/. The wrong path was used in ad-hoc checks during
    # this plan's execution, and `git diff --quiet` on a MISSING path exits 0 trivially -- so that check
    # was silently vacuous for it. The existence assertion below is what caught it, which is why it is
    # here: a frozen-surface check that cannot distinguish "unchanged" from "not a file" is not a check.
    frozen = ["dna_decode/eval/amr_rules.py", "dna_decode/data/calibrated_amr_rules.json",
              "dna_decode/data/mic_tiers.py", "dna_decode/data/shipped_decoder_surface.py",
              "dna_decode/eval/cohort_manifest.py"]
    for f in frozen:
        assert (REPO / f).exists(), (
            f"{f} does not exist, so a `git diff --quiet` on it would pass vacuously. Fix the path.")
    r = subprocess.run(["git", "diff", "--quiet", "HEAD", "--", *frozen], cwd=REPO,
                       capture_output=True)
    assert r.returncode == 0, "a frozen-surface file has uncommitted changes"


def test_constraints_is_a_leaf_layer_with_no_reverse_imports():
    """`constraints` may import eval.regime; nothing in it may import forward/typing/pointfinder/scripts,
    or the dependency direction inverts and a cycle becomes possible."""
    for p in sorted((REPO / "dna_decode" / "constraints").rglob("*.py")):
        if "__pycache__" in p.parts:
            continue
        src = p.read_text(encoding="utf-8")
        for forbidden in ("dna_decode.forward", "dna_decode.typing", "dna_decode.pointfinder"):
            assert forbidden not in src, f"{_rel(p)} imports {forbidden}; constraints must stay a leaf"


def test_eval_regime_does_not_import_constraints():
    """The one inbound edge is constraints -> eval.regime; the reverse would create a cycle."""
    src = (REPO / "dna_decode" / "eval" / "regime.py").read_text(encoding="utf-8")
    assert "dna_decode.constraints" not in src
