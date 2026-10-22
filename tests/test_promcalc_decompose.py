"""Guards for stage B of the promcalc chain (`scripts/promcalc_decompose.py`).

**Why this file exists.** Stage B's module docstring asserts that it imports NOTHING from
`dna_decode` and tells the reader to *keep it that way* -- and until 2026-10-22 nothing enforced it
(`grep -rln promcalc_decompose tests/` returned no matches). A stated invariant with no mechanism is
the `--help` over-claim pattern this repo has repeatedly had to retract.

The boundary is load-bearing rather than decorative: stage B runs under a SEPARATE interpreter because
`promotercalculator` resolves sklearn to 1.9.1 in a fresh environment while the project is pinned at
1.8.0, and `uv.lock` is a frozen reproducibility unit so the dependency cannot be added to it. If
stage B ever imported `dna_decode`, that separation would silently stop being possible.

These guards are AST-only and deliberately so. They pin the import boundary, the laziness of the
third-party import and the two 2026-10-22 version breaks; they do NOT exercise the decomposition
itself, which needs the real package and is covered at the real surface by the bounded `--smoke` run
instead. Stated here so the coverage is not over-read.
"""
from __future__ import annotations

import ast
import importlib.util
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
STAGE_B = REPO / "scripts" / "promcalc_decompose.py"


def _tree() -> ast.Module:
    return ast.parse(STAGE_B.read_text(encoding="utf-8"))


def _load_stage_b(tag: str):
    """Import stage B under a private name.

    This is itself a real check: the project environment has no `promotercalculator`, so this only
    succeeds while the third-party import stays lazy.
    """
    spec = importlib.util.spec_from_file_location(f"_promcalc_stage_b_{tag}", STAGE_B)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _module_level_imports(tree: ast.Module) -> list[str]:
    """Only the imports at module scope -- the ones that run on a plain `import`."""
    out: list[str] = []
    for node in tree.body:
        if isinstance(node, ast.Import):
            out += [a.name for a in node.names]
        elif isinstance(node, ast.ImportFrom):
            out.append(node.module or "")
    return out


def _all_imports(tree: ast.Module) -> list[str]:
    """Every import at ANY nesting depth, so a function-local import cannot hide."""
    out: list[str] = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            out += [a.name for a in node.names]
        elif isinstance(node, ast.ImportFrom):
            out.append(node.module or "")
    return out


# --- the stated-but-unenforced import boundary --------------------------------------------------

def test_module_level_imports_are_STDLIB_ONLY():
    mods = _module_level_imports(_tree())
    offenders = sorted({m for m in mods if m.split(".")[0] not in sys.stdlib_module_names})
    assert not offenders, (
        f"stage B gained a non-stdlib module-level import: {offenders}. It must be importable in "
        f"ANY interpreter, including one holding no project dependencies, or the separate-environment "
        f"boundary stops working."
    )


def test_the_module_level_import_set_is_NOT_EMPTY():
    """Guard the guard: an empty set would make the stdlib-only assertion vacuous.

    A parse failure, or a restructure that moved every import inside a function, would otherwise read
    as "no third-party imports" -- a green test over a file it never actually inspected.
    """
    assert _module_level_imports(_tree()), "no module-level imports found; the stdlib guard is vacuous"


def test_dna_decode_is_never_imported_AT_ANY_DEPTH():
    """The docstring's stated contract, now enforced.

    Checked by AST rather than by grepping for the string, because this module legitimately MENTIONS
    `dna_decode` twice -- once in prose and once inside a default path -- and neither is an import.
    A regex cannot tell those apart, and that ambiguity is exactly why this is an AST test.
    """
    offenders = sorted({m for m in _all_imports(_tree()) if m.split(".")[0] == "dna_decode"})
    assert not offenders, (
        f"stage B imports {offenders}. Its docstring states it imports NOTHING from dna_decode, and "
        f"the separate-interpreter boundary depends on that staying true."
    )


def _lazy_imports_of(tree: ast.Module, prefix: str) -> int:
    """How many imports of `prefix` sit inside a function (i.e. are lazy)."""
    n = 0
    for fn in ast.walk(tree):
        if not isinstance(fn, (ast.FunctionDef, ast.AsyncFunctionDef)):
            continue
        for node in ast.walk(fn):
            if isinstance(node, ast.ImportFrom) and (node.module or "").split(".")[0] == prefix:
                n += 1
            elif isinstance(node, ast.Import):
                n += sum(1 for a in node.names if a.name.split(".")[0] == prefix)
    return n


def test_promoter_calculator_is_imported_LAZILY_inside_a_function():
    """The third-party import must not run at module-import time.

    That laziness is what lets the PROJECT environment import stage B to test it at all -- the project
    env does not have `promotercalculator` installed and must not need it.
    """
    tree = _tree()
    assert _lazy_imports_of(tree, "promoter_calculator"), (
        "no function-local `promoter_calculator` import found; if it moved to module level, stage B "
        "can no longer be imported in an environment without the package installed."
    )
    at_module_level = [m for m in _module_level_imports(tree)
                       if m.split(".")[0] == "promoter_calculator"]
    assert not at_module_level, (
        f"`promoter_calculator` is imported at MODULE level ({at_module_level}); it must stay lazy."
    )


def test_stage_b_imports_cleanly_in_THIS_interpreter():
    """The real consequence of the guards above, exercised rather than asserted."""
    mod = _load_stage_b("import")
    assert hasattr(mod, "decompose_one")
    assert hasattr(mod, "SITE_FIELDS")


def test_the_install_name_differs_from_the_import_name_and_that_is_RECORDED():
    """`pip install promotercalculator`, `import promoter_calculator`.

    A naive install of the import name 404s on PyPI. The trap is recorded next to the code so the next
    reader provisioning the interpreter does not lose time to it.
    """
    src = STAGE_B.read_text(encoding="utf-8")
    assert "promotercalculator" in src, (
        "the DISTRIBUTION name `promotercalculator` is not mentioned in stage B; without it a reader "
        "sees only the import name, and installing that 404s."
    )
    assert "promoter_calculator" in src


# --- the two 2026-10-22 version breaks ----------------------------------------------------------

def test_the_class_import_TOLERATES_BOTH_package_layouts():
    """`promotercalculator` 1.2.4 stopped re-exporting the CLASS at the package root.

    Its `__init__.py` exports only the lowercase `promoter_calculator` FUNCTION from `wrapper`, so the
    recorded `from promoter_calculator import Promoter_Calculator` raises ImportError on EVERY tile.
    Measured: 40 errored tiles and zero rows in 24 s, caught by the bounded smoke gate rather than
    after ~5 h of compute.

    The class itself is UNCHANGED -- 1.2.4's own `wrapper.py` constructs it with the same signature and
    reads the same `output()['Forward_Predictions_per_TSS']` key -- so this is an import-PATH fix, not
    a model swap, which is what lets stage C's frozen reconcile anchors stay the right acceptance test.
    """
    src = STAGE_B.read_text(encoding="utf-8")
    assert "_load_calculator" in src, "the version-tolerant class resolver is gone"
    assert "promoter_calculator.promoter_calculator import Promoter_Calculator" in src, (
        "the >=1.2.4 submodule fallback is missing; stage B will fail on every tile against a current "
        "install of the package."
    )

    fn = next((n for n in _tree().body
               if isinstance(n, ast.FunctionDef) and n.name == "_load_calculator"), None)
    assert fn is not None, "_load_calculator is not a module-level function"

    # The recorded path must be attempted FIRST, or an older install silently changes which code
    # object it runs.
    order = [(node.module or "") for node in ast.walk(fn) if isinstance(node, ast.ImportFrom)]
    assert order, "_load_calculator contains no imports"
    assert order[0] == "promoter_calculator", (
        f"the first attempted import is {order[0]!r}; the RECORDED top-level path must be tried first "
        f"so an older install keeps its original import."
    )


def test_per_site_field_access_does_NOT_use_dot_get():
    """Each per-TSS value is a `PromoCalcResults` object in 1.2.4, not a dict.

    It defines `__getitem__`, `__setitem__` and `keys()` but NOT `get` -- so the original
    `d.get(field, "nan")` raises AttributeError per site, and because `decompose_one` catches broadly
    the whole tile is then counted as FAILED. Subscripting works on both shapes, which is why
    `_site_field` subscripts and catches instead.
    """
    src = STAGE_B.read_text(encoding="utf-8")
    assert "_site_field" in src, "the shape-tolerant field accessor is gone"

    fn = next((n for n in ast.walk(_tree())
               if isinstance(n, ast.FunctionDef) and n.name == "decompose_one"), None)
    assert fn is not None

    bad = [
        node for node in ast.walk(fn)
        if isinstance(node, ast.Call)
        and isinstance(node.func, ast.Attribute)
        and node.func.attr == "get"
        and isinstance(node.func.value, ast.Name)
        and node.func.value.id == "d"
    ]
    assert not bad, (
        f"{len(bad)} `d.get(...)` call(s) remain in decompose_one. PromoCalcResults has no `.get`, so "
        f"each one raises AttributeError and silently converts every site into a failed tile."
    )


def test_site_field_handles_both_a_dict_and_an_object_without_get():
    """`_site_field` exercised against BOTH real shapes, not merely asserted about."""
    site_field = _load_stage_b("sitefield")._site_field

    # shape 1: the recorded dict
    assert site_field({"dG_total": -1.5}, "dG_total") == -1.5
    assert site_field({}, "dG_total", "nan") == "nan"

    # shape 2: the 1.2.4 object -- subscriptable, and DELIBERATELY without `.get`
    class PromoCalcResultsLike:
        __annotations__ = {"dG_total": float}

        def __init__(self, **kw):
            self.__dict__.update(kw)

        def __getitem__(self, key):
            if key in self.__annotations__:
                return getattr(self, key)
            raise KeyError(key)

        def keys(self):
            return self.__annotations__.keys()

    obj = PromoCalcResultsLike(dG_total=-2.25)
    assert not hasattr(obj, "get"), "the fixture must reproduce the missing-`.get` shape"
    assert site_field(obj, "dG_total") == -2.25
    assert site_field(obj, "nope", "nan") == "nan"


def test_int_start_parses_the_LIST_form_that_1_2_4_returns():
    """`hex35_position` comes back as a real list `[25, 31]`, not the string the docstring describes.

    `_int_start` stringifies before parsing, so it handles both -- measured on real output, 0 of 2,920
    rows failed to parse. Pinned because a tidy-up that stopped stringifying would break the live shape.
    """
    mod = _load_stage_b("intstart")
    assert mod._int_start([25, 31]) == 25        # the 1.2.4 shape
    assert mod._int_start("[25, 31]") == 25      # the shape the docstring records
    assert mod._int_start(None) == -1            # unparseable -> sentinel, never a crash
    assert mod._int_start("") == -1


def test_shard_paths_never_mixes_smoke_output_into_a_full_run():
    """A bare `decomposition*` glob also matches `decomposition_smoke.csv.gz`.

    That silently merged a 20-tile smoke output into a full substrate once, and was harmless only
    because the smoke tiles were a prefix of the real ones. The 2026-10-22 run wrote a 40-tile smoke
    into the SAME directory as the full output, so the separation is live, not hypothetical.
    """
    src = STAGE_B.read_text(encoding="utf-8")
    assert 'glob("decomposition_smoke*.csv.gz")' in src
    assert 'glob("decomposition.part-*.csv.gz")' in src
    assert 'glob("decomposition*' not in src, (
        "a bare `decomposition*` glob is back; it matches decomposition_smoke*.csv.gz and will merge "
        "smoke rows into a full substrate."
    )


def test_SITE_FIELDS_matches_the_header_the_real_run_wrote():
    """The committed tuple against the header measured from the 2026-10-22 smoke output.

    A mismatch means the installed package returns different fields -- i.e. version drift of a kind
    that would invalidate stage C's reconcile anchors. Pinned as data so the two cannot drift apart.
    """
    measured = (
        "tile_id", "tss",
        "dG_total", "dG_10", "dG_35", "dG_ext10", "dG_spacer", "dG_disc", "dG_ITR", "dG_UP", "dG_bind",
        "Tx_rate", "hex35_start", "hex10_start", "spacer_len",
    )
    assert tuple(_load_stage_b("fields").SITE_FIELDS) == measured
