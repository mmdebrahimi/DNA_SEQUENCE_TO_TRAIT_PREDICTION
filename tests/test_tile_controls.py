"""Guards for the non-tile control cohort loader.

The load-bearing test is `test_pos_control_at_150_RAISES_rather_than_returning_empty`: an empty list
and a refusal are the same value to a careless caller, and they mean opposite things. Everything else
pins the parse contract against a synthetic fixture plus the two measured real-data counts.
"""

from __future__ import annotations

import gzip
from pathlib import Path

import pytest

from dna_decode.glm.tile_controls import (
    EXPECTED_AT_150,
    ControlCohortError,
    control_counts,
    load_control_tiles,
)

# EIGHT columns, because the loader now routes through genomewide.load_peak_tiles, which
# enforces MIN_COLUMNS = 8 and requires variant / expn_med / peak_start / peak_end / strand /
# RNA_exp_sum_1 / RNA_exp_sum_2. The five-column fixture these tests used predates the routing
# and would fail at the owner's header check. Widening the FIXTURE is the fix; the assertions
# below are the contract and are unchanged.
HDR = ["variant", "expn_med", "category", "RNA_exp_sum_1", "RNA_exp_sum_2",
       "peak_start", "peak_end", "strand"]


def _write(tmp_path: Path, rows: list[list[str]]) -> Path:
    p = tmp_path / "GSE144621_peak_tile_expression_formatted_std.txt.gz"
    with gzip.open(p, "wt") as fh:
        fh.write("\t".join(HDR) + "\n")
        for r in rows:
            fh.write("\t".join(r) + "\n")
    return p


def _row(seq: str, cat: str, e: str = "0.5", peak: str = "100") -> list[str]:
    """One substrate row. `peak` varies the peak id, which the owner composes from
    peak_start/peak_end/strand; controls drop it, so it only has to parse."""
    return [seq, e, cat, "1.0", "2.0", peak, str(int(peak) + 150), "+"]


def test_filters_on_category_and_length(tmp_path):
    rows = [
        _row("A" * 150, "random"),
        _row("C" * 150, "neg_control"),
        _row("G" * 150, "tile"),          # the thing being controlled FOR
        _row("T" * 100, "random"),        # wrong geometry
    ]
    _write(tmp_path, rows)
    got = load_control_tiles(tmp_path)
    assert control_counts(got) == {"random": 1, "neg_control": 1}


def test_target_column_is_expn_med_and_reps_are_carried(tmp_path):
    _write(tmp_path, [["A" * 150, "7.25", "random", "3.5", "4.5", "100", "250", "+"]])
    # single-category request on purpose: the DEFAULT asks for neg_control too, and this fixture has
    # none, so the default would (correctly) refuse. That refusal is pinned in its own test below.
    (c,) = load_control_tiles(tmp_path, categories=("random",))
    assert c.expression == 7.25 and c.rep1 == 3.5 and c.rep2 == 4.5


def test_non_acgt_is_dropped_like_load_peak_tiles_does(tmp_path):
    _write(tmp_path, [_row("N" + "A" * 149, "random"), _row("A" * 150, "random")])
    assert len(load_control_tiles(tmp_path, categories=("random",))) == 1


def test_pos_control_at_150_RAISES_rather_than_returning_empty(tmp_path):
    """The geometry-incompatibility guard, and the reason this module exists in this shape.

    `pos_control` rows are 48-75 bp in the real substrate. Returning `[]` would read downstream as
    "the encoder could not separate these", a claim about the encoder, when the fact is that the
    category does not exist at this geometry. Padding is refused: it was measured to perturb this
    encoder's pooled embedding by max |delta| 1.232.
    """
    _write(tmp_path, [_row("A" * 60, "pos_control"), _row("A" * 150, "random")])
    with pytest.raises(ControlCohortError) as ei:
        load_control_tiles(tmp_path, categories=("random", "pos_control"))
    assert "pos_control" in str(ei.value)


def test_a_category_absent_entirely_also_raises(tmp_path):
    _write(tmp_path, [_row("A" * 150, "random")])
    with pytest.raises(ControlCohortError):
        load_control_tiles(tmp_path, categories=("random", "neg_control"))


def test_missing_required_column_raises(tmp_path):
    # Still deliberately malformed. It now fails at the OWNER's header check (MIN_COLUMNS = 8)
    # rather than at a local column lookup, and the error is translated to ControlCohortError.
    p = tmp_path / "GSE144621_peak_tile_expression_formatted_std.txt.gz"
    with gzip.open(p, "wt") as fh:
        fh.write("variant\tcategory\n" + "A" * 150 + "\trandom\n")
    with pytest.raises(ControlCohortError):
        load_control_tiles(tmp_path)


def test_allow_download_is_refused_not_honoured(tmp_path):
    _write(tmp_path, [_row("A" * 150, "random")])
    with pytest.raises(ControlCohortError):
        load_control_tiles(tmp_path, allow_download=True)


def test_absent_substrate_raises_rather_than_returning_empty(tmp_path):
    with pytest.raises(ControlCohortError):
        load_control_tiles(tmp_path / "nope")


def test_the_module_does_not_PARSE_the_substrate_itself():
    """The structural proof that the parser is gone, asserted at the SYNTAX level.

    A substring search for "gzip" was specifically REJECTED and must not be reinstated: the module
    docstring explains that the fragment files are space-delimited while the peak-tile file is
    tab-delimited, and the word survives there, so a text search would go red on correct code. That
    is the same weakness the ownership guard's own `\\bGSE\\d+` pattern has -- a text search over
    source cannot tell an executable reference from a comment. This parses instead, following the
    convention in tests/test_amr_uncounted.py.
    """
    import ast

    src = Path("dna_decode/glm/tile_controls.py").read_text(encoding="utf-8")
    tree = ast.parse(src)

    imported = set()
    for n in ast.walk(tree):
        if isinstance(n, ast.Import):
            imported.update(a.name.split(".")[0] for a in n.names)
        elif isinstance(n, ast.ImportFrom) and n.module:
            imported.add(n.module.split(".")[0])
    assert "gzip" not in imported, sorted(imported)

    attrs = {n.value.id for n in ast.walk(tree)
             if isinstance(n, ast.Attribute) and isinstance(n.value, ast.Name)}
    assert "gzip" not in attrs, sorted(attrs)

    # and it DOES call the owner -- the positive half, so the test cannot pass on an empty module
    assert "genomewide" in {m.split(".")[0] for m in
                            {n.module for n in ast.walk(tree)
                             if isinstance(n, ast.ImportFrom) and n.module}} or \
           any("genomewide" in (n.module or "") for n in ast.walk(tree)
               if isinstance(n, ast.ImportFrom))


def test_passing_a_FILE_path_now_raises_because_the_owner_resolves_the_filename(tmp_path):
    """Pins a deliberate signature narrowing rather than leaving it as undocumented drift.

    The previous version accepted either a directory or the file itself. The owner joins a cache
    DIRECTORY with the substrate's filename, so a file path is no longer meaningful.
    """
    p = _write(tmp_path, [_row("A" * 150, "random"), _row("C" * 150, "neg_control")])
    assert p.is_file()
    with pytest.raises(ControlCohortError) as ei:
        load_control_tiles(p)
    assert "DIRECTORY" in str(ei.value)


REAL = Path("data/mpra/GSE144621_peak_tile_expression_formatted_std.txt.gz")


@pytest.mark.skipif(not REAL.is_file(), reason="gitignored substrate not cached")
def test_real_substrate_control_counts_are_exactly_the_measured_ones():
    """936 random + 470 neg_control at 150 bp, measured 2026-11-01. These are the 1,406 sequences
    every arm on this rung has discarded via `category == 'tile'`."""
    got = control_counts(load_control_tiles(REAL.parent))
    assert got == EXPECTED_AT_150, got
    assert sum(got.values()) == 1406


# --------------------------------------------------------------------------------------------------
# gaps found by /test-epilogue after the owner-routing correction
# --------------------------------------------------------------------------------------------------
def test_the_refusal_is_about_GEOMETRY_not_about_the_CATEGORY_NAME(tmp_path):
    """Non-vacuity for the load-bearing refusal, and it is the half that was missing.

    `test_pos_control_at_150_RAISES_rather_than_returning_empty` proves the refusal fires. It does
    NOT prove the refusal is keyed on geometry: a hardcoded `if "pos_control" in categories: raise`
    would pass it identically. This is the positive half -- the SAME category at a geometry it does
    exist at must LOAD, which it can only do if `length` is a real filter rather than a name ban.

    It also pins the escape hatch the error message advertises ("request a geometry-compatible
    category"): a caller who narrows `length` instead gets data, not a refusal.
    """
    _write(tmp_path, [_row("A" * 60, "pos_control"), _row("C" * 150, "random")])
    (c,) = load_control_tiles(tmp_path, categories=("pos_control",), length=60)
    assert c.category == "pos_control" and len(c.seq) == 60
    # and the same request at the arm's native geometry still refuses
    with pytest.raises(ControlCohortError):
        load_control_tiles(tmp_path, categories=("pos_control",), length=150)


def test_the_length_filter_is_EXACT_not_a_minimum(tmp_path):
    """A `>=` instead of `==` would silently admit the fragment-length rows the arm's geometry
    excludes, and the real-data count assertion is the only thing that would ever notice."""
    _write(tmp_path, [_row("A" * 151, "random"), _row("C" * 149, "random"), _row("G" * 150, "random")])
    (c,) = load_control_tiles(tmp_path, categories=("random",))
    assert c.seq == "G" * 150


def test_a_str_cache_dir_is_accepted_like_a_Path(tmp_path):
    """The annotation is `Path | str`; only `Path` was exercised. `Path(cache_dir)` is the one line
    that makes the `str` half true, and a future refactor that forwards `cache_dir` straight to the
    owner would break it without any existing test noticing."""
    _write(tmp_path, [_row("A" * 150, "random")])
    assert len(load_control_tiles(str(tmp_path), categories=("random",))) == 1


def test_the_measured_counts_constant_COVERS_the_default_categories():
    """Offline guard for an assertion that is otherwise only checked when gitignored data is present.

    `test_real_substrate_control_counts_are_exactly_the_measured_ones` compares
    `control_counts(...)` to `EXPECTED_AT_150` by equality, so widening `CONTROL_CATEGORIES` without
    extending `EXPECTED_AT_150` makes that test fail -- but ONLY on a host that has the substrate
    cached, which is not the default. This catches the same drift with no data at all.
    """
    from dna_decode.glm.tile_controls import CONTROL_CATEGORIES

    assert set(CONTROL_CATEGORIES) == set(EXPECTED_AT_150), (
        sorted(CONTROL_CATEGORIES), sorted(EXPECTED_AT_150))
    assert sum(EXPECTED_AT_150.values()) == 1406
    # pos_control must never be a default: it is geometry-incompatible at 150 bp, so defaulting to
    # it would make the loader refuse every call.
    assert "pos_control" not in CONTROL_CATEGORIES
    assert "tile" not in CONTROL_CATEGORIES, "tile is the thing being controlled FOR"


def test_ControlSeq_carries_no_peak_field_so_nobody_can_split_controls_on_it():
    """Pins the stated reason `ControlSeq` is not the owner's `Tile`.

    The docstring's argument is that `Tile.peak` is the SPLIT UNIT for tiles and is meaningless for
    a control, so returning `Tile` would "invite someone to split controls on it". Nothing asserted
    that, and the cheapest future tidy-up -- returning the owner's type directly -- would read as a
    simplification while silently reintroducing the field.
    """
    from dna_decode.glm.genomewide import Tile
    from dna_decode.glm.tile_controls import ControlSeq

    assert ControlSeq is not Tile
    fields = set(ControlSeq.__dataclass_fields__)
    assert "peak" not in fields, fields
    assert fields == {"seq", "expression", "category", "rep1", "rep2"}, fields
    # frozen, so a caller cannot mutate a control row in place and have it diverge from its count
    assert ControlSeq.__dataclass_params__.frozen


def test_the_owners_error_is_TRANSLATED_but_its_message_and_cause_survive(tmp_path):
    """The translation comment promises a stable public failure type; this pins that it is a
    translation and not a replacement. A bare `raise ControlCohortError("substrate unavailable")`
    would satisfy `pytest.raises(ControlCohortError)` while discarding which file was missing --
    the only information a caller can act on."""
    from dna_decode.glm.genomewide import GenomeWideDataError

    with pytest.raises(ControlCohortError) as ei:
        load_control_tiles(tmp_path / "absent")
    assert isinstance(ei.value.__cause__, GenomeWideDataError)
    assert str(ei.value) == str(ei.value.__cause__)
    assert "GSE144621_peak_tile_expression_formatted_std.txt.gz" in str(ei.value)


def test_an_EMPTY_category_request_returns_EMPTY_the_one_hole_in_the_refusal(tmp_path):
    """RECORDED, not claimed closed -- the same shape as
    `test_the_REMAINING_hole_is_recorded_not_claimed_closed` in the ownership guard.

    The refusal is built from `per_cat`, which is keyed by the REQUESTED categories, so an empty
    request has nothing to find empty and the function returns `[]`. That is precisely the value the
    module's docstring says must never be returned, because downstream it reads as "the control found
    nothing" -- a claim about the encoder -- when the truth is that nothing was asked for.

    No call site passes an empty tuple today, so this is a latent hole rather than a live defect, and
    it is pinned here so a future session finds it as a known fact instead of rediscovering it from
    a wrong number. Closing it is a one-line guard on `categories` and would make this test fail
    loudly, which is the intended signal.
    """
    _write(tmp_path, [_row("A" * 150, "random")])
    assert load_control_tiles(tmp_path, categories=()) == []


def test_control_counts_omits_absent_categories_rather_than_reporting_zero():
    """`control_counts` builds its keys from the rows present, so a category with no survivors is
    ABSENT from the result, never `0`. That asymmetry against the loader's internal `per_cat` is why
    the real-data test can compare by equality; a version that seeded zeros would make
    `== EXPECTED_AT_150` fail for a reason unrelated to the counts."""
    from dna_decode.glm.tile_controls import ControlSeq

    rows = [ControlSeq(seq="A", expression=0.0, category="random") for _ in range(3)]
    rows.append(ControlSeq(seq="C", expression=0.0, category="neg_control"))
    assert control_counts(rows) == {"random": 3, "neg_control": 1}
    assert control_counts([]) == {}


def test_every_CALL_SITE_survives_the_directory_only_narrowing():
    """The narrowing pinned by `test_passing_a_FILE_path_now_raises...` is a BREAKING change, and
    nothing checked that the existing callers still satisfy it.

    All three consumers are scripts that need the real cached substrate, so they cannot be executed
    here; what is checkable offline is the shape of the argument they pass. A call site handing over
    a `.gz` filename now raises at run time, and the only symptom would be a script that refuses on
    a host where the substrate IS present -- which reads as missing data.
    """
    import ast

    sites = 0
    for rel in ("scripts/promoteratlas_control_separation.py",
                "scripts/promoteratlas_nonlinear_bound.py",
                "scripts/promoteratlas_subspace_removal.py"):
        src = Path(rel).read_text(encoding="utf-8")
        tree = ast.parse(src)
        # resolve the module-level constants each call site passes
        consts = {t.id: n.value for n in ast.walk(tree) if isinstance(n, ast.Assign)
                  for t in n.targets if isinstance(t, ast.Name)}
        for n in ast.walk(tree):
            if (isinstance(n, ast.Call) and isinstance(n.func, ast.Name)
                    and n.func.id == "load_control_tiles"):
                sites += 1
                assert n.args, f"{rel}: load_control_tiles called with no cache_dir"
                arg = n.args[0]
                assert isinstance(arg, ast.Name), f"{rel}: expected a named constant, got {arg!r}"
                src_expr = ast.unparse(consts[arg.id])
                assert ".gz" not in src_expr, f"{rel}: {arg.id} = {src_expr} is a FILE path"
                assert src_expr.endswith('"data/mpra"') or src_expr.endswith("'data/mpra'"), (
                    f"{rel}: {arg.id} = {src_expr}; expected the cache DIRECTORY")
    assert sites == 3, f"expected 3 call sites, found {sites}; the registry above is stale"
