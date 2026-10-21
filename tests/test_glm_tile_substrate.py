"""The tile chain's two substrate dependencies, routed through `external_root`.

All offline and all hermetic: an "unreachable" directory is produced by monkeypatching
`genomewide._reachable`, never by naming a drive letter. That seam exists because of a real defect one
module over -- `test_describe_names_the_missing_requirement` once asserted wording against whatever
volumes happened to be attached, so it passed on a host carrying a name-matching cache and failed the
moment that volume detached. A unit test must not be a statement about the operator's hardware.

THREE LOAD-BEARING TESTS, and the non-vacuity ones are the point:

* `test_a_REACHABLE_cache_does_not_trip_the_guard` -- the guard's whole safety argument is that the
  `configured`-status path is bit-for-bit unchanged, so on the host that HAS the drive nothing moved.
  A guard that refused everything would satisfy the refusal test and break the fragment-rung arms.
* `test_allow_discovered_is_never_passed_True_anywhere_in_the_tile_chain` -- an AST pin, because the
  defect this guards is a later convenience flip, not today's code. Auto-adopting a discovered root
  converts a loud missing-ROOT into a quiet missing-DATASET.
* `test_the_in_dir_check_comes_AFTER_the_tile_load` -- also AST, and it protects a MEASUREMENT rather
  than a behaviour. The two dependencies have different wall classes; getting past the tile load and
  failing only on promcalc is what tells them apart, and a well-meaning "validate all inputs up front"
  refactor would silently destroy that diagnosis while every test still passed.
"""
from __future__ import annotations

import ast
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from dna_decode.glm import genomewide as G  # noqa: E402

REPO = Path(__file__).resolve().parents[1]
PRED_SRC = REPO / "scripts" / "glm_tile_predictions.py"
GENOMEWIDE_SRC = REPO / "dna_decode" / "glm" / "genomewide.py"


# ------------------------------------------------------------------ the reachability guard
def test_an_unreachable_cache_raises_a_DIAGNOSABLE_error_not_a_PermissionError(monkeypatch, tmp_path):
    """Measured 2026-10-21: this used to be an unhandled `WinError 21` from the `mkdir` below."""
    monkeypatch.setattr(G, "_reachable", lambda d: False)
    with pytest.raises(G.GenomeWideDataError) as ei:
        G.ensure_cached(tmp_path / "nowhere" / "mpra", which="peak_tile", allow_download=False)
    msg = str(ei.value)
    assert "not attached" in msg
    assert "ABSENT" in msg
    assert "nowhere" in msg.replace(chr(92), "/")


def test_the_error_names_the_candidate_subdirectory_to_pass(monkeypatch, tmp_path):
    """A candidate is a ROOT that HOLDS `mpra`; passing the root itself would be wrong."""
    monkeypatch.setattr(G, "_reachable", lambda d: False)
    with pytest.raises(G.GenomeWideDataError) as ei:
        G.ensure_cached(tmp_path / "mpra", which="peak_tile", allow_download=False)
    assert "<candidate>/mpra" in str(ei.value)


def test_the_error_says_retrying_will_not_help(monkeypatch, tmp_path):
    """An absent VOLUME is not a download problem, and the two have different remedies."""
    monkeypatch.setattr(G, "_reachable", lambda d: False)
    with pytest.raises(G.GenomeWideDataError) as ei:
        G.ensure_cached(tmp_path / "mpra", which="peak_tile", allow_download=False)
    assert "NOT a download problem" in str(ei.value)


def test_a_REACHABLE_cache_does_not_trip_the_guard(tmp_path):
    """NON-VACUITY. A guard that refused everything would pass the refusal tests and break the host
    that has the drive -- `load_fragments` and `load_peak_tiles` both come through here."""
    with pytest.raises(G.GenomeWideDataError) as ei:
        G.ensure_cached(tmp_path, which="peak_tile", allow_download=False)
    msg = str(ei.value)
    assert "missing and downloading is disabled" in msg      # the PRE-EXISTING error
    assert "not attached" not in msg                          # the new guard did NOT fire


def test_a_present_file_is_returned_unchanged(tmp_path):
    """The fully-happy path: reachable directory, file present, no download, no guard."""
    (tmp_path / G.PEAK_TILE_FILE).write_bytes(b"never opened by ensure_cached")
    assert G.ensure_cached(tmp_path, which="peak_tile", allow_download=False) == \
        tmp_path / G.PEAK_TILE_FILE


def test_a_missing_leaf_under_a_reachable_parent_is_still_created(tmp_path):
    """The parent exists, the leaf does not -- `mkdir` could always create this, so it must still."""
    leaf = tmp_path / "mpra"
    with pytest.raises(G.GenomeWideDataError, match="downloading is disabled"):
        G.ensure_cached(leaf, which="peak_tile", allow_download=False)
    assert leaf.is_dir()


def test_reachable_walks_up_to_an_existing_ancestor(tmp_path):
    assert G._reachable(tmp_path) is True
    assert G._reachable(tmp_path / "a" / "b" / "c") is True


def test_reachable_is_false_when_no_ancestor_exists(monkeypatch, tmp_path):
    """Hermetic stand-in for a drive with no media: nothing on the chain resolves."""
    monkeypatch.setattr(Path, "exists", lambda self: False)
    assert G._reachable(tmp_path / "a") is False


def test_reachable_treats_an_OSError_as_unreachable(monkeypatch, tmp_path):
    """Fail-safe direction: an undecidable stat refuses rather than proceeding to mkdir."""
    def _boom(self):
        raise OSError("device not ready")
    monkeypatch.setattr(Path, "exists", _boom)
    assert G._reachable(tmp_path) is False


def test_the_substrate_report_is_fail_soft(monkeypatch, tmp_path):
    """A drive scan that raises must not replace a data error with an import error."""
    import dna_decode.data.external_root as ER
    monkeypatch.setattr(ER, "describe", lambda *a, **k: (_ for _ in ()).throw(RuntimeError("scan")))
    out = G._substrate_report(tmp_path / "mpra")
    assert "ABSENT" in out and "probe unavailable" in out


def test_an_unknown_file_key_still_fails_before_any_path_work():
    with pytest.raises(G.GenomeWideDataError, match="unknown file key"):
        G.ensure_cached("D:/whatever", which="not_a_key", allow_download=False)


def test_the_guard_also_covers_the_FRAGMENT_keys(monkeypatch, tmp_path):
    """`load_fragments` comes through the same `ensure_cached`, so scoping the guard to `peak_tile`
    would leave the fragment arms with the raw WinError this was written to replace."""
    monkeypatch.setattr(G, "_reachable", lambda d: False)
    for key in sorted(G.FRAG_FILES):
        with pytest.raises(G.GenomeWideDataError, match="not attached"):
            G.ensure_cached(tmp_path / "mpra", which=key, allow_download=False)


def test_the_substrate_report_asks_for_the_mpra_SUBDIRECTORY(monkeypatch, tmp_path):
    """`required=("mpra",)` is what makes the report DIAGNOSTIC rather than decorative.

    Dropped, every discovered candidate renders as a bare "name match" and the line that actually
    tells the operator what to do -- `MISSING : mpra  <- do not point the project here` -- never
    prints. The report would then name a same-named root with the wrong CONTENTS as if it were
    usable, which is the silent-substitution failure `external_root` exists to prevent.
    """
    seen: dict = {}
    import dna_decode.data.external_root as ER

    def _spy(configured, **kw):
        seen["configured"] = configured
        seen["kw"] = kw
        return "REPORT-BODY"

    monkeypatch.setattr(ER, "describe", _spy)
    assert G._substrate_report(tmp_path / "mpra") == "REPORT-BODY"
    assert seen["kw"].get("required") == ("mpra",)
    assert Path(seen["configured"]) == tmp_path / "mpra"


def test_the_substrate_report_is_embedded_verbatim_in_the_guard_error(monkeypatch, tmp_path):
    """Success branch. Only the fail-soft branch was pinned, so a report that silently returned
    nothing useful -- or was built but never interpolated -- would have gone unnoticed."""
    import dna_decode.data.external_root as ER
    monkeypatch.setattr(ER, "describe", lambda *a, **k: "SENTINEL-REPORT-LINE")
    monkeypatch.setattr(G, "_reachable", lambda d: False)
    with pytest.raises(G.GenomeWideDataError) as ei:
        G.ensure_cached(tmp_path / "mpra", which="peak_tile", allow_download=False)
    assert "SENTINEL-REPORT-LINE" in str(ei.value)


# ------------------------------------------------------------------ the script's refusals
def test_require_external_dir_returns_a_configured_path(tmp_path):
    import glm_tile_predictions as pred
    assert pred.require_external_dir(str(tmp_path), required=(), flag="--in-dir") == tmp_path


def test_require_external_dir_refuses_an_absent_path_with_the_report(tmp_path):
    import glm_tile_predictions as pred
    with pytest.raises(SystemExit) as ei:
        pred.require_external_dir(str(tmp_path / "gone"), required=("promcalc",), flag="--in-dir")
    msg = str(ei.value)
    assert "--in-dir" in msg
    assert "ABSENT" in msg
    assert "NOT substituted automatically" in msg


def test_require_external_dir_names_the_flag_it_was_given(tmp_path):
    """The refusal must tell the operator WHICH argument to pass, not just that something is missing."""
    import glm_tile_predictions as pred
    with pytest.raises(SystemExit) as ei:
        pred.require_external_dir(str(tmp_path / "gone"), required=(), flag="--mpra-cache")
    assert "--mpra-cache" in str(ei.value)


def test_the_script_exposes_an_mpra_cache_flag():
    import glm_tile_predictions as pred
    assert "--mpra-cache" in Path(pred.__file__).read_text(encoding="utf-8")


def test_require_external_dir_refuses_EVERY_non_configured_status(monkeypatch, tmp_path):
    """Hermetic, and it covers the statuses this host does not happen to produce.

    `resolve_external_root` can return `absent`, `absent-candidates-unfit` or `discovered`; which one
    the live probe yields depends on which volumes are attached, so an assertion driven by the real
    resolver tests the operator's hardware. `configured` is the ONLY status that may pass -- a
    `discovered` one especially must not, since the whole point of `allow_discovered=False` is that a
    right NAME is not evidence of right CONTENTS.
    """
    import glm_tile_predictions as pred
    import dna_decode.data.external_root as ER
    monkeypatch.setattr(ER, "describe", lambda *a, **k: "REPORT")

    for status in ("absent", "absent-candidates-unfit", "discovered"):
        monkeypatch.setattr(ER, "resolve_external_root",
                            lambda c, **k: (Path(c), status))
        with pytest.raises(SystemExit) as ei:
            pred.require_external_dir(str(tmp_path), required=("promcalc",), flag="--in-dir")
        assert status in str(ei.value), f"the refusal must name the status it saw ({status})"


def test_require_external_dir_refuses_a_FILE_at_the_configured_path(tmp_path):
    """A directory flag handed a file resolves to a non-`configured` status, so it must refuse rather
    than return a path that every downstream `iterdir`/glob would fail on."""
    import glm_tile_predictions as pred
    f = tmp_path / "promcalc"
    f.write_text("", encoding="utf-8")
    with pytest.raises(SystemExit):
        pred.require_external_dir(str(f), required=(), flag="--in-dir")


def test_the_mpra_cache_flag_is_REGISTERED_not_merely_mentioned():
    """The substring test above passes on a comment, and the script has one at line ~339.

    A flag that appears only in prose parses as an unrecognised argument, so pin the `add_argument`
    registration -- and `default=None`, which is what makes `a.mpra_cache or DEFAULT_CACHE` fall back
    to the published path instead of to an empty string.
    """
    adds = [c for c in _calls(PRED_SRC, "add_argument")
            if c.args and isinstance(c.args[0], ast.Constant) and c.args[0].value == "--mpra-cache"]
    assert len(adds) == 1, "--mpra-cache is not registered exactly once with argparse"
    default = next((k.value for k in adds[0].keywords if k.arg == "default"), None)
    assert isinstance(default, ast.Constant) and default.value is None


def test_the_tile_load_consumes_the_mpra_cache_FLAG_VALUE():
    """A flag that parses and is then ignored is worse than no flag: it reads as configurable while
    the loader keeps taking the hardcoded `DEFAULT_CACHE`. The sibling AST pin only requires SOME
    argument, so it would pass on `load_peak_tiles(cache_dir=DEFAULT_CACHE)`."""
    tree = ast.parse(PRED_SRC.read_text(encoding="utf-8"))
    from_flag = {
        t.id
        for n in ast.walk(tree) if isinstance(n, ast.Assign)
        for t in n.targets if isinstance(t, ast.Name)
        if any(isinstance(s, ast.Attribute) and s.attr == "mpra_cache"
               for s in ast.walk(n.value))
    }
    assert from_flag, "nothing is assigned from a.mpra_cache"
    passed = set()
    for c in _calls(PRED_SRC, "load_peak_tiles"):
        for kw in c.keywords:
            if kw.arg == "cache_dir" and isinstance(kw.value, ast.Name):
                passed.add(kw.value.id)
    assert passed & from_flag, f"load_peak_tiles(cache_dir=) gets {passed}, none from the flag"


def test_the_presented_refusal_suppresses_the_exception_CHAIN():
    """`from None` is deliberate. Without it the operator gets "During handling of the above
    exception" plus the loader's traceback -- which is exactly the "this tool is broken" reading the
    typed error was introduced to replace, restored one layer up."""
    tree = ast.parse(PRED_SRC.read_text(encoding="utf-8"))
    handlers = [h for n in ast.walk(tree) if isinstance(n, ast.Try) for h in n.handlers
                if isinstance(h.type, ast.Name) and h.type.id == "GenomeWideDataError"]
    assert handlers, "the GenomeWideDataError handler is gone"
    raises = [n for h in handlers for n in ast.walk(h) if isinstance(n, ast.Raise)]
    assert raises, "the handler no longer raises"
    for r in raises:
        assert isinstance(r.cause, ast.Constant) and r.cause.value is None, \
            "raise ... from None was dropped; the traceback chain is back"


# ------------------------------------------------------------------ AST pins
def _calls(src: Path, func: str) -> list[ast.Call]:
    tree = ast.parse(src.read_text(encoding="utf-8"))
    out = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Call):
            name = node.func.attr if isinstance(node.func, ast.Attribute) else \
                getattr(node.func, "id", None)
            if name == func:
                out.append(node)
    return out


def test_load_peak_tiles_is_never_called_with_no_arguments_in_the_tile_chain():
    """The single line that caused the crash: `load_peak_tiles()` bare, taking the hardcoded D: path."""
    calls = _calls(PRED_SRC, "load_peak_tiles")
    assert calls, "expected at least one call; did the script move?"
    for c in calls:
        assert c.args or c.keywords, "load_peak_tiles() with no arguments takes the hardcoded D: path"


def test_allow_discovered_is_never_passed_True_anywhere_in_the_tile_chain():
    """Pinned on SOURCE SHAPE, because the defect being guarded is a FUTURE convenience flip.

    Auto-adopting a discovered root is the silent-substitution failure `external_root` exists to
    prevent -- a same-named root with different contents turns a loud missing-ROOT into a quiet
    missing-DATASET.
    """
    for src in (PRED_SRC, GENOMEWIDE_SRC):
        for c in _calls(src, "resolve_external_root"):
            for kw in c.keywords:
                if kw.arg == "allow_discovered":
                    assert not (isinstance(kw.value, ast.Constant) and kw.value.value is True), \
                        f"{src.name} passes allow_discovered=True"


def test_the_in_dir_check_comes_AFTER_the_tile_load():
    """Protects a MEASUREMENT, not a behaviour.

    The two substrate dependencies have DIFFERENT wall classes. Getting past the tile load and failing
    only on promcalc is the measurement that separates them; hoisting both checks to argument-parsing
    time would collapse that into one undifferentiated "substrate missing" while every other test here
    still passed.
    """
    tile = _calls(PRED_SRC, "load_peak_tiles")
    indir = [c for c in _calls(PRED_SRC, "require_external_dir")
             if any(isinstance(k.value, ast.Constant) and k.value.value == "--in-dir"
                    for k in c.keywords)]
    assert tile and indir
    assert min(c.lineno for c in tile) < min(c.lineno for c in indir)


def test_the_reachability_check_precedes_the_mkdir():
    """Ordering IS the fix: `mkdir` on an absent volume is what raised the raw WinError."""
    tree = ast.parse(GENOMEWIDE_SRC.read_text(encoding="utf-8"))
    fn = next(n for n in ast.walk(tree)
              if isinstance(n, ast.FunctionDef) and n.name == "ensure_cached")
    reach = [n.lineno for n in ast.walk(fn)
             if isinstance(n, ast.Call) and getattr(n.func, "id", None) == "_reachable"]
    mkdir = [n.lineno for n in ast.walk(fn)
             if isinstance(n, ast.Call) and getattr(n.func, "attr", None) == "mkdir"]
    assert reach and mkdir and min(reach) < min(mkdir)


def test_the_loader_raises_typed_and_the_script_presents(monkeypatch, tmp_path):
    """Layering pin: a library that called SystemExit would be unusable outside a CLI, and a script
    that let the exception escape shows a traceback -- which reads as "this tool is broken"."""
    src = PRED_SRC.read_text(encoding="utf-8")
    assert "except GenomeWideDataError" in src
    gsrc = GENOMEWIDE_SRC.read_text(encoding="utf-8")
    fn_start = gsrc.index("def ensure_cached")
    assert "SystemExit" not in gsrc[fn_start:fn_start + 3000]
