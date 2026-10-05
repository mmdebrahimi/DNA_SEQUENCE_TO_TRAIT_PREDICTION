"""Every disclosure layer `trust_block` attaches must be able to REACH a human, not just the JSON.

THE BUG CLASS. `trust_block` attaches five namespace-separate layers (`DISCLOSURE_LAYERS`). Each exists
because some caveat is decision-relevant at call time. A layer that is attached but never rendered is a
silent half-disclosure: the machine-readable record is honest and the human output is not. The CLI's own
comment names the standard -- "a block carried only in JSON is not a disclosure" -- and two layers had
already violated it: `organism_scope` (fixed 2026-09-05) and `prospective` (fixed here).

`prospective` is the costly one to have missed. It is the project's STRONGEST tier -- an isolate
post-dating the lock is leakage-free BY CONSTRUCTION -- and both E. coli cells currently sit at
`superseded_by_surface_change` after the v2 gentamicin lock revised `amr_rules.py`. Their post-lock
numbers describe the RETIRED rule and are withheld. A reader who remembers "this cell has prospective
validation" was, at the CLI, relying on evidence for a rule that no longer ships.
"""
from __future__ import annotations

import io
import sys
import tempfile
from contextlib import redirect_stdout
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from dna_decode.amr.cli import main as amr_main  # noqa: E402
from dna_decode.data.trust_surface import (DISCLOSURE_LAYERS, prospective_one_line,  # noqa: E402
                                           species_composition_one_line, trust_block)

_HEADER = ("Protein id\tContig id\tStart\tStop\tStrand\tElement symbol\tElement name\tScope\tType\t"
           "Subtype\tClass\tSubclass\tMethod\tTarget length\tReference sequence length\t"
           "% Coverage of reference\t% Identity to reference\tAlignment length\t"
           "Closest reference accession\tClosest reference name\tHMM accession\tHMM description")


def _run_dir(tmp: Path, rows) -> Path:
    cells_rows = []
    for sym, cls, sub, meth in rows:
        cells = [""] * 22
        cells[5] = sym; cells[10] = cls; cells[11] = sub; cells[12] = meth
        cells_rows.append("\t".join(cells))
    d = tmp / "run"; d.mkdir()
    (d / "main.tsv").write_text("\n".join([_HEADER, *cells_rows]) + "\n", encoding="utf-8")
    return d


def _invoke(argv):
    buf = io.StringIO()
    with redirect_stdout(buf):
        rc = amr_main(argv)
    return rc, buf.getvalue()


_RMT_ROW = [("rmtB", "AMINOGLYCOSIDE", "AMINOGLYCOSIDE", "EXACTX")]


# --- the real surface -----------------------------------------------------------------------------

def test_prospective_status_reaches_the_printed_output():
    with tempfile.TemporaryDirectory() as td:
        rd = _run_dir(Path(td), _RMT_ROW)
        rc, out = _invoke(["--drug", "gentamicin", "--amrfinder-run", str(rd),
                           "--organism", "Escherichia_coli_Shigella"])
    assert rc == 0
    assert "prospective-lock" in out, f"prospective layer never printed:\n{out}"
    assert "SUPERSEDED" in out


def test_the_superseded_line_says_the_numbers_are_withheld_and_the_clock_restarts():
    """Naming the state without its consequence would leave a reader thinking the evidence still counts."""
    with tempfile.TemporaryDirectory() as td:
        rd = _run_dir(Path(td), _RMT_ROW)
        _rc, out = _invoke(["--drug", "gentamicin", "--amrfinder-run", str(rd),
                            "--organism", "Escherichia_coli_Shigella"])
    assert "WITHHELD" in out
    assert "clock restarts" in out
    assert "amr_rules.py" in out          # names WHICH file drifted


def test_it_is_silent_for_a_cell_with_no_prospective_evidence():
    """FIXTURE DERIVED, not hardcoded (fixed 2026-10-04).

    This test named Klebsiella x meropenem because that cell had no prospective artifact when it was
    written -- then one was committed on 2026-10-03 (`prospective_lock_validation_Klebsiella_meropenem_
    2026-10-03.json`, 16 post-lock isolates) and the test began failing, asserting silence about a cell
    that had since acquired exactly the evidence it was checking for the absence of. A hardcoded
    negative fixture goes stale the moment the project does the work.

    So the cell is now SELECTED by checking which cells actually lack an artifact. If a future accrual
    covers this one too, the fixture moves by itself; if every cell gains one, the test SKIPS rather
    than passing vacuously -- silence with nothing to be silent about proves nothing.
    """
    import glob
    covered = {Path(p).name.split("prospective_lock_validation_")[1].rsplit("_", 1)[0]
               for p in glob.glob("wiki/prospective_lock_validation_*.json")}
    candidates = [("meropenem", "Klebsiella_pneumoniae", "Klebsiella_meropenem",
                   ("blaKPC-2", "BETA-LACTAM", "CARBAPENEM", "EXACTX")),
                  ("tetracycline", "Klebsiella_pneumoniae", "Klebsiella_tetracycline",
                   ("tet(A)", "TETRACYCLINE", "TETRACYCLINE", "EXACTX")),
                  ("ceftriaxone", "Klebsiella_pneumoniae", "Klebsiella_ceftriaxone",
                   ("blaCTX-M-15", "BETA-LACTAM", "CEPHALOSPORIN", "EXACTX"))]
    pick = next((c for c in candidates if c[2] not in covered), None)
    if pick is None:
        pytest.skip("every candidate cell now has a prospective artifact -- nothing to assert silence on")
    drug, organism, _key, row = pick
    with tempfile.TemporaryDirectory() as td:
        rd = _run_dir(Path(td), [row])
        _rc, out = _invoke(["--drug", drug, "--amrfinder-run", str(rd), "--organism", organism])
    assert "prospective-lock" not in out, f"{organism} x {drug} has no artifact yet printed a block"


def test_the_call_is_unchanged_by_the_disclosure():
    """L2 qualifies a call; it must never alter one."""
    with tempfile.TemporaryDirectory() as td:
        rd = _run_dir(Path(td), _RMT_ROW)
        _rc, out = _invoke(["--drug", "gentamicin", "--amrfinder-run", str(rd),
                            "--organism", "Escherichia_coli_Shigella", "--json-only"])
    import json
    assert json.loads(out)["prediction"] == "R"


# --- the renderer's states ------------------------------------------------------------------------

def test_scored_reports_numbers_and_names_the_construction():
    line = prospective_one_line({"prospective": {
        "status": "scored", "n_scored": 61, "acc": 0.967, "sens": 0.917, "spec": 1.0,
        "powering": "POWERED"}})
    assert "61" in line and "0.967" in line
    assert "BY CONSTRUCTION" in line


def test_an_underpowered_score_is_flagged_not_quoted_flat():
    line = prospective_one_line({"prospective": {
        "status": "scored", "n_scored": 8, "acc": 0.5, "powering": "UNDERPOWERED"}})
    assert "UNDERPOWERED" in line and "indicative" in line


def test_lock_unverified_withholds_the_numbers():
    line = prospective_one_line({"prospective": {"status": "lock_unverified",
                                                 "generated": "2026-08-24"}})
    assert "WITHHELD" in line


def test_not_accrued_is_silent_and_silence_never_means_validated():
    """An armed lock with nothing accrued has no call-time consequence; printing it every call is noise."""
    assert prospective_one_line({"prospective": {"status": "not_accrued"}}) is None
    assert prospective_one_line({}) is None
    assert prospective_one_line({"prospective": None}) is None


# --- the structural guard -------------------------------------------------------------------------

_RENDERED = {"lineage", "source_concentration", "prospective", "organism_scope", "error_rates",
             # species_composition added 2026-10-05. It has a renderer
             # (`species_composition_one_line`) because the reachability guard in
             # test_evidence_surface_reachable.py required one: the plan that added it to the
             # report card left "should it reach a CALL?" open, and that guard had already
             # answered it -- a layer nobody can see from the tool is not a disclosure.
             "species_composition"}


def test_the_known_layer_set_has_not_grown_unnoticed():
    """A NEW layer added to DISCLOSURE_LAYERS without a renderer is the exact bug this file exists for.

    `doubt_layer` is deliberately NOT in `_RENDERED`: see the tripwire below.
    """
    assert set(DISCLOSURE_LAYERS) == _RENDERED | {"doubt_layer"}, (
        "DISCLOSURE_LAYERS changed -- wire a renderer for the new layer into dna_decode/amr/cli.py "
        "and add it to _RENDERED, or justify it here the way doubt_layer is justified.")


def test_error_rates_renders_and_leads_with_the_dangerous_direction():
    """Added 2026-09-10, and this file's tripwire is what required the renderer rather than a card-only
    block. VME must come FIRST in the line: it is the direction that leaves an infection untreated, and
    burying it behind ME would reproduce the ordering problem the layer exists to fix."""
    from dna_decode.data.trust_surface import error_rates_one_line
    line = error_rates_one_line({"error_rates": {
        "vme": 0.533, "vme_n_resistant": 30, "vme_ci": [0.361, 0.698],
        "me": 0.1, "me_n_susceptible": 30, "me_ci": [0.035, 0.256]}})
    assert line and line.index("VME") < line.index("ME " if "ME " in line else "ME")
    assert "0.533" in line and "30 resistant" in line
    assert "reports susceptible" in line, "the line must say what the number MEANS, not just name it"


def test_error_rates_is_silent_without_a_confusion_matrix():
    """A cell with no counts has no rate. Printing a bare label would imply a measurement."""
    from dna_decode.data.trust_surface import error_rates_one_line
    assert error_rates_one_line({}) is None
    assert error_rates_one_line({"error_rates": None}) is None
    assert error_rates_one_line({"error_rates": {"vme": None, "me": 0.1}}) is None


def test_doubt_layer_is_inert_today_and_this_trips_the_moment_it_is_not():
    """TRIPWIRE, not a renderer.

    `doubt_layer` (the per-cell determinant-COMPLETENESS screen -- distinct from the record's top-level
    per-call `doubt`) is attached to the badge and has no renderer. Writing one now would be dead code:
    measured across every deployed cell, NO cell carries a strong completeness signal, so the line could
    never fire. The honest move is to leave it unwritten and make its absence DETECTABLE the moment it
    starts mattering, rather than ship speculative output.

    If this fails, a cell has acquired a strong completeness signal that a CLI reader cannot see. Write
    the renderer then.
    """
    offenders = []
    for drug in ("ciprofloxacin", "gentamicin", "ceftriaxone", "tetracycline", "meropenem"):
        for org in ("Escherichia_coli_Shigella", "Klebsiella_pneumoniae", "Salmonella_enterica"):
            dl = (trust_block(drug, org) or {}).get("doubt_layer") or {}
            if dl.get("n_strong_completeness_signals"):
                offenders.append((drug, org, dl.get("strong_families")))
    assert not offenders, (
        f"doubt_layer now carries a strong completeness signal {offenders} but has NO renderer, so it "
        "reaches the JSON and never a human. Write a doubt_layer one-liner and wire it into the CLI.")


def test_EVERY_render_loop_in_the_cli_renders_EVERY_rendered_layer():
    """TWO render loops exist -- the target-site (viral) path and the bacterial path -- and wiring only
    one leaves the layer invisible on the other. That is exactly what happened when
    `species_composition` was added 2026-10-05: the viral loop got it, the bacterial loop did not, and
    the reachability guards all passed because `trust_block` carried the block. Only running the real
    bacterial CLI showed nothing printed.

    Same class as the documented salmserovar defect where fixing one of two allele-selection points left
    the bug live. Counting the loops and asserting each is complete is what makes it not recur.
    """
    import re
    from pathlib import Path
    src = (Path(__file__).resolve().parent.parent / "dna_decode" / "amr" / "cli.py").read_text(
        encoding="utf-8")
    loops = re.findall(r"for _extra in \((.*?)\):", src, re.S)
    assert len(loops) >= 2, f"expected >=2 render loops, found {len(loops)} -- has the CLI been restructured?"
    # Renderer names do NOT uniformly derive from layer names (`source_concentration` ->
    # `concentration_one_line`), so the mapping is explicit rather than computed -- a derived name made
    # this guard fail on a correctly-wired layer the first time it ran.
    renderer = {"lineage": "lineage_one_line",
                "source_concentration": "concentration_one_line",
                "prospective": "prospective_one_line",
                "error_rates": "error_rates_one_line",
                "organism_scope": "organism_scope_one_line",
                "species_composition": "species_composition_one_line"}
    assert set(renderer) == _RENDERED, "renderer map drifted from _RENDERED"
    for i, body in enumerate(loops):
        for layer in sorted(_RENDERED):
            fn = renderer[layer]
            assert fn in body, f"render loop #{i} does not call {fn} -- layer {layer!r} is invisible there"


# --- species_composition: the real surface, and the renderer's states ------------------------------
#
# Added by /test-epilogue 2026-10-05. The layer shipped with a STRUCTURAL guard only -- the test above
# asserts `species_composition_one_line` is *named* in both render loops. That cannot see whether the
# function produces the right string, or any string at all: a renderer that returned None for every
# input would satisfy it. The sibling layers (`prospective_one_line`, `concentration_one_line`) each
# have both a real-CLI test and per-state renderer tests; this one had neither.


def _card_species_block(organism: str, drug: str) -> dict | None:
    """The committed card's species block for one cell, or None. Fixture DERIVED, never hardcoded."""
    import json
    p = Path(__file__).resolve().parent.parent / "wiki" / "decoder_validation_report_card.json"
    if not p.exists():
        return None
    for c in json.loads(p.read_text(encoding="utf-8")).get("cells", []):
        if (c.get("organism"), c.get("drug")) == (organism, drug):
            return c.get("species_composition")
    return None


def _species_line(out: str) -> str | None:
    return next((ln for ln in out.splitlines() if "species composition:" in ln), None)


def test_species_composition_reaches_the_printed_output():
    """THE anchor for this layer. `klebsiella x meropenem` is the cell the whole audit exists for: its
    published sens 0.467 comes from a cohort that is 38% K. aerogenes scored as K. pneumoniae.

    If this ever stops printing, the caller reading 0.467 is back to not being able to see that -- and
    the structural guard will still pass, because the renderer is still *named* in the loop.
    """
    if not _card_species_block("klebsiella", "meropenem"):
        pytest.skip("the meropenem cell carries no species block (card not built?)")
    with tempfile.TemporaryDirectory() as td:
        rd = _run_dir(Path(td), [("blaKPC-2", "BETA-LACTAM", "CARBAPENEM", "EXACTX")])
        rc, out = _invoke(["--drug", "meropenem", "--amrfinder-run", str(rd),
                           "--organism", "Klebsiella_pneumoniae"])
    assert rc == 0
    line = _species_line(out)
    assert line is not None, f"species layer never printed:\n{out}"
    assert "MIXED" in line
    assert "Klebsiella aerogenes" in line, "the line must name WHAT the cohort actually is"
    # the on-scored-species figure must ship WITH its disclaimer, never as a replacement metric
    assert "user authority call" in line


def test_a_cohort_with_incomplete_cached_runs_prints_composition_but_NO_outcome_FIGURE():
    """NON-VACUITY for the withheld branch, on real data. Three Klebsiella cohorts have partial cached
    determinant runs, so their per-species outcomes are withheld -- and the danger is a line that still
    quotes an outcome number for them. Attributing outcomes on a partial run set is a rate over a
    silently-shrunken denominator, which is the failure the withholding exists to prevent.
    """
    blk = _card_species_block("klebsiella", "ceftriaxone")
    if not blk or "outcome_attribution_withheld" not in blk:
        pytest.skip("no cohort with a withheld cross-tab on the committed card")
    with tempfile.TemporaryDirectory() as td:
        rd = _run_dir(Path(td), [("blaCTX-M-15", "BETA-LACTAM", "CEPHALOSPORIN", "EXACTX")])
        rc, out = _invoke(["--drug", "ceftriaxone", "--amrfinder-run", str(rd),
                           "--organism", "Klebsiella_pneumoniae"])
    assert rc == 0
    line = _species_line(out)
    assert line is not None, f"species layer never printed:\n{out}"
    assert "WITHHELD" in line
    assert "sens" not in line, f"a withheld cohort must not quote an outcome figure: {line}"


def test_species_composition_is_SILENT_for_a_cohort_that_matches_its_scored_organism():
    """Pairs with the two tests above: without this, a renderer that printed on every call would pass
    them both. A reassuring line for a clean cohort is noise, and invites reading the presence of a
    species line as meaningful when it is unconditional.
    """
    if _card_species_block("escherichia_coli_shigella", "ciprofloxacin") is not None:
        pytest.skip("the E. coli cipro cohort now carries a species block -- pick a clean cell instead")
    with tempfile.TemporaryDirectory() as td:
        rd = _run_dir(Path(td), [("gyrA_S83L", "QUINOLONE", "FLUOROQUINOLONE", "POINTX")])
        rc, out = _invoke(["--drug", "ciprofloxacin", "--amrfinder-run", str(rd),
                           "--organism", "Escherichia_coli_Shigella"])
    assert rc == 0
    assert _species_line(out) is None, f"a cohort matching its scored organism printed a block:\n{out}"


def test_an_unresolved_composition_says_UNKNOWN_rather_than_reading_as_clean():
    """The refusal state. A cohort too thin on resolved species to judge must not be silent -- silence
    is what a CLEAN cohort produces, so reusing it here would make "we could not tell" and "we checked
    and it is fine" indistinguishable at the only surface a human reads.
    """
    line = species_composition_one_line({"species_composition": {
        "status": "insufficient_resolution", "reason": "40/60 assemblies had no readable ORGANISM line"}})
    assert line is not None
    assert "UNRESOLVED" in line and "not as clean" in line
    assert "40/60" in line, "the refusal must carry its own reason"
    # and it must not invent composition figures it does not have
    assert "MIXED" not in line and "%" not in line


def test_the_renderer_is_silent_on_every_state_it_has_nothing_to_say_about():
    assert species_composition_one_line({}) is None
    assert species_composition_one_line({"species_composition": None}) is None
    assert species_composition_one_line({"species_composition": {"status": "matches"}}) is None
    # a future status this renderer does not know must be silent, never half-rendered
    assert species_composition_one_line({"species_composition": {"status": "something_new"}}) is None


def test_an_undefined_on_species_sensitivity_is_omitted_not_printed_as_None():
    """`sens` is None when the scored species' subset has no resistant isolates at all (tp+fn == 0).
    Printing `sens None` to a human reads as a measured value of nothing; the clause must be dropped.
    """
    line = species_composition_one_line({"species_composition": {
        "status": "measured", "scored_as": "Klebsiella_pneumoniae", "n_total": 60, "n_off_species": 23,
        "off_species_fraction": 0.3833, "species_counts": {"Klebsiella pneumoniae": 37},
        "on_scored_species_only": {"n": 37, "fn": 0, "sens": None}}})
    assert line is not None and "MIXED" in line
    assert "None" not in line, f"an undefined sensitivity leaked into the human line: {line}"
    assert "sens" not in line
