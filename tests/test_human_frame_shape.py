"""Guards for `scripts/human_frame_shape.py`.

Two failure modes are worth more attention than the rest.

`--select all` is 7.68 GB against 33.8 GB free on a 93.4%-full disk, so the argv is pinned
rather than trusted. And every wall -- fetch failure, MD5 failure, join rate below threshold --
must exit non-zero with its own code and report NO shape, because a partial shape comparison
reads exactly like a complete one.

The artifact must also never carry the bacterial 0.941 in a human field. Both
`wiki/eukaryotic_arm_result_*.json` currently carry it as context, and repeating that
substitution one level down is how a borrowed ceiling becomes a measured one.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from dna_decode.glm.active_frame import shape_summary
from scripts.human_frame_shape import (
    CELLS,
    EXIT_FETCH,
    EXIT_JOIN_RATE,
    EXIT_MD5,
    NOT_A_REPLICATE,
    REPLICATE_PATTERNS,
    cell_shape,
    fetch_argv,
    main,
    scan_replicate_columns,
)

HUMAN_DIR = Path("data/mpra_human")


# ----------------------------------------------------------------------------------------
# the fetch argv is pinned
# ----------------------------------------------------------------------------------------

def test_fetch_argv_selects_tables_and_overrides_dest():
    argv = fetch_argv("data/mpra_human")
    assert "--select" in argv and argv[argv.index("--select") + 1] == "tables"
    assert "--dest" in argv and argv[argv.index("--dest") + 1] == "data/mpra_human"
    assert "--yes" in argv
    assert "scripts/fetch_lentimpra.py" in argv


def test_fetch_argv_never_requests_the_full_deposit():
    """`all` is 7.68 GB; the disk has 33.8 GB free at 93.4% full."""
    assert "all" not in fetch_argv("x")


def test_fetch_argv_never_points_at_the_detached_D_drive():
    argv = " ".join(fetch_argv("data/mpra_human"))
    assert "D:" not in argv and "dna_decode_cache" not in argv


def test_the_fetcher_really_accepts_the_flags_this_script_passes():
    """Verified against the real parser, not inherited from a reference."""
    src = Path("scripts/fetch_lentimpra.py").read_text(encoding="utf-8")
    for flag in ("--dest", "--select", "--yes"):
        assert f'"{flag}"' in src, f"{flag} is not an argument of the real fetcher"
    assert '"tables"' in src or "tables" in src


# ----------------------------------------------------------------------------------------
# every wall exits with its own code and reports NO shape
# ----------------------------------------------------------------------------------------

def _fake_fetch(returncode, blob):
    def _inner(dest, *, timeout=1800):
        return {"argv": fetch_argv(dest), "returncode": returncode, "ok": returncode == 0,
                "md5_failure_detected": ("md5" in blob.lower() and "mismatch" in blob.lower()),
                "wall_seconds": 0.0, "tail": blob.splitlines()}
    return _inner


def test_a_plain_fetch_failure_exits_with_the_fetch_code_and_withholds_the_shape(
        monkeypatch, tmp_path):
    import scripts.human_frame_shape as m

    monkeypatch.setattr(m, "run_fetch", _fake_fetch(2, "connection reset by peer"))
    out = tmp_path / "a.json"
    rc = main(["--dest", str(tmp_path / "dep"), "--out", str(out)])
    assert rc == EXIT_FETCH
    d = json.loads(out.read_text(encoding="utf-8"))
    assert d["status"] == "INDETERMINATE_FETCH_FAILED"
    assert d["human_reproducibility_curve"] is None
    assert d["partial_measurement_withheld"] is True
    assert "cells" not in d, "a failed fetch must not report per-cell shape"


def test_an_md5_failure_exits_with_its_OWN_code_distinct_from_a_fetch_failure(
        monkeypatch, tmp_path):
    import scripts.human_frame_shape as m

    monkeypatch.setattr(m, "run_fetch",
                        _fake_fetch(1, "ERROR: md5 mismatch for TableS6.txt"))
    out = tmp_path / "b.json"
    rc = main(["--dest", str(tmp_path / "dep"), "--out", str(out)])
    assert rc == EXIT_MD5
    assert rc != EXIT_FETCH, "the two walls must be distinguishable"
    d = json.loads(out.read_text(encoding="utf-8"))
    assert d["status"] == "INDETERMINATE_MD5_VERIFICATION_FAILED"
    assert d["human_ceiling"] is None


def test_no_cell_loading_exits_with_the_join_rate_code_and_withholds_the_shape(
        monkeypatch, tmp_path):
    import scripts.human_frame_shape as m

    monkeypatch.setattr(m, "run_fetch", _fake_fetch(0, "ok"))
    monkeypatch.setattr(m, "cell_shape", lambda c, d, **k: {
        "cell": c, "loaded": False, "refusal": "join rate 0.0000 below 0.95"})
    out = tmp_path / "c.json"
    rc = main(["--dest", str(tmp_path / "dep"), "--out", str(out)])
    assert rc == EXIT_JOIN_RATE
    d = json.loads(out.read_text(encoding="utf-8"))
    assert d["status"] == "INDETERMINATE_JOIN_RATE_BELOW_THRESHOLD"
    assert d["human_reproducibility_curve"] is None
    assert d["human_ceiling"] is None
    assert "bacterial_comparison" not in d


def test_the_three_exit_codes_are_mutually_distinct():
    assert len({EXIT_FETCH, EXIT_MD5, EXIT_JOIN_RATE}) == 3


# ----------------------------------------------------------------------------------------
# the replicate scan, and what it refuses to count
# ----------------------------------------------------------------------------------------

def test_the_scan_reports_absence_rather_than_assuming_it(tmp_path):
    d = tmp_path / "dep"
    d.mkdir()
    (d / "TableS6.txt").write_text("name\tactivity\tfold\nA\t1.0\t3\n", encoding="utf-8")
    out = scan_replicate_columns(d)
    assert out["exists"] is True
    assert out["has_replicates"] is False
    assert out["replicate_like_columns"] == []
    assert "no header matched" in out["reason"]
    assert out["files"][0]["columns"] == ["name", "activity", "fold"]


def test_fold_is_NOT_counted_as_a_replicate():
    """`fold` is the published 10-fold CV assignment; counting it would invent a ceiling.

    The first version of this test probed `cv_fold_1` and FAILED its own vacuity guard: no
    replicate pattern matches that string, so excluding it removes nothing. The exclusion only
    bites where a name BOTH looks like a replicate AND contains `fold`, so that is what is
    tested here -- and `scan_replicate_columns` now reports `fold_exclusion_was_inert` so a run
    on real headers shows whether the control did any work at all.
    """
    name = "fold_rep1"
    assert any(p.search(name) for p in REPLICATE_PATTERNS), (
        "vacuity guard: the test means nothing unless a replicate pattern would have matched"
    )
    assert NOT_A_REPLICATE.search(name), "the exclusion must fire on this name"
    # and the plain CV column must not be mistaken for a replicate in the first place
    assert not any(p.search("cv_fold_1") for p in REPLICATE_PATTERNS)


def test_the_scan_reports_whether_the_fold_exclusion_did_any_work(tmp_path):
    """A filter that removes nothing is not a control, so inertness must be visible."""
    d = tmp_path / "dep"
    d.mkdir()
    (d / "plain.tsv").write_text("name\tactivity\tcv_fold_1\nA\t1\t3\n", encoding="utf-8")
    inert = scan_replicate_columns(d)
    assert inert["fold_exclusion_was_inert"] is True
    assert inert["excluded_as_fold"] == []
    assert "removed NOTHING" in inert["what_inert_means_here"]

    (d / "active.tsv").write_text("name\tfold_rep1\tfold_rep2\nA\t1\t2\n", encoding="utf-8")
    live = scan_replicate_columns(d)
    assert live["fold_exclusion_was_inert"] is False
    assert live["n_excluded_as_fold"] == 2
    assert live["has_replicates"] is False, (
        "a table whose only replicate-shaped columns are fold columns must not claim replicates"
    )


def test_the_scan_finds_genuine_replicate_columns_when_present(tmp_path):
    d = tmp_path / "dep"
    d.mkdir()
    (d / "T.tsv").write_text("name\trep1\trep2\tactivity\nA\t1\t2\t1.5\n", encoding="utf-8")
    out = scan_replicate_columns(d)
    assert out["has_replicates"] is True
    assert any("rep1" in c for c in out["replicate_like_columns"])
    assert out["reason"] is None


def test_the_scan_reports_a_missing_directory_without_raising(tmp_path):
    out = scan_replicate_columns(tmp_path / "nope")
    assert out["exists"] is False
    assert out["has_replicates"] is False
    assert "does not exist" in out["reason"]


# ----------------------------------------------------------------------------------------
# the artifact must not borrow the bacterial ceiling, and must not regrow the retracted field
# ----------------------------------------------------------------------------------------

def test_the_shape_keys_are_exactly_shape_summarys_key_set():
    """So a threshold-Spearman field cannot reappear on the human side under any name."""
    expected = set(shape_summary([1.0, 2.0, 3.0, 4.0]).keys())
    assert "modal_bin_share" in expected
    assert not any("threshold" in k.lower() for k in expected)


def test_a_withheld_run_never_writes_0_941_into_a_human_field(monkeypatch, tmp_path):
    import scripts.human_frame_shape as m

    monkeypatch.setattr(m, "run_fetch", _fake_fetch(2, "network down"))
    out = tmp_path / "d.json"
    main(["--dest", str(tmp_path / "dep"), "--out", str(out)])
    blob = out.read_text(encoding="utf-8")
    assert "0.941" not in blob, "the bacterial ceiling must never appear as a human number"


def test_cells_includes_wtc11_so_its_refusal_is_recorded_rather_than_skipped():
    assert "WTC11" in CELLS, (
        "WTC11 refuses at join rate 0.0000; omitting it would hide a known refusal instead of "
        "recording it"
    )


def test_cell_shape_records_a_refusal_instead_of_raising(monkeypatch):
    """A refusal is data; an exception escaping here would lose the WTC11 observation."""
    from dna_decode.glm import lentimpra as L

    def boom(cell, data_dir, *, min_join_rate=0.95):
        raise L.LentiMpraDataError("join rate 0.0000 below 0.95")

    monkeypatch.setattr(L, "load_cell", boom)
    out = cell_shape("WTC11", "whatever")
    assert out["loaded"] is False
    assert "0.0000" in out["refusal"]
    assert out["refusal_is_expected_for_WTC11"] is True


# ----------------------------------------------------------------------------------------
# real data, gated
# ----------------------------------------------------------------------------------------

@pytest.mark.skipif(not HUMAN_DIR.exists(),
                    reason=f"absent human deposit: {HUMAN_DIR} (run the script to fetch it)")
def test_real_deposit_scan_enumerates_its_headers():
    out = scan_replicate_columns(HUMAN_DIR)
    assert out["exists"] is True
    assert out["n_files"] > 0
    # whichever way it lands, the verdict must be stated rather than left implicit
    assert out["has_replicates"] in (True, False)
    if not out["has_replicates"]:
        assert out["reason"]
