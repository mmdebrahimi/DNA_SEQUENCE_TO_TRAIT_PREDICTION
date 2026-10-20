"""Guards for the sourmash distance backend wired into `compute_lineage_metrics`.

The ten committed Mash cells and a sourmash re-derivation share the same
(organism, drug) cell key. A single sidecar would let one backend silently
overwrite the other's numbers, so the namespace separation is the load-bearing
property here -- not the new backend itself.
"""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
SCRIPT = REPO / "scripts" / "compute_lineage_metrics.py"


def _mod():
    spec = importlib.util.spec_from_file_location("_clm", SCRIPT)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def test_each_backend_owns_a_DISTINCT_sidecar():
    m = _mod()
    paths = [m.sidecar_path_for(b) for b in m.SIDECAR_BY_BACKEND]
    assert len(set(paths)) == len(paths), "two backends share a sidecar; cells would collide"


def test_mash_keeps_the_COMMITTED_path_so_the_default_run_is_unchanged():
    m = _mod()
    assert m.sidecar_path_for("mash") == Path("wiki/provdisjoint_lineage_metrics.json")


def test_sourmash_does_NOT_write_the_committed_mash_sidecar():
    m = _mod()
    assert m.sidecar_path_for("sourmash") != m.sidecar_path_for("mash")


def test_an_unknown_backend_RAISES_rather_than_guessing_a_path():
    m = _mod()
    with pytest.raises(ValueError, match="unknown distance backend"):
        m.sidecar_path_for("mashh")


def test_pointing_sourmash_at_the_mash_sidecar_is_REFUSED(tmp_path, capsys):
    """The cross-backend write is the one that destroys committed numbers."""
    m = _mod()
    rc = m.main([
        "--distance-backend", "sourmash",
        "--sidecar", str(m.sidecar_path_for("mash")),
        "--cohort-dir", str(tmp_path / "nonexistent"),
    ])
    assert rc == 2, "a cross-backend write must fail closed"
    assert "REFUSING" in capsys.readouterr().err


def test_the_refusal_is_NOT_vacuous_the_backends_OWN_sidecar_is_accepted(tmp_path):
    """Proves the guard discriminates: the same call shape passes on the owned path.

    It must NOT return 2 (the refusal code). Any other outcome is fine -- the
    cohort dir is absent, so the run has nothing to do.
    """
    m = _mod()
    rc = m.main([
        "--distance-backend", "sourmash",
        "--sidecar", str(tmp_path / "own.json"),
        "--cohort-dir", str(tmp_path / "nonexistent"),
    ])
    assert rc != 2, "the owned-path call was refused; the guard is over-broad"


def test_the_committed_mash_sidecar_still_holds_its_TEN_complete_cells():
    """Regression floor: the sourmash wiring must not have touched the Mash artifact."""
    art = REPO / "wiki" / "provdisjoint_lineage_metrics.json"
    if not art.is_file():
        pytest.skip("committed Mash sidecar absent on this host")
    cells = json.loads(art.read_text(encoding="utf-8"))["cells"]
    assert len(cells) == 10
    assert all(c["partial"] is False for c in cells), "a Mash cell lost its lineage tier"


def test_an_EARLY_RETURN_cell_does_not_claim_a_backend_that_never_RAN():
    """A reconcile failure returns before the distance producer is reached.

    Stamping the requested backend onto such a cell attributes the cell to a tool
    that computed nothing -- false provenance, which is exactly what the
    `searched_roots` discipline exists to prevent elsewhere in this repo.
    """
    m = _mod()
    cell = m.build_cell(organism="Klebsiella", drug="meropenem", cohort="c",
                        raw={}, raw_reconciled=False, partial=False,
                        n_genomes_missing=0, threshold_results=None)
    assert "distance_backend" not in cell


def test_main_REFUSES_to_report_success_when_the_backend_ran_zero_times(tmp_path, capsys):
    """Exit 0 on an all-early-return run is how a never-exercised path escapes review."""
    m = _mod()
    (tmp_path / "data" / "raw").mkdir(parents=True)
    rc = m.main([
        "--distance-backend", "sourmash",
        "--sidecar", str(tmp_path / "own.json"),
        "--cohort-dir", str(tmp_path / "absent_provdisjoint_x"),
    ])
    assert rc == 4, "a run that exercised the backend zero times must not exit 0"
    assert "ZERO times" in capsys.readouterr().err
