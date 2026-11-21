"""The human natural set must refuse rather than degrade, and its split must be the publisher's.

These tests SKIP when the lentiMPRA deposit is absent, following this repo's pattern for gitignored
data, so an absent substrate does not add to the standing failure set. The pure-logic guards are
tested WITHOUT the deposit by monkeypatching the loader, so the length and disjointness rails are
covered even on a host that has never fetched it.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from dna_decode.glm import human_natural as HN
from dna_decode.glm.human_natural import (ADAPTOR_BP, ELEMENT_BP, FULL_BP, HumanNaturalError,
                                          fit_eval_by_fold, human_natural_set)

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "mpra_human"
HAVE_DATA = (DATA / "Table S3 - large-scale lib design.xlsx").is_file()
needs_data = pytest.mark.skipif(
    not HAVE_DATA,
    reason="lentiMPRA deposit absent at data/mpra_human; fetch with scripts/fetch_lentimpra.py")


# --------------------------------------------------------------------------- pure geometry

def test_the_geometry_constants_match_the_publishers_sequence():
    """'15nt 5' adaptor - 200nt element - 15nt 3' adaptor' = 230 bp."""
    assert FULL_BP == 230
    assert ADAPTOR_BP == 15
    assert ELEMENT_BP == 200
    assert FULL_BP - 2 * ADAPTOR_BP == ELEMENT_BP


# --------------------------------------------------------------------------- refusals, no deposit

def test_an_absent_data_dir_refuses_and_names_the_fix(tmp_path):
    with pytest.raises(HumanNaturalError) as e:
        human_natural_set("HepG2", tmp_path / "nope")
    msg = str(e.value)
    assert "fetch_lentimpra.py" in msg, "the failure must be actionable"
    # It must also warn about the trap that cost this plan a wrong premise.
    assert "DEFAULT_DIR" in msg and "detached" in msg


def test_an_unknown_cell_refuses_before_touching_disk():
    with pytest.raises(HumanNaturalError) as e:
        human_natural_set("HeLa", "/definitely/not/a/path")
    assert "unknown cell" in str(e.value)


class _Rec:
    """Minimal stand-in for LentiRecord."""

    def __init__(self, name, seq, fold):
        self.name, self.seq, self.fold = name, seq, fold
        self.activity, self.cell = 0.0, "HepG2"


class _St:
    def __init__(self, folds):
        self.fold_sizes = {f: 1 for f in folds}
        self.n_folds = len(folds)
        self.n_designed = self.n_measured = self.n_joined = len(folds)
        self.join_rate = 1.0
        self.sequence_column = "stub"


def _patch_loader(monkeypatch, records, folds):
    monkeypatch.setattr(HN, "load_cell", lambda cell, p, **kw: (records, _St(folds)))
    monkeypatch.setattr(HN, "_require_dir", lambda d: Path("."))


def test_an_off_length_sequence_REFUSES(monkeypatch):
    """NON-VACUITY on the length rail: a 150 bp row must raise, not be silently kept.

    This is the bacterial length, so it is exactly the mistake that would go unnoticed — a natural
    set at the wrong geometry produces a distinguishability number that looks fine.
    """
    recs = [_Rec("a", "ACGT" * 57 + "AC", 1)]  # 230
    assert len(recs[0].seq) == FULL_BP
    _patch_loader(monkeypatch, recs, [1])
    seqs, _ = human_natural_set("HepG2", "ignored")
    assert len(seqs) == 1  # passes on a well-formed row

    short = [_Rec("b", "ACGT" * 37 + "AC", 1)]  # 150
    assert len(short[0].seq) == 150
    _patch_loader(monkeypatch, short, [1])
    with pytest.raises(HumanNaturalError) as e:
        human_natural_set("HepG2", "ignored")
    assert "150 bp" in str(e.value) and "230 bp" in str(e.value)


def test_a_non_ACGT_sequence_REFUSES(monkeypatch):
    bad = [_Rec("n", "N" + "ACGT" * 57 + "A", 1)]
    assert len(bad[0].seq) == FULL_BP
    _patch_loader(monkeypatch, bad, [1])
    with pytest.raises(HumanNaturalError) as e:
        human_natural_set("HepG2", "ignored")
    assert "alphabet_valid_rate" in str(e.value)


def test_the_alphabet_check_reads_RATES_not_an_invented_boolean(monkeypatch):
    """Regression pin. The first draft gated on a non-existent `ok` key, which defaulted to False
    and made a clean set raise. Assert a clean set PASSES."""
    recs = [_Rec(str(i), "ACGT" * 57 + "AC", 1) for i in range(20)]
    _patch_loader(monkeypatch, recs, [1])
    seqs, st = human_natural_set("HepG2", "ignored")
    assert len(seqs) == 20
    assert st["alphabet_and_length"]["alphabet_valid_rate"] == 1.0
    assert st["alphabet_and_length"]["exact_length_rate"] == 1.0


def test_eval_folds_consuming_every_fold_REFUSES(monkeypatch):
    recs = [_Rec(str(i), "ACGT" * 57 + "AC", f) for i, f in enumerate([1, 2])]
    _patch_loader(monkeypatch, recs, [1, 2])
    with pytest.raises(HumanNaturalError) as e:
        fit_eval_by_fold("HepG2", "ignored", eval_folds={1, 2})
    assert "nothing to fit" in str(e.value)


def test_an_unavailable_eval_fold_REFUSES(monkeypatch):
    recs = [_Rec("a", "ACGT" * 57 + "AC", 1)]
    _patch_loader(monkeypatch, recs, [1])
    with pytest.raises(HumanNaturalError) as e:
        fit_eval_by_fold("HepG2", "ignored", eval_folds={99})
    assert "not in the publisher's folds" in str(e.value)


# --------------------------------------------------------------------------- real deposit

@needs_data
def test_the_real_deposit_loads_at_the_publishers_geometry():
    seqs, st = human_natural_set("HepG2", DATA, limit=200)
    assert st["bp"] == FULL_BP
    assert all(len(s) == FULL_BP for s in seqs)
    assert st["geometry"] == "full_230nt_as_shipped"
    assert st["join_rate"] == 1.0
    assert st["folds_available"] == list(range(1, 11)), "the publisher ships 10 folds"


@needs_data
def test_stripping_adaptors_yields_the_200nt_element():
    full, _ = human_natural_set("HepG2", DATA, limit=10)
    strip, st = human_natural_set("HepG2", DATA, limit=10, strip_adaptors=True)
    assert st["bp"] == ELEMENT_BP
    assert all(len(s) == ELEMENT_BP for s in strip)
    # The element must be the INTERIOR of the full sequence, not a re-read.
    assert strip[0] == full[0][ADAPTOR_BP:FULL_BP - ADAPTOR_BP]


@needs_data
def test_fold_selection_round_trips():
    a, sa = human_natural_set("HepG2", DATA, folds={3})
    b, sb = human_natural_set("HepG2", DATA, folds={3, 4})
    assert sa["folds_selected"] == [3] and sb["folds_selected"] == [3, 4]
    assert sb["n"] > sa["n"], "adding a fold must add records"


@needs_data
def test_the_split_is_disjoint_on_CONTENT_not_just_on_fold_ids():
    fit, ev, st = fit_eval_by_fold("HepG2", DATA, eval_folds={1},
                                   limit_fit=600, limit_eval=300)
    assert st["eval_folds"] == [1]
    assert 1 not in st["fit_folds"]
    assert st["content_overlap"] == 0
    assert not (set(fit) & set(ev))
    assert st["split_is_the_publishers"] is True


@needs_data
def test_WTC11_still_refuses_with_the_loaders_own_reason():
    """A sudden WTC11 pass would mean the join-rate guard had broken, not that data improved."""
    with pytest.raises(HumanNaturalError) as e:
        human_natural_set("WTC11", DATA)
    msg = str(e.value)
    assert "WTC11" in msg
    assert "0.95" in msg or "join" in msg.lower()


@needs_data
def test_both_usable_cells_load():
    for cell in ("HepG2", "K562"):
        seqs, st = human_natural_set(cell, DATA, limit=50)
        assert len(seqs) == 50 and st["cell"] == cell
