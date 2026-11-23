"""Guards for the ruled coordinate-overlap purge.

The load-bearing pair is directional and must stay that way:
  * `test_purged_counts_equal_what_was_RULED` -- the purge removes exactly 56 / 287 / 15.
  * `test_the_UNPURGED_path_is_untouched` -- `load_cell` still returns every joined record, so the
    sha256 regression proving HepG2 and K562 survived the WTC11 rescue byte-identically still holds.
Only one of those failing is the interesting case: a purge that works but silently became the
default would pass the first and fail the second, and that is the hazard this module is about.

Every refusal is tested by CONSTRUCTING the bad input, because none of them can occur in the
committed list -- a guard that cannot be shown to fire is not a guard.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from dna_decode.glm.fold_leak_purge import (  # noqa: E402
    DEFAULT_LIST,
    FoldLeakPurgeError,
    apply_purge,
    purge_names,
    purge_summary,
)
from dna_decode.glm.lentimpra import S3_FILE, S6_FILE, load_cell  # noqa: E402

pytest.importorskip("openpyxl", reason="openpyxl not installed")

#: The RULING, verbatim. Hardcoded here deliberately: a test that read these from the same file the
#: code reads would compare the artifact to itself and pass however wrong the artifact became.
RULED = {"HepG2": 56, "K562": 287, "WTC11": 15}
#: The NARROWER set the rescore actually dropped. Kept so the two can never be conflated in prose.
RESCORED = {"HepG2": 27, "K562": 131}


def _resolve_data_dir():
    from dna_decode.glm.lentimpra import DEFAULT_DIR
    for cand in (DEFAULT_DIR, REPO / "data" / "mpra_human"):
        if (cand / S3_FILE).is_file() and (cand / S6_FILE).is_file():
            return cand
    return DEFAULT_DIR


DATA_DIR = _resolve_data_dir()
HAVE_DATA = (DATA_DIR / S3_FILE).is_file() and (DATA_DIR / S6_FILE).is_file()
needs_data = pytest.mark.skipif(not HAVE_DATA, reason="lentiMPRA tables absent (gitignored)")


# ------------------------------------------------------------------ the list itself

def test_the_committed_list_exists_and_verifies():
    assert DEFAULT_LIST.is_file(), "the purge is defined by a committed list, not a recomputation"
    for cell, n in RULED.items():
        assert len(purge_names(cell)) == n


@pytest.mark.parametrize("cell,n", sorted(RULED.items()))
def test_every_cell_matches_the_ruling(cell, n):
    assert len(purge_names(cell)) == n


def test_the_list_records_that_it_is_NOT_the_rescored_set():
    """The 0.0014 delta describes the narrower set. Conflating them would over-claim the rescue."""
    s = purge_summary()
    assert s["NOT_the_rescored_set"], "the list must say so in its own text"
    for cell, exp in RESCORED.items():
        comp = s["cells"][cell]["composition"]
        assert comp["touched_by_a_TRAINING_fold"] + comp["held_out_vs_held_out_ONLY"] == RULED[cell]
        assert comp["held_out_vs_held_out_ONLY"] > 0, (
            "if this were 0 the ruled and rescored sets would coincide and the warning would be "
            "vacuous -- it is not")
        # The rescore's set was fold 9 ONLY against generator-or-oracle. This purge's
        # training-touched count covers fold 9 AND fold 10, so it must be strictly larger -- that
        # inequality IS the relationship, and asserting mere inequality would be vacuous.
        assert comp["touched_by_a_TRAINING_fold"] > exp, (
            f"{cell}: training-touched {comp['touched_by_a_TRAINING_fold']} should exceed the "
            f"rescore's fold-9-only {exp}, because this set also covers fold 10")


def test_the_purge_touches_only_held_out_folds():
    """Dropping a TRAINING element would change what a model learned -- a larger change."""
    for cell in RULED:
        folds = {int(k) for k in purge_summary()["cells"][cell]["per_fold"]}
        assert folds <= {9, 10}, f"{cell} purges a non-held-out fold: {folds}"


def test_the_authority_is_recorded_not_assumed():
    a = purge_summary()["authority"] or ""
    assert "USER RULING" in a
    assert "relay" in a.lower(), "the decision reached this session via a peer; that must be visible"


# ------------------------------------------------------------------ refusals, each constructed

def test_an_absent_list_REFUSES_instead_of_becoming_a_no_op(tmp_path):
    with pytest.raises(FoldLeakPurgeError, match="absent"):
        purge_names("HepG2", tmp_path / "nope.json")


def test_an_unknown_cell_REFUSES_instead_of_returning_an_empty_set():
    with pytest.raises(FoldLeakPurgeError, match="no purge entry"):
        purge_names("Jurkat")


def test_an_EDITED_list_fails_its_own_digest(tmp_path):
    """Proves the digest guard FIRES. Without this the check could be silently inert."""
    d = json.loads(DEFAULT_LIST.read_text(encoding="utf-8"))
    d["cells"]["HepG2"]["names"] = d["cells"]["HepG2"]["names"][:-1] + ["NOT_A_REAL_ELEMENT"]
    p = tmp_path / "tampered.json"
    p.write_text(json.dumps(d), encoding="utf-8")
    with pytest.raises(FoldLeakPurgeError, match="does not match its own digest"):
        purge_names("HepG2", p)


def test_a_list_whose_count_disagrees_with_the_RULING_refuses(tmp_path):
    import hashlib
    d = json.loads(DEFAULT_LIST.read_text(encoding="utf-8"))
    names = d["cells"]["WTC11"]["names"][:-1]          # one fewer than ruled
    d["cells"]["WTC11"]["names"] = names
    d["cells"]["WTC11"]["n_purged"] = len(names)       # self-consistent, still wrong vs the ruling
    d["sha256_of_names_only"]["WTC11"] = hashlib.sha256("\n".join(names).encode()).hexdigest()
    p = tmp_path / "short.json"
    p.write_text(json.dumps(d), encoding="utf-8")
    with pytest.raises(FoldLeakPurgeError, match="RULING"):
        purge_names("WTC11", p)


def test_a_list_missing_its_digest_refuses(tmp_path):
    d = json.loads(DEFAULT_LIST.read_text(encoding="utf-8"))
    d.pop("sha256_of_names_only")
    p = tmp_path / "nodigest.json"
    p.write_text(json.dumps(d), encoding="utf-8")
    with pytest.raises(FoldLeakPurgeError, match="no sha256"):
        purge_names("HepG2", p)


# ------------------------------------------------------------------ against the real cohort

@needs_data
@pytest.mark.parametrize("cell,n", sorted(RULED.items()))
def test_purged_counts_equal_what_was_RULED(cell, n):
    recs, _ = load_cell(cell, DATA_DIR)
    kept, dropped = apply_purge(recs, cell)
    assert dropped == n
    assert len(kept) == len(recs) - n
    assert {r.name for r in kept}.isdisjoint(purge_names(cell))


@needs_data
@pytest.mark.parametrize("cell", sorted(RULED))
def test_the_UNPURGED_path_is_untouched(cell):
    """The purge is OPT-IN. If this fails the purge became a default and every existing caller's
    numbers moved silently -- the exact failure mode the module documents."""
    recs, st = load_cell(cell, DATA_DIR)
    assert len(recs) == st.n_joined == st.n_measured
    assert st.join_rate == pytest.approx(1.0, abs=1e-9)


@needs_data
def test_every_purged_name_is_actually_in_its_cohort():
    """A name in the list that no cohort contains would mean the list and the data disagree."""
    for cell in RULED:
        recs, _ = load_cell(cell, DATA_DIR)
        have = {r.name for r in recs}
        missing = purge_names(cell) - have
        assert not missing, f"{cell}: {len(missing)} purged names absent from the cohort"
