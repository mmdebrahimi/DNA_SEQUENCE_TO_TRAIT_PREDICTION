"""Guards for the semantic (right-QUANTITY, not just right-VALUE) contract-number check.

The load-bearing tests here are the ones that pin the bug this script was BORN with: prose words were
matched by raw substring while the field side was already token-aware, so `Jaccard` labelled a number
as an accuracy, `alphamissense` as a sensitivity, `corroborating` as a correlation and `discovery` as
a coverage. That asymmetry produced 12 flags of which 10 were the checker's own fault -- the same
strict-on-one-side/loose-on-the-other defect the parent audit had. A checker must be held to the
standard it enforces, so these cases are pinned in both directions.
"""

from __future__ import annotations

import json

import pytest

from scripts.contract_number_semantics import (
    LABEL_VOCAB,
    ST_CONFIRMED,
    ST_MISMATCH,
    audit_semantics,
    label_for,
    numeric_leaves,
    path_is_consistent,
    value_matches,
)


def _label(prose: str, token: str) -> str | None:
    """Label the single occurrence of `token` in `prose`, with the real inter-number scoping."""
    from scripts.contract_number_semantics import DECIMAL_RE

    ms = list(DECIMAL_RE.finditer(prose))
    idx = next(i for i, m in enumerate(ms) if m.group(1) == token)
    m = ms[idx]
    prev_end = ms[idx - 1].end(1) if idx else 0
    next_start = ms[idx + 1].start(1) if idx + 1 < len(ms) else None
    return label_for(prose, m.start(1), m.end(1), prev_end, next_start)


# --- the substring bug, pinned case by case -----------------------------------------------------

@pytest.mark.parametrize("prose,token,forbidden", [
    ("mean per-genome Jaccard 0.9967 across genomes", "0.9967", "acc"),       # J-acc-ard
    ("`alphamissense` PTEN 0.539 on human assays", "0.539", "sens"),          # alphamis-sens-e
    ("no-call delta EXACTLY 0.0000, corroborating the rule", "0.0000", "corr"),  # corr-oborating
    ("NOT the discovery +0.155 -- the gain shrank", "0.155", "coverage"),     # dis-cov-ery
    ("no mapping on this substrate (max 0.21)", "0.21", "rate"),              # subst-rate
])
def test_prose_labels_are_word_boundary_guarded(prose, token, forbidden):
    """Each of these returned `forbidden` before the boundary guard, and every one was wrong."""
    assert _label(prose, token) != forbidden


def test_jaccard_labels_as_jaccard_not_accuracy():
    """Non-vacuity for the guard above: the word IS still recognised as its own quantity."""
    assert _label("mean per-genome Jaccard 0.9967 across genomes", "0.9967") == "jaccard"


# --- leading-vs-trailing attribution -------------------------------------------------------------

def test_leading_label_wins_so_each_number_gets_its_own_quantity():
    """`sens 1.00 / spec 0.939` mis-assigned 1.00 to spec under nearest-wins (3 chars vs 6)."""
    prose = "carriership sens 1.00 / spec 0.939 / PPV 0.941, zygosity"
    assert _label(prose, "1.00") == "sens"
    assert _label(prose, "0.939") == "spec"
    assert _label(prose, "0.941") == "prec"


def test_a_label_between_two_numbers_is_not_stolen_by_the_earlier_one():
    """Trailing labels are not used: `2.44 (AMRFinder 2.46), Jaccard 0.4126` called 2.46 a Jaccard.

    2.46 is a mean-genes-per-genome figure. It must come back unlabelled (not measured) rather than
    carrying the next number's quantity.
    """
    prose = "aminoglycoside 9.16 -> 2.44 (AMRFinder 2.46), Jaccard 0.4126 -> 0.7786"
    assert _label(prose, "2.46") is None
    assert _label(prose, "0.4126") == "jaccard"


# --- field-path consistency ----------------------------------------------------------------------

@pytest.mark.parametrize("path,canon,expected", [
    ("metrics.sens", "sens", True),
    ("metrics.spec", "sens", False),
    ("confusion.accession_count", "acc", False),   # accession must not satisfy acc
    ("cohort.specimen_type", "spec", False),       # specimen must not satisfy spec
    ("summary.mean_abs_spearman", "corr", True),
    ("ours.accuracy", "acc", True),
])
def test_path_consistency_is_token_aware(path, canon, expected):
    assert path_is_consistent(path, canon) is expected


# --- value matching + leaf extraction ------------------------------------------------------------

def test_value_matches_at_the_tokens_own_precision():
    assert value_matches("0.993", 0.9928)      # contracts quote rounded figures
    assert value_matches("1.0", 1.0)
    assert not value_matches("0.993", 0.98)


def test_numeric_leaves_excludes_booleans():
    """bool is a subclass of int, so `verified: true` would otherwise enter the pool as 1.0 and
    could 'confirm' a prose 1.0."""
    leaves = dict((p, v) for p, v in numeric_leaves({"verified": True, "n": 3, "m": {"acc": 0.5}}))
    assert "verified" not in leaves
    assert leaves == {"n": 3.0, "m.acc": 0.5}


def test_numeric_leaves_indexes_lists():
    assert ("r[1].x", 2.0) in numeric_leaves({"r": [{"x": 1.0}, {"x": 2.0}]})


# --- non-vacuity: a real mismatch IS caught ------------------------------------------------------

def test_a_planted_mismatch_is_caught_and_a_planted_match_is_not(tmp_path, monkeypatch):
    """The whole point: prose claiming `sens` against an artifact carrying only `spec` must flag.

    Without this the suite could pass on a checker that never flags anything -- which is exactly the
    state the first run was in (0 of 180 measurable, 0 mismatches, and it meant nothing).
    """
    import scripts.contract_number_semantics as mod

    art = tmp_path / "wiki" / "planted_2026-09-28.json"
    art.parent.mkdir(parents=True)
    art.write_text(json.dumps({"metrics": {"spec": 0.993}}), encoding="utf-8")

    class Cell:
        cell_id = "planted:cell"
        claim = "sens 0.993 (wiki/planted_2026-09-28.json)"
        validation_slice = label_provenance = demotion_rule = claim_status = ""

    monkeypatch.setattr(mod, "REPO", tmp_path)
    monkeypatch.setattr(mod, "_load_artifact_text", lambda a: art.read_text(encoding="utf-8"))
    monkeypatch.setattr(mod, "_cited_artifacts", lambda blob: ["wiki/planted_2026-09-28.json"])

    import dna_decode.data.cell_registry as reg
    monkeypatch.setattr(reg, "cells", lambda: [Cell()])

    rep = mod.audit_semantics()
    row = next(r for r in rep["rows"] if r["number"] == "0.993")
    assert row["claimed_label"] == "sens"
    assert row["status"] == ST_MISMATCH, row
    assert any("spec" in p for p in row["competing_paths"])

    # ...and the same number at a label-CONSISTENT path must NOT flag.
    art.write_text(json.dumps({"metrics": {"sens": 0.993}}), encoding="utf-8")
    rep2 = mod.audit_semantics()
    assert next(r for r in rep2["rows"] if r["number"] == "0.993")["status"] == ST_CONFIRMED


# --- integration invariants on the live registry --------------------------------------------------

def test_live_run_is_measurable_and_self_consistent():
    """A zero-measurable run is a plumbing signature, not a clean bill of health -- the first version
    of this script reported exactly that (a JSON loader bug meant no artifact ever parsed)."""
    rep = audit_semantics()
    assert rep["n_numbers"] > 100
    assert rep["n_measurable"] > 0, "zero measurable = the check ran on nothing"
    assert sum(rep["status_counts"].values()) == rep["n_numbers"], "statuses must account for every number"
    assert rep["n_mismatch"] == len(rep["mismatches"])
    for r in rep["mismatches"]:
        assert r["claimed_label"] in LABEL_VOCAB
        assert r["competing_paths"], "a mismatch must name where the value actually sits"


# --- metric bindings: the provenance-pointer path (added 2026-09-28) -----------------------------

def _cell(**kw):
    """A throwaway contract-like object carrying only what audit_semantics reads."""
    class C:
        cell_id = kw.get("cell_id", "planted:cell")
        claim = kw.get("claim", "")
        validation_slice = label_provenance = demotion_rule = claim_status = ""
        metric_bindings = kw.get("metric_bindings", ())
    return C()


def _run_with(monkeypatch, tmp_path, cell, artifact_obj, artifact_name="planted_2026-09-28.json"):
    import scripts.contract_number_semantics as mod
    art = tmp_path / "wiki" / artifact_name
    art.parent.mkdir(parents=True, exist_ok=True)
    art.write_text(json.dumps(artifact_obj), encoding="utf-8")
    monkeypatch.setattr(mod, "REPO", tmp_path)
    monkeypatch.setattr(mod, "_load_artifact_text", lambda a: art.read_text(encoding="utf-8"))
    monkeypatch.setattr(mod, "_cited_artifacts", lambda blob: [f"wiki/{artifact_name}"])
    import dna_decode.data.cell_registry as reg
    monkeypatch.setattr(reg, "cells", lambda: [cell])
    return mod.audit_semantics()


def test_a_binding_pointing_at_the_WRONG_quantity_is_caught(monkeypatch, tmp_path):
    """THE cheapest test that this design is not ceremony: a declaration must not be able to launder a
    number into the wrong quantity. A binding claiming `sens` while pointing at a `spec` path must FAIL,
    otherwise the binding becomes a second truth surface and the audit stops auditing."""
    from dna_decode.data.cell_registry import ContractNumberBinding
    from scripts.contract_number_semantics import ST_MISMATCH

    cell = _cell(claim="sens 0.993 (wiki/planted_2026-09-28.json)",
                 metric_bindings=(ContractNumberBinding(
                     "0.993", "sens", "wiki/planted_2026-09-28.json", "metrics.spec"),))
    rep = _run_with(monkeypatch, tmp_path, cell, {"metrics": {"spec": 0.993}})
    row = next(r for r in rep["rows"] if r["number"] == "0.993")
    assert row["status"] == ST_MISMATCH, row
    assert row["from_binding"] is True
    assert any("not consistent" in p for p in row["competing_paths"]), row["competing_paths"]


def test_a_binding_whose_token_is_absent_from_prose_is_a_DEFECT(monkeypatch, tmp_path):
    """A binding could otherwise describe a number no reader ever sees — auditing the declaration
    instead of the prose, which is exactly the inversion this design exists to prevent."""
    from dna_decode.data.cell_registry import ContractNumberBinding
    from scripts.contract_number_semantics import ST_BOUND_NO_PROSE

    cell = _cell(claim="the prose mentions no decimal at all here",
                 metric_bindings=(ContractNumberBinding(
                     "0.993", "sens", "wiki/planted_2026-09-28.json", "metrics.sens"),))
    rep = _run_with(monkeypatch, tmp_path, cell, {"metrics": {"sens": 0.993}})
    row = next(r for r in rep["rows"] if r["number"] == "0.993")
    assert row["status"] == ST_BOUND_NO_PROSE, row
    assert rep["n_binding_defects"] == 1


def test_a_binding_whose_path_does_not_resolve_is_a_DEFECT(monkeypatch, tmp_path):
    from dna_decode.data.cell_registry import ContractNumberBinding
    from scripts.contract_number_semantics import ST_BOUND_UNRESOLVED

    cell = _cell(claim="sens 0.993 (wiki/planted_2026-09-28.json)",
                 metric_bindings=(ContractNumberBinding(
                     "0.993", "sens", "wiki/planted_2026-09-28.json", "metrics.nonexistent"),))
    rep = _run_with(monkeypatch, tmp_path, cell, {"metrics": {"sens": 0.993}})
    row = next(r for r in rep["rows"] if r["number"] == "0.993")
    assert row["status"] == ST_BOUND_UNRESOLVED, row
    assert rep["n_binding_defects"] == 1


def test_a_correct_binding_confirms_and_is_reported_SEPARATELY_from_the_heuristic(monkeypatch, tmp_path):
    """The split is load-bearing: a rising measurability figure must not be readable as either 'the audit
    got stronger' or 'an author declared more' indifferently."""
    from dna_decode.data.cell_registry import ContractNumberBinding
    from scripts.contract_number_semantics import ST_BOUND

    cell = _cell(claim="the value 0.993 with no quantity word nearby",
                 metric_bindings=(ContractNumberBinding(
                     "0.993", "sens", "wiki/planted_2026-09-28.json", "metrics.sens"),))
    rep = _run_with(monkeypatch, tmp_path, cell, {"metrics": {"sens": 0.993}})
    row = next(r for r in rep["rows"] if r["number"] == "0.993")
    assert row["status"] == ST_BOUND
    assert rep["measurable_by_binding"] == 1
    assert rep["measurable_by_heuristic"] == 0, "this number is unlabeled to the heuristic"


def test_the_live_registry_bindings_convert_previously_unmeasurable_numbers():
    """Non-vacuity on real data, and the REVERT CRITERION for this whole schema change: the three declared
    bindings must actually raise measurability. Baseline before bindings was 16 measurable of 180 (0.0889),
    all by heuristic, with 147 label_unlabeled. If a future change drops these conversions to zero the
    field is ceremony and should be removed.

    HONEST LIMIT, visible in these numbers: the gain is 1:1 per binding (3 declared -> 3 converted), so
    coverage is LINEAR in authoring effort. This is a tool for the numbers you care most about, NOT a fix
    for the remaining ~144 unlabeled ones.
    """
    rep = audit_semantics()
    assert rep["n_bindings_declared"] == 3, rep["n_bindings_declared"]
    assert rep["measurable_by_binding"] == 3, rep["measurable_by_binding"]
    assert rep["measurable_by_heuristic"] == 16, rep["measurable_by_heuristic"]
    assert rep["n_measurable"] == 19 > 16
    assert rep["n_binding_defects"] == 0, rep["binding_defects"]
