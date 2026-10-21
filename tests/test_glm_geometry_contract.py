"""The frozen-geometry comparability contract.

WHY (2026-10-20): the sampling-sweep plan said "all other flags identical to Step 3" and
"recorded geometry", and nothing checked it. The reproduction anchor was launched without
`--bins`, taking the runner's DEFAULT of 10 against a baseline recorded at 6. The `kmer` arm
(which does not use positional bins) reproduced BIT-IDENTICALLY, while the `both`-mode Markov
null moved 0.6679 -> 0.7475 -- AWAY from the pre-registered comparator -- so a 25-minute run
could not be compared to the thing it existed to be compared to.

Same shape as the dead temperature knob one layer out: a run whose validity rests on a
variable being held, with nothing asserting it was held.
"""
from __future__ import annotations

import importlib.util
import json
import pathlib

import pytest

from dna_decode.glm.knob_contract import (
    COMPARISON_GEOMETRY_FIELDS,
    GeometryMismatchError,
    assert_geometry_comparable,
    geometry_diff,
    geometry_of,
)

ROOT = pathlib.Path(__file__).resolve().parent.parent

# A minimal artifact carrying every frozen field, so a test that changes ONE field is
# changing only that field.
_BASE = {
    "window_length_bp": 150,
    "kmer_k": 4,
    "positional_bins": 6,
    "markov_order": 3,
    "n_generated_per_candidate": 120,
    "natural_set_arg": "promoter",
    "feature_modes": "kmer,both",
    "conditional": True,
    "context_bp": 300,
}


def _with(**over):
    return {**_BASE, **over}


# ------------------------------------------------------------------ the failure that bit
def test_the_bins_mismatch_that_VOIDED_a_real_run_is_REFUSED():
    with pytest.raises(GeometryMismatchError, match="positional_bins"):
        assert_geometry_comparable(_BASE, _with(positional_bins=10))


def test_the_refusal_names_BOTH_values_so_the_reader_need_not_diff_by_hand():
    with pytest.raises(GeometryMismatchError) as e:
        assert_geometry_comparable(_BASE, _with(positional_bins=10))
    msg = str(e.value)
    assert "recorded=6" in msg and "this run=10" in msg
    # and it must carry the measured consequence, not just "they differ"
    assert "0.6679" in msg and "0.7475" in msg


def test_the_REAL_committed_artifacts_reproduce_the_diagnosis():
    """Pin the diagnosis to the actual files, not to a fixture written to agree with it."""
    rec = ROOT / "wiki" / "glm_generator_falsifier_REAL_2026-10-07.json"
    bad = ROOT / "wiki" / "glm_sweep_anchor_bins10_UNCOMPARABLE_2026-10-20.json"
    if not (rec.is_file() and bad.is_file()):
        pytest.skip("committed artifacts absent")
    d = geometry_diff(json.loads(rec.read_text(encoding="utf-8")),
                      json.loads(bad.read_text(encoding="utf-8")))
    assert d == {"positional_bins": (6, 10)}, d


# ------------------------------------------- non-vacuity: it must ACCEPT the right runs
# A guard that refuses everything is exactly as useless as one that refuses nothing, and far
# more annoying. These are the tests that make the refusals above mean something.
def test_an_IDENTICAL_geometry_passes():
    assert geometry_diff(_BASE, dict(_BASE)) == {}
    assert_geometry_comparable(_BASE, dict(_BASE))


def test_a_DIFFERING_SAMPLING_setting_is_EXEMPT_because_it_is_the_SWEPT_variable():
    """The whole point of the sweep is to move sampling; refusing on it refuses every point."""
    a = {**_BASE, "sampling": {"temperature": 1.0, "top_k": 0}}
    b = {**_BASE, "sampling": {"temperature": 0.6, "top_k": 0}}
    assert_geometry_comparable(a, b)


def test_a_DIFFERING_DEVICE_is_EXEMPT_on_measured_grounds():
    """Excluded on evidence: the kmer arm reproduced bit-identically across GTX 860M -> CPU."""
    a = {**_BASE, "device": {"device_name": "NVIDIA GeForce GTX 860M", "cuda": True}}
    b = {**_BASE, "device": {"device_name": None, "cuda": False}}
    assert_geometry_comparable(a, b)


@pytest.mark.parametrize("noise_field", ["date", "results", "honest_limits", "record"])
def test_outputs_and_metadata_are_EXEMPT(noise_field):
    assert_geometry_comparable({**_BASE, noise_field: "x"}, {**_BASE, noise_field: "y"})


# ------------------------------------------------------------------ every frozen field bites
@pytest.mark.parametrize("field,other", [
    ("window_length_bp", 200),
    ("kmer_k", 5),
    ("positional_bins", 10),
    ("markov_order", 2),
    ("n_generated_per_candidate", 400),
    ("natural_set_arg", "genomic"),
    ("feature_modes", "kmer"),
    ("conditional", False),
    ("context_bp", 500),
])
def test_EVERY_declared_frozen_field_is_actually_checked(field, other):
    """Membership in the tuple is a claim; this is what makes it one.

    Without this a field could be listed and silently skipped by the diff -- the
    "declared but inert" pattern this whole module exists to prevent.
    """
    with pytest.raises(GeometryMismatchError, match=field):
        assert_geometry_comparable(_BASE, _with(**{field: other}))


def test_the_fixture_covers_the_declared_field_set_exactly():
    """If a field joins the tuple and not `_BASE`, the test above goes vacuous for it."""
    assert set(_BASE) == set(COMPARISON_GEOMETRY_FIELDS)


# ------------------------------------------------------------- a MISSING field is not a match
def test_an_ABSENT_field_is_REPORTED_not_DEFAULTED():
    """Defaulting a missing field is how two runs come to look comparable when one never
    recorded the value -- silent agreement, which is worse than a loud mismatch."""
    thin = {k: v for k, v in _BASE.items() if k != "positional_bins"}
    d = geometry_diff(_BASE, thin)
    assert d == {"positional_bins": (6, "<absent>")}, d
    with pytest.raises(GeometryMismatchError, match="absent"):
        assert_geometry_comparable(_BASE, thin)


def test_geometry_of_does_not_invent_keys():
    assert geometry_of({}) == {}
    assert set(geometry_of(_BASE)) == set(COMPARISON_GEOMETRY_FIELDS)


# ------------------------------------------------------- ordering inside the real runner
def _load_runner():
    spec = importlib.util.spec_from_file_location(
        "_glm_falsifier_geom", str(ROOT / "scripts" / "glm_generator_falsifier.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


class _CorpusReached(Exception):
    """Raised from the first expensive step the guard must precede."""


def _runner(tmp_path, monkeypatch):
    mod = _load_runner()
    ref = tmp_path / "GCF_000005845.2"
    ref.mkdir()
    (ref / "genome.fna").write_text(">stub\nACGT\n", encoding="utf-8")
    monkeypatch.setattr(mod, "REFSEQ", ref)

    def _tripwire(*_a, **_k):
        raise _CorpusReached("build_natural ran")

    monkeypatch.setattr(mod, "build_natural", _tripwire)
    return mod


def _baseline_file(tmp_path, **over):
    p = tmp_path / "baseline.json"
    p.write_text(json.dumps(_with(**over)), encoding="utf-8")
    return str(p)


_ARGS = ["--natural-set", "promoter", "--conditional", "--context", "300",
         "--length", "150", "--k", "4", "--markov-order", "3", "--n", "120",
         "--modes", "kmer,both"]


def test_the_runner_REFUSES_a_drifted_geometry_BEFORE_building_a_corpus(tmp_path, monkeypatch):
    mod = _runner(tmp_path, monkeypatch)
    with pytest.raises(GeometryMismatchError):
        mod.main([*_ARGS, "--bins", "10", "--compare-to", _baseline_file(tmp_path)])


def test_and_that_ordering_proof_is_NOT_vacuous(tmp_path, monkeypatch):
    """Matching geometry must REACH the corpus build, so the guard is the only difference."""
    mod = _runner(tmp_path, monkeypatch)
    with pytest.raises(_CorpusReached):
        mod.main([*_ARGS, "--bins", "6", "--compare-to", _baseline_file(tmp_path)])


def test_WITHOUT_compare_to_the_runner_does_NOT_refuse(tmp_path, monkeypatch):
    """Opt-in: the guard must not break every existing invocation that cites no baseline."""
    mod = _runner(tmp_path, monkeypatch)
    with pytest.raises(_CorpusReached):
        mod.main([*_ARGS, "--bins", "10"])


# ===================================================================================== #
# Step 1 of the pre-grid-defects plan: the geometry is ONE expression, not two.
# ===================================================================================== #
# The original defect was not a wrong value, it was a SECOND description of the same run.
# The artifact recorded `positional_bins` as `a.bins`; the --compare-to projection recorded
# `None` unless the mode was positional/both. A `--modes kmer` run compared against its OWN
# kmer-only baseline was therefore refused. Tests that only check values would pass again the
# moment someone re-adds a parallel literal that happens to agree, so the pin is on the SOURCE.
import ast as _ast
import types as _types


def _runner_source_tree():
    return _ast.parse((ROOT / "scripts" / "glm_generator_falsifier.py").read_text(encoding="utf-8"))


def _artifact_dict_node(tree):
    """The `art = {...}` dict literal inside `main()`."""
    for node in _ast.walk(tree):
        if (isinstance(node, _ast.Assign)
                and any(isinstance(t, _ast.Name) and t.id == "art" for t in node.targets)
                and isinstance(node.value, _ast.Dict)):
            return node.value
    raise AssertionError("could not find the `art = {...}` artifact dict in the runner")


def _ns(**over):
    base = dict(length=150, k=4, bins=10, markov_order=3, n=120,
                natural_set="promoter", modes="kmer,both", conditional=True, context=300)
    base.update(over)
    return _types.SimpleNamespace(**base)


def test_geometry_from_args_covers_the_declared_field_set_EXACTLY():
    """A field in the tuple but missing here would silently compare as `<absent>`."""
    mod = _load_runner()
    assert set(mod.geometry_from_args(_ns())) == set(COMPARISON_GEOMETRY_FIELDS)


def test_the_MEASURED_false_refusal_is_gone():
    """A kmer-only run against its own kmer-only baseline used to refuse with (10, None)."""
    mod = _load_runner()
    ns = _ns(modes="kmer")
    recorded = mod.geometry_from_args(ns)          # what the artifact writes for this run
    assert geometry_diff(recorded, mod.geometry_from_args(ns)) == {}


def test_the_REAL_bins_mismatch_still_refuses_after_the_refactor():
    """Fixing the false refusal must not cost the true one."""
    rec = ROOT / "wiki" / "glm_generator_falsifier_REAL_2026-10-07.json"
    if not rec.is_file():
        pytest.skip("recorded baseline absent")
    mod = _load_runner()
    recorded = json.loads(rec.read_text(encoding="utf-8"))
    assert geometry_diff(recorded, mod.geometry_from_args(_ns(bins=6))) == {}
    assert geometry_diff(recorded, mod.geometry_from_args(_ns(bins=10))) == {"positional_bins": (6, 10)}


def test_the_ARTIFACT_restates_NONE_of_the_frozen_fields_as_its_own_literal():
    """The pin that would have caught the original divergence.

    Two parallel literals can agree today and drift tomorrow, so this asserts the artifact
    carries NO literal key among the nine — they can only arrive through the shared spread.
    It also closes the duplicate-key trap: a restated key placed AFTER the spread silently
    WINS over it, so an agreeing duplicate would re-create the second expression invisibly.
    """
    art = _artifact_dict_node(_runner_source_tree())
    restated = [k.value for k in art.keys
                if isinstance(k, _ast.Constant) and k.value in COMPARISON_GEOMETRY_FIELDS]
    assert restated == [], f"artifact restates frozen geometry fields: {restated}"


def test_and_the_artifact_DOES_take_the_shared_spread():
    """Non-vacuity for the test above: absent this, deleting the fields entirely would pass."""
    art = _artifact_dict_node(_runner_source_tree())
    spreads = [v for k, v in zip(art.keys, art.values) if k is None]      # `**x` => key None
    calls = [s for s in spreads
             if isinstance(s, _ast.Call) and getattr(s.func, "id", None) == "geometry_from_args"]
    assert len(calls) == 1, f"expected exactly one **geometry_from_args(...) spread, got {len(calls)}"


def test_BOTH_consumers_read_the_same_function():
    """The guard and the artifact must both call it — one call site would mean one is stale."""
    tree = _runner_source_tree()
    n = sum(1 for node in _ast.walk(tree)
            if isinstance(node, _ast.Call) and getattr(node.func, "id", None) == "geometry_from_args")
    assert n == 2, f"expected 2 call sites (--compare-to guard + artifact), found {n}"
