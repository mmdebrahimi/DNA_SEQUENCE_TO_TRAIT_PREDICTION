"""Guards for the G-D fine-tune split (`scripts/glm_finetune_window_split.py`).

The one failure that matters is LEAKAGE, and it is the failure that would look like SUCCESS: a LoRA
run trained on even one of the falsifier's eval windows would score better on the frozen comparison
set because a generator reproduces what it memorised. So the test that earns its keep is not "the
split runs" -- it is that `verify()` REFUSES a leaky split. A disjointness check only ever observed
passing proves nothing.

Everything here is synthetic. The real corpus needs the MG1655 reference on disk, so the one test
that touches it skips when absent rather than failing for a data reason.
"""
from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "scripts" / "glm_finetune_window_split.py"
REFSEQ = ROOT / "data" / "cache" / "refseq" / "GCF_000005845.2"


def _mod():
    spec = importlib.util.spec_from_file_location("_glm_split", SRC)
    assert spec and spec.loader
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def _baseline(n_fit=600, n_eval=600, n_windows=4340, overlapping=2708,
              strand=None) -> dict:
    return {
        "natural_sets": {"promoter": {
            "n_cds_total": n_windows, "n_windows": n_windows,
            "n_overlapping_cds": overlapping,
            "strand_counts": strand or {"+": 2118, "-": 2222},
        }},
        "fit_eval_split": {"n_fit": n_fit, "n_eval": n_eval},
    }


def _stats(n_windows=4340, overlapping=2708, strand=None) -> dict:
    return {"n_cds_total": n_windows, "n_windows": n_windows, "n_skipped_edge": 0,
            "n_skipped_ambiguous": 0, "n_overlapping_cds": overlapping,
            "strand_counts": strand or {"+": 2118, "-": 2222}}


# --- THE defect the script exists to prevent ----------------------------------------------------

def test_a_LEAKY_train_set_is_REFUSED_even_though_every_count_looks_right():
    """The non-vacuity test. Counts all correct, fingerprint matching -- and one eval window has
    found its way into train. That must refuse, because training on it would read as a WIN."""
    m = _mod()
    fit = [f"FIT{i}" for i in range(600)]
    evl = [f"EVAL{i}" for i in range(600)]
    s = {"fit": fit, "eval": evl, "unused": [], "train": fit + ["EVAL17"],
         "eval_seqs": set(evl)}
    fail, _ = m.verify(s, _stats(), _baseline())
    assert fail, "a train set containing an eval window must be refused"
    assert any("LEAKAGE" in f for f in fail)
    assert any("memorised" in f for f in fail), "the refusal should say WHY it would look like a win"


def test_a_CLEAN_split_passes_so_the_guard_is_not_merely_strict():
    m = _mod()
    fit = [f"FIT{i}" for i in range(600)]
    evl = [f"EVAL{i}" for i in range(600)]
    unused = [f"TAIL{i}" for i in range(3140)]
    s = {"fit": fit, "eval": evl, "unused": unused, "train": fit + unused,
         "eval_seqs": set(evl)}
    fail, notes = m.verify(s, _stats(), _baseline())
    assert fail == [], fail
    assert any("disjoint" in n for n in notes)
    assert any("6.2x" in n or "x the fit half" in n for n in notes)


# --- the split derives train BY EXCLUSION, not by index arithmetic ------------------------------

def test_split_excludes_an_eval_duplicate_from_train_rather_than_trusting_indices():
    """Index arithmetic says fit+unused; content says a window equal to an eval window is out.

    This is not hypothetical -- on the real corpus it removes 3 windows, so the exclusion is live.
    """
    m = _mod()
    wins = [f"W{i}" for i in range(m.LIMIT_NATURAL)] + ["W700", "TAIL_UNIQUE"]
    s = m.split(wins)
    assert "W700" in {m._seq(w) for w in s["eval"]}, "fixture must place W700 in the eval half"
    assert "W700" not in {m._seq(w) for w in s["train"]}, "the eval duplicate must be excluded"
    assert "TAIL_UNIQUE" in {m._seq(w) for w in s["train"]}
    assert len(s["train"]) == len(s["fit"]) + len(s["unused"]) - 1


def test_split_reproduces_the_falsifiers_halves_exactly():
    m = _mod()
    wins = [f"W{i}" for i in range(4340)]
    s = m.split(wins)
    assert len(s["fit"]) == 600 and len(s["eval"]) == 600
    assert len(s["unused"]) == 4340 - m.LIMIT_NATURAL
    assert s["fit"] == wins[:600] and s["eval"] == wins[600:1200]


def test_the_geometry_is_the_LAUNCHERS_not_the_parsers_default():
    """`run_glm_anchor_detached.bat` passes --limit-natural 1200; the parser's default is 2000.
    Reproducing the default would build a different split and prove disjointness of the wrong set."""
    m = _mod()
    assert m.LIMIT_NATURAL == 1200, "must match the launcher's invocation, not argparse"
    assert m.LENGTH_BP == 150 and m.SEED == 7 and m.KIND == "promoter"


# --- corpus drift refuses, because it moves the eval set under the comparison --------------------

@pytest.mark.parametrize("bad", [
    {"n_windows": 4000},
    {"overlapping": 2000},
    {"strand": {"+": 1, "-": 2}},
])
def test_corpus_drift_REFUSES_rather_than_warns(bad):
    m = _mod()
    fit = [f"FIT{i}" for i in range(600)]
    evl = [f"EVAL{i}" for i in range(600)]
    s = {"fit": fit, "eval": evl, "unused": [], "train": fit, "eval_seqs": set(evl)}
    fail, _ = m.verify(s, _stats(**bad), _baseline())
    assert any("drift" in f for f in fail), f"drift {bad} should refuse, got {fail}"


def test_a_moved_fit_eval_size_is_refused():
    m = _mod()
    fit = [f"FIT{i}" for i in range(500)]
    evl = [f"EVAL{i}" for i in range(600)]
    s = {"fit": fit, "eval": evl, "unused": [], "train": fit, "eval_seqs": set(evl)}
    fail, _ = m.verify(s, _stats(), _baseline())
    assert any("fit half is 500" in f for f in fail)


# --- it must not reimplement the window loader (R21: the falsifier owns this resource) -----------

def test_it_imports_build_natural_instead_of_reimplementing_it():
    """A second reader of the same files is the recorded way a solved trap gets re-derived."""
    src = SRC.read_text(encoding="utf-8")
    assert "import glm_generator_falsifier" in src
    assert "gf.build_natural" in src
    for own in ("extract_upstream_windows(", "random_genomic_windows(", "parse_gff3("):
        assert own not in src, f"{own} would be a second loader of the falsifier's own corpus"


def test_corpus_stats_dataclass_is_normalised_to_a_dict():
    """`build_natural` returns a CorpusStats dataclass; comparing it to the JSON fingerprint
    attribute-wise would read None and compare equal to nothing."""
    import dataclasses

    m = _mod()

    @dataclasses.dataclass
    class S:
        n_windows: int = 4340
        strand_counts: dict | None = None

    out = m._as_dict(S(strand_counts={"+": 1}))
    assert isinstance(out, dict) and out["n_windows"] == 4340
    assert m._as_dict({"a": 1}) == {"a": 1}


# --- the real surface, once (R3), skipped when the reference is absent ---------------------------

@pytest.mark.skipif(not (REFSEQ / "genome.fna").exists(), reason="MG1655 reference absent")
def test_on_the_REAL_corpus_the_split_is_clean_and_matches_the_recorded_fingerprint():
    m = _mod()
    import json
    if not m.BASELINE.exists():
        pytest.skip("baseline artifact absent")
    wins, stats = m.load_windows()
    s = m.split(wins)
    fail, _ = m.verify(s, stats, json.loads(m.BASELINE.read_text(encoding="utf-8")))
    assert fail == [], fail
    assert len(s["eval"]) == 600
    assert len(s["train"]) > len(s["fit"]), "the unused tail should enlarge the train set"
