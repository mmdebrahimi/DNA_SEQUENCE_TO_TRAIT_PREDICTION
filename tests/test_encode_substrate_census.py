"""Offline tests for the ENCODE substrate census probe.

NO NETWORK. Almost every test drives the PURE helpers against trimmed fixtures, which is why the
probe splits summarising from fetching in the first place. The one exception drives `main` with the
module's `get` monkeypatched, to exercise a refusal that lives in `census_one` rather than in a pure
helper; it still makes no request. The live run is manual and its artifact is the deliverable.

The load-bearing test here is `test_missing_total_REFUSES`: a census that silently reports zero is
indistinguishable from a broken query, and zero is the answer that would quietly kill the goal this
probe exists to price. So the refusal is asserted to FIRE, not merely to exist.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts import encode_substrate_census as census  # noqa: E402
from scripts.encode_substrate_census import (  # noqa: E402
    CensusError, census_verdict, estimate_volume, facet_terms, read_total, summarise_sizes,
)

# Two trimmed `@graph`-bearing payloads in the real response shape.
PAYLOAD = {
    "total": 3582,
    "@graph": [],
    "facets": [
        {"field": "biosample_ontology.term_name", "title": "Biosample",
         "terms": [{"key": "K562", "doc_count": 120},
                   {"key": "HepG2", "doc_count": 90},
                   {"key": "liver", "doc_count": 40}]},
        {"field": "assembly", "title": "Genome assembly",
         "terms": [{"key": "GRCh38", "doc_count": 3000}, {"key": "hg19", "doc_count": 582}]},
    ],
}
NO_TOTAL = {"@graph": [], "facets": []}


# --------------------------------------------------------------------------- refusals
def test_missing_total_REFUSES():
    """The non-vacuity test: a response with no `total` must raise, never read as zero."""
    with pytest.raises(CensusError) as e:
        read_total(NO_TOTAL, what="probe")
    assert "no `total`" in str(e.value)


def test_non_int_total_refuses():
    with pytest.raises(CensusError):
        read_total({"total": "3582"}, what="probe")


def test_present_total_is_returned():
    assert read_total(PAYLOAD, what="probe") == 3582


# --------------------------------------------------------------------------- facets
def test_facet_terms_reads_the_tissue_axis():
    terms = facet_terms(PAYLOAD, "biosample_ontology.term_name")
    assert len(terms) == 3
    assert {t["key"] for t in terms} == {"K562", "HepG2", "liver"}


def test_absent_facet_is_empty_not_an_error():
    assert facet_terms(PAYLOAD, "does.not.exist") == []


def test_assembly_facet_separate_from_tissue_facet():
    """Two facets must not be confused for one another."""
    assert [t["key"] for t in facet_terms(PAYLOAD, "assembly")] == ["GRCh38", "hg19"]


# --------------------------------------------------------------------------- byte arithmetic
def test_summarise_sizes_median_and_spread():
    s = summarise_sizes([10, 20, 30, 40, 50])
    assert s["n_sampled"] == 5
    assert s["median_bytes"] == 30
    assert s["min_bytes"] == 10 and s["max_bytes"] == 50


def test_summarise_sizes_drops_invalid_and_nonpositive():
    s = summarise_sizes([0, -5, None, "x", 100, 200])  # type: ignore[list-item]
    assert s["n_sampled"] == 2
    assert s["median_bytes"] == 150


def test_summarise_sizes_empty_reports_no_median():
    s = summarise_sizes([])
    assert s["n_sampled"] == 0 and s["median_bytes"] is None


def test_estimate_volume_is_count_times_sampled_median():
    stats = summarise_sizes([100, 200, 300])
    vol = estimate_volume(10, stats)
    assert vol["estimated_total_bytes"] == 200 * 10


def test_estimate_volume_labels_itself_an_estimate():
    """An estimate presented as a measurement is the failure this wording prevents."""
    vol = estimate_volume(10, summarise_sizes([100, 200, 300]))
    assert "SAMPLED" in vol["method"] and "estimate" in vol["method"]


def test_estimate_volume_unavailable_without_a_median():
    vol = estimate_volume(10, summarise_sizes([]))
    assert vol["estimated_total_bytes"] is None and "unavailable" in vol["method"]


def test_estimate_volume_unavailable_with_zero_files():
    vol = estimate_volume(0, summarise_sizes([100, 200]))
    assert vol["estimated_total_bytes"] is None


# --------------------------------------------------------------------------- the rollup
def _row(assay, n_exp, n_files, tissues, gib):
    return {"assay_title": assay, "n_experiments": n_exp, "n_bigwig_files": n_files,
            "tissue_terms_sample": tissues,
            "volume": {"estimated_total_gib": gib}}


def test_rollup_keeps_files_and_experiments_as_separate_denominators():
    """The distinction this probe exists to enforce: files are not experiments."""
    per = {"a": _row("a", 100, 247, ["K562"], 10.0),
           "b": _row("b", 100, 253, ["HepG2"], 5.0)}
    r = census_verdict(per)
    assert r["total_experiments"] == 200
    assert r["total_bigwig_files"] == 500
    assert r["file_to_experiment_ratio"] == 2.5
    assert r["total_experiments"] != r["total_bigwig_files"]


def test_rollup_dedups_the_tissue_union():
    per = {"a": _row("a", 10, 10, ["K562", "liver"], 1.0),
           "b": _row("b", 10, 10, ["K562", "lung"], 1.0)}
    assert census_verdict(per)["distinct_tissue_terms_union_of_samples"] == 3


def test_rollup_sums_volume_only_over_available_estimates():
    per = {"a": _row("a", 10, 10, ["x"], 2.5), "b": _row("b", 10, 10, ["y"], None)}
    assert census_verdict(per)["estimated_total_gib_all_assays"] == 2.5


def test_rollup_carries_an_explicit_non_claim():
    per = {"a": _row("a", 10, 10, ["x"], 1.0)}
    note = census_verdict(per)["what_this_does_NOT_say"]
    assert "SUFFICES" in note and "count" in note


def test_rollup_ratio_is_none_rather_than_a_zero_division():
    per = {"a": _row("a", 0, 0, [], None)}
    assert census_verdict(per)["file_to_experiment_ratio"] is None


# ------------------------------------------------------- the GiB figure a reader actually sees
# `estimated_total_bytes` is asserted above; `estimated_total_gib` was not, and it is the number
# that reaches the console line, the rollup sum and therefore the reader. The rollup test supplies
# its own gib via `_row`, so the bytes -> GiB conversion itself had no coverage.
def test_gib_is_BINARY_not_decimal():
    """A GiB label over a decimal-GB divisor would read ~7% high and nothing would catch it.

    4 x 2**30 bytes is 4.0 GiB and 4.3 GB. Only the binary divisor gives 4.0, so this one
    assertion separates the two conventions the field name commits to.
    """
    vol = estimate_volume(4, summarise_sizes([2**30]))
    assert vol["estimated_total_gib"] == 4.0


def test_gib_rounds_to_one_decimal_on_a_non_round_input():
    """A round input cannot show rounding. 3 x 1.5e9 B is 4.190951... GiB -> 4.2."""
    vol = estimate_volume(3, summarise_sizes([1_500_000_000]))
    assert vol["estimated_total_bytes"] == 4_500_000_000
    assert vol["estimated_total_gib"] == 4.2


# ------------------------------------------------------- the spread the module claims to report
def test_quartiles_are_ordered_within_the_min_max_envelope():
    """The module's stated purpose for these fields is that the estimate's spread be VISIBLE.

    Asserted as an ORDERING, not as pinned indices: `clean[n // 4]` is a crude order statistic, and
    freezing its index would break a legitimate move to `statistics.quantiles`. The ordering is what
    any correct implementation must satisfy.

    Measured, the skew is not uniform: `clean[n // 4]` sits one rank ABOVE nearest-rank only when n
    is divisible by four (n=4, 8, 200) and equals nearest-rank otherwise (n=5, 6, 7, 9). n=8 below
    is therefore a divisible-by-four case, as is the live sample size. Left alone in the source --
    the volume estimate uses only the median, so the skew reaches no shipped figure.
    """
    s = summarise_sizes([10, 20, 30, 40, 50, 60, 70, 80])
    assert s["min_bytes"] <= s["p25_bytes"] <= s["median_bytes"] <= s["p75_bytes"] <= s["max_bytes"]


def test_quartiles_land_near_the_true_quartiles_at_two_hundred_values():
    """200 is the size the live probe draws today (`BIGWIG_SAMPLE`), but the literal is deliberate.

    The bands below are derived for n=200 specifically, so reading the size from the constant would
    leave the test half-coupled: retuning `BIGWIG_SAMPLE` would move the inputs and not the expected
    values, and the failure would look like a quartile defect rather than a stale test.

    A tolerant band rather than an equality: interpolated quantiles would give 50.75 / 150.25 and
    this index arithmetic gives 51 / 151. Both are acceptable; a field reporting the median or the
    max instead of a quartile is not, and that is what this excludes.
    """
    s = summarise_sizes(list(range(1, 201)))
    assert s["n_sampled"] == 200
    assert s["p25_bytes"] in (50, 51)
    assert s["p75_bytes"] in (150, 151)


# ------------------------------------------------------- the SECOND refusal, driven through main
def test_zero_experiments_REFUSES_AND_WRITES_NO_ARTIFACT(tmp_path, monkeypatch):
    """`read_total`'s refusal is covered above; `census_one`'s zero-experiment refusal was not.

    A zero here means the assay_title is wrong or the substrate is absent, and the probe exists to
    price a goal that a silent zero would quietly kill. So this asserts the whole path: the guard
    fires, main maps it to exit 3, and NO artifact is left behind for a later reader to mistake for
    a measurement. Patching `get` keeps it offline, matching the sibling W0 probe's test.
    """
    out = tmp_path / "census.json"
    monkeypatch.setattr(census, "get", lambda params: {"total": 0, "@graph": [], "facets": []})

    assert census.main(["--out", str(out), "--assays", "DNase-seq"]) == 3
    assert not out.exists(), "a refused census must leave no artifact behind"
