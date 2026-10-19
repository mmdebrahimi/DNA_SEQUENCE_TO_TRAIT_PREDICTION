"""Tests for the sourmash-vs-Mash crosscheck runner.

All offline. The clause evaluator is pure over partition-derived counts, so the
interesting cases (a flipped TRIP, an off-by-one near the bar) are synthetic and
need neither sourmash nor a drive.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest

from dna_decode.eval import rejection_gates
from scripts.sourmash_vs_mash_crosscheck import (
    INSUFFICIENT_DATA,
    OK,
    BarRefused,
    eol_renderings,
    evaluate_clauses,
    load_bar,
    main,
    mash_reference,
    near_bar_labels,
    read_selected,
    genome_roots,
    resolve_genomes,
    trip_labels,
    verify_derived_from,
)

REPO = Path(__file__).resolve().parent.parent
BAR = REPO / "wiki" / "sourmash_vs_mash_crosscheck_bar.json"
TRIP_LABEL = "Klebsiella x ciprofloxacin @ 0.005"


def _bar() -> dict:
    return json.loads(BAR.read_text(encoding="utf-8"))


def _identical_observation(bar: dict) -> dict[str, dict]:
    """What a perfectly agreeing sourmash run would produce."""
    return {
        label: {"minority": row["minority"], "verdict": row["mash_verdict"]}
        for label, row in mash_reference(bar).items()
    }


# --------------------------------------------------------------------------- #
# clause evaluation
# --------------------------------------------------------------------------- #

def test_identical_partitions_satisfy_every_clause():
    bar = _bar()
    outcomes = evaluate_clauses(bar, _identical_observation(bar))
    assert [o.met for o in outcomes] == [True, True, True]
    assert all(not o.failures for o in outcomes)


def test_a_flipped_trip_fails_clause_3_and_clause_1_together():
    """The honest version of 'clause 3 bites'.

    The plan expected a case failing clause 3 while satisfying clause 1. That is
    UNCONSTRUCTIBLE: clause 1 is verdict identity on ALL pairs, so a TRIP flipping
    to PASS violates it too. Clause 3 is logically IMPLIED by clause 1 -- recorded
    in the bar's own `clause_logic` block -- and is retained as a severity marker.
    This test pins the implication rather than asserting an impossible independence.
    """
    bar = _bar()
    observed = _identical_observation(bar)
    # A clustering that splits the TRIP cell's minority class into 4 lineages.
    observed[TRIP_LABEL] = {"minority": 4, "verdict": rejection_gates.PASS}
    outcomes = {o.clause_id: o for o in evaluate_clauses(bar, observed)}
    assert outcomes[3].met is False, "the TRIP no longer reproduces"
    assert outcomes[1].met is False, "verdict identity must fail too (clause 3 ⊂ 1)"
    assert outcomes[2].met is False, "and the near-bar count moved 2 -> 4"
    assert any("MOST DANGEROUS" in str(f.get("severity", ""))
               for f in outcomes[3].failures)


def test_clause_3_is_implied_by_clause_1_over_every_single_pair_flip():
    """Exhaustive: no single-pair flip can fail 3 while 1 holds."""
    bar = _bar()
    base = _identical_observation(bar)
    for label in base:
        observed = dict(base)
        observed[label] = {"minority": 99, "verdict": rejection_gates.PASS}
        outcomes = {o.clause_id: o for o in evaluate_clauses(bar, observed)}
        if outcomes[3].met is False:
            assert outcomes[1].met is False, (
                f"{label}: clause 3 failed while clause 1 held -- the implication "
                "recorded in the bar is wrong"
            )


def test_an_off_by_one_near_the_bar_fails_clause_2_while_clause_1_holds():
    """THIS is the pair that proves clause 2 is not redundant with clause 1.

    Moving a near-bar count from 4 to 5 keeps the G8 verdict at PASS, so verdict
    identity still holds -- and clause 2 catches it. Clause 2, not clause 3, is the
    clause that adds discriminating power over clause 1.
    """
    bar = _bar()
    observed = _identical_observation(bar)
    label = "Escherichia_coli_Shigella x ciprofloxacin @ 0.005"
    assert observed[label]["minority"] == 4
    observed[label] = {"minority": 5, "verdict": rejection_gates.PASS}
    outcomes = {o.clause_id: o for o in evaluate_clauses(bar, observed)}
    assert outcomes[1].met is True, "the verdict did not change, so clause 1 holds"
    assert outcomes[2].met is False, "but the exact count moved, so clause 2 fails"


def test_a_far_from_bar_count_change_is_not_gated_by_clause_2():
    """Clause 2 is scoped to the near-bar set by construction."""
    bar = _bar()
    observed = _identical_observation(bar)
    far = next(lbl for lbl, row in mash_reference(bar).items()
               if lbl not in near_bar_labels(bar) and row["minority"] > 10)
    observed[far] = {"minority": observed[far]["minority"] + 3,
                     "verdict": rejection_gates.PASS}
    outcomes = {o.clause_id: o for o in evaluate_clauses(bar, observed)}
    assert outcomes[1].met is True
    assert outcomes[2].met is True


def test_a_missing_pair_is_a_failure_not_a_pass():
    """An unobserved pair must never count as agreement."""
    bar = _bar()
    observed = _identical_observation(bar)
    del observed[TRIP_LABEL]
    outcomes = {o.clause_id: o for o in evaluate_clauses(bar, observed)}
    assert outcomes[1].met is False
    assert outcomes[3].met is False


def test_the_bar_exposes_exactly_one_trip_and_three_near_bar_pairs():
    bar = _bar()
    assert trip_labels(bar) == {TRIP_LABEL}
    assert len(near_bar_labels(bar)) == 3
    assert TRIP_LABEL in near_bar_labels(bar)


# --------------------------------------------------------------------------- #
# the bar must be a commitment
# --------------------------------------------------------------------------- #

def test_live_bar_loads():
    bar = load_bar(BAR)
    assert bar["frozen_before_any_agreement_number_existed"] is True


def test_a_changed_lineage_table_is_refused_rather_than_scored(tmp_path):
    """A bar derived from a since-changed table is not a commitment about it."""
    lineage_rel = "wiki/provdisjoint_lineage_metrics.json"
    (tmp_path / "wiki").mkdir()
    (tmp_path / lineage_rel).write_text('{"cells": []}', encoding="utf-8")
    bar = _bar()
    bar_path = tmp_path / "bar.json"
    bar_path.write_text(json.dumps(bar), encoding="utf-8")
    with pytest.raises(BarRefused, match="CHANGED"):
        load_bar(bar_path, repo=tmp_path)


def test_a_bar_not_claiming_to_predate_the_number_is_refused(tmp_path):
    bar = _bar()
    bar["frozen_before_any_agreement_number_existed"] = False
    bar_path = tmp_path / "bar.json"
    bar_path.write_text(json.dumps(bar), encoding="utf-8")
    with pytest.raises(BarRefused, match="predate"):
        load_bar(bar_path, repo=REPO)


def test_a_missing_derived_from_table_is_refused(tmp_path):
    bar = _bar()
    bar_path = tmp_path / "bar.json"
    bar_path.write_text(json.dumps(bar), encoding="utf-8")
    with pytest.raises(BarRefused, match="missing"):
        load_bar(bar_path, repo=tmp_path)


# --------------------------------------------------------------------------- #
# cohort resolution
# --------------------------------------------------------------------------- #

def test_read_selected_reads_the_real_headerless_format(tmp_path):
    """`selected.tsv` is bare `<accession>\t<label>` with NO header, exact `R`.

    That contract is asserted by the FROZEN `eval/cohort_manifest.py`, by
    `scripts/compute_lineage_metrics.py` (which produced the Mash numbers this bar is
    derived from), and by `tests/test_ktype_kaptive_concordance.py`. An earlier
    version of this test supplied an `assembly_accession\tphenotype` header -- a file
    shape that does not occur anywhere in this repo -- so it only ever validated the
    reader against its own assumption.
    """
    p = tmp_path / "selected.tsv"
    p.write_text("GCA_1\tR\nGCA_2\tS\nGCA_3\tR\n", encoding="utf-8")
    assert read_selected(p) == {"GCA_1": 1, "GCA_2": 0, "GCA_3": 1}


def test_read_selected_does_not_eat_the_first_accession(tmp_path):
    """A header-assuming reader drops row 1 of a headerless file.

    `tests/test_ktype_kaptive_concordance.py` records that exact bug by name -- "a
    DictReader silently eats the first genome" -- and this script had it.
    """
    p = tmp_path / "selected.tsv"
    p.write_text("GCA_FIRST\tR\nGCA_2\tS\n", encoding="utf-8")
    labels = read_selected(p)
    assert "GCA_FIRST" in labels, "the first accession was consumed as a header row"
    assert labels["GCA_FIRST"] == 1


def test_a_headerless_cohort_does_not_parse_to_all_susceptible(tmp_path):
    """The manufactured-negative regression, asserted as its consequence.

    Parsing the real format with a header + named-phenotype-column reader found no
    phenotype column, so every isolate fell through to susceptible. Zero R means
    minority 0, G8 TRIPs on all 20 pairs, and the run reports its own documented
    exit-1 "a real finding: sourmash does NOT substitute for Mash" -- off a label
    parser, not off a partition. Measured before the fix: the R class vanished
    entirely (2 of 4 -> 0 of 3) and the verdict became TRIP.
    """
    p = tmp_path / "selected.tsv"
    p.write_text(
        "".join("GCA_%d\t%s\n" % (i, "R" if i % 2 else "S") for i in range(1, 11)),
        encoding="utf-8",
    )
    labels = read_selected(p)
    assert len(labels) == 10, "no row may be consumed as a header"
    assert sum(labels.values()) == 5, "the R class must survive the parse"
    minority = min(sum(labels.values()), len(labels) - sum(labels.values()))
    assert rejection_gates.g8_dedup_collapses_balance(
        {"min_effective_lineages": minority}
    ).verdict != rejection_gates.TRIP, "a balanced cohort must not read as collapsed"


@pytest.mark.parametrize("content", [
    "GCA_1\tR\nGCA_2\tS\n",
    "GCA_1\tR\nGCA_2\tS\n\nGCA_3\tX\n",      # blank line + unknown label
    "GCA_1\tR\n",                                  # one row
    "",                                              # empty file
    "GCA_1\tR\textra\nGCA_2\tS\textra\n",      # a third column
    "GCA_1\tresistant\nGCA_2\tS\n",              # a spelled-out label
])
def test_read_selected_agrees_with_the_producer_of_the_mash_numbers(tmp_path, content):
    """THE consistency gate.

    The bar's reference numbers were produced by
    `scripts/compute_lineage_metrics.read_selected`. Labels fix the R/S split, the
    split fixes the effective-lineage counts, and those counts are exactly what
    clauses 1-3 compare. If the two readers disagree on ANY file shape, part of the
    measured "agreement" is a parser difference wearing the costume of a partition
    difference -- the same discipline as computing verdicts THROUGH the imported gate
    rather than restating `MIN_EFFECTIVE_LINEAGES`.

    Asserts EQUALITY rather than specific labels for the odd shapes on purpose:
    pinning a label for a three-column row would assert that the producer's quirk is
    correct, when what must hold is that both sides answer identically.
    """
    from scripts.compute_lineage_metrics import read_selected as producer

    p = tmp_path / "selected.tsv"
    p.write_text(content, encoding="utf-8")
    assert read_selected(p) == producer(p)


def test_read_selected_on_an_empty_file_returns_nothing(tmp_path):
    p = tmp_path / "selected.tsv"
    p.write_text("", encoding="utf-8")
    assert read_selected(p) == {}


def test_resolve_genomes_reports_what_is_missing(tmp_path):
    cache = tmp_path / "refseq"
    (cache / "GCA_1").mkdir(parents=True)
    (cache / "GCA_1" / "genome.fna").write_text(">c\nACGT\n", encoding="utf-8")
    found, missing = resolve_genomes({"GCA_1": 1, "GCA_2": 0}, cache)
    assert set(found) == {"GCA_1"}
    assert missing == ["GCA_2"]


def test_resolve_genomes_with_no_cache_reports_everything_missing():
    found, missing = resolve_genomes({"GCA_1": 1}, None)
    assert found == {} and missing == ["GCA_1"]


# --------------------------------------------------------------------------- #
# the PER-COHORT assembly layout (regression, 2026-10-19)
#
# The ten reference cohorts keep their assemblies at
# `data/raw/<cohort>/refseq/<acc>/`, not in a shared cache. Measured on two
# independent hosts: ZERO of the 490 unique accessions have a directory in a
# global refseq cache. Resolving only against that cache reports
# INSUFFICIENT_DATA on every cohort WITH THE DATA FULLY PRESENT -- the same
# verdict the bar pre-commits for absent genomes, so the bug cannot be told
# apart from the honest outcome by reading the artifact.
# --------------------------------------------------------------------------- #


def test_resolve_genomes_finds_the_per_cohort_layout(tmp_path):
    """THE REGRESSION: assemblies under `<cohort>/refseq/<acc>/` must resolve."""
    cohort_dir = tmp_path / "klebsiella_provdisjoint_meropenem"
    acc = "GCA_000567225.2"
    (cohort_dir / "refseq" / acc).mkdir(parents=True)
    (cohort_dir / "refseq" / acc / "genome.fna").write_text(">c" + chr(10) + "ACGT" + chr(10), encoding="utf-8")
    # No global cache at all -- the per-cohort root alone must be enough.
    found, missing = resolve_genomes({acc: 1}, None, cohort_dir=cohort_dir)
    assert set(found) == {acc}
    assert missing == []


def test_the_per_cohort_lookup_is_not_vacuous(tmp_path):
    """NON-VACUITY: an accession absent from BOTH roots is still reported missing.

    Without this, a lookup that accepted anything would satisfy the test above.
    """
    cohort_dir = tmp_path / "cohort"
    (cohort_dir / "refseq" / "GCA_present.1").mkdir(parents=True)
    (cohort_dir / "refseq" / "GCA_present.1" / "genome.fna").write_text(">c" + chr(10) + "ACGT" + chr(10), encoding="utf-8")
    cache = tmp_path / "refseq_cache"
    cache.mkdir()
    found, missing = resolve_genomes(
        {"GCA_present.1": 1, "GCA_absent.1": 0}, cache, cohort_dir=cohort_dir
    )
    assert set(found) == {"GCA_present.1"}
    assert missing == ["GCA_absent.1"]


def test_a_directory_without_a_fasta_does_not_resolve(tmp_path):
    """An EMPTY accession dir is not a genome -- 214 of this host's 680 cache
    entries are exactly that, and a directory count read as a genome count is
    what overstated local coverage in the first place."""
    cohort_dir = tmp_path / "cohort"
    (cohort_dir / "refseq" / "GCA_1.1").mkdir(parents=True)  # dir, no ">c" + chr(10) + "ACGT" + chr(10)
    found, missing = resolve_genomes({"GCA_1.1": 1}, None, cohort_dir=cohort_dir)
    assert found == {} and missing == ["GCA_1.1"]


def test_genome_roots_searches_the_cohort_before_the_shared_cache(tmp_path):
    """Priority order is load-bearing: the per-cohort root is where the data
    measurably is, so it must be consulted first."""
    cohort_dir = tmp_path / "cohort"
    cache = tmp_path / "refseq_cache"
    assert genome_roots(cache, cohort_dir) == [cohort_dir / "refseq", cache]
    assert genome_roots(None, cohort_dir) == [cohort_dir / "refseq"]
    assert genome_roots(cache, None) == [cache]
    assert genome_roots(None, None) == []


# --------------------------------------------------------------------------- #
# the runner's exit contract
# --------------------------------------------------------------------------- #

def test_an_empty_cohort_root_exits_2_and_withholds_the_number(tmp_path):
    stem = tmp_path / "out"
    code = main(["--cohort-root", str(tmp_path / "nothing"), "--out-stem", str(stem)])
    assert code == 2
    payload = json.loads(Path(f"{stem}.json").read_text(encoding="utf-8"))
    assert payload["status"] == INSUFFICIENT_DATA
    assert payload["agreement_number_withheld"] is True
    assert payload["why_withheld"]
    # No clause may read as met.
    assert all(c["met"] is None for c in payload["clauses"])


def test_insufficient_data_never_writes_a_bar_met_verdict(tmp_path):
    stem = tmp_path / "out"
    main(["--cohort-root", str(tmp_path / "nothing"), "--out-stem", str(stem)])
    text = Path(f"{stem}.md").read_text(encoding="utf-8")
    assert "INSUFFICIENT_DATA" in text
    assert "| 1 | verdict_identity | — |" in text


def test_a_partial_cohort_is_insufficient_not_scored(tmp_path, monkeypatch):
    """A cohort with a deliberate hole must withhold, not score on what resolved."""
    import scripts.sourmash_vs_mash_crosscheck as mod

    lineage = json.loads(
        (REPO / "wiki" / "provdisjoint_lineage_metrics.json").read_text(encoding="utf-8")
    )
    cohort = lineage["cells"][0]["cohort"]
    root = tmp_path / "raw" / cohort
    root.mkdir(parents=True)
    (root / "selected.tsv").write_text(
        "GCA_present\tR\nGCA_absent\tS\n", encoding="utf-8",
    )
    cache = tmp_path / "refseq"
    (cache / "GCA_present").mkdir(parents=True)
    (cache / "GCA_present" / "genome.fna").write_text(">c\nACGT\n", encoding="utf-8")
    monkeypatch.setattr(mod, "find_refseq_cache", lambda repo=REPO: cache)

    stem = tmp_path / "out"
    code = main(["--cohort-root", str(tmp_path / "raw"), "--out-stem", str(stem)])
    assert code == 2
    payload = json.loads(Path(f"{stem}.json").read_text(encoding="utf-8"))
    cell = next(c for c in payload["cells"] if c["cohort"] == cohort)
    assert cell["status"] == INSUFFICIENT_DATA
    assert cell["n_genomes_missing"] == 1
    assert payload["agreement_number_withheld"] is True


def test_sourmash_absent_exits_3_without_a_verdict(tmp_path, monkeypatch):
    """Exit 3 must be reachable, and must not write a bar-met verdict."""
    import scripts.sourmash_vs_mash_crosscheck as mod
    from dna_decode.eval import sourmash_distance as sd

    lineage = json.loads(
        (REPO / "wiki" / "provdisjoint_lineage_metrics.json").read_text(encoding="utf-8")
    )
    cohort = lineage["cells"][0]["cohort"]
    root = tmp_path / "raw" / cohort
    root.mkdir(parents=True)
    # Headerless, and BOTH label classes present: a single-class parse is now refused
    # as INSUFFICIENT_DATA (exit 2) before the distance call is reached, so a one-row
    # fixture could no longer exercise exit 3.
    (root / "selected.tsv").write_text(
        "GCA_present\tR\nGCA_other\tS\n", encoding="utf-8")
    cache = tmp_path / "refseq"
    for acc in ("GCA_present", "GCA_other"):
        (cache / acc).mkdir(parents=True)
        (cache / acc / "genome.fna").write_text(">c\nACGT\n", encoding="utf-8")
    monkeypatch.setattr(mod, "find_refseq_cache", lambda repo=REPO: cache)

    def boom(*a, **k):
        raise sd.SourmashUnavailable("simulated absence")

    monkeypatch.setattr(sd, "compute_sourmash_distances", boom)
    stem = tmp_path / "out"
    assert main(["--cohort-root", str(tmp_path / "raw"), "--out-stem", str(stem)]) == 3
    assert not Path(f"{stem}.json").exists(), "exit 3 must not write a verdict"


def test_the_live_artifact_is_SELF_CONSISTENT_in_whatever_state_it_records():
    """The committed artifact's status, exit code, blocked set and clauses must AGREE.

    REWRITTEN 2026-10-19. The previous version asserted
    `status == INSUFFICIENT_DATA and len(blocked) == 10`, with the reason
    "all ten cohorts are blocked while D: is offline" -- a TRANSIENT state pinned as a
    test. The assemblies arrived the same day (pulled from server2 over SSH), the run
    scored OK on 20/20 pairs, and the assertion inverted. That is the exact trap this
    repo already recorded elsewhere: a test that asserted two pgx regions were
    unmeasurable, which flipped when both were fetched hours later. TEST THE BEHAVIOUR.

    The durable invariant is INTERNAL AGREEMENT, which holds in every state and cannot
    go stale: a verdict must never be reported alongside a blocked cohort, and an
    INSUFFICIENT_DATA artifact must never carry a met clause. Both directions are
    checked, so this cannot pass by the artifact simply being empty.
    """
    path = REPO / "wiki" / "sourmash_vs_mash_crosscheck_2026-10-19.json"
    payload = json.loads(path.read_text(encoding="utf-8"))
    status, code = payload["status"], payload["exit_code"]
    blocked, clauses = payload["blocked"], payload["clauses"]
    met = [c["met"] for c in clauses]

    assert status in {OK, INSUFFICIENT_DATA, "CLAUSE_FAILED"}, status
    assert clauses, "an artifact with no clauses records nothing"

    if status == INSUFFICIENT_DATA:
        assert code == 2
        assert blocked, "INSUFFICIENT_DATA with nothing blocked is a contradiction"
        assert all(m is None for m in met), "a blocked run must not evaluate a clause"
        assert payload["agreement_number_withheld"] is True
    elif status == OK:
        assert code == 0
        assert not blocked, "OK with a blocked cohort is a contradiction"
        assert all(m is True for m in met), "OK requires every clause met"
        assert payload["agreement_number_withheld"] is False
    else:
        assert code == 1
        assert any(m is False for m in met), "CLAUSE_FAILED needs a failed clause"

    # Provenance, only when the artifact is new enough to carry it: a root named as a
    # genome source must actually have supplied genomes. Guarded on PRESENCE so an
    # artifact written before the field existed is not retroactively failed -- the
    # schema-drift version of the transient-state trap this test was rewritten to avoid.
    if status == OK and "genome_sources" in payload:
        sources = payload["genome_sources"] or {}
        assert sum(sources.values()) > 0, "an OK run scored genomes from somewhere"


# --------------------------------------------------------------------------- #
# a degenerate label parse is refused, not scored
# --------------------------------------------------------------------------- #

def test_every_reference_cohort_carries_both_label_classes():
    """Non-vacuity for the guard below: single-class is impossible for these cohorts.

    Measured across all 20 reference pairs: the minimum effective-lineage counts are
    N_R 2 and N_S 7. So a raw parse yielding one class cannot be a collapsed cohort;
    it is a parse that did not match the file.
    """
    rows = _bar()["mash_reference_table"]
    assert len(rows) == 20
    assert all(r["effective_lineage_N_R"] >= 1 and r["effective_lineage_N_S"] >= 1
               for r in rows)


def test_a_single_class_label_parse_is_refused_rather_than_scored(tmp_path, monkeypatch):
    """Defence in depth: an all-R (or all-S) parse must withhold, not produce TRIPs.

    Same shape as the TB caller refusing a VCF that holds data rows yet parses to
    zero calls instead of printing a confident S -- both states otherwise look
    identical downstream, which is what made the label-parse bug invisible.
    """
    import scripts.sourmash_vs_mash_crosscheck as mod

    lineage = json.loads(
        (REPO / "wiki" / "provdisjoint_lineage_metrics.json").read_text(encoding="utf-8")
    )
    cohort = lineage["cells"][0]["cohort"]
    root = tmp_path / "raw" / cohort
    root.mkdir(parents=True)
    root.joinpath("selected.tsv").write_text(
        "GCA_a\tR\nGCA_b\tR\n", encoding="utf-8")  # every isolate one class
    cache = tmp_path / "refseq"
    for acc in ("GCA_a", "GCA_b"):
        (cache / acc).mkdir(parents=True)
        (cache / acc / "genome.fna").write_text(">c\nACGT\n", encoding="utf-8")
    monkeypatch.setattr(mod, "find_refseq_cache", lambda repo=REPO: cache)

    stem = tmp_path / "out"
    assert main(["--cohort-root", str(tmp_path / "raw"), "--out-stem", str(stem)]) == 2
    payload = json.loads(Path(f"{stem}.json").read_text(encoding="utf-8"))
    cell = next(c for c in payload["cells"] if c["cohort"] == cohort)
    assert cell["status"] == INSUFFICIENT_DATA
    assert "single class" in cell["reason"]
    assert payload["agreement_number_withheld"] is True
    assert all(c["met"] is None for c in payload["clauses"]), (
        "no clause may read as met off a degenerate parse"
    )


# --------------------------------------------------------------------------- #
# the advertised exit codes must be reachable
# --------------------------------------------------------------------------- #

def test_a_refused_bar_exits_4_and_writes_no_artifact(tmp_path):
    """Exit 4 is advertised in the module docstring; nothing exercised it.

    `load_bar`'s refusals were tested directly, but not that `main` converts one into
    exit 4 WITHOUT writing a result -- the same pairing already pinned for exit 3,
    where a written artifact would leave a verdict-shaped file behind a refusal.
    """
    bar = _bar()
    bar["frozen_before_any_agreement_number_existed"] = False
    bar_path = tmp_path / "bar.json"
    bar_path.write_text(json.dumps(bar), encoding="utf-8")

    stem = tmp_path / "out"
    assert main(["--bar", str(bar_path), "--out-stem", str(stem)]) == 4
    assert not Path(f"{stem}.json").exists(), "exit 4 must not write a verdict"
    assert not Path(f"{stem}.md").exists()


# --------------------------------------------------------------------------- #
# EOL-tolerant derived-from verification
#
# The pinned digest was taken on a CRLF checkout and
# `wiki/provdisjoint_lineage_metrics.json` carries NO .gitattributes `eol` rule
# (measured: `text: unspecified`, `eol: unspecified`, 715 CRLF). So on any LF
# checkout -- CI, Linux, a clone with different core.autocrlf -- the raw digest is
# 42f944cc... against a pinned 785dcc8a..., and the bar used to refuse with "the
# table has CHANGED" when not one character had.
#
# The fix is in the SCORER, never in the frozen bar: content is verified across EOL
# renderings. These tests pin both halves -- that the LF case is accepted, and that
# a genuinely changed table is STILL refused. The second half is what makes the
# first half safe; without it the fix would be indistinguishable from deleting the
# check.
# --------------------------------------------------------------------------- #

CRLF_B = bytes((13, 10))
LF_B = bytes((10,))


def _sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def test_eol_renderings_puts_raw_first_and_never_duplicates():
    crlf_doc = b"a" + CRLF_B + b"b" + CRLF_B
    labels = [lbl for lbl, _ in eol_renderings(crlf_doc)]
    assert labels[0] == "raw"
    assert labels == ["raw", "lf"], "a CRLF doc needs no separate crlf rendering"
    digests = [_sha(d) for _, d in eol_renderings(crlf_doc)]
    assert len(digests) == len(set(digests))


def test_eol_renderings_of_an_lf_doc_offers_the_crlf_form():
    lf_doc = b"a" + LF_B + b"b" + LF_B
    rendered = dict(eol_renderings(lf_doc))
    assert set(rendered) == {"raw", "crlf"}
    assert rendered["crlf"] == b"a" + CRLF_B + b"b" + CRLF_B


def test_eol_renderings_never_produces_crcrlf():
    """Expanding an already-CRLF doc directly gives CRCRLF; it must normalise first."""
    crlf_doc = b"x" + CRLF_B + b"y" + CRLF_B
    for _, data in eol_renderings(crlf_doc):
        assert bytes((13, 13, 10)) not in data


def test_a_doc_with_no_line_endings_has_exactly_one_rendering():
    assert [lbl for lbl, _ in eol_renderings(b"nolines")] == ["raw"]


def test_verify_derived_from_matches_raw_exactly():
    doc = b"a" + CRLF_B + b"b"
    assert verify_derived_from(doc, _sha(doc)) == "raw"


def test_verify_derived_from_accepts_an_lf_checkout_of_a_crlf_pin():
    """THE REGRESSION. Pin taken on CRLF; the file arrives LF-normalised."""
    crlf_doc = b"a" + CRLF_B + b"b" + CRLF_B
    pin = _sha(crlf_doc)
    lf_doc = crlf_doc.replace(CRLF_B, LF_B)
    assert _sha(lf_doc) != pin, "fixture must actually differ at the byte level"
    assert verify_derived_from(lf_doc, pin) == "crlf"


def test_verify_derived_from_accepts_a_crlf_checkout_of_an_lf_pin():
    """Both directions, because a pin may itself have been taken on either host."""
    lf_doc = b"a" + LF_B + b"b" + LF_B
    pin = _sha(lf_doc)
    crlf_doc = lf_doc.replace(LF_B, CRLF_B)
    assert _sha(crlf_doc) != pin
    assert verify_derived_from(crlf_doc, pin) == "lf"


def test_verify_derived_from_refuses_changed_content_under_every_rendering():
    """NON-VACUITY: EOL tolerance must not become content tolerance."""
    pin = _sha(b"a" + CRLF_B + b"b" + CRLF_B)
    changed = b"a" + CRLF_B + b"CHANGED" + CRLF_B
    assert verify_derived_from(changed, pin) is None
    assert verify_derived_from(changed.replace(CRLF_B, LF_B), pin) is None


def test_verify_derived_from_refuses_a_single_byte_change():
    doc = b"a" + CRLF_B + b"b" + CRLF_B
    pin = _sha(doc)
    assert verify_derived_from(doc.replace(b"b", b"c"), pin) is None


def test_verify_derived_from_refuses_a_missing_pin():
    """A bar with no pinned digest is not a commitment; absence must not read as a pass."""
    assert verify_derived_from(b"anything", None) is None
    assert verify_derived_from(b"anything", "") is None


def test_the_real_committed_table_verifies_under_both_eol_forms():
    bar = json.loads(BAR.read_text(encoding="utf-8"))
    pin = bar["derived_from"]["sha256"]
    raw = (REPO / bar["derived_from"]["path"]).read_bytes()
    assert verify_derived_from(raw, pin) is not None
    assert verify_derived_from(raw.replace(CRLF_B, LF_B), pin) is not None


def test_load_bar_succeeds_on_a_simulated_lf_checkout(tmp_path):
    bar = json.loads(BAR.read_text(encoding="utf-8"))
    src = REPO / bar["derived_from"]["path"]
    (tmp_path / "wiki").mkdir()
    (tmp_path / "wiki" / BAR.name).write_bytes(BAR.read_bytes())
    (tmp_path / "wiki" / src.name).write_bytes(src.read_bytes().replace(CRLF_B, LF_B))

    loaded = load_bar(tmp_path / "wiki" / BAR.name, repo=tmp_path)
    assert loaded["_derived_from_eol_match"] == "crlf"
    assert [c["id"] for c in loaded["clauses"]] == [1, 2, 3]


def test_load_bar_on_an_lf_checkout_still_refuses_a_tampered_table(tmp_path):
    """The end-to-end non-vacuity half: EOL tolerance did not disable the check."""
    bar = json.loads(BAR.read_text(encoding="utf-8"))
    src = REPO / bar["derived_from"]["path"]
    (tmp_path / "wiki").mkdir()
    (tmp_path / "wiki" / BAR.name).write_bytes(BAR.read_bytes())
    table = json.loads(src.read_text(encoding="utf-8"))
    table["cells"][0]["cohort"] = "TAMPERED"
    (tmp_path / "wiki" / src.name).write_text(
        json.dumps(table, indent=2) + chr(10), encoding="utf-8"
    )
    with pytest.raises(BarRefused, match="has CHANGED"):
        load_bar(tmp_path / "wiki" / BAR.name, repo=tmp_path)


def test_load_bar_reports_an_eol_match_on_this_working_copy():
    assert load_bar(BAR)["_derived_from_eol_match"] in {"raw", "lf", "crlf"}


def test_load_bar_does_not_modify_the_frozen_bar_on_disk():
    """The eol match is attached to the RETURNED dict, never written back."""
    before = hashlib.sha256(BAR.read_bytes()).hexdigest()
    load_bar(BAR)
    assert hashlib.sha256(BAR.read_bytes()).hexdigest() == before
