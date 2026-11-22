"""The training corpus must refuse to leak, must deduplicate by sequence, and must demand a domain label.

The leak gate is the reason this module exists and it gets the most tests, because its failure mode is
the worst kind available: a generator trained on the windows it is later scored against produces
excellent numbers, raises nothing, and looks like a result. Every other defect in this arm announced
itself. That one would not.

`test_withholding_one_seed_would_still_leak_into_the_others` is the load-bearing measurement. The
falsifier reshuffles the corpus per seed, so each seed has a DIFFERENT eval half; withholding only the
seed the committed anchor uses would leave most other seeds' scoring windows sitting in training. The
union over seeds is therefore the correct withhold set, and the test measures how much bigger it is
rather than asserting that it matters.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from dna_decode.glm.training_corpus import (
    VALID_DOMAINS,
    CorpusEntry,
    CorpusSpecError,
    assert_one_genome_per_lineage,
    build_training_corpus,
    domain_balance_report,
    falsifier_eval_windows,
)

ROOT = Path(__file__).resolve().parents[1]
REF = ROOT / "data" / "cache" / "refseq" / "GCF_000005845.2"
CACHE = ROOT / "data" / "refseq_cache"

needs_ref = pytest.mark.skipif(not (REF / "genome.fna").exists(),
                               reason="local reference absent on this host")


def _entry(acc: str, ref: Path, *, domain="bacteria", species="E. coli", lineage=None) -> CorpusEntry:
    return CorpusEntry(accession=acc, species=species, domain=domain,
                       lineage_key=lineage or acc, reference_dir=str(ref))


# ------------------------------------------------------------------ the domain label

def test_a_missing_or_unknown_domain_is_refused():
    """The label is required, not an annotation: an unlabelled genome cannot be conditioned on later."""
    with pytest.raises(CorpusSpecError, match="not one of"):
        CorpusEntry("GCF_1", "x", "", "L1", "/tmp/x")
    with pytest.raises(CorpusSpecError, match="not one of"):
        CorpusEntry("GCF_1", "x", "prokaryote", "L1", "/tmp/x")


def test_the_valid_domains_span_all_life_not_just_bacteria():
    """The generator's scope is all life forms; a bacteria-only enum would re-introduce the drift."""
    assert {"bacteria", "archaea", "fungi", "plant", "animal"} <= VALID_DOMAINS


def test_every_valid_domain_actually_constructs():
    for d in sorted(VALID_DOMAINS):
        assert CorpusEntry("GCF_1", "x", d, "L1", "/tmp/x").domain == d


# ------------------------------------------------------------------ lineage keys

def test_a_missing_lineage_key_is_refused():
    with pytest.raises(CorpusSpecError, match="lineage_key"):
        CorpusEntry("GCF_1", "x", "bacteria", "", "/tmp/x")


def test_two_genomes_of_one_lineage_are_refused_not_silently_collapsed():
    """Silently dropping one would hide a selection rule that was never applied."""
    es = [_entry("GCF_1", REF, lineage="L1"), _entry("GCF_2", REF, lineage="L1")]
    with pytest.raises(CorpusSpecError, match="more than once"):
        assert_one_genome_per_lineage(es)


def test_distinct_lineages_pass_the_check():
    assert_one_genome_per_lineage([_entry("GCF_1", REF, lineage="L1"),
                                   _entry("GCF_2", REF, lineage="L2")]) is None


# ------------------------------------------------------------------ THE LEAK GATE

def test_building_without_a_withhold_set_is_REFUSED():
    """There is no safe default. A leaked generator produces beautiful numbers and no error."""
    with pytest.raises(CorpusSpecError, match="no withhold set"):
        build_training_corpus([_entry("GCF_1", REF)])


def test_the_refusal_names_the_explicit_override():
    """A refusal that does not say how to proceed deliberately gets worked around badly."""
    with pytest.raises(CorpusSpecError, match="allow_falsifier_overlap"):
        build_training_corpus([_entry("GCF_1", REF)])


@needs_ref
def test_withheld_sequences_do_not_appear_in_the_corpus():
    hold = falsifier_eval_windows(REF, length=150, limit=1200, seeds=(7,), kind="promoter")
    wins, man = build_training_corpus([_entry("GCF_ref", REF)], withhold=hold)
    assert hold, "the withhold set is empty, so this test would pass vacuously"
    assert not (set(wins) & hold)
    assert man.n_windows_withheld_as_falsifier_eval > 0


@needs_ref
def test_the_leak_gate_is_NOT_vacuous_on_real_data():
    """It must actually remove something. A gate that withholds zero windows is not protecting anything."""
    hold = falsifier_eval_windows(REF, length=150, limit=1200, seeds=(7,), kind="promoter")
    _wins, man = build_training_corpus([_entry("GCF_ref", REF)], withhold=hold)
    assert man.n_windows_withheld_as_falsifier_eval >= 500, man.as_dict()
    assert man.falsifier_leak_gate["enforced"] is True
    assert man.falsifier_leak_gate["n_sequences_withheld"] == len(hold)


@needs_ref
def test_withholding_one_seed_would_still_leak_into_the_others():
    """THE measurement behind using a union over seeds rather than one seed's eval half.

    The falsifier reshuffles per seed, so each seed has a different eval half. Withholding only the
    committed anchor's half would leave most other seeds' scoring windows in training -- and nothing
    would report it.
    """
    one = falsifier_eval_windows(REF, length=150, limit=1200, seeds=(7,), kind="promoter")
    many = falsifier_eval_windows(REF, length=150, limit=1200, seeds=tuple(range(10)),
                                  kind="promoter")
    assert len(many) > len(one) * 2, (len(one), len(many))
    leaked = many - one
    assert len(leaked) > 0
    # and those leaked windows WOULD have landed in training under the one-seed gate
    wins, _man = build_training_corpus([_entry("GCF_ref", REF)], withhold=one)
    assert set(wins) & leaked, "the one-seed gate must demonstrably let other seeds' windows through"


@needs_ref
def test_the_explicit_override_is_honoured_and_recorded():
    """A caller may say in writing that overlap is acceptable; the manifest must then say so too."""
    wins, man = build_training_corpus([_entry("GCF_ref", REF)], allow_falsifier_overlap=True)
    assert wins
    assert man.falsifier_leak_gate["enforced"] is False


# ------------------------------------------------------------------ sequence deduplication

@needs_ref
def test_the_same_genome_twice_contributes_nothing_the_second_time():
    """Measured necessity: 12 E. coli genomes carry only 2.72x one genome's distinct content, so
    without window-level dedup a corpus looks multiples larger than it is."""
    hold = falsifier_eval_windows(REF, length=150, limit=1200, seeds=(7,), kind="promoter")
    es = [_entry("A", REF, lineage="L1"), _entry("B", REF, lineage="L2")]
    wins, man = build_training_corpus(es, withhold=hold)
    assert man.entries[0]["n_windows_kept"] > 0
    assert man.entries[1]["n_windows_kept"] == 0, man.entries[1]
    assert man.entries[1]["n_windows_duplicate_dropped"] > 0
    assert len(wins) == len(set(wins)), "the corpus must contain no duplicate sequence"


@needs_ref
def test_deduplication_is_by_sequence_not_by_coordinate():
    """The same upstream region sits at different coordinates in different assemblies, and a model
    only ever sees the sequence."""
    hold = falsifier_eval_windows(REF, length=150, limit=1200, seeds=(7,), kind="promoter")
    wins, _ = build_training_corpus([_entry("A", REF)], withhold=hold)
    assert len(wins) == len(set(wins))


@needs_ref
def test_the_manifest_arithmetic_reconciles():
    """kept + dropped + withheld must equal the raw total, or a window went missing unaccounted for."""
    hold = falsifier_eval_windows(REF, length=150, limit=1200, seeds=(7,), kind="promoter")
    es = [_entry("A", REF, lineage="L1"), _entry("B", REF, lineage="L2")]
    _wins, man = build_training_corpus(es, withhold=hold)
    d = man.as_dict()
    assert (d["n_windows_kept"] + d["n_windows_duplicate_dropped"]
            + d["n_windows_withheld_as_falsifier_eval"]) == d["n_windows_raw_total"]
    assert d["corpus_bp"] == d["n_windows_kept"] * 150
    assert sum(man.by_domain.values()) == d["n_windows_kept"]


# ------------------------------------------------------------------ determinism

@needs_ref
def test_the_corpus_order_is_deterministic():
    """Two runs must be comparable; the trainer owns shuffling, not the builder."""
    hold = falsifier_eval_windows(REF, length=150, limit=1200, seeds=(7,), kind="promoter")
    a, _ = build_training_corpus([_entry("A", REF)], withhold=hold)
    b, _ = build_training_corpus([_entry("A", REF)], withhold=hold)
    assert a == b


# ------------------------------------------------------------------ the balance report

def test_a_single_domain_corpus_is_flagged_as_such():
    """A 95%-one-domain corpus trains a single-domain model whatever its labels say."""
    from dna_decode.glm.training_corpus import CorpusManifest

    man = CorpusManifest(window_length_bp=150)
    man.by_domain = {"bacteria": 1000}
    r = domain_balance_report(man)
    assert r["single_domain_in_practice"] is True
    assert r["largest_share"] == 1.0
    assert "NOT corrected here" in r["reading"]


def test_a_mixed_corpus_is_not_flagged():
    from dna_decode.glm.training_corpus import CorpusManifest

    man = CorpusManifest(window_length_bp=150)
    man.by_domain = {"bacteria": 500, "plant": 300, "animal": 200}
    r = domain_balance_report(man)
    assert r["single_domain_in_practice"] is False
    assert r["n_domains"] == 3
    assert sum(r["shares"].values()) == pytest.approx(1.0)


def test_an_empty_corpus_reports_empty_rather_than_dividing_by_zero():
    from dna_decode.glm.training_corpus import CorpusManifest

    assert domain_balance_report(CorpusManifest())["status"] == "empty"


# ------------------------------------------------------------------ single-sourcing

def test_the_builder_does_not_reimplement_the_window_definition():
    """A second description of a promoter window is the hazard this module is written to avoid: the
    falsifier would score one definition while the trainer learned another."""
    src = (ROOT / "dna_decode" / "glm" / "training_corpus.py").read_text(encoding="utf-8")
    assert "build_natural_set" in src
    assert "never a second definition" in src
    # the module must not grow its own slicing or revcomp logic
    assert "revcomp" not in src
    assert "_upstream_slice" not in src
