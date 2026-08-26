"""The colour-cell substrate screen: pure classifier + the trust-surface corrections it forced.

WHY THIS EXISTS. The animal colour/plumage family reached 19 CLI cells, all KNOWLEDGE_BASELINE, before
anyone asked whether they COULD be validated. The screen (scripts/colour_cell_substrate_screen.py) derives
that per-cell from the committed catalogs and found two walls: 40 of 65 loci (62%) record no causal variant
at all, and 14 of the 25 that DO are indel/structural -- off-panel for any SNP array.

These tests pin the classifier (a text heuristic, so it needs anchoring) + the corrections that landed.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

# The pure logic moved to `dna_decode/pigment/substrate_screen.py` (the freeze module cannot import from
# `scripts/`); the script keeps the self-check anchor + the catalog-gap record + artifact writing.
import colour_cell_substrate_screen as screen_cli  # noqa: E402
from colour_cell_substrate_screen import _CATALOG_GAPS, _DOG_TRUTH, self_check  # noqa: E402
from dna_decode.pigment.substrate_screen import (  # noqa: E402
    _loci_of, _source_of, classify_variant, collect, snv_panel_scorable, summarise,
    trait_for_species, verdicts,
)

SCREEN_MD = ROOT / "wiki" / "colour_cell_substrate_screen_2026-08-26.md"
NEG_MAP = ROOT / "wiki" / "negative_results_map_2026-06-13.md"


# ------------------------------------------------------------------------ the classifier

@pytest.mark.parametrize("text,expect", [
    ("MC1R p.Arg306Ter = recessive red `e`", "SNV"),
    ("TYR c.604C>G p.His202Asp, recessive", "SNV"),
    ("MLPH c.-22G>A `d`", "SNV"),                      # negative-offset promoter coordinate
    ("CBD103 beta-defensin c.67_69delGGT", "INDEL"),
    ("ASIP non-agouti black c.181_184delTTCA", "INDEL"),
    ("exon-2 frameshift", "INDEL"),
    ("ASIP A^y/a^t SINE insertion + coding", "STRUCTURAL"),
    ("a 190kb duplication", "STRUCTURAL"),
    ("rabbit C (TYR): C full > chinchilla > Himalayan > c albino", "UNRECORDED"),
    ("mouse a locus (ASIP): A agouti > at > a non-agouti", "UNRECORDED"),
])
def test_classify_variant_cases(text, expect):
    assert classify_variant(text) == expect


def test_structural_beats_indel_because_a_sine_insertion_also_matches_ins():
    """ORDER IS LOAD-BEARING. 'SINE insertion' contains 'insertion'; if INDEL were tested first, every
    structural variant would be mis-filed as an indel and the STRUCTURAL count would read zero."""
    assert classify_variant("ASIP SINE insertion") == "STRUCTURAL"
    assert classify_variant("a 190kb duplication of ASIP") == "STRUCTURAL"


def test_a_frameshift_that_also_cites_a_point_coordinate_is_still_an_indel():
    """Real provenance strings mix notations in one sentence; the coarser class must win."""
    assert classify_variant("MLPH c.667_668insC frameshift (see also c.-22G>A)") == "INDEL"


def test_snv_panel_scorability_is_tri_state_not_boolean():
    """UNRECORDED must be None, never False -- 'we did not write the variant down' is NOT the same claim
    as 'a SNP panel cannot carry it', and collapsing them would fabricate evidence about the substrate."""
    assert snv_panel_scorable("SNV") is True
    assert snv_panel_scorable("INDEL") is False
    assert snv_panel_scorable("STRUCTURAL") is False
    assert snv_panel_scorable("UNRECORDED") is None


def test_empty_source_is_unrecorded_not_a_crash():
    assert classify_variant("") == "UNRECORDED"


# ------------------------------------------------------------------------ the summary verdicts

def _rows(*classes):
    return [{"locus": str(i), "variant_class": c} for i, c in enumerate(classes)]


@pytest.mark.parametrize("classes,verdict", [
    (("UNRECORDED", "UNRECORDED"), "UNSCREENABLE_NO_CAUSAL_VARIANTS_RECORDED"),
    (("SNV", "SNV"), "FULLY_SNV_TRACTABLE"),
    (("SNV", "UNRECORDED"), "SNV_TRACTABLE_WHERE_RECORDED"),
    (("INDEL", "STRUCTURAL"), "NO_LOCUS_SNV_TRACTABLE"),
    (("SNV", "INDEL"), "PARTIALLY_SNV_TRACTABLE"),
])
def test_summarise_verdicts(classes, verdict):
    assert summarise(_rows(*classes))["verdict"] == verdict


def test_blocked_count_excludes_unrecorded():
    """A cell with no recorded variants has ZERO blocked loci -- it is unscreenable, not blocked. Counting
    unrecorded as blocked would overstate the substrate wall, which is the whole error this screen exists
    to avoid making."""
    s = summarise(_rows("UNRECORDED", "UNRECORDED", "SNV"))
    assert s["n_snv_panel_blocked"] == 0


# ------------------------------------------------------------------------ real catalogs

def test_the_screen_runs_on_every_committed_colour_catalog():
    """No network, no D: — the catalogs are committed Python."""
    data = collect()
    assert len(data) >= 19, f"expected >=19 colour cells, got {len(data)}: {sorted(data)}"
    assert "dog" in data and "rabbit" in data


def test_self_check_against_the_dog_catalog_passes():
    """The classifier is a text heuristic, so it is anchored on the ONE cell with measured ground truth.
    This check EARNED ITS KEEP on the first run: it flagged dog A as UNRECORDED against an expectation of
    STRUCTURAL, and the classifier was RIGHT -- the expectation had encoded the literature, not the
    catalog (the dog ASIP entry never names the SINE)."""
    assert self_check(collect()) == []


def test_the_dog_asip_catalog_gap_is_recorded_rather_than_papered_over():
    """The gap between what the measured artifact knows and what the catalog records is itself a finding:
    even the most-developed colour cell omits one of its five causal variants."""
    assert "dog/A" in _CATALOG_GAPS
    from dna_decode.pigment import dog_coat
    src = getattr(dog_coat.LOCI["A"], "source", "")
    assert "SINE" not in src, "the catalog now records the SINE — retire the gap entry deliberately"


def test_the_two_headline_counts_are_reproducible_from_the_catalogs():
    """Pins the memo's numbers to the code. If a catalog gains a causal variant these MOVE — update the
    memo in the same commit rather than loosening the test."""
    data = collect()
    tot = {}
    for rows in data.values():
        for r in rows:
            tot[r["variant_class"]] = tot.get(r["variant_class"], 0) + 1
    assert sum(tot.values()) == 65, f"loci total moved: {tot}"
    assert tot.get("UNRECORDED") == 40
    assert tot.get("INDEL", 0) + tot.get("STRUCTURAL", 0) == 14


# ------------------------------------------------------------------------ trust-surface corrections

def test_coatcolor_reports_its_measured_result_instead_of_calling_it_pending():
    """UNDER-CLAIM regression. The contract framed a run that HAPPENED (2026-07-30) as 'the v0.1 measured
    tier', and never reported black 0.994. Under-claiming is as much a trust-surface falsehood as
    over-claiming."""
    from dna_decode.cli import TRAITS
    v = TRAITS["coatcolor"]["validation"]
    assert "the v0.1 measured tier" not in v
    assert "0.994" in v and "UNSCORABLE ON THAT SUBSTRATE" in v
    assert (ROOT / "wiki" / "dog_coat_darwins_ark_measured_2026-07-30.md").exists()


def test_the_two_decoder_side_gates_are_in_the_negative_results_map():
    """G1-G8 all gate the LABEL; a curated-catalog cell can fail before a label is even relevant."""
    t = NEG_MAP.read_text(encoding="utf-8", errors="replace")
    assert "| G9 |" in t and "| G10 |" in t
    assert "## The 10 rejection gates" in t, "the gate-count heading must track the table"
    assert "screen it against G1–G10" in t, "the how-to-use line must reference every gate"


def _flat(p: Path) -> str:
    """Markdown prose is hard-wrapped, so a phrase spanning a line break defeats a naive `in` check.
    Collapse whitespace before matching -- reflowing the MEMO to suit a test would be backwards."""
    return " ".join(p.read_text(encoding="utf-8", errors="replace").split())


@pytest.mark.skipif(not SCREEN_MD.exists(), reason="screen memo absent")
def test_the_memo_does_not_claim_unrecorded_loci_are_evidence_about_substrate():
    t = _flat(SCREEN_MD)
    assert "statement about the catalog, never evidence about the substrate" in t
    assert "unvalidatable as written" in t.lower()


# --------------------------------------------------------------- CLI trust-surface strings (Step 4)
# A user reading "deterministic curated OMIA epistatic rule" learns nothing about the fact that when this
# family was measured, only the eumelanin default call survived -- or that 7 cells cannot be validated at
# all as written. These assert against the LIVE derivation, not a copy of the prose.

def _colour_traits() -> dict:
    from dna_decode.cli import TRAITS
    return {k: v for k, v in TRAITS.items() if k.endswith("color") or k == "plumage"}


def test_every_colour_trait_cites_the_screen_artifact():
    missing = [k for k, v in _colour_traits().items()
               if "colour_cell_substrate_screen_2026-08-26" not in v["validation"]]
    assert not missing, f"colour traits with no screen verdict on the trust surface: {sorted(missing)}"


def test_the_seven_unscreenable_traits_say_unvalidatable_as_written():
    """The strongest claim the screen supports, and the one a user most needs to see."""
    unscreenable = {t for t, v in verdicts().items()
                    if v == "UNSCREENABLE_NO_CAUSAL_VARIANTS_RECORDED"}
    assert len(unscreenable) == 7
    traits = _colour_traits()
    for t in unscreenable:
        assert "UNVALIDATABLE AS WRITTEN" in traits[t]["validation"], t
        assert "NO COHORT WOULD HELP" in traits[t]["validation"], t


def test_the_two_clear_traits_are_not_described_as_gated():
    """donkey + roe deer pass G9 and G10. Saying otherwise would overstate the wall."""
    for t in ("donkeycolor", "roedeercolor"):
        v = _colour_traits()[t]["validation"]
        assert "screened CLEAR" in v
        assert "UNVALIDATABLE AS WRITTEN" not in v


def test_each_colour_trait_verdict_matches_the_live_screen():
    """A wording-only change is trivial to fake green, so tie the strings to the derivation: a cell whose
    screen verdict is not UNSCREENABLE must NOT claim to be, and vice versa."""
    live = verdicts()
    for t, v in _colour_traits().items():
        claims_unscreenable = "UNVALIDATABLE AS WRITTEN" in v["validation"]
        is_unscreenable = live[t] == "UNSCREENABLE_NO_CAUSAL_VARIANTS_RECORDED"
        assert claims_unscreenable == is_unscreenable, (
            f"{t}: screen says {live[t]} but its CLI string "
            f"{'claims' if claims_unscreenable else 'does not claim'} unvalidatable")


def test_the_trait_and_species_key_spaces_agree():
    """The screen keys by species, the CLI by trait; `trait_for_species` bridges them. If a new colour
    cell forgets the mapping it surfaces here rather than silently dropping out of the screen."""
    assert {trait_for_species(sp) for sp in collect()} == set(_colour_traits())


# ------------------------------------------------------- the readers `collect` is built on
# `_loci_of` and `_source_of` are private, but they are the ONLY place a catalog's text reaches the
# classifier. A silent shrink there re-classifies real variants as UNRECORDED, which is precisely the
# direction that would OVERSTATE the curation wall this family was frozen for. The headline-count test
# nets that in aggregate; these pin the mechanism, so a regression names its own cause.

def test_loci_of_reads_both_catalog_shapes():
    """Two conventions coexist: the 14 mammal cells expose `.loci` on a MammalCatalog instance, the 5
    bird/dog/cat/horse cells expose a module-level `LOCI` dict. Reading only one silently drops 5 or 14
    cells from the screen."""
    from dna_decode.pigment import dog_coat
    from dna_decode.pigment.mammal_color import MAMMAL_CATALOGS

    assert _loci_of(MAMMAL_CATALOGS["rabbit"]), "the `.loci` attribute shape returned nothing"
    assert _loci_of(dog_coat), "the module-level LOCI shape returned nothing"
    assert set(_loci_of(dog_coat)) >= set(_DOG_TRUTH)


def test_loci_of_on_something_with_neither_shape_is_empty_not_a_crash():
    """A catalog that grows a third convention must drop out visibly (the freeze roster guard fails),
    never take the screen down with it."""
    assert _loci_of(object()) == {}


def test_loci_of_ignores_a_non_dict_loci_attribute():
    class Weird:
        loci = ["A", "B"]
        LOCI = {"K": "real"}

    assert _loci_of(Weird()) == {"K": "real"}


def test_source_of_reads_every_string_field_not_a_hand_listed_few():
    """REGRESSION on the bug the module docstring records: an early reader looked for `notes` while the
    dataclass field is `note` (SINGULAR), and silently returned less text. A hand-listed field tuple
    cannot survive a rename; iterating `dataclasses.fields` can. This asserts the mechanism by giving a
    locus a field NO reader would have thought to list."""
    import dataclasses

    @dataclasses.dataclass
    class Locus:
        gene: str
        note: str
        some_field_added_next_year: str
        not_a_string: int = 7

    src = _source_of(Locus(gene="MLPH", note="c.-22G>A",
                           some_field_added_next_year="frameshift"))
    assert "MLPH" in src and "c.-22G>A" in src
    assert "frameshift" in src, "a newly added string field was dropped — the hand-listed bug is back"
    assert "7" not in src


def test_source_of_skips_empty_strings_so_they_cannot_pad_the_evidence():
    import dataclasses

    @dataclasses.dataclass
    class Locus:
        gene: str = "TYR"
        note: str = ""

    assert _source_of(Locus()) == "TYR"


def test_source_of_falls_back_to_named_attributes_on_a_non_dataclass():
    """Not every future catalog must be a dataclass; the fallback keeps such a locus screenable."""
    class Plain:
        gene = "CBD103"
        source = "c.67_69delGGT"
        note = ""

    src = _source_of(Plain())
    assert "CBD103" in src and "c.67_69delGGT" in src


def test_a_locus_carrying_no_text_at_all_is_unrecorded_end_to_end():
    """The two readers and the classifier must compose to UNRECORDED rather than to a crash or a guess."""
    class Blank:
        pass

    assert classify_variant(_source_of(Blank())) == "UNRECORDED"


# ------------------------------------------------------------------------ collect / verdicts contracts

def test_collect_emits_the_row_shape_its_three_consumers_read():
    """`self_check`, the CLI table, and the JSON artifact all index these keys by name."""
    rows = collect()["dog"]
    for r in rows:
        assert set(r) == {"locus", "gene", "variant_class", "snv_panel_scorable", "source"}
        assert r["variant_class"] in {"SNV", "INDEL", "STRUCTURAL", "UNRECORDED"}
        assert r["snv_panel_scorable"] is snv_panel_scorable(r["variant_class"])
        assert len(r["source"]) <= 220, "the artifact-size truncation moved"


def test_collect_omits_a_catalog_with_no_loci_rather_than_emitting_an_empty_cell():
    """An empty cell would count toward `n_cells` and dilute every headline the freeze rests on."""
    assert all(rows for rows in collect().values())


def test_verdicts_is_pure_given_data_and_does_not_re_read_the_catalogs():
    """THE single derivation every consumer reads. If a future edit made it ignore its argument and always
    call `collect()`, the freeze module and the contract guards would silently stop being testable against
    anything but the live catalogs."""
    out = verdicts({"unicorn": [{"variant_class": "SNV"}],
                    "dog": [{"variant_class": "UNRECORDED"}]})
    assert out == {"unicorncolor": "FULLY_SNV_TRACTABLE",
                   "coatcolor": "UNSCREENABLE_NO_CAUSAL_VARIANTS_RECORDED"}


def test_verdicts_is_keyed_by_trait_not_species():
    """It is joined against `cli.TRAITS` and the cell registry, whose keys are TRAITS. Keying by species
    would make `dog` and `chicken` silently miss every join."""
    live = verdicts()
    assert "coatcolor" in live and "plumage" in live
    assert "dog" not in live and "chicken" not in live


def test_trait_for_species_names_the_two_exceptions_and_defaults_for_the_rest():
    assert trait_for_species("dog") == "coatcolor"
    assert trait_for_species("chicken") == "plumage"
    assert trait_for_species("rabbit") == "rabbitcolor"


def test_snv_panel_scorable_on_an_unknown_class_is_none_not_false():
    """Same tri-state discipline as UNRECORDED: an unrecognised class is unknowable, and calling it False
    would assert a substrate limit nobody derived."""
    assert snv_panel_scorable("SOMETHING_NEW") is None


# ------------------------------------------------------------------- the self-check anchor can FAIL
# `test_self_check_against_the_dog_catalog_passes` proves the anchor is satisfied. These prove it is not
# VACUOUS -- an anchor for a text heuristic that could never fail would certify the heuristic against
# nothing, the same defect `test_the_freeze_guard_is_not_vacuous` exists to rule out for the roster guard.

def test_the_self_check_reports_a_mismatch_rather_than_agreeing_with_itself():
    wrong = [{"locus": k, "variant_class": "SNV"} for k in _DOG_TRUTH]
    fails = self_check({"dog": wrong})
    assert fails, "the anchor accepted a classifier that got the dog catalog wrong"
    assert any("classifier said SNV" in f and "records INDEL" in f for f in fails), fails


def test_the_self_check_reports_a_missing_locus_distinctly_from_a_wrong_one():
    """A locus vanishing from the catalog is a different failure from a mis-classified one, and the
    message must say which -- a screen that reads zero loci would otherwise look merely 'mismatched'."""
    fails = self_check({})
    assert len(fails) == len(_DOG_TRUTH)
    assert all("absent from the screen" in f for f in fails), fails


def test_the_self_check_anchors_on_a_mixture_of_classes():
    """Anchoring on five loci that were all SNV would not exercise the ordering that makes the classifier
    correct. The dog cell is the anchor precisely because it spans SNV, INDEL and UNRECORDED."""
    assert set(_DOG_TRUTH.values()) >= {"SNV", "INDEL", "UNRECORDED"}


# ------------------------------------------------------------------------ the script CLI (argparse)
# `--self-check` was converted from `"--self-check" in sys.argv` to argparse so the repo's
# `test_every_flag_named_in_the_docs_is_declared_somewhere` guard (a static scan for `add_argument`) can
# SEE it. That guard proves the flag is declared; these prove it is also WIRED -- per the repo rule that
# an advertised command is a promise, validated through its real parser.

def test_self_check_flag_runs_the_anchor_and_writes_no_artifact(tmp_path, monkeypatch, capsys):
    untouched = tmp_path / "must-not-be-written.json"
    monkeypatch.setattr(screen_cli, "OUT", untouched)
    assert screen_cli.main(["--self-check"]) == 0
    assert not untouched.exists(), "--self-check wrote the artifact; it is meant to check and exit"
    assert "self-check: PASS" in capsys.readouterr().out


def test_the_flag_is_parsed_by_argparse_not_by_inspecting_the_real_sys_argv(tmp_path, monkeypatch):
    """The conversion's actual behaviour change: argv is an ARGUMENT. Under the old `in sys.argv` form
    this test's own pytest command line would leak into the decision."""
    monkeypatch.setattr(screen_cli, "OUT", tmp_path / "out.json")
    monkeypatch.setattr(sys, "argv", ["prog", "--self-check"])
    assert screen_cli.main([]) == 0
    assert (tmp_path / "out.json").exists(), "an empty argv must take the artifact path, not the flag path"


def test_an_unknown_flag_is_rejected_rather_than_silently_ignored():
    with pytest.raises(SystemExit) as e:
        screen_cli.main(["--no-such-flag"])
    assert e.value.code == 2


def test_the_artifact_totals_agree_with_the_pinned_headline_counts(tmp_path, monkeypatch):
    """The script recomputes the totals by counting rows, INDEPENDENTLY of `summarise` — which is what
    `EXPECTED_TOTALS` is pinned from. Two computations of one headline can disagree; this is the only
    place they meet."""
    import json

    from dna_decode.data.colour_cell_freeze import EXPECTED_TOTALS

    out = tmp_path / "screen.json"
    monkeypatch.setattr(screen_cli, "OUT", out)
    assert screen_cli.main([]) == 0

    report = json.loads(out.read_text(encoding="utf-8"))
    assert report["_schema"] == "colour-cell-substrate-screen-v1"
    assert report["self_check_failures"] == []
    for k, v in EXPECTED_TOTALS.items():
        assert report["totals"][k] == v, f"{k}: artifact says {report['totals'][k]}, pinned {v}"


def test_the_artifact_records_its_scope_caveat_and_the_catalog_gap(tmp_path, monkeypatch):
    """The caveat is the artifact's load-bearing honesty: UNRECORDED is a finding about the CATALOG, not
    evidence about the substrate. A machine-readable report that dropped it would be read as the stronger
    claim."""
    import json

    out = tmp_path / "screen.json"
    monkeypatch.setattr(screen_cli, "OUT", out)
    screen_cli.main([])
    report = json.loads(out.read_text(encoding="utf-8"))
    assert "NOT evidence about the" in report["honest_scope"]
    assert "dog/A" in report["catalog_gaps_vs_measured_artifact"]
    assert report["ground_truth"].endswith("dog_coat_darwins_ark_measured_2026-07-30.md")


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-q"]))
