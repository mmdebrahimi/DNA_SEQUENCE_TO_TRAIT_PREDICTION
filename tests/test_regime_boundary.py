"""Pin the g->p regime boundary, including the SCOPE that keeps getting lost.

The boundary has been mis-stated three times, always by compressing a scoped negative into a general
one. These tests make each of those three compressions fail loudly.

Offline, pure, no fixtures.
"""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

import pytest  # noqa: E402

from dna_decode.eval.regime import (CLOSED_NEGATIVE, LOSES_TO_CATALOG, OPEN,  # noqa: E402
                                    ORGANISM_TRANSFER_LEVELS, ORGANISM_TRANSFER_UNMEASURED, REGIMES,
                                    REQUIRES_DECONFOUNDING, SPLIT_UNITS, WORKS, Regime,
                                    _transfer_rank, classify_regime,
                                    organism_transfer_is_unmeasured, screen_proposal,
                                    split_units_for)


# --- the three compressions that have actually happened ---

def test_organism_level_gp_is_NOT_a_closed_negative_in_general():
    """COMPRESSION 1, made three times. A clean organism-level positive exists (yeast segregant cross,
    12/12 traits at r 0.46-0.80), so a blanket 'organism g->p is closed' is false."""
    assert screen_proposal("constructed", "organism", "supervised").verdict == WORKS
    assert not screen_proposal("constructed", "organism", "supervised").refused


def test_the_negative_is_zero_shot_scoped_not_learning_scoped():
    """COMPRESSION 2. The 0-for-5 is ZERO-SHOT-only; the shipped architecture pairs the deterministic
    catalogue with a SUPERVISED complement. Refusing a supervised proposal would re-commit the error."""
    zero = screen_proposal("natural", "organism", "zero_shot")
    sup = screen_proposal("natural", "organism", "supervised")
    assert zero.verdict == CLOSED_NEGATIVE and zero.refused
    assert sup.verdict == REQUIRES_DECONFOUNDING and not sup.refused, (
        "a supervised natural-population proposal must get CONDITIONS, not a refusal")
    assert any("WITHIN-GROUP" in c for c in sup.conditions)


def test_the_discriminating_variable_is_population_design_not_organism_complexity():
    """COMPRESSION 3. Holding endpoint and method fixed, flipping ONLY the population design flips the
    verdict — which is what 'the discriminator is population design' means operationally."""
    natural = screen_proposal("natural", "organism", "zero_shot").verdict
    constructed = screen_proposal("constructed", "organism", "supervised").verdict
    assert natural == CLOSED_NEGATIVE and constructed == WORKS


# --- the refusal is narrow ---

def test_only_the_measured_dead_regime_is_refused():
    """A boundary that refused broadly would block real work. Exactly one regime refuses."""
    refused = [r for r in REGIMES
               if screen_proposal(r.population, r.endpoint, r.method).refused]
    assert len(refused) == 1, f"expected exactly one refusing regime, got {[r.key for r in refused]}"
    assert refused[0].key == "natural_organism_zeroshot"


def test_an_unscreened_combination_is_open_not_endorsed():
    """Absence of a measurement must not read as promise."""
    res = screen_proposal("constructed", "organism_condition_switch", "zero_shot")
    assert res.verdict == OPEN and not res.refused
    assert "not the same as promising" in res.reason


def test_the_condition_switch_cell_is_open_not_solved():
    """The one genuinely open cell — it must not be reported as working."""
    res = screen_proposal("constructed", "organism_condition_switch", "supervised")
    assert res.verdict == OPEN
    assert "silent" in res.evidence.lower()


# --- the catalog rule, and its inversion ---

def test_a_curated_catalog_beats_a_learned_scorer_wherever_one_exists():
    res = screen_proposal("constructed", "molecular", "supervised", curated_catalog_exists=True)
    assert res.verdict == LOSES_TO_CATALOG
    assert "0.454" in res.evidence, "the BELOW-chance measurement must ship with the verdict"


def test_the_catalog_rule_does_not_fire_on_the_catalog_itself():
    """A deterministic catalog proposal must not be told it loses to a catalog."""
    res = screen_proposal("natural", "molecular", "deterministic_catalog", curated_catalog_exists=True)
    assert res.verdict != LOSES_TO_CATALOG


def test_without_a_catalog_the_working_molecular_regime_still_works():
    assert screen_proposal("constructed", "molecular", "supervised").verdict == WORKS


# --- hygiene ---

def test_scale_is_never_offered_as_the_remedy_for_the_closed_negative():
    """The specific wrong next step. It must be named as wrong in the output, not just omitted."""
    res = screen_proposal("natural", "organism", "zero_shot")
    assert any("do NOT re-run this at larger scale" in c for c in res.conditions)


def test_every_regime_cites_an_artifact_that_exists():
    """A claim whose evidence file is missing is not a claim. Same rule as the doc-citation guard."""
    missing = [r.key for r in REGIMES if not (ROOT / r.artifact).exists()]
    assert not missing, f"regimes cite artifacts that do not exist: {missing}"


def test_bad_axis_values_are_refused_not_guessed():
    for bad in (("eukaryote", "organism", "zero_shot"),
                ("natural", "phenotype", "zero_shot"),
                ("natural", "organism", "finetuned")):
        assert screen_proposal(*bad).verdict == "UNKNOWN"


def test_classify_regime_is_exact_not_fuzzy():
    assert classify_regime("NATURAL", "Organism", "ZERO_SHOT").key == "natural_organism_zeroshot"
    assert classify_regime("natural", "organism", "deterministic_catalog") is None


def test_the_catalog_result_is_about_REPLACING_a_catalog_not_complementing_it():
    """MY OWN OVER-COMPRESSION, inside the module that exists to prevent it (fixed 2026-09-30).

    `screen_proposal` short-circuited to LOSES_TO_CATALOG for ANY learned method whenever a catalog
    existed. That is the measured truth about REPLACING a catalog (ESM 0.454 vs catalog 0.926) -- but it
    was also being returned for a blind-spot COMPLEMENT, which is the one shape now measured to WORK
    (supervised genotype tokens, 0.8142 leave-one-study-out on the catalog-NEGATIVE subset).
    """
    replace = screen_proposal("natural", "molecular", "supervised",
                              curated_catalog_exists=True, target="replace")
    assert replace.verdict == LOSES_TO_CATALOG, "replacing a catalog is still measured to lose"

    complement = screen_proposal("natural", "molecular", "supervised",
                                 curated_catalog_exists=True, target="blind_spot_complement")
    assert complement.verdict == WORKS
    assert complement.regime == "natural_molecular_supervised_blindspot"
    assert "blind" in complement.evidence.lower() or "blind" in complement.reason.lower()


def test_a_ZERO_SHOT_complement_still_loses_because_that_is_what_was_measured():
    """The scoping must not become a blanket exemption for anything labelled 'complement'. Zero-shot ESM2
    was measured ON the blind spot and FAILED there (0.4485 against a >=0.65 bar), so it must still lose
    -- reached via the regime match rather than the short-circuit, but the verdict is the same."""
    r = screen_proposal("natural", "molecular", "zero_shot",
                        curated_catalog_exists=True, target="blind_spot_complement")
    assert r.verdict == LOSES_TO_CATALOG


def test_an_unknown_target_is_refused_rather_than_silently_treated_as_replace():
    """A typo'd target must not quietly inherit the refusing default."""
    r = screen_proposal("natural", "molecular", "supervised",
                        curated_catalog_exists=True, target="complement")   # not the real token
    assert r.verdict == "UNKNOWN" and "target" in r.reason


def test_the_supervised_blindspot_row_cites_an_artifact_that_exists():
    """Same rule the map already enforces elsewhere: a regime cannot certify itself from prose."""
    from pathlib import Path
    r = classify_regime("natural", "molecular", "supervised")
    assert r is not None
    assert (Path(__file__).resolve().parent.parent / r.artifact).exists(), \
        f"cited artifact missing: {r.artifact}"


# --- the transfer axis (F1 Step 1): WHICH SPLIT produced each number -----------------------------

def test_every_regime_declares_a_split_unit_and_it_is_REQUIRED():
    """A row that could omit its split would ship silently unmeasured -- the failure the axis exists
    to prevent. `split_unit` is keyword-only AND has no default, so omitting it is a TypeError."""
    for r in REGIMES:
        assert r.split_unit <= SPLIT_UNITS, f"{r.key} declares a split unit outside the vocabulary"
        assert r.split_unit, f"{r.key} declares an EMPTY split_unit; use {{'none'}} to say so explicitly"

    with pytest.raises(TypeError):
        Regime("x", "natural", "organism", "zero_shot", OPEN, "ev", "wiki/x.md")   # no split_unit


def test_organism_transfer_defaults_to_unmeasured_which_is_both_failclosed_and_true():
    r = Regime("x", "natural", "organism", "zero_shot", OPEN, "ev", "wiki/x.md",
               split_unit=frozenset({"none"}))
    assert r.organism_transfer == ORGANISM_TRANSFER_UNMEASURED
    for row in REGIMES:
        assert row.organism_transfer in ORGANISM_TRANSFER_LEVELS


def test_split_units_are_UNORDERED_so_the_single_ladder_cannot_reappear():
    """THE DESIGN POINT. The first draft ranked position < study < protein < organism < clade on one
    scale, which would let this module announce that a cross-protein split beats a leave-one-study
    split -- a relation nobody measured. Ranking is defined ONLY over the narrow organism-transfer
    field, so asking it to rank a split unit must fail."""
    assert _transfer_rank("unmeasured") < _transfer_rank("held_out_organism") < \
        _transfer_rank("held_out_clade")
    for incomparable in ("position", "study", "protein", "condition", "none"):
        assert incomparable in SPLIT_UNITS
        with pytest.raises(ValueError):
            _transfer_rank(incomparable)


def test_organism_transfer_is_unmeasured_is_a_TRIPWIRE_not_a_frozen_count():
    """UPDATED DELIBERATELY 2026-10-06, which is what this tripwire asked for.

    It previously asserted that EVERY row was unmeasured, and it FIRED when
    `curated_catalog_cross_organism` landed with held-out-organism evidence from the transfer ladder --
    the event its own message named as "the headline F1 was built to surface". The relation is now
    "exactly the rows that are measured are measured": the carrier set is named explicitly, so the NEXT
    regime to acquire this evidence still trips the wire rather than slipping in.
    """
    unmeasured = organism_transfer_is_unmeasured()
    carriers = {r.key for r in REGIMES if r.organism_transfer != "unmeasured"}
    assert carriers == {"curated_catalog_cross_organism"}, (
        "the set of regimes carrying held-out-ORGANISM evidence changed: %s. That is the headline F1 was "
        "built to surface -- update the plan's claim and this test deliberately, together." % carriers)
    assert set(unmeasured) == {r.key for r in REGIMES} - carriers
    assert len(unmeasured) == 8
    assert split_units_for("natural_molecular_supervised_blindspot") == frozenset({"study"})
    with pytest.raises(KeyError):
        split_units_for("no_such_regime")


def test_the_transfer_claim_is_AUGMENT_ONLY_verified_by_diff():
    """Proved, not promised: across every regime x target x catalog state, claiming organism transfer
    may only APPEND conditions. Non-vacuous -- at least one case must actually gain one."""
    gained = 0
    for r in REGIMES:
        for target in ("replace", "blind_spot_complement"):
            for cat in (False, True):
                base = screen_proposal(r.population, r.endpoint, r.method,
                                       curated_catalog_exists=cat, target=target).as_dict()
                claimed = screen_proposal(r.population, r.endpoint, r.method,
                                          curated_catalog_exists=cat, target=target,
                                          claims_organism_transfer="held_out_clade").as_dict()
                for k in base:
                    if k == "conditions":
                        if base[k] != claimed[k]:
                            gained += 1
                            assert claimed[k][:len(base[k])] == base[k], \
                                "conditions must be APPENDED, never replaced"
                        continue
                    assert base[k] == claimed[k], f"{k} changed for {r.key} -- not augment-only"
    assert gained > 0, "no case gained a condition; the guard would pass vacuously"


def test_a_supported_transfer_claim_adds_nothing():
    """The condition fires on an UNSUPPORTED claim only. Claiming `unmeasured` matches every row
    today, so it must be silent -- otherwise the layer would nag on a truthful claim."""
    r = REGIMES[0]
    a = screen_proposal(r.population, r.endpoint, r.method).as_dict()
    b = screen_proposal(r.population, r.endpoint, r.method,
                        claims_organism_transfer=ORGANISM_TRANSFER_UNMEASURED).as_dict()
    assert a == b


def test_a_bad_organism_transfer_claim_is_refused_not_guessed():
    res = screen_proposal("natural", "organism", "zero_shot",
                          claims_organism_transfer="held_out_planet")
    assert res.verdict == "UNKNOWN" and "claims_organism_transfer" in res.reason


def test_a_bad_transfer_claim_is_refused_on_the_CATALOG_path_TOO_not_just_the_regime_path():
    """The ordering the source comments on ('validated BEFORE the catalog short-circuit so a bad value
    cannot slip through on either path') but which the test above cannot reach -- it screens with no
    catalog, so it only ever exercises the regime-match path.

    The catalog branch returns early WITH a conditions list built from the same claim, so if the
    validation ever moved below it a typo'd level would reach `_transfer_rank` and raise ValueError, or
    worse be rank-compared as if it were a real level. Non-vacuous: the good value on the identical
    call is NOT refused, so this is testing the validation and not the branch.
    """
    bad = screen_proposal("natural", "molecular", "supervised", curated_catalog_exists=True,
                          target="replace", claims_organism_transfer="held_out_planet")
    assert bad.verdict == "UNKNOWN" and "claims_organism_transfer" in bad.reason
    assert bad.regime is None, "a refused input must not be attributed to a regime"

    good = screen_proposal("natural", "molecular", "supervised", curated_catalog_exists=True,
                           target="replace", claims_organism_transfer="held_out_clade")
    assert good.verdict == LOSES_TO_CATALOG
    assert any("NOT supported by the cited evidence" in c for c in good.conditions)


def test_an_UNSCREENED_combination_gains_no_transfer_condition_because_there_is_nothing_to_compare():
    """Scope pin. `_transfer_gap_conditions` needs a matched regime to compare the claim against, so
    the OPEN fall-through (no regime at all) appends nothing -- a reader must not assume the layer
    always speaks up. The absence is correct here: nothing was claimed ABOUT anything."""
    res = screen_proposal("constructed", "organism_condition_switch", "zero_shot",
                          claims_organism_transfer="held_out_clade")
    assert res.verdict == OPEN and res.regime is None
    assert not any("organism_transfer" in c for c in res.conditions), res.conditions
    assert res.conditions == ["measure a de-confounded baseline first"]


def test_as_dict_carries_both_new_fields_json_serialisably():
    """scripts/regime_map.py json.dumps() this; a raw frozenset would raise at write time."""
    import json
    d = REGIMES[0].as_dict()
    assert isinstance(d["split_unit"], list) and d["organism_transfer"] in ORGANISM_TRANSFER_LEVELS
    json.dumps(d)


def test_the_map_script_RENDERS_the_new_fields_and_the_headline(monkeypatch, tmp_path, capsys):
    """A field carried only into learned_regime_map.json is NOT a disclosure -- this repo has already
    paid for that twice. WIKI is redirected to tmp_path so the test cannot dirty the tracked artifact
    (the date-stamped-rewrite trap recorded against the pgx report-card test)."""
    import scripts.regime_map as rm

    monkeypatch.setattr(rm, "WIKI", tmp_path)
    # main() parses sys.argv, which under pytest carries pytest's own flags -> argparse SystemExit(2).
    monkeypatch.setattr(sys, "argv", ["regime_map.py"])
    rc = rm.main()
    out = capsys.readouterr().out

    assert rc == 0, "all-unmeasured is a truthful state, not a tool failure"
    assert "split unit" in out and "org transfer" in out
    assert f"{len(organism_transfer_is_unmeasured())} of {len(REGIMES)} regimes carry NO " \
        "held-out-ORGANISM transfer evidence." in out
    assert "UNORDERED" in out
    assert (tmp_path / "learned_regime_map.json").exists()


def _run_map(monkeypatch, tmp_path, argv):
    """Drive `scripts/regime_map.py` through its real argparse. WIKI is redirected so no tracked
    artifact can be rewritten, and sys.argv is replaced because main() parses it directly (under
    pytest it otherwise carries pytest's own flags -> SystemExit(2))."""
    import scripts.regime_map as rm
    monkeypatch.setattr(rm, "WIKI", tmp_path)
    monkeypatch.setattr(sys, "argv", ["regime_map.py", *argv])
    return rm.main()


def test_the_screen_subcommand_exits_1_on_a_refusal_and_0_on_a_live_regime(monkeypatch, tmp_path,
                                                                          capsys):
    """The whole `--screen` branch was unexercised: it is the invocation a human actually types, it
    returns a DIFFERENT exit code from the map path, and it writes no artifact."""
    import json

    rc = _run_map(monkeypatch, tmp_path, ["--screen", "natural", "organism", "zero_shot"])
    out = capsys.readouterr().out
    assert rc == 1, "the one measured dead regime must exit non-zero when screened"
    assert "REFUSED" in out and CLOSED_NEGATIVE in out
    payload = json.loads(out.split("\n\n")[0])
    assert payload["regime"] == "natural_organism_zeroshot" and payload["refused"] is True
    assert not list(tmp_path.iterdir()), "--screen must not write the map artifact"

    rc_ok = _run_map(monkeypatch, tmp_path, ["--screen", "constructed", "molecular", "supervised"])
    assert rc_ok == 0 and "not refused" in capsys.readouterr().out


def test_screening_a_TYPO_exits_2_distinctly_from_a_clean_pass_and_a_refusal(monkeypatch, tmp_path,
                                                                             capsys):
    """FIXED 2026-10-03 (this test previously recorded the gap). `--screen` returned
    `1 if res.refused else 0`, and refused is true only for CLOSED_NEGATIVE -- so an invalid axis
    yielded UNKNOWN and exit 0, the same code as a successfully screened WORKS regime, and
    `regime_map.py --screen ... && proceed` proceeded on a typo. UNKNOWN now exits 2, leaving all
    three outcomes distinguishable."""
    rc = _run_map(monkeypatch, tmp_path, ["--screen", "natural", "organism", "finetuned"])
    out = capsys.readouterr().out
    assert rc == 2, "an invalid axis must not share an exit code with a clean screen"
    assert '"verdict": "UNKNOWN"' in out

    assert _run_map(monkeypatch, tmp_path, ["--screen", "constructed", "molecular",
                                            "supervised"]) == 0
    capsys.readouterr()
    assert _run_map(monkeypatch, tmp_path, ["--screen", "natural", "organism", "zero_shot"]) == 1
    capsys.readouterr()


def test_the_map_REFUSES_to_certify_itself_when_a_cited_artifact_is_missing(monkeypatch, tmp_path,
                                                                           capsys):
    """The script's ONLY legitimate non-zero exit on the map path, and it was untested. A regime whose
    evidence file does not resolve is a memory, not a regime -- the failure mode the module exists to
    prevent -- so it must break the exit code even though the json is still written."""
    import json

    import dna_decode.eval.regime as reg

    phantom = Regime("phantom", "natural", "organism", "supervised", OPEN, "ev",
                     "wiki/this_artifact_does_not_exist.md", split_unit=frozenset({"none"}))
    monkeypatch.setattr(reg, "REGIMES", REGIMES + (phantom,))
    rc = _run_map(monkeypatch, tmp_path, [])
    out = capsys.readouterr().out

    assert rc == 1, "a missing cited artifact must not exit 0"
    assert "ARTIFACT MISSING" in out and "REFUSING to certify" in out
    written = json.loads((tmp_path / "learned_regime_map.json").read_text(encoding="utf-8"))
    assert written["artifacts_missing"] == ["phantom"], \
        "the artifact must RECORD the break, not merely exit on it"


def test_the_missing_artifact_guard_is_non_vacuous_on_the_real_map(monkeypatch, tmp_path, capsys):
    """Non-vacuity for the test above: with the real REGIMES every artifact resolves, so the same code
    path reports nothing missing and exits 0. If this ever fails, a regime's evidence has gone."""
    import json

    rc = _run_map(monkeypatch, tmp_path, [])
    assert rc == 0
    assert "ARTIFACT MISSING" not in capsys.readouterr().out
    written = json.loads((tmp_path / "learned_regime_map.json").read_text(encoding="utf-8"))
    assert written["artifacts_missing"] == []
    assert len(written["regimes"]) == len(REGIMES)
    assert all(r["artifact_exists"] for r in written["regimes"])


# ===================================================================================================
# Transfer-ladder Step 9 (2026-10-06): the first row to carry held-out-ORGANISM evidence, and the
# retirement of a restated fact that became false the moment it landed.
# ===================================================================================================
def test_the_unmeasured_set_SHRANK_by_exactly_one_and_is_not_empty():
    """The event `organism_transfer_is_unmeasured`'s own docstring names as the one worth noticing.

    Asserting 'dropped by one' rather than 'is empty' matters: the other 8 regimes are GENUINELY
    unmeasured and an empty return would be the over-claim this field exists to prevent.
    """
    from dna_decode.eval.regime import REGIMES, organism_transfer_is_unmeasured

    unmeasured = organism_transfer_is_unmeasured()
    assert len(REGIMES) == 9
    assert len(unmeasured) == 8, unmeasured
    assert "curated_catalog_cross_organism" not in unmeasured


def test_the_cross_organism_row_is_shaped_as_held_out_organism():
    from dna_decode.eval.regime import REGIMES

    r = next(x for x in REGIMES if x.key == "curated_catalog_cross_organism")
    assert r.organism_transfer == "held_out_organism"
    assert r.split_unit == frozenset({"organism"})
    assert r.population == "constructed" and r.endpoint == "organism"


def test_the_new_rows_verdict_is_OPEN_because_the_LADDER_was_indeterminate():
    """Set FROM the measurement, not chosen: the cross-organism number exists and clears its null ~4x
    (so not CLOSED_NEGATIVE), but the ladder returned INDETERMINATE_INSUFFICIENT_RUNGS (so not WORKS)."""
    from dna_decode.eval.regime import OPEN, REGIMES

    r = next(x for x in REGIMES if x.key == "curated_catalog_cross_organism")
    assert r.verdict == OPEN
    assert "INDETERMINATE_INSUFFICIENT_RUNGS" in r.evidence
    assert "0.1633" in r.evidence and "0.0446" in r.evidence
    assert "coverage_lift, NOT AUROC" in r.note, "the note must steer a reader off AUROC"


def test_the_new_row_warns_against_reading_auroc_and_names_the_phrasing_gap():
    from dna_decode.eval.regime import REGIMES

    r = next(x for x in REGIMES if x.key == "curated_catalog_cross_organism")
    assert "tie-mass" in r.note
    assert "PHRASING gap" in r.note
    assert "EQUIDISTANT" in r.note, "the eukaryote tie must travel with the row"


def test_the_gap_condition_is_DERIVED_and_names_the_carrier():
    """It used to end with the literal 'no regime in this map carries one today', which became FALSE the
    moment a row carried that evidence. A restated fact goes stale silently."""
    from dna_decode.eval import regime as rg

    conds = rg._transfer_gap_conditions(rg.REGIMES[0], "held_out_organism")
    joined = " ".join(conds)
    assert "no regime in this map carries one today" not in joined
    assert "curated_catalog_cross_organism" in joined
    assert "regime(s) that DO carry held-out-organism evidence" in joined


def test_the_fallback_sentence_still_exists_for_an_empty_table():
    """The honest wording must survive for the case where nothing carries the evidence -- otherwise the
    derivation would simply have replaced one unconditional claim with another."""
    from dna_decode.eval import regime as rg

    src = (rg.__file__ and open(rg.__file__, encoding="utf-8").read()) or ""
    assert "no regime in\n" in src or "no regime in " in src, \
        "the empty-table branch must still be reachable in source"


def test_a_claim_of_held_out_CLADE_still_gets_a_gap_condition_against_the_new_row():
    """The ladder gives held_out_ORGANISM, which is strictly below held_out_CLADE -- the ordering must
    still bite rather than being satisfied by the new row."""
    from dna_decode.eval import regime as rg

    r = next(x for x in rg.REGIMES if x.key == "curated_catalog_cross_organism")
    conds = rg._transfer_gap_conditions(r, "held_out_clade")
    assert conds, "a clade claim against an organism-level row must still be flagged"
    assert "held_out_organism" in " ".join(conds)
    # and a claim it DOES support produces nothing
    assert rg._transfer_gap_conditions(r, "held_out_organism") == []


def test_the_cited_artifact_exists_so_regime_map_can_certify_it():
    from pathlib import Path

    from dna_decode.eval.regime import REGIMES

    root = Path(__file__).resolve().parents[1]
    r = next(x for x in REGIMES if x.key == "curated_catalog_cross_organism")
    assert (root / r.artifact).exists(), r.artifact
