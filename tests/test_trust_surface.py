"""Pins the inline trust-surface: every decoder cell resolves to its HONEST validation tier, no tier is
fabricated, and the independence flag is consistent. Card-dependent tiers skip if a card JSON is absent;
the structural invariants (always-a-dict / no-fabrication / genus-normalization / UNKNOWN) always run.
"""
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
from dna_decode.data import trust_surface as ts  # noqa: E402

_WIKI = REPO / "wiki"
_REQUIRED_KEYS = {"tier", "independent", "headline", "metric", "n", "cell", "source_card", "caveat"}


def _has(card: str) -> bool:
    return (_WIKI / card).exists()


# --- structural invariants (always run) ---

def test_always_returns_full_dict():
    for drug, org in [("ciprofloxacin", "Escherichia"), ("zzz_nonsense", None), ("rifampicin", None)]:
        b = ts.trust_block(drug, org)
        assert _REQUIRED_KEYS <= set(b), f"missing keys for {drug}"
        assert isinstance(b["caveat"], str) and b["caveat"]


def test_unknown_drug_is_unknown_and_never_fabricates():
    b = ts.lookup_trust("not_a_real_drug_xyz", "Nowhere")
    assert b["tier"] == ts.UNKNOWN
    assert b["metric"] is None and b["independent"] is False


def test_independent_flag_is_consistent():
    # only the two independent tiers may carry independent=True
    for drug, org in [("efavirenz", None), ("rifampicin", None), ("ciprofloxacin", "Escherichia"),
                      ("nirmatrelvir", None), ("fluconazole", "Candida_auris"), ("oxacillin", "Staphylococcus_aureus")]:
        b = ts.trust_block(drug, org)
        assert b["independent"] == (b["tier"] in (ts.INDEPENDENT_WETLAB, ts.INDEPENDENT_MEASURED))


def test_genus_normalization_collapses_organism_spellings():
    a = ts.trust_block("ciprofloxacin", "Escherichia")
    b = ts.trust_block("ciprofloxacin", "Escherichia_coli_Shigella")
    assert a["tier"] == b["tier"] and a["cell"] == b["cell"]


def test_one_line_is_ascii_safe():
    s = ts.one_line(ts.trust_block("rifampicin"))
    assert isinstance(s, str) and "validation:" in s
    assert "—" not in s  # no em-dash (cp1252-console trap)


# --- per-tier pins (skip if the backing card is absent) ---

@pytest.mark.skipif(not _has("hiv_decoder_report_card.json"), reason="hiv card absent")
def test_hiv_is_free_wetlab_independent():
    b = ts.trust_block("efavirenz")
    assert b["tier"] == ts.INDEPENDENT_WETLAB and b["independent"] is True
    assert b["metric"] is not None and "wetlab" in b["caveat"].lower() or "wet-lab" in b["caveat"].lower()


@pytest.mark.skipif(not _has("tb_report_card.json"), reason="tb card absent")
def test_tb_is_independent_measured():
    b = ts.trust_block("rifampicin", "Mycobacterium_tuberculosis")
    assert b["tier"] == ts.INDEPENDENT_MEASURED and b["independent"] is True
    assert b["source_card"].endswith("tb_report_card.md")


@pytest.mark.skipif(not _has("amr_portal_independent_report_card.json"), reason="amr portal card absent")
def test_ecoli_cipro_is_independent_measured():
    b = ts.trust_block("ciprofloxacin", "Escherichia")
    assert b["tier"] == ts.INDEPENDENT_MEASURED and b["independent"] is True
    assert b["metric"] is not None


def test_sarscov2_is_in_distribution():
    b = ts.trust_block("nirmatrelvir")
    assert b["tier"] == ts.IN_DISTRIBUTION and b["independent"] is False


def test_fungal_is_no_free_source():
    b = ts.trust_block("fluconazole", "Candida_auris")
    assert b["tier"] == ts.NO_FREE_PHENOTYPE_SOURCE and b["independent"] is False


def test_oxacillin_is_label_confounded():
    b = ts.trust_block("oxacillin", "Staphylococcus_aureus")
    assert b["tier"] == ts.LABEL_CONFOUNDED


def test_meropenem_acinetobacter_abstains():
    b = ts.trust_block("meropenem", "Acinetobacter")
    assert b["tier"] == ts.ABSTAINS_BY_DESIGN


if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_") and callable(v)]
    for fn in fns:
        try:
            fn(); print(f"PASS {fn.__name__}")
        except Exception as e:  # pragma: no cover
            print(f"FAIL {fn.__name__}: {e}")
    print(f"\n{len(fns)} tests")


# --- prospective-regression annotation (2026-08-24) -------------------------------------------------
# A cell's badge tier resolves from whichever card has the best evidence, so a prospective regression
# recorded on the validation card would never reach E. coli x gentamicin -- that cell resolves at the
# AMR-Portal card and quotes acc 0.987, while its own post-lock cohort measured sens 0.429. The
# annotation is therefore CROSS-CUTTING (applied in trust_block, after any card wins).


def _fake_card(tmp_path, monkeypatch, *, regression=True):
    import json as _json
    import dna_decode.data.trust_surface as ts
    cells = [{"organism": "escherichia_coli_shigella", "drug": "gentamicin", "state": "SCORED",
              "acc": 0.95, "n": 60,
              "prospective": {"status": "scored", "sens": 0.429, "n_scored": 62,
                              "lock_date": "2026-06-13", "regression": regression}},
             {"organism": "escherichia_coli_shigella", "drug": "ciprofloxacin", "state": "SCORED",
              "acc": 0.95, "n": 60,
              "prospective": {"status": "scored", "sens": 0.917, "n_scored": 61,
                              "lock_date": "2026-06-13", "regression": False}}]
    if regression:
        cells[0]["prospective_regression"] = True
        cells[0]["deployment_caveat"] = "prospective sens 0.429 -- the frozen rule under-calls this cell"
    # distinct filename per variant + an explicit cache clear: `_load` is @lru_cache'd, so rewriting the
    # same path silently returns the FIRST variant and the fixture would test nothing.
    card = tmp_path / f"card_{'reg' if regression else 'clean'}.json"
    card.write_text(_json.dumps({"cells": cells}), encoding="utf-8")
    ts._load.cache_clear()
    monkeypatch.setattr(ts, "_card_path",
                        lambda name: card if name == "decoder_validation_report_card.json"
                        else tmp_path / "absent.json")
    return ts


def test_prospective_regression_is_attached_without_changing_tier_or_metric(tmp_path, monkeypatch):
    ts = _fake_card(tmp_path, monkeypatch)
    raw = ts.lookup_trust("gentamicin", "Escherichia_coli_Shigella")
    badge = ts.trust_block("gentamicin", "Escherichia_coli_Shigella")
    assert badge["tier"] == raw["tier"] and badge["metric"] == raw["metric"]   # augment, never demote
    assert badge["prospective_regression"]["sens"] == 0.429
    assert "PROSPECTIVE REGRESSION" in badge["caveat"] and "UNDER-CALLS" in badge["caveat"]
    assert "PROSPECTIVE REGRESSION" in ts.one_line(badge)


def test_a_healthy_cell_gets_no_annotation(tmp_path, monkeypatch):
    """Non-vacuity: the note must be selective, or it means nothing."""
    ts = _fake_card(tmp_path, monkeypatch)
    badge = ts.trust_block("ciprofloxacin", "Escherichia_coli_Shigella")
    assert badge.get("prospective_regression") is None
    assert "PROSPECTIVE REGRESSION" not in badge["caveat"]

    ts2 = _fake_card(tmp_path, monkeypatch, regression=False)
    assert ts2.trust_block("gentamicin", "Escherichia_coli_Shigella").get("prospective_regression") is None


def test_prospective_regressions_listing_and_missing_card_are_safe(tmp_path, monkeypatch):
    ts = _fake_card(tmp_path, monkeypatch)
    assert [(r["organism"], r["drug"]) for r in ts.prospective_regressions()] == [
        ("escherichia_coli_shigella", "gentamicin")]

    import dna_decode.data.trust_surface as ts2
    monkeypatch.setattr(ts2, "_card_path", lambda name: tmp_path / "nope.json")
    ts2._load.cache_clear()
    assert ts2.prospective_regressions() == []                      # absent card -> silent, not a crash
    assert ts2.prospective_regression_for("gentamicin", "Escherichia_coli_Shigella") is None
    ts2._load.cache_clear()                                         # leave no stale card for other tests


def test_dna_decode_list_surfaces_the_regression_on_the_authoritative_surface(capsys, monkeypatch):
    """`dna-decode list` is the project's stated authoritative support surface, and its per-trait
    validation strings are STATIC (they quote in-distribution "gent 0.945"). Without this the one place a
    user checks before trusting a call would never mention that a post-lock cohort contradicts it."""
    import dna_decode.cli as uni
    import dna_decode.data.trust_surface as ts
    monkeypatch.setattr(ts, "prospective_regressions",
                        lambda: [{"organism": "escherichia_coli_shigella", "drug": "gentamicin",
                                  "sens": 0.429, "n": 62, "lock_date": "2026-06-13", "caveat": "x"}])
    assert uni.main(["list"]) == 0
    out = capsys.readouterr().out
    assert "PROSPECTIVE REGRESSION" in out and "0.429" in out
    assert "IN-DISTRIBUTION" in out                       # names WHY the quoted accuracy is not it
