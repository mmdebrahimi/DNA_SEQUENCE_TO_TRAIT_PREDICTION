"""Offline tests for the stress-axis feasibility probe.

Pure/contract tests -- no feba.db, no cobra solves. These pin the SHAPE of the finding (an exchange is
not a stress) so a future contributor cannot quietly reintroduce the medium-swap contract for stress.
"""
from __future__ import annotations

from scripts.fba_stress_feasibility_probe import (
    CANDIDATE_EXCHANGES,
    MOLECULAR_TARGETS,
    STRESS_TOLERANCE,
    probe_growth_effect,
    target_in_model,
)


class _Rxn:
    def __init__(self, rid):
        self.id = rid


class _Gene:
    def __init__(self, gid, name=None):
        self.id = gid
        self.name = name


class _FakeModel:
    """Stand-in whose growth responds to the medium via an injected table."""

    def __init__(self, growth_by_medium=None, reactions=(), genes=()):
        self.id = "fake"
        self.medium = {}
        self._growth = growth_by_medium or {}
        self.reactions = [_Rxn(r) for r in reactions]
        self.genes = [_Gene(g) if isinstance(g, str) else _Gene(*g) for g in genes]

    def __enter__(self):
        self._saved = dict(self.medium)
        return self

    def __exit__(self, *a):
        self.medium = self._saved
        return False


def _patch_growth(monkeypatch, table):
    import scripts.fba_stress_feasibility_probe as mod
    monkeypatch.setattr(mod, "wildtype_growth", lambda m: table[tuple(sorted(m.medium))])


def test_a_supplement_that_raises_growth_is_not_a_stress(monkeypatch):
    """Acetate: real behaviour. Opening the exchange ADDS a carbon source, so growth goes UP."""
    _patch_growth(monkeypatch, {(): 0.877, ("EX_ac_e",): 1.1215})
    r = probe_growth_effect(_FakeModel(), "EX_ac_e", baseline=0.877)
    assert r["delta"] > 0
    assert r["reduces_growth"] is False


def test_an_inert_exchange_is_not_a_stress(monkeypatch):
    """The metal ions: growth is unchanged, so the model cannot express their toxicity."""
    _patch_growth(monkeypatch, {(): 0.877, ("EX_cu2_e",): 0.877})
    r = probe_growth_effect(_FakeModel(), "EX_cu2_e", baseline=0.877)
    assert r["delta"] == 0.0
    assert r["reduces_growth"] is False


def test_a_genuine_growth_reduction_would_register_as_a_stress(monkeypatch):
    """The gate is not rigged to always say NO -- a real reduction passes."""
    _patch_growth(monkeypatch, {(): 0.877, ("EX_x_e",): 0.400})
    r = probe_growth_effect(_FakeModel(), "EX_x_e", baseline=0.877)
    assert r["reduces_growth"] is True


def test_a_reduction_below_tolerance_does_not_count(monkeypatch):
    _patch_growth(monkeypatch, {(): 0.877, ("EX_x_e",): 0.877 - STRESS_TOLERANCE / 10})
    r = probe_growth_effect(_FakeModel(), "EX_x_e", baseline=0.877)
    assert r["reduces_growth"] is False


def test_metabolic_target_is_found_by_reaction_or_gene():
    m = _FakeModel(reactions=("UAGCVT", "ALAALAr"), genes=(("b3189", "murA"), ("b4053", "alr")))
    assert target_in_model(m, ["murA", "UAGCVT"]) == ["rxn:UAGCVT", "gene:murA"] or \
           set(target_in_model(m, ["murA", "UAGCVT"])) == {"rxn:UAGCVT", "gene:murA"}


def test_non_metabolic_targets_are_absent_by_construction():
    """Ribosome/gyrase targets carry an EMPTY id list -- outside a metabolic model by construction."""
    m = _FakeModel(reactions=("UAGCVT",), genes=(("b3189", "murA"),))
    for ab in ("Chloramphenicol", "Tetracycline hydrochloride", "Nalidixic acid sodium salt"):
        _, ids = MOLECULAR_TARGETS[ab]
        assert ids == []
        assert target_in_model(m, ids) == []


def test_fluoride_is_not_mapped_to_iron():
    """Regression on a wrong mapping I drafted: EX_fe2_e is FERROUS IRON, not fluoride.

    Mapping it would have manufactured a data point out of an unrelated metabolite.
    """
    assert "sodium fluoride" not in CANDIDATE_EXCHANGES
    assert "EX_fe2_e" not in CANDIDATE_EXCHANGES.values()


def test_candidate_exchanges_are_all_metabolites_not_antibiotics():
    """No antibiotic may acquire a medium mapping -- that is the contract error this probe found."""
    antibiotic_words = ("mycin", "cillin", "cycline", "cephalo", "bacitracin", "fusidic", "nalidixic")
    for name in CANDIDATE_EXCHANGES:
        low = name.lower()
        assert not any(w in low for w in antibiotic_words), f"{name} is an antibiotic, not a nutrient"
