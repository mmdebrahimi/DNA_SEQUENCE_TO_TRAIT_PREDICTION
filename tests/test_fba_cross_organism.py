"""Offline tests for the cross-organism axis.

The equivalence of `load_org_records` to the pinned `load_records` is the load-bearing check: the whole
P. putida result rests on a loader I wrote, so it must reproduce the audited one exactly on E. coli. The
real-data version of that check needs feba.db, so it is skipped when the DB is absent; these pure tests
pin the CONTRACT that makes the equivalence possible.
"""
from __future__ import annotations

import sqlite3

import pytest

from scripts.fba_cross_organism import load_org_records, org_conditions


class _Met:
    def __init__(self, name):
        self.name = name


class _Rxn:
    """Exchange stand-in. `build_exchange_name_index` reads the exchanged METABOLITE's name, so a
    reaction id alone is not enough -- the fake has to carry one."""

    def __init__(self, rid, met_name):
        self.id = rid
        self.metabolites = [_Met(met_name)]


class _FakeModel:
    def __init__(self, exchanges: dict):
        self.id = "fake"
        self.reactions = [_Rxn(rid, nm) for rid, nm in exchanges.items()]
        self.exchanges = list(self.reactions)


def _db() -> sqlite3.Connection:
    """Two organisms, deliberately overlapping on experiment names, to catch cross-organism bleed."""
    c = sqlite3.connect(":memory:")
    c.execute("CREATE TABLE Experiment (orgId TEXT, expName TEXT, expGroup TEXT, condition_1 TEXT)")
    c.execute("CREATE TABLE Gene (orgId TEXT, locusId TEXT, sysName TEXT)")
    c.execute("CREATE TABLE GeneFitness (orgId TEXT, locusId TEXT, expName TEXT, fit REAL)")
    rows = [
        ("Putida", "e1", "carbon source", "D-Glucose"),
        ("Putida", "e2", "carbon source", "D-Glucose"),   # replicate
        ("Putida", "e3", "carbon source", "Ethanol"),
        ("Putida", "e4", "nitrogen source", "Ammonium"),  # different group -- must be excluded
        ("Keio", "e1", "carbon source", "D-Glucose"),     # SAME expName, different organism
    ]
    c.executemany("INSERT INTO Experiment VALUES (?,?,?,?)", rows)
    c.executemany("INSERT INTO Gene VALUES (?,?,?)",
                  [("Putida", "L1", "PP_0001"), ("Putida", "L2", "PP_0002"),
                   ("Keio", "K1", "b0001")])
    c.executemany("INSERT INTO GeneFitness VALUES (?,?,?,?)", [
        ("Putida", "L1", "e1", -3.0), ("Putida", "L1", "e2", -1.0),  # mean -2.0 (NOT < -2.0)
        ("Putida", "L1", "e3", 0.5),
        ("Putida", "L2", "e1", -5.0), ("Putida", "L2", "e2", -5.0), ("Putida", "L2", "e3", 0.1),
        ("Keio", "K1", "e1", -9.0),                                  # must never reach Putida records
    ])
    return c


CONDS = {"D-Glucose": "EX_glc__D_e", "Ethanol": "EX_etoh_e"}


def test_replicates_are_averaged_not_taken_singly():
    """L1: fit -3.0 and -1.0 average to exactly -2.0, which is NOT < -2.0 -> dispensable."""
    recs = {r.gene_id: r for r in load_org_records(_db(), "Putida", "carbon source", CONDS)}
    assert recs["PP_0001"].experimental["D-Glucose"] is False
    assert recs["PP_0002"].experimental["D-Glucose"] is True


def test_another_organism_sharing_an_expName_cannot_bleed_in():
    """`Keio` reuses expName 'e1'. A join that forgets orgId would import its fitness value."""
    recs = load_org_records(_db(), "Putida", "carbon source", CONDS)
    assert {r.gene_id for r in recs} == {"PP_0001", "PP_0002"}
    assert all(not r.gene_id.startswith("b") for r in recs)


def test_a_different_expGroup_is_excluded():
    recs = load_org_records(_db(), "Putida", "carbon source", CONDS)
    assert all(set(r.experimental) == set(CONDS) for r in recs)


def test_a_partial_row_is_dropped_not_defaulted():
    """A gene missing a condition would look dispensable there. It must be dropped entirely."""
    c = _db()
    c.execute("INSERT INTO Gene VALUES ('Putida','L3','PP_0003')")
    c.execute("INSERT INTO GeneFitness VALUES ('Putida','L3','e1',-9.0)")  # glucose only
    recs = {r.gene_id for r in load_org_records(c, "Putida", "carbon source", CONDS)}
    assert "PP_0003" not in recs


def test_gene_filter_restricts_to_model_genes():
    recs = load_org_records(_db(), "Putida", "carbon source", CONDS, gene_filter={"PP_0002"})
    assert {r.gene_id for r in recs} == {"PP_0002"}


def test_two_sided_flag_is_recomputed_from_the_calls():
    recs = {r.gene_id: r for r in load_org_records(_db(), "Putida", "carbon source", CONDS)}
    r2 = recs["PP_0002"]
    assert r2.experimental == {"D-Glucose": True, "Ethanol": False}
    assert r2.conditionally_essential is True


def test_org_conditions_only_returns_sources_the_model_can_represent():
    c = _db()
    m = _FakeModel({"EX_glc__D_e": "D-Glucose"})  # model lacks an ethanol exchange
    got = org_conditions(c, m, "Putida", "carbon source")
    assert "Ethanol" not in got


@pytest.mark.skipif(
    not __import__("pathlib").Path("D:/dna_decode_cache/fitness_browser/feba.db").exists(),
    reason="feba.db not attached")
def test_generalized_loader_equals_the_pinned_loader_on_ecoli():
    """The load-bearing check: my loader must reproduce the audited one EXACTLY on E. coli.

    Measured 2026-08-20: 1,339 genes both ways, identical gene set, 0 differing experimental calls.
    """
    from dna_decode.fba.fitness_browser import (
        ESSENTIAL_FITNESS, carbon_conditions, load_records, open_db,
    )
    from dna_decode.fba.model import load_model

    m = load_model()
    conn = open_db()
    conds = carbon_conditions(conn, m)
    gf = {g.id for g in m.genes}
    pinned = {r.gene_id: tuple(sorted(r.experimental.items()))
              for r in load_records(conn, conds, gene_filter=gf, threshold=ESSENTIAL_FITNESS)}
    mine = {r.gene_id: tuple(sorted(r.experimental.items()))
            for r in load_org_records(conn, "Keio", "carbon source", conds,
                                      gene_filter=gf, threshold=ESSENTIAL_FITNESS)}
    assert set(pinned) == set(mine)
    assert [g for g in pinned if pinned[g] != mine[g]] == []
