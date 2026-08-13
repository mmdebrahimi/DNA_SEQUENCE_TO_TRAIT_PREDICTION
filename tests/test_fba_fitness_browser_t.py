"""The Fitness Browser t-statistic axis — in-memory SQLite, no 7.4 GB feba.db.

Context: a published honest-limits sentence claimed "the t-statistic that the loader reads is not used".
The loader did not read it at all — the SQL selected only `fit`. The column exists on the live db
(`['orgId','locusId','expName','fit','t']`), so the axis was one word away the whole time. These tests
pin that it is now selected, that the default path is unchanged, and that it cannot be silently dropped
again.
"""
from __future__ import annotations

import sqlite3

from dna_decode.fba import fitness_browser as fb
from dna_decode.fba.fitness_browser import load_records, mean_t_matrix

CONDS = {"glucose": "EX_glc__D_e", "acetate": "EX_ac_e"}


def _db(rows):
    """rows: (sysName, cond, fit, t). Builds the real 5-column GeneFitness shape."""
    conn = sqlite3.connect(":memory:")
    conn.executescript(
        "CREATE TABLE Experiment (orgId TEXT, expName TEXT, expGroup TEXT, condition_1 TEXT);"
        "CREATE TABLE Gene (orgId TEXT, locusId TEXT, sysName TEXT);"
        "CREATE TABLE GeneFitness (orgId TEXT, locusId TEXT, expName TEXT, fit REAL, t REAL);")
    exps, genes = {}, {}
    for i, (sysname, cond, fit, t) in enumerate(rows):
        if cond not in exps:
            exps[cond] = f"exp{len(exps)}"
            conn.execute("INSERT INTO Experiment VALUES (?,?,?,?)",
                         (fb.ORG_ID, exps[cond], "carbon source", cond))
        if sysname not in genes:
            genes[sysname] = f"loc{len(genes)}"
            conn.execute("INSERT INTO Gene VALUES (?,?,?)", (fb.ORG_ID, genes[sysname], sysname))
        conn.execute("INSERT INTO GeneFitness VALUES (?,?,?,?,?)",
                     (fb.ORG_ID, genes[sysname], exps[cond], fit, t))
    conn.commit()
    return conn


def test_the_loader_selects_the_t_column():
    """Guard against a future edit silently dropping `f.t` again -- the exact regression that made a
    published sentence false."""
    import inspect

    src = inspect.getsource(load_records)
    assert "f.t" in src, "load_records no longer selects the t-statistic"


def test_default_path_is_byte_identical_to_the_no_t_behaviour():
    """min_abs_t=None must not move a single record, or every shipped 25-source number is in play."""
    conn = _db([
        ("b0001", "glucose", -3.0, 8.0), ("b0001", "acetate", 0.1, 0.4),
        ("b0002", "glucose", 0.2, 5.0), ("b0002", "acetate", -4.0, 9.0),
    ])
    recs = load_records(conn, CONDS)
    assert {r.gene_id for r in recs} == {"b0001", "b0002"}
    assert all(r.conditionally_essential for r in recs)
    assert recs == load_records(conn, CONDS, min_abs_t=None)


def test_a_low_confidence_measurement_drops_the_gene_when_a_t_bar_is_set():
    """b0001's acetate call rests on |t|=0.4 -- noise. With a bar it leaves; without one it stays,
    which is precisely the sensitivity the 217-gene set has never been tested for."""
    conn = _db([
        ("b0001", "glucose", -3.0, 8.0), ("b0001", "acetate", 0.1, 0.4),
        ("b0002", "glucose", 0.2, 5.0), ("b0002", "acetate", -4.0, 9.0),
    ])
    assert {r.gene_id for r in load_records(conn, CONDS)} == {"b0001", "b0002"}
    kept = {r.gene_id for r in load_records(conn, CONDS, min_abs_t=3.0)}
    assert kept == {"b0002"}


def test_the_t_bar_uses_absolute_value_so_strong_negatives_survive():
    """A strongly-depleted gene has a large NEGATIVE t; a naive `t >= bar` would drop exactly the
    genes the analysis is about."""
    conn = _db([("b0003", "glucose", -5.0, -12.0), ("b0003", "acetate", 0.0, 7.0)])
    assert {r.gene_id for r in load_records(conn, CONDS, min_abs_t=5.0)} == {"b0003"}


def test_replicates_are_averaged_on_the_t_axis_the_same_way_as_fitness():
    conn = _db([
        ("b0001", "glucose", -3.0, 8.0), ("b0001", "glucose", -1.0, 4.0),
        ("b0001", "acetate", 0.1, 1.0),
    ])
    got = mean_t_matrix(conn, CONDS, {"b0001"})
    assert got["glucose"]["b0001"] == 6.0            # (8 + 4) / 2
    assert got["acetate"]["b0001"] == 1.0


def test_a_gene_missing_t_entirely_is_dropped_under_a_bar_not_admitted_by_default():
    """Fail-closed: an unprovable confidence is not a passing one."""
    conn = _db([("b0004", "glucose", -3.0, None), ("b0004", "acetate", 0.1, None)])
    assert {r.gene_id for r in load_records(conn, CONDS)} == {"b0004"}      # default unaffected
    assert load_records(conn, CONDS, min_abs_t=1.0) == []                   # bar -> dropped
