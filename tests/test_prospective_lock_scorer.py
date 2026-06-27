"""Regression tests for the prospective-lock live scorer (_predict_eligible) fix (2026-06-27).

Two latent bugs (never bit because the post-lock cohort is empty/accruing):
  1. `_run_dir(own_runs, gca)` called with swapped/short args vs `_run_dir(acc, own, reuse_glob)`.
     Fixed to use `ensure_run`'s returned main.tsv path directly.
  2. `y` stored as the raw "R"/"S" string while `independent_cohort_validate._conf` needs int 1/0.

Mocks `ensure_run` + `call_resistance` (lazy-imported inside the function) so no network/Docker fires.
"""
from __future__ import annotations

import scripts.prospective_lock_validate as plv


class _Args:
    drug = "ciprofloxacin"
    organism = "Escherichia_coli_Shigella"
    amrfinder_organism = "Escherichia"


def test_predict_eligible_empty_returns_empty():
    assert plv._predict_eligible([], _Args()) == []


def test_predict_eligible_label_to_int_and_ensure_run_path(tmp_path, monkeypatch):
    import dna_decode.eval.amr_rules as amr
    import scripts.organism_drug_validate as odv

    fake_main = tmp_path / "main.tsv"
    fake_main.write_text("Element symbol\tClass\n", encoding="utf-8")
    monkeypatch.setattr(odv, "ensure_run", lambda *a, **k: fake_main)
    monkeypatch.setattr(amr, "call_resistance", lambda mt, drug, organism=None: {"prediction": "R"})

    recs = plv._predict_eligible([{"gca": "GCA_1.1", "label": "R"},
                                  {"gca": "GCA_2.2", "label": "S"}], _Args())
    assert recs == [{"gca": "GCA_1.1", "prediction": "R", "y": 1},
                    {"gca": "GCA_2.2", "prediction": "R", "y": 0}]
    assert all(isinstance(r["y"], int) for r in recs)  # the bug-2 fix: y is int, not "R"/"S"


def test_predict_eligible_skips_when_ensure_run_none(monkeypatch):
    import scripts.organism_drug_validate as odv
    monkeypatch.setattr(odv, "ensure_run", lambda *a, **k: None)
    assert plv._predict_eligible([{"gca": "GCA_x.1", "label": "R"}], _Args()) == []


def test_score_consumes_int_y_end_to_end(tmp_path, monkeypatch):
    """The y=int records must flow cleanly into _score -> _conf (the original crash path)."""
    import dna_decode.eval.amr_rules as amr
    import scripts.organism_drug_validate as odv

    fake_main = tmp_path / "main.tsv"
    fake_main.write_text("Element symbol\tClass\n", encoding="utf-8")
    monkeypatch.setattr(odv, "ensure_run", lambda *a, **k: fake_main)
    # one true-positive (R/R) + one false-positive (R-pred / S-label)
    monkeypatch.setattr(amr, "call_resistance", lambda mt, drug, organism=None: {"prediction": "R"})

    recs = plv._predict_eligible([{"gca": "g1", "label": "R"}, {"gca": "g2", "label": "S"}], _Args())
    conf = plv.conf_from_records(recs)
    assert conf["tp"] == 1 and conf["fp"] == 1
