"""Offline tests for the SARS-CoV-2 Mpro CoV-RDB validation pure logic (operator-aware fold + scoring)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from scripts.sarscov2_mpro_validate import binarize_fold, catalog_predict, score  # noqa: E402


def test_binarize_fold_operator_aware():
    # '=' -> threshold compare
    assert binarize_fold("25", "=", 2.5) == "R"
    assert binarize_fold("1.1", "=", 2.5) == "S"
    # '>' lower bound: R only if the bound already clears threshold; else can't call S (the MIC lesson)
    assert binarize_fold("428", ">", 2.5) == "R"
    assert binarize_fold("2", ">", 2.5) is None          # '>2' with cutoff 2.5 -> ambiguous, NOT S
    # '<' upper bound: S only if the bound is below threshold; else can't call R
    assert binarize_fold("1", "<", 2.5) == "S"
    assert binarize_fold("5", "<", 2.5) is None          # '<5' with cutoff 2.5 -> ambiguous, NOT R
    # junk
    assert binarize_fold("NULL", "=", 2.5) is None and binarize_fold(None, "=", 2.5) is None


def test_catalog_predict_position_mutant():
    assert catalog_predict({(166, "V")}) == "R"          # E166V catalogued
    assert catalog_predict({(132, "H")}) == "S"          # Omicron P132H not catalogued
    assert catalog_predict(set()) == "S"
    assert catalog_predict({(132, "H"), (166, "V")}) == "R"   # any catalog hit -> R


def test_score_confusion_and_rates():
    recs = [
        {"fold": "25", "cmp": "=", "mpro_posmut": {(166, "V")}},     # pred R, pheno R -> TP
        {"fold": "1.1", "cmp": "=", "mpro_posmut": {(21, "I")}},     # pred R (21I catalogued), pheno S -> FP
        {"fold": "1.0", "cmp": "=", "mpro_posmut": {(132, "H")}},    # pred S, pheno S -> TN
        {"fold": "30", "cmp": "=", "mpro_posmut": {(999, "Z")}},     # pred S (uncatalogued), pheno R -> FN
        {"fold": "2", "cmp": ">", "mpro_posmut": {(166, "V")}},      # ambiguous pheno -> not scored
    ]
    res = score(recs, 2.5)
    assert (res["tp"], res["fp"], res["tn"], res["fn"]) == (1, 1, 1, 1)
    assert res["n_scored"] == 4 and res["n_records"] == 5     # the '>2' record dropped from scoring
    assert res["sens"] == 0.5 and res["spec"] == 0.5
    assert "166V" in res["per_mutation_fold"]                # per-mutation fold table populated


if __name__ == "__main__":
    import pytest
    sys.exit(pytest.main([__file__, "-q"]))
