"""Pure-helper tests for the naive-AMRFinder comparators (Oxford + provenance-disjoint).

The comparators score a NAIVE baseline ("any AMRFinder determinant of the drug class -> R")
against the FROZEN call_resistance rule to validate the wrapper-vs-underlying-tool rail. These
tests pin the pure parsing/prediction helpers; the full scoring run needs cached AMRFinder data
(gitignored) and is exercised manually.
"""
from __future__ import annotations

from pathlib import Path

import scripts.naive_baseline_provdisjoint as pd
import scripts.oxford_naive_baseline as ox

_HEADER = ("Name\tProtein identifier\tContig id\tStart\tStop\tStrand\tElement symbol\t"
           "Sequence name\tScope\tElement type\tElement subtype\tClass\tSubclass\tMethod")


def _write_main(tmp: Path, rows: list[str]) -> Path:
    p = tmp / "main.tsv"
    p.write_text(_HEADER + "\n" + "\n".join(rows) + ("\n" if rows else ""), encoding="utf-8")
    return p


def _row(symbol: str, cls: str) -> str:
    cells = [""] * 14
    cells[6] = symbol
    cells[11] = cls
    return "\t".join(cells)


def test_naive_predict_quinolone_class_is_R_for_cipro(tmp_path):
    mt = _write_main(tmp_path, [_row("gyrA_S83L", "QUINOLONE")])
    assert pd.naive_predict(mt, "ciprofloxacin") == "R"


def test_naive_predict_betalactam_only_is_S_for_cipro(tmp_path):
    # a beta-lactam determinant must NOT trigger a ciprofloxacin R call
    mt = _write_main(tmp_path, [_row("blaTEM-1", "BETA-LACTAM")])
    assert pd.naive_predict(mt, "ciprofloxacin") == "S"


def test_naive_predict_betalactam_is_R_for_ceftriaxone(tmp_path):
    # the intrinsic-gene over-call: ANY beta-lactam class -> naive cef R (the curated rule refines this away)
    mt = _write_main(tmp_path, [_row("blaEC", "BETA-LACTAM")])
    assert pd.naive_predict(mt, "ceftriaxone") == "R"


def test_naive_predict_empty_is_S(tmp_path):
    mt = _write_main(tmp_path, [])
    assert pd.naive_predict(mt, "gentamicin") == "S"


def test_read_selected_parses_acc_label(tmp_path):
    slug = tmp_path
    (slug / "selected.tsv").write_text("GCA_1.1\tR\nGCA_2.2\tS\n\nGCA_3.3\tX\n", encoding="utf-8")
    sel = pd.read_selected(slug)
    assert sel == [("GCA_1.1", "R"), ("GCA_2.2", "S")]  # blank + invalid-label rows dropped


def test_cm_tuple_order(tmp_path):
    conf = {"tp": 1, "fp": 2, "tn": 3, "fn": 4}
    assert pd._cm(conf) == (1, 2, 3, 4)


def test_oxford_class_col_index():
    header = "\t".join(["Name", "Element symbol", "Class", "Subclass"])
    assert ox._class_col_index(header) == 2


def test_oxford_naive_predict_class_match(tmp_path):
    amr = {"g1": ["g1\tx\tQUINOLONE".replace("x", "x")]}  # 3 cols: Name, sym, Class
    header = "Name\tElement symbol\tClass"
    ci = ox._class_col_index(header)
    assert ox.naive_predict("g1", "ciprofloxacin", amr, ci) == "R"
    assert ox.naive_predict("absent", "ciprofloxacin", amr, ci) == "S"
