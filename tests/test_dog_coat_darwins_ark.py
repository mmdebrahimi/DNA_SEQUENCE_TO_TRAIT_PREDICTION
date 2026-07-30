"""Offline tests for the Darwin's Ark coat-colour phenotype ingestion (multi-hot schema, verified 2026-07-30).

Synthetic multi-hot fixture (no 2.67 GB genotype data / no network). Pins: read_phenotypes builds per-dog
colour SETS from the 0/1 presence columns, scoring_target classifies base/multi/abstain, and phenotype_summary
counts the directly-scorable single-base dogs. Runnable via pytest OR standalone.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from scripts.dog_coat_darwins_ark_validate import (  # noqa: E402
    phenotype_summary,
    read_phenotypes,
    scoring_target,
)

_HDR = ("dog_id\tQ243_black_coat_color\tQ243_liver_or_brown_coat_color\tQ243_white_coat_color"
        "\tQ243_red_coat_color\tQ243_yellow_coat_color\tQ243_grey_or_blue_coat_color"
        "\tQ243_tan_coat_color\tQ243_cream_coat_color\tnumber_of_colors_in_coat\tsingle_color_in_coat")
# dog_id, black, liver, white, red, yellow, grey, tan, cream, n_colors, single
_ROWS = [
    "1\t1\t0\t0\t0\t0\t0\t0\t0\t1\t1",     # single black -> ('base','black')
    "2\t0\t0\t1\t0\t0\t0\t0\t0\t1\t1",     # single white -> ('abstain', ...)
    "3\t1\t0\t1\t0\t0\t0\t1\t0\t3\t0",     # black+tan+white -> ('multi', {black,tan})  (tan-points)
    "4\t0\t0\t0\t1\t0\t0\t0\t0\t1\t1",     # single red -> ('base','red/yellow')
    "5\t0\t1\t0\t0\t0\t0\t0\t0\t1\t1",     # single liver -> ('base','brown/liver')
]


def _fixture(tmp_path):
    p = tmp_path / "coat.tsv"
    p.write_text(_HDR + "\n" + "\n".join(_ROWS) + "\n", encoding="utf-8")
    return p


def test_read_phenotypes_builds_color_sets(tmp_path):
    ph = {p["dog_id"]: p for p in read_phenotypes(_fixture(tmp_path))}
    assert ph["1"]["colors"] == {"black"} and ph["1"]["single"] is True
    assert ph["3"]["colors"] == {"black", "tan", "white"} and ph["3"]["n_colors"] == 3
    assert ph["5"]["colors"] == {"brown/liver"}


def test_scoring_target_classification(tmp_path):
    ph = {p["dog_id"]: p for p in read_phenotypes(_fixture(tmp_path))}
    assert scoring_target(ph["1"]) == ("base", "black")
    assert scoring_target(ph["2"])[0] == "abstain"
    assert scoring_target(ph["3"]) == ("multi", frozenset({"black", "tan"}))   # white stripped -> tan-points
    assert scoring_target(ph["4"]) == ("base", "red/yellow")
    assert scoring_target(ph["5"]) == ("base", "brown/liver")


def test_phenotype_summary_counts(tmp_path):
    s = phenotype_summary(read_phenotypes(_fixture(tmp_path)))
    assert s["n"] == 5 and s["single_colour"] == 4
    assert s["directly_scorable_single_base"] == 3        # dogs 1(black),4(red),5(liver); 2=abstain, 3=multi
    assert s["single_base_colour_counts"]["black"] == 1
    assert s["target_kinds"].get("abstain") == 1 and s["target_kinds"].get("multi") == 1


if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_") and callable(v)]
    import tempfile
    for fn in fns:
        fn(Path(tempfile.mkdtemp())); print(f"PASS {fn.__name__}")
    print(f"\n{len(fns)} passed")
