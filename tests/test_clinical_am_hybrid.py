"""Offline tests for the AM/hybrid clinical-AUROC extension (scripts/clinical_am_hybrid_auroc.py).

Pure file-parsing helpers only — no AlphaMissense download, no ESM2/ProSST model, no network. The AUROC math
+ join logic are covered by tests/test_clinical_variant_effect.py (reused here).
"""
from __future__ import annotations

import gzip

from scripts.clinical_am_hybrid_auroc import build_am_filter, load_am, GENE_UNIPROT


def _write_am_gz(path, rows):
    with gzip.open(path, "wt", encoding="utf-8") as f:
        f.write("# Copyright DeepMind — CC BY-NC-SA\n")
        f.write("uniprot_id\tprotein_variant\tam_pathogenicity\tam_class\n")
        for up, var, path_score, cls in rows:
            f.write(f"{up}\t{var}\t{path_score}\t{cls}\n")


def test_build_am_filter_keeps_only_requested_uniprots(tmp_path):
    gz = tmp_path / "am.tsv.gz"
    out = tmp_path / "filtered.tsv"
    _write_am_gz(gz, [
        ("P04637", "R248W", "0.99", "pathogenic"),
        ("P04637", "P72R", "0.05", "benign"),
        ("Q99999", "A1V", "0.5", "ambiguous"),   # not requested -> dropped
        ("P43246", "G674R", "0.97", "pathogenic"),
    ])
    n = build_am_filter({"P04637", "P43246"}, am_gz=gz, out=out)
    assert n == 3   # the Q99999 row is excluded
    text = out.read_text(encoding="utf-8")
    assert "Q99999" not in text
    assert "R248W" in text and "G674R" in text


def test_build_am_filter_skips_header_and_comment(tmp_path):
    gz = tmp_path / "am.tsv.gz"
    out = tmp_path / "filtered.tsv"
    _write_am_gz(gz, [("P04637", "R248W", "0.99", "pathogenic")])
    build_am_filter({"P04637"}, am_gz=gz, out=out)
    # the header line 'uniprot_id\t...' must NOT be written (its first field != a requested uniprot)
    assert "uniprot_id" not in out.read_text(encoding="utf-8")


def test_load_am_parses_variants_and_scores(tmp_path):
    out = tmp_path / "filtered.tsv"
    out.write_text(
        "P04637\tR248W\t0.99\tpathogenic\n"
        "P04637\tP72R\t0.05\tbenign\n"
        "P43246\tG674R\t0.97\tpathogenic\n",
        encoding="utf-8")
    am = load_am("P04637", filtered=out)
    assert am[("R", 248, "W")] == 0.99
    assert am[("P", 72, "R")] == 0.05
    assert ("G", 674, "R") not in am   # that row is P43246, not P04637
    assert load_am("P43246", filtered=out)[("G", 674, "R")] == 0.97


def test_load_am_skips_malformed_rows(tmp_path):
    out = tmp_path / "filtered.tsv"
    out.write_text(
        "P04637\tR248W\t0.99\tpathogenic\n"
        "P04637\tBADVARIANT\t0.5\tambiguous\n"      # pos not integer -> skipped
        "P04637\tX\tnotafloat\tambiguous\n",         # malformed -> skipped
        encoding="utf-8")
    am = load_am("P04637", filtered=out)
    assert am == {("R", 248, "W"): 0.99}


def test_load_am_missing_file_returns_empty(tmp_path):
    assert load_am("P04637", filtered=tmp_path / "nope.tsv") == {}


def test_gene_uniprot_map_is_the_auroc_viable_pair():
    assert GENE_UNIPROT == {"TP53": "P04637", "MSH2": "P43246"}
