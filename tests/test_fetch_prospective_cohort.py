"""Offline tests for the prospective-cohort fetch script (pure parsers + funnel; no network)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from scripts.fetch_prospective_cohort import (  # noqa: E402
    build_cohort_rows, cells_by_taxon, parse_bvbrc_amr, parse_bvbrc_genomes,
    parse_datasets_report, phenotype_to_label, write_cohort_tsv,
)


def test_phenotype_to_label():
    assert phenotype_to_label("Resistant") == "R"
    assert phenotype_to_label("Susceptible") == "S"
    assert phenotype_to_label("Intermediate") is None
    assert phenotype_to_label("Non-susceptible") is None
    assert phenotype_to_label("") is None and phenotype_to_label(None) is None


def test_cells_by_taxon_covers_scored_grid():
    c = cells_by_taxon()
    assert c[562] >= {"ciprofloxacin", "ceftriaxone", "gentamicin", "tetracycline"}   # E. coli
    assert c[573] >= {"ciprofloxacin", "meropenem", "tetracycline"}                   # Klebsiella
    assert "ciprofloxacin" in c[194]                                                  # Campylobacter


def test_parse_bvbrc_genomes_keeps_only_public_with_assembly():
    rows = [
        {"genome_id": "562.1", "assembly_accession": "GCA_001.1", "public": True},
        {"genome_id": "562.2", "assembly_accession": "GCF_002.1", "public": True},
        {"genome_id": "562.3", "assembly_accession": "", "public": True},          # no assembly -> drop
        {"genome_id": "562.4", "assembly_accession": "GCA_004.1", "public": False}, # not public -> drop
    ]
    assert parse_bvbrc_genomes(rows) == {"562.1": "GCA_001.1", "562.2": "GCF_002.1"}


def test_parse_bvbrc_amr_excludes_computational_and_offtarget_and_ambiguous():
    rows = [
        {"genome_id": "562.1", "antibiotic": "ciprofloxacin", "resistant_phenotype": "Resistant",
         "laboratory_typing_method": "Disk diffusion"},                              # keep -> R
        {"genome_id": "562.2", "antibiotic": "ciprofloxacin", "resistant_phenotype": "Susceptible",
         "laboratory_typing_method": "Computational Prediction"},                    # circular -> drop
        {"genome_id": "562.3", "antibiotic": "colistin", "resistant_phenotype": "Resistant",
         "laboratory_typing_method": "MIC"},                                         # off-target drug -> drop
        {"genome_id": "562.4", "antibiotic": "ciprofloxacin", "resistant_phenotype": "Intermediate",
         "laboratory_typing_method": "MIC"},                                         # ambiguous -> drop
    ]
    got = parse_bvbrc_amr(rows, {"ciprofloxacin"})
    assert got == [{"genome_id": "562.1", "drug": "ciprofloxacin", "label": "R", "method": "Disk diffusion"}]


def test_parse_datasets_report():
    rep = {"reports": [{"assembly_info": {"release_date": "2026-07-01T00:00:00Z",
                                          "biosample": {"accession": "SAMN123"}, "assembly_status": "current"}}]}
    out = parse_datasets_report(rep)
    assert out == {"release_date": "2026-07-01", "biosample": "SAMN123", "status": "current"}
    assert parse_datasets_report({}) == {"release_date": "", "biosample": "", "status": ""}


def test_build_cohort_rows_funnel():
    amr = [
        {"genome_id": "g1", "drug": "ciprofloxacin", "label": "R"},   # post-lock -> eligible
        {"genome_id": "g2", "drug": "ceftriaxone", "label": "S"},     # pre-lock -> excluded
        {"genome_id": "g3", "drug": "tetracycline", "label": "R"},    # undatable -> excluded (fail-closed)
        {"genome_id": "g4", "drug": "gentamicin", "label": "R"},      # no GCA -> dropped before resolve
    ]
    gid_to_gca = {"g1": "GCA_1.1", "g2": "GCA_2.1", "g3": "GCA_3.1"}
    gca_release = {
        "GCA_1.1": {"release_date": "2026-07-01", "biosample": "SAMN1"},
        "GCA_2.1": {"release_date": "2026-05-01", "biosample": "SAMN2"},
        "GCA_3.1": {"release_date": "", "biosample": ""},
    }
    rows, stats = build_cohort_rows(amr, gid_to_gca, gca_release, "2026-06-13")
    assert rows == [{"biosample": "SAMN1", "first_public_date": "2026-07-01", "gca": "GCA_1.1",
                     "drug": "ciprofloxacin", "label": "R"}]
    assert stats["amr_records"] == 4 and stats["with_gca"] == 3 and stats["eligible"] == 1
    assert stats["excluded_pre_or_undatable"] == 1   # g2 (pre-lock); g3 fails at resolve (no date), not here


def test_write_cohort_tsv_format(tmp_path):
    rows = [{"biosample": "SAMN1", "first_public_date": "2026-07-01", "gca": "GCA_1.1",
             "drug": "ciprofloxacin", "label": "R"}]
    p = tmp_path / "prospective_cohort.tsv"
    write_cohort_tsv(rows, p)
    lines = p.read_text(encoding="utf-8").splitlines()
    assert lines[0] == "biosample\tfirst_public_date\tgca\tdrug\tlabel"
    assert lines[1] == "SAMN1\t2026-07-01\tGCA_1.1\tciprofloxacin\tR"


def test_write_cohort_tsv_header_only_when_empty(tmp_path):
    p = tmp_path / "empty.tsv"
    write_cohort_tsv([], p)
    assert p.read_text(encoding="utf-8").splitlines() == ["biosample\tfirst_public_date\tgca\tdrug\tlabel"]


if __name__ == "__main__":
    import pytest
    sys.exit(pytest.main([__file__, "-q"]))
