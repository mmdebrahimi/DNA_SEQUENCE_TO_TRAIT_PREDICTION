"""EMBL-EBI AMR Portal header dialect + the fail-closed guard.

The load-bearing test is `test_the_portal_dialect_silently_flips_a_real_determinant_to_S`: it pins the
MEASURED defect (gyrA_T86I scoring SUSCEPTIBLE for ciprofloxacin) through the DEPLOYED `call_resistance`,
so the guard can never be weakened into a no-op without a red test.
"""
from __future__ import annotations

import gzip
from pathlib import Path

import pytest

from dna_decode.data.amr_portal import (
    DIALECT_AGNOSTIC_COLUMNS,
    REPO_SPACE_COLUMNS,
    PortalHeaderIncompatible,
    assert_repo_compatible,
    genome_url,
    header_dialect,
    normalize_header,
    prepare,
    read_text,
)
from dna_decode.eval.amr_rules import call_resistance

# A REAL portal record, transcribed from
# genomes/GCA/001/416/GCA_001416835.1/GCA_001416835.1_amrfinderplus.tsv.gz (Campylobacter jejuni).
# Trimmed to the columns the rules read; the dialect and the determinant symbols are verbatim.
_PORTAL_COLS = ("Protein_id\tContig_id\tStart\tStop\tStrand\tElement_symbol\tElement_name\tScope\t"
                "Type\tSubtype\tClass\tSubclass\tMethod")
_PORTAL_ROWS = [
    "p1\tc1\t1\t2\t+\tgyrA_T86I\tquinolone resistance gyrA T86I\tcore\tAMR\tPOINT\tQUINOLONE\tQUINOLONE\tPOINTP",
    "p2\tc1\t3\t4\t+\ttet(O)\ttetracycline resistance Tet(O)\tcore\tAMR\tAMR\tTETRACYCLINE\tTETRACYCLINE\tBLASTP",
    "p3\tc1\t5\t6\t-\taph(2'')-If\taminoglycoside phosphotransferase\tcore\tAMR\tAMR\tAMINOGLYCOSIDE\t"
    "AMIKACIN/GENTAMICIN/KANAMYCIN/TOBRAMYCIN\tEXACTP",
]
PORTAL_TSV = _PORTAL_COLS + "\n" + "\n".join(_PORTAL_ROWS) + "\n"
REPO_TSV = normalize_header(PORTAL_TSV)


# ------------------------------------------------------------------ the measured defect

def test_the_portal_dialect_silently_flips_a_real_determinant_to_S(tmp_path):
    """MEASURED on real portal data: as-delivered, gyrA_T86I scores SUSCEPTIBLE for ciprofloxacin.

    gyrA_T86I is the determinant this project's own Campylobacter cipro cell is overwhelmingly built
    on, so this is not an edge case. It must fail loudly if the guard is ever weakened.
    """
    asis = tmp_path / "asis.tsv"
    asis.write_text(PORTAL_TSV, encoding="utf-8")
    norm = tmp_path / "norm.tsv"
    norm.write_text(REPO_TSV, encoding="utf-8")

    a = call_resistance(asis, "ciprofloxacin", organism="Campylobacter")
    b = call_resistance(norm, "ciprofloxacin", organism="Campylobacter")
    assert a["prediction"] == "S", "fixture no longer reproduces the as-delivered defect"
    assert b["prediction"] == "R", "normalized header should recover the gyrA_T86I call"


def test_tetracycline_also_flips_because_its_rule_is_symbol_keyed(tmp_path):
    asis = tmp_path / "a.tsv"; asis.write_text(PORTAL_TSV, encoding="utf-8")
    norm = tmp_path / "n.tsv"; norm.write_text(REPO_TSV, encoding="utf-8")
    assert call_resistance(asis, "tetracycline", organism="Campylobacter")["prediction"] == "S"
    assert call_resistance(norm, "tetracycline", organism="Campylobacter")["prediction"] == "R"


def test_gentamicin_is_RIGHT_BY_LUCK_in_both_dialects(tmp_path):
    """The dangerous part, pinned: a Subclass-keyed rule is unaffected, so spot-checking gentamicin
    alone makes the dialect bug invisible."""
    asis = tmp_path / "a.tsv"; asis.write_text(PORTAL_TSV, encoding="utf-8")
    norm = tmp_path / "n.tsv"; norm.write_text(REPO_TSV, encoding="utf-8")
    assert call_resistance(asis, "gentamicin", organism="Campylobacter")["prediction"] == \
           call_resistance(norm, "gentamicin", organism="Campylobacter")["prediction"]
    assert "Subclass" in DIALECT_AGNOSTIC_COLUMNS


# ------------------------------------------------------------------ dialect detection

def test_detects_each_dialect():
    assert header_dialect(PORTAL_TSV) == "portal"
    assert header_dialect(REPO_TSV) == "repo"
    assert header_dialect("alpha\tbeta\ngamma\tdelta\n") == "unknown"
    assert header_dialect("") == "unknown"


def test_guard_refuses_the_portal_dialect_and_names_the_affected_rules():
    with pytest.raises(PortalHeaderIncompatible) as e:
        assert_repo_compatible(PORTAL_TSV)
    msg = str(e.value)
    assert "Element_symbol" in msg
    assert "qrdr_point" in msg and "gene_prefixes" in msg
    assert "Subclass-keyed rules are unaffected" in msg      # the spot-check warning must survive


def test_guard_refuses_an_unrecognised_schema_rather_than_scoring_it():
    with pytest.raises(PortalHeaderIncompatible, match="either dialect"):
        assert_repo_compatible("foo\tbar\n1\t2\n")


def test_guard_accepts_the_repo_dialect():
    assert_repo_compatible(REPO_TSV)        # must not raise


# ------------------------------------------------------------------ normalization is header-only

def test_normalize_touches_ONLY_the_header_and_preserves_symbol_underscores():
    """Rewriting body underscores would corrupt `gyrA_T86I` -- exactly the QRDR point mutations the
    ciprofloxacin rule counts. This is the subtlest requirement in the module."""
    out = normalize_header(PORTAL_TSV)
    assert "Element symbol" in out.splitlines()[0]
    assert "Element_symbol" not in out.splitlines()[0]
    assert "gyrA_T86I" in out, "body underscore was rewritten -- the determinant symbol is corrupted"
    assert "rplV" not in out or "_" in out.split("\n", 1)[1]


def test_normalize_is_idempotent_on_the_repo_dialect():
    assert normalize_header(REPO_TSV) == REPO_TSV


def test_normalize_handles_empty_input():
    assert normalize_header("") == ""


# ------------------------------------------------------------------ prepare()

def test_prepare_normalizes_and_writes_a_repo_readable_file(tmp_path):
    src = tmp_path / "p.tsv"; src.write_text(PORTAL_TSV, encoding="utf-8")
    out = prepare(src, tmp_path / "sub" / "out.tsv")
    assert out.exists()
    assert_repo_compatible(out.read_text(encoding="utf-8"))
    assert call_resistance(out, "ciprofloxacin", organism="Campylobacter")["prediction"] == "R"


def test_prepare_passes_through_a_repo_dialect_file_unchanged(tmp_path):
    src = tmp_path / "r.tsv"; src.write_text(REPO_TSV, encoding="utf-8")
    out = prepare(src, tmp_path / "o.tsv")
    assert out.read_text(encoding="utf-8") == REPO_TSV


def test_prepare_refuses_an_unrecognised_schema(tmp_path):
    src = tmp_path / "x.tsv"; src.write_text("a\tb\n1\t2\n", encoding="utf-8")
    with pytest.raises(PortalHeaderIncompatible):
        prepare(src, tmp_path / "o.tsv")


def test_read_text_transparently_gunzips(tmp_path):
    gz = tmp_path / "g.tsv.gz"
    with gzip.open(gz, "wt", encoding="utf-8") as fh:
        fh.write(PORTAL_TSV)
    assert read_text(gz) == PORTAL_TSV
    assert header_dialect(read_text(gz)) == "portal"


# ------------------------------------------------------------------ URL layout

def test_genome_url_matches_the_documented_three_three_split():
    # verified live: this exact URL served the real file used for the fixture above
    assert genome_url("GCA_001416835.1") == (
        "https://ftp.ebi.ac.uk/pub/databases/amr_portal/genomes/GCA/001/416/"
        "GCA_001416835.1/GCA_001416835.1_amrfinderplus.tsv.gz")


def test_genome_url_handles_an_ERZ_accession():
    u = genome_url("ERZ25223465")
    assert "/ERZ/252/234/ERZ25223465/" in u


def test_genome_url_refuses_a_non_assembly_accession():
    for bad in ("SAMEA1028830", "ERS022205", "GCF_000005845.2", "GCA_12.1"):
        with pytest.raises(ValueError):
            genome_url(bad)


def test_repo_space_columns_are_actually_what_amr_rules_reads():
    """Guard against this list drifting away from the parser it describes."""
    src = Path("dna_decode/eval/amr_rules.py").read_text(encoding="utf-8")
    assert '"Element symbol"' in src
    assert "Element symbol" in REPO_SPACE_COLUMNS
