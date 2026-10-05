"""Tests for the single GenBank `  ORGANISM` reader (Step 1).

Stdlib only, no network. The real-data tests skip when the cached assemblies are absent, because the
`data/` tree is a junction to D: and may not be mounted.
"""
from __future__ import annotations

from pathlib import Path

import pytest

from dna_decode.data.genbank_organism import (
    DEFAULT_MAX_HEADER_LINES,
    genus_of,
    organism_from_genbank,
    species_of,
)

REPO = Path(__file__).resolve().parent.parent
MEROP = REPO / "data" / "raw" / "klebsiella_provdisjoint_meropenem" / "refseq"


def _gbk(acc: str) -> Path:
    return MEROP / acc / "annotations.gbk"


# --------------------------------------------------------------------------- parsing

def test_reads_a_real_organism_line(tmp_path):
    p = tmp_path / "a.gbk"
    p.write_text("LOCUS       x\nDEFINITION  y\n  ORGANISM  Klebsiella aerogenes MGH 77\n",
                 encoding="utf-8")
    assert organism_from_genbank(p) == "Klebsiella aerogenes MGH 77"


def test_a_missing_file_returns_None_rather_than_raising():
    """One bad genome must not abort a 2,000-genome enumeration."""
    assert organism_from_genbank(Path("/nope/absent.gbk")) is None


def test_an_empty_file_returns_None(tmp_path):
    p = tmp_path / "e.gbk"; p.write_text("", encoding="utf-8")
    assert organism_from_genbank(p) is None


def test_a_file_with_no_organism_line_returns_None(tmp_path):
    p = tmp_path / "n.gbk"
    p.write_text("LOCUS x\nDEFINITION y\nSOURCE z\n", encoding="utf-8")
    assert organism_from_genbank(p) is None


def test_an_empty_organism_value_returns_None_not_empty_string(tmp_path):
    """`None` and `''` mean different things to a caller counting unresolved genomes."""
    p = tmp_path / "b.gbk"; p.write_text("  ORGANISM  \n", encoding="utf-8")
    assert organism_from_genbank(p) is None


def test_the_header_line_budget_actually_binds(tmp_path):
    """A budget that never stops the scan is not a budget. The ORGANISM line is planted PAST it."""
    p = tmp_path / "late.gbk"
    p.write_text("filler\n" * (DEFAULT_MAX_HEADER_LINES + 5) + "  ORGANISM  Escherichia coli\n",
                 encoding="utf-8")
    assert organism_from_genbank(p) is None
    # ...and raising the budget finds it, proving the None above was the budget and not a parse failure
    assert organism_from_genbank(p, max_header_lines=10_000) == "Escherichia coli"


def test_accepts_a_str_path_too(tmp_path):
    p = tmp_path / "s.gbk"; p.write_text("  ORGANISM  Salmonella enterica\n", encoding="utf-8")
    assert organism_from_genbank(str(p)) == "Salmonella enterica"  # type: ignore[arg-type]


# --------------------------------------------------------------------------- genus vs species

def test_genus_and_species_split_the_two_levels_a_genus_named_cohort_conflates():
    """The whole point: `Klebsiella aerogenes` is the SAME genus and a DIFFERENT species from
    `Klebsiella pneumoniae`, which is exactly the case a `klebsiella_*` cohort name hides."""
    aer, pne = "Klebsiella aerogenes MGH 77", "Klebsiella pneumoniae BIDMC 53"
    assert genus_of(aer) == genus_of(pne) == "Klebsiella"
    assert species_of(aer) == "Klebsiella aerogenes"
    assert species_of(pne) == "Klebsiella pneumoniae"
    assert species_of(aer) != species_of(pne)


def test_species_of_refuses_a_genus_only_name():
    """Returning `'Klebsiella'` would let a caller mistake a genus for a resolved species."""
    assert species_of("Klebsiella") is None
    assert genus_of("Klebsiella") == "Klebsiella"


def test_both_handle_empty_input():
    for f in (genus_of, species_of):
        assert f("") is None


# --------------------------------------------------------------------------- the relocation

def test_the_script_reexport_is_the_SAME_object_not_a_copy():
    """A second implementation would drift. The script must alias, not re-define."""
    from scripts.build_identify_reference import organism_from_genbank as from_script
    assert from_script is organism_from_genbank


def test_the_relocated_function_is_byte_identical_on_REAL_assemblies():
    """BEHAVIOUR-PRESERVATION PIN for the move. A refactor that changes a value is a defect, not a move.

    Compares against the pre-move implementation, reconstructed inline here rather than imported -- the
    original body is gone, so the only way to assert equivalence is to re-run the old algorithm.
    Deliberately includes a *K. aerogenes* genome, because that species is the finding this whole plan
    exists to measure and a parser that mangled its name would invalidate the audit.
    """
    if not MEROP.is_dir():
        pytest.skip("cached meropenem assemblies absent (data/ junction may be unmounted)")

    def pre_move(gbk: Path, max_header_lines: int = 60) -> str | None:
        try:
            with gbk.open(encoding="utf-8", errors="replace") as fh:
                for i, line in enumerate(fh):
                    if line.startswith("  ORGANISM"):
                        got = line.split("ORGANISM", 1)[1].strip()
                        return got or None
                    if i > max_header_lines:
                        return None
        except OSError:
            return None
        return None

    # one K. pneumoniae (a true positive), two K. aerogenes (false negatives) -- the cases that matter
    accs = ["GCA_000567045.1", "GCA_000692195.1", "GCA_003950575.1"]
    checked = 0
    for acc in accs:
        p = _gbk(acc)
        if not p.exists():
            continue
        assert organism_from_genbank(p) == pre_move(p), acc
        checked += 1
    assert checked >= 2, f"pin is vacuous: only {checked} real assemblies were comparable"


def test_the_real_aerogenes_genome_resolves_to_aerogenes():
    """Anchors the plan's premise on real data: this accession is a cohort false negative and it is
    NOT K. pneumoniae. If this ever flips, the species finding needs re-deriving, not patching."""
    p = _gbk("GCA_000692195.1")
    if not p.exists():
        pytest.skip("cached assembly absent")
    got = organism_from_genbank(p)
    assert got is not None and got.startswith("Klebsiella aerogenes"), got
    assert species_of(got) == "Klebsiella aerogenes"
