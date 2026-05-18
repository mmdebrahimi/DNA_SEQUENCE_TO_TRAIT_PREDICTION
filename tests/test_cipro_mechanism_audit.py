"""Tests for scripts/cipro_mechanism_audit.py — pure-logic functions only.

Pins the per-symbol classification + synonymous-SNP filter + AMRFinder
parsing dedupe logic that drives the EP1 mechanism-audit verdict. The
AMRFinder Docker invocation + cohort iteration are orchestration (skipped).
Tests cover: symbol → mechanism class mapping (with tolerant prefix match),
synonymous-SNP detection, class filter (MULTIDRUG kept for regulatory hits),
and the main.tsv POINTX → kind=mutation routing + dedupe across both TSVs.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from scripts.cipro_mechanism_audit import (
    CIPRO_LOCI_BY_MECHANISM,
    CIPRO_RELEVANT_AMR_CLASSES,
    QUINOLONE_CLASSES,
    _classify_symbol,
    _is_synonymous_point,
    _parse_amrfinder_outputs,
)


# ---- CIPRO_LOCI_BY_MECHANISM contract ---------------------------------------


def test_locus_catalog_includes_qrdr_textbook_loci():
    qrdr = CIPRO_LOCI_BY_MECHANISM["QRDR_target_alteration"]
    assert "gyrA" in qrdr
    assert "gyrB" in qrdr
    assert "parC" in qrdr
    assert "parE" in qrdr


def test_locus_catalog_includes_acrR_in_regulatory():
    # acrR added 2026-05-17 to fix the missed-regulatory-frameshift bug
    assert "acrR" in CIPRO_LOCI_BY_MECHANISM["regulatory"]


def test_locus_catalog_includes_plasmid_qnr_family():
    plasmid = CIPRO_LOCI_BY_MECHANISM["plasmid_protect_modify"]
    assert "qnrA" in plasmid
    assert "qnrB" in plasmid
    assert "qnrS" in plasmid


def test_cipro_relevant_amr_classes_keeps_multidrug():
    # Regulatory mutations (marR_V84WfsTer, acrR_S30HfsTer) come through with
    # AMRFinder Class=MULTIDRUG, not QUINOLONE. Filter MUST keep MULTIDRUG.
    assert "MULTIDRUG" in CIPRO_RELEVANT_AMR_CLASSES
    assert "QUINOLONE" in CIPRO_RELEVANT_AMR_CLASSES
    assert "FLUOROQUINOLONE" in CIPRO_RELEVANT_AMR_CLASSES


def test_quinolone_classes_are_subset_of_cipro_relevant():
    assert QUINOLONE_CLASSES <= CIPRO_RELEVANT_AMR_CLASSES


# ---- _classify_symbol --------------------------------------------------------


def test_classify_symbol_gyrA_to_qrdr():
    assert _classify_symbol("gyrA") == "QRDR_target_alteration"
    assert _classify_symbol("gyrA_S83L") == "QRDR_target_alteration"
    assert _classify_symbol("gyrA_D87N") == "QRDR_target_alteration"


def test_classify_symbol_parC_to_qrdr():
    assert _classify_symbol("parC_S80I") == "QRDR_target_alteration"
    assert _classify_symbol("parC_E84V") == "QRDR_target_alteration"


def test_classify_symbol_acrR_to_regulatory():
    # The acrR addition (2026-05-17) catches regulatory frameshifts
    assert _classify_symbol("acrR_S30HfsTer59") == "regulatory"


def test_classify_symbol_marR_to_regulatory():
    assert _classify_symbol("marR_V84WfsTer14") == "regulatory"


def test_classify_symbol_qnrB_tolerant_prefix():
    # qnrB19 not in catalog explicitly; prefix-match to qnrB should fire
    assert _classify_symbol("qnrB19") == "plasmid_protect_modify"


def test_classify_symbol_empty_returns_empty():
    assert _classify_symbol("") == ""


def test_classify_symbol_unknown_returns_empty():
    assert _classify_symbol("randomGene") == ""


def test_classify_symbol_strips_mutation_suffix():
    # "gyrA_S83L" -> split on _ -> "gyrA" -> matches catalog
    assert _classify_symbol("gyrA_S83L") == "QRDR_target_alteration"


# ---- _is_synonymous_point ---------------------------------------------------


def test_synonymous_same_aa_before_and_after():
    # G141G = synonymous (first AA == last AA)
    assert _is_synonymous_point("ompF_G141G") is True
    assert _is_synonymous_point("acrB_R620R") is True
    assert _is_synonymous_point("soxR_G121G") is True


def test_synonymous_different_aa_is_not_synonymous():
    # P226L = non-synonymous
    assert _is_synonymous_point("ompF_P226L") is False
    assert _is_synonymous_point("gyrA_S83L") is False
    assert _is_synonymous_point("parC_E84V") is False


def test_synonymous_no_underscore_returns_false():
    # Bare gene symbol without mutation suffix
    assert _is_synonymous_point("gyrA") is False
    assert _is_synonymous_point("acrR") is False


def test_synonymous_frameshift_indicator_is_not_synonymous():
    # Frameshift like S30HfsTer59 - first AA S != last char 9 (digit) -> not syn
    assert _is_synonymous_point("acrR_S30HfsTer59") is False


def test_synonymous_short_mutation_returns_false():
    # Too-short mutation suffix
    assert _is_synonymous_point("gyrA_X") is False


# ---- _parse_amrfinder_outputs ----------------------------------------------


@pytest.fixture
def temp_tsv_dir(tmp_path: Path) -> Path:
    """Create empty main.tsv + mutations.tsv files in a temp dir."""
    main = tmp_path / "main.tsv"
    mut = tmp_path / "mutations.tsv"
    main.write_text("", encoding="utf-8")
    mut.write_text("", encoding="utf-8")
    return tmp_path


def _write_tsv(path: Path, rows: list[dict]) -> None:
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    df = pd.DataFrame(rows)
    df.to_csv(path, sep="\t", index=False)


def test_parse_main_tsv_pointx_routed_to_kind_mutation(tmp_path: Path):
    """main.tsv POINTX rows should be kind=mutation, not kind=acquired."""
    main = tmp_path / "main.tsv"
    mut = tmp_path / "mutations.tsv"
    _write_tsv(main, [
        {"Gene symbol": "gyrA_S83L", "Class": "QUINOLONE", "Method": "POINTX",
         "Subclass": "", "Scope": "", "Element type": "AMR",
         "% Identity to reference sequence": "", "% Coverage of reference sequence": ""},
    ])
    _write_tsv(mut, [])

    result = _parse_amrfinder_outputs(main, mut)
    assert len(result["hits"]) == 1
    assert result["hits"][0]["kind"] == "mutation"  # POINTX → mutation
    assert result["hits"][0]["symbol"] == "gyrA_S83L"
    assert result["hits"][0]["mechanism"] == "QRDR_target_alteration"


def test_parse_main_tsv_allelex_routed_to_kind_acquired(tmp_path: Path):
    """main.tsv ALLELEX (acquired-gene Method) -> kind=acquired."""
    main = tmp_path / "main.tsv"
    mut = tmp_path / "mutations.tsv"
    _write_tsv(main, [
        {"Gene symbol": "blaCMY-2", "Class": "BETA-LACTAM", "Method": "ALLELEX",
         "Subclass": "", "Scope": "", "Element type": "AMR",
         "% Identity to reference sequence": "", "% Coverage of reference sequence": ""},
    ])
    _write_tsv(mut, [])

    result = _parse_amrfinder_outputs(main, mut)
    assert len(result["hits"]) == 1
    assert result["hits"][0]["kind"] == "acquired"
    assert result["hits"][0]["symbol"] == "blaCMY-2"


def test_parse_dedupes_pointx_across_main_and_mutations(tmp_path: Path):
    """Same symbol in BOTH main.tsv (POINTX) AND mutations.tsv -> single hit."""
    main = tmp_path / "main.tsv"
    mut = tmp_path / "mutations.tsv"
    _write_tsv(main, [
        {"Gene symbol": "gyrA_S83L", "Class": "QUINOLONE", "Method": "POINTX",
         "Subclass": "", "Scope": "", "Element type": "AMR",
         "% Identity to reference sequence": "", "% Coverage of reference sequence": ""},
    ])
    _write_tsv(mut, [
        {"Gene symbol": "gyrA_S83L", "Class": "QUINOLONE", "Method": "POINTX",
         "Subclass": "", "Scope": "", "Element type": "AMR",
         "% Identity to reference sequence": "", "% Coverage of reference sequence": ""},
    ])

    result = _parse_amrfinder_outputs(main, mut)
    gyrA_hits = [h for h in result["hits"] if h["symbol"] == "gyrA_S83L"]
    assert len(gyrA_hits) == 1  # deduped


def test_parse_mutations_synonymous_dropped(tmp_path: Path):
    """Synonymous SNPs (G141G) in mutations.tsv should be filtered out."""
    main = tmp_path / "main.tsv"
    mut = tmp_path / "mutations.tsv"
    _write_tsv(main, [])
    _write_tsv(mut, [
        {"Gene symbol": "ompF_G141G", "Class": "QUINOLONE", "Method": "POINTX",
         "Subclass": "", "Scope": "", "Element type": "AMR",
         "% Identity to reference sequence": "", "% Coverage of reference sequence": ""},
    ])

    result = _parse_amrfinder_outputs(main, mut)
    # G141G is synonymous; should be dropped
    g141g_hits = [h for h in result["hits"] if h["symbol"] == "ompF_G141G"]
    assert len(g141g_hits) == 0


def test_parse_mutations_non_cipro_class_dropped(tmp_path: Path):
    """Mutations with non-CIPRO_RELEVANT class (e.g., TETRACYCLINE) dropped."""
    main = tmp_path / "main.tsv"
    mut = tmp_path / "mutations.tsv"
    _write_tsv(main, [])
    _write_tsv(mut, [
        {"Gene symbol": "acrB_R620R", "Class": "TETRACYCLINE", "Method": "POINTX",
         "Subclass": "", "Scope": "", "Element type": "AMR",
         "% Identity to reference sequence": "", "% Coverage of reference sequence": ""},
    ])

    result = _parse_amrfinder_outputs(main, mut)
    acrB_hits = [h for h in result["hits"] if h["symbol"] == "acrB_R620R"]
    assert len(acrB_hits) == 0  # TETRACYCLINE not in CIPRO_RELEVANT_AMR_CLASSES


def test_parse_mutations_multidrug_regulatory_kept(tmp_path: Path):
    """marR_V84WfsTer (MULTIDRUG class, regulatory mech) MUST be kept."""
    main = tmp_path / "main.tsv"
    mut = tmp_path / "mutations.tsv"
    _write_tsv(main, [])
    _write_tsv(mut, [
        {"Gene symbol": "marR_V84WfsTer14", "Class": "MULTIDRUG", "Method": "POINTX",
         "Subclass": "", "Scope": "", "Element type": "AMR",
         "% Identity to reference sequence": "", "% Coverage of reference sequence": ""},
    ])

    result = _parse_amrfinder_outputs(main, mut)
    marR_hits = [h for h in result["hits"] if h["symbol"] == "marR_V84WfsTer14"]
    assert len(marR_hits) == 1
    assert marR_hits[0]["mechanism"] == "regulatory"


def test_parse_aggregates_mechanisms_present(tmp_path: Path):
    main = tmp_path / "main.tsv"
    mut = tmp_path / "mutations.tsv"
    _write_tsv(main, [
        {"Gene symbol": "gyrA_S83L", "Class": "QUINOLONE", "Method": "POINTX",
         "Subclass": "", "Scope": "", "Element type": "AMR",
         "% Identity to reference sequence": "", "% Coverage of reference sequence": ""},
        {"Gene symbol": "qnrB19", "Class": "QUINOLONE", "Method": "ALLELEX",
         "Subclass": "", "Scope": "", "Element type": "AMR",
         "% Identity to reference sequence": "", "% Coverage of reference sequence": ""},
    ])
    _write_tsv(mut, [])

    result = _parse_amrfinder_outputs(main, mut)
    assert "QRDR_target_alteration" in result["mechanisms_present"]
    assert "plasmid_protect_modify" in result["mechanisms_present"]


def test_parse_empty_tsvs_yields_no_mechanism(tmp_path: Path):
    main = tmp_path / "main.tsv"
    mut = tmp_path / "mutations.tsv"
    main.write_text("", encoding="utf-8")
    mut.write_text("", encoding="utf-8")

    result = _parse_amrfinder_outputs(main, mut)
    assert result["n_hits"] == 0
    assert result["mechanisms_present"] == []
    assert result["primary_mechanism_class"] == "NO_MECHANISM"
