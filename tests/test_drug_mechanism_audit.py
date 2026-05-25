"""Tests for scripts/drug_mechanism_audit.py pure-logic helpers.

AMRFinder + Docker not available locally — these tests pin the pure-logic
parser + verdict computation against synthetic TSVs.
"""
from __future__ import annotations

from pathlib import Path

import pytest


# ---- _is_synonymous_point ----


def test_is_synonymous_point_detects_synonymous():
    from scripts.drug_mechanism_audit import _is_synonymous_point
    assert _is_synonymous_point("gyrA_G141G") is True
    assert _is_synonymous_point("parC_S80S") is True


def test_is_synonymous_point_rejects_nonsynonymous():
    from scripts.drug_mechanism_audit import _is_synonymous_point
    assert _is_synonymous_point("gyrA_S83L") is False
    assert _is_synonymous_point("parC_S80I") is False


def test_is_synonymous_point_handles_unparseable():
    from scripts.drug_mechanism_audit import _is_synonymous_point
    assert _is_synonymous_point("gyrA") is False
    assert _is_synonymous_point("") is False
    assert _is_synonymous_point("no_underscore_format") is False
    assert _is_synonymous_point("gyrA_") is False


# ---- compute_verdict ----


def test_compute_verdict_primary_dominant_for_cipro():
    """≥ 70% of R have a primary mechanism (QRDR_target_alteration) -> PRIMARY_DOMINANT."""
    from scripts.drug_mechanism_audit import compute_verdict
    r_strains = [
        {"mechanisms_present": ["QRDR_target_alteration"]},
        {"mechanisms_present": ["QRDR_target_alteration", "plasmid_protect_modify"]},
        {"mechanisms_present": ["QRDR_target_alteration"]},
        {"mechanisms_present": ["efflux"]},  # no primary
    ]
    verdict, breakdown = compute_verdict(r_strains, drug="ciprofloxacin")
    assert verdict == "PRIMARY_DOMINANT"
    assert breakdown["n_R_strains"] == 4
    assert breakdown["n_R_with_primary_mechanism"] == 3


def test_compute_verdict_mixed_when_below_primary_bar():
    """50-69% any mechanism -> MIXED."""
    from scripts.drug_mechanism_audit import compute_verdict
    r_strains = [
        {"mechanisms_present": ["QRDR_target_alteration"]},
        {"mechanisms_present": ["efflux"]},
        {"mechanisms_present": ["regulatory"]},
        {"mechanisms_present": []},
        {"mechanisms_present": []},
    ]
    verdict, breakdown = compute_verdict(r_strains, drug="ciprofloxacin")
    assert verdict == "MIXED_MECHANISMS"
    assert breakdown["n_R_with_primary_mechanism"] == 1  # 20%, below 70%
    assert breakdown["n_R_with_any_mechanism"] == 3  # 60%


def test_compute_verdict_mostly_unknown_when_below_mixed_bar():
    """< 50% any mechanism -> MOSTLY_UNKNOWN."""
    from scripts.drug_mechanism_audit import compute_verdict
    r_strains = [
        {"mechanisms_present": ["QRDR_target_alteration"]},
        {"mechanisms_present": []},
        {"mechanisms_present": []},
        {"mechanisms_present": []},
    ]
    verdict, _ = compute_verdict(r_strains, drug="ciprofloxacin")
    assert verdict == "MOSTLY_UNKNOWN"


def test_compute_verdict_empty_r_set():
    from scripts.drug_mechanism_audit import compute_verdict
    verdict, breakdown = compute_verdict([], drug="ciprofloxacin")
    assert verdict == "EMPTY_R_SET"
    assert breakdown == {}


# ---- parse_amrfinder_outputs_for_drug ----


def _write_tsv(path: Path, headers: list[str], rows: list[list[str]]) -> None:
    """Write a tab-separated TSV with the given headers + rows."""
    lines = ["\t".join(headers)]
    for r in rows:
        lines.append("\t".join(r))
    path.write_text("\n".join(lines), encoding="utf-8")


def test_parse_amrfinder_recovers_cipro_qrdr_mutations(tmp_path: Path):
    """gyrA + parC POINT mutations should classify under QRDR_target_alteration."""
    from scripts.drug_mechanism_audit import parse_amrfinder_outputs_for_drug
    main_tsv = tmp_path / "main.tsv"
    mut_tsv = tmp_path / "mutations.tsv"
    headers = [
        "Gene symbol", "Class", "Subclass", "Method", "Scope",
        "Element type", "% Identity to reference sequence", "% Coverage of reference sequence",
    ]
    _write_tsv(main_tsv, headers, [
        ["gyrA_S83L", "QUINOLONE", "QUINOLONE", "POINTX", "core", "POINT", "100.0", "100.0"],
        ["parC_S80I", "QUINOLONE", "QUINOLONE", "POINTX", "core", "POINT", "100.0", "100.0"],
    ])
    _write_tsv(mut_tsv, headers, [])
    parsed = parse_amrfinder_outputs_for_drug(main_tsv, mut_tsv, drug="ciprofloxacin")
    assert "QRDR_target_alteration" in parsed["mechanisms_present"]
    assert parsed["primary_mechanism_class"] == "QRDR_target_alteration"
    assert parsed["n_hits"] == 2


def test_parse_amrfinder_filters_synonymous(tmp_path: Path):
    """Synonymous POINT mutations should not appear."""
    from scripts.drug_mechanism_audit import parse_amrfinder_outputs_for_drug
    main_tsv = tmp_path / "main.tsv"
    mut_tsv = tmp_path / "mutations.tsv"
    headers = [
        "Gene symbol", "Class", "Subclass", "Method", "Scope",
        "Element type", "% Identity to reference sequence", "% Coverage of reference sequence",
    ]
    _write_tsv(main_tsv, headers, [
        ["gyrA_G141G", "QUINOLONE", "QUINOLONE", "POINTX", "core", "POINT", "100.0", "100.0"],
    ])
    _write_tsv(mut_tsv, headers, [])
    parsed = parse_amrfinder_outputs_for_drug(main_tsv, mut_tsv, drug="ciprofloxacin")
    assert parsed["n_hits"] == 0
    assert parsed["primary_mechanism_class"] == "NO_MECHANISM"


def test_parse_amrfinder_filters_non_cipro_classes(tmp_path: Path):
    """TETRACYCLINE class mutation should not appear in cipro audit."""
    from scripts.drug_mechanism_audit import parse_amrfinder_outputs_for_drug
    main_tsv = tmp_path / "main.tsv"
    mut_tsv = tmp_path / "mutations.tsv"
    headers = [
        "Gene symbol", "Class", "Subclass", "Method", "Scope",
        "Element type", "% Identity to reference sequence", "% Coverage of reference sequence",
    ]
    _write_tsv(main_tsv, headers, [
        ["tetA", "TETRACYCLINE", "TETRACYCLINE", "EXACTX", "core", "AMR", "100.0", "100.0"],
    ])
    _write_tsv(mut_tsv, headers, [])
    parsed = parse_amrfinder_outputs_for_drug(main_tsv, mut_tsv, drug="ciprofloxacin")
    assert parsed["n_hits"] == 0


def test_parse_amrfinder_picks_up_tet_when_drug_is_tet(tmp_path: Path):
    """tetA hit in TETRACYCLINE class should classify when drug=tetracycline."""
    from scripts.drug_mechanism_audit import parse_amrfinder_outputs_for_drug
    main_tsv = tmp_path / "main.tsv"
    mut_tsv = tmp_path / "mutations.tsv"
    headers = [
        "Gene symbol", "Class", "Subclass", "Method", "Scope",
        "Element type", "% Identity to reference sequence", "% Coverage of reference sequence",
    ]
    _write_tsv(main_tsv, headers, [
        ["tetA", "TETRACYCLINE", "TETRACYCLINE", "EXACTX", "core", "AMR", "100.0", "100.0"],
    ])
    _write_tsv(mut_tsv, headers, [])
    parsed = parse_amrfinder_outputs_for_drug(main_tsv, mut_tsv, drug="tetracycline")
    # tetA should classify under efflux per mic_tiers (or another tet mechanism)
    assert parsed["n_hits"] == 1
    assert parsed["mechanisms_present"]  # non-empty


def test_parse_amrfinder_dedupes_pointx_across_main_and_mutations(tmp_path: Path):
    """gyrA_S83L appearing in BOTH main.tsv (POINTX) AND mutations.tsv -> single hit."""
    from scripts.drug_mechanism_audit import parse_amrfinder_outputs_for_drug
    main_tsv = tmp_path / "main.tsv"
    mut_tsv = tmp_path / "mutations.tsv"
    headers = [
        "Gene symbol", "Class", "Subclass", "Method", "Scope",
        "Element type", "% Identity to reference sequence", "% Coverage of reference sequence",
    ]
    _write_tsv(main_tsv, headers, [
        ["gyrA_S83L", "QUINOLONE", "QUINOLONE", "POINTX", "core", "POINT", "100.0", "100.0"],
    ])
    _write_tsv(mut_tsv, headers, [
        ["gyrA_S83L", "QUINOLONE", "QUINOLONE", "POINT", "core", "POINT", "100.0", "100.0"],
    ])
    parsed = parse_amrfinder_outputs_for_drug(main_tsv, mut_tsv, drug="ciprofloxacin")
    assert parsed["n_hits"] == 1  # deduped
    assert parsed["mechanisms_present"] == ["QRDR_target_alteration"]


def test_parse_amrfinder_handles_missing_files(tmp_path: Path):
    """Both TSVs missing -> empty result (don't crash)."""
    from scripts.drug_mechanism_audit import parse_amrfinder_outputs_for_drug
    parsed = parse_amrfinder_outputs_for_drug(
        tmp_path / "missing_main.tsv",
        tmp_path / "missing_mut.tsv",
        drug="ciprofloxacin",
    )
    assert parsed["n_hits"] == 0
    assert parsed["mechanisms_present"] == []
    assert parsed["primary_mechanism_class"] == "NO_MECHANISM"


def test_parse_amrfinder_recovers_cef_beta_lactamase(tmp_path: Path):
    """CTX-M acquired-gene hit in BETA-LACTAM class should classify under cef."""
    from scripts.drug_mechanism_audit import parse_amrfinder_outputs_for_drug
    main_tsv = tmp_path / "main.tsv"
    mut_tsv = tmp_path / "mutations.tsv"
    headers = [
        "Gene symbol", "Class", "Subclass", "Method", "Scope",
        "Element type", "% Identity to reference sequence", "% Coverage of reference sequence",
    ]
    _write_tsv(main_tsv, headers, [
        ["blaCTX-M-15", "BETA-LACTAM", "CEPHALOSPORIN", "EXACTX", "core", "AMR", "100.0", "100.0"],
    ])
    _write_tsv(mut_tsv, headers, [])
    parsed = parse_amrfinder_outputs_for_drug(main_tsv, mut_tsv, drug="ceftriaxone")
    assert parsed["n_hits"] >= 1
    assert parsed["mechanisms_present"]  # non-empty
