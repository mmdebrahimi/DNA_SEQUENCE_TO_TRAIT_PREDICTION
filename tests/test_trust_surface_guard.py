"""Pins the trust-surface NAMESPACE GUARD (the productization-move-3 hardening): an HIV/TB/SARS card lends
its evidence ONLY when the requested organism is absent (drug-only) or normalizes into that namespace. A
contradictory organism is REFUSED (tier UNKNOWN + reason='namespace_mismatch'), with the rejected candidate
exposed in evidence_cell so the borrowing is auditable. Closes the fail-open where e.g. rifampicin on an
E. coli call would wear TB's independence badge.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from dna_decode.data import trust_surface as ts  # noqa: E402


def test_contradictory_organism_is_refused_not_borrowed():
    # rifampicin is a TB cell; on an E. coli call it must NOT lend TB's INDEPENDENT_MEASURED badge
    b = ts.lookup_trust("rifampicin", "Escherichia")
    assert b["tier"] == ts.UNKNOWN
    assert b["reason"] == "namespace_mismatch"
    assert b["independent"] is False
    assert b["metric"] is None
    assert b["evidence_cell"] == "M.tuberculosis|rifampicin"   # the rejected candidate is auditable
    assert b["requested_cell"].lower().startswith("escherichia")


def test_drug_only_lookup_still_allowed_when_organism_absent():
    # absent organism (drug-only) is DIFFERENT from a contradictory one -> the namespace may lend evidence
    assert ts.lookup_trust("rifampicin", None)["tier"] == ts.INDEPENDENT_MEASURED
    assert ts.lookup_trust("efavirenz", None)["tier"] == ts.INDEPENDENT_WETLAB
    assert ts.lookup_trust("nirmatrelvir", None)["tier"] == ts.IN_DISTRIBUTION


def test_compatible_organism_fires_the_namespace():
    assert ts.lookup_trust("rifampicin", "Mycobacterium_tuberculosis")["tier"] == ts.INDEPENDENT_MEASURED
    assert ts.lookup_trust("efavirenz", "HIV-1")["tier"] == ts.INDEPENDENT_WETLAB


def test_hiv_and_sars_mismatch_on_bacterial_organism():
    for drug, ev in [("efavirenz", "HIV-1|efavirenz"), ("nirmatrelvir", "SARS-CoV-2|nirmatrelvir")]:
        b = ts.lookup_trust(drug, "Klebsiella")
        assert b["tier"] == ts.UNKNOWN and b["reason"] == "namespace_mismatch"
        assert b["evidence_cell"] == ev


def test_bacterial_match_unchanged_and_carries_cells():
    b = ts.lookup_trust("ciprofloxacin", "Escherichia")
    assert b["tier"] == ts.INDEPENDENT_MEASURED and b["reason"] is None
    assert b["requested_cell"] and b["evidence_cell"]   # both present for a legitimate match


def test_unknown_drug_has_no_evidence_cell():
    b = ts.lookup_trust("totally_made_up", "Escherichia")
    assert b["tier"] == ts.UNKNOWN and b["reason"] is None and b["evidence_cell"] is None


if __name__ == "__main__":
    import pytest
    raise SystemExit(pytest.main([__file__, "-q"]))
