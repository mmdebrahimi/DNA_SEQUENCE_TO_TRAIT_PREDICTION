"""Offline test for the AMR Portal TB-cohort pivot (per-isolate RIF/INH with discordant-repeat blanking)."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from scripts.amr_portal_tb_cohort import pivot_tb  # noqa: E402


def test_pivot_tb_per_isolate_and_discordant():
    rows = [
        ("SAMEA1", "ERS1", "GCA_1", "rifampin", "resistant"),
        ("SAMEA1", "ERS1", "GCA_1", "isoniazid", "susceptible"),
        ("SAMEA2", "ERS2", "GCA_2", "rifampicin", "susceptible"),
        ("SAMEA2", "ERS2", "GCA_2", "rifampicin", "resistant"),   # discordant repeat -> blank rif
        ("SAMEA3", "", "", "ethambutol", "resistant"),            # non-TB-drug -> ignored
    ]
    iso = pivot_tb(rows)
    assert iso["SAMEA1"]["rif"] == "R" and iso["SAMEA1"]["inh"] == "S"
    assert iso["SAMEA1"]["assembly"] == "GCA_1"
    assert iso["SAMEA2"]["rif"] == ""                              # discordant -> blanked
    assert "SAMEA3" not in iso                                     # no TB-drug label


if __name__ == "__main__":
    import pytest
    sys.exit(pytest.main([__file__, "-q"]))
