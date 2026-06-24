"""Offline tests for the TB BioSample-independence check pure logic (composition + overlap)."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from scripts.tb_independence_biosample_check import composition, overlap, probe_ncbi_side  # noqa: E402


def test_composition_splits_ena_vs_ncbi_side():
    rows = [
        {"leaked": "0", "assembly": "GCA_1.1", "biosample": "SAMEA1", "sra": "ERS1"},   # ENA-side
        {"leaked": "0", "assembly": "GCA_2.1", "biosample": "SAMN2", "sra": "SRS2"},     # NCBI-side
        {"leaked": "1", "assembly": "GCA_3.1", "biosample": "SAMEA3", "sra": "ERS3"},    # leaked -> excluded
    ]
    c = composition(rows)
    assert c["ena_side_already_biosample_grade"] == 1
    assert c["ncbi_side_cross_archive_candidates"] == 1


def test_overlap_case_insensitive():
    assert overlap({"err999", "samn5"}, {"ERR999", "ERS1"}) == {"ERR999"}
    assert overlap({"samn5"}, {"ERR999"}) == set()


def test_probe_uses_injected_fetch_and_flags_overlap():
    rows = [{"leaked": "0", "assembly": "GCA_1.1", "biosample": "SAMN_clean", "sra": "SRS1"},
            {"leaked": "0", "assembly": "GCA_2.1", "biosample": "SAMN_leak", "sra": "SRS2"}]
    cryptic = {"ERR_LEAK"}
    def fake_fetch(url):
        if "SAMN_leak" in url:
            return "run_accession\tsample_accession\tsecondary_sample_accession\nERR_LEAK\tSAMN_leak\tERS_x\n"
        return "run_accession\tsample_accession\tsecondary_sample_accession\nSRR9\tSAMN_clean\t\n"
    p = probe_ncbi_side(rows, cryptic, sample=2, fetch=fake_fetch)
    assert p["n_probed"] == 2 and p["n_ena_mirrored"] == 2
    assert p["n_cross_archive_overlap"] == 1 and p["hits"][0]["biosample"] == "SAMN_leak"


if __name__ == "__main__":
    import pytest
    sys.exit(pytest.main([__file__, "-q"]))
