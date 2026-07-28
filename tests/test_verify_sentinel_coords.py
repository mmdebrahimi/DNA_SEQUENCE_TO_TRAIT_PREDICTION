"""Offline tests for the sentinel-coord verifier (anti-fabrication rail). Mock the Ensembl fetch -- NO network.
Runnable via pytest OR standalone."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from scripts.verify_sentinel_coords import verify_sentinel  # noqa: E402


def _mapping(chrom, start, allele_string, assembly="GRCh38", coord="chromosome"):
    return {"assembly_name": assembly, "coord_system": coord, "seq_region_name": chrom,
            "start": start, "end": start, "allele_string": allele_string}


def _resp(name, mappings):
    return {"name": name, "mappings": mappings}


def test_ok_match():
    f = lambda rs: _resp("rs28371686", [_mapping("10", 94981301, "C/A")])
    r = verify_sentinel("rs28371686", "10", 94981301, "C", "A", fetch=f)
    assert r.ok and r.status == "OK"


def test_pos_mismatch_fails():
    f = lambda rs: _resp("rsX", [_mapping("10", 99999999, "C/A")])
    r = verify_sentinel("rsX", "10", 94981301, "C", "A", fetch=f)
    assert r.status == "MISMATCH" and "pos" in r.detail


def test_ref_mismatch_fails():
    f = lambda rs: _resp("rsX", [_mapping("10", 94981301, "G/A")])   # ref G != expected C
    r = verify_sentinel("rsX", "10", 94981301, "C", "A", fetch=f)
    assert r.status == "MISMATCH" and "ref" in r.detail


def test_no_grch38_chromosome_mapping_fails():
    # only a GRCh37 mapping + a scaffold -> no qualifying GRCh38 chromosome mapping
    f = lambda rs: _resp("rsX", [_mapping("10", 94981301, "C/A", assembly="GRCh37"),
                                 _mapping("HSCHR10", 1, "C/A", coord="scaffold")])
    r = verify_sentinel("rsX", "10", 94981301, "C", "A", fetch=f)
    assert r.status == "MISMATCH" and "no GRCh38" in r.detail


def test_multi_mapping_fails():
    f = lambda rs: _resp("rsX", [_mapping("10", 94981301, "C/A"), _mapping("10", 94981302, "C/A")])
    r = verify_sentinel("rsX", "10", 94981301, "C", "A", fetch=f)
    assert r.status == "MISMATCH" and "multi-map" in r.detail


def test_wrong_chrom_fails():
    f = lambda rs: _resp("rsX", [_mapping("6", 94981301, "C/A")])   # mapping on chr6, expected chr10
    r = verify_sentinel("rsX", "10", 94981301, "C", "A", fetch=f)
    assert r.status == "MISMATCH"


def test_merged_synonym_surfaced_but_ok():
    f = lambda rs: _resp("rs9999999", [_mapping("10", 94981301, "C/A")])  # current name != queried
    r = verify_sentinel("rsOLD", "10", 94981301, "C", "A", fetch=f)
    assert r.ok and r.merged_into == "rs9999999"


def test_alt_warning_when_alt_absent():
    f = lambda rs: _resp("rsX", [_mapping("10", 94981301, "C/A")])   # ALT set = {A}; expected alt T
    r = verify_sentinel("rsX", "10", 94981301, "C", "T", fetch=f)
    assert r.ok and r.alt_warning is True


def test_wildcard_alt_never_warns():
    f = lambda rs: _resp("rsX", [_mapping("10", 94981301, "C/A")])
    r = verify_sentinel("rsX", "10", 94981301, "C", "*", fetch=f)
    assert r.ok and r.alt_warning is False


def test_unreachable_is_unverified_not_silent_pass():
    def boom(rs):
        raise OSError("network down")
    r = verify_sentinel("rsX", "10", 94981301, "C", "A", fetch=boom)
    assert r.status == "UNVERIFIED" and not r.ok


def test_chr_prefix_normalized():
    f = lambda rs: _resp("rsX", [_mapping("chr10", 94981301, "C/A")])
    r = verify_sentinel("rsX", "chr10", 94981301, "C", "A", fetch=f)
    assert r.ok


if __name__ == "__main__":
    fns = [(k, v) for k, v in sorted(globals().items()) if k.startswith("test_") and callable(v)]
    for k, fn in fns:
        fn(); print(f"PASS {k}")
    print(f"\n{len(fns)} passed")
