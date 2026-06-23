"""Offline test for the TB lineage-collapse post-processor pure helper (_pass_mask faithfulness)."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from scripts.tb_independent_lineage_collapse import _pass_mask  # noqa: E402


def test_pass_mask_filter_only():
    vcf = ("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n"
           "NC_000962.3\t2155168\t.\tC\tG\t60\t.\tQNAME=x\tGT\t1/1\n")
    body = [l for l in _pass_mask(vcf).splitlines() if not l.startswith("#")][0].split("\t")
    assert body[6] == "PASS" and body[1] == "2155168" and body[4] == "G"   # FILTER->PASS; rest intact


if __name__ == "__main__":
    import pytest
    sys.exit(pytest.main([__file__, "-q"]))
