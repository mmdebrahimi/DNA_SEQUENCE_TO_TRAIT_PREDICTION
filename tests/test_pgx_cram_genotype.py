"""Offline tests for the read-level PGx CRAM site-genotyper (pure parsing/calling logic; mock docker_run)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from scripts.pgx_cram_genotype import call_site, parse_mpileup, genotype_sites, sites_from_gene  # noqa: E402


def test_call_site_het():
    d, rc, ac, frac, gt, present = call_site("atTTTAAAtaTTAaTttA", "T", "A")  # mix of A + T
    assert gt == "0/1" and present and 0.15 <= frac <= 0.85 and d == rc + ac


def test_call_site_hom_alt():
    d, rc, ac, frac, gt, present = call_site("AAAAAAAAAAAA", "T", "A")
    assert gt == "1/1" and present and frac == 1.0


def test_call_site_hom_ref():
    d, rc, ac, frac, gt, present = call_site("TTTTTTTTTTTT", "T", "A")
    assert gt == "0/0" and not present and frac == 0.0


def test_call_site_uncallable_low_depth():
    d, rc, ac, frac, gt, present = call_site("AT", "T", "A")  # depth 2 < MIN_DEPTH
    assert gt == "UNCALLABLE" and not present


def test_parse_mpileup():
    raw = "chr6\t18133845\tN\t36\tatTTTAAA\tqual\nchr6\t18133846\tN\t10\t....\tq"
    p = parse_mpileup(raw)
    assert 18133845 in p and p[18133845][1] == "atTTTAAA"


def test_genotype_sites_with_mock():
    sites = [{"label": "*6", "chrom": "6", "pos": 18133845, "ref": "T", "alt": "A"}]

    class _Out:
        returncode = 0
        stdout = "chr6\t18133845\tN\t36\t" + "A" * 16 + "T" * 20 + "\tq"
        stderr = ""

    def fake_run(img, cmd, **kw):
        return _Out()

    calls = genotype_sites("NA18603", "http://x.cram", sites, docker_run=fake_run)
    assert len(calls) == 1 and calls[0].alt_present and calls[0].genotype == "0/1"
    assert calls[0].chrom == "chr6"  # normalized


def test_sites_from_gene_tpmt():
    sites = sites_from_gene("tpmt")
    assert len(sites) == 10 and all({"label", "chrom", "pos", "ref", "alt"} <= set(s) for s in sites)
    assert any(s["pos"] == 18133845 for s in sites)  # the *6 site


if __name__ == "__main__":
    fns = [(k, v) for k, v in sorted(globals().items()) if k.startswith("test_") and callable(v)]
    for k, fn in fns:
        fn(); print(f"PASS {k}")
    print(f"\n{len(fns)} passed")
