"""Offline pins for the captured summary-stat indexes + the GWAS single-SNP miner filter (no network)."""
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from scripts.mine_gwas_single_snp_candidates import OR_MAX, OR_MIN, P_MAX, _RSID  # noqa: E402

_DIR = Path(__file__).resolve().parent.parent / "data" / "summary_stat_sources"


def _rows(name, cols):
    p = _DIR / name
    if not p.exists():
        import pytest
        pytest.skip(f"{name} not captured yet")
    with p.open(encoding="utf-8") as f:
        rdr = csv.DictReader(f, delimiter="\t")
        assert all(c in rdr.fieldnames for c in cols)
        return list(rdr)


def test_pgs_index_captured():
    r = _rows("pgs_catalog_index.tsv", ["pgs_id", "trait", "n_variants"])
    assert len(r) > 3000 and r[0]["pgs_id"].startswith("PGS")


def test_panukbb_and_finngen_indexes_captured():
    assert len(_rows("pan_ukbb_phenotype_index.tsv", ["phenocode", "description"])) > 5000
    assert len(_rows("finngen_r12_endpoint_index.tsv", ["phenocode", "phenotype"])) > 2000


def test_gwas_candidates_are_clean_single_snp_strong():
    r = _rows("gwas_single_snp_candidates.tsv", ["rsid", "trait", "or_or_beta", "p_value"])
    assert 0 < len(r) <= 400
    for row in r:
        assert _RSID.match(row["rsid"])                          # single rsID
        assert OR_MIN <= float(row["or_or_beta"]) <= OR_MAX      # OR-artifact cap held (verify-in-batch fix)
        assert float(row["p_value"]) <= P_MAX                    # genome-wide significant


def test_filter_constants():
    assert P_MAX == 5e-8 and OR_MIN == 2.0 and OR_MAX == 20.0
