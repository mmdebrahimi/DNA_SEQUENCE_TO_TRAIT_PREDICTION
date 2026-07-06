"""Offline pins for the PGP-UK PGx cohort-card builder (`scripts/pgx_pgp_uk_cohort_card.py`)."""
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "scripts"))

import pgx_pgp_uk_cohort_card as m  # noqa: E402


def _mk(sample, cyp2c19_dp, cyp2c19_ab, vkorc1_gt):
    return {"sample_id": sample, "results": {
        "cyp2c19": {"diplotype": cyp2c19_dp, "phenotype_abbrev": cyp2c19_ab},
        "cyp2c9": {"diplotype": "*1/*1", "phenotype_abbrev": "NM"},
        "cyp3a5": {"diplotype": "*1/*1", "phenotype_abbrev": "NM"},
        "tpmt": {"diplotype": "*1/*1", "phenotype_abbrev": "NM"},
        "cyp2b6": {"diplotype": "*1/*1", "phenotype_abbrev": "NM"},
        "vkorc1": {"genotype": vkorc1_gt}, "slco1b1": {"genotype": "T/T"}}}


def test_cohort_distribution_and_allele_counts():
    results = [_mk("A", "*1/*2", "IM", "G/G"), _mk("B", "*1/*17", "RM", "A/A"), _mk("C", "*1/*1", "NM", "G/A")]
    card = m.build_card(results)
    assert card["n_individuals"] == 3
    c19 = card["per_gene"]["cyp2c19"]
    assert c19["phenotypes"] == {"IM": 1, "RM": 1, "NM": 1}
    # observed allele counts: *1 appears in all 3 diplotypes' first slot + *1/*1 twice -> *1 x4, *2 x1, *17 x1
    assert c19["observed_allele_counts"]["*1"] == 4
    assert c19["observed_allele_counts"]["*2"] == 1 and c19["observed_allele_counts"]["*17"] == 1
    assert card["per_gene"]["vkorc1"]["genotypes"] == {"G/G": 1, "A/A": 1, "G/A": 1}


def test_render_md_has_per_individual_rows():
    results = [_mk("A", "*1/*2", "IM", "G/G")]
    for r in results:
        m._RESULTS_CACHE[r["sample_id"]] = r
    md = m.render_md(m.build_card(results))
    assert "A" in md and "*1/*2 IM" in md and "N = 1" in md


if __name__ == "__main__":
    import pytest
    raise SystemExit(pytest.main([__file__, "-q"]))
