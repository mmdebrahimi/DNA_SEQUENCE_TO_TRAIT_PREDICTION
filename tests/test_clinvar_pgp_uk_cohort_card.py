"""Offline pins for the ClinVar PGP-UK cohort-card builder (`scripts/clinvar_pgp_uk_cohort_card.py`)."""
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "scripts"))

import clinvar_pgp_uk_cohort_card as m  # noqa: E402


def _mk(sample, n_path, n_benign, benign_by_gene, path_hits=None):
    return {"sample_id": sample, "n_pathogenic": n_path, "n_benign": n_benign,
            "benign_by_gene": benign_by_gene, "pathogenic_hits": path_hits or []}


def test_cohort_aggregates_benign_and_pathogenic():
    results = [_mk("A", 0, 100, {"BRCA1": 60, "TTN": 40}),
               _mk("B", 1, 80, {"BRCA1": 50, "RYR2": 30},
                   [{"gene": "MSH2", "chrom": "2", "pos": "47", "ref": "A", "alt": "T",
                     "significance": "Pathogenic", "stars": 3, "disease": "Lynch"}])]
    card = m.build_card(results, "test panel")
    assert card["n_individuals"] == 2
    assert card["cohort_n_pathogenic"] == 1
    assert card["benign_carrier_load_range"] == [80, 100]
    # cohort per-gene benign = summed across samples
    assert card["top_genes_by_benign_coverage"]["BRCA1"] == 110
    assert len(card["cohort_pathogenic_hits"]) == 1 and card["cohort_pathogenic_hits"][0]["gene"] == "MSH2"


def test_render_md_surfaces_pathogenic_detail():
    results = [_mk("A", 1, 10, {"BRCA1": 10},
                   [{"gene": "MSH2", "chrom": "2", "pos": "47", "ref": "A", "alt": "T",
                     "significance": "Pathogenic", "stars": 3, "disease": "Lynch"}])]
    md = m.render_md(m.build_card(results, "p"))
    assert "MSH2" in md and "Lynch" in md and "N = 1" in md


def test_render_md_zero_pathogenic_states_expected():
    results = [_mk("A", 0, 500, {"TTN": 500})]
    md = m.render_md(m.build_card(results, "p"))
    assert "reportable-pathogenic findings: 0" in md and "expected" in md.lower()


if __name__ == "__main__":
    import pytest
    raise SystemExit(pytest.main([__file__, "-q"]))
