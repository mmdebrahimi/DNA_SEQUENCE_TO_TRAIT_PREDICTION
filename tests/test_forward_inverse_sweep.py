"""Tests for the cross-protein inverse boundary map.

Pins the two decompositions that make the sweep honest: the Q1/Q2 split (does inverse design work at all
vs does the LEARNED oracle earn its keep), and the PAIRED per-split test.
"""
from __future__ import annotations

import pytest

from scripts.forward_inverse_sweep import ASSAYS, MATERIAL_MARGIN, build_candidates


def test_assay_table_spans_multiple_kingdoms():
    """A one-organism sweep cannot find a boundary. Pin the tree-of-life spread."""
    orgs = {a["organism"] for a in ASSAYS}
    assert {"E. coli", "human", "yeast", "Arabidopsis"} <= orgs


def test_material_margin_is_shared_with_the_blatem_falsifier_not_retuned():
    """A per-protein bar would let the sweep manufacture whatever verdict it likes."""
    assert MATERIAL_MARGIN == 0.25


def test_candidate_space_is_labelled_and_differs_with_and_without_a_cds():
    """Only blaTEM has a committed CDS. The two candidate spaces answer DIFFERENT questions
    (genome-edit-reachable vs protein-level), so the space must be reported, never implied."""
    target = "MK"
    dms = {"M1L": -1.0, "M1W": -2.0, "K2R": -0.5}
    with_cds, space_a = build_candidates(target, dms, "ATGAAA")
    no_cds, space_b = build_candidates(target, dms, None)
    assert space_a == "single_nt_accessible_from_real_CDS"
    assert space_b == "all_DMS_variants_protein_level_no_CDS"
    # M1W needs 2 nt changes from ATG -> reachable only in the unrestricted space
    assert "M1W" not in {c.mutant for c in with_cds}
    assert "M1W" in {c.mutant for c in no_cds}


def test_build_candidates_drops_variants_that_disagree_with_the_reference():
    """Coordinate guard: a mutant whose stated WT is not the reference residue is refused, not coerced."""
    cands, _ = build_candidates("MK", {"A1L": -1.0, "M1L": -1.0}, None)
    assert {c.mutant for c in cands} == {"M1L"}


def test_build_candidates_excludes_nonsense():
    cands, _ = build_candidates("MK", {"M1*": -3.0, "M1L": -1.0}, None)
    assert {c.mutant for c in cands} == {"M1L"}


@pytest.mark.parametrize("field", ["beats_null", "beats_blosum", "paired_vs_null", "paired_vs_blosum",
                                   "esm2_err_over_span", "magnitude_certifiable", "forward_rank_on_pool"])
def test_scored_rows_carry_the_decomposed_verdict_fields(field):
    """REGRESSION PIN. A single collapsed verdict MISLABELLED CcdB as 'no discriminating power' when its
    inverse in fact beat the null by +77% -- it merely failed to beat BLOSUM. Q1 and Q2 must stay split,
    and the cross-protein error must stay span-normalized."""
    import json
    from pathlib import Path

    art = Path(__file__).resolve().parent.parent / "wiki" / "forward_inverse_sweep_2026-07-17.json"
    if not art.exists():
        pytest.skip("sweep artifact not present")
    scored = [r for r in json.loads(art.read_text(encoding="utf-8"))["assays"] if r["status"] == "SCORED"]
    assert scored
    for r in scored:
        assert field in r, f"{r['dms_id']} missing {field}"


def test_the_committed_sweep_records_the_rank_falsification():
    """The pre-registered question's answer: PTEN (0.5185) and RL40A (0.5190) have near-identical forward
    rank ON THE SAME POOL and OPPOSITE Q2 verdicts -> inverse utility does not track rank.

    HISTORY (kept deliberately): the first version of this finding used CcdB as the counterexample. CcdB
    turned out to be a CENSORED assay (79.3% of variants tied at the -2.00 ceiling) and is now excluded,
    so that pairing was an artifact. The finding SURVIVED re-anchoring to RL40A -- an even tighter pair
    (0.0005 apart) on a non-degenerate assay. Pin it so a re-run that quietly loses it is caught."""
    import json
    from pathlib import Path

    art = Path(__file__).resolve().parent.parent / "wiki" / "forward_inverse_sweep_2026-07-17.json"
    if not art.exists():
        pytest.skip("sweep artifact not present")
    rows = {r["dms_id"]: r for r in json.loads(art.read_text(encoding="utf-8"))["assays"]}
    pten, rl40a = rows["PTEN_HUMAN_Mighell_2018"], rows["RL40A_YEAST_Mavor_2016"]
    assert pten["status"] == rl40a["status"] == "SCORED"
    # near-identical forward rank, measured on each sweep's own pool...
    assert abs(pten["forward_rank_on_pool"] - rl40a["forward_rank_on_pool"]) < 0.01
    # ...opposite utility verdict
    assert pten["beats_blosum"] is True
    assert rl40a["beats_blosum"] is False
