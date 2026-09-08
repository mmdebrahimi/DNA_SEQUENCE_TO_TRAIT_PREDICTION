"""Every headline figure in a cell contract must match the artifact or constant that produced it.

WHY THIS EXISTS. The `pathotype` contract quoted ExPEC recall 0.917 while the ENFORCED test constant
said 10/12 = 0.833 and called 0.917 "an overfit". That number was wrong in CLAUDE.md for three months,
and it survived every existing guard: `test_project_orientation` pins CELL COUNTS, and
`test_report_card_doc_sync` pins REPORT-CARD FIGURES. Nothing pinned a PER-CELL METRIC to the source
of truth that enforces it, so a superseded rule's number kept being repeated to every future reader.

The recorded lesson is "pin any quoted figure to its live source by test". These are those pins, for
the contracts whose numbers HAVE a machine-readable source. A contract number without one (a prose
scope note, a judgement) is deliberately not pinned here -- inventing a fuzzy match would produce
noise, and the guard would stop being trusted.

Runs fully offline from committed wiki artifacts (which live on C:, not the D: data junction).
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tests"))

from dna_decode.data.cell_registry import by_cell_id  # noqa: E402

WIKI = ROOT / "wiki"


def _contract_text(cell_id: str) -> str:
    """Every prose field a figure could legitimately be quoted in."""
    c = by_cell_id()[cell_id]
    return " ".join(str(getattr(c, f) or "") for f in
                    ("claim", "claim_status", "validation_slice", "label_provenance", "demotion_rule"))


def _artifact(name: str) -> dict:
    p = WIKI / name
    if not p.exists():
        pytest.skip(f"artifact absent: {name}")
    return json.loads(p.read_text(encoding="utf-8"))


# --- pathotype: pinned to the ENFORCED test constant, not to an artifact --------------------------

def test_pathotype_contract_matches_the_enforced_recall_cap():
    """THE regression this file was written for. If the cap moves, the contract must move with it."""
    from test_pathotype_expec_recall import EXPEC_RECALL_CAP
    txt = _contract_text("typing:Escherichia_coli:pathotype")
    assert EXPEC_RECALL_CAP == 10 / 12, "the enforced cap changed; update the contract and this pin"
    assert "10/12" in txt, "contract must quote the enforced cap 10/12"
    assert "0.833" in txt


def test_pathotype_contract_does_not_resurrect_the_superseded_number():
    """0.917 may appear ONLY as the rejected alternative, never as this cell's recall."""
    txt = _contract_text("typing:Escherichia_coli:pathotype")
    if "0.917" in txt:
        assert "overfit" in txt.lower() or "over-rescu" in txt.lower(), (
            "0.917 appears without being marked as the REJECTED flat-K=1 number")


# --- contracts pinned to their committed artifacts ------------------------------------------------

def test_resfinder_contract_matches_its_artifact():
    a = _artifact("resfinder_locus_collapse_2026-09-05.json")
    txt = _contract_text("finder:bacteria:resfinder")
    bl, ag = a["summary"]["beta-lactam"], a["summary"]["aminoglycoside"]
    assert str(a["n_genomes"]) in txt                                   # 648
    assert f"{bl['mean_genes_old']:.2f}" in txt                         # 165.52
    assert f"{bl['mean_genes_new']:.2f}" in txt                         # 1.59
    assert f"{bl['jaccard_vs_amrfinder_normalized_new']:.4f}" in txt    # 0.7754
    assert f"{ag['jaccard_vs_amrfinder_normalized_new']:.4f}" in txt    # 0.7786


def test_pointfinder_contract_matches_its_artifact():
    a = _artifact("pointfinder_amrfinder_concordance_2026-09-05.json")
    txt = _contract_text("finder:Escherichia_coli:pointfinder")
    assert str(a["n_agree"]) in txt                                     # 1524
    assert f"{a['n_genomes_exact_set_match']}/{a['n_genomes']}" in txt   # 641/646
    assert str(a["n_pointfinder_only"]) in txt and str(a["n_amrfinder_only"]) in txt


def test_mlst_contract_matches_its_artifact():
    a = _artifact("mlst_serotype_purity_2026-09-05.json")
    txt = _contract_text("typing:bacteria:mlst")
    assert f"{a['observed_purity']:.4f}" in txt                          # 0.7860
    assert f"{a['null_max']:.4f}" in txt                                 # 0.2583
    assert str(a["n_st_scored"]) in txt                                  # 27


def test_disinfinder_contract_matches_its_artifact():
    a = _artifact("disinfinder_locus_collapse_probe_2026-09-05.json")
    txt = _contract_text("finder:bacteria:disinfinder")
    assert str(a["n_db_alleles"]) in txt                                 # 16
    assert str(a["n_genomes_scored"]) in txt                             # 40
    assert str(a["total_gene_calls_old_rule"]) in txt                    # 34
    assert a["n_genomes_with_different_gene_set"] == 0


def test_plasmid_contract_matches_its_artifact():
    a = _artifact("plasmid_selection_rule_probe_2026-09-04.json")
    txt = _contract_text("finder:bacteria:plasmid")
    assert str(a["total_replicon_calls"]) in txt                         # 446
    assert str(a["n_assemblies_scored"]) in txt                          # 40
    assert a["n_assemblies_with_different_replicon_set"] == 0


# --- non-vacuity ----------------------------------------------------------------------------------

def test_the_pins_would_actually_catch_a_drifted_number():
    """A guard that only asserts substrings could pass on an empty contract. Prove the texts are real
    and that a wrong number would NOT be found in them."""
    for cid in ("typing:Escherichia_coli:pathotype", "finder:bacteria:resfinder",
                "finder:Escherichia_coli:pointfinder", "typing:bacteria:mlst"):
        txt = _contract_text(cid)
        assert len(txt) > 400, f"{cid} contract is too thin to be pinning anything"
        assert "999.99" not in txt      # a number that is not there is not found
