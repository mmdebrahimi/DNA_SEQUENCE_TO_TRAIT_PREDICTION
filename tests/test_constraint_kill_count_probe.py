"""Guards for the kill-count probe (plan Step 8, 2026-10-01)."""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
_spec = importlib.util.spec_from_file_location(
    "constraint_kill_count_probe", REPO / "scripts" / "constraint_kill_count_probe.py")
P = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(P)

ART = REPO / "wiki" / "constraint_kill_count_2026-10-01.json"


def _art():
    if not ART.exists():
        pytest.skip("probe artifact absent")
    return json.loads(ART.read_text(encoding="utf-8"))


# ---------------------------------------------------------------- offline unit coverage

def test_read_fasta_strips_headers_and_uppercases(tmp_path):
    f = tmp_path / "x.fna"
    f.write_text(">h\nacg\ntga\n", encoding="utf-8")
    assert P.read_fasta(f) == "ACGTGA"


def test_the_expectation_is_stated_BEFORE_the_run_and_the_correction_is_recorded():
    """A prediction written after seeing the result is not a prediction."""
    assert "predicted INERT" in P.PRE_STATED_EXPECTATION
    assert "NOT to be tuned until something fails" in P.PRE_STATED_EXPECTATION
    assert "refuted 5 of 8" in P.FIRST_RUN_CORRECTION
    assert "EXPECTATION was wrong" in P.FIRST_RUN_CORRECTION


def test_only_genuinely_complete_genes_are_declared_complete():
    """Asserting completeness for a polyprotein extract would manufacture a false refutation."""
    assert P.COMPLETE_GENES == {"cauris_erg11_wt", "cauris_erg11_y132f", "cauris_erg11_k143r"}
    for name in ("hiv1_rt_hxb2", "hiv1_pr_hxb2", "sarscov2_mpro"):
        assert name not in P.COMPLETE_GENES, f"{name} is a polyprotein extract"


def test_the_probe_refuses_rather_than_reporting_when_no_reference_is_readable(monkeypatch, tmp_path):
    monkeypatch.setattr(P, "REFERENCE_CDS", {"nope": "data/does_not_exist.fna"})
    monkeypatch.setattr(P, "MYCOPLASMA_FIXTURE", "data/also_missing.fna")
    assert P.main(["--out", str(tmp_path / "o.json")]) == 2


# ---------------------------------------------------------------- real-data assertions

def test_the_clade_demo_is_REAL_not_synthetic():
    """The plan allowed a synthetic fallback; the fetch succeeded, so this must say so honestly."""
    d = _art()["clade_demo"]
    assert d["source"] == "real" and d["label"] == "real_organism_cds"
    assert "Mycoplasma genitalium" in d["organism"]
    assert "L43967.2" in d["provenance"]
    assert (REPO / d["fixture"]).exists()


def test_the_clade_table_flips_the_verdict_on_the_real_cds():
    d = _art()["clade_demo"]
    assert d["verdict_under_standard_table"] == "refuted"
    assert d["verdict_under_clade_table"] == "satisfied"
    assert d["n_internal_stops_under_table1"] >= 1


def test_the_demo_gene_is_not_a_cherry_pick():
    """A single gene proves little; the census over the whole G37 CDS set is what makes it general."""
    c = _art()["clade_demo"]["census"]
    assert c["n_clean_cds"] > 100
    assert c["n_false_nonsense_under_table1"] / c["n_clean_cds"] > 0.5
    assert "efetch" in _art()["clade_demo"]["reproduce_census"]


def test_the_table4_translation_passes_a_biology_sanity_check():
    """Independent evidence the clade reading is correct, not just internally consistent: the protein
    opens with a canonical P-loop nucleotide-binding motif, as a uridine kinase should."""
    d = _art()["clade_demo"]
    assert d["protein_prefix_under_table4"].startswith("M")
    assert "SGGSCSGKTT" in d["protein_prefix_under_table4"]
    assert "P-loop" in d["biology_sanity_check"]


def test_every_inert_constraint_is_acknowledged_with_a_reason():
    """The whole point of the layer: a filter that removes nothing may not pose as a control."""
    kc = _art()["kill_counts"]["constraints"]
    for key, c in kc.items():
        if c["status"] == "inert":
            assert c["acknowledged_reason"], f"{key} is inert without a reason"


def test_inert_and_not_evaluated_are_reported_separately():
    kc = _art()["kill_counts"]
    statuses = {c["status"] for c in kc["constraints"].values()}
    assert "inert" in statuses
    assert kc["summary"]["n_inert"] >= 1
    assert "not evidence of protection" in kc["reading"].lower().replace("--", "")


def test_exactly_the_predicted_law_is_active():
    """The pre-stated prediction was that only `no_internal_stop_codon` would be active, and only via
    the Mycoplasma CDS scored under the wrong code. That is a falsifiable claim; this checks it held."""
    kc = _art()["kill_counts"]["constraints"]
    active = sorted(k for k, c in kc.items() if c["status"] == "active")
    assert active == ["no_internal_stop_codon"], active
    assert kc["no_internal_stop_codon"]["n_refuted"] == 1


def test_the_reachability_law_abstains_rather_than_claiming_a_clean_pass():
    """The probe supplies no substitution, so this law could not be evaluated at all -- it must read
    not_evaluated, not inert and certainly not satisfied."""
    kc = _art()["kill_counts"]["constraints"]["substitution_reachable_by_single_nt"]
    assert kc["n_evaluated"] == 0 and kc["status"] == "not_evaluated"
    assert kc["n_inapplicable"] > 0


def test_the_honest_scope_names_the_census_limit():
    scope = " ".join(_art()["honest_scope"]).lower()
    assert "not across all mycoplasma" in scope or "not across all" in scope
    assert "inert does not mean protective" in scope
