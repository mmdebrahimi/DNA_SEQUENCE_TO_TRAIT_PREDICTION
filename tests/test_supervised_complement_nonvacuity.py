"""Guards for the per-class non-vacuity scan of the HIV supervised blind-spot complement (2026-09-30).

The real scan needs gitignored Stanford data; these pin the PURE parts plus the guard whose absence
produced two confident wrong design conclusions in a row.
"""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
_spec = importlib.util.spec_from_file_location(
    "supervised_complement_nonvacuity", REPO / "scripts" / "supervised_complement_nonvacuity.py")
S = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(S)


def test_the_partition_guard_refuses_every_one_bucket_shape():
    """THE GUARD BOTH EARLIER SCANS LACKED. One reported '0 catalog-susceptible isolates' after swallowing
    an error on every row; the other put EVERY isolate in the blind spot via a token-format mismatch and
    called the full cohort's AUROC a blind-spot result. Both failures have this ONE signature."""
    assert S.partition_is_non_vacuous(500, 500)
    assert S.partition_is_non_vacuous(2339, 1877)          # the real NNRTI partition
    for bad in ((0, 900), (900, 0), (0, 0), (999, 1), (1, 999)):
        assert not S.partition_is_non_vacuous(*bad), f"{bad} must refuse"


def test_enrichment_refuses_a_zero_base_rate_instead_of_dividing():
    """A blind spot with no resistant isolate has an UNDEFINED enrichment, not an infinite one. Refusing is
    what keeps 'the threshold found nothing' distinct from 'there was nothing to find'."""
    with pytest.raises(S.ScanRefused):
        S.enrichment(0.5, 0.0)
    assert S.enrichment(0.729, 0.047) == pytest.approx(15.51, abs=0.02)


def test_self_check_passes():
    S.self_check()


def test_the_two_dataset_variants_are_distinct_and_named():
    """`.Full.txt` is the model's own TRAINING set and `.txt` is what the rest of the repo validates on.
    Conflating them would report a maximally in-sample number as if it were the repo's usual one."""
    assert S.DATASETS["training"] == ".Full.txt"
    assert S.DATASETS["validation"] == ".txt"
    assert set(S.CASES) == {"NNRTI", "PI", "INSTI"}


def test_every_class_scanned_is_one_the_complement_actually_supports():
    """A scan of a class the complement does not cover would report a shipped model's absence as a result."""
    from dna_decode.data import hiv_supervised_complement as C
    for cls in S.CASES:
        assert cls in C.SUPPORTED_CLASSES


def test_the_scanned_drug_routes_to_the_gene_the_scan_assumes():
    """The blind-spot partition calls the DEPLOYED catalog, which is keyed by gene; a wrong gene key would
    silently score an empty observed set and call every isolate susceptible."""
    from dna_decode.data.hiv_amr import gene_for_hiv_drug
    assert gene_for_hiv_drug("efavirenz") == "RT"
    assert gene_for_hiv_drug("lopinavir") == "PR"
    assert gene_for_hiv_drug("raltegravir") == "IN"


def test_a_refused_cell_carries_no_metrics():
    """`refused` entries must never be readable as a measured zero."""
    art = REPO / "wiki" / "supervised_complement_nonvacuity_2026-09-30.json"
    if not art.exists():
        pytest.skip("artifact not present (gitignored data absent at build time)")
    doc = json.loads(art.read_text(encoding="utf-8"))
    for r in doc.get("refused", []):
        assert "refused" in r and not any(k in r for k in ("precision", "recall", "enrichment"))


def test_the_committed_artifact_keeps_PI_metric_free_and_labels_the_in_sample_auroc():
    """PI's arm NEVER FIRES on either dataset (2-3 truly-resistant isolates in a 614-903 blind spot). It
    must therefore carry a verdict and NO precision/recall/enrichment -- a number there would assert
    call-time value the evidence does not support. And the in-sample AUROC must stay named as such: it is
    0.94-0.99 while the deployable claim is the leave-study-out 0.81/0.89."""
    art = REPO / "wiki" / "supervised_complement_nonvacuity_2026-09-30.json"
    if not art.exists():
        pytest.skip("artifact not present (gitignored data absent at build time)")
    doc = json.loads(art.read_text(encoding="utf-8"))
    pi = [c for c in doc["cells"] if c["class"] == "PI"]
    assert pi, "PI must appear -- a silent arm still has to be reported"
    for c in pi:
        assert c["verdict"] == "THRESHOLD_NEVER_FIRES" and c["n_flagged"] == 0
        assert not any(k in c for k in ("precision", "recall", "enrichment"))
        assert c["risk_max"] < c["threshold"]
    for c in doc["cells"]:
        assert "blindspot_auroc" not in c, "an unqualified AUROC key would read as performance"
        if any("auroc" in k for k in c):
            assert any("IN_SAMPLE" in k for k in c if "auroc" in k)


def test_INSTI_is_non_vacuous_in_the_committed_artifact():
    """The question this scan existed to answer: does the complement's shipped threshold separate anything
    on a class other than NNRTI? For INSTI, yes -- and it must be recorded as beating its base rate."""
    art = REPO / "wiki" / "supervised_complement_nonvacuity_2026-09-30.json"
    if not art.exists():
        pytest.skip("artifact not present (gitignored data absent at build time)")
    doc = json.loads(art.read_text(encoding="utf-8"))
    insti = [c for c in doc["cells"] if c["class"] == "INSTI"]
    assert insti
    for c in insti:
        assert c["verdict"] == "NON_VACUOUS"
        assert c["precision"] > c["base_rate"] and c["enrichment"] > 1.0
