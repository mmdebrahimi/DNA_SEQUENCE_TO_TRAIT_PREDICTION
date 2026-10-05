"""Tests for the provenance-disjoint species-composition audit runner (Step 3).

Offline. The fixture cohorts are hand-built in tmp_path; the real-data tests skip when `data/` is absent
(it is a junction to D: and may be unmounted).
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from scripts import provdisjoint_species_audit as aud

REPO = Path(__file__).resolve().parent.parent
MEROP = REPO / "data" / "raw" / "klebsiella_provdisjoint_meropenem"

MAIN_HEADER = ("Protein id\tContig id\tStart\tStop\tElement symbol\tElement name\tClass\tSubclass\t"
               "% Identity to reference\tType\n")
KPC_ROW = ("NA\tc1\t1\t900\tblaKPC-2\tcarbapenemase\tBETA-LACTAM\tCARBAPENEM\t100.00\tAMR\n")
AMPC_ROW = ("NA\tc1\t1\t900\tampC-Kaer\tcephalosporinase\tBETA-LACTAM\tCEPHALOSPORIN\t99.00\tAMR\n")


def _cohort(tmp: Path, name: str, rows: dict[str, tuple[str, str, str]],
            *, runs_for: list[str] | None = None) -> Path:
    """rows = {accession: (R|S, organism_name, main_tsv_body)}. `runs_for` limits cached runs."""
    d = tmp / "data" / "raw" / name
    (d / "refseq").mkdir(parents=True)
    (d / "amrfinder_runs").mkdir(parents=True)
    lines = []
    for acc, (rs, org, body) in rows.items():
        lines.append(f"{acc}\t{rs}")
        g = d / "refseq" / acc
        g.mkdir(parents=True)
        (g / "annotations.gbk").write_text(f"LOCUS x\n  ORGANISM  {org}\n", encoding="utf-8")
        if runs_for is None or acc in runs_for:
            r = d / "amrfinder_runs" / acc
            r.mkdir(parents=True)
            (r / "main.tsv").write_text(MAIN_HEADER + body, encoding="utf-8")
    (d / "selected.tsv").write_text("\n".join(lines) + "\n", encoding="utf-8")
    return d


def _artifact(tmp: Path, slug: str, drug_abbr: str, metrics: dict, *,
              amrfinder_organism: str = "Klebsiella_pneumoniae",
              organism: str = "Klebsiella", registry_organism: str = "Klebsiella") -> Path:
    w = tmp / "wiki"
    w.mkdir(parents=True, exist_ok=True)
    p = w / f"provenance_disjoint_validation_{slug}_{drug_abbr}_2026-01-01.json"
    p.write_text(json.dumps({
        "_schema": "provenance-disjoint-validation-v1", "organism": organism,
        "registry_organism": registry_organism, "amrfinder_organism": amrfinder_organism,
        "drug": "meropenem", "metrics": metrics}), encoding="utf-8")
    return w


# --------------------------------------------------------------------------- discovery

def test_discovery_requires_both_labels_and_assemblies(tmp_path):
    raw = tmp_path / "data" / "raw"
    (raw / "x_provdisjoint_meropenem").mkdir(parents=True)          # neither -> excluded
    good = _cohort(tmp_path, "klebsiella_provdisjoint_meropenem",
                   {"A": ("R", "Klebsiella pneumoniae", KPC_ROW)})
    found = [p.name for p in aud.provdisjoint_cohorts(raw)]
    assert found == [good.name]


# --------------------------------------------------------------------------- reconcile gate

def test_reconcile_pass_yields_a_crosstab(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    d = _cohort(tmp_path, "klebsiella_provdisjoint_meropenem", {
        "A": ("R", "Klebsiella pneumoniae", KPC_ROW),     # pred R, truth R -> tp
        "B": ("R", "Klebsiella aerogenes", AMPC_ROW),     # pred S, truth R -> fn
        "C": ("S", "Klebsiella pneumoniae", AMPC_ROW),    # pred S, truth S -> tn
    })
    w = _artifact(tmp_path, "klebsiella", "merop",
                  {"n_scored": 3, "tp": 1, "fp": 0, "tn": 1, "fn": 1, "sens": 0.5, "spec": 1.0})
    row = aud.audit_cohort(d, wiki=w)
    assert row["reconciled"] is True
    assert row["outcome_crosstab"]["totals"] == {"tp": 1, "fn": 1, "tn": 1}
    assert row["false_negative_anatomy"]["n_false_negatives"] == 1
    assert row["false_negative_anatomy"]["n_fn_off_species"] == 1


def test_reconcile_FAILURE_withholds_the_crosstab_and_says_why(tmp_path, monkeypatch):
    """NON-VACUITY: the committed artifact disagrees by one cell, so the gate must bite."""
    monkeypatch.chdir(tmp_path)
    d = _cohort(tmp_path, "klebsiella_provdisjoint_meropenem", {
        "A": ("R", "Klebsiella pneumoniae", KPC_ROW),
    })
    w = _artifact(tmp_path, "klebsiella", "merop",
                  {"n_scored": 1, "tp": 99, "fp": 0, "tn": 0, "fn": 0})   # deliberately wrong
    row = aud.audit_cohort(d, wiki=w)
    assert row["reconciled"] is False
    assert row["outcome_crosstab"] is None
    assert "RECONCILE FAILED" in row["crosstab_withheld_reason"]
    # composition is unaffected -- it does not depend on predictions
    assert row["composition"]["verdict"] in ("MATCHES_SCORED_ORGANISM", "MIXED_SPECIES")


def test_partial_cached_runs_report_composition_and_NO_crosstab(tmp_path, monkeypatch):
    """Attributing outcomes on a partial run set is the silent-denominator failure."""
    monkeypatch.chdir(tmp_path)
    d = _cohort(tmp_path, "klebsiella_provdisjoint_meropenem", {
        "A": ("R", "Klebsiella pneumoniae", KPC_ROW),
        "B": ("S", "Klebsiella aerogenes", AMPC_ROW),
    }, runs_for=["A"])                                  # only 1 of 2 cached
    w = _artifact(tmp_path, "klebsiella", "merop", {"n_scored": 2, "tp": 1, "fp": 0, "tn": 1, "fn": 0})
    row = aud.audit_cohort(d, wiki=w)
    assert row["cached_amrfinder_runs"] == {"have": 1, "total": 2}
    assert row["reconciled"] is False
    assert row["outcome_crosstab"] is None
    assert "incomplete" in row["crosstab_withheld_reason"]
    assert row["composition"]["n_total"] == 2          # composition still measured


def test_a_missing_artifact_raises_rather_than_inventing_one(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    d = _cohort(tmp_path, "klebsiella_provdisjoint_meropenem",
                {"A": ("R", "Klebsiella pneumoniae", KPC_ROW)})
    (tmp_path / "wiki").mkdir(exist_ok=True)
    with pytest.raises(aud.AuditError, match="no committed provdisjoint artifact"):
        aud.audit_cohort(d, wiki=tmp_path / "wiki")


def test_the_scoring_identity_comes_from_the_ARTIFACT_not_the_directory_name(tmp_path, monkeypatch):
    """The directory name being untrustworthy at species level is the finding under test, so the
    audit must never read the organism out of it."""
    monkeypatch.chdir(tmp_path)
    d = _cohort(tmp_path, "klebsiella_provdisjoint_meropenem",
                {"A": ("R", "Klebsiella pneumoniae", KPC_ROW)})
    w = _artifact(tmp_path, "klebsiella", "merop",
                  {"n_scored": 1, "tp": 1, "fp": 0, "tn": 0, "fn": 0},
                  amrfinder_organism="Klebsiella_oxytoca")      # artifact disagrees with dir name
    row = aud.audit_cohort(d, wiki=w)
    assert row["amrfinder_organism"] == "Klebsiella_oxytoca"
    # ...and K. pneumoniae is therefore off-species for THIS artifact
    assert row["composition"]["n_same_genus_other_species"] == 1


# --------------------------------------------------------------------------- output shape

def test_main_writes_both_artifacts_and_exits_0(tmp_path, monkeypatch, capsys):
    monkeypatch.chdir(tmp_path)
    _cohort(tmp_path, "klebsiella_provdisjoint_meropenem",
            {"A": ("R", "Klebsiella pneumoniae", KPC_ROW)})
    w = _artifact(tmp_path, "klebsiella", "merop", {"n_scored": 1, "tp": 1, "fp": 0, "tn": 0, "fn": 0})
    rc = aud.main(["--raw", str(tmp_path / "data" / "raw"), "--wiki", str(w), "--out-date", "T"])
    assert rc == 0
    j = json.loads((w / "provdisjoint_species_audit_T.json").read_text())
    assert j["schema"] == aud.SCHEMA
    assert isinstance(j["cells"], list) and j["cells"][0]["organism"] == "Klebsiella"
    assert (w / "provdisjoint_species_audit_T.md").exists()


def test_cells_is_a_LIST_with_organism_and_drug_matching_the_sibling_sidecar(tmp_path, monkeypatch):
    """`provdisjoint_source_concentration.json` keys cells as a LIST of {organism, drug, ...} and its
    loader re-derives canonical_cell_key. A dict keyed by a stringified tuple would be a second
    convention for the same join."""
    monkeypatch.chdir(tmp_path)
    _cohort(tmp_path, "klebsiella_provdisjoint_meropenem",
            {"A": ("R", "Klebsiella pneumoniae", KPC_ROW)})
    w = _artifact(tmp_path, "klebsiella", "merop", {"n_scored": 1, "tp": 1, "fp": 0, "tn": 0, "fn": 0})
    aud.main(["--raw", str(tmp_path / "data" / "raw"), "--wiki", str(w), "--out-date", "T"])
    cells = json.loads((w / "provdisjoint_species_audit_T.json").read_text())["cells"]
    assert isinstance(cells, list)
    for c in cells:
        assert "organism" in c and "drug" in c
        from dna_decode.data.cell_key import canonical_cell_key
        assert canonical_cell_key(c["organism"], c["drug"])   # joinable


def test_an_unknown_cohort_name_exits_2(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    _cohort(tmp_path, "klebsiella_provdisjoint_meropenem",
            {"A": ("R", "Klebsiella pneumoniae", KPC_ROW)})
    w = _artifact(tmp_path, "klebsiella", "merop", {"n_scored": 1, "tp": 1, "fp": 0, "tn": 0, "fn": 0})
    assert aud.main(["--raw", str(tmp_path / "data" / "raw"), "--wiki", str(w),
                     "--cohort", "nope"]) == 2


# --------------------------------------------------------------------------- real data

def _committed() -> dict | None:
    got = sorted((REPO / "wiki").glob("provdisjoint_species_audit_*.json"))
    return json.loads(got[-1].read_text(encoding="utf-8")) if got else None


def test_the_committed_artifact_is_parseable_and_well_shaped():
    """A hand-unreadable or malformed artifact is cited as evidence and nobody notices -- this repo has
    shipped an invalid-JSON artifact before."""
    d = _committed()
    if d is None:
        pytest.skip("no committed species-audit artifact yet")
    assert d["schema"] == aud.SCHEMA
    assert isinstance(d["cells"], list) and d["cells"]
    for c in d["cells"]:
        for k in ("cohort_dir", "organism", "drug", "amrfinder_organism", "composition",
                  "cached_amrfinder_runs", "reconciled"):
            assert k in c, f"{c.get('cohort_dir')}: missing {k}"
        comp = c["composition"]
        assert comp["verdict"] in ("MATCHES_SCORED_ORGANISM", "MIXED_SPECIES",
                                   "INSUFFICIENT_RESOLUTION")
        if comp["verdict"] != "INSUFFICIENT_RESOLUTION":
            # the four buckets must PARTITION -- a denominator that does not sum is the failure mode
            assert (comp["n_matches_expected"] + comp["n_same_genus_other_species"]
                    + comp["n_other_genus"] + comp["n_unresolved"]) == comp["n_total"]
        if not c["reconciled"]:
            assert c["outcome_crosstab"] is None
            assert c.get("crosstab_withheld_reason"), "a withheld cross-tab must say why"


def test_the_committed_artifact_records_the_meropenem_finding():
    """Guards the headline the memo quotes: 16/16 FN off-species, and zero FN on the scored species."""
    d = _committed()
    if d is None:
        pytest.skip("no committed species-audit artifact yet")
    m = [c for c in d["cells"] if c["drug"] == "meropenem" and c["organism"] == "Klebsiella"]
    if not m:
        pytest.skip("meropenem cohort not in the committed artifact")
    c = m[0]
    assert c["composition"]["verdict"] == "MIXED_SPECIES"
    assert c["false_negative_anatomy"] == {
        "n_false_negatives": 16, "n_fn_off_species": 16, "n_fn_expected_species": 0,
        "fn_species_counts": {"Klebsiella aerogenes": 16}}
    assert c["outcome_crosstab"]["by_species"]["Klebsiella pneumoniae"]["fn"] == 0


def test_every_cohort_resolved_its_species_so_the_bar_never_bound():
    """The memo states all 580 assemblies resolved an ORGANISM line. If that stops holding, the
    MAX_UNRESOLVED_FRACTION bar starts mattering and the memo's claim needs re-deriving."""
    d = _committed()
    if d is None:
        pytest.skip("no committed species-audit artifact yet")
    total = sum(c["composition"]["n_total"] for c in d["cells"])
    unresolved = sum(c["composition"].get("n_unresolved", 0) for c in d["cells"])
    assert unresolved == 0, f"{unresolved} of {total} assemblies no longer resolve a species"


def test_the_real_meropenem_cohort_reconciles_to_14_3_27_16():
    """THE anchor. If this stops reproducing, the audit's attribution is void and must be re-derived
    rather than patched -- the committed artifact is a frozen unit."""
    if not MEROP.is_dir():
        pytest.skip("cached meropenem cohort absent (data/ junction may be unmounted)")
    row = aud.audit_cohort(MEROP)
    assert row["reconciled"] is True
    assert row["reconciled_confusion"]["tp"] == 14
    assert row["reconciled_confusion"]["fp"] == 3
    assert row["reconciled_confusion"]["tn"] == 27
    assert row["reconciled_confusion"]["fn"] == 16


def test_the_real_meropenem_cohort_is_MIXED_and_every_FN_is_off_species():
    """The measured finding: all 16 false negatives are K. aerogenes, zero are K. pneumoniae."""
    if not MEROP.is_dir():
        pytest.skip("cached meropenem cohort absent")
    row = aud.audit_cohort(MEROP)
    assert row["composition"]["verdict"] == "MIXED_SPECIES"
    fna = row["false_negative_anatomy"]
    assert fna["n_false_negatives"] == 16
    assert fna["n_fn_off_species"] == 16
    assert fna["n_fn_expected_species"] == 0
    assert set(fna["fn_species_counts"]) == {"Klebsiella aerogenes"}
    # and on the species the rule was validated on, there are NO misses at all
    by_sp = row["outcome_crosstab"]["by_species"]
    assert by_sp["Klebsiella pneumoniae"]["fn"] == 0
