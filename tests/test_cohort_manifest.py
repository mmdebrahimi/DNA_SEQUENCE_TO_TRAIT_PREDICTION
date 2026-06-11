"""Pin the accession-manifest registry (dna_decode/eval/cohort_manifest.py).

The leakage-safety property is load-bearing for the whole provenance-disjoint tool: exact-self exclusion
(not substring) + fail-closed on incomplete load. These tests pin both.
"""
from __future__ import annotations

from dna_decode.eval import cohort_manifest as cm


def _mk_cohort(d, name, accs):
    p = d / name
    p.mkdir(parents=True)
    (p / "selected.tsv").write_text("".join(f"{a}\tR\n" for a in accs), encoding="utf-8")


def test_build_manifest_scans_selected_tsv(tmp_path):
    raw = tmp_path / "raw"
    _mk_cohort(raw, "klebsiella_cipro", ["GCA_1.1", "GCA_2.1"])
    _mk_cohort(raw, "klebsiella_provdisjoint_ciprofloxacin", ["GCA_9.1"])
    m = cm.build_manifest(data_raw=str(raw), data_processed=str(tmp_path / "none"))
    assert not m.incomplete
    names = set(m.cohort_names())
    assert names == {"klebsiella_cipro", "klebsiella_provdisjoint_ciprofloxacin"}
    cal = next(c for c in m.cohorts if c.name == "klebsiella_cipro")
    assert cal.accessions == {"GCA_1.1", "GCA_2.1"} and cal.source == "selected_tsv"


def test_prior_accessions_exact_self_excludes_only_named_cohort(tmp_path):
    """C1: exclude ONLY the exact current cohort — every OTHER cohort (incl. a prior provdisjoint run for
    the same organism) is excluded from the available pool, forcing fresh accessions."""
    raw = tmp_path / "raw"
    _mk_cohort(raw, "klebsiella_cipro", ["GCA_cal.1"])
    _mk_cohort(raw, "klebsiella_provdisjoint_ciprofloxacin_v1", ["GCA_prior.1"])
    _mk_cohort(raw, "klebsiella_provdisjoint_ciprofloxacin_v2", ["GCA_self.1"])
    m = cm.build_manifest(data_raw=str(raw), data_processed=str(tmp_path / "none"))
    excl = cm.prior_accessions(m, exclude_cohort="klebsiella_provdisjoint_ciprofloxacin_v2")
    # the prior provdisjoint cohort IS excluded (in the set); only the exact self is not
    assert "GCA_prior.1" in excl and "GCA_cal.1" in excl
    assert "GCA_self.1" not in excl


def test_substring_collision_does_not_over_exclude(tmp_path):
    """Exact identity must not match a different cohort that merely shares a name prefix."""
    raw = tmp_path / "raw"
    _mk_cohort(raw, "kleb_provdisjoint_cipro", ["GCA_a.1"])
    _mk_cohort(raw, "kleb_provdisjoint_cipro_extra", ["GCA_b.1"])
    m = cm.build_manifest(data_raw=str(raw), data_processed=str(tmp_path / "none"))
    excl = cm.prior_accessions(m, exclude_cohort="kleb_provdisjoint_cipro")
    assert "GCA_b.1" in excl and "GCA_a.1" not in excl  # only the exact name excluded


def test_malformed_parquet_marks_incomplete(tmp_path):
    """C1 fail-closed signal: a parquet that won't load sets incomplete=True (caller must refuse)."""
    raw = tmp_path / "raw"; _mk_cohort(raw, "x_cipro", ["GCA_1.1"])
    proc = tmp_path / "proc"; proc.mkdir()
    (proc / "broken_cohort.parquet").write_bytes(b"NOT A PARQUET FILE")
    m = cm.build_manifest(data_raw=str(raw), data_processed=str(proc))
    assert m.incomplete is True
    assert any("broken_cohort" in w for w in m.warnings)


def test_empty_dirs_complete(tmp_path):
    m = cm.build_manifest(data_raw=str(tmp_path / "a"), data_processed=str(tmp_path / "b"))
    assert m.incomplete is False and m.cohorts == []
