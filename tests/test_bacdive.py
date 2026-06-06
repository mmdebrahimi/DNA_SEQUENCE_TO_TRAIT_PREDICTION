"""Tests for the BacDive carbon-utilization loader (dna_decode/data/bacdive).

Pins: tolerant column mapping, multi-vocabulary binarization, ambiguous-token
drop, organism filter, majority-vote conflict resolution, per-source census +
the >=min_strains accessor. Synthetic fixtures only — no network, no BacDive API.
Runnable via pytest OR standalone.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import pandas as pd  # noqa: E402

from dna_decode.data.bacdive import (  # noqa: E402
    binarize_utilization,
    census_carbon_sources,
    get_binary_labels,
    get_carbon_source_list,
    load_bacdive_carbon,
)


def _write(tmp_path, rows, header="strain_id,carbon_source,utilization,organism,assembly_accession"):
    p = tmp_path / "export.csv"
    p.write_text(header + "\n" + "\n".join(rows), encoding="utf-8")
    return p


def test_binarize_vocabularies():
    for tok in ("+", "positive", "YES", "growth", "1", "utilized"):
        assert binarize_utilization(tok) == 1
    for tok in ("-", "negative", "no", "no_growth", "0", "FALSE"):
        assert binarize_utilization(tok) == 0


def test_binarize_unknown_raises():
    try:
        binarize_utilization("maybe")
    except ValueError:
        return
    raise AssertionError("expected ValueError on unknown token")


def test_load_filters_organism(tmp_path):
    rows = [
        "562.1,glucose,+,Escherichia coli,GCF_1",
        "999.1,glucose,+,Salmonella enterica,GCF_2",
    ]
    t = load_bacdive_carbon(_write(tmp_path, rows))
    assert set(t["strain_id"]) == {"562.1"}


def test_load_drops_ambiguous(tmp_path):
    rows = [
        "562.1,glucose,+,Escherichia coli,GCF_1",
        "562.2,glucose,weak,Escherichia coli,GCF_2",   # ambiguous → dropped
        "562.3,glucose,?,Escherichia coli,GCF_3",       # ambiguous → dropped
        "562.4,glucose,-,Escherichia coli,GCF_4",
    ]
    t = load_bacdive_carbon(_write(tmp_path, rows))
    assert len(t) == 2
    assert set(t["strain_id"]) == {"562.1", "562.4"}


def test_tolerant_column_mapping(tmp_path):
    # bacdive_id / substrate / value / species / ncbi_accession aliases
    p = tmp_path / "alt.tsv"
    p.write_text(
        "bacdive_id\tsubstrate\tvalue\tspecies\tncbi_accession\n"
        "562.1\tD-glucose\tpositive\tEscherichia coli\tGCF_1\n"
        "562.2\tD-glucose\tnegative\tEscherichia coli\tGCF_2\n",
        encoding="utf-8",
    )
    t = load_bacdive_carbon(p)
    assert len(t) == 2
    assert set(t["carbon_source"]) == {"d-glucose"}
    assert t.iloc[0]["assembly_accession"] == "GCF_1"


def test_get_binary_labels_majority_vote(tmp_path):
    rows = [
        "562.1,glucose,+,Escherichia coli,GCF_1",
        "562.1,glucose,+,Escherichia coli,GCF_1",
        "562.1,glucose,-,Escherichia coli,GCF_1",   # 2 pos / 1 neg → 1
        "562.2,glucose,-,Escherichia coli,GCF_2",
    ]
    labels = get_binary_labels(load_bacdive_carbon(_write(tmp_path, rows)), "glucose")
    assert labels == {"562.1": 1, "562.2": 0}


def test_carbon_source_list_min_strains(tmp_path):
    rows = []
    for i in range(60):
        rows.append(f"562.{i},glucose,{'+' if i % 2 else '-'},Escherichia coli,GCF_{i}")
    for i in range(10):
        rows.append(f"562.{i},sorbitol,+,Escherichia coli,GCF_{i}")
    t = load_bacdive_carbon(_write(tmp_path, rows))
    assert get_carbon_source_list(t, min_strains=50) == ["glucose"]
    assert set(get_carbon_source_list(t, min_strains=5)) == {"glucose", "sorbitol"}


def test_census_counts_and_accession(tmp_path):
    rows = [
        "562.1,glucose,+,Escherichia coli,GCF_1",
        "562.2,glucose,-,Escherichia coli,",          # no accession
        "562.3,glucose,+,Escherichia coli,GCF_3",
    ]
    cen = census_carbon_sources(load_bacdive_carbon(_write(tmp_path, rows)))
    assert len(cen) == 1
    c = cen[0]
    assert c.carbon_source == "glucose"
    assert c.n_strains == 3
    assert c.n_positive == 2 and c.n_negative == 1
    assert c.n_with_accession == 2          # 562.1 + 562.3
    assert abs(c.minority_fraction - (1 / 3)) < 1e-9


if __name__ == "__main__":
    import traceback
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_") and callable(v)]
    import tempfile
    passed = 0
    for fn in fns:
        try:
            import inspect
            if "tmp_path" in inspect.signature(fn).parameters:
                with tempfile.TemporaryDirectory() as d:
                    fn(Path(d))
            else:
                fn()
            passed += 1
            print(f"PASS {fn.__name__}")
        except Exception:
            print(f"FAIL {fn.__name__}")
            traceback.print_exc()
    print(f"\n{passed}/{len(fns)} passed")
    sys.exit(0 if passed == len(fns) else 1)
