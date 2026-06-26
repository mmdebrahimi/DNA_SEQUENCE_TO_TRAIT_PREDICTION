"""Functional-token featurizer — synthetic AMRFinder-cache fixtures (no Docker, no real genomes)."""
from __future__ import annotations

import pandas as pd

from dna_decode.eval.functional_tokens import build_feature_matrix, strain_functional_tokens

DRUG = "ciprofloxacin"


def _write_run(root, accession, *, mutations=(), main=()):
    d = root / accession
    d.mkdir(parents=True, exist_ok=True)
    if mutations:
        (d / "mutations.tsv").write_text("Element symbol\n" + "\n".join(mutations), encoding="utf-8")
    if main:
        (d / "main.tsv").write_text("Element symbol\n" + "\n".join(main), encoding="utf-8")


def _row(strain_id, accession):
    return pd.Series({"strain_id": strain_id, "assembly_accession": accession})


def test_qrdr_allele_and_mech_rollup(tmp_path):
    _write_run(tmp_path, "GCA_1", mutations=["gyrA_S83L"])
    toks = strain_functional_tokens(_row("s1", "GCA_1"), tmp_path, DRUG)
    assert "gyrA_S83L" in toks
    assert "mech:qrdr" in toks
    assert "mech:plasmid" not in toks


def test_plasmid_token_and_mech_rollup(tmp_path):
    _write_run(tmp_path, "GCA_2", main=["qnrS"])
    toks = strain_functional_tokens(_row("s2", "GCA_2"), tmp_path, DRUG)
    assert any(t.startswith("plasmid:") for t in toks)
    assert "mech:plasmid" in toks


def test_synonymous_dropped(tmp_path):
    _write_run(tmp_path, "GCA_3", mutations=["gyrA_S83S"])  # synonymous -> not a token
    toks = strain_functional_tokens(_row("s3", "GCA_3"), tmp_path, DRUG)
    assert "gyrA_S83S" not in toks
    assert "mech:qrdr" not in toks  # no non-synonymous QRDR allele


def test_no_lineage_identity_token(tmp_path):
    """The alphabet is mechanism-only -- no mlst:<ST> identity token (the /probe Open-Q2 honesty point)."""
    _write_run(tmp_path, "GCA_4", mutations=["gyrA_S83L"])
    toks = strain_functional_tokens(_row("s4", "GCA_4"), tmp_path, DRUG)
    assert not any(t.startswith("mlst:") for t in toks)


def test_none_when_cache_absent(tmp_path):
    assert strain_functional_tokens(_row("s5", "GCA_MISSING"), tmp_path, DRUG) is None


def test_build_feature_matrix_determinism_and_dropped(tmp_path):
    _write_run(tmp_path, "GCA_A", mutations=["gyrA_S83L"])
    _write_run(tmp_path, "GCA_B", main=["qnrB"])
    df = pd.DataFrame([
        {"strain_id": "a", "assembly_accession": "GCA_A"},
        {"strain_id": "b", "assembly_accession": "GCA_B"},
        {"strain_id": "c", "assembly_accession": "GCA_MISSING"},  # dropped
    ])
    X, vocab, ids, dropped = build_feature_matrix(df, tmp_path, DRUG)
    assert ids == ["a", "b"]
    assert dropped == ["c"]
    assert X.shape == (2, len(vocab))
    assert vocab == sorted(vocab)  # deterministic column order
    # rebuild -> identical
    X2, vocab2, ids2, _ = build_feature_matrix(df, tmp_path, DRUG)
    assert vocab2 == vocab and ids2 == ids and (X2 == X).all()
