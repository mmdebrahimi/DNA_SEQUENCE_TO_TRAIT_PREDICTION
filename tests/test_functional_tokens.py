"""Drug-general functional-token featurizer — synthetic AMRFinder-cache fixtures (no Docker/genomes)."""
from __future__ import annotations

import pandas as pd

from dna_decode.eval.functional_tokens import build_feature_matrix, strain_functional_tokens

CIPRO = "ciprofloxacin"
TET = "tetracycline"


def _write_run(root, accession, *, mutations=(), main=()):
    d = root / accession
    d.mkdir(parents=True, exist_ok=True)
    if mutations:
        (d / "mutations.tsv").write_text("Element symbol\n" + "\n".join(mutations), encoding="utf-8")
    if main:
        (d / "main.tsv").write_text("Element symbol\n" + "\n".join(main), encoding="utf-8")


def _row(strain_id, accession):
    return pd.Series({"strain_id": strain_id, "assembly_accession": accession})


def test_cipro_qrdr_allele_and_mech_rollup(tmp_path):
    _write_run(tmp_path, "GCA_1", mutations=["gyrA_S83L"])
    toks = strain_functional_tokens(_row("s1", "GCA_1"), tmp_path, CIPRO)
    assert "gyrA_S83L" in toks
    assert "mech:QRDR_target_alteration" in toks


def test_cipro_plasmid_qnr_gene(tmp_path):
    _write_run(tmp_path, "GCA_2", main=["qnrS"])
    toks = strain_functional_tokens(_row("s2", "GCA_2"), tmp_path, CIPRO)
    assert "gene:qnrS" in toks
    assert "mech:plasmid_protect_modify" in toks


def test_tet_efflux_and_ribosomal_genes(tmp_path):
    """The drug-general alphabet picks up tet determinants in AMRFinder's PARENTHESISED form (the bug the
    pilot caught: AMRFinder reports `tet(A)`, the catalog has `tetA`)."""
    _write_run(tmp_path, "GCA_T", main=["tet(A)", "tet(M)"])
    toks = strain_functional_tokens(_row("t", "GCA_T"), tmp_path, TET)
    assert "gene:tetA" in toks and "mech:tet_efflux" in toks
    assert "gene:tetM" in toks and "mech:tet_ribosomal_protection" in toks
    assert not any("gyrA" in t for t in toks)  # a cipro QRDR allele is NOT a tet token


def test_synonymous_point_dropped(tmp_path):
    _write_run(tmp_path, "GCA_3", mutations=["gyrA_S83S"])
    toks = strain_functional_tokens(_row("s3", "GCA_3"), tmp_path, CIPRO)
    assert "gyrA_S83S" not in toks


def test_no_lineage_identity_token(tmp_path):
    _write_run(tmp_path, "GCA_4", mutations=["gyrA_S83L"])
    toks = strain_functional_tokens(_row("s4", "GCA_4"), tmp_path, CIPRO)
    assert not any(t.startswith("mlst:") for t in toks)


def test_none_when_cache_absent(tmp_path):
    assert strain_functional_tokens(_row("s5", "GCA_MISSING"), tmp_path, CIPRO) is None


def test_build_feature_matrix_determinism_and_dropped(tmp_path):
    _write_run(tmp_path, "GCA_A", mutations=["gyrA_S83L"])
    _write_run(tmp_path, "GCA_B", main=["qnrB"])
    df = pd.DataFrame([
        {"strain_id": "a", "assembly_accession": "GCA_A"},
        {"strain_id": "b", "assembly_accession": "GCA_B"},
        {"strain_id": "c", "assembly_accession": "GCA_MISSING"},  # dropped
    ])
    X, vocab, ids, dropped = build_feature_matrix(df, tmp_path, CIPRO)
    assert ids == ["a", "b"]
    assert dropped == ["c"]
    assert X.shape == (2, len(vocab))
    assert vocab == sorted(vocab)
    X2, vocab2, ids2, _ = build_feature_matrix(df, tmp_path, CIPRO)
    assert vocab2 == vocab and ids2 == ids and (X2 == X).all()
