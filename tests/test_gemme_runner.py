"""Offline tests for the finalized GEMME runner (dna_decode/forward/gemme_scorer).

No Docker, no ColabFold — pins the two PURE helpers that turn real GEMME I/O into a variant table:
the a3m->aligned-FASTA conversion and the GEMME normPred-matrix parser. Format fixtures are taken
verbatim from the real TEM-1 run (Spearman 0.7191).
"""
from __future__ import annotations

from dna_decode.forward.gemme_scorer import a3m_to_aligned_fasta, _parse_gemme_matrix, run_gemme, GemmeUnavailable
import pytest


def test_a3m_to_aligned_fasta_drops_lowercase_inserts(tmp_path):
    # a3m: query all-uppercase (match cols); homolog has a lowercase INSERT 'q' that must be dropped,
    # so every emitted row is query length.
    a3m = tmp_path / "m.a3m"
    a3m.write_text(">query\nMKV-L\n>homolog\nMKqVAL\n", encoding="utf-8")  # homolog: M K (q insert) V A L
    out = tmp_path / "ali.fasta"
    qseq, n = a3m_to_aligned_fasta(a3m, out)
    assert n == 2
    assert qseq == "MKVL"                       # gaps stripped from the returned query sequence
    lines = out.read_text(encoding="utf-8").splitlines()
    seqs = [lines[i] for i in range(1, len(lines), 2)]
    assert all(len(s) == 5 for s in seqs)        # query match length = MKV-L = 5 cols
    assert seqs[1] == "MKVAL"                     # lowercase 'q' insert dropped -> query-length row


def test_parse_gemme_matrix_real_format(tmp_path):
    # verbatim GEMME normPred format: header "V1".."Vn", rows start with a LOWERCASE aa label.
    # 3-position query 'MK C'; GEMME native scale is negative=deleterious (already higher=preserved).
    pred = tmp_path / "q_normPred_evolCombi.txt"
    pred.write_text(
        '"V1" "V2" "V3"\n'
        '"a" -0.5 -1.2 NA\n'
        '"c" -2.0 -0.1 -0.9\n'
        '"m" NA -3.0 -1.1\n',
        encoding="utf-8")
    tbl = _parse_gemme_matrix(pred, query_seq="MKC")
    # WT at pos1=M,2=K,3=C. Cell (aa='A',pos1)=-0.5 -> 'M1A': -0.5 (M!=A). (aa='A',pos3)=NA skipped.
    assert tbl["M1A"] == -0.5
    assert tbl["K2A"] == -1.2
    assert tbl["M1C"] == -2.0
    assert tbl["C3M"] == -1.1
    # WT-to-WT is never emitted (pos3 aa='C' would be C3C -> skipped)
    assert "C3C" not in tbl
    # higher = preserved: the least-negative score is the most preserved (M1A -0.5 > K2A -1.2)
    assert tbl["M1A"] > tbl["K2A"]


def test_parse_gemme_matrix_raises_on_empty(tmp_path):
    empty = tmp_path / "e.txt"
    empty.write_text("\n", encoding="utf-8")
    with pytest.raises(ValueError):
        _parse_gemme_matrix(empty, query_seq="MK")


def test_run_gemme_raises_without_docker(monkeypatch, tmp_path):
    import dna_decode.forward.gemme_scorer as g
    monkeypatch.setattr(g.shutil, "which", lambda name: None)   # no docker on PATH
    (tmp_path / "m.a3m").write_text(">query\nMK\n", encoding="utf-8")
    with pytest.raises(GemmeUnavailable):
        run_gemme(tmp_path / "m.a3m", "MK", work_dir=tmp_path)
