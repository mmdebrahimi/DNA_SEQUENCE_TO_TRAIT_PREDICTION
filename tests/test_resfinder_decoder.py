"""Acquired-AMR-gene (ResFinder) decoder — offline units + real-BLAST e2e on synthetic fixtures."""
import json
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import pytest  # noqa: E402

from dna_decode.resfinder.cli import main  # noqa: E402
from dna_decode.resfinder.runner import call_resistance_genes, gene_of  # noqa: E402

_HAS_BLAST = bool(shutil.which("blastn", path="C:/Users/Farshad/ncbi-blast/bin") or shutil.which("blastn")) \
    and bool(shutil.which("makeblastdb", path="C:/Users/Farshad/ncbi-blast/bin") or shutil.which("makeblastdb"))


def _varied(n, seed=7):
    b = "ACGT"; x = seed; o = []
    for _ in range(n):
        x = (1103515245 * x + 12345) & 0x7FFFFFFF
        o.append(b[(x >> 24) & 3])
    return "".join(o)


def test_gene_of_parses_real_resfinder_headers():
    assert gene_of("blaNDM-19_1_MF370080") == "blaNDM-19"
    assert gene_of("blaOXA-368_1_KT736121") == "blaOXA-368"
    assert gene_of("aac(6')-Ib_2_M23634") == "aac(6')-Ib"
    assert gene_of("aac(6')-30-aac(6')-Ib'_1_AJ584652") == "aac(6')-30-aac(6')-Ib'"


def test_offline_safe(monkeypatch, tmp_path):
    import dna_decode.typing.blast_caller as bc
    monkeypatch.setattr(bc, "find_blastn", lambda: None)
    db = tmp_path / "bl.fsa"; db.write_text(">blaNDM-1_1_X\nACGT\n", encoding="utf-8")
    g = tmp_path / "g.fna"; g.write_text(">c\nACGT\n", encoding="utf-8")
    r = call_resistance_genes(g, db, drug_class="beta-lactam")
    assert r["status"] == "unavailable" and r["genes"] == []


def test_cli_missing_fasta_exit2():
    assert main(["/nonexistent.fna"]) == 2


def test_cli_no_db_dir_unavailable(tmp_path):
    g = tmp_path / "g.fna"; g.write_text(">c\nACGT\n", encoding="utf-8")
    out = tmp_path / "o.json"
    rc = main([str(g), "--db-dir", str(tmp_path / "nope"), "--out", str(out), "--json-only"])
    assert rc == 3
    assert json.loads(out.read_text())["status"] == "unavailable"


@pytest.mark.skipif(not _HAS_BLAST, reason="BLAST+ not installed")
def test_end_to_end_real_blast_detects_gene(tmp_path):
    gene = _varied(600, 1)            # blaNDM-1
    absent = _varied(600, 2)          # blaKPC-2 (not in assembly)
    db = tmp_path / "beta-lactam.fsa"
    db.write_text(f">blaNDM-1_1_ACC\n{gene}\n>blaKPC-2_1_ACC\n{absent}\n", encoding="utf-8")
    flank = _varied(400, 9)
    asm = tmp_path / "asm.fna"
    asm.write_text(f">contig1\n{flank}{gene}{flank}\n", encoding="utf-8")
    res = call_resistance_genes(asm, db, drug_class="beta-lactam")
    assert res["status"] == "ok", res
    names = {g["gene"] for g in res["genes"]}
    assert "blaNDM-1" in names and "blaKPC-2" not in names
    assert res["genes"][0]["drug_class"] == "beta-lactam"
    # CLI dir-scan path (single class file) returns 0 + independent-baseline flag
    out = tmp_path / "o.json"
    assert main([str(asm), "--db-dir", str(tmp_path), "--out", str(out), "--json-only"]) == 0
    rec = json.loads(out.read_text())
    assert "blaNDM-1" in rec["genes"] and rec["caller"]["caller_is_independent_baseline"] is True


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-q"]))
