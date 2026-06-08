"""Biocide/disinfectant resistance (DisinFinder) decoder — offline units + real-BLAST e2e."""
import json
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import pytest  # noqa: E402

from dna_decode.disinfinder.cli import main  # noqa: E402
from dna_decode.disinfinder.runner import call_disinfectant_genes  # noqa: E402

_HAS_BLAST = bool(shutil.which("blastn", path="C:/Users/Farshad/ncbi-blast/bin") or shutil.which("blastn")) \
    and bool(shutil.which("makeblastdb", path="C:/Users/Farshad/ncbi-blast/bin") or shutil.which("makeblastdb"))


def _varied(n, seed):
    b = "ACGT"; x = seed; o = []
    for _ in range(n):
        x = (1103515245 * x + 12345) & 0x7FFFFFFF; o.append(b[(x >> 24) & 3])
    return "".join(o)


def test_offline_safe(monkeypatch, tmp_path):
    import dna_decode.typing.blast_caller as bc
    monkeypatch.setattr(bc, "find_blastn", lambda: None)
    db = tmp_path / "d.fsa"; db.write_text(">qacA_1_X\nACGT\n", encoding="utf-8")
    g = tmp_path / "g.fna"; g.write_text(">c\nACGT\n", encoding="utf-8")
    r = call_disinfectant_genes(g, db)
    assert r["status"] == "unavailable" and r["genes"] == []


def test_cli_missing_fasta_exit2():
    assert main(["/nonexistent.fna"]) == 2


@pytest.mark.skipif(not _HAS_BLAST, reason="BLAST+ not installed")
def test_e2e_detects_qacA(tmp_path):
    gene = _varied(600, 1)
    absent = _varied(600, 2)
    db = tmp_path / "disinfectants.fsa"
    db.write_text(f">qacA_1_AB566410\n{gene}\n>formA_1_X73835\n{absent}\n", encoding="utf-8")
    asm = tmp_path / "asm.fna"
    asm.write_text(f">c\n{_varied(300,7)}{gene}{_varied(300,8)}\n", encoding="utf-8")
    out = tmp_path / "d.json"
    rc = main([str(asm), "--db", str(db), "--out", str(out), "--json-only"])
    assert rc == 0
    rec = json.loads(out.read_text())
    assert "qacA" in rec["genes"] and "formA" not in rec["genes"]
    assert rec["schema"] == "disinfectant-gene-call-v0"


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-q"]))
