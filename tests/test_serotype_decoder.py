"""E. coli serotype decoder — offline units + real-BLAST e2e on synthetic fixtures. No network/real DB."""
import json
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import pytest  # noqa: E402

from dna_decode.serotype.cli import main  # noqa: E402
from dna_decode.serotype.runner import antigen_of, call_serotype, gene_of  # noqa: E402

_HAS_BLAST = bool(shutil.which("blastn", path="C:/Users/Farshad/ncbi-blast/bin") or shutil.which("blastn")) \
    and bool(shutil.which("makeblastdb", path="C:/Users/Farshad/ncbi-blast/bin") or shutil.which("makeblastdb"))


def _varied(n, seed=7):
    b = "ACGT"; x = seed; o = []
    for _ in range(n):
        x = (1103515245 * x + 12345) & 0x7FFFFFFF
        o.append(b[(x >> 24) & 3])
    return "".join(o)


def test_antigen_and_gene_parse_real_headers():
    assert antigen_of("wzx_1_GU299791_O1") == "O1"
    assert antigen_of("fliC_307_AY249994_H9") == "H9"
    assert antigen_of("wzm_19_CP011331_O8") == "O8"
    assert antigen_of("gnd_5_X_O157") == "O157"
    assert antigen_of("foo_1_X_notantigen") is None
    assert gene_of("wzx_1_GU299791_O1") == "wzx"


def test_offline_safe(monkeypatch, tmp_path):
    import dna_decode.typing.blast_caller as bc
    monkeypatch.setattr(bc, "find_blastn", lambda: None)
    db = tmp_path / "s.fsa"; db.write_text(">wzx_1_X_O1\nACGT\n", encoding="utf-8")
    g = tmp_path / "g.fna"; g.write_text(">c\nACGT\n", encoding="utf-8")
    r = call_serotype(g, db)
    assert r["status"] == "unavailable" and r["serotype"] is None


def test_cli_missing_fasta_exit2():
    assert main(["/nonexistent.fna"]) == 2


@pytest.mark.skipif(not _HAS_BLAST, reason="BLAST+ not installed")
def test_end_to_end_real_blast_calls_O_and_H(tmp_path):
    o_allele = _varied(600, 1)      # wzx O25
    h_allele = _varied(600, 2)      # fliC H4
    absent = _varied(600, 3)        # wzx O1 (not in assembly)
    db = tmp_path / "sero.fsa"
    db.write_text(f">wzx_1_ACC_O25\n{o_allele}\n>fliC_1_ACC_H4\n{h_allele}\n>wzx_2_ACC_O1\n{absent}\n",
                  encoding="utf-8")
    flank = _varied(400, 42)
    asm = tmp_path / "asm.fna"
    asm.write_text(f">contig1\n{flank}{o_allele}{flank}{h_allele}{flank}\n", encoding="utf-8")
    res = call_serotype(asm, db)
    assert res["status"] == "ok", res
    assert res["o_antigen"] == "O25" and res["h_antigen"] == "H4", res
    assert res["serotype"] == "O25:H4"
    # CLI returns 0
    assert main([str(asm), "--db", str(db), "--json-only"]) == 0


@pytest.mark.skipif(not _HAS_BLAST, reason="BLAST+ not installed")
def test_partial_serotype_O_only(tmp_path):
    o_allele = _varied(600, 5)
    db = tmp_path / "sero.fsa"
    db.write_text(f">wzy_1_ACC_O104\n{o_allele}\n>fliC_1_ACC_H7\n{_varied(600, 9)}\n", encoding="utf-8")
    asm = tmp_path / "asm.fna"
    asm.write_text(f">c\n{_varied(300, 1)}{o_allele}{_varied(300, 2)}\n", encoding="utf-8")  # only O present
    res = call_serotype(asm, db)
    assert res["o_antigen"] == "O104" and res["h_antigen"] is None
    assert res["serotype"] == "O104:H?"


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-q"]))
