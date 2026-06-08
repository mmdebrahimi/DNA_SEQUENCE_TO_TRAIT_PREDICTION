"""Plasmid replicon decoder — offline units + real-BLAST end-to-end on synthetic fixtures.

Pure parsing/aggregation + the offline-safe degrade path need no BLAST. The end-to-end test builds a tiny
synthetic replicon DB + an assembly embedding one replicon allele, runs real blastn (skipif absent), and
asserts the replicon is called. No network, no real PlasmidFinder DB download.
"""
import json
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import pytest  # noqa: E402

from dna_decode.plasmid.cli import main  # noqa: E402
from dna_decode.plasmid.runner import call_replicons, replicon_family  # noqa: E402


def _varied_seq(n: int, seed: int = 7) -> str:
    """Deterministic non-repetitive ACGT (low-complexity is dust-masked by blastn)."""
    b = "ACGT"; x = seed; out = []
    for _ in range(n):
        x = (1103515245 * x + 12345) & 0x7FFFFFFF
        out.append(b[(x >> 24) & 3])   # HIGH bits — LCG low bits have period <=4 (low-complexity, dust-masked)
    return "".join(out)

_HAS_BLAST = bool(shutil.which("blastn", path="C:/Users/Farshad/ncbi-blast/bin") or shutil.which("blastn")) \
    and bool(shutil.which("makeblastdb", path="C:/Users/Farshad/ncbi-blast/bin") or shutil.which("makeblastdb"))


def test_replicon_family_parses_real_plasmidfinder_headers():
    # exact header shapes from the real enterobacteriales.fsa
    assert replicon_family("IncFIA_1__AP001918") == "IncFIA"
    assert replicon_family("IncHI1B(R27)_1_R27_AF250878") == "IncHI1B(R27)"
    assert replicon_family("IncHI2_1__BX664015") == "IncHI2"
    assert replicon_family("IncI1-I(Alpha)_1__AP005147") == "IncI1-I(Alpha)"
    assert replicon_family("IncB/O/K/Z_2__GU256641") == "IncB/O/K/Z"
    assert replicon_family("pKPC-CAV1321_1__CP011611") == "pKPC-CAV1321"


def test_call_replicons_offline_safe_no_blast(tmp_path, monkeypatch):
    # force the no-blastn path -> unavailable, never raises
    import dna_decode.plasmid.runner as r
    monkeypatch.setattr(r, "find_blastn", lambda: None)
    db = tmp_path / "db.fsa"; db.write_text(">IncX3_1__JN247852\nACGT\n", encoding="utf-8")
    g = tmp_path / "g.fna"; g.write_text(">c\nACGT\n", encoding="utf-8")
    res = call_replicons(g, db)
    assert res["status"] == "unavailable" and res["replicons"] == [] and "blastn" in res["reason"]


def test_cli_offline_safe_exit3(tmp_path, monkeypatch, capsys):
    import dna_decode.plasmid.runner as r
    monkeypatch.setattr(r, "find_blastn", lambda: None)
    db = tmp_path / "db.fsa"; db.write_text(">IncX3_1__JN247852\nACGT\n", encoding="utf-8")
    g = tmp_path / "g.fna"; g.write_text(">c\nACGT\n", encoding="utf-8")
    out = tmp_path / "o.json"
    rc = main([str(g), "--db", str(db), "--out", str(out), "--json-only"])
    assert rc == 3
    rec = json.loads(out.read_text())
    assert rec["status"] == "unavailable" and rec["schema"] == "plasmid-replicon-call-v0"


def test_cli_missing_fasta_exit2(capsys):
    assert main(["/nonexistent/x.fna"]) == 2


@pytest.mark.skipif(not _HAS_BLAST, reason="BLAST+ not installed")
def test_end_to_end_real_blast_calls_planted_replicon(tmp_path):
    # synthetic replicon allele (300bp, non-repetitive) + an assembly that contains it verbatim in flanks
    allele = _varied_seq(300, seed=11)
    absent = _varied_seq(300, seed=99)
    db = tmp_path / "repdb.fsa"
    db.write_text(f">IncX3_1__SYNTH\n{allele}\n>IncFIA_1__OTHER\n{absent}\n", encoding="utf-8")
    flank = _varied_seq(400, seed=42)
    asm = tmp_path / "asm.fna"
    asm.write_text(f">contig1\n{flank}{allele}{flank}\n", encoding="utf-8")
    res = call_replicons(asm, db)
    assert res["status"] == "ok", res
    reps = {r["replicon"] for r in res["replicons"]}
    assert "IncX3" in reps, res                      # planted replicon called
    assert "IncFIA" not in reps                       # absent allele not called
    # CLI end-to-end returns 0 on a successful call
    rc = main([str(asm), "--db", str(db), "--json-only"])
    assert rc == 0


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-q"]))
