"""MLST decoder — pure profiles/lookup core + real-BLAST exact-allele -> profile -> ST e2e (synthetic scheme)."""
import json
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import pytest  # noqa: E402

from dna_decode.mlst.cli import main  # noqa: E402
from dna_decode.mlst.core import allele_number, lookup_st, parse_profiles  # noqa: E402
from dna_decode.mlst.runner import call_mlst  # noqa: E402

_HAS_BLAST = bool(shutil.which("blastn", path="C:/Users/Farshad/ncbi-blast/bin") or shutil.which("blastn")) \
    and bool(shutil.which("makeblastdb", path="C:/Users/Farshad/ncbi-blast/bin") or shutil.which("makeblastdb"))

_PROFILES = ("ST\tadk\tfumC\tclonal_complex\n"
             "1\t4\t2\tCC1\n"
             "2\t5\t3\t\n")


def _varied(n, seed):
    b = "ACGT"; x = seed; o = []
    for _ in range(n):
        x = (1103515245 * x + 12345) & 0x7FFFFFFF; o.append(b[(x >> 24) & 3])
    return "".join(o)


def test_parse_profiles_and_lookup():
    loci, st_of, cc_of = parse_profiles(_PROFILES)
    assert loci == ["adk", "fumC"]                       # clonal_complex excluded from loci
    assert st_of[(4, 2)] == "1" and st_of[(5, 3)] == "2"
    assert cc_of[(4, 2)] == "CC1"
    assert lookup_st({"adk": 4, "fumC": 2}, loci, st_of) == "1"
    assert lookup_st({"adk": 9, "fumC": 9}, loci, st_of) is None   # novel
    assert lookup_st({"adk": 4, "fumC": None}, loci, st_of) is None  # incomplete


def test_allele_number():
    assert allele_number("adk_4") == ("adk", 4)
    assert allele_number("fumC_127") == ("fumC", 127)
    assert allele_number("nounderscore") is None


def test_offline_safe(monkeypatch, tmp_path):
    import dna_decode.typing.blast_caller as bc
    monkeypatch.setattr(bc, "find_blastn", lambda: None)
    (tmp_path / "profiles.tsv").write_text(_PROFILES, encoding="utf-8")
    (tmp_path / "adk.fasta").write_text(">adk_4\nACGT\n", encoding="utf-8")
    res = call_mlst(tmp_path / "g.fna", {"adk": tmp_path / "adk.fasta"}, tmp_path / "profiles.tsv")
    assert res["status"] == "unavailable"


def test_cli_no_db_exit3(tmp_path):
    g = tmp_path / "g.fna"; g.write_text(">c\nACGT\n", encoding="utf-8")
    out = tmp_path / "m.json"
    rc = main([str(g), "--db-dir", str(tmp_path / "nope"), "--out", str(out), "--json-only"])
    assert rc == 3 and json.loads(out.read_text())["status"] == "unavailable"


@pytest.mark.skipif(not _HAS_BLAST, reason="BLAST+ not installed")
def test_e2e_calls_ST1(tmp_path):
    # synthetic 2-locus scheme; genome carries adk allele#4 + fumC allele#2 -> profile (4,2) -> ST1
    adk4 = _varied(500, 40); adk5 = _varied(500, 50)
    fumc2 = _varied(500, 20); fumc3 = _varied(500, 30)
    db = tmp_path / "scheme"; db.mkdir()
    (db / "adk.fasta").write_text(f">adk_4\n{adk4}\n>adk_5\n{adk5}\n", encoding="utf-8")
    (db / "fumC.fasta").write_text(f">fumC_2\n{fumc2}\n>fumC_3\n{fumc3}\n", encoding="utf-8")
    (db / "profiles.tsv").write_text(_PROFILES, encoding="utf-8")
    flank = _varied(300, 7)
    asm = tmp_path / "asm.fna"
    asm.write_text(f">c1\n{flank}{adk4}{flank}{fumc2}{flank}\n", encoding="utf-8")
    out = tmp_path / "m.json"
    rc = main([str(asm), "--db-dir", str(db), "--loci", "adk,fumC", "--out", str(out), "--json-only"])
    assert rc == 0
    rec = json.loads(out.read_text())
    assert rec["profile"] == {"adk": 4, "fumC": 2}, rec
    assert rec["sequence_type"] == "1" and rec["clonal_complex"] == "CC1"
    assert rec["complete"] is True and rec["novel"] is False


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-q"]))
