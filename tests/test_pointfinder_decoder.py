"""PointFinder chromosomal point-mutation decoder — overview parse + codon-map + real-BLAST S83L e2e."""
import json
import shutil
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import pytest  # noqa: E402

from dna_decode.pointfinder.cli import main  # noqa: E402
from dna_decode.pointfinder.runner import call_point_mutations, parse_overview  # noqa: E402
from dna_decode.typing.codon_map import subject_aa_by_codon, translate  # noqa: E402

_HAS_BLAST = bool(shutil.which("blastn", path="C:/Users/Farshad/ncbi-blast/bin") or shutil.which("blastn")) \
    and bool(shutil.which("makeblastdb", path="C:/Users/Farshad/ncbi-blast/bin") or shutil.which("makeblastdb"))

# varied, stop-free codons so blastn aligns cleanly (low-complexity is dust-masked)
_CYC = ["GCT", "GAT", "TTT", "GGT", "CAT", "ATT", "AAA", "CTT", "CCT", "CAA",
        "CGT", "ACT", "GTT", "TGG", "GTG", "AAT", "GAA", "GGA", "CTG", "TAT"]


def _cds(n_codons=200, at=None):
    """In-frame CDS: codon1=ATG, others cycled (stop-free); `at` = {1based_pos: codon} overrides."""
    codons = ["ATG"]
    for i in range(1, n_codons):
        codons.append(_CYC[i % len(_CYC)])
    for pos, cod in (at or {}).items():
        codons[pos - 1] = cod
    return "".join(codons)


def test_parse_overview(tmp_path):
    ov = tmp_path / "resistens-overview.txt"
    ov.write_text(
        "#Gene_ID\tGene_name\tCodon_pos\tRef_nuc\tRef_codon\tRes_codon\tResistance\tPMID\tMechanism\tNotes\tRequired_mut\n"
        "gyrA\tgyrA\t83\tTCG\tS\tL,W,A\tNalidixic acid,Ciprofloxacin\t123\tTarget mod\t\t\n"
        "gyrA\tgyrA\t87\tGAC\tD\tN\tCiprofloxacin\t456\tTarget mod\t\tgyrA_S83L\n",
        encoding="utf-8")
    o = parse_overview(ov)
    assert o[("gyrA", 83)]["ref_aa"] == "S"
    assert o[("gyrA", 83)]["res_aas"] == {"L", "W", "A"}
    assert "Ciprofloxacin" in o[("gyrA", 83)]["resistances"]
    assert o[("gyrA", 87)]["required"] == {"gyrA_S83L"}


def test_subject_aa_by_codon_ungapped():
    ref = _cds(100, at={83: "TCG"})        # S at 83
    # subject identical except codon 83 -> CTG (L)
    subj = _cds(100, at={83: "CTG"})
    m = subject_aa_by_codon(ref, subj, 1, 100)
    assert m[83] == "L" and translate(ref)[82] == "S"


def test_offline_safe(monkeypatch, tmp_path):
    import dna_decode.pointfinder.runner as r
    monkeypatch.setattr(r, "find_blastn", lambda: None)
    res = call_point_mutations(tmp_path / "g.fna", {"gyrA": tmp_path / "gyrA.fsa"}, {})
    assert res["status"] == "unavailable" and res["mutations"] == []


def test_cli_missing_fasta_exit2():
    assert main(["/nonexistent.fna"]) == 2


@pytest.mark.skipif(not _HAS_BLAST, reason="BLAST+ not installed")
def test_e2e_calls_S83L(tmp_path):
    ref_cds = _cds(200, at={83: "TCG"})        # gyrA ref: Ser83
    db = tmp_path / "escherichia_coli"; db.mkdir()
    (db / "gyrA.fsa").write_text(f">gyrA_1_REF\n{ref_cds}\n", encoding="utf-8")
    (db / "resistens-overview.txt").write_text(
        "#Gene_ID\tGene_name\tCodon_pos\tRef_nuc\tRef_codon\tRes_codon\tResistance\tPMID\tMechanism\tNotes\tRequired_mut\n"
        "gyrA\tgyrA\t83\tTCG\tS\tL\tCiprofloxacin\t1\tTarget mod\t\t\n", encoding="utf-8")
    flank = "ACGTACGT" * 60
    # resistant genome: codon 83 TCG(S) -> CTG(L)
    res_cds = _cds(200, at={83: "CTG"})
    (tmp_path / "R.fna").write_text(f">contig1\n{flank}{res_cds}{flank}\n", encoding="utf-8")
    (tmp_path / "S.fna").write_text(f">contig1\n{flank}{ref_cds}{flank}\n", encoding="utf-8")

    out = tmp_path / "r.json"
    rc = main([str(tmp_path / "R.fna"), "--db-dir", str(db), "--genes", "gyrA",
               "--out", str(out), "--json-only"])
    assert rc == 0
    rec = json.loads(out.read_text())
    assert "gyrA S83L".split()[1] in rec["mutations"] or "S83L" in rec["mutations"], rec
    assert "Ciprofloxacin" in rec["resistances"]
    # wild-type genome -> no mutation
    out2 = tmp_path / "s.json"
    main([str(tmp_path / "S.fna"), "--db-dir", str(db), "--genes", "gyrA", "--out", str(out2), "--json-only"])
    rec2 = json.loads(out2.read_text())
    assert rec2["mutations"] == [] and "gyrA" in rec2["genes_aligned"]


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-q"]))
