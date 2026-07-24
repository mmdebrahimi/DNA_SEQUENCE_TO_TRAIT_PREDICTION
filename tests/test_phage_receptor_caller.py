"""Offline tests for the phage receptor caller (dna_decode/phage/receptor_caller.py).

No real BLAST: the db-fasta relabelling + exclusion, the leave-one-out accounting (via a stubbed
call), the manifest->receptor loader (via the catalog lineage lookup), and the no-blastn INDETERMINATE
degrade are all exercised deterministically.
"""
from __future__ import annotations

from pathlib import Path

import dna_decode.phage.receptor_caller as caller


def _write_fna(p: Path, seqs: dict[str, str]) -> None:
    with open(p, "w", encoding="utf-8") as fh:
        for sid, seq in seqs.items():
            fh.write(f">{sid} desc\n{seq}\n")


def test_labeled_db_fasta_relabels_and_excludes(tmp_path):
    a = tmp_path / "A.fna"; _write_fna(a, {"contig1": "ACGT", "contig2": "TTTT"})
    b = tmp_path / "B.fna"; _write_fna(b, {"x": "GGGG"})
    out = tmp_path / "db.fna"
    id_to_label = caller._write_labeled_db_fasta({"A": str(a), "B": str(b)}, out, exclude="B")
    # B excluded; A's two contigs both map back to label "A"
    assert set(id_to_label.values()) == {"A"}
    assert id_to_label["A"] == "A" and id_to_label["A__1"] == "A"
    assert "B" not in id_to_label


def test_no_blastn_degrades_to_indeterminate(tmp_path, monkeypatch):
    monkeypatch.setattr(caller, "find_blastn", lambda: None)
    q = tmp_path / "q.fna"; _write_fna(q, {"q": "ACGT"})
    call = caller.call_receptor(str(q), {"A": str(q)}, {"A": "OmpC"})
    assert call.status == "INDETERMINATE" and call.predicted_receptor is None
    assert "blastn" in call.reason.lower()


def test_leave_one_out_accounting(tmp_path, monkeypatch):
    # 4 phages: 2 OmpC, 2 BtuB. Stub call_receptor: a phage's nearest is the OTHER same-receptor one
    # for OmpC (correct), but BtuB phages get mispredicted as OmpC (wrong) -> tests per-receptor split.
    refs = {p: str(tmp_path / f"{p}.fna") for p in ("o1", "o2", "b1", "b2")}
    receptors = {"o1": "OmpC", "o2": "OmpC", "b1": "BtuB", "b2": "BtuB"}
    for p in refs:
        _write_fna(Path(refs[p]), {p: "ACGT"})

    def stub(query, r, rec, *, exclude=None, blastn_bin=None):
        pred = "OmpC"  # everything predicted OmpC
        return caller.ReceptorCall("CALLED", pred, "x", 99.0, 100.0)
    monkeypatch.setattr(caller, "call_receptor", stub)

    res = caller.leave_one_out(refs, receptors)
    assert res.n_total == 4 and res.n_called == 4
    assert res.n_correct == 2                       # the 2 OmpC are right; 2 BtuB wrong
    assert res.accuracy == 0.5
    assert res.per_receptor["OmpC"] == [2, 2]       # 2 correct / 2 called
    assert res.per_receptor["BtuB"] == [0, 2]       # 0 correct / 2 called


def test_load_manifest_attaches_receptor_via_lineage(tmp_path):
    gdir = tmp_path / "genomes"; gdir.mkdir()
    # two genomes: one Tequatrovirus (OmpC), one uncatalogued genus (dropped)
    for acc in ("MZ000001.1", "MZ000002.1", "MZ000003.1"):
        _write_fna(gdir / f"{acc}.fna", {acc: "ACGT"})
    man = tmp_path / "manifest.tsv"
    man.write_text(
        "accession\tphage_name\tgenus\tfamily\n"
        "MZ000001.1\tSomeQuintavirus\tVequintavirus\tunclassified\n"   # clade-conserved -> ECA
        "MZ000002.1\tSomeT4like\tTequatrovirus\tStraboviridae\n"       # RBP-variable -> dropped
        "MZ000003.1\tMystery\tNovelgenus\tNovelviridae\n",             # uncatalogued -> dropped
        encoding="utf-8")
    refs, receptors = caller._load_manifest(str(man), str(gdir))
    assert set(refs) == {"MZ000001.1"}              # RBP-variable + uncatalogued dropped, not guessed
    assert receptors["MZ000001.1"] == "ECA"
