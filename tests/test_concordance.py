"""AMR cross-tool concordance — pure-core units + main.tsv parse + real-BLAST gene-set e2e."""
import json
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import pytest  # noqa: E402

from dna_decode.concordance.cli import main  # noqa: E402
from dna_decode.concordance.core import (  # noqa: E402
    amr_acquired_genes_from_main,
    compare,
    family_normalize,
)

_HAS_BLAST = bool(shutil.which("blastn", path="C:/Users/Farshad/ncbi-blast/bin") or shutil.which("blastn")) \
    and bool(shutil.which("makeblastdb", path="C:/Users/Farshad/ncbi-blast/bin") or shutil.which("makeblastdb"))


def test_family_normalize_keeps_distinct_genes_merges_variants():
    assert family_normalize("blaNDM-19") == "blandm"
    assert family_normalize("blaNDM-1") == "blandm"          # same family as NDM-19
    assert family_normalize("blaCTX-M-15") == "blactx-m"
    assert family_normalize("sul1") == "sul1"                # NOT merged with sul2 (no hyphen-number)
    assert family_normalize("sul2") == "sul2"
    assert family_normalize("qnrS1") == "qnrs1"
    assert family_normalize("tet(A)") == "tet(a)"


def test_compare_buckets_and_agreement():
    amr = ["blaCTX-M-15", "sul1", "aac(6')-Ib-cr"]
    res = ["blaCTX-M-27", "sul2", "tet(A)"]   # CTX-M family matches; sul1!=sul2; rest disjoint
    c = compare(amr, res)
    assert "blactx-m" in c["both"]
    assert "sul1" in c["amr_only"] and "sul2" in c["resfinder_only"]
    assert c["n_both"] == 1
    assert 0.0 < c["agreement"] < 1.0


def test_compare_empty():
    c = compare([], [])
    assert c["agreement"] is None and c["both"] == []


def test_amr_acquired_genes_from_main_excludes_point(tmp_path):
    tsv = tmp_path / "main.tsv"
    tsv.write_text(
        "Gene symbol\tElement type\tMethod\n"
        "blaCTX-M-15\tAMR\tBLASTX\n"
        "gyrA_S83L\tAMR\tPOINTX\n"          # point mutation -> excluded
        "sul1\tAMR\tEXACTX\n"
        "stx2\tVIRULENCE\tEXACTX\n",        # not AMR -> excluded
        encoding="utf-8")
    genes = amr_acquired_genes_from_main(str(tsv))
    assert genes == {"blaCTX-M-15", "sul1"}


def test_cli_gene_set_mode(tmp_path, capsys):
    out = tmp_path / "c.json"
    rc = main(["--amr-genes", "blaCTX-M-15,sul1", "--resfinder-genes", "blaCTX-M-27_1_X,tet(A)_6_Y",
               "--out", str(out), "--json-only"])
    assert rc == 0
    rec = json.loads(out.read_text())
    assert "blactx-m" in rec["both"] and rec["schema"] == "amr-concordance-v0"


@pytest.mark.skipif(not _HAS_BLAST, reason="BLAST+ not installed")
def test_cli_real_blast_with_amrfinder_run(tmp_path):
    def varied(n, seed):
        b = "ACGT"; x = seed; o = []
        for _ in range(n):
            x = (1103515245 * x + 12345) & 0x7FFFFFFF; o.append(b[(x >> 24) & 3])
        return "".join(o)
    gene = varied(600, 1)
    db = tmp_path / "beta-lactam.fsa"
    db.write_text(f">blaNDM-1_1_ACC\n{gene}\n", encoding="utf-8")
    asm = tmp_path / "asm.fna"
    asm.write_text(f">c\n{varied(300,9)}{gene}{varied(300,8)}\n", encoding="utf-8")
    run = tmp_path / "amrrun"; run.mkdir()
    (run / "main.tsv").write_text("Gene symbol\tElement type\tMethod\nblaNDM-1\tAMR\tBLASTX\n", encoding="utf-8")
    out = tmp_path / "c.json"
    rc = main([str(asm), "--amrfinder-run", str(run), "--resfinder-db-dir", str(tmp_path),
               "--out", str(out), "--json-only"])
    assert rc == 0
    rec = json.loads(out.read_text())
    assert "blandm" in rec["both"] and rec["agreement"] == 1.0


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-q"]))
