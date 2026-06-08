"""Co-localization — engine positions-mode + pure core + real-BLAST same/different-contig e2e."""
import json
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import pytest  # noqa: E402

from dna_decode.colocalization.cli import main  # noqa: E402
from dna_decode.colocalization.core import colocalize  # noqa: E402
from dna_decode.typing.blast_caller import call_alleles  # noqa: E402

_HAS_BLAST = bool(shutil.which("blastn", path="C:/Users/Farshad/ncbi-blast/bin") or shutil.which("blastn")) \
    and bool(shutil.which("makeblastdb", path="C:/Users/Farshad/ncbi-blast/bin") or shutil.which("makeblastdb"))


def _varied(n, seed):
    b = "ACGT"; x = seed; o = []
    for _ in range(n):
        x = (1103515245 * x + 12345) & 0x7FFFFFFF; o.append(b[(x >> 24) & 3])
    return "".join(o)


def test_core_links_same_contig_gene_to_replicon():
    genes = [{"gene": "blaNDM-1", "contig": "c1"}, {"gene": "sul1", "contig": "c2"}]
    replicons = [{"replicon": "IncX3", "contig": "c1"}]
    co = colocalize(genes, replicons)
    by = {g["gene"]: g for g in co["gene_calls"]}
    assert by["blaNDM-1"]["plasmid_borne"] is True and by["blaNDM-1"]["replicons_on_contig"] == ["IncX3"]
    assert by["sul1"]["plasmid_borne"] is False         # different contig, no replicon
    assert co["summary"]["n_plasmid_borne"] == 1 and co["summary"]["n_chromosomal_or_unplaced"] == 1


def test_core_empty():
    co = colocalize([], [])
    assert co["summary"]["n_genes"] == 0 and co["gene_calls"] == []


def test_cli_missing_fasta_exit2():
    assert main(["/nonexistent.fna"]) == 2


@pytest.mark.skipif(not _HAS_BLAST, reason="BLAST+ not installed")
def test_engine_positions_mode_returns_contig(tmp_path):
    allele = _varied(400, 1)
    db = tmp_path / "a.fsa"; db.write_text(f">x_1_ACC\n{allele}\n", encoding="utf-8")
    asm = tmp_path / "g.fna"; asm.write_text(f">myContig\n{_varied(200,3)}{allele}{_varied(200,4)}\n",
                                             encoding="utf-8")
    res = call_alleles(asm, db, identity_threshold=90, coverage_threshold=60, with_positions=True)
    hit = res["per_allele"]["x_1_ACC"]
    assert hit["called"] and hit["contig"] == "myContig"
    assert isinstance(hit["sstart"], int) and isinstance(hit["send"], int)


@pytest.mark.skipif(not _HAS_BLAST, reason="BLAST+ not installed")
def test_e2e_gene_on_plasmid_contig_is_plasmid_borne(tmp_path):
    gene = _varied(600, 1)       # blaNDM-1
    replicon = _varied(600, 2)   # IncX3
    chrom_gene = _varied(600, 3)  # sul1 on a replicon-free contig
    res_dir = tmp_path / "resfinder_db"; res_dir.mkdir()
    (res_dir / "beta-lactam.fsa").write_text(f">blaNDM-1_1_ACC\n{gene}\n", encoding="utf-8")
    (res_dir / "sulphonamide.fsa").write_text(f">sul1_1_ACC\n{chrom_gene}\n", encoding="utf-8")
    pdb = tmp_path / "plasmid.fsa"; pdb.write_text(f">IncX3_1__ACC\n{replicon}\n", encoding="utf-8")
    # contig 'plasmid1' carries the replicon + blaNDM-1; contig 'chrom1' carries only sul1
    asm = tmp_path / "asm.fna"
    asm.write_text(f">plasmid1\n{_varied(150,7)}{replicon}{_varied(150,8)}{gene}{_varied(150,9)}\n"
                   f">chrom1\n{_varied(150,5)}{chrom_gene}{_varied(150,6)}\n", encoding="utf-8")
    out = tmp_path / "co.json"
    rc = main([str(asm), "--plasmid-db", str(pdb), "--resfinder-db-dir", str(res_dir),
               "--out", str(out), "--json-only"])
    assert rc == 0
    rec = json.loads(out.read_text())
    by = {g["gene"]: g for g in rec["gene_calls"]}
    assert by["blaNDM-1"]["plasmid_borne"] is True and "IncX3" in by["blaNDM-1"]["replicons_on_contig"]
    assert by["sul1"]["plasmid_borne"] is False
    assert rec["summary"]["n_plasmid_borne"] == 1


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-q"]))
