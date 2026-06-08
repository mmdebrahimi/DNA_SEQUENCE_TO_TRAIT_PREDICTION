"""Offline unit tests for the G2 dry-manifest gate (scripts/g2_dry_manifest).

Pure-logic on synthetic fixtures — no 19GB download, no GPU, no network. Validates the laptop-written gate
the workhorse will run on real Arabidopsis data: accession join, empirical pseudogenome-pattern resolution,
phenotype-agnostic window-table generation from a GFF, N-fraction QC, group-label presence, GREEN/RED verdict.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from scripts.g2_dry_manifest import (  # noqa: E402
    build_window_table,
    dry_manifest_report,
    group_labels_present,
    intersect_accessions,
    n_fraction,
    resolve_pseudogenome_pattern,
    window_n_fractions,
)


def test_resolve_pseudogenome_pattern_empirical():
    accs = {"10000", "100000", "5832"}
    files = ["pseudo_10000.fasta.gz", "pseudo_100000.fasta.gz", "pseudo_5832.fasta.gz", "README.txt"]
    tmpl, stats = resolve_pseudogenome_pattern(accs, files)
    assert tmpl == "pseudo_{id}.fasta.gz", tmpl
    assert stats["n_matched"] == 3
    # a different real-world scheme must also resolve (NEVER assume the pseudo{id} form)
    tmpl2, _ = resolve_pseudogenome_pattern({"6909"}, ["6909_pseudogenome.fa"])
    assert tmpl2 == "{id}_pseudogenome.fa", tmpl2


def test_intersect_accessions_reports_drops():
    pheno = {"a1", "a2", "a3", "a4"}
    pseudo = {"a1", "a2", "a3"}        # a4 has no pseudogenome
    snp = {"a1", "a2", "a4"}           # a3 not in SNP matrix
    r = intersect_accessions(pheno, pseudo, snp)
    assert r["analysis_n"] == 2 and set(r["analysis_ids"]) == {"a1", "a2"}
    assert r["dropped_no_pseudogenome"] == ["a4"]
    assert r["dropped_no_snp"] == ["a3"]


def test_build_window_table_agnostic_gene_flanks(tmp_path):
    gff = tmp_path / "ann.gff3"
    gff.write_text(
        "Chr1\tphytozome\tgene\t2000\t2300\t.\t+\t.\tID=AT1G001;Name=GENEA\n"
        "Chr1\tphytozome\tmRNA\t2000\t2300\t.\t+\t.\tID=mrna1;Parent=AT1G001\n"   # ignored (not gene)
        "Chr2\tphytozome\tgene\t500\t560\t.\t-\t.\tID=AT2G002\n",
        encoding="utf-8")
    ws = build_window_table(gff, flank=100, window=200, stride=100,
                            chrom_len={"Chr1": 10000, "Chr2": 10000})
    # gene A interval = 1900..2400 (300+2*100) tiled at window200/stride100; gene B = 400..660
    assert {w.gene_id for w in ws} == {"AT1G001", "AT2G002"}
    a = [w for w in ws if w.gene_id == "AT1G001"]
    assert a[0].start == 1900 and a[0].end == 2099           # first window
    assert all(w.end - w.start + 1 <= 200 for w in ws)       # window-size bound
    assert a[-1].end <= 2400                                  # clipped to interval end


def test_build_window_table_clips_to_chrom_len(tmp_path):
    gff = tmp_path / "ann.gff3"
    gff.write_text("Chr1\tx\tgene\t9900\t9950\t.\t+\t.\tID=EDGE\n", encoding="utf-8")
    ws = build_window_table(gff, flank=500, window=200, stride=200, chrom_len={"Chr1": 10000})
    assert ws and max(w.end for w in ws) <= 10000            # never past chromosome end
    assert min(w.start for w in ws) >= 1


def test_n_fraction_and_window_qc(tmp_path):
    assert n_fraction("ACGT") == 0.0
    assert n_fraction("ACNN") == 0.5
    assert n_fraction("") == 1.0
    gff = tmp_path / "g.gff3"
    gff.write_text("Chr1\tx\tgene\t5\t8\t.\t+\t.\tID=G1\n", encoding="utf-8")
    ws = build_window_table(gff, flank=2, window=4, stride=4, chrom_len={"Chr1": 20})
    clean = {"Chr1": "ACGT" * 5}                              # 20bp, no Ns
    qc = window_n_fractions(ws, clean, max_n=0.10)
    assert qc["n_flagged"] == 0 and qc["mean_n_fraction"] == 0.0
    dirty = {"Chr1": "NNNN" + "ACGT" * 4}                     # first window all-N
    qc2 = window_n_fractions(ws, dirty, max_n=0.10)
    assert qc2["n_flagged"] >= 1


def test_group_labels_present():
    ids = ["a1", "a2", "a3"]
    g = group_labels_present(ids, {"a1": "germany", "a2": "germany", "a3": "iberia"})
    assert g["n_missing"] == 0 and g["n_groups"] == 2 and g["groups"]["germany"] == 2
    g2 = group_labels_present(ids, {"a1": "germany"})
    assert g2["n_missing"] == 2 and sorted(g2["missing"]) == ["a2", "a3"]


def test_dry_manifest_verdict_green_and_red():
    inter = {"analysis_n": 150, "analysis_ids": [f"a{i}" for i in range(150)]}
    pstats = {"template": "pseudo_{id}.fa", "n_matched": 150}
    ws = build_window_table.__wrapped__ if hasattr(build_window_table, "__wrapped__") else None  # noqa
    # build a small real window list instead of mocking
    from scripts.g2_dry_manifest import Window
    wlist = [Window("Chr1", 1, 512, "G1", 0)]
    groups = {"n_missing": 0, "n_groups": 9, "groups": {}, "missing": []}
    green = dry_manifest_report(intersect=inter, pattern_stats=pstats, window_table=wlist,
                               n_qc=None, groups=groups)
    assert green["verdict"] == "GREEN" and green["red_checks"] == []
    # RED: too few accessions + missing groups
    red = dry_manifest_report(intersect={"analysis_n": 40, "analysis_ids": ["a"]},
                              pattern_stats={"template": None, "n_matched": 0},
                              window_table=wlist,
                              n_qc={"n_flagged": 3, "flagged": ["w"], "mean_n_fraction": 0.2,
                                    "n_windows": 1, "max_n_threshold": 0.1},
                              groups={"n_missing": 5, "n_groups": 1, "groups": {}, "missing": ["x"]})
    assert red["verdict"] == "RED"
    assert set(red["red_checks"]) >= {"accession_intersection", "pseudogenome_pattern",
                                      "n_fraction_qc", "group_labels"}


if __name__ == "__main__":
    import pytest
    sys.exit(pytest.main([__file__, "-q"]))
