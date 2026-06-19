"""Single-genome CLI tests — render + offline end-to-end (no live tools)."""
from __future__ import annotations

from pathlib import Path

from dna_decode.genome_map import (
    TIER_DETERMINANT_PHENOTYPE,
    TIER_HOMOLOGY_HYPOTHESIS,
    TIER_UNKNOWN,
)
from dna_decode.genome_map.gate import evaluate_gate
from scripts.genome_map import main, render_genome_summary_md
from scripts.genome_map_spike import run_genome_map_for


def _offline_gff(tmp_path: Path) -> Path:
    gff = tmp_path / "g.gff3"
    gff.write_text(
        "##gff-version 3\n"
        "##sequence-region contig_1 1 5000\n"
        "contig_1\tBakta\tCDS\t100\t900\t.\t+\t0\tID=c1;gene=gyrA;locus_tag=T1;product=DNA gyrase subunit A\n"
        "contig_1\tBakta\tCDS\t1000\t1500\t.\t+\t0\tID=c2;locus_tag=T2;product=putative transporter\n"
        "contig_1\tBakta\tCDS\t2000\t2500\t.\t+\t0\tID=c3;locus_tag=T3;product=hypothetical protein\n"
        "##FASTA\n>contig_1\nACGT\n",
        encoding="utf-8",
    )
    return gff


def test_render_summary_contains_sections(tmp_path: Path):
    gff = _offline_gff(tmp_path)
    main_tsv = tmp_path / "main.tsv"
    main_tsv.write_text(
        "Protein id\tContig id\tStart\tStop\tElement symbol\tElement name\tClass\tSubclass\tMethod\n"
        "T1\tcontig_1\t100\t900\tgyrA_S83L\tgyrase\tQUINOLONE\tQUINOLONE\tPOINTX\n",
        encoding="utf-8",
    )
    gm = run_genome_map_for("S1", "Escherichia", gff, main_tsv, ["ciprofloxacin"])
    md = render_genome_summary_md(gm, evaluate_gate(gm), generated="2026-06-19")
    assert "Genome map — S1" in md
    assert "unknown_under_bakta_db_light" in md
    assert "Determinant-phenotype features" in md
    assert "Honesty gate" in md


def test_cli_offline_no_amrfinder(tmp_path: Path):
    gff = _offline_gff(tmp_path)
    out = tmp_path / "out"
    rc = main([
        "--gff", str(gff), "--no-amrfinder", "--sample-id", "S2",
        "--organism", "Escherichia", "--out-dir", str(out),
    ])
    assert rc == 0
    import json
    gm = json.loads((out / "genome_map_S2.json").read_text())
    assert gm["degraded_coverage"] is True
    # tiers still assigned from the GFF; no determinant phenotype
    assert gm["metrics"]["determinant_phenotype_feature_count"] == 0
    tiers = [f["primary_tier"] for f in gm["features"]]
    assert TIER_HOMOLOGY_HYPOTHESIS in tiers  # putative transporter
    assert TIER_UNKNOWN in tiers              # hypothetical protein
    assert (out / "genome_map_S2.md").exists()
    assert (out / "genome_map_S2_table.json").exists()


def test_cli_requires_an_input():
    assert main(["--no-amrfinder"]) == 2  # neither --genome-fasta nor --gff


def test_render_degraded_banner(tmp_path: Path):
    gff = _offline_gff(tmp_path)
    from dna_decode.genome_map import ingest
    from dna_decode.genome_map.build_map import build_genome_map

    features = ingest.load_genome_gff(gff)
    gm = build_genome_map("S3", "Escherichia", features, [],
                          {"n_main_rows": 0, "n_high_confidence_join": 0,
                           "n_symbol_fallback": 0, "n_unjoined": 0}, degraded=True)
    md = render_genome_summary_md(gm, evaluate_gate(gm))
    assert "DEGRADED COVERAGE" in md
