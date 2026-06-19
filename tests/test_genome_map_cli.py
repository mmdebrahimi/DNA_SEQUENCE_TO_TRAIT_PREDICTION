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


def test_cli_offline_overlay_status(tmp_path: Path):
    import json
    gff = _offline_gff(tmp_path)
    out = tmp_path / "out"
    main(["--gff", str(gff), "--no-amrfinder", "--sample-id", "S", "--out-dir", str(out)])
    gm = json.loads((out / "genome_map_S.json").read_text())
    assert gm["overlay_status"] == "OFFLINE_NO_AMRFINDER"


def test_cli_gff_only_no_fasta_status(tmp_path: Path):
    import json
    gff = _offline_gff(tmp_path)
    out = tmp_path / "out"
    # only --gff (no --genome-fasta, no --no-amrfinder) -> nothing to scan -> GFF_ONLY_NO_FASTA
    main(["--gff", str(gff), "--sample-id", "S", "--out-dir", str(out)])
    gm = json.loads((out / "genome_map_S.json").read_text())
    assert gm["overlay_status"] == "GFF_ONLY_NO_FASTA"
    assert gm["degraded_coverage"] is True


def test_cli_live_amrfinder_failure_blocks_by_default(tmp_path: Path, monkeypatch):
    # LIVE mode (a FASTA supplied) + AMRFinder failure -> BLOCKED (exit 3), no map emitted
    from dna_decode.genome_map import amrfinder

    def _boom(*a, **k):
        raise RuntimeError("docker down")

    monkeypatch.setattr(amrfinder, "run_amrfinder", _boom)
    gff = _offline_gff(tmp_path)
    fasta = tmp_path / "g.fna"
    fasta.write_text(">c\nACGT\n", encoding="utf-8")
    out = tmp_path / "out"
    rc = main(["--genome-fasta", str(fasta), "--gff", str(gff), "--sample-id", "S",
               "--out-dir", str(out)])
    assert rc == 3
    assert not (out / "genome_map_S.json").exists()  # fail-closed: no silent no-determinant map


def test_cli_live_amrfinder_failure_allow_degraded(tmp_path: Path, monkeypatch):
    import json
    from dna_decode.genome_map import amrfinder

    def _boom(*a, **k):
        raise RuntimeError("docker down")

    monkeypatch.setattr(amrfinder, "run_amrfinder", _boom)
    gff = _offline_gff(tmp_path)
    fasta = tmp_path / "g.fna"
    fasta.write_text(">c\nACGT\n", encoding="utf-8")
    out = tmp_path / "out"
    rc = main(["--genome-fasta", str(fasta), "--gff", str(gff), "--allow-degraded",
               "--sample-id", "S", "--out-dir", str(out)])
    assert rc == 0
    gm = json.loads((out / "genome_map_S.json").read_text())
    assert gm["overlay_status"] == "DEGRADED_USER_ACCEPTED"
    assert gm["degraded_coverage"] is True
    assert gm["metrics"]["determinant_phenotype_feature_count"] == 0


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
