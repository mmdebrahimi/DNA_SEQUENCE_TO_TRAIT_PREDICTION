"""Step 7 tests — spike aggregation + verdict rendering (synthetic, no live tools)."""
from __future__ import annotations

from pathlib import Path

import pandas as pd

from dna_decode.genome_map import (
    TIER_DETERMINANT_PHENOTYPE,
    TIER_HOMOLOGY_HYPOTHESIS,
    TIER_UNKNOWN,
)
from dna_decode.genome_map.gate import VERDICT_GO, VERDICT_NO_GO, evaluate_gate
from scripts.genome_map_spike import (
    VERDICT_BLOCKED,
    bakta_contig_lengths,
    fasta_contig_lengths,
    render_verdict_md,
    run_genome_map_for,
    summarize_spike,
)


def _map_with_features(features, *, all_symbol_fallback=False, acc="GCA_X"):
    return {
        "genome_accession": acc,
        "amrfinder_organism": "Escherichia",
        "features": features,
        "metrics": {
            "total_features": len(features),
            "per_tier_counts": {"determinant-phenotype": 0, "curated-molecular-function": 0,
                                "homology-only-hypothesis": 0, "unknown": 0},
            "unknown_under_bakta_db_light": 0.1,
            "determinant_phenotype_feature_count": sum(
                1 for f in features if f["primary_tier"] == TIER_DETERMINANT_PHENOTYPE),
            "join_quality": {"n_main_rows": 1, "n_high_confidence_join": 0 if all_symbol_fallback else 1,
                             "n_symbol_fallback": 1 if all_symbol_fallback else 0, "n_unjoined": 0},
            "all_joins_symbol_fallback": all_symbol_fallback,
            "genome_level_calls": {"ciprofloxacin": {"prediction": "R"}},
        },
    }


def _feat(i, tier, raw_product="", phenotype=None):
    return {"feature_index": i, "primary_tier": tier, "raw_product": raw_product,
            "raw_gene_symbol": "", "classification_reason": "", "phenotype": phenotype or []}


def _passing_map(acc):
    feats = [
        _feat(0, TIER_HOMOLOGY_HYPOTHESIS, "putative a"),
        _feat(1, TIER_HOMOLOGY_HYPOTHESIS, "putative b"),
        _feat(2, TIER_DETERMINANT_PHENOTYPE, "DNA gyrase subunit A",
              phenotype=[{"drug": "ciprofloxacin", "phenotype": "R"}]),
    ]
    return _map_with_features(feats, acc=acc)


# ---- summarize_spike ----


def test_summarize_go():
    e1 = {"accession": "G1", "status": "OK", "genome_map": _passing_map("G1")}
    e1["gate_result"] = evaluate_gate(e1["genome_map"])
    e2 = {"accession": "G2", "status": "OK", "genome_map": _passing_map("G2")}
    e2["gate_result"] = evaluate_gate(e2["genome_map"])
    agg = summarize_spike([e1, e2])
    assert agg["verdict"] == VERDICT_GO


def test_summarize_blocked_when_any_genome_blocked():
    e1 = {"accession": "G1", "status": "OK", "genome_map": _passing_map("G1")}
    e1["gate_result"] = evaluate_gate(e1["genome_map"])
    e2 = {"accession": "G2", "status": "BAKTA_ANNOTATION_BLOCKED"}
    agg = summarize_spike([e1, e2])
    assert agg["verdict"] == VERDICT_BLOCKED
    assert any("G2" in r for r in agg["reasons"])


def test_summarize_no_go_when_all_symbol_fallback():
    e1 = {"accession": "G1", "status": "OK", "genome_map": _passing_map("G1")}
    e1["gate_result"] = evaluate_gate(e1["genome_map"])
    trap = _map_with_features(
        [_feat(i, TIER_HOMOLOGY_HYPOTHESIS, "putative x") for i in range(3)],
        all_symbol_fallback=True, acc="G2")
    e2 = {"accession": "G2", "status": "OK", "genome_map": trap}
    e2["gate_result"] = evaluate_gate(trap)
    agg = summarize_spike([e1, e2])
    assert agg["verdict"] == VERDICT_NO_GO


# ---- render ----


def test_render_verdict_md_contains_sections():
    e1 = {"accession": "G1", "label": "test", "status": "OK", "genome_map": _passing_map("G1")}
    e1["gate_result"] = evaluate_gate(e1["genome_map"])
    agg = summarize_spike([e1])
    md = render_verdict_md([e1], agg, generated="2026-06-18")
    assert "GO/NO-GO verdict" in md
    assert "unknown_under_bakta_db_light" in md
    assert "G1" in md and "G2" in md
    assert "join quality" in md


def test_render_blocked_genome_no_fake_map():
    e = {"accession": "G2", "label": "x", "status": "BAKTA_ANNOTATION_BLOCKED"}
    agg = summarize_spike([e])
    md = render_verdict_md([e], agg, generated="2026-06-18")
    assert "BAKTA_ANNOTATION_BLOCKED" in md
    assert agg["verdict"] == VERDICT_BLOCKED


# ---- contig-length helpers ----


def test_fasta_contig_lengths(tmp_path: Path):
    fa = tmp_path / "g.fna"
    fa.write_text(">CP012345.1 chromosome\nACGTACGTAC\nGGGG\n>CP012346.1\nACGT\n", encoding="utf-8")
    lens = fasta_contig_lengths(fa)
    assert lens["CP012345.1"] == 14
    assert lens["CP012346.1"] == 4


def test_bakta_contig_lengths(tmp_path: Path):
    gff = tmp_path / "g.gff3"
    gff.write_text(
        "##gff-version 3\n"
        "##sequence-region contig_1 1 14\n"
        "##sequence-region contig_2 1 4\n"
        "contig_1\tBakta\tCDS\t1\t9\t.\t+\t0\tID=c1\n"
        "##FASTA\n>contig_1\nACGT\n",
        encoding="utf-8",
    )
    lens = bakta_contig_lengths(gff)
    assert lens == {"contig_1": 14, "contig_2": 4}


# ---- run_genome_map_for end-to-end on synthetic files ----


def test_run_genome_map_for_offline_files(tmp_path: Path):
    gff = tmp_path / "g.gff3"
    gff.write_text(
        "##gff-version 3\n"
        "##sequence-region contig_1 1 5000\n"
        "contig_1\tBakta\tCDS\t100\t900\t.\t+\t0\tID=cds1;gene=gyrA;locus_tag=T1;product=DNA gyrase subunit A\n"
        "contig_1\tBakta\tCDS\t1000\t1500\t.\t+\t0\tID=cds2;locus_tag=T2;product=putative transporter\n"
        "contig_1\tBakta\tCDS\t2000\t2500\t.\t+\t0\tID=cds3;locus_tag=T3;product=hypothetical protein\n"
        "##FASTA\n>contig_1\nACGT\n",
        encoding="utf-8",
    )
    fasta = tmp_path / "genome.fna"
    fasta.write_text(">CP012345.1\n" + "A" * 5000 + "\n", encoding="utf-8")
    main_tsv = tmp_path / "main.tsv"
    main_tsv.write_text(
        "Protein identifier\tContig id\tStart\tStop\tElement symbol\tElement name\tClass\tSubclass\tMethod\n"
        "T1\tCP012345.1\t100\t900\tgyrA_S83L\tgyrase\tQUINOLONE\tQUINOLONE\tPOINTX\n",
        encoding="utf-8",
    )
    gm = run_genome_map_for("GCA_T", "Escherichia", gff, main_tsv,
                            drugs=["ciprofloxacin"], fasta_path=fasta)
    # gyrA feature joined by protein-id (T1) -> determinant-phenotype
    f0 = gm["features"][0]
    assert f0["primary_tier"] == TIER_DETERMINANT_PHENOTYPE
    assert f0["phenotype"]
    # putative -> homology; hypothetical -> unknown
    assert gm["features"][1]["primary_tier"] == TIER_HOMOLOGY_HYPOTHESIS
    assert gm["features"][2]["primary_tier"] == TIER_UNKNOWN
    assert gm["metrics"]["join_quality"]["n_high_confidence_join"] == 1
