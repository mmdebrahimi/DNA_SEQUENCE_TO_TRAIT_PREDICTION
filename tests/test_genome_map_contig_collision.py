"""Contig length-reconciliation diagnostics (the deferred honesty-visibility item).

`contig_collision_count` exposes WHY some determinant hits fall to symbol-fallback: the
AMRFinder(FASTA-named) vs Bakta-renamed contigs are reconciled by UNIQUE length, and
length-ambiguous contigs can't reconcile -> their hits go to symbol-fallback (already
surfaced via n_symbol_fallback). This makes the REASON auditable. Pure + offline.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from dna_decode.genome_map.phenotype_overlay import (  # noqa: E402
    build_contig_name_map, contig_collision_count,
)


def test_no_collision_all_reconcile():
    fasta = {"CP1": 100, "CP2": 200, "CP3": 300}
    bakta = {"contig_1": 100, "contig_2": 200, "contig_3": 300}
    cc = contig_collision_count(fasta, bakta)
    assert cc["n_reconciled"] == 3
    assert cc["n_fasta_ambiguous_contigs"] == 0 and cc["n_bakta_ambiguous_contigs"] == 0


def test_length_collision_is_counted_and_unreconciled():
    # CP2 + CP3 share length 200 -> ambiguous on the FASTA side -> cannot reconcile by length.
    fasta = {"CP1": 100, "CP2": 200, "CP3": 200}
    bakta = {"contig_1": 100, "contig_2": 200, "contig_3": 200}
    cc = contig_collision_count(fasta, bakta)
    assert cc["n_fasta_ambiguous_contigs"] == 2  # CP2 + CP3
    assert cc["n_bakta_ambiguous_contigs"] == 2
    assert cc["n_reconciled"] == 1               # only the unique-length CP1 reconciles
    # consistency with the reconciler itself
    assert len(build_contig_name_map(fasta, bakta)) == 1


def test_renderer_surfaces_contig_reconciliation():
    import scripts.genome_map as gm
    genome_map = {
        "genome_accession": "GCA_x", "amrfinder_organism": "Escherichia",
        "metrics": {
            "total_features": 5, "per_tier_counts": {"unknown": 1},
            "unknown_under_bakta_db_light": 0.2,
            "join_quality": {"n_main_rows": 3, "n_high_confidence_join": 1,
                             "n_symbol_fallback": 2, "n_unjoined": 0},
            "contig_reconciliation": {"n_fasta_contigs": 3, "n_bakta_contigs": 3,
                                      "n_reconciled": 1, "n_fasta_ambiguous_contigs": 2,
                                      "n_bakta_ambiguous_contigs": 2},
            "determinant_phenotype_feature_count": 0, "determinant_phenotype_features": [],
            "genome_level_calls": {},
        },
    }
    gate_result = {"g1_features": [], "g1_demote_count": 0, "g1_surface_count": 0,
                   "g2_spotcheck": {"pass": True, "violations": []},
                   "all_joins_symbol_fallback": False}
    md = gm.render_genome_summary_md(genome_map, gate_result, generated="2026-06-28")
    assert "contig length-reconciliation" in md
    assert "length-ambiguous=4" in md  # 2 fasta + 2 bakta
    assert "reconciled=1/3" in md
