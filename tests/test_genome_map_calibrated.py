"""C1-B: per-feature drug labels track the DEPLOYED (calibrated) verdict end-to-end.

These run through the REAL `call_resistance(organism=...)` calibrated-registry path
(no Docker, no genomes — synthetic GFF + main.tsv on disk), proving the verdict-
derived per-feature labels mirror the deployed call even when the calibrated rule
diverges from the default rule. The divergence the C1 fix closes.
"""
from __future__ import annotations

from pathlib import Path

from dna_decode.genome_map import TIER_DETERMINANT_PHENOTYPE
from scripts.genome_map_spike import run_genome_map_for


def _write_pair(tmp_path: Path, locus: str, contig: str, start: int, stop: int,
                symbol: str, cls: str, sub: str, method: str) -> tuple[Path, Path]:
    """Write a 1-CDS GFF + a 1-row AMRFinder main.tsv that protein-id-joins to it."""
    gff = tmp_path / "g.gff3"
    gff.write_text(
        "##gff-version 3\n"
        f"##sequence-region {contig} 1 100000\n"
        f"{contig}\tBakta\tCDS\t{start}\t{stop}\t.\t+\t0\tID=cds-1;locus_tag={locus};product=hypothetical protein\n"
        "##FASTA\n>" + contig + "\nACGT\n",
        encoding="utf-8",
    )
    main = tmp_path / "main.tsv"
    main.write_text(
        "Protein id\tContig id\tStart\tStop\tElement symbol\tElement name\tClass\tSubclass\tMethod\n"
        f"{locus}\t{contig}\t{start}\t{stop}\t{symbol}\tname\t{cls}\t{sub}\t{method}\n",
        encoding="utf-8",
    )
    return gff, main


def test_salmonella_cipro_broad_counts_qnr(tmp_path: Path):
    # Salmonella|ciprofloxacin is CALIBRATED counter=broad/threshold=1 -> a qnr (QUINOLONE
    # class, NOT a QRDR point) IS counted by the deployed call -> the per-feature label
    # must show ciprofloxacin (the DEFAULT qrdr_point rule would have excluded it).
    gff, main = _write_pair(tmp_path, "TAG_Q", "contig_1", 100, 900,
                            "qnrB19", "QUINOLONE", "QUINOLONE", "EXACTX")
    gm = run_genome_map_for("S1", "Salmonella", gff, main, drugs=["ciprofloxacin"])
    f0 = gm["features"][0]
    assert f0["primary_tier"] == TIER_DETERMINANT_PHENOTYPE
    cip = [p for p in f0["phenotype"] if p.get("drug") == "ciprofloxacin"]
    assert cip, "Salmonella broad-counted qnr must be drug-labelled ciprofloxacin"
    assert cip[0]["drug_rule_counted"] is True
    # mirrors the deployed call: genome-level prediction is R (threshold 1 met)
    assert gm["metrics"]["genome_level_calls"]["ciprofloxacin"]["prediction"] == "R"


def test_acinetobacter_meropenem_expression_floor_abstains(tmp_path: Path):
    # Acinetobacter|meropenem is EXPRESSION_FLOOR -> the deployed call ABSTAINs
    # (determinants=[]) -> a carbapenemase must NOT be R-labelled; it surfaces as an
    # explicit ABSTAIN entry (AC8), never a forced meropenem-R call.
    gff, main = _write_pair(tmp_path, "TAG_C", "contig_1", 100, 1300,
                            "blaOXA-23", "CARBAPENEM", "CARBAPENEM", "EXACTX")
    gm = run_genome_map_for("A1", "Acinetobacter_baumannii", gff, main, drugs=["meropenem"])
    f0 = gm["features"][0]
    assert f0["primary_tier"] == TIER_DETERMINANT_PHENOTYPE
    mero = [p for p in f0["phenotype"] if p.get("drug") == "meropenem"]
    assert mero, "an abstaining drug's relevant determinant should carry an ABSTAIN entry (AC8)"
    assert mero[0]["phenotype"] == "ABSTAIN"
    assert mero[0]["drug_rule_counted"] is False  # NOT a forced resistance call
    assert gm["metrics"]["genome_level_calls"]["meropenem"]["prediction"] == "ABSTAIN"


def test_escherichia_default_qrdr_still_excludes_qnr(tmp_path: Path):
    # Control: E. coli (NO registry entry -> DEFAULT qrdr_point cipro rule) -> a qnr is
    # NOT counted -> NOT ciprofloxacin-labelled (DETERMINANT_PRESENT). Proves the per-drug
    # label genuinely tracks the deployed rule, default OR calibrated.
    gff, main = _write_pair(tmp_path, "TAG_Q", "contig_1", 100, 900,
                            "qnrB19", "QUINOLONE", "QUINOLONE", "EXACTX")
    gm = run_genome_map_for("E1", "Escherichia", gff, main, drugs=["ciprofloxacin"])
    f0 = gm["features"][0]
    assert all(p.get("drug") is None for p in f0["phenotype"])  # DETERMINANT_PRESENT, not cipro
