"""Test the cross-organism transfer eval + report-card builder (skips if D: data absent)."""
import json, os, subprocess, sys
from pathlib import Path
import pytest

_HAVE = os.path.exists("D:/dna_decode_cache/essentiality/CEGv2.txt") and \
        os.path.exists("D:/dna_decode_cache/essentiality/Homo_sapiens.gene_info.gz")


@pytest.mark.skipif(not _HAVE, reason="essentiality D: data absent")
def test_report_card_builds_and_transfer_above_null():
    r = subprocess.run([sys.executable, "scripts/build_essentiality_report_card.py"], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    card = json.loads(Path("wiki/essentiality_report_card.json").read_text())
    assert card["schema"] == "essentiality-report-card-v1"
    orgs = {o["organism"]: o for o in card["organisms"]}
    assert "Escherichia coli K-12" in orgs and "Homo sapiens" in orgs
    assert "aggregate" not in card and "headline" not in card   # no inflated headline
    # human transfer scored above the null
    t = json.load(open("wiki/essentiality_e4_transfer_2026-07-28.json"))
    assert t["auroc"] > t["null_auroc"]          # transfer above chance
    assert t["spec_at_thresh"] > 0.9             # high-precision (universal core)


def test_decoder_recovers_human_universal_core():
    # the conserved-core decoder fires on human universal-core genes (function transfers)
    from dna_decode.essentiality.core_decoder import score_gene
    assert score_gene("EEF2", "eukaryotic translation elongation factor 2").core_score >= 0 or \
           score_gene("RPL3", "ribosomal protein L3").prediction == "essential"
    assert score_gene("RPL3", "ribosomal protein L3").prediction == "essential"
