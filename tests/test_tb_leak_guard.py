"""Step 6 — frozen-surface leak guard for the TB cell.

HONEST SCOPE: this proves (a) the TB rule is NOT fitted on CRyPTIC phenotype labels, (b) the rule sources
are checksum/version-pinned, and (c) the frozen E. coli surface is byte-untouched. It CANNOT prove
biological independence — a CRyPTIC-scored number stays an in-distribution KNOWLEDGE_BASELINE; real
independence is the separate post-2023 gold-set arm (Step 7).
"""
from __future__ import annotations

import hashlib
from pathlib import Path

import pytest

from dna_decode.data import tb_lineage_barcode, tb_who_catalogue

ROOT = Path(__file__).resolve().parents[1]

# Frozen surface, byte-frozen at commit b3761c8 (2026-06-13). The TB work must never edit these.
FROZEN_SHA = {
    "dna_decode/eval/amr_rules.py":
        "a983bf28e4ff4f89034b152404e49ec3aa6b3907ac81c670b14e60e7cfe1fad4",
    "dna_decode/data/calibrated_amr_rules.json":
        "d442b768979f35caa9b33875bf289a21a736fec68f5d85f1483578c7cab352ad",
}

# TB rule-construction modules — must NOT read CRyPTIC phenotype columns.
TB_RULE_MODULES = [
    "dna_decode/organism_rules/tb_vcf.py",
    "dna_decode/organism_rules/tb_amr.py",
    "dna_decode/organism_rules/tb_lineage.py",
    "dna_decode/data/tb_who_catalogue.py",
    "dna_decode/data/tb_lineage_barcode.py",
]
FORBIDDEN_PHENOTYPE_TOKENS = ["BINARY_PHENOTYPE", "PHENOTYPE_QUALITY", "_MIC"]


@pytest.mark.parametrize("rel,want", FROZEN_SHA.items())
def test_frozen_surface_byte_untouched(rel, want):
    got = hashlib.sha256((ROOT / rel).read_bytes()).hexdigest()
    assert got == want, f"FROZEN SURFACE EDITED: {rel} (TB work must never touch it)"


@pytest.mark.parametrize("rel", TB_RULE_MODULES)
def test_tb_rule_modules_do_not_read_cryptic_phenotype(rel):
    src = (ROOT / rel).read_text(encoding="utf-8")
    hits = [t for t in FORBIDDEN_PHENOTYPE_TOKENS if t in src]
    assert not hits, f"{rel} references CRyPTIC phenotype columns {hits} during rule construction"


def test_who_catalogue_is_checksum_pinned():
    assert (tb_who_catalogue.CAT_DIR / tb_who_catalogue.CHECKSUMS_FILE).exists()
    if tb_who_catalogue.catalogue_available():
        tb_who_catalogue.verify_pins()  # raises on drift


def test_barcode_is_checksum_pinned():
    assert (tb_lineage_barcode.BARCODE_DIR / tb_lineage_barcode.CHECKSUMS_FILE).exists()
    if tb_lineage_barcode.barcode_available():
        tb_lineage_barcode.verify_pin()  # raises on drift
