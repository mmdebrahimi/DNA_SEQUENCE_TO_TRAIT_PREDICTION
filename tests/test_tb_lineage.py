"""Step 4 — TB lineage-barcode caller + clusters map."""
from __future__ import annotations

import pytest

from dna_decode.data import tb_lineage_barcode as bc
from dna_decode.data.tb_lineage_barcode import BarcodeSNP
from dna_decode.organism_rules import tb_lineage
from dna_decode.organism_rules.tb_vcf import VariantCall

# synthetic barcode: nested lineage4 path + a separate lineage2 marker
BARCODE = [
    BarcodeSNP(100, "lineage4", "A"),
    BarcodeSNP(200, "lineage4.2", "C"),
    BarcodeSNP(300, "lineage4.2.2", "G"),
    BarcodeSNP(400, "lineage2", "T"),
]


def _call(pos, alt):
    return VariantCall(pos=pos, ref="X", alt=alt, gt="1/1")


def test_assign_deepest_supported_sublineage():
    calls = {100: _call(100, "A"), 200: _call(200, "C"), 300: _call(300, "G")}
    assert tb_lineage.assign_lineage(calls, BARCODE) == "lineage4.2.2"


def test_assign_shallower_when_only_top_carried():
    calls = {100: _call(100, "A")}
    assert tb_lineage.assign_lineage(calls, BARCODE) == "lineage4"


def test_wrong_allele_does_not_support():
    calls = {100: _call(100, "G")}  # barcode wants A at 100
    assert tb_lineage.assign_lineage(calls, BARCODE) == tb_lineage.UNASSIGNED


def test_no_barcode_hit_is_unassigned():
    assert tb_lineage.assign_lineage({}, BARCODE) == tb_lineage.UNASSIGNED


def test_clusters_same_lineage_share_id_unassigned_singleton():
    calls_by_strain = {
        "s1": {100: _call(100, "A"), 200: _call(200, "C")},   # lineage4.2
        "s2": {100: _call(100, "A"), 200: _call(200, "C")},   # lineage4.2 (same)
        "s3": {400: _call(400, "T")},                          # lineage2
        "s4": {},                                              # UNASSIGNED
        "s5": {},                                              # UNASSIGNED
    }
    clusters = tb_lineage.lineage_clusters(calls_by_strain, BARCODE)
    assert clusters["s1"] == clusters["s2"]            # same sublineage -> one vote
    assert clusters["s1"] != clusters["s3"]            # different lineage
    assert clusters["s4"] != clusters["s5"]            # UNASSIGNED never merged
    assert len({clusters[s] for s in calls_by_strain}) == 4  # {4.2, 2, s4, s5}


# ---- real pinned barcode (skips if absent) -------------------------------------------------------

_REAL = pytest.mark.skipif(not bc.barcode_available(), reason="pinned barcode not present")


@_REAL
def test_real_barcode_pins_and_loads_1111():
    bc.verify_pin()
    snps = bc.load_barcode()
    assert len(snps) == 1111
    assert all(s.pos > 0 and s.allele and s.lineage for s in snps)


@_REAL
def test_real_barcode_assigns_a_deep_sublineage():
    snps = bc.load_barcode()
    # carry every derived allele of one deep lineage4 path -> assignment is that lineage (most dots present)
    deep = max((s.lineage for s in snps if s.lineage.startswith("lineage4")), key=lambda L: L.count("."))
    calls = {s.pos: _call(s.pos, s.allele) for s in snps if s.lineage == deep}
    assert tb_lineage.assign_lineage(calls, snps) == deep
