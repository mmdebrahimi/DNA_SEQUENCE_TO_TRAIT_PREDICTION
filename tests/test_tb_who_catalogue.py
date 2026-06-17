"""Step 2 — WHO catalogue loader + pin-verify + RIF/INH grade-1/2 join."""
from __future__ import annotations

import hashlib

import pytest

from dna_decode.data import tb_who_catalogue as cat

# ---- synthetic fixture (always runs; no gitignored data needed) ----------------------------------

_MASTER_COLS = ["drug", "gene", "mutation", "variant", "tier", cat.GRADE_COL]
_MASTER_ROWS = [
    ["Rifampicin", "rpoB", "p.Ser450Leu", "rpoB_p.Ser450Leu", "1", "1) Assoc w R"],
    ["Isoniazid", "katG", "p.Ser315Thr", "katG_p.Ser315Thr", "1", "1) Assoc w R"],
    ["Isoniazid", "inhA", "c.-15C>T", "inhA_c.-15C>T", "2", "2) Assoc w R - Interim"],
    ["Isoniazid", "ahpC", "p.Foo", "ahpC_x", "2", "3) Uncertain significance"],  # NOT grade-1/2
]
_COORDS_ROWS = [
    ["rpoB_p.Ser450Leu", cat.CHROM, "761155", "C", "T"],
    ["katG_p.Ser315Thr", cat.CHROM, "2155168", "C", "G"],
    ["inhA_c.-15C>T", cat.CHROM, "1673425", "C", "T"],
]


def _tsv(rows):
    return "\n".join("\t".join(r) for r in rows) + "\n"


@pytest.fixture
def cat_dir(tmp_path):
    m = _tsv([_MASTER_COLS] + _MASTER_ROWS)
    c = _tsv([["variant", "chromosome", "position", "reference_nucleotide",
               "alternative_nucleotide"]] + _COORDS_ROWS)
    (tmp_path / cat.MASTER_FILE).write_text(m, encoding="utf-8")
    (tmp_path / cat.COORDS_FILE).write_text(c, encoding="utf-8")
    lines = ["# pinned_commit: deadbeef"]
    for fname in (cat.MASTER_FILE, cat.COORDS_FILE):
        raw = (tmp_path / fname).read_bytes()  # hash on-disk bytes (Windows may CRLF-translate)
        lines.append(f"{hashlib.sha256(raw).hexdigest()}  {len(raw)}  {fname}")
    (tmp_path / cat.CHECKSUMS_FILE).write_text("\n".join(lines) + "\n", encoding="utf-8")
    return tmp_path


def test_verify_pins_passes_on_matching_files(cat_dir):
    assert cat.verify_pins(cat_dir) == {cat.MASTER_FILE: True, cat.COORDS_FILE: True}


def test_verify_pins_raises_on_drift(cat_dir):
    (cat_dir / cat.MASTER_FILE).write_text("tampered\n", encoding="utf-8")
    with pytest.raises(cat.CataloguePinError):
        cat.verify_pins(cat_dir)


def test_join_rif_resolves_coords(cat_dir):
    dets = cat.load_determinants("rifampicin", cat_dir)
    assert len(dets) == 1
    d = dets[0]
    assert (d.gene, d.variant, d.pos, d.ref, d.alt) == ("rpoB", "rpoB_p.Ser450Leu", 761155, "C", "T")
    assert cat.is_grade_1_2(d.grade)


def test_inh_scope_is_all_grade12_loci_not_just_katg(cat_dir):
    # Ratified A: all grade-1/2 INH loci. katG + inhA included; the grade-3 ahpC row excluded.
    dets = cat.load_determinants("isoniazid", cat_dir)
    genes = {d.gene for d in dets}
    assert genes == {"katG", "inhA"}
    assert all(cat.is_grade_1_2(d.grade) for d in dets)


def test_grade_1_2_count_excludes_uncertain(cat_dir):
    assert cat.grade_1_2_count(cat_dir) == 3  # rpoB + katG + inhA; ahpC (grade 3) excluded


# ---- real pinned catalogue (skips when the gitignored files are absent) --------------------------

_REAL = pytest.mark.skipif(not cat.catalogue_available(),
                           reason="pinned WHO catalogue not present (gitignored)")


@_REAL
def test_real_catalogue_pins_verify():
    cat.verify_pins()


@_REAL
def test_real_rif_includes_s450l_aligned():
    dets = cat.load_determinants("rifampicin")
    s450l = [d for d in dets if d.variant == "rpoB_p.Ser450Leu"]
    assert s450l and any(d.pos == 761155 and d.ref == "C" and d.alt == "T" for d in s450l)
    assert all(d.gene == "rpoB" for d in dets)  # RIF grade-1/2 is rpoB


@_REAL
def test_real_inh_spans_multiple_loci_and_has_s315t():
    dets = cat.load_determinants("isoniazid")
    genes = {d.gene for d in dets}
    assert "katG" in genes and len(genes) > 1  # ratified A — not katG-only
    assert any(d.variant == "katG_p.Ser315Thr" and d.pos == 2155168 for d in dets)


@_REAL
def test_real_grade_1_2_count_is_438():
    assert cat.grade_1_2_count() == 438
