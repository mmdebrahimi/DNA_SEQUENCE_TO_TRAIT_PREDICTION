"""Shared sentinel-withhold layer (Wave 1.5 refactor) -- proves the brainstorm's grounded findings are fixed:

- Issue 1: the COMPOUND caller (TPMT's path) now applies the non-core sentinel withhold (was a no-op).
- Issue 2: the CYP2C19 *35-shared-with-*2 rule is generalized to `SentinelVariant.accounted_by_core` DATA,
           not hardcoded in the generic caller.
- Issue 3: exact-ALT matching -> a benign DIFFERENT ALT at a sentinel coordinate does NOT false-withhold.
- Issue 4: two sentinels at the SAME (chrom,pos) with different ALTs are both counted (no overwrite).

Offline, no network/Docker. Runnable via pytest OR standalone.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from dna_decode.pgx.caller import (  # noqa: E402
    DiplotypeResult, VariantCall, apply_sentinel_withhold, _scan_sentinel_counts)
from dna_decode.pgx.cyp2c19_catalog import SentinelVariant  # noqa: E402
from dna_decode.pgx.compound_caller import assemble_compound_diplotype  # noqa: E402
from dna_decode.pgx import tpmt_catalog as tp  # noqa: E402

_HDR = ("##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")


def _vcf(tmp_path, rows, name="s.vcf"):
    """rows: (chrom, pos, rsid, ref, alt, gt)."""
    lines = [_HDR.rstrip("\n")]
    for chrom, pos, rsid, ref, alt, gt in rows:
        lines.append(f"{chrom}\t{pos}\t{rsid}\t{ref}\t{alt}\t.\tPASS\t.\tGT\t{gt}")
    p = tmp_path / name
    p.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return p


def _res():
    return DiplotypeResult("ok", "*1/*1", "*1", "*1", "Normal Metabolizer", "unphased",
                           core_proxy_diplotype="*1/*1")


# --- Issue 2: accounted_by_core generalizes the *35/*2 rule (no hardcode) --------------------------------

def test_accounted_by_core_subtracts_core_signal():
    s = SentinelVariant("rsX", "10", 1, "A", "G", "*35", "shared with *2", accounted_by_core="*2")
    star2 = VariantCall("*2", "rs2", 1, True, "0/1", False, 1, False, None)  # one *2 copy
    # n == core copies -> fully accounted -> NOT a non-core hit
    r = apply_sentinel_withhold(_res(), {"*2": star2}, {"rsX": 1}, [s])
    assert r.phenotype_status == "ok" and r.phenotype == "Normal Metabolizer"
    # n IN EXCESS of core copies -> a genuine *35 signal -> withhold
    r2 = apply_sentinel_withhold(_res(), {"*2": star2}, {"rsX": 2}, [s])
    assert r2.phenotype_status == "phenotype_withheld" and r2.phenotype is None


def test_no_accounted_by_core_withholds_on_any_hit():
    s = SentinelVariant("rsY", "6", 2, "A", "G", "*8", "plain non-core")  # accounted_by_core defaults None
    r = apply_sentinel_withhold(_res(), {}, {"rsY": 1}, [s])
    assert r.phenotype_status == "phenotype_withheld" and "non_core_allele_sentinel" in r.flags


def test_cyp2c19_star35_sentinel_carries_the_data():
    from dna_decode.pgx.cyp2c19_catalog import SENTINELS
    star35 = [s for s in SENTINELS if s.implies == "*35"][0]
    assert star35.accounted_by_core == "*2"   # the rule is now DATA, not a hardcode in caller.py


# --- Issue 4: same-coordinate different-ALT sentinels are both counted -----------------------------------

def test_scan_same_pos_two_alts_no_overwrite(tmp_path):
    a = SentinelVariant("rsA", "6", 100, "C", "T", "*A", "alt T")
    b = SentinelVariant("rsB", "6", 100, "C", "G", "*B", "alt G")
    v = _vcf(tmp_path, [("6", 100, "multi", "C", "T,G", "1/2")])   # one T copy + one G copy
    counts = _scan_sentinel_counts(v, sentinels=[a, b])
    assert counts["rsA"] == 1 and counts["rsB"] == 1   # neither overwrote the other


# --- Issue 1: the COMPOUND caller now withholds (was a no-op) --------------------------------------------

def _synth_tpmt_sentinel():
    # a synthetic TPMT non-core sentinel on chr6 (coordinate distinct from the two component SNPs)
    return SentinelVariant("rs_synth8", "6", 18140000, "A", "G", "*8", "synthetic *8 non-core test site")


def test_compound_path_withholds_with_sentinel(tmp_path):
    s = _synth_tpmt_sentinel()
    # component sites present as reference (0/0) -> a clean *1/*1 compound call ...
    rows = [("6", 18138997, "rs1800460", "C", "T", "0/0"),
            ("6", 18130687, "rs1142345", "T", "C", "0/0"),
            ("6", 18140000, "rs_synth8", "A", "G", "0/1")]   # ... + the non-core sentinel present
    v = _vcf(tmp_path, rows)
    res = assemble_compound_diplotype(v, tp.COMPONENTS, tp.COMPOUND_RULES,
                                      reference_allele=tp.REFERENCE_ALLELE,
                                      phenotype_fn=tp.diplotype_phenotype, gene=tp.GENE, sentinels=[s])
    assert res.phenotype_status == "phenotype_withheld" and res.phenotype is None
    assert any(h["implies"] == "*8" for h in res.sentinel_hits)


def test_compound_path_empty_sentinels_is_noop(tmp_path):
    rows = [("6", 18138997, "rs1800460", "C", "T", "0/0"),
            ("6", 18130687, "rs1142345", "T", "C", "0/0")]
    v = _vcf(tmp_path, rows)
    res = assemble_compound_diplotype(v, tp.COMPONENTS, tp.COMPOUND_RULES,
                                      reference_allele=tp.REFERENCE_ALLELE,
                                      phenotype_fn=tp.diplotype_phenotype, gene=tp.GENE, sentinels=[])
    assert res.phenotype_status == "ok" and res.phenotype is not None   # unchanged behavior


# --- Issue 3: exact ALT -> a benign different ALT at the coord does NOT false-withhold -------------------

def test_exact_alt_no_false_withhold(tmp_path):
    s = _synth_tpmt_sentinel()   # expects ALT "G" at 6:18140000
    rows = [("6", 18138997, "rs1800460", "C", "T", "0/0"),
            ("6", 18130687, "rs1142345", "T", "C", "0/0"),
            ("6", 18140000, "rs_benign", "A", "C", "0/1")]   # a DIFFERENT (benign) ALT at the same coord
    v = _vcf(tmp_path, rows)
    res = assemble_compound_diplotype(v, tp.COMPONENTS, tp.COMPOUND_RULES,
                                      reference_allele=tp.REFERENCE_ALLELE,
                                      phenotype_fn=tp.diplotype_phenotype, gene=tp.GENE, sentinels=[s])
    assert res.phenotype_status == "ok" and res.phenotype is not None   # exact-ALT: NOT withheld


def test_tpmt_populated_sentinels_withhold_real_allele(tmp_path):
    """Uses the REAL populated tpmt_catalog.SENTINELS (10 rows, Ensembl-verified) via the runner path."""
    from dna_decode.pgx.runner import call_tpmt
    assert len(tp.SENTINELS) == 10 and all(s.accounted_by_core is None for s in tp.SENTINELS)
    # a real TPMT*8 carrier: rs56161402 (6:18130762 C>T) het, component sites reference -> WITHHELD not *1
    rows = [("6", 18138997, "rs1800460", "C", "T", "0/0"),   # *3B component = ref
            ("6", 18130687, "rs1142345", "T", "C", "0/0"),   # *3C component = ref
            ("6", 18130762, "rs56161402", "C", "T", "0/1")]  # *8 non-core sentinel present
    v = _vcf(tmp_path, rows)
    rec = call_tpmt(v)
    assert rec["phenotype_status"] == "phenotype_withheld" and rec["phenotype"] is None


def test_cyp2b6_sentinels_withhold_and_no_false_withhold(tmp_path):
    """Uses REAL populated cyp2b6_catalog.SENTINELS (3 distinctive-SNP rows) via the runner."""
    from dna_decode.pgx.runner import call_cyp2b6
    from dna_decode.pgx import cyp2b6_catalog as c6
    assert len(c6.SENTINELS) == 3 and all(s.accounted_by_core is None for s in c6.SENTINELS)
    # a real *2 carrier (rs8192709 19:40991369 C>T) -> WITHHELD
    v = _vcf(tmp_path, [("19", 41006936, "rs3745274", "G", "T", "0/0"),      # core *6 site = ref
                        ("19", 40991369, "rs8192709", "C", "T", "0/1")], "c6a.vcf")  # *2 sentinel present
    assert call_cyp2b6(v)["phenotype_status"] == "phenotype_withheld"
    # a real core *6 carrier (516T at rs3745274) with NO distinctive non-core SNP -> NOT withheld (trap check)
    v2 = _vcf(tmp_path, [("19", 41006936, "rs3745274", "G", "T", "1/1")], "c6b.vcf")
    rec = call_cyp2b6(v2)
    assert rec["phenotype_status"] != "phenotype_withheld" and rec["phenotype"] is not None


def test_cyp2c8_sentinels_withhold(tmp_path):
    """Uses REAL populated cyp2c8_catalog.SENTINELS (4 PharmVar-sourced rows) via the runner."""
    from dna_decode.pgx.runner import call_cyp2c8
    from dna_decode.pgx import cyp2c8_catalog as c8
    assert len(c8.SENTINELS) == 4 and all(s.accounted_by_core is None for s in c8.SENTINELS)
    # a real *15 carrier (rs41286886 10:95064901 C>T) -> WITHHELD
    v = _vcf(tmp_path, [("10", 95058349, "rs11572103", "T", "A", "0/0"),   # *2 core = ref
                        ("10", 95064901, "rs41286886", "C", "T", "0/1")], "c8.vcf")  # *15 sentinel present
    assert call_cyp2c8(v)["phenotype_status"] == "phenotype_withheld"


def test_cyp3a5_sentinels_populated_and_fire_synthetic(tmp_path):
    """CYP3A5 is UNDERPOWERED on the GeT-RM cohort (0 non-core carriers); pin that the populated sentinels
    are correct + fire on a synthetic *8 carrier (the real-cohort withhold is validation-deferred)."""
    from dna_decode.pgx.runner import call_cyp3a5
    from dna_decode.pgx import cyp3a5_catalog as c3
    assert len(c3.SENTINELS) == 4 and all(s.accounted_by_core is None for s in c3.SENTINELS)
    v = _vcf(tmp_path, [("7", 99672916, "rs776746", "T", "C", "0/0"),      # *3 core = ref
                        ("7", 99676198, "rs55817950", "G", "A", "0/1")], "c3.vcf")  # *8 sentinel present
    assert call_cyp3a5(v)["phenotype_status"] == "phenotype_withheld"


if __name__ == "__main__":
    import tempfile
    fns = [(k, v) for k, v in sorted(globals().items()) if k.startswith("test_") and callable(v)]
    for k, fn in fns:
        if "tmp_path" in fn.__code__.co_varnames[:fn.__code__.co_argcount]:
            with tempfile.TemporaryDirectory() as d:
                fn(Path(d))
        else:
            fn()
        print(f"PASS {k}")
    print(f"\n{len(fns)} passed")
