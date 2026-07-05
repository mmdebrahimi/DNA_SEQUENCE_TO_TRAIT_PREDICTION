"""Pins the CYP2C8 PGx cell -- star-allele CALLING only (NO CPIC metabolizer phenotype; substrate-dependent).

Third CYP2C-cluster gene (joins C9 + C19 on chr10). Independently validatable vs the GeT-RM consensus
(CYP2C8_getrm_ngs). These tests pin: the Ensembl-verified GRCh38 coordinates, the calling logic on
synthetic VCFs, and the load-bearing HONESTY invariant (has_cpic_phenotype is False; no PM/IM/NM ever
emitted for CYP2C8).
"""
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from dna_decode.pgx import PGX_GENES  # noqa: E402
from dna_decode.pgx import cyp2c8_catalog as c8  # noqa: E402
from dna_decode.pgx.caller import call_diplotype  # noqa: E402
from dna_decode.pgx.runner import call_cyp2c8  # noqa: E402

_HEADER = ("##fileformat=VCFv4.2\n"
           "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")
# Ensembl-verified GRCh38 forward-strand coords (CYP2C8 minus strand; coding change is the revcomp).
_C8 = {
    "*2": (95058349, "T", "A", "rs11572103"),
    "*3": (95067273, "C", "T", "rs11572080"),
    "*4": (95058362, "G", "C", "rs1058930"),
}


def _c8_vcf(tmp_path, rows, name="c8.vcf"):
    lines = [_HEADER.rstrip("\n")]
    for star, gt in rows:
        pos, ref, alt, rsid = _C8[star]
        lines.append(f"chr10\t{pos}\t{rsid}\t{ref}\t{alt}\t.\tPASS\t.\tGT\t{gt}")
    p = tmp_path / name
    p.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return p


def _call_c8(vcf, sample=None):
    return call_diplotype(vcf, sample=sample, defining=c8.CORE_DEFINING, sentinels=c8.SENTINELS,
                          reference_allele=c8.REFERENCE_ALLELE, phenotype_fn=c8.diplotype_phenotype,
                          gene=c8.GENE)


# --- catalog / provenance ---
def test_c8_in_pgx_genes():
    assert "cyp2c8" in PGX_GENES


def test_c8_coords_grounded_ensembl():
    """GRCh38 forward-strand coords VERIFIED via Ensembl REST; a drift here is a real provenance break."""
    by = {d.star: d for d in c8.CORE_DEFINING}
    assert (by["*2"].pos, by["*2"].ref, by["*2"].alt, by["*2"].rsid) == (95058349, "T", "A", "rs11572103")
    assert (by["*3"].pos, by["*3"].ref, by["*3"].alt, by["*3"].rsid) == (95067273, "C", "T", "rs11572080")
    assert (by["*4"].pos, by["*4"].ref, by["*4"].alt, by["*4"].rsid) == (95058362, "G", "C", "rs1058930")
    assert all(d.chrom == "10" for d in c8.CORE_DEFINING)


def test_c8_no_cpic_phenotype_invariant():
    """LOAD-BEARING honesty: CYP2C8 function is substrate-dependent -> NO CPIC metabolizer phenotype."""
    assert c8.HAS_CPIC_PHENOTYPE is False
    # the function seam never returns a PM/IM/NM metabolizer label
    for a1 in ("*1", "*2", "*3", "*4"):
        for a2 in ("*1", "*2", "*3", "*4"):
            out = c8.diplotype_function(a1, a2)
            assert "Metabolizer" not in out
            assert out  # non-empty descriptor


def test_c8_unknown_allele_indeterminate():
    assert c8.diplotype_function("*2", "*99") == "Indeterminate"


# --- caller (VCF) ---
def test_c8_star2_het(tmp_path):
    r = _call_c8(_c8_vcf(tmp_path, [("*2", "0/1")]))
    assert r.diplotype == "*1/*2"


def test_c8_star3_hom(tmp_path):
    r = _call_c8(_c8_vcf(tmp_path, [("*3", "1/1")]))
    assert r.diplotype == "*3/*3"


def test_c8_star4_het(tmp_path):
    r = _call_c8(_c8_vcf(tmp_path, [("*4", "0/1")]))
    assert r.diplotype == "*1/*4"


def test_c8_all_ref_star1(tmp_path):
    r = _call_c8(_c8_vcf(tmp_path, [("*2", "0/0"), ("*3", "0/0"), ("*4", "0/0")]))
    assert r.diplotype == "*1/*1"


def test_c8_compound_het_unphased(tmp_path):
    # *2 het + *4 het, unphased -> trans assumption -> *2/*4
    r = _call_c8(_c8_vcf(tmp_path, [("*2", "0/1"), ("*4", "0/1")]))
    assert set(r.diplotype.split("/")) == {"*2", "*4"}


# --- runner record (honest shape) ---
def test_c8_runner_record_shape(tmp_path):
    rec = call_cyp2c8(_c8_vcf(tmp_path, [("*2", "0/1")]), sample_id="C8_1")
    assert rec["gene"] == "CYP2C8"
    assert rec["trait"] == "pgx_star_allele_diplotype"
    assert rec["diplotype"] == "*1/*2"
    assert rec["has_cpic_phenotype"] is False
    assert rec["phenotype"] is None                       # never a CPIC phenotype
    assert rec["phenotype_abbrev"] is None
    assert rec["function_annotation"]                     # substrate-dependent descriptor present
    assert rec["caller"]["calling_independently_validatable"] is True
    assert rec["caller"]["phenotype_is_faithful_to_cpic"] is False


def test_c8_cli_routing(tmp_path, capsys):
    from dna_decode.pgx.cli import main
    rc = main([str(_c8_vcf(tmp_path, [("*3", "1/1")])), "--gene", "cyp2c8", "--json-only"])
    assert rc == 0                                        # a valid *3/*3 call -> exit 0
    out = capsys.readouterr().out
    assert "CYP2C8" in out
    assert "Metabolizer" not in out                       # honesty: no phenotype rendered


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
