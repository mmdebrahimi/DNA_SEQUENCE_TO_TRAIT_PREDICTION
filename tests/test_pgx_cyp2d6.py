"""Pins the CYP2D6 PGx cell — the last major pharmacogene, SNP surface (chr22).

CYP2D6 is the canonical STRUCTURAL gene (deletions *5, duplications *xN, CYP2D6-CYP2D7 hybrids
*13/*36/*68). A full structural caller is infeasible from a phased SNP VCF, so this cell is a SNP-surface
star-allele caller validated on the SNP-decodable GeT-RM subset — 46/47 core-comparable diplotype
concordance (wiki/pgx_getrm_concordance_cyp2d6_2026-07-06; the single miss is a diagnosed structural
confound). Structural alleles are EXCLUDED (cnv_hybrid_unassessed), NOT withheld.

The LOAD-BEARING tests pin the PRIORITY-ORDERED per-haplotype resolver: CYP2D6 alleles share a SNP
background (*2's 2851/486 is carried by *4/*17/*29/*35/*41; *4 carries *10's 100C>T), so the most-specific
defining SNP must win over the shared background (1846 *4 before 100 *10; every allele-specific SNP before
the 2851 *2 background). The real-data 46/47 number is a committed artifact (VCF is gitignored).
"""
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from dna_decode.pgx import PGX_GENES  # noqa: E402
from dna_decode.pgx import cyp2d6_catalog as c2d6  # noqa: E402
from dna_decode.pgx.cyp2d6_caller import assemble_cyp2d6_diplotype  # noqa: E402
from dna_decode.pgx.runner import call_cyp2d6  # noqa: E402

_HEADER = ("##fileformat=VCFv4.2\n"
           "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")
# component tag -> (pos, ref, alt, rsid): NCBI-verified GRCh38 forward-strand (chr22 minus strand); indels
# use the EXACT 1000G left-anchored representation.
_C = {d.star: (d.pos, d.ref, d.alt, d.rsid) for d in c2d6.COMPONENTS}


def _vcf(tmp_path, rows, name="c2d6.vcf"):
    lines = [_HEADER.rstrip("\n")]
    for tag, gt in rows:
        pos, ref, alt, rsid = _C[tag]
        lines.append(f"chr22\t{pos}\t{rsid}\t{ref}\t{alt}\t.\tPASS\t.\tGT\t{gt}")
    p = tmp_path / name
    p.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return p


def _call(vcf, sample=None):
    return assemble_cyp2d6_diplotype(vcf, c2d6.COMPONENTS, c2d6.STAR_PRIORITY,
                                     reference_allele=c2d6.REFERENCE_ALLELE,
                                     phenotype_fn=c2d6.diplotype_phenotype, gene=c2d6.GENE, sample=sample)


# --- catalog / provenance ---
def test_c2d6_in_pgx_genes():
    assert "cyp2d6" in PGX_GENES


def test_c2d6_coords_grounded_ncbi():
    assert _C["1846A"] == (42128945, "C", "T", "rs3892097")     # *4
    assert _C["100T"] == (42130692, "G", "A", "rs1065852")      # *10 / *4 bkg
    assert _C["2851T"] == (42127941, "G", "A", "rs16947")       # *2 background
    assert _C["486T"] == (42126611, "C", "G", "rs1135840")      # background co-marker
    assert _C["2549delA"] == (42128241, "CT", "C", "rs35742686")  # *3 indel (1000G left-anchored)
    assert _C["1707delT"] == (42129083, "CA", "C", "rs5030655")   # *6 indel
    assert _C["2615del"] == (42128173, "CCTT", "C", "rs5030656")  # *9 indel (K281del)
    assert all(d.chrom == "22" for d in c2d6.COMPONENTS)


@pytest.mark.parametrize("a1,a2,pheno", [
    ("*1", "*1", "Normal Metabolizer"),        # AS 2.0
    ("*1", "*2", "Normal Metabolizer"),        # AS 2.0
    ("*2", "*35", "Normal Metabolizer"),       # AS 2.0
    ("*1", "*4", "Intermediate Metabolizer"),  # AS 1.0
    ("*1", "*10", "Normal Metabolizer"),       # AS 1.25 -> NM
    ("*4", "*10", "Intermediate Metabolizer"), # AS 0.25
    ("*4", "*41", "Intermediate Metabolizer"), # AS 0.25
    ("*4", "*4", "Poor Metabolizer"),          # AS 0.0
    ("*3", "*6", "Poor Metabolizer"),          # AS 0.0
    ("*10", "*10", "Intermediate Metabolizer"),# AS 0.5
    ("*17", "*41", "Intermediate Metabolizer"),# AS 0.75
])
def test_c2d6_activity_score_phenotype(a1, a2, pheno):
    assert c2d6.diplotype_phenotype(a1, a2) == pheno
    assert c2d6.diplotype_phenotype(a2, a1) == pheno


def test_c2d6_unknown_allele_indeterminate():
    assert c2d6.diplotype_phenotype("*4", "*99") == "Indeterminate"


# --- LOAD-BEARING: priority resolution over the shared SNP background ---
def test_c2d6_star4_beats_star10_on_shared_100T(tmp_path):
    # *4 haplotype carries BOTH 1846 (*4) and 100 (*10). Priority: 1846 wins -> *4, not *10.
    # hap0 = {1846A, 100T} (*4); hap1 = {100T} (*10)  -> *4/*10
    r = _call(_vcf(tmp_path, [("1846A", "1|0"), ("100T", "1|1")]))
    assert r.diplotype == "*4/*10"


def test_c2d6_star4_hom(tmp_path):
    r = _call(_vcf(tmp_path, [("1846A", "1|1"), ("100T", "1|1")]))
    assert r.diplotype == "*4/*4" and r.phenotype == "Poor Metabolizer"


def test_c2d6_star17_beats_star2_background(tmp_path):
    # *17 haplotype carries 1023 (*17) + 2851 (*2 bkg) + 486. Priority: 1023 wins -> *17, not *2.
    r = _call(_vcf(tmp_path, [("1023T", "1|0"), ("2851T", "1|0"), ("486T", "1|0")]))
    assert r.diplotype == "*1/*17" and r.phenotype == "Normal Metabolizer"


def test_c2d6_star35_beats_star2_background(tmp_path):
    r = _call(_vcf(tmp_path, [("31A", "1|0"), ("2851T", "1|0"), ("486T", "1|0")]))
    assert r.diplotype == "*1/*35"


def test_c2d6_star2_background_alone(tmp_path):
    # 2851 + 486 with NO allele-specific SNP -> *2
    r = _call(_vcf(tmp_path, [("2851T", "1|0"), ("486T", "1|0")]))
    assert r.diplotype == "*1/*2"


def test_c2d6_star10_alone_no_1846(tmp_path):
    # 100 WITHOUT 1846 -> *10 (not *4)
    r = _call(_vcf(tmp_path, [("100T", "1|0"), ("486T", "1|0")]))
    assert r.diplotype == "*1/*10"


def test_c2d6_486T_alone_is_reference(tmp_path):
    # the background co-marker alone (no 2851, no specific) is NOT a star -> *1/*1
    r = _call(_vcf(tmp_path, [("486T", "1|1")]))
    assert r.diplotype == "*1/*1"


def test_c2d6_indel_star3(tmp_path):
    r = _call(_vcf(tmp_path, [("2549delA", "1|0")]))
    assert r.diplotype == "*1/*3" and r.phenotype == "Intermediate Metabolizer"


def test_c2d6_indel_star9_het(tmp_path):
    r = _call(_vcf(tmp_path, [("2615del", "0|1")]))
    assert r.diplotype == "*1/*9"


def test_c2d6_two_specific_on_one_haplotype_flagged(tmp_path):
    # a data anomaly: two allele-specific SNPs in cis -> priority still resolves, but it is FLAGGED
    r = _call(_vcf(tmp_path, [("1846A", "1|0"), ("1023T", "1|0")]))
    assert any(f.startswith("multi_specific_haplotype") for f in r.flags)


def test_c2d6_all_ref(tmp_path):
    r = _call(_vcf(tmp_path, [("1846A", "0|0"), ("2851T", "0|0")]))
    assert r.diplotype == "*1/*1" and r.phenotype == "Normal Metabolizer"


# --- runner record: the load-bearing cnv_hybrid_unassessed honesty field ---
def test_c2d6_runner_record(tmp_path):
    rec = call_cyp2d6(_vcf(tmp_path, [("1846A", "1|0"), ("100T", "1|1")]), sample_id="T1")
    assert rec["gene"] == "CYP2D6"
    assert rec["diplotype"] == "*4/*10"
    assert rec["activity_score"] == 0.25
    assert rec["phenotype_abbrev"] == "IM"
    assert rec["cnv_hybrid_unassessed"] is True          # structural surface honesty
    assert rec["caller"]["cnv_hybrid_unassessed"] is True
    assert rec["caller"]["is_priority_resolver"] is True
    assert rec["caller"]["phenotype_is_faithful_to_cpic"] is True


def test_c2d6_cli_routing(tmp_path, capsys):
    from dna_decode.pgx.cli import main
    rc = main([str(_vcf(tmp_path, [("1846A", "1|1")])), "--gene", "cyp2d6", "--json-only"])
    assert rc == 0
    out = capsys.readouterr().out
    assert "CYP2D6" in out and "cnv_hybrid_unassessed" in out


# --- honest tiered truth classification (structural/ambiguous EXCLUDED; never scored as match/miss) ---
def test_c2d6_truth_classifier_tiers():
    from scripts.pgx_getrm_concordance import _classify_cyp2d6_truth as clf
    assert clf("*1/*4")[0] == "core_snp"
    assert clf("*2/*3")[0] == "core_snp"
    assert clf("*4/*5")[0] == "structural"        # *5 gene deletion
    assert clf("*2x2/*71")[0] == "structural"     # *xN duplication
    assert clf("*1/(*68)+*4")[0] == "structural"  # *68 hybrid (inside a paren+compound)
    assert clf("*2 (*35)/*9")[0] == "ambiguous"   # parenthetical alternative annotation
    assert clf("*1/*21")[0] == "noncore_snp"      # non-core SNP allele
    assert clf("*1/*14")[0] == "noncore_snp"


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
