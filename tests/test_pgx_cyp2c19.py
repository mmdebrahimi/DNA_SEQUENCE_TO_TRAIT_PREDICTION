"""Pins the human-PGx CYP2C19 cell (non-frozen). Pure-stdlib VCF -> diplotype -> CPIC phenotype.

The diplotype->phenotype expectations are the CPIC standardized consensus (Caudle 2020) -- the same
table the GeT-RM consensus diplotypes are interpreted against -- so reproducing them from a synthetic VCF
encoding a real defining genotype validates the deterministic CALLING + ASSIGNMENT logic end to end. The
full live GeT-RM (Coriell) x 1000 Genomes cohort-concordance number is the P3 follow-up (needs the VCF
fetch); these tests pin the logic that run depends on.
"""
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from dna_decode.pgx.caller import call_diplotype  # noqa: E402
from dna_decode.pgx.cyp2c19_catalog import (  # noqa: E402
    ALLELE_FUNCTION,
    CORE_DEFINING,
    diplotype_phenotype,
)

_HEADER = (
    "##fileformat=VCFv4.2\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
)
# GRCh38 chr10 positions (grounded) for the three core defining SNPs + the two v0.1 sentinels.
_POS = {"*2": (94781859, "G", "A", "rs4244285"),
        "*3": (94780653, "G", "A", "rs4986893"),
        "*17": (94761900, "C", "T", "rs12248560"),
        "rs28399504": (94762706, "A", "G", "rs28399504"),   # *4 sentinel
        "rs12769205": (94775367, "A", "G", "rs12769205")}   # *35 sentinel


def _vcf(tmp_path, rows, name="s.vcf"):
    """rows: list of (key, GT). key is a core star (*2/*3/*17) or a sentinel rsid. One record per row."""
    lines = [_HEADER.rstrip("\n")]
    for key, gt in rows:
        pos, ref, alt, rsid = _POS[key]
        lines.append(f"chr10\t{pos}\t{rsid}\t{ref}\t{alt}\t.\tPASS\t.\tGT\t{gt}")
    p = tmp_path / name
    p.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return p


# --- catalog (pure) ---------------------------------------------------------

def test_allele_functions_grounded():
    assert ALLELE_FUNCTION == {"*1": "normal", "*2": "none", "*3": "none", "*17": "increased"}


def test_defining_coordinates_grounded():
    by = {d.star: d for d in CORE_DEFINING}
    assert (by["*2"].pos, by["*2"].rsid) == (94781859, "rs4244285")
    assert (by["*3"].pos, by["*3"].rsid) == (94780653, "rs4986893")
    assert (by["*17"].pos, by["*17"].rsid) == (94761900, "rs12248560")


@pytest.mark.parametrize("a1,a2,expected", [
    ("*1", "*1", "Normal Metabolizer"),
    ("*1", "*17", "Rapid Metabolizer"),
    ("*17", "*17", "Ultrarapid Metabolizer"),
    ("*1", "*2", "Intermediate Metabolizer"),
    ("*2", "*17", "Intermediate Metabolizer"),   # provisional IM per CPIC
    ("*2", "*2", "Poor Metabolizer"),
    ("*2", "*3", "Poor Metabolizer"),
])
def test_cpic_diplotype_phenotype(a1, a2, expected):
    assert diplotype_phenotype(a1, a2) == expected
    assert diplotype_phenotype(a2, a1) == expected  # order-independent


def test_noncore_allele_indeterminate():
    assert diplotype_phenotype("*2", "*35") == "Indeterminate"


# --- caller (VCF end to end) ------------------------------------------------

def test_getrm_consensus_star2_homozygous_poor_metabolizer(tmp_path):
    """A *2/*2 sample (rs4244285 hom-alt) -- a canonical GeT-RM-consensus genotype -> Poor Metabolizer."""
    vcf = _vcf(tmp_path, [("*2", "1/1")])
    r = call_diplotype(vcf)
    assert r.status == "ok"
    assert r.diplotype == "*2/*2"
    assert r.phenotype == "Poor Metabolizer"


def test_star1_star2_het_intermediate(tmp_path):
    vcf = _vcf(tmp_path, [("*2", "0/1")])
    r = call_diplotype(vcf)
    assert r.diplotype == "*1/*2"
    assert r.phenotype == "Intermediate Metabolizer"


def test_all_reference_is_star1_star1_normal(tmp_path):
    """Defining sites present but homozygous-reference -> *1/*1 Normal Metabolizer."""
    vcf = _vcf(tmp_path, [("*2", "0/0"), ("*3", "0/0"), ("*17", "0/0")])
    r = call_diplotype(vcf)
    assert r.diplotype == "*1/*1"
    assert r.phenotype == "Normal Metabolizer"


def test_star1_star17_rapid(tmp_path):
    vcf = _vcf(tmp_path, [("*17", "0|1")])
    r = call_diplotype(vcf)
    assert r.diplotype == "*1/*17"
    assert r.phenotype == "Rapid Metabolizer"


def test_phased_two_hets_trans_resolved(tmp_path):
    """Phased *2 on hap0, *17 on hap1 -> exact *2/*17 (IM), phasing='phased'."""
    vcf = _vcf(tmp_path, [("*2", "1|0"), ("*17", "0|1")])
    r = call_diplotype(vcf)
    assert r.phasing == "phased"
    assert set((r.allele1, r.allele2)) == {"*2", "*17"}
    assert r.phenotype == "Intermediate Metabolizer"


def test_unphased_two_hets_flags_trans_assumption(tmp_path):
    vcf = _vcf(tmp_path, [("*2", "0/1"), ("*17", "0/1")])
    r = call_diplotype(vcf)
    assert r.phasing == "unphased"
    assert "unphased_trans_assumption" in r.flags
    assert set((r.allele1, r.allele2)) == {"*2", "*17"}


def test_absent_sites_flagged_not_silent(tmp_path):
    """Only *2 record present (het); *3/*17 absent -> *1/*2 with an explicit assumed-reference flag."""
    vcf = _vcf(tmp_path, [("*2", "0/1")])
    r = call_diplotype(vcf)
    assert r.diplotype == "*1/*2"
    assert any(f.startswith("assumed_reference_at_uncalled_sites") for f in r.flags)


def test_no_call_flagged(tmp_path):
    vcf = _vcf(tmp_path, [("*2", "./."), ("*3", "0/0"), ("*17", "0/0")])
    r = call_diplotype(vcf)
    assert any(f.startswith("no_call_at") for f in r.flags)


def test_no_defining_site_is_no_input(tmp_path):
    p = tmp_path / "empty.vcf"
    p.write_text(_HEADER + "chr10\t11111\trsX\tA\tG\t.\tPASS\t.\tGT\t0/1\n", encoding="utf-8")
    r = call_diplotype(p)
    assert r.status == "no_input"


def test_multiallelic_site(tmp_path):
    """rs4244285 as a multiallelic record (ALT 'C,A'); the *2 ALT is index 2 -> hom *2/*2."""
    pos, ref, _alt, rsid = _POS["*2"]
    p = tmp_path / "ma.vcf"
    p.write_text(_HEADER + f"chr10\t{pos}\t{rsid}\t{ref}\tC,A\t.\tPASS\t.\tGT\t2/2\n", encoding="utf-8")
    r = call_diplotype(p)
    assert r.diplotype == "*2/*2"
    assert r.phenotype == "Poor Metabolizer"


# --- runner + CLI smoke -----------------------------------------------------

def test_runner_record_shape(tmp_path):
    from dna_decode.pgx.runner import call_cyp2c19
    vcf = _vcf(tmp_path, [("*2", "1/1")])
    rec = call_cyp2c19(vcf, sample_id="T1")
    assert rec["schema"] == "pgx-diplotype-call-v0"
    assert rec["sample_id"] == "T1"
    assert rec["diplotype"] == "*2/*2"
    assert rec["phenotype_abbrev"] == "PM"
    assert rec["phenotype_status"] == "ok"
    # v0.1 provenance: no overclaim of achieved independence
    assert rec["caller"]["calling_independently_validatable"] is True
    assert "pending" in rec["caller"]["independent_validation_status"]
    assert rec["caller"]["phenotype_is_faithful_to_cpic"] is True
    assert rec["caller"]["is_core_marker_proxy"] is True


# --- v0.1 sentinel layer + phase ambiguity + missing-sample -----------------

def test_sentinel_star4_withholds(tmp_path):
    """rs28399504 (*4) ALT alongside a *17 call -> *4b, not *17 -> phenotype WITHHELD (not increased-func)."""
    from dna_decode.pgx.runner import call_cyp2c19
    vcf = _vcf(tmp_path, [("*17", "0/1"), ("rs28399504", "0/1")])
    r = call_diplotype(vcf)
    assert r.phenotype_status == "phenotype_withheld"
    assert r.phenotype is None
    assert any(h["implies"] == "*4" for h in r.sentinel_hits)
    assert r.core_proxy_diplotype == "*1/*17"   # the proxy call is still visible
    # runner mirrors + CLI exits nonzero
    rec = call_cyp2c19(vcf)
    assert rec["phenotype_status"] == "phenotype_withheld" and rec["phenotype"] is None


def test_sentinel_star35_withholds(tmp_path):
    """rs12769205 ALT WITHOUT rs4244285 (*2) -> *35 mis-called *1 -> WITHHELD. (Core sites present as ref,
    as in a real VCF -- a *35 sample is hom-ref at the *2/*3/*17 positions.)"""
    vcf = _vcf(tmp_path, [("*2", "0/0"), ("*3", "0/0"), ("*17", "0/0"), ("rs12769205", "0/1")])
    r = call_diplotype(vcf)
    assert r.phenotype_status == "phenotype_withheld"
    assert any(h["implies"] == "*35" for h in r.sentinel_hits)


def test_na19122_style_star2_star35_withholds(tmp_path):
    """*2/*35: rs4244285 het + rs12769205 hom -> the *35 excess copy fires the sentinel -> WITHHELD
    (the real NA19122 case: GeT-RM *2/*35, v0 core proxy would say *1/*2)."""
    vcf = _vcf(tmp_path, [("*2", "1|0"), ("rs12769205", "1|1")])
    r = call_diplotype(vcf)
    assert r.phenotype_status == "phenotype_withheld"
    assert r.core_proxy_diplotype == "*1/*2"


def test_star2_alone_does_not_false_fire_star35(tmp_path):
    """*1/*2 with rs12769205 het (the *2 haplotype carries it) -> excess 0 -> NOT a *35 hit -> clean call."""
    vcf = _vcf(tmp_path, [("*2", "0/1"), ("rs12769205", "0/1")])
    r = call_diplotype(vcf)
    assert r.phenotype_status == "ok"
    assert r.diplotype == "*1/*2" and r.phenotype == "Intermediate Metabolizer"


def test_phase_ambiguous_unphased_two_hets(tmp_path):
    """Unphased *2 + *17 hets -> trans *2/*17 IM kept, but flagged phase_ambiguous w/ low confidence."""
    vcf = _vcf(tmp_path, [("*2", "0/1"), ("*17", "0/1")])
    r = call_diplotype(vcf)
    assert r.phenotype_status == "phase_ambiguous"
    assert r.phenotype_confidence == "low"
    assert r.phenotype == "Intermediate Metabolizer"     # standard trans call kept
    assert r.alternate_phenotype is not None


def test_missing_sample_raises(tmp_path):
    vcf = _vcf(tmp_path, [("*2", "1/1")], name="multi.vcf")
    with pytest.raises(ValueError, match="not found in VCF header"):
        call_diplotype(vcf, sample="NOSUCHSAMPLE")


def test_cli_smoke(tmp_path, capsys):
    from dna_decode.pgx.cli import main
    vcf = _vcf(tmp_path, [("*17", "1/1")])
    rc = main([str(vcf), "--sample-id", "UM1"])
    out = capsys.readouterr().out
    assert rc == 0
    assert "*17/*17" in out and "Ultrarapid" in out


def test_dispatch_via_dna_decode(tmp_path, capsys):
    from dna_decode.cli import main as decode_main
    vcf = _vcf(tmp_path, [("*2", "1/1")])
    rc = decode_main(["pgx", str(vcf), "--json-only"])
    assert rc == 0
    assert "pgx-diplotype-call-v0" in capsys.readouterr().out


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
