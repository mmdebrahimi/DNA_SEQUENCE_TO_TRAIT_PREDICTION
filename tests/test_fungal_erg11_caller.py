"""G0 validation for the fungal ERG11 BLAST caller (scripts/fungal_erg11_caller).

Exercises the REAL makeblastdb + blastn pipeline + codon-mapping against a KNOWN planted mutation
(synthetic in-frame ERG11-like CDS with Y at codon 132; a 'resistant genome' carrying TTC=F at codon 132).
Skips cleanly if BLAST+ is not installed (offline-safe, like the pathotype vf_diff tests).

This validates the MACHINERY (frame + position-mapping + catalog wiring). Real-C.auris-reference +
real-resistant-genome validation is EP-7 step-3 (G0-completion), not this test.
"""
import shutil
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import pytest  # noqa: E402

from scripts.fungal_erg11_caller import call_erg11  # noqa: E402

_HAS_BLAST = bool(shutil.which("blastn", path="C:/Users/Farshad/ncbi-blast/bin") or shutil.which("blastn")) \
    and bool(shutil.which("makeblastdb", path="C:/Users/Farshad/ncbi-blast/bin") or shutil.which("makeblastdb"))

# varied, stop-free codons to build a non-repetitive CDS blastn can align cleanly
_CYC = ["GCT", "GAT", "TTT", "GGT", "CAT", "ATT", "AAA", "CTT", "CCT", "CAA",
        "CGT", "TCT", "ACT", "GTT", "TGG", "GTG", "AAT", "GAA", "GGA", "CTG"]


def _synthetic_cds(n_codons: int = 450) -> str:
    """In-frame CDS: codon1=ATG, codon132=TAC(Y), rest cycled (stop-free)."""
    codons = ["ATG"]
    for i in range(1, n_codons):
        codons.append("TAC" if i == 131 else _CYC[i % len(_CYC)])  # i=131 → 1-based pos 132
    return "".join(codons)


def _write(p: Path, header: str, seq: str) -> str:
    p.write_text(f">{header}\n{seq}\n", encoding="utf-8")
    return str(p)


@pytest.mark.skipif(not _HAS_BLAST, reason="BLAST+ (blastn/makeblastdb) not installed")
def test_caller_detects_planted_Y132F():
    cds = _synthetic_cds()
    flank = "ACGT" * 60
    # resistant genome: codon 132 (nt 394..396, 1-based) TAC -> TTC (Y->F)
    pos = 3 * 132 - 3  # 0-based start of codon 132
    r_cds = cds[:pos] + "TTC" + cds[pos + 3:]
    with tempfile.TemporaryDirectory() as td:
        ref = _write(Path(td) / "erg11_ref.fna", "ERG11_ref_CDS", cds)
        r_genome = _write(Path(td) / "R.fna", "contig1", flank + r_cds + flank)
        s_genome = _write(Path(td) / "S.fna", "contig1", flank + cds + flank)
        rc = call_erg11(r_genome, ref, "fluconazole")
        sc = call_erg11(s_genome, ref, "fluconazole")
    assert rc.prediction == "R", rc
    assert "ERG11:Y132F" in rc.determinants, rc
    assert sc.prediction == "S", sc            # wild-type → susceptible
    assert sc.undetectable_mechanisms          # S surfaces efflux/aneuploidy blind spots


_REF = Path(__file__).resolve().parent.parent / "data" / "fungal_ref" / "Cauris_ERG11_cds.fna"
_REAL = [  # (committed public GenBank isolate allele, truth, expected call)
    ("Cauris_ERG11_PV630306_WT.fna", "WT", "S"),
    ("Cauris_ERG11_PV630305_Y132F.fna", "Y132F", "R"),
    ("Cauris_ERG11_PV630302_K143R.fna", "K143R", "R"),
]


@pytest.mark.skipif(not _HAS_BLAST, reason="BLAST+ (blastn/makeblastdb) not installed")
@pytest.mark.skipif(not _REF.exists(), reason="real C. auris ERG11 reference fixture absent")
@pytest.mark.parametrize("fname,truth,expected", _REAL)
def test_caller_on_real_cauris_alleles(fname, truth, expected):
    """G0-COMPLETION: real C. auris ERG11 reference (RefSeq XM_029033208.2, numbering Y132/K143/V125
    confirmed) vs real GenBank isolate alleles carrying DOCUMENTED Y132F / K143R. Validates that the
    catalog numbering matches the real reference on real mutations — not just a planted synthetic one."""
    genome = str(_REF.parent / fname)
    c = call_erg11(genome, str(_REF), "fluconazole")
    assert c.prediction == expected, c
    if truth != "WT":
        assert f"ERG11:{truth}" in c.determinants, c
    else:
        assert c.undetectable_mechanisms  # S surfaces efflux/aneuploidy blind spots


def test_caller_indeterminate_without_blast(monkeypatch=None):
    # force the no-BLAST path → INDETERMINATE (offline-safe contract)
    import scripts.fungal_erg11_caller as m
    orig = m._find
    m._find = lambda tool: None
    try:
        with tempfile.TemporaryDirectory() as td:
            g = _write(Path(td) / "g.fna", "c", "ACGT" * 50)
            r = _write(Path(td) / "r.fna", "ref", "ATG" + "GCT" * 50)
            c = call_erg11(g, r, "fluconazole")
        assert c.prediction == "INDETERMINATE" and "blastn" in c.caveat
    finally:
        m._find = orig


if __name__ == "__main__":
    if not _HAS_BLAST:
        print("SKIP — BLAST+ not installed")
    test_caller_indeterminate_without_blast()
    print("PASS test_caller_indeterminate_without_blast")
    if _HAS_BLAST:
        test_caller_detects_planted_Y132F()
        print("PASS test_caller_detects_planted_Y132F")
