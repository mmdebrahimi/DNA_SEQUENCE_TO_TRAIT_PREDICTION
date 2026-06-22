"""Validation for the HIV PI/INSTI/CAI genome-mode references + the generic RT/PR/IN/CA BLAST caller.

Layers (mirroring tests/test_hiv_rt_caller.py for the RT cell):
  1. REFERENCE SELF-CHECK (no BLAST) — the integrity gate: translate each committed HXB2 CDS reference
     (PR/IN/CA) and assert its WT residue matches the catalog `cls.wt` at EVERY catalogued major position.
     A frame/coordinate error or an isolate that differs from consensus-B at a major position fails here.
  2. OFFLINE DEGRADE (no BLAST) — force the no-BLAST path -> INDETERMINATE.
  3. PLANTED-MUTATION CALL (BLAST-gated) — BLAST the real reference vs a genome built from it with a known
     DRM planted (PI V82A / INSTI Q148H / CAI M66I), validating frame+numbering+catalog wiring per gene.
"""
import shutil
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import pytest  # noqa: E402

from dna_decode.data.hiv_amr import CAI_CLASS, INSTI_CLASS, PI_CLASS  # noqa: E402
from scripts.fungal_erg11_caller import _read_single_fasta, _translate  # noqa: E402
from scripts.hiv_rt_caller import call_hiv_target  # noqa: E402

_HAS_BLAST = bool(shutil.which("blastn", path="C:/Users/Farshad/ncbi-blast/bin") or shutil.which("blastn")) \
    and bool(shutil.which("makeblastdb", path="C:/Users/Farshad/ncbi-blast/bin") or shutil.which("makeblastdb"))
_REFDIR = Path(__file__).resolve().parent.parent / "data" / "hiv_ref"

# (class, ref filename, gene, n_codons, drug, planted (pos,new_codon), expected determinant)
_CASES = [
    (PI_CLASS, "HIV1_PR_HXB2_cds.fna", "PR", 99, "lopinavir", (82, "GCC"), "PR:V82A"),     # V82A
    (INSTI_CLASS, "HIV1_IN_HXB2_cds.fna", "IN", 288, "dolutegravir", (148, "CAC"), "IN:Q148H"),  # Q148H
    (CAI_CLASS, "HIV1_CA_HXB2_cds.fna", "CA", 231, "lenacapavir", (66, "ATC"), "CA:M66I"),  # M66I
]


def _ref(name: str) -> Path:
    return _REFDIR / name


def _mutate_codon(seq: str, pos_1based: int, new_codon: str) -> str:
    i = 3 * (pos_1based - 1)
    return seq[:i] + new_codon + seq[i + 3:]


@pytest.mark.parametrize("cls,fname,gene,ncod,_drug,_planted,_det", _CASES)
def test_reference_wt_matches_catalog(cls, fname, gene, ncod, _drug, _planted, _det):
    """Integrity gate — committed reference WT == catalog WT at every catalogued major position."""
    assert _ref(fname).exists(), f"committed reference missing: {_ref(fname)}"
    prot = _translate(_read_single_fasta(str(_ref(fname))))
    assert len(prot) == ncod, f"{gene}: expected {ncod} codons, got {len(prot)}"
    assert "*" not in prot[:-1], f"{gene}: internal stop -> out of frame"
    mismatches = {p: (exp, prot[p - 1]) for p, exp in cls.wt.items() if prot[p - 1] != exp}
    assert not mismatches, f"{gene} reference WT != catalog WT (pos: (expected, got)): {mismatches}"


@pytest.mark.parametrize("cls,fname,gene,ncod,drug,_planted,_det", _CASES)
def test_caller_indeterminate_without_blast(cls, fname, gene, ncod, drug, _planted, _det):
    import scripts.fungal_erg11_caller as m
    orig = m._find
    m._find = lambda tool: None
    try:
        with tempfile.TemporaryDirectory() as td:
            g = Path(td) / "g.fna"
            g.write_text(">c\n" + "ACGT" * 50 + "\n", encoding="utf-8")
            c = call_hiv_target(str(g), str(_ref(fname)), drug, gene)
        assert c.prediction == "INDETERMINATE" and "blastn" in c.caveat
    finally:
        m._find = orig


@pytest.mark.skipif(not _HAS_BLAST, reason="BLAST+ (blastn/makeblastdb) not installed")
@pytest.mark.parametrize("cls,fname,gene,ncod,drug,planted,det", _CASES)
def test_caller_detects_planted_mutation(cls, fname, gene, ncod, drug, planted, det):
    ref_seq = _read_single_fasta(str(_ref(fname)))
    flank = "ACGT" * 60
    pos, new_codon = planted
    r_cds = _mutate_codon(ref_seq, pos, new_codon)
    with tempfile.TemporaryDirectory() as td:
        r_genome = Path(td) / "R.fna"
        r_genome.write_text(">contig1\n" + flank + r_cds + flank + "\n", encoding="utf-8")
        s_genome = Path(td) / "S.fna"
        s_genome.write_text(">contig1\n" + flank + ref_seq + flank + "\n", encoding="utf-8")
        rc = call_hiv_target(str(r_genome), str(_ref(fname)), drug, gene)
        sc = call_hiv_target(str(s_genome), str(_ref(fname)), drug, gene)
    assert rc.prediction == "R", rc
    assert det in rc.determinants, rc
    assert sc.prediction == "S", sc
    assert sc.undetectable_mechanisms


if __name__ == "__main__":
    for case in _CASES:
        test_reference_wt_matches_catalog(*case)
        test_caller_indeterminate_without_blast(*case)
    print("PASS no-BLAST tests (reference self-checks + offline degrade)")
    if _HAS_BLAST:
        for case in _CASES:
            test_caller_detects_planted_mutation(*case)
        print("PASS BLAST planted-mutation tests")
    else:
        print("SKIP — BLAST+ not installed")
