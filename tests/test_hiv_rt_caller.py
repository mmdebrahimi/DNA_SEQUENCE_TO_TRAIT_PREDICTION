"""Validation for the HIV-1 RT BLAST genome-mode caller (scripts/hiv_rt_caller).

Three layers:
  1. REFERENCE SELF-CHECK (no BLAST) — the integrity gate: translate the committed HXB2 RT CDS reference
     and assert its WT residue matches the Stanford-sourced catalog (`_RT_WT` + `NRTI_RT_WT`) at EVERY
     catalogued NNRTI+NRTI DRM position. If the reference is out of frame, wrong-coordinate, or differs
     from consensus-B at a DRM position, a substitution call would be wrong — this test fails loudly first.
  2. OFFLINE DEGRADE (no BLAST) — force the no-BLAST path → INDETERMINATE (offline-safe contract).
  3. PLANTED-MUTATION CALL (BLAST-gated) — BLAST the real HXB2 RT reference vs a genome built from that
     reference with a known DRM planted (NNRTI K103N / NRTI M184V), validating frame+numbering+catalog
     wiring on the REAL reference (the HIV analogue of the fungal G0-completion test).

Skips the BLAST-dependent tests cleanly when BLAST+ is absent (like the fungal/influenza caller tests).
"""
import shutil
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import pytest  # noqa: E402

from dna_decode.data.hiv_amr import _RT_WT, NRTI_RT_WT  # noqa: E402
from scripts.fungal_erg11_caller import _read_single_fasta, _translate  # noqa: E402
from scripts.hiv_rt_caller import call_hiv_rt  # noqa: E402

_HAS_BLAST = bool(shutil.which("blastn", path="C:/Users/Farshad/ncbi-blast/bin") or shutil.which("blastn")) \
    and bool(shutil.which("makeblastdb", path="C:/Users/Farshad/ncbi-blast/bin") or shutil.which("makeblastdb"))

_REF = Path(__file__).resolve().parent.parent / "data" / "hiv_ref" / "HIV1_RT_HXB2_cds.fna"


def _write(p: Path, header: str, seq: str) -> str:
    p.write_text(f">{header}\n{seq}\n", encoding="utf-8")
    return str(p)


def _mutate_codon(seq: str, pos_1based: int, new_codon: str) -> str:
    i = 3 * (pos_1based - 1)
    return seq[:i] + new_codon + seq[i + 3:]


def test_reference_exists():
    assert _REF.exists(), f"committed HXB2 RT reference missing at {_REF}"


def test_reference_wt_matches_catalog_at_every_drm_position():
    """The integrity gate — committed reference WT must equal the catalog WT at every DRM position."""
    prot = _translate(_read_single_fasta(str(_REF)))
    assert len(prot) == 560, f"HXB2 RT p66 should be 560 codons, got {len(prot)}"
    wt = {**NRTI_RT_WT, **_RT_WT}
    mismatches = {pos: (exp, prot[pos - 1]) for pos, exp in wt.items() if prot[pos - 1] != exp}
    assert not mismatches, f"reference WT != catalog WT (pos: (expected, got)): {mismatches}"


def test_reference_in_frame_no_internal_stops():
    prot = _translate(_read_single_fasta(str(_REF)))
    assert "*" not in prot[:-1], "internal stop codon → reference is out of frame"


def test_caller_indeterminate_without_blast():
    """Force the no-BLAST path → INDETERMINATE (offline-safe contract)."""
    import scripts.fungal_erg11_caller as m
    orig = m._find
    m._find = lambda tool: None
    try:
        with tempfile.TemporaryDirectory() as td:
            g = _write(Path(td) / "g.fna", "c", "ACGT" * 50)
            c = call_hiv_rt(g, str(_REF), "efavirenz")
        assert c.prediction == "INDETERMINATE" and "blastn" in c.caveat
    finally:
        m._find = orig


@pytest.mark.skipif(not _HAS_BLAST, reason="BLAST+ (blastn/makeblastdb) not installed")
@pytest.mark.skipif(not _REF.exists(), reason="HXB2 RT reference fixture absent")
def test_caller_detects_planted_nnrti_K103N():
    ref_seq = _read_single_fasta(str(_REF))
    flank = "ACGT" * 60
    r_cds = _mutate_codon(ref_seq, 103, "AAC")   # K(AAA) -> N(AAC)
    with tempfile.TemporaryDirectory() as td:
        r_genome = _write(Path(td) / "R.fna", "contig1", flank + r_cds + flank)
        s_genome = _write(Path(td) / "S.fna", "contig1", flank + ref_seq + flank)
        rc = call_hiv_rt(r_genome, str(_REF), "efavirenz", is_nrti=False)
        sc = call_hiv_rt(s_genome, str(_REF), "efavirenz", is_nrti=False)
    assert rc.prediction == "R", rc
    assert "RT:K103N" in rc.determinants, rc
    assert sc.prediction == "S", sc
    assert sc.undetectable_mechanisms          # S surfaces NNRTI blind spots


@pytest.mark.skipif(not _HAS_BLAST, reason="BLAST+ (blastn/makeblastdb) not installed")
@pytest.mark.skipif(not _REF.exists(), reason="HXB2 RT reference fixture absent")
def test_caller_detects_planted_nrti_M184V():
    ref_seq = _read_single_fasta(str(_REF))
    flank = "ACGT" * 60
    r_cds = _mutate_codon(ref_seq, 184, "GTG")   # M(ATG) -> V(GTG); 184 is an NRTI major position
    with tempfile.TemporaryDirectory() as td:
        r_genome = _write(Path(td) / "R.fna", "contig1", flank + r_cds + flank)
        s_genome = _write(Path(td) / "S.fna", "contig1", flank + ref_seq + flank)
        rc = call_hiv_rt(r_genome, str(_REF), "lamivudine", is_nrti=True)
        sc = call_hiv_rt(s_genome, str(_REF), "lamivudine", is_nrti=True)
    assert rc.prediction == "R", rc
    assert "RT:M184V" in rc.determinants, rc
    assert sc.prediction == "S", sc
    assert sc.undetectable_mechanisms          # S surfaces NRTI blind spots


if __name__ == "__main__":
    test_reference_exists()
    test_reference_wt_matches_catalog_at_every_drm_position()
    test_reference_in_frame_no_internal_stops()
    test_caller_indeterminate_without_blast()
    print("PASS no-BLAST tests (reference self-check + offline degrade)")
    if _HAS_BLAST:
        test_caller_detects_planted_nnrti_K103N()
        test_caller_detects_planted_nrti_M184V()
        print("PASS BLAST planted-mutation tests")
    else:
        print("SKIP — BLAST+ not installed")
