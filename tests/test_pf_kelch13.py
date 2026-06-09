"""G0 validation for the P. falciparum K13 caller (scripts/pf_kelch13_caller) — the protozoan kingdom jump.

Exercises the REAL makeblastdb + blastn pipeline + codon-mapping against a KNOWN planted mutation: a
synthetic in-frame K13-like CDS (726 codons; codon 580 = TGT = Cys = WT) and a 'resistant genome' carrying
TAC = Tyr at codon 580 (the canonical C580Y artemisinin marker). Skips cleanly if BLAST+ is absent
(offline-safe, like the fungal + pathotype tests). Validates the MACHINERY + catalog wiring; real-3D7-K13
reference + real C580Y isolate validation is G0-completion (the data/antimalarial_ref/ fixtures, optional).
"""
import shutil
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import pytest  # noqa: E402

from dna_decode.data.antimalarial_amr import (  # noqa: E402
    call_from_observed_substitutions, is_resistance_mutation, supported_antimalarial_drugs,
)
from scripts.pf_kelch13_caller import call_kelch13  # noqa: E402

_HAS_BLAST = bool(shutil.which("blastn", path="C:/Users/Farshad/ncbi-blast/bin") or shutil.which("blastn")) \
    and bool(shutil.which("makeblastdb", path="C:/Users/Farshad/ncbi-blast/bin") or shutil.which("makeblastdb"))

_CYC = ["GCT", "GAT", "TTT", "GGT", "CAT", "ATT", "AAA", "CTT", "CCT", "CAA",
        "CGT", "TCT", "ACT", "GTT", "TGG", "GTG", "AAT", "GAA", "GGA", "CTG"]


def _synthetic_k13_cds(n_codons: int = 726) -> str:
    """In-frame CDS: codon1=ATG, codon580=TGT(Cys, WT), rest cycled (stop-free)."""
    codons = ["ATG"]
    for i in range(1, n_codons):
        codons.append("TGT" if i == 579 else _CYC[i % len(_CYC)])  # i=579 → 1-based pos 580
    return "".join(codons)


def _write(p: Path, header: str, seq: str) -> str:
    p.write_text(f">{header}\n{seq}\n", encoding="utf-8")
    return str(p)


# ---------- pure catalog (no BLAST) ----------
def test_catalog_c580y_is_validated_marker():
    assert is_resistance_mutation("artemisinin", "K13", "C580Y")
    assert is_resistance_mutation("artemisinin", "K13", "R561H")    # African-emergence marker
    assert not is_resistance_mutation("artemisinin", "K13", "C580C")  # synonymous-ish non-marker
    assert "artemisinin" in supported_antimalarial_drugs()


def test_call_from_observed_R_and_S():
    r = call_from_observed_substitutions("artemisinin", {"K13": {"C580Y"}})
    assert r.prediction == "R" and "K13:C580Y" in r.determinants
    s = call_from_observed_substitutions("artemisinin", {"K13": set()})
    assert s.prediction == "S" and s.undetectable_mechanisms        # S surfaces non-K13/partner blind spots


def test_call_unknown_drug_indeterminate():
    c = call_from_observed_substitutions("chloroquine", {"K13": {"C580Y"}})
    assert c.prediction == "INDETERMINATE"


# ---------- real BLAST machinery (skip if no BLAST+) ----------
@pytest.mark.skipif(not _HAS_BLAST, reason="BLAST+ (blastn/makeblastdb) not installed")
def test_caller_detects_planted_C580Y():
    cds = _synthetic_k13_cds()
    flank = "ACGT" * 60
    pos = 3 * 580 - 3  # 0-based start of codon 580
    r_cds = cds[:pos] + "TAC" + cds[pos + 3:]      # TGT(Cys) -> TAC(Tyr) = C580Y
    with tempfile.TemporaryDirectory() as td:
        ref = _write(Path(td) / "k13_ref.fna", "K13_ref_CDS", cds)
        r_genome = _write(Path(td) / "R.fna", "contig1", flank + r_cds + flank)
        s_genome = _write(Path(td) / "S.fna", "contig1", flank + cds + flank)
        rc = call_kelch13(r_genome, ref, "artemisinin")
        sc = call_kelch13(s_genome, ref, "artemisinin")
    assert rc.prediction == "R", rc
    assert "K13:C580Y" in rc.determinants, rc
    assert sc.prediction == "S", sc
    assert sc.undetectable_mechanisms


# ---------- G0-COMPLETION: real 3D7 K13 reference numbering ----------
_REAL_REF = Path(__file__).resolve().parent.parent / "data" / "antimalarial_ref" / "Pf3D7_K13_cds.fna"


@pytest.mark.skipif(not _HAS_BLAST, reason="BLAST+ (blastn/makeblastdb) not installed")
@pytest.mark.skipif(not _REAL_REF.exists(), reason="real 3D7 K13 reference fixture absent")
def test_real_3d7_reference_numbering_C580Y():
    """Validates the catalog numbering matches the REAL 3D7 K13 reference (NCBI XM_001350122.1, 726aa,
    WT Cys@580) — not just a synthetic CDS. WT real ref → S; the real ref with codon 580 TGT→TAC → C580Y/R.
    Proves C580 sits at the catalogued position on the actual reference (the G0-completion numbering check)."""
    ref_seq = "".join(l.strip() for l in _REAL_REF.read_text().splitlines() if not l.startswith(">")).upper()
    flank = "ACGT" * 60
    pos = 3 * 580 - 3
    assert ref_seq[pos:pos + 3] in ("TGT", "TGC"), f"real ref codon580 not Cys: {ref_seq[pos:pos+3]}"
    r_seq = ref_seq[:pos] + "TAC" + ref_seq[pos + 3:]   # C580Y on the REAL reference
    with tempfile.TemporaryDirectory() as td:
        ref = _write(Path(td) / "ref.fna", "K13ref", ref_seq)
        rg = _write(Path(td) / "R.fna", "contig1", flank + r_seq + flank)
        sg = _write(Path(td) / "S.fna", "contig1", flank + ref_seq + flank)
        rc = call_kelch13(rg, ref, "artemisinin")
        sc = call_kelch13(sg, ref, "artemisinin")
    assert rc.prediction == "R" and "K13:C580Y" in rc.determinants, rc
    assert sc.prediction == "S", sc


# ---------- offline-safe contract ----------
def test_caller_indeterminate_without_blast():
    import scripts.fungal_erg11_caller as m
    orig = m._find
    m._find = lambda tool: None
    try:
        with tempfile.TemporaryDirectory() as td:
            g = _write(Path(td) / "g.fna", "c", "ACGT" * 50)
            r = _write(Path(td) / "r.fna", "ref", "ATG" + "GCT" * 50)
            c = call_kelch13(g, r, "artemisinin")
        assert c.prediction == "INDETERMINATE" and "blastn" in c.caveat
    finally:
        m._find = orig
