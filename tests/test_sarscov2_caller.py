"""Reference-integrity gate + offline-degrade for the SARS-CoV-2 Mpro caller (mirrors test_hiv_rt_caller).

THE load-bearing test: the committed Mpro reference (Wuhan-Hu-1 NC_045512.2 nsp5) MUST translate to the Mpro
WT at every catalogued position (`MPRO_WT`) AND show the catalytic dyad H41/C145 + the key nirmatrelvir
position E166 — a frame/coordinate error fails loudly before any call. Real-blastn tests skip when BLAST+
is absent (offline-safe)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from dna_decode.data.sarscov2_amr import MPRO_WT  # noqa: E402

_REF = Path(__file__).resolve().parent.parent / "data" / "sarscov2_ref" / "SARSCoV2_Mpro_NC045512_cds.fna"

_CODON = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M", "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T", "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*", "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K", "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W", "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R", "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


def _translate_ref() -> str:
    seq = "".join(l.strip() for l in _REF.read_text().splitlines() if not l.startswith(">")).upper()
    return "".join(_CODON.get(seq[i:i + 3], "X") for i in range(0, len(seq) - 2, 3))


def test_reference_exists_and_is_clean_orf():
    assert _REF.exists(), f"committed Mpro reference missing: {_REF}"
    prot = _translate_ref()
    assert len(prot) == 306, f"Mpro should be 306 aa, got {len(prot)}"
    assert prot[:6] == "SGFRKM", f"Mpro N-terminus should be SGFRKM, got {prot[:6]}"
    assert "*" not in prot, "no internal stop codon expected in the in-frame Mpro CDS"


def test_reference_catalytic_dyad_and_e166():
    """Coordinate sanity: the Mpro catalytic dyad His41/Cys145 + the key nirmatrelvir position Glu166."""
    prot = _translate_ref()
    assert prot[40] == "H", "Mpro His41 (catalytic) mismatch -> wrong frame/coordinates"
    assert prot[144] == "C", "Mpro Cys145 (catalytic) mismatch -> wrong frame/coordinates"
    assert prot[165] == "E", "Mpro Glu166 (nirmatrelvir contact) mismatch -> wrong frame/coordinates"
    assert prot[49] == "L", "Mpro Leu50 (L50F resistance) WT mismatch"


def test_reference_matches_catalog_wt_at_every_position():
    """Integrity gate: the committed reference WT == MPRO_WT at every catalogued position (else fail loud)."""
    prot = _translate_ref()
    mism = {p: (wt, prot[p - 1]) for p, wt in MPRO_WT.items() if prot[p - 1] != wt}
    assert not mism, f"reference vs catalog WT mismatch at {mism}"


def test_caller_offline_indeterminate(monkeypatch):
    """Absent BLAST+ -> INDETERMINATE (offline-safe), not a crash."""
    import scripts.sarscov2_caller as sc
    monkeypatch.setattr(sc, "observed_substitutions", lambda *a, **k: None)
    c = sc.call_sarscov2_target("x.fna", str(_REF), "nirmatrelvir")
    assert c.prediction == "INDETERMINATE" and "blast" in c.caveat.lower()


if __name__ == "__main__":
    import pytest
    sys.exit(pytest.main([__file__, "-q"]))
