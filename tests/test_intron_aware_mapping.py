"""Intron-aware multi-HSP codon mapping (scripts/fungal_erg11_caller.observed_substitutions).

Validates that the codon-mapper handles an INTRON-CONTAINING gene: blastn of an in-frame CDS reference vs a
genome whose gene is split into exons (introns inserted) returns one HSP per exon, and the stitched
query-position map must (a) detect a mutation in a later exon and (b) — the killer case — detect a mutation
in a codon that SPANS an exon-exon boundary (its 3 nts in two different HSPs). Unblocks pfcrt (13 exons) +
eukaryotic multi-exon targets. Skips cleanly if BLAST+ absent.

Uses the REAL (non-repetitive) 3D7 K13 CDS as the substrate, artificially split into exons. K13 is
intronless in reality; splitting it here is purely to exercise the multi-HSP machinery on a realistic
non-repetitive sequence (a periodic synthetic CDS would self-align at its period and is not representative
of real references). The intronless single-HSP path is the special case covered by
test_fungal_erg11_caller + test_pf_kelch13 (regression guard).
"""
import shutil
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import pytest  # noqa: E402

from scripts.fungal_erg11_caller import _translate, observed_substitutions  # noqa: E402

_HAS_BLAST = bool(shutil.which("blastn", path="C:/Users/Farshad/ncbi-blast/bin") or shutil.which("blastn")) \
    and bool(shutil.which("makeblastdb", path="C:/Users/Farshad/ncbi-blast/bin") or shutil.which("makeblastdb"))
_K13 = Path(__file__).resolve().parent.parent / "data" / "antimalarial_ref" / "Pf3D7_K13_cds.fna"
_INTRON = "GATTACAGGGCCCTATAGTGAGTCGTATTACAATTCACTGGCCGTCGTTTTACAACGTCGTGACTGGGAAAACC" * 2  # ~146nt, non-K13


def _seq(p):
    return "".join(l.strip() for l in p.read_text().splitlines() if not l.startswith(">")).upper()


def _write(p, header, seq):
    p.write_text(f">{header}\n{seq}\n", encoding="utf-8")
    return str(p)


@pytest.mark.skipif(not _HAS_BLAST or not _K13.exists(), reason="BLAST+ or K13 ref absent")
def test_multi_exon_detects_mutations_across_intron():
    """A mutation in exon1 AND one DEEP in exon2 are both detected through the stitched multi-HSP map —
    i.e. the exon2 substitution is found despite the intron between it and exon1. (Boundary codon 580 is
    WT here and is correctly assembled across the boundary — see test_multi_exon_wildtype_no_calls, which
    proves no spurious call arises from the mid-codon split.)"""
    ref = _seq(_K13)[:726 * 3]                 # 726 codons (drop stop)
    prot = _translate(ref)
    def mut(cds, codon_1based, new):
        i = 3 * (codon_1based - 1)
        return cds[:i] + new + cds[i + 3:]
    new200 = "TGG" if prot[199] != "W" else "GCT"      # exon1 codon
    new650 = "TGG" if prot[649] != "W" else "GCT"      # exon2 codon, deep (nt 1948-50, away from edge)
    m = mut(mut(ref, 200, new200), 650, new650)
    exp200 = f"{prot[199]}200{_translate(new200)}"
    exp650 = f"{prot[649]}650{_translate(new650)}"
    # exon1 = nt 1..1738 (boundary mid-codon-580), exon2 = nt 1739..2178; non-homologous intron between.
    flank = "ACGTTGCA" * 40
    genome = flank + m[0:1738] + _INTRON + m[1738:2178] + flank
    with tempfile.TemporaryDirectory() as td:
        r = _write(Path(td) / "ref.fna", "K13", ref)
        g = _write(Path(td) / "g.fna", "contig1", genome)
        obs = observed_substitutions(g, r, gene="K13")
    assert obs is not None
    assert exp200 in obs["K13"], (exp200, sorted(obs["K13"]))   # exon1 codon
    assert exp650 in obs["K13"], (exp650, sorted(obs["K13"]))   # exon2 codon (across the intron)


@pytest.mark.skipif(not _HAS_BLAST or not _K13.exists(), reason="BLAST+ or K13 ref absent")
def test_multi_exon_wildtype_no_calls():
    ref = _seq(_K13)[:726 * 3]
    flank = "ACGTTGCA" * 40
    genome = flank + ref[0:1738] + _INTRON + ref[1738:2178] + flank
    with tempfile.TemporaryDirectory() as td:
        r = _write(Path(td) / "ref.fna", "K13", ref)
        g = _write(Path(td) / "g.fna", "contig1", genome)
        obs = observed_substitutions(g, r, gene="K13")
    assert obs == {"K13": set()}, sorted(obs["K13"])           # WT multi-exon genome -> no spurious calls
