"""Reference-integrity + genome-mode tests for the HCMV caller (scripts/hcmv_caller.py).

THE load-bearing test: the committed Merlin NC_006273.2 UL97/UL54/UL56 CDS references MUST translate to the
catalog WT at every catalogued RESISTANCE position (zero tolerance) and at benign positions except the ONE
documented Merlin-vs-AD169 strain difference (UL54:897, a benign polymorphic site: Merlin L897 vs catalog
S897). A frame/coordinate error or a numbering offset fails LOUDLY here before any genome call — the same
reference-integrity discipline as the HIV/SARS-CoV-2 cells.
"""
from __future__ import annotations

import re
import shutil
from pathlib import Path

import pytest

from dna_decode.data.hcmv_amr import RESISTANCE_BY_GENE, BENIGN_BY_GENE

_REF_DIR = Path(__file__).resolve().parent.parent / "data" / "hcmv_ref"
_POINT = re.compile(r"^([ACDEFGHIKLMNPQRSTVWY])(\d+)([ACDEFGHIKLMNPQRSTVWY])$")
_EXP_LEN = {"UL97": 707, "UL54": 1242, "UL56": 850}
# documented Merlin-vs-AD169 strain difference at a BENIGN polymorphic site (does not affect any resistance call)
_MERLIN_WT_EXCEPTIONS = {("UL54", 897)}   # catalog S897 (AD169) vs Merlin L897

_CODON = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M', 'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*', 'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K', 'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W', 'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'}


def _translate(gene: str) -> str:
    seq = "".join(l.strip() for l in (_REF_DIR / f"HCMV_{gene}_Merlin_cds.fna").read_text().splitlines()
                  if not l.startswith(">")).upper().replace("U", "T")
    return "".join(_CODON.get(seq[i:i + 3], "X") for i in range(0, len(seq) - 2, 3))


def _catalog_wt(gene: str) -> dict[int, str]:
    wt: dict[int, str] = {}
    for m in list(RESISTANCE_BY_GENE[gene]) + list(BENIGN_BY_GENE[gene]):
        mm = _POINT.match(m)
        if mm:
            wt[int(mm.group(2))] = mm.group(1)
    return wt


@pytest.mark.parametrize("gene", ["UL97", "UL54", "UL56"])
def test_reference_is_clean_orf(gene):
    assert (_REF_DIR / f"HCMV_{gene}_Merlin_cds.fna").exists(), f"missing committed reference for {gene}"
    prot = _translate(gene).rstrip("*")
    assert len(prot) == _EXP_LEN[gene], f"{gene} should be {_EXP_LEN[gene]} aa, got {len(prot)}"
    assert "*" not in prot, f"no internal stop codon expected in the in-frame {gene} CDS"


def test_reference_wt_matches_catalog_at_every_resistance_position():
    # ZERO tolerance at resistance positions -- a resistance call depends on the WT being right here.
    for gene in ("UL97", "UL54", "UL56"):
        prot = _translate(gene)
        res_positions = {int(_POINT.match(m).group(2)): _POINT.match(m).group(1)
                         for m in RESISTANCE_BY_GENE[gene] if _POINT.match(m)}
        for pos, wt in res_positions.items():
            assert prot[pos - 1] == wt, f"{gene} resistance pos {pos}: catalog WT {wt}, Merlin {prot[pos-1]}"


def test_reference_wt_matches_benign_except_documented_merlin_polymorphism():
    mismatches = []
    for gene in ("UL97", "UL54", "UL56"):
        prot = _translate(gene)
        for pos, wt in _catalog_wt(gene).items():
            if prot[pos - 1] != wt and (gene, pos) not in _MERLIN_WT_EXCEPTIONS:
                mismatches.append((gene, pos, wt, prot[pos - 1]))
    assert not mismatches, f"undocumented Merlin-vs-catalog WT mismatches: {mismatches}"


def test_documented_exception_is_real_and_benign():
    # the pinned exception must actually be a Merlin/L vs catalog/S benign-position difference (not a silent pass)
    prot = _translate("UL54")
    assert prot[896] == "L", "UL54:897 exception expects Merlin L897"
    assert "S897P" in BENIGN_BY_GENE["UL54"], "UL54:897 must be a BENIGN entry (not resistance)"


# ---- genome-mode functional test (needs BLAST+; skips cleanly without it) ----

def _has_blast():
    return shutil.which("blastn") is not None and shutil.which("makeblastdb") is not None


def test_genome_mode_degrades_without_blast(tmp_path, monkeypatch):
    # force the no-BLAST path -> INDETERMINATE (never a silent wrong call)
    import scripts.hcmv_caller as hc
    monkeypatch.setattr(hc, "observed_substitutions", lambda *a, **k: None)
    g = tmp_path / "g.fna"
    g.write_text(">c\nACGT\n")
    call = hc.call_hcmv_target(str(g), "ganciclovir")
    assert call.prediction == "INDETERMINATE" and "BLAST" in call.caveat


@pytest.mark.skipif(not _has_blast(), reason="BLAST+ not installed")
def test_genome_mode_calls_planted_mutation(tmp_path):
    # a synthetic 'genome' = the UL97 CDS with a planted M460V -> R for ganciclovir
    from scripts.hcmv_caller import call_hcmv_target, default_ref_for
    ref = "".join(l.strip() for l in default_ref_for("UL97").read_text().splitlines() if not l.startswith(">"))
    # M460V: codon 460 (nt 1378..1380) Met(ATG) -> Val(GTG)
    i = (460 - 1) * 3
    assert ref[i:i + 3].upper() == "ATG"
    mut = ref[:i] + "GTG" + ref[i + 3:]
    g = tmp_path / "genome.fna"
    g.write_text(f">syn_UL97_M460V\n{mut}\n")
    call = call_hcmv_target(str(g), "ganciclovir")
    assert call.prediction == "R" and "UL97:M460V" in call.determinants
