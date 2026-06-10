"""Unit tests for the expression_context detector (ISAba1-upstream-of-blaOXA-51 junction).

Builds synthetic assembly contigs from the committed refs (ISAba1 + OXA-51-family) and asserts the PRIMARY
frozen rule (same-contig + strand-aware upstream proximity, NO IS-orientation). Real-blastn cases skipif no
BLAST+; the offline-degrade case runs without BLAST. The primary path must IGNORE orientation (the
orientation refinement is default-off and separately re-falsified).
"""
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import pytest  # noqa: E402

from dna_decode.eval.expression_context import detect_is_upstream_junction  # noqa: E402

_REF_DIR = Path(__file__).resolve().parent.parent / "data" / "isaba1_ref"
_IS_REF = _REF_DIR / "ISAba1_ref.fna"
_OXA_REF = _REF_DIR / "OXA51fam_ref.fna"
_HAS_BLAST = bool(shutil.which("blastn", path="C:/Users/Farshad/ncbi-blast/bin") or shutil.which("blastn")) \
    and bool(shutil.which("makeblastdb", path="C:/Users/Farshad/ncbi-blast/bin") or shutil.which("makeblastdb"))
_blast_or_skip = pytest.mark.skipif(not (_HAS_BLAST and _IS_REF.exists() and _OXA_REF.exists()),
                                    reason="BLAST+ or ISAba1/OXA refs absent")


def _seq(fasta: Path) -> str:
    return "".join(l.strip() for l in fasta.read_text().splitlines() if not l.startswith(">")).upper()


def _rc(s: str) -> str:
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def _write_contigs(path: Path, contigs: dict[str, str]) -> None:
    path.write_text("".join(f">{name}\n{seq}\n" for name, seq in contigs.items()), encoding="utf-8")


def test_offline_degrade_no_raise(tmp_path):
    """blastn forced absent -> status unavailable, signal False, never raises."""
    g = tmp_path / "g.fna"
    g.write_text(">c1\n" + "ACGT" * 50 + "\n", encoding="utf-8")
    out = detect_is_upstream_junction(g, is_ref=_IS_REF, target_ref=_OXA_REF, blastn_bin="/nonexistent/blastn")
    assert out["status"] == "unavailable" and out["signal"] is False


@_blast_or_skip
def test_junction_positive_same_contig(tmp_path):
    """ISAba1 placed 200 bp upstream of OXA on the SAME contig (+ strand) -> signal True."""
    isa, oxa = _seq(_IS_REF), _seq(_OXA_REF)
    spacer = "AT" * 100                                   # 200 bp upstream gap (< 400)
    flank = "GC" * 80
    contig = flank + isa + spacer + oxa + flank
    g = tmp_path / "pos.fna"
    _write_contigs(g, {"contig1": contig})
    out = detect_is_upstream_junction(g, is_ref=_IS_REF, target_ref=_OXA_REF)
    assert out["status"] == "ok" and out["signal"] is True
    assert out["evidence"]["junction"]["contig"] == "contig1"
    assert out["evidence"]["raw_hits"]["is"] and out["evidence"]["raw_hits"]["target"]


@_blast_or_skip
def test_junction_negative_too_far(tmp_path):
    """ISAba1 ~1 kb upstream of OXA (> 400 bp window) -> signal False."""
    isa, oxa = _seq(_IS_REF), _seq(_OXA_REF)
    contig = ("GC" * 80) + isa + ("AT" * 500) + oxa + ("GC" * 80)   # 1000 bp gap
    g = tmp_path / "far.fna"
    _write_contigs(g, {"contig1": contig})
    out = detect_is_upstream_junction(g, is_ref=_IS_REF, target_ref=_OXA_REF)
    assert out["status"] == "ok" and out["signal"] is False


@_blast_or_skip
def test_different_contig_no_signal(tmp_path):
    """ISAba1 and OXA on DIFFERENT contigs -> signal False (the same-contig test)."""
    isa, oxa = _seq(_IS_REF), _seq(_OXA_REF)
    g = tmp_path / "split.fna"
    _write_contigs(g, {"contig_is": ("GC" * 80) + isa + ("GC" * 80),
                       "contig_oxa": ("AT" * 80) + oxa + ("AT" * 80)})
    out = detect_is_upstream_junction(g, is_ref=_IS_REF, target_ref=_OXA_REF)
    assert out["status"] == "ok" and out["signal"] is False


@_blast_or_skip
def test_no_is_element_no_signal(tmp_path):
    """OXA present, no ISAba1 anywhere -> signal False (and n_is_hits == 0)."""
    oxa = _seq(_OXA_REF)
    g = tmp_path / "oxa_only.fna"
    _write_contigs(g, {"contig1": ("GC" * 200) + oxa + ("GC" * 200)})
    out = detect_is_upstream_junction(g, is_ref=_IS_REF, target_ref=_OXA_REF)
    assert out["status"] == "ok" and out["signal"] is False and out["evidence"]["n_is_hits"] == 0


@_blast_or_skip
def test_multicopy_is_no_truncation(tmp_path):
    """Many ISAba1 copies, only one upstream of OXA -> signal True (proves no -max_target_seqs truncation)."""
    isa, oxa = _seq(_IS_REF), _seq(_OXA_REF)
    decoys = "".join(isa + ("CT" * 300) for _ in range(8))          # 8 ISAba1 copies far from OXA
    contig = decoys + ("GC" * 80) + isa + ("AT" * 100) + oxa + ("GC" * 80)   # the 9th is 200 bp upstream
    g = tmp_path / "multi.fna"
    _write_contigs(g, {"contig1": contig})
    out = detect_is_upstream_junction(g, is_ref=_IS_REF, target_ref=_OXA_REF)
    assert out["status"] == "ok" and out["signal"] is True
    assert out["evidence"]["n_is_hits"] >= 9                         # all copies enumerated, not capped at 5


@_blast_or_skip
def test_primary_ignores_orientation(tmp_path):
    """Default path must NOT enforce IS-orientation: ISAba1 on the opposite strand, still upstream -> True."""
    isa, oxa = _seq(_IS_REF), _seq(_OXA_REF)
    contig = ("GC" * 80) + _rc(isa) + ("AT" * 100) + oxa + ("GC" * 80)   # IS reverse-complemented
    g = tmp_path / "revis.fna"
    _write_contigs(g, {"contig1": contig})
    primary = detect_is_upstream_junction(g, is_ref=_IS_REF, target_ref=_OXA_REF)
    assert primary["status"] == "ok" and primary["signal"] is True       # primary ignores orientation
    refined = detect_is_upstream_junction(g, is_ref=_IS_REF, target_ref=_OXA_REF, is_orientation=True)
    assert refined["evidence"]["is_orientation_enforced"] is True        # refinement path is reachable + flagged


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-q"]))
