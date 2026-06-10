"""Regression tests for the gated EXPRESSION_FLOOR -> R expression-context override in call_resistance.

Pins the default-OFF contract (the shipped registry must not change Acinetobacter|meropenem behavior) and
the opt-in override behavior via a test registry. The override fires ONLY with: explicit genome_fasta +
enabled:true registry block + a detected junction. Junction cases need real BLAST (skipif).
"""
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import pytest  # noqa: E402

from dna_decode.eval.amr_rules import call_resistance, load_calibrated_registry  # noqa: E402

_REF_DIR = Path(__file__).resolve().parent.parent / "data" / "isaba1_ref"
_IS_REF = _REF_DIR / "ISAba1_ref.fna"
_OXA_REF = _REF_DIR / "OXA51fam_ref.fna"
_HAS_BLAST = bool(shutil.which("blastn", path="C:/Users/Farshad/ncbi-blast/bin") or shutil.which("blastn")) \
    and bool(shutil.which("makeblastdb", path="C:/Users/Farshad/ncbi-blast/bin") or shutil.which("makeblastdb"))
_blast_or_skip = pytest.mark.skipif(not (_HAS_BLAST and _IS_REF.exists() and _OXA_REF.exists()),
                                    reason="BLAST+ or ISAba1/OXA refs absent")


def _ec_block(enabled: bool) -> dict:
    return {"enabled": enabled, "experimental": True, "is_ref": str(_IS_REF), "target_ref": str(_OXA_REF),
            "upstream_bp": 400}


def _registry(enabled: bool) -> dict:
    return {"rules": {"Acinetobacter|meropenem": {
        "verdict": "EXPRESSION_FLOOR", "counter": "broad", "threshold": 1,
        "intrinsic_families_excluded": ["blaADC", "blaOXA-51-family"], "loo_balanced_accuracy": 0.13,
        "n": 30, "n_R": 15, "n_S": 15, "expression_context": _ec_block(enabled)}}}


def _main_tsv(tmp_path: Path) -> Path:
    """Minimal AMRFinder main.tsv — the EXPRESSION_FLOOR branch doesn't read determinants, but the file
    must exist (call_resistance returns INDETERMINATE on a missing main.tsv)."""
    p = tmp_path / "main.tsv"
    p.write_text("Protein identifier\tGene symbol\tClass\tSubclass\n", encoding="utf-8")
    return p


def _seq(fasta: Path) -> str:
    return "".join(l.strip() for l in fasta.read_text().splitlines() if not l.startswith(">")).upper()


def _junction_genome(tmp_path: Path) -> Path:
    g = tmp_path / "pos.fna"
    g.write_text(f">c1\n{'GC'*80}{_seq(_IS_REF)}{'AT'*100}{_seq(_OXA_REF)}{'GC'*80}\n", encoding="utf-8")
    return g


def _no_junction_genome(tmp_path: Path) -> Path:
    g = tmp_path / "neg.fna"
    g.write_text(f">c1\n{'GC'*200}{_seq(_OXA_REF)}{'GC'*200}\n", encoding="utf-8")   # OXA only, no ISAba1
    return g


def test_shipped_registry_acinetobacter_is_disabled():
    """The SHIPPED registry must keep Acinetobacter|meropenem expression_context OFF (default-off contract)."""
    reg = load_calibrated_registry()
    ec = reg["rules"]["Acinetobacter|meropenem"].get("expression_context", {})
    assert ec.get("enabled") is False and ec.get("experimental") is True


def test_default_off_genome_present_stays_abstain(tmp_path):
    """enabled:false + genome present -> ABSTAIN unchanged (no override)."""
    g = _no_junction_genome(tmp_path) if False else (tmp_path / "any.fna")
    g.write_text(">c1\nACGTACGTACGT\n", encoding="utf-8")
    out = call_resistance(_main_tsv(tmp_path), "meropenem", organism="Acinetobacter_baumannii",
                          registry=_registry(enabled=False), genome_fasta=g)
    assert out["prediction"] == "ABSTAIN" and "expression_context" not in out


def test_genome_absent_stays_abstain(tmp_path):
    """enabled:true but NO genome supplied -> ABSTAIN (override needs an explicit genome path)."""
    out = call_resistance(_main_tsv(tmp_path), "meropenem", organism="Acinetobacter_baumannii",
                          registry=_registry(enabled=True), genome_fasta=None)
    assert out["prediction"] == "ABSTAIN" and "expression_context" not in out


@_blast_or_skip
def test_enabled_junction_positive_upgrades_to_R(tmp_path):
    """enabled:true + junction-positive genome -> R via expression_context_v1."""
    out = call_resistance(_main_tsv(tmp_path), "meropenem", organism="Acinetobacter_baumannii",
                          registry=_registry(enabled=True), genome_fasta=_junction_genome(tmp_path))
    assert out["prediction"] == "R"
    assert out["rule"].startswith("expression_context_v1")
    assert out["expression_context"]["junction"] is not None


@_blast_or_skip
def test_enabled_no_junction_stays_abstain(tmp_path):
    """enabled:true + genome WITHOUT a junction -> ABSTAIN (no S over-claim)."""
    out = call_resistance(_main_tsv(tmp_path), "meropenem", organism="Acinetobacter_baumannii",
                          registry=_registry(enabled=True), genome_fasta=_no_junction_genome(tmp_path))
    assert out["prediction"] == "ABSTAIN"


def test_calibrated_organism_unaffected_by_genome_param(tmp_path):
    """A CALIBRATED organism (Campylobacter) path is unaffected by the new genome_fasta param."""
    reg = {"rules": {"Campylobacter|ciprofloxacin": {
        "verdict": "CALIBRATED", "counter": "qrdr_point", "threshold": 1,
        "intrinsic_families_excluded": [], "loo_balanced_accuracy": 1.0, "n": 30, "n_R": 15, "n_S": 15}}}
    g = tmp_path / "x.fna"; g.write_text(">c1\nACGTACGT\n", encoding="utf-8")
    out = call_resistance(_main_tsv(tmp_path), "ciprofloxacin", organism="Campylobacter",
                          registry=reg, genome_fasta=g)
    assert out["prediction"] in ("R", "S") and "expression_context" not in out


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-q"]))
