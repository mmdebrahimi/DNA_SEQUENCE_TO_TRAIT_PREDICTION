"""G1 machinery test for scripts/build_fungal_cohort.

Builds a 4-isolate cohort from the 3 committed REAL C. auris ERG11 alleles (data/fungal_ref/):
  iso1 = Y132F allele, MIC 64  -> true R, ERG11-R  => TP
  iso2 = K143R allele, MIC 128 -> true R, ERG11-R  => TP
  iso3 = WT allele,    MIC 4   -> true S, ERG11-S  => TN
  iso4 = WT allele,    MIC 256 -> true R, ERG11-S  => FN (simulated efflux/aneuploidy discordance)

Exercises the full pipeline on REAL sequence via REAL blastn: confusion matrix, within-clade buckets,
the efflux-discordance set, and the DOCUMENTED_FAILURE_MODE verdict branch (every FN is efflux-discordant).
Skips cleanly without BLAST+ (offline-safe). Validates MACHINERY, not the real-cohort G1 outcome.

Also unit-tests the MIC breakpoint + MIC-string parser (no BLAST needed).
"""
import shutil
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import pytest  # noqa: E402

from dna_decode.data.fungal_amr import mic_to_phenotype  # noqa: E402
from scripts.build_fungal_cohort import (  # noqa: E402
    _parse_mic,
    build_cohort_report,
)

_FREF_DIR = Path(__file__).resolve().parent.parent / "data" / "fungal_ref"
_REF = _FREF_DIR / "Cauris_ERG11_cds.fna"
_WT = _FREF_DIR / "Cauris_ERG11_PV630306_WT.fna"
_Y132F = _FREF_DIR / "Cauris_ERG11_PV630305_Y132F.fna"
_K143R = _FREF_DIR / "Cauris_ERG11_PV630302_K143R.fna"

_HAS_BLAST = bool(shutil.which("blastn", path="C:/Users/Farshad/ncbi-blast/bin") or shutil.which("blastn")) \
    and bool(shutil.which("makeblastdb", path="C:/Users/Farshad/ncbi-blast/bin") or shutil.which("makeblastdb"))
_HAS_FIXTURES = _REF.exists() and _WT.exists() and _Y132F.exists() and _K143R.exists()


def test_mic_to_phenotype_breakpoint():
    assert mic_to_phenotype("fluconazole", 32) == "R"   # >= 32 = R
    assert mic_to_phenotype("fluconazole", 31.9) == "S"
    assert mic_to_phenotype("fluconazole", 256) == "R"
    assert mic_to_phenotype("caspofungin", 2) == "R"
    assert mic_to_phenotype("amphotericin_b", 99) is None   # no tentative breakpoint configured


def test_parse_mic_tolerates_inequalities():
    assert _parse_mic(">256") == 256.0
    assert _parse_mic("<=1") == 1.0
    assert _parse_mic(" 64.0 ") == 64.0


@pytest.mark.skipif(not _HAS_BLAST, reason="BLAST+ not installed")
@pytest.mark.skipif(not _HAS_FIXTURES, reason="committed fungal_ref fixtures absent")
def test_cohort_pipeline_on_real_alleles():
    table = (
        "isolate_id\tfluconazole_mic\tclade\tgenome_fasta\n"
        f"iso1_Y132F\t64\tI\t{_Y132F}\n"
        f"iso2_K143R\t128\tI\t{_K143R}\n"
        f"iso3_WT\t4\tIII\t{_WT}\n"
        f"iso4_efflux\t256\tIII\t{_WT}\n"   # WT genome, MIC-R => simulated efflux/aneuploidy FN
    )
    with tempfile.TemporaryDirectory() as td:
        tpath = Path(td) / "labels.tsv"
        tpath.write_text(table, encoding="utf-8")
        rep = build_cohort_report(tpath, _REF, Path(td), "fluconazole")

    assert rep.n_total == 4
    assert rep.n_scored == 4
    assert (rep.tp, rep.tn, rep.fp, rep.fn) == (2, 1, 0, 1), vars(rep)
    assert rep.sensitivity == pytest.approx(2 / 3)
    assert rep.specificity == 1.0
    assert rep.accuracy == 0.75
    # the single FN must be the simulated-efflux isolate, surfaced in the discordance set
    assert rep.discordant_efflux == ["iso4_efflux"], rep.discordant_efflux
    # every FN is efflux-discordant => documented-failure-mode branch (not a caller defect)
    assert rep.verdict() == "DOCUMENTED_FAILURE_MODE", rep.verdict()
    # within-clade de-confounding: clade I = 2 TP, clade III = 1 TN + 1 FN
    assert rep.by_clade["I"] == {"tp": 2, "tn": 0, "fp": 0, "fn": 0}
    assert rep.by_clade["III"] == {"tp": 0, "tn": 1, "fp": 0, "fn": 1}


@pytest.mark.skipif(not _HAS_BLAST, reason="BLAST+ not installed")
@pytest.mark.skipif(not _HAS_FIXTURES, reason="committed fungal_ref fixtures absent")
def test_cohort_pass_verdict_all_concordant():
    """3 concordant isolates (no efflux FN) => clean PASS branch."""
    table = (
        "isolate_id\tfluconazole_mic\tclade\tgenome_fasta\n"
        f"a\t64\tI\t{_Y132F}\n"
        f"b\t128\tI\t{_K143R}\n"
        f"c\t4\tIII\t{_WT}\n"
    )
    with tempfile.TemporaryDirectory() as td:
        tpath = Path(td) / "labels.tsv"
        tpath.write_text(table, encoding="utf-8")
        rep = build_cohort_report(tpath, _REF, Path(td), "fluconazole")
    assert (rep.tp, rep.tn, rep.fp, rep.fn) == (2, 1, 0, 0)
    assert rep.accuracy == 1.0 and rep.sensitivity == 1.0
    assert rep.verdict() == "PASS"


@pytest.mark.skipif(not _HAS_BLAST, reason="BLAST+ not installed")
@pytest.mark.skipif(not _HAS_FIXTURES, reason="committed fungal_ref fixtures absent")
def test_cohort_label_limited_failure_verdict():
    """High sens + low spec where FP isolates CARRY the determinant at reduced-susceptibility MICs
    (>= breakpoint/4) => LABEL_LIMITED_FAILURE (the 'suspect the label' mode), not bare FAIL."""
    table = (
        "isolate_id\tfluconazole_mic\tclade\tgenome_fasta\n"
        f"tp\t64\tI\t{_K143R}\n"          # R, K143R -> TP
        f"fp1\t16\tIII\t{_Y132F}\n"       # true S (MIC16) but carries Y132F -> FP, reduced-suscept
        f"fp2\t16\tIII\t{_K143R}\n"       # true S (MIC16) but carries K143R -> FP, reduced-suscept
        f"tn\t4\tIII\t{_WT}\n"            # true S, WT -> TN
    )
    with tempfile.TemporaryDirectory() as td:
        tpath = Path(td) / "labels.tsv"
        tpath.write_text(table, encoding="utf-8")
        rep = build_cohort_report(tpath, _REF, Path(td), "fluconazole")
    assert rep.sensitivity == 1.0 and rep.specificity is not None and rep.specificity < 0.5, vars(rep)
    assert (rep.tp, rep.fp, rep.tn) == (1, 2, 1), vars(rep)
    assert rep.verdict() == "LABEL_LIMITED_FAILURE", rep.verdict()


if __name__ == "__main__":
    test_mic_to_phenotype_breakpoint()
    test_parse_mic_tolerates_inequalities()
    print("PASS unit tests (breakpoint + parser)")
    if _HAS_BLAST and _HAS_FIXTURES:
        test_cohort_pipeline_on_real_alleles()
        test_cohort_pass_verdict_all_concordant()
        print("PASS cohort-pipeline tests (real alleles)")
    else:
        print("SKIP cohort-pipeline (no BLAST or fixtures)")
