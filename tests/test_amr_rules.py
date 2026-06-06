"""Tests for the deterministic AMR caller (dna_decode/eval/amr_rules).

Pure-logic tests on synthetic main.tsv files + (when the committed AMRFinder cache is present) a cohort
op-char regression pin. Runnable via pytest OR standalone. No Docker / no AMRFinder run.
"""
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from dna_decode.eval.amr_rules import call_resistance, evaluate_cohort

_HEADER = ("Protein id\tContig id\tStart\tStop\tStrand\tElement symbol\tElement name\tScope\tType\t"
           "Subtype\tClass\tSubclass\tMethod\tTarget length\tReference sequence length\t"
           "% Coverage of reference\t% Identity to reference\tAlignment length\tClosest reference accession\t"
           "Closest reference name\tHMM accession\tHMM description")


def _write_main(tmp: Path, rows: list[tuple[str, str, str]]) -> Path:
    """rows = [(element_symbol, class, subclass)]. Returns the main.tsv path."""
    lines = [_HEADER]
    for sym, cls, sub in rows:
        cells = [""] * 22
        cells[5] = sym; cells[10] = cls; cells[11] = sub
        lines.append("\t".join(cells))
    p = tmp / "main.tsv"
    p.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return p


def test_single_determinant_is_susceptible_default_threshold():
    # default threshold=2: a single quinolone determinant → S (boundary), not R (the v1 over-call fix).
    with tempfile.TemporaryDirectory() as td:
        m = _write_main(Path(td), [("qnrS1", "QUINOLONE", "QUINOLONE")])
        c = call_resistance(m, "ciprofloxacin")
    assert c["prediction"] == "S" and c["n_determinants"] == 1 and c["confidence"] == "MODERATE"


def test_resistant_on_two_determinants():
    with tempfile.TemporaryDirectory() as td:
        m = _write_main(Path(td), [("gyrA_S83L", "QUINOLONE", "QUINOLONE"),
                                   ("parC_S80I", "QUINOLONE", "QUINOLONE")])
        c = call_resistance(m, "ciprofloxacin")
    assert c["prediction"] == "R" and c["n_determinants"] == 2


def test_threshold_1_for_acquired_gene_drugs():
    # acquired-gene drugs (one determinant = resistance) override threshold=1.
    with tempfile.TemporaryDirectory() as td:
        m = _write_main(Path(td), [("qnrS1", "QUINOLONE", "QUINOLONE")])
        c = call_resistance(m, "ciprofloxacin", resistance_threshold=1)
    assert c["prediction"] == "R"


def test_multiclass_subclass_matches():
    # AMRFinder emits '/'-joined multi-class strings; a quinolone token anywhere must be COUNTED as a
    # drug-relevant determinant (matching logic, independent of the R/S threshold).
    with tempfile.TemporaryDirectory() as td:
        m = _write_main(Path(td), [("aac(6')-Ib-cr", "AMINOGLYCOSIDE", "AMIKACIN/KANAMYCIN/QUINOLONE/TOBRAMYCIN")])
        c = call_resistance(m, "ciprofloxacin")
    assert c["n_determinants"] == 1 and c["determinants"][0]["symbol"] == "aac(6')-Ib-cr"


def test_susceptible_when_no_quinolone_determinant():
    with tempfile.TemporaryDirectory() as td:
        m = _write_main(Path(td), [("blaTEM-1", "BETA-LACTAM", "BETA-LACTAM"),
                                   ("tet(A)", "TETRACYCLINE", "TETRACYCLINE")])
        c = call_resistance(m, "ciprofloxacin")
    assert c["prediction"] == "S" and c["n_determinants"] == 0 and c["confidence"] == "HIGH"


def test_indeterminate_when_main_missing():
    c = call_resistance(Path("/nonexistent/main.tsv"), "ciprofloxacin")
    assert c["prediction"] == "INDETERMINATE"


def test_cohort_opchars_regression():
    """If the committed AMRFinder cache + cohort are present, pin the cipro deterministic-rule op-chars
    (acc 0.85 / sens 0.96 / spec 0.75 on N=147). Skips cleanly if data/ is absent (gitignored)."""
    root = Path(__file__).resolve().parent.parent
    runs = root / "data/amrfinder_runs"
    cohort_p = root / "data/processed/stage2_n150_cipro_cohort.parquet"
    if not runs.exists() or not cohort_p.exists():
        print("SKIP cohort op-chars (data/ not present — gitignored)")
        return
    from dna_decode.data.cohort import load_cohort
    c = load_cohort(str(cohort_p))
    pairs = [(s.assembly_accession, int(s.ast_labels["ciprofloxacin"]))
             for s in c.strains if "ciprofloxacin" in s.ast_labels]
    r = evaluate_cohort(runs, pairs, "ciprofloxacin")
    if r["n"] < 140:   # cache may be partial on some machines
        print(f"SKIP cohort op-chars (only {r['n']} cached)")
        return
    # tiered rule (threshold=2): acc 0.939 / sens 0.931 / spec 0.947 on cipro N=147.
    assert r["accuracy"] >= 0.90, r
    assert r["sensitivity"] >= 0.88, r
    assert r["specificity"] >= 0.90, r
    print(f"cohort op-chars OK: {r}")


if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_") and callable(v)]
    for fn in fns:
        fn(); print(f"PASS {fn.__name__}")
    print(f"\n{len(fns)} passed")
