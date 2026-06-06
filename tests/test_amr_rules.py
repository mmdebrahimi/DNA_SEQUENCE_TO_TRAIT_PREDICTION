"""Tests for the deterministic AMR caller (dna_decode/eval/amr_rules).

Pure-logic tests on synthetic main.tsv files + (when the committed AMRFinder cache is present) a cohort
op-char regression pin. Runnable via pytest OR standalone. No Docker / no AMRFinder run.
"""
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from dna_decode.eval.amr_rules import DRUG_RULE, call_resistance, evaluate_cohort, rule_for

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


def test_per_drug_rule_config():
    # the validated per-drug config: cipro needs 2 (QRDR), cef/tet need 1; cef has subclass refinement.
    assert rule_for("ciprofloxacin")["threshold"] == 2
    assert rule_for("ceftriaxone")["threshold"] == 1
    assert rule_for("tetracycline")["threshold"] == 1
    assert rule_for("ceftriaxone")["subclass_any"] is not None
    assert rule_for("ciprofloxacin")["subclass_any"] is None
    # unconfigured drug falls back to threshold 2, no refinement
    assert rule_for("nonexistent-drug")["threshold"] == 2


def test_cef_extended_spectrum_excludes_narrow_beta_lactam():
    # blaTEM-1 (Subclass BETA-LACTAM) is ampicillin-R, NOT ceftriaxone-R → must NOT be counted for cef.
    with tempfile.TemporaryDirectory() as td:
        m = _write_main(Path(td), [("blaTEM-1", "BETA-LACTAM", "BETA-LACTAM")])
        c = call_resistance(m, "ceftriaxone")
    assert c["prediction"] == "S" and c["n_determinants"] == 0


def test_cef_extended_spectrum_counts_esbl_at_threshold_1():
    # blaCTX-M-15 (Subclass CEPHALOSPORIN) is a real 3rd-gen-cephalosporin ESBL → cef R at threshold 1.
    with tempfile.TemporaryDirectory() as td:
        m = _write_main(Path(td), [("blaCTX-M-15", "BETA-LACTAM", "CEPHALOSPORIN"),
                                   ("blaTEM-1", "BETA-LACTAM", "BETA-LACTAM")])
        c = call_resistance(m, "ceftriaxone")
    assert c["prediction"] == "R" and c["n_determinants"] == 1
    assert c["determinants"][0]["symbol"] == "blaCTX-M-15"   # narrow blaTEM excluded from report


def test_cef_carbapenemase_counts():
    with tempfile.TemporaryDirectory() as td:
        m = _write_main(Path(td), [("blaNDM-5", "BETA-LACTAM", "CARBAPENEM")])
        c = call_resistance(m, "ceftriaxone")
    assert c["prediction"] == "R" and c["n_determinants"] == 1


def test_tet_single_acquired_gene_is_resistant():
    # tetracycline is acquired-gene (one tet(A) = R) → threshold 1 by default.
    with tempfile.TemporaryDirectory() as td:
        m = _write_main(Path(td), [("tet(A)", "TETRACYCLINE", "TETRACYCLINE")])
        c = call_resistance(m, "tetracycline")
    assert c["prediction"] == "R" and c["n_determinants"] == 1


def test_gent_subclass_excludes_non_gentamicin_aminoglycoside():
    # aph(3')-Ia (Subclass KANAMYCIN/NEOMYCIN) + aadA (STREPTOMYCIN) confer aminoglycoside-R but NOT
    # gentamicin-R → must NOT be counted for gentamicin.
    with tempfile.TemporaryDirectory() as td:
        m = _write_main(Path(td), [("aph(3')-Ia", "AMINOGLYCOSIDE", "KANAMYCIN/NEOMYCIN"),
                                   ("aadA1", "AMINOGLYCOSIDE", "STREPTOMYCIN")])
        c = call_resistance(m, "gentamicin")
    assert c["prediction"] == "S" and c["n_determinants"] == 0


def test_gent_subclass_counts_aac3_at_threshold_1():
    # aac(3)-IIa (Subclass GENTAMICIN) is the gentamicin-conferring family → R at threshold 1.
    with tempfile.TemporaryDirectory() as td:
        m = _write_main(Path(td), [("aac(3)-IIa", "AMINOGLYCOSIDE", "GENTAMICIN"),
                                   ("aph(3')-Ia", "AMINOGLYCOSIDE", "KANAMYCIN/NEOMYCIN")])
        c = call_resistance(m, "gentamicin")
    assert c["prediction"] == "R" and c["n_determinants"] == 1
    assert c["determinants"][0]["symbol"] == "aac(3)-IIa"


def _gent_cohort_pairs(root):
    from dna_decode.data.cohort import load_cohort
    import glob
    runs = root / "data/amrfinder_runs"
    pooled = {}
    for cp in glob.glob(str(root / "data/processed/*.parquet")):
        try:
            c = load_cohort(cp)
        except Exception:
            continue
        for s in c.strains:
            k = [k for k in s.ast_labels if k.lower() == "gentamicin"]
            if k and (runs / s.assembly_accession / "main.tsv").exists():
                pooled[s.assembly_accession] = int(s.ast_labels[k[0]])
    return list(pooled.items())


def test_gent_cohort_opchars_regression():
    """Pin pooled gentamicin op-chars (N=128: acc 0.945 with GENTAMICIN-subclass refinement).
    Skips cleanly if data/ absent."""
    root = Path(__file__).resolve().parent.parent
    if not (root / "data/amrfinder_runs").exists():
        print("SKIP gent cohort op-chars (data/ not present)")
        return
    pairs = _gent_cohort_pairs(root)
    if len(pairs) < 100:
        print(f"SKIP gent cohort op-chars (only {len(pairs)} cached)")
        return
    r = evaluate_cohort(root / "data/amrfinder_runs", pairs, "gentamicin")
    assert r["accuracy"] >= 0.90, r        # validated 0.945
    assert r["specificity"] >= 0.90, r     # validated 0.96 (the refinement's payoff)
    print(f"gent cohort op-chars OK: {r}")


def _cef_cohort_pairs(root):
    from dna_decode.data.cohort import load_cohort
    cohort_p = root / "data/processed/gate_b_cohort.parquet"
    if not cohort_p.exists():
        return None
    c = load_cohort(str(cohort_p))
    return [(s.assembly_accession, int(s.ast_labels[k[0]]))
            for s in c.strains
            for k in [[k for k in s.ast_labels if k.lower() == "ceftriaxone"]] if k]


def test_cef_cohort_opchars_regression():
    """Pin the cef extended-spectrum deterministic-rule op-chars on gate_b_cohort (N=60: acc 0.933).
    Skips cleanly if data/ is absent (gitignored)."""
    root = Path(__file__).resolve().parent.parent
    runs = root / "data/amrfinder_runs"
    if not runs.exists():
        print("SKIP cef cohort op-chars (data/ not present)")
        return
    pairs = _cef_cohort_pairs(root)
    if not pairs:
        print("SKIP cef cohort op-chars (cohort absent)")
        return
    r = evaluate_cohort(runs, pairs, "ceftriaxone")
    if r["n"] < 50:
        print(f"SKIP cef cohort op-chars (only {r['n']} cached)")
        return
    assert r["accuracy"] >= 0.88, r        # validated 0.933
    assert r["sensitivity"] >= 0.90, r     # validated 0.962
    assert r["specificity"] >= 0.85, r     # validated 0.912
    print(f"cef cohort op-chars OK: {r}")


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
