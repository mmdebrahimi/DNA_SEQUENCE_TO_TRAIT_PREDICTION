"""Tests for the deterministic AMR caller (dna_decode/eval/amr_rules).

Pure-logic tests on synthetic main.tsv files + (when the committed AMRFinder cache is present) a cohort
op-char regression pin. Runnable via pytest OR standalone. No Docker / no AMRFinder run.
"""
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from dna_decode.eval.amr_rules import (
    FN_UNDETECTED_MECHANISM, FP_DETERMINANT_NO_PHENOTYPE, UNDETECTABLE_MECHANISMS,
    DRUG_RULE, call_resistance, discordance_bucket, evaluate_cohort, qrdr_point_count, rule_for,
)

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


def test_single_qrdr_mutation_is_susceptible_default_threshold():
    # cipro counter='qrdr_point', threshold=2: a single QRDR point mutation → S (boundary), not R.
    with tempfile.TemporaryDirectory() as td:
        m = _write_main_with_method(Path(td), [("gyrA_S83L", "QUINOLONE", "QUINOLONE", "POINTX")])
        c = call_resistance(m, "ciprofloxacin")
    assert c["prediction"] == "S" and c["n_determinants"] == 1 and c["confidence"] == "MODERATE"


def test_resistant_on_two_qrdr_mutations():
    # gyrA + parC point mutations (the canonical clinical FQ-R combo) → R at threshold 2.
    with tempfile.TemporaryDirectory() as td:
        m = _write_main_with_method(Path(td), [("gyrA_S83L", "QUINOLONE", "QUINOLONE", "POINTX"),
                                               ("parC_S80I", "QUINOLONE", "QUINOLONE", "POINTX")])
        c = call_resistance(m, "ciprofloxacin")
    assert c["prediction"] == "R" and c["n_determinants"] == 2


def test_cipro_qrdr_excludes_intrinsic_efflux_and_acquired():
    # the cross-organism fix: intrinsic oqxAB efflux + acquired qnr are NOT counted (only QRDR POINT).
    with tempfile.TemporaryDirectory() as td:
        m = _write_main_with_method(Path(td), [
            ("gyrA_S83L", "QUINOLONE", "QUINOLONE", "POINTX"),       # QRDR → counts
            ("oqxB", "QUINOLONE", "NITROFURAN/PHENICOL/QUINOLONE", "EXACTX"),  # intrinsic efflux → excluded
            ("qnrS1", "QUINOLONE", "QUINOLONE", "EXACTX"),           # acquired gene → excluded
        ])
        c = call_resistance(m, "ciprofloxacin")
    # only 1 QRDR point mutation → below threshold 2 → S (broad rule would have called R on 3 determinants)
    assert c["prediction"] == "S" and c["n_determinants"] == 1


def test_threshold_1_override_still_works():
    # explicit threshold override still applies (single QRDR mutation + threshold=1 → R).
    with tempfile.TemporaryDirectory() as td:
        m = _write_main_with_method(Path(td), [("gyrA_S83L", "QUINOLONE", "QUINOLONE", "POINTX")])
        c = call_resistance(m, "ciprofloxacin", resistance_threshold=1)
    assert c["prediction"] == "R"


def test_broad_matcher_still_multiclass(check=None):
    # the broad determinant matcher (used by naive baseline / cross-source) still does multiclass matching;
    # cipro's CALL path no longer uses it (qrdr_point), but the matcher itself is unchanged.
    from dna_decode.eval.amr_rules import cipro_determinants_from_main
    with tempfile.TemporaryDirectory() as td:
        m = _write_main(Path(td), [("aac(6')-Ib-cr", "AMINOGLYCOSIDE", "AMIKACIN/KANAMYCIN/QUINOLONE/TOBRAMYCIN")])
        dets = cipro_determinants_from_main(m, "ciprofloxacin")
    assert len(dets) == 1 and dets[0]["symbol"] == "aac(6')-Ib-cr"


def test_susceptible_when_no_quinolone_determinant():
    with tempfile.TemporaryDirectory() as td:
        m = _write_main(Path(td), [("blaTEM-1", "BETA-LACTAM", "BETA-LACTAM"),
                                   ("tet(A)", "TETRACYCLINE", "TETRACYCLINE")])
        c = call_resistance(m, "ciprofloxacin")
    assert c["prediction"] == "S" and c["n_determinants"] == 0 and c["confidence"] == "HIGH"


def test_indeterminate_when_main_missing():
    c = call_resistance(Path("/nonexistent/main.tsv"), "ciprofloxacin")
    assert c["prediction"] == "INDETERMINATE"


def test_susceptible_call_carries_blind_spots():
    # an S call must surface the mechanisms the curated-determinant rule cannot see (honest negative).
    with tempfile.TemporaryDirectory() as td:
        m = _write_main(Path(td), [("blaTEM-1", "BETA-LACTAM", "BETA-LACTAM")])  # not cef-relevant
        c = call_resistance(m, "ceftriaxone")
    assert c["prediction"] == "S"
    assert c["undetectable_mechanisms"] == UNDETECTABLE_MECHANISMS  # efflux/porin_loss/regulatory
    assert "rule out resistance" in c["caveat"]


def test_resistant_call_has_no_blind_spot_list():
    with tempfile.TemporaryDirectory() as td:
        m = _write_main(Path(td), [("blaCTX-M-15", "BETA-LACTAM", "CEPHALOSPORIN")])
        c = call_resistance(m, "ceftriaxone")
    assert c["prediction"] == "R" and c["undetectable_mechanisms"] == []


def test_discordance_bucket_taxonomy():
    assert discordance_bucket("S", 1) == FN_UNDETECTED_MECHANISM        # R phenotype missed
    assert discordance_bucket("R", 0) == FP_DETERMINANT_NO_PHENOTYPE    # called R, susceptible
    assert discordance_bucket("R", 1) is None                          # concordant
    assert discordance_bucket("S", 0) is None                          # concordant
    assert discordance_bucket("INDETERMINATE", 1) is None


def test_evaluate_cohort_emits_discordance_breakdown():
    with tempfile.TemporaryDirectory() as td:
        # one true-R strain with NO cef determinant (FN) + one true-S strain WITH an ESBL (FP)
        import os
        fn_dir = Path(td) / "FN"; fn_dir.mkdir()
        _write_main(fn_dir, [("blaTEM-1", "BETA-LACTAM", "BETA-LACTAM")])  # no extended-spectrum → S
        fp_dir = Path(td) / "FP"; fp_dir.mkdir()
        _write_main(fp_dir, [("blaCTX-M-15", "BETA-LACTAM", "CEPHALOSPORIN")])  # ESBL → R
        os.rename(fn_dir / "main.tsv", fn_dir / "main.tsv")  # no-op; main.tsv already there
        r = evaluate_cohort(Path(td), [("FN", 1), ("FP", 0)], "ceftriaxone")
    assert r["discordance"][FN_UNDETECTED_MECHANISM] == 1
    assert r["discordance"][FP_DETERMINANT_NO_PHENOTYPE] == 1


def _write_main_with_method(tmp, rows):
    """rows = [(element_symbol, class, subclass, method)]. Sets the Method column (index 12)."""
    lines = [_HEADER]
    for sym, cls, sub, meth in rows:
        cells = [""] * 22
        cells[5] = sym; cells[10] = cls; cells[11] = sub; cells[12] = meth
        lines.append("\t".join(cells))
    p = tmp / "main.tsv"
    p.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return p


def test_qrdr_point_count_counts_target_mutations_only():
    with tempfile.TemporaryDirectory() as td:
        m = _write_main_with_method(Path(td), [
            ("gyrA_S83L", "QUINOLONE", "QUINOLONE", "POINTX"),       # QRDR point → counts
            ("parC_S80I", "QUINOLONE", "QUINOLONE", "POINTX"),       # QRDR point → counts
            ("oqxB", "QUINOLONE", "NITROFURAN/PHENICOL/QUINOLONE", "EXACTX"),  # intrinsic efflux → excluded
            ("qnrS1", "QUINOLONE", "QUINOLONE", "EXACTX"),           # acquired gene (not POINT) → excluded
        ])
        assert qrdr_point_count(m) == 2   # only gyrA + parC point mutations


def test_qrdr_point_count_excludes_intrinsic_efflux():
    # a susceptible-shaped Klebsiella strain: intrinsic oqxAB only, no QRDR point mutation → 0.
    with tempfile.TemporaryDirectory() as td:
        m = _write_main_with_method(Path(td), [
            ("oqxA", "QUINOLONE", "NITROFURAN/PHENICOL/QUINOLONE", "EXACTX"),
            ("oqxB", "QUINOLONE", "NITROFURAN/PHENICOL/QUINOLONE", "EXACTX"),
        ])
        assert qrdr_point_count(m) == 0


def test_qrdr_point_count_missing_file():
    assert qrdr_point_count(Path("/nonexistent/main.tsv")) == 0


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


def test_meropenem_carbapenemase_counts_esbl_excluded():
    # meropenem rule: a CARBAPENEM-subclass carbapenemase (blaKPC/NDM/OXA-48) → R; ESBL (CEPHALOSPORIN)
    # raises MIC but doesn't hydrolyze carbapenems → excluded.
    with tempfile.TemporaryDirectory() as td:
        m = _write_main(Path(td), [("blaKPC-2", "BETA-LACTAM", "CARBAPENEM"),
                                   ("blaCTX-M-15", "BETA-LACTAM", "CEPHALOSPORIN")])
        c = call_resistance(m, "meropenem")
    assert c["prediction"] == "R" and c["n_determinants"] == 1
    assert c["determinants"][0]["symbol"] == "blaKPC-2"


def test_meropenem_esbl_only_is_susceptible():
    with tempfile.TemporaryDirectory() as td:
        m = _write_main(Path(td), [("blaCTX-M-15", "BETA-LACTAM", "CEPHALOSPORIN")])  # ESBL, no carbapenemase
        c = call_resistance(m, "meropenem")
    assert c["prediction"] == "S" and c["n_determinants"] == 0


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
    # QRDR-POINT rule (cross-organism, global as of 2026-06-07): E.coli N=147 acc 0.925/sens 0.875/spec 0.973.
    assert r["accuracy"] >= 0.90, r
    assert r["sensitivity"] >= 0.85, r
    assert r["specificity"] >= 0.90, r
    print(f"cohort op-chars OK: {r}")


if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_") and callable(v)]
    for fn in fns:
        fn(); print(f"PASS {fn.__name__}")
    print(f"\n{len(fns)} passed")
