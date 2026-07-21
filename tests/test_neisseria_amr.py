"""Pin the N. gonorrhoeae cipro curated rule (gyrA QRDR Ser91/Asp95 -> R; parC accessory-only)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from dna_decode.organism_rules.neisseria_amr import (  # noqa: E402
    call_ng_amr, call_ng_azithromycin, call_ng_ceftriaxone, call_ng_cefixime,
    call_ng_ciprofloxacin, call_ng_gentamicin, call_ng_penicillin, call_ng_tetracycline,
)


def test_ng_tet_v02_tetM_highlevel_only():
    # v0.2 (NCBI-PD-validated): tet(M) high-level TRNG ONLY. The v0.1 rpsJ/mtrR promotion over-called
    # (spec 0.0 on NCBI-PD: rpsJ_V57M in 34/34 R AND 10/26 S). Narrowed to the clean plasmid marker.
    assert call_ng_tetracycline(["tet(M)"])["prediction"] == "R"
    assert call_ng_tetracycline(["tet(M)_1"])["prediction"] == "R"       # startswith tet(M)
    # v0.2: the chromosomal markers no longer confer R on their own (they over-called)
    assert call_ng_tetracycline(["rpsJ_V57M"])["prediction"] == "S"      # was R in v0.1 -> now S (over-call fix)
    assert call_ng_tetracycline(["mtrR_A-53del"])["prediction"] == "S"   # was R in v0.1 -> now S
    assert call_ng_tetracycline([])["prediction"] == "S"
    assert call_ng_tetracycline(["gyrA_S91F"])["prediction"] == "S"      # unrelated gene -> S
    assert call_ng_tetracycline(["tet(M)", "rpsJ_V57M"])["prediction"] == "R"


def test_gyrA_qrdr_confers_R():
    for sym in ("gyrA_S91F", "gyrA_S91Y", "gyrA_D95N", "gyrA_D95A", "gyrA_D95G"):
        r = call_ng_ciprofloxacin([sym])
        assert r["prediction"] == "R", sym
        assert r["matched_gyrA_qrdr"] == [sym]


def test_no_qrdr_is_S():
    assert call_ng_ciprofloxacin([])["prediction"] == "S"
    assert call_ng_ciprofloxacin(["mtrR", "penA_1", "porB"])["prediction"] == "S"       # irrelevant genes
    assert call_ng_ciprofloxacin(["gyrA_A75S"])["prediction"] == "S"                    # gyrA but NOT a QRDR codon


def test_parC_is_accessory_only():
    r = call_ng_ciprofloxacin(["parC_S87R"])
    assert r["prediction"] == "S"                          # parC alone does NOT flip the binary call
    assert r["accessory_parC_qrdr"] == ["parC_S87R"]
    r2 = call_ng_ciprofloxacin(["gyrA_S91F", "parC_S87R"])
    assert r2["prediction"] == "R" and r2["accessory_parC_qrdr"] == ["parC_S87R"]


def test_rule_is_nonfrozen_scoped():
    r = call_ng_ciprofloxacin(["gyrA_S91F"])
    assert r["rule_status"] == "CURATED_NONFROZEN" and r["rule_scope"] == "scorer_local"


# --- Extension 2026-07-20: azithromycin / ceftriaxone / cefixime / penicillin / gentamicin --------

def test_ng_azithromycin_23S_primary():
    for mut in ("23S_A2045G", "23S_C2597T", "23S_A2059G", "23S_C2611T"):   # WHO + E. coli coords both accepted
        assert call_ng_azithromycin([mut])["prediction"] == "R", mut
    # mtrR efflux is accessory-only (low-level, over-calls) -> does NOT flip the binary call
    r = call_ng_azithromycin(["mtrR_mosaic", "mtr_promoter_a-57del"])
    assert r["prediction"] == "S" and r["accessory_mtr"]
    assert call_ng_azithromycin([])["prediction"] == "S"
    assert call_ng_azithromycin(["23S_C2597T", "mtrR_mosaic"])["prediction"] == "R"


def test_ng_penicillin_v02_specific_determinants():
    # v0.2 (NCBI-PD-validated): blaTEM (PPNG penicillinase) OR ponA_L421P (PBP1). The v0.1 penA-point/mtrR
    # promotion over-called (spec 0.0 on NCBI-PD: those markers are near-universal). Narrowed to specifics.
    assert call_ng_penicillin(["blaTEM-1"])["prediction"] == "R"
    assert call_ng_penicillin(["blaTEM-135"])["prediction"] == "R"
    assert call_ng_penicillin(["ponA_L421P"])["prediction"] == "R"        # chromosomal PBP1 -> R
    # v0.2: penA-point + mtrR no longer confer R on their own (they over-called; lineage-linked not causal)
    assert call_ng_penicillin(["penA_G545S"])["prediction"] == "S"        # was R in v0.1 -> now S (over-call fix)
    assert call_ng_penicillin(["mtrR_A-53del"])["prediction"] == "S"      # was R in v0.1 -> now S
    assert call_ng_penicillin([])["prediction"] == "S"
    assert call_ng_penicillin(["gyrA_S91F"])["prediction"] == "S"         # unrelated -> S


def test_ng_cefixime_v01_mosaic34_core():
    # v0.1 (AR-Bank-validated): cefixime-R requires the mosaic penA-34 CORE {I312M,V316T,N512Y,G545S}
    # (>=3 of 4), NOT any penA ESC point. Fixes v0's spec 0.0 (partial-mosaic S isolates were all FP).
    # R: full mosaic-34 quartet (the AR-Bank R signature)
    R = call_ng_cefixime(["penA_I312M", "penA_V316T", "penA_N512Y", "penA_G545S", "penA_A510V", "penA_F504L"])
    assert R["prediction"] == "R" and len(R["matched_penA_mosaic34_core"]) == 4
    # S: the partial-mosaic reduced-susceptibility signature (A510V/F504L/A516G) -> S (was v0's FP)
    S = call_ng_cefixime(["penA_A510V", "penA_F504L", "penA_A516G", "penA_D346DD"])
    assert S["prediction"] == "S", "partial mosaic must NOT be called cefixime-R (v0 spec-0.0 bug)"
    # a single shared marker alone is NOT enough (< 3 core)
    assert call_ng_cefixime(["penA_A510V"])["prediction"] == "S"
    assert call_ng_cefixime([])["prediction"] == "S"
    acc = call_ng_cefixime(["ponA_L421P", "porB1b_G120K", "mtrR_A-53del"])   # accessory only -> S
    assert acc["prediction"] == "S" and acc["accessory_ponA_porB_mtr"]


def test_ng_ceftriaxone_v01_A501_specific():
    # v0.1: ceftriaxone-R requires the SPECIFIC high-level penA Ala501 marker (A501P/T/V), NOT any penA point.
    for sym in ("penA_A501P", "penA_A501T", "penA_A501V"):
        assert call_ng_ceftriaxone([sym])["prediction"] == "R", sym
    # the reduced-susceptibility mosaic markers (A510V/F504L/G545S/etc) -> ceftriaxone S (they raise cefixime,
    # not ceftriaxone, MIC) -- this is what lifted spec 0.0 -> correct on the AR Bank
    for sym in ("penA_A510V", "penA_F504L", "penA_G545S", "penA_I312M", "penA_V316T", "penA_N512Y"):
        assert call_ng_ceftriaxone([sym])["prediction"] == "S", sym
    assert call_ng_ceftriaxone([])["prediction"] == "S"
    # a real AR#0165-style full-mosaic-minus-A501 vector -> S
    assert call_ng_ceftriaxone(["penA_A510V", "penA_F504L", "penA_G545S", "ponA_L421P", "mtrR_A-53del"])["prediction"] == "S"


def test_ng_real_amrfinder_symbols_AR0165():
    """R3 real-surface pin: the ACTUAL AMRFinder -O Neisseria_gonorrhoeae output for AR Bank #0165
    (GCA_042036815.1). gyrA S91F/D95G -> cipro R; penA mosaic-34 core -> cefixime R; ponA_L421P ->
    penicillin R (v0.2); rpsJ_V57M but NO tet(M) -> tetracycline S (v0.2 over-call fix); gentamicin abstains."""
    syms = ["ponA_L421P", "pbp2", "penA_A510V", "penA_F504L", "penA_G545S", "penA_I312M",
            "penA_N512Y", "penA_V316T", "mtrR_A-53del", "folP_R228S", "rpsJ_V57M",
            "porB1b_A121N", "porB1b_G120K", "gyrA_D95G", "gyrA_S91F", "parC_S87R"]
    assert call_ng_amr("ciprofloxacin", syms)["prediction"] == "R"     # gyrA QRDR
    assert call_ng_amr("cefixime", syms)["prediction"] == "R"          # mosaic-34 core (I312M/V316T/N512Y/G545S)
    assert call_ng_amr("ceftriaxone", syms)["prediction"] == "S"       # v0.1: A510V (not A501) -> ceftriaxone S
    assert call_ng_amr("azithromycin", syms)["prediction"] == "S"      # no 23S mutation on this isolate
    assert call_ng_amr("penicillin", syms)["prediction"] == "R"        # v0.2: ponA_L421P (chromosomal PBP1) -> R
    assert call_ng_amr("tetracycline", syms)["prediction"] == "S"      # v0.2: rpsJ_V57M but no tet(M) -> S
    assert call_ng_amr("gentamicin", syms)["prediction"] == "INDETERMINATE"


def test_ng_gentamicin_abstains():
    r = call_ng_gentamicin(["anything"])
    assert r["prediction"] == "INDETERMINATE"
    assert r["rule_status"] == "ABSTAIN_NO_DETERMINANT"


def test_ng_dispatch():
    assert call_ng_amr("ceftriaxone", ["penA_A501P"])["prediction"] == "R"   # v0.1: A501-class -> ceftriaxone R
    assert call_ng_amr("ceftriaxone", ["penA_G545S"])["prediction"] == "S"   # v0.1: non-A501 mosaic -> S
    assert call_ng_amr("azithromycin", ["23S_A2045G"])["prediction"] == "R"
    assert call_ng_amr("ciprofloxacin", ["gyrA_S91F"])["prediction"] == "R"
    assert call_ng_amr("gentamicin", [])["prediction"] == "INDETERMINATE"
    assert call_ng_amr("meropenem", [])["prediction"] == "INDETERMINATE"   # unsupported drug -> abstain
    # case-insensitive dispatch; v0.1 cefixime needs the mosaic-34 core (>=3 of I312M/V316T/N512Y/G545S)
    assert call_ng_amr("CEFIXIME", ["penA_I312M", "penA_V316T", "penA_N512Y"])["prediction"] == "R"


def test_ng_new_rules_nonfrozen_scoped():
    for fn, sym in ((call_ng_azithromycin, "23S_A2045G"), (call_ng_penicillin, "blaTEM-1"),
                    (call_ng_ceftriaxone, "penA_G545S"), (call_ng_cefixime, "penA_G545S")):
        r = fn([sym])
        assert r["rule_status"] == "CURATED_NONFROZEN" and r["rule_scope"] == "scorer_local"
