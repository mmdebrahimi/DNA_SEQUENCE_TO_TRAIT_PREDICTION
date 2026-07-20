"""Pin the N. gonorrhoeae cipro curated rule (gyrA QRDR Ser91/Asp95 -> R; parC accessory-only)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from dna_decode.organism_rules.neisseria_amr import (  # noqa: E402
    call_ng_amr, call_ng_azithromycin, call_ng_ceftriaxone, call_ng_cefixime,
    call_ng_ciprofloxacin, call_ng_gentamicin, call_ng_penicillin, call_ng_tetracycline,
)


def test_ng_tet_tetM_only():
    # tet(M) (high-level) -> R; rpsJ V57M is accessory-only (low-level, over-calls -> excluded from the call)
    assert call_ng_tetracycline(["tet(M)"])["prediction"] == "R"
    assert call_ng_tetracycline(["tet(M)_1"])["prediction"] == "R"       # startswith tet(M)
    r = call_ng_tetracycline(["rpsJ_V57M"])
    assert r["prediction"] == "S" and r["accessory_rpsJ_V57M"] == ["rpsJ_V57M"]
    assert call_ng_tetracycline([])["prediction"] == "S"
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


def test_ng_penicillin_blaTEM_primary():
    assert call_ng_penicillin(["blaTEM-1"])["prediction"] == "R"
    assert call_ng_penicillin(["blaTEM-135"])["prediction"] == "R"
    # chromosomal penA/mtrR/ponA/porB are accessory-only
    r = call_ng_penicillin(["penA_A501V", "mtrR_G45D", "ponA_L421P"])
    assert r["prediction"] == "S" and r["accessory_chromosomal"]
    assert call_ng_penicillin([])["prediction"] == "S"


def test_ng_esc_penA_mosaic():
    # AMRFinder emits only resistance-curated penA point mutations (all Subclass=CEPHALOSPORIN), incl the
    # real mosaic positions 504/510 an earlier hard-coded codon set missed -> match ANY penA point.
    for sym in ("penA_G545S", "penA_I312M", "penA_V316T", "penA_F504L", "penA_A510V", "penA_N512Y",
                "penA_mosaic_60.001"):
        assert call_ng_ceftriaxone([sym])["prediction"] == "R", sym
        assert call_ng_cefixime([sym])["prediction"] == "R", sym
    assert call_ng_ceftriaxone([])["prediction"] == "S"
    # ponA/porB/mtrR are accessory (do NOT flip the call alone)
    r = call_ng_ceftriaxone(["ponA_L421P", "porB1b_G120K", "mtrR_A-53del"])
    assert r["prediction"] == "S" and r["accessory_ponA_porB_mtr"]


def test_ng_real_amrfinder_symbols_AR0165():
    """R3 real-surface pin: the ACTUAL AMRFinder -O Neisseria_gonorrhoeae output for AR Bank #0165
    (GCA_042036815.1). gyrA S91F/D95G -> cipro R; penA mosaic set -> ESC R; no 23S -> azithro S;
    no blaTEM -> penicillin S (chromosomal accessory only); gentamicin abstains."""
    syms = ["ponA_L421P", "pbp2", "penA_A510V", "penA_F504L", "penA_G545S", "penA_I312M",
            "penA_N512Y", "penA_V316T", "mtrR_A-53del", "folP_R228S", "rpsJ_V57M",
            "porB1b_A121N", "porB1b_G120K", "gyrA_D95G", "gyrA_S91F", "parC_S87R"]
    assert call_ng_amr("ciprofloxacin", syms)["prediction"] == "R"     # gyrA QRDR
    assert call_ng_amr("ceftriaxone", syms)["prediction"] == "R"       # penA mosaic
    assert call_ng_amr("cefixime", syms)["prediction"] == "R"
    assert call_ng_amr("azithromycin", syms)["prediction"] == "S"      # no 23S mutation on this isolate
    assert call_ng_amr("penicillin", syms)["prediction"] == "S"        # no blaTEM (chromosomal accessory)
    assert call_ng_amr("tetracycline", syms)["prediction"] == "S"      # rpsJ V57M accessory-only, no tet(M)
    assert call_ng_amr("gentamicin", syms)["prediction"] == "INDETERMINATE"


def test_ng_gentamicin_abstains():
    r = call_ng_gentamicin(["anything"])
    assert r["prediction"] == "INDETERMINATE"
    assert r["rule_status"] == "ABSTAIN_NO_DETERMINANT"


def test_ng_dispatch():
    assert call_ng_amr("ceftriaxone", ["penA_G545S"])["prediction"] == "R"
    assert call_ng_amr("azithromycin", ["23S_A2045G"])["prediction"] == "R"
    assert call_ng_amr("ciprofloxacin", ["gyrA_S91F"])["prediction"] == "R"
    assert call_ng_amr("gentamicin", [])["prediction"] == "INDETERMINATE"
    assert call_ng_amr("meropenem", [])["prediction"] == "INDETERMINATE"   # unsupported drug -> abstain
    assert call_ng_amr("CEFIXIME", ["penA_A501P"])["prediction"] == "R"     # case-insensitive


def test_ng_new_rules_nonfrozen_scoped():
    for fn, sym in ((call_ng_azithromycin, "23S_A2045G"), (call_ng_penicillin, "blaTEM-1"),
                    (call_ng_ceftriaxone, "penA_G545S"), (call_ng_cefixime, "penA_G545S")):
        r = fn([sym])
        assert r["rule_status"] == "CURATED_NONFROZEN" and r["rule_scope"] == "scorer_local"
