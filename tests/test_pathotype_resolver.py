"""Decision-table tests for the v0 pathotype compatibility resolver.

Exhaustive over the 11-class surface + abstention, on synthetic cluster profiles
(no genome, no I/O). Runnable via pytest OR standalone (`python tests/test_pathotype_resolver.py`)
so it satisfies a test-exit-0 gate even without pytest installed.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from dna_decode.pathotype.resolve import resolve_call


def p(*clusters):
    return {c: True for c in clusters}


def test_ehec():
    r = resolve_call(p("STX2", "LEE"))
    assert r["primary"] == "EHEC_COMPATIBLE" and r["confidence_tier"] == "CONFIDENT"
    assert r["external_validity"] == "scope_limited"  # clean EHEC is scope-limited in v0


def test_stec_non_lee():
    assert resolve_call(p("STX1"))["primary"] == "STEC_NON_LEE"


def test_tepec_supported():
    r = resolve_call(p("LEE", "BFP_EAF"))
    assert r["primary"] == "tEPEC_COMPATIBLE" and r["external_validity"] == "supported"


def test_aepec():
    assert resolve_call(p("LEE"))["primary"] == "aEPEC_COMPATIBLE"


def test_etec_lt_and_st():
    assert resolve_call(p("LT"))["primary"] == "ETEC_COMPATIBLE"
    assert resolve_call(p("ST"))["primary"] == "ETEC_COMPATIBLE"


def test_eaec_confident():
    assert resolve_call(p("EAEC_REG", "AAF_I"))["primary"] == "EAEC_COMPATIBLE"
    assert resolve_call(p("AAF_II", "EAEC_TRANSPORT"))["primary"] == "EAEC_COMPATIBLE"


def test_eaec_ambiguous_regulator_alone():
    r = resolve_call(p("EAEC_REG"))
    assert r["primary"] == "AMBIGUOUS" and "EAEC" in r["reason"]


def test_upec():
    r = resolve_call(p("P_FIMBRIAE", "HEMOLYSIN"))
    assert r["primary"] == "UPEC_COMPATIBLE" and r["external_validity"] == "supported"


def test_upec_needs_two_strong():
    # single strong ExPEC marker -> AMBIGUOUS, not UPEC
    assert resolve_call(p("P_FIMBRIAE"))["primary"] == "AMBIGUOUS"


def test_hybrid_two_dec_modules():
    r = resolve_call(p("LEE", "LT"))  # EPEC + ETEC
    assert r["primary"] == "HYBRID"
    assert set(r["secondary"]) == {"aEPEC_COMPATIBLE", "ETEC_COMPATIBLE"}


def test_hybrid_stx_and_etec():
    r = resolve_call(p("STX2", "LEE", "ST"))  # EHEC + ETEC modules
    assert r["primary"] == "HYBRID"


def test_unclassified_eiec():
    assert resolve_call(p("EIEC_FLAG"))["primary"] == "UNCLASSIFIED"


def test_unclassified_daec_afa_only():
    # afa/dra as the ONLY ExPEC evidence -> DAEC flag -> UNCLASSIFIED (per contract)
    assert resolve_call(p("AFA_DRA"))["primary"] == "UNCLASSIFIED"


def test_commensal_low_marker_burden():
    assert resolve_call({})["primary"] == "COMMENSAL_LOW_MARKER_BURDEN"


def test_low_qc_blocks_commensal():
    r = resolve_call({}, qc_pass=False)
    assert r["primary"] == "AMBIGUOUS_LOW_QC" and r["confidence_tier"] == "AMBIGUOUS"


def test_partial_primary_is_ambiguous():
    r = resolve_call({}, partial_clusters=frozenset({"LT"}))
    assert r["primary"] == "AMBIGUOUS"


def test_epec_with_expec_secondary():
    r = resolve_call(p("LEE", "P_FIMBRIAE", "HEMOLYSIN"))  # EPEC primary + UPEC secondary
    assert r["primary"] == "aEPEC_COMPATIBLE" and r["secondary"] == ["UPEC_COMPATIBLE"]


def test_every_call_has_provenance_fields():
    for prof in [p("STX2", "LEE"), p("LEE"), p("LT"), p("P_FIMBRIAE", "HEMOLYSIN"), {}]:
        r = resolve_call(prof)
        for k in ("primary", "secondary", "confidence_tier", "rule_id", "rule_version",
                  "reason", "external_validity"):
            assert k in r, f"missing {k}"


def _run_standalone() -> int:
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_") and callable(v)]
    failed = 0
    for fn in fns:
        try:
            fn()
            print(f"PASS {fn.__name__}")
        except AssertionError as e:
            failed += 1
            print(f"FAIL {fn.__name__}: {e}")
    print(f"\n{len(fns) - failed}/{len(fns)} passed")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(_run_standalone())
