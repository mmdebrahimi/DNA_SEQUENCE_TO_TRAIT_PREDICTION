"""Pins the CYP2D6 hybrid IDENTITY classifier (Phase B) + the full-N falsifier verdict.

Validated at full N (wiki/cyp2d6_psv_phaseb_falsifier.json): non-hybrid specificity 1.0 (the confound test),
*68 4/4, *36 6/8. These tests pin the pure three-level classifier logic (directional 5'-3' -> *68/*13;
exon-9-tip dip -> *36; else unresolved; callability gate -> not_callable) + the committed GO verdict.
"""
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from dna_decode.pgx.cyp2d6_hybrid_identity import MIN_CALLABLE_PSV, classify_hybrid_identity  # noqa: E402


def _profile(five, three, ex9=None):
    """A synthetic region D6-fraction profile."""
    p = {"exon1": five, "upstream_exon1": five, "intron1": five,
         "downstream_exon9": ex9 if ex9 is not None else three, "intron6": three, "exon6": three,
         "intron2": 0.5, "exon2": 0.5}
    return p


def test_star68_directional_5p_high():
    r = classify_hybrid_identity(_profile(0.66, 0.39), n_callable=115)
    assert r["call"] == "*68" and r["resolved"] is True and r["confidence"] == "high"


def test_star13_directional_5p_low_and_unpowered_note():
    r = classify_hybrid_identity(_profile(0.43, 0.71), n_callable=115)
    assert r["call"] == "*13" and r["resolved"] is True
    assert "UNPOWERED" in r["description"]                       # honest: *13 is n=1


def test_star36_exon9_tip_dip():
    # flat 5'-3' but a sharp downstream_exon9 dip -> *36
    r = classify_hybrid_identity(_profile(0.55, 0.60, ex9=0.40), n_callable=115)
    assert r["call"] == "*36" and r["resolved"] is True


def test_flat_profile_is_unresolved_not_forced():
    # a hybrid-present sample whose profile is flat/ambiguous -> NEVER forced to a specific allele
    r = classify_hybrid_identity(_profile(0.55, 0.56, ex9=0.55), n_callable=115)
    assert r["call"] == "hybrid_present_identity_unresolved" and "resolved" not in r


def test_callability_gate():
    r = classify_hybrid_identity(_profile(0.66, 0.39), n_callable=MIN_CALLABLE_PSV - 1)
    assert r["call"] == "evidence_not_callable"


def test_specificity_high_no_false_call_on_near_flat():
    # a slightly noisy but non-directional profile (a pure dup's flat-elevated shape) -> not resolved
    r = classify_hybrid_identity(_profile(0.60, 0.58, ex9=0.60), n_callable=200)
    assert not r.get("resolved")


def test_committed_full_n_falsifier_is_go():
    import json
    p = REPO / "wiki" / "cyp2d6_psv_phaseb_falsifier.json"
    if not p.exists():
        import pytest
        pytest.skip("falsifier artifact not present")
    d = json.loads(p.read_text(encoding="utf-8"))
    assert d["verdict"] == "GO_BUILD_CLASSIFIER"
    assert d["nonhybrid_specificity"] == 1.0                     # the confound test
    assert d["per_allele"]["*68"]["rate"] == 1.0                 # *68 clean at full N
    assert (d["per_allele"]["*36"]["rate"] or 0) >= 0.6          # *36 partial (subtle conversions)


if __name__ == "__main__":
    import pytest
    raise SystemExit(pytest.main([__file__, "-q"]))
