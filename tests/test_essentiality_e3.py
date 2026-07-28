"""Guard: the E3 learned complement beats the conserved-core baseline (skips if D: data absent)."""
import json, os
import pytest
_HAVE = os.path.exists("D:/dna_decode_cache/essentiality/goodall_TableS1_essential.xlsx")


@pytest.mark.skipif(not _HAVE, reason="essentiality D: data absent")
def test_e3_learned_beats_conserved_core():
    r = json.load(open("wiki/essentiality_e3_learned_2026-07-28.json"))
    assert r["learned_auroc"] > r["conserved_core_auroc"]           # learned earns keep
    assert r["auroc_lift"] > 0.02                                    # meaningful lift
    assert r["tail_recovery"]["learned_auroc_on_tail"] > 0.55        # recovers the core-missed tail
    assert "EARNS_KEEP" in r["verdict"]


def test_human_e3_learned_beats_transfer():
    import json, os
    if not os.path.exists("wiki/essentiality_e3_human_2026-07-28.json"):
        import pytest; pytest.skip("human E3 artifact absent")
    r = json.load(open("wiki/essentiality_e3_human_2026-07-28.json"))
    assert r["learned_gbm_auroc"] > r["conserved_core_auroc"]
    assert r["auroc_lift"] > 0.1                        # human lift is large (E.coli-tuned core transfers poorly)
    assert "EARNS_KEEP" in r["verdict"]
