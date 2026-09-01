"""Pin the v2 `symbol_rescue` — the first WIDENING refinement in the rule engine.

Both prior refinements (`subclass_any`, `gene_prefixes`) NARROW and compose as AND, which is why the
gentamicin `rmt` gap was not expressible as a config change. A widening refinement is the dangerous
kind, so its containment is pinned here rather than trusted.
"""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from dna_decode.eval.amr_rules import (DRUG_RULE, call_resistance,  # noqa: E402
                                       cipro_determinants_from_main, rule_for)

_HDR = ["Element symbol", "Element name", "Class", "Subclass", "% Identity to reference",
        "Element type", "Method"]


def _tsv(tmp_path, rows):
    p = tmp_path / "main.tsv"
    with p.open("w", encoding="utf-8", newline="") as fh:
        fh.write("\t".join(_HDR) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(h, "")) for h in _HDR) + "\n")
    return p


def _row(sym, cls, sub):
    return {"Element symbol": sym, "Element name": sym, "Class": cls, "Subclass": sub,
            "% Identity to reference": "100", "Element type": "AMR", "Method": "EXACTX"}


def test_symbol_rescue_cannot_escape_the_class_gate(tmp_path):
    """THE safety property. The rescue is evaluated INSIDE the drug-class gate, so it can only re-admit
    a determinant already relevant to the drug. A row matching the pattern but belonging to an unrelated
    class must stay excluded -- otherwise a widening refinement could pull in anything."""
    p = _tsv(tmp_path, [_row("rmtB1", "BETA-LACTAM", "CEPHALOSPORIN")])
    dets = cipro_determinants_from_main(p, "gentamicin",
                                        subclass_any=frozenset({"GENTAMICIN"}),
                                        symbol_rescue=r"^(rmt[A-H]\d*|npmA\d*)$")
    assert dets == [], "the rescue admitted an out-of-class determinant"


def test_the_rescue_admits_the_real_gap(tmp_path):
    """rmt* is filed under the GENERIC AMINOGLYCOSIDE subclass -- in gentamicin's class set, but filtered
    out by the GENTAMICIN subclass refinement. That is the whole gap."""
    p = _tsv(tmp_path, [_row("rmtB1", "AMINOGLYCOSIDE", "AMINOGLYCOSIDE")])
    without = cipro_determinants_from_main(p, "gentamicin", subclass_any=frozenset({"GENTAMICIN"}))
    with_r = cipro_determinants_from_main(p, "gentamicin", subclass_any=frozenset({"GENTAMICIN"}),
                                          symbol_rescue=r"^(rmt[A-H]\d*|npmA\d*)$")
    assert without == [] and len(with_r) == 1


def test_no_rescue_leaves_every_other_drug_byte_identical(tmp_path):
    """Default None ⇒ unchanged behaviour. Five of six drugs must not set it."""
    with_rescue = [d for d, cfg in DRUG_RULE.items() if cfg.get("symbol_rescue")]
    assert with_rescue == ["gentamicin"], f"unexpected drugs carry a rescue: {with_rescue}"


def test_the_deployed_rule_calls_an_rmt_carrier_resistant(tmp_path):
    """End-to-end through the real call_resistance, not the matcher."""
    p = _tsv(tmp_path, [_row("rmtB1", "AMINOGLYCOSIDE", "AMINOGLYCOSIDE")])
    assert call_resistance(p, "gentamicin")["prediction"] == "R"


def test_armA_was_already_counted_and_is_NOT_in_the_rescue():
    """Including armA would overstate the change: AMRFinder files it under Subclass GENTAMICIN, so the
    v1 rule always counted it. The gap was rmt*/npmA only."""
    import re
    pat = re.compile(rule_for("gentamicin")["symbol_rescue"], re.I)
    assert not pat.match("armA"), "armA is in the rescue -- it was already counted by subclass"
    assert pat.match("rmtB1") and pat.match("npmA")


def test_an_unrelated_aminoglycoside_gene_is_still_excluded(tmp_path):
    """aph/aadA are streptomycin/kanamycin determinants that do NOT confer gentamicin resistance. The
    v1 rule excluded them deliberately and v2 must not quietly re-admit them."""
    p = _tsv(tmp_path, [_row("aph(6)-Id", "AMINOGLYCOSIDE", "STREPTOMYCIN"),
                        _row("aadA5", "AMINOGLYCOSIDE", "STREPTOMYCIN")])
    assert call_resistance(p, "gentamicin")["prediction"] == "S"
