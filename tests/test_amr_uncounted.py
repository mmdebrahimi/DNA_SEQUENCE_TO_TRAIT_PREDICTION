"""The uncounted-class-determinant DISCLOSURE (dna_decode/amr/uncounted.py).

Surfaces determinants of a drug-relevant AMRFinder CLASS that the frozen rule did NOT count, at the point
of use. It exists because of a MEASURED blind spot: prospective E. coli x gentamicin scored sens 0.429 and
24 of the 28 false negatives carried a 16S rRNA methyltransferase (rmtE1/armA), which AMRFinder files under
the generic `AMINOGLYCOSIDE` subclass where a `Subclass=GENTAMICIN` rule cannot see it
(`wiki/prospective_lock_first_accrual_2026-08-24.md`).

INVARIANT: this is a disclosure, never a rule. It must never change a prediction, and it lives outside the
sha256-pinned frozen surface.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from dna_decode.amr.uncounted import (  # noqa: E402
    disclosure_provenance,
    measured_gap_misses,
    parse_main_tsv_rows,
    primary_mechanism_misses,
    render_note,
    uncounted_class_determinants,
)

_HDR = "Protein id\tContig id\tStart\tStop\tStrand\tElement symbol\tName\tScope\tElement type\tElement subtype\tClass\tSubclass"


def _tsv(*rows: tuple[str, str, str]) -> str:
    out = [_HDR]
    for sym, cls, sub in rows:
        out.append("\t".join(["NA", "c1", "1", "9", "+", sym, "n", "core", "AMR", "AMR", cls, sub]))
    return "\n".join(out) + "\n"


def test_parse_main_tsv_is_header_driven_not_positional():
    """AMRFinder's column set has shifted across versions; positional parsing silently mis-indexes."""
    rows = parse_main_tsv_rows(_tsv(("rmtE1", "AMINOGLYCOSIDE", "AMINOGLYCOSIDE"),
                                    ("aac(3)-VIa", "AMINOGLYCOSIDE", "GENTAMICIN")))
    assert [r["symbol"] for r in rows] == ["rmtE1", "aac(3)-VIa"]
    assert rows[0]["subclass"] == "AMINOGLYCOSIDE" and rows[1]["subclass"] == "GENTAMICIN"
    assert parse_main_tsv_rows("") == []
    assert parse_main_tsv_rows("no\theader\twe\tknow\n") == []      # unknown header -> no junk rows


def test_uncounted_excludes_what_the_rule_already_counted():
    rows = parse_main_tsv_rows(_tsv(("aac(3)-VIa", "AMINOGLYCOSIDE", "GENTAMICIN"),
                                    ("rmtE1", "AMINOGLYCOSIDE", "AMINOGLYCOSIDE"),
                                    ("aph(3')-Ia", "AMINOGLYCOSIDE", "KANAMYCIN")))
    counted = [{"symbol": "aac(3)-VIa"}]
    got = {u["symbol"] for u in uncounted_class_determinants(rows, "gentamicin", counted)}
    assert got == {"rmtE1", "aph(3')-Ia"}          # the counted one is not re-reported


def test_a_determinant_of_an_unrelated_class_is_not_disclosed():
    """Only classes the drug's own catalog calls relevant -- otherwise the note becomes noise."""
    rows = parse_main_tsv_rows(_tsv(("blaTEM-1", "BETA-LACTAM", "BETA-LACTAM"),
                                    ("rmtE1", "AMINOGLYCOSIDE", "AMINOGLYCOSIDE")))
    got = {u["symbol"] for u in uncounted_class_determinants(rows, "gentamicin", [])}
    assert got == {"rmtE1"}                        # the beta-lactam is irrelevant to gentamicin


def test_primary_mechanism_misses_excludes_deliberate_exclusions():
    """The subclass test is load-bearing and cannot be replaced by the mechanism test alone.

    `classify_gene_symbol('gentamicin', "aph(3')-Ia")` returns `aminoglycoside_modifying_enzymes`, which IS
    a primary gentamicin mechanism -- so a mechanism-only rule would flag a KANAMYCIN gene as a miss. What
    separates them is that aph(3')-Ia carries a specific, different drug subclass while rmtE1 carries the
    bare CLASS name and is therefore invisible to any subclass-matching rule.
    """
    rows = parse_main_tsv_rows(_tsv(("rmtE1", "AMINOGLYCOSIDE", "AMINOGLYCOSIDE"),
                                    ("aph(3')-Ia", "AMINOGLYCOSIDE", "KANAMYCIN"),
                                    ("aadA5", "AMINOGLYCOSIDE", "STREPTOMYCIN")))
    u = uncounted_class_determinants(rows, "gentamicin", [])
    assert len(u) == 3                                   # all three are class-relevant + uncounted
    got = {m["symbol"]: m["mechanism"] for m in primary_mechanism_misses(u, "gentamicin")}
    assert got == {"rmtE1": "16S_rRNA_methyltransferase"}


def test_measured_gap_gating_is_evidence_backed_not_inferred():
    """Only (drug, mechanism) pairs with a MEASURED gap reach the human note.

    Both broader triggers were rejected BY MEASUREMENT: class-relevant fired on 70% of gentamicin calls,
    and primary-mechanism still fired on 48% of CEFTRIAXONE calls (a narrow-spectrum blaTEM-1 also files
    under the generic BETA-LACTAM subclass, and excluding it from ceftriaxone is correct).
    """
    rows = parse_main_tsv_rows(_tsv(("rmtE1", "AMINOGLYCOSIDE", "AMINOGLYCOSIDE")))
    u = rows and uncounted_class_determinants(rows, "gentamicin", [])
    assert [m["symbol"] for m in measured_gap_misses(u, "gentamicin")] == ["rmtE1"]
    assert measured_gap_misses(u, "gentamicin")[0]["evidence"]

    # a beta-lactamase under the generic BETA-LACTAM subclass IS a primary-mechanism miss, but ceftriaxone
    # has no MEASURED gap on record -> it must not reach the human note
    cef = parse_main_tsv_rows(_tsv(("blaTEM-1", "BETA-LACTAM", "BETA-LACTAM")))
    ucef = uncounted_class_determinants(cef, "ceftriaxone", [])
    assert measured_gap_misses(ucef, "ceftriaxone") == []


def test_disclosure_provenance_names_and_hashes_the_catalog():
    """Without recording WHICH catalog spoke, two decoder generations could silently diverge."""
    p = disclosure_provenance()
    assert p["disclosure_catalog"] == "dna_decode/data/mic_tiers.py"
    assert len(p["disclosure_catalog_sha256"]) == 64
    assert p["alters_frozen_decision"] is False


def test_the_note_never_cites_another_drugs_measurement():
    """REGRESSION, twice over. v1 printed the gentamicin/rmt paragraph on EVERY drug. v2 narrowed the
    trigger but still hardcoded the gentamicin sentence, so a ciprofloxacin call flagging `qnrB19` cited
    a gentamicin measurement. Evidence is now per-(drug, mechanism) and cannot travel."""
    gent = parse_main_tsv_rows(_tsv(("rmtE1", "AMINOGLYCOSIDE", "AMINOGLYCOSIDE")))
    note = render_note(uncounted_class_determinants(gent, "gentamicin", []), "gentamicin", "S")
    assert "rmtE1" in note and "0.429" in note and "UNRELIABLE" in note

    cip = parse_main_tsv_rows(_tsv(("qnrB19", "QUINOLONE", "QUINOLONE")))
    cnote = render_note(uncounted_class_determinants(cip, "ciprofloxacin", []), "ciprofloxacin", "S")
    assert cnote == "", "ciprofloxacin has no measured gap on record; the note must stay silent"

    # an R call still discloses, but must not tell the reader to distrust an S it did not make
    rnote = render_note(uncounted_class_determinants(gent, "gentamicin", []), "gentamicin", "R")
    assert "does NOT change the call" in rnote and "UNRELIABLE" not in rnote

    assert render_note([], "gentamicin", "S") == ""


def test_disclosure_never_changes_a_prediction_and_never_touches_the_frozen_surface():
    """The two invariants that make this safe to ship on a LOCKED decoder."""
    import hashlib
    import json
    root = Path(__file__).resolve().parent.parent
    manifest = json.loads((root / "wiki" / "prospective_lock_manifest_2026-08-31.json")
                          .read_text(encoding="utf-8"))
    for rel, want in manifest["surface_sha256"].items():
        got = hashlib.sha256((root / rel).read_bytes()).hexdigest()
        assert got == want, f"{rel} drifted from the prospective lock"

    # and the disclosure module must not be able to reach the decision path at all. Checked by AST, not
    # by substring: the module's docstring legitimately DISCUSSES call_resistance in prose, so a naive
    # `"call_resistance" not in src` fails on its own documentation.
    import ast
    tree = ast.parse((root / "dna_decode" / "amr" / "uncounted.py").read_text(encoding="utf-8"))
    called = {n.func.id for n in ast.walk(tree)
              if isinstance(n, ast.Call) and isinstance(n.func, ast.Name)}
    imported = {a.name for n in ast.walk(tree) if isinstance(n, ast.ImportFrom)
                for a in n.names if n.module != "__future__"}
    assert "call_resistance" not in called | imported, \
        "the disclosure layer must never invoke or import the deployed rule"

    # It MAY read frozen CATALOGS (pure lookups, no decision); it may NOT touch the frozen DECISION
    # module. Asserted by module, not by a literal name list -- an earlier version whitelisted exactly
    # {"amrfinder_classes_for"} and had to be widened by hand the moment a second catalog lookup was
    # added, which is the drift this guard is supposed to catch rather than commit.
    frozen_decision_modules = {"dna_decode.eval.amr_rules"}
    modules = {n.module for n in ast.walk(tree) if isinstance(n, ast.ImportFrom) and n.module}
    assert not (modules & frozen_decision_modules), \
        f"disclosure imports the frozen DECISION path: {modules & frozen_decision_modules}"
    catalog_only = {m for m in modules if m.startswith("dna_decode.")} <= {"dna_decode.data.mic_tiers"}
    assert catalog_only, f"unexpected dna_decode import in the disclosure layer: {modules}"


@pytest.mark.parametrize("drug,sym,cls,sub,expect", [
    ("gentamicin", "rmtE1", "AMINOGLYCOSIDE", "AMINOGLYCOSIDE", True),
    ("gentamicin", "aadA5", "AMINOGLYCOSIDE", "STREPTOMYCIN", True),
    ("ciprofloxacin", "qnrB19", "QUINOLONE", "QUINOLONE", True),
    ("ciprofloxacin", "rmtE1", "AMINOGLYCOSIDE", "AMINOGLYCOSIDE", False),
])
def test_relevance_is_per_drug(drug, sym, cls, sub, expect):
    rows = parse_main_tsv_rows(_tsv((sym, cls, sub)))
    assert bool(uncounted_class_determinants(rows, drug, [])) is expect


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-q"]))
