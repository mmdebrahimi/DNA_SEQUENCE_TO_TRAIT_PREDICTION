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
    has_16s_methyltransferase,
    parse_main_tsv_rows,
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


def test_16s_methyltransferase_detection():
    assert has_16s_methyltransferase([{"symbol": "rmtE1"}])
    assert has_16s_methyltransferase([{"symbol": "armA"}])
    assert has_16s_methyltransferase([{"symbol": "npmA"}])
    assert not has_16s_methyltransferase([{"symbol": "aph(3')-Ia"}, {"symbol": "aadA5"}])


def test_the_specific_warning_fires_only_when_a_methyltransferase_is_present():
    """REGRESSION: the first version printed the gentamicin/rmt paragraph unconditionally, so a
    ciprofloxacin call that correctly flagged `qnrB19` also carried an irrelevant gentamicin lecture --
    noise that trains a reader to skip the note entirely."""
    with_rmt = render_note([{"symbol": "rmtE1", "class": "AMINOGLYCOSIDE", "subclass": "AMINOGLYCOSIDE"}],
                           "gentamicin")
    assert "16S rRNA METHYLTRANSFERASE" in with_rmt and "UNRELIABLE" in with_rmt

    without = render_note([{"symbol": "qnrB19", "class": "QUINOLONE", "subclass": "QUINOLONE"}],
                          "ciprofloxacin")
    assert "qnrB19" in without
    assert "16S" not in without and "gentamicin" not in without    # no cross-drug lecture
    assert "DELIBERATE" in without                                 # the honest generic framing survives

    assert render_note([], "gentamicin") == ""                     # nothing to disclose -> silent


def test_disclosure_never_changes_a_prediction_and_never_touches_the_frozen_surface():
    """The two invariants that make this safe to ship on a LOCKED decoder."""
    import hashlib
    import json
    root = Path(__file__).resolve().parent.parent
    manifest = json.loads((root / "wiki" / "prospective_lock_manifest_2026-06-22.json")
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
    # the ONE frozen thing it may read is the per-drug class catalog (read-only, no decision)
    assert imported <= {"amrfinder_classes_for"}, f"unexpected frozen-surface import: {imported}"


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
