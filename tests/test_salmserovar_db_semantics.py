"""The O-axis root cause is a LOST HEADER QUALIFIER, not a missing or defective sequence.

Our antigen DB is SeqSero2's, re-headered into `O__<antigen>__<index>` -- a schema with no field for the
semantic qualifier SeqSero2 attaches (`not_in_3,10`, `wzy_partial`, `wbaV`). The caller then reads every
entry as a standalone positive antigen allele, which produces BOTH measured symptoms from ONE cause:
a differential marker asserted as a positive `1,3,19` call, and the O9-vs-O9,46 discriminator flattened
so plain O9 is unemittable.

These tests pin the evidence that survives without Docker: the two flagged alleles are exactly the
lengths of SeqSero2's qualified entries, and our schema provably cannot express a qualifier.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from dna_decode.salmserovar.runner import parse_axis_antigen, parse_ss2_key  # noqa: E402

DB = ROOT / "data" / "salmserovar_db" / "salmonella_antigens.fasta"
MEMO = ROOT / "wiki" / "salmserovar_db_semantics_2026-09-09.md"

# SeqSero2 1.3.2 entries whose NAME declares special handling, with their lengths.
SS2_QUALIFIED = {
    "O-1,3,19_not_in_3,10": 130,
    "O-9,46_wzy_partial": 216,
    "O-3,10_not_in_1,3,19": 1519,
    "O-9,46_wbaV": 1002,
    "O-9,46,27_partial_wzy": 1019,
}


def _o_lengths() -> dict[str, int]:
    seqs, name, buf = {}, None, []
    for ln in DB.read_text(encoding="utf-8", errors="replace").splitlines():
        if ln.startswith(">"):
            if name:
                seqs[name] = "".join(buf)
            name, buf = ln[1:].strip(), []
        else:
            buf.append(ln.strip())
    if name:
        seqs[name] = "".join(buf)
    return {k: len(v) for k, v in seqs.items()
            if (parse_axis_antigen(k) or ("", ""))[0] == "O"}


def test_the_header_schema_now_carries_the_seqsero2_role():
    """THE FIX, pinned. The schema was `O__<antigen>__<index>` -- three fields, none carrying
    semantics -- which is why the qualifier was lost and the caller read six special-purpose entries
    as ordinary positive alleles. Every O header now carries the ORIGINAL SeqSero2 header as a 4th
    field, and it round-trips exactly (it embeds `__` itself, so it is re-joined, not indexed)."""
    if not DB.exists():
        pytest.skip("antigen DB absent")
    seen = set()
    for k in _o_lengths():
        parts = k.split("__")
        assert len(parts) >= 4, f"{k}: no role field -- the DB is a pre-2026-09-09 build"
        assert parts[0] == "O"
        assert parts[2].isdigit(), f"{k}: third field must remain the source index"
        key = parse_ss2_key(k)
        assert key and key.startswith(("O-", "O:")), f"{k}: role field is not a SeqSero2 O header"
        seen.add(key)
    # Non-vacuity: a schema that merely HAS a 4th field proves nothing. Every one of SeqSero2's six
    # qualified entries -- the ones whose loss caused the defect -- must be recoverable by name.
    for qualified in SS2_QUALIFIED:
        assert any(s.startswith(qualified + "__") for s in seen), f"{qualified} not recoverable"


def test_the_two_flagged_alleles_have_the_lengths_of_seqsero2_qualified_entries():
    """Byte-identity was verified against the live container; length is the part that survives offline.
    If either changes, the DB was rebuilt and the whole diagnosis must be re-derived."""
    if not DB.exists():
        pytest.skip("antigen DB absent")
    by_key = {parse_ss2_key(k): v for k, v in _o_lengths().items()}
    assert by_key.get("O-1,3,19_not_in_3,10__130") == SS2_QUALIFIED["O-1,3,19_not_in_3,10"] == 130
    assert by_key.get("O-9,46_wzy_partial__216") == SS2_QUALIFIED["O-9,46_wzy_partial"] == 216


def test_plain_O9_is_absent_from_BOTH_databases():
    """The correction: SeqSero2 has no plain O9 allele either. O9 is called by wbaV differential logic,
    so 'add a missing O9 sequence' was the wrong fix to reach for."""
    if not DB.exists():
        pytest.skip("antigen DB absent")
    ours = {(parse_axis_antigen(k) or ("", ""))[1] for k in _o_lengths()}
    assert "9" not in ours
    assert not any(k.startswith("O-9_") for k in SS2_QUALIFIED)
    assert "O-9,46_wbaV" in SS2_QUALIFIED      # the discriminator SeqSero2 uses instead


def test_the_memo_corrects_the_external_wall_claim():
    if not MEMO.exists():
        pytest.skip("memo absent")
    t = MEMO.read_text(encoding="utf-8").lower()
    assert "not external" in t
    assert "third time this session i called a wall too early" in t
    assert "deliberately not implemented" in t


def test_the_memo_keeps_the_qualifier_meanings_as_a_reading_not_a_fact():
    """The meanings are inferred from self-documenting names; inferring a rule from a name is exactly
    how the previous wrong causes happened, so the memo must say to read SeqSero2's source first."""
    # Assert on phrases that cannot straddle a markdown line-wrap: the first version checked
    # "did not read seqsero2's source", which the memo happens to wrap mid-phrase.
    t = " ".join(MEMO.read_text(encoding="utf-8").lower().split())
    assert "read seqsero2's source to confirm" in t
    assert "before implementing, read that" in t
