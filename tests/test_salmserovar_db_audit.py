"""The serovar O-axis gap is reference-DB content, and the audit that proves it must not over-flag.

Two things these tests protect. First, the DIAGNOSIS: the shipped DB has no plain O9 allele and carries a
130 bp `1,3,19` fragment that fires on 60% of a cohort where its antigen is 2% -- neither fixable by a
threshold. Second, the CHECK ITSELF: a first version used a flat hit-rate bar and flagged `O__4__231`,
whose 32% hit-rate exactly matches O=4 being 32% of the cohort. That was the check being wrong, not the
allele, so promiscuity is now measured against each antigen's OWN prevalence.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from salmserovar_db_audit import (  # noqa: E402
    PROMISCUITY_FLOOR, PROMISCUITY_OVER_PREVALENCE, SHORT_ALLELE_BP, antigen_of, axis_of, load_fasta,
)

DB = ROOT / "data" / "salmserovar_db" / "salmonella_antigens.fasta"
ARTIFACT = ROOT / "wiki" / "salmserovar_db_audit_2026-09-09.json"


@pytest.fixture(scope="module")
def art() -> dict:
    if not ARTIFACT.exists():
        pytest.skip("db audit artifact absent")
    return json.loads(ARTIFACT.read_text(encoding="utf-8"))


# --- the DB defects, checked against the live DB not just the artifact ----------------------------

def test_the_db_still_has_no_plain_O9_allele():
    """The structural defect. If this starts failing, a plain O9 allele was added -- re-run the
    SeqSero2 comparison, because the 9,46 over-call should collapse with it."""
    if not DB.exists():
        pytest.skip("antigen DB absent")
    o_antigens = {antigen_of(k) for k in load_fasta(DB) if axis_of(k) == "O"}
    assert "9" not in o_antigens, (
        "a plain O9 allele now exists -- the missing-antigen defect may be fixed; re-measure the "
        "O axis against SeqSero2 and update the diagnosis")
    assert "9,46" in o_antigens          # non-vacuous: the parser does find O antigens


def test_the_promiscuous_allele_is_still_the_shortest_and_still_present():
    if not DB.exists():
        pytest.skip("antigen DB absent")
    lens = {k: len(v) for k, v in load_fasta(DB).items() if axis_of(k) == "O"}
    # The header gained a 4th field (the SeqSero2 role) on 2026-09-09, so match on the stable
    # `<axis>__<antigen>__<index>` prefix rather than the whole id. The AUDIT ARTIFACT this file's
    # other tests read still carries the pre-role ids -- it is a record of a run, not of the schema.
    hits = [k for k in lens if k == "O__1,3,19__126" or k.startswith("O__1,3,19__126__")]
    assert len(hits) == 1, f"expected exactly one 1,3,19 marker entry, got {hits}"
    assert lens[hits[0]] == min(lens.values()) == 130


# --- the audit must not over-flag -----------------------------------------------------------------

def test_promiscuity_is_measured_against_prevalence_not_a_flat_rate(art):
    """The correction: O__4__231 hits 32% and O=4 IS 32% of the cohort, so it must NOT be flagged."""
    flagged = {p["allele"] for p in art["promiscuous_alleles"]}
    assert "O__1,3,19__126" in flagged
    assert "O__4__231" not in flagged, (
        "a correctly-prevalent allele is being flagged as promiscuous -- the check regressed to a "
        "flat hit-rate bar")


def test_every_promiscuity_flag_carries_its_prevalence_baseline(art):
    """A hit-rate without the baseline it is judged against is not evidence."""
    for p in art["promiscuous_alleles"]:
        assert "antigen_true_prevalence" in p and "times_over_prevalence" in p
        assert p["times_over_prevalence"] > PROMISCUITY_OVER_PREVALENCE
        assert p["fraction"] >= PROMISCUITY_FLOOR


def test_the_missing_antigen_check_is_data_driven(art):
    """It lists only antigens the REFERENCE actually emitted here -- not a curated 'should exist' list,
    which would be asserting biology from memory."""
    for m in art["missing_antigens"]:
        assert m["reference_called_it_on"] > 0
    assert any(m["antigen"] == "9" for m in art["missing_antigens"])


def test_the_audit_is_non_vacuous(art):
    assert art["n_o_alleles"] > 20
    assert art["verdict"] == "DB_CONTENT_DEFECTS_FOUND"


# --- honesty rails --------------------------------------------------------------------------------

def test_the_audit_declares_it_does_not_fix_anything(art):
    """Fixing needs sourced sequences; inventing them is the fabrication hazard."""
    assert "does not fix" in art["what_this_does_not_do"].lower()
    assert "fabrication hazard" in art["what_this_does_not_do"].lower()


def test_deleting_the_bad_allele_is_not_proposed(art):
    """It is the DB's ONLY 1,3,19 entry, so removal trades false positives for losing E4 entirely."""
    joined = " ".join(art["honest_limits"]).lower()
    assert "only" in joined and "1,3,19" in joined


def test_absence_of_a_flag_is_not_claimed_as_completeness(art):
    joined = " ".join(art["honest_limits"]).lower()
    assert "not proof of completeness" in joined


def test_the_tripwires_are_declared_as_tripwires_not_derived_constants(art):
    joined = " ".join(art["honest_limits"]).lower()
    assert "tripwire" in joined
    assert SHORT_ALLELE_BP == 300
