"""Regression tests pinning the locked sections of the two-machine contract.

Dogfoods the contract's own `Contract Locks` rule: every locked parameter in
plans/Two_Machine_Operating_Contract.md MUST have a same-commit regression
test or KNOWN_DIVERGENCE_TARGETS marker. This file is the test side.

If a future edit silently weakens any of the 5 locked sections, these tests
fail before merge.
"""
from __future__ import annotations

from pathlib import Path

import pytest


CONTRACT_PATH = Path("plans/Two_Machine_Operating_Contract.md")


@pytest.fixture(scope="module")
def contract_text() -> str:
    """Contract markdown with whitespace collapsed.

    Collapses all runs of whitespace (incl. newlines + indentation from
    markdown line-wrapping) to single spaces so test assertions can use
    natural English phrases without worrying about where the markdown wrapped.
    """
    if not CONTRACT_PATH.exists():
        pytest.fail(f"contract document missing at {CONTRACT_PATH}")
    raw = CONTRACT_PATH.read_text(encoding="utf-8")
    return " ".join(raw.split())


def test_two_lanes_locked(contract_text: str):
    """§1: the two-lane split is the load-bearing structural lock."""
    # Both lane names must appear, with the right machine assignment
    assert "Execution" in contract_text and "Codex" in contract_text and "Precision 7780" in contract_text, \
        "execution lane (Codex / Precision 7780) must remain locked in §1"
    assert "Discovery + planning + contract-lock conversion" in contract_text, \
        "discovery lane name must include 'contract-lock conversion' (the Issue 2 amendment)"
    assert "GTX 860M" in contract_text, \
        "Claude side must be identified by the GTX 860M machine spec"
    # Hard rule
    assert "Discovery lane does NOT start implementation" in contract_text, \
        "hard rule against discovery-lane implementation must remain"


def test_handoff_gate_5_checks(contract_text: str):
    """§2.1: all 5 handoff-gate checks must be present + numbered."""
    required_phrases = [
        "Sender push status",
        "Receiver pull",
        "Sync check",
        "Locked-parameter coverage",
        "Contract tests run",
    ]
    for phrase in required_phrases:
        assert phrase in contract_text, f"handoff gate check '{phrase}' is missing from §2.1"
    # Stop-condition on drift
    assert "hard pause" in contract_text.lower(), \
        "stop condition on drift must remain 'hard pause' per §2.2"


def test_contract_locks_section_spec(contract_text: str):
    """§4: the Contract Locks section spec must require enforcement targets."""
    assert "## Contract Locks" in contract_text, \
        "§4 must reference the literal '## Contract Locks' section heading"
    assert "Regression test" in contract_text and "Sync marker" in contract_text and "Not lock-bearing" in contract_text, \
        "§4 must enumerate the 3 valid enforcement targets"
    assert "Empty enforcement target = unacceptable" in contract_text, \
        "§4 must prohibit empty enforcement targets"
    assert "Single-commit rule" in contract_text, \
        "§4 must lock the single-commit rule for parameter + enforcement target"


def test_five_falsification_triggers(contract_text: str):
    """§6: all 5 falsification triggers must be present."""
    # Each trigger has a distinct keyword
    required_triggers = [
        "producing commit is not on `origin/main`",
        "without a same-commit regression test or divergence marker",
        "reports drift after handoff but work continues",
        "more than 24 hr",
        "regression test has to be added LATER",
    ]
    for trigger in required_triggers:
        assert trigger in contract_text, f"falsification trigger '{trigger}' is missing from §6"


def test_single_commit_rule_documented(contract_text: str):
    """§4 + §2.1.4: the single-commit rule must be explicit."""
    assert "SAME commit" in contract_text, \
        "§2.1 item 4 must require the regression test/marker in the SAME commit as the lock"
    # Historical exemption — but only as historical, not future norm
    assert "grandfathered" in contract_text.lower(), \
        "historical split (cd76d6c + 716d214) must be explicitly grandfathered, not normalized"


def test_source_of_truth_ownership_model(contract_text: str):
    """§3: ownership rules must be source-of-truth-based, not author-based."""
    assert "source-of-truth" in contract_text.lower(), \
        "§3 must use the source-of-truth ownership model (Issue 3 amendment)"
    assert "NOT author ownership" in contract_text, \
        "§3 must explicitly reject author-based ownership"
    assert "provisional" in contract_text.lower(), \
        "§3 must define 'provisional' status for cross-lane proposals"


def test_handoff_manifest_required_for_non_repo_artifacts(contract_text: str):
    """§5: handoff manifest for large/gitignored artifacts."""
    assert "Handoff Manifest" in contract_text, \
        "§5 must define the Handoff Manifest format"
    assert "24-hr alarm" in contract_text or "24 hr" in contract_text, \
        "§5 must define the 24-hr alarm for uncommitted contract-bearing artifacts"


def test_first_production_test_cef_closeout(contract_text: str):
    """§7: the first production test of this contract is cef closeout."""
    assert "cef audit-aware closeout" in contract_text, \
        "§7 must name cef audit-aware closeout as the first production test"
    assert "10/10 markers" in contract_text or "10 markers" in contract_text, \
        "§7 must specify the expected 10-marker post-closeout state"
