"""A caller's module docstring must not describe the rule it was fixed for.

TWICE NOW the selection rule changed and the prose describing it did not, leaving the module docstring --
the first thing a reader sees -- describing the DEFECT as if it were the design:

  * `serotype/runner.py` said the call was "best-coverage O + best-coverage H" for a day after
    coverage-only selection was replaced by identity-primary (measured H accuracy 0.770 -> 0.926).
  * `resfinder/runner.py` said "a gene is CALLED when its best allele clears thresholds" after the
    locus collapse landed -- that sentence IS the 165-beta-lactam-genes-per-genome behaviour.

Both are fixed. These pins keep the superseded wording from coming back, the same way
`test_hla_help_does_not_advertise_the_demoted_tags` pins a corrected help string. They check the
DOCSTRING only; the behaviour itself is pinned by the selection and locus-collapse tests.

`disinfinder` is deliberately asserted to KEEP the allele-name wording: that caller was measured and
left unchanged (16 alleles, no dense family, 0 multi-reported), so its description is accurate and
"fixing" it would make the docs wrong.
"""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import dna_decode.disinfinder.runner as disinfinder_runner  # noqa: E402
import dna_decode.resfinder.runner as resfinder_runner  # noqa: E402
import dna_decode.salmserovar.runner as salmserovar_runner  # noqa: E402
import dna_decode.serotype.runner as serotype_runner  # noqa: E402


def _doc(mod) -> str:
    return (mod.__doc__ or "").lower()


def test_serotype_does_not_describe_coverage_only_selection():
    d = _doc(serotype_runner)
    assert "best-coverage" not in d, "module docstring describes the superseded coverage-only rule"
    assert "identity-primary" in d
    assert "cross-hybridize" in d, "the docstring should say WHY identity-primary is load-bearing"


def test_resfinder_describes_the_locus_collapse_not_per_allele_calling():
    d = _doc(resfinder_runner)
    assert "locus" in d, "module docstring omits the load-bearing locus collapse"
    assert "a gene is called when its best allele clears thresholds" not in d


def test_salmserovar_already_described_its_rule_correctly():
    """The control: this sibling was updated when it was fixed, which is why it is not in the fix list."""
    d = _doc(salmserovar_runner._best_per_axis)
    assert "identity" in d
    assert "cross-hybridiz" in d or "wrong antigen" in d


def test_disinfinder_KEEPS_the_allele_wording_because_it_was_measured_inert():
    """Not an oversight. Propagating the resfinder fix here without measuring would have been the
    error the plasmid/pneumoserotype probes exist to prevent."""
    d = _doc(disinfinder_runner)
    assert "best" in d and "allele" in d
    assert "locus" not in d, (
        "disinfinder was measured INERT (16 alleles, 0 multi-reported) and left unchanged; if the "
        "caller has since been given a locus collapse, update this test and its cell contract")
