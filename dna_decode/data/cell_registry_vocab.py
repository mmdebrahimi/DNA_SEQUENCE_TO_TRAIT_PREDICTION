"""Abstention vocabulary — a shared cross-cell GROUPING, NOT a confidence scale (Evidence-Contract Registry v0).

The suite expresses "this cell is not a clean SCORED number" in ~10 heterogeneous in-tree strings
(`ABSTAIN`, `SUSPEND`, `phenotype_withheld`, `O?`, …) scattered across AMR / PGx / genome-map / typing
code. This module collapses those NATIVE terms into ONE controlled enum so a reader can GROUP cells by
why-they-abstain across the whole surface.

HONESTY RAIL (load-bearing, from the pre-exec /brainstorm): this is a VOCABULARY, never a probability or a
ranked scale. The enum has NO numeric value, NO ordering semantics, and members are deliberately string-valued.
Two cells sharing a vocab value are "abstaining for the same KIND of reason" — they are NOT "equally confident".
A cell keeps its own `native_abstention` string; the vocab is only for cross-cell grouping.
"""
from __future__ import annotations

from enum import Enum


class AbstentionVocab(str, Enum):
    """Controlled cross-cell abstention vocabulary. String-valued; NO numeric scale (guardrail)."""

    SCORED = "scored"                      # a real validated number exists for this cell
    ABSTAIN_BY_DESIGN = "abstain_by_design"  # the engine deliberately declines (e.g. intrinsic-gene organism)
    UNDERPOWERED = "underpowered"          # too few labels of one class to score honestly
    NOT_CENSUSED = "not_censused"          # shipped but never run against a label set
    WITHHELD_NONCORE = "withheld_noncore"  # a non-core determinant was seen -> phenotype withheld (PGx sentinel)
    NO_FREE_SOURCE = "no_free_source"      # no free isolate-level phenotype label source exists at all
    LABEL_CONFOUNDED = "label_confounded"  # the phenotype LABEL is an unreliable surrogate (oxacillin/mecA)
    SUSPEND = "suspend"                    # the audit/merge gate suspended the cell (dirty labels)
    GATE_BLOCKED = "gate_blocked"          # a prevent-wrong-inference gate said NO-GO (genome-map)


# The 10 heterogeneous in-tree native terms -> the controlled vocab. Adding a new native term is a one-line
# edit here; tests/test_cell_registry.py asserts every value is a real AbstentionVocab member and that the
# vocab carries no numeric scale.
NATIVE_TO_VOCAB: dict[str, AbstentionVocab] = {
    "SCORED": AbstentionVocab.SCORED,
    "ABSTAIN": AbstentionVocab.ABSTAIN_BY_DESIGN,
    "ABSTAINS_BY_DESIGN": AbstentionVocab.ABSTAIN_BY_DESIGN,
    "UNDERPOWERED": AbstentionVocab.UNDERPOWERED,
    "NOT_CENSUSED": AbstentionVocab.NOT_CENSUSED,
    "NO_GO": AbstentionVocab.GATE_BLOCKED,
    "SUSPEND": AbstentionVocab.SUSPEND,
    "NO_FREE_PHENOTYPE": AbstentionVocab.NO_FREE_SOURCE,
    "LABEL_CONFOUNDED": AbstentionVocab.LABEL_CONFOUNDED,
    "phenotype_withheld": AbstentionVocab.WITHHELD_NONCORE,
    "O?": AbstentionVocab.ABSTAIN_BY_DESIGN,
}


def to_vocab(native: str) -> AbstentionVocab:
    """Collapse a native in-tree abstention term to the controlled vocab. Raises KeyError on an unknown term."""
    return NATIVE_TO_VOCAB[native]
