"""Apply the ruled coordinate-overlap purge to a human lentiMPRA cohort. OPT-IN.

WHAT WAS RULED AND BY WHOM. The publishers assign folds randomly over ELEMENTS rather than by
chromosome, so genomically overlapping elements can take different fold roles while sharing most of
their sequence (measured 2026-11-23: 625 such pairs across the three cells, 215 sharing >=50 bp).
Two remedies were costed and the user chose PURGE -- drop the touched HELD-OUT elements, keep the
publishers' folds -- because a chromosome-held-out rebuild would permanently forfeit the only
external comparison this substrate has (the published per-fold Enformer / Sei / MPRAnn / MPRALegNet
scores are defined on THEIR folds). The ruling reached this session through a peer relay, which is
recorded in the list artifact rather than assumed.

WHY THIS IS OPT-IN AND load_cell IS UNTOUCHED. `lentimpra.load_cell` still returns every joined
record, so the sha256 regression that proves HepG2 and K562 are byte-identical across the WTC11
rescue still holds. A consumer that wants the purged cohort asks for it by name. Making the purge
the default would silently change what every existing caller measures, which is the failure this
module's own subject matter is about.

THE RULED SET IS WIDER THAN THE RESCORED SET, and they must not be quoted interchangeably. The
rescore that showed the exposure immaterial (worst delta 0.0014, about 0.05x the matched-n sd)
dropped only elements a TRAINING fold touched: 27 HepG2 and 131 K562. The ruled set also drops
elements whose only conflict is fold 9 against fold 10 -- 4, 20 and 2 of them. So the 0.0014 figure
describes the narrower set and says nothing about these 26.

FOLD 10 IS NOT SPENT BY ANY OF THIS. The list was generated from names, folds and hg38 coordinates
only; no activity was ever scored against fold 10, so its FINAL_ONLY reservation is intact.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

#: The committed, auditable list. A names-only artifact rather than a recomputation, because the
#: generating sweep needs a 56 MB workbook parse and a cohort definition must not depend on a
#: re-derivation that could silently drift from what was approved.
DEFAULT_LIST = (Path(__file__).resolve().parents[2] / "wiki"
                / "human_fold_leak_purge_list_2026-11-23.json")


class FoldLeakPurgeError(ValueError):
    """Refuse rather than apply a purge that does not match what was ruled."""


def _payload(path: Path | str | None = None) -> dict:
    p = Path(path) if path is not None else DEFAULT_LIST
    if not p.is_file():
        raise FoldLeakPurgeError(
            "%s is absent. The purge is defined by a COMMITTED list, not by recomputation, so a "
            "missing list is a refusal rather than a silent no-op -- a no-op would hand back the "
            "unpurged cohort while the caller believed it was purged." % p)
    return json.loads(p.read_text(encoding="utf-8"))


def purge_names(cell: str, path: Path | str | None = None) -> frozenset[str]:
    """The element names ruled out for `cell`, verified against the list's own sha256.

    The digest is checked because this list is the one thing standing between a caller and a cohort
    nobody approved; an edited list that still parses would otherwise apply silently.
    """
    d = _payload(path)
    cells = d.get("cells") or {}
    if cell not in cells:
        raise FoldLeakPurgeError(
            "no purge entry for cell %r; the list covers %s. A cell absent from the list is NOT the "
            "same as a cell with nothing to purge, and returning an empty set would conflate them."
            % (cell, sorted(cells)))
    names = cells[cell]["names"]
    want = (d.get("sha256_of_names_only") or {}).get(cell)
    if not want:
        raise FoldLeakPurgeError(
            "the list carries no sha256 for %r, so its contents cannot be verified" % cell)
    got = hashlib.sha256("\n".join(names).encode()).hexdigest()
    if got != want:
        raise FoldLeakPurgeError(
            "purge list for %r does not match its own digest (%s vs %s). It has been edited since "
            "it was generated; refusing rather than applying an unapproved set."
            % (cell, got[:16], want[:16]))
    declared, actual = cells[cell]["n_purged"], len(names)
    if declared != actual:
        raise FoldLeakPurgeError(
            "purge list for %r declares %d names and carries %d" % (cell, declared, actual))
    ruled = cells[cell].get("ruled_count")
    if ruled is not None and ruled != actual:
        raise FoldLeakPurgeError(
            "purge list for %r carries %d names but the RULING was %d" % (cell, actual, ruled))
    return frozenset(names)


def apply_purge(records, cell: str, path: Path | str | None = None):
    """`(kept, n_dropped)` for an iterable of `LentiRecord`.

    Returns the count rather than only the filtered list so a caller can assert it matched the
    ruling. A purge that silently removed a different number than approved is the whole hazard.
    """
    drop = purge_names(cell, path)
    recs = list(records)
    kept = [r for r in recs if r.name not in drop]
    n_dropped = len(recs) - len(kept)
    if n_dropped != len(drop):
        raise FoldLeakPurgeError(
            "%s: the list names %d elements but only %d were present in the cohort. The list and "
            "the cohort disagree, so neither the purged nor the unpurged set is trustworthy here."
            % (cell, len(drop), n_dropped))
    return kept, n_dropped


def purge_summary(path: Path | str | None = None) -> dict:
    """Per-cell counts and composition, for an artifact or a report. Reads, never recomputes."""
    d = _payload(path)
    return {
        "date": d.get("date"),
        "authority": d.get("authority"),
        "cells": {c: {"n_purged": v["n_purged"], "per_fold": v["per_fold"],
                      "composition": v["composition"]}
                  for c, v in (d.get("cells") or {}).items()},
        "NOT_the_rescored_set": d.get("NOT_the_rescored_set"),
    }
