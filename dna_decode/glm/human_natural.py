"""The HUMAN natural reference corpus for the distinguishability falsifier.

WHY THIS MODULE EXISTS AS A THIRD PATH. `scripts/glm_generator_falsifier.py:build_natural` accepts
`kind` in {promoter, genomic}, but BOTH branches go through `load_fasta(genome)` -- so the natural
set is organism-locked by the genome file rather than by `kind`. A human natural set needs sequences
taken directly, with no genome and no GFF3. This module supplies them and touches neither
`dna_decode/glm/corpus.py` nor the falsifier CLI, so the bacterial path stays byte-identical and
nothing here can alter a committed bacterial number.

THE SPLIT IS THE PUBLISHER'S. `load_cell` carries the authors' own 10-fold assignment through to
each record, so `fit_eval_by_fold` partitions on that rather than on a split invented here. A Markov
null fitted on an invented split and scored on its complement would be measuring our own
partitioning as much as the sequence.

GEOMETRY. The publisher's sequence is `15nt 5' adaptor - 200nt element - 15nt 3' adaptor` = 230 bp.
`wiki/human_falsifier_bar_2026-11-21.json` fixes the FULL 230 bp as primary, because that is the
geometry the measured activity labels and the 2026-11-01 encoder arm used, with the 200 nt element as
a secondary row -- constant flanks are trivially learnable and could move distinguishability for a
reason that is not about regulatory structure.

REFUSALS, not degradations. An absent directory, an unusable cell or an off-length sequence raises,
carrying the loader's own reason string. A natural set that silently drops malformed rows would make
a distinguishability number that looks fine and means something else.
"""

from __future__ import annotations

from pathlib import Path

from dna_decode.glm.falsifier import alphabet_and_length_ok
from dna_decode.glm.lentimpra import CELL_TYPES, LentiMpraDataError, load_cell

#: The publisher's full element length, adaptors included.
FULL_BP = 230
#: Adaptor length at EACH end. 230 - 2*15 = 200, the element proper.
ADAPTOR_BP = 15
ELEMENT_BP = FULL_BP - 2 * ADAPTOR_BP


class HumanNaturalError(RuntimeError):
    """Raised when a human natural set cannot be built as specified."""


def _require_dir(data_dir: Path | str) -> Path:
    p = Path(data_dir)
    if not p.is_dir():
        raise HumanNaturalError(
            "data_dir %s does not exist. The lentiMPRA deposit is expected at data/mpra_human; "
            "fetch it with `uv run python scripts/fetch_lentimpra.py --select tables`, which "
            "verifies every byte against the publisher's MD5. NOTE: lentimpra.DEFAULT_DIR points at "
            "D:/dna_decode_cache/mpra_human on a drive that may be detached -- pass data_dir "
            "EXPLICITLY rather than relying on it." % p)
    return p


def human_natural_set(cell: str, data_dir: Path | str, *, folds: set[int] | None = None,
                      limit: int | None = None, strip_adaptors: bool = False
                      ) -> tuple[list[str], dict]:
    """Return (sequences, stats) for one cell type.

    `folds` selects the publisher's fold ids; None takes all. `strip_adaptors` trims 15 bp from each
    end to yield the 200 nt element (the bar's SECONDARY geometry). `limit` truncates AFTER fold
    selection and is applied deterministically by record order, never by sampling, so a limited run
    is reproducible.
    """
    if cell not in CELL_TYPES:
        raise HumanNaturalError("unknown cell %r; have %s" % (cell, sorted(CELL_TYPES)))
    p = _require_dir(data_dir)

    try:
        records, st = load_cell(cell, p)
    except LentiMpraDataError as e:
        # Carry the loader's own reason: it is more specific than anything restated here, and for
        # WTC11 it names the join-rate bar rather than a generic failure.
        raise HumanNaturalError("lentiMPRA refused %s: %s" % (cell, e)) from e

    want = sorted({int(f) for f in folds}) if folds else None
    chosen = [r for r in records if want is None or int(r.fold) in want]
    if not chosen:
        raise HumanNaturalError(
            "no records for cell %s with folds=%s; available folds are %s"
            % (cell, want, sorted(st.fold_sizes)))

    expect = ELEMENT_BP if strip_adaptors else FULL_BP
    seqs: list[str] = []
    for r in chosen:
        s = r.seq.strip().upper()
        if len(s) != FULL_BP:
            raise HumanNaturalError(
                "record %s on cell %s is %d bp, expected the publisher's %d bp; a natural set with "
                "off-length rows would produce a distinguishability number that looks fine and "
                "measures something else" % (r.name, cell, len(s), FULL_BP))
        seqs.append(s[ADAPTOR_BP:FULL_BP - ADAPTOR_BP] if strip_adaptors else s)

    if limit is not None:
        seqs = seqs[:limit]

    # alphabet_and_length_ok returns RATES, not a boolean -- it has no `ok` key, and an earlier
    # draft of this module invented one, which made a clean 500-sequence set look malformed. Both
    # rates must be exactly 1.0: a natural set with any invalid or off-length row would produce a
    # distinguishability number that looks fine and measures something else.
    ok = alphabet_and_length_ok(seqs, expect)
    if ok["alphabet_valid_rate"] != 1.0 or ok["exact_length_rate"] != 1.0:
        raise HumanNaturalError(
            "cell %s at %d bp: alphabet_valid_rate=%s exact_length_rate=%s, both must be 1.0 (%s)"
            % (cell, expect, ok["alphabet_valid_rate"], ok["exact_length_rate"], ok))

    stats = {
        "cell": cell,
        "data_dir": str(p).replace("\\", "/"),
        "n": len(seqs),
        "bp": expect,
        "geometry": "element_200nt_adaptors_stripped" if strip_adaptors else "full_230nt_as_shipped",
        "folds_selected": want,
        "folds_available": sorted(st.fold_sizes),
        "n_folds_available": st.n_folds,
        "n_designed": st.n_designed,
        "n_measured": st.n_measured,
        "n_joined": st.n_joined,
        "join_rate": st.join_rate,
        "sequence_column": st.sequence_column,
        "limit_applied": limit,
        "alphabet_and_length": ok,
    }
    return seqs, stats


def fit_eval_by_fold(cell: str, data_dir: Path | str, *, eval_folds: set[int],
                     limit_fit: int | None = None, limit_eval: int | None = None,
                     strip_adaptors: bool = False
                     ) -> tuple[list[str], list[str], dict]:
    """Split on the PUBLISHER's fold assignment: `eval_folds` is held out, the rest is fit.

    Returns (fit, eval, stats). REFUSES if the two sides overlap or if either is empty -- a Markov
    null fitted and scored on the same sequences would report a number about memorisation.
    """
    _, all_stats = human_natural_set(cell, data_dir, strip_adaptors=strip_adaptors)
    available = set(all_stats["folds_available"])
    want_eval = {int(f) for f in eval_folds}
    missing = want_eval - available
    if missing:
        raise HumanNaturalError(
            "eval_folds %s are not in the publisher's folds %s" % (sorted(missing), sorted(available)))
    want_fit = available - want_eval
    if not want_fit:
        raise HumanNaturalError(
            "eval_folds %s consumes every available fold, leaving nothing to fit the null on"
            % sorted(want_eval))

    fit, fit_st = human_natural_set(cell, data_dir, folds=want_fit, limit=limit_fit,
                                    strip_adaptors=strip_adaptors)
    ev, ev_st = human_natural_set(cell, data_dir, folds=want_eval, limit=limit_eval,
                                  strip_adaptors=strip_adaptors)
    if not fit or not ev:
        raise HumanNaturalError("empty side after split: fit=%d eval=%d" % (len(fit), len(ev)))

    # Disjointness is guaranteed by the fold partition, but it is asserted on CONTENT because a
    # fold-id bug would otherwise be invisible and would inflate every downstream number.
    overlap = set(fit) & set(ev)
    if overlap:
        raise HumanNaturalError(
            "fit and eval share %d sequences despite a disjoint fold partition -- a fold-id bug"
            % len(overlap))

    stats = {
        "cell": cell,
        "geometry": fit_st["geometry"],
        "bp": fit_st["bp"],
        "fit_folds": sorted(want_fit),
        "eval_folds": sorted(want_eval),
        "n_fit": len(fit),
        "n_eval": len(ev),
        "split_is_the_publishers": True,
        "content_overlap": 0,
        "folds_available": sorted(available),
        "join_rate": fit_st["join_rate"],
    }
    return fit, ev, stats
