"""Load the Agarwal 2025 human lentiMPRA panel: 230 bp sequences with measured log2(RNA/DNA).

SUBSTRATE, measured 2026-10-27 from the verified deposit (Zenodo 10558183, CC-BY-4.0):

| cell  | designed (S3) | measured (S6) | JOINED        | usable |
|-------|---------------|---------------|---------------|--------|
| HepG2 | 164,307       | 61,463        | 61,463 (1.000)| yes    |
| K562  | 243,778       | 98,332        | 98,332 (1.000)| yes    |
| WTC11 |  75,542       | 23,064        | 23,064 (1.000)| yes    |

So the usable labelled set is **182,859 pairs** over three cell types, every sequence exactly 230 bp,
with the authors' own **10-fold assignment** carried through.

WHY WTC11 NEEDED A NAME RULE, corrected 2026-11-23 (it previously joined at 0.000 and refused, on a
MISSTATED mechanism -- the old text said S6 keys this cell by Ensembl gene id while S3 uses element
names, "different entity types, so the intersection is empty"). Measured, with the two shipped cells
reproducing 1.0000 as controls first:

  * **S6 carries no orientation marker for ANY cell.** Every measured id is bare -- HepG2
    `DNasePeakNoPromoter1`, K562 `peak1`, WTC11 `seq1000` / `ENSG00000000003`.
  * **The difference is on the DESIGN side.** S3 HepG2 and S3 K562 name the forward strand BARE and
    put the reverse under `_Reversed:`, so a bare S6 id matches at identity. **S3 WTC11 is the only
    sheet that suffixes the forward strand too** (`seq1000_F` / `seq1000_R`, neither bare), so
    nothing matches and the rate is exactly 0.
  * **ENSG was a red herring.** It is 4,397 of WTC11's 23,064 measured ids (19%); the other 81% are
    bare `seq...` ids failing for the identical suffix reason. The old "different entity types"
    account could not have explained those, and S3 WTC11 in fact holds 15,000 ENSG-named design rows
    (7,500 genes x 2 orientations), so the two sides were never disjoint entity types at all.

WHICH STRAND A BARE S6 ID DENOTES -- settled on SEQUENCE IDENTITY, not on the naming analogy, because
`_F` and `_R` are different sequences for **23,064 of 23,064** measured rows while the join rate
reads 1.000 either way. That is the shape where a cohort looks complete and is silently half wrong.
Three independent lines agree on `_F`:

  1. The 4,397 bare-ENSG ids WTC11 measured also exist as BARE names in the HepG2/K562 design
     sheets, which the shipped cells pair by identity. That shared-promoter sequence equals WTC11's
     `_F` for **4,397 of 4,397** and `_R` for **0** (0 neither, 0 absent).
  2. The `_F`/`_R` inserts (`seq[15:215]`, the adaptors being constant) are exact reverse
     complements for **37,621 of 37,621** pairs -- so they are true orientation pairs of one 200 nt
     element, not unrelated oligos.
  3. Table S7, the authors' own joint library, takes **18,946 `_F` rows and zero `_R`**.

The rule is therefore `identity, else sid + "_F"`, applied as an EXPLICIT two-step lookup. It is not
a stripped-suffix pool: stripping `_F`/`_R` into one dict makes `seq1000_F` and `seq1000_R` collide
on one key and lets insertion order choose the sequence, which would score the same 1.0000 while
pairing half the cohort with the wrong strand. `load_cell` REFUSES if an id resolves both ways to
different sequences, and `LentiStats.name_rule` records which rule was used.

WHAT THIS STILL DOES NOT ESTABLISH. S6 holds one row per element, which does not by itself separate
"the forward-strand measurement" from "a value aggregated over both orientations". The pairing above
is the construct convention the two shipped cells already use, so it is consistent across all three
cells; but until the paper's Methods or the deposit's model code states how the S6 value was built,
say "one value per element, paired with the forward-designed oligo, as for HepG2/K562" and do NOT
write that S6 *is* the forward measurement.

The join rate stays ENFORCED at MIN_JOIN_RATE. This repo has paid for a plausible-looking partial
cohort before (the sigma70 rung's naive key join returned 3,821 of 10,898, which is why
`expression.load_pairs` enforces it), and Table S7 rescues only **699 of 23,064 (3.03%)** -- the
route a previous correction reached for, and 3% of a cohort looks exactly like a cohort.

TWO PARSING TRAPS, both measured rather than assumed:

1. **The header row MOVES between sheets of the same workbook** -- S3 has it at row 4 (HepG2),
   row 2 (K562) and row 1 (WTC11), because some sheets carry a free-text design note above it. A
   fixed `min_row` reads a 327-character description as a column name and then finds no columns. So
   the header is DISCOVERED by looking for the tokens, and the row it was found at is REPORTED.
2. **The 230 nt includes the cloning adapters** -- the column is named
   `230nt sequence (15nt 5' adaptor - ...)`. The biological insert is shorter than 230. This module
   returns the sequence AS PUBLISHED and records the verbatim column name in `LentiStats`, so a
   caller trimming adapters does so deliberately. It does not trim silently, because a silent trim
   would change what every downstream number is about.

PUBLISHED INCUMBENTS live in the same workbook (`Benchmarking results`): per-fold figures for
Biochemical / MPRAnn / Sei / Enformer / MPRALegNet. **The sheet does not state the METRIC**, so
`published_incumbents()` returns the numbers with `metric="UNSPECIFIED_IN_SOURCE"` and callers must
not compare a Spearman of their own against them until that is resolved -- a readout-type mismatch
is not a comparison.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

#: The deposit lands here, namespace-separate from the E. coli `mpra/` directory so a human table
#: can never be mistaken for a GSE144621 one (the shared-key silent-overwrite trap).
DEFAULT_DIR = Path("D:/dna_decode_cache/mpra_human")

S3_FILE = "Table S3 - large-scale lib design.xlsx"
S6_FILE = "Table S6 - folds and performance.xlsx"

#: Mirrors `expression.load_pairs`. A cohort below this is a key problem, not a small cohort.
MIN_JOIN_RATE = 0.95

#: The strand a bare S6 id denotes, for sheets that suffix BOTH strands (only S3 WTC11 does). See
#: the module docstring for the three measurements that settle this on `_F`. It is a module constant
#: rather than a literal at the call site so the value, and the reasoning it points at, has exactly
#: one home -- and so a future sheet that inverts the convention is a one-line change with a
#: docstring beside it rather than a silent edit inside a join loop.
ORIENTATION_SUFFIX = "_F"

#: Name rules `load_cell` may use, reported in `LentiStats.name_rule`. The pairing a cohort was
#: built with changes what every downstream number is ABOUT, so it is recorded, never implicit --
#: the same discipline as `sequence_column` recording that the 230 nt includes the adaptors.
RULE_IDENTITY = "identity"
RULE_IDENTITY_THEN_FORWARD = "identity_then_forward_suffix"

#: Measured: all 483,627 designed sequences across all three sheets are exactly this long.
EXPECTED_LEN = 230

CELL_TYPES = ("HepG2", "K562", "WTC11")

#: Row kinds of the `Benchmarking results` sheet, as the parser sees them. `ROW_SKIPPED` is the
#: interesting one: those rows are invisible to `published_incumbents` and visible to a human, which
#: is why a claim that "the sheet never states a metric" cannot be made from the parsed subset.
ROW_CELL_TYPE = "cell-type"
ROW_METHOD_HEADER = "method-header"
ROW_FOLD = "fold"
ROW_SKIPPED = "skipped"


class LentiMpraDataError(ValueError):
    """Refuse rather than return a plausible-looking partial cohort."""


#: Parsed-workbook cache, keyed by (kind, cell, data_dir). Table S3 is 56 MB and openpyxl takes
#: tens of seconds per read_only pass over it, so an uncached loader makes any caller that touches
#: two cell types -- or a test module that touches several -- pay minutes. Measured: the first
#: version of `tests/test_lentimpra.py` exceeded a 600 s timeout purely by re-parsing.
#:
#: Values are stored as TUPLES so a cached payload cannot be mutated by one caller and observed
#: changed by the next; each public call rebuilds its own dict from them.
_CACHE: dict = {}


def clear_cache() -> None:
    """Drop the parsed-workbook cache. For tests that swap the data directory under the loader."""
    _CACHE.clear()


@dataclass(frozen=True)
class LentiRecord:
    name: str
    seq: str
    activity: float          # observed log2(RNA/DNA)
    fold: int                # the authors' published fold, 1..10
    cell: str


@dataclass
class LentiStats:
    cell: str = ""
    n_designed: int = 0
    n_measured: int = 0
    n_joined: int = 0
    join_rate: float = 0.0
    s3_sheet: str = ""
    s3_header_row: int = 0
    s6_sheet: str = ""
    s6_header_row: int = 0
    sequence_column: str = ""     # verbatim, because it records that adapters are included
    length_min: int = 0
    length_max: int = 0
    n_folds: int = 0
    name_rule: str = RULE_IDENTITY   # which key rule joined this cell; see ORIENTATION_SUFFIX
    n_joined_identity: int = 0       # rows paired on a bare id, as the two shipped cells are
    n_joined_forward: int = 0        # rows paired via `sid + ORIENTATION_SUFFIX`
    fold_sizes: dict = field(default_factory=dict)

    def as_dict(self) -> dict:
        d = {k: getattr(self, k) for k in (
            "cell", "n_designed", "n_measured", "n_joined", "s3_sheet", "s3_header_row",
            "s6_sheet", "s6_header_row", "sequence_column", "length_min", "length_max", "n_folds",
            "name_rule", "n_joined_identity", "n_joined_forward")}
        d["join_rate"] = round(self.join_rate, 4)
        d["fold_sizes"] = dict(sorted(self.fold_sizes.items()))
        return d


def _find_header(ws, tokens, limit: int = 8):
    """First row (1-based) containing every token, plus its cells. The row MOVES; see module doc."""
    for i, row in enumerate(ws.iter_rows(min_row=1, max_row=limit, values_only=True), start=1):
        cells = ["" if c is None else str(c).strip() for c in row]
        low = [c.lower() for c in cells]
        if all(any(tok in c for c in low) for tok in tokens):
            return i, cells
    raise LentiMpraDataError(
        "no header row containing %s in the first %d rows; the sheet layout changed and a fixed "
        "row index would read a description string as a column name" % (list(tokens), limit))


def _sheet_for(wb, cell: str, *, suffix: bool):
    want = cell.lower()
    for s in wb.sheetnames:
        sl = s.lower()
        if (sl.endswith(want) if suffix else sl.startswith(want)):
            return s
    raise LentiMpraDataError("no sheet for cell %r among %s" % (cell, wb.sheetnames))


def load_design(cell: str, data_dir: Path | str = DEFAULT_DIR) -> tuple[dict, LentiStats]:
    """`{name: seq}` from Table S3, for one cell type."""
    import openpyxl

    st = LentiStats(cell=cell)
    path = Path(data_dir) / S3_FILE
    if not path.is_file():
        raise LentiMpraDataError(
            "%s is absent. Fetch it with `uv run python scripts/fetch_lentimpra.py "
            "--select tables`, which verifies every byte against the publisher's MD5." % path)

    key = ("design", cell, str(path))
    if key in _CACHE:
        pairs, sheet, hrow, seq_col, lmin, lmax = _CACHE[key]
        st.s3_sheet, st.s3_header_row, st.sequence_column = sheet, hrow, seq_col
        st.n_designed, st.length_min, st.length_max = len(pairs), lmin, lmax
        return dict(pairs), st

    wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
    try:
        st.s3_sheet = _sheet_for(wb, cell, suffix=False)
        ws = wb[st.s3_sheet]
        st.s3_header_row, hdr = _find_header(ws, ("name", "sequence"))
        j_name = next(j for j, h in enumerate(hdr) if h.lower() == "name")
        j_seq = next(j for j, h in enumerate(hdr) if "sequence" in h.lower())
        st.sequence_column = hdr[j_seq]
        out, lens = {}, []
        for row in ws.iter_rows(min_row=st.s3_header_row + 1, values_only=True):
            if j_name >= len(row) or j_seq >= len(row):
                continue
            nm, sq = row[j_name], row[j_seq]
            if nm is None or sq is None:
                continue
            nm, sq = str(nm).strip(), str(sq).strip().upper()
            if not nm or not re.fullmatch(r"[ACGTN]+", sq):
                continue
            out[nm] = sq
            lens.append(len(sq))
        if not out:
            raise LentiMpraDataError(
                "%s/%s parsed ZERO sequences. A zero-row parse that is then joined reads as 'no "
                "overlap' rather than 'no data', so this refuses instead." % (path.name, st.s3_sheet))
        st.n_designed = len(out)
        st.length_min, st.length_max = min(lens), max(lens)
        _CACHE[key] = (tuple(out.items()), st.s3_sheet, st.s3_header_row,
                       st.sequence_column, st.length_min, st.length_max)
        return out, st
    finally:
        wb.close()


def load_activity(cell: str, data_dir: Path | str = DEFAULT_DIR) -> tuple[dict, LentiStats]:
    """`{sequence_id: (activity, fold)}` from Table S6, for one cell type."""
    import openpyxl

    st = LentiStats(cell=cell)
    path = Path(data_dir) / S6_FILE
    if not path.is_file():
        raise LentiMpraDataError("%s is absent; see scripts/fetch_lentimpra.py" % path)

    key = ("activity", cell, str(path))
    if key in _CACHE:
        pairs, sheet, hrow = _CACHE[key]
        st.s6_sheet, st.s6_header_row, st.n_measured = sheet, hrow, len(pairs)
        return dict(pairs), st

    wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
    try:
        st.s6_sheet = _sheet_for(wb, cell, suffix=True)
        ws = wb[st.s6_sheet]
        st.s6_header_row, hdr = _find_header(ws, ("sequence id", "observed"))
        low = [h.lower() for h in hdr]
        j_id = next(j for j, h in enumerate(low) if "sequence id" in h)
        j_y = next(j for j, h in enumerate(low) if "observed" in h)
        j_f = next((j for j, h in enumerate(low) if h == "fold"), None)
        if j_f is None:
            raise LentiMpraDataError(
                "%s has no `Fold` column; the published fold assignment is what makes this cohort "
                "comparable to the published benchmark, so it is required. NOTE the scope: the "
                "publishers assign folds RANDOMLY OVER ELEMENTS, not by chromosome, so the split "
                "is leakage-safe at SEQUENCE IDENTITY and NOT at coordinate overlap -- measured "
                "2026-11-23, 625 genomically overlapping pairs straddle a fold-role boundary "
                "across the three cells. Do not call these folds leakage-safe unqualified."
                % st.s6_sheet)
        out = {}
        for row in ws.iter_rows(min_row=st.s6_header_row + 1, values_only=True):
            if max(j_id, j_y, j_f) >= len(row):
                continue
            sid, y, f = row[j_id], row[j_y], row[j_f]
            if sid is None or y is None or f is None:
                continue
            try:
                out[str(sid).strip()] = (float(y), int(f))
            except (TypeError, ValueError):
                continue
        if not out:
            raise LentiMpraDataError("%s parsed ZERO activities" % st.s6_sheet)
        st.n_measured = len(out)
        _CACHE[key] = (tuple(out.items()), st.s6_sheet, st.s6_header_row)
        return out, st
    finally:
        wb.close()


def load_cell(cell: str, data_dir: Path | str = DEFAULT_DIR, *,
              min_join_rate: float = MIN_JOIN_RATE) -> tuple[list, LentiStats]:
    """Join design to activity for one cell type. REFUSES below `min_join_rate`.

    The join rate is computed over the MEASURED side: the design library is deliberately larger than
    the measured set (not every designed element yields a usable readout), so measured-coverage is
    the meaningful denominator and design-coverage is reported but never gated on.
    """
    if cell not in CELL_TYPES:
        raise LentiMpraDataError("unknown cell %r; have %s" % (cell, list(CELL_TYPES)))
    seqs, st3 = load_design(cell, data_dir)
    acts, st6 = load_activity(cell, data_dir)

    st = LentiStats(
        cell=cell, n_designed=st3.n_designed, n_measured=st6.n_measured,
        s3_sheet=st3.s3_sheet, s3_header_row=st3.s3_header_row,
        s6_sheet=st6.s6_sheet, s6_header_row=st6.s6_header_row,
        sequence_column=st3.sequence_column,
        length_min=st3.length_min, length_max=st3.length_max,
    )

    records, lens = [], []
    for sid, (y, fold) in acts.items():
        # Identity FIRST, so the two cells that already join at 1.0000 take an untouched path and
        # the suffix rule can only ever ADD rows. Then the explicit forward-strand name. These are
        # two named lookups rather than one lookup into a suffix-stripped pool, because stripping
        # collides `X_F` and `X_R` onto one key and lets insertion order pick the sequence -- which
        # scores an identical join rate while pairing half the cohort with the wrong strand.
        sq = seqs.get(sid)
        alt = seqs.get(sid + ORIENTATION_SUFFIX)
        if sq is not None and alt is not None and sq != alt:
            raise LentiMpraDataError(
                "%s: id %r resolves BOTH as a bare name and as %r, to DIFFERENT sequences. The "
                "name rule cannot decide which construct was measured, and either choice would "
                "look like a complete join, so this refuses." % (cell, sid, sid + ORIENTATION_SUFFIX))
        if sq is not None:
            st.n_joined_identity += 1
        elif alt is not None:
            sq = alt
            st.n_joined_forward += 1
        else:
            continue
        records.append(LentiRecord(name=sid, seq=sq, activity=y, fold=fold, cell=cell))
        lens.append(len(sq))

    st.name_rule = RULE_IDENTITY_THEN_FORWARD if st.n_joined_forward else RULE_IDENTITY
    st.n_joined = len(records)
    st.join_rate = len(records) / st.n_measured if st.n_measured else 0.0
    for r in records:
        st.fold_sizes[r.fold] = st.fold_sizes.get(r.fold, 0) + 1
    st.n_folds = len(st.fold_sizes)
    if lens:
        st.length_min, st.length_max = min(lens), max(lens)

    if st.join_rate < min_join_rate:
        raise LentiMpraDataError(
            "%s joins at %.4f (%d of %d measured rows), below the %.2f bar. This is a KEY problem, "
            "not a small cohort. S6 ids are bare for every cell; a sheet joins at identity only if "
            "it names the forward strand bare, and `identity, else sid + %r` is the only rule "
            "applied -- so a rate this low means this sheet uses a THIRD naming convention that "
            "neither rule covers. Diagnose the convention; do not relax the bar, because returning "
            "the overlap would hand back a %.1f%%-complete cohort that looks complete."
            % (cell, st.join_rate, st.n_joined, st.n_measured, min_join_rate, ORIENTATION_SUFFIX,
               100 * st.join_rate))

    if st.length_min != st.length_max:
        raise LentiMpraDataError(
            "ragged sequences (%d-%d bp) where the panel is documented uniform; a downstream "
            "fixed-width model would silently pad or truncate" % (st.length_min, st.length_max))
    if st.length_max != EXPECTED_LEN:
        raise LentiMpraDataError(
            "sequences are %d bp, expected %d. The deposit may have changed; refusing rather than "
            "proceeding on a different substrate." % (st.length_max, EXPECTED_LEN))
    return records, st


def usable_cells(data_dir: Path | str = DEFAULT_DIR) -> dict:
    """Which cell types clear the join bar. Reports the failures rather than hiding them."""
    out = {}
    for cell in CELL_TYPES:
        try:
            recs, st = load_cell(cell, data_dir)
            out[cell] = {"usable": True, "n": len(recs), "stats": st.as_dict()}
        except LentiMpraDataError as exc:
            out[cell] = {"usable": False, "reason": str(exc)[:220]}
    return out


def _walk_rows(ws):
    """Yield `(index, kind, first, cells, cur, methods)` for EVERY row of a benchmarking sheet.

    THE SINGLE DEFINITION OF THE CONSUMPTION RULE. A row's kind is not a function of the row alone:
    `cur` and `methods` are carried across rows, the cell-type branch resets `methods`, and the fold
    branch is reachable only once both are set. So the state transitions live here, once, and both
    `published_incumbents` and `unparsed_rows` are pure consumers of this stream.

    Extracting only a predicate would have left each caller carrying its own copy of those
    transitions -- the same duplication one level down, and able to drift in exactly the way that
    makes one consumer disagree with the other about which rows are discarded.

    `index` is 1-based to match what a spreadsheet shows a human.
    """
    cur, methods = None, None
    for index, row in enumerate(ws.iter_rows(values_only=True), start=1):
        cells = list(row)
        first = "" if cells[0] is None else str(cells[0]).strip()
        if first in CELL_TYPES:
            cur, methods = first, None
            yield index, ROW_CELL_TYPE, first, cells, cur, methods
            continue
        if first.lower() == "test fold":
            methods = ["" if c is None else str(c).strip() for c in cells[1:]]
            yield index, ROW_METHOD_HEADER, first, cells, cur, methods
            continue
        if cur and methods and isinstance(cells[0], (int, float)):
            yield index, ROW_FOLD, first, cells, cur, methods
            continue
        yield index, ROW_SKIPPED, first, cells, cur, methods


def unparsed_rows(data_dir: Path | str = DEFAULT_DIR,
                  sheet: str = "Benchmarking results") -> list[dict]:
    """Every row `published_incumbents` does NOT consume, as `{"row": <1-based>, "cells": [...]}`.

    This exists because the refusal's own next check was unrunnable as prose. The metric, if the
    workbook states it anywhere, is in a caption or units row -- and those are the rows the parser
    throws away without a record, so no amount of reading the parsed output can find them.

    Raises `LentiMpraDataError` when the workbook is absent, matching `published_incumbents` rather
    than returning `[]`: an empty list would read as "there are no unparsed rows" when the true
    state is "could not look", and that is the confusion this whole correction is about.
    """
    import openpyxl

    path = Path(data_dir) / S6_FILE
    if not path.is_file():
        raise LentiMpraDataError("%s is absent; see scripts/fetch_lentimpra.py" % path)
    wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
    try:
        ws = wb[sheet]
        out = []
        for index, kind, _first, cells, _cur, _methods in _walk_rows(ws):
            if kind != ROW_SKIPPED:
                continue
            values = [c for c in cells if c is not None and str(c).strip() != ""]
            if values:
                out.append({"row": index, "cells": values})
        return out
    finally:
        wb.close()


def published_incumbents(data_dir: Path | str = DEFAULT_DIR) -> dict:
    """Per-fold published figures for the five benchmarked methods, from Table S6.

    `metric` is **UNSPECIFIED_IN_SOURCE**, and the reason is a READOUT-TYPE MISMATCH: comparing a
    Spearman of one's own against an unknown readout is not a comparison, so the field refuses
    instead of guessing the likely answer. That rationale is unaffected by everything below.

    SCOPE OF THE SILENCE, stated precisely because the earlier wording over-claimed it. The rows
    this function CONSUMES are three shapes only -- a cell-type marker, a `test fold` header
    carrying method names, and numeric fold rows -- and no row it reads states a coefficient. Every
    other row is skipped with no record, so the SHEET has not been shown to be silent; only the
    parsed subset has. A caption or units row is exactly where such a statement would sit, and it
    is exactly what this parser discards. Use `unparsed_rows()` to look at them.
    """
    import openpyxl

    path = Path(data_dir) / S6_FILE
    if not path.is_file():
        raise LentiMpraDataError("%s is absent; see scripts/fetch_lentimpra.py" % path)
    wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
    try:
        ws = wb["Benchmarking results"]
        blocks = {}
        # The consumption rule lives in `_walk_rows` and nowhere else; this is a pure consumer of
        # it, which is what keeps `unparsed_rows` exactly complementary to this function.
        for _index, kind, first, cells, cur, methods in _walk_rows(ws):
            if kind == ROW_CELL_TYPE:
                blocks[first] = {}
            elif kind == ROW_FOLD:
                fold = int(cells[0])
                for m, v in zip(methods, cells[1:]):
                    if m and isinstance(v, (int, float)):
                        blocks[cur].setdefault(m, {})[fold] = float(v)
        return {"metric": "UNSPECIFIED_IN_SOURCE", "source": S6_FILE,
                "sheet": "Benchmarking results", "by_cell": blocks,
                # Resolution ATTEMPTED 2026-10-28 against the paper's own text and REFUSED. Recorded
                # here so the next session does not repeat the search and does not mistake the
                # corroboration below for the answer. Memo: wiki/lentimpra_metric_resolution_2026-10-28.md
                "metric_resolution_attempt": {
                    "verdict": "UNRESOLVED_TABLE_IDENTITY_UNVERIFIED",
                    "paper_reports_both": (
                        "the publication's MPRALegNet panels carry 'corresponding Pearson (r) and "
                        "Spearman (rho) correlation values', so the paper has no single global "
                        "benchmarking metric to copy"
                    ),
                    "table_identity_mismatch": (
                        "the paper's IN-TEXT 'Supplementary Table 6' is a biochemical-feature / "
                        "bigWig manifest (it is cited for '1,506 bigWig files for HepG2 cells'), "
                        "whereas this S6_FILE is 'Table S6 - folds and performance.xlsx' carrying a "
                        "'Benchmarking results' sheet the paper never names. Do NOT cite 'the "
                        "paper's Table S6' as the provenance of these numbers."
                    ),
                    "corroboration_NOT_the_answer": (
                        "the paper states MPRALegNet reached 'Pearson r ~ 0.83' COMBINING folds, and "
                        "the HepG2 MPRALegNet fold-mean here is 0.8333. Suggestive of Pearson; it is "
                        "a different aggregation and says nothing about this sheet."
                    ),
                    "decisive_next_check": (
                        "call unparsed_rows() and read the rows this parser DISCARDS. The header "
                        "row cannot settle it: the header is 'test fold | <method names>' and "
                        "carries no coefficient, so the earlier instruction to read it could never "
                        "have worked. A caption or units row is where a metric statement would "
                        "sit, and those are precisely the skipped rows. Running it needs the "
                        "substrate volume attached; that is a prerequisite, not the reason the "
                        "previous check was wrong."
                    ),
                }}
    finally:
        wb.close()
