"""Genome-wide MPRA: the oracle's honest-limit substrate, and the first NOISE CEILING we can derive.

GSE108535 (the designed sigma70 grid) gave the oracle rho ~0.59 on held-out elements, with the standing
caveat that *a fully-crossed designed grid is not a natural promoter distribution*. **GSE144621 is the
substrate that tests it**: 321,123 sheared genomic fragments from *E. coli* MG1655 measured in **two media
(LB and M9)**, plus 46,713 peak tiles.

Every fact below was MEASURED on the real files, not taken from the paper or from the sibling dataset.

**FOUR TRAPS, and two of them contradict what the GSE108535 loader taught.**

1. **THE DELIMITER IS MIXED WITHIN ONE GEO SERIES.** The two `frag-rLP5_*_expression` files are
   SPACE-delimited; `peak_tile_expression_formatted_std` is TAB-delimited. Assuming either one globally
   parses the other to ZERO rows, and a zero-row parse that is then averaged reads as "no signal" rather
   than "no data". So the delimiter is **sniffed per file** and the parse REFUSES on an implausible column
   count rather than silently yielding nothing.
2. **`RNA_exp_ave` IS ALREADY DNA-NORMALISED — do NOT divide by `DNA_ave`.** Measured:
   `spearman(RNA_exp_ave, DNA_ave) = -0.039`, i.e. no relationship. Raw RNA counts would track DNA
   abundance strongly. Dividing again would double-normalise and inject the DNA column's noise.
3. **A RANDOM SPLIT LEAKS, and this is measured not assumed: 91.5% of consecutive fragments OVERLAP** (a
   sheared library tiles the genome). Two fragments sharing 200 bp are near-duplicates, so the honest split
   is CONTIGUOUS GENOMIC BLOCKS — the bacterial analogue of a held-out chromosome. `position_blocked_split`
   is the only split offered here; there is deliberately no random option to reach for.
4. **FRAGMENTS ARE VARIABLE LENGTH (48-475 bp, median 244, 280 distinct lengths in a 40k sample).** The
   designed grid was a fixed 150 bp, so the oracle's winning `positional_onehot` feature **does not apply
   without choosing an anchor**, and there is no obviously correct anchor (the promoter sits somewhere
   inside the fragment). k-mer features are length-invariant and are therefore the honest default here;
   a positional feature needs a defended anchor and is NOT provided rather than silently anchored at 0.

**THE NOISE CEILING — the thing GSE108535 could not give us.** This assay ships two replicates, so the
assay's own reproducibility is computable, and it BOUNDS any model: a sequence model cannot correlate with
a measurement better than the measurement correlates with itself. Measured on 297,599 shared fragments:
single-replicate agreement **0.4578 (LB) / 0.4774 (M9)**; by Spearman-Brown the 2-replicate mean has
reliability **0.628 / 0.646**, so the ceiling for predicting `RNA_exp_ave` from sequence is
**sqrt(0.628) = 0.7925**. Read every number here against 0.79, never against 1.0. (GSE108535's two
published tables correlate at exactly 1.0 — the same numbers under a monotone transform, not replicates —
so no ceiling was derivable there and its 0.59 stands uncorrected.)
"""
from __future__ import annotations

import gzip
import math
from dataclasses import dataclass, field
from pathlib import Path

GEO_BASE = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE144nnn/GSE144621/suppl"
FRAG_FILES = {
    "LB": "GSE144621_U00096.2_frag-rLP5_LB_expression.txt.gz",
    "M9": "GSE144621_U00096.2_frag-rLP5_M9_expression.txt.gz",
}
PEAK_TILE_FILE = "GSE144621_peak_tile_expression_formatted_std.txt.gz"

DEFAULT_CACHE = Path("D:/dna_decode_cache/mpra")

#: Measured column counts. A parse yielding fewer means the delimiter guess was wrong.
MIN_COLUMNS = 8


class GenomeWideDataError(RuntimeError):
    """The genome-wide MPRA substrate could not be loaded in a trustworthy state."""


@dataclass
class Fragment:
    """One measured genomic fragment."""

    seq: str
    expression: float          # RNA_exp_ave -- ALREADY DNA-normalised (trap 2)
    rep1: float
    rep2: float
    start: int
    end: int
    strand: str


@dataclass
class LoadStats:
    """What the parse produced. Reported so a degraded load is visible rather than inferred."""

    path: str = ""
    delimiter: str = ""
    n_header_columns: int = 0
    n_lines: int = 0
    n_parsed: int = 0
    n_field_count_mismatch: int = 0
    n_non_acgt: int = 0
    n_unparseable_value: int = 0
    length_min: int = 0
    length_median: int = 0
    length_max: int = 0
    n_distinct_lengths: int = 0

    def as_dict(self) -> dict:
        return {k: getattr(self, k) for k in (
            "path", "delimiter", "n_header_columns", "n_lines", "n_parsed",
            "n_field_count_mismatch", "n_non_acgt", "n_unparseable_value",
            "length_min", "length_median", "length_max", "n_distinct_lengths")}


def sniff_delimiter(first_line: str) -> str:
    """Decide tab vs whitespace from the HEADER, per file.

    Load-bearing: within GSE144621 the frag files are space-delimited and the peak-tile file is
    tab-delimited. A global assumption parses one of them to zero rows.
    """
    return "\t" if "\t" in first_line else None      # None => str.split() on any whitespace


def _split(line: str, delim: str | None) -> list[str]:
    return line.rstrip("\n").split(delim) if delim else line.split()


def ensure_cached(cache_dir: Path | str = DEFAULT_CACHE, *, which: str = "LB",
                  allow_download: bool = True) -> Path:
    """Make sure one GSE144621 file is present; download it if not. Returns the file path."""
    name = FRAG_FILES.get(which, PEAK_TILE_FILE if which == "peak_tile" else None)
    if name is None:
        raise GenomeWideDataError(f"unknown file key {which!r}; have {sorted(FRAG_FILES)} + 'peak_tile'")
    d = Path(cache_dir)
    d.mkdir(parents=True, exist_ok=True)
    p = d / name
    if not p.exists():
        if not allow_download:
            raise GenomeWideDataError(f"{p} missing and downloading is disabled")
        import urllib.request
        tmp = p.with_suffix(p.suffix + ".part")
        urllib.request.urlretrieve(f"{GEO_BASE}/{name}", tmp)   # noqa: S310 - fixed NCBI URL
        tmp.replace(p)
    return p


def load_fragments(cache_dir: Path | str = DEFAULT_CACHE, *, which: str = "LB",
                   allow_download: bool = True, limit: int | None = None,
                   ) -> tuple[list[Fragment], LoadStats]:
    """Load (sequence, expression, both replicates, coordinates) for one medium. Fails loudly."""
    p = ensure_cached(cache_dir, which=which, allow_download=allow_download)
    st = LoadStats(path=str(p))
    out: list[Fragment] = []
    lengths: list[int] = []
    with gzip.open(p, "rt", errors="replace") as fh:
        first = next(fh)
        delim = sniff_delimiter(first)
        st.delimiter = "tab" if delim == "\t" else "whitespace"
        hdr = _split(first, delim)
        st.n_header_columns = len(hdr)
        if len(hdr) < MIN_COLUMNS:
            raise GenomeWideDataError(
                f"header parsed to {len(hdr)} columns under {st.delimiter} delimiter; expected "
                f">= {MIN_COLUMNS}. The delimiter sniff is wrong for {p.name} -- fix sniff_delimiter "
                f"rather than lowering MIN_COLUMNS."
            )
        ix = {c: i for i, c in enumerate(hdr)}
        need = ("fragment", "RNA_exp_ave", "RNA_exp_1", "RNA_exp_2", "start", "end", "strand")
        missing = [c for c in need if c not in ix]
        if missing:
            raise GenomeWideDataError(f"{p.name} lacks required columns {missing}; has {hdr}")
        for line in fh:
            st.n_lines += 1
            q = _split(line, delim)
            if len(q) != len(hdr):
                st.n_field_count_mismatch += 1
                continue
            seq = q[ix["fragment"]].upper()
            if set(seq) - set("ACGT"):
                st.n_non_acgt += 1
                continue
            try:
                out.append(Fragment(
                    seq=seq,
                    expression=float(q[ix["RNA_exp_ave"]]),
                    rep1=float(q[ix["RNA_exp_1"]]),
                    rep2=float(q[ix["RNA_exp_2"]]),
                    start=int(q[ix["start"]]), end=int(q[ix["end"]]),
                    strand=q[ix["strand"]],
                ))
            except ValueError:
                st.n_unparseable_value += 1
                continue
            lengths.append(len(seq))
            if limit and len(out) >= limit:
                break
    st.n_parsed = len(out)
    if not out:
        raise GenomeWideDataError(
            f"{p.name} parsed ZERO fragments under the {st.delimiter} delimiter. A zero-row parse that is "
            f"then averaged reads as 'no signal' rather than 'no data', so this refuses instead."
        )
    lengths.sort()
    st.length_min, st.length_max = lengths[0], lengths[-1]
    st.length_median = lengths[len(lengths) // 2]
    st.n_distinct_lengths = len(set(lengths))
    return out, st


# ---------------------------------------------------------------------------------------------------
# the noise ceiling -- what GSE108535 could not give us
# ---------------------------------------------------------------------------------------------------
def spearman(a: list[float], b: list[float]) -> float:
    from scipy.stats import spearmanr
    r = spearmanr(a, b).statistic
    return float(r) if r == r else 0.0


def noise_ceiling(frags: list[Fragment]) -> dict:
    """Replicate reproducibility, and the ceiling it imposes on ANY sequence model.

    A model cannot correlate with a measurement better than the measurement correlates with itself.
    `RNA_exp_ave` is the mean of two replicates, so its reliability is the Spearman-Brown step-up of the
    single-replicate agreement, and the ceiling on predicting it is the square root of that.
    """
    r11 = spearman([f.rep1 for f in frags], [f.rep2 for f in frags])
    rel = 2 * r11 / (1 + r11) if r11 > -1 else 0.0
    return {
        "n": len(frags),
        "single_replicate_agreement": round(r11, 4),
        "reliability_of_2rep_mean_spearman_brown": round(rel, 4),
        "ceiling_for_a_sequence_model": round(math.sqrt(rel), 4) if rel > 0 else 0.0,
        "note": ("read every model number against the ceiling, never against 1.0; GSE108535 offered no "
                 "ceiling because its two published tables correlate at exactly 1.0 (one monotone "
                 "transform of the other, not replicates)"),
    }


def condition_effect(a: list[Fragment], b: list[Fragment], *,
                     coordinate_key: bool = True) -> dict:
    """Does GROWTH CONDITION carry signal beyond noise? The bacterial analogue of cell-state conditioning.

    **The comparison must be NOISE-MATCHED, and getting that wrong inverts the answer.** Comparing
    `RNA_exp_ave` across media is average-vs-average, while replicate agreement is single-vs-single;
    averaging two replicates cuts noise by ~sqrt(2), so the cross-condition figure can come out HIGHER than
    within-condition reproducibility and appear to show that condition does not matter. Measured: 0.5025
    (ave-vs-ave) against 0.458/0.477 (single-vs-single) -- the naive reading is exactly backwards.

    So this reports BOTH: the single-vs-single cross-condition mean (matched to within-condition) and the
    disattenuated true correlation.
    """
    # KEY ON COORDINATES, not sequence alone. Measured on the real files: 340 LB and 302 M9 sequences appear
    # more than once, and EVERY one of them maps to more than one distinct coordinate triple (e.g. a 232-bp
    # sequence at both 211303 and 2604594 -- a genuine repeated genomic element). A sequence-only key
    # collapses those to one arbitrary record.
    #
    # HONEST MAGNITUDE: this does NOT move the published condition numbers. Re-derived under both keys, the
    # noise-matched gap is -0.0614 either way (delta -0.00002; n 297,599 -> 297,868). The reason to fix it is
    # PROSPECTIVE: a collapsed repeated element is block-assigned by whichever coordinate survived, so the two
    # loci can land on opposite sides of a position-blocked split. That is a split-integrity bug waiting to
    # happen, which is a better argument for the key than the published-number one.
    key = (lambda f: (f.start, f.end, f.strand, f.seq)) if coordinate_key else (lambda f: f.seq)
    ka: dict = {}
    kb: dict = {}
    collisions_a = collisions_b = 0
    for f in a:
        k = key(f)
        if k in ka:
            collisions_a += 1
        ka[k] = f
    for f in b:
        k = key(f)
        if k in kb:
            collisions_b += 1
        kb[k] = f
    shared = sorted(set(ka) & set(kb))
    if len(shared) < 100:
        return {"status": "insufficient_shared_fragments", "n_shared": len(shared)}
    A, B = [ka[s] for s in shared], [kb[s] for s in shared]

    within_a = spearman([f.rep1 for f in A], [f.rep2 for f in A])
    within_b = spearman([f.rep1 for f in B], [f.rep2 for f in B])
    cross = [spearman([getattr(f, x) for f in A], [getattr(f, y) for f in B])
             for x in ("rep1", "rep2") for y in ("rep1", "rep2")]
    cross_mean = sum(cross) / len(cross)
    within_mean = (within_a + within_b) / 2

    rel_a = 2 * within_a / (1 + within_a)
    rel_b = 2 * within_b / (1 + within_b)
    ave_ave = spearman([f.expression for f in A], [f.expression for f in B])
    irrelevant_ceiling = math.sqrt(rel_a * rel_b) if rel_a > 0 and rel_b > 0 else 0.0
    return {
        "n_shared": len(shared),
        "key": "coordinate" if coordinate_key else "sequence_only",
        "n_key_collisions_dropped": {"a": collisions_a, "b": collisions_b},
        "within_condition_single_vs_single": {"a": round(within_a, 4), "b": round(within_b, 4),
                                              "mean": round(within_mean, 4)},
        "cross_condition_single_vs_single": {"pairs": [round(c, 4) for c in cross],
                                            "mean": round(cross_mean, 4)},
        "noise_matched_gap": round(cross_mean - within_mean, 4),
        "unmatched_ave_vs_ave": round(ave_ave, 4),
        "ceiling_if_condition_were_irrelevant": round(irrelevant_ceiling, 4),
        "disattenuated_true_cross_condition": (round(ave_ave / irrelevant_ceiling, 4)
                                               if irrelevant_ceiling > 0 else None),
        "verdict": ("CONDITION_CARRIES_SIGNAL" if cross_mean < within_mean
                    else "CONDITION_NOT_DISTINGUISHABLE_FROM_NOISE"),
    }


@dataclass
class Tile:
    """One peak tile: a 150 bp window inside a called expression peak.

    **This is the LIKE-FOR-LIKE comparator to the designed grid, and the frag library is not.** Measured on
    the frag file: 96.7% of random genomic fragments sit below 2x the median expression, tightly packed at a
    ~0.87 baseline with a thin tail to 400. So a random-fragment library asks *"is there any promoter in
    this arbitrary piece of genome?"* -- a rare-event DETECTION task -- while the designed grid asked *"how
    strong is this promoter?"*. Peak tiles are selected from called peaks, so they ask the grid's question
    on real genomic sequence, and they are a fixed 150 bp exactly like the grid, which means the grid's
    winning positional feature applies here and the comparison is finally fair.
    """

    seq: str
    expression: float
    peak: str
    strand: str
    category: str
    active: str
    rep1: float = 0.0
    rep2: float = 0.0


def load_peak_tiles(cache_dir: Path | str = DEFAULT_CACHE, *, allow_download: bool = True,
                    ) -> tuple[list[Tile], LoadStats]:
    """Load the peak-tile table. **TAB-delimited, unlike its space-delimited siblings in the same series.**"""
    p = ensure_cached(cache_dir, which="peak_tile", allow_download=allow_download)
    st = LoadStats(path=str(p))
    out: list[Tile] = []
    lengths: list[int] = []
    with gzip.open(p, "rt", errors="replace") as fh:
        first = next(fh)
        delim = sniff_delimiter(first)
        st.delimiter = "tab" if delim == "\t" else "whitespace"
        hdr = _split(first, delim)
        st.n_header_columns = len(hdr)
        if len(hdr) < MIN_COLUMNS:
            raise GenomeWideDataError(
                f"header parsed to {len(hdr)} columns under {st.delimiter}; the delimiter sniff is wrong "
                f"for {p.name}"
            )
        ix = {c: i for i, c in enumerate(hdr)}
        need = ("variant", "expn_med", "peak_start", "peak_end", "strand",
                "RNA_exp_sum_1", "RNA_exp_sum_2")
        missing = [c for c in need if c not in ix]
        if missing:
            raise GenomeWideDataError(f"{p.name} lacks required columns {missing}; has {hdr}")
        for line in fh:
            st.n_lines += 1
            q = _split(line, delim)
            if len(q) != len(hdr):
                st.n_field_count_mismatch += 1
                continue
            seq = q[ix["variant"]].upper()
            if set(seq) - set("ACGT"):
                st.n_non_acgt += 1
                continue
            try:
                e = float(q[ix["expn_med"]])
            except ValueError:
                st.n_unparseable_value += 1
                continue
            try:
                r1 = float(q[ix["RNA_exp_sum_1"]]); r2 = float(q[ix["RNA_exp_sum_2"]])
            except ValueError:
                st.n_unparseable_value += 1
                continue
            out.append(Tile(seq=seq, expression=e, rep1=r1, rep2=r2,
                            # the PEAK is the split unit: tiles within one peak overlap heavily
                            peak=f"{q[ix['peak_start']]}_{q[ix['peak_end']]}_{q[ix['strand']]}",
                            strand=q[ix["strand"]],
                            category=q[ix["category"]] if "category" in ix else "",
                            active=q[ix["active"]] if "active" in ix else ""))
            lengths.append(len(seq))
    st.n_parsed = len(out)
    if not out:
        raise GenomeWideDataError(f"{p.name} parsed ZERO tiles under the {st.delimiter} delimiter")
    lengths.sort()
    st.length_min, st.length_max = lengths[0], lengths[-1]
    st.length_median = lengths[len(lengths) // 2]
    st.n_distinct_lengths = len(set(lengths))
    return out, st


def leave_peak_out(tiles: list[Tile], *, frac: float = 0.25, seed: int = 0):
    """Hold out WHOLE PEAKS -- the direct analogue of the designed grid's leave-element-out.

    Tiles inside one peak are offset by tens of bp and overlap heavily, so splitting tiles at random puts
    near-identical sequences on both sides. Membership is a pure function of the peak id, so disjointness is
    BY CONSTRUCTION.
    """
    peaks = sorted({t.peak for t in tiles})
    if len(peaks) < 2:
        raise ValueError(f"need >= 2 peaks to split on, got {len(peaks)}")
    import random as _r
    ps = list(peaks)
    _r.Random(seed).shuffle(ps)
    k = max(1, int(round(len(ps) * frac)))
    held = set(ps[:k])
    train = [t for t in tiles if t.peak not in held]
    test = [t for t in tiles if t.peak in held]
    return train, test, len(held)


def positional_onehot(seqs: list[str]) -> list[list[float]]:
    """Flat per-position one-hot. **Valid for peak TILES only** -- they are a fixed 150 bp.

    Deliberately NOT applicable to the frag library (48-475 bp, 280 distinct lengths), and it REFUSES a
    ragged input rather than padding, because a silent pad would produce a plausible meaningless number.
    """
    if not seqs:
        return []
    L = len(seqs[0])
    if any(len(s) != L for s in seqs):
        raise GenomeWideDataError(
            "positional_onehot requires a FIXED length; this input is ragged. Peak tiles are 150 bp and "
            "qualify; sheared fragments are 48-475 bp and do not -- use kmer_features for those rather "
            "than padding, which would invent an anchor."
        )
    idx = {b: i for i, b in enumerate("ACGT")}
    out = []
    for s in seqs:
        v = [0.0] * (L * 4)
        for i, b in enumerate(s):
            j = idx.get(b)
            if j is not None:
                v[i * 4 + j] = 1.0
        out.append(v)
    return out


# ---------------------------------------------------------------------------------------------------
# the split -- contiguous genomic blocks, because 91.5% of fragments overlap
# ---------------------------------------------------------------------------------------------------
def position_blocked_split(frags: list[Fragment], *, n_blocks: int = 10, held_blocks: int = 2,
                           seed: int = 0) -> tuple[list[Fragment], list[Fragment], list[int]]:
    """Hold out CONTIGUOUS GENOMIC BLOCKS. The only split offered here, deliberately.

    A sheared library tiles the genome: 91.5% of consecutive fragments overlap, so two fragments can share
    hundreds of bases and a random split puts near-duplicates on both sides. Blocking by position is the
    bacterial analogue of a held-out chromosome. Disjointness is BY CONSTRUCTION -- block membership is a
    pure function of a fragment's midpoint.

    A fragment STRADDLING a held-block boundary is assigned by midpoint, so a held-out fragment can still
    share a little sequence with a training fragment across the seam. With 10 blocks over a 4.64 Mb genome
    the seams are 2 of ~4,600 fragment-widths, which is reported rather than hidden.
    """
    if not frags:
        raise ValueError("no fragments")
    if not 1 <= held_blocks < n_blocks:
        raise ValueError(f"held_blocks must be in [1, {n_blocks - 1}], got {held_blocks}")
    lo = min(f.start for f in frags)
    hi = max(f.end for f in frags)
    width = (hi - lo) / n_blocks

    def block_of(f: Fragment) -> int:
        mid = (f.start + f.end) / 2
        return min(n_blocks - 1, int((mid - lo) / width))

    import random as _r
    bs = list(range(n_blocks))
    _r.Random(seed).shuffle(bs)
    held = set(bs[:held_blocks])
    train = [f for f in frags if block_of(f) not in held]
    test = [f for f in frags if block_of(f) in held]
    return train, test, sorted(held)


# ---------------------------------------------------------------------------------------------------
# features -- length-invariant only, because fragments are 48-475 bp
# ---------------------------------------------------------------------------------------------------
def kmer_features(seqs: list[str], k: int = 4) -> list[list[float]]:
    """Normalised k-mer frequencies. Length-INVARIANT, which is why this is the default here.

    `positional_onehot` (the designed grid's winner) is deliberately NOT offered: these fragments have 280
    distinct lengths, so it would need an anchor, and no anchor is obviously right because the promoter
    sits somewhere inside the fragment. Silently anchoring at position 0 would produce a plausible,
    meaningless number.
    """
    from itertools import product
    idx = {"".join(p): i for i, p in enumerate(product("ACGT", repeat=k))}
    out = []
    for s in seqs:
        v = [0.0] * len(idx)
        n = 0
        for i in range(len(s) - k + 1):
            j = idx.get(s[i:i + k])
            if j is not None:
                v[j] += 1.0
                n += 1
        if n:
            v = [x / n for x in v]
        out.append(v)
    return out


def gc_features(seqs: list[str]) -> list[list[float]]:
    """GC content alone. The mandatory crude baseline."""
    return [[(s.count("G") + s.count("C")) / len(s)] for s in seqs]
