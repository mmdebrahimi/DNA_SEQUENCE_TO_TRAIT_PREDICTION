"""The NATURAL reference corpus — real regulatory sequence, extracted from a real genome.

Every falsifier claim about a generator needs something real to be compared against, and this supplies it
without a download: the committed MG1655 reference plus its GFF3 give thousands of genuine *E. coli*
promoter regions, which are the windows immediately upstream of annotated CDS starts.

**Strand is the whole difficulty and it is where this kind of code silently goes wrong.** "Upstream" means
*5-prime of the start codon*, which is to the LEFT in genome coordinates for a plus-strand gene and to the
RIGHT for a minus-strand gene — and the minus-strand window must then be reverse-complemented, or the
extracted "promoter" is the reverse complement of a real one. A promoter set that is 50% reverse-complement
noise still looks like DNA, still has plausible GC content, and would quietly destroy every downstream
comparison, so `test_glm_corpus.py` pins both strands against hand-computed slices of a synthetic genome
and against the real reference.

**GFF3 is 1-BASED INCLUSIVE; Python slices are 0-based half-open.** That off-by-one is the second trap, so
the conversion happens in exactly one place here (`_upstream_slice`) rather than at each call site.

Honest scope: an upstream window is a *proxy* for a promoter, not a mapped transcription start site. Some
windows overlap the preceding gene (operons; tight genomes), which is biology, not a bug — so overlap is
MEASURED and reported (`n_overlapping_cds`) and optionally excluded, never silently ignored.
"""
from __future__ import annotations

from dataclasses import dataclass, field

_COMPLEMENT = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")

_DNA = frozenset("ACGT")


def revcomp(seq: str) -> str:
    """Reverse complement. Pure."""
    return seq.translate(_COMPLEMENT)[::-1]


@dataclass(frozen=True)
class UpstreamWindow:
    """One extracted upstream (candidate promoter) region."""

    gene_id: str
    seqid: str
    strand: str
    #: 0-based half-open span on the FORWARD strand, before any reverse-complement.
    start: int
    end: int
    seq: str
    overlaps_cds: bool = False

    @property
    def length(self) -> int:
        return len(self.seq)


@dataclass
class CorpusStats:
    """What the extraction actually produced — reported, so a degraded corpus is visible."""

    n_cds_total: int = 0
    n_windows: int = 0
    n_skipped_edge: int = 0
    n_skipped_ambiguous: int = 0
    n_overlapping_cds: int = 0
    strand_counts: dict = field(default_factory=dict)

    def as_dict(self) -> dict:
        return {
            "n_cds_total": self.n_cds_total,
            "n_windows": self.n_windows,
            "n_skipped_edge": self.n_skipped_edge,
            "n_skipped_ambiguous": self.n_skipped_ambiguous,
            "n_overlapping_cds": self.n_overlapping_cds,
            "strand_counts": dict(self.strand_counts),
        }


def _upstream_slice(start_1b: int, end_1b: int, strand: str, length: int) -> tuple[int, int]:
    """PURE: GFF3 1-based-inclusive CDS bounds -> the 0-based half-open upstream window.

    Plus strand: the window ends where the CDS begins, so it is `[start-1-length, start-1)`.
    Minus strand: the gene reads right-to-left, so 5-prime is PAST `end` — `[end, end+length)`.
    Returns forward-strand coordinates; the caller reverse-complements for the minus strand.
    """
    if strand == "+":
        s = start_1b - 1 - length
        return s, start_1b - 1
    if strand == "-":
        return end_1b, end_1b + length
    raise ValueError(f"strand must be '+' or '-', got {strand!r}")


def as_rows(table) -> list[dict]:
    """Normalise `parse_gff3` output to a list of dicts.

    `AnnotationTable` is an aliased **pandas DataFrame**, and `list(df)` yields COLUMN NAMES rather than
    rows — a silent shape error that surfaces only as `'str' object has no attribute 'get'` further down.
    Accepting both shapes here means no call site has to remember which it was handed.
    """
    if hasattr(table, "to_dict") and hasattr(table, "columns"):
        return table.to_dict("records")
    return list(table)


def extract_upstream_windows(
    genome: dict[str, str],
    rows,
    *,
    length: int = 150,
    feature_type: str = "CDS",
    exclude_overlapping: bool = False,
    require_unambiguous: bool = True,
) -> tuple[list[UpstreamWindow], CorpusStats]:
    """Pull the `length`-bp region 5-prime of every `feature_type` start.

    `rows` is `dna_decode.data.annotations.parse_gff3` output (1-based inclusive, with `strand`) — either
    the DataFrame itself or a list of dicts; `as_rows` normalises both.
    A window running off a contig end is SKIPPED and counted, never truncated — a short window would
    silently change the length distribution the discriminator keys on.
    """
    if length <= 0:
        raise ValueError("length must be positive")
    stats = CorpusStats()
    rows = as_rows(rows)
    cds = [r for r in rows if r.get("type") == feature_type and r.get("strand") in ("+", "-")]
    stats.n_cds_total = len(cds)

    # Interval index per contig, for the overlap measurement. Built once, not per window.
    spans: dict[str, list[tuple[int, int]]] = {}
    for r in cds:
        spans.setdefault(r["seqid"], []).append((r["start"] - 1, r["end"]))
    for v in spans.values():
        v.sort()

    out: list[UpstreamWindow] = []
    for r in cds:
        seqid = r["seqid"]
        contig = genome.get(seqid)
        if contig is None:
            continue
        s, e = _upstream_slice(r["start"], r["end"], r["strand"], length)
        if s < 0 or e > len(contig):
            stats.n_skipped_edge += 1
            continue
        fwd = contig[s:e].upper()
        if require_unambiguous and set(fwd) - _DNA:
            stats.n_skipped_ambiguous += 1
            continue
        seq = fwd if r["strand"] == "+" else revcomp(fwd)
        # Overlap with ANY other CDS. Common and real in a compact genome (operons), so it is counted.
        ov = False
        for (cs, ce) in spans.get(seqid, ()):
            if cs < e and s < ce and not (cs == r["start"] - 1 and ce == r["end"]):
                ov = True
                break
        if ov:
            stats.n_overlapping_cds += 1
            if exclude_overlapping:
                continue
        stats.strand_counts[r["strand"]] = stats.strand_counts.get(r["strand"], 0) + 1
        out.append(UpstreamWindow(
            gene_id=r.get("gene_id") or r.get("locus_tag") or f"{seqid}:{s}",
            seqid=seqid, strand=r["strand"], start=s, end=e, seq=seq, overlaps_cds=ov,
        ))
    stats.n_windows = len(out)
    return out, stats


def window_context(genome: dict[str, str], w: UpstreamWindow, ctx: int) -> str | None:
    """The `ctx` bases immediately PRECEDING `w` in the window's own orientation — a generation prompt.

    This is what makes a CONDITIONAL comparison possible, and that turned out to be necessary rather than
    optional. Measured: against unaligned random genomic windows an order-3 Markov null scores 0.5086
    against a 0.5255 ceiling — i.e. BELOW its own noise floor, so that reference has no discriminating
    power at all. Only the ALIGNED promoter set separates (0.6560 vs 0.5292), because its positional
    structure is anchored to the CDS start. An unconditional sample cannot be scored against an aligned
    set fairly, so the meaningful task is *given the genomic context, continue into the promoter*.

    Returns None when the context would run off the contig — skipped, never truncated, so prompt length
    stays constant across the set.
    """
    contig = genome.get(w.seqid)
    if contig is None or ctx <= 0:
        return None
    if w.strand == "+":
        s = w.start - ctx
        if s < 0:
            return None
        out = contig[s:w.start]
    else:
        # the window was reverse-complemented, so its 5-prime neighbour is to the RIGHT and also flips
        e = w.end + ctx
        if e > len(contig):
            return None
        out = revcomp(contig[w.end:e])
    out = out.upper()
    return None if set(out) - _DNA else out


def random_genomic_windows(
    genome: dict[str, str],
    *,
    length: int = 150,
    n: int = 2000,
    seed: int = 0,
    require_unambiguous: bool = True,
) -> tuple[list[UpstreamWindow], CorpusStats]:
    """Uniformly-sampled windows from anywhere in the genome — the reference for an UNCONDITIONAL generator.

    **Why this exists, and it is a correctness issue, not a nicety.** An unprompted sample from a genome
    language model is drawn from the model's whole-genome prior, and a prokaryotic genome is mostly CODING
    (the first real GENERator sample began `ATGCTGCAA…` — a start codon). Scoring such a sample against a
    PROMOTER set would therefore largely measure "coding vs non-coding", not whether the generator produces
    realistic DNA. The matching reference for unconditional generation is the genome's own distribution.

    The promoter set (`extract_upstream_windows`) stays the right reference for a *conditional* promoter
    generator; these two are not interchangeable and the artifact records which was used.
    """
    import random as _random
    if length <= 0:
        raise ValueError("length must be positive")
    rng = _random.Random(seed)
    stats = CorpusStats()
    names = [c for c, s in genome.items() if len(s) > length]
    if not names:
        return [], stats
    weights = [len(genome[c]) for c in names]
    out: list[UpstreamWindow] = []
    attempts = 0
    while len(out) < n and attempts < n * 20:
        attempts += 1
        c = rng.choices(names, weights=weights, k=1)[0]
        contig = genome[c]
        s = rng.randrange(0, len(contig) - length)
        seq = contig[s:s + length].upper()
        if require_unambiguous and set(seq) - _DNA:
            stats.n_skipped_ambiguous += 1
            continue
        strand = rng.choice("+-")
        if strand == "-":
            seq = revcomp(seq)
        stats.strand_counts[strand] = stats.strand_counts.get(strand, 0) + 1
        out.append(UpstreamWindow(gene_id=f"{c}:{s}", seqid=c, strand=strand,
                                  start=s, end=s + length, seq=seq))
    stats.n_windows = len(out)
    stats.n_cds_total = 0
    return out, stats


def load_fasta(path) -> dict[str, str]:
    """Minimal FASTA reader -> {first_token_of_header: sequence}. No biopython needed."""
    from pathlib import Path
    seqs: dict[str, str] = {}
    name: str | None = None
    buf: list[str] = []
    for line in Path(path).read_text().splitlines():
        if line.startswith(">"):
            if name is not None:
                seqs[name] = "".join(buf)
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip())
    if name is not None:
        seqs[name] = "".join(buf)
    return seqs


def gc_fraction(seq: str) -> float:
    """GC content of a sequence. Returns 0.0 for an empty string rather than raising."""
    if not seq:
        return 0.0
    return (seq.count("G") + seq.count("C") + seq.count("g") + seq.count("c")) / len(seq)
