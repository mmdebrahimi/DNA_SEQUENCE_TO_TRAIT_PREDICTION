"""Measure the window-ordering divergence instead of asserting it. Read-only, offline.

**Why this is its own script rather than an extension of `frame_plan_feasibility_probe.py`, and the
reason is mechanical.** That script fixes its output path as a module constant, includes a
`wall_seconds` field in its record, and writes the file unconditionally. So any run of a modified
version rewrites `wiki/frame_plan_feasibility_2026-11-01.json`, and the timing field guarantees the
rewrite differs even when every input is identical. That artifact is committed and cited by both the
frame-conditioning plan and its review, so extending the script and promising to leave the artifact
untouched was a contradiction. A separate script makes the promise true by construction, where a flag
on the existing script would depend on every future caller remembering to pass it.

**The accepted cost is duplication.** Both scripts carry a corpus-window replica and two copies can
drift into producing different window sets under the same name. If that drift becomes real, the fix
on record is to extract the replica into the package rather than to re-merge the scripts.

**The question.** The upstream extractor threads a running `last_end` through gene features in FILE
order, gating each window on `start - last_end >= min_intergenic_length`. Whether that order is
coordinate-ascending therefore decides which windows get emitted at all. Measured: the cached GFF3's
coding-sequence rows are NOT ascending by start on either strand. So "replicate the upstream function
faithfully" and "sort first" are different programs, and this script measures how different.

Emits `wiki/frame_plan_window_ordering_2026-11-01.json`. A symmetric difference of ZERO is a
legitimate and informative outcome -- it would mean the divergence is not real -- and is reported as
a finding rather than treated as a failed run.
"""

from __future__ import annotations

import json
import pathlib
import sys
import time

REPO = pathlib.Path(__file__).resolve().parents[1]
FNA = REPO / "data/cache/refseq/GCF_000005845.2/genome.fna"
GFF = REPO / "data/cache/refseq/GCF_000005845.2/annotations.gff3"
SUB = REPO / "data/mpra"
OUT = REPO / "wiki/frame_plan_window_ordering_2026-11-01.json"

RL, MIN_IG, TILE_LEN = 200, 5, 150
COMP = str.maketrans("ACGT", "TGCA")


def rc(s: str) -> str:
    return s.translate(COMP)[::-1]


def corpus_windows(fwd, rev, genome_len, *, regulatory_region_length=RL,
                   min_intergenic_length=MIN_IG):
    """Replicate `extract_promoter_regions`, threading `last_end` in the ORDER GIVEN.

    `fwd` and `rev` are sequences of (start0, end0_exclusive). The caller decides their order, which
    is the whole point: the upstream function iterates features in file order, so the order is an
    input to the result rather than a presentation detail.

    Two upstream quirks are reproduced deliberately, not corrected. A gene skipped by the
    `start < regulatory_region_length` guard does NOT update the running `last_end`, because the
    `continue` precedes the assignment. And the reverse-strand branch compares against the previous
    gene's START. Both look like defects; the question this answers is what the corpus CONTAINED,
    not what it should have contained.
    """
    windows = []
    last_end = 0
    for st, en in fwd:
        if st < regulatory_region_length:
            continue
        if st - last_end >= min_intergenic_length:
            windows.append((st - regulatory_region_length, st, "+"))
        last_end = en
    last_end_rev = genome_len
    for st, en in reversed(rev):
        if genome_len - en < regulatory_region_length:
            continue
        if last_end_rev - en >= min_intergenic_length:
            windows.append((en, en + regulatory_region_length, "-"))
        last_end_rev = st
    return windows


def coverage(windows, genome_len) -> bytearray:
    cov = bytearray(genome_len)
    for a, b, _s in windows:
        for p in range(max(a, 0), min(b, genome_len)):
            cov[p] = 1
    return cov


def main() -> int:
    t0 = time.time()
    for p in (FNA, GFF):
        if not p.is_file():
            print("REFUSED: %s absent (gitignored substrate not cached)" % p)
            return 3

    from dna_decode.data.annotations import parse_gff3

    tab = parse_gff3(GFF)
    cds = tab[(tab["type"] == "CDS") & (tab["locus_tag"] != "")]

    # FILE order: exactly the order the parser returned, which is the order the upstream function
    # would have iterated had it read this file.
    fwd_file = [(int(r.start) - 1, int(r.end)) for r in cds[cds["strand"] == "+"].itertuples()]
    rev_file = [(int(r.start) - 1, int(r.end)) for r in cds[cds["strand"] == "-"].itertuples()]
    fwd_sorted, rev_sorted = sorted(fwd_file), sorted(rev_file)

    order_is_ascending = {
        "+": all(b >= a for a, b in zip([s for s, _ in fwd_file], [s for s, _ in fwd_file][1:])),
        "-": all(b >= a for a, b in zip([s for s, _ in rev_file], [s for s, _ in rev_file][1:])),
    }

    chunks = []
    contig = None
    for line in FNA.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            contig = line[1:].split()[0]
        else:
            chunks.append(line.strip())
    G = "".join(chunks).upper()

    w_file = corpus_windows(fwd_file, rev_file, len(G))
    w_sorted = corpus_windows(fwd_sorted, rev_sorted, len(G))
    s_file, s_sorted = set(w_file), set(w_sorted)
    sym = s_file ^ s_sorted

    cov_file = coverage(w_file, len(G))
    cov_sorted = coverage(w_sorted, len(G))

    # ---- how many uniquely-located tiles change their in-window status? -----------------------
    from dna_decode.glm import genomewide as gw

    rows, _ = gw.load_peak_tiles(SUB, allow_download=False)
    tiles = [t for t in rows if t.category == "tile" and len(t.seq) == TILE_LEN]
    want = {}
    for t in tiles:
        want.setdefault(t.seq, 0)
        want[t.seq] += 1
    hits = {}
    for strand, S in (("+", G), ("-", rc(G))):
        n = len(S)
        for p in range(0, n - TILE_LEN + 1):
            k = S[p:p + TILE_LEN]
            if k in want:
                hits.setdefault(k, []).append(p if strand == "+" else n - p - TILE_LEN)
    uniq = {k: v[0] for k, v in hits.items() if len(v) == 1}

    changed = in_file = in_sorted = 0
    for pos in uniq.values():
        a = sum(cov_file[pos:pos + TILE_LEN]) > 0
        b = sum(cov_sorted[pos:pos + TILE_LEN]) > 0
        in_file += a
        in_sorted += b
        if a != b:
            changed += 1

    divergent = len(sym) > 0 or changed > 0
    rec = {
        "schema": "frame-plan-window-ordering-v1",
        "date": "2026-11-01",
        "contig": contig,
        "genome_length": len(G),
        "regulatory_region_length": RL,
        "min_intergenic_length": MIN_IG,
        "cds_rows_with_locus_tag": int(len(cds)),
        "file_order_is_ascending_by_start": order_is_ascending,
        "file_order": {
            "n_windows": len(w_file),
            "n_distinct_windows": len(s_file),
            "bases_covered": sum(cov_file),
            "tiles_in_window": in_file,
        },
        "sorted_order": {
            "n_windows": len(w_sorted),
            "n_distinct_windows": len(s_sorted),
            "bases_covered": sum(cov_sorted),
            "tiles_in_window": in_sorted,
        },
        "symmetric_difference_windows": len(sym),
        "tiles_changing_in_window_status": changed,
        "uniquely_located_tiles": len(uniq),
        "verdict": "ORDERING_DIVERGES" if divergent else "ORDERING_IMMATERIAL",
        "reading": (
            "A non-zero symmetric difference means 'replicate the upstream function faithfully' and "
            "'sort first' are different programs, so a plan cannot claim both. A ZERO difference "
            "would mean the divergence is not real and the ordering parameter is unnecessary "
            "configuration -- a finding, not a failed run."
        ),
        "honest_limits": [
            "The upstream code reads GenBank, not this GFF3, so FILE order here is this file's "
            "order and not the order that built the published corpus. This bounds the SHAPE of the "
            "divergence, not its exact magnitude upstream.",
            "Both quirks of the upstream function are reproduced rather than corrected, which is "
            "deliberate; a corrected replica would measure a rule nobody trained on.",
            "One organism, one annotation build.",
        ],
        "wall_seconds": round(time.time() - t0, 1),
    }
    OUT.write_text(json.dumps(rec, indent=2, sort_keys=True), encoding="utf-8")
    print("file order ascending by start: %s" % order_is_ascending)
    print("windows  file=%d  sorted=%d" % (len(w_file), len(w_sorted)))
    print("bases    file=%d  sorted=%d" % (sum(cov_file), sum(cov_sorted)))
    print("symmetric difference: %d windows" % len(sym))
    print("tiles changing in-window status: %d of %d uniquely located" % (changed, len(uniq)))
    print("VERDICT: %s" % rec["verdict"])
    print("-> %s" % OUT)
    return 0


if __name__ == "__main__":
    sys.exit(main())
