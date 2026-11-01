"""Feasibility + review measurements for the frame-conditioning plan. Read-only, offline, ~30 s.

Emits `wiki/frame_plan_feasibility_2026-11-01.json`. Every number the plan or its review cites comes
from HERE rather than from prose, because this project has a standing rule that a plan's numeric
baseline must be sourced from an artifact -- an unsourced figure reads as authoritative and
propagates.

Five measurements, three of which changed the plan:

  1. STATED-COORDINATE AGREEMENT. The substrate's reference is U00096.2; the cached annotation is
     NC_000913.3. If the stated `start` column agreed, the whole locate step would be unnecessary.
  2. LOCATE RATE by direct sequence search over both strands, at full scale rather than a spot check.
  3. CDS FILE ORDER. The upstream extractor threads a running `last_end` through features in FILE
     order, so whether that order is coordinate-ascending decides whether "replicate it faithfully"
     and "sort first" are the same program. They are not.
  4. CORPUS-WINDOW REPLICA + overlap distribution. Whether `in_window_any` is a usable feature at
     all, and whether it is dominated by trivial 1-bp overlaps.
  5. PER-STRATUM EXPRESSION SPREAD. Whether the two strata a stratified comparison would use are
     exchangeable in outcome distribution. They are not, which is a live confound for that cell.

Nothing here reads the substrate's `start`/`end` columns for anything except measurement (1).
"""

from __future__ import annotations

import collections
import gzip
import json
import pathlib
import statistics
import sys
import time

REPO = pathlib.Path(__file__).resolve().parents[1]
SUB = REPO / "data/mpra/GSE144621_peak_tile_expression_formatted_std.txt.gz"
FNA = REPO / "data/cache/refseq/GCF_000005845.2/genome.fna"
GFF = REPO / "data/cache/refseq/GCF_000005845.2/annotations.gff3"
OUT = REPO / "wiki/frame_plan_feasibility_2026-11-01.json"

RL, MIN_IG, TILE_LEN = 200, 5, 150
COMP = str.maketrans("ACGT", "TGCA")


def rc(s: str) -> str:
    return s.translate(COMP)[::-1]


def main() -> int:
    t0 = time.time()
    for p in (SUB, FNA, GFF):
        if not p.is_file():
            print("REFUSED: %s absent (gitignored substrate not cached)" % p)
            return 3

    # ---- substrate, through the OWNER ---------------------------------------------------------
    # genomewide.load_peak_tiles is the registered owner of GSE144621 (see
    # tests/test_substrate_loader_ownership.py). It is called rather than re-parsed: a second reader
    # does not inherit `sniff_delimiter`, and the peak-tile file is TAB-delimited while its fragment
    # siblings are SPACE-delimited, so a hardcoded split parses one of them to zero rows.
    from dna_decode.glm import genomewide as gw

    all_rows, _ = gw.load_peak_tiles(SUB.parent, allow_download=False)
    tiles = [t for t in all_rows if t.category == "tile" and len(t.seq) == TILE_LEN]

    # the stated-coordinate column is NOT on the Tile dataclass, so measurement (1) needs the raw
    # column. Read it separately and use it for nothing else.
    stated: dict[str, int] = {}
    with gzip.open(SUB, "rt", errors="replace") as fh:
        hdr = next(fh).rstrip("\n").split("\t")
        ix = {c: i for i, c in enumerate(hdr)}
        for line in fh:
            q = line.rstrip("\n").split("\t")
            if len(q) != len(hdr) or q[ix["category"]] != "tile":
                continue
            try:
                stated[q[ix["variant"]].upper()] = int(q[ix["start"]])
            except ValueError:
                pass

    # ---- genome -------------------------------------------------------------------------------
    contig, chunks = None, []
    for line in FNA.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            contig = line[1:].split()[0]
        else:
            chunks.append(line.strip())
    G = "".join(chunks).upper()

    # ---- (2) locate -----------------------------------------------------------------------------
    want: dict[str, list] = {}
    for t in tiles:
        want.setdefault(t.seq, []).append(t)
    hits: dict[str, list[int]] = {}
    for strand, S in (("+", G), ("-", rc(G))):
        n = len(S)
        for p in range(0, n - TILE_LEN + 1):
            k = S[p:p + TILE_LEN]
            if k in want:
                hits.setdefault(k, []).append(p if strand == "+" else n - p - TILE_LEN)
    uniq = {k: v[0] for k, v in hits.items() if len(v) == 1}

    # ---- (1) stated-coordinate agreement, on the uniquely-located set --------------------------
    checked = agree = 0
    offsets = collections.Counter()
    for k, pos in uniq.items():
        st = stated.get(k)
        if st is None:
            continue
        checked += 1
        offsets[pos - st] += 1
        if pos == st or pos == st - 1:
            agree += 1

    # ---- (3) CDS file order ---------------------------------------------------------------------
    from dna_decode.data.annotations import parse_gff3

    tab = parse_gff3(GFF)
    cds = tab[(tab["type"] == "CDS") & (tab["locus_tag"] != "")]
    order = {}
    for strand in ("+", "-"):
        starts = list(cds[cds["strand"] == strand]["start"])
        order[strand] = {
            "n_rows": len(starts),
            "file_order_ascending_by_start": bool(
                all(b >= a for a, b in zip(starts, starts[1:]))
            ),
        }

    # ---- (4) corpus-window replica --------------------------------------------------------------
    # SORTED ascending, which DIVERGES from the upstream extractor's reliance on file order. The
    # divergence is deliberate and recorded: measurement (3) shows file order is NOT ascending, so
    # threading `last_end` through it would produce a window set matching neither the corpus nor a
    # coherent rule. See the plan's Step 1 and the review's critical issue C2.
    fwd = sorted((int(r.start) - 1, int(r.end)) for r in cds[cds["strand"] == "+"].itertuples())
    rev = sorted((int(r.start) - 1, int(r.end)) for r in cds[cds["strand"] == "-"].itertuples())
    windows: list[tuple[int, int]] = []
    last_end = 0
    for st, en in fwd:
        if st < RL:
            continue
        if st - last_end >= MIN_IG:
            windows.append((st - RL, st))
        last_end = en
    last_end_rev = len(G)
    for st, en in reversed(rev):
        if len(G) - en < RL:
            continue
        if last_end_rev - en >= MIN_IG:
            windows.append((en, en + RL))
        last_end_rev = st

    cover = bytearray(len(G))
    for a, b in windows:
        for p in range(max(a, 0), min(b, len(G))):
            cover[p] = 1

    frac_hist = collections.Counter()
    in_any = in_full = trivial = 0
    inw: list[float] = []
    outw: list[float] = []
    for k, pos in uniq.items():
        c = sum(cover[pos:pos + TILE_LEN])
        f = c / TILE_LEN
        frac_hist[round(f, 1)] += 1
        e = want[k][0].expression
        if c > 0:
            in_any += 1
            inw.append(e)
            if f <= 0.2:
                trivial += 1
        else:
            outw.append(e)
        if c == TILE_LEN:
            in_full += 1

    def spread(xs: list[float]) -> dict:
        xs = sorted(xs)
        n = len(xs)
        q1, q3 = xs[n // 4], xs[3 * n // 4]
        return {
            "n": n,
            "median": round(statistics.median(xs), 4),
            "q1": round(q1, 4),
            "q3": round(q3, 4),
            "iqr": round(q3 - q1, 4),
            "sd": round(statistics.pstdev(xs), 3),
            "max": round(max(xs), 1),
        }

    rec = {
        "schema": "frame-plan-feasibility-v1",
        "date": "2026-11-01",
        "contig": contig,
        "genome_length": len(G),
        "regulatory_region_length": RL,
        "min_intergenic_length": MIN_IG,
        "substrate_owner_called": "dna_decode.glm.genomewide.load_peak_tiles",
        "rows_from_owner": len(all_rows),
        "tiles_at_150bp": len(tiles),
        "distinct_tile_sequences": len(want),
        "rows_sharing_a_sequence": len(tiles) - len(want),
        "locate": {
            "located_distinct": sum(1 for k in want if k in hits),
            "uniquely_located": len(uniq),
            "ambiguous_multi_hit": sum(1 for k in want if len(hits.get(k, [])) > 1),
            "not_found": sum(1 for k in want if k not in hits),
            "unique_rate": round(len(uniq) / len(want), 4),
        },
        "stated_coordinate_agreement": {
            "checked": checked,
            "agree": agree,
            "rate": round(agree / checked, 4) if checked else None,
            "top_offsets": offsets.most_common(8),
            "reading": (
                "an agreement rate at this level is a COINCIDENCE rate, and the top offsets form "
                "FAMILIES rather than one constant, which is cumulative indel drift between "
                "assembly versions -- not a shift that could be corrected. Locate by sequence."
            ),
        },
        "cds_file_order": {
            "per_strand": order,
            "why_it_matters": (
                "the upstream extractor threads a running `last_end` through features in FILE "
                "order. If that order is not coordinate-ascending, then 'replicate faithfully' and "
                "'sort first' are DIFFERENT programs producing different window sets. The replica "
                "below sorts, which is a recorded DIVERGENCE, not a reproduction."
            ),
        },
        "corpus_windows": {
            "n_windows": len(windows),
            "bases_covered": sum(cover),
            "genome_fraction": round(sum(cover) / len(G), 4),
            "ordering": "sorted-ascending (diverges from upstream file-order threading)",
        },
        "overlap": {
            "in_window_any_frac_gt_0": in_any,
            "in_window_full": in_full,
            "trivial_overlap_le_0_2": trivial,
            "trivial_share_of_any": round(trivial / in_any, 4) if in_any else None,
            "fraction_histogram_rounded_1dp": {str(k): v for k, v in sorted(frac_hist.items())},
            "reading": (
                "the distribution is bimodal-ish with substantial mass in between, so a BINARY "
                "indicator discards a genuinely graded variable. It is NOT dominated by trivial "
                "overlaps, which refutes that specific worry. Note also that `in_window_any` has no "
                "stated threshold: at `>0` it fires on the count above, and a small epsilon would "
                "lower it, so the threshold must be named."
            ),
        },
        "per_stratum_expression": {
            "in_window": spread(inw),
            "out_of_window": spread(outw),
            "all": spread(inw + outw),
            "reading": (
                "the two strata are NOT exchangeable in outcome distribution: the out-of-window IQR "
                "is materially narrower while its sd is LARGER (outlier-driven), so the shapes "
                "differ, not just the scales. A rank correlation on a compressed outcome attenuates "
                "for BOTH channels, so a smaller within-stratum margin is expected for reasons "
                "unrelated to frame absorption. A size-matched RANDOM stratification is the needed "
                "null, and each stratum's spread must be reported beside its margin."
            ),
        },
        "wall_seconds": round(time.time() - t0, 1),
    }
    OUT.write_text(json.dumps(rec, indent=2, sort_keys=True), encoding="utf-8")
    print("locate: %d unique of %d distinct (%.4f); ambiguous %d; not found %d"
          % (len(uniq), len(want), rec["locate"]["unique_rate"],
             rec["locate"]["ambiguous_multi_hit"], rec["locate"]["not_found"]))
    print("stated-coordinate agreement: %d/%d" % (agree, checked))
    print("cds file order ascending: %s" % {k: v["file_order_ascending_by_start"]
                                            for k, v in order.items()})
    print("windows: %d covering %d bp (%.4f of genome)"
          % (len(windows), sum(cover), rec["corpus_windows"]["genome_fraction"]))
    print("in_window_any %d / full %d / trivial %d" % (in_any, in_full, trivial))
    print("IQR in-window %.4f vs out-of-window %.4f"
          % (rec["per_stratum_expression"]["in_window"]["iqr"],
             rec["per_stratum_expression"]["out_of_window"]["iqr"]))
    print("-> %s" % OUT)
    return 0


if __name__ == "__main__":
    sys.exit(main())
