"""W0 probe: pin the S3 coordinate schema, re-verify the oracle footprint, freeze the tier set.

Step 1 of the locus-level oracle-overlap plan. READ-ONLY on every input. It exists because the
exact-sequence contamination result (100 of 182,859 records byte-identical to a candidate oracle's
training splits) states twice in its own `honest_limits` that it bounds VERBATIM reuse and that the
locus-level question needs a coordinate intersect -- and that intersect cannot be built until the
S3 coordinate schema is pinned by measurement rather than assumed.

WHAT THIS REFUSES ON (non-zero exit), because a silent zero reads as independence:
  * a cell sheet whose header row or coordinate columns cannot be located
  * start/end semantics that do not resolve to a single convention
  * any footprint figure disagreeing with the committed sidecar
  * a footprint BED that is not actually merged (unsorted / overlapping / abutting runs)
  * an empty tier set or outcome vocabulary -- a guard looping over an empty list passes
    vacuously forever, which is the recorded vacuous-statistic shape

THE ASSEMBLY RAIL. No coordinate is ever derived from, or cross-checked against, an element ID.
A peer session on the oracle side derived windows as `pos +/- 100` from element IDs, verified the
rule on ONE row where it held, and a 400-row batch check then held for 1 and failed for 399 with a
constant ~64,620 bp offset -- because the IDs carry hg19 while the GRCh38 files are lifted over,
and the single confirming row was chr1's start where the offset is zero. Our side carries the
mirror exposure (HepG2's sheet is hg19-native), so this probe reads coordinates from NAMED
ASSEMBLY COLUMNS only and reports the hg38 availability per cell instead of filling it from hg19.
"""

from __future__ import annotations

import argparse
import gzip
import json
import sys
from collections import defaultdict
from pathlib import Path

SCHEMA = "oracle-locus-w0-probe-v1"

#: Frozen AFTER denominators and nulls are known, BEFORE any observed overlap count exists.
#: `identical_interval` is the >=200 bp rung for a 200 nt element; `same_interval_different_sequence`
#: is a separate CLASS, not a tier -- the exact-sequence hash join was blind to it by construction.
TIERS: tuple[tuple[str, int], ...] = (
    (">=1bp", 1),
    (">=50bp", 50),
    (">=100bp", 100),
    ("identical_interval", 200),
)

OUTCOMES: tuple[str, ...] = (
    "LOCUS_OVERLAP_EXCEEDS_NULL",
    "LOCUS_OVERLAP_AT_NULL",
    "LOCUS_OVERLAP_BELOW_NULL",
    "REFUSED_POSITIVE_CONTROL_FAILED",
)

#: The element the oracle's footprint was built from, and ours. Both 200 nt of biology.
ELEMENT_BP = 200
#: The peer's own denominator, carried so the comparison is like-for-like rather than re-based.
GENOME_BP = 3_100_000_000

#: Published split rule, quoted from the oracle's own notebook. Gives a plumbing check that is
#: TRUE BY CONSTRUCTION -- it detects a broken parse and NOTHING ELSE, and is blind to both
#: assembly and off-by-one error. Never report it as evidence.
SPLIT_CHROMS = {
    "val": {"chr19", "chr21", "chrX"},
    "test": {"chr13", "chr7"},
}


class ProbeRefusal(RuntimeError):
    """Raised instead of returning a number that would read as a clean result."""


# --------------------------------------------------------------------------------------------
# chromosome naming -- ONE normaliser, both directions, because the three cells disagree
# --------------------------------------------------------------------------------------------

def to_chr(token: str) -> str:
    """`1` -> `chr1`, `chr1` -> `chr1`. WTC11 writes every chromosome BARE."""
    t = str(token).strip()
    return t if t.startswith("chr") else "chr" + t


def strip_chr(token: str) -> str:
    """Inverse of `to_chr`, so a join can be made in either convention."""
    t = str(token).strip()
    return t[3:] if t.startswith("chr") else t


# --------------------------------------------------------------------------------------------
# S3 schema -- measured per sheet, never inherited
# --------------------------------------------------------------------------------------------

def find_header(rows: list[tuple], want: str = "name") -> int:
    """Index of the header row. The header MOVES between sheets (measured 4 / 2 / 1)."""
    for i, r in enumerate(rows):
        if r and str(r[0] or "").strip().lower() == want:
            return i
    raise ProbeRefusal("no header row whose first cell is %r in the first %d rows"
                       % (want, len(rows)))


def classify_columns(header: tuple) -> dict:
    """Which coordinate columns exist, and which assembly each belongs to."""
    cols = {str(c).strip(): i for i, c in enumerate(header) if c is not None}
    out: dict = {"columns": list(cols), "assemblies": {}}
    for asm in ("hg19", "hg38"):
        trio = {k: cols.get("%s.%s" % (k, asm)) for k in ("chr", "start", "stop")}
        if all(v is not None for v in trio.values()):
            out["assemblies"][asm] = trio
    if "hg38" not in out["assemblies"]:
        raise ProbeRefusal("no complete hg38 chr/start/stop column trio; columns were %r"
                           % (list(cols),))
    return out


def probe_sheet(ws, sheet: str, sample: int) -> dict:
    """Pin one cell sheet's schema and measure its hg38 availability."""
    head = [r for r in ws.iter_rows(min_row=1, max_row=8, values_only=True)]
    hi = find_header(head)
    schema = classify_columns(head[hi])

    idx = schema["assemblies"]["hg38"]
    bare = 0
    n = 0
    na_hg38 = 0
    widths: dict[int, int] = defaultdict(int)
    suffix_f = 0
    for r in ws.iter_rows(min_row=hi + 2, max_row=hi + 1 + sample, values_only=True):
        if r is None or r[0] is None:
            continue
        n += 1
        if str(r[0]).strip().endswith("_F"):
            suffix_f += 1
        c, s, e = (r[idx["chr"]], r[idx["start"]], r[idx["stop"]])
        if c is None or s is None or e is None:
            na_hg38 += 1
            continue
        if not str(c).strip().startswith("chr"):
            bare += 1
        widths[int(e) - int(s)] += 1

    if not widths:
        raise ProbeRefusal("sheet %r yielded no parseable hg38 interval in %d sampled rows"
                           % (sheet, sample))

    # DERIVE the convention rather than assume it: the elements are 200 nt of biology, so a
    # stop-start of exactly 200 is half-open [start, stop) and 201 would be 1-based inclusive.
    modal = max(widths.items(), key=lambda kv: kv[1])[0]
    if modal == ELEMENT_BP:
        convention = "half_open_0_based"
    elif modal == ELEMENT_BP + 1:
        convention = "closed_1_based"
    else:
        raise ProbeRefusal(
            "sheet %r modal width is %d, neither %d (half-open) nor %d (closed) -- the "
            "start/end convention does not resolve and must not be guessed"
            % (sheet, modal, ELEMENT_BP, ELEMENT_BP + 1))

    return {
        "sheet": sheet,
        "header_row_1based": hi + 1,
        "columns": schema["columns"],
        "assemblies_present": sorted(schema["assemblies"]),
        "hg38_column_index": idx,
        "is_hg19_native": "hg19" in schema["assemblies"],
        "n_sampled": n,
        "n_hg38_absent_in_sample": na_hg38,
        "chromosome_convention": ("bare" if bare > n / 2 else "chr_prefixed"),
        "n_bare_chromosome_in_sample": bare,
        "n_names_ending_underscore_F": suffix_f,
        "interval_convention": convention,
        "interval_convention_basis": (
            "modal stop-start is %d against a %d nt element" % (modal, ELEMENT_BP)),
        "width_census_sample": dict(sorted(widths.items(), key=lambda kv: -kv[1])[:6]),
    }


# --------------------------------------------------------------------------------------------
# footprint re-verification + the per-tier null
# --------------------------------------------------------------------------------------------

def read_footprint(path: Path) -> tuple[dict, dict]:
    """Per-chromosome runs plus parse accounting. Skips `track` / `browser` / `#` lines."""
    by_chrom: dict[str, list[tuple[int, int]]] = defaultdict(list)
    acc = {"track_or_comment": 0, "malformed": 0, "blank": 0}
    with gzip.open(path, "rt") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                acc["blank"] += 1
                continue
            if line.startswith(("track", "browser", "#")):
                acc["track_or_comment"] += 1
                continue
            f = line.split("\t")
            if len(f) < 3:
                acc["malformed"] += 1
                continue
            by_chrom[f[0]].append((int(f[1]), int(f[2])))
    return by_chrom, acc


def merge_integrity(by_chrom: dict) -> dict:
    unsorted_ = overlapping = abutting = 0
    for v in by_chrom.values():
        for (s1, e1), (s2, _e2) in zip(v, v[1:]):
            if s2 < s1:
                unsorted_ += 1
            elif s2 < e1:
                overlapping += 1
            elif s2 == e1:
                abutting += 1
    return {"unsorted": unsorted_, "overlapping": overlapping, "abutting": abutting}


def start_window_union(by_chrom: dict, k: int, element_bp: int = ELEMENT_BP) -> int:
    """EXACT count of element starts overlapping any run by >= k bp.

    An element of length L starting at s overlaps run [a, b) by >= k bp iff
    s in [a - L + k, b - k]. The measure of the UNION of those per-run start-windows IS the
    number of qualifying starts.

    Summing the windows instead would DOUBLE-COUNT wherever two merged runs sit closer than the
    element length, and variant-centred tiling around clustered variants is exactly where that
    happens -- measured 1.04x inflation at k=1 on the real train footprint.
    """
    union = 0
    for v in by_chrom.values():
        wins = sorted((s - element_bp + k, e - k + 1) for s, e in v)
        cur_s = cur_e = None
        for lo, hi in wins:
            if hi <= lo:
                continue
            if cur_s is None:
                cur_s, cur_e = lo, hi
            elif lo <= cur_e:
                cur_e = max(cur_e, hi)
            else:
                union += cur_e - cur_s
                cur_s, cur_e = lo, hi
        if cur_s is not None:
            union += cur_e - cur_s
    return union


def naive_window_sum(by_chrom: dict, k: int, element_bp: int = ELEMENT_BP) -> int:
    """The WRONG computation, kept so the artifact can state the inflation it avoids."""
    tot = 0
    for v in by_chrom.values():
        for s, e in v:
            lo, hi = s - element_bp + k, e - k + 1
            if hi > lo:
                tot += hi - lo
    return tot


def verify_split(split: str, path: Path, sidecar: dict | None) -> dict:
    by_chrom, acc = read_footprint(path)
    n_runs = sum(len(v) for v in by_chrom.values())
    merged_bp = sum(e - s for v in by_chrom.values() for s, e in v)
    integ = merge_integrity(by_chrom)
    if any(integ.values()):
        raise ProbeRefusal(
            "%s footprint is NOT a genuine merge (%r); every null here assumes merged runs"
            % (split, integ))

    rec = {
        "split": split,
        "n_runs": n_runs,
        "merged_bp": merged_bp,
        "base_fraction_pct": round(100.0 * merged_bp / GENOME_BP, 4),
        "n_chromosomes": len(by_chrom),
        "chromosomes": sorted(by_chrom),
        "parse_accounting": acc,
        "merge_integrity": integ,
    }

    if sidecar is not None:
        # ONLY the quantities a MERGED bed can verify. The sidecar's `intervals` is the
        # PRE-MERGE count and is NOT the merged run count -- equating them was this probe's own
        # first defect, and it raised a false disagreement (sidecar 1,542,544 vs measured
        # 247,767) that read exactly like a corrupted footprint. Confirmed pre-merge by the
        # sidecar's own identity: 1,542,544 x 200 bp / 54,388,968 bp = 5.672268, matching its
        # `self_overlap_compression` 5.672267949632727 to six decimal places.
        for key, got in (("merged_covered_bases", merged_bp), ("chromosomes", len(by_chrom))):
            want = sidecar.get(key)
            if want is not None and int(want) != got:
                raise ProbeRefusal(
                    "%s footprint disagrees with the committed sidecar on %s: sidecar %s, "
                    "measured %d" % (split, key, want, got))
        rec["sidecar_agrees_on_verifiable_fields"] = ["merged_covered_bases", "chromosomes"]

        # Cross-check the pre-merge count through the compression identity rather than
        # discarding it: a merged BED cannot reproduce it directly, so this is the only
        # available check and it is reported as derived, not measured.
        iv, comp = sidecar.get("intervals"), sidecar.get("self_overlap_compression")
        if iv and comp:
            implied = iv * ELEMENT_BP / merged_bp
            rec["sidecar_premerge_cross_check"] = {
                "intervals_premerge": iv,
                "implied_compression": round(implied, 6),
                "sidecar_compression": comp,
                "agrees": abs(implied - comp) < 1e-6,
                "note": ("`intervals` is PRE-MERGE and not verifiable from a merged BED. This "
                         "identity is the available cross-check; it assumes every pre-merge "
                         "interval is exactly %d bp." % ELEMENT_BP),
            }
            if not rec["sidecar_premerge_cross_check"]["agrees"]:
                raise ProbeRefusal(
                    "%s sidecar is internally inconsistent: intervals x %d / merged_bp = %.6f "
                    "but self_overlap_compression says %s"
                    % (split, ELEMENT_BP, implied, comp))

    # PLUMBING GATE, true by construction from the published split rule. Detects a broken parse
    # and nothing else; blind to assembly and off-by-one. NEVER evidence.
    if split in SPLIT_CHROMS:
        rec["plumbing_gate_chromosomes_match_split_rule"] = (
            set(by_chrom) == SPLIT_CHROMS[split])
        rec["plumbing_gate_is_true_by_construction_never_evidence"] = True

    if split == "train":
        nulls = {}
        for name, k in TIERS:
            exact = start_window_union(by_chrom, k)
            naive = naive_window_sum(by_chrom, k)
            nulls[name] = {
                "k_bp": k,
                "exact_qualifying_starts": exact,
                "exact_pct_of_genome": round(100.0 * exact / GENOME_BP, 4),
                "naive_sum_pct": round(100.0 * naive / GENOME_BP, 4),
                "naive_inflation_x": round(naive / exact, 4) if exact else None,
            }
        rec["random_element_null_per_tier"] = nulls
        gaps = sorted(s2 - e1 for v in by_chrom.values()
                      for (_s1, e1), (s2, _e2) in zip(v, v[1:]))
        if gaps:
            rec["inter_run_gaps"] = {
                "n": len(gaps),
                "median_bp": gaps[len(gaps) // 2],
                "p10_bp": gaps[len(gaps) // 10],
                "p25_bp": gaps[len(gaps) // 4],
                "frac_under_element_bp": round(
                    sum(1 for g in gaps if g < ELEMENT_BP) / len(gaps), 4),
            }
    return rec


# --------------------------------------------------------------------------------------------

def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--s3", default="data/mpra_human/Table S3 - large-scale lib design.xlsx")
    ap.add_argument("--footprint-dir", default=(
        "C:/Users/Laptop Portable/pythonprojects/dna-generator-research/data/"
        "oracle_contamination_set"))
    ap.add_argument("--sample", type=int, default=4000,
                    help="rows per sheet for the schema probe; a header read plus a bounded "
                         "sample, never a full 56 MB parse")
    ap.add_argument("--out", default="wiki/oracle_locus_w0_probe_2026-11-23.json")
    ap.add_argument("--self-check", action="store_true",
                    help="run the pure-logic checks and exit; no file or network access")
    a = ap.parse_args(argv)

    if a.self_check:
        return _self_check()

    if not TIERS or not OUTCOMES:
        raise ProbeRefusal("tier set or outcome vocabulary is EMPTY; a guard looping over an "
                           "empty list passes vacuously forever")

    import openpyxl
    wb = openpyxl.load_workbook(a.s3, read_only=True)
    cells = {}
    for sn in wb.sheetnames:
        key = sn.split()[0]
        cells[key] = probe_sheet(wb[sn], sn, a.sample)
    wb.close()

    fdir = Path(a.footprint_dir)
    side = {}
    sp = fdir / "oracle_locus_footprint.json"
    if sp.exists():
        side = json.loads(sp.read_text(encoding="utf-8"))

    splits = {}
    for split in ("train", "val", "test"):
        p = fdir / ("oracle_%s_footprint_GRCh38.bed.gz" % split)
        if not p.exists():
            raise ProbeRefusal("footprint missing: %s" % p)
        splits[split] = verify_split(split, p, side if split == "train" else None)

    rec = {
        "schema": SCHEMA,
        "date": "2026-11-23",
        "question": ("can the locus-level intersect be built -- is the S3 coordinate schema "
                     "pinnable, is the oracle footprint what its sidecar says, and what is the "
                     "random-element null at each tier?"),
        "our_side_s3_schema": cells,
        "oracle_side_footprint": splits,
        "frozen_tier_set": [{"name": n, "min_overlap_bp": k} for n, k in TIERS],
        "frozen_outcome_vocabulary": list(OUTCOMES),
        "separate_class_not_a_tier": {
            "name": "same_interval_different_sequence",
            "why": ("a ref/alt pair at one locus shares the interval and differs in sequence. The "
                    "exact-sequence hash join was blind to this BY CONSTRUCTION, so it is a class "
                    "of its own and must never be folded into the exact-match count."),
        },
        "assembly_rail": (
            "coordinates are read from NAMED assembly columns only. No coordinate is derived from "
            "or cross-checked against an element ID, and an absent hg38 coordinate is reported "
            "missing rather than filled from hg19."),
        "honest_limits": [
            "Schema figures come from a header read plus a bounded row sample, not a full parse, "
            "so per-cell absent-hg38 counts are SAMPLE rates and not cohort totals.",
            "The nulls assume UNIFORM element placement over a 3.1 Gbp genome. They are a "
            "reference for reading an observed rate, never a prediction about our corpus, which "
            "is cCRE/DNase/TSS-selected and non-uniform by construction.",
            "The chromosome-geography checks are TRUE BY CONSTRUCTION from the published split "
            "rule. They detect a broken parse and nothing else, and are blind to assembly and "
            "off-by-one error. They are labelled plumbing gates and are not evidence.",
            "This probe establishes FEASIBILITY and denominators only. It measures no overlap, "
            "and a pinned schema is not a validated join.",
        ],
    }

    Path(a.out).write_text(json.dumps(rec, indent=2) + "\n", encoding="utf-8")

    print("=== our side: S3 coordinate schema ===")
    for k, c in cells.items():
        print("  %-7s header r%d | %s | %s | %s | hg38-absent %d/%d | _F names %d"
              % (k, c["header_row_1based"], "+".join(c["assemblies_present"]),
                 c["chromosome_convention"], c["interval_convention"],
                 c["n_hg38_absent_in_sample"], c["n_sampled"],
                 c["n_names_ending_underscore_F"]))
    print()
    print("=== oracle side: footprint ===")
    for s, r in splits.items():
        print("  %-6s %7d runs | %10d bp (%.4f%%) | %2d chroms | merge %s"
              % (s, r["n_runs"], r["merged_bp"], r["base_fraction_pct"],
                 r["n_chromosomes"],
                 "clean" if not any(r["merge_integrity"].values()) else "DIRTY"))
    print()
    print("=== random-element null per tier (exact start-window union) ===")
    for n, d in splits["train"]["random_element_null_per_tier"].items():
        print("  %-20s %.4f%%   (naive sum %.4f%%, inflation %sx)"
              % (n, d["exact_pct_of_genome"], d["naive_sum_pct"], d["naive_inflation_x"]))
    print()
    print("-> %s" % a.out)
    return 0


def _self_check() -> int:
    """Pure-logic checks with a hand-computable answer. No files, no network."""
    assert to_chr("1") == "chr1" and to_chr("chr1") == "chr1"
    assert strip_chr("chr1") == "1" and strip_chr("1") == "1"

    # two runs 100 bp apart -- CLOSER than the 200 bp element, so the naive sum double-counts
    by = {"chr1": [(1000, 1200), (1300, 1500)]}
    for k in (1, 50, 100, 200):
        u = start_window_union(by, k)
        nv = naive_window_sum(by, k)
        assert u <= nv, (k, u, nv)
    # Hand-computed, and the first version of this assertion was WRONG -- recorded because the
    # arithmetic is the whole point of the function. At k=1 each run's start-window is
    # [a-200+1, b-1+1): run (1000,1200) -> [801,1200) = 399 starts, run (1300,1500) -> [1101,1500)
    # = 399. They OVERLAP from 1101, so the naive sum is 798 while the true union is [801,1500) =
    # 699 -- a 1.14x inflation on a pair only 100 bp apart.
    assert naive_window_sum(by, 1) == 798, naive_window_sum(by, 1)
    assert start_window_union(by, 1) == 699, start_window_union(by, 1)
    # at k=200 only a start exactly at `a` covers a full 200 bp of a 200 bp run
    assert start_window_union(by, 200) == 2, start_window_union(by, 200)

    # far-apart runs must NOT be double-counted away: union == sum
    far = {"chr1": [(0, 200), (10_000, 10_200)]}
    assert start_window_union(far, 1) == naive_window_sum(far, 1)

    # abutting runs are not a merge
    assert merge_integrity({"c": [(0, 100), (100, 200)]})["abutting"] == 1
    assert merge_integrity({"c": [(0, 100), (50, 200)]})["overlapping"] == 1
    assert merge_integrity({"c": [(100, 200), (0, 50)]})["unsorted"] == 1
    assert not any(merge_integrity({"c": [(0, 100), (200, 300)]}).values())

    # the frozen lists must be non-empty, or every downstream guard is vacuous
    assert TIERS and OUTCOMES
    print("self-check OK: %d tiers, %d outcomes" % (len(TIERS), len(OUTCOMES)))
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except ProbeRefusal as exc:
        print("REFUSED: %s" % exc, file=sys.stderr)
        sys.exit(2)
