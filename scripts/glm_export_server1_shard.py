"""Export a SELF-CONTAINED generator-training shard for an off-repo GPU box.

WHY A SHARD AND NOT THE WORKBOOK. The obvious handoff is "copy the lentiMPRA deposit to the other
machine and let its script load it". That is worse on four counts, three of them silent:

1. LEAKAGE SURFACE. The deposit carries folds 1-10. The generator trains on 1-6; 7-8 are the
   oracle's and 9-10 are held out, with fold 10 reserved FINAL_ONLY
   (`wiki/human_final_shard_reservation_2026-11-22.json`). Shipping the workbook puts the
   final-only shard on a third machine running a live session. This export contains folds 1-6 and
   NOTHING ELSE, so the reserved shard physically cannot be touched there.
2. DATA VINTAGE. The remote copy is on a different disk and may be a different download. A
   mismatched vintage changes the quintile boundaries and therefore the conditioning LABEL, which
   is the one quantity that must be identical across machines or the two runs are not comparable.
   The manifest pins the md5 of every source file and of the shard itself.
3. RE-IMPLEMENTED PARSING. `lentimpra.load_cell` is ~200 lines carrying two measured traps (the
   header row MOVES between sheets; the 230 nt INCLUDES cloning adaptors) plus a join-rate gate.
   Re-implementing that standalone is the "re-implemented harness must be gated against its
   runner" hazard with no runner to gate against. The remote script parses NOTHING: it reads two
   columns of text.
4. COST. Table S3 is 56 MB and openpyxl takes tens of seconds per read_only pass.

WHAT THE REMOTE SCRIPT THEREFORE NEVER DOES: derive a quintile boundary, strip an adaptor, or
decide a fold. Those are decided here, once, and travel as data.

THE N FILTER IS NOT COSMETIC. `lentimpra` admits `[ACGTN]+` by design, but
`dna_decode.glm.diffusion.encode` REFUSES a non-ACGT character -- correctly, since an N mapped to a
base teaches the model that an ambiguity code is real sequence. Unfiltered, a single N would abort
a remote run hours in. They are dropped here and the count is reported.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from dna_decode.glm.lentimpra import S3_FILE, S6_FILE, load_cell  # noqa: E402

ADAPTOR = 15            # 5' and 3', per the publishers' own column name
ELEMENT_LEN = 200       # the frozen bar's primary geometry
GENERATOR_FOLDS = (1, 2, 3, 4, 5, 6)
N_BINS = 5


def md5(path: Path) -> str:
    h = hashlib.md5()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def equal_count_bins(recs, n_bins: int):
    """Equal-COUNT strata by activity rank, mirroring `glm_human_stratified_null.quintiles`.

    Equal-count, NOT equal-width: activity is skewed, so equal-width bins leave the top stratum
    nearly empty. The remote job receives the bin as a COLUMN, so this rule never has to be
    reproduced anywhere else -- which is the point of exporting a shard at all.
    """
    s = sorted(recs, key=lambda r: r.activity)
    step = len(s) / n_bins
    out = []
    for i in range(n_bins):
        for r in s[int(i * step):int((i + 1) * step)]:
            out.append((i, r))
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data-dir", required=True,
                    help="lentiMPRA deposit directory. REQUIRED: lentimpra.DEFAULT_DIR points at "
                         "D:/dna_decode_cache/mpra_human, which does not exist on this host, and a "
                         "bare load reports absent-data for data that is present.")
    ap.add_argument("--cell", default="K562", choices=("K562", "HepG2"))
    ap.add_argument("--out-dir", required=True)
    a = ap.parse_args(argv)

    data_dir, out_dir = Path(a.data_dir), Path(a.out_dir)
    if not data_dir.is_dir():
        print("REFUSE: --data-dir %s is not a directory" % data_dir, file=sys.stderr)
        return 2
    out_dir.mkdir(parents=True, exist_ok=True)

    recs, st = load_cell(a.cell, data_dir)
    print("%s: joined %d at %d-%d bp (join rate %.4f)"
          % (a.cell, len(recs), st.length_min, st.length_max, st.join_rate))

    train = [r for r in recs if r.fold in GENERATOR_FOLDS]
    print("generator folds %s: %d records" % (list(GENERATOR_FOLDS), len(train)))

    kept, n_dropped_n, n_dropped_len = [], 0, 0
    for r in train:
        el = r.seq[ADAPTOR:len(r.seq) - ADAPTOR].upper()
        if len(el) != ELEMENT_LEN:
            n_dropped_len += 1
            continue
        if set(el) - set("ACGT"):
            n_dropped_n += 1
            continue
        kept.append((el, r))
    if n_dropped_len:
        print("REFUSE: %d elements are not %d bp after trimming %d bp from each end"
              % (n_dropped_len, ELEMENT_LEN, ADAPTOR), file=sys.stderr)
        return 3
    print("dropped %d with a non-ACGT base in the element; kept %d" % (n_dropped_n, len(kept)))
    if not kept:
        print("REFUSE: zero usable training elements", file=sys.stderr)
        return 4

    binned = equal_count_bins([r for _el, r in kept], N_BINS)
    el_by_name = {r.name: el for el, r in kept}

    shard = out_dir / ("train_%s_folds1-6.tsv" % a.cell)
    per_bin: dict[int, int] = {}
    with shard.open("w", encoding="utf-8", newline="\n") as fh:
        fh.write("bin\tsequence\n")
        for b, r in binned:
            fh.write("%d\t%s\n" % (b, el_by_name[r.name]))
            per_bin[b] = per_bin.get(b, 0) + 1

    ranges = {}
    for b in range(N_BINS):
        ys = sorted(r.activity for bb, r in binned if bb == b)
        ranges[b] = {"n": len(ys), "activity_min": round(ys[0], 6),
                     "activity_max": round(ys[-1], 6)}

    manifest_name = "shard_manifest_%s.json" % a.cell
    man = {
        "record": "glm-server1-training-shard-v1",
        "date": "2026-11-22",
        "cell": a.cell,
        "bar": "wiki/human_design_bar_2026-11-22.json",
        "geometry": "element_200nt_adaptors_stripped",
        "adaptor_bp_trimmed_each_end": ADAPTOR,
        "folds_included": list(GENERATOR_FOLDS),
        "folds_DELIBERATELY_ABSENT": {
            "7_8": "the oracle's training folds -- the judge must stay independent of the generator",
            "9": "WORKING_HELD_OUT",
            "10": "FINAL_ONLY, reserved by wiki/human_final_shard_reservation_2026-11-22.json",
            "why_absent_rather_than_present_and_unused": (
                "a fold that is not on the machine cannot be read by a script that is edited later. "
                "This is the leakage control, not a size optimisation."),
        },
        "conditioning": {
            "n_bins": N_BINS,
            "rule": "equal-COUNT activity strata by rank, mirroring "
                    "glm_human_stratified_null.quintiles",
            "bin_0_is": "lowest activity",
            "bin_4_is": "highest activity",
            "why_binned_not_raw": (
                "raw log2(RNA/DNA) is a sufficient RETRIEVAL KEY for 1.2-3.5% of records (726 K562 / "
                "1291 HepG2 have a unique value; 12.3 bits of target against 15.8 needed to name one "
                "sequence). Five bins carry 2.3 bits with a largest bin of 40,987, so retrieval is "
                "impossible by BOUND rather than unlikely by argument. Free, because the bar's "
                "decision quantity is already quintile monotonicity."),
            "per_bin": ranges,
        },
        "counts": {"joined_all_folds": len(recs), "generator_folds": len(train),
                   "kept": len(kept), "dropped_non_acgt": n_dropped_n,
                   "per_bin": dict(sorted(per_bin.items()))},
        "shard": {"file": shard.name, "md5": md5(shard), "bytes": shard.stat().st_size,
                  "format": "TSV, header `bin<TAB>sequence`, LF endings, bin in 0..4, "
                            "sequence exactly 200 ACGT"},
        "source_md5": {S3_FILE: md5(data_dir / S3_FILE), S6_FILE: md5(data_dir / S6_FILE)},
        "loader_stats": st.as_dict(),
        "honest_limits": [
            "Equal-count bins are computed over folds 1-6 ONLY. They are the generator's own "
            "conditioning scale and are NOT the scoring strata; the scorer re-derives its strata on "
            "the fold it scores. Reusing these boundaries at scoring time would be fitting the "
            "readout to the training set.",
            "The shard carries no fold column: every row is a generator fold by construction, so a "
            "column would only create the chance of filtering on it wrongly.",
        ],
    }
    (out_dir / manifest_name).write_text(json.dumps(man, indent=2), encoding="utf-8")

    print("-> %s  (%d bytes, md5 %s)" % (shard, shard.stat().st_size, man["shard"]["md5"]))
    print("-> %s" % (out_dir / manifest_name))
    print("per-bin: %s" % man["counts"]["per_bin"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
