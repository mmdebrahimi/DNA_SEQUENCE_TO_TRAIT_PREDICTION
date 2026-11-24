"""After removing the draft E. coli, how much verbatim DNA does a held-out MG1655 window still
share with the training corpus?

WHY THIS EXISTS. `scripts/glm_build_representative_corpus.py` removes the draft E. coli and the
representative manifest then asserts the leak gate is "STRICTLY STRONGER ... this corpus contains NO
E. coli, so it cannot hold a falsifier eval window under any scoring seed". That is a NAME-keyed
claim: it removes a SPECIES. What actually bounds a near-ortholog leak is a SEQUENCE-keyed fact --
how much verbatim DNA remains shared -- because identical or near-identical DNA under a different
organism's name is exactly what a name-keyed check cannot see.

The committed corpus's figure is known: held-out MG1655 windows share a median 55.0 bp verbatim with
it, max 148 of 150, 53.3% at or above 30 bp (`wiki/glm_copying_penalty_2026-11-22.json`). That is the
near-ortholog leak, measured. This probe asks what the same windows reach against the representative
corpus. If the removal worked it should fall toward the ~15 bp chance level that artifact also
records for uniform random sequence.

PAIRED BY CONSTRUCTION. The same subsampled windows are profiled against BOTH corpora, so the drop
is a paired per-window quantity and not two independent draws. Two medians cannot say whether a drop
is uniform or carried by a handful of windows; the paired difference can.

The committed arm doubles as a sanity check against the published 55.0. It is NOT a bit-level
reproduction: the published run's own subsample is not recoverable from its artifact, so a few bp of
difference there is sampling, not disagreement. Said here rather than discovered later.

`longest_shared_substring` is reused from scripts/glm_copying_penalty.py rather than reimplemented --
that is the part which is easy to get subtly wrong, and that module closes the non-ACGT separator
trap this probe would otherwise have to re-close. Its `profile` wrapper is deliberately NOT reused,
because it sorts its lengths before computing stats and so destroys the pairing.

Read-only. No network, no GPU, no model.
"""
from __future__ import annotations

import argparse
import json
import random
import statistics
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from dna_decode.glm.corpus import build_natural_set  # noqa: E402
from scripts.glm_copying_penalty import LONG_BP, SEP, longest_shared_substring  # noqa: E402

MG1655 = ROOT / "data" / "cache" / "refseq" / "GCF_000005845.2"
COMMITTED = ROOT / "wiki" / "glm_training_corpus_2026-11-21.txt"
REPRESENTATIVE = ROOT / "wiki" / "glm_training_corpus_representative_2026-11-22.txt"
PUBLISHED = ROOT / "wiki" / "glm_copying_penalty_2026-11-22.json"
OUT = ROOT / "wiki" / "glm_residual_leak_2026-11-24.json"
LENGTH = 150


def load_corpus(path: Path) -> list[str]:
    """The ONE way this probe reads a corpus file, matching the scorer's own reader."""
    return [ln.strip().upper() for ln in path.open(encoding="utf-8")
            if ln.strip() and not ln.startswith("#")]


def stats(vals: list[int]) -> dict:
    return {
        "n": len(vals),
        "median_bp": statistics.median(vals),
        "mean_bp": round(statistics.fmean(vals), 2),
        "max_bp": max(vals),
        "min_bp": min(vals),
        "long_bp_threshold": LONG_BP,
        "frac_at_or_above_long": round(sum(v >= LONG_BP for v in vals) / len(vals), 4),
    }


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--n", type=int, default=120,
                    help="windows to profile; 0 = all held-out windows, which costs hours")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", type=Path, default=OUT)
    a = ap.parse_args(argv)
    t0 = time.time()

    for p in (MG1655 / "annotations.gff3", COMMITTED, REPRESENTATIVE):
        if not p.exists():
            print(f"REFUSING: missing {p}", file=sys.stderr)
            return 2

    committed = load_corpus(COMMITTED)
    representative = load_corpus(REPRESENTATIVE)
    cset = set(committed)

    wins, _ = build_natural_set(MG1655, length=LENGTH, limit=None, seed=0, kind="promoter")
    mg = sorted({w.seq.upper() for w in wins})
    held = [s for s in mg if s not in cset]
    print(f"MG1655 promoter windows: {len(mg)} distinct | "
          f"held out of the committed corpus: {len(held)}")
    if not held:
        print("REFUSING: no held-out windows, so this cannot be the corpus the set was built "
              "against", file=sys.stderr)
        return 3

    # The representative corpus must be a SUBSET of the committed one, or the two arms are not
    # paired and no window-level difference is interpretable. Asserted, not assumed.
    extra = set(representative) - cset
    if extra:
        print(f"REFUSING: the representative corpus holds {len(extra)} windows absent from the "
              "committed one, so it is not a subset and the arms are not paired.", file=sys.stderr)
        return 4

    pool = held if a.n <= 0 else random.Random(a.seed).sample(held, min(a.n, len(held)))
    print(f"profiling {len(pool)} windows against both corpora (paired)\n")

    per_window: dict[str, list[int]] = {}
    for name, corpus in (("committed", committed), ("representative", representative)):
        hay = SEP.join(corpus) + SEP
        print(f"  {name}: {len(corpus)} windows")
        vals: list[int] = []
        for i, seq in enumerate(pool):
            vals.append(longest_shared_substring(seq, hay))
            if (i + 1) % 25 == 0:
                print(f"    {name}: {i + 1}/{len(pool)} ({time.time() - t0:.0f}s)", flush=True)
        per_window[name] = vals
        st = stats(vals)
        print(f"  -> median {st['median_bp']} bp | max {st['max_bp']} | "
              f"frac>=30bp {st['frac_at_or_above_long']}\n")

    drops = [c - r for c, r in zip(per_window["committed"], per_window["representative"])]
    paired = {
        "definition": "committed minus representative, per window; positive = the removal CUT the share",
        "median_drop_bp": statistics.median(drops),
        "mean_drop_bp": round(statistics.fmean(drops), 2),
        "max_drop_bp": max(drops),
        "n_windows_cut": sum(d > 0 for d in drops),
        "n_windows_unchanged": sum(d == 0 for d in drops),
        "n_windows_increased": sum(d < 0 for d in drops),
        "why_increase_would_be_a_defect": ("the representative corpus is a verified subset, so no "
                                           "window can gain a longer match; a nonzero count here "
                                           "means the pairing or the subset check is wrong"),
    }
    if paired["n_windows_increased"]:
        print(f"REFUSING: {paired['n_windows_increased']} windows share MORE with a subset corpus, "
              "which is impossible; the pairing is wrong.", file=sys.stderr)
        return 5

    pub = json.loads(PUBLISHED.read_text(encoding="utf-8")) if PUBLISHED.exists() else {}
    rec = {
        "record": "glm-residual-leak-v1",
        "date": time.strftime("%Y-%m-%d"),
        "question": ("after removing the draft E. coli, how much verbatim DNA does a held-out "
                     "MG1655 window still share with the training corpus"),
        "why_sequence_keyed": ("removing a species is a NAME-keyed operation; a near-ortholog leak "
                               "is a SEQUENCE-keyed fact. The representative manifest asserts the "
                               "gate is stronger because no E. coli remains. That is true, and it "
                               "is not the same claim as bounding residual verbatim overlap."),
        "n_mg1655_distinct": len(mg),
        "n_held_out_of_committed": len(held),
        "n_profiled": len(pool),
        "subsample_seed": a.seed,
        "paired": True,
        "representative_is_verified_subset_of_committed": True,
        "arms": {name: stats(vals) for name, vals in per_window.items()},
        "paired_drop": paired,
        "per_window": per_window,
        "published_committed_for_sanity": {
            "source": PUBLISHED.name,
            "real_held_out": (pub.get("real_held_out") or {}),
            "uniform_chance_median_bp": (pub.get("uniform_chance") or {}).get("median_bp"),
            "note": ("the published run's subsample is not recoverable from its artifact, so the "
                     "committed arm here is a sanity check on the same population, NOT a bit-level "
                     "reproduction"),
        },
        "honest_limits": [
            "A longest-verbatim-substring statistic bounds EXACT shared DNA. A near-ortholog that "
            "diverges every few bases can carry real homology while scoring low here, so a fall to "
            "chance level bounds the verbatim leak and does NOT prove the absence of homology.",
            "The default n=120 is a subsample of the held-out population, chosen to match the "
            "published figure's n rather than to be exhaustive.",
        ],
        "wall_seconds": round(time.time() - t0, 1),
    }
    a.out.write_text(json.dumps(rec, indent=2), encoding="utf-8")
    print(f"-> {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
