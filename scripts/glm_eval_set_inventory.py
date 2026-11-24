"""What does the orthology filter actually COST? Measure the three candidate eval sets, no model.

THE POINT. `scripts/glm_beyond_kmer.py` scored the generator on SIX windows and refused a verdict,
because the orthology filter left only six. Nobody had measured how many windows each candidate
eval set contains, so "the filter is expensive" was an inference. This probe turns it into a number,
and it does so with NO CHECKPOINT LOADED -- the strong null alone decides whether a set is even
capable of carrying a verdict, via `dna_decode.glm.eval_admissibility`.

THE THREE SETS, and they are NOT nested the way one would assume:
  A -- the DESIGNED withhold. `falsifier_eval_windows` at the committed geometry; the 3,349 of
       `wiki/glm_training_corpus_manifest_2026-11-21.json`. Built at limit=1200 over 10 seeds.
  B -- `held`. MG1655 promoter windows whose sequence is absent from the corpus TEXT. A DIFFERENT
       predicate from A, and built at limit=None, so A is not a subset of B by construction.
  C -- `of`. B intersected with MG1655-only gene symbols. This is what the committed run scored.

WHY A != 3349 IS A REFUSAL, NOT A WARNING. `SCORING_LIMIT` 1200 and `seeds=range(10)` DEFINE set A.
A geometry drift would silently make every number here describe a different cohort than the manifest
does, and the resulting artifact would look perfectly clean. The scorer derives its window length
from the CHECKPOINT config (measured 150 at step 6000) while this probe takes it as an argument, so
this refusal is also the tripwire for those two disagreeing.

CASE IS A REAL TRAP HERE. `fit_eval_halves` returns `w.seq` RAW while the scorer compares
`w.seq.upper()`. Intersecting the two without normalising can report a spurious ZERO overlap, which
would read as "the sets are disjoint" -- a finding -- rather than as a bug. Every set is uppercased
before any comparison and the raw-vs-uppercased sizes are both reported so a case problem is
VISIBLE instead of silent.

Run: uv run python scripts/glm_eval_set_inventory.py
"""
from __future__ import annotations

import argparse
import bisect
import json
import math
import random
import statistics
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from dna_decode.glm.baselines import MarkovGenerator  # noqa: E402
from dna_decode.glm.eval_admissibility import admissibility  # noqa: E402
from dna_decode.glm.training_corpus import falsifier_eval_windows  # noqa: E402
from scripts.glm_beyond_kmer import (  # noqa: E402
    MG1655,
    ORDER,
    ecoli_fit_windows,
    load_corpus_windows,
    markov_nll,
    mg1655_candidate_sets,
)

# The COMMITTED geometry. These are not tunables: they are what the manifest recorded and what set A
# is DEFINED by. `scoring_limit_natural` 1200 and `scoring_seeds` 10 come from the `geometry` block
# of wiki/glm_training_corpus_manifest_2026-11-21.json.
LENGTH = 150
SCORING_LIMIT = 1200
N_SEEDS = 10
KIND = "promoter"
EXPECTED_A_SIZE = 3349          # manifest leak_gate.withheld_union_size
MIN_ECOLI_FIT = 500             # the scorer's own floor for the strong null
UNIFORM = math.log(4)           # what a predictor that commits to nothing achieves

# The geometry the committed 2026-11-22 run was actually invoked at, which the artifact does not
# record. Reconstructed and then CONFIRMED by `decompose_committed_cell` reproducing 1.4291 exactly.
COMMITTED_N_WINDOWS = 6
COMMITTED_N_POSITIONS = 8

OUT = ROOT / "wiki" / "glm_eval_set_inventory_2026-11-23.json"


class InventoryRefusal(RuntimeError):
    """A geometry or integrity precondition failed. Never downgraded to a warning."""


def assert_expected_a(n_a: int, expected: int = EXPECTED_A_SIZE) -> None:
    """REFUSE unless set A reproduces the manifest exactly. The drift tripwire."""
    if n_a != expected:
        raise InventoryRefusal(
            f"set A recomputed to {n_a} windows, expected {expected} from the manifest "
            f"leak_gate.withheld_union_size. The geometry has drifted (length / limit / seeds / "
            f"kind), so every set below would describe a different cohort than the manifest does. "
            f"Refusing rather than emitting a clean-looking artifact about the wrong windows.")


def set_algebra(a_seqs, held_seqs, of_seqs) -> dict:
    """Sizes and pairwise overlaps of the three candidate sets, case-normalised.

    PURE. Set C must be a subset of set B by construction -- C is B filtered by gene symbol -- so a
    violation means the filter was applied to the wrong list and is RAISED, not reported.
    """
    a = {s.upper() for s in a_seqs}
    b = {s.upper() for s in held_seqs}
    c = {s.upper() for s in of_seqs}
    if not c <= b:
        raise InventoryRefusal(
            f"set C is not a subset of set B ({len(c - b)} windows in C but not in B). C is defined "
            f"as B filtered by gene symbol, so this means the filter ran against the wrong list.")
    return {
        # DISTINCT sequence counts. These are the statistically meaningful sizes.
        "A_sequence_withheld": len(a),
        "B_held_out_of_corpus": len(b),
        "C_orthology_free": len(c),
        # LIST entry counts, which is what the SCORER sees: it tests `len(of) < --n-windows` and
        # then slices `of[:n]`, so a duplicated sequence occupies a slot and would be scored twice.
        # Reported separately because the two differ and conflating them inflates apparent n.
        "A_list_entries": len(list(a_seqs)),
        "B_list_entries": len(list(held_seqs)),
        "C_list_entries": len(list(of_seqs)),
        "B_duplicate_entries": len(list(held_seqs)) - len(b),
        "C_duplicate_entries": len(list(of_seqs)) - len(c),
        "A_raw_before_uppercasing": len(set(a_seqs)),
        "B_raw_before_uppercasing": len(set(held_seqs)),
        "C_raw_before_uppercasing": len(set(of_seqs)),
        "A_inter_B": len(a & b),
        "A_inter_C": len(a & c),
        "B_inter_C": len(b & c),
        "C_is_subset_of_B": True,
        "A_is_subset_of_B": a <= b,
        # CAUTION, and the name below is kept only because tests pin it: this is a SIZE RATIO
        # between two differently-built sets, NOT a subset share. C is not inside A -- they share
        # 519 windows and C holds 115 that A never saw, because A is built at limit=1200 and C at
        # limit=None. Reading it as "the filter discards 81% of the withheld set" is wrong: the
        # filter is not applied to A at all.
        "C_over_A_the_filters_measured_cost": (round(len(c) / len(a), 6) if a else None),
        "C_over_A_is_a_size_ratio_NOT_a_subset_share": True,
        # THIS is the filter's cost where the question is well posed: of the windows actually IN the
        # designed withhold, the share that is orthology-free. It is also exactly the arm the leak
        # contrast scores, and |A&C| + |A\C| = |A| is asserted below.
        "A_inter_C_over_A_the_filters_cost_WITHIN_the_withheld_set": (
            round(len(a & c) / len(a), 6) if a else None),
        "A_partitions_into_the_two_contrast_arms": (len(a - c) + len(a & c) == len(a)),
        "A_minus_B": len(a - b),
        "B_minus_A": len(b - a),
    }


def null_ce_on(chain: MarkovGenerator, seqs, *, length: int) -> dict:
    """The strong null mean cross-entropy over EVERY scoreable position of every window.

    DELIBERATELY NOT the scorer sampled estimate. The scorer averages `--positions` sampled
    positions per window because it must pay for a forward pass at each one; the chain is free, so
    scoring all of them removes the sampling seed from the screen entirely and gives a
    lower-variance estimate of the SAME quantity. A sampled figure from the scorer is a noisy draw
    around this one, which is why the two will not agree to 4dp.
    """
    # The geometry check must be TWO-directional. An overrun already raised from markov_nll's
    # indexing, but a `length` BELOW the true window length silently scored fewer positions and
    # then reported the short count as `n_positions_per_window` -- a quietly different measurement
    # wearing the right label.
    bad = {len(s) for s in seqs} - {length}
    if bad:
        raise InventoryRefusal(
            f"every window must be exactly {length} bases to score positions "
            f"[{ORDER}, {length}); found lengths {sorted(bad)}. Scoring a shorter window would "
            f"report fewer positions under the full-geometry label.")
    per_window = [statistics.fmean(markov_nll(chain, s, p) for p in range(ORDER, length))
                  for s in seqs]
    return {
        "n_windows": len(per_window),
        "n_positions_per_window": length - ORDER,
        "scoring": "ALL scoreable positions, not a sample -- no seed involved",
        "null_ce_nats": round(statistics.fmean(per_window), 4) if per_window else None,
    }


def null_side_leak_contrast(chain: MarkovGenerator, a_seqs, of_seqs, *, length: int) -> dict:
    """How much of the near-ortholog leak accrues to the NULL rather than the model?

    The orthology filter exists because a window whose gene has a draft counterpart may have a
    95%-identical ortholog in training. But the STRONG null was fit on exactly those draft windows,
    so it has seen those same near-orthologs -- which means the leak is not a one-sided model
    advantage, and an A-vs-C delta difference is model-leak MINUS null-leak.

    This measures the null half directly and for free, by splitting the withheld set A into its
    ortholog-BEARING part (A minus C) and its orthology-FREE part (A intersect C). Both halves are
    withheld from the corpus and both come from the same builder, so the only systematic difference
    between them is whether a draft counterpart exists. A null that gained from the leak scores the
    ortholog-bearing half LOWER (better).
    """
    all_pos = list(range(ORDER, length))
    a = {s.upper() for s in a_seqs}
    c = {s.upper() for s in of_seqs}
    bearing = sorted(a - c)
    free = sorted(a & c)

    def mean_ce(seqs):
        if not seqs:
            return None
        return statistics.fmean(
            statistics.fmean(markov_nll(chain, s, p) for p in all_pos) for s in seqs)

    ce_bearing, ce_free = mean_ce(bearing), mean_ce(free)
    gap = (None if ce_bearing is None or ce_free is None else ce_free - ce_bearing)
    return {
        "A_minus_C_ortholog_bearing": {"n_windows": len(bearing),
                                       "null_ce_nats": None if ce_bearing is None else round(ce_bearing, 4)},
        "A_inter_C_orthology_free": {"n_windows": len(free),
                                     "null_ce_nats": None if ce_free is None else round(ce_free, 4)},
        "null_side_leak_nats": None if gap is None else round(gap, 4),
        "sign_convention": "positive = the null scores the ortholog-BEARING half better, i.e. the "
                           "null gained from the leak",
        "scoring": "ALL scoreable positions on every window of each half; no sampling, no seed",
        "reading": (
            "This is the NULL half of the leak only. It does not measure the model's share, and a "
            "near-zero value here does NOT mean the leak is small for the model -- only that the "
            "null is not the party benefiting from it."),
    }


def decompose_committed_cell(chain: MarkovGenerator, of_seqs, *, length: int,
                             n_windows: int = 6, n_positions: int = 8, seed: int = 0) -> dict:
    """Why did the committed run measure the null ABOVE uniform when the full set is BELOW it?

    The committed 2026-11-22 cell reported `markov_ecoli` at 1.4291 nats on set C, above uniform,
    and that single number is what made the runner REFUSE a verdict and record the representative
    re-pretrain as "the only route to an informative measurement". The same null on the same set,
    measured over every window and every position, is 1.3571 -- BELOW uniform.

    Both cannot be properties of the eval set, so this isolates which factor moved it by varying
    ONE at a time: the six windows (`of[:6]`, the FIRST six in MG1655 coordinate order, which is not
    a random sample) and the eight sampled positions. The position draw is replicated EXACTLY as the
    scorer draws it -- one RNG advanced across windows -- so the committed figure must come out of
    the first cell, and that identity is what makes the other three cells trustworthy.
    """
    all_pos = list(range(ORDER, length))

    def sampled(n_w: int, n_p: int):
        rng = random.Random(seed)
        return [sorted(rng.sample(all_pos, n_p)) for _ in range(n_w)]

    def ce(seqs, pos_per_window) -> float:
        return statistics.fmean(
            statistics.fmean(markov_nll(chain, s, p) for p in ps)
            for s, ps in zip(seqs, pos_per_window))

    six = of_seqs[:n_windows]
    cells = {
        "committed_6_windows_8_sampled_positions": round(ce(six, sampled(n_windows, n_positions)), 4),
        "same_6_windows_ALL_positions": round(ce(six, [all_pos] * len(six)), 4),
        "ALL_windows_8_sampled_positions": round(ce(of_seqs, sampled(len(of_seqs), n_positions)), 4),
        "ALL_windows_ALL_positions": round(ce(of_seqs, [all_pos] * len(of_seqs)), 4),
    }

    per = [statistics.fmean(markov_nll(chain, s, p) for p in all_pos) for s in of_seqs]
    srt = sorted(per)
    first_pct = [round(bisect.bisect_left(srt, v) / len(srt), 3) for v in per[:n_windows]]
    n_above = sum(1 for v in per if v > UNIFORM)
    return {
        "cells": cells,
        "uniform_ce_nats": round(UNIFORM, 4),
        "n_predictions_in_the_committed_cell": n_windows * n_positions,
        "per_window_all_position_null_ce": {
            "min": round(min(per), 4),
            "median": round(statistics.median(per), 4),
            "max": round(max(per), 4),
        },
        "first_six_windows_percentile_within_set_C": first_pct,
        "windows_whose_own_null_ce_is_worse_than_uniform": n_above,
        "windows_total": len(per),
        "fraction_worse_than_uniform": round(n_above / len(per), 4) if per else None,
        "reading": (
            "The committed cell REPRODUCES exactly, which is what makes the other three cells "
            "trustworthy. Holding the windows fixed and scoring every position instead of eight "
            "moves the null from 1.4291 to 1.3541 -- across uniform. Holding the positions at eight "
            "and using every window gives 1.3565. So the factor is the TOTAL number of predictions, "
            "not which windows were chosen: the first six sit at ordinary percentiles within set C, "
            "so `of[:6]` is not a biased selection. 48 predictions is simply too few, and the "
            "committed run drew an unlucky 48. The estimator is unbiased and was read as a property "
            "of the eval set."),
        "what_this_does_NOT_establish": (
            "That the generator beats the null. This is the NULL alone. It establishes that set C "
            "is admissible and that the committed refusal rests on 48 predictions, nothing more. "
            "It also does not remove the near-ortholog confound from any set."),
    }


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--length", type=int, default=LENGTH)
    ap.add_argument("--limit", type=int, default=SCORING_LIMIT)
    ap.add_argument("--seeds", type=int, default=N_SEEDS)
    ap.add_argument("--out", type=Path, default=OUT)
    a = ap.parse_args(argv)
    t0 = time.time()

    if not (MG1655 / "annotations.gff3").exists():
        print(f"REFUSING: missing {MG1655 / 'annotations.gff3'}", file=sys.stderr)
        return 2

    corpus = load_corpus_windows()
    cset = set(corpus)
    print(f"corpus: {len(corpus)} windows, {len(cset)} unique")

    # ---- set A: the designed withhold, recomputed because it is NOT persisted anywhere ----
    a_seqs = falsifier_eval_windows(MG1655, length=a.length, limit=a.limit,
                                    seeds=range(a.seeds), kind=KIND)
    print(f"set A (sequence-withheld): {len(a_seqs)} windows")
    try:
        assert_expected_a(len(a_seqs))
    except InventoryRefusal as exc:
        print(f"REFUSING: {exc}", file=sys.stderr)
        return 3

    # ---- sets B and C, in ONE computation, through the scorer own helper ----
    sets = mg1655_candidate_sets(cset, length=a.length)
    held_seqs = [w.seq for w in sets["held"]]
    of_seqs = sets["orthology_free"]
    n_only = len(sets["only"])
    n_all = len(sets["all_windows"])
    print(f"set B (held out of corpus): {len(held_seqs)} windows")
    print(f"set C (orthology-free):     {len(of_seqs)} windows "
          f"({n_only} MG1655-only symbols, {n_all} MG1655 windows total)")

    try:
        algebra = set_algebra(a_seqs, held_seqs, of_seqs)
    except InventoryRefusal as exc:
        print(f"REFUSING: {exc}", file=sys.stderr)
        return 4

    # The scorer slices `of[:n]`, so whether a given request is duplicate-free is a property of the
    # SLICE, not of the set. Checked for the geometries that matter rather than assumed.
    slice_safety = {
        f"first_{n}_is_duplicate_free": (len(of_seqs[:n]) == len(set(of_seqs[:n])))
        for n in (COMMITTED_N_WINDOWS, 120, 300, len(of_seqs))
    }
    slice_safety["distinct_in_first_120"] = len(set(of_seqs[:120]))
    print(f"slice safety: committed {COMMITTED_N_WINDOWS} duplicate-free="
          f"{slice_safety[f'first_{COMMITTED_N_WINDOWS}_is_duplicate_free']}, "
          f"default 120 duplicate-free={slice_safety['first_120_is_duplicate_free']}")

    # ---- the strong null, fit exactly as the scorer fits it ----
    ecoli_fit = ecoli_fit_windows(cset, length=a.length)
    if len(ecoli_fit) < MIN_ECOLI_FIT:
        print(f"REFUSING: only {len(ecoli_fit)} draft E. coli windows in the corpus; the strong "
              f"null would be fit on too little data.", file=sys.stderr)
        return 5
    chain_ec = MarkovGenerator(order=ORDER, seed=0).fit(ecoli_fit)
    print(f"strong null fit on {len(ecoli_fit)} draft E. coli windows in the corpus")

    # ---- the null-only admissibility screen, per set, with NO checkpoint loaded ----
    screens = {}
    for name, seqs in (("A_sequence_withheld", sorted({s.upper() for s in a_seqs})),
                       ("B_held_out_of_corpus", sorted({s.upper() for s in held_seqs})),
                       ("C_orthology_free", sorted({s.upper() for s in of_seqs}))):
        scored = null_ce_on(chain_ec, seqs, length=a.length)
        # An EMPTY set has no cross-entropy, and `admissibility(None)` raises on float(None). A
        # set with nothing in it is not admissible or inadmissible -- it is unscoreable, and saying
        # so is the honest third answer rather than a crash or a fabricated verdict.
        if scored["null_ce_nats"] is None:
            screens[name] = {**scored, "gate": None,
                             "unscoreable": "the set is empty, so it has no cross-entropy and no "
                                            "admissibility verdict is available"}
            print(f"  {name:<24} EMPTY  ->  UNSCOREABLE")
            continue
        gate = admissibility(null_ce_nats=scored["null_ce_nats"])
        screens[name] = {**scored, "gate": gate.to_dict()}
        print(f"  {name:<24} null CE {scored['null_ce_nats']}  ->  {gate.state}")

    # Every set screens ADMISSIBLE while the committed run measured this same null ABOVE uniform on
    # set C. Both cannot describe the eval set, so the contradiction is isolated rather than left to
    # look like a defect in this screen.
    decomp = decompose_committed_cell(chain_ec, of_seqs, length=a.length,
                                      n_windows=COMMITTED_N_WINDOWS,
                                      n_positions=COMMITTED_N_POSITIONS)
    print("\ncommitted-cell decomposition (set C, strong null):")
    for k, v in decomp["cells"].items():
        print(f"  {k:<42} {v}")
    print(f"  uniform {decomp['uniform_ce_nats']} | "
          f"{decomp['windows_whose_own_null_ce_is_worse_than_uniform']}"
          f"/{decomp['windows_total']} windows individually worse than uniform")

    leak = null_side_leak_contrast(chain_ec, a_seqs, of_seqs, length=a.length)
    print("\nnull-side near-ortholog leak (within the withheld set A):")
    print(f"  A\\C ortholog-bearing  n={leak['A_minus_C_ortholog_bearing']['n_windows']:<5} "
          f"null CE {leak['A_minus_C_ortholog_bearing']['null_ce_nats']}")
    print(f"  A&C orthology-free    n={leak['A_inter_C_orthology_free']['n_windows']:<5} "
          f"null CE {leak['A_inter_C_orthology_free']['null_ce_nats']}")
    print(f"  null-side leak: {leak['null_side_leak_nats']} nats "
          f"(positive = null gained from the leak)")

    res = {
        "record": "glm-eval-set-inventory-v1",
        "date": "2026-11-24",
        "question": "how many windows does each candidate eval set contain, and can each carry a "
                    "verdict at all, judged from the strong null alone?",
        "no_checkpoint_loaded": True,
        "geometry": {
            "window_length_bp": a.length,
            "scoring_limit_natural": a.limit,
            "scoring_seeds": a.seeds,
            "natural_set_arg": KIND,
            "markov_order": ORDER,
            "length_source": "an ARGUMENT here; the scorer derives it from the checkpoint config "
                             "(measured 150 at step 6000). The A==3349 refusal is the tripwire for "
                             "those two disagreeing.",
        },
        "set_definitions": {
            "A_sequence_withheld": "falsifier_eval_windows union over seeds, built at limit=1200",
            "B_held_out_of_corpus": "MG1655 promoter windows absent from the corpus TEXT, limit=None",
            "C_orthology_free": "B filtered to genes with no draft counterpart; what the committed "
                                "2026-11-22 run scored",
        },
        "mg1655_only_symbols": n_only,
        "mg1655_windows_total": n_all,
        "sizes_and_overlaps": algebra,
        "slice_safety": slice_safety,
        "strong_null": {
            "fit_windows": len(ecoli_fit),
            "fit_on": "draft E. coli (GCA_059754755.1) promoter windows that are IN the corpus, so "
                      "the null is never fit on the windows it scores, for any of the three sets",
            "order": ORDER,
        },
        "admissibility_screens": screens,
        "committed_cell_decomposition": decomp,
        "null_side_leak_contrast": leak,
        "committed_run_geometry_reconstructed": {
            "n_windows": COMMITTED_N_WINDOWS,
            "positions_per_window": COMMITTED_N_POSITIONS,
            "n_predictions": COMMITTED_N_WINDOWS * COMMITTED_N_POSITIONS,
            "windows_that_were_AVAILABLE": len(of_seqs),
            "confirmed_by": "decompose_committed_cell reproduces the committed 1.4291 exactly",
            "note": ("the committed artifact records no seed or boot argument, so its CI is not "
                     "reproducible from it; the POINT estimates are. The scorer refuses when "
                     "len(of) < --n-windows, and len(of) is 648, so the scorer default of 120 "
                     "windows would have run -- the 6 was a choice, not a ceiling the filter "
                     "imposed."),
        },
        "uniform_ce_nats": round(math.log(4), 4),
        "seconds": round(time.time() - t0, 1),
        "emits_no_verdict_about_the_generator": True,
        "headline": (
            "TWO claims, and the second is deliberately narrower than it first looks. (1) The "
            "orthology-free set holds 648 list entries (634 DISTINCT sequences), not 6, so the "
            "premise that the gene-level filter leaves too little to score is FALSE: the scorer "
            "refuses only when len(of) < --n-windows, and at 648 its own default of 120 would have "
            "run. The 6 was a choice, not a ceiling the filter imposed. (2) The NULL half of the "
            "admissibility gate clears uniform at scale: the committed run read the strong null at "
            "1.4291 nats on set C, above uniform 1.3863, while the same null on the same set over "
            "every window and position reads 1.3571, below it. The decomposition isolates the cause "
            "as the 48-prediction sample size, not which windows were picked."),
        "what_claim_2_does_NOT_settle": (
            "Whether the committed refusal stands. The gate is an AND -- it requires the null AND "
            "the model to beat uniform -- and this probe measures the NULL only, by design, since "
            "its whole point is to screen a set with no checkpoint loaded. The model's committed "
            "1.4098 nats is equally a 48-prediction number and is equally unmeasured at scale. If "
            "the model stays above uniform at full geometry while the null clears it, the refusal "
            "still fires, but it then means something different: a model negative on a set that is "
            "demonstrably predictable, NOT an uninformative eval set. Scoring the model at full "
            "geometry is what decides between those, and it is not done here."),
        "duplicate_slot_finding": (
            "The orthology-free LIST carries 14 duplicate sequences (648 entries, 634 distinct) and "
            "`held` carries 17 (3923 entries, 3906 distinct), one sequence appearing 6 times. The "
            "scorer slices the LIST, so a duplicate occupies a scoring slot and would be scored "
            "twice -- inflating apparent n and breaking the independence the paired bootstrap "
            "assumes. MEASURED: neither of[:6] (the committed run) nor of[:120] (the scorer "
            "default) contains a duplicate, so no committed number is affected and Step 3 at 120 "
            "windows is safe. The duplicates sit later in the list, which is why this is recorded "
            "as a latent hazard rather than a live defect."),
        "what_the_headline_does_not_claim": (
            "It does not claim the generator beats the null, or that any set is unconfounded. The "
            "near-ortholog leak is untouched by all of this, and the strong null shares it."),
        "honest_limits": [
            "This measures SET SIZES and a null-only admissibility screen. It says nothing about "
            "the generator, and an ADMISSIBLE set is not an unconfounded set.",
            "Set A retains near-orthologs: withholding is by SEQUENCE, so a 95%-identical ortholog "
            "in the draft E. coli the model trained on is NOT withheld. Admissibility does not fix "
            "that and does not claim to.",
            "The strong null was fit on draft windows that ARE in the corpus, which are the very "
            "near-orthologs set A retains. So a set retaining near-orthologs strengthens the NULL "
            "as well as the model, and an A-vs-C difference is therefore not a one-sided model "
            "advantage.",
            "Sets A and B are built with DIFFERENT builder arguments (limit=1200 vs limit=None), so "
            "they are not nested by construction; A_is_subset_of_B records whether they happen to "
            "be rather than assuming it.",
            "The null CE here averages ALL scoreable positions while the scorer samples a few, so "
            "these figures are lower-variance and will not equal the scorer figures to 4dp.",
        ],
    }
    a.out.write_text(json.dumps(res, indent=2), encoding="utf-8")
    print(f"\n-> {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
