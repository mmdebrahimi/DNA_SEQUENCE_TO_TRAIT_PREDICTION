"""Where do the ~0.03 nats of promoter signal live -- in CONTEXT, in POSITION, or in COMPOSITION?

WHY THIS EXISTS. The mission's highest-budget gap (`generator-beats-its-baseline-panel`, strategy
budget 4) requires the generator to beat a MANDATORY five-baseline panel: mutation-count, BLOSUM62,
one-hot ridge, site-independent, shuffled null. Two of those five -- **site-independent** and
**shuffled null** -- are well defined for a generative per-position cross-entropy readout and have
NEVER been measured on it. The committed scorer compares the generator against order-3 Markov chains
and uniform only. So the panel has a hole exactly where it is cheapest to fill.

WHAT THIS ADDS, AND WHY IT IS NOT A RELABELLING. The order-3 chain uses three bases of CONTEXT and
knows nothing about WHERE in the window it sits. A site-independent model is its exact complement:
it knows the position and nothing about context. Promoters are positionally anchored (-10, -35,
spacer), so these two nulls are not nested and neither dominates by construction. Adding them
decomposes the ~0.03 nats of available signal into three orthogonal sources:

    context-only   = markov (order 3)          -- committed, re-derived here as the reconcile gate
    position-only  = site_independent          -- NEW
    composition    = shuffled_null (order 0)   -- NEW, the floor both of the above must clear

ZERO TORCH, BY CONSTRUCTION. The generator's per-window cross-entropy is already PERSISTED in the
`per_window` block of all six artifacts in `wiki/glm_repretrain_scored_2026-11-25/`, together with
the `run_args` that generated the windows and positions. So the new columns are computed on the
IDENTICAL positions and paired against the stored generator numbers without loading a checkpoint.
The script asserts `torch` never enters `sys.modules`.

THE RECONCILE GATE IS THE WHOLE BASIS OF TRUST. Before any new column is added, every window's
`markov_ecoli`, `markov_pooled` and `uniform` are RECOMPUTED from the recorded seeds and compared to
the persisted values at 6 decimals -- the precision the artifacts store. That one check validates
the corpus read, the eval-set reconstruction, the window selection, the position sampling and the
chain fit TOGETHER. A mismatch means the reconstruction is wrong and the script REFUSES (exit 2)
rather than reporting a number against windows it cannot prove are the scored ones.

WHAT THIS DELIBERATELY DOES NOT DO. It does NOT write `wiki/glm_generator_baseline_panel.json`.
That path is the mission's live gating predicate for the budget-4 gap AND the self-firing expiry
argument of `north-star-step4-animals-deferred` in `soraya_runs/constraint-expiry.json` AND one of
the five non-vacuity fixtures in `tests/test_mission_gap_predicates.py`. Creating it would close a
mission gap, re-open the north star's step 4, and flip a test fixture -- three authority calls, none
of them mine. This writes a dated artifact and records the definability split as a PROPOSAL.

HONEST LIMITS, stated here rather than discovered later.
* Two of the five mandated baselines (mutation-count, BLOSUM62) are NOT defined on this readout and
  are recorded as such with reasons, not silently dropped. One (one-hot ridge) is defined on the
  ACTIVITY criterion and was measured there; it is cross-referenced, not re-measured.
* `shuffled_null` here is the FIT-SIDE reading: a site-independent model whose positional
  information has been destroyed, which collapses to the global base composition. The other reading
  -- the GENERATOR scored on shuffled eval windows -- needs torch and is recorded as un-run.
* Both new nulls are fit on the same windows as their Markov counterparts, with the SAME add-one
  smoothing, so a difference between them is not a smoothing artifact.
* These are infilling cross-entropies. They are necessary and diagnostic, NOT generation quality.
"""
from __future__ import annotations

import argparse
import json
import math
import random
import statistics
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

# The scorer's OWN functions, imported rather than reimplemented so the reconcile gate exercises the
# real machinery. `dna_decode.glm.diffusion` imports torch lazily (inside functions), so this import
# does not load it -- asserted below, because that property is load-bearing and could regress.
from dna_decode.glm.baselines import MarkovGenerator  # noqa: E402
from scripts.glm_beyond_kmer import (  # noqa: E402
    BASES,
    ORDER,
    EvalSetUnavailable,
    ecoli_fit_windows,
    load_corpus_windows,
    markov_nll,
    mg1655_candidate_sets,
    select_eval_windows,
)

if "torch" in sys.modules:  # pragma: no cover - a guard on an import-time property
    raise RuntimeError(
        "torch was imported; this script must stay torch-free so it can run beside the test suite")

MG1655_DIR = ROOT / "data" / "cache" / "refseq" / "GCF_000005845.2"
SCORED_DIR = ROOT / "wiki" / "glm_repretrain_scored_2026-11-25"
OUT = ROOT / "wiki" / "glm_baseline_panel_ce_2026-11-28.json"
LENGTH = 150          # the committed window geometry; asserted against the corpus below
SMOOTH = 1.0          # add-one, matching MarkovGenerator's smoothing exactly
RECONCILE_DECIMALS = 6
SCORING_LIMIT = 1200  # mirrors glm_beyond_kmer's set-A geometry
N_SEEDS = 10

#: FROZEN BEFORE ANY NUMBER EXISTED. All four outcomes are named, including the two that are worse
#: for the project than the status quo.
PREREGISTERED = {
    "question": ("does the generator beat the POSITION-only and COMPOSITION-only nulls it has "
                 "never been scored against, and where do the ~0.03 available nats live?"),
    "outcomes": {
        "GENERATOR_CLEARS_FULL_PANEL": (
            "the generator beats markov_ecoli AND site_independent AND shuffled_null, each with a "
            "paired 95% CI excluding zero, on the criterion's own eval set"),
        "NEW_NULL_DOMINATES_GENERATOR": (
            "at least one NEW null (site_independent or shuffled_null) BEATS the generator with a "
            "paired CI excluding zero -- strictly worse than the existing negative, because those "
            "nulls carry no context at all"),
        "GENERATOR_BEATS_ONLY_THE_WEAK_NULLS": (
            "the generator beats shuffled_null and/or site_independent but still loses to "
            "markov_ecoli -- the existing negative, now with its floor located"),
        "UNRESOLVED": "every relevant paired CI spans zero",
    },
    "resolution_rule": ("a paired bootstrap 95% CI excluding zero, INHERITED verbatim from the "
                        "committed scorer's `resolved` field rather than a new margin invented here"),
    "secondary_prediction": (
        "site_independent should be WEAKER than markov_ecoli on promoters: "
        "wiki/glm_frame_headroom_probe_2026-11-25.json measured the full-position headroom over a "
        "homogeneous chain at only +0.0042 nats. An EXACT tie between them is a misalignment "
        "signature, not a finding, and is reported as such."),
    "why_not_a_new_margin": (
        "the arm bars in this project already resolve on CI-excludes-zero; inventing a nats margin "
        "here would let the panel clear a bar its own siblings never had to."),
}

#: The five mandated baselines, each classified for THIS readout with a reason. A baseline that is
#: undefined is recorded as undefined -- never silently dropped, and never coerced into a number.
MANDATED_PANEL = {
    "site_independent": {
        "status": "MEASURED_HERE",
        "definition": ("per-position base frequencies, add-one smoothed, fit on the same windows "
                       "as the Markov null it is compared against"),
    },
    "shuffled_null": {
        "status": "MEASURED_HERE",
        "definition": ("a site-independent model with its positional information destroyed, which "
                       "collapses to the global base composition (order 0)"),
        "variant_not_run": "the generator scored on column-shuffled EVAL windows; needs torch",
    },
    "one_hot_ridge": {
        "status": "DEFINED_ON_A_DIFFERENT_CRITERION",
        "why": ("a supervised regressor needs a continuous target; masked-nucleotide cross-entropy "
                "has none. It IS defined on the ACTIVITY criterion and was measured there."),
        "measured_at": "wiki/glm_onehot_ridge_backrun_2026-11-28.json",
    },
    "mutation_count": {
        "status": "UNDEFINED_ON_THIS_READOUT",
        "why": ("counting mutations needs a variant and a reference; the eval windows are natural "
                "held-out sequence, not variants of anything."),
    },
    "blosum62": {
        "status": "UNDEFINED_ON_THIS_READOUT",
        "why": ("BLOSUM62 scores AMINO-ACID substitutions; promoter windows are non-coding, so "
                "there is no residue to substitute."),
    },
}


class PanelError(RuntimeError):
    pass


def compute_verdict(cells: dict, arm: str) -> dict:
    """The frozen bar, applied by CODE to one arm's cells. Ordered; first match wins.

    Reading the numbers myself and then naming an outcome is how a bar gets re-fitted to its
    result. This maps cells -> outcome mechanically, so the verdict is whatever the frozen rule
    says it is.
    """
    def rel(null: str) -> str:
        c = cells[f"{arm}|{null}"]
        if not c["resolved"]:
            return "tie"
        return "gen" if c["paired_delta_nats"] > 0 else "null"

    comp, site, ctx = rel("shuffled_null_ecoli"), rel("site_independent_ecoli"), rel("markov_ecoli")
    if "null" in (comp, site):
        v = "NEW_NULL_DOMINATES_GENERATOR"
    elif comp == site == ctx == "gen":
        v = "GENERATOR_CLEARS_FULL_PANEL"
    elif "gen" in (comp, site):
        v = "GENERATOR_BEATS_ONLY_THE_WEAK_NULLS"
    else:
        v = "UNRESOLVED"
    return {
        "arm": arm,
        "VERDICT": v,
        "vs_shuffled_null": comp,
        "vs_site_independent": site,
        "vs_markov_ecoli": ctx,
        "reading": PREREGISTERED["outcomes"][v],
    }


# --------------------------------------------------------------------------- the two new nulls


def fit_site_independent(seqs, *, length=LENGTH, smooth=SMOOTH):
    """Per-position base probabilities. The POSITION-only null: no context, position fully known."""
    if not seqs:
        raise PanelError("site-independent fit needs at least one sequence")
    counts = [[smooth] * 4 for _ in range(length)]
    for s in seqs:
        if len(s) != length:
            raise PanelError(f"expected length {length}, got {len(s)}")
        for i, ch in enumerate(s):
            j = BASES.find(ch)
            if j < 0:
                raise PanelError(f"non-ACGT base {ch!r} at {i}")
            counts[i][j] += 1.0
    return [[c / sum(row) for c in row] for row in counts]


def fit_composition(seqs, *, smooth=SMOOTH):
    """Global base probabilities -- the SHUFFLED null, and the floor the other two must clear.

    This is what a site-independent model becomes once its columns are shuffled: every position
    shares one distribution, so all positional information is gone and only composition remains.
    """
    if not seqs:
        raise PanelError("composition fit needs at least one sequence")
    counts = [smooth] * 4
    for s in seqs:
        for ch in s:
            j = BASES.find(ch)
            if j < 0:
                raise PanelError(f"non-ACGT base {ch!r}")
            counts[j] += 1.0
    tot = sum(counts)
    return [c / tot for c in counts]


def site_nll(table, seq, pos):
    return -math.log(table[pos][BASES.index(seq[pos])])


def comp_nll(probs, seq, pos):
    return -math.log(probs[BASES.index(seq[pos])])


# --------------------------------------------------------------------------- paired statistics


def paired(better, worse, *, seed, boot):
    """`better` minus `worse`, bootstrapped. Positive delta = `better` has the LOWER cross-entropy.

    Mirrors glm_beyond_kmer.summarise exactly: same difference orientation, same percentile
    interval, same `resolved` rule.
    """
    if len(better) != len(worse):
        raise PanelError("paired arrays differ in length")
    d = [w - b for b, w in zip(better, worse)]
    r = random.Random(seed)
    boots = []
    for _ in range(boot):
        s = [d[r.randrange(len(d))] for _ in range(len(d))]
        boots.append(statistics.fmean(s))
    boots.sort()
    lo, hi = boots[int(0.025 * len(boots))], boots[int(0.975 * len(boots))]
    return {
        "n_windows": len(d),
        "left_ce_nats": round(statistics.fmean(better), 4),
        "right_ce_nats": round(statistics.fmean(worse), 4),
        "paired_delta_nats": round(statistics.fmean(d), 4),
        "ci95": [round(lo, 4), round(hi, 4)],
        "resolved": bool(lo > 0 or hi < 0),
        "n_windows_left_better": sum(1 for x in d if x > 0),
    }


# --------------------------------------------------------------------------- reconstruction


def rebuild_positions(n_windows, positions, seed, length=LENGTH):
    """The scorer's position draw, replicated exactly.

    The RNG is created ONCE and shared across windows, so the draws form a single dependent stream.
    Re-seeding per window would silently give every window identical positions and would still
    reconcile on window 0 alone -- which is why the gate checks every window, not a sample.
    """
    rng = random.Random(seed)
    return [sorted(rng.sample(range(ORDER, length), positions)) for _ in range(n_windows)]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--scored-dir", type=Path, default=SCORED_DIR)
    ap.add_argument("--out", type=Path, default=OUT)
    ap.add_argument("--boot", type=int, default=2000)
    a = ap.parse_args()

    corpus = load_corpus_windows()
    if not corpus:
        print("REFUSING: empty corpus", file=sys.stderr)
        return 2
    if {len(w) for w in corpus} != {LENGTH}:
        print(f"REFUSING: corpus windows are not all length {LENGTH}", file=sys.stderr)
        return 2
    cset = set(corpus)

    ecoli_fit = ecoli_fit_windows(cset, length=LENGTH)
    if len(ecoli_fit) < 500:
        print(f"REFUSING: only {len(ecoli_fit)} draft E. coli windows in corpus", file=sys.stderr)
        return 2
    print(f"corpus {len(corpus)} windows | E. coli fit set {len(ecoli_fit)}", flush=True)

    # The two committed chains, and the two NEW nulls, each fit on BOTH the pooled corpus and the
    # E. coli subset -- a 2x2 of {pooled, ecoli} x {context, position} plus the composition floor.
    chain_pooled = MarkovGenerator(order=ORDER, seed=0).fit(corpus)
    chain_ecoli = MarkovGenerator(order=ORDER, seed=0).fit(ecoli_fit)
    si_pooled = fit_site_independent(corpus)
    si_ecoli = fit_site_independent(ecoli_fit)
    comp_pooled = fit_composition(corpus)
    comp_ecoli = fit_composition(ecoli_fit)
    print("composition (ecoli fit): "
          + ", ".join(f"{b} {p:.4f}" for b, p in zip(BASES, comp_ecoli)), flush=True)

    sets = mg1655_candidate_sets(cset, length=LENGTH)
    of = sets["orthology_free"]

    arts = sorted(p for p in a.scored_dir.glob("*.json") if p.name != "_kernel_manifest.json")
    if not arts:
        print(f"REFUSING: no scored artifacts in {a.scored_dir}", file=sys.stderr)
        return 2
    loaded = []
    needs_a = False
    for p in arts:
        d = json.loads(p.read_text(encoding="utf-8"))
        if not d.get("per_window") or not d.get("run_args"):
            print(f"REFUSING: {p.name} lacks per_window/run_args", file=sys.stderr)
            return 2
        loaded.append((p, d))
        if d["run_args"]["eval_set"] != "orthology_free":
            needs_a = True

    # Set A is only needed by the two non-default sets and costs ten natural-set builds.
    withheld = set()
    if needs_a:
        from dna_decode.glm.training_corpus import falsifier_eval_windows
        print("building set A (ten natural-set builds)...", flush=True)
        withheld = falsifier_eval_windows(MG1655_DIR, length=LENGTH, limit=SCORING_LIMIT,
                                          seeds=range(N_SEEDS), kind="promoter")
        print(f"set A: {len(withheld)} windows", flush=True)

    results = {}
    recon = {}

    for p, d in loaded:
        ra = d["run_args"]
        rows = d["per_window"]
        try:
            seqs, prov = select_eval_windows(
                ra["eval_set"], of_seqs=of, withheld_seqs=withheld,
                n_windows=ra["n_windows"], selection_seed=ra["selection_seed"])
        except EvalSetUnavailable as exc:
            print(f"REFUSING: {p.name}: {exc}", file=sys.stderr)
            return 2
        if len(seqs) != len(rows):
            print(f"REFUSING: {p.name}: rebuilt {len(seqs)} windows vs {len(rows)} persisted",
                  file=sys.stderr)
            return 2

        posns = rebuild_positions(len(seqs), ra["positions"], ra["seed"])

        # ---- RECONCILE GATE: every window, every committed column, before any new number ----
        bad = []
        for wi, (seq, pos, row) in enumerate(zip(seqs, posns, rows)):
            got = {
                "markov_pooled": statistics.fmean(markov_nll(chain_pooled, seq, q) for q in pos),
                "markov_ecoli": statistics.fmean(markov_nll(chain_ecoli, seq, q) for q in pos),
                "uniform": math.log(4),
            }
            for k, v in got.items():
                if k in row and round(v, RECONCILE_DECIMALS) != round(row[k], RECONCILE_DECIMALS):
                    bad.append({"window": wi, "column": k,
                                "recomputed": round(v, RECONCILE_DECIMALS),
                                "persisted": round(row[k], RECONCILE_DECIMALS)})
        recon[p.name] = {"n_windows": len(rows), "n_mismatches": len(bad),
                         "first_mismatches": bad[:5], "pass": not bad,
                         "decimals": RECONCILE_DECIMALS}
        if bad:
            print(f"RECONCILE FAILED on {p.name}: {len(bad)} mismatches; first {bad[0]}",
                  file=sys.stderr)
            print("REFUSING: the reconstruction is not the scored one, so no new column is "
                  "comparable to the persisted generator numbers.", file=sys.stderr)
            return 2
        print(f"  reconcile OK {p.name} ({len(rows)} windows)", flush=True)

        # ---- the NEW columns, on the SAME positions ----
        cols = {k: [] for k in ("site_independent_ecoli", "site_independent_pooled",
                                "shuffled_null_ecoli", "shuffled_null_pooled")}
        for seq, pos in zip(seqs, posns):
            cols["site_independent_ecoli"].append(
                statistics.fmean(site_nll(si_ecoli, seq, q) for q in pos))
            cols["site_independent_pooled"].append(
                statistics.fmean(site_nll(si_pooled, seq, q) for q in pos))
            cols["shuffled_null_ecoli"].append(
                statistics.fmean(comp_nll(comp_ecoli, seq, q) for q in pos))
            cols["shuffled_null_pooled"].append(
                statistics.fmean(comp_nll(comp_pooled, seq, q) for q in pos))

        for k in ("markov_ecoli", "markov_pooled", "uniform"):
            cols[k] = [r[k] for r in rows]
        for arm in ("left_only", "bidirectional"):
            cols[arm] = [r[arm] for r in rows]

        cells = {}
        nulls = ("markov_ecoli", "site_independent_ecoli", "shuffled_null_ecoli",
                 "markov_pooled", "site_independent_pooled", "shuffled_null_pooled", "uniform")
        for arm in ("left_only", "bidirectional"):
            for null in nulls:
                cells[f"{arm}|{null}"] = paired(cols[arm], cols[null],
                                                seed=ra["seed"] + 1, boot=a.boot)
        # The null-vs-null contrasts are where the decomposition actually lives.
        for left, right in (("site_independent_ecoli", "shuffled_null_ecoli"),
                            ("markov_ecoli", "shuffled_null_ecoli"),
                            ("markov_ecoli", "site_independent_ecoli")):
            cells[f"{left}|{right}"] = paired(cols[left], cols[right],
                                              seed=ra["seed"] + 1, boot=a.boot)

        results[p.name] = {
            "eval_set": ra["eval_set"],
            "checkpoint_step": d.get("checkpoint_step"),
            "ckpt": ra.get("ckpt"),
            "provenance": prov,
            "mean_ce_nats": {k: round(statistics.fmean(v), 4) for k, v in cols.items()},
            # The nats each null captures below uniform -- the decomposition, stated directly so a
            # reader does not have to subtract 1.3863 by hand to see where the signal lives.
            "captures_below_uniform_nats": {
                k: round(math.log(4) - statistics.fmean(cols[k]), 4)
                for k in ("shuffled_null_ecoli", "site_independent_ecoli", "markov_ecoli",
                          "left_only", "bidirectional")},
            "verdicts": {arm: compute_verdict(cells, arm)
                         for arm in ("left_only", "bidirectional")},
            "cells": cells,
            "per_window": [{k: round(cols[k][i], 6) for k in sorted(cols)}
                           for i in range(len(rows))],
        }
        print(f"    {ra['eval_set']}: site-indep(ecoli) "
              f"{statistics.fmean(cols['site_independent_ecoli']):.4f} | "
              f"shuffled {statistics.fmean(cols['shuffled_null_ecoli']):.4f} | "
              f"markov_ecoli {statistics.fmean(cols['markov_ecoli']):.4f} | "
              f"left_only {statistics.fmean(cols['left_only']):.4f}", flush=True)

    out = {
        "record": "glm-baseline-panel-ce-v1",
        "date": "2026-11-28",
        "question": PREREGISTERED["question"],
        "PREREGISTERED": PREREGISTERED,
        "mandated_panel": MANDATED_PANEL,
        "torch_free": True,
        "generator_ce_source": ("the persisted `per_window` block of the six artifacts in "
                                "wiki/glm_repretrain_scored_2026-11-25/ -- no checkpoint loaded"),
        "reconcile_gate": recon,
        # The headline is restricted to the CRITERION's own eval set, and to the PRIMARY arm, so a
        # more favourable set or the information-asymmetric arm cannot stand in for it.
        "headline": {
            "eval_set": "withheld_orthology_free",
            "primary_arm": "left_only",
            "why_primary": ("left_only matches the information budget -- both the generator and the "
                            "nulls see only the unmasked left prefix. bidirectional lets the "
                            "generator see both flanks while every null does not, so it is the "
                            "more favourable framing and is reported second."),
            "verdicts": {n: v["verdicts"]["left_only"]["VERDICT"]
                         for n, v in results.items()
                         if v["eval_set"] == "withheld_orthology_free"},
            "secondary_arm_verdicts": {n: v["verdicts"]["bidirectional"]["VERDICT"]
                                       for n, v in results.items()
                                       if v["eval_set"] == "withheld_orthology_free"},
        },
        "arms": results,
        "does_NOT_write_the_mission_gate": {
            "path": "wiki/glm_generator_baseline_panel.json",
            "why": ("that path is the live `file-exists` predicate of the budget-4 mission gap, the "
                    "self-firing expiry argument of `north-star-step4-animals-deferred`, and one of "
                    "five non-vacuity fixtures in tests/test_mission_gap_predicates.py. Writing it "
                    "would close a mission gap, re-open the north star's step 4 and flip a test "
                    "fixture. All three are authority calls."),
            "proposal": ("if the 2-of-5 definability split above is ratified, the gate can be met "
                         "by this artifact plus a random-proposal null; until then the gap stands "
                         "open."),
        },
        "honest_limits": [
            "Infilling cross-entropy is NOT generation quality; necessary and diagnostic only.",
            "Two mandated baselines are undefined on this readout and are recorded as undefined.",
            "`shuffled_null` is the FIT-SIDE reading (composition); the generator-on-shuffled-eval "
            "variant needs torch and was not run.",
            "The scored checkpoints are the v5 arms (step 6000), not the committed checkpoint, so "
            "these are the ARMS' generator numbers -- the same ones the arms ledger records.",
            "Both new nulls use add-one smoothing, matching the chain, so a difference between them "
            "is not a smoothing artifact.",
        ],
    }
    a.out.write_text(json.dumps(out, indent=1) + "\n", encoding="utf-8")
    print(f"\nwrote {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
