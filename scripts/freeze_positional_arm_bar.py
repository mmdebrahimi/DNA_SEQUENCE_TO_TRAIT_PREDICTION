"""Freeze the positional arm's acceptance bar BEFORE any score exists. Refuses on anchor drift.

A bar frozen after a run is not a bar. This script writes
`wiki/positional_arm_acceptance_bar.json` from two sources and nothing else: the Step 2 feasibility
artifact (for figures that had to be MEASURED) and the committed comparator plus the arms ledger
(for figures that are INHERITED). No threshold here is re-chosen for this arm.

**`verify()` re-reads its anchors and REFUSES to write if they have drifted.** The atlas median, the
GC median and the arm's `planned` status are all read live, because a bar anchored to retyped
numbers is anchored to nothing -- and because the arm moving to `measured` before its bar is frozen
would mean a number already exists.

**What is INHERITED, not chosen here:**

- margin **0.02**, verbatim from the tile-rung arms, so the bars stay commensurable. It is the same
  margin the 499bp fragment arm, the hybrid/stack arms and the measured atlas arm were graded
  against, and the inheritance is the reason it is defensible.
- the comparator: the committed atlas per-seed values in
  `data/promoteratlas/local_reconcile_2026-11-01.json` (`verdict: RECONCILED`, `mismatched: []`),
  so the incumbent is never re-run and can never be retro-edited by this arm.
- the split: `leave_peak_out(frac=0.25)` over seeds 0-9, with the committed `n_test_per_seed`.
- the scorer: `atlas_embed.ridge_spearman`, unchanged.
- the shuffled-label null and its bar, in the same FORM as the fragment-pyramid arm: on failure the
  response is MORE SEEDS, never a weaker control.

**What is MEASURED and inherited from the probe rather than guessed:** the parity tolerance (derived
from the achieved deviation by a stated rule), the autocorrelation lag set, the per-rung RAM
verdicts and the throughput.

**Multiplicity is corrected in the bar, not afterwards.** Six rungs against one incumbent is a
garden of forking paths if the best is reported. Two CO-PRIMARIES test DIFFERENT hypotheses --
`bins10_mean_oriented` asks whether a position-resolved reduction beats a position-averaged one,
`autocorr` asks whether periodic SPACING carries signal -- and each is judged at a Bonferroni
two-comparison level. The four secondaries are judged at the level for every rung measured and are
explicitly INELIGIBLE to be reported as the arm's result.
"""

from __future__ import annotations

import argparse
import json
import math
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
OUT = REPO / "wiki/positional_arm_acceptance_bar.json"
FEASIBILITY = REPO / "wiki/positional_arm_feasibility_2026-11-01.json"
RECONCILE = REPO / "data/promoteratlas/local_reconcile_2026-11-01.json"
LEDGER = REPO / "soraya_runs/arms-ledger.json"

CRITERION_PREFIX = "ceiling:"
ARM = "positional-non-mean-pooled-promoteratlas-150bp"
INCUMBENT_ARM = "promoteratlas-base-embeddings-frozen-150bp-padded-to-200"

#: INHERITED verbatim from the tile-rung arms. Never re-chosen here.
MARGIN = 0.02

#: The anchors `verify()` refuses on. Values are what the committed artifacts said when this bar
#: was designed; a drift means the bar would be describing a different rung.
EXPECTED_ATLAS_MEDIAN = 0.5076
EXPECTED_GC_MEDIAN = 0.3112
EXPECTED_CONV = 0.3234
EXPECTED_N_SEEDS = 10
EXPECTED_N_TILES = 44_106

FAMILY_ALPHA = 0.05

CO_PRIMARIES = ("bins10_mean_oriented", "autocorr")
SECONDARIES = ("bins10_mean_asgiven", "bins10_max_oriented", "bins2_mean_oriented",
               "bins5_mean_oriented")

#: The CLOSED outcome set. The reachable set must EQUAL this, which is what catches a dead branch.
OUTCOMES = (
    "POSITIONAL_ADDS",
    "POSITIONAL_NO_GAIN",
    "POSITIONAL_AMBIGUOUS",
    "INDETERMINATE_PARITY_FAIL",
    "INDETERMINATE_SPLIT_MISMATCH",
    "INDETERMINATE_NULL_FAIL",
)

#: Validity FIRST, so a broken run can never read as a finding. Split equality is checked before
#: everything, because every other number is meaningless if the test sets differ.
CHECK_ORDER = (
    ("split_equality", "INDETERMINATE_SPLIT_MISMATCH"),
    ("parity", "INDETERMINATE_PARITY_FAIL"),
    ("shuffled_label_null", "INDETERMINATE_NULL_FAIL"),
    ("co_primary_decision", ("POSITIONAL_ADDS", "POSITIONAL_NO_GAIN", "POSITIONAL_AMBIGUOUS")),
)


class BarRefused(RuntimeError):
    """An anchor drifted, or the arm is no longer in a state where a bar may be frozen."""


# ---------------------------------------------------------------------------------------------------
# PURE helpers
# ---------------------------------------------------------------------------------------------------
def bonferroni_percentiles(n_comparisons: int, *, family_alpha: float = FAMILY_ALPHA):
    """Bootstrap percentiles for a two-sided interval at the Bonferroni-corrected level.

    `family_alpha / n` per comparison, so the interval is `1 - alpha_each` wide. Returned as
    percentiles because that is what `np.percentile` consumes in the scoring runner.
    """
    n = int(n_comparisons)
    if n < 1:
        raise ValueError("n_comparisons must be >= 1, got %r" % (n_comparisons,))
    if not 0.0 < float(family_alpha) < 1.0:
        raise ValueError("family_alpha must be in (0, 1), got %r" % (family_alpha,))
    alpha_each = float(family_alpha) / n
    return (100.0 * alpha_each / 2.0, 100.0 * (1.0 - alpha_each / 2.0))


def derive_parity_tolerance(measured_deviation: float) -> float:
    """Smallest power of ten at least 10x the ACHIEVED deviation.

    Stated as a rule so the tolerance is derived from the measurement rather than preferred: the
    probe measured 4.352e-07 (float32 round-off), 10x is 4.352e-06, and the smallest power of ten
    at or above that is 1e-05. Tight enough that a real forward-pass drift cannot hide inside it,
    loose enough that float32 summation order cannot trip it.
    """
    m = float(measured_deviation)
    if not math.isfinite(m) or m <= 0:
        raise ValueError("measured_deviation must be finite and > 0, got %r" % (measured_deviation,))
    return float(10.0 ** math.ceil(math.log10(m * 10.0)))


def reachable_outcomes(check_order=CHECK_ORDER) -> set[str]:
    """Every outcome the check order can actually return.

    Compared for EQUALITY against the closed set, not containment: containment would pass a bar
    with an unreachable outcome, which is a dead validity branch nobody would notice.
    """
    out: set[str] = set()
    for _gate, result in check_order:
        if isinstance(result, str):
            out.add(result)
        else:
            out.update(result)
    return out


def classify(co_primary_cis: dict, *, margin: float = MARGIN) -> tuple[str, str]:
    """The arm's verdict from the co-primary intervals alone. Pure, ordered, first match wins.

    Each value is `(ci_lo, ci_hi)` for the PAIRED gain of that rung over the committed atlas
    per-seed values, at the two-comparison corrected level. Secondaries are not consulted: they are
    ineligible as the arm's result by construction, not by convention.
    """
    missing = [k for k in CO_PRIMARIES if k not in co_primary_cis]
    if missing:
        raise ValueError(
            "classify needs an interval for every co-primary; missing %r. Returning a verdict on a "
            "subset would silently drop a declared hypothesis." % missing
        )
    adds = [k for k in CO_PRIMARIES if co_primary_cis[k][0] > margin]
    if adds:
        return ("POSITIONAL_ADDS",
                "the corrected interval for %s lies entirely above the inherited margin %+.4f"
                % (" and ".join(adds), margin))
    if all(co_primary_cis[k][1] < margin for k in CO_PRIMARIES):
        return ("POSITIONAL_NO_GAIN",
                "both co-primary intervals lie entirely BELOW the margin %+.4f" % margin)
    return ("POSITIONAL_AMBIGUOUS",
            "at least one co-primary interval contains the margin %+.4f, so no difference is "
            "detectable at the corrected level" % margin)


# ---------------------------------------------------------------------------------------------------
# verify -- re-read the anchors, refuse on drift
# ---------------------------------------------------------------------------------------------------
def _ceiling_arms(ledger: dict) -> dict:
    crits = [k for k in ledger.get("criteria", {}) if k.startswith(CRITERION_PREFIX)]
    if len(crits) != 1:
        raise BarRefused(
            "expected exactly one '%s' criterion in the ledger, found %d: %r"
            % (CRITERION_PREFIX, len(crits), crits)
        )
    return ledger["criteria"][crits[0]], crits[0]


def verify(*, reconcile_path: Path = RECONCILE, ledger_path: Path = LEDGER,
           feasibility_path: Path = FEASIBILITY) -> dict:
    """Re-read every anchor live. Raises `BarRefused` on any drift. Returns what it verified."""
    for p in (reconcile_path, ledger_path, feasibility_path):
        if not Path(p).is_file():
            raise BarRefused("anchor source absent: %s" % p)

    rec = json.loads(Path(reconcile_path).read_text(encoding="utf-8"))
    if rec.get("verdict") != "RECONCILED":
        raise BarRefused(
            "the comparator's verdict is %r, not RECONCILED. A bar anchored to an unreconciled "
            "comparator would grade this arm against numbers nobody has checked."
            % rec.get("verdict")
        )
    if rec.get("mismatched"):
        raise BarRefused("the comparator reports mismatched seeds %r" % rec["mismatched"])
    if rec.get("atlas_median") != EXPECTED_ATLAS_MEDIAN:
        raise BarRefused(
            "atlas median drifted: comparator says %r, this bar was designed against %r"
            % (rec.get("atlas_median"), EXPECTED_ATLAS_MEDIAN)
        )
    if rec.get("gc_median") != EXPECTED_GC_MEDIAN:
        raise BarRefused(
            "GC median drifted: comparator says %r, this bar was designed against %r"
            % (rec.get("gc_median"), EXPECTED_GC_MEDIAN)
        )
    if int(rec.get("n_seeds", -1)) != EXPECTED_N_SEEDS:
        raise BarRefused("comparator has %r seeds, expected %d" % (rec.get("n_seeds"),
                                                                   EXPECTED_N_SEEDS))
    if int(rec.get("n_tiles", -1)) != EXPECTED_N_TILES:
        raise BarRefused("comparator has %r tiles, expected %d" % (rec.get("n_tiles"),
                                                                   EXPECTED_N_TILES))
    for key in ("atlas_per_seed", "gc_per_seed", "n_test_per_seed"):
        if len(rec.get(key, {})) != EXPECTED_N_SEEDS:
            raise BarRefused("comparator's %s has %d entries, expected %d"
                             % (key, len(rec.get(key, {})), EXPECTED_N_SEEDS))

    ledger = json.loads(Path(ledger_path).read_text(encoding="utf-8"))
    arms, criterion = _ceiling_arms(ledger)
    if ARM not in arms:
        raise BarRefused(
            "the arm %r is NOT declared on %r. A bar for an undeclared arm would grade work that "
            "the ledger does not know exists." % (ARM, criterion)
        )
    status = arms[ARM].get("status")
    if status != "planned":
        raise BarRefused(
            "the arm's status is %r, not 'planned'. A bar may only be frozen while no number "
            "exists; %s means one already does, and freezing now would be freezing after the fact."
            % (status, status)
        )
    if arms[ARM].get("metric") is not None:
        raise BarRefused("the arm already carries metric %r" % arms[ARM].get("metric"))
    if INCUMBENT_ARM not in arms or arms[INCUMBENT_ARM].get("metric") != EXPECTED_ATLAS_MEDIAN:
        raise BarRefused(
            "the incumbent %r reads %r in the ledger, expected %r"
            % (INCUMBENT_ARM, arms.get(INCUMBENT_ARM, {}).get("metric"), EXPECTED_ATLAS_MEDIAN)
        )
    if arms.get("conv-3width-13k-params-150bp", {}).get("metric") != EXPECTED_CONV:
        raise BarRefused(
            "conv anchor reads %r, expected %r"
            % (arms.get("conv-3width-13k-params-150bp", {}).get("metric"), EXPECTED_CONV)
        )
    if arms.get("gc-single-feature", {}).get("metric") != EXPECTED_GC_MEDIAN:
        raise BarRefused(
            "gc anchor reads %r, expected %r"
            % (arms.get("gc-single-feature", {}).get("metric"), EXPECTED_GC_MEDIAN)
        )

    feas = json.loads(Path(feasibility_path).read_text(encoding="utf-8"))
    if feas.get("verdict") != "FEASIBLE" or not feas["parity"]["passed"]:
        raise BarRefused(
            "feasibility verdict is %r with parity passed=%r; a bar frozen on a failed probe would "
            "grade a run whose representation is not the incumbent's"
            % (feas.get("verdict"), feas["parity"]["passed"])
        )

    extra = reachable_outcomes() ^ set(OUTCOMES)
    if extra:
        raise BarRefused(
            "the reachable outcome set and the closed set differ by %r. Containment is not enough: "
            "an unreachable outcome is a dead validity branch and an unlisted one is an "
            "ungradeable result." % sorted(extra)
        )

    return {"reconcile": rec, "arms": arms, "criterion": criterion, "feasibility": feas}


# ---------------------------------------------------------------------------------------------------
# build
# ---------------------------------------------------------------------------------------------------
def build_bar(verified: dict) -> dict:
    rec, feas, criterion = verified["reconcile"], verified["feasibility"], verified["criterion"]
    measured_dev = float(feas["parity"]["max_abs_deviation"])
    parity_tol = derive_parity_tolerance(measured_dev)
    lo2, hi2 = bonferroni_percentiles(len(CO_PRIMARIES))
    n_all = len(CO_PRIMARIES) + len(SECONDARIES)
    loN, hiN = bonferroni_percentiles(n_all)

    return {
        "schema": "positional-arm-acceptance-bar-v1",
        "frozen": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "frozen_before_any_score_exists": True,
        "arm": ARM,
        "criterion": criterion,
        "question": (
            "Does a POSITION-RESOLVED reduction of the SAME frozen PromoterAtlas latent beat the "
            "position-AVERAGED mean-pool incumbent (median Spearman 0.5076) by the rung's inherited "
            "margin of 0.02?"
        ),
        "why_this_is_not_a_re_read_of_the_measured_arm": (
            "The incumbent arm reduces the latent with z[:, :, :L].mean(dim=2). This arm reduces the "
            "SAME forward pass differently and never re-reads the incumbent's number, which is read "
            "from its committed per-seed table. Under R6 it is a distinct arm, declared before this "
            "bar existed."
        ),
        "scoping_corrections_fixed_in_advance": [
            "This arm CANNOT test the -35/-10 spacing claim its own ledger rationale asserts. "
            "genomewide.Tile exposes seq, expression, peak, strand, category, active, rep1, rep2 -- "
            "there is no TSS, no promoter register and no peak-relative offset, so a bin index is a "
            "position in an ARBITRARY 150 bp window, not a promoter coordinate. What is testable is "
            "narrower: does a position-resolved reduction beat a position-averaged one.",
            "Because of that, the readout matching the stated spacing mechanism is NOT bins. Spacing "
            "is translation-invariant, so the anchor-free autocorrelation is a declared CO-PRIMARY "
            "testing a DIFFERENT hypothesis -- not a second shot at the same one.",
        ],
        "inherited_not_chosen": {
            "margin": {
                "value": MARGIN,
                "derivation": ("INHERITED verbatim from the tile-rung arms so the bars stay "
                               "commensurable; the same margin the 499bp fragment arm, the "
                               "hybrid/stack arms and the measured atlas arm were graded against."),
            },
            "comparator": {
                "source": "data/promoteratlas/local_reconcile_2026-11-01.json",
                "verdict": rec["verdict"],
                "atlas_median": rec["atlas_median"],
                "atlas_per_seed": rec["atlas_per_seed"],
                "gc_median": rec["gc_median"],
                "gc_per_seed": rec["gc_per_seed"],
                "must_be_recomputed_in_run": False,
                "why": ("the incumbent is never re-run by this arm, so this arm cannot move the "
                        "number it is being graded against"),
            },
            "split": {
                "spec": "genomewide.leave_peak_out(tiles, frac=0.25, seed=s) for s in 0..9",
                "frac": rec["frac"], "n_seeds": rec["n_seeds"], "seed_range": "0..9",
                "n_test_per_seed": rec["n_test_per_seed"],
                "leakage_unit": "the called expression PEAK, not the tile",
            },
            "scorer": {
                "spec": "dna_decode.glm.atlas_embed.ridge_spearman, unchanged",
                "protocol": ("StandardScaler on train, Ridge(alpha=1.0), fit on log1p(expression), "
                             "score Spearman against RAW expression"),
            },
        },
        "measured_and_inherited_from_the_probe": {
            "source": "wiki/positional_arm_feasibility_2026-11-01.json",
            "parity_max_abs_deviation_achieved": measured_dev,
            "parity_tolerance": parity_tol,
            "parity_tolerance_rule": ("the smallest power of ten at least 10x the achieved "
                                      "deviation; derived from the measurement, not preferred"),
            "lag_set": feas["lag_set"]["lags"],
            "lag_set_name": feas["lag_set"]["chosen"],
            "lag_set_rule": feas["lag_set"]["rule"],
            "checkpoint_sha256": feas["checkpoint"]["sha256"],
            "throughput_seq_per_second": feas["throughput"]["seq_per_second"],
            "free_ram_bytes_at_freeze": feas["resources"]["free_ram_bytes"],
            "ram_headroom_factor": feas["resources"]["headroom_factor"],
            "rungs_fitting_in_ram": {r["rung"]: r["fits_in_ram"] for r in feas["resources"]["rungs"]},
        },
        "multiplicity": {
            "family_alpha": FAMILY_ALPHA,
            "co_primaries": list(CO_PRIMARIES),
            "co_primary_hypotheses": {
                "bins10_mean_oriented": ("does a POSITION-RESOLVED reduction beat a "
                                         "position-AVERAGED one"),
                "autocorr": ("does periodic SPACING carry signal, measured anchor-free so it does "
                             "not need a promoter register"),
            },
            "co_primary_n_comparisons": len(CO_PRIMARIES),
            "co_primary_alpha_each": FAMILY_ALPHA / len(CO_PRIMARIES),
            "co_primary_ci_percentiles": [lo2, hi2],
            "secondaries": list(SECONDARIES),
            "secondary_n_comparisons": n_all,
            "secondary_alpha_each": FAMILY_ALPHA / n_all,
            "secondary_ci_percentiles": [loN, hiN],
            "secondaries_are_ineligible_as_the_arms_result": True,
            "why_ineligible": ("six rungs against one incumbent is a garden of forking paths if the "
                               "best is reported; classify() consults the co-primaries ONLY, so "
                               "ineligibility is enforced by construction rather than by convention"),
            "divisor_discrepancy_recorded": (
                "The plan's prose named a '5-comparison corrected level' for the secondaries while "
                "listing FOUR of them. Neither 5 nor 4 is the family-wise count over everything "
                "measured, so this bar uses %d = every rung on the ladder. That is strictly more "
                "conservative than the plan's 5, it is internally consistent, and it is "
                "outcome-neutral: a secondary can never be the arm's result, so its interval width "
                "moves no verdict. Recorded rather than silently applied." % n_all
            ),
            "autocorr_needs_no_orientation_variant": (
                "an autocorrelation is invariant under reversal of its series -- the set of position "
                "pairs separated by k is unchanged when the series is read backwards -- so the "
                "autocorr co-primary is IMMUNE to the unresolved orientation convention that the "
                "bins rungs carry a 1-bit inference about"
            ),
        },
        "validity_controls": {
            "split_equality": {
                "spec": ("recompute leave_peak_out(frac=0.25) for seeds 0-9 and assert "
                         "n_test_per_seed matches the committed values for ALL ten"),
                "on_failure": "INDETERMINATE_SPLIT_MISMATCH, and stop",
                "why_first": ("every other number is meaningless if the test sets differ; the arm "
                              "would not be scored on the incumbent's held-out peaks"),
            },
            "parity": {
                "spec": ("derive the global mean from the STORED 10-bin means, score it through "
                         "ridge_spearman on those same splits, and assert it reproduces the "
                         "committed atlas per-seed values"),
                "tolerance": parity_tol,
                "on_failure": "INDETERMINATE_PARITY_FAIL, and stop",
                "why": ("this is what makes every later number comparable to 0.5076 rather than "
                        "merely adjacent to it"),
            },
            "shuffled_label_null": {
                "spec": "refit the ridge on shuffled tile labels, same splits and seeds",
                "bar": "|median gain over GC| < %.2f" % MARGIN,
                "derivation": ("same FORM and bar as the fragment-pyramid arm, whose own 10-seed "
                               "null landed at +0.0517 and voided its verdict until a 270-seed "
                               "pooled null cleared the gate"),
                "on_failure": "INDETERMINATE_NULL_FAIL, and SPEND MORE SEEDS -- never re-specify",
            },
        },
        "outcomes_defined_before_running": {
            "POSITIONAL_ADDS": ("the corrected interval for at least one CO-PRIMARY's paired gain "
                                "over the committed atlas per-seed values lies entirely above "
                                "+%.2f" % MARGIN),
            "POSITIONAL_NO_GAIN": ("BOTH co-primary intervals lie entirely below +%.2f" % MARGIN),
            "POSITIONAL_AMBIGUOUS": ("at least one co-primary interval contains +%.2f, so no "
                                     "difference is detectable at the corrected level" % MARGIN),
            "INDETERMINATE_SPLIT_MISMATCH": "the recomputed splits differ from the committed ones",
            "INDETERMINATE_PARITY_FAIL": ("the global mean derived from the stored bin means does "
                                          "not reproduce the committed embedding"),
            "INDETERMINATE_NULL_FAIL": "the shuffled-label control cleared the bar; nothing is graded",
        },
        "verdict_check_order": [g for g, _ in CHECK_ORDER],
        "check_order_rationale": (
            "Validity FIRST, so a broken run can never read as a finding, and split equality before "
            "everything else because every other number is meaningless if the test sets differ."
        ),
        "closed_outcome_set": list(OUTCOMES),
        "reachable_outcome_set": sorted(reachable_outcomes()),
        "reachable_equals_closed": sorted(reachable_outcomes()) == sorted(OUTCOMES),
        "pre_commitments": [
            "The margin 0.02 is INHERITED and will NOT be re-specified after a number exists.",
            "If the shuffled-label null fails, the response is MORE SEEDS, never a weaker control.",
            "The measured sibling's 0.5076 is NEVER retro-edited; it is read from its committed "
            "per-seed table and this arm cannot move it.",
            "A secondary rung is NEVER reported as the arm's result, however it scores. classify() "
            "cannot see the secondaries at all.",
            "The verdict field emitted by the run is whatever this bar mechanically returns. A "
            "withdrawn reading goes in a namespace-separate block, leaving the verdict unedited.",
            "A miss does NOT kill the fine-tuned variant or the four shipped expression heads. Each "
            "is a NEW arm to declare under R6, never a re-read of this one.",
            "A rung whose design matrix fails the RAM pre-check is recorded as SKIPPED with its "
            "byte estimate, never silently dropped and never quietly retried at a smaller width.",
        ],
        "honest_limits_fixed_in_advance": [
            "No promoter anchor exists on this substrate, so a NO_GAIN says that position bins "
            "relative to an ARBITRARY window add nothing -- NOT that positional structure is absent. "
            "The negative must never be over-read as closing positional readouts generally.",
            "The orientation convention of `seq` is UNRESOLVED and is inferred from a performance "
            "difference between two rungs. That is a 1-bit inference, not established biology: only "
            "576 of 22,285 minus sequences (2.6%) have their reverse complement among the plus "
            "sequences, and ZERO of those share a peak, so plus and minus tiles are distinct "
            "genomic windows rather than two orientations of one element.",
            "689 duplicate sequences are known on this substrate and are UNCHANGED by this arm. A "
            "train/test-straddling duplicate affects this arm exactly as it affects the incumbent, "
            "so the PAIRED comparison is the right shape while the absolute levels inherit the issue.",
            "The encoder is FROZEN and its hyperparameters were chosen elsewhere; this is not a "
            "tuned positional model and a NO_GAIN does not bound what a fine-tuned one could do.",
            "Ten seeds is the inherited pairing. The fragment-pyramid arm needed 270 pooled seeds "
            "before its null cleared, so ten may be too few here for the same reason.",
            "The rung ceiling is 0.941 and replicate-derived. No tile-rung number may be read as a "
            "fraction of it without saying so.",
        ],
        "forecast": {
            "recorded_before_any_number_exists": True,
            "predicted_outcome": "POSITIONAL_NO_GAIN",
            "p_positional_adds": 0.25,
            "p_positional_no_gain": 0.45,
            "p_positional_ambiguous": 0.30,
            "p_autocorr_is_the_one_that_clears_given_adds": 0.4,
            "rationale": (
                "Against: every sequence-statistics arm on this rung clusters 0.2691-0.3234 while "
                "the mean-pooled atlas sits at 0.5076, and a ridge over 1,280 or 2,560 correlated "
                "columns on ~33,000 training rows has far more variance to pay for than the 128-dim "
                "pool. Without a promoter register the bins are positions in an arbitrary window, "
                "which is the scoping correction above. For: mean-pooling provably discards "
                "information the latent carries, the incumbent's own recorded limit says so, and "
                "the ridge can ignore useless columns at alpha=1.0 more cheaply than it can "
                "recover discarded ones."
            ),
            "my_forecast_record_in_this_project": (
                "3 of 7. The fragment-null P1 was falsified twice, P2 mislocated its interval, and "
                "the atlas arm's own bar put p=0.55 on LOSES_ALONE and p=0.08 on what actually "
                "happened. DO NOT defer to these priors."
            ),
        },
        "reproduce": [
            "uv run python -m scripts.positional_arm_probe",
            "uv run python scripts/freeze_positional_arm_bar.py",
            "uv run python -m scripts.positional_arm_extract",
            "uv run python -m scripts.positional_arm_score",
        ],
    }


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--print-only", action="store_true", help="verify and render, but do not write")
    args = ap.parse_args(argv)

    try:
        verified = verify()
    except BarRefused as e:
        print("REFUSED TO FREEZE: %s" % e, flush=True)
        return 2

    bar = build_bar(verified)
    if args.print_only:
        print(json.dumps(bar, indent=2, sort_keys=True))
        return 0

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(bar, indent=2, sort_keys=True), encoding="utf-8")

    m = bar["multiplicity"]
    print("anchors      atlas %.4f / gc %.4f / conv %.4f -- all re-read, none drifted"
          % (EXPECTED_ATLAS_MEDIAN, EXPECTED_GC_MEDIAN, EXPECTED_CONV))
    print("arm status   planned (a bar may only be frozen while no number exists)")
    print("margin       %+.2f INHERITED from the tile-rung arms" % MARGIN)
    print("co-primaries %r at alpha/%d -> CI percentiles [%.4f, %.4f]"
          % (list(CO_PRIMARIES), m["co_primary_n_comparisons"], *m["co_primary_ci_percentiles"]))
    print("secondaries  %d rungs at alpha/%d -> CI percentiles [%.4f, %.4f] (INELIGIBLE as result)"
          % (len(SECONDARIES), m["secondary_n_comparisons"], *m["secondary_ci_percentiles"]))
    print("parity tol   %.0e (derived from the achieved %.3e by the stated rule)"
          % (bar["measured_and_inherited_from_the_probe"]["parity_tolerance"],
             bar["measured_and_inherited_from_the_probe"]["parity_max_abs_deviation_achieved"]))
    print("lag set      %s (%d lags)" % (bar["measured_and_inherited_from_the_probe"]["lag_set_name"],
                                         len(bar["measured_and_inherited_from_the_probe"]["lag_set"])))
    print("outcomes     %d closed, reachable set EQUALS closed set: %s"
          % (len(OUTCOMES), bar["reachable_equals_closed"]))
    print("forecast     %s (p=%.2f) -- record 3 of 7, do not defer"
          % (bar["forecast"]["predicted_outcome"], bar["forecast"]["p_positional_no_gain"]))
    print("-> %s" % OUT)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
