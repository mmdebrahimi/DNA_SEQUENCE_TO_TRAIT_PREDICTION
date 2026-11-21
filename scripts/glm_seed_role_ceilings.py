"""Measure the per-seed NOISE FLOOR the committed paired grid never recorded.

    uv run python scripts/glm_seed_role_ceilings.py            # read-only; writes one artifact

**Why this probe exists.** Every distinguishability number in the generator arm is only meaningful against
the discriminator's own noise floor -- the AUROC it reaches separating natural sequence from MORE natural
sequence. `wiki/glm_generator_paired_sweep_2026-11-20.json` records, per cell,
`{model, epoch, seed, both, kmer, null{kmer,both}}` and **no ceiling field at all**. So the grid can say
the generator loses to the Markov null, and cannot say whether the null itself is above the floor. The
per-seed margin has therefore never existed, and one seed's margin was read as if it were typical.

**What a margin means here.** `margin = null - ceiling`. A POSITIVE margin is a winnable region: the null
is distinguishable above the noise floor, so a generator that beat the null would be making a real claim. A
NEGATIVE margin means the null sits BELOW its own noise floor -- in that mode there is nothing to win,
because the target a generator is being asked to beat is itself inside the measurement noise. **Nothing
clamps at zero.** A negative margin is the finding, not an error.

**The refusal.** This script will not write unless it reproduces the one committed ceiling pair at seed 7
(0.5447 kmer / 0.5865 both, from `wiki/glm_generator_falsifier_REAL_2026-10-07.json`). A probe that cannot
reproduce the single ceiling the repo already has, has no standing to report nine new ones.

Reads the natural corpus through `dna_decode.glm.corpus.build_natural_set` -- the same definition the
runner scores on, not a second description of it. No network, no GPU, no torch, no D: volume.
"""

from __future__ import annotations

import json
import statistics
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
REFSEQ = ROOT / "data" / "cache" / "refseq" / "GCF_000005845.2"
GRID = ROOT / "wiki" / "glm_generator_paired_sweep_2026-11-20.json"
OUT = ROOT / "wiki" / "glm_seed_role_ceilings_2026-11-21.json"

# The frozen geometry of the committed grid. Changing any of these makes the join invalid.
LENGTH = 150
LIMIT = 1200
KIND = "promoter"
KMER_K = 4
POSITIONAL_BINS = 6
MODES = ("kmer", "both")

# The committed seed-7 reproduction gate.
COMMITTED = {"kmer": 0.5447, "both": 0.5865}
ANCHOR_SEED = 7
GATE_TOL = 5e-4          # the source artifact rounds to 4 dp


def margin(null: float, ceiling: float) -> float:
    """The winnable region at one seed in one mode. MAY BE NEGATIVE -- deliberately unclamped.

    Kept as a named pure function because clamping it is the one edit that would silently convert the
    finding into its opposite: a `max(0, ...)` here turns "the null is below its own noise floor, so this
    mode has no winnable region" into "the margin is zero", which reads as a near-miss.
    """
    return null - ceiling


def seeds_and_nulls_from_grid(grid: dict) -> dict[int, dict[str, float]]:
    """Read the seed list and per-seed nulls FROM the committed grid, never from a literal.

    Restating the ten seeds here would let a grid built on different seeds be silently scored against the
    wrong ones. The nulls are also verified to be seed-determined rather than model-determined: all four
    models in a seed share one null pair (same corpus seed, same split, same Markov fit), so collapsing to
    one null per seed is exact and not an average. If that ever stops holding, this refuses.
    """
    per: dict[int, set] = {}
    for c in grid["grid"]:
        per.setdefault(c["seed"], set()).add((c["null"]["kmer"], c["null"]["both"]))
    out: dict[int, dict[str, float]] = {}
    for seed, pairs in sorted(per.items()):
        if len(pairs) != 1:
            raise SystemExit(
                f"REFUSE: seed {seed} carries {len(pairs)} distinct null pairs {sorted(pairs)}. The null "
                "is supposed to be a property of the seed; collapsing to one per seed would average over "
                "a real difference.")
        kmer, both = next(iter(pairs))
        out[seed] = {"kmer": kmer, "both": both}
    return out


def ceilings_for_seed(seed: int) -> dict[str, dict]:
    """The discriminator's own noise floor at one corpus seed, per feature mode."""
    from dna_decode.glm.corpus import build_natural_set, fit_eval_halves
    from dna_decode.glm.falsifier import self_control_ceiling

    wins, _stats = build_natural_set(REFSEQ, length=LENGTH, limit=LIMIT, seed=seed, kind=KIND)
    _fit, ev = fit_eval_halves(wins)
    out = {}
    for mode in MODES:
        kw = {"k": KMER_K, "mode": mode}
        if mode == "both":
            kw["bins"] = POSITIONAL_BINS
        ceil_max, vals = self_control_ceiling(ev, **kw)
        out[mode] = {"ceiling": ceil_max, "values": vals, "n_eval": len(ev)}
    return out


def build_rows(nulls_by_seed: dict[int, dict[str, float]],
               ceilings_by_seed: dict[int, dict[str, dict]]) -> list[dict]:
    rows = []
    for seed in sorted(nulls_by_seed):
        for mode in MODES:
            null = nulls_by_seed[seed][mode]
            c = ceilings_by_seed[seed][mode]
            rows.append({
                "seed": seed, "mode": mode,
                "null": null,
                "ceiling": round(c["ceiling"], 4) if c["ceiling"] is not None else None,
                "margin": round(margin(null, c["ceiling"]), 4) if c["ceiling"] is not None else None,
                "ceiling_values": [round(v, 4) for v in c["values"]],
                "n_eval": c["n_eval"],
            })
    return rows


def summarise(rows: list[dict]) -> dict:
    """Per-mode mean/SD of each column, plus the count of seeds with a winnable region."""
    out = {}
    for mode in MODES:
        sel = [r for r in rows if r["mode"] == mode and r["margin"] is not None]
        nulls = [r["null"] for r in sel]
        ceils = [r["ceiling"] for r in sel]
        margins = [r["margin"] for r in sel]
        out[mode] = {
            "n_seeds": len(sel),
            "null_mean": round(statistics.fmean(nulls), 4),
            "null_sd": round(statistics.stdev(nulls), 4) if len(nulls) > 1 else None,
            "null_min": min(nulls), "null_max": max(nulls),
            "ceiling_mean": round(statistics.fmean(ceils), 4),
            "ceiling_sd": round(statistics.stdev(ceils), 4) if len(ceils) > 1 else None,
            "margin_mean": round(statistics.fmean(margins), 4),
            "margin_sd": round(statistics.stdev(margins), 4) if len(margins) > 1 else None,
            "margin_min": min(margins), "margin_max": max(margins),
            "seeds_with_positive_margin": sum(1 for m in margins if m > 0),
            "seeds_with_negative_margin": sum(1 for m in margins if m < 0),
        }
    return out


def verdict_for_mode(s: dict) -> str:
    """Ordered, first match wins. Describes the MODE's winnable region, never a generator's quality."""
    if s["margin_mean"] <= 0:
        return "NO_WINNABLE_REGION_ON_AVERAGE"
    if s["seeds_with_negative_margin"] > 0:
        return "WINNABLE_REGION_SEED_DEPENDENT"
    return "WINNABLE_REGION_AT_EVERY_SEED"


def main() -> int:
    if not (REFSEQ / "genome.fna").exists():
        print(f"REFUSE: no local reference at {REFSEQ}. data/ is gitignored; see "
              "scripts/glm_generator_falsifier.py docstring for the downloader.")
        return 2
    if not GRID.exists():
        print(f"REFUSE: committed grid absent at {GRID}")
        return 2

    grid = json.loads(GRID.read_text(encoding="utf-8"))
    nulls_by_seed = seeds_and_nulls_from_grid(grid)
    print(f"grid: {len(grid['grid'])} cells, {len(nulls_by_seed)} seeds {sorted(nulls_by_seed)}")

    ceilings_by_seed = {}
    for seed in sorted(nulls_by_seed):
        ceilings_by_seed[seed] = ceilings_for_seed(seed)
        c = ceilings_by_seed[seed]
        print(f"  seed {seed}: ceiling kmer {c['kmer']['ceiling']:.4f}  both {c['both']['ceiling']:.4f}")

    # THE REFUSAL -- before any write.
    if ANCHOR_SEED not in ceilings_by_seed:
        print(f"REFUSE: anchor seed {ANCHOR_SEED} is not in the grid seed list")
        return 3
    for mode, expected in COMMITTED.items():
        got = ceilings_by_seed[ANCHOR_SEED][mode]["ceiling"]
        if got is None or abs(got - expected) > GATE_TOL:
            print(f"REFUSE: seed-{ANCHOR_SEED} {mode} ceiling is {got}, committed value is {expected} "
                  f"(tol {GATE_TOL}). A probe that cannot reproduce the one committed ceiling has no "
                  "standing to report the others.")
            return 3
    print(f"gate: seed-{ANCHOR_SEED} ceilings reproduce the committed "
          f"{COMMITTED['kmer']} / {COMMITTED['both']}")

    rows = build_rows(nulls_by_seed, ceilings_by_seed)
    summary = summarise(rows)
    rec = {
        "record": "glm-seed-role-ceilings-v1",
        "date": date.today().isoformat(),
        "question": ("At each of the committed grid seeds, is the Markov null actually ABOVE the "
                     "discriminator own noise floor -- i.e. is there a winnable region at all?"),
        "geometry": {"window_length_bp": LENGTH, "limit_natural": LIMIT, "natural_set_arg": KIND,
                     "kmer_k": KMER_K, "positional_bins": POSITIONAL_BINS, "feature_modes": list(MODES)},
        "source_grid": GRID.name,
        "anchor_gate": {"seed": ANCHOR_SEED, "committed": COMMITTED, "tol": GATE_TOL, "passed": True},
        "why_this_was_missing": (
            "The committed grid cell schema is {model, epoch, seed, both, kmer, null{kmer,both}} with no "
            "ceiling field, so the per-seed margin could not be computed from it. The one committed "
            "ceiling pair (0.5447/0.5865) comes from a DIFFERENT artifact at a single seed, and was read "
            "as though it applied at every seed."),
        "margin_definition": ("null - ceiling. Positive = the null is distinguishable above the noise "
                              "floor, so beating it would be a real claim. Negative = the null sits below "
                              "its own noise floor and the mode has no winnable region. NOT clamped."),
        "rows": rows,
        "summary": summary,
        "verdict_by_mode": {m: verdict_for_mode(summary[m]) for m in MODES},
        "honest_limits": [
            "The nulls are READ from the committed grid, not re-run, so this inherits that grid "
            "generation. The ceilings are freshly computed here.",
            "The ceiling is the MAX over 5 internal self-control seeds (the falsifier own rule), so it "
            "is a conservative floor by construction -- a mean-based floor would be lower and every "
            "margin correspondingly more positive.",
            "This bounds the INSTRUMENT. It says nothing about any generator quality, and a mode with "
            "no winnable region does not mean the generator is good.",
            "One geometry (150bp promoter windows, limit 1200) and one organism. The margin is a property "
            "of a (substrate, geometry, feature-mode) triple, not of discriminators in general.",
        ],
    }
    OUT.write_text(json.dumps(rec, indent=2), encoding="utf-8")

    print()
    for mode in MODES:
        s = summary[mode]
        print(f"[{mode}] null {s['null_mean']:.4f} (sd {s['null_sd']:.4f}, "
              f"range {s['null_min']}-{s['null_max']}) vs ceiling {s['ceiling_mean']:.4f} "
              f"-> margin {s['margin_mean']:+.4f} (sd {s['margin_sd']:.4f})")
        print(f"         positive at {s['seeds_with_positive_margin']}/{s['n_seeds']} seeds, "
              f"negative at {s['seeds_with_negative_margin']}/{s['n_seeds']} "
              f"-> {verdict_for_mode(s)}")
    print(f"\n-> {OUT.relative_to(ROOT).as_posix()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
