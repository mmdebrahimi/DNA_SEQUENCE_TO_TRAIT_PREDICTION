"""Score a trained generator with the COMMITTED falsifier — the north-star number, not a proxy.

    uv run python scripts/glm_score_generator.py --ckpt data/glm_diffusion_ckpt/latest.pt

**This runs `falsifier.run_falsifier`, it does not re-implement it.** Every scoring claim in this arm is
produced by that function, and a second scorer written here would be a second definition of the number
the whole project turns on. So the generator's sequences are handed to the committed path with the
committed geometry, and what comes back is directly comparable to every recorded result.

**What a score here means, stated because the direction is counter-intuitive.** The falsifier reports
DISTINGUISHABILITY from natural sequence, where **LOWER IS BETTER** — it is the AUROC of a discriminator
trying to tell generated sequence from real. Three reference points are reported beside it, and the
generator is only interesting relative to them:

  * the **Markov null** at the same seed, fitted on the same fit half. This is the thing to beat, and
    beating it means scoring BELOW it.
  * the **perfect-generator floor** — the same eval half scored against real held-out natural windows at
    the same class balance. That is what a generator emitting genuine promoter sequence would earn, and
    nothing can honestly go below it.
  * the **uniform control**, which bounds the instrument's power rather than the target.

**Scored in POSITIONAL mode as well as the committed two, deliberately.** Measured earlier today: a
k-mer discriminator is saturated by an order-3 Markov chain BY CONSTRUCTION, because the chain's
sufficient statistic IS the k-mer conditional distribution. So a generator can clear the kmer view by
matching local composition alone. The positional margin is the one no Markov order erodes, and it is
where a generator that has learned WHERE motifs sit will separate from one that has only learned what
they are made of. Reporting only the committed modes would hide exactly the axis this arm cares about.

**One interpretive guard this script enforces in its output.** If the model was trained on a corpus whose
species composition does not match the scoring organism, the score is not a clean statement about the
architecture. The manifest records the composition; this reports the scoring organism's share of training
and REFUSES to print a headline verdict when that share is small, because a poor score there has an
obvious confound that would otherwise get attributed to the model.

No network, no money. CPU.
"""

from __future__ import annotations

import argparse
import json
import statistics
import time
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
REF = ROOT / "data" / "cache" / "refseq" / "GCF_000005845.2"
MANIFEST = ROOT / "wiki" / "glm_training_corpus_manifest_2026-11-21.json"
CKPT = ROOT / "data" / "glm_diffusion_ckpt" / "latest.pt"
OUT_DIR = ROOT / "wiki"

# The committed kernel geometry. Changing any of these breaks comparability with every anchor.
LENGTH = 150
LIMIT = 1200
N_GENERATED = 120
KIND = "promoter"
KMER_K = 4
POSITIONAL_BINS = 6
MARKOV_ORDER = 3
MODES = ("kmer", "positional", "both")
SEEDS = tuple(range(10))

# Below this share of training data, a poor score has an obvious confound.
MIN_SCORING_ORGANISM_SHARE = 0.25


def training_composition(manifest: Path, organism_substring: str = "Escherichia") -> dict:
    """What fraction of the training corpus is the organism we are about to score against?

    The confound this exists to surface is real and was live in this arm: a corpus built for lineage
    diversity across 18 species is 5.4% E. coli, while the falsifier scores only E. coli -- so a model
    can train 95% on organisms it is never graded on and the resulting poor score looks like an
    architecture result.
    """
    if not manifest.exists():
        return {"status": "no-manifest"}
    rec = json.loads(manifest.read_text(encoding="utf-8"))
    man = rec.get("manifest") or {}
    entries = man.get("entries") or []
    total = man.get("n_windows_kept") or sum(e.get("n_windows_kept", 0) for e in entries)
    hit = sum(e.get("n_windows_kept", 0) for e in entries
              if organism_substring.lower() in (e.get("species") or "").lower())
    return {
        "status": "ok",
        "n_windows_total": total,
        "n_windows_scoring_organism": hit,
        "share": round(hit / total, 4) if total else None,
        "organism_substring": organism_substring,
        "leak_gate": man.get("falsifier_leak_gate"),
    }


def corpus_with_spares(seed: int):
    """The committed draw plus sequence-disjoint spares, for the perfect-generator floor."""
    from dna_decode.glm.corpus import build_natural_set, fit_eval_halves

    wins, _ = build_natural_set(REF, length=LENGTH, limit=None, seed=seed, kind=KIND)
    committed = wins[:LIMIT]
    fit, ev = fit_eval_halves(committed)
    used = set(fit) | set(ev)
    spare = []
    for w in wins[LIMIT:]:
        if len(spare) >= N_GENERATED:
            break
        if w.seq in used:
            continue
        spare.append(w.seq)
        used.add(w.seq)
    return fit, ev, spare


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--ckpt", default=str(CKPT))
    p.add_argument("--manifest", default=str(MANIFEST))
    p.add_argument("--samples-file", default=None,
                   help="score sequences from a file instead of generating from a checkpoint")
    p.add_argument("--label", default="diffusion")
    p.add_argument("--sampling-steps", type=int, default=64)
    p.add_argument("--seeds", type=int, default=len(SEEDS))
    p.add_argument("--scoring-organism", default="Escherichia")
    p.add_argument("--out", default=None)
    a = p.parse_args(argv)

    from dna_decode.glm.baselines import MarkovGenerator, UniformGenerator
    from dna_decode.glm.falsifier import heldout_auroc, run_falsifier

    t0 = time.time()
    comp = training_composition(Path(a.manifest), a.scoring_organism)
    if comp.get("status") == "ok":
        print(f"training composition: {comp['n_windows_scoring_organism']:,} of "
              f"{comp['n_windows_total']:,} windows are {a.scoring_organism} "
              f"({comp['share']*100:.1f}%)")
        gate = comp.get("leak_gate") or {}
        print(f"leak gate: {'ENFORCED' if gate.get('enforced') else 'NOT ENFORCED'}")
        if not gate.get("enforced"):
            print("REFUSE: the training corpus manifest does not record an enforced leak gate, so a "
                  "score against the falsifier's own windows would be meaningless.")
            return 3

    # --- obtain the generated sequences -------------------------------------------------------
    if a.samples_file:
        gen = [ln.strip() for ln in Path(a.samples_file).read_text(encoding="utf-8").splitlines()
               if ln.strip()]
        source = {"kind": "file", "path": a.samples_file}
    else:
        import torch

        from dna_decode.glm.diffusion import DiffusionConfig, build_model, sample

        ck = Path(a.ckpt)
        if not ck.exists():
            print(f"REFUSE: no checkpoint at {ck}")
            return 2
        blob = torch.load(ck, map_location="cpu", weights_only=False)
        cfg = DiffusionConfig(**{k: v for k, v in blob["config"].items()
                                 if k in DiffusionConfig.__dataclass_fields__})
        model = build_model(cfg)
        model.load_state_dict(blob["model"])
        print(f"checkpoint: step {blob['step']}, "
              f"{sum(q.numel() for q in model.parameters()):,} params")
        gen = sample(model, N_GENERATED, cfg, steps=a.sampling_steps, seed=0)
        source = {"kind": "checkpoint", "path": str(ck), "step": int(blob["step"]),
                  "sampling_steps": a.sampling_steps}

    bad = [s for s in gen if len(s) != LENGTH or set(s) - set("ACGT")]
    if bad:
        print(f"REFUSE: {len(bad)} of {len(gen)} sequences are malformed")
        return 3
    if len(gen) < N_GENERATED:
        print(f"REFUSE: only {len(gen)} sequences, need {N_GENERATED} for the committed geometry")
        return 3
    gen = gen[:N_GENERATED]
    print(f"scoring {len(gen)} sequences of {LENGTH} bp across {a.seeds} seeds, "
          f"modes {list(MODES)}\n")

    rows = []
    for s in SEEDS[:a.seeds]:
        fit, ev, spare = corpus_with_spares(s)
        null_seqs = MarkovGenerator(order=MARKOV_ORDER, seed=s).fit(fit).generate(N_GENERATED, LENGTH)
        uni = UniformGenerator(seed=s).generate(N_GENERATED, LENGTH)
        for mode in MODES:
            kw = {"k": KMER_K, "mode": mode, "bins": POSITIONAL_BINS}
            r = run_falsifier(ev, gen, candidate_name=a.label, null_sequences=null_seqs,
                              uniform_sequences=uni, length=LENGTH, seed=s, **kw)
            floor = heldout_auroc(ev, spare, k=KMER_K, seed=s, mode=mode, bins=POSITIONAL_BINS)
            rows.append({
                "seed": s, "mode": mode,
                "generator": r.distinguishability,
                "markov_null": r.null_distinguishability,
                "perfect_floor": round(floor, 4) if floor is not None else None,
                "uniform_control": r.control_uniform,
                "self_control_max": r.control_self_max,
                "verdict": r.verdict,
                "beats_null": (r.distinguishability is not None
                               and r.null_distinguishability is not None
                               and r.distinguishability < r.null_distinguishability),
            })

    summary = {}
    for mode in MODES:
        sel = [r for r in rows if r["mode"] == mode and r["generator"] is not None]
        if not sel:
            continue
        g = [r["generator"] for r in sel]
        n = [r["markov_null"] for r in sel if r["markov_null"] is not None]
        f = [r["perfect_floor"] for r in sel if r["perfect_floor"] is not None]
        deltas = [r["generator"] - r["markov_null"] for r in sel
                  if r["markov_null"] is not None]
        summary[mode] = {
            "n_seeds": len(sel),
            "generator_mean": round(statistics.fmean(g), 4),
            "markov_null_mean": round(statistics.fmean(n), 4) if n else None,
            "perfect_floor_mean": round(statistics.fmean(f), 4) if f else None,
            "generator_minus_null_mean": round(statistics.fmean(deltas), 4) if deltas else None,
            "generator_minus_null_sd": (round(statistics.stdev(deltas), 4)
                                        if len(deltas) > 1 else None),
            "seeds_beating_null": sum(1 for r in sel if r["beats_null"]),
            "generator_above_floor": (round(statistics.fmean(g) - statistics.fmean(f), 4)
                                      if f else None),
        }
        e = summary[mode]
        print(f"[{mode:10s}] generator {e['generator_mean']:.4f}  vs null "
              f"{e['markov_null_mean']:.4f}  vs floor {e['perfect_floor_mean']:.4f}   "
              f"delta {e['generator_minus_null_mean']:+.4f}   beats null at "
              f"{e['seeds_beating_null']}/{e['n_seeds']} seeds")

    # --- the interpretive guard ---------------------------------------------------------------
    share = comp.get("share")
    confounded = share is not None and share < MIN_SCORING_ORGANISM_SHARE
    verdict = ("CONFOUNDED_BY_TRAINING_COMPOSITION" if confounded else
               "READABLE")
    print()
    if confounded:
        print(f"VERDICT WITHHELD: {a.scoring_organism} is only {share*100:.1f}% of the training "
              f"corpus (floor {MIN_SCORING_ORGANISM_SHARE*100:.0f}%), so a poor score here has an "
              "obvious confound and must NOT be read as a statement about the architecture. "
              "Fine-tune on the scoring organism first, then re-score.")
    else:
        print(f"Scores are readable: {a.scoring_organism} is {share*100:.1f}% of training.")

    rec = {
        "record": "glm-generator-score-v1",
        "date": date.today().isoformat(),
        "label": a.label,
        "source": source,
        "geometry": {"window_length_bp": LENGTH, "limit_natural": LIMIT,
                     "n_generated_per_candidate": N_GENERATED, "natural_set_arg": KIND,
                     "kmer_k": KMER_K, "positional_bins": POSITIONAL_BINS,
                     "markov_order": MARKOV_ORDER, "feature_modes": list(MODES),
                     "seeds": list(SEEDS[:a.seeds])},
        "scored_by": "dna_decode.glm.falsifier.run_falsifier (the committed scorer, not a copy)",
        "direction": "distinguishability from natural sequence; LOWER IS BETTER",
        "training_composition": comp,
        "interpretive_guard": {
            "verdict": verdict,
            "scoring_organism_share": share,
            "floor": MIN_SCORING_ORGANISM_SHARE,
            "why": ("a model trained mostly on organisms it is never scored against will score poorly "
                    "for a reason unrelated to its architecture; attributing that to the model is the "
                    "mistake this guard exists to prevent"),
        },
        "rows": rows,
        "summary": summary,
        "power_gate_finding": {
            "what": ("positional-ONLY mode FAILS the committed falsifier's power gate, so "
                     "run_falsifier returns INSTRUMENT_UNDERPOWERED and emits NO number for it"),
            "measured_uniform_control_over_5_seeds": {"kmer": 0.9317, "positional": 0.8611,
                                                      "both": 0.9494, "MIN_POWER_AUROC": 0.90},
            "why": ("the gate is calibrated on UNIFORM sequence, which is an easy target that a "
                    "256-dimensional k-mer histogram crushes. Positional features are WEAKER against "
                    "uniform and STRONGER against the Markov null -- the target that actually matters. "
                    "So the gate rejects the mode that is best at the hard task because it is worse at "
                    "the easy one."),
            "consequence": ("`both` passes at 0.9494 AND carries the positional block, so it is the "
                            "practical route to positional signal through the committed scorer. This "
                            "WEAKENS the case for scoring on positional alone, whose advantage over "
                            "`both` was never resolvable (-0.0364, interval crossing zero) and was "
                            "feature economy rather than score."),
            "not_resolved_here": ("whether MIN_POWER_AUROC or the uniform control is the right "
                                  "calibration for positional features is an instrument decision, not "
                                  "this script's to make"),
        },
        "honest_limits": [
            "A low distinguishability means the discriminator could not separate generated from "
            "natural sequence ON THE FEATURES IT SEES. All three modes see LOCAL composition and "
            "positional composition; none sees long-range organisation, so clearing them is NECESSARY "
            "and not sufficient for functional realism.",
            "The generator produces one sample set per run at a single sampling seed; the per-seed "
            "variation reported here is the SCORING seed, not generator sampling variance. A full "
            "account would resample the generator too.",
            "The perfect-generator floor is empirical at one draw per seed and carries its own "
            "sampling variation.",
            "positional mode is NOT in the committed grid, which scores kmer and both. Its numbers "
            "here are new cells, not a re-reading of any recorded result.",
            "Beating the Markov null on these features is the arm's stated bar. It is not evidence "
            "that generated sequence is biologically functional, which needs a wet-lab or an "
            "oracle-scored design task.",
        ],
        "wall_seconds": round(time.time() - t0, 1),
    }
    out = Path(a.out) if a.out else OUT_DIR / f"glm_generator_score_{a.label}_{date.today()}.json"
    out = out if out.is_absolute() else (ROOT / out)
    out.write_text(json.dumps(rec, indent=2), encoding="utf-8")
    # Display relative when possible, absolute otherwise. An earlier version called `relative_to`
    # unconditionally, which raises on a path outside the repo -- AFTER the artifact was written and
    # every measurement made, so the run exited non-zero having actually succeeded. Same shape as the
    # `t0` NameError earlier today: a display-layer bug at the very end of a long job.
    try:
        shown = out.relative_to(ROOT).as_posix()
    except ValueError:
        shown = str(out)
    print(f"\n{rec['wall_seconds']}s -> {shown}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
