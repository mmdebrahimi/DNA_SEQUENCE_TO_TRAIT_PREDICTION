"""Train the masked-diffusion DNA generator on the leak-gated corpus, checkpointing as it goes.

    uv run python scripts/glm_train_diffusion.py --smoke          # seconds, proves the pipeline
    uv run python scripts/glm_train_diffusion.py --steps 20000    # the real de-risk run
    uv run python scripts/glm_train_diffusion.py --resume         # continue from the last checkpoint

**It refuses a corpus whose leak gate was not enforced.** The corpus manifest records whether the
falsifier's scoring windows were withheld. If they were not, training proceeds happily and the resulting
generator scores beautifully against sequences it was trained on, with no error anywhere. That is the one
failure in this arm that manufactures success, so the check is a refusal at startup rather than a note in
a report, and `--allow-leaky-corpus` exists only so that overriding it has to be written down.

**Checkpointing is frequent and resuming is real, because the intended host is a shared gaming box.**
Every `--ckpt-every` steps the model, optimizer, step count and config are written atomically, so a run
killed for a game loses minutes rather than hours. `--resume` restores all four; the config is compared
against the checkpoint's and a mismatch refuses, because silently resuming a different architecture from
an old optimizer state produces a model that trains without obviously breaking.

**A held-out split comes from the TRAINING corpus, not from the falsifier's windows.** Those are already
withheld and must stay unseen -- using them to watch for overfitting would be the leak by another route.
So a fraction of the training corpus is held back purely to report eval loss.

**Tokens per parameter is reported, not enforced.** The corpus supports a given model size only under a
scaling heuristic imported from language models and not validated for DNA. The number is printed so the
ratio is visible; it is not a gate, because a de-risk run that is deliberately over-parameterised is a
legitimate thing to do and the point here is to prove the pipeline end to end.

CPU by default. No network, no download, no money.
"""

from __future__ import annotations

import argparse
import json
import math
import time
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
CORPUS = ROOT / "wiki" / "glm_training_corpus_2026-11-21.txt"
MANIFEST = ROOT / "wiki" / "glm_training_corpus_manifest_2026-11-21.json"
CKPT_DIR = ROOT / "data" / "glm_diffusion_ckpt"
OUT_DIR = ROOT / "wiki"

EVAL_FRACTION = 0.02
TOKENS_PER_PARAM_HEURISTIC = 20


def load_corpus(path: Path, manifest: Path, *, allow_leaky: bool) -> tuple[list[str], dict]:
    """Read the corpus and REFUSE unless its manifest says the leak gate was enforced."""
    if not path.exists():
        raise SystemExit(f"REFUSE: corpus not found at {path}. Build it with "
                         "scripts/glm_build_corpus.py")
    man = {}
    if manifest.exists():
        man = json.loads(manifest.read_text(encoding="utf-8"))
    gate = (man.get("manifest") or {}).get("falsifier_leak_gate") or {}
    enforced = bool(gate.get("enforced"))
    if not enforced and not allow_leaky:
        raise SystemExit(
            "REFUSE: the corpus manifest does not record an enforced leak gate "
            f"({manifest.name}). Training on the windows the falsifier scores against produces a "
            "generator with excellent numbers and no error anywhere, which is the one failure in this "
            "arm that manufactures success. Rebuild the corpus with the gate on, or pass "
            "--allow-leaky-corpus to state in writing that this run is not scoreable.")
    seqs = [ln.strip() for ln in path.read_text(encoding="utf-8").splitlines() if ln.strip()]
    if not seqs:
        raise SystemExit(f"REFUSE: corpus at {path} is empty")
    return seqs, man


def split_train_eval(seqs: list[str], fraction: float, seed: int) -> tuple[list[str], list[str]]:
    """Deterministic held-out split for eval loss only.

    Taken from the TRAINING corpus, never from the falsifier's withheld windows -- those must stay
    unseen, and using them to watch for overfitting would reintroduce the leak by another route.
    """
    import random

    idx = list(range(len(seqs)))
    random.Random(seed).shuffle(idx)
    n_eval = max(1, int(len(seqs) * fraction))
    ev = [seqs[i] for i in idx[:n_eval]]
    tr = [seqs[i] for i in idx[n_eval:]]
    return tr, ev


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--corpus", default=str(CORPUS))
    p.add_argument("--manifest", default=str(MANIFEST))
    p.add_argument("--steps", type=int, default=20000)
    p.add_argument("--batch", type=int, default=64)
    p.add_argument("--lr", type=float, default=3e-4)
    p.add_argument("--d-model", type=int, default=128)
    p.add_argument("--n-layers", type=int, default=4)
    p.add_argument("--n-heads", type=int, default=4)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--ckpt-every", type=int, default=500)
    p.add_argument("--eval-every", type=int, default=250)
    p.add_argument("--ckpt-dir", default=str(CKPT_DIR))
    p.add_argument("--resume", action="store_true")
    p.add_argument("--smoke", action="store_true",
                   help="tiny model, tiny corpus, few steps: proves the pipeline in seconds")
    p.add_argument("--allow-leaky-corpus", action="store_true")
    p.add_argument("--n-samples", type=int, default=120,
                   help="sequences to generate at the end, in the falsifier's geometry")
    a = p.parse_args(argv)

    import torch

    from dna_decode.glm.diffusion import DiffusionConfig, build_model, encode, loss_fn, sample

    seqs, man = load_corpus(Path(a.corpus), Path(a.manifest),
                            allow_leaky=a.allow_leaky_corpus)
    length = len(seqs[0])
    if a.smoke:
        # --smoke OVERRIDES other flags, and says so. Silently clobbering an explicit --d-model is how
        # a verification run ends up testing something other than what was asked for: a config-mismatch
        # check written as `--smoke --d-model 64` quietly got d_model 32 and therefore no mismatch,
        # reporting a pass for a guard it never exercised.
        overrides = {"corpus_windows": 512, "steps": 60, "batch": 16, "ckpt_every": 30,
                     "eval_every": 20, "d_model": 32, "n_layers": 2, "n_heads": 2, "n_samples": 8}
        changed = {k: (getattr(a, k), v) for k, v in overrides.items()
                   if hasattr(a, k) and getattr(a, k) != v}
        seqs = seqs[:512]
        a.steps, a.batch, a.ckpt_every, a.eval_every = 60, 16, 30, 20
        a.d_model, a.n_layers, a.n_heads, a.n_samples = 32, 2, 2, 8
        if changed:
            print("--smoke overrides: "
                  + ", ".join(f"{k} {was}->{now}" for k, (was, now) in sorted(changed.items())))

    cfg = DiffusionConfig(length=length, d_model=a.d_model, n_layers=a.n_layers,
                          n_heads=a.n_heads, d_ff=4 * a.d_model, seed=a.seed)
    train, ev = split_train_eval(seqs, EVAL_FRACTION, a.seed)
    model = build_model(cfg)
    n_params = sum(q.numel() for q in model.parameters())
    corpus_tokens = len(train) * length

    print(f"corpus: {len(seqs):,} windows x {length} bp "
          f"({len(train):,} train / {len(ev):,} eval)")
    print(f"leak gate: {'ENFORCED' if (man.get('manifest') or {}).get('falsifier_leak_gate', {}).get('enforced') else 'NOT ENFORCED (overridden)'}")
    print(f"model: {n_params:,} params  (d_model={cfg.d_model} layers={cfg.n_layers} "
          f"heads={cfg.n_heads})")
    print(f"tokens/param: {corpus_tokens / n_params:.1f} "
          f"(heuristic target {TOKENS_PER_PARAM_HEURISTIC}; imported from language models, "
          f"NOT validated for DNA)")

    opt = torch.optim.AdamW(model.parameters(), lr=a.lr, weight_decay=0.01)
    ckpt_dir = Path(a.ckpt_dir)
    ckpt_dir.mkdir(parents=True, exist_ok=True)
    ckpt_path = ckpt_dir / "latest.pt"
    start = 0
    if a.resume and ckpt_path.exists():
        blob = torch.load(ckpt_path, map_location="cpu", weights_only=False)
        if blob["config"] != cfg.as_dict():
            raise SystemExit(
                "REFUSE: the checkpoint's config differs from this run's. Resuming a different "
                "architecture from an old optimizer state produces a model that trains without "
                "obviously breaking.\n"
                f"  checkpoint: {blob['config']}\n  this run:   {cfg.as_dict()}")
        model.load_state_dict(blob["model"])
        opt.load_state_dict(blob["optimizer"])
        start = int(blob["step"])
        print(f"resumed from step {start}")

    x_train = encode(train, length)
    x_eval = encode(ev, length)
    g = torch.Generator().manual_seed(a.seed)
    t0 = time.time()
    history = []

    def save(step: int) -> None:
        tmp = ckpt_path.with_suffix(".tmp")
        torch.save({"model": model.state_dict(), "optimizer": opt.state_dict(),
                    "step": step, "config": cfg.as_dict()}, tmp)
        tmp.replace(ckpt_path)

    model.train()
    for step in range(start, a.steps):
        idx = torch.randint(0, x_train.shape[0], (a.batch,), generator=g)
        loss = loss_fn(model, x_train[idx], generator=g)
        opt.zero_grad(set_to_none=True)
        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
        opt.step()

        if (step + 1) % a.eval_every == 0 or step == a.steps - 1:
            model.eval()
            with torch.no_grad():
                ge = torch.Generator().manual_seed(12345)
                ev_loss = float(loss_fn(model, x_eval[:min(256, x_eval.shape[0])], generator=ge))
            model.train()
            row = {"step": step + 1, "train_loss": round(float(loss), 4),
                   "eval_loss": round(ev_loss, 4),
                   "seconds": round(time.time() - t0, 1)}
            history.append(row)
            print(f"  step {step+1:6d}  train {row['train_loss']:8.4f}  eval {row['eval_loss']:8.4f}"
                  f"  {row['seconds']:7.1f}s")
        if (step + 1) % a.ckpt_every == 0:
            save(step + 1)
    save(a.steps)

    print(f"\ngenerating {a.n_samples} sequences at {length} bp...")
    gen = sample(model, a.n_samples, cfg, steps=64, seed=a.seed)
    bad = [s for s in gen if len(s) != length or set(s) - set("ACGT")]
    if bad:
        print(f"REFUSE: {len(bad)} generated sequences are malformed")
        return 3
    gc = sum(s.count("G") + s.count("C") for s in gen) / (len(gen) * length)
    gc_nat = sum(s.count("G") + s.count("C") for s in train[:2000]) / (min(2000, len(train)) * length)
    print(f"  GC generated {gc:.4f} vs corpus {gc_nat:.4f}")

    samples_path = OUT_DIR / ("glm_diffusion_samples_smoke.txt" if a.smoke
                              else "glm_diffusion_samples_2026-11-21.txt")
    samples_path.write_text("\n".join(gen) + "\n", encoding="utf-8")

    rec = {
        "record": "glm-diffusion-train-v1",
        "date": date.today().isoformat(),
        "mode": "smoke" if a.smoke else "full",
        "config": cfg.as_dict(),
        "n_params": n_params,
        "corpus": {"file": Path(a.corpus).name, "n_windows": len(seqs),
                   "n_train": len(train), "n_eval": len(ev), "length_bp": length,
                   "tokens_per_param": round(corpus_tokens / n_params, 2),
                   "leak_gate": (man.get("manifest") or {}).get("falsifier_leak_gate")},
        "optim": {"steps": a.steps, "batch": a.batch, "lr": a.lr, "optimizer": "AdamW",
                  "weight_decay": 0.01, "grad_clip": 1.0},
        "history": history,
        "final": history[-1] if history else None,
        "samples_file": samples_path.name,
        "gc": {"generated": round(gc, 4), "corpus_subset": round(gc_nat, 4)},
        "honest_limits": [
            "GC agreement is a necessary sanity check and nothing more -- a model that matched GC and "
            "nothing else would pass it. The real score is the falsifier's held-out AUROC against the "
            "withheld natural windows, which is a separate run.",
            "Eval loss is computed on a held-out slice of the TRAINING corpus, not on the falsifier's "
            "windows. It tracks overfitting to the training distribution and is not a generator score.",
            "Tokens per parameter is reported against a heuristic imported from language-model scaling "
            "and NOT validated for DNA. It is not enforced.",
            "One architecture, one corpus, one seed. Nothing here establishes that diffusion beats an "
            "autoregressive model on this substrate; it establishes that the pipeline runs end to end "
            "and produces well-formed sequence in the geometry the falsifier scores.",
        ],
        "wall_seconds": round(time.time() - t0, 1),
    }
    out = OUT_DIR / ("glm_diffusion_train_smoke.json" if a.smoke
                     else "glm_diffusion_train_2026-11-21.json")
    out.write_text(json.dumps(rec, indent=2), encoding="utf-8")
    print(f"\n{rec['wall_seconds']}s -> {out.relative_to(ROOT).as_posix()}")
    print(f"              -> {samples_path.relative_to(ROOT).as_posix()}")
    print(f"              -> {ckpt_path.relative_to(ROOT).as_posix()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
