"""Measure, on a real Kaggle T4, what size of Enformer-shaped DNA model is trainable here.

WHY THIS EXISTS. The project's most ambitious goal was given up on 2026-05-12: "a tissue-specific
gene-regulatory-element activity predictor using a transformer-based DNA model trained on
ENCODE / Roadmap Epigenomics multi-tissue data (Enformer-shaped)". The cut was on LOCAL compute
(GTX 860M, 4 GiB). Constraint `local-gpu-vram-blocks-dna-lm-finetune` recorded it with a
`gpu-vram-gte 8` expiry test, and that test CANNOT FIRE: it is excluded from the default scan by a
`resolved` flag whose `last_checked` is empty, and forced in it returns `unknown` ("no GPU visible")
because the compute that actually arrived is CLOUD, not a new local card.

So the open question is NOT "is there an 8 GB GPU" -- there is, and the route is demonstrated. It is
"what does the available cloud envelope actually fit, and how fast does it run". NOBODY HAS MEASURED
IT. Constraint `encode-roadmap-regulatory-goal-unsized-on-cloud` expires on
`wiki/encode_roadmap_cloud_sizing.md`, and this kernel produces the numbers that memo needs.

WHAT IT MEASURES, AND WHY IT USES RANDOM TENSORS. The question is COMPUTE CAPACITY, so the inputs are
random one-hot DNA of the right SHAPE. That is deliberate: it isolates the compute question from the
separate and still-unmeasured question of whether ENCODE / Roadmap data is reachable on this host. A
sizing run that also had to fetch multi-tissue tracks would conflate two walls and settle neither.

Per configuration it runs a real forward AND backward pass and records peak VRAM and milliseconds per
step, sweeping context length x width x depth across three memory regimes: fp32, AMP fp16 (the T4 has
tensor cores at sm_75), and AMP fp16 PLUS gradient checkpointing. Checkpointing is included
deliberately -- it is THE standard lever for long context, and omitting it would under-report
feasibility, which is the wrong direction of error for a question about re-opening a goal.

OOM is a RESULT, not a crash: each config is caught, the cache is emptied, and the sweep continues,
so the output is an envelope rather than a single pass/fail.

WHAT IT DOES NOT SETTLE, stated here so the artifact cannot be over-read:
  * whether ENCODE / Roadmap data is reachable here -- not touched, separate constraint;
  * whether the measured envelope SUFFICES for the scientific goal -- that needs a sourced training
    budget and a target metric, neither of which this kernel has. It reports what fits and how fast.
    Reading "it fits" as "the goal is feasible" is exactly the over-read to avoid.
  * anything about accuracy. No model here is trained to convergence on anything.
"""
import json
import os
import platform
import sys
import time
from pathlib import Path

OUT = Path("/kaggle/working")
BIN_BP = 128          # Enformer-style binning: this many bp per output bin
SESSION_HOURS = 9.0   # Kaggle GPU sessions cap near here
WARMUP, MEASURE = 2, 3


def build_model(torch, nn, *, channels: int, depth: int, heads: int, n_tracks: int):
    """A deliberately plain Enformer-SHAPED net: conv tower to bins, transformer, per-bin head.

    Not a reimplementation of Enformer and not trying to be -- the point is the memory/throughput
    SHAPE of that architecture class (bp-resolution conv stem, attention over bins, multi-track
    head), which is what determines whether it fits.
    """
    n_halvings = 7  # 2**7 = 128 bp per bin

    class Res(nn.Module):
        def __init__(self, c):
            super().__init__()
            self.b = nn.Sequential(nn.BatchNorm1d(c), nn.GELU(), nn.Conv1d(c, c, 5, padding=2))

        def forward(self, x):
            return x + self.b(x)

    layers = [nn.Conv1d(4, channels, 15, padding=7)]
    for _ in range(n_halvings):
        layers += [Res(channels), nn.MaxPool1d(2)]
    tower = nn.Sequential(*layers)

    block = nn.TransformerEncoderLayer(d_model=channels, nhead=heads,
                                       dim_feedforward=channels * 2, batch_first=True,
                                       norm_first=True, dropout=0.0)
    trunk = nn.TransformerEncoder(block, num_layers=depth)
    head = nn.Linear(channels, n_tracks)

    class Net(nn.Module):
        def __init__(self):
            super().__init__()
            self.tower, self.trunk, self.head = tower, trunk, head

        def forward(self, x, *, checkpoint=False):
            h = self.tower(x).transpose(1, 2)            # B, L_bins, C
            if checkpoint:
                from torch.utils.checkpoint import checkpoint as ck
                for layer in self.trunk.layers:
                    h = ck(layer, h, use_reentrant=False)
            else:
                h = self.trunk(h)
            return self.head(h)

    return Net()


def try_config(torch, nn, *, context_bp: int, channels: int, depth: int, heads: int,
               n_tracks: int, mode: str) -> dict:
    """Run a real fwd+bwd. Returns a record; OOM is reported, never raised."""
    rec = {"context_bp": context_bp, "channels": channels, "depth": depth,
           "n_tracks": n_tracks, "bins": context_bp // BIN_BP, "mode": mode}
    amp = mode in ("amp", "amp+ckpt")
    ckpt = mode == "amp+ckpt"
    try:
        torch.cuda.empty_cache()
        torch.cuda.reset_peak_memory_stats()
        model = build_model(torch, nn, channels=channels, depth=depth, heads=heads,
                            n_tracks=n_tracks).cuda()
        rec["params_m"] = round(sum(p.numel() for p in model.parameters()) / 1e6, 2)
        opt = torch.optim.Adam(model.parameters(), lr=1e-4)
        # The newer torch.amp API: torch.cuda.amp.* is deprecated in 2.4+ and warns loudly on the
        # 2.11 image the sibling kernel reported. Fall back only if this build predates it.
        try:
            scaler = torch.amp.GradScaler("cuda", enabled=amp)
            autocast = lambda: torch.amp.autocast("cuda", enabled=amp)  # noqa: E731
        except (AttributeError, TypeError):
            scaler = torch.cuda.amp.GradScaler(enabled=amp)
            autocast = lambda: torch.cuda.amp.autocast(enabled=amp)  # noqa: E731
        x = torch.randn(1, 4, context_bp, device="cuda")
        y = torch.randn(1, context_bp // BIN_BP, n_tracks, device="cuda")

        def step():
            opt.zero_grad(set_to_none=True)
            with autocast():
                loss = ((model(x, checkpoint=ckpt) - y) ** 2).mean()
            scaler.scale(loss).backward()
            scaler.step(opt)
            scaler.update()

        for _ in range(WARMUP):
            step()
        torch.cuda.synchronize()
        t0 = time.time()
        for _ in range(MEASURE):
            step()
        torch.cuda.synchronize()
        ms = (time.time() - t0) / MEASURE * 1000.0

        rec["fits"] = True
        rec["ms_per_step"] = round(ms, 1)
        rec["peak_vram_gib"] = round(torch.cuda.max_memory_allocated() / 2**30, 2)
        rec["steps_per_session"] = int(SESSION_HOURS * 3600 * 1000 / ms)
        del model, opt, x, y
    except torch.cuda.OutOfMemoryError:
        rec["fits"] = False
        rec["reason"] = "CUDA OOM"
    except RuntimeError as e:
        msg = str(e)
        rec["fits"] = False
        rec["reason"] = ("CUDA OOM" if "out of memory" in msg.lower() else f"RuntimeError: {msg[:160]}")
    finally:
        torch.cuda.empty_cache()
    return rec


def main() -> int:
    t0 = time.time()
    OUT.mkdir(parents=True, exist_ok=True)
    import torch
    from torch import nn

    if not torch.cuda.is_available():
        print("REFUSING: no CUDA device. This measurement is meaningless on CPU.", file=sys.stderr)
        return 2

    name = torch.cuda.get_device_name(0)
    cap = torch.cuda.get_device_capability(0)
    total = torch.cuda.get_device_properties(0).total_memory / 2**30
    env = {"gpu": name, "capability": f"sm_{cap[0]}{cap[1]}", "total_vram_gib": round(total, 2),
           "torch": torch.__version__, "python": sys.version.split()[0],
           "platform": platform.platform(), "session_hours_assumed": SESSION_HOURS,
           "bin_bp": BIN_BP}
    # The documented trap: `enable_gpu` alone can provision a P100 (sm_60). Record what we GOT, and
    # say so loudly rather than silently sizing against the wrong card.
    env["is_t4"] = "T4" in name
    env["card_note"] = ("T4 as pinned" if env["is_t4"] else
                        f"NOT a T4 -- got {name}. Numbers describe THIS card, not a T4.")
    print("ENV", json.dumps(env), flush=True)

    # Sweep. Context grows toward Enformer's 196,608 bp; width/depth grow toward its ~250M class.
    contexts = [16_384, 32_768, 65_536, 131_072, 196_608]
    shapes = [(384, 4), (768, 8), (1024, 11)]
    n_tracks = 5313          # Enformer's output-track count, the realistic head width
    results = []
    for channels, depth in shapes:
        heads = 8 if channels % 8 == 0 else 4
        for ctx in contexts:
            for mode in ("fp32", "amp", "amp+ckpt"):
                r = try_config(torch, nn, context_bp=ctx, channels=channels, depth=depth,
                               heads=heads, n_tracks=n_tracks, mode=mode)
                results.append(r)
                tag = (f"{r['peak_vram_gib']} GiB, {r['ms_per_step']} ms/step, "
                       f"{r['steps_per_session']:,} steps/{SESSION_HOURS:g}h"
                       if r["fits"] else r["reason"])
                print(f"  C{channels} D{depth} ctx{ctx:>7} {mode:9s} -> "
                      f"{'FITS ' if r['fits'] else 'no   '} {tag}", flush=True)

    fits = [r for r in results if r["fits"]]
    rec = {
        "record": "encode-roadmap-cloud-sizing-v1",
        "env": env,
        "wall_seconds": round(time.time() - t0, 1),
        "n_configs": len(results),
        "n_fitting": len(fits),
        "results": results,
    }
    if fits:
        widest = max(fits, key=lambda r: (r["context_bp"], r["params_m"]))
        fastest_big = max(fits, key=lambda r: (r["context_bp"], r["steps_per_session"]))
        rec["largest_context_that_fits_bp"] = widest["context_bp"]
        rec["largest_fitting_config"] = widest
        rec["best_throughput_at_largest_context"] = fastest_big
        rec["enformer_context_bp"] = 196_608
        rec["reaches_enformer_context"] = widest["context_bp"] >= 196_608
    rec["WHAT_THIS_DOES_NOT_SETTLE"] = [
        "Whether ENCODE / Roadmap data is reachable on this host. Not touched here: the inputs are "
        "random tensors of the right SHAPE, deliberately, so a compute answer is not entangled with "
        "a data-availability answer.",
        "Whether the measured envelope SUFFICES for the goal. That needs a sourced training budget "
        "and a target metric, neither of which this kernel has. It reports what FITS and HOW FAST. "
        "Reading 'it fits' as 'the goal is feasible' is the over-read this list exists to block.",
        "Anything about accuracy. Nothing here is trained to convergence on anything, and the loss "
        "is against random targets.",
        "Multi-GPU or paid-tier scaling. One free T4, one session.",
    ]
    (OUT / "encode_sizing.json").write_text(json.dumps(rec, indent=2), encoding="utf-8")
    print("\n=== summary ===", flush=True)
    print(f"card: {env['gpu']} ({env['total_vram_gib']} GiB, {env['capability']}) — "
          f"{env['card_note']}", flush=True)
    print(f"configs tried {len(results)}, fitting {len(fits)}", flush=True)
    if fits:
        w = rec["largest_fitting_config"]
        print(f"largest context that fits: {w['context_bp']:,} bp "
              f"({w['params_m']}M params, {w['mode']}, {w['peak_vram_gib']} GiB, "
              f"{w['steps_per_session']:,} steps per {SESSION_HOURS:g} h)", flush=True)
        print(f"reaches Enformer's 196,608 bp context: {rec['reaches_enformer_context']}", flush=True)
    else:
        print("NOTHING FIT — that is itself the answer for this envelope.", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
