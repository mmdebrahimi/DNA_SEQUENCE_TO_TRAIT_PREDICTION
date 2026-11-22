"""STANDALONE activity-conditioned masked-diffusion DNA generator. Trains, samples, scores NOTHING.

Target host: a shared gaming box with a GTX 1070 (Pascal, sm_61, 8 GB). Deps: torch only.

===============================================================================================
THIS SCRIPT DELIBERATELY DOES NOT SCORE. Read this before adding a metric to it.
===============================================================================================
The judge is an oracle trained on folds 7-8 and the decision quantities are frozen in
`wiki/human_design_bar_2026-11-22.json`. If this box scored its own output with its own re-derived
oracle, every number it produced would be a different generator seen through a different judge, and
the Kaggle run and these runs would not be comparable. So this emits SAMPLES and a MANIFEST; the
repo host scores every run with the one frozen oracle on the one reserved shard.

It is also why this script never sees folds 7-10 -- see the shard manifest's
`folds_DELIBERATELY_ABSENT`. A fold that is not on the machine cannot be read by a later edit.

===============================================================================================
WHAT IT MIRRORS, AND THE ONE THING IT FIXES
===============================================================================================
The formulation is copied from `dna_decode/glm/diffusion.py` so the two machines train the same
objective: absorbing-state masked diffusion, t ~ U[eps,1], mask each position i.i.d. with prob t,
predict the ORIGINAL base at masked positions, cross-entropy weighted by 1/t. The model is NOT fed
t -- corruption level is implicit in the mask density (the MDLM time-independent parameterisation),
and that matches the reference exactly.

ONE REAL FIX, and it applies to the Kaggle kernel too: the reference `sample()` builds
`torch.Generator(device="cpu")` and then calls `torch.rand(x.shape, generator=g)` with no device, so
on CUDA it throws a device mismatch. The shipped sampler is CPU-only by construction. Here the
generator and every draw live on the model's device. Consequence to accept: CUDA and CPU RNG streams
differ, so a given seed does NOT reproduce across devices. The seed reproduces on the SAME device,
which is what seed-spread needs.

ONE ADDITION: conditioning. An `nn.Embedding(n_bins, d_model)` added to every position, which is
the whole conditioning mechanism. The bin arrives as a COLUMN in the shard -- this script never
derives a quintile boundary, strips an adaptor or decides a fold.

===============================================================================================
--shuffle-labels IS NOT A DEBUG FLAG. It is the arm that makes a positive result readable.
===============================================================================================
It permutes the bin column across training rows, preserving the label MARGINAL while destroying the
sequence-label association. A model trained this way must show ~zero target monotonicity. Measuring
against it is "differs from RANDOM", not "differs from zero" -- the difference that caught a
near-harmless-looking prompt leak on the E. coli arm. Run it at MATCHED capacity and matched steps,
or it compares two things at once.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import random
import signal
import sys
import time
from pathlib import Path

BASES = "ACGT"
MASK_ID = 4
VOCAB = 5
ELEMENT_LEN = 200
N_BINS = 5

#: Named points on the model-size curve. Same `--steps` for every size ON PURPOSE: a learning curve
#: that varies capacity AND token budget at once cannot attribute a difference to either. Fixed
#: steps at fixed batch = fixed tokens seen, so size is the only moving part.
SIZES = {
    "s":  dict(d_model=128, n_layers=4,  n_heads=4, d_ff=512),
    "m":  dict(d_model=192, n_layers=6,  n_heads=6, d_ff=768),
    "l":  dict(d_model=256, n_layers=8,  n_heads=8, d_ff=1024),
    "xl": dict(d_model=384, n_layers=10, n_heads=8, d_ff=1536),
}

_STOP = {"flag": False, "why": ""}


def _on_signal(signum, _frame):
    _STOP["flag"] = True
    _STOP["why"] = "signal:%d" % signum
    print("\n[yield] signal %d -- will checkpoint and exit 0 at the next step boundary" % signum,
          flush=True)


def n_params_estimate(d_model: int, n_layers: int, d_ff: int, length: int = ELEMENT_LEN) -> int:
    emb = (VOCAB + length + N_BINS) * d_model
    per_layer = 4 * d_model * d_model + 2 * d_model * d_ff
    return emb + n_layers * per_layer + d_model * len(BASES)


def md5(path: Path) -> str:
    h = hashlib.md5()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load_shard(path: Path, expect_md5: str | None):
    """Read `bin<TAB>sequence`. Verifies length and alphabet per row and REFUSES on the first bad one.

    A silent coercion here is the whole ballgame: an N mapped to a base teaches the model that an
    ambiguity code is real sequence, and a short row silently padded teaches it that elements end
    early. The export already filtered both; this re-checks because the file crossed a machine.
    """
    if not path.is_file():
        raise SystemExit("REFUSE: shard %s is absent" % path)
    got = md5(path)
    if expect_md5 and got != expect_md5:
        raise SystemExit(
            "REFUSE: shard md5 %s != expected %s.\nThe conditioning LABEL is the one quantity that "
            "must be identical across machines; a different shard vintage silently changes the "
            "quintile boundaries and makes these runs incomparable to the main run." % (got, expect_md5))
    bins, seqs = [], []
    with path.open("r", encoding="utf-8") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        if header[:2] != ["bin", "sequence"]:
            raise SystemExit("REFUSE: shard header is %r, expected ['bin', 'sequence']" % header[:2])
        for ln, line in enumerate(fh, start=2):
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) != 2:
                raise SystemExit("REFUSE: shard line %d has %d fields" % (ln, len(parts)))
            b, s = parts
            s = s.upper()
            if len(s) != ELEMENT_LEN:
                raise SystemExit("REFUSE: shard line %d sequence is %d bp, expected %d"
                                 % (ln, len(s), ELEMENT_LEN))
            if set(s) - set(BASES):
                raise SystemExit("REFUSE: shard line %d has a non-ACGT base" % ln)
            bi = int(b)
            if not 0 <= bi < N_BINS:
                raise SystemExit("REFUSE: shard line %d bin %d outside 0..%d" % (ln, bi, N_BINS - 1))
            bins.append(bi)
            seqs.append(s)
    if not seqs:
        raise SystemExit("REFUSE: shard parsed zero rows")
    return bins, seqs, got


def assert_pascal_ok(torch, device: str) -> dict:
    """Refuse loudly if the wheel has no kernel for this GPU, instead of discovering it mid-run.

    A torch built for CUDA 12.8+ dropped Pascal; the symptom is `no kernel image is available for
    execution on the device`, which can surface late. sm_61 has no usable tensor cores and fp16 is
    not faster than fp32 there, so this runs fp32 and does not enable AMP.
    """
    info = {"torch": torch.__version__, "cuda_available": bool(torch.cuda.is_available())}
    if device == "cpu":
        info["device"] = "cpu"
        return info
    if not torch.cuda.is_available():
        raise SystemExit("REFUSE: --device cuda but torch.cuda.is_available() is False. Install a "
                         "cu126-or-older wheel; CUDA 12.8+ builds dropped Pascal (sm_61).")
    cap = torch.cuda.get_device_capability(0)
    name = torch.cuda.get_device_name(0)
    archs = torch.cuda.get_arch_list()
    sm = "sm_%d%d" % cap
    info.update(device="cuda", gpu=name, capability=sm, arch_list=archs,
                total_vram_gb=round(torch.cuda.get_device_properties(0).total_memory / 1e9, 2))
    if sm not in archs and "compute_%d%d" % cap not in archs:
        raise SystemExit(
            "REFUSE: this torch (%s) was compiled for %s and carries no kernel for %s (%s).\n"
            "Install a cu126 or older wheel -- CUDA 12.8+ builds dropped Pascal."
            % (torch.__version__, archs, sm, name))
    print("[gpu] %s  %s  %.1f GB  torch %s  archs=%s"
          % (name, sm, info["total_vram_gb"], torch.__version__, archs), flush=True)
    return info


def build_model(torch, cfg: dict, seed: int):
    from torch import nn
    torch.manual_seed(seed)

    class CondMaskedDiffusionDNA(nn.Module):
        def __init__(self) -> None:
            super().__init__()
            d = cfg["d_model"]
            self.tok = nn.Embedding(VOCAB, d)
            self.pos = nn.Embedding(ELEMENT_LEN, d)
            self.cond = nn.Embedding(N_BINS, d)
            layer = nn.TransformerEncoderLayer(
                d_model=d, nhead=cfg["n_heads"], dim_feedforward=cfg["d_ff"],
                dropout=cfg["dropout"], batch_first=True, norm_first=True, activation="gelu")
            self.enc = nn.TransformerEncoder(layer, num_layers=cfg["n_layers"])
            self.norm = nn.LayerNorm(d)
            self.head = nn.Linear(d, len(BASES))

        def forward(self, x, c):
            idx = torch.arange(x.shape[1], device=x.device)
            h = self.tok(x) + self.pos(idx)[None, :, :] + self.cond(c)[:, None, :]
            return self.head(self.norm(self.enc(h)))

    return CondMaskedDiffusionDNA()


def loss_fn(torch, model, x0, c, eps: float = 1e-3):
    """1/t-weighted CE over masked positions only. Mirrors `diffusion.loss_fn` exactly."""
    import torch.nn.functional as F
    b = x0.shape[0]
    t = torch.rand(b, device=x0.device) * (1.0 - eps) + eps
    mask = torch.rand(x0.shape, device=x0.device) < t[:, None].expand_as(x0)
    xt = x0.masked_fill(mask, MASK_ID)
    logits = model(xt, c)
    if not bool(mask.any()):
        return logits.sum() * 0.0
    ce = F.cross_entropy(logits[mask], x0[mask], reduction="none")
    w = (1.0 / t)[:, None].expand_as(mask)[mask]
    return (ce * w).sum() / mask.sum()


def sample(torch, model, cfg: dict, bins, steps: int, seed: int, temperature: float = 1.0):
    """Reverse the masking process from a fully masked start, conditioned on `bins`.

    Device-correct, unlike the reference: the generator and every draw live on the model's device.
    """
    dev = next(model.parameters()).device
    model.eval()
    g = torch.Generator(device=dev.type).manual_seed(seed)
    n = len(bins)
    c = torch.tensor(bins, dtype=torch.long, device=dev)
    x = torch.full((n, ELEMENT_LEN), MASK_ID, dtype=torch.long, device=dev)
    ts = torch.linspace(1.0, 0.0, steps + 1)
    with torch.no_grad():
        for i in range(steps):
            t, s = float(ts[i]), float(ts[i + 1])
            masked = x == MASK_ID
            if not bool(masked.any()):
                break
            logits = model(x, c) / max(temperature, 1e-6)
            probs = torch.softmax(logits, dim=-1)
            draw = torch.multinomial(probs.reshape(-1, len(BASES)), 1,
                                     generator=g).reshape(x.shape)
            if s <= 0.0:
                reveal = masked
            else:
                p = (t - s) / max(t, 1e-9)
                reveal = masked & (torch.rand(x.shape, device=dev, generator=g) < p)
            x = torch.where(reveal, draw, x)
    if bool((x == MASK_ID).any()):
        raise RuntimeError("%d positions are still MASK after sampling; a leftover absorbing state "
                           "would be decoded as a base and scored as real sequence."
                           % int((x == MASK_ID).sum()))
    model.train()
    return ["".join(BASES[int(v)] for v in row) for row in x.cpu()]


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--shard", required=True,
                    help="train_<cell>_folds1-6.tsv. REQUIRED and explicit: there is no default "
                         "data path, because a defaulted one reports absent for data that is present.")
    ap.add_argument("--shard-md5", default=None,
                    help="expected md5 from the shard manifest. Omitting it skips the vintage check.")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--size", required=True, choices=sorted(SIZES))
    ap.add_argument("--seed", type=int, required=True)
    ap.add_argument("--cell", default="K562")
    ap.add_argument("--steps", type=int, default=20000,
                    help="SAME for every size, so the curve isolates capacity (default 20000)")
    ap.add_argument("--batch", type=int, default=64)
    ap.add_argument("--lr", type=float, default=3e-4)
    ap.add_argument("--dropout", type=float, default=0.1)
    ap.add_argument("--val-frac", type=float, default=0.02,
                    help="held-out slice of folds 1-6 for eval loss ONLY. Folds 9/10 are not on "
                         "this machine, so eval loss here is in-corpus and is NOT a design number.")
    ap.add_argument("--shuffle-labels", action="store_true",
                    help="THE NO-INFORMATION ARM. Permutes the bin column across rows, preserving "
                         "the marginal, destroying the association. Not a debug flag.")
    ap.add_argument("--sample-per-bin", type=int, default=240)
    ap.add_argument("--sample-steps", type=int, default=64)
    ap.add_argument("--ckpt-every", type=int, default=250)
    ap.add_argument("--log-every", type=int, default=100)
    ap.add_argument("--resume", default="auto", choices=("auto", "never"))
    ap.add_argument("--device", default="cuda", choices=("cuda", "cpu"))
    ap.add_argument("--on-pause", default="exit", choices=("exit", "wait"),
                    help="what to do when a PAUSE file appears in --out-dir. `exit` (default) "
                         "checkpoints and exits 0, releasing the GPU completely -- the right "
                         "behaviour on a shared gaming box. `wait` offloads to CPU and polls.")
    ap.add_argument("--gpu-sleep-ms", type=float, default=0.0,
                    help="sleep between steps to leave the GPU responsive while someone games")
    a = ap.parse_args(argv)

    try:
        import torch
    except ImportError:
        print("REFUSE: torch is not installed. See the install line in the handoff memo; it must "
              "be a cu126-or-older wheel because CUDA 12.8+ dropped Pascal (sm_61).",
              file=sys.stderr)
        return 2

    out_dir = Path(a.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    tag = "%s_%s_seed%d%s" % (a.cell, a.size, a.seed, "_shuffled" if a.shuffle_labels else "")
    ckpt_path = out_dir / ("ckpt_%s.pt" % tag)
    pause_file = out_dir / "PAUSE"

    signal.signal(signal.SIGINT, _on_signal)
    signal.signal(signal.SIGTERM, _on_signal)

    env = assert_pascal_ok(torch, a.device)
    dev = torch.device(a.device)
    torch.manual_seed(a.seed)

    bins, seqs, shard_md5 = load_shard(Path(a.shard), a.shard_md5)
    print("[data] %d rows from %s (md5 %s)" % (len(seqs), Path(a.shard).name, shard_md5), flush=True)

    if a.shuffle_labels:
        rng = random.Random(a.seed ^ 0x5EED)
        bins = bins[:]
        rng.shuffle(bins)
        print("[arm] NO-INFORMATION: bin column permuted (marginal preserved)", flush=True)

    lut = {c: i for i, c in enumerate(BASES)}
    x_all = torch.tensor([[lut[ch] for ch in s] for s in seqs], dtype=torch.long)
    c_all = torch.tensor(bins, dtype=torch.long)

    order = list(range(len(seqs)))
    random.Random(12345).shuffle(order)       # split seed FIXED across runs: same val slice always
    n_val = max(1, int(len(order) * a.val_frac))
    val_idx = torch.tensor(order[:n_val], dtype=torch.long)
    trn_idx = torch.tensor(order[n_val:], dtype=torch.long)
    x_trn, c_trn = x_all[trn_idx].to(dev), c_all[trn_idx].to(dev)
    x_val, c_val = x_all[val_idx][:512].to(dev), c_all[val_idx][:512].to(dev)
    print("[split] train %d / val %d (val seed fixed at 12345 across every run)"
          % (x_trn.shape[0], x_val.shape[0]), flush=True)

    cfg = dict(SIZES[a.size], dropout=a.dropout)
    model = build_model(torch, cfg, a.seed).to(dev)
    n_real = sum(p.numel() for p in model.parameters())
    print("[model] size=%s %s  params=%d (estimate %d)"
          % (a.size, cfg, n_real, n_params_estimate(cfg["d_model"], cfg["n_layers"], cfg["d_ff"])),
          flush=True)
    opt = torch.optim.AdamW(model.parameters(), lr=a.lr, weight_decay=0.01)

    step0, history = 0, []
    if a.resume == "auto" and ckpt_path.is_file():
        ck = torch.load(ckpt_path, map_location=dev, weights_only=False)
        if ck.get("cfg") != cfg or ck.get("shard_md5") != shard_md5:
            print("REFUSE: checkpoint %s was written under a different config or shard "
                  "(cfg match=%s, shard match=%s). Delete it or change --out-dir rather than "
                  "resuming a different experiment into this one."
                  % (ckpt_path.name, ck.get("cfg") == cfg, ck.get("shard_md5") == shard_md5),
                  file=sys.stderr)
            return 3
        model.load_state_dict(ck["model"])
        opt.load_state_dict(ck["opt"])
        step0, history = ck["step"], ck.get("history", [])
        print("[resume] from step %d" % step0, flush=True)

    def save(step: int) -> None:
        tmp = ckpt_path.with_suffix(".pt.tmp")
        torch.save({"model": model.state_dict(), "opt": opt.state_dict(), "step": step,
                    "cfg": cfg, "shard_md5": shard_md5, "history": history,
                    "args": vars(a)}, tmp)
        os.replace(tmp, ckpt_path)      # atomic: a torn checkpoint is worse than none

    def eval_loss() -> float:
        model.eval()
        with torch.no_grad():
            torch.manual_seed(999)      # fixed corruption draw, so eval loss is comparable per step
            v = float(loss_fn(torch, model, x_val, c_val))
        model.train()
        return v

    g_cpu = random.Random(a.seed)
    n_trn = x_trn.shape[0]
    model.train()
    t_start, t_mark, step = time.time(), time.time(), step0
    print("[train] %d -> %d steps, batch %d, fp32 (no AMP: sm_61 has no usable tensor cores)"
          % (step0, a.steps, a.batch), flush=True)

    while step < a.steps:
        if _STOP["flag"]:
            save(step)
            print("[yield] checkpointed at step %d (%s). Relaunch the identical command; "
                  "--resume auto picks up." % (step, _STOP["why"]), flush=True)
            return 0
        if pause_file.exists():
            save(step)
            if a.on_pause == "exit":
                print("[yield] PAUSE file present -- checkpointed at step %d and releasing the GPU. "
                      "Delete %s and relaunch the identical command to continue."
                      % (step, pause_file), flush=True)
                return 0
            print("[yield] PAUSE present -- offloading to CPU and polling every 15 s", flush=True)
            model.to("cpu")
            torch.cuda.empty_cache() if dev.type == "cuda" else None
            while pause_file.exists() and not _STOP["flag"]:
                time.sleep(15)
            model.to(dev)
            print("[resume] PAUSE cleared at step %d" % step, flush=True)
            continue

        idx = torch.tensor([g_cpu.randrange(n_trn) for _ in range(a.batch)],
                           dtype=torch.long, device=dev)
        loss = loss_fn(torch, model, x_trn[idx], c_trn[idx])
        opt.zero_grad(set_to_none=True)
        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
        opt.step()
        step += 1

        if a.gpu_sleep_ms > 0:
            time.sleep(a.gpu_sleep_ms / 1000.0)

        if step % a.log_every == 0:
            dt = time.time() - t_mark
            t_mark = time.time()
            sps = a.log_every / dt if dt > 0 else float("nan")
            ev = eval_loss()
            row = {"step": step, "train_loss": round(float(loss), 4), "eval_loss": round(ev, 4),
                   "steps_per_s": round(sps, 2)}
            history.append(row)
            vram = (torch.cuda.max_memory_allocated() / 1e9) if dev.type == "cuda" else 0.0
            print("  step %6d  train %8.4f  eval %8.4f  %6.2f step/s  peak %.2f GB  eta %s"
                  % (step, row["train_loss"], ev, sps, vram,
                     time.strftime("%H:%M:%S", time.gmtime((a.steps - step) / max(sps, 1e-9)))),
                  flush=True)
        if step % a.ckpt_every == 0:
            save(step)

    save(step)
    train_s = time.time() - t_start
    print("[train] done %d steps in %.1f s" % (step - step0, train_s), flush=True)

    targets = [b for b in range(N_BINS) for _ in range(a.sample_per_bin)]
    print("[sample] %d sequences (%d per bin x %d bins), %d reverse steps"
          % (len(targets), a.sample_per_bin, N_BINS, a.sample_steps), flush=True)
    gen = sample(torch, model, cfg, targets, steps=a.sample_steps, seed=a.seed)

    sfile = out_dir / ("samples_%s.tsv" % tag)
    with sfile.open("w", encoding="utf-8", newline="\n") as fh:
        fh.write("target_bin\tsequence\n")
        for b, s in zip(targets, gen):
            fh.write("%d\t%s\n" % (b, s))

    distinct = len(set(gen))
    man = {
        "record": "glm-human-diffusion-server1-run-v1",
        "tag": tag, "cell": a.cell, "size": a.size, "seed": a.seed,
        "arm": "no_information_shuffled_labels" if a.shuffle_labels else "conditioned",
        "scored_here": False,
        "why_not_scored_here": (
            "the judge is an oracle trained on folds 7-8 and the decision quantities are frozen in "
            "wiki/human_design_bar_2026-11-22.json. Scoring on this box with a re-derived oracle "
            "would make every run a different generator seen through a different judge."),
        "config": cfg, "n_params": n_real,
        "training": {"steps_target": a.steps, "steps_done": step, "batch": a.batch, "lr": a.lr,
                     "tokens_seen": step * a.batch * ELEMENT_LEN,
                     "train_seconds_this_process": round(train_s, 1),
                     "precision": "fp32", "amp": False,
                     "why_fp32": "sm_61 has no usable tensor cores; fp16 is not faster there"},
        "data": {"shard": Path(a.shard).name, "shard_md5": shard_md5,
                 "n_rows": len(seqs), "n_train": int(x_trn.shape[0]), "n_val": int(x_val.shape[0]),
                 "val_split_seed": 12345},
        "samples": {"file": sfile.name, "per_bin": a.sample_per_bin, "n_bins": N_BINS,
                    "total": len(gen), "sample_steps": a.sample_steps,
                    "n_distinct": distinct,
                    "distinct_fraction": round(distinct / len(gen), 4)},
        "final_losses": history[-1] if history else None,
        "history": history,
        "environment": env,
        "honest_limits": [
            "Eval loss is a held-out slice of folds 1-6, not fold 9. It is an in-corpus loss and "
            "is NOT a design number; the bar's decision quantities are monotonicity and lift, "
            "neither of which this script computes.",
            "CUDA and CPU RNG streams differ, so a seed reproduces on the SAME device only.",
            "distinct_fraction is a cheap degeneracy tripwire, not the copying control. The real "
            "copying check is longest shared substring against the training folds, run on the host.",
        ],
    }
    (out_dir / ("run_manifest_%s.json" % tag)).write_text(json.dumps(man, indent=2),
                                                          encoding="utf-8")
    print("-> %s" % sfile, flush=True)
    print("-> %s" % (out_dir / ("run_manifest_%s.json" % tag)), flush=True)
    print("[done] distinct %d/%d (%.3f)  final %s"
          % (distinct, len(gen), distinct / len(gen), man["final_losses"]), flush=True)
    if distinct < len(gen) * 0.5:
        print("[WARN] under half the samples are distinct -- degenerate output. Score it anyway, "
              "but read the copying and FLOOR controls before reading any monotonicity.", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
