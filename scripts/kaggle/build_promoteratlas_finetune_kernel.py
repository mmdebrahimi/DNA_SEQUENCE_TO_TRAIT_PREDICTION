#!/usr/bin/env python
"""Build the Kaggle T4 kernel for the `finetuned-promoteratlas-encoder-150bp` arm.

WHY THIS KERNEL EXISTS AT ALL. The arm's frozen bar routed it to `local_cpu`, and that choice was
never about compute: its stated decisive ground was retrievability -- *"a Kaggle chain requires an
invocation to retrieve its output and nothing will invoke me inside this window"* -- where the
window was a 216.0 h authorised unattended absence. That absence ended, so the route was re-decided
to `kaggle_t4` on explicit user authority (2026-11-14), with a diff proving the change moved no
predicate: `check_order`, the FINETUNE_ADDS/TIES/NO_GAIN decision, every threshold and all of
`frozen_hyperparameters` are byte-identical. Re-routing also removes the failure mode that already
destroyed the arm's only launch -- this host recorded THREE unexpected shutdowns in ~42 h and the
launch died in the third (`wiki/finetune_arm_launch_death_2026-11-14.json`).

IT CALLS `run()`, NOT `main()`, AND THAT IS THE DESIGN. `finetune_arm_run.main()` refuses any route
but `local_cpu`. That refusal is CORRECT and is deliberately left alone: it is a guard on the LOCAL
launch surface, and local is no longer the chosen route. The route gate is not in `run()`, so the
kernel calls `run()` with the same bar, the same tiles and the same trajectory semantics. Nothing
is patched, nothing is forked.

WHY THE SOURCES ARE EMBEDDED RATHER THAN CLONED. `fragment_pyramid_kernel.py` argues, correctly,
that a clone at a pinned SHA beats an embedded tarball -- but that needs the commit to EXIST on the
remote. It does not: this working tree is 60+ commits ahead of origin and nothing has been pushed
this whole arc, including the bar this kernel implements. So this follows the other documented
pattern (`build_glm_sampling_sweep_kernel.py`, `build_promoteratlas_tile_arm_kernel.py`) and embeds
sources to stay version-exact against a tree newer than what is published.

THE EMBED SET IS RESOLVED BY AST, NOT HAND-LISTED. A regex sweep for imports over these modules
returns `from the` -- a docstring line -- and misses nothing only by luck. A hand-listed set is also
the exact shape of this repo's recorded `hardcoded_exclusion_list_undercovers` failure: it drifts
and fails by PASSING. `closure()` walks real `ast.Import`/`ast.ImportFrom` nodes from the entry
modules and keeps whatever resolves inside the repo, so a new import is carried automatically.

THE PROMOTERATLAS CLONE AND THE COHORT NEED NO CODE CHANGES. `atlas_embed.DEFAULT_PA_REPO` reads
`$PA_REPO` and `positional_arm_probe.SUBSTRATE_CACHE` reads `$MPRA_CACHE`, so the kernel points both
at Kaggle paths by environment alone. The PromoterAtlas repo is public and MIT, so cloning it is a
read of someone else's published artifact, not a push of ours.

ORDER IS SMOKE THEN FULL, AND THE PER-EPOCH COST IS MEASURED. The local probe's own lesson was that
a guessed rate under-sizes a run: its figure counted training only and omitted the validation pass,
a measured +28%. So the kernel TIMES one real epoch on the actual device and sizes the horizon from
that, rather than carrying the local 1568.1 s/epoch onto a T4 where it is meaningless. The smoke
epoch runs in a THROWAWAY checkpoint directory -- see the hazard note below.

A WALL-CLOCK STOP STILL YIELDS A VERDICT, which is what makes a 12 h session cap survivable.
`run()` orders its phases budget -> nulls -> extend precisely because the bar makes the shuffled
null a GATE: a run cut at the ceiling after training everything to the frozen 7-epoch budget and
completing the nulls is gradeable, and only the ineligible 8-to-38 diagnostic is lost. An early stop
stamps BOUNDED_BY_WALL_CLOCK and never COMPLETE.

HAZARD, RECORDED BECAUSE IT IS SILENT. Do NOT point the smoke run at the real checkpoint directory.
The resume contract READS any existing per-seed checkpoint and never recomputes it ("Already at or
past this depth: READ, never recomputed"), so a one-epoch smoke under smoke hyperparameters would be
silently adopted as the real run's trajectory prefix. The kernel uses a separate throwaway dir and
asserts the real one is untouched by the smoke.

THE BUILDER VALIDATES THE EMITTED KERNEL, not just itself. On 2026-10-31 a builder that checked only
its own syntax shipped an unterminated string literal to Kaggle twice, because a collapsed escape
broke the EMITTED body. `ast.parse(body)` runs before the write and `--self-test` proves that check
is non-vacuous.

PUSH (free tier -> no money gate; T4 pinned because `enable_gpu` alone can provision a P100):
    uv run python scripts/kaggle/build_promoteratlas_finetune_kernel.py --out /tmp/ft_kernel.py
    uv run python scripts/kaggle_push_poll.py push /tmp/ft_kernel.py pa-finetune-arm-v1 --gpu
    uv run python scripts/kaggle_push_poll.py status pa-finetune-arm-v1
    uv run python scripts/kaggle_push_poll.py pull   pa-finetune-arm-v1 ./kaggle_out

CHAINING ACROSS SESSIONS, if one 12 h session is not enough. Pull the output, upload the checkpoint
directory as a dataset, and attach it on the next push; the kernel seeds its real checkpoint dir
from `$FT_CKPT_IN` (default `/kaggle/input`) before running, so finished seeds are READ and never
recomputed. The resume property is by construction -- each epoch's shuffle is seeded by
`(seed, epoch)`, never by the accumulated global RNG -- so a resumed run reproduces an
uninterrupted one exactly.
"""
from __future__ import annotations

import argparse
import ast
import base64
import gzip
import io
import json
import tarfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]

#: Entry modules. `finetune_arm_run` is the trainer; `positional_arm_probe` supplies
#: `load_arm_cohort`; `genomewide` supplies the downloader and `spearman`. Their transitive
#: in-repo imports are resolved from these by AST.
ENTRIES = (
    "scripts.finetune_arm_run",
    "scripts.finetune_arm_probe",
    "scripts.positional_arm_probe",
    "dna_decode.glm.genomewide",
)

BAR_REL = "wiki/finetune_arm_acceptance_bar.json"
PA_URL = "https://github.com/LucasCoppens/PromoterAtlas.git"


# ---------------------------------------------------------------------------------------------------
# the embed set, resolved by AST
# ---------------------------------------------------------------------------------------------------
def _module_path(mod: str) -> Path | None:
    """Resolve a dotted module name to a file inside this repo, or None if it is external."""
    p = REPO / (mod.replace(".", "/") + ".py")
    if p.is_file():
        return p
    pkg = REPO / mod.replace(".", "/") / "__init__.py"
    return pkg if pkg.is_file() else None


def _imports_of(path: Path) -> set[str]:
    """Dotted module names imported by `path`, read from the AST.

    A docstring that happens to begin 'from the ...' is an `ast.Expr`, never an import, so it
    cannot land here -- which is the whole reason this is not a regex.
    """
    found: set[str] = set()
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for a in node.names:
                found.add(a.name)
        elif isinstance(node, ast.ImportFrom):
            if node.level:                      # relative import inside a package
                parts = path.relative_to(REPO).with_suffix("").parts[:-1]
                base = ".".join(parts[: len(parts) - node.level + 1])
                found.add(f"{base}.{node.module}" if node.module else base)
            elif node.module:
                found.add(node.module)
                # `from pkg.mod import name` may name a SUBMODULE, not an attribute.
                for a in node.names:
                    found.add(f"{node.module}.{a.name}")
    return found


def closure(entries=ENTRIES) -> dict[str, Path]:
    """Every in-repo module reachable from `entries`, keyed by dotted name."""
    out: dict[str, Path] = {}
    queue = list(entries)
    while queue:
        mod = queue.pop()
        if mod in out:
            continue
        p = _module_path(mod)
        if p is None:
            continue                            # external (stdlib / numpy / torch) -- not embedded
        out[mod] = p
        queue.extend(_imports_of(p))
    return out


def _package_inits(mods: dict[str, Path]) -> dict[str, Path]:
    """`__init__.py` for every package that holds an embedded module, so imports resolve."""
    extra: dict[str, Path] = {}
    for mod in mods:
        parts = mod.split(".")
        for i in range(1, len(parts)):
            pkg = ".".join(parts[:i])
            ini = REPO / pkg.replace(".", "/") / "__init__.py"
            if ini.is_file() and pkg not in mods:
                extra[pkg] = ini
    return extra


def build_payload() -> tuple[str, dict[str, Path], str]:
    """(base64 tar.gz of the sources, the embed set, base64 gzip of the bar)."""
    mods = closure()
    mods.update(_package_inits(mods))
    if not mods:
        raise SystemExit("REFUSED: the import closure is empty; the resolver is broken")

    buf = io.BytesIO()
    # Deterministic tar: sorted members, zeroed mtime, fixed uid/gid.
    with tarfile.open(fileobj=buf, mode="w") as tf:
        for mod in sorted(mods):
            src = mods[mod]
            rel = src.relative_to(REPO).as_posix()
            data = src.read_bytes()
            ti = tarfile.TarInfo(rel)
            ti.size = len(data)
            ti.mtime = 0
            ti.uid = ti.gid = 0
            ti.uname = ti.gname = ""
            tf.addfile(ti, io.BytesIO(data))
    sources_b64 = base64.b64encode(gzip.compress(buf.getvalue(), mtime=0)).decode("ascii")

    bar_path = REPO / BAR_REL
    bar = json.loads(bar_path.read_text(encoding="utf-8"))
    if bar["route"]["chosen"] != "kaggle_t4":
        raise SystemExit(
            "REFUSED: the frozen bar routes this arm to %r, not kaggle_t4. Building a Kaggle "
            "kernel for a non-Kaggle route would ship a kernel whose own bar disowns it."
            % bar["route"]["chosen"]
        )
    bar_b64 = base64.b64encode(
        gzip.compress(bar_path.read_bytes(), mtime=0)
    ).decode("ascii")
    return sources_b64, mods, bar_b64


# ---------------------------------------------------------------------------------------------------
# the emitted kernel
# ---------------------------------------------------------------------------------------------------
KERNEL = '''#!/usr/bin/env python
"""GENERATED by scripts/kaggle/build_promoteratlas_finetune_kernel.py. Do not hand-edit.

The `finetuned-promoteratlas-encoder-150bp` arm, on a Kaggle T4. Calls the repo's own
`finetune_arm_run.run()` with the frozen bar -- not `main()`, whose route gate guards the LOCAL
launch surface and correctly refuses a non-local route.
"""
import base64
import gzip
import io
import json
import os
import shutil
import subprocess
import sys
import tarfile
import time
from pathlib import Path

WORK = Path(os.environ.get("FT_WORK", "/kaggle/working"))
SRC = WORK / "src"
CKPT = WORK / "finetune_checkpoints"
SMOKE_CKPT = WORK / "smoke_checkpoints_THROWAWAY"
CKPT_IN = Path(os.environ.get("FT_CKPT_IN", "/kaggle/input"))

#: Kaggle's GPU session cap is 12 h. Leave margin for clone + download + write so the artifact is
#: always flushed; `run()` stamps BOUNDED_BY_WALL_CLOCK if it stops at this, never COMPLETE.
WALL_HOURS = float(os.environ.get("FT_WALL_HOURS", "10.5"))
SMOKE_ONLY = os.environ.get("FT_SMOKE_ONLY", "") == "1"

PA_URL = {pa_url!r}
SOURCES_B64 = "{sources_b64}"
BAR_B64 = "{bar_b64}"
EMBEDDED = {embedded!r}


def sh(cmd, cwd=None, check=True):
    print("$ " + " ".join(cmd), flush=True)
    return subprocess.run(cmd, cwd=cwd, check=check)


def device_report():
    import torch
    print("torch %s  cuda=%s" % (torch.__version__, torch.cuda.is_available()), flush=True)
    if torch.cuda.is_available():
        name = torch.cuda.get_device_name(0)
        cap = torch.cuda.get_device_capability(0)
        print("device: %s  capability sm_%d%d" % (name, cap[0], cap[1]), flush=True)
        if cap[0] < 7:
            print("WARNING: compute capability < 7.0 (a P100 was provisioned, not a T4). This arm "
                  "is fp32 either way, so the run is slow-but-correct.", flush=True)
    else:
        print("WARNING: no CUDA device. The whole point of this route is the GPU; a CPU session "
              "will not finish inside the session cap.", flush=True)


def unpack():
    SRC.mkdir(parents=True, exist_ok=True)
    raw = gzip.decompress(base64.b64decode(SOURCES_B64))
    with tarfile.open(fileobj=io.BytesIO(raw), mode="r") as tf:
        tf.extractall(SRC)
    got = sorted(p.relative_to(SRC).as_posix() for p in SRC.rglob("*.py"))
    print("embedded %d module(s): %s" % (len(got), ", ".join(got)), flush=True)
    missing = [m for m in EMBEDDED if not (SRC / EMBEDDED[m]).is_file()]
    if missing:
        raise SystemExit("REFUSED: embedded payload is missing %s" % missing)
    sys.path.insert(0, str(SRC))

    bar_path = SRC / "wiki/finetune_arm_acceptance_bar.json"
    bar_path.parent.mkdir(parents=True, exist_ok=True)
    bar_path.write_bytes(gzip.decompress(base64.b64decode(BAR_B64)))
    return json.loads(bar_path.read_text(encoding="utf-8"))


def seed_checkpoints():
    """Copy any attached prior checkpoints in, so finished seeds are READ not recomputed."""
    CKPT.mkdir(parents=True, exist_ok=True)
    if not CKPT_IN.is_dir():
        return 0
    n = 0
    for p in CKPT_IN.rglob("*"):
        if p.is_file() and p.parent.name == "finetune_checkpoints":
            shutil.copyfile(p, CKPT / p.name)
            n += 1
    if n:
        print("seeded %d prior checkpoint file(s) from %s" % (n, CKPT_IN), flush=True)
    return n


def load_cohort():
    """Download the GEO cohort, then load it through the repo's own refusing loader.

    `load_arm_cohort` calls `load_peak_tiles(..., allow_download=False)` on purpose, so the fetch
    happens FIRST and separately. The bounded retry mirrors a measured 2026-11-01 failure: a run
    died with ContentTooShortError "got only 0 out of 9604523 bytes" while the run before it
    downloaded the same file fine -- a transient NCBI blip, not a code defect.
    """
    cache = Path(os.environ["MPRA_CACHE"])
    cache.mkdir(parents=True, exist_ok=True)
    from dna_decode.glm import genomewide as gw
    for attempt in range(1, 5):
        try:
            gw.load_peak_tiles(cache, allow_download=True)
            break
        except Exception as exc:
            print("cohort fetch %d/4 failed: %s: %s" % (attempt, type(exc).__name__, exc),
                  flush=True)
            if attempt == 4:
                raise
            for junk in list(cache.glob("*.tmp")) + list(cache.glob("*.part")):
                try:
                    junk.unlink()
                    print("   removed partial %s" % junk.name, flush=True)
                except OSError:
                    pass
            time.sleep(20 * attempt)

    from scripts.positional_arm_probe import load_arm_cohort
    tiles, st = load_arm_cohort(cache)
    print("COHORT   %d tiles x %d bp from %s"
          % (len(tiles), len(tiles[0].seq), getattr(st, "path", "?")), flush=True)
    return tiles


def measure_seconds_per_epoch(bar, tiles):
    """Time ONE real epoch on this device, in a THROWAWAY checkpoint dir.

    Never the real dir: the resume contract reads any existing per-seed checkpoint and never
    recomputes it, so a smoke epoch written there would become the real run's trajectory prefix.
    """
    from scripts.finetune_arm_run import train_one_seed
    hp = bar["frozen_hyperparameters"]
    if SMOKE_CKPT.exists():
        shutil.rmtree(SMOKE_CKPT)
    SMOKE_CKPT.mkdir(parents=True, exist_ok=True)
    before = sorted(p.name for p in CKPT.glob("*")) if CKPT.is_dir() else []

    t0 = time.time()
    train_one_seed(tiles, seed=0, horizon_epochs=1, lr=hp["learning_rate"],
                   batch_size=hp["batch_size"], split_frac=bar["comparator"]["split_frac"],
                   ckpt_dir=SMOKE_CKPT, build=None, kind="smoke", verbose=True)
    spe = time.time() - t0

    after = sorted(p.name for p in CKPT.glob("*")) if CKPT.is_dir() else []
    if before != after:
        raise SystemExit("REFUSED: the smoke run touched the REAL checkpoint dir (%s -> %s). "
                         "That would have become the real trajectory's prefix." % (before, after))
    shutil.rmtree(SMOKE_CKPT, ignore_errors=True)
    print("MEASURED %.1f s/epoch on this device (one real epoch, throwaway ckpt)" % spe, flush=True)
    return spe


def main():
    device_report()
    WORK.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("PA_REPO", str(WORK / "PromoterAtlas"))
    os.environ.setdefault("MPRA_CACHE", str(WORK / "mpra"))

    bar = unpack()
    print("BAR route=%s  horizon=%s  seeds=%s  budget=%s epochs"
          % (bar["route"]["chosen"], bar["horizon"]["diagnostic_horizon_epochs"],
             bar["horizon"]["n_seeds"], bar["frozen_hyperparameters"]["n_epochs_budget"]),
          flush=True)
    if bar["route"]["chosen"] != "kaggle_t4":
        raise SystemExit("REFUSED: embedded bar routes to %r" % bar["route"]["chosen"])

    pa = Path(os.environ["PA_REPO"])
    if not pa.exists():
        sh(["git", "clone", "--depth", "1", PA_URL, str(pa)])
    sh(["git", "rev-parse", "HEAD"], cwd=str(pa), check=False)

    tiles = load_cohort()
    seed_checkpoints()
    spe = measure_seconds_per_epoch(bar, tiles)

    n_seeds = bar["horizon"]["n_seeds"]
    horizon = bar["horizon"]["diagnostic_horizon_epochs"]
    print("PROJECTION %.2f h for %d seeds x %d epochs at the measured rate; wall ceiling %.2f h"
          % (spe * n_seeds * horizon / 3600.0, n_seeds, horizon, WALL_HOURS), flush=True)
    if SMOKE_ONLY:
        print("FT_SMOKE_ONLY=1 -- stopping after the measurement. Nothing was scored.", flush=True)
        (WORK / "finetune_arm_smoke.json").write_text(
            json.dumps({{"schema": "finetune-arm-kaggle-smoke-v1",
                         "seconds_per_epoch_measured": spe,
                         "projected_hours_at_horizon": spe * n_seeds * horizon / 3600.0,
                         "wall_ceiling_hours": WALL_HOURS, "scored": False}}, indent=2),
            encoding="utf-8")
        return 0

    from scripts.finetune_arm_run import run
    rec = run(bar=bar, tiles=tiles, ckpt_dir=CKPT, build=None,
              wall_ceiling_hours=WALL_HOURS, seconds_per_epoch=spe, verbose=True)

    out = WORK / "finetune_arm_run_kaggle.json"
    out.write_text(json.dumps(rec, indent=2, sort_keys=True), encoding="utf-8")
    print()
    print("STATUS   %s" % rec.get("status"), flush=True)
    print("         %s" % rec.get("why_status"), flush=True)
    print("WROTE    %s" % out, flush=True)
    if rec.get("bounded_by_wall_clock"):
        print("BOUNDED BY WALL CLOCK -- not COMPLETE. The frozen-budget and null phases run FIRST, "
              "so a verdict may still be gradeable; the 8-to-38 diagnostic is what is short. "
              "Re-push with the checkpoints attached to extend.", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
'''


def emit(sources_b64: str, mods: dict[str, Path], bar_b64: str) -> str:
    embedded = {m: p.relative_to(REPO).as_posix() for m, p in sorted(mods.items())}
    body = KERNEL.format(pa_url=PA_URL, sources_b64=sources_b64, bar_b64=bar_b64,
                         embedded=embedded)
    # Validate the EMITTED body, not just this builder. A collapsed escape has shipped twice.
    ast.parse(body)
    return body


def self_test() -> int:
    """Prove the closure is real, the payload round-trips, and the ast check is non-vacuous."""
    mods = closure()
    assert "scripts.finetune_arm_run" in mods, "entry module missing from its own closure"
    # The closure must find modules the entries import only TRANSITIVELY.
    assert "scripts.positional_arm_score" in mods, "transitive import not resolved"
    assert "dna_decode.glm.genomewide" in mods, "genomewide not resolved"
    # External modules must NOT be embedded.
    for ext in ("numpy", "torch", "json", "pathlib"):
        assert ext not in mods, "external module %r was embedded" % ext
    # A docstring opening 'from the' must not be mistaken for an import.
    assert not any(m.startswith("the") for m in mods), "a prose line was parsed as an import"

    sources_b64, mods2, bar_b64 = build_payload()
    raw = gzip.decompress(base64.b64decode(sources_b64))
    with tarfile.open(fileobj=io.BytesIO(raw), mode="r") as tf:
        names = sorted(tf.getnames())
    for m, p in mods2.items():
        rel = p.relative_to(REPO).as_posix()
        assert rel in names, "%s (%s) missing from the tar" % (m, rel)

    body = emit(sources_b64, mods2, bar_b64)
    ast.parse(body)

    # Non-vacuity: a broken body MUST fail the same check that passed above.
    broken = body.replace('SOURCES_B64 = "', 'SOURCES_B64 = "\n', 1)
    try:
        ast.parse(broken)
    except SyntaxError:
        pass
    else:
        raise AssertionError("ast.parse accepted a deliberately broken body; the check is vacuous")

    print("self-test OK")
    print("  closure: %d module(s)" % len(mods2))
    for m in sorted(mods2):
        print("    %s -> %s" % (m, mods2[m].relative_to(REPO).as_posix()))
    print("  tar members: %d | kernel bytes: %d" % (len(names), len(body)))
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", help="write the kernel here")
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args(argv)

    if args.self_test:
        return self_test()
    if not args.out:
        ap.error("--out is required (or pass --self-test)")

    sources_b64, mods, bar_b64 = build_payload()
    body = emit(sources_b64, mods, bar_b64)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(body, encoding="utf-8")
    print("WROTE %s (%d bytes)" % (out, len(body)))
    print("embedded %d module(s) + the frozen bar" % len(mods))
    for m in sorted(mods):
        print("  %s" % mods[m].relative_to(REPO).as_posix())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
