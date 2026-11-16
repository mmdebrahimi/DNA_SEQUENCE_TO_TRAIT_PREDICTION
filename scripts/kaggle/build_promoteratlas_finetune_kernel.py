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

A KAGGLE SCRIPT KERNEL RECEIVES NO ENVIRONMENT VARIABLES, so the knobs are BAKED at build time.
The metadata `scripts/kaggle_push_poll.py` writes carries `id`, `title`, `code_file`, `language`,
`kernel_type`, `is_private`, `enable_gpu`, `enable_internet`, `machine_shape` and an optional
`dataset_sources` -- there is no env field anywhere in it. An `os.environ.get` in the emitted kernel
can therefore only ever return its own default on Kaggle, which left smoke-first sizing and
cross-session chaining unreachable. `--smoke`, `--wall-hours`, `--ckpt-dataset` and
`--horizon-epochs` bake literals that WIN over the env reads; the env reads stay live underneath so
a local dry-run can still override.

A BOUNDED CHAIN, which is what `--horizon-epochs` exists for. At the T4-measured 1276.04 s/epoch the
bar's full 10-seed x 38-epoch horizon projects to 134.69 h, and even the gradeable prefix -- the
frozen 7-epoch budget across 10 seeds plus the 21 null epochs -- costs 32.3 h. No session is that
long, and `plan_run` REFUSES whenever the ceiling leaves fewer than `frozen_epochs` per seed. That
refusal is right for ONE session and wrong for a chain: `run()` calls `plan_run` only when
`horizon_epochs is None`, so supplying it hands planning to the caller and records
"horizon supplied by the caller, not planned from the ceiling". The session then self-bounds on
`advance()`'s `time.time() >= deadline` check, so it writes its artifact gracefully instead of being
killed, and -- because `_write_curve` fires INSIDE the per-epoch loop -- leaves every completed epoch
on disk for the next session to READ. The default stays unset so an un-flagged build still refuses.

SMOKE FIRST, THEN SIZE. The local 1568.1 s/epoch was measured on 4 CPU threads and does not transfer
to a T4 in either direction, so the horizon is sized from a measurement on the device that actually
arrives -- the same discipline that rejected a guessed local rate:
    uv run python scripts/kaggle/build_promoteratlas_finetune_kernel.py --smoke --out /tmp/ft_smoke.py
    uv run python scripts/kaggle_push_poll.py push /tmp/ft_smoke.py pa-finetune-arm-v1 --gpu
    uv run python scripts/kaggle_push_poll.py status pa-finetune-arm-v1
    uv run python scripts/kaggle_push_poll.py pull   pa-finetune-arm-v1 ./kaggle_out

THE FULL RUN (free tier -> no money gate; T4 pinned because `enable_gpu` alone can provision a P100):
    uv run python scripts/kaggle/build_promoteratlas_finetune_kernel.py --out /tmp/ft_kernel.py
    uv run python scripts/kaggle_push_poll.py push /tmp/ft_kernel.py pa-finetune-arm-v1 --gpu

CHAINING ACROSS SESSIONS, if one session is not enough. Pull the output, upload the checkpoint
directory as a dataset, and attach it on the next push with BOTH flags -- `--ckpt-dataset <slug>` so
the kernel knows where to look, and `kaggle_push_poll.py push --dataset <slug>` so Kaggle mounts it:
    ... --ckpt-dataset ft-ckpt-round1 --out /tmp/ft_chain.py
    uv run python scripts/kaggle_push_poll.py push /tmp/ft_chain.py pa-finetune-arm-v1 --gpu \
        --dataset ft-ckpt-round1
The kernel seeds its real checkpoint dir from that mount before running, so finished seeds are READ
and never recomputed. The resume property is by construction -- each epoch's shuffle is seeded by
`(seed, epoch)`, never by the accumulated global RNG -- so a resumed run reproduces an uninterrupted
one exactly.
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

#: BAKED AT BUILD TIME, because a Kaggle SCRIPT KERNEL RECEIVES NO ENVIRONMENT VARIABLES. The
#: metadata scripts/kaggle_push_poll.py writes carries id, title, code_file, language, kernel_type,
#: is_private, enable_gpu, enable_internet, machine_shape and an optional dataset_sources -- there is
#: no env field anywhere in it, so on Kaggle an os.environ read can only ever return its own default.
#: A baked literal WINS; the env read survives below it only so a LOCAL dry-run can still override.
#: None means "not baked", and that is what keeps the env path live rather than dead code.
BAKED_SMOKE_ONLY = {baked_smoke_only!r}
BAKED_WALL_HOURS = {baked_wall_hours!r}
BAKED_CKPT_IN = {baked_ckpt_in!r}
BAKED_HORIZON_EPOCHS = {baked_horizon_epochs!r}

WORK = Path(os.environ.get("FT_WORK", "/kaggle/working"))
SRC = WORK / "src"
CKPT = WORK / "finetune_checkpoints"
SMOKE_CKPT = WORK / "smoke_checkpoints_THROWAWAY"
CKPT_IN = (Path(BAKED_CKPT_IN) if BAKED_CKPT_IN is not None
           else Path(os.environ.get("FT_CKPT_IN", "/kaggle/input")))

#: The root Kaggle mounts attached datasets somewhere under. Used as the FALLBACK search root when a
#: baked mount path does not exist, which is the measured case rather than the hypothetical one.
INPUT_ROOT = Path(os.environ.get("FT_INPUT_ROOT", "/kaggle/input"))


#: 10.5 h is DERIVED, not picked. Kaggle's batch GPU session cap is 12.0 h, and this kernel spends
#: wall clock on either side of training: the shallow PromoterAtlas clone, the GSE144621 GEO download
#: (bounded 4-attempt retry, after a measured ContentTooShortError), the cohort build, and the final
#: artifact write. 12.0 - 10.5 = 1.5 h of margin pays for those, so the artifact is always flushed
#: inside the session. `run()` stamps BOUNDED_BY_WALL_CLOCK when it stops at this ceiling, never
#: COMPLETE -- the frozen-budget and null phases run FIRST, so such a stop is still gradeable.
WALL_HOURS = (float(BAKED_WALL_HOURS) if BAKED_WALL_HOURS is not None
              else float(os.environ.get("FT_WALL_HOURS", "10.5")))
SMOKE_ONLY = (bool(BAKED_SMOKE_ONLY) if BAKED_SMOKE_ONLY is not None
              else os.environ.get("FT_SMOKE_ONLY", "") == "1")

#: None means "let `run()` plan the horizon from the ceiling". That is the DEFAULT and it is the
#: right default: at the measured T4 rate `plan_run` REFUSES every session-sized ceiling, and that
#: refusal is correct for a single unchained session. A baked integer instead BYPASSES planning --
#: `run()` takes the horizon as given and `wall_ceiling_hours` alone decides when the session stops,
#: which is what makes a chained sequence of self-bounding sessions possible. The bypass is recorded
#: in the run artifact under `horizon_divergence` rather than left for a reader to infer.
HORIZON_EPOCHS = (int(BAKED_HORIZON_EPOCHS) if BAKED_HORIZON_EPOCHS is not None
                  else (int(os.environ["FT_HORIZON_EPOCHS"])
                        if os.environ.get("FT_HORIZON_EPOCHS") else None))

PA_URL = {pa_url!r}
SOURCES_B64 = "{sources_b64}"
BAR_B64 = "{bar_b64}"
EMBEDDED = {embedded!r}


def sh(cmd, cwd=None, check=True):
    print("$ " + " ".join(cmd), flush=True)
    return subprocess.run(cmd, cwd=cwd, check=check)


def device_report():
    """Print the device AND return it, so a measured rate is never reported without its device.

    A seconds-per-epoch figure is a property of the hardware that produced it. The T4 pin can still
    yield a P100 (CC 6.0), and the arm is fp32 either way, so a P100 is slow-but-correct -- the risk
    is a rate attributed to the wrong device and then used to size the horizon.
    """
    import torch
    rep = {{"torch": torch.__version__, "cuda_available": bool(torch.cuda.is_available()),
           "device_name": None, "compute_capability": None, "t4_pin_honoured": None}}
    print("torch %s  cuda=%s" % (torch.__version__, torch.cuda.is_available()), flush=True)
    if torch.cuda.is_available():
        name = torch.cuda.get_device_name(0)
        cap = torch.cuda.get_device_capability(0)
        rep["device_name"] = name
        rep["compute_capability"] = "sm_%d%d" % (cap[0], cap[1])
        rep["t4_pin_honoured"] = bool(cap[0] >= 7)
        print("device: %s  capability sm_%d%d" % (name, cap[0], cap[1]), flush=True)
        if cap[0] < 7:
            print("WARNING: compute capability < 7.0 (a P100 was provisioned, not a T4). This arm "
                  "is fp32 either way, so the run is slow-but-correct.", flush=True)
    else:
        print("WARNING: no CUDA device. The whole point of this route is the GPU; a CPU session "
              "will not finish inside the session cap.", flush=True)
    return rep


#: Every bar path this kernel dereferences, checked UP FRONT. A missing key used to surface as a
#: bare KeyError from inside a call chain 42 s in, after a 25 s clone and the cohort build; now it
#: costs ~17 s and names the path. Keep this list in step with the dereferences below.
REQUIRED_BAR_PATHS = [
    ("route", "chosen"),
    ("horizon", "diagnostic_horizon_epochs"),
    ("horizon", "n_seeds"),
    ("frozen_hyperparameters", "n_epochs_budget"),
    ("frozen_hyperparameters", "learning_rate"),
    ("frozen_hyperparameters", "batch_size"),
    ("inherited_not_chosen", "comparator", "split_frac"),
]


def require_bar_paths(bar):
    """Refuse on a bar whose SHAPE this kernel cannot read, naming the path rather than KeyError-ing.

    This is the guard the first Kaggle run earned: `bar["comparator"]["split_frac"]` was one level
    too shallow, every local check was static, and the cost of finding out was a whole session.
    """
    missing = []
    for path in REQUIRED_BAR_PATHS:
        node = bar
        for key in path:
            if not isinstance(node, dict) or key not in node:
                missing.append(".".join(path))
                break
            node = node[key]
    if missing:
        raise SystemExit("REFUSED: the embedded bar is missing %d required path(s): %s. "
                         "The kernel dereferences these and would otherwise fail mid-run."
                         % (len(missing), ", ".join(missing)))
    print("BAR SHAPE ok -- all %d required path(s) present" % len(REQUIRED_BAR_PATHS), flush=True)


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
    """Copy any attached prior checkpoints in, so finished seeds are READ not recomputed.

    A zero count means two OPPOSITE things and the original `if n:` print collapsed them. With no
    --ckpt-dataset, nothing is attached and zero is CORRECT, so silence is right. With a dataset
    baked, zero is a CONTRADICTION -- the point of the flag is that finished seeds are read rather
    than recomputed, and a chained session that seeds nothing silently recomputes every seed while
    looking entirely normal. So the baked case prints its count UNCONDITIONALLY, a positive signal
    rather than the absence of a complaint, and REFUSES on zero.

    SELECTION IS BY FILENAME AND THE SEARCH ROOT IS NOT PREDICTED. Both halves are measured
    corrections, both made on 2026-11-15, and together they are why this function no longer asserts
    anything about Kaggle's layout:

      * The directory level inside the dataset does not survive. `kaggle datasets create` skips
        directories unless --dir-mode is given; --dir-mode zip uploads `finetune_checkpoints.zip`;
        Kaggle then extracts that archive and FLATTENS it. A dataset built from
        `finetune_checkpoints/seed_99_curve.json` stores exactly `seed_99_curve.json` at its root.
        So a rule keyed on `p.parent.name` could never seed a correctly-built dataset.
      * The mount path is not /kaggle/input/<slug>. A dataset attached via `dataset_sources` and
        reporting status=ready did not appear there at all; /kaggle/input held only `datasets/`.

    Hence: the baked path is a HINT, and when it is absent the search falls back to /kaggle/input and
    recurses. The filenames identify a checkpoint; the path never did. The fallback widens WHERE we
    look and never WHETHER a zero is a fault -- the refusal below still fires.
    """
    # The SINGLE definition of which filenames are real prior work, imported from the module that
    # PRODUCES those names. Step 7's uploader imports the same object, so the two ends of the round
    # trip cannot disagree about what counts as a checkpoint.
    from scripts.finetune_arm_run import CKPT_NAME_RE

    CKPT.mkdir(parents=True, exist_ok=True)

    root = None
    if CKPT_IN.is_dir():
        root = CKPT_IN
    elif BAKED_CKPT_IN is not None and INPUT_ROOT.is_dir():
        root = INPUT_ROOT
        print("the baked mount %s does not exist; searching %s recursively instead, because "
              "Kaggle's mount layout is not ours to predict" % (CKPT_IN, INPUT_ROOT), flush=True)

    n = 0
    taken = {{}}
    where = set()
    if root is not None:
        for p in sorted(root.rglob("*")):
            if not p.is_file() or not CKPT_NAME_RE.match(p.name):
                continue
            if p.name in taken:
                raise SystemExit(
                    "REFUSING to run: two different files under the checkpoint search root are "
                    "both named %r (%s and %s). The copy into the real checkpoint directory is "
                    "FLAT, so this would silently keep one curve for that seed and discard the "
                    "other." % (p.name, taken[p.name], p))
            taken[p.name] = p
            where.add(str(p.parent))
            shutil.copyfile(p, CKPT / p.name)
            n += 1

    if BAKED_CKPT_IN is None:
        # Nothing was attached and nothing was asked for; zero is the right answer. Stay quiet.
        if n:
            print("seeded %d prior checkpoint file(s) from %s" % (n, CKPT_IN), flush=True)
        return n

    print("seeded %d prior checkpoint file(s) for the BAKED dataset %s (searched %s); found in: %s"
          % (n, BAKED_CKPT_IN, root if root is not None else "<nothing searchable>",
             ", ".join(sorted(where)) if where else "<none>"), flush=True)
    if n == 0:
        found = []
        probe = root if root is not None else INPUT_ROOT
        if probe.is_dir():
            for p in sorted(probe.rglob("*"))[:40]:
                found.append("%s%s" % (p.relative_to(probe).as_posix(),
                                       "/" if p.is_dir() else ""))
            if not found:
                found.append("<EMPTY>")
        else:
            found.append("%s does not exist either" % probe)
        raise SystemExit(
            "REFUSING to run: --ckpt-dataset baked %s and ZERO checkpoint files were found under "
            "%s. A chained session that seeds nothing RECOMPUTES EVERY SEED while looking normal, "
            "so this refuses rather than burn a GPU session on finished work. A file is taken only "
            "when its NAME matches %s -- the `<kind>_<seed>.pt` and `<kind>_<seed>_curve.json` "
            "files written flat by _ckpt_path and _curve_path, for kind 'seed' or 'null' only. "
            "`smoke_*` is deliberately NOT accepted, so a search root holding only smoke output "
            "seeds zero BY DESIGN. Two causes remain and they need opposite fixes: the dataset was "
            "never attached (compare the listing against the names above), or it is genuinely "
            "empty. Present under %s: %s"
            % (BAKED_CKPT_IN, probe, CKPT_NAME_RE.pattern, probe, ", ".join(found)))
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
    # `split_frac` lives under `inherited_not_chosen`, which is where finetune_arm_run.run() reads
    # it (`inh = bar["inherited_not_chosen"]; split_frac = float(inh["comparator"]["split_frac"])`).
    # The smoke MUST use the same value as the real run or it times a different split.
    inh = bar["inherited_not_chosen"]
    # `train_one_seed` calls build() UNCONDITIONALLY and has no default, while `run()` resolves one
    # when build is None (`from scripts.finetune_arm_probe import _build_trainable`). The smoke calls
    # train_one_seed DIRECTLY, so it must do that resolution itself -- and it must resolve the SAME
    # builder, because a smoke that times a different model is not measuring the run it is sizing.
    from scripts.finetune_arm_probe import _build_trainable
    train_one_seed(tiles, seed=0, horizon_epochs=1, lr=hp["learning_rate"],
                   batch_size=hp["batch_size"],
                   split_frac=float(inh["comparator"]["split_frac"]),
                   ckpt_dir=SMOKE_CKPT, build=_build_trainable, kind="smoke", verbose=True)
    spe = time.time() - t0

    after = sorted(p.name for p in CKPT.glob("*")) if CKPT.is_dir() else []
    if before != after:
        raise SystemExit("REFUSED: the smoke run touched the REAL checkpoint dir (%s -> %s). "
                         "That would have become the real trajectory's prefix." % (before, after))
    shutil.rmtree(SMOKE_CKPT, ignore_errors=True)
    print("MEASURED %.1f s/epoch on this device (one real epoch, throwaway ckpt)" % spe, flush=True)
    return spe


def main():
    dev = device_report()
    print("KNOBS smoke_only=%s wall_hours=%.2f ckpt_in=%s (baked: smoke=%r wall=%r ckpt=%r)"
          % (SMOKE_ONLY, WALL_HOURS, CKPT_IN, BAKED_SMOKE_ONLY, BAKED_WALL_HOURS, BAKED_CKPT_IN),
          flush=True)
    WORK.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("PA_REPO", str(WORK / "PromoterAtlas"))
    os.environ.setdefault("MPRA_CACHE", str(WORK / "mpra"))

    bar = unpack()
    require_bar_paths(bar)
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
        print("SMOKE_ONLY -- stopping after the measurement. Nothing was scored and no run "
              "artifact was written.", flush=True)
        projected = spe * n_seeds * horizon / 3600.0
        (WORK / "finetune_arm_smoke.json").write_text(
            json.dumps({{"schema": "finetune-arm-kaggle-smoke-v1",
                         "device": dev,
                         "seconds_per_epoch_measured": spe,
                         "n_seeds": n_seeds,
                         "horizon_epochs": horizon,
                         "projected_hours_at_horizon": projected,
                         "wall_ceiling_hours": WALL_HOURS,
                         "fits_one_session": bool(projected <= WALL_HOURS),
                         "rate_is_for_this_device_only":
                             "seconds_per_epoch is a property of device above, not of the T4 pin. "
                             "If t4_pin_honoured is false a P100 was provisioned and this is a P100 "
                             "rate.",
                         "scored": False}}, indent=2),
            encoding="utf-8")
        return 0

    from scripts.finetune_arm_run import run
    rec = run(bar=bar, tiles=tiles, ckpt_dir=CKPT, build=None,
              horizon_epochs=HORIZON_EPOCHS,
              wall_ceiling_hours=WALL_HOURS, seconds_per_epoch=spe, verbose=True)

    out = WORK / "finetune_arm_run_kaggle.json"
    # Additive and NAMESPACE-SEPARATE. `rec` is what the scorer reads, so the divergence note goes
    # beside it under its own key and never into a field the scorer consumes.
    payload = dict(rec)
    if HORIZON_EPOCHS is not None:
        payload["horizon_divergence"] = {{
            "schema": "finetune-arm-horizon-divergence-v1",
            "horizon_epochs_supplied": int(HORIZON_EPOCHS),
            "planning_bypassed": True,
            "what_was_bypassed":
                "run() calls plan_run only when horizon_epochs is None. This session supplied it, "
                "so no ceiling-derived horizon was computed and wall_ceiling_hours alone decided "
                "when the session stopped.",
            "plan_run_would_have_refused_saying":
                "after reserving N epochs for K null replicates, the C h ceiling leaves H epochs "
                "per seed across S seeds -- fewer than the frozen budget of B. The arm's own "
                "verdict would not fit, so this run must not start.",
            "why_that_refusal_does_not_apply":
                "plan_run asks whether ONE session can hold the whole run. It cannot, which is why "
                "this is a CHAINED sequence: each session self-bounds at wall_ceiling_hours, writes "
                "every completed epoch to the checkpoint directory, and the next session READS "
                "those curves instead of recomputing them. The verdict is assembled across "
                "sessions, so it does not have to fit inside one.",
            "why_a_bounded_stop_is_still_gradeable":
                "run() orders the phases budget -> nulls -> diagnostic extension precisely so a "
                "stop after the first two yields a verdict. A stop INSIDE the budget phase does "
                "not, and the scorer says so by returning INDETERMINATE_BUDGET_EXHAUSTED rather "
                "than a substantive verdict.",
            "wall_ceiling_hours": WALL_HOURS,
            "seconds_per_epoch_measured": spe,
        }}
    out.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
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


def emit(sources_b64: str, mods: dict[str, Path], bar_b64: str, *,
         smoke: bool = False, wall_hours: float | None = None,
         ckpt_dataset: str | None = None, horizon_epochs: int | None = None) -> str:
    """Render the kernel, baking any knob the build set.

    `None` for a knob means "not baked", which leaves the emitted kernel's env read live. That
    distinction is the whole point: Kaggle supplies no environment, so a knob that is not baked is
    pinned to its default there, and a knob that IS baked must not be overridable by an env var that
    can never be set.
    """
    embedded = {m: p.relative_to(REPO).as_posix() for m, p in sorted(mods.items())}
    body = KERNEL.format(pa_url=PA_URL, sources_b64=sources_b64, bar_b64=bar_b64,
                         embedded=embedded,
                         baked_smoke_only=(True if smoke else None),
                         baked_wall_hours=(None if wall_hours is None else float(wall_hours)),
                         baked_ckpt_in=(None if ckpt_dataset is None
                                        else "/kaggle/input/" + ckpt_dataset),
                         baked_horizon_epochs=(None if horizon_epochs is None
                                               else int(horizon_epochs)))
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

    # The flags must reach the EMITTED body, and a non-baked build must NOT look baked.
    smoke_body = emit(sources_b64, mods2, bar_b64, smoke=True, wall_hours=3.5,
                      ckpt_dataset="ft-ckpt-round1", horizon_epochs=38)
    ast.parse(smoke_body)
    assert "BAKED_SMOKE_ONLY = True" in smoke_body, "--smoke did not reach the emitted kernel"
    assert "BAKED_WALL_HOURS = 3.5" in smoke_body, "--wall-hours did not reach the emitted kernel"
    assert "BAKED_CKPT_IN = '/kaggle/input/ft-ckpt-round1'" in smoke_body, \
        "--ckpt-dataset did not reach the emitted kernel"
    assert "BAKED_HORIZON_EPOCHS = 38" in smoke_body, \
        "--horizon-epochs did not reach the emitted kernel"
    assert "BAKED_SMOKE_ONLY = None" in body, "a build without --smoke must not bake a smoke flag"
    assert "BAKED_WALL_HOURS = None" in body, "a build without --wall-hours must not bake one"
    assert "BAKED_HORIZON_EPOCHS = None" in body, (
        "a build without --horizon-epochs must not bake one -- an un-flagged build has to reach "
        "plan_run and still refuse")

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
    ap.add_argument("--smoke", action="store_true",
                    help="bake SMOKE_ONLY: measure ONE real epoch on the provisioned device in a "
                         "throwaway checkpoint dir, write the smoke artifact, score nothing")
    ap.add_argument("--wall-hours", type=float, default=None,
                    help="bake the wall ceiling in hours; unset leaves the kernel's own derived "
                         "10.5 h (12.0 h session cap minus 1.5 h for clone + GEO + artifact write)")
    ap.add_argument("--ckpt-dataset", default=None, metavar="SLUG",
                    help="bake the checkpoint input to /kaggle/input/SLUG so a chained session reads "
                         "finished seeds; pass the same slug to kaggle_push_poll.py push --dataset")
    ap.add_argument("--horizon-epochs", type=int, default=None, metavar="N",
                    help="bake the per-seed horizon, BYPASSING plan_run so wall-hours alone bounds "
                         "the session. Required for chaining: at the measured T4 rate plan_run "
                         "refuses every session-sized ceiling. Unset leaves planning in charge, "
                         "which is correct for a single unchained session.")
    args = ap.parse_args(argv)

    if args.self_test:
        return self_test()
    if not args.out:
        ap.error("--out is required (or pass --self-test)")

    sources_b64, mods, bar_b64 = build_payload()
    body = emit(sources_b64, mods, bar_b64, smoke=args.smoke, wall_hours=args.wall_hours,
                ckpt_dataset=args.ckpt_dataset, horizon_epochs=args.horizon_epochs)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(body, encoding="utf-8")
    print("WROTE %s (%d bytes)" % (out, len(body)))
    print("baked: smoke=%s wall_hours=%s ckpt_dataset=%s horizon_epochs=%s"
          % (args.smoke, args.wall_hours, args.ckpt_dataset, args.horizon_epochs))
    print("embedded %d module(s) + the frozen bar" % len(mods))
    for m in sorted(mods):
        print("  %s" % mods[m].relative_to(REPO).as_posix())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
