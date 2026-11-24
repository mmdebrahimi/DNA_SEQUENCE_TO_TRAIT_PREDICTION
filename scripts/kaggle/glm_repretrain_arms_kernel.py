"""Kaggle CPU kernel: run BOTH re-pretrain arms to 6000 steps and emit the two checkpoints.

WHY KAGGLE AND WHY CPU. The laptop this normally runs on produced four kernel bugchecks in about
100 minutes under exactly this load (Event 41/6008, codes 0x3B / 0xEF / 0xA / 0x139), killing the
arm three times -- once at step 3500 of 6000 after 3 h 13 m. And `scripts/glm_train_diffusion.py`
never touches cuda, so this is a pure CPU torch job: a GPU would idle while spending the 30 h/week
accelerator quota. CPU kernels also cap near 12 h rather than 9 h, which is the margin this needs.

WHY THIS VERSION LOOKS DIFFERENT FROM v1-v4. Those four each tested ONE MORE PREDICTION about how
Kaggle differs from the local venv, and each cost a push-and-wait cycle to disprove:
  v1  the corpora are tracked        -> they are gitignored BY NAME (.gitignore:187-188)
  v2  the dataset was attached late  -> an unready dataset mounts as NOTHING, silently
  v3  it mounts at /kaggle/input/<slug> -> it mounted NESTED under /kaggle/input/datasets/
  v4  the trainer needs only torch   -> `from dna_decode.glm.diffusion import ...` sits at LINE 129
                                        INSIDE main(), and importing it runs dna_decode/glm/__init__
                                        which pulls in Bio (biopython) and pandas
So this version stops predicting. It verifies the dependency chain by IMPORTING it, then runs the
REAL configuration for a few dozen steps before committing to 6000 -- which exercises imports,
corpus load, the training step, the final save, sample generation and artifact writing, all in
minutes. A guard belongs before the cost, not after it.

IT ALSO REFUSES ON TIME. No Kaggle step rate has ever been measured for this job, and a run cut at
the cap short of step 6000 yields NOTHING usable -- the grader refuses any checkpoint whose internal
step is not 6000, because a partial one scores perfectly cleanly and would read as the finished arm.
So the pre-flight measures the rate from a two-point slope (which cancels fixed overhead) and
projects; if 6000 steps will not fit with margin, it stops and REPORTS THE MEASURED RATE rather than
silently picking a workaround. Splitting the arms or chaining --resume both change what the noise
floor controls for, so that is an authority call, not a kernel decision.

THE ONE SCIENTIFIC CONSTRAINT, and why both arms share one kernel. Run 2 is the NOISE FLOOR for
run 1. Its whole job is to let any environment difference against the already-trained committed
checkpoint cancel out of the comparison, which only holds if both new runs share that environment.
Two kernels would be two containers and could land on different base images. One kernel makes it
true by construction, and the ENV record below makes the claim checkable afterwards.
"""
import hashlib
import json
import os
import platform
import shutil
import subprocess
import sys
import time
from pathlib import Path

REPO = "https://github.com/mmdebrahimi/DNA_SEQUENCE_TO_TRAIT_PREDICTION.git"
SHA = "b10b54750d2ad93cc3261971fde89f8bf8268637"
SRC = Path("/tmp/dna_decode")
OUT = Path("/kaggle/working")
STEPS = 6000
CAP_HOURS = 11.5          # the kernel is killed near 12 h; leave room to write output
FIT_MARGIN = 1.25         # projected time must clear the cap by this factor

#: Third-party packages the import chain actually needs, walked with ast over
#: dna_decode/__init__.py, dna_decode/glm/__init__.py and dna_decode/glm/diffusion.py.
#: `Bio` is the one v4 would never have reached: it arrives via dna_decode.glm.corpus.
NEEDED = {"torch": "torch", "pandas": "pandas", "Bio": "biopython"}

#: There is deliberately NO mount-path constant. `/kaggle/input/<slug>` is a convention, not a
#: guarantee (v3 proved it), so the corpora are LOCATED BY NAME and verified by sha256.
#: (label, corpus, manifest, seed, ckpt_dir, tag, expected lines, expected sha256)
ARMS = [
    ("run1_representative", "wiki/glm_training_corpus_representative_2026-11-22.txt",
     "wiki/glm_training_corpus_representative_manifest_2026-11-22.json", 0,
     "data/glm_ckpt_representative", "representative", 60462,
     "228d705e77c9206a1b4d58adf7bf0a7563b8775eb2996aabec91953c9e01ba58"),
    ("run2_reseed_seed1", "wiki/glm_training_corpus_2026-11-21.txt",
     "wiki/glm_training_corpus_manifest_2026-11-21.json", 1,
     "data/glm_ckpt_reseed_control", "reseed_seed1", 63925,
     "f1479ef524aa463c5c7512841ba2c306fdb774d7dba073baf52e08a3a5a33333"),
]


def sh(*cmd, cwd=None):
    print("+", " ".join(str(c) for c in cmd), flush=True)
    subprocess.run([str(c) for c in cmd], cwd=cwd, check=True)


def arm_env(threads: int) -> dict:
    """PYTHONPATH is load-bearing: v4 died because the clone root was not importable.

    Locally `dna_decode` resolves only because `uv sync` installs the project into the venv. The
    Kaggle image has no such install, and running `python scripts/x.py` puts `scripts/` on sys.path
    rather than the repo root, so the package is invisible without this.
    """
    return dict(os.environ, PYTHONPATH=str(SRC), OMP_NUM_THREADS=str(threads),
                MKL_NUM_THREADS=str(threads), PYTHONUNBUFFERED="1")


def run_arms(steps: int, threads: int, *, suffix: str, n_samples: int | None,
             logdir: Path) -> dict:
    """Launch both arms in PARALLEL and wait. Returns {label: {rc, wall_s, ckpt}}.

    Parallel in the pre-flight too, deliberately: the measured rate must reflect the same core
    contention the real run will face, or the projection is of a different experiment.
    """
    procs, out = {}, {}
    for label, corpus, man, seed, ckpt_dir, tag, _w, _s in ARMS:
        cd = f"{ckpt_dir}{suffix}"
        (SRC / cd).mkdir(parents=True, exist_ok=True)
        cmd = [sys.executable, "-u", "scripts/glm_train_diffusion.py",
               "--corpus", corpus, "--manifest", man,
               "--steps", str(steps), "--batch", "64", "--lr", "3e-4",
               "--seed", str(seed), "--ckpt-dir", cd, "--tag", f"{tag}{suffix}"]
        if n_samples is not None:
            cmd += ["--n-samples", str(n_samples)]
        log = logdir / f"{label}{suffix}.log"
        fh = log.open("w", encoding="utf-8")
        procs[label] = (subprocess.Popen(cmd, cwd=SRC, env=arm_env(threads),
                                         stdout=fh, stderr=subprocess.STDOUT),
                        fh, cd, log, time.time())
    for label, (proc, fh, cd, log, t_start) in procs.items():
        proc.wait()
        fh.close()
        out[label] = {"rc": proc.returncode, "wall_s": round(time.time() - t_start, 1),
                      "ckpt": (SRC / cd / "latest.pt").exists(), "log": str(log),
                      "tail": "\n".join(log.read_text(encoding="utf-8", errors="replace")
                                        .splitlines()[-12:])}
    return out


def main() -> int:
    t0 = time.time()
    cores = os.cpu_count() or 2
    threads = max(1, cores // 2)
    env_rec = {"python": sys.version.split()[0], "platform": platform.platform(),
               "machine": platform.machine(), "cpu_count": cores,
               "threads_per_arm": threads, "pinned_commit": SHA}
    OUT.mkdir(parents=True, exist_ok=True)

    # ---- code, pinned to a COMMIT (a branch does not pin what runs) ----
    if SRC.exists():
        shutil.rmtree(SRC)
    sh("git", "clone", "--quiet", REPO, SRC)
    sh("git", "checkout", "--quiet", SHA, cwd=SRC)
    head = subprocess.run(["git", "rev-parse", "HEAD"], cwd=SRC, capture_output=True,
                          text=True).stdout.strip()
    if head != SHA:
        print(f"REFUSING: checkout landed on {head}, not the pinned {SHA}", file=sys.stderr)
        return 2

    # ---- corpora: LOCATE BY NAME under /kaggle/input, then verify the BYTES ----
    root = Path("/kaggle/input")
    if not root.exists():
        print("REFUSING: /kaggle/input does not exist, so no dataset is attached.", file=sys.stderr)
        return 3
    found: dict[str, list[Path]] = {}
    for p in root.rglob("*"):
        if p.is_file():
            found.setdefault(p.name, []).append(p)
    print(f"/kaggle/input top level: {sorted(q.name for q in root.iterdir())}", flush=True)
    print(f"/kaggle/input holds {sum(len(v) for v in found.values())} files", flush=True)

    for label, corpus, man, _sd, _cd, _tg, want, want_sha in ARMS:
        for rel in (corpus, man):
            hits = found.get(Path(rel).name, [])
            if len(hits) != 1:
                # Ambiguity is REFUSED, not resolved by sort order: choosing between two files of
                # the same name is how the wrong corpus gets trained on silently, and the arm and
                # its noise floor differ ONLY by which corpus they saw.
                print(f"REFUSING: {Path(rel).name} found {len(hits)} times under /kaggle/input"
                      + (f": {[str(h) for h in hits]}" if hits else
                         ". The corpora are gitignored (.gitignore:187-188) so the clone cannot "
                         "supply them; the dataset must be attached AND `kaggle datasets status` "
                         "must read `ready` at push time."), file=sys.stderr)
                return 3
            dst = SRC / rel
            dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(hits[0], dst)
        p = SRC / corpus
        got_sha = hashlib.sha256(p.read_bytes()).hexdigest()
        n = sum(1 for _ in p.open(encoding="utf-8"))
        print(f"{label}: {Path(corpus).name} -> {n} lines, sha256 {got_sha[:16]}...", flush=True)
        if got_sha != want_sha:
            print(f"REFUSING: {corpus} sha256 {got_sha} != expected {want_sha}", file=sys.stderr)
            return 3
        if n != want:
            print(f"REFUSING: {corpus} has {n} lines, expected {want}", file=sys.stderr)
            return 3

    # ---- dependencies: IMPORT them, do not assume them. v4's whole failure mode. ----
    probe = ("import json,importlib;"
             "print(json.dumps({m: getattr(importlib.import_module(m),'__version__','?') "
             "for m in " + repr(sorted(NEEDED)) + "}))")
    r = subprocess.run([sys.executable, "-c", probe], cwd=SRC, env=arm_env(threads),
                       capture_output=True, text=True)
    if r.returncode != 0:
        missing = [pip for mod, pip in NEEDED.items() if f"'{mod}'" in r.stderr or mod in r.stderr]
        print(f"import probe failed; installing {missing or list(NEEDED.values())}", flush=True)
        subprocess.run([sys.executable, "-m", "pip", "install", "--quiet",
                        *(missing or NEEDED.values())], check=False)
        r = subprocess.run([sys.executable, "-c", probe], cwd=SRC, env=arm_env(threads),
                           capture_output=True, text=True)
        if r.returncode != 0:
            print(f"REFUSING: dependency chain still unimportable.\n{r.stderr[-1500:]}",
                  file=sys.stderr)
            return 4
    env_rec["deps"] = json.loads(r.stdout.strip().splitlines()[-1])
    print("ENV", json.dumps(env_rec), flush=True)

    # ---- PRE-FLIGHT: the REAL configuration, briefly, at two step counts ----
    # Two points so the slope cancels fixed overhead (torch import, corpus load, sample gen).
    # Throwaway ckpt dirs and distinct tags so nothing here can contaminate the real arms.
    pre = {}
    for steps, suffix in ((20, "_pre20"), (120, "_pre120")):
        print(f"--- pre-flight: {steps} steps, both arms in parallel ---", flush=True)
        pre[steps] = run_arms(steps, threads, suffix=suffix, n_samples=8, logdir=OUT)
        for label, rec in pre[steps].items():
            print(f"  {label}: rc={rec['rc']} wall={rec['wall_s']}s ckpt={rec['ckpt']}", flush=True)
            if rec["rc"] != 0 or not rec["ckpt"]:
                print(f"REFUSING: pre-flight failed for {label} at {steps} steps "
                      f"(rc={rec['rc']}, checkpoint written={rec['ckpt']}). The real run would "
                      f"have failed the same way after hours.\n--- tail ---\n{rec['tail']}",
                      file=sys.stderr)
                (OUT / "preflight-failure.json").write_text(
                    json.dumps({"env": env_rec, "preflight": pre}, indent=2, default=str),
                    encoding="utf-8")
                return 5

    # ---- PROJECT, and refuse on time rather than discover it at the cap ----
    proj = {}
    for label in pre[20]:
        w20, w120 = pre[20][label]["wall_s"], pre[120][label]["wall_s"]
        slope = (w120 - w20) / 100.0
        fixed = max(0.0, w20 - 20 * slope)
        total_h = (fixed + STEPS * slope) / 3600.0
        proj[label] = {"wall_20s": w20, "wall_120s": w120, "s_per_step": round(slope, 4),
                       "fixed_overhead_s": round(fixed, 1),
                       "projected_hours": round(total_h, 2)}
        print(f"  {label}: {slope:.3f} s/step, overhead {fixed:.0f}s -> "
              f"{total_h:.2f} h for {STEPS} steps", flush=True)
    worst = max(p["projected_hours"] for p in proj.values())
    need = worst * FIT_MARGIN
    print(f"worst projection {worst:.2f} h; with {FIT_MARGIN}x margin {need:.2f} h "
          f"vs cap {CAP_HOURS} h", flush=True)
    if need > CAP_HOURS:
        print(f"REFUSING: {STEPS} steps will not fit. Measured rate is in the artifact. "
              f"Splitting the arms across kernels or chaining --resume BOTH change what the noise "
              f"floor controls for, so that is an authority call, not a kernel decision.",
              file=sys.stderr)
        (OUT / "manifest.json").write_text(json.dumps(
            {"env": env_rec, "preflight": pre, "projection": proj,
             "VERDICT": "REFUSED_WILL_NOT_FIT_IN_CAP",
             "cap_hours": CAP_HOURS, "margin": FIT_MARGIN}, indent=2, default=str),
            encoding="utf-8")
        return 6

    # ---- the real arms ----
    print(f"--- real run: {STEPS} steps per arm, both in parallel ---", flush=True)
    real = run_arms(STEPS, threads, suffix="", n_samples=None, logdir=OUT)

    manifest = {"env": env_rec, "preflight": pre, "projection": proj,
                "wall_hours": round((time.time() - t0) / 3600.0, 3),
                "steps_requested": STEPS, "arms": {}}
    import torch
    for (label, _c, _m, _sd, ckpt_dir, tag, _w, _s) in ARMS:
        rec = dict(real[label]); rec.pop("tail", None)
        rec["ckpt_dir"], rec["tag"] = ckpt_dir, tag
        src = SRC / ckpt_dir / "latest.pt"
        if src.exists():
            dst = OUT / f"{label}_latest.pt"
            shutil.copy2(src, dst)
            blob = torch.load(dst, map_location="cpu", weights_only=False)
            rec["internal_step"] = int(blob.get("step", -1))
            rec["finished_at_requested_steps"] = rec["internal_step"] == STEPS
            rec["config"] = blob.get("config")
            rec["bytes"] = dst.stat().st_size
            rec["sha256"] = hashlib.sha256(dst.read_bytes()).hexdigest()
            rec["output_file"] = dst.name
        else:
            rec["internal_step"] = None
            rec["finished_at_requested_steps"] = False
        for extra in sorted((SRC / "wiki").glob(f"*{tag}*.json")):
            shutil.copy2(extra, OUT / extra.name)
        manifest["arms"][label] = rec
        print(f"{label}: rc={rec['rc']} step={rec['internal_step']} "
              f"finished={rec['finished_at_requested_steps']}", flush=True)

    ok = all(v["finished_at_requested_steps"] for v in manifest["arms"].values())
    manifest["BOTH_ARMS_FINISHED"] = ok
    manifest["what_to_do_next"] = (
        "pull both *_latest.pt, score each on orthology_free (plus ortholog_bearing and "
        "withheld_orthology_free for the leak-identification layer) with scripts/glm_beyond_kmer.py "
        "passing an EXPLICIT --out, then grade with scripts/glm_repretrain_grade.py"
        if ok else
        "DO NOT SCORE. An arm did not reach step 6000, and the grader refuses those by design "
        "because a partial checkpoint scores cleanly and would read as the finished arm.")
    (OUT / "manifest.json").write_text(json.dumps(manifest, indent=2, default=str),
                                       encoding="utf-8")
    print("BOTH_ARMS_FINISHED =", ok, flush=True)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
