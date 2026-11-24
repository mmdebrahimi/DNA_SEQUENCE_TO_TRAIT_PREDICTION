"""Kaggle CPU kernel: run BOTH re-pretrain arms to 6000 steps and emit the two checkpoints.

WHY KAGGLE AND WHY CPU. The laptop this normally runs on produced four kernel bugchecks in about
100 minutes under exactly this load (Event 41/6008, codes 0x3B / 0xEF / 0xA / 0x139), killing the
arm three times -- once at step 3500 of 6000 after 3 h 13 m. And `scripts/glm_train_diffusion.py`
contains ZERO references to cuda / device / .to(), so this is a pure CPU torch job: a GPU would be
idle and would spend the 30 h/week accelerator quota for nothing. CPU kernels also get a 12 h cap
rather than 9 h, which is the margin this needs.

THE ONE SCIENTIFIC CONSTRAINT, and it is why both arms run in ONE kernel rather than two.
Run 2 is the NOISE FLOOR for run 1. Its entire job is to let any environment difference against the
already-trained committed checkpoint cancel out of the comparison, which only works if both new runs
share that environment exactly. Two kernels would be two containers and, worse, could land on two
different base images. One kernel makes "same machine, same torch build" true by construction, and
this script records torch version / platform / core count so the claim is checkable afterwards.

The code is pinned to a COMMIT, not a branch, because a pin selects the code that runs and a branch
does not. Output: the two `latest.pt` files, both logs, and a manifest carrying each checkpoint's
INTERNAL step plus its sha256 -- the grader refuses anything whose step is not 6000, since a
partially-trained checkpoint scores perfectly cleanly and would otherwise read as the finished arm.
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

#: THE CORPORA DO NOT COME FROM THE CLONE. v1 of this kernel assumed they were tracked and refused
#: at 11.9 s with `missing corpus wiki/glm_training_corpus_representative_2026-11-22.txt` -- they are
#: DELIBERATELY gitignored BY NAME at .gitignore:187-188 (9.2 MB and 9.7 MB), so no clone can ever
#: carry them. They ride in as a private Kaggle dataset instead, and the clone supplies CODE only.
#: The guard firing in seconds rather than after six hours is why it sits before the cost.
DATASET = Path("/kaggle/input/glm-repretrain-corpora")

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


def main() -> int:
    t0 = time.time()
    cores = os.cpu_count() or 2
    import torch
    env_rec = {
        "torch": torch.__version__,
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "machine": platform.machine(),
        "cpu_count": cores,
        "pinned_commit": SHA,
    }
    print("ENV", json.dumps(env_rec), flush=True)

    # ---- code, pinned to a commit ----
    if SRC.exists():
        shutil.rmtree(SRC)
    sh("git", "clone", "--quiet", REPO, SRC)
    sh("git", "checkout", "--quiet", SHA, cwd=SRC)
    head = subprocess.run(["git", "rev-parse", "HEAD"], cwd=SRC, capture_output=True,
                           text=True).stdout.strip()
    if head != SHA:
        print(f"REFUSING: checkout landed on {head}, not the pinned {SHA}", file=sys.stderr)
        return 2

    # ---- stage the corpora from the dataset, then verify BEFORE spending hours ----
    if not DATASET.exists():
        print(f"REFUSING: the corpora dataset is not mounted at {DATASET}. The corpora are "
              f"gitignored, so the clone cannot supply them; attach the dataset to the kernel.",
              file=sys.stderr)
        return 3
    print("dataset contents:", sorted(p.name for p in DATASET.iterdir()), flush=True)

    for label, corpus, man, _seed, _cd, _tag, want, want_sha in ARMS:
        for rel in (corpus, man):
            src, dst = DATASET / Path(rel).name, SRC / rel
            if not src.exists():
                print(f"REFUSING: {src.name} absent from the dataset", file=sys.stderr)
                return 3
            dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(src, dst)
        p = SRC / corpus
        # sha256 FIRST: it pins the exact bytes, where a line count only pins the shape.
        got_sha = hashlib.sha256(p.read_bytes()).hexdigest()
        n = sum(1 for _ in p.open(encoding="utf-8"))
        print(f"{label}: {Path(corpus).name} -> {n} lines, sha256 {got_sha[:16]}... "
              f"(expect {want} lines, {want_sha[:16]}...)", flush=True)
        if got_sha != want_sha:
            print(f"REFUSING: {corpus} sha256 is {got_sha}, expected {want_sha}. The dataset does "
                  f"not hold the corpus this arm was specified against.", file=sys.stderr)
            return 3
        if n != want:
            print(f"REFUSING: {corpus} has {n} lines, expected {want}. Wrong corpus.",
                  file=sys.stderr)
            return 3

    # ---- launch both arms in parallel, partitioning the cores ----
    threads = max(1, cores // 2)
    print(f"cores={cores} -> OMP_NUM_THREADS={threads} per arm, both in parallel", flush=True)
    procs, logs = {}, {}
    for label, corpus, man, seed, ckpt_dir, tag, _want, _sha in ARMS:
        (SRC / ckpt_dir).mkdir(parents=True, exist_ok=True)
        log = OUT / f"{label}.log"
        logs[label] = log
        env = dict(os.environ,
                   OMP_NUM_THREADS=str(threads), MKL_NUM_THREADS=str(threads),
                   PYTHONUNBUFFERED="1")
        fh = log.open("w", encoding="utf-8")
        procs[label] = (subprocess.Popen(
            [sys.executable, "-u", "scripts/glm_train_diffusion.py",
             "--corpus", corpus, "--manifest", man,
             "--steps", str(STEPS), "--batch", "64", "--lr", "3e-4",
             "--seed", str(seed), "--ckpt-dir", ckpt_dir, "--tag", tag],
            cwd=SRC, env=env, stdout=fh, stderr=subprocess.STDOUT), fh, ckpt_dir, tag)
        print(f"launched {label} pid={procs[label][0].pid}", flush=True)

    # ---- wait, reporting progress from the logs rather than guessing ----
    while True:
        alive = {k: v for k, v in procs.items() if v[0].poll() is None}
        el = (time.time() - t0) / 3600.0
        tails = {}
        for label, log in logs.items():
            try:
                ls = [l for l in log.read_text(encoding="utf-8", errors="replace").splitlines()
                      if " step " in l]
                tails[label] = ls[-1].strip() if ls else "(no step line yet)"
            except OSError:
                tails[label] = "(log unreadable)"
        print(f"[{el:5.2f} h] alive={list(alive)} | " +
              " | ".join(f"{k}: {v}" for k, v in tails.items()), flush=True)
        if not alive:
            break
        if el > CAP_HOURS:
            print(f"CAP: {el:.2f} h exceeds {CAP_HOURS} h; terminating so output can be written. "
                  f"Any checkpoint below step {STEPS} is UNUSABLE and the grader will refuse it.",
                  file=sys.stderr, flush=True)
            for v in alive.values():
                v[0].terminate()
            break
        time.sleep(300)

    # ---- collect, and report each checkpoint's INTERNAL step ----
    manifest = {"env": env_rec, "wall_hours": round((time.time() - t0) / 3600.0, 3),
                "steps_requested": STEPS, "arms": {}}
    for label, (proc, fh, ckpt_dir, tag) in procs.items():
        fh.close()
        rc = proc.poll()
        rec = {"returncode": rc, "ckpt_dir": ckpt_dir, "tag": tag}
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
            rec["note"] = "no checkpoint written -- the arm did not reach the first save"
        # the artifact the trainer writes for itself, if it got that far. Keyed on the --tag it was
        # given, NOT on the checkpoint directory -- `glm_ckpt_reseed_control` would yield "control"
        # while the tag is "reseed_seed1", so the glob would silently match nothing.
        for extra in sorted((SRC / "wiki").glob(f"*{tag}*.json")):
            shutil.copy2(extra, OUT / extra.name)
            rec.setdefault("extra_artifacts", []).append(extra.name)
        manifest["arms"][label] = rec
        print(f"{label}: rc={rc} step={rec['internal_step']} "
              f"finished={rec['finished_at_requested_steps']}", flush=True)

    ok = all(v["finished_at_requested_steps"] for v in manifest["arms"].values())
    manifest["BOTH_ARMS_FINISHED"] = ok
    manifest["what_to_do_next"] = (
        "pull both *_latest.pt, score each on orthology_free (and, for the leak-identification "
        "layer, on ortholog_bearing and withheld_orthology_free) with scripts/glm_beyond_kmer.py "
        "passing an EXPLICIT --out, then grade with scripts/glm_repretrain_grade.py"
        if ok else
        "DO NOT SCORE. At least one arm did not reach step 6000; the grader refuses those by "
        "design because a partial checkpoint scores cleanly and would read as the finished arm.")
    (OUT / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print("MANIFEST", json.dumps(manifest["arms"], default=str), flush=True)
    print("BOTH_ARMS_FINISHED =", ok, flush=True)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
