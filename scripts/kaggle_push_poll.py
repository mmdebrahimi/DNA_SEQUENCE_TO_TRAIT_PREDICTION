#!/usr/bin/env python
"""Push a self-contained Kaggle SCRIPT kernel + poll its status / pull its output.

Reuses the authenticated kaggle module (auth already works on this host as emanueleebrahimi). A script
kernel runs one .py file headless on Kaggle's CPU/GPU with internet, up to 12 h (CPU). Free tier -> no
money gate.

  push  : uv run python scripts/kaggle_push_poll.py push  <code.py> <slug> [--gpu] [--no-internet]
  status: uv run python scripts/kaggle_push_poll.py status <slug>
  pull  : uv run python scripts/kaggle_push_poll.py pull   <slug> <dest_dir>
"""
import argparse
import json
import shutil
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))


def _api():
    from kaggle.api.kaggle_api_extended import KaggleApi
    api = KaggleApi()
    api.authenticate()
    return api


def push(code_py: str, slug: str, gpu: bool, internet: bool, dataset: str | None = None) -> int:
    api = _api()
    user = api.get_config_value("username") or "emanueleebrahimi"
    kid = f"{user}/{slug}"
    with tempfile.TemporaryDirectory() as td:
        tdp = Path(td)
        code_name = Path(code_py).name
        shutil.copyfile(code_py, tdp / code_name)
        meta = {
            "id": kid, "title": slug, "code_file": code_name,
            "language": "python", "kernel_type": "script",
            "is_private": True, "enable_gpu": bool(gpu), "enable_internet": bool(internet),
        }
        if dataset:
            # attach a dataset (e.g. the precomputed GEMME tables) at /kaggle/input/<slug>/
            meta["dataset_sources"] = [dataset if "/" in dataset else f"{user}/{dataset}"]
        if gpu:
            # enable_gpu alone provisions a Tesla P100 (CC 6.0) that Kaggle's current torch (CC 7.0+)
            # cannot run ("no kernel image available"). Pin the T4. See memory
            # reference_kaggle_headless_gpu_kernels.
            meta["machine_shape"] = "NvidiaTeslaT4"
        (tdp / "kernel-metadata.json").write_text(json.dumps(meta, indent=2))
        print(f"pushing {kid} (gpu={gpu} internet={internet}, code={code_name})...")
        res = api.kernels_push(str(tdp))
        print("push result:", getattr(res, "ref", res), getattr(res, "url", ""))
        print(f"watch: https://www.kaggle.com/code/{kid}")
    return 0


def status(slug: str) -> int:
    api = _api()
    user = api.get_config_value("username") or "emanueleebrahimi"
    st = api.kernels_status(f"{user}/{slug}")
    print("status:", getattr(st, "status", st), "| failureMessage:", getattr(st, "failureMessage", ""))
    return 0


def pull(slug: str, dest: str) -> int:
    api = _api()
    user = api.get_config_value("username") or "emanueleebrahimi"
    Path(dest).mkdir(parents=True, exist_ok=True)
    api.kernels_output(f"{user}/{slug}", path=dest)
    print(f"pulled output of {user}/{slug} -> {dest}")
    for p in sorted(Path(dest).iterdir()):
        print("  ", p.name)
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("push"); p.add_argument("code_py"); p.add_argument("slug")
    p.add_argument("--gpu", action="store_true"); p.add_argument("--no-internet", action="store_true")
    p.add_argument("--dataset", default=None, help="attach a Kaggle dataset slug at /kaggle/input/<slug>/")
    s = sub.add_parser("status"); s.add_argument("slug")
    q = sub.add_parser("pull"); q.add_argument("slug"); q.add_argument("dest")
    a = ap.parse_args()
    if a.cmd == "push":
        return push(a.code_py, a.slug, a.gpu, not a.no_internet, a.dataset)
    if a.cmd == "status":
        return status(a.slug)
    if a.cmd == "pull":
        return pull(a.slug, a.dest)
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
