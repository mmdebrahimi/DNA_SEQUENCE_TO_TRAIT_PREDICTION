"""MEASURE the storage a Kaggle session actually gets, instead of reading it off a docs page.

WHY THIS EXISTS. The sizing memo's one open question was whether a 301.7-517.1 GiB training tensor
fits free-tier cloud storage. Seven fetch attempts across five Kaggle-owned URLs returned only HTML
titles, because those docs pages are client-rendered -- so the number was recorded as UNVERIFIED
rather than recalled from memory. That framed it as a documentation problem. It is not: the limit is
a property of the running session, so it can be MEASURED, and a measurement of the real surface beats
a published figure anyway (the published number is a policy, the measured one is what the job gets).

CPU-ONLY ON PURPOSE. Disk needs no accelerator, and the weekly GPU allowance is 30 h -- spending any
of it to read `shutil.disk_usage` would be waste. `kaggle quota` already answers the accelerator half
authoritatively from the API (GPU 30.00 h / TPU 20.00 h per week, measured 2026-11-25).

WHAT THIS CANNOT MEASURE, stated so the artifact does not over-claim: the per-dataset and per-account
UPLOAD ceilings are server-side policy and are not visible from inside a kernel. What IS measurable
is every mount a session can read or write, which is what the decision actually turns on -- a tensor
that cannot be staged anywhere a session can reach is blocked regardless of what the upload policy
permits.
"""
import json
import os
import shutil
import subprocess

GIB = 2 ** 30
TENSOR_GIB = {           # from wiki/encode_roadmap_training_budget_sources.json, sourced inputs
    "896 bins float16": 34021 * 896 * 5313 * 2 / GIB,
    "896 bins float32": 34021 * 896 * 5313 * 4 / GIB,
    "1536 bins float16": 34021 * 1536 * 5313 * 2 / GIB,
}

PATHS = ["/", "/kaggle", "/kaggle/working", "/kaggle/input", "/kaggle/temp", "/tmp", "/root",
         "/opt/conda", "/dev/shm"]

out = {"record": "kaggle-session-storage-measured-v1", "paths": {}, "notes": []}

for p in PATHS:
    row = {"exists": os.path.isdir(p)}
    if row["exists"]:
        try:
            u = shutil.disk_usage(p)
            row.update(total_gib=round(u.total / GIB, 2),
                       used_gib=round(u.used / GIB, 2),
                       free_gib=round(u.free / GIB, 2))
        except OSError as e:
            row["error"] = str(e)
        row["writable"] = os.access(p, os.W_OK)
    out["paths"][p] = row

# NESTED DATASET MOUNTS (added 2026-11-25 after this probe's first run misled me).
# `shutil.disk_usage("/kaggle/input")` reports the device backing that DIRECTORY -- the same ~20 GiB
# /dev/loop2 as /kaggle/working -- while each attached dataset is mounted SEPARATELY beneath it via
# gcsfuse and reports ~1.0 PiB. Reading the per-path figure alone would have manufactured a false
# ~19.5 GiB cap on attached data. Enumerate the nested mounts so the number cannot be misread.
out["dataset_mounts"] = {}
_base = "/kaggle/input"
if os.path.isdir(_base):
    for root, dirs, _ in os.walk(_base):
        if os.path.ismount(root) and root != _base:
            try:
                u = shutil.disk_usage(root)
                out["dataset_mounts"][root] = {"total_gib": round(u.total / GIB, 2),
                                               "free_gib": round(u.free / GIB, 2),
                                               "writable": os.access(root, os.W_OK)}
            except OSError as e:
                out["dataset_mounts"][root] = {"error": str(e)}
            dirs[:] = []            # a mount's interior is the dataset, not more mounts
if not out["dataset_mounts"]:
    out["notes"].append("no nested dataset mount found -- attach a dataset or the per-path "
                        "/kaggle/input figure will be the loop device, NOT the dataset capacity")

# Distinct filesystems, so two mounts on one device are not counted twice.
try:
    out["df"] = subprocess.run(["df", "-h"], capture_output=True, text=True, timeout=60).stdout
except Exception as e:                                           # noqa: BLE001
    out["df"] = f"unavailable: {e}"
try:
    out["mounts_kaggle"] = [ln for ln in open("/proc/mounts", encoding="utf-8")
                            if "kaggle" in ln or " /tmp " in ln or ln.split()[1] == "/"]
except Exception as e:                                           # noqa: BLE001
    out["mounts_kaggle"] = [f"unavailable: {e}"]

out["env_kaggle"] = {k: v for k, v in sorted(os.environ.items())
                     if "KAGGLE" in k.upper() and "KEY" not in k.upper()
                     and "TOKEN" not in k.upper() and "SECRET" not in k.upper()}

# The actual question. Compare the tensor against the LARGEST writable free space found.
writable = {p: r["free_gib"] for p, r in out["paths"].items()
            if r.get("writable") and isinstance(r.get("free_gib"), (int, float))}
best_p, best_free = (max(writable.items(), key=lambda kv: kv[1]) if writable else (None, 0.0))
out["verdict"] = {
    "largest_writable_mount": best_p,
    "largest_writable_free_gib": best_free,
    "training_tensor_gib": {k: round(v, 1) for k, v in TENSOR_GIB.items()},
    "fits": {k: bool(best_free >= v) for k, v in TENSOR_GIB.items()},
    "shortfall_gib_smallest_variant": round(min(TENSOR_GIB.values()) - best_free, 1),
    "what_this_does_NOT_measure": (
        "the per-dataset and per-account UPLOAD ceilings, which are server-side policy and invisible "
        "from inside a kernel. This measures what a running session can read and write."),
}

print(json.dumps(out, indent=2))
with open("/kaggle/working/kaggle_limits_measured.json", "w", encoding="utf-8") as fh:
    json.dump(out, fh, indent=2)
print("\nSMALLEST TENSOR VARIANT: %.1f GiB | LARGEST WRITABLE FREE: %.1f GiB (%s) | FITS: %s"
      % (min(TENSOR_GIB.values()), best_free, best_p,
         best_free >= min(TENSOR_GIB.values())))
