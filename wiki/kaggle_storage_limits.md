# Kaggle session storage and accelerator limits — MEASURED

**Date:** 2026-11-25
**Expires:** constraint `regulatory-training-tensor-vs-unverified-cloud-storage`
**Method:** `scripts/kaggle/kaggle_limits_probe.py` run as a CPU-only kernel
(`emanueleebrahimi/kaggle-limits-probe-v1`) + `kaggle quota` from the API.
**Raw:** `wiki/kaggle_limits_measured_2026-11-25.json`

---

## DECISION: storage does not block the regulatory goal

The sizing memo left one question open — whether a **301.7–603.3 GiB** training tensor fits free-tier
cloud storage. It does, with large margin, and the answer came from **measuring a running session**
rather than from a docs page.

| surface | measured | backing device |
|---|---|---|
| ephemeral scratch (`/`, `/tmp`, `/root`, `/kaggle`) | **8.0 TiB total, 2.5 TiB free** | `overlay` |
| **attached datasets** (`/kaggle/input/datasets/<owner>/<name>`) | **1.0 PiB reported, read-only** | `gcsfuse` (`kkb-production-mirror-us-west1`) |
| persistent output (`/kaggle/working`) | **19.52 GiB total, 19.5 free** | `/dev/loop2` |
| shared memory (`/dev/shm`) | 14.0 GiB | `shm` |
| **GPU allowance** | **30.00 h/week** (29.94 remaining, refreshes 2026-11-28) | `kaggle quota` |
| **TPU allowance** | **20.00 h/week** | `kaggle quota` |

**Reading it.** The tensor fits ephemeral scratch roughly **four to eight times over**. More
importantly, an attached dataset is a **gcsfuse mount streamed from cloud storage**, so it consumes
**no local disk at all** — a large dataset does not have to fit the session's filesystem to be read
by it.

**So the binding limits are not storage.** They are (a) **persistent output at 19.5 GiB**, which
constrains what a run may *keep*, not what it may read, and (b) **30 GPU hours per week**, which at
the measured 31,509 steps per 9-hour session is about **3.3 sessions per week**. For a 134.63M-param
model a checkpoint is roughly 0.5–1.6 GiB, so 19.5 GiB holds a dozen or more.

---

## CORRECTION — my own probe's single number was misleading

**`/kaggle/input` reported 19.52 GiB, and taking that at face value would have manufactured a false
cap.** `shutil.disk_usage` on a mount *point* reports the device backing that **directory**, and
`/kaggle/input` is on the same 20 GiB `/dev/loop2` as `/kaggle/working`. The datasets themselves are
mounted **nested beneath it**, one gcsfuse mount per dataset, and those report 1.0 PiB.

So a probe that printed only `disk_usage("/kaggle/input")` would have concluded that attached data is
capped at ~19.5 GiB — which would have **re-registered a storage wall that does not exist**, days
after retracting a different false wall in the same memo. What disambiguated it was `df -h` and
`/proc/mounts`, captured in the same run precisely because a single figure is not a measurement of a
layered filesystem.

The probe has since been corrected to enumerate nested mounts under `/kaggle/input` directly, so a
future run does not mislead. **This run's verdict rests on the `df` and `/proc/mounts` evidence, not
on the per-path figure.**

---

## Honest limits

1. **The dataset UPLOAD ceiling is still unmeasured.** Per-dataset and per-account upload limits are
   server-side policy and are not visible from inside a kernel. What is established is that a
   *mounted* dataset costs no local disk — not that an arbitrarily large one can be uploaded. The
   largest dataset in this account today is 186 MB, which is evidence about what has been uploaded,
   not about the ceiling.
2. **Scratch is EPHEMERAL.** The 2.5 TiB disappears when the session ends. Data must arrive each run
   — as an attached dataset (free, gcsfuse) or by download (needs internet enabled and spends
   wall-clock against the 9-hour session).
3. **2.5 TiB free of 8.0 TiB is one observation**, on one container in `us-west1-c`, with 70% of the
   overlay already used by the image. It is not a guaranteed allocation and could differ per session.
4. **Reported capacity is not throughput.** A 1.0 PiB gcsfuse mount says nothing about read speed,
   and for a 302 GiB tensor streamed per epoch, **bandwidth is the untested question** — plausibly
   the real constraint, and not one this probe touched.
5. Measured on CPU deliberately, so none of the 30-hour GPU allowance was spent establishing it.
