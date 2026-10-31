#!/usr/bin/env python
"""Build the Kaggle T4 kernel for family G-D: the LoRA fine-tune of GENERator-v2-prokaryote-1.2b.

WHY A BUILDER, not a git clone. Same reason as the sampling-sweep sibling: the repo is private and
the published wheel does not ship `scripts/`, so the exact import closure is embedded as a base64
tar.gz and token-substituted into the kernel body. Token replacement, NOT str.format -- the body is
full of f-strings and escaping every brace is how the sweep builder broke the first time.

WHAT THE RUN DOES, and the ordering is the whole design:

    1  unpack sources, pip install peft, IMPORT EVERY MODULE          (seconds)
    2  fetch MG1655, rebuild the disjoint split, assert it reproduces  (~1 min, ~5 MB)
    3  download the model                                             (~4.8 GB)
    4  ONE LoRA training step as a shape check                        (seconds)
    5  ARM A: zero-shot baseline on THIS device                       (~120 generations)
    6  LoRA fine-tune, 3 epochs                                       (the actual lever)
    7  ARM B: fine-tuned, same eval windows, same seed                (~120 generations)
    8  verdict artifact

Every cheap guard sits BEFORE the expensive thing it protects. That is deliberate and it is this
project's own recorded lesson twice over: the sweep kernel once died on a missing import AFTER
paying for the full 4.8 GB download, and the local anchor runs died mid-generation with nothing
banked. Here an import error, a corpus drift, or a broken LoRA wiring all surface in the first two
minutes.

WHY THE BASELINE IS RE-MEASURED INSTEAD OF READ OFF THE 2026-10-07 ARTIFACT. That number (0.7706)
was produced on CPU; this is a T4, and the project has already measured that CUDA and CPU yield
DIFFERENT sampled sequences for this generator. So the comparison that carries the verdict is
internal to this run -- same device, dtype, eval windows, seed, and code path, with the weights as
the ONLY difference. The recorded number is a cross-check, not the comparator. Pre-registered in
`wiki/glm_generator_finetune_acceptance_bar.json`.

BOTH ARMS GO THROUGH THE SAME CLI. The fine-tuned model is merged and saved, then scored by invoking
`scripts/glm_generator_falsifier.py` exactly as the baseline arm does, with `--model <path>`. Two
different code paths for the two arms would make the delta a measurement of the code.

    uv run python scripts/kaggle/build_glm_finetune_kernel.py
    uv run python scripts/kaggle_push_poll.py push \
        scripts/kaggle/glm_finetune_kernel.py glm-generator-finetune --gpu
"""
from __future__ import annotations

import base64
import io
import tarfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent.parent

#: The import closure of `scripts/glm_generator_falsifier.py` PLUS the split script. Kept VERBATIM.
#: `glm_finetune_window_split.py` is included rather than reimplemented because it is the module
#: that OWNS the disjointness proof -- a second reader of the same windows is exactly how a solved
#: trap gets re-derived, and the guard would live in the copy that was not called.
FILES = [
    "dna_decode/__init__.py",
    "dna_decode/data/__init__.py",
    "dna_decode/data/annotations.py",
    "dna_decode/data/refseq.py",
    "dna_decode/glm/__init__.py",
    "dna_decode/glm/baselines.py",
    "dna_decode/glm/corpus.py",
    "dna_decode/glm/edit.py",
    "dna_decode/glm/falsifier.py",
    "dna_decode/glm/generate.py",
    "dna_decode/glm/knob_contract.py",
    "dna_decode/glm/sweep_gates.py",
    "scripts/glm_generator_falsifier.py",
    "scripts/glm_finetune_window_split.py",
    # NOT optional, and a .py-only payload is what broke the first two pushes. The split script
    # verifies the rebuilt corpus fingerprint AGAINST this baseline and REFUSES (rc=3,
    # "baseline ... absent; the fingerprint cannot be checked") when it is missing -- which is the
    # correct behaviour: a fingerprint it cannot compare is not a fingerprint it has checked.
    # Measured 2026-10-31 from the diagnostic's own traceback.
    "wiki/glm_generator_falsifier_REAL_2026-10-07.json",
]

KERNEL = r'''"""Kaggle T4 kernel: family G-D, LoRA fine-tune of GENERator-v2-prokaryote-1.2b-base.

THE QUESTION. Zero-shot, this generator is WORSE than a 3-mer Markov chain at telling its output
from real E. coli promoters: 0.7706 vs the null's 0.6679 at feature mode `both` (lower is better).
The 2026-10-31 sampling sweep established that NO temperature setting closes that gap. Fine-tuning
the weights is the remaining lever, and the family's refined goal names the bar verbatim: bring
distinguishability below 0.6679 on a comparison set the model was never trained on.

THE BAR IS FROZEN before this ran, in wiki/glm_generator_finetune_acceptance_bar.json, with all
four outcomes named (beats null / improves-but-short / no improvement / instrument invalid). This
kernel does not decide anything; it measures and stamps.

DISJOINTNESS IS THE INTEGRITY CRUX. A generator reproduces what it memorised, so training on even
one eval window would look like a WIN. The train set is derived BY EXCLUSION from the falsifier's
eval half by scripts/glm_finetune_window_split.py, which asserts zero overlap on sequence content.
This kernel re-derives it and REFUSES unless it reproduces the committed manifest exactly.

Everything is printed so the log itself is the audit trail.
"""
import base64
import io
import json
import os
import subprocess
import sys
import tarfile
import time
import traceback

os.environ.setdefault("PYTHONUTF8", "1")
OUT = "/kaggle/working"
REPO = OUT + "/repo"
HF = OUT + "/hf_cache"
#: The merged fine-tuned model goes to /tmp, NOT /kaggle/working, on purpose: at fp32 it is ~4.8 GB
#: and /kaggle/working is collected as the kernel OUTPUT. A 4.8 GB artifact would bloat the output
#: past anything retrievable and buy nothing -- the adapter is what is worth keeping, and it is tiny.
MERGED = "/tmp/glm_finetuned_merged"
ADAPTER = OUT + "/lora_adapter"

MODEL_ID = "GenerTeam/GENERator-v2-prokaryote-1.2b-base"

# ---- the recorded CPU baseline, used as a CROSS-CHECK only, never as the comparator ----
RECORDED = {
    "artifact": "wiki/glm_generator_falsifier_REAL_2026-10-07.json",
    "n": 120, "limit_natural": 1200, "bins": 6, "length": 150,
    "context": 300, "conditional": True, "natural_set": "promoter",
    "hf": {"kmer": 0.7237, "both": 0.7706},
    "markov_null": {"kmer": 0.5990, "both": 0.6679},
    "uniform": {"kmer": 0.9408, "both": 0.9641},
    "device_note": "CPU / GTX 860M host at fp32. CUDA yields DIFFERENT sampled sequences for this "
                   "generator, so this is NOT the comparator for a T4 run.",
}

# ---- the frozen bar (mirrors wiki/glm_generator_finetune_acceptance_bar.json) ----
BAR_BOTH = 0.6679
BAR_KMER = 0.5990

# ---- the frozen LoRA method ----
LORA = {"r": 16, "lora_alpha": 32, "lora_dropout": 0.05,
        "target_modules": ["q_proj", "k_proj", "v_proj", "o_proj"]}
LR = 1e-4
EPOCHS = 3
BATCH = 8
DTYPE = "float32"

# ---- what the committed split manifest says; a drift here REFUSES ----
EXPECT_SPLIT = {"corpus_total": 4340, "falsifier_fit": 600, "falsifier_eval": 600,
                "finetune_train": 3737}
EXPECT_FINGERPRINT = {"n_cds_total": 4340, "n_windows": 4340, "n_overlapping_cds": 2708,
                      "strand_counts": {"+": 2118, "-": 2222}}

#: DIAGNOSTIC MODE. When true the run STOPS after step 2 -- imports, pip, device report, reference
#: fetch and the disjoint-split check -- and never downloads the model. It exists to separate "my
#: code is broken in the Kaggle environment" from "no GPU was allocated to this session", which the
#: 2026-10-31 failure could not distinguish: the kernel ERRORed with an EMPTY failureMessage and
#: ZERO output files, so there was no log to read. A CPU diagnostic costs no GPU quota and answers
#: it in about two minutes.
DIAG = __DIAG__

PAYLOAD = "__PAYLOAD__"


def sh(cmd):
    print("$ " + cmd, flush=True)
    r = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if r.stdout:
        print(r.stdout[-4000:], flush=True)
    if r.returncode and r.stderr:
        print("STDERR:", r.stderr[-3000:], flush=True)
    return r.returncode


def unpack():
    os.makedirs(REPO, exist_ok=True)
    with tarfile.open(fileobj=io.BytesIO(base64.b64decode(PAYLOAD.encode())), mode="r:gz") as tf:
        tf.extractall(REPO)
    print("unpacked sources:", flush=True)
    for root, _d, fs in os.walk(REPO):
        for f in sorted(fs):
            p = os.path.join(root, f)
            print("   ", os.path.relpath(p, REPO), os.path.getsize(p), "B", flush=True)
    sys.path.insert(0, REPO)
    sys.path.insert(0, REPO + "/scripts")


def step1_imports():
    """STEP 1 -- every import, BEFORE the 4.8 GB download. The sweep kernel once died here, after."""
    print("\n=== STEP 1: imports (cheapest guard, runs first) ===", flush=True)
    sh(sys.executable + " -m pip install -q peft")
    import peft
    import torch
    import transformers
    print("peft %s  transformers %s  torch %s" % (peft.__version__, transformers.__version__,
                                                  torch.__version__), flush=True)
    import glm_generator_falsifier as cli          # noqa: F401
    import glm_finetune_window_split as split      # noqa: F401
    from dna_decode.glm import falsifier, generate, knob_contract, sweep_gates  # noqa: F401
    print("every module imported OK", flush=True)
    return cli, split


def device_report():
    import torch
    if not torch.cuda.is_available():
        print("!! NO CUDA. This kernel is pinned to NvidiaTeslaT4; a CPU session means the "
              "machine_shape was not honoured.", flush=True)
        return {"cuda": False}
    name = torch.cuda.get_device_name(0)
    cc = torch.cuda.get_device_capability(0)
    vram = torch.cuda.get_device_properties(0).total_memory / 2 ** 30
    d = {"cuda": True, "device_name": name, "compute_capability": "%d.%d" % cc,
         "vram_gib": round(vram, 2), "torch": torch.__version__, "dtype_used": DTYPE,
         "note": "fp32 for BOTH arms. CC>=7.5 would make fp16 fast, but fp16 sampling can change "
                 "generated sequences and the weights must be the only difference between arms."}
    print("device: %s  CC %s  %.1f GiB" % (name, d["compute_capability"], vram), flush=True)
    return d


def step2_corpus(split):
    """STEP 2 -- rebuild the disjoint split and REFUSE on any drift. Before the model download."""
    print("\n=== STEP 2: reference + disjoint split ===", flush=True)
    from dna_decode.data.refseq import download_genome
    dest = REPO + "/data/cache/refseq"
    os.makedirs(dest, exist_ok=True)
    t0 = time.time()
    download_genome("GCF_000005845.2", dest)
    print("reference fetched in %.0fs" % (time.time() - t0), flush=True)

    # Run the split script's OWN entry point and read its manifest back, rather than calling into
    # its internals. Its `main()` is what also runs `verify()` against the baseline fingerprint,
    # so invoking it gets the drift check for free instead of reimplementing it here.
    rc = subprocess.run([sys.executable, REPO + "/scripts/glm_finetune_window_split.py"],
                        cwd=REPO, capture_output=True, text=True)
    print(rc.stdout[-3000:], flush=True)
    if rc.returncode:
        print("STDERR:", rc.stderr[-3000:], flush=True)
        raise SystemExit("split script failed rc=%d" % rc.returncode)
    man = json.load(open(REPO + "/wiki/glm_finetune_window_split.json"))

    counts = man["counts"]
    fp = man["corpus_fingerprint"]
    print("counts:", json.dumps(counts), flush=True)
    print("fingerprint:", json.dumps(fp), flush=True)
    print("overlap:", man["disjointness"]["train_vs_eval_overlap"], flush=True)

    for k, v in EXPECT_SPLIT.items():
        if counts.get(k) != v:
            raise SystemExit("SPLIT DRIFT: counts[%s]=%r, expected %r. The eval set has moved "
                             "underneath the comparison; refusing." % (k, counts.get(k), v))
    for k, v in EXPECT_FINGERPRINT.items():
        if fp.get(k) != v:
            raise SystemExit("CORPUS DRIFT: fingerprint[%s]=%r, expected %r. NCBI re-annotates "
                             "under the same accession; refusing." % (k, fp.get(k), v))
    if man["disjointness"]["train_vs_eval_overlap"] != 0:
        raise SystemExit("LEAKAGE: train/eval overlap is nonzero; refusing.")
    print("split verified: %d train / %d held-out / 0 overlap"
          % (counts["finetune_train"], counts["falsifier_eval"]), flush=True)
    return man


def train_sequences(split):
    """The 3,737 train windows as plain strings, from the module that OWNS the split.

    `load_windows` -> `split` -> `_seq` is that module's real API. Reimplementing the exclusion
    here is exactly the second-reader mistake its own docstring warns about: the disjointness
    guard lives in `split()`, and a copy of the index arithmetic would not carry it.
    """
    wins, stats = split.load_windows()
    s = split.split(wins)
    seqs = [split._seq(w) for w in s["train"]]
    evalset = s["eval_seqs"]
    bad = [q for q in seqs if q in evalset]
    if bad:
        raise SystemExit("LEAKAGE: %d train sequences are in the eval set; refusing." % len(bad))
    print("train sequences re-derived: %d, eval held out: %d, overlap 0 (re-asserted here)"
          % (len(seqs), len(evalset)), flush=True)
    return seqs, stats


def step3_model():
    print("\n=== STEP 3: model download ===", flush=True)
    import torch
    from transformers import AutoModelForCausalLM, AutoTokenizer
    os.makedirs(HF, exist_ok=True)
    t0 = time.time()
    tok = AutoTokenizer.from_pretrained(MODEL_ID, cache_dir=HF, trust_remote_code=True)
    model = AutoModelForCausalLM.from_pretrained(
        MODEL_ID, cache_dir=HF, dtype=getattr(torch, DTYPE), trust_remote_code=True)
    print("model loaded in %.0fs" % (time.time() - t0), flush=True)
    # The remote modeling code registers its bp buffers by RE-FETCHING its own tokenizer with no
    # cache_dir; when that fails it only PRINTS a warning and generate() later raises. Wire it
    # explicitly, exactly as dna_decode/glm/generate.py does.
    setup = getattr(model, "setup_tokenizer", None)
    if callable(setup) and not hasattr(model, "_bp_base_index"):
        setup(tok)
        print("setup_tokenizer wired explicitly", flush=True)
    return tok, model


def step4_lora_shapecheck(tok, model, seqs):
    """STEP 4 -- ONE training step, to prove the wiring before paying for 3 epochs."""
    print("\n=== STEP 4: LoRA shape check (one step) ===", flush=True)
    import torch
    from peft import LoraConfig, get_peft_model
    cfg = LoraConfig(task_type="CAUSAL_LM", **LORA)
    peft_model = get_peft_model(model, cfg)
    peft_model.print_trainable_parameters()
    dev = "cuda" if torch.cuda.is_available() else "cpu"
    peft_model.to(dev)
    batch = tok(seqs[:2], return_tensors="pt", padding=True, truncation=True, max_length=256)
    batch = {k: v.to(dev) for k, v in batch.items()}
    batch["labels"] = batch["input_ids"].clone()
    peft_model.train()
    out = peft_model(**batch)
    loss = out.loss
    print("one forward: loss=%.4f" % float(loss), flush=True)
    loss.backward()
    g = [p for p in peft_model.parameters() if p.requires_grad and p.grad is not None]
    if not g:
        raise SystemExit("LORA DEAD: no trainable parameter received a gradient. The adapter is "
                         "not attached to anything -- refusing before spending 3 epochs.")
    print("gradients reached %d trainable tensors -> adapter is live" % len(g), flush=True)
    peft_model.zero_grad(set_to_none=True)
    return peft_model


def run_arm(cli, label, model_path, device, outfile):
    """Score one arm through the SAME falsifier CLI. Only `--model` differs between arms."""
    print("\n=== ARM %s: %s ===" % (label, model_path), flush=True)
    argv = [
        "--model", model_path, "--device", device, "--conditional",
        "--context", str(RECORDED["context"]),
        "--natural-set", RECORDED["natural_set"],
        "--n", str(RECORDED["n"]), "--length", str(RECORDED["length"]),
        "--limit-natural", str(RECORDED["limit_natural"]),
        "--modes", "kmer,both", "--markov-order", "3",
        "--bins", str(RECORDED["bins"]),
        "--cache-dir", HF,
        "--out", outfile,
    ]
    print("argv: " + " ".join(argv), flush=True)
    t0 = time.time()
    rc = cli.main(argv)
    print("arm %s rc=%d in %.0f min" % (label, rc, (time.time() - t0) / 60.0), flush=True)
    d = json.load(open(outfile))
    got = {}
    for r in d.get("results", []):
        if r.get("candidate", "").startswith("hf:") or "GENER" in r.get("candidate", ""):
            got[r["feature_mode"]] = r["distinguishability_auroc_LOWER_IS_BETTER"]
        if r.get("candidate") == "markov-k3":
            got.setdefault("_null", {})[r["feature_mode"]] = \
                r["distinguishability_auroc_LOWER_IS_BETTER"]
    print("arm %s -> %s" % (label, json.dumps(got)), flush=True)
    return d, got


def step6_finetune(peft_model, tok, seqs):
    print("\n=== STEP 6: LoRA fine-tune, %d epochs on %d windows ===" % (EPOCHS, len(seqs)),
          flush=True)
    import torch
    dev = "cuda" if torch.cuda.is_available() else "cpu"
    opt = torch.optim.AdamW([p for p in peft_model.parameters() if p.requires_grad], lr=LR)
    peft_model.train()
    curve = []
    t0 = time.time()
    for ep in range(EPOCHS):
        tot, nb = 0.0, 0
        for i in range(0, len(seqs), BATCH):
            chunk = seqs[i:i + BATCH]
            b = tok(chunk, return_tensors="pt", padding=True, truncation=True, max_length=256)
            b = {k: v.to(dev) for k, v in b.items()}
            b["labels"] = b["input_ids"].clone()
            loss = peft_model(**b).loss
            loss.backward()
            opt.step()
            opt.zero_grad(set_to_none=True)
            tot += float(loss); nb += 1
            if nb % 50 == 0:
                print("  ep%d %d/%d  loss %.4f" % (ep, i, len(seqs), tot / nb), flush=True)
        curve.append({"epoch": ep, "mean_loss": round(tot / max(nb, 1), 5)})
        print("epoch %d mean loss %.4f  (%.1f min elapsed)"
              % (ep, tot / max(nb, 1), (time.time() - t0) / 60.0), flush=True)
    print("loss curve:", json.dumps(curve), flush=True)
    return curve


def verdict(base, ft, null_ok):
    """The FROZEN reading. Ordered; first match wins. No threshold is computed here."""
    if not null_ok:
        return "D_INSTRUMENT_INVALID", "the Markov null did not reproduce"
    b, f = base.get("both"), ft.get("both")
    if b is None or f is None:
        return "D_INSTRUMENT_INVALID", "an arm produced no mode-both number"
    if f < BAR_BOTH:
        return "A_BEATS_NULL", "fine-tuned %.4f < null %.4f" % (f, BAR_BOTH)
    if f < b:
        return "B_IMPROVES_BUT_BAR_UNMET", \
            "fine-tuned %.4f < baseline %.4f but >= null %.4f" % (f, b, BAR_BOTH)
    return "C_NO_IMPROVEMENT", "fine-tuned %.4f >= baseline %.4f" % (f, b)


def main():
    t_run = time.time()
    unpack()
    cli, split = step1_imports()
    dev = device_report()
    device = "cuda" if dev.get("cuda") else "cpu"
    man = step2_corpus(split)

    seqs, _stats = train_sequences(split)
    print("train sequences: %d, first len %d" % (len(seqs), len(seqs[0])), flush=True)

    if DIAG:
        json.dump({"record": "glm-finetune-diagnostic-v1",
                   "reached": "step2-complete",
                   "device": dev, "counts": man["counts"],
                   "fingerprint": man["corpus_fingerprint"],
                   "n_train_sequences": len(seqs),
                   "note": "DIAGNOSTIC ONLY. Stops before the model download. A COMPLETE status "
                           "here means the sources, pip install, reference fetch and disjoint "
                           "split all work in the Kaggle environment, so a GPU run that ERRORs "
                           "with no log is environmental (allocation / quota), not a code fault."},
                  open(OUT + "/glm_finetune_diagnostic.json", "w"), indent=2)
        print("")
        print("=== DIAGNOSTIC COMPLETE: stopping before the model download ===", flush=True)
        print("wrote " + OUT + "/glm_finetune_diagnostic.json", flush=True)
        return

    tok, model = step3_model()

    # ARM A first, on the PRISTINE model, before any adapter touches it.
    base_json, base = run_arm(cli, "A-zeroshot", MODEL_ID, device, OUT + "/arm_a_zeroshot.json")

    peft_model = step4_lora_shapecheck(tok, model, seqs)
    curve = step6_finetune(peft_model, tok, seqs)

    print("\n=== merging adapter + saving ===", flush=True)
    peft_model.save_pretrained(ADAPTER)
    merged = peft_model.merge_and_unload()
    merged.save_pretrained(MERGED, safe_serialization=True)
    tok.save_pretrained(MERGED)
    print("merged model at %s" % MERGED, flush=True)

    ft_json, ft = run_arm(cli, "B-finetuned", MERGED, device, OUT + "/arm_b_finetuned.json")

    nb = base.get("_null", {}); nf = ft.get("_null", {})
    null_ok = (abs(nb.get("both", 0) - BAR_BOTH) < 1e-6 and
               abs(nf.get("both", 0) - BAR_BOTH) < 1e-6 and
               abs(nb.get("kmer", 0) - BAR_KMER) < 1e-6 and
               abs(nf.get("kmer", 0) - BAR_KMER) < 1e-6)
    print("null reproduced in both arms: %s  (A=%s B=%s)"
          % (null_ok, json.dumps(nb), json.dumps(nf)), flush=True)

    v, why = verdict(base, ft, null_ok)
    rec = {
        "record": "glm-generator-finetune-verdict-v1",
        "date": time.strftime("%Y-%m-%d"),
        "acceptance_bar": "wiki/glm_generator_finetune_acceptance_bar.json",
        "verdict": v,
        "why": why,
        "bar": {"mode": "both", "threshold": BAR_BOTH, "direction": "lower is better"},
        "arm_a_zeroshot_on_device": base,
        "arm_b_finetuned": ft,
        "delta_both": (None if (base.get("both") is None or ft.get("both") is None)
                       else round(ft["both"] - base["both"], 4)),
        "markov_null_reproduced": null_ok,
        "recorded_cpu_crosscheck": RECORDED,
        "lora": LORA, "lr": LR, "epochs": EPOCHS, "batch": BATCH, "dtype": DTYPE,
        "loss_curve": curve,
        "split": {"counts": man["counts"], "fingerprint": man["corpus_fingerprint"],
                  "train_vs_eval_overlap": man["disjointness"]["train_vs_eval_overlap"]},
        "device": dev,
        "wall_minutes": round((time.time() - t_run) / 60.0, 1),
        "honest_limits": [
            "The comparator is the ON-DEVICE zero-shot arm measured in THIS run, not the recorded "
            "CPU 0.7706 -- CUDA and CPU yield different sampled sequences for this generator.",
            "Clearing the bar would be NECESSARY, not sufficient: a k-mer + positional "
            "discriminator measures LOCAL composition, and published audits name long-range "
            "organisation as the real failure mode.",
            "n=120 generated sequences per arm, one seed. This is the recorded geometry, not a "
            "powered comparison of two generators.",
            "3,737 windows (~560 kb) is a small corpus for a 1.2B model; the ledger rates 'is that "
            "enough' at LOW confidence. A miss does not separate 'fine-tuning does not work' from "
            "'not enough data'.",
        ],
    }
    json.dump(rec, open(OUT + "/glm_generator_finetune_verdict.json", "w"), indent=2)
    print("\n==== VERDICT: %s ====" % v, flush=True)
    print(why, flush=True)
    print("wrote " + OUT + "/glm_generator_finetune_verdict.json", flush=True)


if __name__ == "__main__":
    # WHY THE DIAGNOSTIC SWALLOWS ITS OWN FAILURE AND EXITS 0. Kaggle exposes a kernel's output
    # files ONLY when its status is COMPLETE. A kernel that ERRORs returns zero files and an EMPTY
    # failureMessage through the API, so its traceback is unreachable without a browser -- which is
    # why the 2026-10-31 failures were indistinguishable from an infrastructure problem and cost a
    # wrong diagnosis before a right one. So in DIAG mode the run always COMPLETES and reports what
    # happened IN AN ARTIFACT. A diagnostic that cannot report its own failure is useless.
    if DIAG:
        try:
            main()
        except BaseException:
            tb = traceback.format_exc()
            print(tb, flush=True)
            try:
                json.dump({"record": "glm-finetune-diagnostic-v1",
                           "reached": "FAILED",
                           "traceback": tb.splitlines()[-40:],
                           "note": "The diagnostic itself failed. Because it exits 0 anyway, this "
                                   "artifact is retrievable -- an ERROR status would have hidden "
                                   "the traceback entirely."},
                          open(OUT + "/glm_finetune_diagnostic.json", "w"), indent=2)
            except BaseException:
                pass
        raise SystemExit(0)
    try:
        main()
    except SystemExit as e:
        print("REFUSED: %s" % e, flush=True)
        raise
    except Exception:
        traceback.print_exc()
        raise
'''


def main(argv=None) -> int:
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--diagnostic", action="store_true",
                    help="emit the CPU diagnostic variant that stops before the model download")
    a = ap.parse_args(argv)
    buf = io.BytesIO()
    with tarfile.open(fileobj=buf, mode="w:gz") as tf:
        for rel in FILES:
            p = ROOT / rel
            if not p.is_file():
                raise SystemExit("missing source for the payload: %s" % rel)
            tf.add(p, arcname=rel)
    payload = base64.b64encode(buf.getvalue()).decode()
    body = KERNEL.replace("__DIAG__", "True" if a.diagnostic else "False")
    body = body.replace("__PAYLOAD__", payload)
    # THE EMITTED KERNEL MUST PARSE, and this guard exists because its absence cost two Kaggle
    # sessions on 2026-10-31. The builder's OWN syntax says NOTHING about the kernel's: the body is
    # a raw triple-quoted string, so a literal newline inside a `print("...")` is perfectly legal
    # *in the builder* and produces an unterminated string literal *in the emitted file*. That is
    # exactly what happened -- the builder parsed, the kernel did not, and Kaggle reported ERROR
    # with an EMPTY failureMessage and ZERO output files, which is indistinguishable from an
    # infrastructure failure. Checking the artifact rather than its generator is the lesson.
    import ast
    try:
        ast.parse(body)
    except SyntaxError as e:
        raise SystemExit("REFUSING to write a kernel that does not parse: line %d: %s -- %s"
                         % (e.lineno, e.msg, (e.text or "").rstrip()))
    name = "glm_finetune_diag_kernel.py" if a.diagnostic else "glm_finetune_kernel.py"
    out = ROOT / "scripts" / "kaggle" / name
    out.write_text(body, encoding="utf-8")
    print("wrote %s  (%d KB, payload %d KB, %d files)"
          % (out, len(body) // 1024, len(payload) // 1024, len(FILES)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
