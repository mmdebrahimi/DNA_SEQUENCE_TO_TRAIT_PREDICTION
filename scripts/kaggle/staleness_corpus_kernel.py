"""Kaggle kernel: run the benchmark-validated staleness auditor over the full CLAUDE.md corpus.

Reads the 110 (claim, artifact) pairs from the attached `dna-staleness-corpus` dataset and emits one
verdict per pair. The SYSTEM PROMPT and the user-turn shape are byte-identical to the ones the 10-item
benchmark scored (5/5 TP, 1 FP) -- if they drifted, the measured performance would not transfer, and the
whole point of running the benchmark first was to earn the right to trust this run.

Output is a FLAGGING pass, not a verdict on the docs: every `stale` flag is adjudicated by hand against
its artifact before any documentation changes. At the benchmark's measured 1-in-5 false-positive rate on
capability-shaped claims, that adjudication is the product, not an optional step.
"""
import json
import os

os.environ.setdefault("PYTHONUTF8", "1")
import torch
from transformers import AutoModelForCausalLM, AutoTokenizer

MODEL = "Qwen/Qwen3-8B"
MAX_NEW_TOKENS = 2500
# Total prompt+generation ceiling. RAISED 4200 -> 5500 after measuring how often it actually binds:
# at 4200 it clipped 107 of 110 items to a median 1292 generated tokens, against the 2500 every measured
# run used. That is not an outlier guard, it is a different configuration -- and 1292 sits right on the
# v1 cap of 1200 that produced a truncated, unparseable answer, so it would have re-created a failure
# mode already measured and fixed. At 5500 it binds on 9/110 (8%), which is the intended shape: clip the
# genuine outliers, leave the typical item on the measured 2500.
# Prompt tokens measured over the real corpus: median 2906, p90 3000, max 3003 (chars/3.3 estimate).
# STILL NOT YET VERIFIED end-to-end: no clean 110/110 run has completed with this value.
TOTAL_TOKEN_BUDGET = 5500

SYSTEM = """You audit documentation claims for STALENESS.

You are given a CLAIM made in a project's documentation, and an EXCERPT of the ARTIFACT that claim is
about. Decide whether the artifact still SUPPORTS the claim.

Answer exactly one of:
  stale     - the artifact CONTRADICTS the claim (e.g. the claim says work is "deferred"/"pending"/
              "blocked" but the artifact shows it was DONE, or the claim states a number/status the
              artifact disagrees with)
  supported - the artifact is consistent with the claim, OR the text merely MENTIONS a deferral that
              refers to something ELSE (a different sub-task), OR the text is itself a CORRECTION
              explaining that something was previously mis-stated
  unclear   - the excerpt genuinely does not settle it

CRITICAL: prose often mentions "deferred" or "pending" NEAR a file path without the claim being about that
file. Do not flag on keyword proximity. Flag only when the artifact's CONTENT contradicts the claim's
SUBSTANCE. A sentence that says "X was deferred, and here is the correction recording that it shipped" is
`supported` -- it is accurate text about a past error, not a stale claim.

ARTIFACT FACTS are evidence too, but they are the LAST rule, not the first. Apply them ONLY after every
exception above has been checked, and ONLY to a CAPABILITY claim.

A CAPABILITY claim asserts something is NOT BUILT ("X is deferred", "the browser is not implemented",
"this is still pending"). For those, and only those, ARTIFACT FACTS showing implemented code REFUTE the
claim -> `stale`, even if the excerpt's prose repeats "deferred" (a module often describes itself as "the
deferred X" because it is the thing that finally implemented X).

Implemented code does NOT make a claim stale in any of these cases -- each is `supported`:
  - FINDING claims. "X is infeasible", "the substrate does not clear the bar", "the cohort is too small".
    The script is the INSTRUMENT that produced the finding; its existence is what makes the finding
    trustworthy. Code existing SUPPORTS a finding claim. This is the single most common error.
  - CORRECTION text. "X was deferred, and here is the correction recording that it shipped." Already
    stated above; it OUTRANKS this rule, it is not overridden by it.
  - HISTORICAL / past-tense claims. "the gates this family needed and did not have", "until 2026-08-25
    this line said X". A past-tense statement about what was missing is not refuted by the thing now
    existing.

When you cannot tell whether a claim is CAPABILITY or FINDING, answer `supported`. A missed stale claim
costs one unnoticed line; a false flag costs a human the time to check it.

Reply with STRICT JSON only, no prose outside it:
{"verdict": "stale|supported|unclear", "evidence": "<one sentence quoting what decided it>"}"""

# DISCOVER the input rather than trusting a hardcoded mount path. The first run died with
# FileNotFoundError on /kaggle/input/dna-staleness-corpus/... -- a dataset that did not attach where
# expected looks exactly like a model failure in the results, which is the failure mode the benchmark
# avoided by inlining. Search, and if nothing is found, say so loudly instead of dying on a path.
import glob
cands = sorted(glob.glob("/kaggle/input/**/staleness_corpus_items.json", recursive=True))
print("available under /kaggle/input:", sorted(glob.glob("/kaggle/input/*")), flush=True)
if not cands:
    raise SystemExit("corpus dataset did not attach -- no staleness_corpus_items.json under /kaggle/input")
print("using:", cands[0], flush=True)
items = json.load(open(cands[0], encoding="utf-8"))
print(f"loaded {len(items)} pairs", flush=True)
print("device:", torch.cuda.get_device_name(0) if torch.cuda.is_available() else "CPU", flush=True)

tok = AutoTokenizer.from_pretrained(MODEL)
model = AutoModelForCausalLM.from_pretrained(MODEL, torch_dtype=torch.bfloat16, device_map="auto")
model.eval()

out = []
for n, it in enumerate(items, 1):
    user = (f"CLAIM (from project documentation):\n{it['claim']}\n\n"
            f"ARTIFACT the claim is about: {it['artifact']}\n"
            f"ARTIFACT FACTS: {it.get('facts','')}\n"
            f"ARTIFACT EXCERPT:\n{it['artifact_text']}\n\n"
            f"Does the artifact still support the claim?")
    msgs = [{"role": "system", "content": SYSTEM}, {"role": "user", "content": user}]
    text = tok.apply_chat_template(msgs, tokenize=False, add_generation_prompt=True)
    ids = tok([text], return_tensors="pt").to(model.device)
    # Bound the TOTAL sequence, not the input. The OOM is a per-item PEAK (one 3.94 GiB allocation), so
    # empty_cache() between items does not prevent it -- measured: the cache-clearing version still died
    # at the same item. Truncating the excerpt does prevent it but costs half the recall. So cap the
    # generation for the few longest prompts instead, which trades reasoning length (cheap: the P4 fix
    # showed ~1600-3000 tokens is typical) for input length (expensive: it is the evidence). The budget
    # must be set so this clips OUTLIERS only -- see the TOTAL_TOKEN_BUDGET note; a budget that binds on
    # nearly every item silently changes the configuration rather than guarding it.
    n_in = ids.input_ids.shape[-1]
    budget = max(600, min(MAX_NEW_TOKENS, TOTAL_TOKEN_BUDGET - n_in))
    with torch.no_grad():
        gen = model.generate(**ids, max_new_tokens=budget, do_sample=False)
    raw = tok.decode(gen[0][ids.input_ids.shape[-1]:], skip_special_tokens=True)
    # Free the KV cache per item rather than truncating the excerpt. The 6000->3000 cap "fixed" the
    # OOM by removing the evidence the model needs -- measured at HALF the recall (0.667 -> 0.333 on
    # the same prompt). Treat the memory symptom in memory.
    del gen, ids
    if torch.cuda.is_available():
        torch.cuda.empty_cache()
    out.append({"item_id": it["item_id"], "artifact": it["artifact"], "raw": raw})
    if n % 10 == 0 or n == len(items):
        # checkpoint as we go: a 110-item run that dies at item 90 must not lose 90 items of work
        json.dump(out, open("/kaggle/working/results.json", "w", encoding="utf-8"), indent=2)
        print(f"[{n}/{len(items)}] checkpointed", flush=True)

json.dump(out, open("/kaggle/working/results.json", "w", encoding="utf-8"), indent=2)
print(f"wrote {len(out)} verdicts to /kaggle/working/results.json", flush=True)
