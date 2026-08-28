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

ALSO CRITICAL: weigh ARTIFACT FACTS as evidence, not just the excerpt text. If a claim says work is
"deferred"/"pending"/"not built" and the ARTIFACT FACTS show implemented code exists, the claim is `stale`
-- EVEN IF the excerpt's own prose repeats the word "deferred". A module frequently describes ITSELF as
"the deferred X" because it is the thing that finally implemented X; that phrasing is about what it
implements, and its existence refutes the claim that X is still deferred.

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
    with torch.no_grad():
        gen = model.generate(**ids, max_new_tokens=MAX_NEW_TOKENS, do_sample=False)
    raw = tok.decode(gen[0][ids.input_ids.shape[-1]:], skip_special_tokens=True)
    out.append({"item_id": it["item_id"], "artifact": it["artifact"], "raw": raw})
    if n % 10 == 0 or n == len(items):
        # checkpoint as we go: a 110-item run that dies at item 90 must not lose 90 items of work
        json.dump(out, open("/kaggle/working/results.json", "w", encoding="utf-8"), indent=2)
        print(f"[{n}/{len(items)}] checkpointed", flush=True)

json.dump(out, open("/kaggle/working/results.json", "w", encoding="utf-8"), indent=2)
print(f"wrote {len(out)} verdicts to /kaggle/working/results.json", flush=True)
