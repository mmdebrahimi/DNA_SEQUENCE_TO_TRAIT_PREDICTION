"""Does relaxing the O-antigen coverage cut rescue more than it breaks? Decision rule fixed in advance.

WHY THE ORDER MATTERS HERE MORE THAN USUAL. Three causal claims about this cell were asserted and then
measured wrong in a row: phase-2 flagellin as the dominant defect, then DB coverage as the real
priority, then (nearly) "so lower the threshold". Each was reached by reasoning from the previous
result rather than by measuring. This run therefore fixes its verdict rule BEFORE reading any number,
so the third claim cannot be rationalised into existence the way the first two were.

WHAT IS ALREADY KNOWN. 21 of 200 isolates abstain with an unresolved O antigen. Re-probing 21 of them
permissively found 14 DO hit an O allele below the deployed 90 identity / 80 coverage cut, all 14
naming the CORRECT O group, at identity median 99.8 and coverage median 58.4 (max 78.9). So the
information is present. What is NOT known -- and is the entire point of this run -- is what admitting
it COSTS on the isolates that currently resolve correctly.

THE TRADE IS REAL AND IT IS NOT FREE. Abstention and error are different failures: this cell abstains
on 29.5% and is explicitly not deployed as a drop-in caller, so turning silence into confident error is
the wrong direction to be wrong in. A relaxed cut can (a) rescue a correct call from an abstention,
(b) leave a call unchanged, (c) turn a correct call into a wrong one, or (d) turn one wrong call into a
different wrong call. Only (a) is a gain, only (c) is a real loss.

PRE-REGISTERED DECISION RULE (frozen in `PREREGISTERED` below, written before the sweep ran):
  ADOPT     if rescued_correct >= 7 (a third of the 21) AND newly_wrong <= 2 AND net >= +5
  REJECT    if newly_wrong > 2, regardless of how much is rescued
  NO_CHANGE otherwise -- the honest default, since the deployed cut stays and the cell keeps abstaining

The asymmetry (7 rescues demanded, 2 new errors tolerated) is deliberate: an abstention is recoverable
by a human, a confident wrong serovar is not.

Also reports SELECTIVE-CLASSIFICATION metrics -- coverage, accuracy-conditional-on-coverage, and the
forced-call baseline -- defined here rather than after the fact, so "high accuracy after abstaining"
cannot become a slogan. A caller that abstains on 90% can post any accuracy.

Needs blastn + the cached assemblies. Writes wiki/salmserovar_threshold_tradeoff_<date>.json.
"""
from __future__ import annotations

import argparse
import collections
import json
import sys
from datetime import date as _date
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from dna_decode.salmserovar.equivalence import equivalent, load_formula_index  # noqa: E402
from dna_decode.salmserovar.runner import call_serovar  # noqa: E402

BLASTN = "C:/Users/Farshad/ncbi-blast/bin/blastn.exe"
DEPLOYED = {"identity_threshold": 90.0, "coverage_threshold": 80.0}

# FROZEN BEFORE THE SWEEP RAN. Do not tune these to the result.
PREREGISTERED = {
    "registered_before_run": True,
    "adopt_if": {"rescued_correct_min": 7, "newly_wrong_max": 2, "net_gain_min": 5},
    "reject_if": "newly_wrong > 2 regardless of rescues",
    "default": "NO_CHANGE",
    "asymmetry_rationale": ("an abstention is recoverable by a human; a confident wrong serovar is not, "
                            "so the bar demands many rescues and tolerates few new errors"),
}


def outcome(pred, label, idx) -> str:
    if not pred:
        return "abstain"
    ok, _ = equivalent(pred, label, idx)
    return "correct" if ok else "wrong"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--baseline", type=Path,
                    default=Path("D:/dna_decode_cache/salm_asm/results.jsonl"))
    ap.add_argument("--asm-root", type=Path, default=Path("D:/dna_decode_cache/salm_asm"))
    ap.add_argument("--db-dir", type=Path, default=ROOT / "data" / "salmserovar_db")
    ap.add_argument("--relaxed-coverage", type=float, default=40.0)
    ap.add_argument("--relaxed-identity", type=float, default=90.0,
                    help="identity stays at the deployed value -- only COVERAGE is relaxed, since the "
                         "measured sub-threshold hits were near-perfect on identity")
    ap.add_argument("--checkpoint", type=Path,
                    default=Path("D:/dna_decode_cache/salm_asm/relaxed_results.jsonl"))
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--blastn", default=BLASTN)
    ap.add_argument("--out", type=Path, default=ROOT / "wiki" /
                    f"salmserovar_threshold_tradeoff_{_date.today().isoformat()}.json")
    a = ap.parse_args()

    idx = load_formula_index(a.db_dir / "serovar_table.tsv")
    base = [json.loads(ln) for ln in a.baseline.read_text(encoding="utf-8").splitlines() if ln.strip()]
    base = [r for r in base if r.get("status") == "ok"]
    if a.limit:
        base = base[:a.limit]
    print(f"baseline isolates: {len(base)}   deployed cut {DEPLOYED['identity_threshold']}/"
          f"{DEPLOYED['coverage_threshold']} -> relaxed {a.relaxed_identity}/{a.relaxed_coverage}\n")

    done: dict[str, dict] = {}
    if a.checkpoint.exists():
        for ln in a.checkpoint.read_text(encoding="utf-8").splitlines():
            if ln.strip():
                r = json.loads(ln)
                done[r["asm_acc"]] = r

    fh = open(a.checkpoint, "a", encoding="utf-8")
    for n, r in enumerate(base, 1):
        acc = r["asm_acc"]
        if acc in done:
            continue
        rec = {"asm_acc": acc, "label": r.get("label")}
        try:
            fa = next((p for p in (a.asm_root / acc).glob("*.fna")), None)
            if fa is None:
                rec["status"] = "assembly_missing"
            else:
                call = call_serovar(fa, a.db_dir, identity_threshold=a.relaxed_identity,
                                    coverage_threshold=a.relaxed_coverage, blastn_bin=a.blastn)
                rec["status"] = "ok"
                rec["pred_serovar"] = call.get("serovar")
                rec["pred_formula"] = call.get("antigenic_formula") or call.get("formula")
        except Exception as e:                       # noqa: BLE001
            rec["status"] = f"error:{type(e).__name__}"
            rec["error"] = str(e)[:180]
        fh.write(json.dumps(rec) + "\n"); fh.flush()
        done[acc] = rec
        if n % 25 == 0:
            print(f"  [{n}/{len(base)}]", flush=True)
    fh.close()

    # ---- transition matrix: baseline outcome -> relaxed outcome -------------------------------
    trans: collections.Counter = collections.Counter()
    rescued, newly_wrong, examples = [], [], []
    scored = 0
    for r in base:
        acc, lab = r["asm_acc"], r.get("label")
        rel = done.get(acc)
        if not lab or not rel or rel.get("status") != "ok":
            continue
        scored += 1
        b = outcome(r.get("pred_serovar"), lab, idx)
        v = outcome(rel.get("pred_serovar"), lab, idx)
        trans[(b, v)] += 1
        if b == "abstain" and v == "correct":
            rescued.append({"asm": acc, "label": lab, "now": rel.get("pred_serovar")})
        if b == "correct" and v == "wrong":
            newly_wrong.append({"asm": acc, "label": lab, "was": r.get("pred_serovar"),
                                "now": rel.get("pred_serovar")})
        if b != v and len(examples) < 20:
            examples.append({"asm": acc, "label": lab, "baseline": b, "relaxed": v,
                             "was": r.get("pred_serovar"), "now": rel.get("pred_serovar")})

    def profile(get) -> dict:
        c = collections.Counter(get(r) for r in base if r.get("label"))
        tot = sum(c.values())
        resolved = c["correct"] + c["wrong"]
        return {"correct": c["correct"], "wrong": c["wrong"], "abstain": c["abstain"], "n": tot,
                "coverage": (resolved / tot) if tot else None,
                "accuracy_on_covered": (c["correct"] / resolved) if resolved else None,
                "accuracy_forced_call": (c["correct"] / tot) if tot else None}

    b_prof = profile(lambda r: outcome(r.get("pred_serovar"), r["label"], idx))
    r_prof = profile(lambda r: outcome((done.get(r["asm_acc"]) or {}).get("pred_serovar"),
                                       r["label"], idx)
                     if (done.get(r["asm_acc"]) or {}).get("status") == "ok" else "abstain")

    n_res, n_wrong = len(rescued), len(newly_wrong)
    net = n_res - n_wrong
    print(f"scored on {scored} isolates with a label and a relaxed call\n")
    print("  transition (baseline -> relaxed):")
    for (b, v), k in sorted(trans.items(), key=lambda kv: -kv[1]):
        mark = "  <- GAIN" if (b, v) == ("abstain", "correct") else (
            "  <- LOSS" if (b, v) == ("correct", "wrong") else "")
        print(f"    {b:<8} -> {v:<8} {k:>4}{mark}")
    print(f"\n  rescued (abstain->correct): {n_res}")
    print(f"  newly wrong (correct->wrong): {n_wrong}")
    print(f"  net: {net:+d}")
    print("\n  selective-classification profile:")
    for nm, p in (("deployed", b_prof), ("relaxed ", r_prof)):
        f = lambda x: "n/a" if x is None else f"{x:.4f}"                       # noqa: E731
        print(f"    {nm}: coverage={f(p['coverage'])} acc_on_covered={f(p['accuracy_on_covered'])} "
              f"acc_forced={f(p['accuracy_forced_call'])}")

    pr = PREREGISTERED["adopt_if"]
    if n_wrong > pr["newly_wrong_max"]:
        verdict = "REJECT"
        why = (f"{n_wrong} correct calls became wrong, above the pre-registered ceiling of "
               f"{pr['newly_wrong_max']}. Rejected regardless of the {n_res} rescued -- a confident "
               "wrong serovar is not recoverable by a human, an abstention is.")
    elif (n_res >= pr["rescued_correct_min"] and n_wrong <= pr["newly_wrong_max"]
          and net >= pr["net_gain_min"]):
        verdict = "ADOPT"
        why = (f"{n_res} abstentions became correct calls for {n_wrong} new errors (net {net:+d}), "
               "clearing every pre-registered bar.")
    else:
        verdict = "NO_CHANGE"
        why = (f"{n_res} rescued / {n_wrong} newly wrong / net {net:+d} does not clear the "
               f"pre-registered bar (>= {pr['rescued_correct_min']} rescued AND <= "
               f"{pr['newly_wrong_max']} newly wrong AND net >= {pr['net_gain_min']}). The deployed "
               "cut stays; the cell keeps abstaining, which is the safer failure.")
    print(f"\nVERDICT: {verdict}\n  {why}")

    out = {"schema": "salmserovar-threshold-tradeoff-v1", "date": _date.today().isoformat(),
           "question": "does relaxing the O-antigen coverage cut rescue more correct calls than it breaks?",
           "preregistered": PREREGISTERED,
           "deployed_thresholds": DEPLOYED,
           "relaxed_thresholds": {"identity_threshold": a.relaxed_identity,
                                  "coverage_threshold": a.relaxed_coverage},
           "n_scored": scored,
           "transition_matrix": {f"{b}->{v}": k for (b, v), k in trans.items()},
           "rescued_correct": n_res, "newly_wrong": n_wrong, "net": net,
           "rescued": rescued, "newly_wrong_detail": newly_wrong, "sample_changes": examples,
           "selective_classification": {"deployed": b_prof, "relaxed": r_prof},
           "verdict": verdict, "why": why,
           "honest_limits": [
               "Only the COVERAGE cut moves; identity stays at the deployed 90 because the measured "
               "sub-threshold hits were near-perfect on identity. A joint sweep was not run.",
               "One cohort (N=200, reference-lab-filtered) and one antigen-DB build.",
               "Equivalence uses the same notation+table rule as the original validation, applied "
               "identically to both threshold settings, so no leniency can favour either.",
               "A NO_CHANGE or REJECT verdict does NOT mean the O-antigen defect is unfixable -- it "
               "means THIS lever does not pay for itself. The allele-length hypothesis (11 of 14 "
               "sub-threshold hits were one O group) remains untested.",
           ]}
    a.out.write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")
    print(f"\nwrote {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
