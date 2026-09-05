"""Is coverage 40 the right cut, or just the first one tried? Sweep — without selecting on the test set.

WHY THIS EXISTS. The coverage cut was moved 80 -> 40 against a pre-registered bar and it cleared that bar
decisively (36 rescued / 1 new error). But 40 was the ONLY value tried. The shipped artifact says so:
"40 is not shown to be optimal -- it clears a bar". This closes that debt.

THE TRAP THIS IS DESIGNED AGAINST. Sweeping several cuts and keeping the best one on the SAME 200
isolates is selection on the test set: the winner's margin is then part noise, and the reported number
is optimistic by construction. That is the same class of error as the unblinded fix earlier in this
track, which shrank from +0.155 to +0.106 when replicated.

TWO GUARDS, BOTH FIXED BEFORE THE SWEEP RAN:

  1. SPLIT. The cohort is partitioned deterministically by a hash of the accession into SELECT and
     CONFIRM halves. The cut is chosen on SELECT alone; the reported number comes from CONFIRM, which
     the choice never touched.
  2. THE SELECTION RULE IS NOT ARGMAX. Argmax picks the luckiest cut. The rule instead takes the MOST
     CONSERVATIVE (highest-coverage) cut that clears the adopt bar on SELECT -- so ties and noise
     resolve toward caution, and a cut only wins by being sufficient, not by being maximal.

ONE BLASTN PASS, NOT ONE PER CUT. Every cut is re-derived offline from a single permissive blastn run
per isolate, so the candidate cuts are compared on identical alignments and nothing varies but the
threshold. It is also ~N times cheaper.

Needs blastn + cached assemblies. Writes wiki/salmserovar_threshold_sweep_<date>.json.
"""
from __future__ import annotations

import argparse
import collections
import hashlib
import json
import sys
from datetime import date as _date
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from dna_decode.salmserovar.equivalence import equivalent, load_formula_index  # noqa: E402
from dna_decode.salmserovar.runner import (  # noqa: E402
    _best_per_axis, load_serovar_table, parse_axis_antigen,
)
from dna_decode.typing.blast_caller import call_alleles  # noqa: E402

BLASTN = "C:/Users/Farshad/ncbi-blast/bin/blastn.exe"
PROBE = {"identity": 85.0, "coverage": 20.0}          # permissive superset; every cut filters from this
CANDIDATE_COVERAGE = [70.0, 60.0, 50.0, 40.0, 30.0]   # identity held at 90 -- never the failing axis
IDENTITY = 90.0

PREREGISTERED = {
    "registered_before_sweep": True,
    "split": "deterministic md5(accession) parity -> SELECT / CONFIRM; the cut is chosen on SELECT only",
    "selection_rule": ("the MOST CONSERVATIVE (highest-coverage) candidate that clears the adopt bar on "
                       "SELECT -- explicitly NOT argmax, so noise and ties resolve toward caution"),
    "adopt_bar_on_select": {"net_gain_min": 1, "newly_wrong_max": 2},
    "reported_number_comes_from": "CONFIRM only",
    "why": ("choosing the best of several cuts on the same isolates that report the result is selection "
            "on the test set; the split and the non-argmax rule are what make the reported number honest"),
}


def resolve(per_allele: dict, table: dict, ident: float, cov: float) -> str | None:
    """Re-derive the serovar call at one (identity, coverage) cut from a single blastn pass."""
    filtered = {}
    for aid, hit in per_allele.items():
        h = dict(hit)
        h["called"] = (hit["percent_identity"] >= ident and hit["percent_coverage"] >= cov)
        filtered[aid] = h
    ax = _best_per_axis(filtered)
    o = ax.get("O", {}).get("antigen")
    h1 = ax.get("H1", {}).get("antigen")
    h2 = ax.get("H2", {}).get("antigen")
    if not (o and h1):
        return None
    sv = table.get((o, h1, h2 or "-")) or table.get((o, h1, h2 or ""))
    if sv is None:
        cands = {v for (to, th1, _), v in table.items() if to == o and th1 == h1}
        sv = next(iter(cands)) if len(cands) == 1 else None
    return sv


def half(acc: str) -> str:
    return "SELECT" if int(hashlib.md5(acc.encode()).hexdigest(), 16) % 2 == 0 else "CONFIRM"


def score(rows, cov, table, idx) -> dict:
    c = collections.Counter()
    for r in rows:
        pred = resolve(r["per_allele"], table, IDENTITY, cov)
        if not pred:
            c["abstain"] += 1
        else:
            ok, _ = equivalent(pred, r["label"], idx)
            c["correct" if ok else "wrong"] += 1
    n = sum(c.values())
    res = c["correct"] + c["wrong"]
    # A Counter omits zero-count keys, so a cut that never abstains would produce a dict MISSING
    # `abstain` entirely -- which reads downstream as absent data rather than as a zero.
    fixed = {k: c[k] for k in ("correct", "wrong", "abstain")}
    return {"coverage_cut": cov, "n": n, **fixed,
            "coverage_rate": (res / n) if n else None,
            "accuracy_on_covered": (c["correct"] / res) if res else None,
            "correct_per_isolate": (c["correct"] / n) if n else None}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--baseline", type=Path,
                    default=Path("D:/dna_decode_cache/salm_asm/results.jsonl"))
    ap.add_argument("--asm-root", type=Path, default=Path("D:/dna_decode_cache/salm_asm"))
    ap.add_argument("--db-dir", type=Path, default=ROOT / "data" / "salmserovar_db")
    ap.add_argument("--hits-cache", type=Path,
                    default=Path("D:/dna_decode_cache/salm_asm/sweep_hits.jsonl"))
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--blastn", default=BLASTN)
    ap.add_argument("--out", type=Path, default=ROOT / "wiki" /
                    f"salmserovar_threshold_sweep_{_date.today().isoformat()}.json")
    a = ap.parse_args()

    table = load_serovar_table(a.db_dir / "serovar_table.tsv")
    idx = load_formula_index(a.db_dir / "serovar_table.tsv")
    base = [json.loads(ln) for ln in a.baseline.read_text(encoding="utf-8").splitlines() if ln.strip()]
    base = [r for r in base if r.get("status") == "ok" and r.get("label")]
    if a.limit:
        base = base[:a.limit]

    cached: dict[str, dict] = {}
    if a.hits_cache.exists():
        for ln in a.hits_cache.read_text(encoding="utf-8").splitlines():
            if ln.strip():
                r = json.loads(ln)
                cached[r["asm_acc"]] = r
    print(f"{len(base)} labelled isolates | {len(cached)} blastn passes cached")

    fh = open(a.hits_cache, "a", encoding="utf-8")
    for n, r in enumerate(base, 1):
        acc = r["asm_acc"]
        if acc in cached:
            continue
        fa = next((p for p in (a.asm_root / acc).glob("*.fna")), None)
        if fa is None:
            continue
        try:
            res = call_alleles(fa, a.db_dir / "salmonella_antigens.fasta",
                               identity_threshold=PROBE["identity"],
                               coverage_threshold=PROBE["coverage"],
                               blastn_bin=a.blastn, timeout=600)
        except Exception as e:                       # noqa: BLE001
            print(f"  skip {acc}: {type(e).__name__}", flush=True)
            continue
        if res.get("status") != "ok":
            continue
        # keep only antigen-parseable alleles -- the rest can never contribute to a formula
        pa = {k: {"percent_identity": v["percent_identity"],
                  "percent_coverage": v["percent_coverage"], "called": True}
              for k, v in res["per_allele"].items() if parse_axis_antigen(k)}
        rec = {"asm_acc": acc, "label": r["label"], "per_allele": pa}
        fh.write(json.dumps(rec) + "\n"); fh.flush()
        cached[acc] = rec
        if n % 25 == 0:
            print(f"  [{n}/{len(base)}]", flush=True)
    fh.close()

    rows = [cached[r["asm_acc"]] for r in base if r["asm_acc"] in cached]
    sel = [r for r in rows if half(r["asm_acc"]) == "SELECT"]
    con = [r for r in rows if half(r["asm_acc"]) == "CONFIRM"]
    if not sel or not con:
        print("split produced an empty half", file=sys.stderr)
        return 3
    print(f"\nsplit: SELECT {len(sel)} | CONFIRM {len(con)}\n")

    base_sel = score(sel, 80.0, table, idx)
    sweep_sel = [score(sel, c, table, idx) for c in CANDIDATE_COVERAGE]
    print("SELECT half (choice is made here, and only here):")
    print(f"  {'cut':>5} {'correct':>8} {'wrong':>6} {'abstain':>8} {'cov':>7} {'acc|cov':>8} {'net':>5}")
    print(f"  {80.0:>5} {base_sel['correct']:>8} {base_sel['wrong']:>6} {base_sel['abstain']:>8} "
          f"{base_sel['coverage_rate']:>7.3f} {base_sel['accuracy_on_covered']:>8.3f}  (deployed-before)")
    eligible = []
    for s in sweep_sel:
        net = s["correct"] - base_sel["correct"]
        newly_wrong = max(0, s["wrong"] - base_sel["wrong"])
        s["net_vs_80"], s["newly_wrong_vs_80"] = net, newly_wrong
        ok = (net >= PREREGISTERED["adopt_bar_on_select"]["net_gain_min"]
              and newly_wrong <= PREREGISTERED["adopt_bar_on_select"]["newly_wrong_max"])
        s["clears_bar_on_select"] = ok
        if ok:
            eligible.append(s)
        print(f"  {s['coverage_cut']:>5} {s['correct']:>8} {s['wrong']:>6} {s['abstain']:>8} "
              f"{s['coverage_rate']:>7.3f} {s['accuracy_on_covered']:>8.3f} {net:>+5}"
              f"{'  CLEARS' if ok else ''}")

    if not eligible:
        print("\nno candidate clears the bar on SELECT -- keeping the deployed cut", file=sys.stderr)
        chosen = None
    else:
        chosen = max(eligible, key=lambda s: s["coverage_cut"])   # MOST CONSERVATIVE, not argmax
        argmax = max(eligible, key=lambda s: s["correct"])
        print(f"\n  selection rule picks the MOST CONSERVATIVE clearing cut: {chosen['coverage_cut']}")
        if argmax["coverage_cut"] != chosen["coverage_cut"]:
            print(f"  (argmax on SELECT would have picked {argmax['coverage_cut']} -- deliberately NOT used)")

    out = {"schema": "salmserovar-threshold-sweep-v1", "date": _date.today().isoformat(),
           "question": "is coverage 40 the right cut, or just the first one tried?",
           "preregistered": PREREGISTERED, "identity_held_at": IDENTITY,
           "candidates": CANDIDATE_COVERAGE,
           "split": {"select_n": len(sel), "confirm_n": len(con)},
           "select_half": {"baseline_80": base_sel, "sweep": sweep_sel},
           "chosen_cut": (chosen["coverage_cut"] if chosen else None)}

    if chosen:
        c_base = score(con, 80.0, table, idx)
        c_chosen = score(con, chosen["coverage_cut"], table, idx)
        c_deployed = score(con, 40.0, table, idx)
        out["confirm_half"] = {"baseline_80": c_base, "chosen": c_chosen, "deployed_40": c_deployed}
        gain = c_chosen["correct"] - c_base["correct"]
        wrong_delta = c_chosen["wrong"] - c_base["wrong"]
        out["confirm_gain_vs_80"] = gain
        out["confirm_newly_wrong_vs_80"] = wrong_delta
        print(f"\nCONFIRM half (never touched by the choice), n={len(con)}:")
        for nm, s in (("cut 80 (old)", c_base), (f"cut {chosen['coverage_cut']} (chosen)", c_chosen),
                      ("cut 40 (deployed)", c_deployed)):
            print(f"  {nm:<22} correct={s['correct']:>3} wrong={s['wrong']:>3} abstain={s['abstain']:>3} "
                  f"cov={s['coverage_rate']:.3f} acc|cov={s['accuracy_on_covered']:.3f}")
        same = chosen["coverage_cut"] == 40.0
        if gain > 0 and wrong_delta <= 2:
            verdict = "CONFIRMED_KEEP_40" if same else "CONFIRMED_CHANGE_TO_CHOSEN"
            why = (f"the chosen cut {chosen['coverage_cut']} holds up on the untouched half "
                   f"(+{gain} correct, {wrong_delta:+d} wrong vs the old 80). "
                   + ("It IS the deployed value, so nothing changes -- 40 is now shown to be a "
                      "defensible choice, not merely the first one tried."
                      if same else
                      f"It differs from the deployed 40, so the deployed cut should move to "
                      f"{chosen['coverage_cut']}."))
        else:
            verdict = "NOT_CONFIRMED_ON_HELDOUT"
            why = (f"the cut chosen on SELECT does not hold on CONFIRM (+{gain} correct, "
                   f"{wrong_delta:+d} wrong) -- evidence the SELECT margin was partly noise.")
        out["verdict"], out["why"] = verdict, why
        print(f"\nVERDICT: {verdict}\n  {why}")

    out["honest_limits"] = [
        "Identity is held at 90 throughout -- this sweeps ONE axis. A joint identity x coverage grid was "
        "not run.",
        "The split is by ISOLATE, not by lineage or serovar; near-identical genomes could land on both "
        "sides, which would make CONFIRM optimistic.",
        "CONFIRM is ~half the cohort, so its counts are small and its margins correspondingly noisy.",
        "The candidate grid is coarse (5 values). A finer grid might find a better cut; the point here "
        "is whether the deployed value is DEFENSIBLE, not whether it is globally optimal.",
        "All cuts are derived from one permissive blastn pass per isolate, so they share any limitation "
        "of that pass (DB coverage, alignment parameters) equally.",
    ]
    a.out.write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")
    print(f"\nwrote {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
