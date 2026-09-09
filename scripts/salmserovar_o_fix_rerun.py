"""Re-score the serovar caller after the O-procedure port, against a bar frozen BEFORE the code existed.

WHAT CHANGED. `dna_decode/salmserovar/runner.py` now ports SeqSero2 1.3.2's three-branch O decision
procedure instead of picking the best-scoring O allele, and the antigen DB header carries the SeqSero2
role the old build discarded. Nothing else moved: the H1/H2 headers are byte-identical to the previous
build, and the identity/coverage thresholds are unchanged.

WHY SEQSERO2 IS NOT RE-RUN. Its binary, image and database are untouched by this change, so its column
is fixed by construction and is read from the 2026-09-08 checkpoint. That reuse is only legitimate if
the cache IS the run the bar was frozen against, so this script REFUSES unless the cached SeqSero2
tally reproduces the frozen baseline exactly. Reusing a cache from a different run would silently
compare our new numbers against someone else's denominator.

THE BAR IS FROZEN AND IS NOT RE-DERIVED HERE. `wiki/salmserovar_o_fix_acceptance_bar.json` was
registered before any of this code was written; this script reads it and reports ADOPT/REJECT against
it. Accuracy means WET-LAB accuracy on all 200 isolates -- never agreement with SeqSero2, which
becomes circular the moment its logic is ported.
"""
from __future__ import annotations

import argparse
import json
import sys
from datetime import date as _date
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from dna_decode.salmserovar.equivalence import equivalent, load_formula_index  # noqa: E402
from dna_decode.salmserovar.runner import call_serovar  # noqa: E402

BLASTN = "C:/Users/Farshad/ncbi-blast/bin/blastn.exe"
SEROVAR_DB = ROOT / "data" / "salmserovar_db"
BAR = ROOT / "wiki" / "salmserovar_o_fix_acceptance_bar.json"
NON_SPECIFIC = {"", "-", "na", "n/a", "none", "unknown", "undetermined", "pending",
                "not applicable", "not determined", "untypeable", "untypable"}


def score(call: str | None, truth: str, idx: dict) -> str:
    """hit / miss / no_call -- the SAME three-way split the 2026-09-08 run used.

    Abstention is counted separately from error on purpose: pooling them flatters whichever caller
    abstains most, and this change could plausibly trade one for the other."""
    if not call or str(call).strip().lower() in NON_SPECIFIC:
        return "no_call"
    ok, _ = equivalent(str(call), truth, idx)
    return "hit" if ok else "miss"


def tally(calls: list[str | None], truths: list[str], idx: dict) -> dict:
    c = {"hit": 0, "miss": 0, "no_call": 0}
    for call, truth in zip(calls, truths):
        c[score(call, truth, idx)] += 1
    n = len(truths)
    return {**c, "n": n, "accuracy": c["hit"] / n if n else None}


def axis_of(formula: str | None, i: int) -> str | None:
    """Axis `i` of an `O:H1:H2` formula, or None when that axis did not resolve.

    `O?` / `H?` / `-` are the caller's own unresolved markers; treating them as ordinary strings would
    let two abstentions 'agree' with each other and inflate the agreement count."""
    if not formula:
        return None
    parts = formula.split(":")
    if len(parts) <= i:
        return None
    v = parts[i].strip()
    return None if v in ("", "O?", "H?", "-") else v


def axis_agreement(ours: list[str | None], theirs: list[str | None]) -> dict:
    """Per-axis agreement on the FIXED cohort denominator.

    `ours_unresolved` is its OWN bucket rather than a disagreement, and the three buckets always sum
    to the full cohort -- which is what stops a fix from 'improving' agreement by abstaining more."""
    agree = differ = unresolved = 0
    for o, t in zip(ours, theirs):
        if o is None:
            unresolved += 1
        elif o == t:
            agree += 1
        else:
            differ += 1
    return {"agree": agree, "differ": differ, "ours_unresolved": unresolved,
            "n": len(ours)}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cohort", type=Path,
                    default=ROOT / "wiki" / "salmserovar_cohort_2026-09-04.json")
    ap.add_argument("--asm-root", type=Path, default=Path("D:/dna_decode_cache/salm_asm"))
    ap.add_argument("--db-dir", type=Path, default=SEROVAR_DB)
    ap.add_argument("--blastn", default=BLASTN)
    ap.add_argument("--ss2-checkpoint", type=Path,
                    default=Path("D:/dna_decode_cache/ss2_checkpoint/rows.jsonl"))
    ap.add_argument("--checkpoint", type=Path,
                    default=Path("D:/dna_decode_cache/ss2_checkpoint/ofix_rows.jsonl"))
    ap.add_argument("--bar", type=Path, default=BAR)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--out", type=Path, default=ROOT / "wiki" /
                    f"salmserovar_o_fix_result_{_date.today().isoformat()}.json")
    a = ap.parse_args()

    bar = json.loads(a.bar.read_text(encoding="utf-8"))
    base = bar["baseline_2026-09-08"]
    idx = load_formula_index(a.db_dir / "serovar_table.tsv")

    cached = {}
    for ln in a.ss2_checkpoint.read_text(encoding="utf-8").splitlines():
        if ln.strip():
            r = json.loads(ln)
            cached[r["asm_acc"]] = r
    rows = [r for r in cached.values() if r.get("status") == "ok"]
    if a.limit:
        rows = rows[:a.limit]

    truths = [r["truth"] for r in rows]
    ss2_calls = [(r.get("seqsero2") or {}).get("serotype") for r in rows]
    ss2_tally = tally(ss2_calls, truths, idx)

    # --- the cache-identity gate ------------------------------------------------------------------
    # Reusing the SeqSero2 column is only sound if this cache is the run the bar was frozen against.
    if not a.limit:
        want = base["seqsero2_overall"]
        got = {k: ss2_tally[k] for k in ("hit", "miss", "no_call")}
        if len(rows) != base["n_cohort"] or got != {k: want[k] for k in ("hit", "miss", "no_call")}:
            print(f"REFUSING: cached SeqSero2 column does not reproduce the frozen baseline.\n"
                  f"  n={len(rows)} (want {base['n_cohort']}), tally={got} (want "
                  f"{ {k: want[k] for k in ('hit','miss','no_call')} })", file=sys.stderr)
            return 3

    # --- re-run OUR caller only --------------------------------------------------------------------
    a.checkpoint.parent.mkdir(parents=True, exist_ok=True)
    done: dict[str, dict] = {}
    if a.checkpoint.exists():
        for ln in a.checkpoint.read_text(encoding="utf-8").splitlines():
            if ln.strip():
                r = json.loads(ln)
                done[r["asm_acc"]] = r
        print(f"resuming: {len(done)} already re-scored")
    fh = open(a.checkpoint, "a", encoding="utf-8")
    for n, r in enumerate(rows, 1):
        acc = r["asm_acc"]
        if acc in done:
            continue
        fa = next((p for p in (a.asm_root / acc).glob("*.fna")), None) \
            if (a.asm_root / acc).is_dir() else None
        rec = {"asm_acc": acc}
        if fa is None:
            rec["error"] = "assembly_missing"
        else:
            try:
                out = call_serovar(fa, a.db_dir, blastn_bin=a.blastn)
                rec.update({"serovar": out.get("serovar"), "formula": out.get("antigenic_formula"),
                            "o_antigen_rule": out.get("o_antigen_rule"), "status": out.get("status")})
            except Exception as e:                            # noqa: BLE001
                rec["error"] = type(e).__name__
        fh.write(json.dumps(rec) + "\n")
        fh.flush()
        done[acc] = rec
        if n % 25 == 0:
            print(f"  [{n}/{len(rows)}] {acc} -> {rec.get('formula')} "
                  f"({rec.get('o_antigen_rule')})", flush=True)
    fh.close()

    missing = [acc for acc in (r["asm_acc"] for r in rows) if done[acc].get("error")]
    if missing:
        # A row our caller could not score is not an abstention -- it is an absent measurement, and
        # counting it as `no_call` would credit the change with caution it did not exercise.
        print(f"REFUSING: {len(missing)} isolates could not be re-scored "
              f"(e.g. {missing[:3]}). The denominator would not match the frozen bar.",
              file=sys.stderr)
        return 3

    new_serovar = [done[r["asm_acc"]]["serovar"] for r in rows]
    new_formula = [done[r["asm_acc"]]["formula"] for r in rows]
    old_formula = [r.get("ours_formula") for r in rows]
    ours_tally = tally(new_serovar, truths, idx)
    old_tally = tally([r.get("ours") for r in rows], truths, idx)

    axes = {}
    for i, name in enumerate(("O", "H1", "H2")):
        theirs = [(r.get("seqsero2") or {}).get(name) for r in rows]
        theirs = [None if (t in (None, "", "-")) else t for t in theirs]
        axes[name] = {"after": axis_agreement([axis_of(f, i) for f in new_formula], theirs),
                      "before": axis_agreement([axis_of(f, i) for f in old_formula], theirs)}

    # --- the frozen bar, applied ------------------------------------------------------------------
    b_o = base["per_axis_vs_seqsero2"]["O"]
    checks = [
        ("O agreement rises by >= +15 over 159 on the fixed 200 denominator",
         axes["O"]["after"]["agree"] >= b_o["agree"] + 15,
         f"{b_o['agree']} -> {axes['O']['after']['agree']} "
         f"({axes['O']['after']['agree'] - b_o['agree']:+d})"),
    ]
    for name in ("H1", "H2"):
        wb = base["per_axis_vs_seqsero2"][name]
        aft = axes[name]["after"]
        checks.append((f"{name} agreement does not fall below baseline",
                       aft["agree"] >= wb.get("agree", 0),
                       f"{wb.get('agree', 0)} -> {aft['agree']}"))
        checks.append((f"{name} ours_unresolved does not rise above baseline",
                       aft["ours_unresolved"] <= wb.get("ours_unresolved", 0),
                       f"{wb.get('ours_unresolved', 0)} -> {aft['ours_unresolved']}"))
    wo = base["ours_overall"]
    checks.append(("overall wet-lab accuracy does not fall below 0.7050",
                   ours_tally["accuracy"] >= wo["accuracy"],
                   f"{wo['accuracy']:.4f} -> {ours_tally['accuracy']:.4f}"))
    checks.append(("overall no_call does not rise above 20",
                   ours_tally["no_call"] <= wo["no_call"],
                   f"{wo['no_call']} -> {ours_tally['no_call']}"))

    verdict = "ADOPT" if all(ok for _, ok, _ in checks) else "REJECT"

    print(f"\n=== {len(rows)} isolates, wet-lab labels ===")
    print(f"  ours BEFORE  acc {old_tally['accuracy']:.4f}  "
          f"(hit {old_tally['hit']} / miss {old_tally['miss']} / no_call {old_tally['no_call']})")
    print(f"  ours AFTER   acc {ours_tally['accuracy']:.4f}  "
          f"(hit {ours_tally['hit']} / miss {ours_tally['miss']} / no_call {ours_tally['no_call']})")
    print(f"  SeqSero2     acc {ss2_tally['accuracy']:.4f}  "
          f"(hit {ss2_tally['hit']} / miss {ss2_tally['miss']} / no_call {ss2_tally['no_call']})")
    print("\n  per-axis agreement with SeqSero2 (fixed denominator):")
    for name in ("O", "H1", "H2"):
        bf, af = axes[name]["before"], axes[name]["after"]
        print(f"    {name:2s} agree {bf['agree']:3d} -> {af['agree']:3d}   "
              f"differ {bf['differ']:3d} -> {af['differ']:3d}   "
              f"unresolved {bf['ours_unresolved']:3d} -> {af['ours_unresolved']:3d}")
    print("\n  frozen bar:")
    for label, ok, detail in checks:
        print(f"    [{'PASS' if ok else 'FAIL'}] {label}  ({detail})")
    print(f"\nVERDICT: {verdict}")

    rules: dict[str, int] = {}
    for r in rows:
        rule = done[r["asm_acc"]].get("o_antigen_rule") or "?"
        rules[rule] = rules.get(rule, 0) + 1

    # The PER-ISOLATE transition table. A tally cannot distinguish "20 misses became hits" from
    # "20 hits became misses and 20 misses became hits", and the bar's `no_call` ceiling is a NET
    # quantity -- so the net figure alone cannot say whether a rise came from laundering wrong calls
    # into abstentions (which the bar forbids) or from rescuing abstentions into hits (which it
    # wants). Only the transitions separate those, so they are computed here rather than argued.
    trans: dict[str, int] = {}
    for r, call in zip(rows, new_serovar):
        before = score(r.get("ours"), r["truth"], idx)
        after = score(call, r["truth"], idx)
        trans[f"{before}->{after}"] = trans.get(f"{before}->{after}", 0) + 1
    regressions = sum(v for k, v in trans.items() if k.startswith("hit->") and k != "hit->hit")
    print("\n  per-isolate transitions:")
    for k, v in sorted(trans.items(), key=lambda kv: -kv[1]):
        tag = "  <-- REGRESSION" if k.startswith("hit->") and k != "hit->hit" else ""
        print(f"    {k:20s} x{v:3d}{tag}")
    print(f"    regressions (hit -> anything else): {regressions}")

    out = {
        "schema": "salmserovar-o-fix-result-v1", "date": _date.today().isoformat(),
        "bar": str(a.bar.relative_to(ROOT)), "bar_status": bar["status"],
        "change": ("ported SeqSero2 1.3.2's three-branch O decision procedure "
                   "(call_O_and_H_type) and restored the SeqSero2 role to the antigen DB header"),
        "n": len(rows),
        "ours_before": old_tally, "ours_after": ours_tally, "seqsero2": ss2_tally,
        "delta_vs_reference_tool_before": old_tally["accuracy"] - ss2_tally["accuracy"],
        "delta_vs_reference_tool_after": ours_tally["accuracy"] - ss2_tally["accuracy"],
        "per_axis_vs_seqsero2": axes,
        "o_antigen_rule_counts": dict(sorted(rules.items(), key=lambda kv: -kv[1])),
        "transitions": dict(sorted(trans.items(), key=lambda kv: -kv[1])),
        "n_regressions": regressions,
        "accuracy_on_called_before": (old_tally["hit"] / (old_tally["hit"] + old_tally["miss"])
                                      if old_tally["hit"] + old_tally["miss"] else None),
        "accuracy_on_called_after": (ours_tally["hit"] / (ours_tally["hit"] + ours_tally["miss"])
                                     if ours_tally["hit"] + ours_tally["miss"] else None),
        "bar_checks": [{"check": c, "pass": ok, "detail": d} for c, ok, d in checks],
        "verdict": verdict,
        "verdict_is_mechanical": ("computed from the frozen bar with no judgment applied. If this "
                                  "cell is shipped despite a REJECT, that is an OVERRIDE and is "
                                  "recorded as such in the memo -- it is never written back here."),
        "honest_limits": [
            "The O procedure is now a PORT of SeqSero2's, so agreement with SeqSero2 on the O axis is "
            "a COMPATIBILITY measure, not independent corroboration. Only the WET-LAB accuracy is an "
            "independent number, which is why the bar is written on it.",
            "SeqSero2's branch gates are k-mer coverage-of-reference scores; ours are BLAST coverage. "
            "The mapping is argued on the shared axis, not measured -- see "
            "wiki/seqsero2_o_algorithm_reading_2026-09-09.md.",
            "Our O_dict membership bar is the caller's coverage threshold (40), STRICTER than "
            "SeqSero2's 25. That divergence is deliberately left in place: 40 was itself adopted "
            "against a pre-registered bar, and moving two things at once would make neither "
            "measurable.",
            "SeqSero2's column is reused from the 2026-09-08 checkpoint (its binary and DB are "
            "untouched by this change); the run refuses unless that cache reproduces the frozen "
            "baseline tally exactly.",
            "Per-serovar cap 12 flattens prevalence, so these are per-isolate accuracies on a "
            "deliberately diverse mix, NOT population-weighted rates.",
        ],
    }
    a.out.write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {a.out}")
    return 0 if verdict == "ADOPT" else 1


if __name__ == "__main__":
    raise SystemExit(main())
