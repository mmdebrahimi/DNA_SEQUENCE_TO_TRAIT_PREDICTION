"""Does extending the gentamicin rule to 16S rRNA methyltransferases REGRESS the validated cells?

WHAT THIS IS, AND WHAT IT IS NOT. The first prospective accrual (2026-08-24) located a real catalog gap:
E. coli x gentamicin scored sens 0.429, and 24 of the 28 false negatives carried a 16S rRNA
methyltransferase (`rmtE1`/`rmtE`/`armA`) while 0 carried `aac(3)`. The frozen rule is
`subclass_any={"GENTAMICIN"}` and AMRFinder files `rmt` under the generic `AMINOGLYCOSIDE` subclass, so
the rule cannot see them. That diagnosis is not re-litigated here.

This script answers the NEXT question, which is the one that gates any revision: **would adding `rmt`
break the cells that currently work?** A fix that rescues the prospective misses but wrecks the validated
provenance-disjoint cells is not a fix. That check runs entirely on cached AMRFinder output — no Docker,
no network.

IT DOES NOT MEASURE THE RESCUE. The prospective cohort's per-isolate determinant calls are not on disk
(the committed artifact carries only a confusion matrix), so the "does it fix the 28 FN" half cannot be
computed here and is NOT claimed. In-distribution regression is what is measurable today.

THE FROZEN SURFACE IS NOT TOUCHED. `dna_decode/eval/amr_rules.py` is read-only to this script; the
candidate rule is applied scorer-locally, mirroring the `experimental_drug_rules.py` overlay pattern.

VERIFIED, NOT ASSUMED: the frozen rule matches Subclass by TOKEN, not exact equality, so compound
subclasses that name gentamicin (`GENTAMICIN/KANAMYCIN/TOBRAMYCIN`, `APRAMYCIN/GENTAMICIN/TOBRAMYCIN`)
are ALREADY counted. I expected them to be a second blind spot; they are not. Only the generic
`AMINOGLYCOSIDE` filing is invisible.

Run: uv run python scripts/gentamicin_rmt_candidate.py
"""
from __future__ import annotations

import csv
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

# 16S rRNA methyltransferases confer high-level 4,6-deoxystreptamine resistance INCLUDING gentamicin.
# Explicit family regexes rather than a substring: `rmt` alone would also catch unrelated symbols, and
# this project has been bitten by loose prefix matching before.
METHYLTRANSFERASE = re.compile(r"^(rmt[A-H]\d*|armA|npmA\d*)$", re.I)

COHORTS = ("escherichia_coli_shigella_provdisjoint_gentamicin",
           "klebsiella_provdisjoint_gentamicin",
           "klebsiella_gentamicin")


def gene_symbol(row: dict) -> str:
    return (row.get("Element symbol") or row.get("Gene symbol") or "").strip()


def is_methyltransferase(symbol: str) -> bool:
    """PURE. True for a 16S rRNA methyltransferase gene symbol."""
    return bool(METHYLTRANSFERASE.match(symbol.strip()))


_INDEX: dict[str, Path] | None = None


def amrfinder_index() -> dict[str, Path]:
    """accession -> main.tsv, built ONCE across the whole data tree. Cached.

    AMRFinder output is per-GENOME, not per-drug, so cohorts SHARE runs: klebsiella_gentamicin's
    accessions have their output under klebsiella_cipro/amrfinder_runs/. Looking only in a cohort's own
    directory silently loses most of the cohort -- the first version of this script scored 63 of ~148
    isolates for exactly that reason and reported the rest as "missing AMRFinder", which reads as a data
    gap rather than a lookup bug.
    """
    global _INDEX
    if _INDEX is None:
        _INDEX = {}
        for m in (ROOT / "data").rglob("amrfinder_runs/*/main.tsv"):
            _INDEX.setdefault(m.parent.name, m)
    return _INDEX


def read_rows(main_tsv: Path) -> list[dict]:
    if not main_tsv.exists():
        return []
    with open(main_tsv, encoding="utf-8", errors="replace") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def frozen_call(rows: list[dict]) -> bool:
    """Reproduce the FROZEN rule locally: >=1 row whose Class or Subclass carries a GENTAMICIN token.

    Deliberately re-implemented rather than imported so this script cannot mutate the frozen module by
    accident; `main()` asserts it agrees with the real `call_resistance` before any number is reported.
    """
    for r in rows:
        blob = ((r.get("Class") or "") + "|" + (r.get("Subclass") or "")).upper()
        if "GENTAMICIN" in blob:
            return True
    return False


def candidate_call(rows: list[dict]) -> bool:
    """FROZEN rule OR a 16S methyltransferase. Scorer-local; the frozen module is untouched."""
    return frozen_call(rows) or any(is_methyltransferase(gene_symbol(r)) for r in rows)


def score(cohort: str) -> dict | None:
    base = ROOT / "data" / "raw" / cohort
    sel = base / "selected.tsv"
    if not sel.exists():
        return None
    out = {"cohort": cohort, "n": 0, "frozen": {"tp": 0, "fp": 0, "tn": 0, "fn": 0},
           "candidate": {"tp": 0, "fp": 0, "tn": 0, "fn": 0}, "no_amrfinder": 0,
           "carries_mtase": 0, "mtase_genes": {}, "changed": []}
    for line in sel.read_text(encoding="utf-8").splitlines():
        parts = line.split("\t")
        if len(parts) < 2:
            continue
        acc, label = parts[0].strip(), parts[1].strip().upper()
        if label not in ("R", "S"):
            continue
        main = amrfinder_index().get(acc, base / "amrfinder_runs" / acc / "main.tsv")
        rows = read_rows(main)
        if not rows and not main.exists():
            out["no_amrfinder"] += 1
            continue
        out["n"] += 1
        mt = [gene_symbol(r) for r in rows if is_methyltransferase(gene_symbol(r))]
        if mt:
            out["carries_mtase"] += 1
            for g in mt:
                out["mtase_genes"][g] = out["mtase_genes"].get(g, 0) + 1
        for name, call in (("frozen", frozen_call(rows)), ("candidate", candidate_call(rows))):
            key = ("tp" if call else "fn") if label == "R" else ("fp" if call else "tn")
            out[name][key] += 1
        if frozen_call(rows) != candidate_call(rows):
            out["changed"].append({"accession": acc, "label": label, "mtase": mt})
    return out


def _s_labelled_carriers(results: list[dict]) -> int:
    """S-labelled methyltransferase carriers -- the ONLY isolates that can become a new false positive.

    Zero of them makes the specificity comparison vacuous, which is a materially different statement from
    "the candidate is safe". Counted, not assumed.
    """
    n = 0
    for r in results:
        base = ROOT / "data" / "raw" / r["cohort"] / "selected.tsv"
        if not base.exists():
            continue
        for line in base.read_text(encoding="utf-8").splitlines():
            parts = line.split("	")
            if len(parts) < 2 or parts[1].strip().upper() != "S":
                continue
            main = amrfinder_index().get(parts[0].strip())
            if main and any(is_methyltransferase(gene_symbol(x)) for x in read_rows(main)):
                n += 1
    return n


def metrics(c: dict) -> dict:
    tp, fp, tn, fn = c["tp"], c["fp"], c["tn"], c["fn"]
    n = tp + fp + tn + fn
    return {"acc": round((tp + tn) / n, 3) if n else None,
            "sens": round(tp / (tp + fn), 3) if (tp + fn) else None,
            "spec": round(tn / (tn + fp), 3) if (tn + fp) else None}


def main() -> int:
    from dna_decode.eval.amr_rules import call_resistance

    results = [r for r in (score(c) for c in COHORTS) if r]
    if not results:
        print("no cohorts on disk")
        return 1

    # CONTROL: the local re-implementation must agree with the real frozen rule, or nothing below means
    # anything. Checked on every genome that has cached AMRFinder output.
    checked = mismatch = 0
    for cohort in COHORTS:
        base = ROOT / "data" / "raw" / cohort
        if not (base / "selected.tsv").exists():
            continue
        for line in (base / "selected.tsv").read_text(encoding="utf-8").splitlines():
            acc = line.split("\t")[0].strip()
            main = amrfinder_index().get(acc)
            if main is None or not main.exists():
                continue
            real = call_resistance(main, "gentamicin")["prediction"] == "R"
            if real != frozen_call(read_rows(main)):
                mismatch += 1
            checked += 1
    print(f"CONTROL: local frozen re-implementation vs call_resistance -> "
          f"{checked - mismatch}/{checked} agree")
    if mismatch:
        print("  REFUSING to report: the re-implementation does not match the frozen rule.")
        return 2

    print()
    total = {"frozen": {"tp": 0, "fp": 0, "tn": 0, "fn": 0}, "candidate": {"tp": 0, "fp": 0, "tn": 0, "fn": 0}}
    for r in results:
        print(f"{r['cohort']}  (n={r['n']}, missing AMRFinder {r['no_amrfinder']})")
        print(f"  carries a 16S methyltransferase: {r['carries_mtase']}  {r['mtase_genes']}")
        for name in ("frozen", "candidate"):
            m = metrics(r[name])
            c = r[name]
            print(f"  {name:10} acc {m['acc']} sens {m['sens']} spec {m['spec']}   "
                  f"(tp{c['tp']} fp{c['fp']} tn{c['tn']} fn{c['fn']})")
            for k in total[name]:
                total[name][k] += c[k]
        if r["changed"]:
            print(f"  calls CHANGED by the candidate: {len(r['changed'])}")
            for ch in r["changed"][:6]:
                print(f"    {ch['accession']} label={ch['label']} mtase={ch['mtase']}")
    print()
    print("POOLED")
    for name in ("frozen", "candidate"):
        m = metrics(total[name])
        c = total[name]
        print(f"  {name:10} acc {m['acc']} sens {m['sens']} spec {m['spec']}   "
              f"(tp{c['tp']} fp{c['fp']} tn{c['tn']} fn{c['fn']})")
    regressed = total["candidate"]["fp"] > total["frozen"]["fp"]
    # How much does "specificity unchanged" actually say? Only an S-labelled CARRIER can turn into a new
    # false positive. If there are none, the specificity check is VACUOUS -- true, and carrying no
    # information about over-calling. Report that rather than let an unchanged number imply a test.
    s_carriers = sum(1 for r in results for ch in [] ) or _s_labelled_carriers(results)
    print()
    print(f"false positives {total['frozen']['fp']} -> {total['candidate']['fp']}; "
          f"S-labelled methyltransferase carriers in the data: {s_carriers}")
    if s_carriers == 0:
        print("VERDICT: the specificity check is VACUOUS. Every methyltransferase carrier in these")
        print("  cohorts is R-labelled, so the candidate CANNOT produce a false positive here. Its")
        print("  over-calling risk is UNTESTED, not zero -- do not read 'spec unchanged' as evidence.")
    else:
        print(f"VERDICT: candidate {'REGRESSES' if regressed else 'does not regress'} specificity over "
              f"{s_carriers} S-labelled carrier(s).")
    print("ALSO NOT MEASURED: whether it rescues the prospective false negatives -- those isolates'")
    print("determinant calls are not on disk, so that half is unproven and is not claimed.")

    (ROOT / "wiki" / "gentamicin_rmt_candidate.json").write_text(
        json.dumps({"cohorts": results, "pooled": total,
                    "regresses_specificity": regressed,
                    "s_labelled_carriers": s_carriers,
                    "specificity_check_vacuous": s_carriers == 0}, indent=2), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
