"""Held-out confirmation of the identity-primary serotype fix, on isolates it was never chosen from.

WHY THIS RUN EXISTS. The fix was found by inspecting the E. coli discovery cohort's failure pattern and
then applied: H accuracy 0.770 -> 0.926 on those same 150 isolates. That is a strong result and an
UNBLINDED one -- the change was selected after seeing what it would repair, so its effect size is
exactly the kind that shrinks on replication. The discovery memo says so; this run is the check.

DISJOINT BY CONSTRUCTION, two ways: a different cohort seed, AND an explicit exclusion of every
accession already scored in the discovery checkpoint. Overlap is measured and reported, never assumed.

BOTH RULES ARE RUN ON THE SAME GENOMES IN THE SAME PASS, from one blastn result per isolate. Scoring
two orderings over identical hits removes any run-to-run variation from the comparison -- the only
difference between the two numbers is the sort key.

The prediction was written down BEFORE this ran (H gain > +0.05; resolution unchanged), so the run can
fail. It is stamped into the artifact so a reader can check the claim was not adjusted afterwards.
"""
from __future__ import annotations

import argparse
import collections
import json
import sys
import traceback
from datetime import date as _date
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from dna_decode.data.refseq import download_genome, fasta_path  # noqa: E402
from dna_decode.typing.blast_caller import call_alleles  # noqa: E402
from dna_decode.serotype.runner import antigen_of, gene_of  # noqa: E402
from serotype_oh_validate import build_cohort, norm_call  # noqa: E402

BLASTN = "C:/Users/Farshad/ncbi-blast/bin/blastn.exe"

PREREGISTERED = {
    "primary": "H accuracy under identity-primary exceeds coverage-only on held-out isolates",
    "quantitative_bar": 0.05,
    "secondary": "resolution (no_call) unchanged between rules",
    "falsified_if": "H gain <= 0, or resolution changes materially",
    "registered_before_run": True,
}


def call_both_rules(per_allele: dict) -> dict:
    """Score BOTH orderings from ONE blastn result, so only the sort key differs."""
    out = {}
    for name, identity_primary in (("identity_primary", True), ("coverage_primary", False)):
        ag_best: dict[str, dict] = {}
        for allele_id, hit in per_allele.items():
            if not hit["called"]:
                continue
            ag = antigen_of(allele_id)
            if ag is None:
                continue
            key = ((hit["percent_identity"], hit["percent_coverage"]) if identity_primary
                   else (hit["percent_coverage"], hit["percent_identity"]))
            cur = ag_best.get(ag)
            cur_key = ((cur["percent_identity"], cur["percent_coverage"]) if (cur and identity_primary)
                       else ((cur["percent_coverage"], cur["percent_identity"]) if cur else None))
            if cur is None or key > cur_key:
                ag_best[ag] = {"antigen": ag, "gene": gene_of(allele_id),
                               "percent_identity": hit["percent_identity"],
                               "percent_coverage": hit["percent_coverage"]}

        def top(prefix):
            c = [v for k, v in ag_best.items() if k.startswith(prefix)]
            if not c:
                return None
            return max(c, key=(lambda v: (v["percent_identity"], v["percent_coverage"]))
                       if identity_primary else
                       (lambda v: (v["percent_coverage"], v["percent_identity"])))["antigen"]
        out[name] = {"O": top("O"), "H": top("H")}
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--target", type=int, default=150)
    ap.add_argument("--seed", type=int, default=77, help="different from the discovery seed (23)")
    ap.add_argument("--max-rows", type=int, default=250000)
    ap.add_argument("--db", type=Path, default=ROOT / "data" / "serotypefinder_db" / "serotypefinder.fsa")
    ap.add_argument("--asm-dir", type=Path, default=Path("D:/dna_decode_cache/ecoli_sero_heldout"))
    ap.add_argument("--discovery-checkpoint", type=Path,
                    default=Path("D:/dna_decode_cache/ecoli_sero_asm/results.jsonl"))
    ap.add_argument("--checkpoint", type=Path,
                    default=Path("D:/dna_decode_cache/ecoli_sero_heldout/results.jsonl"))
    ap.add_argument("--blastn", default=BLASTN)
    ap.add_argument("--out", type=Path,
                    default=ROOT / "wiki" / f"serotype_heldout_confirm_{_date.today().isoformat()}.json")
    a = ap.parse_args()

    seen: set[str] = set()
    if a.discovery_checkpoint.exists():
        for line in a.discovery_checkpoint.read_text(encoding="utf-8").splitlines():
            if line.strip():
                seen.add(json.loads(line)["asm_acc"])
    print(f"discovery cohort accessions to EXCLUDE: {len(seen)}")

    cohort, meta = build_cohort(a.max_rows, a.target + len(seen), a.seed)
    heldout = [c for c in cohort if c["asm_acc"] not in seen][:a.target]
    overlap = sum(1 for c in cohort if c["asm_acc"] in seen)
    print(f"cohort {len(cohort)} -> held-out {len(heldout)} (overlap with discovery removed: {overlap})")

    a.asm_dir.mkdir(parents=True, exist_ok=True)
    done: dict[str, dict] = {}
    if a.checkpoint.exists():
        for line in a.checkpoint.read_text(encoding="utf-8").splitlines():
            if line.strip():
                r = json.loads(line)
                done[r["asm_acc"]] = r

    fh = open(a.checkpoint, "a", encoding="utf-8")
    for n, iso in enumerate(heldout, 1):
        acc = iso["asm_acc"]
        if acc in done:
            continue
        rec = {"asm_acc": acc, "O_label": iso["O"], "H_label": iso["H"]}
        try:
            download_genome(acc, a.asm_dir)
            fa = fasta_path(acc, a.asm_dir)
            if not Path(fa).exists():
                rec["status"] = "assembly_missing"
            else:
                res = call_alleles(fa, a.db, identity_threshold=85.0, coverage_threshold=60.0,
                                   blastn_bin=a.blastn, timeout=600)
                if res.get("status") != "ok":
                    rec["status"] = "blast_unavailable"
                    rec["reason"] = str(res.get("reason"))[:150]
                else:
                    rec["status"] = "ok"
                    rec["calls"] = call_both_rules(res["per_allele"])
        except Exception as e:                        # noqa: BLE001
            rec["status"] = f"error:{type(e).__name__}"
            rec["error"] = str(e)[:200]
            rec["trace"] = traceback.format_exc()[-300:]
        fh.write(json.dumps(rec) + "\n"); fh.flush()
        done[acc] = rec
        if n % 15 == 0 or rec["status"] != "ok":
            print(f"  [{n}/{len(heldout)}] {acc} {rec['status']}", flush=True)
    fh.close()

    rows = [done[c["asm_acc"]] for c in heldout if c["asm_acc"] in done]
    ok = [r for r in rows if r.get("status") == "ok"]
    stats = {rule: {ax: collections.Counter() for ax in ("O", "H")}
             for rule in ("identity_primary", "coverage_primary")}
    for r in ok:
        for rule in stats:
            for ax in ("O", "H"):
                lab = r.get(f"{ax}_label")
                if not lab:
                    continue
                call = norm_call(r["calls"][rule][ax], ax)
                if call is None:
                    stats[rule][ax]["no_call"] += 1
                elif call == lab:
                    stats[rule][ax]["hit"] += 1
                else:
                    stats[rule][ax]["miss"] += 1

    def acc_of(c):
        s = c["hit"] + c["miss"]
        return (c["hit"] / s) if s else None

    def nocall_of(c):
        t = sum(c.values())
        return (c["no_call"] / t) if t else None

    print(f"\n=== held-out: {len(ok)} isolates called ===")
    for rule in ("coverage_primary", "identity_primary"):
        for ax in ("O", "H"):
            c = stats[rule][ax]
            print(f"  {rule:<18} {ax}: hit={c['hit']:<4} miss={c['miss']:<4} no_call={c['no_call']:<4} "
                  f"acc={acc_of(c) if acc_of(c) is None else round(acc_of(c),4)}")

    h_gain = (acc_of(stats["identity_primary"]["H"]) or 0) - (acc_of(stats["coverage_primary"]["H"]) or 0)
    o_gain = (acc_of(stats["identity_primary"]["O"]) or 0) - (acc_of(stats["coverage_primary"]["O"]) or 0)
    res_delta = abs((nocall_of(stats["identity_primary"]["H"]) or 0)
                    - (nocall_of(stats["coverage_primary"]["H"]) or 0))
    print(f"\n  H gain {h_gain:+.4f}   O gain {o_gain:+.4f}   |no_call delta| {res_delta:.4f}")

    bar = PREREGISTERED["quantitative_bar"]
    if h_gain <= 0:
        verdict, why = "FALSIFIED", (f"identity-primary does NOT beat coverage-only on held-out "
                                     f"isolates (H gain {h_gain:+.4f}); the discovery result does not "
                                     "replicate and the fix should be reconsidered.")
    elif h_gain > bar:
        verdict, why = "CONFIRMED", (f"identity-primary beats coverage-only by {h_gain:+.4f} on "
                                     f"isolates the fix was never chosen from, clearing the "
                                     f"pre-registered +{bar} bar.")
    else:
        verdict, why = "DIRECTIONALLY_CONFIRMED_BUT_SMALLER", (
            f"the gain is positive ({h_gain:+.4f}) but below the pre-registered +{bar} bar, so the "
            "discovery estimate was inflated by the unblinded choice. The fix helps; quote the "
            "held-out number, not the discovery one.")
    print(f"\nVERDICT: {verdict}\n  {why}")

    out = {"schema": "serotype-heldout-confirm-v1", "date": _date.today().isoformat(),
           "preregistered": PREREGISTERED,
           "discovery": {"n": 150, "seed": 23, "H_before": 0.7703, "H_after": 0.9257,
                         "note": "UNBLINDED -- the fix was chosen after seeing this cohort's failures"},
           "heldout": {"n_requested": len(heldout), "n_called": len(ok), "seed": a.seed,
                       "n_excluded_as_discovery_overlap": overlap,
                       "statuses": dict(collections.Counter(r.get("status", "?") for r in rows)),
                       "cohort_meta": meta},
           "coverage_primary": {ax: dict(stats["coverage_primary"][ax]) for ax in ("O", "H")},
           "identity_primary": {ax: dict(stats["identity_primary"][ax]) for ax in ("O", "H")},
           "H_accuracy_coverage_primary": acc_of(stats["coverage_primary"]["H"]),
           "H_accuracy_identity_primary": acc_of(stats["identity_primary"]["H"]),
           "H_gain": h_gain, "O_gain": o_gain, "no_call_delta": res_delta,
           "verdict": verdict, "why": why,
           "honest_limits": [
               "Both rules are scored from ONE blastn pass per isolate, so the only difference is the "
               "sort key -- but that also means shared blast parameters/DB limits affect both equally.",
               "Disjointness is by accession against the discovery checkpoint; two accessions could "
               "still be near-identical genomes, so this is held-out by ISOLATE, not by lineage.",
               "Still no in-silico incumbent for E. coli (PD leaves computed_types null), so the "
               "ABSOLUTE accuracy remains weakly anchored; the RULE COMPARISON is what this run tests.",
               "O-only labels score the O axis alone -- the axes have different denominators.",
           ]}
    a.out.write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")
    print(f"\nwrote {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
