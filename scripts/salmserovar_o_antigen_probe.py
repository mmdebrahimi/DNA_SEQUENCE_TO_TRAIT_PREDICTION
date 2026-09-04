"""Is the O-antigen failure a THRESHOLD problem or a DATABASE-COVERAGE problem? They cost differently.

WHY THIS EXISTS. The abstention anatomy found the O antigen to be the largest cause of no-calls
(21 of 59, 35.6%) and I immediately wrote that its fix is "DB coverage -- data engineering, not a code
change". That was ASSERTED, not measured, one commit after correcting a different assertion made the
same way. So it gets the same treatment.

The two hypotheses cost wildly different amounts and are cleanly separable:

  THRESHOLD  -- the true O allele IS in the DB and blastn DOES hit it, but below the identity/coverage
                cut, so it is discarded. Fix: tune a number. Cheap.
  COVERAGE   -- no allele for that O group exists in the DB at all, so no threshold admits it.
                Fix: extend the wzx/wzy allele set. Data engineering.

The probe re-runs the SAME genomes at a deliberately permissive threshold and asks what appears. A
sub-threshold hit that would have named the right O group is evidence for THRESHOLD; nothing at any
threshold is evidence for COVERAGE.

IT CANNOT BE READ AS A RECOMMENDATION TO LOWER THE THRESHOLD. Finding recoverable hits at 60% identity
says the information is present, not that admitting everything at 60% is safe -- that would trade
abstentions for wrong calls, and abstention is the safer failure. The output is a DIAGNOSIS of where
the information is, and the wrong-call cost of actually moving the cut is explicitly not measured here.

Offline: cached assemblies + the committed antigen DB + native blastn.
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

from dna_decode.typing.blast_caller import call_alleles  # noqa: E402
from dna_decode.salmserovar.runner import parse_axis_antigen  # noqa: E402

BLASTN = "C:/Users/Farshad/ncbi-blast/bin/blastn.exe"
DEPLOYED_ID, DEPLOYED_COV = 90.0, 80.0


def o_groups_in_db(fasta: Path) -> set[str]:
    out = set()
    for line in fasta.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith(">"):
            pa = parse_axis_antigen(line[1:].strip())
            if pa and pa[0] == "O":
                out.add(pa[1])
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--results", type=Path,
                    default=Path("D:/dna_decode_cache/salm_asm/results.jsonl"))
    ap.add_argument("--asm-root", type=Path, default=Path("D:/dna_decode_cache/salm_asm"))
    ap.add_argument("--db-dir", type=Path, default=ROOT / "data" / "salmserovar_db")
    ap.add_argument("--probe-identity", type=float, default=60.0)
    ap.add_argument("--probe-coverage", type=float, default=30.0)
    ap.add_argument("--limit", type=int, default=21)
    ap.add_argument("--blastn", default=BLASTN)
    ap.add_argument("--out", type=Path,
                    default=ROOT / "wiki" / f"salmserovar_o_antigen_probe_{_date.today().isoformat()}.json")
    a = ap.parse_args()

    rows = [json.loads(ln) for ln in a.results.read_text(encoding="utf-8").splitlines() if ln.strip()]
    # the O-unresolved abstentions, exactly as the anatomy classified them
    targets = [r for r in rows
               if r.get("status") == "ok" and not r.get("pred_serovar")
               and (r.get("pred_formula") or "").split(":")[0].strip() in ("O?", "", "-")][:a.limit]
    if not targets:
        print("no O-unresolved abstentions to probe", file=sys.stderr)
        return 3

    db = a.db_dir / "salmonella_antigens.fasta"
    in_db = o_groups_in_db(db)
    print(f"O groups present in the DB: {len(in_db)}")
    print(f"probing {len(targets)} O-unresolved isolates at identity>={a.probe_identity} "
          f"coverage>={a.probe_coverage} (deployed: {DEPLOYED_ID}/{DEPLOYED_COV})\n")

    out_rows, verdicts = [], collections.Counter()
    for n, r in enumerate(targets, 1):
        acc = r["asm_acc"]
        fa = next((p for p in (a.asm_root / acc).glob("*.fna")), None)
        if fa is None:
            verdicts["assembly_missing"] += 1
            continue
        try:
            res = call_alleles(fa, db, identity_threshold=a.probe_identity,
                               coverage_threshold=a.probe_coverage, blastn_bin=a.blastn, timeout=600)
        except Exception as e:                        # noqa: BLE001
            verdicts[f"error:{type(e).__name__}"] += 1
            continue
        if res.get("status") != "ok":
            verdicts["blast_unavailable"] += 1
            continue

        o_hits = []
        for allele_id, hit in res["per_allele"].items():
            pa = parse_axis_antigen(allele_id)
            if not pa or pa[0] != "O":
                continue
            if hit["percent_identity"] >= a.probe_identity and \
               hit["percent_coverage"] >= a.probe_coverage:
                o_hits.append({"antigen": pa[1], "pid": hit["percent_identity"],
                               "cov": hit["percent_coverage"]})
        o_hits.sort(key=lambda h: (-h["pid"], -h["cov"]))
        best = o_hits[0] if o_hits else None
        # would the DEPLOYED thresholds have admitted the best permissive hit?
        sub = bool(best and (best["pid"] < DEPLOYED_ID or best["cov"] < DEPLOYED_COV))
        v = ("no_O_hit_at_any_threshold" if not o_hits
             else ("sub_threshold_hit_exists" if sub else "hit_above_deployed_threshold"))
        verdicts[v] += 1
        out_rows.append({"asm_acc": acc, "label": r.get("label"), "verdict": v,
                         "n_O_hits_permissive": len(o_hits), "best": best})
        if n % 5 == 0:
            print(f"  [{n}/{len(targets)}] {acc} {v}", flush=True)

    n = sum(verdicts[k] for k in ("no_O_hit_at_any_threshold", "sub_threshold_hit_exists",
                                  "hit_above_deployed_threshold"))
    print(f"\n=== {n} isolates probed ===")
    for k, v in verdicts.most_common():
        print(f"  {v:>3}  {k}")

    if n == 0:
        print("\nREFUSING: nothing was probed, so neither hypothesis was tested.", file=sys.stderr)
        return 3

    sub_n = verdicts["sub_threshold_hit_exists"]
    none_n = verdicts["no_O_hit_at_any_threshold"]
    if none_n > sub_n:
        verdict = "COVERAGE_IS_THE_BINDING_CONSTRAINT"
        why = (f"{none_n} of {n} O-unresolved isolates produce NO O-antigen hit even at "
               f"identity>={a.probe_identity}/coverage>={a.probe_coverage}. No threshold admits an "
               "allele that is not in the DB, so the fix is extending the wzx/wzy allele set -- data "
               "engineering, as claimed.")
    else:
        verdict = "THRESHOLD_IS_THE_BINDING_CONSTRAINT"
        why = (f"{sub_n} of {n} O-unresolved isolates DO hit an O allele below the deployed "
               f"{DEPLOYED_ID}/{DEPLOYED_COV} cut, so the information is already in the DB and the "
               "earlier 'DB coverage / data engineering' claim was WRONG. NOTE: this does NOT license "
               "lowering the threshold -- that trades abstentions for wrong calls, and the wrong-call "
               "cost is NOT measured here.")
    print(f"\nVERDICT: {verdict}\n  {why}")

    art = {"schema": "salmserovar-o-antigen-probe-v1", "date": _date.today().isoformat(),
           "question": "is the O-antigen abstention a threshold problem or a DB-coverage problem?",
           "tests_the_claim": ("the 2026-09-04 anatomy asserted 'the real priority is O-antigen DB "
                               "coverage -- data engineering, not a code change' WITHOUT measuring it"),
           "deployed_thresholds": {"identity": DEPLOYED_ID, "coverage": DEPLOYED_COV},
           "probe_thresholds": {"identity": a.probe_identity, "coverage": a.probe_coverage},
           "n_O_groups_in_db": len(in_db), "n_probed": n,
           "verdicts": dict(verdicts), "rows": out_rows,
           "verdict": verdict, "why": why,
           "honest_limits": [
               "A permissive hit means the INFORMATION is present, NOT that admitting it is safe. "
               "Lowering the cut trades abstentions for wrong calls and that cost is NOT measured "
               "here; abstention is the safer failure.",
               "Only the O-unresolved abstentions are probed -- this says nothing about the H axes or "
               "about isolates that already resolve.",
               "One antigen-DB build (62 O alleles) and one cohort. A different DB changes the answer.",
           ]}
    a.out.write_text(json.dumps(art, indent=2) + "\n", encoding="utf-8")
    print(f"\nwrote {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
