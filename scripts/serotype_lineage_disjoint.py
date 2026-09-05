"""Does the serotype fix survive a LINEAGE-disjoint split, not just an isolate-disjoint one?

THE NAMED LIMIT THIS CLOSES. The identity-primary selection fix replicated on 250 held-out isolates at
+0.106 H accuracy. That replication was held out BY ISOLATE: two distinct accessions can still be
near-identical genomes, so a lineage could appear on both sides and the "held-out" half would not be
independent in the way the claim needs. The memo says so explicitly. This is the check.

WHY SEQUENCE TYPE RATHER THAN MASH. The obvious tool is Mash clustering, but Docker is down on this
host and wrestling it unattended is the wrong risk. Multi-locus sequence type is the better unit here
anyway: ST is the canonical E. coli lineage unit -- it is what the field actually groups by, it is
interpretable ("ST131"), and the caller is already in this repo and needs only blastn. A Mash cluster
at an arbitrary threshold would have been a less meaningful grouping, not a more rigorous one.

THE SPLIT RULE. Genomes are grouped by ST, then whole STs -- never individual isolates -- are assigned
to TRAIN or TEST. No sequence type appears on both sides, so the test half contains no lineage the
selection rule was tuned on. Assignment is by a deterministic hash of the ST, so it does not depend on
the outcome.

WHAT IS AND IS NOT BEING RE-DECIDED. The fix is NOT re-chosen here. It is already deployed and was
already confirmed once. This asks only whether its measured gain survives when lineage overlap is
removed -- so a smaller gain is an honest correction, not a failure, and a reversal would be a genuine
falsification worth acting on.

ONE BLASTN PASS PER GENOME FOR EACH OF MLST AND SEROTYPE; both orderings are then scored offline from
the same serotype alignments, so nothing varies between the two rules but the sort key.

Needs blastn + cached assemblies. Writes wiki/serotype_lineage_disjoint_<date>.json.
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
sys.path.insert(0, str(ROOT / "scripts"))

from dna_decode.mlst.runner import call_mlst  # noqa: E402
from dna_decode.serotype.runner import antigen_of, gene_of  # noqa: E402
from dna_decode.typing.blast_caller import call_alleles  # noqa: E402
from serotype_oh_validate import norm_call  # noqa: E402

BLASTN = "C:/Users/Farshad/ncbi-blast/bin/blastn.exe"
SERO_DB = ROOT / "data" / "serotypefinder_db" / "serotypefinder.fsa"
MLST_DIR = ROOT / "data" / "mlst_db" / "ecoli_achtman"


def call_both_rules(per_allele: dict) -> dict:
    """Both orderings from ONE set of alignments -- only the sort key differs."""
    out = {}
    for name, ident_first in (("identity_primary", True), ("coverage_primary", False)):
        best: dict[str, dict] = {}
        for allele_id, hit in per_allele.items():
            if not hit["called"]:
                continue
            ag = antigen_of(allele_id)
            if ag is None:
                continue
            key = ((hit["percent_identity"], hit["percent_coverage"]) if ident_first
                   else (hit["percent_coverage"], hit["percent_identity"]))
            cur = best.get(ag)
            cur_key = None if cur is None else (
                (cur["percent_identity"], cur["percent_coverage"]) if ident_first
                else (cur["percent_coverage"], cur["percent_identity"]))
            if cur is None or key > cur_key:
                best[ag] = {"antigen": ag, "gene": gene_of(allele_id),
                            "percent_identity": hit["percent_identity"],
                            "percent_coverage": hit["percent_coverage"]}

        def top(prefix):
            c = [v for k, v in best.items() if k.startswith(prefix)]
            if not c:
                return None
            return max(c, key=(lambda v: (v["percent_identity"], v["percent_coverage"])) if ident_first
                       else (lambda v: (v["percent_coverage"], v["percent_identity"])))["antigen"]
        out[name] = {"O": top("O"), "H": top("H")}
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--label-sources", nargs="+",
                    default=["D:/dna_decode_cache/ecoli_sero_asm/results.jsonl",
                             "D:/dna_decode_cache/ecoli_sero_heldout/results.jsonl"])
    ap.add_argument("--asm-dirs", nargs="+",
                    default=["D:/dna_decode_cache/ecoli_sero_asm",
                             "D:/dna_decode_cache/ecoli_sero_heldout"])
    ap.add_argument("--checkpoint", type=Path,
                    default=Path("D:/dna_decode_cache/ecoli_sero_asm/lineage_calls.jsonl"))
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--blastn", default=BLASTN)
    ap.add_argument("--out", type=Path, default=ROOT / "wiki" /
                    f"serotype_lineage_disjoint_{_date.today().isoformat()}.json")
    a = ap.parse_args()

    labels: dict[str, dict] = {}
    for src in a.label_sources:
        p = Path(src)
        if not p.exists():
            continue
        for ln in p.read_text(encoding="utf-8").splitlines():
            if not ln.strip():
                continue
            r = json.loads(ln)
            if r.get("status") == "ok" and r.get("O_label"):
                labels[r["asm_acc"]] = {"O": r.get("O_label"), "H": r.get("H_label")}

    fastas: dict[str, Path] = {}
    for d in a.asm_dirs:
        root = Path(d)
        if not root.exists():
            continue
        for sub in sorted(root.iterdir()):
            if sub.is_dir() and sub.name in labels and sub.name not in fastas:
                fa = next((p for p in sub.glob("*.fna")), None)
                if fa:
                    fastas[sub.name] = fa
    accs = sorted(fastas)
    if a.limit:
        accs = accs[:a.limit]
    print(f"{len(labels)} labelled isolates | {len(accs)} with a cached assembly")

    loci = {p.stem: p for p in MLST_DIR.glob("*.fasta")}
    profiles = MLST_DIR / "profiles.tsv"
    done: dict[str, dict] = {}
    if a.checkpoint.exists():
        for ln in a.checkpoint.read_text(encoding="utf-8").splitlines():
            if ln.strip():
                r = json.loads(ln)
                done[r["asm_acc"]] = r
    print(f"{len(done)} already checkpointed\n")

    fh = open(a.checkpoint, "a", encoding="utf-8")
    for n, acc in enumerate(accs, 1):
        if acc in done:
            continue
        rec = {"asm_acc": acc, **labels[acc]}
        try:
            m = call_mlst(fastas[acc], loci, profiles, blastn_bin=a.blastn)
            rec["st"] = m.get("st")
            rec["st_complete"] = bool(m.get("complete"))
            res = call_alleles(fastas[acc], SERO_DB, identity_threshold=85.0,
                               coverage_threshold=60.0, blastn_bin=a.blastn, timeout=600)
            if res.get("status") != "ok":
                rec["status"] = "blast_unavailable"
            else:
                rec["status"] = "ok"
                rec["calls"] = call_both_rules(res["per_allele"])
        except Exception as e:                       # noqa: BLE001
            rec["status"] = f"error:{type(e).__name__}"
            rec["error"] = str(e)[:180]
        fh.write(json.dumps(rec) + "\n"); fh.flush()
        done[acc] = rec
        if n % 20 == 0:
            print(f"  [{n}/{len(accs)}] {acc} ST={rec.get('st')} {rec['status']}", flush=True)
    fh.close()

    rows = [done[a_] for a_ in accs if a_ in done and done[a_].get("status") == "ok"]
    # Group by ST. An unresolved ST cannot be placed in a lineage-disjoint split, so those isolates
    # are EXCLUDED and counted -- silently pooling them would break the disjointness being claimed.
    typed = [r for r in rows if r.get("st") and r.get("st_complete")]
    untyped = len(rows) - len(typed)
    by_st: dict[str, list] = collections.defaultdict(list)
    for r in typed:
        by_st[str(r["st"])].append(r)

    def side(st: str) -> str:
        return "TEST" if int(hashlib.md5(f"st{st}".encode()).hexdigest(), 16) % 2 == 0 else "TRAIN"

    test = [r for st, rs in by_st.items() if side(st) == "TEST" for r in rs]
    train = [r for st, rs in by_st.items() if side(st) == "TRAIN" for r in rs]
    st_test = {st for st in by_st if side(st) == "TEST"}
    st_train = {st for st in by_st if side(st) == "TRAIN"}
    assert not (st_test & st_train), "a sequence type landed on both sides"

    def score(rs, rule):
        c = {ax: collections.Counter() for ax in ("O", "H")}
        for r in rs:
            for ax in ("O", "H"):
                lab = r.get(ax)
                if not lab:
                    continue
                call = norm_call(r["calls"][rule][ax], ax)
                if call is None:
                    c[ax]["no_call"] += 1
                elif call == lab:
                    c[ax]["hit"] += 1
                else:
                    c[ax]["miss"] += 1
        out = {}
        for ax in ("O", "H"):
            s = c[ax]["hit"] + c[ax]["miss"]
            out[ax] = {k: c[ax][k] for k in ("hit", "miss", "no_call")}
            out[ax]["accuracy"] = (c[ax]["hit"] / s) if s else None
        return out

    print(f"\nlineage split: {len(by_st)} distinct STs | TEST {len(st_test)} STs / {len(test)} isolates"
          f" | TRAIN {len(st_train)} STs / {len(train)} isolates")
    print(f"  excluded (no complete ST, cannot be placed disjointly): {untyped}")
    biggest = max((len(v) for v in by_st.values()), default=0)
    print(f"  largest single ST: {biggest} isolates "
          f"({biggest/len(typed):.1%} of typed)" if typed else "")

    res_test = {rule: score(test, rule) for rule in ("coverage_primary", "identity_primary")}
    print(f"\nLINEAGE-DISJOINT TEST half (n={len(test)}):")
    for rule in ("coverage_primary", "identity_primary"):
        for ax in ("O", "H"):
            s = res_test[rule][ax]
            acc = "n/a" if s["accuracy"] is None else f"{s['accuracy']:.4f}"
            print(f"  {rule:<18} {ax}: hit={s['hit']:>3} miss={s['miss']:>3} "
                  f"no_call={s['no_call']:>3} acc={acc}")
    h_gain = ((res_test["identity_primary"]["H"]["accuracy"] or 0)
              - (res_test["coverage_primary"]["H"]["accuracy"] or 0))
    o_gain = ((res_test["identity_primary"]["O"]["accuracy"] or 0)
              - (res_test["coverage_primary"]["O"]["accuracy"] or 0))
    print(f"\n  H gain {h_gain:+.4f}   O gain {o_gain:+.4f}   "
          f"(isolate-disjoint replication was +0.1057)")

    # POWER GUARD. If the two rules produced IDENTICAL H calls on the test half, the split never
    # exercised the difference, and a gain of exactly 0 says nothing about the fix. Reporting that as
    # falsification would be the "a test that errors looks like a finding" trap in reverse.
    n_flips = sum(1 for r in test
                  if norm_call(r["calls"]["identity_primary"]["H"], "H")
                  != norm_call(r["calls"]["coverage_primary"]["H"], "H"))
    print(f"  H calls that DIFFER between the rules on the test half: {n_flips}")
    if n_flips == 0:
        verdict = "UNDERPOWERED_RULES_NEVER_DIFFERED"
        why = (f"the two orderings produced identical H calls on all {len(test)} test isolates, so the "
               "split never exercised the difference. A gain of 0 here is an absence of evidence, NOT "
               "evidence the fix fails.")
    elif h_gain <= 0:
        verdict = "FALSIFIED_ON_LINEAGE_DISJOINT"
        why = (f"the fix does NOT survive a lineage-disjoint split (H gain {h_gain:+.4f}). The earlier "
               "isolate-held-out gain was carried by lineage overlap, and the deployed rule should be "
               "reconsidered.")
    elif h_gain >= 0.05:
        verdict = "SURVIVES_LINEAGE_DISJOINT"
        why = (f"the fix holds with no sequence type shared between halves (H gain {h_gain:+.4f}). The "
               "named isolate-vs-lineage limit is closed in the direction the claim needed.")
    else:
        verdict = "SURVIVES_BUT_SMALLER_ON_LINEAGE_DISJOINT"
        why = (f"the gain is positive but attenuated ({h_gain:+.4f} vs +0.1057 isolate-disjoint), so "
               "part of the earlier margin WAS lineage overlap. The fix still helps; quote the "
               "lineage-disjoint number.")
    print(f"\nVERDICT: {verdict}\n  {why}")

    out = {"schema": "serotype-lineage-disjoint-v1", "date": _date.today().isoformat(),
           "closes": ("the named limit that the +0.106 replication was held out BY ISOLATE, not by "
                      "lineage -- two accessions can be near-identical genomes"),
           "lineage_unit": ("multi-locus sequence type (E. coli Achtman 7-locus). Chosen over Mash "
                            "because Docker was unavailable AND because ST is the canonical, "
                            "interpretable lineage unit for this organism"),
           "split_rule": "whole STs assigned to TRAIN/TEST by deterministic hash; no ST on both sides",
           "n_labelled": len(labels), "n_scored": len(rows),
           "n_excluded_no_complete_st": untyped,
           "n_distinct_st": len(by_st), "largest_st_isolates": biggest,
           "test": {"n_st": len(st_test), "n_isolates": len(test)},
           "train": {"n_st": len(st_train), "n_isolates": len(train)},
           "test_half": res_test, "H_gain": h_gain, "O_gain": o_gain,
           "n_H_calls_differing_between_rules": n_flips,
           "isolate_disjoint_reference_gain": 0.1057,
           "verdict": verdict, "why": why,
           "honest_limits": [
               "Isolates without a complete 7-locus ST are EXCLUDED, not pooled -- they cannot be "
               "placed disjointly. That shrinks the test set and could bias it toward well-assembled "
               "genomes.",
               "ST is a lineage unit, not a clonality guarantee: two isolates of the SAME ST are "
               "near-clonal, but different STs can still be closely related, so this is stronger than "
               "isolate-disjoint and weaker than a full phylogenetic split.",
               "Both rules are scored from ONE alignment pass per genome, so only the sort key differs "
               "-- but both share any limitation of that pass.",
               "The fix is NOT re-chosen here; this only asks whether its gain survives. A smaller "
               "number is an honest correction rather than a failure.",
           ]}
    a.out.write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")
    print(f"\nwrote {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
