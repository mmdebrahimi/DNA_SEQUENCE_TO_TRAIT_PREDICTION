"""Which determinant families does a deployed rule NOT count? The L2 "doubt" screen.

WHY THIS EXISTS. The catalog's failure mode is COMPLETENESS, not accuracy, and it has now been found
twice with the same shape:

  gentamicin  the rule matches AMRFinder `Subclass == GENTAMICIN`; `rmt*` files under the generic
              `AMINOGLYCOSIDE` subclass -> invisible. sens 0.523 vs 0.893; 24 of 31 FN carry `rmt`.
  HIV NNRTI   the catalog is scoped to 8 major positions; 53 resistant isolates carry drivers outside
              them.

Both were invisible until an independent label set arrived. Neither is a model failure. This screen is
the generalisation: **a determinant family present in the data but unrepresentable by the rule**.

WHAT IT IS NOT. It does not predict resistance and it never emits a call. A determinant the rule does not
count is a CANDIDATE gap for human review, nothing more -- and many are excluded ON PURPOSE (blaTEM-1 is
correctly not ceftriaxone-R; `aph`/`aadA` are correctly not gentamicin-R). Ranking is what makes the
output usable: a family carried by many R-labelled isolates and NO S-labelled ones has the `rmt`
signature; one carried by both is probably a correct exclusion.

HOW IT ASKS. It does NOT reimplement the rule. For each distinct determinant it writes a ONE-ROW table
carrying the ORIGINAL header and the row VERBATIM, then asks the deployed `call_resistance` whether that
alone yields R. Re-deriving the rule's logic here would drift from `DRUG_RULE` the moment either changed;
probing the deployed function cannot.

Read-only. No network, no Docker, no model. Frozen surface untouched.

Run: uv run python scripts/determinant_completeness_screen.py [--drug gentamicin] [--limit N]
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
import tempfile
from collections import defaultdict
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))
WIKI = ROOT / "wiki"

REGISTRY_ORGANISM = {"Escherichia_coli_Shigella": "Escherichia_coli_Shigella",
                     "Klebsiella": "Klebsiella_pneumoniae"}


def determinant_key(row: dict) -> tuple[str, str, str]:
    """The identity of a determinant FAMILY. PURE.

    Keyed on (symbol, Class, Subclass) -- NOT on the allele. `rmtB` and `rmtE1` are different alleles of
    one blind family, and collapsing them to a bare prefix would over-merge unrelated genes. The Subclass
    is part of the key because it is precisely what the rule matches on, so two rows sharing a symbol but
    differing in Subclass are genuinely different cases for this question.
    """
    return ((row.get("Element symbol") or "").strip(),
            (row.get("Class") or "").strip().upper(),
            (row.get("Subclass") or "").strip().upper())


def rank_candidates(families: dict, min_carriers: int = 1) -> list[dict]:
    """Rank uncounted families by how much they look like a real gap. PURE.

    The `rmt` signature is: carried by many R-labelled isolates and NO S-labelled ones. A family carried
    by both classes is probably a CORRECT exclusion, so it sorts down. Families with no labelled carriers
    at all sort last -- they are unassessable here, not innocent.
    """
    out = []
    for key, rec in families.items():
        r, s = rec["r_carriers"], rec["s_carriers"]
        labelled = r + s
        # purity is undefined without labels; treat it as 0 so unlabelled families cannot outrank
        # a family with real evidence behind it.
        purity = (r / labelled) if labelled else 0.0
        out.append({"symbol": key[0], "class": key[1], "subclass": key[2],
                    "n_genomes": rec["n_genomes"], "r_carriers": r, "s_carriers": s,
                    "r_purity": round(purity, 3),
                    "signature": ("rmt_like" if r >= 3 and s == 0
                                  else "mixed" if labelled else "unlabelled")})
    out = [c for c in out if c["n_genomes"] >= min_carriers]
    # SIGNATURE FIRST, then volume. Sorting by raw R-count buries the actionable family under prevalent
    # but MIXED ones: the first run ranked rmtE1 (36R/0S -- the known gap) 5th, beneath aph/aadA at
    # 62R/28S, which are CORRECT exclusions. Purity is what separates a gap from a deliberate exclusion,
    # so it leads.
    _RANK = {"rmt_like": 0, "mixed": 1, "unlabelled": 2}
    out.sort(key=lambda c: (_RANK[c["signature"]], -c["r_carriers"], -c["r_purity"],
                            -c["n_genomes"], c["symbol"]))
    return out


def rule_counts_determinant(header: list[str], row: dict, drug: str, organism: str | None) -> bool:
    """Can the DEPLOYED rule represent this determinant at all?

    Writes a table with the original header and the row verbatim, then asks `call_resistance`. Nothing
    here mirrors the rule's logic, so nothing here can drift from it.

    THE ROW IS REPEATED TO THE RULE'S THRESHOLD, and that is load-bearing. A one-row probe against a
    threshold-2 rule can NEVER return R, so the first run flagged every QRDR point mutation as an
    uncounted "gap" -- ciprofloxacin reported 0 of 51 determinants counted, with parC_S80I at 60R/0S on
    top -- when the rule counts them perfectly and simply requires TWO. That conflates "the rule cannot
    REPRESENT this determinant" with "the rule needs more than one of them", and only the first is a
    completeness gap.
    """
    from dna_decode.eval.amr_rules import call_resistance, rule_for

    n_needed = max(1, int(rule_for(drug).get("threshold") or 1))
    with tempfile.TemporaryDirectory() as td:
        p = Path(td) / "main.tsv"
        with p.open("w", encoding="utf-8", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=header, delimiter="	", extrasaction="ignore")
            w.writeheader()
            for _ in range(n_needed):
                w.writerow(row)
        try:
            return call_resistance(p, drug, organism=organism)["prediction"] == "R"
        except Exception:
            return False       # a determinant the rule cannot even evaluate is, for us, not counted


def screen_drug(drug: str, idx: dict, labels: dict, limit: int | None = None) -> dict:
    """Aggregate uncounted drug-relevant determinant families across the cached genomes."""
    from dna_decode.data.mic_tiers import amrfinder_classes_for

    classes = {c.upper() for c in amrfinder_classes_for(drug)}
    counted_cache: dict[tuple, bool] = {}
    families: dict[tuple, dict] = defaultdict(lambda: {"n_genomes": 0, "r_carriers": 0, "s_carriers": 0})
    n_scanned = 0

    # LABELLED GENOMES FIRST. The ranking is driven by R/S carrier counts, so scanning in arbitrary
    # index order makes a truncated run rank on raw prevalence instead -- the first smoke run surfaced
    # `aph`/`aadA` at STREPTOMYCIN subclass (CORRECT exclusions) purely because the first 250 accessions
    # happened to be unlabelled. Order by "has a label for THIS drug" so a capped run sees the evidence.
    ordered = sorted(idx.items(),
                     key=lambda kv: (0 if (labels.get(kv[0], {}).get("calls", {}) or {}).get(drug)
                                     else 1, kv[0]))
    for acc, tsv in ordered[:limit]:
        p = Path(tsv)
        if not p.exists():
            continue
        n_scanned += 1
        lab = (labels.get(acc, {}).get("calls", {}) or {}).get(drug, "")
        org = REGISTRY_ORGANISM.get(labels.get(acc, {}).get("group", ""))
        try:
            with p.open(encoding="utf-8") as fh:
                rdr = csv.DictReader(fh, delimiter="\t")
                header = rdr.fieldnames or []
                seen_here: set[tuple] = set()
                for row in rdr:
                    cls = (row.get("Class") or "").upper()
                    sub = (row.get("Subclass") or "").upper()
                    if not any(c in cls or c in sub for c in classes):
                        continue                      # not drug-relevant at all: out of scope
                    key = determinant_key(row)
                    if key not in counted_cache:
                        counted_cache[key] = rule_counts_determinant(header, row, drug, org)
                    if counted_cache[key]:
                        continue                      # the rule sees it; not a gap
                    if key in seen_here:
                        continue                      # count each genome once per family
                    seen_here.add(key)
                    fam = families[key]
                    fam["n_genomes"] += 1
                    if lab.upper() == "R":
                        fam["r_carriers"] += 1
                    elif lab.upper() == "S":
                        fam["s_carriers"] += 1
        except OSError:
            continue

    return {"drug": drug, "n_genomes_scanned": n_scanned,
            "n_distinct_determinants_probed": len(counted_cache),
            "n_counted_by_rule": sum(1 for v in counted_cache.values() if v),
            "candidates": rank_candidates(families)}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--drug", action="append", help="repeatable; default = every deployed DRUG_RULE drug")
    ap.add_argument("--limit", type=int, default=None, help="cap genomes scanned (smoke runs)")
    args = ap.parse_args()

    from gentamicin_rmt_candidate import amrfinder_index

    from dna_decode.eval.amr_rules import DRUG_RULE

    idx = amrfinder_index()
    if not idx:
        print("no cached AMRFinder output found -- nothing to screen")
        return 1

    lf = WIKI / "unscored_genome_label_census.json"
    labels = {}
    if lf.exists():
        d = json.loads(lf.read_text(encoding="utf-8"))
        labels = d.get("labels") or {}

    drugs = args.drug or sorted(DRUG_RULE)
    out = {"schema": "determinant-completeness-screen-v1", "generated": date.today().isoformat(),
           "n_cached_genomes": len(idx), "n_labelled_genomes": len(labels),
           "contract": ("A DOUBT signal, never a call. An uncounted determinant is a CANDIDATE gap for "
                        "human review; many exclusions are deliberate and correct."),
           "drugs": []}

    for drug in drugs:
        res = screen_drug(drug, idx, labels, args.limit)
        out["drugs"].append(res)
        top = res["candidates"][:5]
        print(f"\n{drug}: scanned {res['n_genomes_scanned']} genomes | "
              f"{res['n_distinct_determinants_probed']} distinct drug-relevant determinants "
              f"({res['n_counted_by_rule']} counted by the rule)")
        if not top:
            print("   no uncounted drug-relevant determinant families")
        for c in top:
            print(f"   {c['signature']:11} {c['symbol']:14} {c['subclass'][:26]:28} "
                  f"genomes={c['n_genomes']:>4} R={c['r_carriers']:>3} S={c['s_carriers']:>3}")

    (WIKI / f"determinant_completeness_screen_{out['generated']}.json").write_text(
        json.dumps(out, indent=2), encoding="utf-8")
    print(f"\nwrote wiki/determinant_completeness_screen_{out['generated']}.json")
    print("DOUBT signal only -- an uncounted family is a candidate for review, NOT a resistance call.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
