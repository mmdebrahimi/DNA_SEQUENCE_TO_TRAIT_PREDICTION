"""Hunt S-labelled `rmt` carriers in NCBI-PD's OWN genotype field -- the deployed v2 rule's untested risk.

WHAT IS UNTESTED. The gentamicin v2 rule (deployed 2026-08-31 under a new lock) widens the frozen
`subclass_any={GENTAMICIN}` with `symbol_rescue=^(rmt[A-H]\\d*|npmA\\d*)$`. Its sensitivity gain is measured
(0.523 -> 0.892 on 131 leakage-gated isolates). Its SPECIFICITY claim is not: no S-labelled `rmt` carrier
exists in any dataset checked, so "specificity unchanged" is an ABSENCE, not a bound. One over-call would
be a real false positive on a deployed rule.

WHY THE PREVIOUS HUNT COULD NOT SETTLE IT -- a STRUCTURAL limit, not an effort one.
`gentamicin_rmt_label_hunt.py` keeps only PD rows whose `asm_acc` is already in the LOCAL AMRFinder cache.
It therefore answers "of MY 109 carriers, which are labelled?" (63 were; 62 R, and the 1 S carries `armA`
only, which the frozen rule already called R). It can never reach a carrier we have not run AMRFinder on,
so a counter-example living outside that cache is INVISIBLE to it -- and that is exactly where one would be.

THE INVERSION. PD metadata ships `AMR_genotypes` -- NCBI's OWN AMRFinder gene calls -- in the same row as
the measured `AST_phenotypes`. So the question can be asked the other way round, over EVERY labelled
isolate rather than over our cache: of all PD isolates with a gentamicin AST label, which carry `rmt`, and
what is their phenotype? That population is orders of magnitude larger and is not cache-bounded.

The rescue regex is IMPORTED from the deployed rule, never restated -- a hunt scored against a re-typed
pattern would not be testing what ships.

Offline-safe: network only. Writes wiki/gentamicin_rmt_specificity_hunt.json.
"""
from __future__ import annotations

import argparse
import json
import re
import sys
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from dna_decode.data.pd_ast import ast_label_for  # noqa: E402
from dna_decode.eval.amr_rules import DRUG_RULE  # noqa: E402

BASE = "https://ftp.ncbi.nlm.nih.gov/pathogen/Results/{group}/latest_snps/Metadata/"

# Groups where 16S methyltransferases actually circulate. Wider than the previous hunt's four, because
# the point is to reach carriers the cache-bounded search could not.
DEFAULT_GROUPS = ("Klebsiella", "Escherichia_coli_Shigella", "Acinetobacter", "Salmonella",
                  "Enterobacter_hormaechei", "Enterobacter_cloacae", "Pseudomonas_aeruginosa",
                  "Serratia", "Citrobacter_freundii", "Providencia", "Morganella")

# THE DEPLOYED PATTERN -- imported, not restated.
RESCUE = DRUG_RULE["gentamicin"].get("symbol_rescue")
if not RESCUE:
    raise SystemExit("gentamicin rule carries no symbol_rescue; this hunt targets the v2 rule.")
RESCUE_RE = re.compile(RESCUE)
ARMA_RE = re.compile(r"^armA\d*$")


def parse_amr_genotypes(field: str | None) -> set[str]:
    """`AMR_genotypes` -> bare gene symbols.

    Two traps, both real in the live field:
      - the whole cell is wrapped in literal double quotes, so the quote rides on the FIRST and LAST
        tokens (the identical trap already fixed for `AST_phenotypes`);
      - each symbol may carry an `=SUFFIX` annotation (`=PARTIAL`, `=POINT`, `=MISTRANSLATION`,
        `=PARTIAL_END_OF_CONTIG`), which must be stripped before matching a symbol regex.
    """
    if not field or field.strip().upper() in {"", "NULL", "NA"}:
        return set()
    out = set()
    for tok in field.strip().strip('"').split(","):
        sym = tok.strip().strip('"').split("=", 1)[0].strip()
        if sym:
            out.add(sym)
    return out


def latest_metadata_url(group: str) -> str:
    base = BASE.format(group=group)
    html = urllib.request.urlopen(base, timeout=180).read().decode("utf8", "replace")
    files = sorted(set(re.findall(r"(PDG[0-9.]+\.metadata\.tsv)", html)))
    if not files:
        raise RuntimeError(f"no metadata TSV listed for {group}")
    return base + files[-1]


def scan_group(group: str, drug: str, row_cap: int | None) -> dict:
    url = latest_metadata_url(group)
    r = urllib.request.urlopen(url, timeout=600)
    cols = r.readline().decode("utf8", "replace").rstrip("\n").split("\t")
    try:
        gi, ai = cols.index("AMR_genotypes"), cols.index("AST_phenotypes")
    except ValueError:
        return {"group": group, "error": "missing AMR_genotypes/AST_phenotypes columns"}
    acc_i = cols.index("asm_acc") if "asm_acc" in cols else None
    idx = {c: i for i, c in enumerate(cols)}

    n_rows = n_labelled = 0
    rmt = {"R": [], "S": [], "I": []}
    arma_only = {"R": [], "S": [], "I": []}
    for line in r:
        n_rows += 1
        if row_cap and n_rows > row_cap:
            break
        f = line.decode("utf8", "replace").rstrip("\n").split("\t")
        if len(f) <= max(gi, ai):
            continue
        label = ast_label_for(f[ai], drug)
        if label not in ("R", "S", "I"):
            continue
        n_labelled += 1
        syms = parse_amr_genotypes(f[gi])
        acc = f[acc_i] if acc_i is not None and len(f) > acc_i else ""
        has_rmt = any(RESCUE_RE.match(s) for s in syms)
        has_arma = any(ARMA_RE.match(s) for s in syms)
        if has_rmt:
            rec = {"acc": acc or "?",
                   "rmt": sorted(s for s in syms if RESCUE_RE.match(s)),
                   "has_gent_aac3": any(s.startswith("aac(3)") for s in syms)}
            for col in ("bioproject_acc", "bioproject_center", "collected_by", "sra_center",
                        "epi_type", "collection_date", "geo_loc_name"):
                j = idx.get(col)
                if j is not None and len(f) > j:
                    rec[col] = f[j]
            rmt[label].append(rec)
        elif has_arma:
            arma_only[label].append(acc or "?")

    return {"group": group, "url": url, "n_rows_scanned": n_rows, "n_labelled": n_labelled,
            "rmt_counts": {k: len(v) for k, v in rmt.items()},
            "rmt_S_records": rmt["S"], "rmt_R_records": rmt["R"], "rmt_I_records": rmt["I"],
            "arma_only_counts": {k: len(v) for k, v in arma_only.items()},
            "arma_only_S_accessions": arma_only["S"][:50],
            "truncated": bool(row_cap and n_rows > row_cap)}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--groups", default=",".join(DEFAULT_GROUPS))
    ap.add_argument("--drug", default="gentamicin")
    ap.add_argument("--row-cap", type=int, default=None,
                    help="smoke only; a capped run CANNOT claim a specificity bound")
    ap.add_argument("--out", type=Path,
                    default=Path(__file__).resolve().parents[1] / "wiki" /
                    "gentamicin_rmt_specificity_hunt.json")
    a = ap.parse_args()

    print(f"deployed rescue pattern (imported): {RESCUE}\n")
    per_group, errors = [], {}
    for g in [x.strip() for x in a.groups.split(",") if x.strip()]:
        try:
            res = scan_group(g, a.drug, a.row_cap)
        except Exception as e:                       # a fetch failure must not read as "no carriers"
            errors[g] = f"{type(e).__name__}: {e}"
            print(f"{g:26} ERROR {errors[g][:70]}", flush=True)
            continue
        if "error" in res:
            errors[g] = res["error"]
            print(f"{g:26} ERROR {res['error']}", flush=True)
            continue
        per_group.append(res)
        c = res["rmt_counts"]
        print(f"{g:26} labelled={res['n_labelled']:6d}  rmt R/S/I = "
              f"{c['R']:4d}/{c['S']:3d}/{c['I']:3d}   armA-only S={res['arma_only_counts']['S']}",
              flush=True)

    tot = {k: sum(g["rmt_counts"][k] for g in per_group) for k in ("R", "S", "I")}
    s_recs = [r for g in per_group for r in g["rmt_S_records"]]
    r_recs = [r for g in per_group for r in g["rmt_R_records"]]
    s_accs = [r["acc"] for r in s_recs]
    n_lab = sum(g["n_labelled"] for g in per_group)
    complete = not errors and not a.row_cap

    print(f"\nTOTAL over {len(per_group)} groups: {n_lab:,} gentamicin-labelled isolates")
    print(f"  rmt carriers  R={tot['R']}  S={tot['S']}  I={tot['I']}")
    if tot["R"] + tot["S"]:
        print(f"  PPV(rmt -> R) = {tot['R']}/{tot['R'] + tot['S']} = "
              f"{tot['R'] / (tot['R'] + tot['S']):.4f}")
    if s_accs:
        print(f"  S-labelled rmt carriers (the sought counter-examples): {s_accs[:20]}")
    else:
        print("  NO S-labelled rmt carrier found.")

    out = {
        "schema": "gentamicin-rmt-specificity-hunt-v1",
        "drug": a.drug,
        "deployed_rescue_pattern": RESCUE,
        "population": "ALL PD isolates carrying a gentamicin AST label -- NOT restricted to the local "
                      "AMRFinder cache, which is what the previous hunt could not escape",
        "n_labelled_isolates": n_lab,
        "rmt_carrier_counts": tot,
        "rmt_S_accessions": s_accs,
        "rmt_S_records": s_recs,
        "rmt_R_records": r_recs,
        "ppv_rmt_to_R": (tot["R"] / (tot["R"] + tot["S"])) if (tot["R"] + tot["S"]) else None,
        "per_group": per_group,
        "errors": errors,
        "complete": complete,
        "honest_limits": [
            "The CARRIER call is NCBI's own AMRFinder run (PD's AMR_genotypes), not ours -- a different "
            "version/DB. That makes it independent of our pipeline, but it is still a tool-derived "
            "feature; only the PHENOTYPE is measured.",
            "Same archive (NCBI-PD) as the earlier hunt, but a DIFFERENT population: every labelled "
            "isolate rather than our 109 cached carriers. It can reach counter-examples the previous "
            "search structurally could not.",
            "A run with errors or --row-cap is NOT a specificity bound; `complete` records which.",
            "An S label is a phenotype claim from a submitting lab; a single S carrier warrants "
            "inspection (co-carriage, method, breakpoint) before it is treated as a rule defect.",
        ],
    }
    a.out.write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")
    print(f"\ncomplete={complete}  wrote {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
