"""Find S-labelled 16S-methyltransferase carriers — the one measurement the candidate rule still needs.

WHY. The `rmt` candidate for gentamicin rescues 1 of 1 rescuable false negative and changes no other call,
but its specificity result is VACUOUS: every methyltransferase carrier in the labelled gentamicin data is
R, so the rule *cannot* produce a false positive there. Its over-calling risk is UNTESTED, not zero.

Only an **S-labelled carrier** can settle that. This repo has 109 genomes with a cached methyltransferase
call — they just aren't labelled for gentamicin locally. NCBI Pathogen Detection carries `AST_phenotypes`
for many of them, free, and is the same source the 10 frozen SCORED cells came from.

So: stream the PD metadata for the relevant organism groups, keep only rows whose `asm_acc` is one of our
carriers, and read the gentamicin call out of `AST_phenotypes`.

REUSES THE FIXED PARSER. `pd_ast.parse_ast_phenotypes` handles the quoted comma-separated field whose
quotes ride on the first and last token — the end-position bug that silently cost the census ~22R+26R on
two drugs. Re-implementing that split here would re-introduce it.

HONEST SCOPE. A carrier that PD has no AST for is simply unmeasured; absence is never read as S. And the
result bounds over-calling only over the carriers PD happens to label — it is a floor on the evidence, not
a population estimate.

Run: uv run python scripts/gentamicin_rmt_label_hunt.py [--groups Klebsiella,...]
"""
from __future__ import annotations

import argparse
import csv
import gzip
import io
import json
import sys
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))

# Where the carriers actually live, so we stream only groups that can contain them.
DEFAULT_GROUPS = ("Klebsiella", "Escherichia_coli_Shigella", "Acinetobacter", "Salmonella")


def carrier_accessions() -> dict[str, list[str]]:
    """accession -> methyltransferase gene symbols, from cached AMRFinder output."""
    from gentamicin_rmt_candidate import (amrfinder_index, gene_symbol,
                                          is_methyltransferase, read_rows)
    out: dict[str, list[str]] = {}
    for acc, main in amrfinder_index().items():
        mt = [gene_symbol(r) for r in read_rows(main) if is_methyltransferase(gene_symbol(r))]
        if mt:
            out[acc] = sorted(set(mt))
    return out


def stream_group(group: str, wanted: set[str], drug: str = "gentamicin") -> dict[str, str]:
    """{asm_acc: R/S/...} for wanted accessions in one PD group. Network."""
    from ncbi_pd_provenance_census import latest_metadata_url
    from dna_decode.data.pd_ast import parse_ast_phenotypes

    url = latest_metadata_url(group)
    found: dict[str, str] = {}
    with urllib.request.urlopen(url, timeout=600) as resp:
        raw = gzip.GzipFile(fileobj=resp) if url.endswith(".gz") else resp
        text = io.TextIOWrapper(raw, encoding="utf-8", errors="replace")
        rd = csv.reader(text, delimiter="\t")
        header = next(rd, None)
        if not header:
            return found
        idx = {c: i for i, c in enumerate(header)}
        ai, pi = idx.get("asm_acc"), idx.get("AST_phenotypes")
        if ai is None or pi is None:
            return found
        for row in rd:
            if len(row) <= max(ai, pi):
                continue
            acc = row[ai].strip()
            if acc not in wanted:
                continue
            call = parse_ast_phenotypes(row[pi], {drug}).get(drug)
            if call:
                found[acc] = call
    return found


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--groups", default=",".join(DEFAULT_GROUPS))
    ap.add_argument("--drug", default="gentamicin")
    args = ap.parse_args()

    carriers = carrier_accessions()
    print(f"methyltransferase carriers with cached AMRFinder: {len(carriers)}")
    if not carriers:
        return 1

    labels: dict[str, str] = {}
    per_group: dict[str, int] = {}
    errors: dict[str, str] = {}
    for g in [x.strip() for x in args.groups.split(",") if x.strip()]:
        try:
            got = stream_group(g, set(carriers), args.drug)
        except Exception as exc:                      # network/format — record, never silently zero
            errors[g] = f"{type(exc).__name__}: {exc}"
            print(f"  {g}: FAILED ({errors[g]})")
            continue
        per_group[g] = len(got)
        labels.update(got)
        print(f"  {g}: {len(got)} carrier(s) with a {args.drug} call")

    r = sorted(a for a, v in labels.items() if v.upper() == "R")
    s = sorted(a for a, v in labels.items() if v.upper() == "S")
    print()
    print(f"{args.drug} calls found for carriers: {len(labels)}  (R {len(r)}, S {len(s)})")

    if errors:
        # A partial sweep cannot bound anything -- say so rather than reporting a floor as a result.
        print(f"INCOMPLETE: {len(errors)} group(s) failed; the S-count is a partial floor, not a bound.")

    if s:
        from gentamicin_rmt_candidate import amrfinder_index, read_rows, frozen_call, candidate_call
        idx = amrfinder_index()
        flips = [a for a in s
                 if not frozen_call(read_rows(idx[a])) and candidate_call(read_rows(idx[a]))]
        print(f"S-labelled carriers the candidate would newly call R (FALSE POSITIVES): "
              f"{len(flips)}/{len(s)}")
        for a in flips[:10]:
            print(f"  {a}  {carriers[a]}")
        print()
        print(f"=> over-calling cost, MEASURED on {len(s)} S-labelled carrier(s): {len(flips)} FP")
    elif errors:
        # Do NOT draw a conclusion from a sweep that failed. Printing "still zero" right after
        # "INCOMPLETE" would be the reassuring-verdict pattern: a confident line that the run above
        # already contradicted.
        print("=> NO CONCLUSION. The sweep did not complete, so 'zero S-labelled carriers' here means")
        print("   'the sweep failed', not 'none exist'. Re-run before reading anything into it.")
    else:
        print("=> still ZERO S-labelled carriers across a COMPLETE sweep. The over-calling risk remains")
        print("   UNTESTED, and that is now a property of the available public labels rather than of")
        print("   this repo's local data.")

    out = {"drug": args.drug, "n_carriers": len(carriers), "n_labelled": len(labels),
           "n_R": len(r), "n_S": len(s), "per_group": per_group, "errors": errors,
           "s_labelled": s, "complete": not errors}
    (ROOT / "wiki" / "gentamicin_rmt_label_hunt.json").write_text(
        json.dumps(out, indent=2), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
