"""Does the shipped PointFinder caller agree with an INDEPENDENT curated caller on QRDR point mutations?

`finder:Escherichia_coli:pointfinder` is the last untested finder cell, registered FAITHFUL_TO_TOOL --
checked against the reference METHOD, never against reality. It is the exact complement of the ResFinder
locus-collapse work, which compared ACQUIRED GENES and could not touch point mutations at all.

WHY THIS ONE MATTERS BEYOND THE CELL. The catalogued genes are gyrA / gyrB / parC / parE -- the QRDR --
and the FROZEN cipro decoder rule (`qrdr_point`) consumes exactly these determinants from AMRFinder. So
an independent re-derivation of those calls touches a determinant class the DEPLOYED surface depends on.
This script does NOT modify the frozen surface; it is read-only.

THE COMPARATOR, AND WHERE IT ACTUALLY LIVES. AMRFinder's point-mutation screen is NOT in `main.tsv` --
measured across all 1818 committed runs, main.tsv `Type` is only AMR/STRESS with ZERO POINT rows. The
screen is in a separate `mutations.tsv`, in `--mutation_all` form, so it contains the WILDTYPE rows too
(a row per screened position, most of them "no mutation here"). Those are filtered out; only genuine
non-wildtype calls are compared.

THE POSITION RESTRICTION, AND WHY IT IS NOT A FREE PASS. PointFinder can only report positions in its
own `resistens-overview.txt`, so scoring it against AMRFinder calls at UNCATALOGUED positions would
charge it with missing something it structurally cannot express. The comparison is therefore restricted
to catalogued (gene, codon) pairs. Having just been burned by a control that excluded nothing, this
script REPORTS how many AMRFinder calls the restriction actually removes and REFUSES to claim the
restriction is meaningful if it removes none.

Offline: cached assemblies + the committed PointFinder DB + native blastn + committed AMRFinder runs.
"""
from __future__ import annotations

import argparse
import collections
import csv
import json
import re
import sys
from datetime import date as _date
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from dna_decode.pointfinder.runner import call_point_mutations, parse_overview  # noqa: E402
from gentamicin_rmt_candidate import amrfinder_index  # noqa: E402

BLASTN = "C:/Users/Farshad/ncbi-blast/bin/blastn.exe"
POINTFINDER_DB = ROOT / "data" / "pointfinder_db" / "escherichia_coli"

# `S83L` -> ('S', 83, 'L'); also handles multi-residue forms like `DVADD871EPEEE`.
_MUT_RE = re.compile(r"^(?P<ref>[A-Za-z*]+)(?P<pos>\d+)(?P<alt>[A-Za-z*]+)$")


def split_symbol(symbol: str) -> tuple[str, str] | None:
    """'gyrA_S83L' -> ('gyrA', 'S83L'). rsplit, because gene names themselves contain underscores
    (e.g. '16S_rrsC_A226G' must give gene '16S_rrsC', not '16S')."""
    if "_" not in symbol:
        return None
    gene, mut = symbol.rsplit("_", 1)
    return (gene.strip(), mut.strip())


def codon_of(mutation: str) -> int | None:
    m = _MUT_RE.match(mutation)
    return int(m.group("pos")) if m else None


def amrfinder_point_calls(mutations_tsv: Path, genes: set[str]) -> list[tuple[str, str, int]]:
    """Non-wildtype point calls in `genes`, as (gene, mutation, codon).

    mutations.tsv is AMRFinder's --mutation_all output: it lists every SCREENED position, marking the
    ones with no mutation '[WILDTYPE]' in Element name. Keeping those would count a screened-and-clean
    position as a called mutation.
    """
    out = []
    if not mutations_tsv.exists():
        return out
    with open(mutations_tsv, encoding="utf-8", errors="replace") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if "[WILDTYPE]" in (row.get("Element name") or ""):
                continue
            parts = split_symbol((row.get("Element symbol") or "").strip())
            if not parts:
                continue
            gene, mut = parts
            if gene not in genes:
                continue
            pos = codon_of(mut)
            if pos is not None:
                out.append((gene, mut, pos))
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--refseq-cache", type=Path, default=Path("D:/dna_decode_cache/refseq"))
    ap.add_argument("--db", type=Path, default=POINTFINDER_DB)
    ap.add_argument("--limit", type=int, default=0, help="0 = all")
    ap.add_argument("--blastn", default=BLASTN)
    ap.add_argument("--checkpoint", type=Path,
                    default=Path("D:/dna_decode_cache/pointfinder_concordance/rows.jsonl"))
    ap.add_argument("--out", type=Path, default=ROOT / "wiki" /
                    f"pointfinder_amrfinder_concordance_{_date.today().isoformat()}.json")
    a = ap.parse_args()

    if not a.db.exists():
        print(f"PointFinder DB absent at {a.db}", file=sys.stderr)
        return 2
    refs = {p.stem: p for p in a.db.glob("*.fsa")}
    overview = parse_overview(a.db / "resistens-overview.txt")
    genes = set(refs)
    # The catalogued (gene, codon) pairs we can actually call: a catalogue entry for a gene with no
    # reference CDS on disk (rpoB, 16S, ampC promoter...) is unreachable for this caller.
    catalogued = {(g, p) for (g, p) in overview if g in genes}
    print(f"PointFinder reference genes: {sorted(genes)}")
    print(f"catalogued (gene,codon) reachable with a reference CDS: {len(catalogued)}"
          f"  (catalogue holds {len(overview)} across all genes)\n")

    idx = amrfinder_index()
    fastas = {}
    for d in sorted(a.refseq_cache.iterdir()) if a.refseq_cache.exists() else []:
        if d.is_dir() and d.name in idx:
            fa = next((p for p in d.glob("*.fna")), None)
            if fa:
                fastas[d.name] = fa
    accs = sorted(fastas)
    if a.limit:
        accs = accs[:a.limit]
    if not accs:
        print("no accession has BOTH a cached assembly and an AMRFinder run", file=sys.stderr)
        return 2
    print(f"{len(accs)} accessions with both a cached assembly and a committed AMRFinder run\n")

    a.checkpoint.parent.mkdir(parents=True, exist_ok=True)
    done: dict[str, dict] = {}
    if a.checkpoint.exists():
        for line in a.checkpoint.read_text(encoding="utf-8").splitlines():
            if line.strip():
                r = json.loads(line)
                done[r["acc"]] = r
        print(f"  resuming: {len(done)} already scored")
    fh = open(a.checkpoint, "a", encoding="utf-8")

    for n, acc in enumerate(accs, 1):
        if acc in done:
            continue
        rec: dict = {"acc": acc}
        try:
            res = call_point_mutations(fastas[acc], refs, overview, blastn_bin=a.blastn, timeout=600)
        except Exception as e:                       # noqa: BLE001
            rec["status"] = f"error:{type(e).__name__}"
            fh.write(json.dumps(rec) + "\n"); fh.flush(); done[acc] = rec
            continue
        if res.get("status") != "ok":
            rec["status"] = res.get("status")
            fh.write(json.dumps(rec) + "\n"); fh.flush(); done[acc] = rec
            continue

        pf = sorted({f"{m['gene']}_{m['mutation']}" for m in res["mutations"]})
        amr_all = amrfinder_point_calls(idx[acc].parent / "mutations.tsv", genes)
        amr_restricted = sorted({f"{g}_{m}" for g, m, p in amr_all if (g, p) in catalogued})
        rec.update({"status": "ok", "pointfinder": pf,
                    "genes_aligned": sorted(res.get("genes_aligned") or []),
                    "amrfinder_restricted": amr_restricted,
                    "n_amrfinder_all": len({f"{g}_{m}" for g, m, _p in amr_all}),
                    "n_amrfinder_restricted": len(amr_restricted)})
        fh.write(json.dumps(rec) + "\n"); fh.flush(); done[acc] = rec
        if n % 50 == 0:
            print(f"  [{n}/{len(accs)}] {acc} pf={len(pf)} amr={len(amr_restricted)}", flush=True)
    fh.close()

    all_ok = [done[x] for x in accs if x in done and done[x].get("status") == "ok"]
    if not all_ok:
        print("\nREFUSING: no genome scored successfully.", file=sys.stderr)
        return 3

    # ABSTENTION IS NOT AGREEMENT. On a fragmented draft assembly the reference genes may not align at
    # all; the caller then returns status 'ok' with zero mutations (the CLI does say "genes aligned:
    # none", so this is visible to a human -- but a genome that aligned NOTHING has not agreed with
    # anything). Counting those as an exact set match would score an abstention as a success, inflating
    # agreement by exactly the genomes where the caller did no work.
    rows = [r for r in all_ok if r.get("genes_aligned")]
    abstained = [r for r in all_ok if not r.get("genes_aligned")]
    # An abstention is only harmless if the comparator ALSO found nothing there; when AMRFinder did
    # find catalogued mutations in a genome PointFinder could not align, those are real missed calls.
    abstained_with_amr_calls = [r for r in abstained if r.get("amrfinder_restricted")]
    if not rows:
        print(f"\nREFUSING: {len(all_ok)} genomes ran but NONE aligned a single reference gene, so "
              "nothing was ever compared.", file=sys.stderr)
        return 3

    exact = both = pf_only = amr_only = 0
    n_pf_total = n_amr_total = n_amr_all_total = 0
    jacc = []
    per_mut = collections.defaultdict(lambda: {"both": 0, "pf_only": 0, "amr_only": 0})
    for r in rows:
        p, m = set(r["pointfinder"]), set(r["amrfinder_restricted"])
        n_pf_total += len(p); n_amr_total += len(m); n_amr_all_total += r["n_amrfinder_all"]
        both += len(p & m); pf_only += len(p - m); amr_only += len(m - p)
        exact += (p == m)
        if p | m:
            jacc.append(len(p & m) / len(p | m))
        for k in p & m:
            per_mut[k]["both"] += 1
        for k in p - m:
            per_mut[k]["pf_only"] += 1
        for k in m - p:
            per_mut[k]["amr_only"] += 1

    # --- NON-VACUITY, two separate ways this could report a meaningless number -------------------
    # (1) If neither caller ever called anything, perfect agreement is a statement about a broken
    #     pipeline. (2) If the position restriction removed nothing, calling it a fairness control
    #     would over-claim it -- exactly the error corrected in the ResFinder write-up.
    removed = n_amr_all_total - n_amr_total
    if n_pf_total == 0 and n_amr_total == 0:
        print(f"\nREFUSING: both callers returned zero mutations across {len(rows)} genomes. "
              "Agreement here is a plumbing result, not a finding.", file=sys.stderr)
        return 3
    restriction_is_live = removed > 0

    mean_j = sum(jacc) / len(jacc) if jacc else None
    n_missed_by_abstention = sum(len(r["amrfinder_restricted"]) for r in abstained_with_amr_calls)
    print(f"\n=== {len(all_ok)} genomes ran ===")
    print(f"  aligned >=1 reference gene (SCORED)        : {len(rows)}")
    print(f"  aligned NOTHING (abstained, NOT scored)    : {len(abstained)} "
          f"({len(abstained)/len(all_ok):.1%})")
    print(f"    of those, AMRFinder DID find catalogued  : {len(abstained_with_amr_calls)} "
          f"genomes / {n_missed_by_abstention} calls")
    print(f"\n=== agreement, on the {len(rows)} SCORED genomes ===")
    print(f"  PointFinder calls          : {n_pf_total}")
    print(f"  AMRFinder calls (all)      : {n_amr_all_total}")
    print(f"  AMRFinder calls (catalogued): {n_amr_total}   "
          f"[position restriction removed {removed}]")
    print(f"  agree / PF-only / AMR-only : {both} / {pf_only} / {amr_only}")
    print(f"  genomes with an EXACT set match: {exact}/{len(rows)} = {exact/len(rows):.4f}")
    print(f"  mean per-genome Jaccard    : {mean_j:.4f}" if mean_j is not None else "")

    print("\nper-mutation concordance (top by frequency):")
    for k, v in sorted(per_mut.items(), key=lambda kv: -(kv[1]["both"] + kv[1]["pf_only"]
                                                         + kv[1]["amr_only"]))[:12]:
        tot = v["both"] + v["pf_only"] + v["amr_only"]
        print(f"   {k:<18s} both={v['both']:4d}  PF-only={v['pf_only']:4d}  "
              f"AMR-only={v['amr_only']:4d}   agreement={v['both']/tot:.3f}")

    if mean_j is not None and mean_j >= 0.95 and exact / len(rows) >= 0.90:
        verdict = "POINTFINDER_AGREES_WITH_AN_INDEPENDENT_CALLER"
        why = (f"on catalogued QRDR positions the two independent callers agree on "
               f"{exact}/{len(rows)} genomes exactly (mean Jaccard {mean_j:.4f}). The cell reproduces "
               "an independent implementation of the same determinant class the FROZEN cipro rule "
               "consumes.")
    elif mean_j is not None and mean_j >= 0.80:
        verdict = "POINTFINDER_MOSTLY_AGREES_WITH_RESIDUAL_DISCORDANCE"
        why = (f"mean per-genome Jaccard {mean_j:.4f} with {pf_only} PointFinder-only and {amr_only} "
               "AMRFinder-only calls. Substantial agreement, but the discordant calls are real and "
               "are enumerated per-mutation above rather than averaged away.")
    else:
        verdict = "POINTFINDER_DISAGREES_WITH_THE_INDEPENDENT_CALLER"
        why = (f"mean per-genome Jaccard {mean_j}. The two callers do NOT agree on catalogued "
               "positions, which is a defect signal in one of them and needs diagnosis before the "
               "cell is relied on.")
    print(f"\nVERDICT: {verdict}\n  {why}")

    out = {"schema": "pointfinder-amrfinder-concordance-v1", "date": _date.today().isoformat(),
           "question": ("does the shipped PointFinder caller reproduce an independent curated caller's "
                        "QRDR point-mutation calls on catalogued positions?"),
           "comparator": ("AMRFinder's committed per-genome mutations.tsv for the SAME accession -- NOT "
                          "main.tsv, which carries zero POINT rows. --mutation_all output, so WILDTYPE "
                          "rows are filtered out before comparison. A TOOL, NOT A WET-LAB LABEL."),
           "position_restriction": {
               "why": ("PointFinder can only report positions in its own resistens-overview.txt, so "
                       "scoring it against AMRFinder calls at uncatalogued positions would charge it "
                       "with missing what it structurally cannot express"),
               "n_catalogued_reachable": len(catalogued),
               "n_amrfinder_calls_removed": removed,
               "restriction_is_live": restriction_is_live,
               "note": ("reported explicitly because a restriction that removes nothing is not a "
                        "fairness control -- that exact over-claim was corrected in the ResFinder "
                        "write-up the day before")},
           "abstention": {
               "why_separated": ("on a fragmented draft assembly the reference genes may not align at "
                                 "all; the caller returns status 'ok' with zero mutations (the CLI does "
                                 "print 'genes aligned: none', so it is visible to a human). A genome "
                                 "that aligned NOTHING has not AGREED with anything -- counting it as "
                                 "an exact set match would score an abstention as a success"),
               "n_genomes_ran": len(all_ok), "n_scored": len(rows), "n_abstained": len(abstained),
               "abstention_rate": len(abstained) / len(all_ok),
               "n_abstained_where_amrfinder_found_catalogued_mutations":
                   len(abstained_with_amr_calls),
               "n_calls_missed_via_abstention": n_missed_by_abstention,
               "abstained_accessions": [r["acc"] for r in abstained][:50]},
           "n_genomes": len(rows),
           "n_pointfinder_calls": n_pf_total, "n_amrfinder_calls_all": n_amr_all_total,
           "n_amrfinder_calls_catalogued": n_amr_total,
           "n_agree": both, "n_pointfinder_only": pf_only, "n_amrfinder_only": amr_only,
           "n_genomes_exact_set_match": exact, "exact_set_match_rate": exact / len(rows),
           "mean_per_genome_jaccard": mean_j,
           "per_mutation": {k: v for k, v in sorted(
               per_mut.items(),
               key=lambda kv: -(kv[1]["both"] + kv[1]["pf_only"] + kv[1]["amr_only"]))},
           "verdict": verdict, "why": why,
           "honest_limits": [
               "The comparator is a TOOL, not a wet-lab label. This measures agreement with an "
               "independent implementation, NOT correctness -- both callers could be wrong together, "
               "and both ultimately derive from the same published QRDR literature. The cell stays "
               "FAITHFUL_TO_TOOL.",
               "Scope is the four genes with a committed reference CDS (gyrA/gyrB/parC/parE). The "
               "PointFinder catalogue also lists rpoB / 16S / 23S / ampC-promoter / pmrAB / folP "
               "positions that have NO reference sequence on disk, so this says nothing about them.",
               "Genomes are whatever this project had cached, drawn from AMR cohorts -- ENRICHED for "
               "resistance and not a random sample of the species.",
               "Agreement on a position both callers catalogue does not validate either one's "
               "CATALOGUE. A resistance position missing from both would be invisible here.",
               "Epistasis ('Required_mut' in the overview) is recorded but NOT enforced by the v0 "
               "caller, so a mutation whose resistance effect depends on a partner is still reported "
               "as conferring resistance.",
               "Genomes where no reference gene aligned are EXCLUDED from the agreement metrics and "
               "reported separately as abstentions. The agreement figures therefore describe the "
               "caller WHEN IT ALIGNS, and the abstention rate is the separate coverage question -- "
               "quoting the agreement alone would overstate the cell.",
               "A discordance at a catalogued POSITION can still be a catalogue difference rather than "
               "a caller defect: PointFinder emits only when the observed residue is in that "
               "position's Res_codon set, so e.g. gyrA_S83I (AMRFinder) is not callable because "
               "PointFinder lists only S83 A/L/V/W. Position-level restriction does not equalize the "
               "two catalogues.",
           ]}
    a.out.write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")
    print(f"\nwrote {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
