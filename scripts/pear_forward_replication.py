"""Run the SHIPPED genome-edit forward path against PEAR's measured blaCTX-M-14 fitness.

This is external replication of the one learned regime that works. The forward cell was validated on
TEM-1 + ampicillin (genome-edit path, Spearman 0.7611 over 1,715 variants). PEAR is the same SHAPE on a
DIFFERENT beta-lactamase (CTX-M-14) with DIFFERENT drugs (cefotaxime, ceftazidime), measured by an
independent lab. Nothing about the forward cell is re-fit here: `predict_genome_edit` is called exactly
as shipped.

COORDINATES, established by measurement rather than assumption. PEAR's `C648T` notation is 1-based on
the authors' own 795-nt reference (`Genotype_barcode_calling/Ref_CTXM.fasta`), NOT on the full CTX-M-14
CDS: all 2,114 variants satisfy ref[pos-1] == wt, while the +81 offset convention used by their Figure 2
axis matches only 452/2114 (chance). The reference is the MATURE protein — the full gene's 81-nt signal
peptide is trimmed — so `aa_pos` here is relative to that reference and is NOT Ambler numbering. Do not
quote a residue number from this script against the literature without converting.

THE SANITY CHECK IS THE POINT. A coordinate error produces plausible garbage rather than an error, so
the run asserts something the data must satisfy independently of whether the predictor works: NONSENSE
variants must be strongly less fit than SILENT ones. If that fails, the mapping is wrong and no
correlation from this script means anything.

Offline for blosum62. Writes wiki/pear_forward_replication_<date>.json.
"""
from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from dna_decode.forward.genome_edit import predict_genome_edit, translate_cds  # noqa: E402

REF = Path("D:/dna_decode_cache/pear/CTXM-14/Genotype_barcode_calling/Ref_CTXM.fasta")
TABLE = Path("D:/dna_decode_cache/pear/extracted/Figure3.A__data.tsv")
VARIANT = re.compile(r"^([ACGT])(\d+)([ACGT])$")


def load_ref(p: Path) -> str:
    return "".join(l.strip() for l in p.read_text().splitlines() if not l.startswith(">"))


def spearman(xs: list[float], ys: list[float]) -> float | None:
    """Rank correlation with MID-RANKS for ties (the documented tie trap)."""
    n = len(xs)
    if n < 3:
        return None

    def ranks(v):
        order = sorted(range(n), key=lambda i: v[i])
        r = [0.0] * n
        i = 0
        while i < n:
            j = i
            while j + 1 < n and v[order[j + 1]] == v[order[i]]:
                j += 1
            mid = (i + j) / 2.0 + 1.0
            for k in range(i, j + 1):
                r[order[k]] = mid
            i = j + 1
        return r

    rx, ry = ranks(xs), ranks(ys)
    mx, my = sum(rx) / n, sum(ry) / n
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    dx = sum((a - mx) ** 2 for a in rx) ** 0.5
    dy = sum((b - my) ** 2 for b in ry) ** 0.5
    return None if dx == 0 or dy == 0 else num / (dx * dy)


def median(v: list[float]) -> float:
    s = sorted(v)
    n = len(s)
    return s[n // 2] if n % 2 else (s[n // 2 - 1] + s[n // 2]) / 2.0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--method", default="blosum62", help="blosum62 | esm2 | prosst | hybrid")
    ap.add_argument("--esm-table", type=Path, default=None, help="precomputed ESM score table (JSON)")
    ap.add_argument("--prosst-table", type=Path, default=None,
                    help="ProSST variant table (JSON, 'wt{pos}alt' -> log-ratio); with --method hybrid "
                         "it is rank-averaged with the ESM table by the shipped rank_average_hybrid")
    ap.add_argument("--out", type=Path, default=None)
    a = ap.parse_args()

    if not REF.is_file() or not TABLE.is_file():
        print(f"missing inputs: {REF} / {TABLE}")
        return 2

    cds = load_ref(REF)
    protein = translate_cds(cds).rstrip("*")
    print(f"reference: {len(cds)} nt -> {len(protein)} aa (mature protein; signal peptide trimmed)")

    # Position keys MUST be ints. JSON round-trips them as strings and `predict_effect` then raises
    # KeyError per variant -- which this script records as a SKIP rather than a crash, so an uncoerced
    # table silently reduces the run to synonymous variants only and reports "missense n=0" instead of
    # failing. Coerce here; the guard below refuses to report if it still scores nothing.
    esm_table = None
    if a.esm_table:
        esm_table = {int(k): v for k, v in json.loads(a.esm_table.read_text()).items()}

    # For the hybrid: convert the ESM POSITION table to a VARIANT table so both modalities are keyed
    # the same way, then hand both to the shipped rank_average_hybrid via predict_effect.
    prosst_table = json.loads(a.prosst_table.read_text()) if a.prosst_table else None

    hybrid_tables = None
    if a.method == "hybrid":
        from dna_decode.forward.variant_effect import esm_pos_table_to_variant_table
        if not (esm_table and a.prosst_table):
            print("method=hybrid needs BOTH --esm-table and --prosst-table")
            return 2
        hybrid_tables = [esm_pos_table_to_variant_table(esm_table, protein), prosst_table]
        inter = set(hybrid_tables[0]) & set(hybrid_tables[1])
        print(f"hybrid tables: esm={len(hybrid_tables[0])} prosst={len(prosst_table)} "
              f"intersection={len(inter)}")
        if len(inter) < 1000:
            print("REFUSING: hybrid intersection too small; the tables are not keyed compatibly.")
            return 3

    rows = list(csv.DictReader(TABLE.open(encoding="utf-8"), delimiter="\t"))
    recs, skipped = [], []
    for r in rows:
        m = VARIANT.match(r["genotype"])
        if not m:
            skipped.append((r["genotype"], "unparsed"))
            continue
        wt, pos, alt = m.group(1), int(m.group(2)), m.group(3)
        try:
            p = predict_genome_edit(cds, pos, wt, alt, protein_seq=protein,
                                    protein="blaCTX-M-14 (mature)", method=a.method,
                                    esm_table=esm_table, prosst_table=prosst_table,
                                    hybrid_tables=hybrid_tables)
        except Exception as e:                      # a REF mismatch must surface, never be swallowed
            skipped.append((r["genotype"], f"{type(e).__name__}: {e}"))
            continue
        try:
            ctx, caz = float(r["CTX"]), float(r["CAZ"])
        except ValueError:
            skipped.append((r["genotype"], "unparseable measurement"))
            continue
        score = None
        if p.protein_prediction is not None:
            score = getattr(p.protein_prediction, "raw_score", None)
        recs.append({"genotype": r["genotype"], "consequence": p.consequence, "aa_mut": p.aa_mutation,
                     "score": score, "ctx": ctx, "caz": caz,
                     "abstain": bool(getattr(p.protein_prediction, "abstain", False))
                     if p.protein_prediction is not None else False})

    by = {}
    for rec in recs:
        by.setdefault(rec["consequence"], []).append(rec)
    print(f"scored {len(recs)} variants, skipped {len(skipped)}")
    for k, v in sorted(by.items()):
        print(f"  {k:9} n={len(v):5d}  median CTX={median([x['ctx'] for x in v]):+.4f}  "
              f"median CAZ={median([x['caz'] for x in v]):+.4f}")

    # ---- the independent pipeline check: nonsense must be less fit than silent -------------------
    sil, non = by.get("silent", []), by.get("nonsense", [])
    check = {"n_silent": len(sil), "n_nonsense": len(non)}
    if sil and non:
        for drug in ("ctx", "caz"):
            ms, mn = median([x[drug] for x in sil]), median([x[drug] for x in non])
            check[f"median_silent_{drug}"] = round(ms, 5)
            check[f"median_nonsense_{drug}"] = round(mn, 5)
            check[f"nonsense_below_silent_{drug}"] = bool(mn < ms)
        check["passed"] = all(check[f"nonsense_below_silent_{d}"] for d in ("ctx", "caz"))
    else:
        check["passed"] = None
        check["note"] = "not evaluable: no silent and/or nonsense variants in the table"
    print(f"\ncoordinate/pipeline sanity check: {check}")

    # ---- the actual correlation, on MISSENSE only ------------------------------------------------
    mis = [r for r in by.get("missense", []) if r["score"] is not None]
    if not mis:
        print(f"REFUSING to report: method={a.method} scored 0 missense variants "
              f"({len(skipped)} skipped, e.g. {skipped[:2]}). That is a table/key problem, not a result.")
        return 3
    corr = {}
    for drug in ("ctx", "caz"):
        corr[drug] = {"n": len(mis),
                      "spearman_score_vs_fitness": (round(spearman([r["score"] for r in mis],
                                                                   [r[drug] for r in mis]), 4)
                                                    if len(mis) >= 3 else None)}
    print(f"\nmissense n={len(mis)}")
    for drug, c in corr.items():
        print(f"  {drug.upper():4} spearman(predicted_score, measured_fitness) = {c['spearman_score_vs_fitness']}")

    out = {"schema": "pear-forward-replication-v1",
           "method": a.method,
           "reference": {"path": str(REF), "n_nt": len(cds), "n_aa": len(protein),
                         "note": "mature protein; the full gene's 81-nt signal peptide is trimmed. "
                                 "aa_pos is NOT Ambler numbering."},
           "coordinate_basis": "1-based on the authors' 795-nt reference; 2114/2114 ref-base matches "
                               "(the +81 Figure-2 axis convention matches only 452/2114)",
           "n_scored": len(recs), "n_skipped": len(skipped), "skipped_examples": skipped[:10],
           "by_consequence": {k: len(v) for k, v in sorted(by.items())},
           "sanity_check": check,
           "correlation": corr,
           "sign_convention": "predict_effect returns a score where MORE NEGATIVE = more damaging; "
                              "measured fitness is relative growth where HIGHER = fitter. A working "
                              "predictor therefore gives a POSITIVE Spearman.",
           "benchmark": {"tem1_ampicillin_genome_edit_spearman": 0.7611,
                         "note": "the forward cell's own validation, on a DIFFERENT beta-lactamase "
                                 "and a DIFFERENT drug. Not a target; a reference point."}}
    dest = a.out or (Path(__file__).resolve().parents[1] / "wiki" /
                     f"pear_forward_replication_{a.method}_2026-09-01.json")
    dest.write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")
    print(f"\nwrote {dest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
