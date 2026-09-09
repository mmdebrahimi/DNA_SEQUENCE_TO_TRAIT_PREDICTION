"""Audit the Salmonella antigen DB itself — the O-axis gap is DB CONTENT, not code.

WHY THIS EXISTS. Scoring the serovar caller against pinned SeqSero2 1.3.2 (2026-09-08) showed the whole
gap is one axis: H1 agrees 200/200, H2 146/146, O only 159/200. Chasing that down found the cause is
neither the selection rule nor the coverage threshold -- BOTH were tested and exonerated -- but two
defects in the reference DB:

  1. MISSING ANTIGEN. There is no plain `O__9__*` allele. O group 9 (D1: Enteritidis, Typhi, Dublin) is
     among the most prevalent in Salmonella, and the reference tool calls it on 17 of our 200 cohort
     genomes. Our caller cannot emit it at all, so those genomes are forced onto the nearest available
     allele, `9,46`. This is structural: no threshold change can fix a missing reference.

  2. NON-DISCRIMINATING ALLELE. `O__1,3,19__126` is 130 bp against an O-allele median of 1155 -- the
     SHORTEST in the DB -- and it aligns at 99.2-100% coverage on 60% of a diverse cohort. A 130 bp
     sequence present in most Salmonella is not typing an antigen. Our caller emits `1,3,19` on 14
     genomes where the reference says 4.

WHAT THIS SCRIPT DOES *NOT* DO. It does not fix either defect. Both fixes need real, sourced reference
sequences -- a genuine O9 wzx/wzy, a full-length 1,3,19 -- and writing biological reference data from
memory is the fabrication hazard this project guards against elsewhere. This audit MEASURES, so the
defects are visible, sized, and cannot silently regress.

Offline: reads the committed DB, plus the cached permissive blastn sweep when present for promiscuity.
"""
from __future__ import annotations

import argparse
import json
import statistics
import sys
from collections import Counter
from datetime import date as _date
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from dna_decode.salmserovar.runner import parse_axis_antigen  # noqa: E402

SHORT_ALLELE_BP = 300          # an O allele this short cannot carry antigen-specific signal
PROMISCUITY_OVER_PREVALENCE = 5.0   # hit-rate this many times the antigen's TRUE prevalence
PROMISCUITY_FLOOR = 0.10            # ...and hitting at least this fraction, so rare antigens
                                    # cannot trip the ratio on one or two hits


def load_fasta(path: Path) -> dict[str, str]:
    seqs: dict[str, str] = {}
    name, buf = None, []
    for ln in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if ln.startswith(">"):
            if name:
                seqs[name] = "".join(buf)
            name, buf = ln[1:].strip(), []
        else:
            buf.append(ln.strip())
    if name:
        seqs[name] = "".join(buf)
    return seqs


def axis_of(allele_id: str) -> str:
    return (parse_axis_antigen(allele_id) or ("", ""))[0]


def antigen_of(allele_id: str) -> str:
    return (parse_axis_antigen(allele_id) or ("", ""))[1]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--db", type=Path,
                    default=ROOT / "data" / "salmserovar_db" / "salmonella_antigens.fasta")
    ap.add_argument("--sweep", type=Path,
                    default=Path("D:/dna_decode_cache/salm_asm/sweep_hits.jsonl"),
                    help="cached permissive blastn pass, for the promiscuity check")
    ap.add_argument("--comparison", type=Path,
                    default=ROOT / "wiki" / "salmserovar_seqsero2_2026-09-08.json")
    ap.add_argument("--out", type=Path, default=ROOT / "wiki" /
                    f"salmserovar_db_audit_{_date.today().isoformat()}.json")
    a = ap.parse_args()

    if not a.db.exists():
        print(f"antigen DB absent: {a.db}", file=sys.stderr)
        return 2
    seqs = load_fasta(a.db)
    o_len = {k: len(v) for k, v in seqs.items() if axis_of(k) == "O"}
    if not o_len:
        print("REFUSING: no O-axis alleles parsed from the DB -- the header shape changed.",
              file=sys.stderr)
        return 3
    median = statistics.median(o_len.values())
    antigens = sorted({antigen_of(k) for k in o_len})

    short = sorted(((k, n) for k, n in o_len.items() if n < SHORT_ALLELE_BP), key=lambda x: x[1])

    # --- promiscuity (needs the cached sweep) -----------------------------------------------------
    # Reference-tool O-antigen prevalence on our own cohort: the baseline every promiscuity claim is
    # measured against. Without it a flat bar mistakes a common antigen for a broken allele.
    ref_prevalence: dict[str, float] = {}
    _ck = Path("D:/dna_decode_cache/ss2_checkpoint/rows.jsonl")
    if _ck.exists():
        _ref: Counter = Counter()
        for _ln in _ck.read_text(encoding="utf-8").splitlines():
            if _ln.strip():
                _o = ((json.loads(_ln).get("seqsero2") or {}).get("O") or "").strip()
                if _o and _o != "-":
                    _ref[_o] += 1
        _tot = sum(_ref.values()) or 1
        ref_prevalence = {k: v / _tot for k, v in _ref.items()}

    promiscuous: list[dict] = []
    n_genomes = 0
    if a.sweep.exists():
        rows = [json.loads(ln) for ln in a.sweep.read_text(encoding="utf-8").splitlines() if ln.strip()]
        n_genomes = len(rows)
        hits: Counter = Counter()
        for r in rows:
            for allele, h in (r.get("per_allele") or {}).items():
                if axis_of(allele) != "O":
                    continue
                if (h.get("percent_identity") or 0) >= 90 and (h.get("percent_coverage") or 0) >= 40:
                    hits[allele] += 1
        # A flat hit-rate bar CANNOT separate promiscuity from prevalence: the allele for a genuinely
        # common antigen SHOULD hit at that antigen's rate. `O__4__231` hits 32% and O=4 really is 32%
        # of this cohort -- flagging it was the CHECK being wrong, not the allele. Compare each allele's
        # hit-rate against its own antigen's reference prevalence instead.
        for allele, c in hits.most_common():
            if not n_genomes:
                continue
            frac = c / n_genomes
            prev = ref_prevalence.get(antigen_of(allele))
            if prev is None:            # the reference never called this antigen here -- no baseline
                continue
            ratio = frac / prev if prev > 0 else float("inf")
            if frac >= PROMISCUITY_FLOOR and ratio > PROMISCUITY_OVER_PREVALENCE:
                promiscuous.append({"allele": allele, "hits": c, "of": n_genomes,
                                    "fraction": frac, "antigen_true_prevalence": prev,
                                    "times_over_prevalence": ratio,
                                    "length_bp": o_len.get(allele)})

    # --- missing antigens the REFERENCE TOOL actually called ---------------------------------------
    # Data-driven, not curated from memory: an antigen the reference emitted on our own cohort, that
    # our DB has no allele for, is a hole we can prove rather than assert.
    missing: list[dict] = []
    if a.comparison.exists():
        ck = Path("D:/dna_decode_cache/ss2_checkpoint/rows.jsonl")
        if ck.exists():
            ref = Counter()
            for ln in ck.read_text(encoding="utf-8").splitlines():
                if not ln.strip():
                    continue
                r = json.loads(ln)
                o = ((r.get("seqsero2") or {}).get("O") or "").strip()
                if o and o != "-":
                    ref[o] += 1
            for o, c in ref.most_common():
                if o not in antigens:
                    missing.append({"antigen": o, "reference_called_it_on": c,
                                    "of_genomes": sum(ref.values())})

    print(f"=== {a.db.name}: {len(o_len)} O-axis alleles, {len(antigens)} distinct antigens ===")
    print(f"  length median {median:.0f} bp (range {min(o_len.values())}-{max(o_len.values())})")
    print(f"\n  ANOMALOUSLY SHORT (<{SHORT_ALLELE_BP} bp): {len(short)}")
    for k, n in short:
        print(f"     {k:<22s} {n:4d} bp")
    if promiscuous:
        print(f"\n  NON-DISCRIMINATING (hit-rate >{PROMISCUITY_OVER_PREVALENCE:.0f}x its antigen's "
              f"true prevalence, on {n_genomes} genomes): "
              f"{len(promiscuous)}")
        for p in promiscuous:
            print(f"     {p['allele']:<22s} {p['hits']:3d}/{p['of']} = {p['fraction']:.0%} vs true "
                  f"prevalence {p['antigen_true_prevalence']:.0%} "
                  f"= {p['times_over_prevalence']:.0f}x over  ({p['length_bp']} bp)")
    else:
        print(f"\n  promiscuity check: {'no cached sweep' if not a.sweep.exists() else 'none flagged'}")
    if missing:
        print(f"\n  MISSING ANTIGENS the reference tool called on our own cohort: {len(missing)}")
        for m in missing:
            print(f"     O={m['antigen']:<8s} reference called it on {m['reference_called_it_on']} "
                  f"of {m['of_genomes']} genomes -- our DB has NO allele for it")

    findings = len(short) + len(promiscuous) + len(missing)
    verdict = "DB_CONTENT_DEFECTS_FOUND" if findings else "DB_CLEAN_ON_THESE_CHECKS"
    why = ((f"{len(missing)} antigen(s) the reference tool calls on our own cohort have NO allele in "
            f"this DB, {len(promiscuous)} allele(s) hit far more often than their antigen's true "
            f"prevalence in the same cohort, and "
            f"{len(short)} are under {SHORT_ALLELE_BP} bp against a {median:.0f} bp "
            "median. These are REFERENCE-DATA defects: no threshold or selection change can fix a "
            "missing allele, and a 130 bp sequence present in most Salmonella is not typing an antigen.")
           if findings else
           "no missing antigen, no promiscuous allele, and no anomalously short allele on these checks.")
    print(f"\nVERDICT: {verdict}\n  {why}")

    out = {"schema": "salmserovar-db-audit-v1", "date": _date.today().isoformat(),
           "db": str(a.db), "n_o_alleles": len(o_len), "n_distinct_o_antigens": len(antigens),
           "o_length_median_bp": median,
           "o_length_range_bp": [min(o_len.values()), max(o_len.values())],
           "thresholds": {"short_allele_bp": SHORT_ALLELE_BP,
                          "promiscuity_over_prevalence": PROMISCUITY_OVER_PREVALENCE,
                          "promiscuity_floor": PROMISCUITY_FLOOR},
           "short_alleles": [{"allele": k, "length_bp": n} for k, n in short],
           "promiscuous_alleles": promiscuous,
           "missing_antigens": missing,
           "distinct_o_antigens": antigens,
           "verdict": verdict, "why": why,
           "what_this_does_not_do": ("It does not fix either defect. Adding a genuine O9 wzx/wzy "
                                     "reference, or replacing the 130 bp 1,3,19 fragment with a "
                                     "full-length allele, requires SOURCED sequences; writing "
                                     "biological reference data from memory is a fabrication hazard."),
           "honest_limits": [
               "The missing-antigen check is data-driven, not curated: it lists only antigens the "
               "REFERENCE TOOL actually emitted on our own 200-genome cohort. An antigen neither the "
               "cohort nor the tool exercised would not be flagged, so absence of a flag is not proof "
               "of completeness.",
               "The promiscuity check needs the cached permissive blastn sweep; without it that "
               "section is skipped rather than reported as clean.",
               "Promiscuity is measured against each antigen's OWN reference prevalence, not a flat "
               "hit-rate: a first version used a flat 25% bar and flagged O__4__231, whose 32% hit "
               "rate exactly matches O=4 being 32% of this cohort -- the check was wrong, not the "
               "allele. The 300 bp / 5x / 10% bars are hygiene tripwires separating clear outliers "
               "from the bulk, NOT derived constants; they flag for inspection, never adjudicate.",
               "Removing the flagged 1,3,19 allele is NOT proposed here: it is the DB's ONLY 1,3,19 "
               "entry, so deleting it would trade ~10 false positives for the loss of all E4 calling, "
               "which the reference says is correct on 4 cohort genomes.",
           ]}
    a.out.write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")
    print(f"\nwrote {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
