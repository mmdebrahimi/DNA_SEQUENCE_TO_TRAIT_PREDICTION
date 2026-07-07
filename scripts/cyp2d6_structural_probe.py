"""CYP2D6 structural surface — read-depth copy-number VALIDATION on real 1000G CRAMs.

Validates the `dna_decode.pgx.cyp2d6_structural` copy-number caller against the GeT-RM structural truth,
using CYP2D6-body / single-copy-control DEPTH RATIOS measured off the 1000G 30x CRAMs (per-sample, via
samtools over remote CRAM with ENA reference auto-fetch — no full-reference download).

Two modes:
  * default: read the COMMITTED ratios TSV (data/pgx_1000g/cyp2d6_structural_ratios.tsv) -> validate +
    emit wiki/cyp2d6_structural_<date>.{md,json}. Offline, reproducible.
  * --compute <cohort.tsv>: (live, needs Docker) compute the ratios first via samtools, write the TSV.

The headline is CN-class concordance (deletion / normal / duplication) on the NON-hybrid truths; hybrid-
containing truths (*13/*36/*68) are MEASURED but EXCLUDED from CN scoring (their copy contribution isn't
determinable from the star label — this surface never resolves hybrid identity).
"""
from __future__ import annotations

import argparse
import csv
import datetime
import json
import statistics
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
from dna_decode.pgx.cyp2d6_structural import (  # noqa: E402
    CONTROL_REGION,
    CYP2D6_REGION,
    NORMAL_BASELINE,
    structural_call,
    truth_copy_number_class,
)

RATIOS_TSV = REPO / "tests" / "data" / "pgx_getrm" / "cyp2d6_structural_ratios.tsv"
_SAMTOOLS_IMAGE = "quay.io/biocontainers/samtools:1.21--h50ea8bc_0"
_ENA_REF = "https://www.ebi.ac.uk/ena/cram/md5/%s"


def _region(r):
    return f"{r[0]}:{r[1]}-{r[2]}"


def compute_ratios(cohort_tsv: Path) -> list[dict]:
    """LIVE: compute CYP2D6/control depth ratios for each (sample, truth, cram_url) via Docker samtools."""
    from tools.docker_runner import run as docker_run
    rows = []
    for rec in csv.DictReader(cohort_tsv.open(encoding="utf-8"), delimiter="\t"):
        cram = rec["cram_url"].strip()
        script = (
            f'd6=$(samtools depth -a -r {_region(CYP2D6_REGION)} "{cram}" 2>/dev/null | '
            'awk "{s+=\\$3;n++} END{if(n)printf \\"%.2f\\",s/n; else print 0}"); '
            f'ctl=$(samtools depth -a -r {_region(CONTROL_REGION)} "{cram}" 2>/dev/null | '
            'awk "{s+=\\$3;n++} END{if(n)printf \\"%.2f\\",s/n; else print 0}"); '
            'echo "$d6 $ctl"')
        out = docker_run(_SAMTOOLS_IMAGE, ["bash", "-c", script], env={"REF_PATH": _ENA_REF},
                         capture_output=True, check=False, timeout=300)
        d6, ctl = (out.stdout or "0 0").strip().split()[:2]
        ratio = round(float(d6) / float(ctl), 3) if float(ctl) > 0 else None
        rows.append({"sample": rec["sample"], "truth": rec["truth"], "d6": d6, "ctl": ctl, "ratio": ratio})
    return rows


def validate(rows: list[dict]) -> dict:
    parsed = [r for r in rows if r.get("ratio") not in (None, "", "NA")]
    normals = [float(r["ratio"]) for r in parsed if truth_copy_number_class(r["truth"]) == "normal_copy_number"]
    baseline_measured = round(statistics.median(normals), 3) if normals else None

    scored, correct = [], 0
    measured_only = []
    dels = {"deletion", "homozygous_deletion"}
    for r in parsed:
        tclass = truth_copy_number_class(r["truth"])
        call = structural_call(float(r["ratio"]))
        rec = {"sample": r["sample"], "truth": r["truth"], "ratio": float(r["ratio"]),
               "copy_number": call.copy_number, "predicted": call.status, "truth_class": tclass}
        if tclass is None:                          # hybrid / ambiguous -> measured, not CN-scored
            measured_only.append(rec)
            continue
        ok = (call.status == tclass) or (call.status in dels and tclass in dels)
        rec["match"] = ok
        correct += ok
        scored.append(rec)

    return {
        "schema": "cyp2d6-structural-cn-v0",
        "analysis_date": datetime.date.today().isoformat(),
        "method": ("CYP2D6-body/control read-depth ratio off 1000G 30x CRAMs (samtools + ENA reference "
                   "auto-fetch, no full-reference download); ratio -> integer copy number -> deletion/normal/"
                   "duplication."),
        "cyp2d6_region": _region(CYP2D6_REGION), "control_region": _region(CONTROL_REGION),
        "normal_baseline_pinned": NORMAL_BASELINE, "normal_baseline_measured": baseline_measured,
        "n_samples_measured": len(parsed), "n_cn_scored": len(scored),
        "n_hybrid_or_ambiguous_excluded": len(measured_only),
        "cn_class_concordance": round(correct / len(scored), 4) if scored else None,
        "cn_class_hits": f"{correct}/{len(scored)}",
        "honesty_tier": ("Real-CRAM read-depth copy-number caller (the structural surface the SNP cell is "
                         "blind to). Resolves *5 deletion + *xN duplication; NEVER resolves hybrid IDENTITY "
                         "(*13/*36/*68 -> excluded from CN scoring, measured only; needs CYP2D6-vs-CYP2D7 PSV "
                         "analysis, Cyrius-class). Coarse integer CN, not a breakpoint caller."),
        "scored_rows": sorted(scored, key=lambda x: x["ratio"]),
        "measured_only_rows": sorted(measured_only, key=lambda x: x["ratio"]),
    }


def _write_reports(rep: dict) -> None:
    stem = f"cyp2d6_structural_{rep['analysis_date']}"
    (REPO / "wiki" / f"{stem}.json").write_text(json.dumps(rep, indent=2), encoding="utf-8")
    L = [f"# CYP2D6 STRUCTURAL surface — read-depth copy-number on real 1000G CRAMs ({rep['analysis_date']})",
         "", f"_{rep['method']}_", "",
         f"- Regions: CYP2D6 `{rep['cyp2d6_region']}` / control `{rep['control_region']}`",
         f"- NORMAL baseline: pinned **{rep['normal_baseline_pinned']}** vs measured "
         f"**{rep['normal_baseline_measured']}** (median of NORMAL-truth ratios)",
         f"- Samples measured: **{rep['n_samples_measured']}**  "
         f"(CN-scored **{rep['n_cn_scored']}**; hybrid/ambiguous excluded **{rep['n_hybrid_or_ambiguous_excluded']}**)",
         f"- **CN-class concordance (deletion/normal/duplication): {rep['cn_class_hits']} "
         f"({rep['cn_class_concordance']})**", "", f"_{rep['honesty_tier']}_", "",
         "## CN-scored samples (non-hybrid truth)", "",
         "| sample | truth | ratio | CN | predicted | truth class | match |", "|---|---|---|---|---|---|---|"]
    for r in rep["scored_rows"]:
        L.append(f"| {r['sample']} | `{r['truth']}` | {r['ratio']:.2f} | {r['copy_number']} | "
                 f"{r['predicted']} | {r['truth_class']} | {'OK' if r['match'] else 'X'} |")
    L += ["", "## Measured-only (hybrid / ambiguous — NOT CN-scored; identity unresolved)", "",
          "| sample | truth | ratio | CN | predicted |", "|---|---|---|---|---|"]
    for r in rep["measured_only_rows"]:
        L.append(f"| {r['sample']} | `{r['truth']}` | {r['ratio']:.2f} | {r['copy_number']} | {r['predicted']} |")
    L.append("")
    (REPO / "wiki" / f"{stem}.md").write_text("\n".join(L), encoding="utf-8")
    print("\n".join(L[:12]))
    print(f"[report -> wiki/{stem}.{{md,json}}]")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="CYP2D6 structural copy-number validation on real 1000G CRAMs.")
    ap.add_argument("--compute", type=Path, default=None,
                    help="cohort TSV (sample/truth/cram_url) -> compute ratios live via Docker, write the TSV")
    ap.add_argument("--ratios", type=Path, default=RATIOS_TSV, help="committed ratios TSV (default)")
    args = ap.parse_args(argv)

    if args.compute:
        rows = compute_ratios(args.compute)
        with args.ratios.open("w", encoding="utf-8", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=["sample", "truth", "d6", "ctl", "ratio"], delimiter="\t")
            w.writeheader()
            for r in rows:
                w.writerow(r)
        print(f"[wrote {len(rows)} ratios -> {args.ratios}]")
    if not args.ratios.exists():
        print(f"ERROR: no ratios TSV at {args.ratios} (run with --compute <cohort.tsv> first)")
        return 2
    rows = list(csv.DictReader(args.ratios.open(encoding="utf-8"), delimiter="\t"))
    rep = validate(rows)
    _write_reports(rep)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
