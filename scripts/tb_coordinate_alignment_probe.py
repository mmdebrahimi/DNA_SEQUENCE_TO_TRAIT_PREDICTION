"""TB coordinate-alignment probe — does the WHO catalogue align to the CRyPTIC VCFs?

The TB-decoder technical plan (`plans/TB_AMR_Decoder_CRyPTIC_Technical_Plan.md`) flags ONE
high-severity, load-bearing assumption that must be *verified, not assumed* before any number:

    CRyPTIC per-isolate VCFs and the WHO M. tuberculosis mutation catalogue v2 (2023) both use
    H37Rv NC_000962.3 — so a catalogue determinant's genomic position/ref/alt should appear
    verbatim as a variant record in a resistant isolate's VCF.

This probe verifies that on the two canonical sentinels (RIF `rpoB S450L`, INH `katG S315T`),
using ONLY data already on disk: the pinned WHO catalogue (`data/raw/who_tb_catalogue/`, see
CHECKSUMS) + the cached CRyPTIC VCFs (`data/raw/cryptic/vcf_cache/`). NO network, NO new download.

It is an ACQUISITION/VERIFICATION probe (sibling of `scripts/cryptic_feasibility_probe.py`), NOT a
plan code module — it does not build the decoder, score a cohort, or claim sens/spec. It de-risks
Stage 1 (catalogue join) + Stage 3 (coordinate-alignment fixture) of the plan.

Run:
  uv run python scripts/tb_coordinate_alignment_probe.py
Exit: 0 = ALIGNED (both sentinels matched), 1 = MISALIGNED/MISSING, 2 = inputs absent.
"""
from __future__ import annotations

import csv
import json
import sys
from datetime import date as _date
from pathlib import Path

WHO_DIR = Path("data/raw/who_tb_catalogue")
MASTER = WHO_DIR / "WHO-UCN-TB-2023.6-eng_catalogue_master_file.txt"
COORDS = WHO_DIR / "WHO-UCN-TB-2023.7-eng_genomic_coordinates.txt"
CHECKSUMS = WHO_DIR / "CHECKSUMS"
VCF_CACHE = Path("data/raw/cryptic/vcf_cache")

# Grade-1/2 = "Associated with resistance" (the rule source). Column "INITIAL CONFIDENCE GRADING"
# in the master file encodes "1) Assoc w R" / "2) Assoc w R - Interim" / "3) Uncertain ..." etc.
GRADE_1_2_PREFIXES = ("1)", "2)")
GRADE_COL = "INITIAL CONFIDENCE GRADING"

# The two canonical sentinels: variant name (the JOIN key shared by master + coords files) + drug.
SENTINELS = [
    {"variant": "rpoB_p.Ser450Leu", "drug": "Rifampicin", "label": "RIF rpoB S450L"},
    {"variant": "katG_p.Ser315Thr", "drug": "Isoniazid", "label": "INH katG S315T"},
]


def _load_master_grades() -> dict[str, list[dict]]:
    """variant -> list of {drug, gene, tier, grade} rows from the grade-bearing master file."""
    out: dict[str, list[dict]] = {}
    with open(MASTER, encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            v = (row.get("variant") or "").strip()
            if not v:
                continue
            out.setdefault(v, []).append({
                "drug": (row.get("drug") or "").strip(),
                "gene": (row.get("gene") or "").strip(),
                "tier": (row.get("tier") or "").strip(),
                "grade": (row.get(GRADE_COL) or "").strip(),
            })
    return out


def _load_coords() -> dict[str, list[dict]]:
    """variant -> list of {chrom, pos, ref, alt} from the genomic-coordinates file."""
    out: dict[str, list[dict]] = {}
    with open(COORDS, encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            v = (row.get("variant") or "").strip()
            if not v:
                continue
            out.setdefault(v, []).append({
                "chrom": (row.get("chromosome") or "").strip(),
                "pos": int(row["position"]),
                "ref": (row.get("reference_nucleotide") or "").strip(),
                "alt": (row.get("alternative_nucleotide") or "").strip(),
            })
    return out


def _is_grade_1_2(rows: list[dict]) -> bool:
    return any(r["grade"].startswith(GRADE_1_2_PREFIXES) for r in rows)


def _scan_vcfs_for(targets: list[dict]) -> list[dict]:
    """Find PASS, non-reference (GT allele >=1) VCF records matching any (pos, ref, alt) target.

    Cached files are DECOMPRESSED VCF text (the feasibility probe gzip-decompresses before caching),
    so they are read as plaintext despite the `.vcf.gz` name.
    """
    by_pos: dict[int, list[dict]] = {}
    for t in targets:
        by_pos.setdefault(t["pos"], []).append(t)
    hits: list[dict] = []
    for vcf in sorted(VCF_CACHE.glob("*")):
        try:
            text = vcf.read_text(encoding="utf-8", errors="replace")
        except OSError:
            continue
        for line in text.splitlines():
            if line.startswith("#") or not line.strip():
                continue
            f = line.split("\t")
            if len(f) < 10:
                continue
            try:
                pos = int(f[1])
            except ValueError:
                continue
            if pos not in by_pos:
                continue
            ref, alt, filt, gt_field = f[3].strip(), f[4].strip(), f[6].strip(), f[9]
            gt = gt_field.split(":")[0]
            non_ref = any(a.isdigit() and int(a) >= 1 for a in gt.replace("|", "/").split("/"))
            for tgt in by_pos[pos]:
                if ref == tgt["ref"] and alt == tgt["alt"] and filt == "PASS" and non_ref:
                    hits.append({"isolate": vcf.name[:48], "pos": pos, "ref": ref, "alt": alt,
                                 "filter": filt, "gt": gt, "variant": tgt["variant"],
                                 "label": tgt["label"]})
    return hits


def main() -> int:
    for p in (MASTER, COORDS):
        if not p.exists():
            print(f"missing {p} — fetch the pinned WHO catalogue first (see {CHECKSUMS}).")
            return 2
    if not any(VCF_CACHE.glob("*")):
        print(f"no cached CRyPTIC VCFs under {VCF_CACHE}.")
        return 2

    master = _load_master_grades()
    coords = _load_coords()
    n_grade12 = sum(1 for v, rows in master.items() if _is_grade_1_2(rows))
    print(f"WHO catalogue: {len(master)} distinct variants, {n_grade12} grade-1/2 (Assoc w R).")
    print(f"genomic-coordinate rows: {sum(len(v) for v in coords.values())} for "
          f"{len(coords)} variants.\n")

    results = []
    all_aligned = True
    for s in SENTINELS:
        v = s["variant"]
        m_rows = master.get(v, [])
        c_rows = coords.get(v, [])
        grade12 = _is_grade_1_2(m_rows)
        targets = [{"variant": v, "label": s["label"], **c} for c in c_rows]
        hits = _scan_vcfs_for(targets)
        aligned = bool(grade12 and c_rows and hits)
        all_aligned = all_aligned and aligned
        status = "ALIGNED" if aligned else "NOT_ALIGNED"
        coord_preview = [(c["chrom"], c["pos"], c["ref"] + ">" + c["alt"]) for c in c_rows][:4]
        print(f"=== {s['label']} ({v}) -> {status} ===")
        print(f"  master grade-1/2: {grade12}  ({sorted({r['grade'] for r in m_rows})})")
        print(f"  catalogue coords: {coord_preview}")
        print(f"  VCF matches (PASS + GT>=1, ref/alt exact): {len(hits)} isolate(s)")
        for h in hits[:3]:
            print(f"    {h['isolate']}  {h['pos']} {h['ref']}>{h['alt']} {h['filter']} GT={h['gt']}")
        print()
        results.append({"sentinel": s["label"], "variant": v, "grade_1_2": grade12,
                        "n_catalogue_coords": len(c_rows), "n_vcf_matches": len(hits),
                        "aligned": aligned,
                        "example_matches": hits[:3]})

    verdict = "ALIGNED" if all_aligned else "MISALIGNED_OR_MISSING"
    pinned = ""
    if CHECKSUMS.exists():
        for ln in CHECKSUMS.read_text(encoding="utf-8").splitlines():
            if ln.startswith("# pinned_commit:"):
                pinned = ln.split(":", 1)[1].strip()

    art = {
        "_schema": "tb-coordinate-alignment-probe-v1",
        "date": _date.today().isoformat(),
        "who_catalogue_pinned_commit": pinned,
        "reference": "NC_000962.3 (H37Rv) — both CRyPTIC VCFs and WHO catalogue",
        "n_catalogue_variants": len(master),
        "n_grade_1_2": n_grade12,
        "n_cached_vcfs": len(list(VCF_CACHE.glob("*"))),
        "verdict": verdict,
        "sentinels": results,
        "honesty": ("Verifies coordinate alignment on 2 canonical sentinels against the cached VCF "
                    "subset. Proves the catalogue<->VCF reference frame matches; does NOT establish "
                    "cohort sens/spec or biological independence (WHO catalogue built partly from "
                    "CRyPTIC -> any CRyPTIC score is a KNOWLEDGE_BASELINE)."),
    }
    out_json = Path(f"wiki/tb_coordinate_alignment_probe_{_date.today().isoformat()}.json")
    out_json.write_text(json.dumps(art, indent=2), encoding="utf-8")

    md = [f"# TB coordinate-alignment probe — {verdict} ({_date.today().isoformat()})", "",
          f"Verifies the TB-decoder plan's load-bearing assumption: the WHO catalogue v2 (pinned "
          f"`{pinned[:12]}`) and the cached CRyPTIC VCFs share reference **NC_000962.3** and align "
          f"at the variant level. Data-only (no network). Sibling of the CRyPTIC feasibility probe.",
          "",
          f"- WHO catalogue: **{len(master)}** distinct variants, **{n_grade12}** grade-1/2 (Assoc w R).",
          f"- Cached CRyPTIC VCFs scanned: **{len(list(VCF_CACHE.glob('*')))}**.",
          f"- Verdict: **{verdict}**.", "",
          "## Sentinels", "",
          "| sentinel | grade-1/2 | catalogue coords | VCF matches (PASS+GT≥1, exact ref/alt) | aligned |",
          "|---|---|---|---|---|"]
    for r in results:
        ex = r["example_matches"][0] if r["example_matches"] else None
        coordstr = f"{ex['pos']} {ex['ref']}>{ex['alt']}" if ex else "—"
        md.append(f"| {r['sentinel']} ({r['variant']}) | {r['grade_1_2']} | "
                  f"{coordstr} | {r['n_vcf_matches']} | **{r['aligned']}** |")
    md += ["", "## Honesty", "", art["honesty"], ""]
    out_md = Path(f"wiki/tb_coordinate_alignment_probe_{_date.today().isoformat()}.md")
    out_md.write_text("\n".join(md), encoding="utf-8")

    print(f"verdict: {verdict}")
    print(f"artifact -> {out_md}")
    print(f"artifact -> {out_json}")
    return 0 if all_aligned else 1


if __name__ == "__main__":
    raise SystemExit(main())
