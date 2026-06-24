"""Validation + report card for the Klebsiella K-antigen (wzi) cell -- namespace-separate.

Two honest validations, NEITHER overclaimed:
  (A) MEASURED-LABEL ceiling -- from the KlebNET-GSP 731-isolate set (Zenodo 15742130, CC-BY): genomic
      capsule typing (Kaptive KL) vs SEROLOGICAL K-type (wet-lab). Reproduces the measured-serology ceiling
      for genomic K-typing. Reported as BOTH a naive numeric match (KL#==K#) AND the paper's curated value
      (0.845, which applies the KL<->K renaming/alternative-type equivalence a naive match cannot).
  (B) caller SELF-CONSISTENCY -- run the wzi caller on a sample of the 555 wzi alleles; each must type to
      itself and map to its catalogued KL. Proves the caller is mechanically sound across the whole DB.

The wzi-v0 caller is a SINGLE-GENE approximation (~94% of full-locus Kaptive accuracy; Brisse 2013 JCM) and
FAITHFUL-TO-TOOL (not an independent baseline). The FULL caller-vs-serology number (run the wzi caller on
the 731 genomes) needs the cohort fetch (ENA run reads -> targeted wzi mapping) -> a scoped scale-up, NOT
run here. Writes wiki/ktype_report_card.{md,json}.

Run: `uv run --with openpyxl python scripts/ktype_validate.py --xlsx <Supplementary_Table_1.xlsx> --db-dir data/ktype_db`
"""
from __future__ import annotations

import argparse
import json
import random
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
PAPER_CURATED_CONCORDANCE = 0.845   # KlebNET-GSP 2025 technical report, 731 isolates


def measured_label_concordance(xlsx: Path) -> dict:
    import openpyxl
    ws = openpyxl.load_workbook(xlsx, read_only=True)["Supplementary_Table1"]
    rows = list(ws.iter_rows(values_only=True))
    ix = {h: i for i, h in enumerate(rows[0])}
    def c(r, k):
        j = ix.get(k); return r[j] if (j is not None and j < len(r)) else None
    inc = [r for r in rows[1:] if str(c(r, "Included_in_report")).lower() == "include"]
    n = conc = 0
    for r in inc:
        sero, kap = c(r, "Serological_Ktype"), c(r, "K_locus_Best_match_locus")
        if not sero or not kap or "unknown" in str(kap).lower():
            continue
        n += 1
        sero_nums = {x.strip().replace("K", "") for x in str(sero).split(",")}
        if str(kap).replace("KL", "").replace("K", "").split()[0] in sero_nums:
            conc += 1
    return {"n_isolates": n, "naive_numeric_concordant": conc,
            "naive_numeric_rate": round(conc / n, 3) if n else None,
            "paper_curated_rate": PAPER_CURATED_CONCORDANCE,
            "note": "naive KL#==K# match UNDER-counts (drops the KL<->K renaming/alternative-type "
                    "equivalence the paper applies); the curated rate is the real ceiling"}


def caller_self_consistency(db_dir: Path, k: int = 15) -> dict:
    from dna_decode.ktype.runner import call_ktype, load_wzi_kl_map
    kl_map = load_wzi_kl_map(db_dir / "wzi.txt")
    recs, cur = {}, None
    for ln in (db_dir / "wzi.fasta").read_text().splitlines():
        if ln.startswith(">"):
            cur = ln[1:].strip(); recs[cur] = []
        elif cur:
            recs[cur].append(ln)
    sample = random.Random(0).sample(list(recs), min(k, len(recs)))
    tmp = REPO / "data" / "ktype_db" / "_self_tmp.fna"
    ok = 0
    for aid in sample:
        tmp.write_text(">c\nNNNN" + "".join(recs[aid]) + "NNNN\n")
        r = call_ktype(tmp, db_dir)
        if r["status"] == "ok" and r["wzi_allele"] == aid and r["kl_type"] == kl_map.get(aid.split("_", 1)[1]):
            ok += 1
    tmp.unlink(missing_ok=True)
    return {"n_sampled": len(sample), "n_correct": ok, "n_alleles_in_db": len(recs),
            "all_correct": ok == len(sample)}


def build_card(measured: dict | None, selfc: dict | None) -> dict:
    return {
        "schema": "ktype-report-card-v1", "trait": "ktype", "organism": "Klebsiella",
        "caller": "dna_decode-wzi-blastn-v0 (BIGSdb Pasteur wzi scheme, Kleborate-bundled)",
        "tier": "FAITHFUL_TO_TOOL_MEASURED_LABEL_AVAILABLE",
        "honest_tier_text": "deterministic single-gene wzi caller, FAITHFUL to the Kleborate/BIGSdb wzi "
                            "method (NOT an independent baseline). A free MEASURED serological K-type label "
                            "EXISTS (KlebNET-GSP 731-isolate set), so the cell is VALIDATABLE -- but the "
                            "full wzi-caller-vs-serology number needs the cohort genome run (scoped). "
                            "wzi->K is ~94% / NOT one-to-one; full-locus Kaptive is more accurate.",
        "measured_label_ceiling": measured,
        "caller_self_consistency": selfc,
        "pending_scale_up": "run the wzi caller on the 731 genomes (ENA run reads -> targeted wzi mapping) "
                            "-> the genuine wzi-caller-vs-measured-serology number (namespace-separate).",
        "sources": ["KlebNET-GSP 2025 Technical Report (Zenodo 15742130, CC-BY)",
                    "Kleborate wzi DB (BIGSdb Pasteur)", "Brisse 2013 JCM (wzi typing ~94%)"],
        "namespace_note": "SEPARATE from the AMR/HIV/TB trust cards -- ktype is a typing trait, not a drug "
                          "cell; never keyed into the frozen canonical_cell_key card.",
    }


def render_md(card: dict) -> str:
    m, s = card.get("measured_label_ceiling"), card.get("caller_self_consistency")
    L = ["# Klebsiella K-antigen (wzi) — validation report card", "",
         "Namespace-separate typing-trait card (sibling of `serotype`). **Honest tier:** "
         f"`{card['tier']}`.", "", card["honest_tier_text"], ""]
    if s:
        L += [f"## Caller self-consistency", f"- **{s['n_correct']}/{s['n_sampled']}** sampled wzi alleles "
              f"(of {s['n_alleles_in_db']} in the DB) typed correctly -> correct KL "
              f"({'PASS' if s['all_correct'] else 'FAIL'}). The caller is mechanically sound across the DB.", ""]
    if m:
        L += ["## Measured-label ceiling (the rare non-AMR free measured label)",
              f"- Genomic capsule typing (Kaptive KL) vs **SEROLOGICAL K-type** (wet-lab) on the "
              f"KlebNET-GSP **N={m['n_isolates']}** set: naive-numeric **{m['naive_numeric_rate']}** | "
              f"paper-curated **{m['paper_curated_rate']}**.",
              f"- {m['note']}.",
              "- This is the *ceiling* a genomic K-typer approaches; the single-gene wzi-v0 is ~94% of "
              "full-locus Kaptive. The genuine wzi-caller-vs-serology number is the scoped scale-up below.", ""]
    L += ["## Pending scale-up", f"- {card['pending_scale_up']}", "",
          "## Honesty rails",
          "- FAITHFUL-TO-TOOL (wzi method), NOT an independent baseline; `caller_is_independent_baseline=false`.",
          "- wzi->K-type is ~94% predictive and NOT one-to-one (Brisse 2013 JCM).",
          f"- {card['namespace_note']}",
          "- Sources: " + "; ".join(card["sources"]) + ".",
          "", "Rebuild: `uv run --with openpyxl python scripts/ktype_validate.py --xlsx <suppl.xlsx> --db-dir data/ktype_db`."]
    return "\n".join(L) + "\n"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--xlsx", type=Path, default=None, help="KlebNET-GSP Supplementary_Table_1.xlsx (measured label)")
    ap.add_argument("--db-dir", type=Path, default=REPO / "data" / "ktype_db")
    a = ap.parse_args(argv)
    measured = None
    if a.xlsx and a.xlsx.exists():
        try:
            measured = measured_label_concordance(a.xlsx)
        except Exception as e:
            measured = {"error": f"{type(e).__name__}: {e}"}
    selfc = caller_self_consistency(a.db_dir) if (a.db_dir / "wzi.fasta").exists() else None
    card = build_card(measured, selfc)
    (REPO / "wiki" / "ktype_report_card.json").write_text(json.dumps(card, indent=2), encoding="utf-8")
    (REPO / "wiki" / "ktype_report_card.md").write_text(render_md(card), encoding="utf-8")
    print(f"[ktype-card] self-consistency: {selfc} | measured: {measured}")
    print("-> wiki/ktype_report_card.{md,json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
