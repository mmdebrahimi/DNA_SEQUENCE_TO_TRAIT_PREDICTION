"""Growth-coupled strain design over a genome-scale model — the DESIGN direction of the FBA cell.

Searches gene knockouts that make producing a target metabolite NECESSARY for growth, and emits a
`.md` + `.json` artifact pair. Every design is a HYPOTHESIS FOR THE BENCH, never a validated strain.

    uv run python scripts/fba_strain_design.py --target succ
    uv run python scripts/fba_strain_design.py --target EX_lac__D_e --organism ecoli --max-knockouts 2
"""
from __future__ import annotations

import argparse
import json
import sys
from datetime import date
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from dna_decode.fba.design import find_coupled_designs  # noqa: E402
from dna_decode.fba.model import load_model, organism_for, resolve_model_id, wildtype_growth  # noqa: E402


def render_md(rec: dict) -> str:
    L = [
        f"# Growth-coupled strain design — {rec['target_name']} ({rec['organism']}, {rec['date']})",
        "",
        f"Condition: **{rec['condition']}**.",
        "",
        f"Model **{rec['model']}** ({rec['organism']}) · target **`{rec['target_reaction']}`** · "
        f"growth floor **{rec['growth_floor_frac']:.0%} of each strain's OWN max growth** "
        f"(wild type max {rec['wildtype_growth_per_h']} /h) · knockout level **{rec['knockout_level']}**.",
        "",
        "**Growth coupling** means every flux distribution that grows at the floor also secretes the "
        "product — production is obligatory, not merely allowed. That is what makes an engineered strain "
        "stable: selection for growth becomes selection for production.",
        "",
        "## Baseline (wild type)",
        "",
        f"- min product flux **{rec['baseline']['min_product_flux']}** · max **{rec['baseline']['max_product_flux']}** "
        f"→ **{rec['baseline']['coupling']}**",
        "",
    ]
    if rec.get("wildtype_already_coupled"):
        L += [
            f"> **The wild type is ALREADY growth-coupled here** (guaranteed "
            f"{rec['baseline']['min_product_flux']}), so `OBLIGATORY` alone is not evidence of a design. "
            f"Every design below must **beat** that floor; **{rec.get('n_inheriting_baseline_only', 0)}** "
            f"knockout sets were coupled only by inheriting it and are NOT counted.",
            "",
        ]
    L += [
        "## Search",
        "",
        f"- {rec['search']}",
        f"- {rec['knockout_level']} candidates scanned **{rec['n_candidates_scanned']}**; non-viable at the growth floor **{rec['n_nonviable_at_floor']}**; "
        f"singles evaluated **{rec['n_singles_evaluated']}**; pairs **{rec['n_pairs_evaluated']}**; "
        f"triples **{rec.get('n_triples_evaluated', 0)}**",
        "",
        f"## Result — {rec['n_coupled_designs']} growth-coupled design(s)",
        "",
    ]
    if rec["designs"]:
        L += [
            "| knockouts | guaranteed product flux | gain over wild type | max product flux | growth (/h) |",
            "|---|---|---|---|---|",
        ]
        for d in rec["designs"]:
            ko = ", ".join(f"`{k}`" for k in d["knockouts"])
            L.append(f"| {ko} | **{d['min_product_flux']}** | +{d.get('improvement_over_wildtype', 0)} | "
                     f"{d['max_product_flux']} | {d['growth_per_h']} |")
        L.append("")
    else:
        L += [
            "**No growth-coupled design found** under this search. The closest candidates (still "
            "uncoupled — the cell can avoid producing) were:",
            "",
            "| knockouts | min product | max product | growth (/h) | coupling |",
            "|---|---|---|---|---|",
        ]
        for d in rec["best_uncoupled"]:
            ko = ", ".join(f"`{k}`" for k in d["knockouts"]) or "(none)"
            L.append(
                f"| {ko} | {d['min_product_flux']} | {d['max_product_flux']} | {d['growth_per_h']} | {d['coupling']} |"
            )
        L += ["", "A negative result here is informative: it bounds what single/paired knockouts can do "
              "for this target on this medium, and points at medium or pathway changes instead.", ""]
    L += ["## Scope", "", f"{rec['scope']}", ""]
    return "\n".join(L)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--target", required=True, help="product: reaction id (EX_succ_e) or short name (succ)")
    ap.add_argument("--organism", default=None, help="ecoli(default) | saureus | salmonella | pputida | yeast")
    ap.add_argument("--growth-frac", type=float, default=0.9,
                    help="growth floor as a fraction of each strain's OWN max growth (default 0.9 -- coupling is defined near-optimal; a low floor leaves slack and finds nothing)")
    ap.add_argument("--max-knockouts", type=int, default=3, choices=[1, 2, 3])
    ap.add_argument("--level", default="reaction", choices=["reaction", "gene"],
                    help="knockout level (default reaction -- GPR isozymes make gene-level unable to disable most reactions)")
    ap.add_argument("--candidates", default=None,
                    help="comma-separated knockout candidates to scope the search to (e.g. a pathway: PFL,LDH_D,ALCD2x,ACKr,PTAr). The unrestricted heuristic pool can MISS a design whose members are individually unremarkable.")
    ap.add_argument("--slug-suffix", default="", help="extra tag for the artifact filename")
    ap.add_argument("--pair-pool", type=int, default=40, help="top-N singles to pair (bounded heuristic)")
    ap.add_argument("--anaerobic", action="store_true",
                    help="close O2 uptake; redox must balance through secreted products (where the "
                         "classic growth-coupled designs live)")
    ap.add_argument("--date", default=str(date.today()))
    ap.add_argument("--out-dir", default=None)
    a = ap.parse_args(argv)

    model_id = resolve_model_id(a.organism)
    model = load_model(organism=a.organism)
    wt = wildtype_growth(model)
    print(f"{model_id} ({organism_for(model_id)}) WT {wt:.4f} /h — searching designs for '{a.target}' ...")

    rec = find_coupled_designs(
        model, a.target, growth_frac=a.growth_frac, max_knockouts=a.max_knockouts,
        wt_growth=wt, pair_pool=a.pair_pool, anaerobic=a.anaerobic, level=a.level,
        gene_ids=[c.strip() for c in a.candidates.split(',')] if a.candidates else None,
    )
    rec.update({
        "record": "fba-strain-design-v1",
        "date": a.date,
        "model": model_id,
        "organism": organism_for(model_id),
        "method": "growth-coupled knockout search (two-sided LP at a fixed growth floor); mechanistic, deterministic",
    })

    outdir = Path(a.out_dir) if a.out_dir else Path(__file__).resolve().parent.parent / "wiki"
    outdir.mkdir(parents=True, exist_ok=True)
    slug = rec["target_reaction"].replace("EX_", "").rstrip("_e") or "target"
    # condition is part of the identity -- an aerobic and an anaerobic run must NEVER overwrite each other
    slug += ("_anaerobic" if a.anaerobic else "_aerobic") + f"_{a.level}" + (f"_{a.slug_suffix}" if a.slug_suffix else "")
    stem = outdir / f"fba_strain_design_{slug}_{a.date}"
    stem.with_suffix(".json").write_text(json.dumps(rec, indent=2), encoding="utf-8")
    stem.with_suffix(".md").write_text(render_md(rec), encoding="utf-8")

    print(f"coupled designs: {rec['n_coupled_designs']}  (singles {rec['n_singles_evaluated']}, "
          f"pairs {rec['n_pairs_evaluated']}, non-viable {rec['n_nonviable_at_floor']})")
    for d in rec["designs"][:5]:
        print(f"  KO {d['knockouts']} -> guaranteed {d['min_product_flux']} , growth {d['growth_per_h']} /h")
    print(f"wrote {stem}.md / .json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
