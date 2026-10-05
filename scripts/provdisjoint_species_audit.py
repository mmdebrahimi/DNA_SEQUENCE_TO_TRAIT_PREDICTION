"""What SPECIES is each provenance-disjoint validation cohort actually made of? (read-only, offline)

    uv run python scripts/provdisjoint_species_audit.py
    uv run python scripts/provdisjoint_species_audit.py --cohort klebsiella_provdisjoint_meropenem

WHY. The `Klebsiella x meropenem` SCORED cell reports acc 0.683 / sens 0.467 with 16 false negatives whose
only recorded explanation -- the frozen rule being "blind to porin-loss-mediated R", stated in
`amr_rules.DRUG_RULE['meropenem']['validated']` -- is measurably FALSE: AMRFinder files ompK35/ompK36
truncations under Subclass=CARBAPENEM and the rule's `subclass_any={'CARBAPENEM'}` therefore counts them.
So the cell's dominant error mode has no explanation, and this script measures what those isolates ARE.

RECONCILE BEFORE ATTRIBUTE. Nothing is decomposed until the cohort's committed confusion matrix is
reproduced EXACTLY from the cached AMRFinder runs through the deployed `call_resistance`. That gate is
`compute_lineage_metrics.reconcile_raw_metrics` (reused, not rebuilt -- it already raises rather than
returns on disagreement and hands back the per-accession predictions a cross-tab needs). An attribution
built on an unreconciled prediction set is how a plausible wrong story gets published.

A COHORT THAT CANNOT RECONCILE STILL GETS A COMPOSITION, AND EXPLICITLY NO CROSS-TAB. Three Klebsiella
cohorts have incomplete cached runs (gentamicin 3/60, tetracycline 33/60, ceftriaxone 54/60). Attributing
outcomes on a partial run set would be the silent-denominator failure again, so `outcome_crosstab` is
`null` with a stated reason. Species composition needs only the assemblies, which ARE complete everywhere.

IT DISCLOSES, IT DOES NOT RE-SCORE. The 10 `provenance_disjoint_validation_*.json` files are named frozen
units of the reproducibility freeze. This script never writes them. Exit 0 always -- a report, not a gate.

THE LABEL IS THE SUBMITTER'S. Species come from each assembly's own GenBank ORGANISM line: strong evidence
about what a cohort CONTAINS, never an independent wet-lab species call.
"""
from __future__ import annotations

import argparse
import json
import sys
from datetime import date
from pathlib import Path

from dna_decode.eval.cohort_species import (
    MAX_UNRESOLVED_FRACTION,
    classify_species,
    compose,
    composition_verdict,
    resolve_cohort_species,
    vocabulary_status,
)
from scripts.compute_lineage_metrics import (
    ReconcileMismatch,
    find_artifact,
    parse_cohort_dir,
    read_selected,
    reconcile_raw_metrics,
)

SCHEMA = "provdisjoint-species-audit-v1"
RAW = Path("data/raw")
WIKI = Path("wiki")


class AuditError(RuntimeError):
    """A cohort could not be audited for a reason the operator must see."""


def provdisjoint_cohorts(raw: Path = RAW) -> list[Path]:
    """Every `*_provdisjoint_*` cohort directory that has both labels and assemblies."""
    return sorted(p for p in raw.glob("*provdisjoint*")
                  if p.is_dir() and (p / "selected.tsv").exists() and (p / "refseq").is_dir())


def runs_completeness(cohort_dir: Path, accessions: list[str]) -> tuple[int, int]:
    """(n_with_cached_main_tsv, n_total) for this cohort's OWN amrfinder_runs."""
    runs = cohort_dir / "amrfinder_runs"
    have = sum(1 for a in accessions if (runs / a / "main.tsv").exists())
    return have, len(accessions)


def audit_cohort(cohort_dir: Path, wiki: Path = WIKI) -> dict:
    """Species composition (+ an outcome cross-tab when the cohort reconciles) for one cohort."""
    slug, drug = parse_cohort_dir(cohort_dir.name)
    art_path = find_artifact(slug, drug, wiki=wiki)
    if art_path is None:
        raise AuditError(f"no committed provdisjoint artifact for {cohort_dir.name}")
    art = json.loads(art_path.read_text(encoding="utf-8"))

    # Read the scoring identity from the ARTIFACT, never from the directory name -- the directory name
    # being untrustworthy at species level is the finding under test.
    organism = art["organism"]
    reg_org = art.get("registry_organism") or organism
    amrfinder_org = art.get("amrfinder_organism") or organism
    a_drug = art["drug"]
    metrics = art.get("metrics", {})

    selected = read_selected(cohort_dir / "selected.tsv")
    accs = sorted(selected)
    resolved = resolve_cohort_species(cohort_dir / "refseq", accs)
    cls = classify_species(resolved, amrfinder_org)
    verdict = composition_verdict(cls)

    row: dict = {
        "cohort_dir": cohort_dir.name,
        "artifact": art_path.name,
        "organism": organism,
        "registry_organism": reg_org,
        "amrfinder_organism": amrfinder_org,
        "drug": a_drug,
        "committed_metrics": {k: metrics.get(k) for k in
                              ("n_scored", "tp", "fp", "tn", "fn", "abstain", "acc", "sens", "spec")},
        "composition": verdict,
        "species_vocabulary_status": {sp: vocabulary_status(sp) for sp in sorted(cls.species_counts)},
        "max_unresolved_fraction": MAX_UNRESOLVED_FRACTION,
    }

    have, total = runs_completeness(cohort_dir, accs)
    row["cached_amrfinder_runs"] = {"have": have, "total": total}

    if have < total:
        row["reconciled"] = False
        row["outcome_crosstab"] = None
        row["crosstab_withheld_reason"] = (
            f"cached AMRFinder runs incomplete ({have}/{total}); attributing outcomes on a partial run "
            f"set would compute a rate over a silently-shrunken denominator. Composition above is "
            f"unaffected -- it needs only the assemblies, which are complete."
        )
        return row

    reuse_glob = f"data/raw/{slug}_*/amrfinder_runs"
    try:
        raw_conf, preds = reconcile_raw_metrics(selected, cohort_dir / "amrfinder_runs", reuse_glob,
                                               a_drug, reg_org, metrics)
    except ReconcileMismatch as e:
        row["reconciled"] = False
        row["outcome_crosstab"] = None
        row["crosstab_withheld_reason"] = f"RECONCILE FAILED -- {e}"
        return row

    row["reconciled"] = True
    row["reconciled_confusion"] = {k: raw_conf.get(k) for k in ("tp", "fp", "tn", "fn", "n_scored")}
    row["outcome_crosstab"] = compose(resolved, preds, selected, confusion=raw_conf)

    # The headline the whole audit exists for: of the misses, how many are off-species?
    fn_accs = [a for a, p in preds.items() if p == "S" and selected.get(a) == 1]
    off = set(verdict.get("off_species_accessions", ()))
    row["false_negative_anatomy"] = {
        "n_false_negatives": len(fn_accs),
        "n_fn_off_species": sum(1 for a in fn_accs if a in off),
        "n_fn_expected_species": sum(1 for a in fn_accs if a not in off),
        "fn_species_counts": _counts([resolved.get(a) for a in fn_accs]),
    }
    return row


def _counts(names: list[str | None]) -> dict[str, int]:
    from collections import Counter

    from dna_decode.data.genbank_organism import species_of
    c: Counter[str] = Counter()
    for n in names:
        c[(species_of(n) if n else None) or "(unresolved)"] += 1
    return dict(c)


def render_md(rows: list[dict], today: str) -> str:
    L = [f"# Species composition of the provenance-disjoint validation cohorts ({today})", "",
         "Read-only audit. **It discloses; it does not re-score.** The 10",
         "`provenance_disjoint_validation_*.json` artifacts are frozen units of the reproducibility",
         "freeze and are untouched here.", "",
         "**Species come from each assembly's own GenBank `ORGANISM` line** — the submitter's assertion,",
         "not a wet-lab identification. Strong evidence about what a cohort *contains*; not an",
         "independent species call.", "",
         "**No cohort is decomposed unless its committed confusion matrix reproduces exactly** from the",
         "cached AMRFinder runs through the deployed `call_resistance`. A cohort with incomplete cached",
         "runs reports composition and **no** outcome cross-tab.", "",
         "| cohort | scored as | composition | off-species | reconciled | FN off-species |",
         "|---|---|---|---|---|---|"]
    for r in rows:
        comp = r["composition"]
        v = comp["verdict"]
        off = (comp.get("n_same_genus_other_species", 0) + comp.get("n_other_genus", 0)
               if "n_same_genus_other_species" in comp else None)
        offs = "n/a" if off is None else str(off)
        fna = r.get("false_negative_anatomy")
        fn_txt = (f"{fna['n_fn_off_species']}/{fna['n_false_negatives']}" if fna else "—")
        L.append(f"| `{r['cohort_dir']}` | `{r['amrfinder_organism']}` | {v} | {offs} | "
                 f"{'yes' if r.get('reconciled') else 'no'} | {fn_txt} |")
    L.append("")
    for r in rows:
        comp = r["composition"]
        L.append(f"## {r['cohort_dir']}  (`{r['drug']}`, scored as `{r['amrfinder_organism']}`)")
        L.append("")
        L.append(f"- verdict: **{comp['verdict']}**")
        if "species_counts" in comp:
            for sp, n in sorted(comp["species_counts"].items(), key=lambda kv: -kv[1]):
                L.append(f"  - {n}x `{sp}`  ({r['species_vocabulary_status'].get(sp, '?')})")
        else:
            L.append(f"  - refused: {comp.get('reason')}")
        if not r.get("reconciled"):
            L.append(f"- outcome cross-tab **withheld**: {r.get('crosstab_withheld_reason')}")
        else:
            fna = r["false_negative_anatomy"]
            L.append(f"- committed confusion reproduced exactly: {r['reconciled_confusion']}")
            L.append(f"- false negatives: **{fna['n_false_negatives']}**, of which "
                     f"**{fna['n_fn_off_species']} are off-species** and "
                     f"{fna['n_fn_expected_species']} are the expected species")
            for sp, n in sorted(fna["fn_species_counts"].items(), key=lambda kv: -kv[1]):
                L.append(f"  - FN: {n}x `{sp}`")
        L.append("")
    return "\n".join(L) + "\n"


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cohort", default=None, help="audit only this cohort directory name")
    ap.add_argument("--raw", type=Path, default=RAW)
    ap.add_argument("--wiki", type=Path, default=WIKI)
    ap.add_argument("--out-date", default=None, help="override the artifact date stamp")
    a = ap.parse_args(argv)

    cohorts = provdisjoint_cohorts(a.raw)
    if a.cohort:
        cohorts = [c for c in cohorts if c.name == a.cohort]
        if not cohorts:
            print(f"no such provdisjoint cohort: {a.cohort}", file=sys.stderr)
            return 2
    if not cohorts:
        print("no provenance-disjoint cohorts found (is data/ mounted?)", file=sys.stderr)
        return 2

    rows: list[dict] = []
    for c in cohorts:
        try:
            rows.append(audit_cohort(c, wiki=a.wiki))
            print(f"  audited {c.name}")
        except AuditError as e:
            print(f"  SKIP {c.name}: {e}", file=sys.stderr)

    today = a.out_date or date.today().isoformat()
    payload = {
        "schema": SCHEMA,
        "date": today,
        "n_cohorts": len(rows),
        "scope": ("read-only species-composition disclosure for the provenance-disjoint arm; does NOT "
                  "re-score any cell and does NOT modify any frozen artifact"),
        "label_provenance": ("each assembly's own GenBank ORGANISM line -- the submitter's assertion, not "
                             "a wet-lab identification (the project's own G1 circular-label class)"),
        # A LIST of cells each carrying `organism` + `drug`, matching the sibling
        # `provdisjoint_source_concentration.json` exactly -- its loader re-derives
        # `canonical_cell_key(organism, drug)` rather than reading a pre-stringified key. Emitting a
        # dict keyed by a stringified tuple would have been a second convention for the same join.
        "cells": rows,
        "complete": all(r["composition"]["verdict"] != "INSUFFICIENT_RESOLUTION" for r in rows),
    }
    (a.wiki / f"provdisjoint_species_audit_{today}.json").write_text(
        json.dumps(payload, indent=2, default=str), encoding="utf-8")
    (a.wiki / f"provdisjoint_species_audit_{today}.md").write_text(
        render_md(rows, today), encoding="utf-8")
    print(f"[species audit -> wiki/provdisjoint_species_audit_{today}.{{md,json}}]  {len(rows)} cohorts")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
