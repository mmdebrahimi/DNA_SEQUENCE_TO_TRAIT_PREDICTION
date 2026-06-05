"""Drug-agnostic AMRFinderPlus mechanism audit.

Generalizes the cipro-only `scripts/cipro_mechanism_audit.py` to any drug
supported by `dna_decode/data/mic_tiers.py` (cipro / ceftriaxone / tetracycline
/ gentamicin). Pulls per-drug catalogs (loci by mechanism, AMRFinder Class
filter, gene-symbol classifier) from `mic_tiers.py` rather than hardcoding.

Why this exists as a separate script:
- `scripts/cipro_mechanism_audit.py` has 80+ dedicated cipro-specific tests
  pinning verdict naming + behavior. Refactoring it to take --drug would
  break those tests + offer no value beyond drift prevention.
- New audits for cef / tet / gent (v0.1+ work) use THIS script.
- Some helpers (`_is_synonymous_point`, `_run_amrfinder`) are duplicated
  from cipro_mechanism_audit.py with a TODO marker; future consolidation
  is drift-prevention only, not load-bearing.

Output:
- wiki/<drug>_mechanism_audit_<DATE>.{md,json}

Exit code: 0 success / 1 cohort error / 2 IO/runtime error.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from collections import defaultdict
from datetime import date as _date
from pathlib import Path

import pandas as pd

from dna_decode.data.cohort import load_cohort
from dna_decode.data.mic_tiers import (
    amrfinder_classes_for,
    classify_gene_symbol,
    loci_by_mechanism_for,
    primary_mechanisms_for,
    supported_drugs,
)
from dna_decode.data.refseq import fasta_path
from tools.docker_runner import DockerRunnerError, run as docker_run


# Hardcoded image + DB path matching cipro_mechanism_audit.py CLAUDE.md pin.
AMRFINDER_IMAGE = "ncbi/amr:4.2.7-2026-03-24.1"
AMRFINDER_DB = str(
    Path(
        os.environ.get(
            "DNA_DECODE_AMRFINDER_DB",
            # machine-agnostic default: ~/dna_decode_stage2/amrfinder_db resolves on BOTH the
            # Farshad laptop and the b0652085 workhorse (was hardcoded to the workhorse path —
            # cross-machine drift that broke the Docker mount on the laptop, 2026-06-05).
            str(Path.home() / "dna_decode_stage2" / "amrfinder_db"),
        )
    )
)


# TODO(consolidation): duplicated from scripts/cipro_mechanism_audit.py.
# Extract to a shared helper module if a third caller appears.
def _is_synonymous_point(sym: str) -> bool:
    """Detect a synonymous POINT mutation like 'gyrA_G141G' (first AA == last AA)."""
    if "_" not in sym:
        return False
    parts = sym.rsplit("_", 1)
    if len(parts) != 2:
        return False
    mut = parts[1]
    if len(mut) < 3:
        return False
    wt = mut[0]
    alt = mut[-1]
    middle = mut[1:-1]
    if not middle.isdigit():
        return False
    return wt == alt


# TODO(consolidation): duplicated from scripts/cipro_mechanism_audit.py.
def _run_amrfinder(fasta: Path, out_dir: Path, timeout_sec: float = 600) -> tuple[Path, Path]:
    """Invoke AMRFinderPlus on a single fasta; return (main_tsv, mutations_tsv)."""
    main_out = out_dir / "main.tsv"
    mut_out = out_dir / "mutations.tsv"
    docker_run(
        AMRFINDER_IMAGE,
        [
            "amrfinder",
            "-n", f"/in/{fasta.name}",
            "-O", "Escherichia",
            "--database", "/db/latest",
            "--mutation_all", "/out/mutations.tsv",
            "-o", "/out/main.tsv",
        ],
        mounts={
            str(fasta.parent): "/in:ro",
            AMRFINDER_DB: "/db:ro",
            str(out_dir): "/out",
        },
        timeout=timeout_sec,
    )
    return main_out, mut_out


def parse_amrfinder_outputs_for_drug(
    main_tsv: Path,
    mut_tsv: Path,
    drug: str,
) -> dict:
    """Drug-agnostic AMRFinder parser.

    Filters main.tsv + mutations.tsv to (a) classes relevant to the drug per
    mic_tiers.amrfinder_classes_for(drug) AND (b) gene symbols that classify
    to a known mechanism via mic_tiers.classify_gene_symbol(drug, symbol).
    Dedupes mutations across main.tsv POINTX rows + mutations.tsv per symbol.

    Returns:
        dict with hits, n_hits, mechanisms_present, mech_hits, primary_mechanism_class.
    """
    relevant_classes = amrfinder_classes_for(drug)
    hits: list[dict] = []
    seen_mutations: set[str] = set()

    if main_tsv.exists() and main_tsv.stat().st_size > 0:
        try:
            main_df = pd.read_csv(main_tsv, sep="\t", dtype=str, keep_default_na=False)
            for _, row in main_df.iterrows():
                sym = row.get("Gene symbol", "") or row.get("Element symbol", "")
                cls = row.get("Class", "")
                method = row.get("Method", "") or ""
                is_point = "POINT" in method.upper()
                kind = "mutation" if is_point else "acquired"
                if kind == "mutation":
                    if _is_synonymous_point(sym):
                        continue
                    if cls.upper() not in relevant_classes:
                        continue
                    if not classify_gene_symbol(drug, sym):
                        continue
                    seen_mutations.add(sym)
                else:
                    # Acquired-gene hits: require the symbol to classify
                    # AND/OR the class to be in the drug's relevant set.
                    sym_match = bool(classify_gene_symbol(drug, sym))
                    cls_match = cls.upper() in relevant_classes
                    if not (sym_match or cls_match):
                        continue
                hits.append({
                    "kind": kind,
                    "symbol": sym,
                    "class": cls,
                    "subclass": row.get("Subclass", ""),
                    "method": method,
                    "scope": row.get("Scope", ""),
                    "type": row.get("Element type", ""),
                    "identity": row.get("% Identity to reference sequence", "") or row.get("% Identity to reference", ""),
                    "coverage": row.get("% Coverage of reference sequence", "") or row.get("% Coverage of reference", ""),
                    "mechanism": classify_gene_symbol(drug, sym),
                })
        except Exception as e:
            print(f"  [parse] WARN main.tsv parse failed: {e!r}")

    if mut_tsv.exists() and mut_tsv.stat().st_size > 0:
        try:
            mut_df = pd.read_csv(mut_tsv, sep="\t", dtype=str, keep_default_na=False)
            for _, row in mut_df.iterrows():
                sym = row.get("Gene symbol", "") or row.get("Element symbol", "")
                cls = row.get("Class", "")
                if _is_synonymous_point(sym):
                    continue
                if cls.upper() not in relevant_classes:
                    continue
                mech = classify_gene_symbol(drug, sym)
                if not mech:
                    continue
                if sym in seen_mutations:
                    continue
                seen_mutations.add(sym)
                hits.append({
                    "kind": "mutation",
                    "symbol": sym,
                    "class": cls,
                    "subclass": row.get("Subclass", ""),
                    "method": row.get("Method", ""),
                    "scope": row.get("Scope", ""),
                    "type": row.get("Element type", ""),
                    "identity": row.get("% Identity to reference sequence", "") or row.get("% Identity to reference", ""),
                    "coverage": row.get("% Coverage of reference sequence", "") or row.get("% Coverage of reference", ""),
                    "mechanism": mech,
                })
        except Exception as e:
            print(f"  [parse] WARN mutations.tsv parse failed: {e!r}")

    # Aggregate per mechanism
    mech_hits: dict[str, list[str]] = defaultdict(list)
    for h in hits:
        if h["mechanism"]:
            mech_hits[h["mechanism"]].append(h["symbol"])
    mechs_found = sorted(mech_hits.keys())

    # Primary mechanism = first primary-set hit if any, else first mechanism found, else NO_MECHANISM.
    primary_set = primary_mechanisms_for(drug)
    primary = ""
    for m in mechs_found:
        if m in primary_set:
            primary = m
            break
    if not primary:
        primary = mechs_found[0] if mechs_found else "NO_MECHANISM"

    return {
        "hits": hits,
        "n_hits": len(hits),
        "mechanisms_present": mechs_found,
        "mech_hits": {k: sorted(set(v)) for k, v in mech_hits.items()},
        "primary_mechanism_class": primary,
    }


def compute_verdict(
    r_strains: list[dict],
    drug: str,
    qrdr_dominant_threshold: float = 0.70,
    mixed_threshold: float = 0.50,
) -> tuple[str, dict[str, int]]:
    """Drug-agnostic verdict: PRIMARY_DOMINANT / MIXED_MECHANISMS / MOSTLY_UNKNOWN.

    PRIMARY_DOMINANT: ≥ 70% of R strains have ≥ 1 primary mechanism (per
        mic_tiers.primary_mechanisms_for(drug)).
    MIXED_MECHANISMS: ≥ 50% of R strains have ≥ 1 primary OR co-resistance
        mechanism but below the PRIMARY_DOMINANT bar.
    MOSTLY_UNKNOWN: < 50% of R strains have any known mechanism.

    Returns (verdict, breakdown_counts).
    """
    primary_set = primary_mechanisms_for(drug)
    n_r = len(r_strains)
    if n_r == 0:
        return "EMPTY_R_SET", {}
    n_with_primary = 0
    n_with_any = 0
    for r in r_strains:
        mechs = set(r.get("mechanisms_present", []))
        if mechs & set(primary_set):
            n_with_primary += 1
        if mechs:
            n_with_any += 1
    breakdown = {
        "n_R_strains": n_r,
        "n_R_with_primary_mechanism": n_with_primary,
        "n_R_with_any_mechanism": n_with_any,
    }
    primary_frac = n_with_primary / n_r
    any_frac = n_with_any / n_r
    if primary_frac >= qrdr_dominant_threshold:
        return "PRIMARY_DOMINANT", breakdown
    if any_frac >= mixed_threshold:
        return "MIXED_MECHANISMS", breakdown
    return "MOSTLY_UNKNOWN", breakdown


def _target_strains_for_drug(cohort, drug: str, *, all_labeled: bool = False):
    drug_lower = drug.lower()
    if all_labeled:
        return [s for s in cohort.strains if drug_lower in s.ast_labels], "labeled"
    pool_ids = cohort.per_drug_strain_ids.get(drug_lower, [])
    if pool_ids:
        pool_set = set(pool_ids)
        return [s for s in cohort.strains if s.strain_id in pool_set], "pool"
    return [s for s in cohort.strains if drug_lower in s.ast_labels], "labeled"


def main(argv: list[str] | None = None) -> int:
    global AMRFINDER_DB, AMRFINDER_IMAGE
    parser = argparse.ArgumentParser(description="Drug-agnostic AMRFinderPlus mechanism audit.")
    parser.add_argument("--drug", required=True, choices=supported_drugs(),
                        help="Drug to audit (per mic_tiers.supported_drugs())")
    parser.add_argument("--cohort", type=Path, required=True,
                        help="Cohort parquet path")
    parser.add_argument("--refseq-cache", type=Path, required=True,
                        help="RefSeq cache root (FASTAs live at <root>/<accession>/genome.fna)")
    parser.add_argument("--out-root", type=Path, default=Path("data/amrfinder_runs"),
                        help="Per-strain AMRFinder output root")
    parser.add_argument("--output", type=Path, default=None,
                        help="Output markdown path; JSON sidecar emitted alongside")
    parser.add_argument("--start-index", type=int, default=0)
    parser.add_argument("--end-index", type=int, default=None)
    parser.add_argument("--per-strain-timeout", type=float, default=600.0)
    parser.add_argument("--skip-existing", action="store_true", default=True)
    parser.add_argument(
        "--all-labeled",
        action="store_true",
        help="Process all labeled strains instead of the saved per-drug in_pool_<drug> subset.",
    )
    parser.add_argument(
        "--amrfinder-db",
        type=Path,
        default=Path(AMRFINDER_DB),
        help="Local AMRFinderPlus DB root; container will read /db/latest.",
    )
    parser.add_argument(
        "--amrfinder-image",
        default=AMRFINDER_IMAGE,
        help="Pinned AMRFinderPlus Docker image tag.",
    )
    args = parser.parse_args(argv)

    drug_lower = args.drug.lower()
    if args.output is None:
        args.output = Path(f"wiki/{drug_lower}_mechanism_audit_{_date.today().isoformat()}.md")

    AMRFINDER_DB = str(args.amrfinder_db)
    AMRFINDER_IMAGE = args.amrfinder_image

    cohort = load_cohort(args.cohort)
    strains, cohort_scope = _target_strains_for_drug(cohort, args.drug, all_labeled=args.all_labeled)
    strains.sort(key=lambda s: s.strain_id)
    end = args.end_index if args.end_index is not None else len(strains)
    target = strains[args.start_index:end]
    print(f"[drug_mech_audit] drug={args.drug} cohort: {len(strains)} {cohort_scope} strains; "
          f"processing indices [{args.start_index}:{end}) = {len(target)} strains")

    args.out_root.mkdir(parents=True, exist_ok=True)
    per_strain_results: list[dict] = []
    t0 = time.time()
    for i, s in enumerate(target):
        idx = args.start_index + i
        fna = fasta_path(s.assembly_accession, args.refseq_cache)
        if not fna.exists():
            print(f"  [{idx}] {s.strain_id} ({s.assembly_accession}): MISSING fasta at {fna}")
            per_strain_results.append({
                "strain_id": s.strain_id,
                "accession": s.assembly_accession,
                "cohort_label": int(s.ast_labels[drug_lower]),
                "status": "MISSING_FASTA",
            })
            continue

        out_dir = args.out_root / s.assembly_accession
        out_dir.mkdir(parents=True, exist_ok=True)
        main_tsv = out_dir / "main.tsv"
        mut_tsv = out_dir / "mutations.tsv"

        if args.skip_existing and main_tsv.exists() and mut_tsv.exists():
            print(f"  [{idx}] {s.strain_id} ({s.assembly_accession}): skip (cached)")
        else:
            elapsed = time.time() - t0
            print(f"  [{idx}] {s.strain_id} ({s.assembly_accession}): running AMRFinder... (t+{elapsed:.0f}s)")
            try:
                _run_amrfinder(fna, out_dir, timeout_sec=args.per_strain_timeout)
            except DockerRunnerError as e:
                print(f"  [{idx}] FAIL: {e}")
                per_strain_results.append({
                    "strain_id": s.strain_id,
                    "accession": s.assembly_accession,
                    "cohort_label": int(s.ast_labels[drug_lower]),
                    "status": "DOCKER_FAIL",
                    "error": str(e),
                })
                continue

        parsed = parse_amrfinder_outputs_for_drug(main_tsv, mut_tsv, args.drug)
        per_strain_results.append({
            "strain_id": s.strain_id,
            "accession": s.assembly_accession,
            "cohort_label": int(s.ast_labels[drug_lower]),
            "mlst": s.mlst,
            "status": "OK",
            **parsed,
        })

    r_strains = [r for r in per_strain_results if r.get("cohort_label") == 1]
    s_strains = [r for r in per_strain_results if r.get("cohort_label") == 0]

    r_mech_counts: dict[str, int] = defaultdict(int)
    s_mech_counts: dict[str, int] = defaultdict(int)
    r_no_known = 0
    s_silent_mech = 0
    for r in r_strains:
        mechs = r.get("mechanisms_present", [])
        if not mechs:
            r_no_known += 1
        for m in mechs:
            r_mech_counts[m] += 1
    for r in s_strains:
        mechs = r.get("mechanisms_present", [])
        if mechs:
            s_silent_mech += 1
        for m in mechs:
            s_mech_counts[m] += 1

    verdict, breakdown = compute_verdict(r_strains, args.drug)
    elapsed_total = time.time() - t0

    print(f"\n=== Mechanism distribution (R set, n={len(r_strains)}) ===")
    for m, n in sorted(r_mech_counts.items(), key=lambda x: -x[1]):
        print(f"  {m:30s} {n}")
    print(f"  (R with no known mechanism: {r_no_known})")
    print(f"\n=== Mechanism distribution (S set, n={len(s_strains)}) ===")
    for m, n in sorted(s_mech_counts.items(), key=lambda x: -x[1]):
        print(f"  {m:30s} {n}")
    print(f"  (S with silent mechanism hit: {s_silent_mech})")
    print(f"\nVerdict: {verdict}")
    print(f"Total wallclock: {elapsed_total:.0f}s")

    json_path = args.output.with_suffix(".json")
    json_path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "cohort_path": str(args.cohort),
        "refseq_cache": str(args.refseq_cache),
        "drug": args.drug,
        "amrfinder_image": AMRFINDER_IMAGE,
        "amrfinder_db_path": AMRFINDER_DB,
        "n_processed": len(per_strain_results),
        "r_mech_counts": dict(r_mech_counts),
        "s_mech_counts": dict(s_mech_counts),
        "r_no_known_mechanism": r_no_known,
        "s_silent_mechanism_hits": s_silent_mech,
        "verdict": verdict,
        "verdict_breakdown": breakdown,
        "elapsed_sec": elapsed_total,
        "per_strain": per_strain_results,
    }
    json_path.write_text(json.dumps(payload, indent=2, default=str), encoding="utf-8")
    print(f"[drug_mech_audit] wrote JSON: {json_path}")

    lines = [
        f"# {args.drug.title()} AMRFinderPlus mechanism audit - N={len(per_strain_results)} cohort ({_date.today().isoformat()})",
        "",
        f"**Drug:** {args.drug}",
        f"**Tool:** AMRFinderPlus `{AMRFINDER_IMAGE}` on `Escherichia` organism mode + `--mutation_all`.",
        f"**Cohort:** `{args.cohort}` ({len(strains)} {cohort_scope} strains for {args.drug}).",
        f"**Per-drug catalogs:** sourced from `dna_decode/data/mic_tiers.py`.",
        "",
        f"## Verdict: **{verdict}**",
        f"- Primary mechanism found in: **{breakdown.get('n_R_with_primary_mechanism', 0)}** / {len(r_strains)} R strains",
        f"- Any known mechanism found in: **{breakdown.get('n_R_with_any_mechanism', 0)}** / {len(r_strains)} R strains",
        f"- R strains with NO known mechanism: **{r_no_known}** / {len(r_strains)}",
        f"- S strains with silent mechanism hit: **{s_silent_mech}** / {len(s_strains)}",
        "",
        "## Mechanism distribution",
        "",
        "| mechanism | R count | S count |",
        "|---|---:|---:|",
    ]
    all_mechs = sorted(set(r_mech_counts) | set(s_mech_counts))
    for m in all_mechs:
        lines.append(f"| {m} | {r_mech_counts.get(m, 0)} | {s_mech_counts.get(m, 0)} |")
    lines.extend([
        "",
        "## Per-strain mechanism table",
        "",
        "| strain_id | accession | label | status | primary mech | mechanisms | mlst |",
        "|---|---|---|---|---|---|---|",
    ])
    for r in sorted(per_strain_results, key=lambda x: (x.get("cohort_label", 0), x["strain_id"])):
        label = "R" if r.get("cohort_label") == 1 else "S"
        prim = r.get("primary_mechanism_class", r.get("status", ""))
        mechs = ",".join(r.get("mechanisms_present", []))
        lines.append(
            f"| {r['strain_id']} | {r['accession']} | {label} | {r.get('status','?')} | {prim} | {mechs} | {r.get('mlst','')} |"
        )
    lines.extend([
        "",
        "## Verdict interpretation",
        "",
        "- **PRIMARY_DOMINANT:** ≥ 70% of R strains carry a primary mechanism for this drug. NT/classifier should find this signal; if it does not, model failure is implicated.",
        "- **MIXED_MECHANISMS:** 50-69% of R strains carry primary OR co-resistance mechanisms. Per-gene attribution is unlikely to converge on one locus.",
        "- **MOSTLY_UNKNOWN:** < 50% of R strains have any known mechanism. AMRFinder alone may miss novel/regulatory mechanisms; pair with Bakta gene-presence + curated baseline.",
        "",
        f"_JSON sidecar: `{json_path}`_",
    ])
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text("\n".join(lines), encoding="utf-8")
    print(f"[drug_mech_audit] wrote packet: {args.output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
