"""Cipro AMRFinderPlus mechanism audit across all N=38 cohort strains.

Experiment 1 from the Open Questions brainstorm. Runs AMRFinderPlus against
each strain's RefSeq fasta, parses both `main.tsv` (acquired-gene hits) and
`mutations.tsv` (POINT mutations including QRDR), classifies each hit into a
mechanism class, and builds a per-strain table.

Cipro mechanism classes (per textbook + Codex synthesis):
- QRDR_target_alteration: gyrA/B + parC/E POINT mutations
- plasmid_protect_modify: qnr* + aac(6')-Ib-cr (and variants)
- efflux: acrAB-TolC + oqxAB + mdfA + mdtK
- porin_loss: ompC/ompF disruptions (rarely reported by AMRFinder)
- regulatory: marR/marA + soxRS

Output: wiki/cipro_mechanism_audit_<date>.md + .json sidecar.

Runtime: ~95 sec/strain on the install-smoke ST131 case; estimate ~1 hr for 38
strains. Use --start-index / --end-index for resumability if interrupted.
"""
from __future__ import annotations

import argparse
import json
import sys
import tempfile
import time
from collections import defaultdict
from datetime import date as _date
from pathlib import Path

import pandas as pd

from dna_decode.data.cohort import load_cohort
from dna_decode.data.refseq import fasta_path
from tools.docker_runner import DockerRunnerError, run as docker_run


AMRFINDER_IMAGE = "ncbi/amr:4.2.7-2026-03-24.1"
AMRFINDER_DB = "C:/Users/Farshad/dna_decode_stage2/amrfinder_db"

CIPRO_LOCI_BY_MECHANISM: dict[str, set[str]] = {
    "QRDR_target_alteration": {"gyrA", "gyrB", "parC", "parE"},
    "plasmid_protect_modify": {
        "qnrA", "qnrB", "qnrC", "qnrD", "qnrS",
        "aac(6')-Ib-cr", "aac(6')-Ib", "aac6-Ib-cr",
    },
    "efflux": {"acrA", "acrB", "tolC", "oqxA", "oqxB", "mdfA", "mdtK"},
    "porin_loss": {"ompC", "ompF"},
    "regulatory": {"marR", "marA", "marB", "soxR", "soxS"},
}
QUINOLONE_CLASSES = {"QUINOLONE", "FLUOROQUINOLONE"}


def _classify_symbol(sym: str) -> str:
    """Return mechanism class for a gene symbol; '' if not cipro-relevant."""
    if not sym:
        return ""
    sym_norm = sym.split("_")[0].strip()  # gyrA_S83L -> gyrA
    for mech, loci in CIPRO_LOCI_BY_MECHANISM.items():
        if sym_norm in loci:
            return mech
        # tolerant prefix match (e.g., qnrB19 -> qnrB)
        for locus in loci:
            if sym_norm.startswith(locus) and sym_norm != locus:
                return mech
    return ""


def _is_synonymous_point(sym: str) -> bool:
    """Detect a synonymous POINT mutation like 'gyrA_G141G' (first AA == last AA).

    Format expected: '<gene>_<wt><pos><alt>'. If parsing fails, returns False
    (don't filter — let the row through).
    """
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


def _parse_amrfinder_outputs(main_tsv: Path, mut_tsv: Path) -> dict:
    """Parse main + mutations TSVs into a structured per-strain summary."""
    hits: list[dict] = []
    if main_tsv.exists() and main_tsv.stat().st_size > 0:
        try:
            main_df = pd.read_csv(main_tsv, sep="\t", dtype=str, keep_default_na=False)
            for _, row in main_df.iterrows():
                sym = row.get("Gene symbol", "") or row.get("Element symbol", "")
                cls = row.get("Class", "")
                hits.append({
                    "kind": "acquired",
                    "symbol": sym,
                    "class": cls,
                    "subclass": row.get("Subclass", ""),
                    "method": row.get("Method", ""),
                    "scope": row.get("Scope", ""),
                    "type": row.get("Element type", ""),
                    "identity": row.get("% Identity to reference sequence", "") or row.get("% Identity to reference", ""),
                    "coverage": row.get("% Coverage of reference sequence", "") or row.get("% Coverage of reference", ""),
                    "mechanism": _classify_symbol(sym),
                })
        except Exception as e:
            print(f"  [parse] WARN main.tsv parse failed: {e!r}")
    if mut_tsv.exists() and mut_tsv.stat().st_size > 0:
        try:
            mut_df = pd.read_csv(mut_tsv, sep="\t", dtype=str, keep_default_na=False)
            for _, row in mut_df.iterrows():
                sym = row.get("Gene symbol", "") or row.get("Element symbol", "")
                cls = row.get("Class", "")
                # mutations.tsv reports EVERY position vs reference (including
                # synonymous + non-cipro classes). Filter to cipro-relevant only.
                if _is_synonymous_point(sym):
                    continue
                if cls.upper() not in QUINOLONE_CLASSES:
                    continue
                mech = _classify_symbol(sym)
                if not mech:
                    continue
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
    quinolone_class_hits: list[str] = []
    for h in hits:
        if h["mechanism"]:
            mech_hits[h["mechanism"]].append(h["symbol"])
        if h["class"].upper() in QUINOLONE_CLASSES and h["symbol"] not in mech_hits.get(h["mechanism"], []):
            quinolone_class_hits.append(h["symbol"])
    mechs_found = sorted(mech_hits.keys())
    primary = mechs_found[0] if mechs_found else ("QUINOLONE_class_unclassified" if quinolone_class_hits else "NO_MECHANISM")
    return {
        "hits": hits,
        "n_hits": len(hits),
        "mechanisms_present": mechs_found,
        "mech_hits": {k: sorted(set(v)) for k, v in mech_hits.items()},
        "quinolone_class_hits": sorted(set(quinolone_class_hits)),
        "primary_mechanism_class": primary,
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cohort", type=Path, default=Path("data/processed/gate_b_n40_cipro_cohort.parquet"))
    parser.add_argument("--refseq-cache", type=Path, default=Path("D:/dna_decode_cache/refseq"))
    parser.add_argument("--out-root", type=Path, default=Path("data/amrfinder_runs"))
    parser.add_argument("--drug", default="ciprofloxacin")
    parser.add_argument("--output", type=Path, default=None)
    parser.add_argument("--start-index", type=int, default=0)
    parser.add_argument("--end-index", type=int, default=None)
    parser.add_argument("--per-strain-timeout", type=float, default=600.0)
    parser.add_argument("--skip-existing", action="store_true", default=True)
    args = parser.parse_args(argv)

    if args.output is None:
        args.output = Path(f"wiki/cipro_mechanism_audit_{_date.today().isoformat()}.md")

    cohort = load_cohort(args.cohort)
    drug_lower = args.drug.lower()
    strains = [s for s in cohort.strains if drug_lower in s.ast_labels]
    strains.sort(key=lambda s: s.strain_id)
    end = args.end_index if args.end_index is not None else len(strains)
    target = strains[args.start_index:end]
    print(f"[mech_audit] cohort: {len(strains)} {args.drug}-labeled strains; processing indices [{args.start_index}:{end}) = {len(target)} strains")

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

        parsed = _parse_amrfinder_outputs(main_tsv, mut_tsv)
        per_strain_results.append({
            "strain_id": s.strain_id,
            "accession": s.assembly_accession,
            "cohort_label": int(s.ast_labels[drug_lower]),
            "mlst": s.mlst,
            "status": "OK",
            **parsed,
        })

    # Aggregate across cohort
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

    qrdr_in_R = sum(1 for r in r_strains if "QRDR_target_alteration" in r.get("mechanisms_present", []))
    plasmid_in_R = sum(1 for r in r_strains if "plasmid_protect_modify" in r.get("mechanisms_present", []))

    verdict = (
        "QRDR_DOMINANT" if qrdr_in_R / max(1, len(r_strains)) >= 0.7 else
        "MIXED_MECHANISMS" if (qrdr_in_R + plasmid_in_R) / max(1, len(r_strains)) >= 0.5 else
        "MOSTLY_UNKNOWN"
    )

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

    # JSON sidecar
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
        "qrdr_in_R_count": qrdr_in_R,
        "plasmid_in_R_count": plasmid_in_R,
        "verdict": verdict,
        "elapsed_sec": elapsed_total,
        "per_strain": per_strain_results,
    }
    json_path.write_text(json.dumps(payload, indent=2, default=str), encoding="utf-8")
    print(f"[mech_audit] wrote JSON: {json_path}")

    # Markdown packet
    lines = [
        f"# Cipro AMRFinderPlus mechanism audit — N={len(per_strain_results)} cohort ({_date.today().isoformat()})",
        "",
        f"**Purpose:** classify the actual cipro-resistance mechanism in each R strain (and check for silent mechanisms in S strains) before judging whether NT's preflight INCONCLUSIVE_MISS reflects TRUE biology vs a model failure.",
        f"**Tool:** AMRFinderPlus {AMRFINDER_IMAGE} on `Escherichia` organism mode + `--mutation_all`.",
        f"**Cohort:** `{args.cohort}` ({len(strains)} {args.drug}-labeled strains).",
        "",
        f"## Verdict: **{verdict}**",
        f"- QRDR mechanism found in: **{qrdr_in_R}** / {len(r_strains)} R strains",
        f"- Plasmid (qnr / aac6-Ib-cr) found in: **{plasmid_in_R}** / {len(r_strains)} R strains",
        f"- R strains with NO known cipro mechanism: **{r_no_known}** / {len(r_strains)}",
        f"- S strains with silent mechanism hit: **{s_silent_mech}** / {len(s_strains)}",
        "",
        "## Mechanism distribution (R set)",
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
        "## How to use",
        "",
        "- If verdict is QRDR_DOMINANT: NT's preflight INCONCLUSIVE_MISS is a MODEL failure (the signal is there in 70%+ of R strains; NT just isn't finding it).",
        "- If verdict is MIXED_MECHANISMS: cipro resistance is heterogeneous; per-gene attribution should NOT expect to converge on one locus. Use per-strain known-mechanism overlap as ground truth for attribution.",
        "- If verdict is MOSTLY_UNKNOWN: the biology is unclear; AMRFinderPlus alone may miss novel/regulatory mechanisms. Need Bakta gene-presence + downstream curated baseline.",
        f"- Silent-S count of {s_silent_mech} suggests label noise upper bound — these strains may be mislabeled or have functional-but-not-clinical-MIC resistance.",
        "",
        f"_JSON sidecar: `{json_path}`_",
    ])
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text("\n".join(lines), encoding="utf-8")
    print(f"[mech_audit] wrote packet: {args.output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
