"""Single-genome map CLI — point the genome-map at ONE genome.

The product's primary surface ("point the tool at ONE microbial genome → an
honest, evidence-tiered per-feature map"). The 3-genome spike
(`scripts/genome_map_spike.py`) validated the pipeline (verdict GO 2026-06-18);
this is the thin, reusable single-genome front door over the SAME tested core
(`run_genome_map_for` + the Step-1 runners + the Step-6 gate).

Modes:
  - Live:    --genome-fasta X.fna [--organism Escherichia]  → Bakta + AMRFinder via Docker.
  - Hybrid:  --genome-fasta X.fna --gff Y.gff3              → skip Bakta, run AMRFinder.
  - Offline: --gff Y.gff3 --no-amrfinder                    → tiers from a provided GFF
             + the determinant cells; degraded-coverage flag set (AC12). No Docker.

Emits a JSON map + a flat feature table + a readable markdown summary, and prints
the per-tier counts + DB-labelled unknown rate + the G1/G2 gate result.

Usage:
    MSYS_NO_PATHCONV=1 uv run python -m scripts.genome_map \
        --genome-fasta D:/dna_decode_cache/refseq/GCA_x/genome.fna \
        --organism Escherichia --sample-id GCA_x --out-dir wiki/genome_map_GCA_x
Exit: 0 = map emitted, 2 = usage/IO error, 3 = annotation BLOCKED (no map).
"""
from __future__ import annotations

import argparse
import json
import sys
from datetime import date as _date
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from dna_decode.genome_map.build_map import build_feature_table
from dna_decode.genome_map.gate import evaluate_gate
from scripts.genome_map_spike import run_genome_map_for


def render_genome_summary_md(genome_map: dict, gate_result: dict, *, generated: str | None = None) -> str:
    """Readable per-genome markdown summary (pure — testable without live tools)."""
    generated = generated or _date.today().isoformat()
    m = genome_map["metrics"]
    lines: list[str] = []
    acc = genome_map.get("genome_accession")
    lines.append(f"# Genome map — {acc} ({generated})")
    lines.append("")
    overlay_status = genome_map.get("overlay_status")
    if overlay_status and overlay_status != "FULL":
        lines.append(f"> **overlay_status: {overlay_status}** — no determinant overlay; "
                     "tiers from the GFF only. Phenotype claims require a full overlay.")
        lines.append("")
    elif genome_map.get("degraded_coverage"):
        lines.append("> **DEGRADED COVERAGE** — no determinant overlay; tiers from the GFF only.")
        lines.append("")
    if overlay_status:
        lines.append(f"- overlay_status: `{overlay_status}`")
    lines.append(f"- organism (-O): `{genome_map.get('amrfinder_organism')}`")
    lines.append(f"- total features: {m['total_features']}")
    lines.append("- per-tier counts: " + ", ".join(f"{k}={v}" for k, v in m["per_tier_counts"].items()))
    lines.append(f"- **`unknown_under_bakta_db_light`: {m['unknown_under_bakta_db_light']:.3f}** "
                 "(db-light coverage caveat is IN the field name — not biological unknown)")
    jq = m["join_quality"]
    lines.append(f"- determinant join quality: n_main_rows={jq.get('n_main_rows')} "
                 f"high_confidence={jq.get('n_high_confidence_join')} "
                 f"symbol_fallback={jq.get('n_symbol_fallback')} unjoined={jq.get('n_unjoined')}")
    lines.append(f"- determinant-phenotype features: {m['determinant_phenotype_feature_count']}")
    gl = ", ".join(f"{d}={v.get('prediction')}" for d, v in m["genome_level_calls"].items())
    lines.append(f"- genome-level R/S calls (separate from features): {gl or '(none)'}")
    lines.append("")
    lines.append("## Determinant-phenotype features (the only tier with a phenotype claim)")
    if m["determinant_phenotype_features"]:
        for f in m["determinant_phenotype_features"]:
            drugs = sorted({(p.get("drug") or p.get("phenotype") or "") for p in f["phenotype"]})
            label = f["raw_gene_symbol"] or f["raw_product"][:50]
            lines.append(f"- `{label}` @ {f['seqid']}:{f['start']}-{f['end']} → "
                         f"[{', '.join(s for s in drugs if s)}]")
    else:
        lines.append("- (none — no curated determinant hit this genome; the MODAL case for an arbitrary genome)")
    lines.append("")
    lines.append("## Honesty gate (informational for a single genome)")
    lines.append(f"- G1 prevent-wrong-inference features: {len(gate_result['g1_features'])} "
                 f"(demote={gate_result['g1_demote_count']}, surface={gate_result['g1_surface_count']})")
    lines.append(f"- G2 phenotype-wall: pass={gate_result['g2_spotcheck']['pass']} "
                 f"(violations={len(gate_result['g2_spotcheck']['violations'])})")
    lines.append(f"- all_joins_symbol_fallback: {gate_result['all_joins_symbol_fallback']}")
    return "\n".join(lines)


def _supported_drugs_default() -> list[str]:
    from dna_decode.data.mic_tiers import supported_drugs

    return list(supported_drugs())


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="genome_map", description=__doc__)
    ap.add_argument("--genome-fasta", type=Path, default=None,
                    help="genome FASTA (live/hybrid mode; runs Bakta if --gff absent)")
    ap.add_argument("--gff", type=Path, default=None,
                    help="precomputed Bakta-compatible GFF3 (skips Bakta; offline-capable)")
    ap.add_argument("--organism", default="Escherichia",
                    help="AMRFinder -O organism; 'none' = generic (no -O)")
    ap.add_argument("--sample-id", default=None, help="label for the map (default: FASTA/GFF stem)")
    ap.add_argument("--drugs", default=None,
                    help="comma-separated drugs for the overlay (default: all mic_tiers-supported)")
    ap.add_argument("--no-amrfinder", action="store_true",
                    help="skip AMRFinder (offline/degraded — tiers from GFF + determinant cells only)")
    ap.add_argument("--allow-degraded", action="store_true",
                    help="in LIVE mode, emit a tiers-only map (overlay_status=DEGRADED_USER_ACCEPTED) "
                         "if AMRFinder fails, instead of exiting BLOCKED. Off by default so a live "
                         "AMRFinder failure does not silently produce a no-determinant map.")
    ap.add_argument("--out-dir", type=Path, default=None,
                    help="output dir (default: wiki/genome_map_<sample>_<date>/)")
    ap.add_argument("--bakta-out", type=Path, default=None)
    ap.add_argument("--amrfinder-out", type=Path, default=None)
    args = ap.parse_args(argv)

    from dna_decode.genome_map import amrfinder, annotate
    from tools.docker_runner import DockerRunnerError

    if not args.genome_fasta and not args.gff:
        print("ERROR: provide --genome-fasta (live) and/or --gff (offline)", file=sys.stderr)
        return 2

    organism = None if str(args.organism).lower() == "none" else args.organism
    drugs = [d.strip() for d in args.drugs.split(",")] if args.drugs else _supported_drugs_default()
    sample_id = args.sample_id or (args.genome_fasta or args.gff).stem
    today = _date.today().isoformat()
    out_dir = args.out_dir or (REPO / "wiki" / f"genome_map_{sample_id}_{today}")
    out_dir.mkdir(parents=True, exist_ok=True)
    base = (args.genome_fasta or args.gff).parent

    # --- resolve the GFF (provided, or annotate the FASTA with Bakta) ---
    gff = args.gff
    if gff is None:
        try:
            gff = annotate.run_bakta(args.genome_fasta, args.bakta_out or base / "bakta", prefix=sample_id)
        except (DockerRunnerError, RuntimeError, OSError) as e:
            print(f"BAKTA_ANNOTATION_BLOCKED: {e}", file=sys.stderr)
            return 3

    # --- AMRFinder + overlay_status ---
    # overlay_status is an explicit top-level field so a degraded map is never
    # mistaken for a full overlay. In LIVE mode (a FASTA was supplied expecting the
    # determinant overlay) an AMRFinder failure is BLOCKED, not a silent degrade —
    # a no-determinant map + exit 0 would be a fail-open. --allow-degraded opts in.
    main_tsv = None
    degraded = False
    if args.no_amrfinder:
        degraded = True
        overlay_status = "OFFLINE_NO_AMRFINDER"      # user asked to skip (intended)
    elif args.genome_fasta is None:
        degraded = True
        overlay_status = "GFF_ONLY_NO_FASTA"          # no FASTA to scan (intended)
    else:
        try:
            main_tsv, _ = amrfinder.run_amrfinder(
                args.genome_fasta, args.amrfinder_out or base / "amrfinder", organism=organism)
            overlay_status = "FULL"
        except (DockerRunnerError, RuntimeError, OSError) as e:
            if not args.allow_degraded:
                print(f"AMRFINDER_BLOCKED: {e}\n"
                      f"(live mode expected the determinant overlay; pass --allow-degraded to emit a "
                      f"tiers-only map instead)", file=sys.stderr)
                return 3
            print(f"WARN AMRFinder failed -> DEGRADED_USER_ACCEPTED (tiers-only map): {e}", file=sys.stderr)
            degraded = True
            overlay_status = "DEGRADED_USER_ACCEPTED"

    if main_tsv is None:
        # offline / degraded: build the map with no determinants (tiers + unknown only).
        from dna_decode.genome_map import ingest
        from dna_decode.genome_map.build_map import build_genome_map

        features = ingest.load_genome_gff(gff)
        gm = build_genome_map(sample_id, organism, features, [],
                              {"n_main_rows": 0, "n_high_confidence_join": 0,
                               "n_symbol_fallback": 0, "n_unjoined": 0},
                              drug_verdicts={}, drugs=drugs, degraded=True)
    else:
        gm = run_genome_map_for(sample_id, organism, gff, main_tsv, drugs,
                                fasta_path=args.genome_fasta)
        gm["degraded_coverage"] = degraded

    gm["overlay_status"] = overlay_status
    gate_result = evaluate_gate(gm)

    (out_dir / f"genome_map_{sample_id}.json").write_text(json.dumps(gm, indent=2), encoding="utf-8")
    (out_dir / f"genome_map_{sample_id}_table.json").write_text(
        json.dumps(build_feature_table(gm), indent=2), encoding="utf-8")
    md = render_genome_summary_md(gm, gate_result, generated=today)
    (out_dir / f"genome_map_{sample_id}.md").write_text(md, encoding="utf-8")

    m = gm["metrics"]
    print(f"Wrote {out_dir} | overlay_status={gm['overlay_status']} "
          f"features={m['total_features']} "
          f"determinant-phenotype={m['determinant_phenotype_feature_count']} "
          f"unknown_under_bakta_db_light={m['unknown_under_bakta_db_light']:.3f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
