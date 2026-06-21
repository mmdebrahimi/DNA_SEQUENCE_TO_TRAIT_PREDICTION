"""Genome-map prototype spike (3 bacterial genomes) + GO/NO-GO verdict (Step 7).

Runs the full pipeline on 3 BACTERIAL prototype genomes (brainstorm catch M1 —
TB is out of the v1 spike):

  annotate (Bakta) + AMRFinder [Step 1] -> manifest [Step 2] -> tier [Step 3]
  -> determinant overlay [Step 4] -> assemble map [Step 5] -> gate [Step 6]

and emits wiki/genome_map_spike_verdict_<date>.md (+ per-genome JSON maps): per
genome the tier counts + DB-labelled unknown rate + join-quality counts + the
G1 prevent-wrong-inference evidence + the G2 spot-check, then the cross-genome
GO/NO-GO verdict. An honest NO-GO / BLOCKED is a valid outcome (a Bakta wedge
-> BAKTA_ANNOTATION_BLOCKED, never a fake map).

The 3 prototype genomes (Open Question A, resolved 2026-06-18):
  (i)   GCA_002180195.1  E. coli ST131, cipro-R (rich AMR overlay)  -O Escherichia
  (ii)  GCA_000417105.1  K. pneumoniae, meropenem/carbapenem-R       -O Klebsiella_pneumoniae
  (iii) Gemmata obscuriglobus (Planctomycete; homology/hypothetical-rich
        honesty stress test)                                          no -O (generic)

The pure aggregation / rendering functions are unit-tested on synthetic inputs;
the live 3-genome run is the manual deliverable.

Usage:
    uv run python -m scripts.genome_map_spike --refseq-cache D:/dna_decode_cache/refseq
Exit: 0 = GO, 1 = NO_GO, 2 = BLOCKED, 3 = usage/IO error.
"""
from __future__ import annotations

import argparse
import json
import sys
from dataclasses import dataclass, field
from datetime import date as _date
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from dna_decode.genome_map import ingest
from dna_decode.genome_map.build_map import build_feature_table, build_genome_map
from dna_decode.genome_map.gate import (
    VERDICT_GO,
    VERDICT_NO_GO,
    aggregate_spike_verdict,
    evaluate_gate,
)
from dna_decode.genome_map.phenotype_overlay import (
    build_contig_name_map,
    join_hits,
    parse_determinant_hits,
)
from dna_decode.genome_map.virulence_overlay import (
    genome_pathotype_call,
    join_virulence,
    parse_virulence_hits,
    virulence_organism_in_scope,
)
from dna_decode.pathotype.detect import parse_fasta
from dna_decode.pathotype.vf_runner import run_canonical_vf

VERDICT_BLOCKED = "BLOCKED"

# The committed v1 VirulenceFinder E. coli allele DB (E. coli / Shigella scope).
VF_DB_DEFAULT = REPO / "data" / "virulencefinder_db" / "virulence_ecoli.fsa"


@dataclass
class GenomeConfig:
    accession: str
    organism: str | None      # AMRFinder -O; None = generic (no -O)
    label: str
    drugs: list[str] = field(default_factory=list)


# The 3 prototype genomes. Gemmata obscuriglobus type-strain accession resolved
# at run time via NCBI Datasets; the homology-heavy honesty stress test.
DEFAULT_GENOMES = [
    GenomeConfig("GCA_002180195.1", "Escherichia", "E. coli ST131 (cipro-R)",
                 drugs=["ciprofloxacin", "ceftriaxone", "tetracycline", "gentamicin"]),
    GenomeConfig("GCA_000417105.1", "Klebsiella_pneumoniae", "K. pneumoniae (carbapenem-R)",
                 drugs=["meropenem", "ciprofloxacin", "ceftriaxone", "gentamicin", "tetracycline"]),
    GenomeConfig("GCF_000171775.1", None, "Gemmata obscuriglobus (homology stress test)",
                 drugs=[]),
]


# ---------- contig-length helpers (for the Step-4 coordinate join) ----------

def fasta_contig_lengths(fasta_path: Path | str) -> dict[str, int]:
    """Map FASTA contig name (header up to first whitespace) -> sequence length."""
    out: dict[str, int] = {}
    name = None
    n = 0
    with open(fasta_path, encoding="utf-8") as f:
        for line in f:
            if line.startswith(">"):
                if name is not None:
                    out[name] = n
                name = line[1:].split()[0]
                n = 0
            else:
                n += len(line.strip())
    if name is not None:
        out[name] = n
    return out


def bakta_contig_lengths(gff_path: Path | str) -> dict[str, int]:
    """Map Bakta seqid -> length from the GFF ``##sequence-region`` pragmas."""
    out: dict[str, int] = {}
    with open(gff_path, encoding="utf-8") as f:
        for line in f:
            if line.startswith("##sequence-region"):
                parts = line.split()
                # ##sequence-region <seqid> <start> <end>
                if len(parts) >= 4:
                    try:
                        out[parts[1]] = int(parts[3]) - int(parts[2]) + 1
                    except ValueError:
                        continue
            elif line.startswith("##FASTA"):
                break
    return out


# ---------- per-genome map (file IO, no Docker) ----------

def _virulence_for_genome(organism, fasta_path, features, contig_map, vf_db, *, enabled):
    """Run the VF overlay (Step 5) → (status, joined, counts, pathotype_call, db_sha).

    status ∈ {FULL, UNAVAILABLE_NO_BLASTN, SKIPPED_NON_ECOLI, SKIPPED_USER}. Reuses the
    SAME contig-name map already built for the AMR join. Fail-soft: a non-ok VF result
    (no blastn / blastn error / no FASTA to scan) yields UNAVAILABLE_NO_BLASTN with no
    tier — the map still emits.
    """
    if not enabled:
        return "SKIPPED_USER", None, None, None, None
    if not virulence_organism_in_scope(organism):
        return "SKIPPED_NON_ECOLI", None, None, None, None
    if fasta_path is None or not Path(fasta_path).exists():
        return "UNAVAILABLE_NO_BLASTN", None, None, None, None   # nothing to scan
    vf_db = vf_db or VF_DB_DEFAULT
    result = run_canonical_vf(str(fasta_path), str(vf_db), all_hits=True)
    if result.get("status") != "ok":
        return "UNAVAILABLE_NO_BLASTN", None, None, None, result.get("db_sha")
    contigs = parse_fasta(fasta_path)
    contig_names = [name for name, _ in contigs]
    vir_hits = parse_virulence_hits(result)
    joined, counts = join_virulence(features, vir_hits, contig_name_map=contig_map,
                                    contig_names=contig_names)
    pathotype_call = genome_pathotype_call(result, contigs)
    return "FULL", joined, counts, pathotype_call, result.get("db_sha")


def run_genome_map_for(
    accession: str,
    organism: str | None,
    gff_path: Path | str,
    main_tsv_path: Path | str,
    drugs: list[str],
    fasta_path: Path | str | None = None,
    *,
    virulence: bool = False,
    vf_db: Path | str | None = None,
) -> dict:
    """Build the genome map from a Bakta GFF + an AMRFinder main.tsv (no live tools).

    Reconciles AMRFinder contig names to Bakta seqids by length (when a FASTA is
    given) so the coordinate join can fire; computes per-drug genome-level
    verdicts via the frozen ``amr_rules.call_resistance`` (organism-aware).

    `virulence=True` (the single-genome product surface) additionally runs the VF
    overlay when the organism is in scope + a FASTA is present + blastn resolves —
    reusing the SAME contig-name map — and stamps `gm["virulence_status"]`. Default
    False keeps the AMR-only spike + direct callers byte-unchanged.
    """
    from dna_decode.eval.amr_rules import call_resistance

    features = ingest.load_genome_gff(gff_path)
    hits = parse_determinant_hits(main_tsv_path)

    contig_map = None
    if fasta_path is not None and Path(fasta_path).exists():
        contig_map = build_contig_name_map(
            fasta_contig_lengths(fasta_path), bakta_contig_lengths(gff_path)
        )
    joined, counts = join_hits(features, hits, contig_name_map=contig_map)

    drug_verdicts = {}
    for d in drugs:
        drug_verdicts[d] = call_resistance(Path(main_tsv_path), d, organism=organism)

    vstatus, vir_joined, vir_counts, pathotype_call, vir_db_sha = _virulence_for_genome(
        organism, fasta_path, features, contig_map, vf_db, enabled=virulence)

    gm = build_genome_map(
        accession, organism, features, joined, counts,
        drug_verdicts=drug_verdicts, drugs=drugs,
        virulence_joined_hits=vir_joined, virulence_join_counts=vir_counts,
        pathotype_call=pathotype_call, virulence_db_sha=vir_db_sha,
    )
    gm["virulence_status"] = vstatus
    return gm


# ---------- aggregation + rendering (pure) ----------

def summarize_spike(entries: list[dict]) -> dict:
    """Aggregate per-genome entries into the spike verdict.

    Each entry: {accession, status, genome_map?, gate_result?}. If ANY genome is
    not OK (BLOCKED), the spike verdict is BLOCKED (honest — the spike could not
    complete). Otherwise apply the cross-genome G1/G2/all-symbol-fallback rule.
    """
    blocked = [e for e in entries if e.get("status") != "OK"]
    if blocked:
        return {
            "verdict": VERDICT_BLOCKED,
            "reasons": [f"{e['accession']}: {e.get('status')}" for e in blocked],
            "spike_g1_pass": None, "spike_g2_pass": None, "any_all_symbol_fallback": None,
        }
    gate_results = [e["gate_result"] for e in entries]
    return aggregate_spike_verdict(gate_results)


def render_verdict_md(entries: list[dict], aggregate: dict, *, generated: str | None = None) -> str:
    """Render the GO/NO-GO verdict markdown from per-genome entries + the aggregate."""
    generated = generated or _date.today().isoformat()
    lines: list[str] = []
    lines.append(f"# Genome-map v1 spike — GO/NO-GO verdict ({generated})")
    lines.append("")
    lines.append(f"**Tiering verdict (Bakta honesty re-tiering): {aggregate['verdict']}**")
    if "overlay_go" in aggregate:
        og = aggregate["overlay_go"]
        lines.append(f"**Overlay-integrity verdict (determinant->feature join): "
                     f"{'GO' if og else 'NOT DEMONSTRATED'}**")
    lines.append("")
    for r in aggregate.get("reasons", []):
        lines.append(f"- {r}")
    if aggregate.get("overlay_reason"):
        lines.append(f"- {aggregate['overlay_reason']}")
    lines.append("")
    lines.append("Honesty contract: phenotype claims appear ONLY behind a high-confidence "
                 "determinant join (symbol-fallback excluded); the unknown rate is DB-labelled "
                 "`unknown_under_bakta_db_light` (db-light = reduced functional coverage, not "
                 "biological unknown).")
    lines.append("")

    for e in entries:
        acc = e["accession"]
        lines.append(f"## {acc} — {e.get('label', '')}")
        if e.get("status") != "OK":
            lines.append(f"- **STATUS: {e.get('status')}** (no map emitted — honest BLOCKED, not a fake map)")
            lines.append("")
            continue
        gm = e["genome_map"]
        gr = e["gate_result"]
        m = gm["metrics"]
        lines.append(f"- organism (-O): `{gm.get('amrfinder_organism')}`")
        lines.append(f"- total features: {m['total_features']}")
        tier_str = ", ".join(f"{k}={v}" for k, v in m["per_tier_counts"].items())
        lines.append(f"- per-tier counts: {tier_str}")
        lines.append(f"- `unknown_under_bakta_db_light`: {m['unknown_under_bakta_db_light']:.3f}")
        jq = m["join_quality"]
        lines.append(f"- join quality: n_main_rows={jq.get('n_main_rows')} "
                     f"high_confidence={jq.get('n_high_confidence_join')} "
                     f"symbol_fallback={jq.get('n_symbol_fallback')} unjoined={jq.get('n_unjoined')}")
        lines.append(f"- all_joins_symbol_fallback: {m['all_joins_symbol_fallback']}")
        lines.append(f"- determinant-phenotype features: {m['determinant_phenotype_feature_count']}")
        gl = ", ".join(f"{d}={v.get('prediction')}" for d, v in m["genome_level_calls"].items())
        lines.append(f"- genome-level R/S calls (separate from features): {gl or '(none)'}")
        lines.append(f"- **G1** prevent-wrong-inference: {len(gr['g1_features'])} features "
                     f"(demote={gr['g1_demote_count']}, surface={gr['g1_surface_count']}) "
                     f"-> g1_pass={gr['g1_pass']}")
        # show up to 5 G1 examples
        for g1f in gr["g1_features"][:5]:
            if g1f["type"] == "demote-homology":
                lines.append(f"    - DEMOTE: raw `{g1f['raw_product']}` -> homology-only-hypothesis "
                             f"({g1f['classification_reason']})")
            else:
                drugs = sorted({(p.get('drug') or p.get('phenotype') or '') for p in g1f.get('phenotype', [])})
                lines.append(f"    - SURFACE: raw `{g1f['raw_product']}` "
                             f"(gene `{g1f.get('raw_gene_symbol')}`) -> determinant-phenotype "
                             f"[{', '.join(s for s in drugs if s)}]")
        lines.append(f"- **G2** no-tier-confusion: pass={gr['g2_spotcheck']['pass']} "
                     f"(violations={len(gr['g2_spotcheck']['violations'])})")
        og = gr.get("overlay_go")
        og_str = "n/a (determinant-free)" if og is None else ("GO" if og else "NOT DEMONSTRATED")
        oe = gr.get("overlay_evidence", {})
        lines.append(f"- **overlay-integrity**: {og_str} "
                     f"(main_rows={oe.get('n_main_rows')}, high_conf_joins={oe.get('n_high_confidence_join')}, "
                     f"surfaced={oe.get('surface_determinant_count')})")
        lines.append(f"- per-genome tiering verdict: {gr['verdict']}")
        lines.append("")
    return "\n".join(lines)


# ---------- live driver ----------

def _resolve_fasta(accession: str, refseq_cache: str) -> Path:
    from dna_decode.data.refseq import download_genome

    cache_dir = download_genome(accession, refseq_cache)
    return Path(cache_dir) / "genome.fna"


def run_spike_live(genomes: list[GenomeConfig], refseq_cache: str, out_dir: Path) -> list[dict]:
    """Run the live pipeline per genome; return entries for summarize/render."""
    from dna_decode.genome_map import amrfinder, annotate
    from tools.docker_runner import DockerRunnerError

    entries: list[dict] = []
    for g in genomes:
        entry = {"accession": g.accession, "label": g.label, "status": "OK"}
        try:
            fasta = _resolve_fasta(g.accession, refseq_cache)
        except Exception as e:  # noqa: BLE001
            entry["status"] = f"GENOME_DOWNLOAD_FAILED: {e}"
            entries.append(entry)
            continue

        base = Path(fasta).parent
        # Bakta
        try:
            gff = annotate.run_bakta(fasta, base / "bakta", prefix=g.accession)
        except (DockerRunnerError, RuntimeError, OSError) as e:
            entry["status"] = "BAKTA_ANNOTATION_BLOCKED"
            entry["detail"] = str(e)
            entries.append(entry)
            continue
        # AMRFinder
        try:
            main_tsv, _ = amrfinder.run_amrfinder(fasta, base / "amrfinder", organism=g.organism)
        except (DockerRunnerError, RuntimeError, OSError) as e:
            entry["status"] = "AMRFINDER_BLOCKED"
            entry["detail"] = str(e)
            entries.append(entry)
            continue

        gm = run_genome_map_for(g.accession, g.organism, gff, main_tsv, g.drugs, fasta_path=fasta)
        entry["genome_map"] = gm
        entry["gate_result"] = evaluate_gate(gm)
        # persist per-genome JSON map + flat table
        out_dir.mkdir(parents=True, exist_ok=True)
        (out_dir / f"genome_map_{g.accession}.json").write_text(json.dumps(gm, indent=2), encoding="utf-8")
        (out_dir / f"genome_map_{g.accession}_table.json").write_text(
            json.dumps(build_feature_table(gm), indent=2), encoding="utf-8")
        entries.append(entry)
    return entries


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="genome_map_spike", description=__doc__)
    ap.add_argument("--refseq-cache", default="D:/dna_decode_cache/refseq")
    ap.add_argument("--out-md", type=Path, default=None,
                    help="verdict markdown (default: wiki/genome_map_spike_verdict_<date>.md)")
    ap.add_argument("--maps-out", type=Path, default=None,
                    help="per-genome JSON map dir (default: wiki/genome_map_spike_<date>/)")
    args = ap.parse_args(argv)

    today = _date.today().isoformat()
    out_md = args.out_md or (REPO / "wiki" / f"genome_map_spike_verdict_{today}.md")
    maps_out = args.maps_out or (REPO / "wiki" / f"genome_map_spike_{today}")

    entries = run_spike_live(DEFAULT_GENOMES, args.refseq_cache, maps_out)
    aggregate = summarize_spike(entries)
    md = render_verdict_md(entries, aggregate, generated=today)
    out_md.parent.mkdir(parents=True, exist_ok=True)
    out_md.write_text(md, encoding="utf-8")
    # also drop a machine-readable verdict sidecar
    out_md.with_suffix(".json").write_text(
        json.dumps({"verdict": aggregate, "generated": today,
                    "genomes": [{"accession": e["accession"], "status": e["status"]} for e in entries]},
                   indent=2), encoding="utf-8")
    print(f"Wrote {out_md} (verdict={aggregate['verdict']})")

    return {VERDICT_GO: 0, VERDICT_NO_GO: 1, VERDICT_BLOCKED: 2}.get(aggregate["verdict"], 3)


if __name__ == "__main__":
    raise SystemExit(main())
