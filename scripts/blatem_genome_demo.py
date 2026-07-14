"""Real blaTEM-1 CDS end-to-end forward demo (R3 real-surface): a REAL nucleotide edit in the actual
blaTEM coding sequence -> codon -> amino-acid change -> the Regime-B ESM2 predictor -> cross-checked against
the wet-lab DMS ampicillin-fitness measurement for that exact variant.

Substrate: PZ538321.1 (E. coli class A beta-lactamase CDS, 861 nt) — VERIFIED to translate byte-identically
to the ProteinGym BLAT_ECOLX_Stiffler_2015 target protein (so the coordinate frame is real, not synthetic).

Enumerates ALL single-nucleotide substitutions of the real CDS (sense region), classifies each
silent/missense/nonsense via predict_genome_edit, and for every MISSENSE edit that realizes a DMS-measured
variant, drives it through the genome->codon->ESM2 path and correlates the ESM score with the measured
fitness. Uses the CACHED ESM table (no model run). This is the honest single-nucleotide-accessible mutation
set — the subset of the full DMS that a real point edit can actually reach.
"""
from __future__ import annotations

import csv
import json
import sys
from datetime import date as _date
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from dna_decode.forward import predict_genome_edit, translate_codon  # noqa: E402
from scripts.tem1_forward_cell import spearman  # noqa: E402

CDS_FASTA = REPO / "data" / "forward_ref" / "blatem_3349172526.fna"   # PZ538321.1
PG = Path("D:/dna_decode_cache/proteingym")
DMS_ID = "BLAT_ECOLX_Stiffler_2015"
ESM_TABLE = Path(f"D:/dna_decode_cache/esm/esm2_t33_650M_UR50D__{DMS_ID}.json")


def load_cds(f: Path) -> str:
    return "".join(ln.strip() for ln in f.read_text(encoding="utf-8").splitlines()
                   if not ln.startswith(">")).upper()


def load_target_and_dms():
    target = next(r["target_seq"] for r in csv.DictReader(open(PG / "pg_reference.csv", encoding="utf-8"))
                  if r["DMS_id"] == DMS_ID)
    dms = {}
    for r in csv.DictReader(open(PG / "pg_dms" / "DMS_ProteinGym_substitutions" / f"{DMS_ID}.csv",
                                 encoding="utf-8")):
        m = (r.get("mutant") or "").strip()
        if ":" not in m:
            try:
                dms[m] = float(r["DMS_score"])
            except (TypeError, ValueError, KeyError):
                pass
    return target, dms


def main() -> int:
    if not CDS_FASTA.exists():
        print(f"ERROR: real CDS {CDS_FASTA} missing", file=sys.stderr)
        return 2
    if not ESM_TABLE.exists():
        print(f"ERROR: cached ESM table {ESM_TABLE} missing (run tem1_forward_cell.py --method esm2 first)",
              file=sys.stderr)
        return 2
    cds = load_cds(CDS_FASTA)
    target, dms = load_target_and_dms()
    esm_table = {int(k): v for k, v in json.loads(ESM_TABLE.read_text(encoding="utf-8")).items()}

    # INTEGRITY: the real CDS must translate to the exact DMS target protein (coordinate frame is real)
    translated = "".join(translate_codon(cds[i:i + 3]) for i in range(0, len(cds) - 2, 3)).rstrip("*")
    integrity_ok = (translated == target)
    n_sense_nt = len(target) * 3   # 286 codons -> 858 nt (exclude the terminal stop codon)

    counts = {"silent": 0, "missense": 0, "nonsense": 0}
    matched = []          # (aa_mutation, esm_score, dms_score)
    esm_scores, dms_scores = [], []
    examples = {"benign": None, "damaging": None, "nonsense": None, "silent": None}

    for nt_pos in range(1, n_sense_nt + 1):
        ref = cds[nt_pos - 1]
        for alt in "ACGT":
            if alt == ref:
                continue
            ge = predict_genome_edit(cds, nt_pos, ref, alt, protein_seq=target,
                                     protein="blaTEM-1 (PZ538321.1)",
                                     phenotype_axis="ampicillin fitness (DMS)",
                                     method="esm2", esm_table=esm_table)
            counts[ge.consequence] += 1
            if ge.consequence == "silent" and examples["silent"] is None:
                examples["silent"] = {"nt_edit": f"{ref}{nt_pos}{alt}", "codon": f"{ge.wt_codon}>{ge.alt_codon}",
                                      "aa_pos": ge.aa_pos, "aa": ge.wt_aa, "consequence": "silent"}
            if ge.consequence == "nonsense" and examples["nonsense"] is None:
                examples["nonsense"] = {"nt_edit": f"{ref}{nt_pos}{alt}", "codon": f"{ge.wt_codon}>{ge.alt_codon}",
                                        "aa_mutation": ge.aa_mutation,
                                        "esm_score": round(ge.protein_prediction.raw_score, 3)}
            if ge.consequence == "missense" and ge.aa_mutation in dms:
                s = ge.protein_prediction.raw_score
                d = dms[ge.aa_mutation]
                matched.append((ge.aa_mutation, round(s, 3), round(d, 3), f"{ref}{nt_pos}{alt}",
                                f"{ge.wt_codon}>{ge.alt_codon}"))
                esm_scores.append(s)
                dms_scores.append(d)

    rho = spearman(esm_scores, dms_scores) if len(esm_scores) >= 3 else float("nan")
    # narrated examples: most-damaging (lowest DMS) + most-benign (highest DMS) matched real edits
    if matched:
        dmg = min(matched, key=lambda t: t[2])
        ben = max(matched, key=lambda t: t[2])
        examples["damaging"] = {"aa_mutation": dmg[0], "esm_score": dmg[1], "dms_score": dmg[2],
                                "nt_edit": dmg[3], "codon": dmg[4]}
        examples["benign"] = {"aa_mutation": ben[0], "esm_score": ben[1], "dms_score": ben[2],
                              "nt_edit": ben[3], "codon": ben[4]}

    res = {
        "demo": "blatem_genome_end_to_end",
        "cds_source": "PZ538321.1 (E. coli class A beta-lactamase, NCBI)",
        "cds_len_nt": len(cds),
        "translation_matches_dms_target": integrity_ok,      # MUST be True — real coordinate frame
        "dms_id": DMS_ID,
        "method": "esm2_zeroshot (cached table)",
        "single_nt_edits_enumerated": sum(counts.values()),
        "consequence_counts": counts,
        "missense_edits_matching_a_dms_variant": len(matched),
        "spearman_esm_vs_measured_dms_on_single_nt_accessible_set": round(rho, 4),
        "examples": examples,
        "honesty": ("Real blaTEM CDS (translation byte-identical to the DMS protein). Single-nucleotide-"
                    "accessible missense subset only (a point edit cannot reach codon-distant AAs — that is "
                    "the honest natural mutation set). ESM2 cached table; cross-checked vs measured DMS."),
        "status": "GENOME_PATH_REAL_SURFACE_VALIDATED" if (integrity_ok and len(matched) >= 100) else "DEGRADED",
    }
    out = REPO / "wiki" / f"blatem_genome_demo_{_date.today().isoformat()}.json"
    out.write_text(json.dumps(res, indent=2), encoding="utf-8")
    print(f"[blatem-genome-demo] CDS=PZ538321.1 len={len(cds)}nt translation_matches_target={integrity_ok}")
    print(f"  single-nt edits: {sum(counts.values())} -> silent {counts['silent']} / missense "
          f"{counts['missense']} / nonsense {counts['nonsense']}")
    print(f"  missense edits realizing a DMS variant: {len(matched)}")
    print(f"  Spearman(ESM2, measured DMS) on the single-nt-accessible set: {res['spearman_esm_vs_measured_dms_on_single_nt_accessible_set']}")
    if examples["damaging"]:
        e = examples["damaging"]; print(f"  most-damaging real edit: {e['nt_edit']} ({e['codon']}) -> {e['aa_mutation']} | ESM {e['esm_score']} / DMS {e['dms_score']}")
        b = examples["benign"];  print(f"  most-benign   real edit: {b['nt_edit']} ({b['codon']}) -> {b['aa_mutation']} | ESM {b['esm_score']} / DMS {b['dms_score']}")
    print(f"  status={res['status']}  artifact -> {out}")
    return 0 if res["status"] == "GENOME_PATH_REAL_SURFACE_VALIDATED" else 1


if __name__ == "__main__":
    raise SystemExit(main())
