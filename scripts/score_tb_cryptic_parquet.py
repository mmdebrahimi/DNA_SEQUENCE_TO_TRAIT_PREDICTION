"""Score the TB cell on the CRyPTIC Zenodo PARQUET dump (the in-distribution baseline, deliverable-a).

The shipped `scripts/score_tb_cryptic.py` consumes per-isolate masked VCFs (the ~1.6 TB regeno cohort that
was never staged). The CRyPTIC Zenodo dump (now on D:) ships the SAME determinant information as
`VARIANTS.parquet` — per-isolate, per-position rows whose `VARIANT` is the GENOMIC nucleotide form
(`761155c>t` = pos 761155, ref C, alt T). That matches the WHO catalogue's genomic `(pos, ref, alt)`
determinants exactly (rpoB S450L -> 761155 C>T), so this is a thin PARQUET->calls ADAPTER over the FROZEN
scorer: `tb_amr.score_drug` (genomic SNV match) + `tb_lineage.lineage_clusters` + `score_tb_cryptic
.score_cohort` are reused unchanged.

HONEST (unchanged): WHO_CATALOGUE_ON_CRYPTIC_KNOWLEDGE_BASELINE — the catalogue was built partly from
CRyPTIC, so this is in-distribution, NOT independent validation. No regeno here => callability is
unassessed (non-match -> S, never ABSTAIN); flagged in the result. CRyPTIC's per-position decomposition
means an MNV determinant is matched base-by-base (the documented masked-VCF FN cause is sidestepped).

Reads VARIANTS.parquet in row-batches (it is 59 M rows / 2.9 GB) and keeps ONLY rows at determinant ∪
barcode positions for major alleles — memory-bounded on a constrained host. Labels = the reuse table's
`{CODE}_BINARY_PHENOTYPE` (measured), quality HIGH.
"""
from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from datetime import date as _date
from pathlib import Path

import pyarrow.parquet as pq

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from dna_decode.data import tb_lineage_barcode, tb_who_catalogue  # noqa: E402
from dna_decode.organism_rules import tb_amr, tb_lineage, tb_vcf  # noqa: E402
from scripts.score_tb_cryptic import score_cohort  # noqa: E402

DEFAULT_DUMP = Path("D:/dna_decode_cache/data files donwload/lotsa SMILES")
DEFAULT_REUSE = Path("D:/dna_decode_cache/data files donwload/CRyPTIC TB MIC compendium/"
                     "CRyPTIC_reuse_table_20240917.csv")
# drug -> CRyPTIC reuse-table phenotype-column code. First-line (RIF/INH) shipped; the rest are the
# SECOND-LINE + new-drug extension (2026-07-14) — each key is a DRUG_CATALOGUE_NAME key AND has a
# `<code>_BINARY_PHENOTYPE` column verified present in CRyPTIC_reuse_table_20240917.csv. The scorer
# (tb_amr.score_drug) is already drug-generic; this map is the only thing that gated second-line.
DRUG_CODE = {
    "rifampicin": "RIF", "isoniazid": "INH",           # first-line (shipped)
    "ethambutol": "EMB",                                 # first-line companion
    "moxifloxacin": "MXF", "levofloxacin": "LEV",        # fluoroquinolones (gyrA/gyrB)
    "amikacin": "AMI", "kanamycin": "KAN",               # aminoglycosides (rrs/eis)
    "ethionamide": "ETH",                                # ethA/inhA
    "bedaquiline": "BDQ", "clofazimine": "CFZ",          # Rv0678/atpE/pepQ (cross-resistant)
    "delamanid": "DLM",                                  # ddn/fbiA-C/fgd1
    "linezolid": "LZD",                                  # rplC/rrl
}

_SNV_RE = re.compile(r"^(\d+)([acgt])>([acgt])$")
_INDEL_RE = re.compile(r"^(\d+)_(del|ins)_([acgt]+)$")


def parse_cryptic_variant(variant: str) -> tuple[int, str, str] | None:
    """CRyPTIC genomic SNV string `761155c>t` -> (pos, REF, ALT) uppercase; None for non-SNV (indels)."""
    m = _SNV_RE.match((variant or "").strip())
    if not m:
        return None
    return int(m.group(1)), m.group(2).upper(), m.group(3).upper()


def who_indel_to_cryptic_string(pos: int, ref: str, alt: str) -> str | None:
    """A WHO VCF-anchored indel determinant -> the EXACT CRyPTIC VARIANT string, or None if not a simple
    indel (complex delins where alt/ref share no clean prefix — needs reference left-alignment, unmapped).

    Conventions VERIFIED against CRyPTIC VARIANTS.parquet (scripts/_tb_indel_probe.py):
      deletion  (alt is prefix of ref): `{pos+len(alt)}_del_{ref[len(alt):]}`   e.g. 761100 CAATTCATGG>C -> 761101_del_aattcatgg
      insertion (ref is prefix of alt): `{pos+len(ref)-1}_ins_{alt[len(ref):]}` e.g. 761102 A>ATTC       -> 761102_ins_ttc
    EXACT match (no positional tolerance) => zero false-positive risk."""
    ref, alt = ref.upper(), alt.upper()
    if len(ref) > len(alt) and ref.startswith(alt):
        return f"{pos + len(alt)}_del_{ref[len(alt):].lower()}"
    if len(alt) > len(ref) and alt.startswith(ref):
        return f"{pos + len(ref) - 1}_ins_{alt[len(ref):].lower()}"
    return None


def indel_determinant_targets(determinants) -> dict[str, object]:
    """{exact CRyPTIC indel VARIANT string -> Determinant} for every simply-mappable indel determinant.
    Complex delins (who_indel_to_cryptic_string -> None) are skipped (reported separately as unmapped)."""
    targets: dict[str, object] = {}
    for d in determinants:
        if len(d.ref) == len(d.alt):
            continue
        s = who_indel_to_cryptic_string(d.pos, d.ref, d.alt)
        if s is not None:
            targets[s] = d
    return targets


def barcode_positions(barcode) -> set[int]:
    """Genomic positions referenced by the lineage barcode (so lineage SNPs survive the position filter)."""
    pos: set[int] = set()
    for item in barcode:
        p = getattr(item, "pos", None)
        if p is None and isinstance(item, dict):
            p = item.get("pos") or item.get("position")
        if isinstance(p, int):
            pos.add(p)
    return pos


def load_calls_by_strain(parquet_path: Path, wanted_pos: set[int], indel_targets: dict[str, object] | None = None,
                         batch_size: int = 1_000_000, progress: bool = True):
    """Stream VARIANTS.parquet once -> ({UNIQUEID: {pos: VariantCall}} for MAJOR-allele SNVs at wanted_pos,
    {UNIQUEID: set(matched indel VARIANT strings)} for MAJOR-allele indels hitting indel_targets)."""
    indel_targets = indel_targets or {}
    pf = pq.ParquetFile(str(parquet_path))
    calls: dict[str, dict[int, tb_vcf.VariantCall]] = {}
    indel_hits: dict[str, set] = {}
    seen = 0
    for batch in pf.iter_batches(batch_size=batch_size,
                                 columns=["VARIANT", "IS_MINOR_ALLELE", "UNIQUEID", "FRS"]):
        d = batch.to_pydict()
        for v, minor, uid, frs in zip(d["VARIANT"], d["IS_MINOR_ALLELE"], d["UNIQUEID"], d["FRS"]):
            if minor:                      # keep major alleles only (mirror masked PASS/GT>=1 major call)
                continue
            vs = (v or "").strip()
            parsed = parse_cryptic_variant(vs)
            if parsed is not None:
                pos, ref, alt = parsed
                if pos in wanted_pos:
                    calls.setdefault(str(uid), {})[pos] = tb_vcf.VariantCall(
                        pos=pos, ref=ref, alt=alt, gt="1/1", frs=(float(frs) if frs is not None else None))
            elif vs in indel_targets:      # exact-match an indel determinant realization
                indel_hits.setdefault(str(uid), set()).add(vs)
        seen += batch.num_rows
        if progress:
            print(f"    ...scanned {seen:,} variant rows; {len(calls):,} isolates with SNV calls; "
                  f"{len(indel_hits):,} with indel-determinant hits")
    return calls, indel_hits


def load_labels(reuse_csv: Path, code: str) -> dict[str, str]:
    """{UNIQUEID: 'R'/'S'} from the reuse table, measured BINARY_PHENOTYPE with quality HIGH."""
    labels: dict[str, str] = {}
    with open(reuse_csv, encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            ph = (r.get(f"{code}_BINARY_PHENOTYPE") or "").strip().upper()
            q = (r.get(f"{code}_PHENOTYPE_QUALITY") or "").strip().upper()
            uid = (r.get("UNIQUEID") or "").strip()
            if uid and ph in ("R", "S") and q == "HIGH":
                labels[uid] = ph
    return labels


def run(drug: str, dump_dir: Path, reuse_csv: Path, max_isolates: int = 0, match_indels: bool = True) -> dict:
    code = DRUG_CODE[drug]
    tb_who_catalogue.verify_pins()
    dets = tb_who_catalogue.load_determinants(drug)
    barcode = tb_lineage_barcode.load_barcode()
    det_pos = {d.pos for d in dets}
    wanted = det_pos | barcode_positions(barcode)
    indel_dets = [d for d in dets if len(d.ref) != len(d.alt)]
    targets = indel_determinant_targets(dets) if match_indels else {}
    n_complex = len(indel_dets) - len({who_indel_to_cryptic_string(d.pos, d.ref, d.alt) for d in indel_dets
                                       if who_indel_to_cryptic_string(d.pos, d.ref, d.alt)})
    print(f"[tb-cryptic-parquet] {drug}: {len(dets)} determinant rows ({len(det_pos)} positions) + "
          f"{len(wanted) - len(det_pos)} barcode positions")
    print(f"[tb-cryptic-parquet] indel determinants: {len(indel_dets)} -> {len(targets)} exact CRyPTIC "
          f"targets ({n_complex} complex delins unmapped); indel matching={'ON' if match_indels else 'OFF'}")

    labels = load_labels(reuse_csv, code)
    print(f"[tb-cryptic-parquet] {len(labels):,} measured HIGH-quality {code} labels")
    calls, indel_hits = load_calls_by_strain(dump_dir / "VARIANTS.parquet", wanted, targets)

    # Restrict to isolates that HAVE a measured label; an absent isolate in `calls` = no in-scope variants
    # (all-reference at determinant+barcode positions) -> empty calls (susceptible-by-absence, callability
    # unassessed). Cap for a smoke subset if requested.
    label_uids = list(labels)
    if max_isolates:
        label_uids = label_uids[:max_isolates]
    cohort_calls = {uid: calls.get(uid, {}) for uid in label_uids}
    cohort_labels = {uid: labels[uid] for uid in label_uids}

    # SNV-only prediction (the prior lower bound), then OR in exact indel-determinant hits.
    snv_preds = {uid: tb_amr.score_drug(drug, c, dets).prediction for uid, c in cohort_calls.items()}
    indel_R = {uid for uid in label_uids if indel_hits.get(uid)}
    preds = {uid: ("R" if (snv_preds[uid] == "R" or uid in indel_R) else snv_preds[uid]) for uid in label_uids}
    # audit the marginal effect: isolates flipped S->R purely by an indel determinant, split by truth
    flips = [uid for uid in label_uids if snv_preds[uid] == "S" and uid in indel_R]
    flip_tp = sum(1 for uid in flips if cohort_labels[uid] == "R")
    flip_fp = sum(1 for uid in flips if cohort_labels[uid] == "S")

    clusters = tb_lineage.lineage_clusters(cohort_calls, barcode)
    res = score_cohort(preds, cohort_labels, clusters, drug=drug,
                       cohort_complete=(max_isolates == 0))
    res["source"] = "CRyPTIC Zenodo VARIANTS.parquet (genomic SNV + exact indel-determinant match)"
    res["callability_assessed"] = False
    res["n_with_inscope_calls"] = sum(1 for c in cohort_calls.values() if c)
    res["indel_matching"] = {
        "enabled": match_indels,
        "indel_determinants": len(indel_dets),
        "exact_cryptic_targets": len(targets),
        "complex_delins_unmapped": n_complex,
        "isolates_with_indel_hit": sum(1 for uid in label_uids if indel_hits.get(uid)),
        "snv_to_R_flips": len(flips),
        "flips_true_positive": flip_tp,
        "flips_false_positive": flip_fp,
    }
    # also report the SNV-only raw confusion for a clean before/after delta
    res["raw_snv_only"] = score_cohort(snv_preds, cohort_labels, clusters, drug=drug,
                                       cohort_complete=(max_isolates == 0))["raw"]
    return res


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--drug", default="rifampicin", choices=list(DRUG_CODE))
    ap.add_argument("--dump-dir", type=Path, default=DEFAULT_DUMP)
    ap.add_argument("--reuse-csv", type=Path, default=DEFAULT_REUSE)
    ap.add_argument("--max", type=int, default=0, help="cap isolates (0=all; >0 => TB_SUBSET_PLUMBING)")
    ap.add_argument("--no-indels", action="store_true", help="SNV-only (the prior lower-bound behavior)")
    a = ap.parse_args(argv)
    if not (a.dump_dir / "VARIANTS.parquet").exists():
        print(f"ERROR: VARIANTS.parquet not found under {a.dump_dir}", file=sys.stderr)
        return 2
    if not a.reuse_csv.exists():
        print(f"ERROR: reuse table not found at {a.reuse_csv}", file=sys.stderr)
        return 2
    res = run(a.drug, a.dump_dir, a.reuse_csv, a.max, match_indels=not a.no_indels)
    code = DRUG_CODE[a.drug]
    out = REPO / "wiki" / f"tb_{code.lower()}_cryptic_parquet_baseline_{_date.today().isoformat()}.json"
    out.write_text(json.dumps(res, indent=2, default=str), encoding="utf-8")
    print(f"status={res['status']}  n_isolates={res['n_isolates']}  "
          f"n_with_inscope_calls={res.get('n_with_inscope_calls')}")
    im = res.get("indel_matching", {})
    print(f"  indel matching: {im.get('isolates_with_indel_hit')} isolates hit | "
          f"S->R flips={im.get('snv_to_R_flips')} (TP={im.get('flips_true_positive')} "
          f"FP={im.get('flips_false_positive')}) | complex-unmapped={im.get('complex_delins_unmapped')}")
    print(f"  raw SNV-only : sens={res['raw_snv_only']['sens']} spec={res['raw_snv_only']['spec']}")
    print(f"  raw SNV+indel: sens={res['raw']['sens']} spec={res['raw']['spec']}")
    lc = res.get("lineage_collapsed")
    if lc:
        print(f"  lineage-collapsed: sens={lc['sens']} spec={lc['spec']} "
              f"(R-lineages={lc['n_clusters_R']} S-lineages={lc['n_clusters_S']} discordant={lc['n_discordant']})")
    print(f"artifact -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
