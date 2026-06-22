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
DRUG_CODE = {"rifampicin": "RIF", "isoniazid": "INH"}

_SNV_RE = re.compile(r"^(\d+)([acgt])>([acgt])$")


def parse_cryptic_variant(variant: str) -> tuple[int, str, str] | None:
    """CRyPTIC genomic SNV string `761155c>t` -> (pos, REF, ALT) uppercase; None for non-SNV (indels)."""
    m = _SNV_RE.match((variant or "").strip())
    if not m:
        return None
    return int(m.group(1)), m.group(2).upper(), m.group(3).upper()


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


def load_calls_by_strain(parquet_path: Path, wanted_pos: set[int], batch_size: int = 1_000_000,
                         progress: bool = True) -> dict[str, dict[int, tb_vcf.VariantCall]]:
    """Stream VARIANTS.parquet -> {UNIQUEID: {pos: VariantCall}} for MAJOR-allele SNVs at wanted_pos."""
    pf = pq.ParquetFile(str(parquet_path))
    calls: dict[str, dict[int, tb_vcf.VariantCall]] = {}
    seen = 0
    for batch in pf.iter_batches(batch_size=batch_size,
                                 columns=["VARIANT", "IS_MINOR_ALLELE", "UNIQUEID", "FRS"]):
        d = batch.to_pydict()
        for v, minor, uid, frs in zip(d["VARIANT"], d["IS_MINOR_ALLELE"], d["UNIQUEID"], d["FRS"]):
            if minor:                      # keep major alleles only (mirror masked PASS/GT>=1 major call)
                continue
            parsed = parse_cryptic_variant(v)
            if parsed is None or parsed[0] not in wanted_pos:
                continue
            pos, ref, alt = parsed
            calls.setdefault(str(uid), {})[pos] = tb_vcf.VariantCall(
                pos=pos, ref=ref, alt=alt, gt="1/1", frs=(float(frs) if frs is not None else None))
        seen += batch.num_rows
        if progress:
            print(f"    ...scanned {seen:,} variant rows; {len(calls):,} isolates with in-scope calls")
    return calls


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


def run(drug: str, dump_dir: Path, reuse_csv: Path, max_isolates: int = 0) -> dict:
    code = DRUG_CODE[drug]
    tb_who_catalogue.verify_pins()
    dets = tb_who_catalogue.load_determinants(drug)
    barcode = tb_lineage_barcode.load_barcode()
    det_pos = {d.pos for d in dets}
    wanted = det_pos | barcode_positions(barcode)
    print(f"[tb-cryptic-parquet] {drug}: {len(dets)} determinant rows ({len(det_pos)} positions) + "
          f"{len(wanted) - len(det_pos)} barcode positions")

    labels = load_labels(reuse_csv, code)
    print(f"[tb-cryptic-parquet] {len(labels):,} measured HIGH-quality {code} labels")
    calls = load_calls_by_strain(dump_dir / "VARIANTS.parquet", wanted)

    # Restrict to isolates that HAVE a measured label; an absent isolate in `calls` = no in-scope variants
    # (all-reference at determinant+barcode positions) -> empty calls (susceptible-by-absence, callability
    # unassessed). Cap for a smoke subset if requested.
    label_uids = list(labels)
    if max_isolates:
        label_uids = label_uids[:max_isolates]
    cohort_calls = {uid: calls.get(uid, {}) for uid in label_uids}
    cohort_labels = {uid: labels[uid] for uid in label_uids}

    preds = {uid: tb_amr.score_drug(drug, c, dets).prediction for uid, c in cohort_calls.items()}
    clusters = tb_lineage.lineage_clusters(cohort_calls, barcode)
    res = score_cohort(preds, cohort_labels, clusters, drug=drug,
                       cohort_complete=(max_isolates == 0))
    res["source"] = "CRyPTIC Zenodo VARIANTS.parquet (genomic-nucleotide variant match)"
    res["callability_assessed"] = False
    res["n_with_inscope_calls"] = sum(1 for c in cohort_calls.values() if c)
    return res


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--drug", default="rifampicin", choices=list(DRUG_CODE))
    ap.add_argument("--dump-dir", type=Path, default=DEFAULT_DUMP)
    ap.add_argument("--reuse-csv", type=Path, default=DEFAULT_REUSE)
    ap.add_argument("--max", type=int, default=0, help="cap isolates (0=all; >0 => TB_SUBSET_PLUMBING)")
    a = ap.parse_args(argv)
    if not (a.dump_dir / "VARIANTS.parquet").exists():
        print(f"ERROR: VARIANTS.parquet not found under {a.dump_dir}", file=sys.stderr)
        return 2
    if not a.reuse_csv.exists():
        print(f"ERROR: reuse table not found at {a.reuse_csv}", file=sys.stderr)
        return 2
    res = run(a.drug, a.dump_dir, a.reuse_csv, a.max)
    code = DRUG_CODE[a.drug]
    out = REPO / "wiki" / f"tb_{code.lower()}_cryptic_parquet_baseline_{_date.today().isoformat()}.json"
    out.write_text(json.dumps(res, indent=2, default=str), encoding="utf-8")
    print(f"status={res['status']}  n_isolates={res['n_isolates']}  "
          f"n_with_inscope_calls={res.get('n_with_inscope_calls')}")
    lc = res.get("lineage_collapsed")
    if lc:
        print(f"  lineage-collapsed: sens={lc['sens']} spec={lc['spec']} "
              f"(R-lineages={lc['n_clusters_R']} S-lineages={lc['n_clusters_S']} discordant={lc['n_discordant']})")
    print(f"  raw: sens={res['raw']['sens']} spec={res['raw']['spec']}")
    print(f"artifact -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
