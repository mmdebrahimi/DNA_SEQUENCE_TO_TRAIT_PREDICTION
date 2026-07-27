"""Build the LOCAL Klebsiella depolymerase->KL-type reference from DpoTropiSearch (fetch-only; user-side).

This helper fetches the DpoTropiSearch training data from Zenodo, extracts the depolymerase enzymatic
domain + its (prophage-LCA-inferred) capsule KL-type, dedupes + caps per KL-type, and writes a compact
local reference the Klebsiella cell reads. **The dna-decode package ships NONE of this data** — you fetch
it here, under the license, onto your own machine.

LICENSE (READ — you must decide your use complies):
  Source: DpoTropiSearch (Concha-Eloko et al., Nat Commun 2025; github.com/conchaeloko/DpoTropiSearch;
  data at Zenodo 10.5281/zenodo.14065540). The Zenodo record declares the data CC-BY-4.0 (open, attribution)
  BUT the code repository carries a "Decapsulate Non-Commercial License v1.1" that explicitly restricts
  "the Software AND Data" to Non-Commercial Use. These conflict. This helper does NOT resolve the conflict
  for you — it prints the notice and requires --accept-license so the choice is yours, informed.

Usage:
  uv run --with biopython python scripts/fetch_dpotropisearch.py --accept-license \
    --out data/kleb_ref/depolymerase_kltype_reference.faa
"""
from __future__ import annotations

import argparse
import subprocess
import sys
import zipfile
from collections import defaultdict
from pathlib import Path

ZENODO_TRAINING = "https://zenodo.org/api/records/14065540/files/Training_data.zip/content"
_NOTICE = (
    "\n=== DpoTropiSearch data license notice ===\n"
    "Source: Concha-Eloko et al., Nat Commun 2025; Zenodo 10.5281/zenodo.14065540.\n"
    "Zenodo record: CC-BY-4.0 (open, attribution). Repo LICENSE: Decapsulate Non-Commercial License v1.1\n"
    "(restricts 'Software AND Data' to Non-Commercial Use). These CONFLICT. You must decide whether YOUR use\n"
    "complies. dna-decode ships none of this data; you are fetching it here onto your own machine.\n"
    "Re-run with --accept-license once you have made that determination.\n"
)


def build_reference(tsv_path: Path, out: Path, min_members: int, cap: int) -> tuple[int, int]:
    """Extract domain_seq -> KL_type_LCA from the DpoTropiSearch training TSV; dedupe + cap per KL-type."""
    import csv
    by_kl: dict[str, list[str]] = defaultdict(list)
    with open(tsv_path, encoding="utf-8", errors="replace") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            kl = (r.get("KL_type_LCA") or "").strip()
            seq = (r.get("domain_seq") or "").strip()
            if not seq or not kl or ";" in kl or "," in kl or kl.upper() in ("NAN", "NA", "NONE", ""):
                continue
            by_kl[kl].append(seq)
    rows = []
    for kl, seqs in by_kl.items():
        uniq = list(dict.fromkeys(seqs))
        if len(uniq) >= min_members:
            rows.extend((kl, s) for s in uniq[:cap])
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", encoding="utf-8") as fh:
        fh.write("; Klebsiella depolymerase-domain -> capsule KL-type reference (user-fetched; NOT shipped).\n")
        fh.write("; Source: DpoTropiSearch, Zenodo 10.5281/zenodo.14065540 (CC-BY-4.0 record / repo "
                 "Decapsulate Non-Commercial License v1.1 — verify your use). Header: >{KL_type}|{i}\n")
        for i, (kl, s) in enumerate(sorted(rows)):
            fh.write(f">{kl}|{i}\n{s}\n")
    return len(rows), len({kl for kl, _ in rows})


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--accept-license", action="store_true",
                    help="affirm you have reviewed the DpoTropiSearch license notice + your use complies")
    ap.add_argument("--out", default="data/kleb_ref/depolymerase_kltype_reference.faa")
    ap.add_argument("--work-dir", default="D:/dna_decode_cache/kleb" if Path("D:/").exists() else ".kleb_fetch",
                    help="scratch dir for the ~56MB Zenodo download + 555MB extract (route off C: if tight)")
    ap.add_argument("--min-members", type=int, default=3)
    ap.add_argument("--cap", type=int, default=8)
    args = ap.parse_args()

    if not args.accept_license:
        print(_NOTICE)
        return 2

    work = Path(args.work_dir); work.mkdir(parents=True, exist_ok=True)
    tsv = next(iter(work.glob("*final_df*.tsv")), None) or (work / "TropiGATv2.final_df_v2.tsv")
    if not tsv.exists():
        zip_path = work / "Training_data.zip"
        if not zip_path.exists():
            print(f"fetching Training_data.zip from Zenodo -> {zip_path} (~56MB) ...")
            r = subprocess.run(["curl", "-s", "-m", "600", "-L", ZENODO_TRAINING, "-o", str(zip_path)],
                               capture_output=True, text=True, timeout=650)
            if r.returncode != 0 or not zip_path.exists() or zip_path.stat().st_size < 1_000_000:
                print(f"fetch failed (rc={r.returncode}); {r.stderr[:200]}", file=sys.stderr); return 1
        print(f"extracting -> {work} (555MB TSV; ensure the drive has room) ...")
        with zipfile.ZipFile(zip_path) as z:
            z.extractall(work)
        tsv = next(iter(work.glob("*final_df*.tsv")), None) or tsv
    if not tsv or not tsv.exists():
        print("could not locate the training TSV after extract", file=sys.stderr); return 1

    n, n_kl = build_reference(tsv, Path(args.out), args.min_members, args.cap)
    print(f"built reference: {n} depolymerase domains / {n_kl} KL-types -> {args.out}")
    print("the Klebsiella cell will now find it via resolve_reference(); "
          "or set DPO_KLEB_REFERENCE to its path.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
