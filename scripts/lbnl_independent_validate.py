"""Independent validation of the phage receptor cell on the LBNL/Arkin-Mutalik dataset.

The genuine INDEPENDENT test the phage cell needed (rows 559-560 were blocked on it): score the
shipped BASEL-2021 genome-homology caller on a DIFFERENT-LAB, MEASURED-receptor phage set — the
Moriniere/Noonan/Arkin/Mutalik "Phage Datasheets" (github.com/mjohnson11/PhageDataSheets;
Table_S1_Phages.tsv = per-phage receptor measured by genome-wide genetic screens on E. coli K-12
BW25113). BW25113 is a K-12 derivative, so the receptor labels are comparable to BASEL's K-12 host
(unlike the 2025 O-antigen-restored set).

LEAKAGE GUARD: the LBNL set INCLUDES the BASEL collection (Bas## phages) — the exact phages the catalogue
is built from. Those are EXCLUDED; only the non-Bas LBNL isolates are scored → genuinely disjoint from
the BASEL-2021 reference.

Two honest numbers are reported:
  - overall accuracy on ALL scoreable non-Bas phages (shows the v0 catalogue's CLASS COVERAGE limit —
    it covers ~5 receptor classes; the RBP-variable classes Tsx/OmpC/FhuA/... are out of v0 scope);
  - covered-subset accuracy: restricted to non-Bas phages whose measured receptor is a class the BASEL
    reference actually contains (the fair test of what the shipped v0 CLAIMS to decode).

Usage:  uv run --with biopython python scripts/lbnl_independent_validate.py --repo <PhageDataSheets/Ecoli_phages>
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

# LBNL BW25113 receptor vocabulary -> our RECEPTOR_CLASSES. Combos take the primary (first) token.
RECEPTOR_MAP = {
    "BtuB": "BtuB", "Tsx": "Tsx", "FhuA": "FhuA", "OmpC": "OmpC", "OmpA": "OmpA", "OmpF": "OmpF",
    "LamB": "LamB", "YncD": "YncD", "TolC": "TolC", "FadL": "FadL", "LptD": "LptD", "LPS": "LPS_core",
}
# excluded (not a single mappable receptor): no infection / unresolved / class we don't model
EXCLUDE = {"Resistant", "NGR", "OmpW", "Unknown", "Unknow", ""}


def map_receptor(raw: str) -> str | None:
    raw = (raw or "").strip()
    primary = raw.split(";")[0].strip()   # FhuA;NupG -> FhuA
    if primary in EXCLUDE or not primary:
        return None
    return RECEPTOR_MAP.get(primary)


def ingest(repo: Path, out_dir: Path) -> Path:
    """Parse Table_S1, convert non-Bas genomes .gbk->.fna, write the independent manifest. Returns manifest path."""
    from Bio import SeqIO
    tsv = repo / "data" / "Table_S1_Phages.tsv"
    gdir = repo / "data" / "phage_genomes"
    out_dir.mkdir(parents=True, exist_ok=True)
    manifest = out_dir.parent / "lbnl_independent_manifest.tsv"
    rows = []
    with open(tsv, "r", encoding="utf-8", errors="replace") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(header)}
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) <= idx["Receptor-binding protein"]:
                continue
            phage = p[idx["Phage"]].strip()
            if not phage or phage.startswith("Bas"):     # LEAKAGE GUARD: drop the BASEL collection
                continue
            receptor = map_receptor(p[idx["BW25113 receptor"]])
            if receptor is None:
                continue
            gbk = gdir / f"{phage}.gbk"
            if not gbk.exists():
                continue
            fna = out_dir / f"{phage}.fna"
            try:
                recs = list(SeqIO.parse(str(gbk), "genbank"))
                if not recs:
                    continue
                with open(fna, "w", encoding="utf-8") as out:
                    for r in recs:
                        out.write(f">{phage}_{r.id}\n{str(r.seq)}\n")
            except Exception:
                continue
            rows.append((phage, receptor, p[idx["BW25113 receptor"]].strip(),
                         p[idx["Genus"]].strip(), p[idx["Receptor-binding protein"]].strip()))
    with open(manifest, "w", encoding="utf-8") as out:
        out.write("accession\treceptor\traw_receptor\tgenus\trbp_cds\n")
        for r in rows:
            out.write("\t".join(r) + "\n")
    print(f"ingested {len(rows)} non-Bas LBNL phages -> {manifest}")
    return manifest


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", required=True, help="path to PhageDataSheets/Ecoli_phages")
    ap.add_argument("--ref-manifest", default="data/phage_ref/basel_manifest.tsv")
    ap.add_argument("--ref-dir", default="data/phage_ref/basel")
    ap.add_argument("--out-dir", default="data/phage_ref/lbnl")
    ap.add_argument("--date", default="2026-07-24")
    args = ap.parse_args()

    from dna_decode.phage.receptor_caller import independent_validate, _load_manifest

    out_dir = Path(args.out_dir)
    manifest = ingest(Path(args.repo), out_dir)

    ref_refs, ref_receptors = _load_manifest(args.ref_manifest, args.ref_dir)
    covered_classes = sorted(set(ref_receptors.values()))
    print(f"BASEL-2021 reference: {len(ref_refs)} phages, covers receptor classes {covered_classes}")

    res = independent_validate(args.ref_manifest, args.ref_dir, str(manifest), str(out_dir))

    # covered-subset = test phages whose TRUE receptor is a class the reference contains
    covered = [p for p in res.predictions if p["true"] in covered_classes]
    cov_called = [p for p in covered if p["status"] == "CALLED"]
    cov_correct = sum(1 for p in cov_called if p["correct"])
    from collections import Counter
    true_dist = dict(Counter(p["true"] for p in res.predictions))

    out = {
        "cell": "phage_receptor_class", "date": args.date,
        "independent_source": "LBNL/Arkin-Mutalik Phage Datasheets (Moriniere et al.; github.com/mjohnson11/PhageDataSheets)",
        "label": "measured receptor via genome-wide genetic screens on E. coli K-12 BW25113",
        "independence": "DIFFERENT LAB, measured labels, non-Bas isolates disjoint from the BASEL-2021 reference "
                        "(Bas## excluded as leakage); K-12 host comparable to BASEL",
        "reference": "BASEL-2021 genome-homology caller (shipped v0)",
        "reference_covered_classes": covered_classes,
        "n_test_total": res.n_total, "n_called": res.n_called, "n_correct": res.n_correct,
        "overall_accuracy_on_called": res.accuracy,
        "test_true_receptor_distribution": true_dist,
        "covered_subset_n": len(covered), "covered_subset_called": len(cov_called),
        "covered_subset_correct": cov_correct,
        "covered_subset_accuracy": (cov_correct / len(cov_called)) if cov_called else None,
        "per_receptor_correct_called": res.per_receptor,
        "honest_reading": "The OVERALL number is dragged down by CLASS COVERAGE: the shipped v0 catalogue "
                          "only contains reference phages for classes " + str(covered_classes) + ", so it "
                          "cannot predict the RBP-variable classes (Tsx/OmpC/FhuA/OmpA/OmpF/...) that dominate "
                          "the independent set — those are exactly the (3) RBP-caller scope. The COVERED-SUBSET "
                          "accuracy is the fair independent test of what v0 CLAIMS to decode.",
    }
    wiki = Path("wiki")
    (wiki / f"phage_independent_result_{args.date}.json").write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps({k: out[k] for k in ("n_test_total", "n_called", "overall_accuracy_on_called",
          "covered_subset_n", "covered_subset_called", "covered_subset_correct", "covered_subset_accuracy",
          "reference_covered_classes", "per_receptor_correct_called")}, indent=2))


if __name__ == "__main__":
    main()
