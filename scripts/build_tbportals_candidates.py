"""Convert a TB Portals data export -> the candidate TSV that `build_tb_goldset_manifest.py` consumes.

TB Portals (NIAID, DUA-gated) is the no-author-contact powered independent TB source (shortlist #3). After
the DUA is signed you DOWNLOAD a de-identified export; the genomic reads are PUBLIC on SRA. This adapter
turns the export's DST table into the per-isolate `(accession, rif_label, inh_label)` candidate TSV, then the
existing path takes over:  build_tbportals_candidates -> validate_tb_goldset_candidates -> build_tb_goldset_manifest
(leakage check) -> fetch SRA reads -> variant-call -> score_tb_independent_goldset.

TWO modes (W0-probe discipline — pin the real schema BEFORE trusting it):
  --probe  <export.csv>            : report detected columns + the drug / result / METHOD vocabularies so you
                                     confirm which methods are phenotypic. Run this FIRST on the real file.
  --build  <export.csv> --out X.tsv: write the candidate TSV (phenotypic-only; long-format auto-pivoted).

CIRCULARITY GUARD (gate G1 — load-bearing): TB Portals carries BOTH phenotypic (measured: MGIT/LJ/proportion)
AND molecular/genotypic (Hain/Xpert/WGS/LPA) DST. ONLY phenotypic clears the non-circular bar — scoring our
WHO-catalogue rule against a *molecular* DST call is rule-vs-rule (the same circularity that made Thorpe a
NO-GO). `--build` KEEPS only rows whose method matches a phenotypic token and DROPS molecular rows (logged).
If no method column is detected, `--build` REFUSES unless `--no-method-column-ok` is passed (you assert the
export is phenotypic-only) — it will NOT silently let molecular DST through.

Stdlib only (csv). Pure helpers (`normalize_result` / `classify_method` / `detect_columns` / `pivot_dst`)
are unit-tested; `main` is the CLI wrapper.
"""
from __future__ import annotations

import argparse
import csv
import sys
from collections import Counter, defaultdict
from pathlib import Path

# --- vocabularies (lowercased substring match) -----------------------------------------------------------
PHENOTYPIC_TOKENS = ("mgit", "bactec", "lj", "lowenstein", "proportion", "phenotyp", "broth", "agar",
                     "microdilution", "absolute concentration", "resistance ratio", "culture", "wayne")
MOLECULAR_TOKENS = ("hain", "genotyp", "molecular", "xpert", "genexpert", "line probe", "lpa", "wgs",
                    "sequenc", "pcr", "mtbdr", "geno type", "whole genome", "ngs")

RESISTANT_TOKENS = ("resist",)               # "Resistant", "R"
SUSCEPTIBLE_TOKENS = ("suscept", "sensitive")  # "Susceptible", "Sensitive", "S"
# Intermediate / Indeterminate / Not tested / Contaminated -> blank (no usable measured call)

RIF_TOKENS = ("rifampicin", "rifampin", "rif")
INH_TOKENS = ("isoniazid", "inh")

# --- TB Portals WIDE format (Patient_Cases spreadsheet; verified vs the 2026 data dictionary) ------------
# Per-drug DST is one column per (method x drug): <method>_<drug>. Method is encoded in the column PREFIX.
# Phenotypic (culture-based, measured -> KEEP): bactec_ (BACTEC MIC), le_ (Lowenstein-Jensen).
# Molecular (genotype-derived -> DROP, gate G1): genexpert_, hain_ (Hain LPA), lpaother_, truenat_.
TBPORTALS_PHENO_COLS = {"RIF": ("bactec_rifampicin", "le_rifampicin"),
                        "INH": ("bactec_isoniazid", "le_isoniazid")}
TBPORTALS_MOLECULAR_COLS = ("genexpert_rifampicin", "hain_rifampicin", "lpaother_rifampicin",
                            "truenat_rifampicin", "genexpert_isoniazid", "hain_isoniazid",
                            "lpaother_isoniazid", "truenat_isoniazid")
TBPORTALS_ID_COLS = ("condition_id", "patient_id", "identifier")
# alias column -> TB Portals NCBI column. ncbi_biosample (SAMN/SAMEA) is the clean cross-archive leakage key.
TBPORTALS_ALIAS_MAP = {"biosample_accession": "ncbi_biosample", "sample_accession": "ncbi_sra"}

ACCESSION_HINTS = ("run_accession", "sra", "ena", "biosample", "sample_accession", "accession",
                   "ncbi", "srr", "err", "drr")
DRUG_HINTS = ("drug", "antibiotic", "agent", "compound")
RESULT_HINTS = ("result", "dst_result", "interpretation", "susceptibility", "phenotype", "outcome")
METHOD_HINTS = ("method", "dst_method", "test_method", "technique", "assay", "platform", "test_type")
ID_HINTS = ("condition_id", "patient_id", "specimen", "sample_id", "isolate", "strain", "case_id", "id")


def _norm(s) -> str:
    return (s or "").strip().lower()


def normalize_result(s) -> str:
    """Map a TB Portals DST result value -> 'R' / 'S' / '' (blank for anything not a clean measured call)."""
    v = _norm(s)
    if not v:
        return ""
    if v in ("r", "s"):
        return v.upper()
    if any(t in v for t in RESISTANT_TOKENS):
        return "R"
    if any(t in v for t in SUSCEPTIBLE_TOKENS):
        return "S"
    return ""   # intermediate / indeterminate / not tested / contaminated / unknown


def classify_method(s) -> str:
    """'phenotypic' | 'molecular' | 'unknown' from a DST method/technique string."""
    v = _norm(s)
    if not v:
        return "unknown"
    # molecular wins ties (a row labelled "WGS-confirmed phenotypic" is genotype-derived -> circular).
    if any(t in v for t in MOLECULAR_TOKENS):
        return "molecular"
    if any(t in v for t in PHENOTYPIC_TOKENS):
        return "phenotypic"
    return "unknown"


def classify_drug(s) -> str:
    """'RIF' | 'INH' | '' from a drug name."""
    v = _norm(s)
    # INH check before RIF is irrelevant (disjoint), but guard 'rif' not matching inside other words.
    if any(t == v or t in v.split() or t in v for t in INH_TOKENS):
        if any(t in v for t in INH_TOKENS):
            return "INH"
    if any(t in v for t in RIF_TOKENS):
        return "RIF"
    if any(t in v for t in INH_TOKENS):
        return "INH"
    return ""


def parse_tbportals_result(cell) -> str:
    """Parse a TB Portals wide DST cell -> 'R' / 'S' / '' (blank).

    Cells are coded R / S / I / Ind / Not Reported; multiple records are aggregated like '{S, R}'.
    A cell carrying BOTH R and S (discordant aggregate) is ambiguous -> blank (no clean measured call)."""
    v = _norm(cell)
    if not v or "not reported" in v or v in ("nan", "none"):
        return ""
    up = v.upper()
    import re as _re
    has_r = bool(_re.search(r"\bR\b", up)) or "RESIST" in up
    has_s = bool(_re.search(r"\bS\b", up)) or "SENSITIV" in up or "SUSCEPT" in up
    if has_r and has_s:
        return ""          # {S, R} discordant aggregate
    if has_r:
        return "R"
    if has_s:
        return "S"
    return ""              # I / Ind / unrecognized


def is_tbportals_wide(header: list[str]) -> bool:
    """True if the export is the TB Portals Patient_Cases WIDE shape (method-prefixed drug columns)."""
    hs = set(header)
    return any(c in hs for cols in TBPORTALS_PHENO_COLS.values() for c in cols)


def pivot_tbportals_wide(rows: list[dict], header: list[str]) -> tuple[list[dict], dict]:
    """TB Portals Patient_Cases wide rows -> one candidate dict per isolate.

    Combines phenotypic methods (bactec + le) per drug: agree -> that call; disagree -> blank (conflict).
    Molecular columns (genexpert/hain/lpaother/truenat) are NEVER read (gate G1). Returns (candidates, stats)."""
    hs = set(header)
    id_col = next((c for c in TBPORTALS_ID_COLS if c in hs), None)
    alias_present = {a: col for a, col in TBPORTALS_ALIAS_MAP.items() if col in hs}
    stats = Counter()
    stats["molecular_cols_ignored"] = sum(1 for c in TBPORTALS_MOLECULAR_COLS if c in hs)
    candidates = []
    for r in rows:
        iso = str(r.get(id_col, "")).strip() if id_col else ""
        if not iso:
            stats["no_isolate_id"] += 1
            continue
        labels = {}
        for drug, cols in TBPORTALS_PHENO_COLS.items():
            calls = {parse_tbportals_result(r.get(c)) for c in cols if c in hs}
            calls.discard("")
            if calls == {"R"}:
                labels[drug] = "R"; stats[f"kept_{drug.lower()}_R"] += 1
            elif calls == {"S"}:
                labels[drug] = "S"; stats[f"kept_{drug.lower()}_S"] += 1
            elif len(calls) > 1:
                labels[drug] = ""; stats[f"conflict_{drug.lower()}"] += 1   # bactec vs le disagree
            else:
                labels[drug] = ""
        if not (labels.get("RIF") or labels.get("INH")):
            stats["no_usable_label"] += 1
            continue
        row = {"strain_id": iso, "masked_vcf": "", "regeno_vcf": "",
               "rif_label": labels.get("RIF", ""), "inh_label": labels.get("INH", ""),
               "run_accession": "", "sample_accession": "", "biosample_accession": ""}
        for alias, col in alias_present.items():
            row[alias] = str(r.get(col, "")).strip()
        candidates.append(row)
    stats["n_isolates_usable"] = len(candidates)
    return candidates, dict(stats)


def detect_columns(header: list[str]) -> dict:
    """Best-effort column detection. Returns dict of role -> column name (or None)."""
    low = {h: _norm(h) for h in header}

    def pick(hints):
        # exact-ish first (hint == column), then substring.
        for h in header:
            if low[h] in hints:
                return h
        for h in header:
            if any(hint in low[h] for hint in hints):
                return h
        return None

    return {"id": pick(ID_HINTS), "accession": pick(ACCESSION_HINTS), "drug": pick(DRUG_HINTS),
            "result": pick(RESULT_HINTS), "method": pick(METHOD_HINTS)}


def _accession_cols(header: list[str]) -> dict:
    """Map our 3 alias columns to whatever accession-bearing columns the export has."""
    low = {h: _norm(h) for h in header}
    out = {}
    for alias, toks in (("run_accession", ("run_accession", "srr", "err", "drr", "run")),
                        ("sample_accession", ("sample_accession", "srs", "ers", "drs")),
                        ("biosample_accession", ("biosample", "samea", "samn", "samd"))):
        for h in header:
            if any(t in low[h] for t in toks):
                out[alias] = h
                break
    return out


def pivot_dst(rows: list[dict], cols: dict, acc_cols: dict, *,
              keep_methods=("phenotypic", "unknown"), method_required=True) -> tuple[list[dict], dict]:
    """Long-format DST rows -> one candidate dict per isolate. Returns (candidates, stats).

    cols: detect_columns() result. acc_cols: alias->column map. keep_methods: which classify_method buckets
    to keep (default phenotypic + unknown; pass ('phenotypic',) to be strict). method_required: if True and
    no method column, the caller should refuse (this fn just records method='unknown' for all rows)."""
    id_col = cols.get("id") or cols.get("accession")
    drug_col, result_col, method_col = cols.get("drug"), cols.get("result"), cols.get("method")
    stats = Counter()
    per_iso: dict[str, dict] = defaultdict(lambda: {"rif": "", "inh": "", "aliases": {}})

    for r in rows:
        iso = str(r.get(id_col, "")).strip() if id_col else ""
        if not iso:
            stats["no_isolate_id"] += 1
            continue
        method = classify_method(r.get(method_col)) if method_col else "unknown"
        if method not in keep_methods:
            stats[f"dropped_method_{method}"] += 1
            continue
        drug = classify_drug(r.get(drug_col)) if drug_col else ""
        if drug not in ("RIF", "INH"):
            stats["not_rif_inh"] += 1
            continue
        res = normalize_result(r.get(result_col)) if result_col else ""
        # record aliases regardless (so even a blank-result isolate keeps its accession if later filled)
        for alias, col in acc_cols.items():
            val = str(r.get(col, "")).strip()
            if val:
                per_iso[iso]["aliases"][alias] = val
        if res in ("R", "S"):
            key = "rif" if drug == "RIF" else "inh"
            prior = per_iso[iso][key]
            if prior and prior != res:
                stats[f"conflict_{key}"] += 1   # discordant repeat tests -> keep first, flag
            elif not prior:
                per_iso[iso][key] = res
                stats[f"kept_{key}_{res}"] += 1
        else:
            stats["blank_result"] += 1

    candidates = []
    for iso, d in per_iso.items():
        if not (d["rif"] or d["inh"]):
            continue   # no usable measured label for either drug
        row = {"strain_id": iso, "masked_vcf": "", "regeno_vcf": "",
               "rif_label": d["rif"], "inh_label": d["inh"]}
        row.update({a: d["aliases"].get(a, "") for a in ("run_accession", "sample_accession",
                                                         "biosample_accession")})
        candidates.append(row)
    stats["n_isolates_usable"] = len(candidates)
    return candidates, dict(stats)


# --- CLI -------------------------------------------------------------------------------------------------
def _read_csv(path: Path) -> tuple[list[dict], list[str]]:
    # tolerate csv or tsv by sniffing the delimiter on the header line
    with open(path, encoding="utf-8-sig") as fh:
        sample = fh.readline()
        delim = "\t" if sample.count("\t") > sample.count(",") else ","
    with open(path, encoding="utf-8-sig", newline="") as fh:
        rd = csv.DictReader(fh, delimiter=delim)
        return list(rd), (rd.fieldnames or [])


def _probe(path: Path) -> int:
    rows, header = _read_csv(path)
    print(f"[probe] {path}  rows={len(rows)}  delimiter-sniffed")
    if is_tbportals_wide(header):
        hs = set(header)
        print("[probe] FORMAT: TB Portals WIDE (Patient_Cases) — method-prefixed drug columns detected.")
        id_col = next((c for c in TBPORTALS_ID_COLS if c in hs), None)
        print(f"[probe] isolate id column: {id_col or '(NONE — cannot key isolates!)'}")
        for drug, pcols in TBPORTALS_PHENO_COLS.items():
            present = [c for c in pcols if c in hs]
            print(f"[probe]   {drug} phenotypic cols present (KEEP): {present or 'NONE'}")
        molc = [c for c in TBPORTALS_MOLECULAR_COLS if c in hs]
        print(f"[probe]   molecular cols present (DROP, gate G1): {molc}")
        alias = {a: col for a, col in TBPORTALS_ALIAS_MAP.items() if col in hs}
        print(f"[probe]   accession columns: {alias or '(NONE — leakage check will be weak)'}")
        cands, stats = pivot_tbportals_wide(rows, header)
        rif = sum(1 for c in cands if c["rif_label"] in ("R", "S"))
        inh = sum(1 for c in cands if c["inh_label"] in ("R", "S"))
        print(f"[probe]   -> {len(cands)} isolates with a phenotypic label  (RIF={rif}  INH={inh})")
        print(f"[probe]   stats: {stats}")
        return 0
    cols = detect_columns(header)
    acc = _accession_cols(header)
    print("[probe] FORMAT: generic long-format (no TB Portals wide columns detected).")
    print(f"[probe] columns ({len(header)}): {header}")
    print(f"[probe] detected roles: {cols}")
    print(f"[probe] accession-alias columns: {acc or '(NONE detected — leakage check will be weak)'}")
    if cols.get("drug"):
        dv = Counter(classify_drug(r.get(cols['drug'])) or f"OTHER:{_norm(r.get(cols['drug']))[:20]}"
                     for r in rows)
        print(f"[probe] drug vocab (mapped): {dict(dv.most_common(12))}")
    if cols.get("result"):
        rv = Counter(normalize_result(r.get(cols['result'])) or f"BLANK:{_norm(r.get(cols['result']))[:20]}"
                     for r in rows)
        print(f"[probe] result vocab (mapped): {dict(rv.most_common(12))}")
    if cols.get("method"):
        mv = Counter(classify_method(r.get(cols['method'])) for r in rows)
        raw = Counter(_norm(r.get(cols['method']))[:24] for r in rows)
        print(f"[probe] method buckets: {dict(mv)}")
        print(f"[probe] raw method values: {dict(raw.most_common(12))}")
        print("[probe] CIRCULARITY NOTE: only 'phenotypic' rows are non-circular; 'molecular' rows are "
              "dropped by --build (gate G1).")
    else:
        print("[probe] NO method column detected — --build will REFUSE unless --no-method-column-ok "
              "(you must assert the export is phenotypic-only, else molecular DST silently leaks in).")
    return 0


def _write_candidates(cands: list[dict], out: Path) -> None:
    fields = ["strain_id", "run_accession", "sample_accession", "biosample_accession",
              "masked_vcf", "regeno_vcf", "rif_label", "inh_label"]
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        for c in cands:
            w.writerow({k: c.get(k, "") for k in fields})


def _build(path: Path, out: Path, *, strict_phenotypic: bool, no_method_ok: bool) -> int:
    rows, header = _read_csv(path)
    # TB Portals wide path: phenotypic columns are hardcoded-known, so the molecular drop is structural.
    if is_tbportals_wide(header):
        cands, stats = pivot_tbportals_wide(rows, header)
        _write_candidates(cands, out)
        rif = sum(1 for c in cands if c["rif_label"] in ("R", "S"))
        inh = sum(1 for c in cands if c["inh_label"] in ("R", "S"))
        print(f"[build] TB Portals WIDE: wrote {len(cands)} isolates -> {out}")
        print(f"[build] measured (phenotypic: bactec+le) labels: RIF={rif}  INH={inh}")
        print(f"[build] molecular DST columns IGNORED by construction (gate G1): "
              f"{stats.get('molecular_cols_ignored', 0)}")
        if not any(a in header for a in TBPORTALS_ALIAS_MAP.values()):
            print("WARNING: no ncbi_biosample/ncbi_sra column — add accessions before the leakage check.",
                  file=sys.stderr)
        print(f"[build] stats: {stats}")
        print("Next:\n"
              f"  uv run python scripts/validate_tb_goldset_candidates.py {out}\n"
              f"  uv run python -m scripts.build_tb_goldset_manifest --candidates {out}")
        return 0
    cols = detect_columns(header)
    acc = _accession_cols(header)
    if cols.get("method") is None and not no_method_ok:
        print("ERROR: no DST-method column detected. Re-run --probe to inspect, then either point at an "
              "export WITH a method column, or pass --no-method-column-ok to assert it is phenotypic-only "
              "(REQUIRED so molecular/genotypic DST can't silently make the gold set circular).",
              file=sys.stderr)
        return 2
    if not acc:
        print("WARNING: no accession-alias column detected — the CRyPTIC leakage check will have nothing to "
              "match on. You must add run/sample/biosample accessions before scoring.", file=sys.stderr)
    keep = ("phenotypic",) if strict_phenotypic else ("phenotypic", "unknown")
    cands, stats = pivot_dst(rows, cols, acc, keep_methods=keep)
    fields = ["strain_id", "run_accession", "sample_accession", "biosample_accession",
              "masked_vcf", "regeno_vcf", "rif_label", "inh_label"]
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        for c in cands:
            w.writerow({k: c.get(k, "") for k in fields})
    rif = sum(1 for c in cands if c["rif_label"] in ("R", "S"))
    inh = sum(1 for c in cands if c["inh_label"] in ("R", "S"))
    print(f"[build] wrote {len(cands)} isolates -> {out}")
    print(f"[build] measured labels: RIF={rif}  INH={inh}  (keep_methods={keep})")
    dropped = {k: v for k, v in stats.items() if k.startswith("dropped_method")}
    if dropped:
        print(f"[build] DROPPED molecular/non-phenotypic rows (gate G1): {dropped}")
    print(f"[build] stats: {stats}")
    print("Next:\n"
          f"  uv run python scripts/validate_tb_goldset_candidates.py {out}\n"
          f"  uv run python -m scripts.build_tb_goldset_manifest --candidates {out}")
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("export", type=Path, help="TB Portals DST export (csv or tsv)")
    ap.add_argument("--probe", action="store_true", help="report schema + vocabularies, write nothing")
    ap.add_argument("--build", action="store_true", help="write the candidate TSV")
    ap.add_argument("--out", type=Path, default=Path("data/raw/tb_goldset/tbportals_candidates.tsv"))
    ap.add_argument("--strict-phenotypic", action="store_true",
                    help="keep ONLY rows whose method is classified phenotypic (drop 'unknown' too)")
    ap.add_argument("--no-method-column-ok", action="store_true",
                    help="proceed without a method column (you assert the export is phenotypic-only)")
    a = ap.parse_args(argv)
    if not a.export.exists():
        print(f"ERROR: export not found: {a.export}", file=sys.stderr)
        return 2
    if a.probe == a.build:
        print("ERROR: pass exactly one of --probe / --build", file=sys.stderr)
        return 2
    if a.probe:
        return _probe(a.export)
    return _build(a.export, a.out, strict_phenotypic=a.strict_phenotypic,
                  no_method_ok=a.no_method_column_ok)


if __name__ == "__main__":
    raise SystemExit(main())
