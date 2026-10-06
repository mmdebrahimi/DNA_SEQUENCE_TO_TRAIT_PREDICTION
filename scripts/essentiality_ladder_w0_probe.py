"""W0 empirical probe — pin every transfer-ladder label file's REAL schema before any parser exists.

THE "BUILD AGAINST REAL DATA" GATE, and it exists because writing parsers from prose failed twice inside
one session on this very plan:

  1. A candidate source was the WRONG SHAPE ENTIRELY. `uniprot_ecoli_p1.tsv` was assumed to be a second
     E. coli product annotation; it is 501 rows of `Gene Names (primary) / Gene Names (ordered locus) /
     Disruption phenotype` with NO product column at all.
  2. A header was NOT ROW 0. `goodall_TableS1_essential.xlsx` row 0 is a TITLE
     ("Table S1. Essentiality classification for genes of the TL data") and the real header is row 1.
     Reading row 0 as the header silently mis-joined to 3936 essential / 372 non-essential against the
     committed 351 / 3432 -- a confident, plausible, completely wrong label set.

Neither failure raised. Both produced a number. So this probe records, per source, the facts a parser
needs -- and `w0_verified` is what `label_sources.parse_or_refuse` reads before it will parse anything.

PURE summarizers are split from the fetch/read so they unit-test offline, mirroring
`scripts/oxford_w0_probe.py`. Writes `wiki/essentiality_ladder_w0_probe_<date>.json`.

    uv run python scripts/essentiality_ladder_w0_probe.py --self-check     # offline, no network
    uv run python scripts/essentiality_ladder_w0_probe.py                  # probe local sources
    uv run python scripts/essentiality_ladder_w0_probe.py --fetch-remote   # + attempt remote fetches

Exit 0 always for the local pass (it is a REPORT). `--require-all` exits 1 when any declared source is
unverified, for use as a gate.
"""
from __future__ import annotations

import argparse
import gzip
import json
import sys
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

CACHE = Path("D:/dna_decode_cache/essentiality")
WIKI = ROOT / "wiki"

# ---------------------------------------------------------------------------------------------------
# Declared sources. `local` names a file expected in CACHE; `url` is the remote to try when absent.
# `strain` is recorded because the KEY SPACE belongs to a strain, which is precisely what the FBA arm's
# "walls" are about (NTML is USA300 JE2/FPR3757 `SAUSA300_*`; PLOS GOLD is PA14-keyed).
# ---------------------------------------------------------------------------------------------------
SOURCES = {
    "goodall_tradis": {
        "local": "goodall_TableS1_essential.xlsx", "kind": "xlsx", "organism": "Escherichia coli",
        "strain": "K-12 MG1655", "key_space": "gene_symbol", "technology": "transposon_insertion",
        "url": None,
        "recorded_in": "scripts/essentiality_e3_learned.py::load_labels",
    },
    "bagel_ceg": {
        "local": "CEGv2.txt", "kind": "tsv", "organism": "Homo sapiens", "strain": "n/a",
        "key_space": "entrez_gene_id", "technology": "crispr_ko", "url": None,
        "recorded_in": "scripts/build_essentiality_report_card.py::human_transfer_auroc",
    },
    "bagel_neg": {
        "local": "NEGv1.txt", "kind": "tsv", "organism": "Homo sapiens", "strain": "n/a",
        "key_space": "entrez_gene_id", "technology": "crispr_ko", "url": None,
        "recorded_in": "scripts/build_essentiality_report_card.py::human_transfer_auroc",
    },
    "sgd_phenotype": {
        "local": "phenotype_data.tab", "kind": "tsv", "organism": "Saccharomyces cerevisiae",
        "strain": "S288C", "key_space": "systematic_orf", "technology": "deletion_collection",
        "url": "http://sgd-archive.yeastgenome.org/curation/literature/phenotype_data.tab",
        "recorded_in": "dna_decode/fba/essentiality_labels.py::parse_sgd_essential",
    },
    "ntml_nebraska": {
        "local": "ntml_essential.tsv", "kind": "tsv", "organism": "Staphylococcus aureus",
        "strain": "USA300 JE2 / FPR3757", "key_space": "sausa300_locus_tag",
        "technology": "transposon_insertion", "url": None,
        "recorded_in": "dna_decode/fba/essentiality_labels.py (LABEL_WALLED: model-gene-join wall only)",
    },
    "plos_gold_pa14": {
        "local": "pcbi.1013945.s011.xlsx", "kind": "xlsx", "organism": "Pseudomonas aeruginosa",
        "strain": "PA14", "key_space": "pa14_locus_tag", "technology": "transposon_insertion",
        "url": None,
        "recorded_in": "dna_decode/fba/essentiality_labels.py (MODEL_WALLED: no GEM; label is PA14-keyed)",
    },
}

# A source whose header looks like THIS is the uniprot trap: real columns, wrong KIND of columns.
WRONG_SHAPE_MARKERS = ("disruption phenotype",)


# ---------------------------------------------------------------------------------------------------
# PURE summarizers — no IO, unit-tested offline
# ---------------------------------------------------------------------------------------------------
def detect_header_row(rows, min_cols: int = 2, scan: int = 10) -> int | None:
    """Index of the first row that is a HEADER, not a title.

    THE GOODALL FIX. A spreadsheet title row is a single populated cell; a header has several. Scanning
    for the first row with >= `min_cols` non-empty cells finds row 1 in the Goodall workbook, where
    assuming row 0 silently mis-joins. Returns None when no row in the first `scan` qualifies -- an
    honest 'cannot tell', never a guess of 0.
    """
    for i, row in enumerate(rows[:scan]):
        if row is None:
            continue
        filled = [c for c in row if c is not None and str(c).strip() != ""]
        if len(filled) >= min_cols:
            return i
    return None


def normalise_header(row) -> list[str]:
    return [("" if c is None else str(c).strip()) for c in row]


def looks_wrong_shape(header: list[str]) -> str | None:
    """Detect the uniprot trap: a table with real columns that are the wrong KIND of columns."""
    low = [h.lower() for h in header]
    for marker in WRONG_SHAPE_MARKERS:
        if any(marker in h for h in low):
            return (f"header contains {marker!r} -- this is a phenotype-DESCRIPTION table, not a "
                    f"label or product table (the uniprot_ecoli_p1.tsv trap)")
    return None


def key_column_candidates(header: list[str]) -> list[int]:
    """Columns plausibly holding the join key. Index-based so a duplicate name cannot collapse them."""
    out = []
    for i, h in enumerate(header):
        hl = h.lower()
        if not hl:
            continue
        if any(t in hl for t in ("gene", "orf", "locus", "systematic", "id", "symbol", "tag")):
            out.append(i)
    return out


def label_column_candidates(header: list[str]) -> list[int]:
    """CORROBORATION ONLY -- never the header guard. Measured weakness: the Goodall title row reads
    "Table S1. ESSENTIALITY classification for genes of the TL data", so this lexical match fires on
    `essential` inside `Essentiality` and nominates the TITLE's column 0 as a label column. So a
    row-0-as-header read would NOT be caught by the label-column check; `detect_header_row` is what
    actually catches it. Pinned by test_the_label_column_detector_is_FOOLED_by_a_descriptive_title.
    """
    out = []
    for i, h in enumerate(header):
        hl = h.lower()
        if any(t in hl for t in ("essential", "non-essential", "nonessential", "unclear",
                                 "inviable", "viable", "phenotype", "call")):
            out.append(i)
    return out


def summarize_table(rows, *, min_cols: int = 2) -> dict:
    """PURE schema summary of a row-of-cells table. The record every parser is written against."""
    hdr_idx = detect_header_row(rows, min_cols=min_cols)
    if hdr_idx is None:
        return {"header_row_index": None, "n_rows_total": len(rows),
                "problem": "no row in the first 10 has >= %d populated cells" % min_cols}
    header = normalise_header(rows[hdr_idx])
    data = [r for r in rows[hdr_idx + 1:] if r is not None and any(
        c is not None and str(c).strip() != "" for c in r)]
    keycols = key_column_candidates(header)
    first_key = keycols[0] if keycols else None
    keys = []
    if first_key is not None:
        for r in data:
            if first_key < len(r) and r[first_key] is not None:
                v = str(r[first_key]).strip()
                if v:
                    keys.append(v)
    out = {
        "header_row_index": hdr_idx,
        "title_rows_skipped": hdr_idx,
        "header": header,
        "n_columns": len(header),
        "n_rows_total": len(rows),
        "n_data_rows": len(data),
        "key_column_candidates": keycols,
        "key_column_chosen": first_key,
        "label_column_candidates": label_column_candidates(header),
        "n_keys": len(keys),
        "n_keys_distinct": len(set(keys)),
        "n_duplicate_keys": len(keys) - len(set(keys)),
        "representative_keys": keys[:5],
    }
    ws = looks_wrong_shape(header)
    if ws:
        out["wrong_shape"] = ws
    return out


def verify(summary: dict, declared: dict) -> tuple[bool, str]:
    """Is this summary enough for a parser to be written against? Fail-closed and specific."""
    if summary.get("problem"):
        return False, summary["problem"]
    if summary.get("wrong_shape"):
        return False, summary["wrong_shape"]
    if summary.get("key_column_chosen") is None:
        return False, "no key-column candidate in the header"
    if summary.get("n_data_rows", 0) == 0:
        return False, "zero data rows below the detected header"
    if not summary.get("label_column_candidates") and declared.get("kind") == "xlsx":
        return False, "no label-column candidate -- a label table must carry a call column"
    return True, "schema pinned: header row %d, key column %d, %d data rows, %d distinct keys" % (
        summary["header_row_index"], summary["key_column_chosen"],
        summary["n_data_rows"], summary["n_keys_distinct"])


# ---------------------------------------------------------------------------------------------------
# Thin IO
# ---------------------------------------------------------------------------------------------------
def read_rows(path: Path, kind: str):
    """-> (rows, sheet_names). Rows are lists of cells; sheet_names is [] for non-xlsx."""
    if kind == "xlsx":
        import openpyxl
        wb = openpyxl.load_workbook(path, read_only=True)
        return list(wb[wb.sheetnames[0]].iter_rows(values_only=True)), list(wb.sheetnames)
    op = gzip.open if str(path).endswith(".gz") else open
    with op(path, "rt", encoding="utf-8", errors="replace") as fh:
        return [ln.rstrip("\n").split("\t") for ln in fh], []


def probe_source(sid: str, declared: dict) -> dict:
    rec = {"source_id": sid, "organism": declared["organism"], "strain": declared["strain"],
           "declared_key_space": declared["key_space"], "technology": declared["technology"],
           "kind": declared["kind"], "url": declared["url"], "recorded_in": declared["recorded_in"]}
    p = CACHE / declared["local"]
    rec["local_path"] = str(p)
    if not p.exists():
        rec.update(fetch_outcome="absent_locally", w0_verified=False,
                   reason=("not in %s; remote is %s" % (CACHE, declared["url"] or "NOT RECORDED -- "
                           "reachability is a HYPOTHESIS until fetched")))
        return rec
    try:
        rows, sheets = read_rows(p, declared["kind"])
    except Exception as e:                                   # report, never abort the sweep
        rec.update(fetch_outcome="read_failed", w0_verified=False,
                   reason="%s: %s" % (type(e).__name__, str(e)[:160]))
        return rec
    rec["fetch_outcome"] = "read_local"
    rec["sheet_names"] = sheets
    summary = summarize_table(rows)
    rec["schema"] = summary
    ok, why = verify(summary, declared)
    rec["w0_verified"] = ok
    rec["reason"] = why
    return rec


def _self_check() -> int:
    """Offline proof the summarizers catch BOTH recorded failures."""
    goodall_like = [
        ["Table S1. Essentiality classification for genes of the TL data", None, None, None, None, None],
        ["Gene", "Insertion Index Score", "Log Likelihood Ratio", "Essential", "Non-essential", "Unclear"],
        ["thrL", "0.39", "31.4", "False", "True", "False"],
        ["thrA", "0.22", "17.2", "False", "True", "False"],
    ]
    s = summarize_table(goodall_like)
    assert s["header_row_index"] == 1, s
    assert s["title_rows_skipped"] == 1
    assert s["n_data_rows"] == 2, s
    assert verify(s, SOURCES["goodall_tradis"])[0]

    uniprot_like = [
        ["Gene Names (primary)", "Gene Names (ordered locus)", "Disruption phenotype"],
        ["ndh", "b1109 JW1095", "DISRUPTION PHENOTYPE: ..."],
    ]
    u = summarize_table(uniprot_like)
    ok, why = verify(u, {"kind": "tsv"})
    assert not ok and "Disruption phenotype".lower() in why.lower(), (ok, why)

    title_only = [["Just a title", None, None]]
    assert detect_header_row(title_only) is None, "a title-only sheet must be 'cannot tell', not 0"
    print("self-check OK: Goodall header-row-1 detected, uniprot wrong-shape rejected, "
          "title-only refuses rather than guessing row 0")
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--self-check", action="store_true")
    ap.add_argument("--fetch-remote", action="store_true", help="attempt remote fetches for absent files")
    ap.add_argument("--require-all", action="store_true", help="exit 1 if any source is unverified")
    ap.add_argument("--no-emit", action="store_true")
    a = ap.parse_args(argv)

    if a.self_check:
        return _self_check()

    records = [probe_source(sid, d) for sid, d in SOURCES.items()]
    if a.fetch_remote:
        for r in records:
            if r["w0_verified"] is False and r.get("fetch_outcome") == "absent_locally" and r["url"]:
                r["remote_attempt"] = "NOT IMPLEMENTED in this step -- a fetch is a separate, " \
                                      "manual, network-gated action; this probe never invents one"

    n_ok = sum(1 for r in records if r["w0_verified"])
    art = {
        "record": "essentiality-ladder-w0-probe-v1",
        "date": str(date.today()),
        "purpose": "pin each label source's REAL schema before any parser is written against it",
        "why": "two prose-written parsers failed silently on this plan: a wrong-shape source "
               "(uniprot_ecoli_p1.tsv, 501 rows of Disruption phenotype, no product column) and a "
               "header that is not row 0 (Goodall row 0 is a title; reading it mis-joined 3936/372 "
               "against the true 351/3432). Neither raised -- both produced a number.",
        "n_sources": len(records), "n_verified": n_ok,
        "sources": records,
        "consumed_by": "dna_decode/essentiality/label_sources.py::parse_or_refuse -- refuses to parse a "
                       "source whose w0_verified is false or absent",
        "honest_limits": [
            "a verified schema says a parser CAN be written, not that the labels are correct",
            "absent_locally is not 'unreachable' -- it means not fetched yet; reachability stays a "
            "HYPOTHESIS until a fetch happens (wiki/essentiality_label_wall_2026-07-28.md records "
            "seven dead routes for essentiality labels generally)",
            "key-column and label-column detection is lexical; it proposes candidates and records the "
            "chosen index, it does not understand the biology",
        ],
    }
    print("W0 probe: %d/%d sources verified" % (n_ok, len(records)))
    for r in records:
        mark = "OK " if r["w0_verified"] else "-- "
        print("  %s%-16s %-26s %s" % (mark, r["source_id"], r["organism"], r["reason"][:90]))
        sc = r.get("schema") or {}
        if sc.get("header_row_index") is not None:
            print("       header_row=%s  title_rows_skipped=%s  key_col=%s  data_rows=%s  distinct=%s"
                  % (sc["header_row_index"], sc["title_rows_skipped"], sc["key_column_chosen"],
                     sc["n_data_rows"], sc["n_keys_distinct"]))
    if not a.no_emit:
        p = WIKI / ("essentiality_ladder_w0_probe_%s.json" % date.today())
        p.write_text(json.dumps(art, indent=2), encoding="utf-8")
        print("[-> %s]" % p.relative_to(ROOT))
    if a.require_all and n_ok != len(records):
        print("REQUIRE-ALL: %d source(s) unverified" % (len(records) - n_ok))
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
