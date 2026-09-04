"""Build a Salmonella serovar validation cohort whose label is WET-LAB, not tool-derived.

WHY THIS CELL AND WHY NOW. `typing:Salmonella:salmserovar` ships at `FAITHFUL_TO_TOOL` -- validated
against the reference method, never against reality. That is 1 of 6 molecular-typing cells in the same
state, and the typing track is 34 of the tool's 115 cells with exactly ONE independently-measured cell in
it. The salmserovar report card already recorded both gates as PASSED and the full-cohort number as
"PENDING -- runnable"; the caller, the real antigen DB and blastn are all present. What was missing is
this: a cohort.

THE ONE THING THAT MATTERS HERE IS THE CIRCULARITY GATE (G1). Salmonella serovar is unusual and
valuable because its gold standard is SLIDE AGGLUTINATION -- a wet-lab antisera reaction, not a
computation. That is the free, independent, isolate-level label the AMR track never had. But public
"serovar" strings are a MIXTURE: some are agglutination results, some are SeqSero2/SISTR output pasted
into a metadata field. Scoring against the latter would be scoring the tool against itself.

NCBI-PD makes the distinction addressable because it carries BOTH fields:
  - `serovar`        -- the SUBMITTER's value (the putative wet-lab label)
  - `computed_types` -- NCBI's own in-silico call (`serotype=...`,`antigen_formula=...`)

THREE FILTERS, EACH DOING DIFFERENT WORK:
  1. REFERENCE-LAB provenance. Restrict to public-health/regulatory labs that run traditional
     serotyping as routine (CDC / PHE / FDA / USDA-FSIS / state health departments). This is the main
     defence and it is a JUDGMENT, stated as one.
  2. TOOL-DISAGREEMENT evidence. Rows where `computed_types` FAILED to type (`I -:-:-`, `-:-:-`) while
     `serovar` is populated PROVE the serovar came from somewhere other than that tool. These are
     counted and reported as positive evidence of independence -- they are NOT used as the cohort,
     because selecting on tool-failure would bias toward hard genomes.
  3. SOURCE DIVERSITY. This project refuses to report a cell whose largest source exceeds
     `MAX_SOURCE_SHARE`; that bar is applied HERE, at construction, rather than discovered afterwards.

WHAT THIS CANNOT DO, stated plainly: it cannot prove any INDIVIDUAL serovar string was produced by
antisera. Residual circularity is bounded, not eliminated -- and the reason the result survives it is
that the scorer reports our caller AND the in-silico comparator against the SAME labels. Contamination
inflates both equally, so the DELTA is robust even where the absolute number is optimistic.

Network-only (streams PD metadata). Writes data/salmserovar_cohort.tsv + a JSON provenance sidecar.
"""
from __future__ import annotations

import argparse
import collections
import json
import random
import sys
import urllib.request
from datetime import date as _date
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from source_diverse_validate import MAX_SOURCE_SHARE  # noqa: E402
from gentamicin_rmt_specificity_hunt import latest_metadata_url  # noqa: E402

# Labs whose routine workflow includes traditional Kauffmann-White slide agglutination. This is the
# load-bearing judgment of the whole script, so it is an explicit, auditable list rather than a regex
# over free text. Matched case-insensitively as a substring of `collected_by` or `sra_center`.
REFERENCE_LABS = (
    "cdc", "centers for disease control",
    "phe", "public health england", "ukhsa",
    "fda", "food and drug administration",
    "usda", "fsis",
    "department of health", "dept of health", "public health laboratory",
    "national reference", "health protection",
)

# `computed_types` values meaning the in-silico caller could not resolve a serovar.
TOOL_FAILED_MARKERS = ("serotype=I -:-:-", "antigen_formula=-:-:-", "serotype=-:-:-")

# Serovar strings that are not a specific serovar and cannot be scored as one.
NON_SPECIFIC = {"", "null", "na", "n/a", "unknown", "not determined", "undetermined",
                "not applicable", "missing", "-", "untypeable", "unnamed",
                "salmonella enterica", "enterica", "i", "ii", "iiia", "iiib", "iv", "v", "vi",
                # In-progress / placeholder values seen in the live field. `pending` is a real
                # example: it survived the first build and would have scored as a serovar miss.
                "pending", "in progress", "tbd", "to be determined", "not typed", "not serotyped",
                "under investigation", "no serovar", "none", "other", "unidentified"}


def norm_serovar(s: str) -> str:
    """Normalise a serovar string for comparison, without collapsing distinct serovars.

    Public strings carry subspecies prefixes and punctuation noise ('Salmonella enterica subsp.
    enterica serovar Typhimurium' / 'serovar Typhimurium' / 'Typhimurium'). Only prefix chrome is
    stripped -- never the serovar token itself.
    """
    t = (s or "").strip().strip('"').lower()
    for cut in ("serovar ", "ser. ", "serotype "):
        if cut in t:
            t = t.split(cut, 1)[1]
    for pre in ("salmonella enterica subsp. enterica ", "salmonella enterica subspecies enterica ",
                "salmonella enterica ", "salmonella "):
        if t.startswith(pre):
            t = t[len(pre):]
    return " ".join(t.replace("_", " ").split()).strip(" .,;")


def parse_computed_types(computed: str) -> dict[str, str]:
    """Parse PD's `computed_types` respecting its QUOTING.

    THE TRAP, and it is the same one the `AST_phenotypes` parser hit. The field looks like
        "serotype=I 4,[5],12:i:-","antigen_formula=4:i:-"
    so it is comma-separated -- but a Salmonella antigenic formula CONTAINS COMMAS. Splitting on a
    bare comma shreds `I 4,[5],12:i:-` into `I 4`, which silently turns monophasic Typhimurium into a
    different, truncated string and manufactures disagreements that are not real. Quoted segments must
    be honoured, so this walks the string and only splits on commas OUTSIDE quotes.
    """
    out: dict[str, str] = {}
    buf, in_q, parts = [], False, []
    for ch in (computed or ""):
        if ch == '"':
            in_q = not in_q
            continue
        if ch == "," and not in_q:
            parts.append("".join(buf)); buf = []
            continue
        buf.append(ch)
    parts.append("".join(buf))
    for p in parts:
        p = p.strip()
        if "=" in p:
            k, v = p.split("=", 1)
            out[k.strip()] = v.strip()
    return out


def parse_computed_serotype(computed: str) -> str | None:
    """The in-silico serotype call, or None when the field carries no `serotype=`."""
    return parse_computed_types(computed).get("serotype") or None


def is_reference_lab(*fields: str) -> bool:
    blob = " ".join((f or "").lower() for f in fields)
    return any(lab in blob for lab in REFERENCE_LABS)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--max-rows", type=int, default=400000, help="rows of PD metadata to stream")
    ap.add_argument("--target", type=int, default=200, help="cohort size")
    ap.add_argument("--per-serovar-cap", type=int, default=12,
                    help="cap per serovar so the cohort is not all Typhimurium/Enteritidis")
    ap.add_argument("--seed", type=int, default=17)
    ap.add_argument("--out-tsv", type=Path, default=ROOT / "data" / "salmserovar_cohort.tsv")
    ap.add_argument("--out-json", type=Path,
                    default=ROOT / "wiki" / f"salmserovar_cohort_{_date.today().isoformat()}.json")
    a = ap.parse_args()

    url = latest_metadata_url("Salmonella")
    print(f"streaming {url}")
    r = urllib.request.urlopen(url, timeout=900)
    cols = r.readline().decode("utf8", "replace").rstrip("\n").split("\t")
    ix = {c: i for i, c in enumerate(cols)}
    need = ("serovar", "computed_types", "asm_acc", "collected_by", "sra_center",
            "bioproject_acc", "collection_date", "biosample_acc")
    if any(k not in ix for k in need):
        print(f"PD schema changed; missing {[k for k in need if k not in ix]}", file=sys.stderr)
        return 2

    funnel = collections.Counter()
    eligible: list[dict] = []
    tool_failed_but_labelled = 0

    for line in r:
        f = line.decode("utf8", "replace").rstrip("\n").split("\t")
        if len(f) < len(cols):
            continue
        funnel["rows"] += 1
        if funnel["rows"] > a.max_rows:
            break

        serovar_raw = f[ix["serovar"]].strip()
        computed = f[ix["computed_types"]].strip()
        asm = f[ix["asm_acc"]].strip()

        sv = norm_serovar(serovar_raw)
        if not sv or sv in NON_SPECIFIC:
            continue
        funnel["has_serovar"] += 1
        if not asm or asm.upper() in ("NULL", "NA"):
            continue
        funnel["has_assembly"] += 1

        # Positive evidence of independence: the tool failed here, the label exists anyway.
        if any(m in computed for m in TOOL_FAILED_MARKERS):
            tool_failed_but_labelled += 1

        if not is_reference_lab(f[ix["collected_by"]], f[ix["sra_center"]]):
            continue
        funnel["reference_lab"] += 1

        eligible.append({
            "asm_acc": asm,
            "biosample": f[ix["biosample_acc"]].strip(),
            "serovar_label": sv,
            "serovar_raw": serovar_raw,
            "computed_serotype": parse_computed_serotype(computed),
            "collected_by": f[ix["collected_by"]].strip(),
            "bioproject": f[ix["bioproject_acc"]].strip() or "NO_BP",
            "collection_date": f[ix["collection_date"]].strip(),
        })

    print(f"\nfunnel: {dict(funnel)}")
    print(f"  tool FAILED to type but a serovar label exists: {tool_failed_but_labelled} "
          f"(positive evidence some labels are not tool-derived)")
    if not eligible:
        print("no eligible isolates", file=sys.stderr)
        return 2

    # ---- select: cap per serovar, then spread across BioProjects -------------------------------
    rng = random.Random(a.seed)
    rng.shuffle(eligible)
    by_sv: collections.Counter = collections.Counter()
    by_bp: collections.Counter = collections.Counter()
    cohort: list[dict] = []
    # Two passes: the first respects a per-source cap so no single BioProject dominates by accident;
    # the second fills any shortfall without it, since a short cohort is worse than a concentrated one
    # PROVIDED the concentration is measured and reported (which it is, below).
    bp_cap = max(2, a.target // 8)
    for allow_bp_overflow in (False, True):
        for e in eligible:
            if len(cohort) >= a.target:
                break
            if by_sv[e["serovar_label"]] >= a.per_serovar_cap:
                continue
            if not allow_bp_overflow and by_bp[e["bioproject"]] >= bp_cap:
                continue
            if any(c["asm_acc"] == e["asm_acc"] for c in cohort):
                continue
            cohort.append(e)
            by_sv[e["serovar_label"]] += 1
            by_bp[e["bioproject"]] += 1

    src = collections.Counter(c["bioproject"] for c in cohort)
    top, ntop = src.most_common(1)[0]
    share = ntop / len(cohort)
    passes = share <= MAX_SOURCE_SHARE

    print(f"\ncohort: {len(cohort)} isolates | {len(by_sv)} distinct serovars | "
          f"{len(src)} BioProjects | largest share {share:.3f} "
          f"{'PASSES' if passes else 'FAILS'} the {MAX_SOURCE_SHARE:.0%} bar")
    print("  top serovars:", ", ".join(f"{k}({v})" for k, v in by_sv.most_common(8)))

    a.out_tsv.parent.mkdir(parents=True, exist_ok=True)
    with open(a.out_tsv, "w", encoding="utf-8", newline="") as fh:
        fh.write("accession\tmeasured_label\n")
        for c in cohort:
            fh.write(f"{c['asm_acc']}\t{c['serovar_label']}\n")

    out = {
        "schema": "salmserovar-cohort-v1",
        "built": _date.today().isoformat(),
        "source": url,
        "label_field": "NCBI-PD `serovar` (submitter-provided)",
        "comparator_field": "NCBI-PD `computed_types` (NCBI in-silico serotype)",
        "funnel": dict(funnel),
        "independence_evidence": {
            "tool_failed_but_serovar_present": tool_failed_but_labelled,
            "argument": "in these rows the in-silico caller returned an unresolved serotype while a "
                        "serovar label exists, so that label cannot have come from that tool.",
        },
        "filters": {
            "reference_labs": list(REFERENCE_LABS),
            "reference_lab_is_a_judgment": True,
            "per_serovar_cap": a.per_serovar_cap,
            "source_diversity_bar": MAX_SOURCE_SHARE,
        },
        "cohort": {
            "n": len(cohort), "n_distinct_serovars": len(by_sv), "n_bioprojects": len(src),
            "largest_source_share": share, "passes_source_diversity_bar": passes,
            "serovar_counts": dict(by_sv.most_common()),
        },
        "isolates": cohort,
        "honest_limits": [
            "It CANNOT be proven that any individual serovar string came from slide agglutination. "
            "The reference-lab filter is a provenance JUDGMENT, not a measurement.",
            "Residual circularity is BOUNDED, not eliminated: the scorer reports our caller AND the "
            "in-silico comparator against the SAME labels, so contamination inflates both equally and "
            "the DELTA between them survives it even where the absolute number is optimistic.",
            "Restricting to reference labs biases toward clinical/regulatory isolates -- serovar "
            "distribution here is NOT the natural distribution.",
            "A per-serovar cap deliberately distorts prevalence to avoid a cohort that is 90% "
            "Typhimurium/Enteritidis; accuracy here is per-isolate on a flattened mix, not "
            "population-weighted accuracy.",
        ],
    }
    a.out_json.write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")
    print(f"\nwrote {a.out_tsv}\nwrote {a.out_json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
