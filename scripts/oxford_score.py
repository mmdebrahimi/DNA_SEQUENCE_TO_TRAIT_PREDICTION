"""Score the FROZEN decoder on the Oxford cohort's OWN AMRFinder output + measured MIC.

The Oxford repo (github.com/samlipworth/ecoli_mic_arg) deposits everything keyed by a
study `guuid`: `main_data.csv` (per-drug log2 MIC interval) + `amrfinder.tsv` (standard
AMRFinderPlus output, Name=guuid). There is NO guuid<->BioSample bridge, but none is
needed: we JOIN the two by guuid and apply the shipped `call_resistance` rule directly.

INDEPENDENCE (honest): this is MORE independent than the NCBI-PD provdisjoint cells — the
MIC is Oxford's own broth microdilution AND the genotype is Oxford's own AMRFinder run
(NOT the frozen pipeline's Docker AMRFinder). The single methodology caveat is the
AMRFinder VERSION/DB difference; the decoder's CORE (the call_resistance DRUG_RULE) is
what's under test, and it operates on standard AMRFinder vocab regardless of version.

MIC: the observed dilution MIC = 2**upper (paper encodes MIC=8 as lower=log2(4), upper=log2(8)).
Fed to the existing classify_tier path via external_mic_labels (strict HIGH_R/HIGH_S primary,
relaxed +DECISIVE). No genome download, no Docker — pure join + the frozen rule.
"""
from __future__ import annotations

import csv
import json
import sys
import tempfile
from datetime import date as _date
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from dna_decode.data.external_mic_labels import build_drug_labels
from dna_decode.data.mic_tiers import breakpoints_for
from dna_decode.eval.amr_rules import call_resistance
from scripts.independent_cohort_validate import _conf


def build_binary_labels(rows, drug: str, col: str) -> dict[str, str]:
    """Standard clinical binary call for CLEAN measured MIC (MIC = 2**upper):
    R if MIC >= CLSI-R breakpoint, S if MIC <= CLSI-S, intermediate excluded.
    Appropriate for gold-standard measured MIC (the tier framework gates NOISY labels)."""
    bp = breakpoints_for(drug)
    r_bp, s_bp = bp["clsi_r"], bp["clsi_s"]
    out = {}
    for r in rows:
        up = (r.get(f"{col}_upper") or "").strip()
        if up in ("", "NA"):
            continue
        mic = 2.0 ** float(up)
        if mic >= r_bp:
            out[r["guuid"].strip()] = "R"
        elif mic <= s_bp:
            out[r["guuid"].strip()] = "S"
    return out

PILOT_DRUGS = {"ciprofloxacin": "Ciprofloxacin", "ceftriaxone": "Ceftriaxone", "gentamicin": "Gentamicin"}
REGISTRY_ORGANISM = "Escherichia_coli_Shigella"
AMR_HEADER = None  # set from the file


def load_mic_labels(main_csv: Path) -> tuple[dict, dict, int]:
    """Per pilot drug -> (tiered build_drug_labels, binary R/S labels). MIC = 2**upper."""
    rows = list(csv.DictReader(open(main_csv, encoding="utf-8")))
    tiered, binary = {}, {}
    for drug, col in PILOT_DRUGS.items():
        iso_mics = {}
        for r in rows:
            up = (r.get(f"{col}_upper") or "").strip()
            if up in ("", "NA"):
                continue
            try:
                iso_mics[r["guuid"].strip()] = [str(2.0 ** float(up))]
            except ValueError:
                continue
        tiered[drug] = build_drug_labels(iso_mics, drug)
        binary[drug] = build_binary_labels(rows, drug, col)
    return tiered, binary, len(rows)


# The paper used an OLDER AMRFinderPlus schema; the FROZEN decoder (amr_rules.py) expects
# the NEWER one (Docker ncbi/amr:4.2.7). Normalize the header labels so the frozen parser
# reads the same columns. Column ORDER + row values are untouched; only labels are renamed.
# (Frozen amr_rules.py stays byte-unchanged — the adaptation lives here, in the adapter.)
_HEADER_RENAME = {
    "Gene symbol": "Element symbol",
    "Sequence name": "Element name",
    "% Identity to reference sequence": "% Identity to reference",
    "% Coverage of reference sequence": "% Coverage of reference",
}


def group_amrfinder(amr_tsv: Path) -> tuple[str, dict[str, list[str]]]:
    """Group raw AMRFinder lines by guuid (col 'Name'). Returns (NORMALIZED header, {guuid:[lines]})."""
    by_guuid: dict[str, list[str]] = {}
    with open(amr_tsv, encoding="utf-8") as fh:
        raw_header = fh.readline().rstrip("\n")
        header = "\t".join(_HEADER_RENAME.get(c, c) for c in raw_header.split("\t"))
        for ln in fh:
            ln = ln.rstrip("\n")
            if not ln:
                continue
            g = ln.split("\t", 1)[0].strip()
            by_guuid.setdefault(g, []).append(ln)
    return header, by_guuid


def main() -> int:
    base = Path("data/raw/oxford")
    mic_labels, binary_labels, n_total = load_mic_labels(base / "main_data.csv")
    header, amr = group_amrfinder(base / "amrfinder.tsv")
    print(f"loaded {n_total} MIC rows; {len(amr)} guuids with AMRFinder hits")

    tmpdir = Path(tempfile.mkdtemp(prefix="oxford_amr_"))
    pred_cache: dict[tuple[str, str], str] = {}

    def predict(guuid: str, drug: str) -> str:
        key = (guuid, drug)
        if key in pred_cache:
            return pred_cache[key]
        p = tmpdir / f"{guuid}.tsv"
        if not p.exists():
            p.write_text(header + "\n" + "\n".join(amr.get(guuid, [])) + "\n", encoding="utf-8")
        pred = call_resistance(p, drug, organism=REGISTRY_ORGANISM)["prediction"]
        pred_cache[key] = pred
        return pred

    artifact_cells = {}
    for drug in PILOT_DRUGS:
        res = mic_labels[drug]
        cells = {}
        for tier_name, labelmap in (("strict", res["strict"]), ("relaxed", res["relaxed"]),
                                    ("binary", binary_labels[drug])):
            pairs = [(predict(g, drug), 1 if rs == "R" else 0) for g, rs in labelmap.items()]
            cells[tier_name] = _conf(pairs)
        artifact_cells[drug] = {"strict": cells["strict"], "relaxed": cells["relaxed"],
                                "binary": cells["binary"], "buckets": res["buckets"], "n_total": res["n_total"]}
        b = cells["binary"]
        print(f"\n{drug}: BINARY n={b['n_scored']} acc={b['acc']} sens={b['sens']} spec={b['spec']} "
              f"(R {b['tp']+b['fn']} / S {b['tn']+b['fp']}; abstain {b['abstain']})")
        s = cells["strict"]
        print(f"  strict n={s['n_scored']} acc={s['acc']} sens={s['sens']} spec={s['spec']}")

    out = {
        "_schema": "external-validation-v1",
        "date": _date.today().isoformat(),
        "cohort": "oxford_lipworth_ecoli_mic_arg",
        "organism": REGISTRY_ORGANISM,
        "drugs": artifact_cells,
        "evidence_tier": "external_clinical",
        "independence_tier": ("Oxford OUH measured MIC (broth microdilution) + Oxford's OWN AMRFinder "
                              "genotype, joined by study guuid. MORE independent than NCBI-PD provdisjoint "
                              "(independent genotype caller too). CAVEAT: AMRFinder version/DB differs from "
                              "the frozen Docker pipeline (header-normalized in the adapter; Gene symbol -> "
                              "Element symbol); the decoder's call_resistance RULE is what is validated."),
        "leakage_status": "DISJOINT_PROJECT_LEVEL",
        "leakage_evidence": ("Oxford BioProjects PRJNA604975 + PRJNA1007570 deposit 0 NCBI assemblies "
                             "(reads-only); the decoder's 1039 tuning/validation accessions are 100% "
                             "GCA_/GCF_ assemblies -> an Oxford isolate (no GCA) cannot be in the tuning "
                             "set. Disjoint by accession-type construction."),
        "primary_metric": "binary (clean measured MIC; strict/relaxed tiers over-exclude gold-standard labels)",
        "source": "github.com/samlipworth/ecoli_mic_arg (PRJNA604975 + PRJNA1007570)",
    }
    outp = Path(f"wiki/external_validation_oxford_{_date.today().isoformat()}.json")
    outp.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(f"\nartifact -> {outp}  (leakage_status={out['leakage_status']})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
