"""Score the FROZEN decoder on the 234-isolate Sci Rep 2023 cohort (2nd external cohort).

Source: PMC9829913 (Sci Rep 2023, BioProject PRJNA854358), open CC-BY supplement
MOESM1_ESM.xlsx (committed-out at data/raw/sci234/, gitignored). Everything is keyed by
study `Key`: "Supplementary data 1" = per-isolate GENOTYPE (QRDR POINT mutations gyrA_p*/
parC_p*/parE_p* + all acquired bla/aac/ant/aph alleles), "Supplementary data 2" = per-isolate
measured MIC (broth microdilution) incl CIPROF / GENTAM / CEFOTA.

Like Oxford (scripts/oxford_score.py): a pure JOIN of the cohort's OWN genotype + measured MIC,
scored by the frozen call_resistance RULE. NO download / assembly / Docker — the supplement
ships the genotype. Independence: independent genotype caller (ResFinder/PointFinder-style) +
independent measured MIC; MORE independent than NCBI-PD provdisjoint. Leakage DISJOINT by
construction: PRJNA854358 deposits 0 NCBI assemblies (reads-only); the decoder's tuning set is
100% GCA/GCF -> no 234-cohort isolate can be in it.

CIPRO is scored here (QRDR-count rule needs no Subclass lookup). GENTAM + CEFOTA(ceftriaxone
proxy) need a fam.tsv-based per-gene Subclass resolver (the gent/cef refinement) -> follow-up.
"""
from __future__ import annotations

import json
import re
import sys
import tempfile
from datetime import date as _date
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from dna_decode.data.mic_tiers import breakpoints_for
from dna_decode.eval.amr_rules import call_resistance
from scripts.independent_cohort_validate import _conf

SUPPL = Path("data/raw/sci234/sci234_supplement_MOESM1.xlsx")
REGISTRY_ORGANISM = "Escherichia_coli_Shigella"
_AMR_HDR = "Element symbol\tMethod\tClass\tSubclass\tElement name\t% Identity to reference"


def _present(v) -> bool:
    return v not in (None, 0, "0", "", "-", "NA", "NaN", False)


def _mic(cell) -> float | None:
    s = str(cell)
    for op in ("<=", ">=", ">", "<", "="):
        s = s.replace(op, "")
    s = s.strip()
    try:
        return float(s)
    except ValueError:
        return None


def load_cohort(suppl: Path = SUPPL):
    import openpyxl
    wb = openpyxl.load_workbook(suppl, read_only=True)
    r1 = wb["Supplementary data 1"].iter_rows(values_only=True)
    h1 = list(next(r1))
    ki = h1.index("Key")
    qrdr = [i for i in range(len(h1)) if re.match(r"(gyrA|gyrB|parC|parE)_p", str(h1[i]) or "")]
    geno_qrdr = {}
    for r in r1:
        k = str(r[ki]).strip()
        geno_qrdr[k] = [str(h1[i]).replace("_p", "_", 1) for i in qrdr if i < len(r) and _present(r[i])]
    r2 = wb["Supplementary data 2"].iter_rows(values_only=True)
    h2 = list(next(r2))
    k2 = h2.index("Key")
    mic_cols = {d: h2.index(col) for d, col in (("ciprofloxacin", "CIPROF"),) if col in h2}
    mic_rows = [(str(r[k2]).strip(), r) for r in r2]
    return geno_qrdr, h2, mic_cols, mic_rows


def cipro_labels(mic_rows, ci) -> dict[str, str]:
    bp = breakpoints_for("ciprofloxacin")
    out = {}
    for k, r in mic_rows:
        m = _mic(r[ci]) if ci < len(r) else None
        if m is None:
            continue
        if m >= bp["clsi_r"]:
            out[k] = "R"
        elif m <= bp["clsi_s"]:
            out[k] = "S"
    return out


def score_cipro(geno_qrdr, labels) -> dict:
    tmp = Path(tempfile.mkdtemp(prefix="sci234_"))

    def predict(k):
        rows = [f"{s}\tPOINTX\tQUINOLONE\tQUINOLONE\t{s}\t100" for s in geno_qrdr.get(k, [])]
        p = tmp / f"{k.replace('/', '_')}.tsv"
        p.write_text(_AMR_HDR + "\n" + "\n".join(rows) + "\n", encoding="utf-8")
        return call_resistance(p, "ciprofloxacin", organism=REGISTRY_ORGANISM)["prediction"]

    pairs = [(predict(k), 1 if rs == "R" else 0) for k, rs in labels.items() if k in geno_qrdr]
    return _conf(pairs)


def main() -> int:
    geno_qrdr, h2, mic_cols, mic_rows = load_cohort()
    labels = cipro_labels(mic_rows, mic_cols["ciprofloxacin"])
    conf = score_cipro(geno_qrdr, labels)
    print(f"SCI234 cipro: n={conf['n_scored']} acc={conf['acc']} sens={conf['sens']} spec={conf['spec']} "
          f"(R{conf['tp']+conf['fn']}/S{conf['tn']+conf['fp']})")
    art = {
        "_schema": "external-validation-v1", "date": _date.today().isoformat(),
        "cohort": "sci234_lanellaecoli_PRJNA854358", "organism": REGISTRY_ORGANISM,
        "drugs": {"ciprofloxacin": conf},
        "evidence_tier": "external_clinical",
        "independence_tier": ("Sci Rep 2023 234-isolate E. coli cohort (Spain/Europe): own measured MIC "
                              "(broth microdilution) + own ResFinder/PointFinder-style genotype, joined by study Key. "
                              "Independent genotype caller; the call_resistance RULE is validated. cipro = QRDR-count "
                              "(no Subclass lookup); gent+cef pending a fam.tsv Subclass resolver."),
        "leakage_status": "DISJOINT_PROJECT_LEVEL",
        "leakage_evidence": ("PRJNA854358 deposits 0 NCBI assemblies (reads-only); decoder tuning set is 100% "
                             "GCA/GCF -> a 234-cohort isolate cannot be in it. Disjoint by accession-type construction."),
        "source": "PMC9829913 / Sci Rep 2023 MOESM1_ESM.xlsx (open CC-BY)",
        "note": "Direct supplement score (no assembly/Docker) — the genotype supplement ships QRDR mutations + acquired alleles.",
    }
    out = Path(f"wiki/external_validation_sci234_{_date.today().isoformat()}.json")
    out.write_text(json.dumps(art, indent=2), encoding="utf-8")
    print(f"artifact -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
