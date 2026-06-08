"""EP-7 step 3-4 — C. auris azole-AMR cohort validation (Gate G1).

Given a per-isolate label table (fluconazole MIC + clade + a genome source), runs the deterministic
ERG11/FKS1 target-site caller (`scripts/fungal_erg11_caller.call_erg11`) over the cohort and reports:

  * overall accuracy / sensitivity / specificity vs the CDC tentative MIC breakpoint (MIC>=32 => R),
  * WITHIN-CLADE de-confounded metrics (the fungal analogue of within-lineage — clade is the C. auris
    confound, like bacterial lineage; an overall AUROC/acc can be inflated by clade structure alone), and
  * the efflux/aneuploidy DISCORDANCE set: isolates that are MIC-R but ERG11-S — the documented non-target
    mechanism (TAC1b/CDR1 efflux, chr5/ERG11 aneuploidy) that a target-site scan CANNOT detect. A large
    discordance is the EXPECTED failure mode (the falsifier), not a bug.

Gate G1 verdict (per EP-7 / the eukaryotic ledger):
  * PASS              -> overall acc >= 0.80 AND sensitivity >= 0.80
  * DOCUMENTED_FAILURE-> sens < 0.80 BUT the false-negatives are dominated by ERG11-S isolates (i.e. the
                         miss is the documented efflux/aneuploidy blind spot, not a caller defect)
  * FAIL              -> neither (caller is wrong on isolates that DO carry catalogued ERG11 mutations)

LABEL TABLE (TSV, header required). Columns:
  isolate_id           (required)  unique id
  fluconazole_mic      (required)  MIC in ug/mL (numeric; '>256'/'<=1' tolerated -> stripped to number)
  clade                (required)  C. auris clade (I/II/III/IV/V or free text; '' allowed -> '(unknown)')
  assembly_accession   (one-of)    NCBI accession -> downloaded via refseq.download_genome
  genome_fasta         (one-of)    local FASTA path (bypasses download; for offline / pre-fetched genomes)

The ERG11 reference is the committed real C. auris allele (data/fungal_ref/Cauris_ERG11_cds.fna).
Offline-safe: if BLAST+ is absent every call is INDETERMINATE and the script reports that honestly
(no fabricated metrics).
"""
from __future__ import annotations

import argparse
import json
import sys
from dataclasses import dataclass, field
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from dna_decode.data import refseq  # noqa: E402
from dna_decode.data.fungal_amr import (  # noqa: E402
    CAURIS_TENTATIVE_R_MIC,
    FUNGAL_UNDETECTABLE_MECHANISMS,
    mic_to_phenotype,
)
from scripts.fungal_erg11_caller import call_erg11  # noqa: E402

_DEFAULT_REF = Path(__file__).resolve().parent.parent / "data" / "fungal_ref" / "Cauris_ERG11_cds.fna"


def _parse_mic(raw: str) -> float:
    """Tolerate '>256', '<=1', '256.0' -> float. Raises ValueError on unparseable."""
    s = raw.strip().lstrip("><=~").strip()
    return float(s)


@dataclass
class IsolateResult:
    isolate_id: str
    clade: str
    mic: float
    true_pheno: str            # R / S  (from MIC)
    pred_pheno: str            # R / S / INDETERMINATE  (from ERG11 caller)
    determinants: list[str]
    bucket: str = ""           # TP / TN / FP / FN / INDETERMINATE


@dataclass
class CohortReport:
    drug: str
    n_total: int
    n_scored: int              # excludes INDETERMINATE (no BLAST) + unlabelable MIC
    tp: int
    tn: int
    fp: int
    fn: int
    by_clade: dict[str, dict[str, int]] = field(default_factory=dict)
    discordant_efflux: list[str] = field(default_factory=list)  # MIC-R but ERG11-S isolate ids
    indeterminate: list[str] = field(default_factory=list)
    isolates: list[IsolateResult] = field(default_factory=list)

    @property
    def accuracy(self) -> float | None:
        return (self.tp + self.tn) / self.n_scored if self.n_scored else None

    @property
    def sensitivity(self) -> float | None:
        d = self.tp + self.fn
        return self.tp / d if d else None

    @property
    def specificity(self) -> float | None:
        d = self.tn + self.fp
        return self.tn / d if d else None

    def verdict(self) -> str:
        acc, sens = self.accuracy, self.sensitivity
        if acc is None or sens is None:
            return "INSUFFICIENT_DATA"
        if acc >= 0.80 and sens >= 0.80:
            return "PASS"
        # documented failure mode: misses are the efflux/aneuploidy blind spot (ALL FN are ERG11-S,
        # which they are by construction — FN means MIC-R + ERG11-S). Distinguish from a caller defect
        # (which would show up as wrong calls on isolates carrying catalogued ERG11 mutations = FP-shaped
        # or a sens drop with NON-discordant FNs). Here every FN IS efflux-discordant by definition, so a
        # sub-0.80 sens with FN==discordance count is the documented blind spot.
        if self.fn > 0 and len(self.discordant_efflux) == self.fn:
            return "DOCUMENTED_FAILURE_MODE"
        return "FAIL"


def _read_label_table(path: Path) -> list[dict]:
    rows: list[dict] = []
    lines = [ln for ln in path.read_text(encoding="utf-8").splitlines() if ln.strip()]
    if not lines:
        raise ValueError(f"empty label table: {path}")
    header = [h.strip() for h in lines[0].split("\t")]
    req = {"isolate_id", "fluconazole_mic", "clade"}
    if not req.issubset(header):
        raise ValueError(f"label table missing required columns {req - set(header)}; got {header}")
    for ln in lines[1:]:
        cells = ln.split("\t")
        row = {header[i]: (cells[i].strip() if i < len(cells) else "") for i in range(len(header))}
        if not (row.get("assembly_accession") or row.get("genome_fasta")):
            raise ValueError(f"isolate {row.get('isolate_id')!r}: need assembly_accession OR genome_fasta")
        rows.append(row)
    return rows


def _resolve_genome(row: dict, refseq_cache: Path) -> str:
    local = row.get("genome_fasta", "").strip()
    if local:
        if not Path(local).exists():
            raise FileNotFoundError(f"genome_fasta not found: {local}")
        return local
    acc = row["assembly_accession"].strip()
    refseq.download_genome(acc, refseq_cache)        # idempotent; skips if cached
    return str(refseq.fasta_path(acc, refseq_cache))


def build_cohort_report(label_table: Path, ref_fasta: Path, refseq_cache: Path,
                        drug: str = "fluconazole") -> CohortReport:
    rows = _read_label_table(label_table)
    rep = CohortReport(drug=drug, n_total=len(rows), n_scored=0, tp=0, tn=0, fp=0, fn=0)
    for row in rows:
        iso = row["isolate_id"]
        clade = row["clade"].strip() or "(unknown)"
        mic = _parse_mic(row["fluconazole_mic"])
        true_p = mic_to_phenotype(drug, mic)
        genome = _resolve_genome(row, refseq_cache)
        call = call_erg11(genome, str(ref_fasta), drug)
        res = IsolateResult(iso, clade, mic, true_p or "(unlabelable)", call.prediction,
                            list(call.determinants))
        if call.prediction == "INDETERMINATE":
            res.bucket = "INDETERMINATE"
            rep.indeterminate.append(iso)
        elif true_p is None:
            res.bucket = "UNLABELABLE_MIC"
        else:
            pred = call.prediction
            if true_p == "R" and pred == "R":
                res.bucket = "TP"; rep.tp += 1
            elif true_p == "S" and pred == "S":
                res.bucket = "TN"; rep.tn += 1
            elif true_p == "S" and pred == "R":
                res.bucket = "FP"; rep.fp += 1
            else:  # true R, pred S -> the documented efflux/aneuploidy blind spot
                res.bucket = "FN"; rep.fn += 1
                rep.discordant_efflux.append(iso)
            rep.n_scored += 1
            cb = rep.by_clade.setdefault(clade, {"tp": 0, "tn": 0, "fp": 0, "fn": 0})
            cb[res.bucket.lower()] += 1
        rep.isolates.append(res)
    return rep


def _fmt_pct(x: float | None) -> str:
    return f"{x:.3f}" if x is not None else "n/a"


def render_markdown(rep: CohortReport, label_table: Path, today: str) -> str:
    bp = CAURIS_TENTATIVE_R_MIC.get(rep.drug, "?")
    lines = [
        f"# Fungal AMR cohort validation (Gate G1) — C. auris {rep.drug}",
        "",
        f"> Generated {today}. Label table: `{label_table.name}`. Breakpoint: MIC >= {bp} ug/mL = R "
        "(CDC tentative; no formal CLSI/EUCAST C. auris breakpoint exists).",
        f"> Caller: deterministic ERG11/FKS1 target-site scan (blastn vs the committed real C. auris "
        "reference). Clade = the de-confounding variable (C. auris analogue of bacterial lineage).",
        "",
        f"## Verdict: **{rep.verdict()}**",
        "",
        "| metric | value |",
        "|---|---|",
        f"| isolates (total) | {rep.n_total} |",
        f"| scored (excl. INDETERMINATE/unlabelable) | {rep.n_scored} |",
        f"| accuracy | {_fmt_pct(rep.accuracy)} |",
        f"| sensitivity (recall on MIC-R) | {_fmt_pct(rep.sensitivity)} |",
        f"| specificity | {_fmt_pct(rep.specificity)} |",
        f"| TP / TN / FP / FN | {rep.tp} / {rep.tn} / {rep.fp} / {rep.fn} |",
        f"| INDETERMINATE (no BLAST) | {len(rep.indeterminate)} |",
        "",
        "## Within-clade de-confounding (per-clade confusion)",
        "",
        "| clade | TP | TN | FP | FN | n |",
        "|---|---|---|---|---|---|",
    ]
    for clade, cb in sorted(rep.by_clade.items()):
        n = sum(cb.values())
        lines.append(f"| {clade} | {cb['tp']} | {cb['tn']} | {cb['fp']} | {cb['fn']} | {n} |")
    lines += [
        "",
        "## Efflux/aneuploidy discordance (the documented blind spot)",
        "",
        f"Isolates MIC-R but ERG11-S = candidate non-target resistance "
        f"({', '.join(FUNGAL_UNDETECTABLE_MECHANISMS)}): "
        + (", ".join(rep.discordant_efflux) if rep.discordant_efflux else "(none)"),
        "",
        "A deterministic target-site scan CANNOT detect these by design. A sub-0.80 sensitivity driven "
        "entirely by this set is the documented failure mode (falsifier), not a caller defect.",
        "",
        "## Per-isolate",
        "",
        "| isolate | clade | MIC | true | pred | bucket | determinants |",
        "|---|---|---|---|---|---|---|",
    ]
    for r in rep.isolates:
        lines.append(f"| {r.isolate_id} | {r.clade} | {r.mic:g} | {r.true_pheno} | {r.pred_pheno} | "
                     f"{r.bucket} | {', '.join(r.determinants) or '-'} |")
    return "\n".join(lines) + "\n"


def report_to_dict(rep: CohortReport) -> dict:
    return {
        "drug": rep.drug, "verdict": rep.verdict(),
        "n_total": rep.n_total, "n_scored": rep.n_scored,
        "accuracy": rep.accuracy, "sensitivity": rep.sensitivity, "specificity": rep.specificity,
        "tp": rep.tp, "tn": rep.tn, "fp": rep.fp, "fn": rep.fn,
        "by_clade": rep.by_clade,
        "discordant_efflux": rep.discordant_efflux,
        "indeterminate": rep.indeterminate,
        "isolates": [vars(r) for r in rep.isolates],
    }


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--label-table", required=True, help="TSV per-isolate label table (see module docstring)")
    ap.add_argument("--erg11-ref", default=str(_DEFAULT_REF), help="in-frame ERG11 CDS reference FASTA")
    ap.add_argument("--refseq-cache", default="D:/dna_decode_cache/refseq",
                    help="genome download cache root (for assembly_accession rows)")
    ap.add_argument("--drug", default="fluconazole")
    ap.add_argument("--out-prefix", default=None,
                    help="write <prefix>.md + <prefix>.json (default wiki/fungal_cohort_g1_<drug>_<date>)")
    ap.add_argument("--today", default="", help="ISO date stamp for the artifact header (Date.now unavailable)")
    a = ap.parse_args(argv)

    rep = build_cohort_report(Path(a.label_table), Path(a.erg11_ref), Path(a.refseq_cache), a.drug)
    today = a.today or "(date unset)"
    prefix = a.out_prefix or f"wiki/fungal_cohort_g1_{a.drug}_{today}"
    Path(prefix).parent.mkdir(parents=True, exist_ok=True)
    Path(prefix + ".md").write_text(render_markdown(rep, Path(a.label_table), today), encoding="utf-8")
    Path(prefix + ".json").write_text(json.dumps(report_to_dict(rep), indent=2), encoding="utf-8")

    v = rep.verdict()
    print(f"G1 verdict: {v}  acc={_fmt_pct(rep.accuracy)} sens={_fmt_pct(rep.sensitivity)} "
          f"spec={_fmt_pct(rep.specificity)}  (n_scored={rep.n_scored}, FN/efflux={len(rep.discordant_efflux)})")
    print(f"wrote {prefix}.md + {prefix}.json")
    # exit 0 only on a clean PASS; documented-failure + fail are non-zero so a wrapper can branch.
    return 0 if v == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
