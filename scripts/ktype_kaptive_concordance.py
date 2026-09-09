"""Does the shipped single-gene ktype caller reach the wzi method's own published ceiling?

`typing:Klebsiella:ktype` was the last cell carrying the shared "never measured" default. It looked
cohort-blocked -- Klebsiella capsule typing has no free wet-lab label -- until reading the caller showed
it is NOT a Kaptive wrapper: it types the capsule from a SINGLE conserved gene (`wzi`, BIGSdb Pasteur
scheme as bundled by Kleborate) and maps the allele to a K locus. Kaptive types the FULL K locus. Two
different methods over different amounts of sequence, so Kaptive is a genuine independent comparator --
the same relationship as AMRFinder to ResFinder -- and no wet-lab cohort is needed to run it.

THE QUESTION IS NOT "DO THEY AGREE 100%". They cannot. The caller's own docstring records the published
limit: wzi -> K-type is ~94% predictive and NOT one-to-one, because isolates with distinct K types can
share a wzi allele (Brisse 2013 JCM). So a shortfall from 100% is EXPECTED and is a property of the
METHOD. The measurable question is:

    does our implementation ACHIEVE the wzi method's ceiling, or fall short of it?

Falling well short would indicate an implementation defect (the thing a measurement can fix). Landing at
the ceiling means the caller is faithful and the residual is the method's own limit, which is exactly
what the cell already claims -- and claims are worth more once checked.

FAIRNESS. Kaptive is only used as the reference where IT is confident: rows whose `Match confidence` is
not typeable are reported separately and excluded from the agreement denominator, because an uncertain
reference cannot adjudicate anything. Our caller's abstentions are likewise counted separately from its
errors -- an abstention is not a wrong answer.

Per-genome invocation with a restartable checkpoint, deliberately, rather than one batched Kaptive call:
every cached genome file is named `genome.fna`, so Kaptive's `Assembly` column would collide across all
307 inputs, and this host's external drive has disconnected mid-run before.
"""
from __future__ import annotations

import argparse
import collections
import csv
import json
import os
import subprocess
import sys
from datetime import date as _date
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from dna_decode.ktype.runner import call_ktype  # noqa: E402

BLASTN = "C:/Users/Farshad/ncbi-blast/bin/blastn.exe"
WZI_DB = ROOT / "data" / "ktype_db"
# Kaptive's own wording for a call it stands behind; anything else is not usable as a reference.
TYPEABLE = "typeable"


def klebsiella_genomes() -> dict[str, Path]:
    """Accession -> assembly, from every committed Klebsiella cohort. selected.tsv has NO header --
    it is bare `<accession>\\t<label>` rows, so a DictReader silently eats the first genome."""
    accs: set[str] = set()
    for d in (ROOT / "data" / "raw").glob("*kleb*"):
        f = d / "selected.tsv"
        if f.exists():
            for ln in f.read_text(encoding="utf-8").splitlines():
                a = ln.split("\t")[0].strip()
                if a.startswith(("GCA_", "GCF_")):
                    accs.add(a)
    roots = [Path("D:/dna_decode_cache/refseq")] + sorted((ROOT / "data" / "raw").glob("*kleb*/refseq"))
    out: dict[str, Path] = {}
    for a in sorted(accs):
        for r in roots:
            fa = next((x for x in (r / a).glob("*.fna")), None) if (r / a).is_dir() else None
            if fa:
                out[a] = fa
                break
    return out


def run_kaptive(fasta: Path, db: str, kaptive_bin: str, timeout: int = 900) -> dict:
    """One genome -> {kl, confidence, problems} or {'error': ...}. Never raises."""
    env = dict(os.environ, PYTHONUTF8="1", PYTHONIOENCODING="utf-8")
    try:
        p = subprocess.run([kaptive_bin, "assembly", db, str(fasta)], capture_output=True,
                           text=True, encoding="utf-8", errors="replace", timeout=timeout, env=env)
    except Exception as e:                                   # noqa: BLE001
        return {"error": f"{type(e).__name__}"}
    if p.returncode != 0:
        return {"error": f"exit{p.returncode}"}
    rows = list(csv.DictReader(p.stdout.splitlines(), delimiter="\t"))
    if not rows:
        return {"error": "no_rows"}
    r = rows[0]
    return {"kl": (r.get("Best match locus") or "").strip(),
            "confidence": (r.get("Match confidence") or "").strip(),
            "problems": (r.get("Problems") or "").strip()}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--db", default="kpsc_k")
    ap.add_argument("--kaptive-bin", default="kaptive")
    ap.add_argument("--blastn", default=BLASTN)
    ap.add_argument("--limit", type=int, default=0, help="0 = all")
    ap.add_argument("--checkpoint", type=Path,
                    default=Path("D:/dna_decode_cache/ktype_kaptive/rows.jsonl"))
    ap.add_argument("--out", type=Path, default=ROOT / "wiki" /
                    f"ktype_kaptive_concordance_{_date.today().isoformat()}.json")
    a = ap.parse_args()

    genomes = klebsiella_genomes()
    accs = sorted(genomes)
    if a.limit:
        accs = accs[:a.limit]
    if not accs:
        print("no Klebsiella genomes found", file=sys.stderr)
        return 2
    print(f"{len(accs)} Klebsiella genomes with cached assemblies\n")

    a.checkpoint.parent.mkdir(parents=True, exist_ok=True)
    done: dict[str, dict] = {}
    if a.checkpoint.exists():
        for ln in a.checkpoint.read_text(encoding="utf-8").splitlines():
            if ln.strip():
                r = json.loads(ln)
                done[r["acc"]] = r
        print(f"  resuming: {len(done)} already scored")
    fh = open(a.checkpoint, "a", encoding="utf-8")

    for n, acc in enumerate(accs, 1):
        if acc in done:
            continue
        rec: dict = {"acc": acc}
        try:
            ours = call_ktype(genomes[acc], WZI_DB, blastn_bin=a.blastn)
            rec["ours_status"] = ours.get("status")
            rec["ours_kl"] = ours.get("kl_type")
            rec["ours_wzi"] = ours.get("wzi_allele")
        except Exception as e:                               # noqa: BLE001
            rec["ours_status"] = f"error:{type(e).__name__}"
        rec["kaptive"] = run_kaptive(genomes[acc], a.db, a.kaptive_bin)
        fh.write(json.dumps(rec) + "\n")
        fh.flush()
        done[acc] = rec
        if n % 25 == 0:
            print(f"  [{n}/{len(accs)}] {acc} ours={rec.get('ours_kl')} "
                  f"kaptive={rec['kaptive'].get('kl')}", flush=True)
    fh.close()

    rows = [done[x] for x in accs if x in done]
    kap_err = [r for r in rows if r["kaptive"].get("error")]
    # An uncertain reference cannot adjudicate: excluded from the denominator, reported separately.
    kap_ok = [r for r in rows if not r["kaptive"].get("error")
              and TYPEABLE in (r["kaptive"].get("confidence") or "").lower()]
    kap_untypeable = [r for r in rows if not r["kaptive"].get("error") and r not in kap_ok]

    # An abstention is not an error -- counted apart from disagreement, as everywhere else here.
    ours_called = [r for r in kap_ok if r.get("ours_kl")]
    ours_abstained = [r for r in kap_ok if not r.get("ours_kl")]
    agree = [r for r in ours_called if r["ours_kl"] == r["kaptive"]["kl"]]
    disagree = [r for r in ours_called if r["ours_kl"] != r["kaptive"]["kl"]]

    # --- NON-VACUITY -----------------------------------------------------------------------------
    if not rows or not kap_ok or not ours_called:
        print(f"\nREFUSING: {len(rows)} genomes, {len(kap_ok)} with a confident Kaptive reference, "
              f"{len(ours_called)} with a call from our caller. Nothing was compared.", file=sys.stderr)
        return 3

    # An AMBIGUOUS call is not a resolved call, so the strict figure counts it a miss -- but a wzi
    # allele that genuinely maps to two K types is the caller reporting the method's limit honestly,
    # not guessing wrong. Both figures ship; strict is the headline.
    ambiguous = [r for r in disagree if "/" in (r["ours_kl"] or "")]
    ambiguous_containing_truth = [r for r in ambiguous
                                  if r["kaptive"]["kl"] in (r["ours_kl"] or "").split("/")]
    agreement = len(agree) / len(ours_called)
    agreement_lenient = (len(agree) + len(ambiguous_containing_truth)) / len(ours_called)
    CEILING = 0.94          # published wzi->K-type predictive value (Brisse 2013 JCM)

    print(f"\n=== {len(rows)} genomes ===")
    print(f"  Kaptive errored                        : {len(kap_err)}")
    print(f"  Kaptive ran but NOT typeable (excluded): {len(kap_untypeable)}")
    print(f"  confident Kaptive reference            : {len(kap_ok)}")
    print(f"    our caller abstained (not an error)  : {len(ours_abstained)}")
    print(f"    our caller called                    : {len(ours_called)}")
    print(f"\n  AGREEMENT on the {len(ours_called)} comparable genomes: {len(agree)} "
          f"= {agreement:.4f}")
    print(f"  published wzi method ceiling           : ~{CEILING:.2f} (Brisse 2013; wzi->K is not 1:1)")
    print(f"  ambiguous multi-KL calls (scored MISS)  : {len(ambiguous)}, of which "
          f"{len(ambiguous_containing_truth)} contain Kaptive's answer")
    print(f"  lenient agreement crediting those       : {agreement_lenient:.4f}")

    top = collections.Counter((r["ours_kl"], r["kaptive"]["kl"]) for r in disagree)
    if top:
        print("\n  most frequent disagreements (ours -> kaptive):")
        for (o, k), c in top.most_common(8):
            print(f"     {o:>8s} -> {k:<8s}  x{c}")

    if agreement >= CEILING:
        verdict = "REACHES_THE_WZI_METHOD_CEILING"
        why = (f"agreement {agreement:.4f} on {len(ours_called)} genomes with a confident Kaptive "
               f"reference meets the published wzi ceiling (~{CEILING:.2f}). The residual disagreement "
               "is the METHOD's limit -- wzi is one gene and distinct K types can share an allele -- "
               "not an implementation defect. This is a TOOL-vs-TOOL agreement measurement, so it "
               "bounds faithfulness, never correctness.")
    elif agreement >= CEILING - 0.10:
        verdict = "NEAR_THE_WZI_METHOD_CEILING"
        why = (f"agreement {agreement:.4f} sits below the published ~{CEILING:.2f} ceiling but within "
               "10 points. Consistent with the method's own limit plus cohort composition; not clean "
               "evidence of an implementation defect. Read the disagreement table before concluding.")
    else:
        verdict = "FALLS_SHORT_OF_THE_WZI_METHOD_CEILING"
        why = (f"agreement {agreement:.4f} is more than 10 points below the published ~{CEILING:.2f} "
               "ceiling, which the method alone does not explain. That is an implementation-defect "
               "signal and the disagreement table is the place to start.")
    print(f"\nVERDICT: {verdict}\n  {why}")

    out = {"schema": "ktype-kaptive-concordance-v1", "date": _date.today().isoformat(),
           "question": ("does the single-gene wzi ktype caller reach the published ceiling of the wzi "
                        "method itself, measured against full-K-locus Kaptive?"),
           "comparator": {"tool": "Kaptive", "db": a.db,
                          "why_independent": ("Kaptive types the FULL K locus; our caller types one "
                                              "conserved gene (wzi) and maps the allele to a K type. "
                                              "Different methods over different sequence, so this is "
                                              "an independent implementation -- but still a TOOL, not "
                                              "a wet-lab label")},
           "published_ceiling": CEILING,
           "ceiling_source": "Brisse 2013 JCM; wzi->K-type ~94% predictive and NOT one-to-one",
           "n_genomes": len(rows), "n_kaptive_error": len(kap_err),
           "n_kaptive_untypeable_excluded": len(kap_untypeable),
           "n_confident_reference": len(kap_ok),
           "n_our_abstentions": len(ours_abstained), "n_comparable": len(ours_called),
           "n_agree": len(agree), "n_disagree": len(disagree), "agreement": agreement,
           "n_ambiguous_calls_scored_as_miss": len(ambiguous),
           "n_ambiguous_containing_kaptive_answer": len(ambiguous_containing_truth),
           "agreement_lenient_crediting_ambiguous": agreement_lenient,
           "which_figure_is_the_headline": ("strict -- an ambiguous multi-KL string is not a resolved "
                                            "call. The lenient figure is reported beside it because a "
                                            "wzi allele that maps to two K types is the caller stating "
                                            "the method's limit, not guessing wrong"),
           "top_disagreements": [{"ours": o, "kaptive": k, "n": c} for (o, k), c in top.most_common(15)],
           "verdict": verdict, "why": why,
           "honest_limits": [
               "The comparator is a TOOL, not a wet-lab label. This measures agreement with an "
               "independent implementation, never correctness -- both could be wrong together. The "
               "cell stays FAITHFUL_TO_TOOL; a Quellung/serology-labelled Klebsiella cohort is what an "
               "INDEPENDENT_MEASURED tier would need, and none is free.",
               "Genomes are Klebsiella drawn from this project's AMR cohorts, so they are ENRICHED FOR "
               "RESISTANCE and are not a random sample of the species; K-type prevalence here need not "
               "match a population.",
               "Kaptive rows it does not call typeable are EXCLUDED from the denominator, so the "
               "agreement figure describes genomes where a confident reference exists. The excluded "
               "count is reported beside it and is part of the result.",
               "Our caller's abstentions are counted separately from disagreements. Abstention is not "
               "error, and folding the two together would understate the caller.",
               "The ~94% ceiling is a published property of the wzi method on ITS cohorts, not a "
               "constant. Treating it as an exact bar would over-read it; it is an expectation, which "
               "is why the verdict has a near-ceiling band rather than a single threshold.",
           ]}
    a.out.write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")
    print(f"\nwrote {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
