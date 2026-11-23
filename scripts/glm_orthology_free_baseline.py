"""Is the real-held-out copying baseline inflated by ORTHOLOGY? The one-sided test.

THE QUESTION. `glm_copying_penalty.py` measured MG1655 held-out windows sharing a median 55 bp
verbatim with the training corpus (max 148 of 150), against generated sequence at 15.5 bp and uniform
chance at 15.0. Two readings of that gap:

 (a) real promoter sequence genuinely carries shared content (conserved motifs, paralogues, repeated
     elements) that a good generator should reproduce -> our generator is DEFICIENT; or
 (b) the overlap exists because the corpus contains a DIFFERENT E. COLI STRAIN
     (`GCA_059754755.1`, draft) while scoring uses MG1655 (`GCF_000005845.2`), so the held-out
     windows are ORTHOLOGUES of training windows -> at-chance is the CORRECT behaviour for a
     non-memorising generator, and reproducing that overlap would BE memorisation.

THE TEST. Restrict the real baseline to MG1655 windows whose gene has NO counterpart in the draft
assembly. Those windows have no orthologue anywhere in the corpus by construction. If their overlap
falls to chance, (b) holds. If it stays high, (a) holds.

IT IS ONE-SIDED, AND THAT IS STATED UP FRONT. The draft's `gene=` coverage is 80.7% (4,038 of 5,003
CDS), so "MG1655-only by symbol" is an UPPER bound on strain-specific -- some of these genes are
present in the draft but unannotated. That contamination can only ADD orthologues to this set, which
inflates its overlap. So a fall to chance is decisive for (b), while a high value is ambiguous
between (a) and residual contamination.

THE JOIN IS THE DOCUMENTED TRAP. `UpstreamWindow` carries `gene_id`, which comes from GFF3 `ID=`/
`Name=` and is STRAIN-UNIQUE by construction (`cds-NP_414542.1`). `CLAUDE.md` records that using it
as a cross-strain identifier causes 0% vocabulary overlap. So the cross-strain comparison is done on
`gene=` (the symbol), reached via a WITHIN-MG1655 join from gene_id -> symbol, which is safe because
both come from the same annotation file.
"""
from __future__ import annotations

import json
import re
import statistics
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from dna_decode.glm.corpus import build_natural_set  # noqa: E402
from scripts.glm_copying_penalty import LONG_BP, SEP, longest_shared_substring  # noqa: E402

MG1655 = ROOT / "data" / "cache" / "refseq" / "GCF_000005845.2"
DRAFT_GFF = ROOT / "data" / "refseq_cache" / "GCA_059754755.1" / "annotations.gff3"
CORPUS = ROOT / "wiki" / "glm_training_corpus_2026-11-21.txt"
OUT = ROOT / "wiki" / "glm_orthology_free_baseline_2026-11-22.json"

#: Measured by `glm_copying_penalty.py` on the same corpus, n=120 each. Carried so the comparison is
#: against the SAME reference rather than a re-derived one.
MEASURED_GENERATED_MEDIAN = 15.5
MEASURED_UNIFORM_MEDIAN = 15.0
MEASURED_REAL_ALL_MEDIAN = 55.0

_ATTR = re.compile(r"(?:^|;)(ID|Name|locus_tag|gene)=([^;]+)")


def cds_attr_rows(path: Path):
    with path.open(encoding="utf-8", errors="replace") as f:
        for ln in f:
            if ln.startswith("#"):
                continue
            p = ln.rstrip("\n").split("\t")
            if len(p) < 9 or p[2] != "CDS":
                continue
            yield dict(_ATTR.findall(p[8]))


def symbols(path: Path) -> set[str]:
    return {a["gene"] for a in cds_attr_rows(path) if a.get("gene")}


def id_to_symbol(path: Path) -> dict[str, str]:
    """Within-strain map. Keyed on every identifier form a parser might have stored."""
    out: dict[str, str] = {}
    for a in cds_attr_rows(path):
        sym = a.get("gene")
        if not sym:
            continue
        for k in ("ID", "Name", "locus_tag"):
            v = a.get(k)
            if v:
                out[v] = sym
                # parse_gff3 may strip a `cds-` / `gene-` prefix; register the bare form too.
                if "-" in v:
                    out[v.split("-", 1)[1]] = sym
    return out


def main() -> int:
    mg_gff = MG1655 / "annotations.gff3"
    for p in (mg_gff, DRAFT_GFF, CORPUS):
        if not p.exists():
            print(f"REFUSING: missing {p}", file=sys.stderr)
            return 2

    mg_syms, draft_syms = symbols(mg_gff), symbols(DRAFT_GFF)
    only = mg_syms - draft_syms
    print(f"MG1655 symbols {len(mg_syms)} | draft {len(draft_syms)} | MG1655-ONLY {len(only)}")

    corpus = [ln.strip().upper() for ln in CORPUS.open(encoding="utf-8")
              if ln.strip() and not ln.startswith("#")]
    cset = set(corpus)
    haystack = SEP.join(corpus)

    wins, _ = build_natural_set(MG1655, length=150, limit=None, seed=0, kind="promoter")
    m = id_to_symbol(mg_gff)
    resolved = sum(1 for w in wins if m.get(w.gene_id))
    print(f"windows {len(wins)} | gene_id->symbol resolved {resolved} "
          f"({resolved / max(len(wins), 1):.1%})")
    if resolved < 0.5 * len(wins):
        print("REFUSING: the within-strain gene_id->symbol join resolved under half the windows, so "
              "the orthology filter would be selecting on join failure rather than on orthology.",
              file=sys.stderr)
        return 3

    # Held out = absent from the corpus verbatim. Orthology-free = gene has no draft counterpart.
    held = [w for w in wins if w.seq.upper() not in cset]
    of = [w for w in held if (m.get(w.gene_id) or "") in only]
    print(f"held-out {len(held)} | orthology-free {len(of)}")
    if len(of) < 30:
        print(f"REFUSING: only {len(of)} orthology-free held-out windows; too few to read a median "
              "against a 120-sequence reference.", file=sys.stderr)
        return 4

    sample = [w.seq.upper() for w in of][:120]
    print(f"profiling {len(sample)} orthology-free windows ...", flush=True)
    t0 = time.time()
    lens = []
    for i, s in enumerate(sample):
        lens.append(longest_shared_substring(s, haystack))
        if (i + 1) % 25 == 0:
            print(f"  {i + 1}/{len(sample)} ({time.time() - t0:.0f}s)", flush=True)
    lens.sort()
    med = statistics.median(lens)

    # One-sided reading -- but the MEDIAN alone cannot carry it. A longest-match distribution is
    # driven by its tail, and this one is BIMODAL: most windows at chance, a minority sharing nearly
    # the whole window. A median-only rule declared DECISIVE on exactly that shape, which is the same
    # defect as a knife-edge threshold -- one statistic deciding a qualitative state. The decisive
    # reading therefore requires the median AND the tail to both sit at chance.
    u, r = MEASURED_UNIFORM_MEDIAN, MEASURED_REAL_ALL_MEDIAN
    tail = sum(1 for x in lens if x >= LONG_BP) / len(lens)
    med_at_chance = (med - u) <= 0.25 * (r - u)
    tail_at_chance = tail <= 0.05
    if med_at_chance and tail_at_chance:
        v = "ORTHOLOGY_EXPLAINS_THE_BASELINE"
        why = ("Orthology-free real sequence falls to chance on BOTH median and tail, so the 55 bp "
               "all-windows baseline was cross-strain orthology rather than a property real "
               "promoters carry. A generator at chance is behaving CORRECTLY against this "
               "reference; reproducing the 55 bp would be memorisation. Decisive for this reading, "
               "because contamination could only have pushed the other way.")
    elif med_at_chance:
        v = "ORTHOLOGY_EXPLAINS_THE_MAJORITY_RESIDUAL_TAIL"
        why = (f"The MEDIAN falls to chance ({med} bp vs {u}) so orthology explains the baseline for "
               f"most windows -- but {tail:.1%} still share >={LONG_BP} bp (mean "
               f"{statistics.fmean(lens):.1f} bp, max {max(lens)}), so the distribution is BIMODAL "
               "and the claim does not extend to all real sequence. Most plausibly residual "
               "UNANNOTATED orthologues, since draft gene= coverage is 80.7% and this set is an "
               "upper bound on strain-specific -- but that is not proven, and a generator producing "
               "0% long shares does not reproduce this subpopulation whatever its cause. Report the "
               "median result and this tail together; neither alone is the finding.")
    elif med >= u + 0.5 * (r - u):
        v = "REAL_SEQUENCE_CARRIES_SHARED_CONTENT"
        why = ("Orthology-free real sequence still shares materially more than chance, so real "
               "promoters do carry shared content a generator could be expected to reproduce. "
               "AMBIGUOUS: residual unannotated orthologues (draft gene= coverage is 80.7%) could "
               "produce this same result, so it does not establish generator deficiency.")
    else:
        v = "INTERMEDIATE_UNRESOLVED"
        why = ("Orthology-free overlap sits between chance and the all-windows baseline without "
               "reaching either band. Partial orthology inflation is indicated; the generator "
               "question is not settled either way by this run.")

    res = {
        "record": "glm-orthology-free-baseline-v1",
        "date": "2026-11-22",
        "question": "is the real-held-out copying baseline inflated by cross-strain orthology?",
        "strains": {"scored": "GCF_000005845.2 (MG1655)",
                    "in_corpus": "GCA_059754755.1 (draft)"},
        "symbols": {"mg1655": len(mg_syms), "draft": len(draft_syms), "mg1655_only": len(only),
                    "draft_gene_coverage_note": "draft gene= coverage 80.7%, so MG1655-only is an "
                                                "UPPER bound on strain-specific; contamination "
                                                "inflates this baseline, making the test one-sided"},
        "windows": {"total": len(wins), "held_out": len(held), "orthology_free": len(of),
                    "profiled": len(sample), "gene_id_resolved_frac": round(resolved / len(wins), 4)},
        "orthology_free_real": {
            "n": len(lens), "median_bp": med,
            "mean_bp": round(statistics.fmean(lens), 2), "max_bp": max(lens),
            "frac_at_or_above_long": round(sum(1 for x in lens if x >= LONG_BP) / len(lens), 4),
            "long_bp_threshold": LONG_BP, "seconds": round(time.time() - t0, 1),
        },
        "reference_medians": {"generated": MEASURED_GENERATED_MEDIAN,
                              "uniform_chance": u, "real_all_windows": r},
        "verdict": v, "verdict_detail": why,
    }
    print("\n" + "=" * 74)
    print(f"orthology-free real : median {med} bp  (mean "
          f"{res['orthology_free_real']['mean_bp']}, max {res['orthology_free_real']['max_bp']}, "
          f">={LONG_BP}bp {res['orthology_free_real']['frac_at_or_above_long']:.1%})")
    print(f"real, ALL windows   : median {r} bp")
    print(f"generated           : median {MEASURED_GENERATED_MEDIAN} bp")
    print(f"uniform chance      : median {u} bp")
    print("=" * 74)
    print(f"VERDICT: {v}\n  {why}")
    OUT.write_text(json.dumps(res, indent=2), encoding="utf-8")
    print(f"\n-> {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
