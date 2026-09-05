"""Is the shipped MLST caller producing REAL lineages, or plausible-looking noise? Measure it.

`typing:bacteria:mlst` ships CLI-routable and is registered FAITHFUL_TO_TOOL -- checked against the
reference METHOD, never against reality. It became cheap to test the moment the serotype lineage-disjoint
run called sequence types on 400 real E. coli genomes that ALSO carry a wet-lab O:H serotype label.

THE TEST, AND WHY IT NEEDS NO CURATED BIOLOGY. E. coli clonal lineages conserve their O:H antigens, so a
working MLST caller must produce sequence types that are SEROTYPE-PURE. A broken one -- wrong allele
matching, a mis-joined profile table, an off-by-one in the profile lookup -- produces sequence types that
carve the cohort arbitrarily and land at chance. The obvious alternative test would score against
remembered associations ("ST131 is O25:H4"), which means asserting biology from memory as a scoring key.
This instead measures per-ST purity against a SHUFFLE NULL and lets the data speak.

WHY THE NULL IS THE LOAD-BEARING PART. Purity alone is trivially gameable in BOTH directions: a caller
that gave every genome its own ST would score 1.000, and one that pooled everything into a single ST
would score the cohort's modal-serotype frequency. The null therefore SHUFFLES the wet-lab serotype
labels across genomes while PRESERVING THE OBSERVED ST PARTITION -- identical group sizes, identical
serotype frequencies, only the association destroyed. Over-splitting inflates the null exactly as much as
it inflates the observation, so it cannot manufacture a result.

Offline: reads the committed lineage-disjoint checkpoint. No blastn, no network, no Docker.
"""
from __future__ import annotations

import argparse
import collections
import json
import random
import sys
from datetime import date as _date
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

MIN_ST_SIZE = 3          # an ST of size 1 is 100% pure BY CONSTRUCTION and carries no information
N_SHUFFLES = 1000


def load_rows(path: Path) -> list[dict]:
    return [json.loads(ln) for ln in path.read_text(encoding="utf-8").splitlines() if ln.strip()]


def serotype_of(row: dict) -> str | None:
    """The WET-LAB label, not either caller's prediction."""
    o, h = (row.get("O") or "").strip(), (row.get("H") or "").strip()
    return f"{o}:{h}" if o and h else None


def weighted_purity(groups: dict[str, list[str]], top_k: int = 1) -> tuple[float, int, int]:
    """Fraction of genomes carrying one of their own ST's `top_k` modal serotypes.

    `top_k=2` exists because a real lineage can be genuinely BIMODAL -- so a low top-1 purity is
    ambiguous between "the caller mis-grouped these" and "this lineage really has two serotypes". It
    is a MEASUREMENT of that ambiguity, not an interpretation of any particular ST.

    Returns (purity, n_genomes_scored, n_sts_scored).
    """
    hit = tot = n_st = 0
    for _st, sers in groups.items():
        if len(sers) < MIN_ST_SIZE:
            continue
        n_st += 1
        tot += len(sers)
        hit += sum(k for _s, k in collections.Counter(sers).most_common(top_k))
    return (hit / tot if tot else 0.0), tot, n_st


def shuffle_null(groups: dict[str, list[str]], n: int, top_k: int = 1,
                 seed: int = 12345) -> list[float]:
    """Destroy the ST<->serotype association, KEEP the partition shape and the serotype frequencies.

    Every reported statistic gets its OWN null computed the same way -- a top-2 or O-axis-only purity
    compared against the top-1 null would be an unfair comparison that manufactures significance.
    """
    rng = random.Random(seed)
    sizes = [(st, len(v)) for st, v in groups.items()]
    pool = [s for v in groups.values() for s in v]
    out = []
    for _ in range(n):
        rng.shuffle(pool)
        i, redrawn = 0, {}
        for st, k in sizes:
            redrawn[st] = pool[i:i + k]
            i += k
        out.append(weighted_purity(redrawn, top_k=top_k)[0])
    return out


def measure(groups: dict[str, list[str]], label: str, top_k: int = 1) -> dict:
    """One purity statistic with its OWN matched null."""
    obs = weighted_purity(groups, top_k=top_k)[0]
    null = shuffle_null(groups, N_SHUFFLES, top_k=top_k)
    nm, nx = sum(null) / len(null), max(null)
    return {"statistic": label, "top_k": top_k, "observed": obs, "null_mean": nm, "null_max": nx,
            "p_value": (sum(1 for v in null if v >= obs) + 1) / (len(null) + 1),
            "exceeds_null_max": obs > nx}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--checkpoint", type=Path,
                    default=Path("D:/dna_decode_cache/ecoli_sero_asm/lineage_calls.jsonl"))
    ap.add_argument("--out", type=Path, default=ROOT / "wiki" /
                    f"mlst_serotype_purity_{_date.today().isoformat()}.json")
    a = ap.parse_args()

    if not a.checkpoint.exists():
        print(f"checkpoint absent: {a.checkpoint}", file=sys.stderr)
        return 2
    rows = load_rows(a.checkpoint)
    typed = [r for r in rows
             if r.get("status") == "ok" and r.get("st") and r.get("st_complete")
             and serotype_of(r)]
    groups: dict[str, list[str]] = collections.defaultdict(list)
    for r in typed:
        groups[str(r["st"])].append(serotype_of(r))

    obs, n_scored, n_st = weighted_purity(groups)

    # --- NON-VACUITY -----------------------------------------------------------------------------
    # Three ways this could report a number that means nothing.
    distinct_serotypes = len({s for v in groups.values() for s in v})
    if n_st == 0 or n_scored == 0:
        print(f"REFUSING: {n_st} sequence types reach size {MIN_ST_SIZE}, {n_scored} genomes scored. "
              "Nothing was measured.", file=sys.stderr)
        return 3
    if distinct_serotypes < 2:
        print(f"REFUSING: the cohort carries {distinct_serotypes} distinct serotype(s), so purity is "
              "1.0 by construction and the null is too. This cohort cannot test the caller.",
              file=sys.stderr)
        return 3

    null = shuffle_null(groups, N_SHUFFLES)
    null_mean = sum(null) / len(null)
    null_max = max(null)
    p = (sum(1 for v in null if v >= obs) + 1) / (len(null) + 1)

    # Two extra statistics that turn INTERPRETATIONS of the low-purity STs into MEASUREMENTS.
    # (1) top-2: a genuinely bimodal lineage is not a mis-grouping. (2) O-axis only: if the residual
    # impurity is concentrated on the H antigen, the O antigen is cleaner than the combined figure
    # suggests -- measurable WITHOUT asserting anything about what an H- label means.
    o_groups = {st: [s.split(":", 1)[0] for s in v] for st, v in groups.items()}
    extra = [measure(groups, "top2_serotype_purity", top_k=2),
             measure(o_groups, "O_antigen_only_purity", top_k=1)]

    sizes = sorted((len(v) for v in groups.values()), reverse=True)
    frac_in_qualifying = n_scored / len(typed)

    print(f"genomes typed with a wet-lab serotype : {len(typed)}")
    print(f"distinct sequence types               : {len(groups)}  (largest {sizes[0]}, "
          f"{sizes[0]/len(typed):.1%} of cohort)")
    print(f"STs reaching size {MIN_ST_SIZE}                   : {n_st}  "
          f"covering {n_scored} genomes ({frac_in_qualifying:.1%})")
    print(f"distinct wet-lab serotypes            : {distinct_serotypes}")
    print(f"\nweighted serotype purity  observed : {obs:.4f}")
    print(f"                          null mean: {null_mean:.4f}   null max: {null_max:.4f}   "
          f"p = {p:.4f}")

    for m in extra:
        print(f"{m['statistic']:26s} observed {m['observed']:.4f}  null mean {m['null_mean']:.4f}  "
              f"null max {m['null_max']:.4f}  p={m['p_value']:.4f}")

    print("\nlargest sequence types, modal WET-LAB serotype (data, not a scoring key):")
    for st, sers in sorted(groups.items(), key=lambda kv: -len(kv[1]))[:10]:
        if len(sers) < MIN_ST_SIZE:
            continue
        ser, k = collections.Counter(sers).most_common(1)[0]
        print(f"   ST{st:<6s} n={len(sers):3d}  modal {ser:<12s} {k}/{len(sers)} = {k/len(sers):.2f}")

    if obs > null_max:
        verdict = "MLST_RECOVERS_REAL_LINEAGES"
        why = (f"sequence-type groups are serotype-pure at {obs:.4f}, above the LARGEST of "
               f"{N_SHUFFLES} shuffles of the same partition ({null_max:.4f}, mean {null_mean:.4f}). "
               "The caller is recovering real clonal structure, not carving the cohort arbitrarily. "
               "This is a coherence check against a wet-lab label, NOT a correctness check against a "
               "reference MLST implementation.")
    elif p <= 0.05:
        verdict = "MLST_RECOVERS_REAL_LINEAGES_MARGINAL"
        why = (f"purity {obs:.4f} beats the shuffle null (mean {null_mean:.4f}) at p={p:.4f}, but does "
               "not exceed the null's maximum, so the margin is modest.")
    else:
        verdict = "MLST_PURITY_INDISTINGUISHABLE_FROM_CHANCE"
        why = (f"purity {obs:.4f} vs null mean {null_mean:.4f} at p={p:.4f}. The sequence types carry "
               "no more serotype information than a random partition of the same shape -- which is "
               "what a BROKEN allele-matching or profile-lookup step would produce.")
    print(f"\nVERDICT: {verdict}\n  {why}")

    out = {"schema": "mlst-serotype-purity-v1", "date": _date.today().isoformat(),
           "question": ("do the shipped MLST caller's sequence types group genomes that share a WET-LAB "
                        "serotype, more than a random partition of the same shape would?"),
           "why_no_curated_biology": ("scoring against remembered ST<->serotype associations would mean "
                                      "asserting biology from memory as the scoring key; per-ST purity "
                                      "against a shuffle null needs none"),
           "null": ("wet-lab serotype labels shuffled across genomes with the OBSERVED ST partition "
                    "held fixed -- identical group sizes and serotype frequencies, association "
                    "destroyed. Over-splitting inflates null and observation equally, so it cannot "
                    "manufacture a result"),
           "min_st_size": MIN_ST_SIZE, "n_shuffles": N_SHUFFLES,
           "n_genomes_typed_with_serotype": len(typed), "n_distinct_st": len(groups),
           "largest_st_size": sizes[0], "largest_st_fraction": sizes[0] / len(typed),
           "n_st_scored": n_st, "n_genomes_scored": n_scored,
           "fraction_of_cohort_in_qualifying_sts": frac_in_qualifying,
           "n_distinct_serotypes": distinct_serotypes,
           "observed_purity": obs, "null_mean": null_mean, "null_max": null_max, "p_value": p,
           "additional_statistics": extra,
           "st_composition_for_low_purity_sts": {
               st: dict(collections.Counter(v).most_common())
               for st, v in groups.items()
               if len(v) >= MIN_ST_SIZE
               and collections.Counter(v).most_common(1)[0][1] / len(v) < 0.70},
           "top_sts": [{"st": st, "n": len(v),
                        "modal_serotype": collections.Counter(v).most_common(1)[0][0],
                        "modal_fraction": collections.Counter(v).most_common(1)[0][1] / len(v)}
                       for st, v in sorted(groups.items(), key=lambda kv: -len(kv[1]))[:10]
                       if len(v) >= MIN_ST_SIZE],
           "verdict": verdict, "why": why,
           "honest_limits": [
               "This is a COHERENCE check, not a correctness check. It shows the sequence types track "
               "real clonal structure; it does NOT show they carry the SAME NUMBERS a reference MLST "
               "implementation would assign. A caller with a systematically shifted profile table would "
               "pass this and still report wrong ST numbers -- that needs the reference tool installed "
               "and pinned locally, which was NOT done. The cell stays FAITHFUL_TO_TOOL.",
               "One organism (E. coli), one scheme (Achtman 7-locus), one cohort.",
               "Genomes within a sequence type are clonal, so the genome-level p-value overstates "
               "independence. The result does not rest on the p-value -- the observation exceeding the "
               "null's MAXIMUM over 1000 shuffles is the claim.",
               "Serotype labels are NCBI-PD submitter strings; E. coli O:H typing is traditionally "
               "slide agglutination, but per-isolate method is not provable from the metadata.",
               "STs below the size floor are excluded, so this says nothing about the caller's "
               "behaviour on singleton or near-singleton lineages.",
           ]}
    a.out.write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")
    print(f"\nwrote {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
