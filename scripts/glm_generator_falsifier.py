"""Score GLM sequence generators against real E. coli regulatory sequence. Read-only; writes an artifact.

    uv run python scripts/glm_generator_falsifier.py                  # baselines only, no download
    uv run python scripts/glm_generator_falsifier.py --model auto     # + the real HF generator

**The split that keeps this honest.** The natural corpus is divided into a FIT half and an EVAL half. The
Markov null is fitted on FIT only and every discrimination is against EVAL only, so no generator is ever
compared against sequences it was trained on. Without that, the null would have memorised its own
comparison set and its distinguishability would be flattered.

Needs nothing but the committed MG1655 reference: `data/cache/refseq/GCF_000005845.2/{genome.fna,
annotations.gff3}`. No network, no GPU, no labels.
"""
from __future__ import annotations

import argparse
import json
import random
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
REFSEQ = ROOT / "data" / "cache" / "refseq" / "GCF_000005845.2"
WIKI = ROOT / "wiki"


def build_natural(length: int, limit: int | None, seed: int):
    from dna_decode.data.annotations import parse_gff3
    from dna_decode.glm.corpus import extract_upstream_windows, load_fasta

    genome = load_fasta(REFSEQ / "genome.fna")
    rows = parse_gff3(REFSEQ / "annotations.gff3")   # a pandas DataFrame; as_rows normalises it
    wins, stats = extract_upstream_windows(genome, rows, length=length)
    seqs = [w.seq for w in wins]
    random.Random(seed).shuffle(seqs)
    if limit:
        seqs = seqs[:limit]
    return seqs, stats


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--length", type=int, default=150, help="upstream window length (bp)")
    p.add_argument("--n", type=int, default=400, help="sequences to generate per candidate")
    p.add_argument("--limit-natural", type=int, default=2000)
    p.add_argument("--k", type=int, default=4, help="k-mer size for the discriminator")
    p.add_argument("--markov-order", type=int, default=3)
    p.add_argument("--seed", type=int, default=7)
    p.add_argument("--modes", default="kmer,positional",
                   help="comma-separated feature modes: kmer | positional | both")
    p.add_argument("--bins", type=int, default=10, help="positional bins")
    p.add_argument("--model", default=None,
                   help="HF model id, or 'auto' for the default GENERator prokaryote checkpoint")
    p.add_argument("--device", default="auto", help="auto | cpu | cuda")
    p.add_argument("--cache-dir", default="D:/hf_cache", help="weights cache (keep off C:)")
    p.add_argument("--out", default=None)
    a = p.parse_args(argv)

    from dna_decode.glm.baselines import MarkovGenerator, UniformGenerator
    from dna_decode.glm.falsifier import run_falsifier
    from dna_decode.glm.generate import device_report

    if not (REFSEQ / "genome.fna").exists():
        print(f"REFUSE: reference not found at {REFSEQ}")
        return 2

    natural, cstats = build_natural(a.length, a.limit_natural, a.seed)
    print(f"natural corpus: {len(natural)} windows of {a.length} bp  "
          f"(from {cstats.n_cds_total} CDS; {cstats.n_overlapping_cds} overlap another CDS, "
          f"{cstats.n_skipped_edge} skipped at contig edge)")
    if len(natural) < 200:
        print("REFUSE: natural corpus too small to split into fit/eval halves")
        return 2

    half = len(natural) // 2
    fit_set, eval_set = natural[:half], natural[half:]
    print(f"fit half: {len(fit_set)}   eval half (all discrimination happens here): {len(eval_set)}")

    dev = device_report()
    print(f"device: {dev['device_name'] or 'cpu'}  dtype={dev['recommended_dtype']}"
          + (f"  [{dev['notes'][0]}]" if dev["notes"] else ""))

    uni = UniformGenerator(seed=a.seed).generate(a.n, a.length)
    markov = MarkovGenerator(order=a.markov_order, seed=a.seed).fit(fit_set)
    null_seqs = markov.generate(a.n, a.length)

    candidates = [(markov.name, null_seqs), ("uniform", uni)]

    if a.model:
        model_id = None if a.model == "auto" else a.model
        from dna_decode.glm.generate import DEFAULT_MODEL, HFGenerator
        gen = HFGenerator(model_id or DEFAULT_MODEL, cache_dir=a.cache_dir, seed=a.seed,
                          device=a.device)
        try:
            print(f"loading {gen.model_id} (dtype {gen.dtype}) ...")
            seqs = gen.generate(a.n, a.length)
            if gen.fallback_reason:
                print(f"  NOTE: {gen.fallback_reason}")
            candidates.append((gen.name, seqs))
        except Exception as e:  # noqa: BLE001 - a load/HW failure must not void the baseline result
            print(f"HF generator UNAVAILABLE ({type(e).__name__}: {str(e)[:160]})")
            print("  -> baselines still scored; the artifact records the generator as unavailable")
            candidates.append((f"hf:{(model_id or DEFAULT_MODEL).split('/')[-1]}", []))

    results = []
    for mode in a.modes.split(","):
        mode = mode.strip()
        print(f"\n{'=' * 78}\nFEATURE MODE: {mode}"
              + (f" (bins={a.bins})" if mode in ("positional", "both") else "")
              + f"\n{'=' * 78}")
        for name, seqs in candidates:
            r = run_falsifier(
                eval_set, seqs, candidate_name=name,
                null_sequences=None if name == markov.name else null_seqs,
                uniform_sequences=uni, k=a.k, length=a.length, seed=a.seed,
                mode=mode, bins=a.bins,
            )
            results.append(r.as_dict())
            print(f"\n[{name}]  verdict={r.verdict}")
            print(f"   distinguishability (lower=better): {r.distinguishability}"
                  f"   null: {r.null_distinguishability}")
            print(f"   controls: self-max {r.control_self_max} {r.control_self_values}"
                  f"   uniform {r.control_uniform}")
            print(f"   GC natural {r.gc_natural} vs candidate {r.gc_candidate}")
            print(f"   {r.reason}")

    today = date.today().isoformat()
    art = {
        "record": "glm-generator-falsifier-v1",
        "date": today,
        "reference": "GCF_000005845.2 (E. coli K-12 MG1655), committed locally",
        "window_length_bp": a.length,
        "kmer_k": a.k,
        "feature_modes": a.modes,
        "positional_bins": a.bins,
        "markov_order": a.markov_order,
        "n_generated_per_candidate": a.n,
        "natural_corpus": cstats.as_dict(),
        "fit_eval_split": {"n_fit": len(fit_set), "n_eval": len(eval_set),
                           "note": "Markov null fitted on FIT only; all discrimination on EVAL only, so "
                                   "no generator is scored against sequences it was trained on"},
        "device": dev,
        "results": results,
        "honest_limits": [
            "An upstream window is a PROXY for a promoter, not a mapped transcription start site.",
            "One organism (E. coli K-12). Says nothing about transfer to another host.",
            "A k-mer discriminator measures LOCAL composition; the published audits that separate "
            "generated from natural sequence name LONG-RANGE organisation as the failure, which k-mer "
            "features cannot see. So a low distinguishability here is necessary, NOT sufficient.",
            "dna_decode/constraints is CDS-shaped and INAPPLICABLE to non-coding sequence; it is not run "
            "and no vacuous pass-rate is reported.",
        ],
    }
    out = Path(a.out) if a.out else WIKI / f"glm_generator_falsifier_{today}.json"
    out.write_text(json.dumps(art, indent=2), encoding="utf-8")
    print(f"\nwrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
