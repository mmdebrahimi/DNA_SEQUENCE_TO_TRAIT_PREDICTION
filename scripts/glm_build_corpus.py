"""Assemble the generator's training corpus from the fetched genomes, with the leak gate enforced.

    uv run python scripts/glm_build_corpus.py [--spec PATH] [--seeds 10] [--out-dir DIR]

Reads the species spec plus the four species already cached, builds every genome's promoter windows
through `corpus.build_natural_set`, deduplicates by SEQUENCE across genomes, withholds everything the
committed falsifier could score against, and writes the corpus with a manifest that records why each
genome is in it.

**The withhold set is the union over SCORING seeds, and that is the whole point of the gate.** The
falsifier reshuffles its corpus per seed, so each seed has a different eval half. Withholding only the
committed anchor's half demonstrably leaves other seeds' scoring windows in training, which was measured
rather than assumed. Train on those and the generator is trained on the exam: the numbers come out
excellent and nothing raises an error.

**The withhold is applied by SEQUENCE across every genome, not just the reference.** A conserved promoter
in a different *E. coli* strain can be byte-identical to one in the reference's eval half, so gating only
the reference genome would leak through its relatives. Checking every window against the withheld set
costs nothing and closes that.

Honest scope: this builds a corpus and reports what is in it. Promoter-mode corpora are SMALL by
construction -- yield scales with gene count, not genome size -- so the number this prints is the real
constraint on model size, not a disappointment to be explained away.
"""

from __future__ import annotations

import argparse
import json
import time
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
CACHE = ROOT / "data" / "refseq_cache"
REF = ROOT / "data" / "cache" / "refseq" / "GCF_000005845.2"
DEFAULT_SPEC = Path("C:/Users/Laptop Portable/pythonprojects/dna-generator-research/handoffs"
                    "/bacterial-expansion-species-spec.json")
OUT_DIR = ROOT / "wiki"

LENGTH = 150
KIND = "promoter"
SCORING_LIMIT = 1200        # the committed falsifier geometry, which defines what must be withheld

# The four species already cached, as one representative each. Lineage keys are recorded in the same
# `phylum;class` form the spec uses so the manifest has one vocabulary.
PRIOR_SPECIES = (
    ("Escherichia coli", "Pseudomonadota;Gammaproteobacteria"),
    ("Klebsiella pneumoniae", "Pseudomonadota;Gammaproteobacteria"),
    ("Campylobacter jejuni", "Campylobacterota;Epsilonproteobacteria"),
    ("Campylobacter coli", "Campylobacterota;Epsilonproteobacteria"),
)


def prior_representatives() -> list[tuple[str, Path, str]]:
    """`(species, dir, lineage_key)` for the already-cached species, lexicographically smallest.

    Same blind rule used everywhere else in this arm. These share two lineage keys between them, so
    only the first representative of each key is taken -- the builder refuses duplicates by design and
    that refusal should not be worked around here.
    """
    out, taken = [], set()
    for d in sorted(CACHE.iterdir()):
        if not d.is_dir() or not ((d / "genome.fna").exists() and (d / "annotations.gff3").exists()):
            continue
        try:
            header = (d / "genome.fna").open(encoding="utf-8", errors="replace").readline()
        except OSError:
            continue
        for species, key in PRIOR_SPECIES:
            if species in header and key not in taken:
                taken.add(key)
                out.append((species, d, key))
                break
    return out


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--spec", default=str(DEFAULT_SPEC))
    p.add_argument("--seeds", type=int, default=10,
                   help="number of scoring seeds whose eval halves must be withheld")
    p.add_argument("--out-dir", default=str(OUT_DIR))
    a = p.parse_args(argv)

    from dna_decode.glm.training_corpus import (
        CorpusEntry,
        build_training_corpus,
        domain_balance_report,
        falsifier_eval_windows,
    )

    t0 = time.time()
    entries: list[CorpusEntry] = []
    seen_keys: set[str] = set()

    spec_path = Path(a.spec)
    spec_name = None
    if spec_path.exists():
        rec = json.loads(spec_path.read_text(encoding="utf-8"))
        spec_name = f"{rec.get('spec')} v{rec.get('version', '?')}"
        for g in rec.get("genomes", []):
            acc = g.get("accession")
            if not acc:
                continue
            d = CACHE / acc
            if not ((d / "genome.fna").exists() and (d / "annotations.gff3").exists()):
                print(f"  skip {acc}: not fetched")
                continue
            key = g.get("lineage_key") or acc
            if key in seen_keys:
                print(f"  skip {acc}: lineage {key} already represented")
                continue
            seen_keys.add(key)
            entries.append(CorpusEntry(
                accession=acc, species=g.get("species") or acc, domain=g.get("domain", "bacteria"),
                lineage_key=key, reference_dir=str(d), source="spec",
                note=g.get("ncbi_organism_name") or ""))

    for species, d, key in prior_representatives():
        if key in seen_keys:
            print(f"  skip {d.name}: lineage {key} already represented by the spec")
            continue
        seen_keys.add(key)
        entries.append(CorpusEntry(accession=d.name, species=species, domain="bacteria",
                                   lineage_key=key, reference_dir=str(d), source="prior_cache"))

    if not entries:
        print("REFUSE: no genomes available to build from.")
        return 2
    print(f"genomes: {len(entries)} across {len({e.lineage_key for e in entries})} lineages")

    print(f"computing the withhold set over {a.seeds} scoring seeds...")
    hold = falsifier_eval_windows(REF, length=LENGTH, limit=SCORING_LIMIT,
                                  seeds=tuple(range(a.seeds)), kind=KIND)
    one = falsifier_eval_windows(REF, length=LENGTH, limit=SCORING_LIMIT, seeds=(7,), kind=KIND)
    print(f"  union over {a.seeds} seeds: {len(hold)} sequences "
          f"({len(one)} for the anchor seed alone -- {len(hold) - len(one)} would have leaked)")

    windows, man = build_training_corpus(entries, length=LENGTH, kind=KIND, withhold=hold)
    bal = domain_balance_report(man)

    out_dir = Path(a.out_dir)
    corpus_path = out_dir / "glm_training_corpus_2026-11-21.txt"
    corpus_path.write_text("\n".join(windows) + "\n", encoding="utf-8")

    rec_out = {
        "record": "glm-training-corpus-v1",
        "date": date.today().isoformat(),
        "spec": spec_name,
        "corpus_file": corpus_path.name,
        "geometry": {"window_length_bp": LENGTH, "natural_set_arg": KIND,
                     "scoring_limit_natural": SCORING_LIMIT, "scoring_seeds": a.seeds},
        "leak_gate": {
            "withheld_union_size": len(hold),
            "anchor_seed_only_size": len(one),
            "would_have_leaked_under_single_seed_gate": len(hold) - len(one),
            "why_the_union": ("the falsifier reshuffles per seed, so each seed has a different eval "
                              "half; gating only the anchor leaves the other seeds' scoring windows "
                              "in training, and a generator trained on them produces excellent "
                              "numbers and no error"),
            "applied_across_all_genomes": ("by SEQUENCE, so a conserved promoter that is "
                                           "byte-identical in a relative is withheld too"),
        },
        "manifest": man.as_dict(),
        "domain_balance": bal,
        "honest_limits": [
            "Promoter-mode corpora are SMALL by construction: yield scales with gene count, not genome "
            "size. The corpus size printed here is the real constraint on model size.",
            "One representative per lineage, so within-lineage pangenome content is deliberately not "
            "collected -- diversity comes from species, not strains.",
            "The withhold set is computed from the reference genome the falsifier scores on. A future "
            "scorer using a different reference or geometry would need its own withhold set, and this "
            "corpus would not be gated against it.",
            "Deduplication is exact-sequence. Near-duplicate promoters differing by a base are kept as "
            "distinct, so the corpus is an upper bound on truly independent content.",
            "All bacteria, so the domain balance is trivially one-sided; the field exists for the "
            "cross-domain phase and is reported rather than corrected.",
        ],
        "wall_seconds": round(time.time() - t0, 1),
    }
    man_path = out_dir / "glm_training_corpus_manifest_2026-11-21.json"
    man_path.write_text(json.dumps(rec_out, indent=2), encoding="utf-8")

    d = man.as_dict()
    print()
    print(f"corpus: {d['n_windows_kept']:,} windows = {d['corpus_bp']/1e6:.2f} Mbp")
    print(f"  raw {d['n_windows_raw_total']:,} -> dropped {d['n_windows_duplicate_dropped']:,} "
          f"duplicates ({d['duplicate_fraction']*100:.1f}%), withheld "
          f"{d['n_windows_withheld_as_falsifier_eval']:,} as scoring windows")
    print(f"  {d['n_lineages']} lineages, domains {d['by_domain']}")
    print(f"  largest domain share {bal['largest_share']} "
          f"(single-domain in practice: {bal['single_domain_in_practice']})")
    print(f"\n{rec_out['wall_seconds']}s -> {corpus_path.relative_to(ROOT).as_posix()}")
    print(f"              -> {man_path.relative_to(ROOT).as_posix()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
