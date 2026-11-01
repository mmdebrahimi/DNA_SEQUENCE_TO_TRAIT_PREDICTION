"""Record what is KNOWN about how the PromoterAtlas pretraining corpus was built, with provenance.

The arm's load-bearing honest limit reads: *"the corpus holds 66 Escherichia/Shigella assemblies and
its card calls it regulatory sequence extracted from GenBank -- IF those are upstream-of-gene
regions, resembling the corpus is itself a proxy for being a real promoter."* The conditional was
the open part. This probe resolves it from the authors' own source.

**It emits every claim with an explicit status and refuses to collapse two of them.** A failed fetch
records `fetch_failed`, never `not_stated`: "the source is silent" and "I could not read the source"
are different facts, and this project has collapsed them before. Measured here -- the PMC route
returned a cookie-enablement page and a summarizer duly reported "NOT STATED", which is exactly that
collapse; it is recorded as `fetch_failed`.

Exits 0 always. This is a probe, not a gate.
"""

from __future__ import annotations

import json
import os
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
PA = Path(os.environ.get("PA_REPO", r"C:/Users/LAPTOP~1/AppData/Local/Temp/pa_probe/PromoterAtlas"))
OUT = REPO / "wiki/promoteratlas_corpus_provenance_2026-11-01.json"


def _grab(rel: str, needles: tuple[str, ...]) -> list[dict]:
    """Return the matching lines of a local source file, with line numbers, verbatim."""
    p = PA / rel
    if not p.is_file():
        return [{"status": "source_absent", "path": rel}]
    hits = []
    for i, line in enumerate(p.read_text(encoding="utf-8", errors="replace").splitlines(), 1):
        if any(n in line for n in needles):
            hits.append({"line": i, "text": line.rstrip()})
    return hits


def build() -> dict:
    claims = []

    # ---- CLAIM 1: the extraction rule, read from the authors' own code -------------------------
    g = _grab(
        "src/promoter_atlas/utils/genomics.py",
        ("def extract_promoter_regions", "regulatory_region_length", "UPSTREAM", "prom_seq ="),
    )
    b = _grab(
        "src/promoter_atlas/io/genbank_io.py",
        ("def extract_promoter_regions_genbank", "upstream_length", "prom_seq ="),
    )
    claims.append({
        "claim": "corpus_extraction_rule",
        "status": "quoted" if g and "status" not in g[0] else "fetch_failed",
        "source_url": "https://github.com/LucasCoppens/PromoterAtlas (MIT) — local clone",
        "source_paths": [
            "src/promoter_atlas/utils/genomics.py",
            "src/promoter_atlas/io/genbank_io.py",
        ],
        "verbatim_lines": {"genomics.py": g, "genbank_io.py": b},
        "reading": (
            "The authors' own code extracts, per annotated gene, the REGULATORY_REGION_LENGTH=200 bp "
            "IMMEDIATELY UPSTREAM of the gene start, strand-aware (the minus strand takes the 200 bp "
            "downstream of the gene end and reverse-complements it), gated by "
            "min_intergenic_length=5 so a gene packed against its neighbour is skipped."
        ),
        "what_this_does_NOT_establish": (
            "That THIS function produced the exact ~9M-sequence h5 the base model was trained on. "
            "scripts/train_base_model.py CONSUMES a prebuilt data/processed/sequence_dataset.h5 via "
            "--data_path and no script in the repo calls extract_promoter_regions to build it. So the "
            "provenance chain is CONVERGENT (same 200 bp geometry as the README's annotation path, "
            "same description as the dataset card) but not proven end to end."
        ),
    })

    # ---- CLAIM 2: corpus size / species count, from the README ---------------------------------
    rd = _grab("README.md", ("9 million regulatory sequences", "~9M regulatory sequences"))
    claims.append({
        "claim": "corpus_size_and_breadth",
        "status": "quoted" if rd and "status" not in rd[0] else "fetch_failed",
        "source_url": "https://github.com/LucasCoppens/PromoterAtlas/blob/main/README.md",
        "verbatim_lines": rd,
    })

    # ---- CLAIM 3: the h5 schema carries NO provenance -------------------------------------------
    ds = _grab(
        "src/promoter_atlas/data/sequence_dataset.py",
        ("self.sequences", "total_sequences", "orig_seq"),
    )
    claims.append({
        "claim": "h5_carries_no_coordinates_or_provenance",
        "status": "quoted" if ds and "status" not in ds[0] else "source_absent",
        "source_paths": ["src/promoter_atlas/data/sequence_dataset.py"],
        "verbatim_lines": ds,
        "reading": (
            "The corpus h5 is ONE dataset 'sequences' of shape (N, seq_length, 4) one-hot plus one "
            "attribute 'total_sequences'. No coordinates, no gene names, no source genome. So "
            "downloading the 2.3 GB h5 could NOT have answered the extraction question that reading "
            "138 lines of the authors' source answered for free."
        ),
    })

    # ---- CLAIM 4: the paper's Methods, which could NOT be read ---------------------------------
    claims.append({
        "claim": "paper_methods_extraction_statement",
        "status": "fetch_failed",
        "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC13377108/",
        "failure_mode": (
            "The PMC route returned a cookie-enablement error page rather than article text. A "
            "summarizer reading that page reported 'NOT STATED', which is the exact collapse this "
            "field exists to prevent: the source was never read, so it is silent about nothing. "
            "Recorded as fetch_failed."
        ),
        "verbatim_quote": "",
        "note": (
            "This claim is NOT load-bearing any more. The extraction rule was recovered from the "
            "authors' source code (claim 1), which is a stronger provenance class than a prose "
            "Methods sentence: it is the executable rule, not a description of one."
        ),
    })

    return {
        "schema": "promoteratlas-corpus-provenance-v1",
        "date": "2026-11-01",
        "question": (
            "Are the PromoterAtlas pretraining corpus sequences upstream-of-gene regions? If so, "
            "'resembles the pretraining corpus' is itself a proxy for 'is a real promoter', which "
            "correlates with measured expression -- the arm's leading alternative explanation for "
            "arms.native150 = 0.5076."
        ),
        "answer": (
            "YES for the authors' own extraction code: 200 bp immediately upstream of each annotated "
            "gene start, strand-aware, min 5 bp intergenic. The arm's honest limit stated this as a "
            "conditional ('IF those are upstream-of-gene regions'); the condition is met by the "
            "repository's own source. The exposure mechanism is therefore LIVE, and the magnitude of "
            "any inflation is a separate empirical question the control-separation experiment asks."
        ),
        "corrections_to_the_plan_that_motivated_this_probe": [
            "The plan recorded that no corpus-construction method is readable in the clone and that "
            "the paper's Methods was 'the one remaining readable source'. That is WRONG. "
            "src/promoter_atlas/utils/genomics.py:6 defines extract_promoter_regions(..., "
            "regulatory_region_length=200) and src/promoter_atlas/io/genbank_io.py:57 defines "
            "extract_promoter_regions_genbank(..., upstream_length=200). The plan's research read "
            "sequence_dataset.py and train_base_model.py and stopped; it never read utils/ or io/.",
            "So the standing lesson fired twice in one direction: reading the source was right, and "
            "the first read was not exhaustive. 'I read the source' is not 'I read the whole source'.",
        ],
        "claims": claims,
    }


def main() -> int:
    rec = build()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(rec, indent=2, sort_keys=True), encoding="utf-8")
    print("answer: %s\n" % rec["answer"])
    for c in rec["claims"]:
        print("  %-42s %s" % (c["claim"], c["status"]))
    print("\n-> %s" % OUT)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
