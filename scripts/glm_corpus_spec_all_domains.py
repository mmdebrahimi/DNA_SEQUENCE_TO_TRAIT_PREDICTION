"""The all-domains corpus SPEC, and the fork that decides how big a download it is.

    uv run python scripts/glm_corpus_spec_all_domains.py      # read-only; writes one artifact

**Why this exists.** The generator's scope is ALL LIFE FORMS. A bacterial framing had crept into my
corpus planning for a bad reason: the locally cached assemblies happen to be bacterial, so the available
substrate quietly became the stated goal. That is scope drift, it was caught, and this file re-does the
arithmetic without it. The project's own north star already said so -- "act in bacteria, then climb to
animals", with AMR as substrate #1 and explicitly not the subject.

**The fork, which is the whole point of this file, and it is not what raw genome sizes suggest.**
Going cross-domain changes the corpus by wildly different factors depending on what the generator is
trained to emit:

  * **Promoter-window mode** -- the geometry the committed falsifier actually scores (150 bp upstream of
    an annotated start). Yield scales with GENE COUNT, and gene counts are remarkably similar across
    domains: measured locally at 1,744-5,501 CDS per bacterial genome, against roughly 20,000 for a
    typical animal. So all-domains buys about **5x per genome**, not more.
  * **Whole-genome mode** -- every window, regulatory or not. Yield scales with GENOME SIZE, which
    differs by three orders of magnitude: measured 1.8-5.7 Mbp per bacterial genome against roughly
    3,100 Mbp for a human. So all-domains buys about **700x per genome**.

The 700x is the tempting number and it is mostly the wrong sequence. A large eukaryotic genome is
dominated by repetitive elements and intergenic distance; coding is a few percent and regulatory regions
a small fraction. Training a promoter generator on it means training mostly on transposons. This is the
caveat already declared in `glm_corpus_distinct_content`: distinct k-mer content measures sequence
NOVELTY, not learnable regulatory signal, and the gap between those two widens enormously once genomes
stop being gene-dense.

**So the decision is not "which genomes" but "which TARGET", and it changes the download by ~70x.**
Both arms are computed below and neither is chosen here, because the target determines what the existing
instrument can score and that is an authority call.

**Provenance rail, following this repo's ECOFF precedent.** Every entry carries a `source`. Bacterial
rows are MEASURED, from artifacts produced locally today, and are the only rows in the headline. Rows
for domains with no local assembly are marked `UNVERIFIED_APPROX` and are confined to a clearly
separated illustrative block -- they are order-of-magnitude figures used to show the SHAPE of the fork,
not inputs to a committed number. A wrong genome size produces a fully self-consistent corpus plan that
is entirely wrong, with nothing downstream able to catch it, which is exactly why the ECOFF catalog
refuses to return a value it cannot source. Verify each against the actual assembly before it is used
for anything but this comparison.

**One closed negative that looks like a blocker and is not.** This repo records a de-confounded
EUKARYOTIC negative -- the Arabidopsis flowering-time arm, where the embedding learned population
structure rather than the trait, and the instruction is not to reopen it. That negative is about
DISCRIMINATIVE zero-shot scoring on natural populations. It is not evidence about GENERATION, which the
project's own orientation notes state as error 0 in bold. A eukaryotic generator is not foreclosed by it,
and citing it as a blocker here would be the exact misreading that file warns about.

No network, no GPU, no download. This writes a plan, not a corpus.
"""

from __future__ import annotations

import json
import statistics
from dataclasses import asdict, dataclass
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
MEASURED_ORG = ROOT / "wiki" / "glm_margin_across_organisms_2026-11-21.json"
MEASURED_DISTINCT = ROOT / "wiki" / "glm_corpus_distinct_content_2026-11-21.json"
OUT = ROOT / "wiki" / "glm_corpus_spec_all_domains_2026-11-21.json"

WINDOW_BP = 150            # the committed falsifier geometry
TOKENS_PER_PARAM = 20      # language-model scaling heuristic, NOT validated for DNA
RUNGS = (1.5e6, 7e6, 20e6, 50e6)

MEASURED = "MEASURED"
UNVERIFIED = "UNVERIFIED_APPROX"


@dataclass(frozen=True)
class DomainEntry:
    """One domain's per-genome yield. `source` decides whether it may enter a headline number."""

    domain: str
    exemplar: str
    genome_mbp: float
    n_genes: int
    source: str
    note: str = ""

    def promoter_mbp(self) -> float:
        """Usable corpus when the target is promoter windows: genes x window length."""
        return self.n_genes * WINDOW_BP / 1e6

    def genome_mbp_usable(self) -> float:
        """Usable corpus when the target is any window: the whole genome."""
        return self.genome_mbp

    def regulatory_fraction(self) -> float:
        """What share of the genome the promoter corpus represents. Collapses in large genomes."""
        return self.promoter_mbp() / self.genome_mbp if self.genome_mbp else 0.0


# --- MEASURED rows. Numbers come from artifacts produced locally; see `load_measured`. -------------
# Populated at runtime so the spec cannot drift from the measurement it claims to rest on.

# --- UNVERIFIED rows. Order-of-magnitude only, to show the SHAPE of the fork. NEVER a headline. ---
ILLUSTRATIVE = (
    DomainEntry("archaea", "a halophile or methanogen", 3.0, 2800, UNVERIFIED,
                "gene-dense like bacteria; expected to behave as a bacterial-like rung"),
    DomainEntry("fungi", "a budding yeast", 12.0, 6000, UNVERIFIED,
                "first rung where introns and chromatin enter; still relatively gene-dense"),
    DomainEntry("plant", "a small model plant", 135.0, 27000, UNVERIFIED,
                "highest gene count of these exemplars, so the best promoter-mode yield per genome"),
    DomainEntry("animal", "a nematode", 100.0, 20000, UNVERIFIED,
                "compact animal genome"),
    DomainEntry("animal", "a mammal", 3100.0, 20000, UNVERIFIED,
                "the extreme case: ~700x a bacterial genome in whole-genome mode and only ~4x in "
                "promoter mode, because gene count barely moves while genome size explodes"),
)


def load_measured() -> list[DomainEntry]:
    """Bacterial rows, read from today's artifacts rather than restated.

    Restating them here is how a spec and its evidence drift apart; the organisms artifact already
    records n_cds_total and the per-genome window count for every substrate it scored.
    """
    if not MEASURED_ORG.exists():
        return []
    rec = json.loads(MEASURED_ORG.read_text(encoding="utf-8"))
    out = []
    for name, v in rec["by_substrate"].items():
        if v.get("status") != "ok":
            continue
        row = next((r for r in rec["rows"]
                    if r.get("substrate") == name and r.get("status") == "ok"), None)
        if row is None:
            continue
        out.append(DomainEntry(
            domain="bacteria",
            exemplar=name,
            genome_mbp=round(row["n_windows_available"] * 0 + v["n_cds_total"] * 0 + 0.0, 3)
            or 0.0,   # replaced below from the distinct-content artifact where available
            n_genes=v["n_cds_total"],
            source=MEASURED,
            note=f"{row['n_windows_available']} usable upstream windows; "
                 f"GC {v['gc_natural_mean']}",
        ))
    return out


def bacterial_genome_sizes() -> dict[str, float]:
    """Mean genome size per bacterial species, from the distinct-content artifact."""
    if not MEASURED_DISTINCT.exists():
        return {}
    rec = json.loads(MEASURED_DISTINCT.read_text(encoding="utf-8"))
    return {org: v["mean_single_genome_distinct_bp_est"] / 1e6
            for org, v in rec["per_species"].items()}


def scaling_table(per_genome_mbp: float, label: str) -> dict:
    """Genomes needed per model rung, at a given per-genome yield."""
    out = {}
    for p in RUNGS:
        need_mbp = p * TOKENS_PER_PARAM / 1e6
        out[f"{p/1e6:g}M_params"] = {
            "corpus_needed_mbp": round(need_mbp, 1),
            "genomes_needed": (round(need_mbp / per_genome_mbp, 1)
                               if per_genome_mbp else None),
        }
    return {"mode": label, "per_genome_mbp": round(per_genome_mbp, 3), "rungs": out}


def main() -> int:
    measured = load_measured()
    sizes = bacterial_genome_sizes()
    # attach real genome sizes to the measured rows
    fixed = []
    for e in measured:
        mbp = 0.0
        for org, v in sizes.items():
            if org.split()[0] in e.exemplar and org.split()[-1] in e.exemplar:
                mbp = v
                break
        if not mbp:
            mbp = statistics.fmean(sizes.values()) if sizes else 0.0
        fixed.append(DomainEntry(e.domain, e.exemplar, round(mbp, 3), e.n_genes,
                                 e.source, e.note))
    measured = fixed

    if not measured:
        print("REFUSE: no MEASURED bacterial rows available; the headline would rest on nothing.")
        return 2

    print(f"window={WINDOW_BP}bp  tokens/param={TOKENS_PER_PARAM} (heuristic, unvalidated for DNA)")
    print()
    print("MEASURED rows (bacteria, from today's artifacts):")
    for e in measured:
        print(f"  {e.exemplar[:44]:44s} {e.n_genes:5d} genes  {e.genome_mbp:6.2f} Mbp  "
              f"promoter {e.promoter_mbp():.2f} Mbp  ({e.regulatory_fraction()*100:5.1f}% of genome)")
    bact_prom = statistics.fmean(e.promoter_mbp() for e in measured)
    bact_gen = statistics.fmean(e.genome_mbp for e in measured)

    print()
    print("ILLUSTRATIVE rows (UNVERIFIED order-of-magnitude; never a headline):")
    for e in ILLUSTRATIVE:
        print(f"  {e.domain:9s} {e.exemplar[:30]:30s} {e.n_genes:6d} genes  {e.genome_mbp:7.1f} Mbp  "
              f"promoter {e.promoter_mbp():5.2f} Mbp  ({e.regulatory_fraction()*100:5.2f}% of genome)")

    euk = [e for e in ILLUSTRATIVE if e.domain in ("plant", "animal", "fungi")]
    euk_prom = statistics.fmean(e.promoter_mbp() for e in euk)
    euk_gen = statistics.fmean(e.genome_mbp for e in euk)

    rec = {
        "record": "glm-corpus-spec-all-domains-v1",
        "date": date.today().isoformat(),
        "status": "PLAN, not a measurement. No corpus is downloaded or built by this script.",
        "scope_correction": (
            "The generator's scope is ALL LIFE FORMS. A bacterial framing had entered my corpus "
            "planning because the locally cached assemblies are bacterial -- the available substrate "
            "quietly became the stated goal. Corrected here. The project's north star already said "
            "bacteria is a stage, not the subject."),
        "the_fork": {
            "question": "Does the generator emit PROMOTER WINDOWS or ANY GENOMIC WINDOW?",
            "why_it_dominates": (
                "Promoter-mode corpus scales with GENE COUNT, which barely moves across domains "
                "(1,744-5,501 CDS measured per bacterial genome against roughly 20,000 for an animal). "
                "Whole-genome corpus scales with GENOME SIZE, which spans three orders of magnitude "
                "(1.8-5.7 Mbp measured against roughly 3,100 Mbp for a mammal). So going cross-domain "
                "buys about 5x per genome in one mode and about 700x in the other."),
            "why_the_700x_is_mostly_the_wrong_sequence": (
                "A large eukaryotic genome is dominated by repetitive elements and intergenic distance; "
                "coding is a few percent and regulatory regions a small fraction. A promoter generator "
                "trained on it would train mostly on transposons. This is the already-declared caveat "
                "that distinct k-mer content measures sequence NOVELTY and not learnable regulatory "
                "signal -- and the gap between those two widens sharply once genomes stop being "
                "gene-dense."),
            "consequence": (
                "The download size differs by roughly 70x between the two modes, and the mode also "
                "decides whether the committed falsifier can score the output at all, since it is "
                "built on 150 bp upstream windows."),
            "not_chosen_here": "the target is an authority call, so both arms are computed and neither "
                               "is selected",
        },
        "measured_rows": [asdict(e) | {"promoter_mbp": round(e.promoter_mbp(), 4),
                                       "regulatory_fraction": round(e.regulatory_fraction(), 5)}
                          for e in measured],
        "illustrative_rows_UNVERIFIED": [
            asdict(e) | {"promoter_mbp": round(e.promoter_mbp(), 4),
                         "regulatory_fraction": round(e.regulatory_fraction(), 5)}
            for e in ILLUSTRATIVE],
        "scaling_from_MEASURED_bacteria": {
            "promoter_mode": scaling_table(bact_prom, "promoter windows, bacteria (MEASURED)"),
            "whole_genome_mode": scaling_table(bact_gen, "whole genome, bacteria (MEASURED)"),
        },
        "scaling_from_ILLUSTRATIVE_eukaryotes": {
            "caveat": "every input here is UNVERIFIED_APPROX; this block shows the SHAPE of the fork "
                      "and must not be quoted as a corpus plan",
            "promoter_mode": scaling_table(euk_prom, "promoter windows, eukaryotes (UNVERIFIED)"),
            "whole_genome_mode": scaling_table(euk_gen, "whole genome, eukaryotes (UNVERIFIED)"),
        },
        "corpus_spec_requirements": {
            "domain_label_per_genome": (
                "REQUIRED, as a first-class field: bacteria | archaea | fungi | plant | animal. "
                "Cross-domain regulatory grammar differs structurally -- bacterial -10/-35 boxes and "
                "operons against eukaryotic TATA, enhancers, chromatin and splicing -- so an "
                "undifferentiated pool asks one model to learn mutually inconsistent grammars with no "
                "signal telling it which applies. Domain conditioning or a staged curriculum needs the "
                "label to exist from the start, and retrofitting it after a corpus is built is the "
                "expensive way to get it."),
            "lineage_dedup": (
                "REQUIRED, and measured as necessary rather than assumed: the local cache holds 400 "
                "Campylobacter jejuni assemblies whose 12th genome adds 0.50% new distinct content. "
                "Deduplicate by LINEAGE, not by accession. Un-deduplicated, a corpus over-represents "
                "whichever lineage happens to be best sampled, and the model memorises it -- the "
                "confound that has already closed four arms in this project."),
            "diversity_over_depth": (
                "At the measured bacterial yield, distinct content comes from adding SPECIES, not "
                "strains. One representative per lineage, many lineages."),
            "scoring_geometry": (
                "Whatever the corpus, the generator must emit something the committed falsifier can "
                "score, or a new instrument is needed before any number means anything. Today that is "
                f"{WINDOW_BP} bp upstream windows."),
            "positional_architecture_is_the_target": (
                "Measured today: a k-mer discriminator is saturated BY CONSTRUCTION because a Markov "
                "chain's sufficient statistic IS the k-mer conditional distribution, while the "
                "positional margin holds at +0.1740 and is invariant to the null's order, replicated "
                "at GC 0.266 and 0.453. So a generator that matches only local composition clears the "
                "k-mer view and is caught flat by the positional one. Train for WHERE motifs sit."),
        },
        "a_closed_negative_that_is_NOT_a_blocker_here": {
            "what": "the Arabidopsis flowering-time arm, a de-confounded EUKARYOTIC negative where the "
                    "embedding learned population structure rather than the trait",
            "why_it_does_not_block_this": (
                "it is about DISCRIMINATIVE zero-shot scoring on natural populations, not about "
                "GENERATION. The project's own orientation notes flag conflating those as error 0. A "
                "eukaryotic generator is not foreclosed by it, and citing it as a blocker would be the "
                "exact misreading that file warns against."),
        },
        "local_substrate_status": {
            "bacterial_genomes_with_annotation": 459,
            "measured_distinct_bacterial_content_mbp": 30.9,
            "eukaryotic_genomes_available_locally": 0,
            "note": ("data/arabidopsis/ holds phenotype tables and the 1001-genomes flowering data, "
                     "not an assembly. So any all-domains corpus requires downloads for every "
                     "non-bacterial domain, in both modes."),
        },
        "honest_limits": [
            "This is a PLAN. Nothing here is a measurement of a eukaryotic corpus, because no "
            "eukaryotic assembly is available locally to measure.",
            "Every eukaryotic figure is UNVERIFIED_APPROX and order-of-magnitude. They are confined to "
            "an illustrative block and excluded from every headline, following the ECOFF precedent in "
            "this repo: a wrong constant produces a fully self-consistent plan that is entirely wrong, "
            "with nothing downstream able to catch it. Verify each against the real assembly before "
            "use.",
            "Gene COUNT is used as a proxy for promoter-corpus yield. That ignores alternative "
            "promoters, bidirectional promoters and unannotated transcription start sites, so it is a "
            "lower bound in eukaryotes and close to exact in gene-dense bacteria.",
            "The 20-tokens-per-parameter ratio is a language-model scaling heuristic, imported and NOT "
            "validated for DNA. A different ratio moves every genome count proportionally.",
            "Genome size is taken as distinct content in whole-genome mode, which OVERSTATES it: "
            "repetitive elements are a large fraction of big genomes and compress heavily under "
            "distinct-k-mer counting. So the 700x figure is itself an upper bound on the real gain.",
            "No claim is made that cross-domain training helps. Mixing grammars could easily hurt, "
            "which is precisely why the domain label and a staged curriculum are specified as "
            "requirements rather than options.",
        ],
    }
    OUT.write_text(json.dumps(rec, indent=2), encoding="utf-8")

    print()
    print("THE FORK, genomes needed per rung:")
    for block, src in (("scaling_from_MEASURED_bacteria", "bacteria MEASURED"),
                       ("scaling_from_ILLUSTRATIVE_eukaryotes", "eukaryotes UNVERIFIED")):
        for mode in ("promoter_mode", "whole_genome_mode"):
            t = rec[block][mode]
            cells = "  ".join(
                f"{k.replace('_params','')}={v['genomes_needed']:g}"
                for k, v in t["rungs"].items())
            print(f"  [{src:22s}] {mode:17s} {t['per_genome_mbp']:7.2f} Mbp/genome -> {cells}")
    print(f"\n-> {OUT.relative_to(ROOT).as_posix()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
