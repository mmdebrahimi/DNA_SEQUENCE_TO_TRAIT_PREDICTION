"""The multi-genome TRAINING corpus for a generator — domain-labelled, deduplicated, and leak-gated.

**This is a driver over `corpus.build_natural_set`, never a second definition of a window.** The single
most expensive mistake available here is to grow a parallel description of what a promoter window is:
the falsifier scores one definition, a trainer would learn another, and the two would differ by an
off-by-one nobody could see. So every window in this module comes out of `build_natural_set`, and the
only things added are which genomes to draw from, how to deduplicate them, and what must be withheld.

Three properties this enforces, each because something measured says it is necessary rather than tidy.

**1. LEAK GATE, and this is the one that would silently invent success.** The committed falsifier scores
a generator against the eval half of `build_natural_set(<reference>, limit=1200, seed=s)`. If a training
corpus contains those same windows, the generator is trained on the exact sequences it is later scored
against -- and the failure is invisible, because a leaked generator produces beautiful numbers and no
error. `build_training_corpus` therefore EXCLUDES the falsifier's scoring windows by default, across
every seed the scoring will use, and refuses to include them unless a caller passes
`allow_falsifier_overlap=True` and thereby says so in writing. The count withheld is recorded in the
manifest, so a run that quietly stopped withholding is visible as a zero.

**2. SEQUENCE deduplication across genomes, measured as necessary.** Distinct-content measurement on the
cached assemblies found that twelve *E. coli* genomes carry only 2.72x the distinct content of one, so
roughly 63% of the raw window total is duplicate. Without window-level deduplication a corpus looks
multiples larger than it is and the model sees the same promoter many times, which is memorisation
dressed as scale. Deduplication is by SEQUENCE, not by coordinate: the same upstream region appears at
different coordinates in different assemblies, and a model only ever sees the sequence.

**3. DOMAIN label as a required field, not an annotation.** Bacterial -10/-35 boxes and operons against
eukaryotic TATA, enhancers, chromatin and splicing are structurally different grammars. An
undifferentiated pool asks one model to learn mutually inconsistent rules with no signal saying which
applies, so domain conditioning or a staged curriculum needs the label to exist from the first byte.
`CorpusEntry` refuses to construct without one, because retrofitting it after a corpus is built is the
expensive way to get it.

Honest scope: this builds a corpus and reports what it contains. It makes no claim that cross-domain
training helps -- mixing grammars could easily hurt, which is exactly why the label is mandatory and the
per-domain counts are reported rather than averaged away.
"""

from __future__ import annotations

from dataclasses import dataclass, field

VALID_DOMAINS = frozenset({"bacteria", "archaea", "fungi", "plant", "animal"})


class CorpusSpecError(ValueError):
    """A corpus specification that cannot be built as written."""


@dataclass(frozen=True)
class CorpusEntry:
    """One genome in the training corpus.

    `lineage_key` is recorded rather than derived: whoever selected the species deduplicated on some
    notion of lineage, and the corpus should carry THAT reason rather than a second guess at it. A
    manifest that records the real selection key can be audited against the selection; one that
    re-derives its own cannot.
    """

    accession: str
    species: str
    domain: str
    lineage_key: str
    reference_dir: str
    source: str = "ncbi_datasets_v2"
    note: str = ""

    def __post_init__(self) -> None:
        if self.domain not in VALID_DOMAINS:
            raise CorpusSpecError(
                f"domain {self.domain!r} for {self.accession} is not one of "
                f"{sorted(VALID_DOMAINS)}. The label is required, not optional: an unlabelled genome "
                "in a cross-domain pool cannot be conditioned on or curriculum-staged later.")
        if not self.lineage_key:
            raise CorpusSpecError(
                f"{self.accession} has no lineage_key. Deduplication is by LINEAGE, not accession -- "
                "the cached set holds 400 Campylobacter jejuni assemblies whose 12th genome adds 0.50% "
                "new content, so an accession-keyed corpus over-represents whichever lineage happens to "
                "be best sampled.")


@dataclass
class CorpusManifest:
    """What the corpus contains and what was withheld. Every count is a check someone can run again."""

    entries: list[dict] = field(default_factory=list)
    n_windows_kept: int = 0
    n_windows_duplicate_dropped: int = 0
    n_windows_withheld_as_falsifier_eval: int = 0
    by_domain: dict[str, int] = field(default_factory=dict)
    by_species: dict[str, int] = field(default_factory=dict)
    lineage_keys: list[str] = field(default_factory=list)
    window_length_bp: int = 0
    falsifier_leak_gate: dict = field(default_factory=dict)

    def as_dict(self) -> dict:
        total_raw = (self.n_windows_kept + self.n_windows_duplicate_dropped
                     + self.n_windows_withheld_as_falsifier_eval)
        return {
            "n_genomes": len(self.entries),
            "window_length_bp": self.window_length_bp,
            "n_windows_kept": self.n_windows_kept,
            "n_windows_duplicate_dropped": self.n_windows_duplicate_dropped,
            "n_windows_withheld_as_falsifier_eval": self.n_windows_withheld_as_falsifier_eval,
            "n_windows_raw_total": total_raw,
            "duplicate_fraction": (round(self.n_windows_duplicate_dropped / total_raw, 4)
                                   if total_raw else 0.0),
            "corpus_bp": self.n_windows_kept * self.window_length_bp,
            "by_domain": dict(sorted(self.by_domain.items())),
            "by_species": dict(sorted(self.by_species.items())),
            "n_lineages": len(set(self.lineage_keys)),
            "falsifier_leak_gate": self.falsifier_leak_gate,
            "entries": self.entries,
        }


def assert_one_genome_per_lineage(entries: list[CorpusEntry]) -> None:
    """Refuse a spec that puts two genomes of one lineage in the corpus.

    Checked rather than silently collapsed, because two entries sharing a lineage key usually means the
    selection rule was not applied, and quietly dropping one would hide that from whoever wrote the spec.
    """
    seen: dict[str, str] = {}
    clashes = []
    for e in entries:
        if e.lineage_key in seen:
            clashes.append((e.lineage_key, seen[e.lineage_key], e.accession))
        else:
            seen[e.lineage_key] = e.accession
    if clashes:
        detail = "; ".join(f"{k}: {a} and {b}" for k, a, b in clashes)
        raise CorpusSpecError(
            f"{len(clashes)} lineage key(s) appear more than once ({detail}). Deduplicate by lineage "
            "before building, or pass distinct keys if these really are different lineages.")


def falsifier_eval_windows(reference_dir, *, length: int, limit: int, seeds, kind: str) -> set[str]:
    """Every window the committed falsifier could score a generator against, across all seeds.

    Built through `build_natural_set` and `fit_eval_halves` -- the same functions the scorer uses -- so
    the withheld set is defined by the scorer rather than by a second description of it. The union over
    seeds matters: the scorer reshuffles per seed, so a window in seed 3's eval half may sit in seed 7's
    fit half, and withholding only one seed's eval set would still leak into the others.
    """
    from dna_decode.glm.corpus import build_natural_set, fit_eval_halves

    out: set[str] = set()
    for s in seeds:
        wins, _ = build_natural_set(reference_dir, length=length, limit=limit, seed=s, kind=kind)
        _fit, ev = fit_eval_halves(wins)
        out.update(ev)
    return out


def build_training_corpus(
    entries: list[CorpusEntry],
    *,
    length: int = 150,
    kind: str = "promoter",
    limit_per_genome: int | None = None,
    withhold: set[str] | None = None,
    allow_falsifier_overlap: bool = False,
) -> tuple[list[str], CorpusManifest]:
    """Build a deduplicated, domain-labelled window corpus from several genomes.

    Returns `(windows, manifest)`. Windows are plain sequences in a deterministic order: genomes in the
    order given, windows in `build_natural_set`'s order for each, first occurrence of a sequence kept.
    Determinism matters because a shuffled corpus would make two runs incomparable for no benefit --
    the trainer shuffles, and it should own that.

    `withhold` is the leak gate. Passing None with `allow_falsifier_overlap=False` is a programming
    error rather than a default, because silently training on scoring windows is the failure this module
    exists to prevent and it produces no error of its own.
    """
    from dna_decode.glm.corpus import build_natural_set

    if withhold is None and not allow_falsifier_overlap:
        raise CorpusSpecError(
            "no withhold set given and allow_falsifier_overlap is False. Pass the falsifier's eval "
            "windows to withhold them, or pass allow_falsifier_overlap=True to state in writing that "
            "this corpus may contain the sequences the generator will be scored against. There is no "
            "safe default here: a leaked generator produces excellent numbers and no error.")
    withhold = withhold or set()
    if entries:
        assert_one_genome_per_lineage(entries)

    man = CorpusManifest(window_length_bp=length)
    man.falsifier_leak_gate = {
        "enforced": not allow_falsifier_overlap,
        "n_sequences_withheld": len(withhold),
        "meaning": ("windows the committed falsifier scores against are excluded from training, so a "
                    "generator cannot be trained on the sequences it is graded on"),
    }
    seen: set[str] = set()
    out: list[str] = []
    for e in entries:
        wins, stats = build_natural_set(e.reference_dir, length=length,
                                        limit=limit_per_genome, seed=0, kind=kind)
        kept = dup = held = 0
        for w in wins:
            seq = w.seq
            if seq in withhold:
                held += 1
                continue
            if seq in seen:
                dup += 1
                continue
            seen.add(seq)
            out.append(seq)
            kept += 1
        man.n_windows_kept += kept
        man.n_windows_duplicate_dropped += dup
        man.n_windows_withheld_as_falsifier_eval += held
        man.by_domain[e.domain] = man.by_domain.get(e.domain, 0) + kept
        man.by_species[e.species] = man.by_species.get(e.species, 0) + kept
        man.lineage_keys.append(e.lineage_key)
        man.entries.append({
            "accession": e.accession, "species": e.species, "domain": e.domain,
            "lineage_key": e.lineage_key, "source": e.source, "note": e.note,
            "n_windows_available": len(wins),
            "n_windows_kept": kept,
            "n_windows_duplicate_dropped": dup,
            "n_windows_withheld": held,
            "n_cds_total": getattr(stats, "n_cds_total", None),
        })
    return out, man


def domain_balance_report(man: CorpusManifest) -> dict:
    """How skewed the corpus is by domain. Reported, never silently rebalanced.

    A corpus that is 95% one domain will train a model that is effectively single-domain whatever its
    label says, and that is a fact the trainer's caller needs rather than something to paper over by
    resampling behind their back.
    """
    total = sum(man.by_domain.values())
    if not total:
        return {"status": "empty"}
    shares = {d: round(n / total, 4) for d, n in sorted(man.by_domain.items())}
    top = max(shares, key=lambda d: shares[d])
    return {
        "n_domains": len(shares),
        "shares": shares,
        "largest_domain": top,
        "largest_share": shares[top],
        "single_domain_in_practice": bool(shares[top] >= 0.95),
        "reading": ("a corpus dominated by one domain trains a single-domain model whatever its labels "
                    "say; this is reported so a caller can stage a curriculum or weight the sampler, "
                    "and is deliberately NOT corrected here"),
    }
