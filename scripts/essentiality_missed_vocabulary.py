"""Why does the conserved-core decoder MISS 83% of human essential genes? A fabrication-free answer.

THE QUESTION THIS SETTLES, AND WHY IT IS CHEAP.
The decoder (`dna_decode/essentiality/core_decoder.py`) scores E. coli essentials at AUROC 0.695 and
human at 0.580. A cross-organism "transfer ladder" was proposed to decide whether that decay is
(A) continuous with phylogenetic distance or (B) a cliff at the eukaryote boundary. But the decay is
not a discrimination failure at all -- it is SILENCE: **83.1% of human essential genes score EXACTLY
zero** (E. coli 62.1%), and in human 11 of 13 core patterns have P(hit | non-essential) = 0.0000, so
the patterns never misfire, they just never fire. AUROC blends coverage and precision and hides both.

So the real question is about the MISSED set, and a third mechanism is possible that neither (A) nor
(B) covers:

  (i)   bacteria-only function      -- correctly absent in human (peptidoglycan, BAM complex)
  (ii)  conserved function, BACTERIA-SPECIFIC NAME -- a catalogue PHRASING gap, not biology
  (iii) host-specific function       -- genuinely outside a bacterial catalogue (proteasome, spliceosome)

(ii) matters out of proportion: if it carries most of the gap, the honest conclusion is "the catalogue
has a phrasing bug", nothing about phylogeny, and a 5-rung ladder would be the WRONG INSTRUMENT -- it
would read a vocabulary artifact as a biological law. Two patterns already hint at it: `DNA polymerase
III|...|replicative` and `cell division protein Fts|divisome|...` both drop to EXACTLY 0.00 hits per
1000 human genes, yet human plainly has replicative polymerases and cell division.

HOW THIS AVOIDS ASSERTING BIOLOGY (the load-bearing design constraint).
The tempting implementation is a hand-written synonym map ("POLD1 means replication"). That is the
fabrication hazard this repo already refuses elsewhere -- `dna_decode/data/ecoff_catalog.py` requires a
`source_url` + `verbatim_quote` and RAISES rather than return a plausible default, because a wrong
curated value yields a fully self-consistent evaluation that is entirely wrong. A synonym map would be
exactly that: an unsourced biological claim doing the adjudicating.

So NO biology is written here. Every word comes from NCBI's own `description` / `name` fields, and the
adjudication is a CORPUS COMPARISON the data performs on itself:

  for each content word W frequent among MISSED human essentials, ask where W lives in E. coli:
    W appears in E. coli descriptions the catalogue CATCHES      -> IN_CATALOGUE_REACH   (mechanism ii)
    W appears in E. coli but only in descriptions it MISSES      -> ECOLI_BUT_UNCAUGHT   (ii-or-i)
    W does not appear in E. coli CDS descriptions at all         -> ABSENT_FROM_ECOLI    (mechanism iii)

`IN_CATALOGUE_REACH` is the decisive bucket: the function is present in E. coli AND the catalogue sees
it THERE, so failing to see it in human is a vocabulary gap, not a missing function. Each bucket is a
statement about WORDS IN TWO ANNOTATION CORPORA, which is checkable, not about biology, which would be
asserted.

HONEST LIMITS, stated because they bound the verdict:
  * A word is not a function. "polymerase" appearing in both corpora does not prove the same complex;
    this measures VOCABULARY OVERLAP, which is precisely the thing (ii) is about, and nothing more.
  * The human label set is BAGEL CEGv2/NEGv1 -- two curated extremes, not a genome-wide screen -- so
    its base rate (0.431) is a construction artifact and is NOT comparable to E. coli's 0.0928.
  * Stopwords are a hand-written list, but they are ENGLISH function words, not biology.
  * This cannot distinguish (i) from (ii) inside `ECOLI_BUT_UNCAUGHT`; it is reported as the honest
    joint bucket rather than split by assertion.

Run: uv run python scripts/essentiality_missed_vocabulary.py [--top 25] [--emit]
Read-only unless --emit. Needs the D: essentiality cache; exits 2 if absent.
"""
from __future__ import annotations

import argparse
import gzip
import json
import re
import sys
from collections import Counter
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

TGT = Path("D:/dna_decode_cache/essentiality")
WIKI = ROOT / "wiki"

# ENGLISH function words + annotation boilerplate. Deliberately contains no biology: removing
# "protein" or "domain" changes which words surface, so they are listed as boilerplate, not as a
# judgement about what is or is not a function.
# PROMOTED to the package (2026-10-06, transfer-ladder Step 4) so the ladder does not import a scripts/
# module and there is exactly ONE definition of this vocabulary. Re-exported here because this script's
# own tests and callers reference the module-level names.
from dna_decode.essentiality.phrasing_floor import (  # noqa: E402
    ROBUST_II,
    STOP,
    WORD,
    content_words,
)


def _need(*names):
    missing = [n for n in names if not (TGT / n).exists()]
    if missing:
        print("MISSING from %s: %s" % (TGT, ", ".join(missing)))
        print("This probe needs the D: essentiality cache (BAGEL + gene_info + E. coli feature table).")
        return False
    return True


def human_rows():
    """(symbol, description) for BAGEL CEGv2 essentials and NEGv1 non-essentials."""
    desc = {}
    with gzip.open(TGT / "Homo_sapiens.gene_info.gz", "rt", encoding="utf-8", errors="replace") as f:
        h = f.readline().rstrip("\n").split("\t")
        gid, sym, dsc, ty = h.index("GeneID"), h.index("Symbol"), h.index("description"), h.index("type_of_gene")
        for line in f:
            p = line.rstrip("\n").split("\t")
            if p[ty] == "protein-coding":
                desc[p[gid]] = (p[sym], p[dsc])

    def load(fn):
        out = []
        for c in (TGT / fn).read_text(encoding="utf-8", errors="replace").splitlines()[1:]:
            parts = c.split("\t")
            if len(parts) >= 3 and parts[2] in desc:
                out.append(desc[parts[2]])
        return out

    return load("CEGv2.txt"), load("NEGv1.txt")


def ecoli_rows():
    """(symbol, product) for every E. coli K-12 CDS in the NCBI feature table."""
    out = []
    with gzip.open(TGT / "ecoli_k12_feature_table.txt.gz", "rt", encoding="utf-8", errors="replace") as f:
        h = f.readline().rstrip("\n").split("\t")
        i, n = h.index("symbol"), h.index("name")
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) > max(i, n) and p[0] == "CDS" and p[i]:
                out.append((p[i], p[n]))
    return out


def classify_word(word, ecoli_caught_words, ecoli_all_words):
    """Where does this word live in the E. coli corpus? A corpus fact, not a biological claim."""
    if word in ecoli_caught_words:
        return "IN_CATALOGUE_REACH"      # mechanism (ii): catalogue sees this word's genes in E. coli
    if word in ecoli_all_words:
        return "ECOLI_BUT_UNCAUGHT"      # (ii)-or-(i): present in E. coli, catalogue blind there too
    return "ABSENT_FROM_ECOLI"           # mechanism (iii): host-specific vocabulary


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--top", type=int, default=25)
    ap.add_argument("--emit", action="store_true", help="write the wiki artifact")
    a = ap.parse_args(argv)

    if not _need("CEGv2.txt", "NEGv1.txt", "Homo_sapiens.gene_info.gz", "ecoli_k12_feature_table.txt.gz"):
        return 2

    from dna_decode.essentiality.core_decoder import _CORE, score_gene

    ceg, neg = human_rows()
    eco = ecoli_rows()

    caught = [(g, d) for g, d in ceg if score_gene(g, d).core_score != 0]
    missed = [(g, d) for g, d in ceg if score_gene(g, d).core_score == 0]
    eco_caught = [(g, p) for g, p in eco if score_gene(g, p).core_score != 0]

    eco_caught_words = set()
    for g, p in eco_caught:
        eco_caught_words |= content_words("%s %s" % (g, p))
    eco_all_words = set()
    for g, p in eco:
        eco_all_words |= content_words("%s %s" % (g, p))

    freq = Counter()
    for g, d in missed:
        freq.update(content_words(d))

    rows = []
    for w, n in freq.most_common(a.top):
        rows.append({"word": w, "n_missed_essentials": n,
                     "share_of_missed": round(n / len(missed), 4),
                     "bucket": classify_word(w, eco_caught_words, eco_all_words)}
                    )

    # Gene-level attribution: a missed gene counts toward (ii) if ANY of its content words is already
    # inside the catalogue's E. coli reach. Deliberately GENEROUS to mechanism (ii) -- if (ii) still
    # fails to dominate under a bar tilted in its favour, that is the stronger finding.
    per_gene = Counter()
    for g, d in missed:
        ws = content_words(d)
        if ws & eco_caught_words:
            per_gene["ii_IN_CATALOGUE_REACH"] += 1
        elif ws & eco_all_words:
            per_gene["i_or_ii_ECOLI_BUT_UNCAUGHT"] += 1
        else:
            per_gene["iii_ABSENT_FROM_ECOLI"] += 1

    n_missed = len(missed)

    # ------------------------------------------------------------------------------------------
    # THREE BARS, because ONE BAR PUBLISHED A FALSE VERDICT and the next two disagreed with it.
    #
    # The first version of this probe used only `loose` and printed
    # PHRASING_GAP_DOMINATES_LADDER_IS_WRONG_INSTRUMENT at 0.5548. That was inflated: `loose` counts a
    # gene whenever ANY content word is shared with a caught E. coli gene, and the words doing the work
    # were `binding` / `small` / `beta` / `alpha` / `cell` / `repeat` -- generic modifiers, not
    # functions. Sharing "alpha" with E. coli says nothing about whether the catalogue could see a gene.
    #
    # `strict` is the honest definition of the catalogue's reach: words its REGEXES key on, taken from
    # the pattern sources themselves. It gives 0.2332 -- a DIFFERENT verdict band from `loose`.
    # `strict_nofrag` additionally drops `cell` and `iii`, which are tokenization debris from splitting
    # `cell division protein Fts` and `DNA polymerase III`. It gives ~0.1996 -- a THIRD band.
    #
    # Three defensible tokenizations spanning two threshold crossings means the quantity is not
    # measurable this way, so the probe REFUSES a verdict instead of picking the flattering bar.
    # ------------------------------------------------------------------------------------------
    pat_words = set()
    for _w, pat in _CORE:
        for tok in re.findall(r"[A-Za-z][A-Za-z-]{2,}", pat.pattern):
            pat_words.add(tok.lower())
    pat_words -= STOP
    FRAGMENTS = {"cell", "iii"}   # debris from splitting multi-word regex phrases

    def frac(reach):
        return sum(1 for _g, d in missed if content_words(d) & reach) / n_missed if n_missed else 0.0

    bars = {
        "loose_any_shared_word": round(frac(eco_caught_words), 4),
        "strict_regex_vocabulary": round(frac(pat_words), 4),
        "strict_minus_regex_fragments": round(frac(pat_words - FRAGMENTS), 4),
    }

    def band(f):
        return ("dominates" if f >= 0.50 else "material" if f >= 0.20 else "host_specific")

    bands = {k: band(v) for k, v in bars.items()}
    frac_ii = bars["strict_regex_vocabulary"]        # the defensible point estimate, not the headline

    if len(set(bands.values())) > 1:
        verdict = "INDETERMINATE_BAR_SENSITIVE_VOCABULARY_CANNOT_ADJUDICATE"
    else:
        verdict = {"dominates": "PHRASING_GAP_DOMINATES_LADDER_IS_WRONG_INSTRUMENT",
                   "material": "PHRASING_GAP_MATERIAL_BUT_NOT_DOMINANT",
                   "host_specific": "GAP_IS_GENUINELY_HOST_SPECIFIC_LADDER_MEASURES_BIOLOGY"}[
            next(iter(set(bands.values())))]

    # What survives EVERY bar, because it does not depend on the reach definition at all:
    #   * a word absent from the E. coli corpus is absent under any definition of reach
    #   * a word that is BOTH in the regex vocabulary AND unambiguously functional is in reach always
    # ROBUST_II now lives in dna_decode/essentiality/phrasing_floor.py, where an IMPORT-TIME gate
    # asserts every word already appears in the decoder's own _CORE patterns.
    robust_ii = sum(1 for _g, d in missed if content_words(d) & ROBUST_II)
    robust_iii = per_gene["iii_ABSENT_FROM_ECOLI"]

    print("human BAGEL CEGv2 essentials: %d  (caught %d / missed %d  -> coverage %.4f)"
          % (len(ceg), len(caught), n_missed, len(caught) / len(ceg)))
    print("E. coli CDS: %d  (catalogue catches %d -> %d distinct content words in reach)"
          % (len(eco), len(eco_caught), len(eco_caught_words)))
    print("\nGENE-LEVEL attribution of the %d MISSED human essentials, under the LOOSE bar" % n_missed)
    print("(loose = any single shared word; it is the INFLATED bar -- see BAR SENSITIVITY below)")
    for k in ("ii_IN_CATALOGUE_REACH", "i_or_ii_ECOLI_BUT_UNCAUGHT", "iii_ABSENT_FROM_ECOLI"):
        v = per_gene[k]
        print("  %-28s %5d  %.4f" % (k, v, v / n_missed if n_missed else 0))
    print("\nBAR SENSITIVITY (the reason there is no single frac_ii):")
    for k, v in bars.items():
        print("  %-32s %.4f  -> %s" % (k, v, bands[k]))
    print("\nVERDICT: %s" % verdict)
    print("ROBUST under every bar:")
    print("  mechanism (ii) floor  (polymerase/helicase/replication/division/topoisomerase/primase): "
          "%d  %.4f" % (robust_ii, robust_ii / n_missed))
    print("  mechanism (iii) floor (no E. coli vocabulary at all):                                  "
          "%d  %.4f" % (robust_iii, robust_iii / n_missed))

    print("\nTop %d content words among MISSED human essentials:" % a.top)
    print("  %-22s %6s %8s  bucket" % ("word", "n", "share"))
    for r in rows:
        print("  %-22s %6d %8.4f  %s" % (r["word"], r["n_missed_essentials"],
                                         r["share_of_missed"], r["bucket"]))

    art = {
        "record": "essentiality-missed-vocabulary-v1",
        "date": str(date.today()),
        "question": "is the conserved-core decoder's human miss a catalogue PHRASING gap (mechanism ii) "
                    "or genuinely host-specific function (mechanism iii)?",
        "method": "corpus comparison only -- no curated biology. Content words of MISSED human BAGEL "
                  "essentials are located in the E. coli CDS corpus: inside the catalogue's reach, "
                  "present-but-uncaught, or absent entirely.",
        "n_human_essential": len(ceg), "n_caught": len(caught), "n_missed": n_missed,
        "human_coverage": round(len(caught) / len(ceg), 4),
        "n_ecoli_cds": len(eco), "n_ecoli_caught": len(eco_caught),
        "gene_level_attribution": dict(per_gene),
        "frac_mechanism_ii_strict": round(frac_ii, 4),
        "bar_sensitivity": bars,
        "bar_sensitivity_bands": bands,
        "verdict": verdict,
        "verdict_rule": ">=0.50 phrasing dominates; 0.20-0.50 material; <0.20 host-specific; "
                        "bands disagreeing across bars -> INDETERMINATE",
        "robust_floors": {
            "mechanism_ii_floor_n": robust_ii,
            "mechanism_ii_floor_frac": round(robust_ii / n_missed, 4) if n_missed else 0.0,
            "mechanism_iii_floor_n": robust_iii,
            "mechanism_iii_floor_frac": round(robust_iii / n_missed, 4) if n_missed else 0.0,
            "why_robust": "an absent word is absent under ANY reach definition; the (ii) floor words are "
                          "both in the regex vocabulary and unambiguously functional, so they are in "
                          "reach under any definition. The MIDDLE is what the bars disagree about.",
        },
        "self_correction": "the first version of this probe used only the loose bar, reported 0.5548 and "
                           "a verdict of PHRASING_GAP_DOMINATES_LADDER_IS_WRONG_INSTRUMENT. That was "
                           "inflated >2x by generic modifiers (binding/small/beta/alpha/cell/repeat) and "
                           "the strict bar gives a different band. The verdict was retracted before "
                           "publication, not after.",
        "top_words": rows,
        "honest_limits": [
            "a WORD is not a function: shared vocabulary does not prove the same complex. This measures "
            "vocabulary overlap, which is what mechanism (ii) IS, and nothing more.",
            "human labels are BAGEL CEGv2/NEGv1, two curated extremes, so base rate 0.431 is a "
            "construction artifact and not comparable to E. coli's genome-wide 0.0928.",
            "ECOLI_BUT_UNCAUGHT cannot be split into (i) bacteria-only vs (ii) phrasing without "
            "asserting biology, so it is reported as the joint bucket.",
            "stopwords are hand-written but are English function words plus annotation boilerplate, "
            "not biological judgements.",
        ],
    }
    if a.emit:
        p = WIKI / ("essentiality_missed_vocabulary_%s.json" % date.today())
        p.write_text(json.dumps(art, indent=2), encoding="utf-8")
        print("\n[-> %s]" % p.relative_to(ROOT))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
