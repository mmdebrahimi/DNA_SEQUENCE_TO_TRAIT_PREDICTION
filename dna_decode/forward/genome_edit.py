"""Genome-level forward edit: a single nucleotide change in a coding sequence (CDS) -> the codon it hits ->
the amino-acid change -> the protein-level phenotype predictor (variant_effect.predict_effect).

This lifts the forward cell's INPUT from a protein point-mutation ('M69L') to a real genome edit
(CDS position + ref base + alt base), classifying it as SILENT (synonymous), NONSENSE (premature stop),
or MISSENSE, and routing missense/nonsense to the Regime-B predictor. Reference base is verified against
the CDS (a mismatch fails LOUDLY — the coordinate-integrity discipline), and if a protein sequence is
supplied the translated codon must match it.
"""
from __future__ import annotations

from dataclasses import dataclass

from .variant_effect import ForwardPrediction, predict_effect

# Standard genetic code (NCBI table 1; internal-codon translation is identical to bacterial table 11 —
# they differ only in alternative START codons, which do not affect a point-substitution's residue call).
_CODON = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M", "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T", "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*", "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K", "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W", "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R", "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


def translate_codon(codon: str) -> str:
    """Standard-code translation of one 3-nt codon -> one-letter AA ('*' = stop). Raises on a bad codon."""
    c = codon.upper()
    if len(c) != 3 or any(b not in "ACGT" for b in c):
        raise ValueError(f"not a valid DNA codon: {codon!r}")
    return _CODON[c]


def parse_hgvs_c(spec: str) -> tuple[int, str, str]:
    """Parse an HGVS coding-DNA substitution -> (1-based CDS position, ref base, alt base). PURE.

    Accepts `c.205G>A` (canonical) and the bare `205G>A`. Deliberately NARROW: only single-base
    substitutions on the coding sequence. Anything else — an indel (`c.205delG`), a genomic/protein
    prefix (`g.` / `p.`), an interval, a `*`/`-` UTR offset, an intronic `+`/`-` offset — is REFUSED
    rather than silently coerced, because every one of those has a different coordinate meaning and a
    wrong coordinate is the failure mode this whole module exists to catch loudly.
    """
    import re
    s = str(spec).strip()
    m = re.fullmatch(r"(?:c\.)?(\d+)([ACGTacgt])>([ACGTacgt])", s)
    if not m:
        raise ValueError(
            f"not an HGVS coding substitution: {spec!r}. Expected `c.<pos><REF>><ALT>` (e.g. c.205G>A). "
            f"Only single-base CDS substitutions are supported — indels, `g.`/`p.` coordinates, and "
            f"intronic/UTR offsets are refused rather than guessed."
        )
    pos = int(m.group(1))
    if pos < 1:
        raise ValueError(f"HGVS position must be 1-based and positive; got {pos}")
    return pos, m.group(2).upper(), m.group(3).upper()


def translate_cds(cds: str) -> str:
    """Translate a full CDS -> protein, dropping ONE trailing stop. PURE.

    A codon containing an ambiguity base (N/R/Y/...) becomes 'X' rather than raising, so a real genome
    with Ns still translates; the edited codon itself is still held to strict ACGT by `cds_point_edit`.
    An INTERNAL stop is left in the returned string as '*' — it is not silently trimmed, so a
    wrong-frame or wrong-strand CDS surfaces as a loud WT-residue mismatch downstream instead of
    quietly producing a truncated protein that verifies against nothing.
    """
    s = "".join(str(cds).split()).upper()
    if len(s) % 3 != 0:
        raise ValueError(
            f"CDS length {len(s)} is not a multiple of 3 — not a coding sequence in frame "
            f"(refusing rather than guessing the reading frame)"
        )
    aas = [_CODON.get(s[i:i + 3], "X") for i in range(0, len(s), 3)]
    if aas and aas[-1] == "*":
        aas.pop()
    return "".join(aas)


def cds_record_key(rec: dict) -> str:
    """The key `annotations.extract_cds_sequences` files a CDS row under. Kept in ONE place so the
    lookup cannot drift from the extractor's own rule."""
    return (rec.get("gene_id") or rec.get("locus_tag")
            or f"{rec.get('seqid')}:{rec.get('start')}-{rec.get('end')}")


def _ambiguous_message(gene: str, field: str, hits: list[dict]) -> str:
    """Explain WHY a gene name matched several CDS rows, and name the field that actually separates them.

    Two genuinely different situations, and telling them apart matters (both occur in the real E. coli
    K-12 MG1655 RefSeq GFF3):

      * ALTERNATIVE PRODUCTS -- one locus, several protein accessions (mrcB: NP_414691.1 + YP_010051172.1).
        `locus_tag` is IDENTICAL across them, so advising a locus_tag would be useless; `gene_id` is what
        separates them.
      * JOINED / MULTI-SEGMENT CDS -- one accession spanning several GFF rows (dnaX YP_009518751.1 covers
        492092-493375 and 493375-493386, the -1 programmed ribosomal frameshift that makes the tau/gamma
        subunits). Such a product is NOT one contiguous CDS, so this path refuses it rather than decoding
        an arbitrary segment. NOTE that `annotations.extract_cds_sequences` keys by gene_id and would
        silently keep only the LAST segment (here a 12-nt "CDS") -- which is exactly why selection fails
        closed here instead of trusting that lookup.
    """
    from collections import Counter
    ids = [str(h.get("gene_id") or "") for h in hits]
    counts = Counter(i for i in ids if i)
    joined = sorted(k for k, v in counts.items() if v > 1)
    single = sorted(k for k, v in counts.items() if v == 1)

    span = lambda h: f"{h.get('seqid')}:{h.get('start')}-{h.get('end')}"          # noqa: E731
    lines = [f"{gene!r} matches {len(hits)} CDS features by {field}:"]
    for h in hits[:6]:
        lines.append(f"    {h.get('gene_id') or '?'}  {span(h)}  (locus_tag {h.get('locus_tag') or '-'})")
    if len(hits) > 6:
        lines.append(f"    ... and {len(hits) - 6} more")

    if field != "gene_id" and single:
        lines.append(f"  -> alternative products of one locus. Re-run with --gene <gene_id>, e.g. {single[0]}.")
    if joined:
        lines.append(f"  -> {', '.join(joined)} is a JOINED multi-segment CDS (e.g. a programmed "
                     f"frameshift); it is not one contiguous coding sequence and this path cannot decode "
                     f"a `c.` coordinate on it. Supply the assembled CDS via --cds-fasta instead.")
    if field == "gene_id" and not joined:
        lines.append("  -> the same gene_id appears on several rows; supply the CDS via --cds-fasta.")
    return "\n".join(lines)


def select_gene_cds(records: list[dict], gene: str) -> dict:
    """Pick the ONE CDS annotation row naming `gene`. PURE (operates on plain dicts, no pandas).

    Matches `gene_symbol` FIRST -- that is the CROSS-STRAIN identifier (`gyrA`) -- then `locus_tag`, then
    `gene_id`. This order is load-bearing: `gene_id` is strain-unique by construction (`gene-b0001`), so
    resolving a user's `gyrA` against it is the documented 0%-overlap trap.

    Raises on 0 matches (listing nearby symbols) and on >1 (a multi-copy gene -- e.g. the real 7-copy
    tandem blaTEM array in the genome-map spike -- where silently taking the first copy would decode an
    arbitrary one of them). Matching is case-insensitive; ties are refused, never guessed.
    """
    want = str(gene).strip().lower()
    cds = [r for r in records if str(r.get("type", "")) == "CDS"]
    for field in ("gene_symbol", "locus_tag", "gene_id"):
        hits = [r for r in cds if str(r.get(field) or "").lower() == want]
        if len(hits) == 1:
            return hits[0]
        if len(hits) > 1:
            raise ValueError(_ambiguous_message(gene, field, hits))
    known = sorted({str(r.get("gene_symbol")) for r in cds if r.get("gene_symbol")})
    near = [s for s in known if want[:3] and s.lower().startswith(want[:3])][:8]
    raise ValueError(
        f"no CDS feature named {gene!r} (searched gene_symbol, locus_tag, gene_id over {len(cds)} CDS "
        f"rows; {len(known)} carry a gene_symbol)."
        + (f" Did you mean: {', '.join(near)}?" if near else
           " Note that RefSeq GFF3 populates gene_symbol for only ~11% of CDSs -- try a locus_tag.")
    )


@dataclass
class GenomeEditPrediction:
    nt_pos: int                 # 1-based CDS position of the edited base
    ref_base: str
    alt_base: str
    aa_pos: int                 # 1-based residue the codon encodes
    wt_aa: str
    alt_aa: str
    wt_codon: str
    alt_codon: str
    consequence: str            # "silent" | "missense" | "nonsense"
    aa_mutation: str | None     # e.g. "M69L" (None for silent)
    protein_prediction: ForwardPrediction | None   # None for silent (no protein change)

    def as_dict(self) -> dict:
        d = {
            "nt_pos": self.nt_pos, "ref_base": self.ref_base, "alt_base": self.alt_base,
            "aa_pos": self.aa_pos, "wt_aa": self.wt_aa, "alt_aa": self.alt_aa,
            "wt_codon": self.wt_codon, "alt_codon": self.alt_codon, "consequence": self.consequence,
            "aa_mutation": self.aa_mutation,
            "protein_prediction": (self.protein_prediction.as_dict() if self.protein_prediction else None),
        }
        return d


def cds_point_edit(cds: str, nt_pos: int, ref_base: str, alt_base: str) -> dict:
    """Resolve a 1-based CDS base substitution to its codon consequence. Verifies ref_base against the CDS."""
    if nt_pos < 1 or nt_pos > len(cds):
        raise ValueError(f"nt_pos {nt_pos} out of range for CDS length {len(cds)}")
    ref_base, alt_base = ref_base.upper(), alt_base.upper()
    if alt_base not in "ACGT" or ref_base not in "ACGT":
        raise ValueError(f"ref/alt must be single DNA bases; got {ref_base!r}->{alt_base!r}")
    idx = nt_pos - 1
    have = cds[idx].upper()
    if have != ref_base:
        raise ValueError(f"REF mismatch at CDS pos {nt_pos}: sequence has {have!r}, edit asserts {ref_base!r} "
                         f"(coordinate error — refusing)")
    codon_no = idx // 3                       # 0-based codon index
    within = idx % 3                          # 0..2 position within codon
    cstart = codon_no * 3
    wt_codon = cds[cstart:cstart + 3].upper()
    if len(wt_codon) != 3:
        raise ValueError(f"edit at nt_pos {nt_pos} falls in an incomplete terminal codon (CDS not a multiple of 3?)")
    alt_codon = wt_codon[:within] + alt_base + wt_codon[within + 1:]
    return {
        "aa_pos": codon_no + 1, "within_codon": within,
        "wt_codon": wt_codon, "alt_codon": alt_codon,
        "wt_aa": translate_codon(wt_codon), "alt_aa": translate_codon(alt_codon),
    }


def predict_genome_edit(cds: str, nt_pos: int, ref_base: str, alt_base: str, *,
                        protein_seq: str | None = None, protein: str = "protein",
                        phenotype_axis: str = "molecular fitness (DMS-measured)",
                        method: str = "blosum62", esm_table: dict | None = None) -> GenomeEditPrediction:
    """Full genome-edit -> phenotype path: CDS base substitution -> codon -> AA change -> Regime-B predictor.

    - SILENT (synonymous): no protein change -> protein_prediction=None (predicted neutral at the protein level).
    - NONSENSE / MISSENSE: build the AA mutation and delegate to predict_effect (which re-verifies the WT
      residue against `protein_seq` if supplied — double coordinate check: CDS ref base AND translated WT AA).
    """
    info = cds_point_edit(cds, nt_pos, ref_base, alt_base)
    aa_pos, wt_aa, alt_aa = info["aa_pos"], info["wt_aa"], info["alt_aa"]

    if wt_aa == alt_aa:
        return GenomeEditPrediction(nt_pos, ref_base.upper(), alt_base.upper(), aa_pos, wt_aa, alt_aa,
                                    info["wt_codon"], info["alt_codon"], "silent", None, None)

    consequence = "nonsense" if alt_aa == "*" else "missense"
    aa_mut = f"{wt_aa}{aa_pos}{'*' if alt_aa == '*' else alt_aa}"
    pred = predict_effect(protein_seq or "", aa_mut, protein=protein, phenotype_axis=phenotype_axis,
                          method=method, esm_table=esm_table)
    return GenomeEditPrediction(nt_pos, ref_base.upper(), alt_base.upper(), aa_pos, wt_aa, alt_aa,
                                info["wt_codon"], info["alt_codon"], consequence, aa_mut, pred)
