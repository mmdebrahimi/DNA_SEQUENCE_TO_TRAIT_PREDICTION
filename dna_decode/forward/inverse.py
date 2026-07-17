"""The INVERSE of the forward cell: given a desired effect, propose the edit — RANK only, never a dose.

`dna_decode.forward` answers *edit -> effect*. This answers *effect -> edit*, which is what a design loop
actually wants. It uses the DMS-validated forward oracle as LABEL-FREE ground truth: no phenotype label is
ever consulted, which is the move that dodges this project's binding constraint (labels, not models).

WHAT IT ANSWERS, EXACTLY -- and the narrowness is the point:

    "propose the edits at the p-th percentile of predicted damage, among the edits reachable here"

It ranks. **It does not dose.** It can say *near the top of the damaging tail*; it can never say
*fold-change 4.2*. That limit is measured, not modest -- see WHY below.

WHY NOT A MAGNITUDE (this shape was tested and rejected, `scripts/forward_inverse_deployable.py`):
hitting a target EFFECT requires a score->effect calibrator fit on the TARGET protein's own DMS -- and if
you have that protein's DMS you already know every variant's effect, so you do not need an inverse. The
calibrator cannot be borrowed from another protein either: the assays share no scale (CcdB's whole measured
range [-9.00,-2.00] lies BELOW TEM-1's minimum -3.56, so a TEM-1 calibrator cannot express a CcdB value at
all -- impossible by construction, not merely inaccurate). And on blaTEM the conformal dosage interval is
informative 0/6 splits: it BRACKETS the target while proving nothing, because split-conformal coverage
holds even for a useless model. The rank question needs NONE of that -- no calibrator, no DMS, no label --
which is why it is the one that ships.

MEASURED (wiki/forward_inverse_{roundtrip,sweep,deployable}_2026-07-1{6,7}):
  * the rank inverse beats an exact no-oracle null on **4/4** usable proteins (E. coli, human, yeast,
    Arabidopsis), ~2-5 percentile points of error with 5 proposals;
  * graded non-circularly against MEASURED wet-lab DMS, never against the model's own re-score;
  * **the learned oracle earns its keep over plain BLOSUM62 on only 3/4** -- so `method="blosum62"` (the
    default: instant, wheel-only, no GPU) is often the RIGHT answer, not a fallback;
  * **utility does NOT track forward rank** (PTEN 0.5185 earns its keep; RL40A 0.5190 does not) -- so a
    good leaderboard Spearman does NOT license skipping the per-protein check.

SCOPE: Regime B (molecular fitness/stability) ONLY. Never clinical resistance -- the same scorer class is
BELOW CHANCE there (ESM2 0.454 vs the curated catalogue's 0.926); use the frozen `dna-decode amr`.
Research use only. Pure-python, wheel-only, offline, deterministic.
"""
from __future__ import annotations

from dataclasses import dataclass, field

from .variant_effect import blosum62_score

# The standard genetic code -- used ONLY to restrict proposals to single-nt-reachable edits when a CDS is
# supplied. Without a CDS the candidate space is every substitution (protein design, not genome editing).
CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M", "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T", "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*", "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K", "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W", "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R", "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}
AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"
BASES = "ACGT"

# Measured evidence that ships with every call, so a proposal cannot be read as more than it is.
# The EXACT proteins the evidence was earned on. A user's protein is almost never one of these, and the
# measured finding is that utility does NOT transfer by rank (PTEN 0.5185 earns its keep; RL40A 0.5190 does
# not) -- so "not gated for YOUR protein" is a real caveat, not boilerplate, and it ships per call.
VALIDATED_PROTEINS = ("TEM-1 beta-lactamase (E. coli)", "PTEN (human)", "RL40A/ubiquitin (yeast)",
                      "SR43C (Arabidopsis)")
EVIDENCE = {
    "validated_on": "ProteinGym DMS (measured wet-lab per-variant fitness)",
    "validated_proteins": list(VALIDATED_PROTEINS),
    # HONESTY: the 4/4 was the LEARNED method (ESM) on 4 hand-picked proteins. The SHIPPED CLI default is
    # blosum62 -- and at scale (N=200 ProteinGym, wiki/proteingym_inverse_sweep_2026-07-17.md) blosum62 is
    # materially better than a random pick on only 13.5% of proteins, has any positive edge on 59%, and is
    # frequently WORSE than guessing. Do NOT let the shipped default inherit the learned method's number.
    "esm_rank_inverse_beats_null": "4/4 hand-picked proteins (LEARNED method; needs a GPU per protein)",
    "shipped_blosum62_default_beats_null_at_scale": ("13.5% materially / 59% any-edge over N=200 ProteinGym "
                                                     "assays -- often WORSE than random. The wheel-only "
                                                     "default is NOT a reliable design tool on its own"),
    "typical_error_top5": "~2-5 percentile points ON THE PROTEINS WHERE IT WORKS (not the average)",
    "artifact": "wiki/proteingym_inverse_sweep_2026-07-17.md (blosum62 N=200) + "
                "wiki/forward_inverse_deployable_2026-07-17.md (esm N=4)",
}
UNSUPPORTED_CLAIMS = (
    "a MAGNITUDE / dose / fold-change for a proposed edit -- the calibrator needed for that requires the "
    "target protein's own DMS (which would make the inverse unnecessary), and the conformal interval is "
    "uninformative even where it brackets",
    "a CLINICAL resistance call -- Regime A; this scorer class is below chance there. Use `dna-decode amr`",
    "an ORGANISM-level phenotype -- Regime C, a closed negative; the forward router abstains",
)


class InverseError(ValueError):
    """Raised on an unusable input -- never a silent wrong proposal."""


@dataclass
class ProposedEdit:
    mutation: str                # e.g. "M69L"
    pos: int
    wt: str
    alt: str
    score: float                 # the oracle's raw score (NOT an effect size)
    score_percentile: float      # where this edit sits in the predicted-damage ordering [0,1]
    codon_change: str | None = None   # e.g. "ATG->CTG" when a CDS was supplied

    def as_dict(self) -> dict:
        return {"mutation": self.mutation, "pos": self.pos, "wt": self.wt, "alt": self.alt,
                "score": round(self.score, 4), "score_percentile": round(self.score_percentile, 4),
                "codon_change": self.codon_change}


@dataclass
class InverseProposal:
    target_percentile: float
    method: str
    candidate_space: str
    n_candidates: int
    proposals: list[ProposedEdit] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)
    gated_for_this_protein: bool = False

    def as_dict(self) -> dict:
        return {
            "tool": "dna_decode.forward.inverse", "version": "v0",
            "regime": "B_molecular", "claim": "predicted-damage RANK only -- never a dose",
            "target_percentile": self.target_percentile, "method": self.method,
            "candidate_space": self.candidate_space, "n_candidates": self.n_candidates,
            "proposals": [p.as_dict() for p in self.proposals],
            "evidence": dict(EVIDENCE),
            "gated_for_this_protein": self.gated_for_this_protein,
            "does_not_support": list(UNSUPPORTED_CLAIMS),
            "notes": self.notes,
            "research_use_only": True,
        }


def enumerate_candidates(protein_seq: str, cds: str | None = None) -> tuple[list[tuple[str, int, str, str, str | None]], str]:
    """Every edit to consider. With a CDS -> only single-nt-reachable ones (real genome editing);
    without -> every substitution (protein design). Returns (candidates, space_label)."""
    seq = protein_seq.strip().upper()
    if not seq or any(c not in AMINO_ACIDS for c in seq):
        raise InverseError("protein_seq must be non-empty and contain only the 20 standard amino acids")

    out: list[tuple[str, int, str, str, str | None]] = []
    if cds is None:
        for i, wt in enumerate(seq):
            for alt in AMINO_ACIDS:
                if alt != wt:
                    out.append((f"{wt}{i + 1}{alt}", i + 1, wt, alt, None))
        return out, "all_substitutions_protein_level_no_CDS"

    cds = cds.strip().upper()
    if len(cds) < len(seq) * 3:
        raise InverseError(f"CDS is too short ({len(cds)} nt) for a {len(seq)}-residue protein")
    for i, wt in enumerate(seq):
        codon = cds[i * 3: i * 3 + 3]
        if CODON_TABLE.get(codon) != wt:
            # Coordinate gate: fail loudly rather than propose an edit in the wrong frame.
            raise InverseError(
                f"CDS codon {i + 1} is {codon!r} (={CODON_TABLE.get(codon)}) but the protein has {wt!r} at "
                f"that position -- the CDS does not translate to this protein; refusing to propose edits")
        for k in range(3):
            for b in BASES:
                if b == codon[k]:
                    continue
                new_codon = codon[:k] + b + codon[k + 1:]
                alt = CODON_TABLE.get(new_codon)
                if alt and alt != "*" and alt != wt:
                    out.append((f"{wt}{i + 1}{alt}", i + 1, wt, alt, f"{codon}->{new_codon}"))
    return out, "single_nt_accessible_from_CDS"


def propose_edits(protein_seq: str, target_percentile: float, *, top_k: int = 5,
                  method: str = "blosum62", esm_table: dict | None = None,
                  cds: str | None = None, diverse: bool = True) -> InverseProposal:
    """Propose the `top_k` edits nearest the `target_percentile` of PREDICTED damage.

    percentile 0.0 = the most damaging predicted edit; 1.0 = the most tolerated. (Scores are ordered
    ascending, and a lower oracle score means a more disruptive substitution.)

    `diverse=True` (DEFAULT) returns at most one edit per residue. This is not cosmetic: BLOSUM62 is
    heavily quantized -- 1,874 real blaTEM single-nt candidates take only SEVEN distinct scores, the
    largest tie group holding 383 -- so a plain window returns k edits from one tie group ordered by
    position, i.e. k shots at the SAME residue (measured: D48Y, D48V, L49H, L49P, N50I). For a loop that
    assays k proposals that wastes k-1 of them. The cost was MEASURED, not assumed
    (`scripts/forward_inverse_deployable.py`): diversity costs ESM ~0 (+/-0.005 percentile points) and
    IMPROVES BLOSUM by up to 7 points (RL40A -0.0701), because forcing distinct residues escapes the tie
    groups. Free-to-beneficial, so it is the default. Pass diverse=False for the plain window.
    """
    if not 0.0 <= target_percentile <= 1.0:
        raise InverseError(f"target_percentile must be in [0,1], got {target_percentile}")
    if top_k < 1:
        raise InverseError("top_k must be >= 1")

    cands, space = enumerate_candidates(protein_seq, cds)
    scored: list[ProposedEdit] = []
    for mut, pos, wt, alt, codon_change in cands:
        if method == "blosum62":
            s = blosum62_score(wt, alt)
        elif method == "esm2":
            if esm_table is None:
                raise InverseError("method='esm2' needs a precomputed esm_table {pos: {aa: logprob}} -- the "
                                   "model runs ONCE per protein; see dna_decode/forward/README.md")
            col = esm_table.get(pos) or esm_table.get(str(pos))
            if not col or alt not in col or wt not in col:
                continue
            s = col[alt] - col[wt]
        else:
            raise InverseError(f"unknown method {method!r}; expected 'blosum62' or 'esm2'")
        scored.append(ProposedEdit(mut, pos, wt, alt, s, 0.0, codon_change))

    if not scored:
        raise InverseError("no scorable candidate edits")
    scored.sort(key=lambda p: p.score)
    n = len(scored)
    for i, p in enumerate(scored):
        p.score_percentile = i / max(1, n - 1)

    idx = min(n - 1, max(0, int(round(target_percentile * (n - 1)))))
    if diverse:
        picks, seen = [], set()
        for off in range(n):
            for j in ((idx,) if off == 0 else (idx - off, idx + off)):
                if 0 <= j < n and scored[j].pos not in seen:
                    seen.add(scored[j].pos)
                    picks.append(scored[j])
                    if len(picks) == top_k:
                        break
            if len(picks) == top_k:
                break
    else:
        lo = max(0, min(n - top_k, idx - top_k // 2))
        picks = scored[lo: lo + top_k]

    n_tied = sum(1 for p in scored if p.score == scored[idx].score)
    notes = [
        f"proposals are the {len(picks)} edits nearest percentile {target_percentile:.2f} of PREDICTED "
        f"damage among {n} candidates -- a RANK, not a dose",
        f"{len({p.pos for p in picks})} distinct residues among {len(picks)} proposals"
        + (" (diverse=True: one edit per residue)" if diverse else " (diverse=False: plain window)"),]
    if n_tied > top_k:
        notes.append(
            f"{n_tied} candidates TIE at this score ({n_tied / n:.0%} of the pool) -- {method} is coarsely "
            f"quantized here, so the pick within the tie group is arbitrary. Diversity is measured "
            f"free-to-beneficial precisely because it escapes ties; treat the k proposals as "
            f"interchangeable draws from the tie group, not a ranking among themselves")
    if method == "blosum62":
        notes.append(
            "NOT GATED FOR YOUR PROTEIN + this is the WEAK default: at scale (N=200 ProteinGym) the "
            "blosum62 rank inverse beats a random pick MATERIALLY on only 13.5% of proteins and is often "
            "WORSE than guessing. The 4/4 headline was ESM (a GPU method) on 4 hand-picked proteins. So "
            "treat these proposals as a cheap FIRST PASS, not a validated ranking, and run "
            "scripts/forward_inverse_deployable.py on a DMS assay for YOUR protein before trusting them.")
    else:
        notes.append(
            "NOT GATED FOR YOUR PROTEIN: the 4/4 ESM evidence was earned on " + ", ".join(VALIDATED_PROTEINS) +
            ". Inverse utility does NOT transfer by rank quality (PTEN forward-rank 0.5185 earns its keep "
            "while RL40A 0.5190 does not — a gap of 0.0005, opposite verdicts), so a good Spearman does not "
            "license skipping the per-protein check (scripts/forward_inverse_deployable.py on your DMS).")
    notes += [
        "measured guidance: assay all k proposals and keep the best. Single-shot (top-1) is ~4x worse than "
        "best-of-5 -- the loop assumes you can measure what it proposes",
    ]
    if method == "blosum62":
        notes.append("method=blosum62 (deterministic, no model). Measured: the learned oracle earns its keep "
                     "on only 3/4 proteins, so this is often the right answer rather than a fallback")
    if cds is None:
        notes.append("no CDS supplied -> candidates are ALL substitutions (protein design). Pass a CDS to "
                     "restrict to single-nt-reachable edits (real genome editing)")
    return InverseProposal(target_percentile, method, space, n, picks, notes)
