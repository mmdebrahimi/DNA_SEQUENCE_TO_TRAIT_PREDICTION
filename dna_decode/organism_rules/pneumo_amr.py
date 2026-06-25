"""S. pneumoniae AMR decoder cell — gene-presence determinant rule (NON-FROZEN).

Organism-routed, deterministic, mirrors the `tb_amr` / TMP-SMX overlay branding. Scope v0 = the CLEAN
GENE-PRESENCE drugs only:
  - Macrolide (erythromycin): R iff erm(B) OR mef(A)/msr(D) present.
  - Tetracycline: R iff tet(M) OR tet(O) present.

DELIBERATELY OUT OF SCOPE (deferred): β-lactams (penicillin/cephalosporins) — pneumococcal β-lactam
resistance is PBP-mutation-MIC (pbp1a/2b/2x → MIC, the CDC regression method), NOT gene-presence, and is
breakpoint-AMBIGUOUS (meningitis/non-meningitis/oral flip the R/S of one MIC). That needs a PBP-typing
engine + a pneumo breakpoint set — a separate, larger build (`wiki/pneumo_amr_cell_go_no_go_2026-06-25.md`).

HONESTY (load-bearing):
  - Branding KNOWLEDGE_BASELINE / organism_routed — never confused with the FROZEN E. coli deployed surface.
  - The determinant catalog (erm/mef/tet) is the GPS/CDC open-source set. For GENE-PRESENCE genes the call
    is plain BLAST presence, which AMRFinder replicates near-identically → the faithful-to-tool gap is small
    (unlike the PBP regression). A truly independent number runs OUR AMRFinder on the assemblies; the
    GPS-determinant baseline is `RULE_STATUS`-tagged accordingly by the scorer.
  - Validated vs WET-LAB measured AST (GPS Poland cohort, disc/agar): macrolide acc 0.961 (sens 0.968 /
    spec 0.954, n=127); tetracycline acc 0.932 (sens 0.964 / spec 0.913, n=74). See
    `scripts/pneumo_amr_validate.py` + `wiki/pneumo_amr_validation_result_2026-06-25.md`.

This module makes PER-ISOLATE calls from a determinant-token set. It does NOT call determinants itself
(that is AMRFinder/blastn, organism `Streptococcus_pneumoniae`) nor compute cohort sens/spec (the scorer).
"""
from __future__ import annotations

from dataclasses import dataclass, field

RULE_STATUS = "KNOWLEDGE_BASELINE"
RULE_SCOPE = "organism_routed"
ORGANISM = "Streptococcus_pneumoniae"

R, S = "R", "S"

# Gene-presence determinant rule: drug -> the token substrings (lowercased) that, if ANY is present in the
# isolate's determinant set, call R. Tokens are matched as substrings so allele suffixes don't matter
# (e.g. "ermB_16_X82819" contains "ermb"; "mefA_10_AF376746" contains "mefa").
GENE_PRESENCE_RULES: dict[str, tuple[str, ...]] = {
    "erythromycin": ("ermb", "mefa", "msrd"),   # macrolide: erm(B) target methylation + mef/msr efflux
    "tetracycline": ("tetm", "teto"),           # ribosomal protection
}
# Macrolide-class synonyms route to the erythromycin rule (the measured/representative macrolide).
DRUG_ALIASES = {"azithromycin": "erythromycin", "clarithromycin": "erythromycin", "macrolide": "erythromycin"}

SUPPORTED_DRUGS = tuple(GENE_PRESENCE_RULES)


@dataclass(frozen=True)
class PneumoDrugCall:
    drug: str
    prediction: str                          # R / S
    organism: str = ORGANISM
    rule_status: str = RULE_STATUS
    rule_scope: str = RULE_SCOPE
    matched_tokens: tuple[str, ...] = ()      # which determinant tokens fired
    rule_tokens: tuple[str, ...] = field(default_factory=tuple)


def resolve_drug(drug: str) -> str | None:
    """Map an input drug name to a supported gene-presence rule key, or None if out of scope."""
    d = (drug or "").strip().lower()
    d = DRUG_ALIASES.get(d, d)
    return d if d in GENE_PRESENCE_RULES else None


def call_drug(drug: str, determinant_tokens) -> PneumoDrugCall | None:
    """Per-isolate gene-presence call. `determinant_tokens` = an iterable of determinant strings (e.g. the
    AMRFinder/GPS determinant calls for this isolate). Returns None if the drug is out of v0 scope
    (β-lactams, FQ-without-determinant, etc.)."""
    key = resolve_drug(drug)
    if key is None:
        return None
    tokens = GENE_PRESENCE_RULES[key]
    hay = " ".join(str(t) for t in determinant_tokens).lower()
    matched = tuple(t for t in tokens if t in hay)
    return PneumoDrugCall(drug=key, prediction=(R if matched else S),
                          matched_tokens=matched, rule_tokens=tokens)
