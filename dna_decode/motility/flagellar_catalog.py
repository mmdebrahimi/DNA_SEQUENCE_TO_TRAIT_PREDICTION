"""E. coli flagellar-motility catalog: gene presence -> can the cell SWIM? (the first NON-metabolic cell.)

THE BIOLOGY (the rule this encodes): flagellar motility needs a COMPLETE, functional flagellar system, and
that system is built as a strict regulatory + structural cascade. Missing ANY of the load-bearing modules ->
no working flagellum -> non-motile. The modules (each a decision gate):

  1. MASTER REGULATOR (flhDC = flhD AND flhC) -- the class-1 switch. Without it NO flagellar gene is
     transcribed at all. This is the flagellar analog of a master gate (the K-12 IS-insertion that makes
     lab strains non-motile hits exactly here -- see the honest scope note).
  2. SIGMA-28 (fliA) -- the class-3 sigma factor; without it the filament + motor genes are not expressed.
  3. FLAGELLIN (fliC OR fljB) -- the filament protein. No flagellin -> no filament.
  4. MOTOR (motA AND motB) -- the stator/torque generator. A flagellum that can't rotate -> non-motile.
  5. BASAL BODY + EXPORT (fliF, fliG, flhA, fliI) -- the MS-ring/switch + the type-III export apparatus that
     secretes the flagellar proteins. No basal body -> nothing to build on.

MOTILE iff all 5 modules are satisfied. CHEMOTAXIS (cheA/W/Y/Z) is reported SEPARATELY and does NOT gate
motility: a che-mutant still SWIMS (it just tumbles randomly / can't chase a gradient) -- gating swimming on
chemotaxis would be a biology error.

VALIDATED against literature-known anchors: E. coli K-12 MG1655 (motile, all modules) + Salmonella
Typhimurium (motile) vs Shigella flexneri (NON-motile -- its flagellar genes are pseudogenes/deleted) + any
flhDC / fliC / motAB knockout (non-motile). KNOWLEDGE_BASELINE tier (a curated catalog vs literature
anchors), NOT a big measured cohort -- no free genome-keyed swim-plate cohort is fetchable.

HONEST SCOPE (load-bearing):
  - v0 = the flagellar SWIMMING decision (motile / non-motile) from gene PRESENCE. Direction only, NOT
    swim speed / rate / gradient performance.
  - Presence-based: it CANNOT see a point mutation / IS-insertion that silently inactivates a gene that is
    still ANNOTATED PRESENT (the classic K-12 flhD IS1 case -> a present-but-dead flhD the rule mis-calls
    motile). A v0.1 sequence-integrity genome-mode is the named follow-on, deliberately not fabricated here.
  - Type-IV-pilus twitching motility + gliding + swarming-specific regulators are OUT of v0 (flagellar swim only).

Pure-python, wheel-only, offline, deterministic. NON-frozen (the frozen AMR decoder surface is untouched).
"""
from __future__ import annotations

from dataclasses import dataclass


class MotilityInputError(ValueError):
    """Raised on malformed input (e.g. no gene evidence supplied)."""


# module name -> (combinator, required genes, role). combinator "all"=every gene, "any"=>=1 gene.
MOTILITY_MODULES: dict[str, tuple[str, tuple[str, ...], str]] = {
    "master_regulator": ("all", ("flhD", "flhC"), "class-1 flagellar master (no flagella without it)"),
    "sigma28": ("all", ("fliA",), "sigma-28; class-3 expression of filament + motor"),
    "flagellin": ("any", ("fliC", "fljB"), "flagellar filament protein"),
    "motor": ("all", ("motA", "motB"), "stator / torque generator (flagellar rotation)"),
    "basal_export": ("all", ("fliF", "fliG", "flhA", "fliI"), "MS-ring/switch + type-III export apparatus"),
}
# reported SEPARATELY -- NOT a motility gate (a che-mutant still swims)
CHEMOTAXIS: tuple[str, tuple[str, ...], str] = (
    "all", ("cheA", "cheW", "cheY", "cheZ"), "directed chemotaxis (not required for swimming)")


def catalog_genes() -> set[str]:
    """Every gene the catalog looks at (motility modules + chemotaxis)."""
    genes: set[str] = set()
    for _, gs, _ in MOTILITY_MODULES.values():
        genes.update(gs)
    genes.update(CHEMOTAXIS[1])
    return genes


def _module_satisfied(combinator: str, genes: tuple[str, ...], present: set[str]) -> bool:
    hits = [g for g in genes if g in present]
    return (len(hits) == len(genes)) if combinator == "all" else (len(hits) >= 1)


@dataclass(frozen=True)
class MotilityCall:
    motile: bool
    verdict: str
    module_status: dict[str, bool]        # module -> satisfied
    missing_modules: tuple[str, ...]
    chemotaxis_competent: bool
    note: str


def call_motility(present_genes) -> MotilityCall:
    """Gene presence -> flagellar-motility call. `present_genes` = an iterable of present gene symbols.

    MOTILE iff all 5 flagellar modules are satisfied. Chemotaxis is reported but does NOT gate motility.
    """
    present = {str(g).strip() for g in present_genes if str(g).strip()}
    if not present:
        raise MotilityInputError("no gene evidence supplied (give present gene symbols)")

    status = {name: _module_satisfied(comb, genes, present)
              for name, (comb, genes, _) in MOTILITY_MODULES.items()}
    missing = tuple(name for name, ok in status.items() if not ok)
    motile = not missing
    chemo = _module_satisfied(*CHEMOTAXIS[:2], present) if motile else False

    if motile:
        verdict = "MOTILE"
        note = ("all 5 flagellar modules present; "
                + ("chemotaxis-competent" if chemo else "swims but chemotaxis incomplete (random tumbling)"))
    else:
        verdict = f"NON-MOTILE (missing module(s): {', '.join(missing)})"
        note = ("presence-based -- a gene that is annotated PRESENT but silently inactivated "
                "(e.g. the K-12 flhD IS-insertion) would be mis-called; sequence-mode is v0.1")
    return MotilityCall(
        motile=motile, verdict=verdict, module_status=status,
        missing_modules=missing, chemotaxis_competent=chemo, note=note,
    )
