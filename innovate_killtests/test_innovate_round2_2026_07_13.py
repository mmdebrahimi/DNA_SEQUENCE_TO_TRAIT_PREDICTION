"""Round-2 kill-tests for /innovate on the genome world model (2026-07-13, second deeper pass).

Round 1 reframed the PRODUCT (federation, interventional lever, self-awareness flag). Round 2 hunts a
concretely novel MECHANISM: the curated catalog is ALREADY an interventional (edit->effect) knowledge
base. Anti-theater contract: each test PASSES iff its OWN candidate claim is FALSIFIED (pytest exit-0 =
KILLED; exit-1 = SURVIVED). Every kill-test is UNIQUE to its claim (no shared falsifier — the
2026-07-13 dogfound lesson). Read-only over committed artifacts; frozen surface untouched.
"""
from __future__ import annotations

import collections
import json
import re
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))


def test_g1r2_generative_catalog_falsified():
    """g1r2-generative-catalog (make the catalog GENERATIVE — predict an unseen edit's effect by
    interpolating across a site's catalogued substitutions) DEAD iff there is NO interpolation
    substrate: fewer than 2 DRM positions carry >=2 catalogued substitutions. Committed: 4 positions
    (106,181,188,190) -> substrate exists -> SURVIVES."""
    from dna_decode.data import hiv_amr as H
    pos = collections.Counter()
    for m in H.NNRTI_RT_MAJOR_DRMS:
        mm = re.match(r"([A-Z])(\d+)([A-Z])", m)
        if mm:
            pos[int(mm.group(2))] += 1
    multi = [p for p, c in pos.items() if c >= 2]
    assert len(multi) < 2, (
        f"NOT falsified: {len(multi)} DRM positions carry >=2 catalogued substitutions ({sorted(multi)}) "
        "-> the catalog HAS interventional density to interpolate across (generative substrate exists)")


def test_g8r2_cross_organism_transfer_falsified():
    """g8r2-cross-organism-transfer (a target-site catalog is protein-family-level, so a validated
    catalog TRANSFERS to a new organism with zero new labels) DEAD iff cross-organism transfer does
    NOT already hold: fewer than 3 DISTINCT-genus organisms share a SCORED cipro gyrA/parC-QRDR cell.
    Committed: 4 genera (E. coli, Klebsiella, Campylobacter, Salmonella) all SCORED -> SURVIVES."""
    from dna_decode.data.shipped_decoder_surface import SHIPPED_DECODER_SURFACE as S
    cipro_orgs = {r[0] for r in S if r[1] == "ciprofloxacin"}
    assert len(cipro_orgs) < 3, (
        f"NOT falsified: {len(cipro_orgs)} distinct organisms share a SCORED cipro-QRDR cell "
        f"({sorted(cipro_orgs)}) -> the same target-site catalog mechanism transfers across genera")


def test_g9r2_interventional_kb_falsified():
    """g9r2-interventional-kb (the catalog wins because each entry is a background-INDEPENDENT
    interventional edit->effect statement, so the world model should BE a catalog-shaped interventional
    KB) DEAD iff catalogued DRMs are actually background-DEPENDENT (observational-like): the catalog
    does NOT generalize across held-out studies -> leave-study-out balanced accuracy < 0.70. Committed:
    0.824 -> generalizes across studies (interventional signature) -> SURVIVES."""
    d = json.loads((REPO / "wiki/hiv_catalog_accessory_extension_2026-07-12.json").read_text(encoding="utf-8"))
    bal_acc = float(d["pooled"]["catalog"]["bal_acc"])
    assert bal_acc < 0.70, (
        f"NOT falsified: catalog leave-study-out balanced accuracy = {bal_acc:.3f} >= 0.70 -> catalogued "
        "DRMs generalize across independent studies (the background-independent interventional signature)")


def test_g5r2_self_supervised_catalog_falsified():
    """g5r2-self-supervised-catalog (relax 'a learned world model needs external phenotype labels' ->
    train a generative effect-predictor SELF-SUPERVISED on the catalog's own edit->effect entries) DEAD
    iff the pooled cross-cell catalog is too small to be a training set: < 100 total catalogued
    edit->effect statements. Committed: 16 HIV NNRTI majors + 438 TB WHO grade-1/2 = 454 (before the
    other viral/bacterial cells) -> SURVIVES."""
    from dna_decode.data import hiv_amr as H
    n_hiv = len(H.NNRTI_RT_MAJOR_DRMS)
    n_tb = 438   # WHO catalogue v2 grade-1/2 (committed CHECKSUMS join count, CLAUDE.md)
    total = n_hiv + n_tb
    assert total < 100, (
        f"NOT falsified: {total} catalogued edit->effect statements exist (>=100) -> the catalog is a "
        "usable self-supervised training set for a generative effect-predictor")


def test_control_genewide_interpolation_falsified():
    """CONTROL (tempting-but-wrong): 'gene-wide interpolation — ANY mutation ANYWHERE in the gene
    predicts effect'. DEAD iff the signal is LOCALIZED to known-functional sites, not gene-wide:
    known-functional enrichment R-over-S > 1.5. Committed: 4.203 (localized) -> this test PASSES ->
    the gene-wide claim is correctly KILLED (round-2 discrimination proof)."""
    d = json.loads((REPO / "wiki/hiv_blindspot_pocket_localization_2026-07-09.json").read_text(encoding="utf-8"))
    enrichment = float(d["enrichment_R_over_S_known_functional"])
    assert enrichment > 1.5, (
        f"NOT falsified: known-functional enrichment R-over-S = {enrichment:.2f} <= 1.5 -> signal is NOT "
        "localized, so gene-wide interpolation would be viable (committed data would have to show this)")
