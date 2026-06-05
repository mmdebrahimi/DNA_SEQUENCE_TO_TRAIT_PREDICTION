"""QRDR-POINT knowledge baseline — the 'best classical' comparator for fluoroquinolone AMR.

For ciprofloxacin, resistance is conferred by SPECIFIC point mutations in the QRDR (gyrA S83L/D87N,
parC S80I/E84V, parE, gyrB) plus acquired plasmid-mediated quinolone genes (qnr*, aac(6')-Ib-cr).
Gene PRESENCE is useless (R and S both carry gyrA); the discriminative feature is the ALLELE/position.
So the knowledge baseline = binary presence of each distinct NON-SYNONYMOUS QRDR point mutation +
each plasmid-quinolone gene, read from the per-strain AMRFinderPlus mutations.tsv / main.tsv.

This is the comparator NT must beat to claim "beats best classical" (per the 2× brainstorm + ledger
Phase2_Decision_Gate). Built from the committed AMRFinder cache (data/amrfinder_runs/<accession>/);
no Docker at feature-build time. Reuses mic_tiers loci catalogs + a synonymous-mutation filter.
"""
from __future__ import annotations

import csv
from pathlib import Path

import numpy as np

from dna_decode.data.mic_tiers import loci_by_mechanism_for


def _is_synonymous(symbol: str) -> bool:
    """True for a synonymous POINT mutation like '16S_A523A' / 'gyrA_G141G' (wt AA == alt AA).
    Format '<gene>_<wt><pos><alt>'; non-parseable → False (treat as non-synonymous / keep)."""
    if "_" not in symbol:
        return False
    tail = symbol.split("_", 1)[1]
    # tail like 'S83L' (wt='S', alt='L') or 'A523A'. Need leading + trailing alpha around digits.
    lead = tail[:1]
    trail = tail[-1:]
    if not (lead.isalpha() and trail.isalpha()):
        return False
    has_digit = any(c.isdigit() for c in tail)
    return has_digit and lead == trail


def _locus_of(symbol: str) -> str:
    """Gene locus prefix before the first '_' (e.g. 'gyrA_S83L' -> 'gyrA')."""
    return symbol.split("_", 1)[0]


def strain_point_features(runs_root: Path, accession: str, drug: str) -> set[str] | None:
    """Set of drug-relevant feature tokens for one strain, or None if its AMRFinder cache is absent.

    Tokens: each distinct NON-SYNONYMOUS QRDR point mutation symbol (from mutations.tsv) + each
    acquired plasmid-quinolone gene symbol (from main.tsv). None = strain not yet audited (skip)."""
    by_mech = loci_by_mechanism_for(drug)
    target_loci = set(by_mech.get("QRDR_target_alteration", set()))
    plasmid_loci = set(by_mech.get("plasmid_protect_modify", set()))
    d = Path(runs_root) / accession
    mut = d / "mutations.tsv"
    main = d / "main.tsv"
    if not mut.exists() and not main.exists():
        return None
    feats: set[str] = set()
    # QRDR point mutations (mutations.tsv): keep non-synonymous on a target locus
    if mut.exists():
        with open(mut, encoding="utf-8") as f:
            for row in csv.DictReader(f, delimiter="\t"):
                sym = (row.get("Element symbol") or "").strip()
                if not sym or _is_synonymous(sym):
                    continue
                if _locus_of(sym) in target_loci:
                    feats.add(sym)              # specific allele, e.g. gyrA_S83L
    # acquired plasmid-quinolone genes (main.tsv): presence by locus
    if main.exists():
        with open(main, encoding="utf-8") as f:
            for row in csv.DictReader(f, delimiter="\t"):
                sym = (row.get("Element symbol") or "").strip()
                loc = _locus_of(sym)
                if loc in plasmid_loci or sym in plasmid_loci:
                    feats.add(f"plasmid:{loc}")
    return feats


def build_point_matrix(runs_root, accessions: list[str], drug: str):
    """Return (X, feature_names, present_mask) aligned to `accessions`.
    X[i] = binary presence vector over the cohort-union feature vocab. present_mask[i]=False if the
    strain has no AMRFinder cache (row is all-zero + flagged so the caller can drop it)."""
    per = [strain_point_features(runs_root, a, drug) for a in accessions]
    vocab = sorted({t for s in per if s for t in s})
    idx = {t: i for i, t in enumerate(vocab)}
    X = np.zeros((len(accessions), len(vocab)), dtype=np.float32)
    present = np.zeros(len(accessions), dtype=bool)
    for i, s in enumerate(per):
        if s is None:
            continue
        present[i] = True
        for t in s:
            X[i, idx[t]] = 1.0
    return X, vocab, present
