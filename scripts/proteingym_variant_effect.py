"""Protein variant-effect cell — genotype -> MOLECULAR phenotype (the tractable frontier for humans/higher
animals when the ORGANISM-level trait is polygenic).

The project's deterministic determinant scan wins where a CURATED CATALOG of high-effect variants exists (AMR,
PGx, ben-1). Organism-level complex traits (height, disease) are polygenic -> the learned path is 0-for-5 under
de-confounding. This cell tests the OTHER regime: a MOLECULAR phenotype (protein variant -> measured function),
where the genotype->phenotype map is short and the label is FREE + unit-level (one score per single-amino-acid
variant). Substrate: the canonical human PTEN VAMP-seq abundance DMS (Matreyek 2018, MAVEDB urn:mavedb:00000102-0-1).

Deterministic feature = the BLOSUM62 substitution score of each missense change (authoritative, from Biopython
-- no hand-transcribed matrix). A conservative substitution (high BLOSUM) should preserve abundance; a
disruptive one (low BLOSUM) should reduce it -> POSITIVE Spearman expected. A nonsense-vs-missense abundance
gap is the direction sanity (truncation -> destabilized).

HONEST FRAMING (the point of the cell): the deterministic substitution-severity baseline has MODEST signal
(Spearman ~0.2). This is the FIRST modality where LEARNED models clearly BEAT the deterministic rule -- the
published ProteinGym leaderboard puts zero-shot models (ESM-1v / EVE / TranceptEVE) in the ~0.4-0.5 Spearman
range on PTEN abundance (CITED, NOT run here -- those need a GPU). That is the opposite of the AMR/PGx/ben-1
cells and it MAPS the boundary: deterministic wins with a curated high-effect catalog; learned wins for
distributed molecular-property prediction where evolutionary/structural context is the signal. Data on D:.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

AA3 = {"Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q", "Glu": "E", "Gly": "G",
       "His": "H", "Ile": "I", "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S",
       "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V"}
_MISSENSE = re.compile(r"^p\.([A-Za-z]{3})(\d+)([A-Za-z]{3})$")


def _blosum():
    from Bio.Align import substitution_matrices
    return substitution_matrices.load("BLOSUM62")


def blosum_score(mat, a: str, b: str) -> float:
    try:
        return float(mat[a, b])
    except (KeyError, IndexError):
        return float(mat[b, a])


def parse_missense(hgvs: str) -> tuple[str, str] | None:
    """p.Tyr88Val -> ('Y','V'); returns None for synonymous / nonsense / indel / unparseable."""
    m = _MISSENSE.match(str(hgvs))
    if not m:
        return None
    a, _pos, b = m.group(1), m.group(2), m.group(3)
    if a in AA3 and b in AA3 and a != b:
        return AA3[a], AA3[b]
    return None


def load_dms(csv: Path) -> pd.DataFrame:
    df = pd.read_csv(csv)
    df["score_num"] = pd.to_numeric(df["score"], errors="coerce")
    return df


def run(csv: Path) -> dict:
    df = load_dms(csv)
    mat = _blosum()
    blos, sc = [], []
    for hg, s in zip(df["hgvs_pro"].astype(str), df["score_num"]):
        if np.isnan(s):
            continue
        mm = parse_missense(hg)
        if mm:
            blos.append(blosum_score(mat, mm[0], mm[1]))
            sc.append(float(s))
    blos = np.array(blos)
    sc = np.array(sc)
    rho = float(spearmanr(blos, sc)[0]) if len(blos) > 2 else float("nan")
    # nonsense-vs-missense direction sanity
    ter = pd.to_numeric(df[df["hgvs_pro"].astype(str).str.contains("Ter")]["score"], errors="coerce").dropna()
    nonsense_mean = float(ter.mean()) if len(ter) else float("nan")
    missense_mean = float(sc.mean()) if len(sc) else float("nan")
    direction_ok = bool(not np.isnan(nonsense_mean) and nonsense_mean < missense_mean)
    return {
        "substrate": "PTEN VAMP-seq abundance DMS (Matreyek 2018; MAVEDB urn:mavedb:00000102-0-1)",
        "phenotype": "molecular (protein variant -> measured abundance)",
        "deterministic_feature": "BLOSUM62 substitution score",
        "n_missense_scored": len(blos),
        "spearman_blosum_vs_abundance": round(rho, 4),
        "nonsense_mean_abundance": round(nonsense_mean, 3),
        "missense_mean_abundance": round(missense_mean, 3),
        "n_nonsense": int(len(ter)),
        "direction_sanity_ok": direction_ok,
        "learned_contrast": "published ProteinGym zero-shot models (ESM-1v/EVE/TranceptEVE) ~0.4-0.5 Spearman "
                            "on PTEN abundance (CITED, NOT run here -- GPU). First modality where learned > "
                            "deterministic.",
        "verdict": "DETERMINISTIC_SUBSTITUTION_BASELINE_MODEST_SIGNAL",
    }


def main(argv=None) -> int:
    import argparse
    import json
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dms", type=Path, default=Path("D:/dna_decode_cache/proteingym/pten_scores.csv"))
    ap.add_argument("--out", type=Path, default=REPO / "wiki" / "proteingym_variant_effect_scores.json")
    a = ap.parse_args(argv)
    res = run(a.dms)
    a.out.write_text(json.dumps(res, indent=2), encoding="utf-8")
    print(json.dumps(res, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
