"""Deployable learned-decoder clinical AUROC — place AlphaMissense AND the shipped ESM2+ProSST hybrid between
the BLOSUM62 FLOOR and the DMS-functional CEILING on the same actionable-human-gene ClinVar joins.

Follow-up to `scripts/clinical_variant_effect_validate.py` (the floor+ceiling cell). That cell established,
on the MaveDB-DMS ⋈ ClinVar path/benign joins: DMS-itself AUROC (fitness-alignment CEILING) TP53 0.996 /
MSH2 0.955, BLOSUM62 FLOOR 0.707 / 0.832 — a large gap = the headroom a LEARNED decoder should capture. This
places the deployable learned decoders in that gap, ON THE SAME joined variant set (so every number is
directly comparable):

  - **AlphaMissense** (Cheng 2023): precomputed pathogenicity in [0,1] for every human missense — free, NO
    GPU, NO network at score time. Higher = more pathogenic. The deployable, no-model-run winner.
  - **ESM2-650M** masked-marginals (the sequence model) — higher logP-delta = more preserved.
  - **ProSST-2048** structure-conditioned (AlphaFold v6 structure tokens) — higher = more preserved.
  - **ESM2+ProSST hybrid** = the SHIPPED `predict_effect(method='hybrid')` decoder (rank-average of the two
    orthogonal modalities; `wiki/forward_modality_hybrid_2026-07-17.md`). This validates the ACTUAL deployed
    forward cell on clinical labels, not just a proxy.

Coordinate-integrity gate: variants whose WT residue does not match the UniProt canonical sequence are
dropped before ESM2/ProSST scoring (reported), never mis-scored.

HONEST RAILS (inherited): in-distribution-clinical NOT held-out (AM/ESM2/ProSST saw these proteins in
training; but the ClinVar labels are independent of the DMS-fitness tuning); single-class genes stay
AUROC-inapplicable. AlphaMissense is CC BY-NC-SA (non-commercial) — research use.

  uv run python scripts/clinical_am_hybrid_auroc.py                 # AM only (fast)
  HF_HOME=D:/hf_cache uv run python scripts/clinical_am_hybrid_auroc.py --hybrid   # + ESM2 + ProSST + hybrid
  uv run python scripts/clinical_am_hybrid_auroc.py --build-am-filter              # rebuild the per-gene AM filter

Frozen AMR surface byte-unchanged (READ-only).
"""
from __future__ import annotations

import argparse
import gzip
import json
import sys
import urllib.request
from datetime import date as _date
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from dna_decode.forward.variant_effect import blosum62_score, rank_average_hybrid  # noqa: E402
from scripts.clinical_variant_effect_validate import (  # noqa: E402
    auroc, _spearman_sign, fetch_mavedb_dms, fetch_clinvar_missense, CLINICAL_GENES, MIN_PER_CLASS,
)

AM_GZ = Path("D:/dna_decode_cache/alphamissense/AlphaMissense_aa_substitutions.tsv.gz")
AM_FILTERED = Path("D:/dna_decode_cache/alphamissense/am_clinical_filtered.tsv")
AF_DIR = Path("D:/dna_decode_cache/alphafold")
ESM_DIR = Path("D:/dna_decode_cache/esm")
# UniProt accessions for the AUROC-viable clinical genes (numbering matches the MaveDB DMS + ClinVar, offset 0).
GENE_UNIPROT = {"TP53": "P04637", "MSH2": "P43246"}


# ------------------------------------------------------------------ AlphaMissense -------------------------

def build_am_filter(uniprots: set[str], am_gz: Path = AM_GZ, out: Path = AM_FILTERED) -> int:
    """Stream the ~5.6 GB AlphaMissense gz once and write ONLY the rows for `uniprots` to a small plain TSV."""
    n = 0
    out.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(am_gz, "rt") as fin, out.open("w", encoding="utf-8") as fout:
        for ln in fin:
            if ln.startswith("#"):
                continue
            if ln.split("\t", 1)[0] in uniprots:
                fout.write(ln)
                n += 1
    return n


def load_am(uniprot: str, filtered: Path = AM_FILTERED) -> dict[tuple[str, int, str], float]:
    """{(wt,pos,alt): am_pathogenicity} for one UniProt from the filtered TSV. Higher = more pathogenic."""
    out: dict[tuple[str, int, str], float] = {}
    if not filtered.exists():
        return out
    for ln in filtered.open(encoding="utf-8"):
        p = ln.rstrip("\n").split("\t")
        if len(p) < 3 or p[0] != uniprot:
            continue
        var = p[1]
        if len(var) < 3 or not var[1:-1].isdigit():
            continue
        try:
            out[(var[0], int(var[1:-1]), var[-1])] = float(p[2])
        except ValueError:
            continue
    return out


# ------------------------------------------------------------------ sequence + structure + models --------

_SEQ_CACHE: dict[str, str] = {}


def uniprot_seq(uniprot: str) -> str | None:
    if uniprot in _SEQ_CACHE:
        return _SEQ_CACHE[uniprot]
    try:
        with urllib.request.urlopen(f"https://rest.uniprot.org/uniprotkb/{uniprot}.fasta", timeout=60) as r:
            fa = r.read().decode()
        seq = "".join(l.strip() for l in fa.splitlines() if not l.startswith(">"))
        _SEQ_CACHE[uniprot] = seq
        return seq
    except Exception:  # noqa: BLE001
        return None


def af_pdb(uniprot: str) -> Path | None:
    """Fetch the CURRENT AlphaFold PDB (v6; the shipped helper hardcodes v4 which now 404s). Cached to D:."""
    AF_DIR.mkdir(parents=True, exist_ok=True)
    dest = AF_DIR / f"AF-{uniprot}-F1-v6.pdb"
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    for v in ("v6", "v5", "v4"):
        try:
            urllib.request.urlretrieve(
                f"https://alphafold.ebi.ac.uk/files/AF-{uniprot}-F1-model_{v}.pdb", dest)
            if dest.stat().st_size > 0:
                return dest
        except Exception:  # noqa: BLE001
            continue
    return None


def struct_tokens(uniprot: str) -> list[int] | None:
    cache = AF_DIR / f"{uniprot}_prosst_structok.json"
    if cache.exists():
        return json.loads(cache.read_text())
    pdb = af_pdb(uniprot)
    if pdb is None:
        return None
    try:
        from dna_decode.forward.prosst_scorer import quantize_structure
        toks = quantize_structure(pdb)
    except Exception:  # noqa: BLE001
        return None
    cache.write_text(json.dumps(toks))
    return toks


def esm2_var_table(uniprot: str, seq: str, positions: list[int]) -> dict[tuple[str, int, str], float] | None:
    """ESM2-650M masked-marginal delta {(wt,pos,alt): logP(alt)-logP(wt)} at ONLY `positions` (fast). Cached."""
    cache = ESM_DIR / f"clinical_{uniprot}_esm2_650m.json"
    have: dict[str, float] = {}
    if cache.exists():
        have = json.loads(cache.read_text())
    need = [p for p in positions if f"@{p}" not in have.get("_positions_done", [])]
    if need:
        try:
            from dna_decode.forward.esm_scorer import esm2_logp_table
            from dna_decode.forward.variant_effect import esm_pos_table_to_variant_table
        except Exception:  # noqa: BLE001
            return None
        ESM_DIR.mkdir(parents=True, exist_ok=True)
        # CHUNK + persist per chunk so a mid-run kill (documented D:-USB-hiccup / memory failure mode) only
        # loses the current chunk — a rerun resumes from `_positions_done`. Restartable, per the project lesson.
        CHUNK = 24
        for i in range(0, len(need), CHUNK):
            sub = need[i:i + CHUNK]
            try:
                pos_table = esm2_logp_table(seq, positions=sub)
            except Exception:  # noqa: BLE001
                break  # keep whatever chunks already persisted; caller handles a short table
            for k, v in esm_pos_table_to_variant_table(pos_table, seq).items():
                have[k] = v
            have["_positions_done"] = sorted(set(have.get("_positions_done", [])) | {f"@{p}" for p in sub})
            cache.write_text(json.dumps(have))
    out: dict[tuple[str, int, str], float] = {}
    for k, v in have.items():
        if k.startswith("_") or len(k) < 3 or not k[1:-1].isdigit():
            continue
        out[(k[0], int(k[1:-1]), k[-1])] = v
    return out


def prosst_var_table(seq: str, toks: list[int], variants: list[tuple[str, int, str]]
                     ) -> dict[tuple[str, int, str], float] | None:
    from dna_decode.forward.prosst_scorer import prosst_variant_table
    muts = [f"{w}{p}{a}" for (w, p, a) in variants]
    try:
        tbl = prosst_variant_table(seq, muts, structure_tokens=toks)
    except Exception:  # noqa: BLE001
        return None
    out: dict[tuple[str, int, str], float] = {}
    for k, v in tbl.items():
        if len(k) >= 3 and k[1:-1].isdigit():
            out[(k[0], int(k[1:-1]), k[-1])] = v
    return out


# ------------------------------------------------------------------ per-gene scoring ---------------------

def score_gene(gene: str, meta: dict, *, with_hybrid: bool = False) -> dict:
    uniprot = GENE_UNIPROT[gene]
    dms = fetch_mavedb_dms(meta["urn"])
    clin = fetch_clinvar_missense(gene, use_cache=True)
    am = load_am(uniprot)
    rec = {"gene": gene, "uniprot": uniprot, "urn": meta["urn"],
           "n_dms": len(dms), "n_clinvar": len(clin), "n_am": len(am)}

    # AM standalone on the full ClinVar set (AM's real deployment — no DMS needed).
    clin_am = sorted(set(clin) & set(am))
    p_full = sum(1 for k in clin_am if clin[k] == "PATH")
    rec["am_full_clinvar"] = {"n": len(clin_am), "n_path": p_full, "n_benign": len(clin_am) - p_full}
    if p_full >= MIN_PER_CLASS and (len(clin_am) - p_full) >= MIN_PER_CLASS:
        labels = [clin[k] == "PATH" for k in clin_am]
        rec["am_full_clinvar"]["am_auroc"] = round(auroc(labels, [am[k] for k in clin_am]), 4)

    # Comparable set: MaveDB ∩ ClinVar ∩ AM — floor/ceiling/AM on identical variants.
    shared = sorted(set(dms) & set(clin) & set(am))
    n_path = sum(1 for k in shared if clin[k] == "PATH")
    n_benign = len(shared) - n_path
    comp = {"n_joined": len(shared), "n_path": n_path, "n_benign": n_benign,
            "auroc_applicable": n_path >= MIN_PER_CLASS and n_benign >= MIN_PER_CLASS}
    if comp["auroc_applicable"]:
        labels = [clin[k] == "PATH" for k in shared]
        blos = [blosum62_score(k[0], k[2]) for k in shared]
        dvals = [dms[k] for k in shared]
        dsign = _spearman_sign(dvals, blos) or 1.0
        comp["dms_ceiling_auroc"] = round(auroc(labels, [-(dsign * v) for v in dvals]), 4)
        comp["blosum_floor_auroc"] = round(auroc(labels, [-v for v in blos]), 4)
        comp["am_auroc"] = round(auroc(labels, [am[k] for k in shared]), 4)

        if with_hybrid:
            seq = uniprot_seq(uniprot)
            toks = struct_tokens(uniprot) if seq else None
            if not seq or not toks:
                comp["hybrid_note"] = "ESM2/ProSST unavailable (seq/structure) — see Kaggle follow-up"
            else:
                # coordinate-integrity gate: keep only WT-consistent variants
                ok = [k for k in shared if 1 <= k[1] <= len(seq) and seq[k[1] - 1] == k[0]]
                comp["n_wt_consistent"] = len(ok)
                positions = sorted({k[1] for k in ok})
                esm = esm2_var_table(uniprot, seq, positions)
                pros = prosst_var_table(seq, toks, ok)
                if esm and pros:
                    both = [k for k in ok if k in esm and k in pros]
                    lab2 = [clin[k] == "PATH" for k in both]
                    if sum(lab2) >= MIN_PER_CLASS and (len(lab2) - sum(lab2)) >= MIN_PER_CLASS:
                        esm_s = {f"{k[0]}{k[1]}{k[2]}": esm[k] for k in both}
                        pro_s = {f"{k[0]}{k[1]}{k[2]}": pros[k] for k in both}
                        hyb = rank_average_hybrid([esm_s, pro_s])  # higher = preserved
                        comp["n_hybrid"] = len(both)
                        comp["esm2_auroc"] = round(auroc(lab2, [-esm[k] for k in both]), 4)
                        comp["prosst_auroc"] = round(auroc(lab2, [-pros[k] for k in both]), 4)
                        comp["hybrid_auroc"] = round(
                            auroc(lab2, [-hyb[f"{k[0]}{k[1]}{k[2]}"] for k in both]), 4)
                    else:
                        comp["hybrid_note"] = f"post-WT-gate class balance too small (n={len(both)})"
                else:
                    comp["hybrid_note"] = "ESM2 or ProSST scoring returned empty"
    rec["comparable_set"] = comp
    return rec


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--build-am-filter", action="store_true", help="(re)build the per-gene AM filter from the gz")
    ap.add_argument("--hybrid", action="store_true", help="also run ESM2-650M + ProSST + the hybrid (slow, CPU)")
    a = ap.parse_args()

    if a.build_am_filter or not AM_FILTERED.exists():
        print(f"building AM filter for {sorted(set(GENE_UNIPROT.values()))} ...", flush=True)
        print(f"  wrote {build_am_filter(set(GENE_UNIPROT.values()))} AM rows -> {AM_FILTERED}")

    results = []
    for gene in GENE_UNIPROT:
        print(f"[{gene}] scoring{' AM+ESM2+ProSST+hybrid' if a.hybrid else ' AM'} ...", flush=True)
        rec = score_gene(gene, CLINICAL_GENES[gene], with_hybrid=a.hybrid)
        results.append(rec)
        c = rec["comparable_set"]
        if c.get("auroc_applicable"):
            print(f"  comparable n={c['n_joined']} ({c['n_path']}P/{c['n_benign']}B): "
                  f"FLOOR blosum={c['blosum_floor_auroc']}  AM={c['am_auroc']}  "
                  f"CEILING dms={c['dms_ceiling_auroc']}")
            if "hybrid_auroc" in c:
                print(f"    LEARNED (n_hybrid={c['n_hybrid']}): ESM2={c['esm2_auroc']}  "
                      f"ProSST={c['prosst_auroc']}  HYBRID={c['hybrid_auroc']}")
            elif "hybrid_note" in c:
                print(f"    hybrid: {c['hybrid_note']}")
        f = rec["am_full_clinvar"]
        if "am_auroc" in f:
            print(f"  AM standalone on FULL ClinVar (n={f['n']}, {f['n_path']}P/{f['n_benign']}B): {f['am_auroc']}")

    art = {"_schema": "clinical-am-hybrid-auroc-v1", "date": _date.today().isoformat(),
           "question": "Place the deployable LEARNED decoder (AlphaMissense; ESM2+ProSST hybrid) between the "
                       "BLOSUM floor and the DMS ceiling on the actionable-human-gene ClinVar joins.",
           "am_source": "AlphaMissense_aa_substitutions (Cheng 2023; CC BY-NC-SA; precomputed, no GPU) via D: cache",
           "hybrid_source": "shipped predict_effect(method='hybrid') = rank-avg(ESM2-650M, ProSST-2048 on AlphaFold v6)",
           "tier": "in_distribution_clinical (proteins seen in AM/ESM2/ProSST training; ClinVar labels "
                   "independent of the DMS-fitness tuning)",
           "hybrid_ran": a.hybrid, "results": results, "frozen_surface_changed": False}
    out = Path(f"wiki/clinical_am_hybrid_auroc_{_date.today().isoformat()}.json")
    out.write_text(json.dumps(art, indent=2), encoding="utf-8")
    print(f"\nartifact: {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
