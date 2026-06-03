"""v0 deterministic virulence-marker presence screen (EP-4) — confound-immune.

The learned k-mer/NT track was confounded (study==class) and its top features were
rare CTAG-motif k-mers = assembler batch, not biology (see
pathotype_model_interpret_confound_2026-06-02.json). This is the ledger-LOCKED v0
alternative: search the 24 cached genomes for SPECIFIC known virulence marker genes
from the VirulenceFinder E. coli DB. Whole-gene presence (~1-3 kb) is robust to the
rare-k-mer assembly artifact, and looking for KNOWN biology sidesteps the
"learn whatever separates the two studies" failure mode.

Detection: pure-Python k-mer seeding (k=15, no BLAST). For each marker allele, the
fraction of its k=15 seeds (forward OR reverse-complement) found in the genome's
forward 15-mer set; a gene FAMILY is called present at the max coverage over its
alleles, present if >= COV_THRESH. Genomes are the same 24 cached ENA assemblies as
every other arm (load_strains), so results are directly comparable.

CONFOUND-IMMUNE TEST: ExPEC=Salipante (y=0), EPEC=Hazen (y=1). EPEC is DEFINED by LEE
(eae) +/- bfpA; ExPEC/UPEC by papC/sfa/afa + iron (iutA/fyuA) + toxins (hlyA/cnf1) +
capsule (kpsMII). If those known markers track the labels, the contrast is real
pathotype biology, not assembler batch.
"""
from __future__ import annotations
import json, re, sys
from collections import defaultdict
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
import numpy as np
from sklearn.metrics import roc_auc_score

from scripts.pathotype_kmer_bakeoff import load_strains

DB = REPO / "data/virulencefinder_db/virulence_ecoli.fsa"
K = 15
COV_THRESH = 0.80
_COMP = str.maketrans("ACGTacgt", "TGCAtgca")


def revcomp(s: str) -> str:
    return s.translate(_COMP)[::-1]


# allele-name -> (gene_family, pathotype_group). Matched by startswith (lowercased).
# group keys: LEE (EPEC/EHEC shared), bfp (tEPEC), stx (EHEC), LT/ST (ETEC),
# EAEC, ExPEC_adhesin, ExPEC_iron, ExPEC_toxin, ExPEC_capsule, ExPEC_other.
FAMILY_RULES = [
    ("eae",  "eae",   "LEE"), ("tir", "tir", "LEE"), ("espa", "espA", "LEE"),
    ("espb", "espB", "LEE"), ("espf", "espF", "LEE"),
    ("bfpa", "bfpA",  "bfp"),
    ("stx1", "stx1",  "stx"), ("stx2", "stx2", "stx"),
    ("eltiab","LT",   "ETEC"), ("eltiiab","LT-II","ETEC"), ("elti-","LT","ETEC"),
    ("estap","STa",   "ETEC"), ("estah","STa","ETEC"), ("estb","STb","ETEC"),
    ("aggr","aggR",   "EAEC"), ("aata","aatA","EAEC"), ("aaic","aaiC","EAEC"),
    ("aap","aap",     "EAEC"), ("agg","AAF","EAEC"),
    ("papc","papC",   "ExPEC_adhesin"), ("papa","papA","ExPEC_adhesin"),
    ("sfa","sfa",     "ExPEC_adhesin"), ("foc","foc","ExPEC_adhesin"),
    ("afa","afa",     "ExPEC_adhesin"), ("dra","dra","ExPEC_adhesin"),
    ("iha","iha",     "ExPEC_adhesin"),
    ("iuta","iutA",   "ExPEC_iron"), ("fyua","fyuA","ExPEC_iron"),
    ("irp2","irp2",   "ExPEC_iron"), ("irea","ireA","ExPEC_iron"),
    ("chua","chuA",   "ExPEC_iron"), ("sita","sitA","ExPEC_iron"),
    ("hlya","hlyA",   "ExPEC_toxin"), ("hlyf","hlyF","ExPEC_toxin"),
    ("cnf","cnf1",    "ExPEC_toxin"), ("sat","sat","ExPEC_toxin"),
    ("vat","vat",     "ExPEC_toxin"), ("usp","usp","ExPEC_toxin"),
    ("kpsmii","kpsMII","ExPEC_capsule"), ("kpsm","kpsM","ExPEC_capsule"),
    ("iss","iss",     "ExPEC_other"), ("ompt","ompT","ExPEC_other"),
    ("tsh","tsh",     "ExPEC_other"),
]


def family_of(allele_name: str):
    a = allele_name.lower()
    for prefix, fam, grp in FAMILY_RULES:
        if a.startswith(prefix):
            return fam, grp
    return None, None


def load_db():
    """gene_family -> list[allele_seq]; family -> group."""
    fam_seqs = defaultdict(list)
    fam_group = {}
    name = None
    buf = []
    for line in DB.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith(">"):
            if name is not None:
                fam, grp = family_of(name)
                if fam:
                    fam_seqs[fam].append("".join(buf))
                    fam_group[fam] = grp
            name = line[1:].split(":")[0]
            buf = []
        else:
            buf.append(line.strip())
    if name is not None:
        fam, grp = family_of(name)
        if fam:
            fam_seqs[fam].append("".join(buf)); fam_group[fam] = grp
    return fam_seqs, fam_group


def genome_kmer_set(seq: str, k: int) -> set:
    s = seq.upper()
    out = set()
    for i in range(len(s) - k + 1):
        km = s[i:i + k]
        if "N" not in km:
            out.add(km)
    return out


def allele_coverage(allele: str, gset: set, k: int) -> float:
    a = allele.upper()
    if len(a) < k:
        return 0.0
    def cov(seq):
        kms = [seq[i:i+k] for i in range(len(seq) - k + 1)]
        kms = [x for x in kms if "N" not in x]
        if not kms:
            return 0.0
        return sum(1 for x in kms if x in gset) / len(kms)
    return max(cov(a), cov(revcomp(a)))


def main() -> int:
    seqs, labels, ids = load_strains()
    fam_seqs, fam_group = load_db()
    families = sorted(fam_seqs, key=lambda f: (fam_group[f], f))
    print(f"[v0] {len(ids)} genomes; {len(families)} marker families "
          f"({sum(len(v) for v in fam_seqs.values())} alleles)")

    # per-genome best coverage per family
    cov = {sid: {} for sid in ids}
    for sid in ids:
        gset = genome_kmer_set(seqs[sid], K)
        for fam in families:
            best = max(allele_coverage(al, gset, K) for al in fam_seqs[fam])
            cov[sid][fam] = round(best, 3)
        present = [f for f in families if cov[sid][f] >= COV_THRESH]
        print(f"[v0] {sid.split('|')[0]} (y={labels[sid]}): {len(present)} markers >= {COV_THRESH}", flush=True)

    def present(sid, fam):
        return cov[sid].get(fam, 0) >= COV_THRESH

    # --- confound-immune biology test: do known markers track the labels? ---
    y = np.array([labels[s] for s in ids])  # 0=ExPEC(Salipante), 1=EPEC(Hazen)
    def group_count(sid, grp):
        return sum(1 for f in families if fam_group[f] == grp and present(sid, f))

    rows = []
    for sid in ids:
        rows.append({
            "id": sid.split("|")[0], "y": labels[sid],
            "eae": present(sid, "eae"), "bfpA": present(sid, "bfpA"),
            "stx": group_count(sid, "stx") > 0,
            "LEE_genes": group_count(sid, "LEE"),
            "ExPEC_adhesin": group_count(sid, "ExPEC_adhesin"),
            "ExPEC_iron": group_count(sid, "ExPEC_iron"),
            "ExPEC_toxin": group_count(sid, "ExPEC_toxin"),
            "ExPEC_capsule": group_count(sid, "ExPEC_capsule"),
        })

    # eae presence as an EPEC predictor; ExPEC virulence-score as an ExPEC predictor
    eae_vec = np.array([1 if r["eae"] else 0 for r in rows])
    expec_score = np.array([r["ExPEC_adhesin"] + r["ExPEC_iron"] + r["ExPEC_toxin"] + r["ExPEC_capsule"]
                            for r in rows], dtype=float)
    # EPEC (y=1) should have eae; ExPEC (y=0) should have high expec_score
    eae_in_epec = int(eae_vec[y == 1].sum()); eae_in_expec = int(eae_vec[y == 0].sum())
    auroc_eae = float(roc_auc_score(y, eae_vec)) if len(set(y)) == 2 else float("nan")
    auroc_expec = float(roc_auc_score(1 - y, expec_score)) if len(set(y)) == 2 else float("nan")  # predict ExPEC(y=0)

    print("\n[v0] === confound-immune biology test ===")
    print(f"[v0] eae present in EPEC(Hazen): {eae_in_epec}/12 | in ExPEC(Salipante): {eae_in_expec}/12")
    print(f"[v0] eae-presence -> EPEC AUROC = {auroc_eae:.3f}")
    print(f"[v0] ExPEC-virulence-score -> ExPEC AUROC = {auroc_expec:.3f}")
    print(f"[v0] mean ExPEC-score: ExPEC(Salipante)={expec_score[y==0].mean():.2f} EPEC(Hazen)={expec_score[y==1].mean():.2f}")

    verdict = ("BIOLOGY_TRACKS_LABELS: known markers separate the classes -> the ExPEC/EPEC labels reflect "
               "real pathotype gene content (NOT just assembler batch); v0 marker resolver is the right path."
               if (auroc_eae >= 0.8 or auroc_expec >= 0.8) else
               "MARKERS_DO_NOT_TRACK: known virulence markers do NOT cleanly separate the labels -> either the "
               "labels are weak, the marker set/threshold needs work, or signal really is non-canonical. Investigate before trusting.")
    print(f"[v0] VERDICT: {verdict}")

    res = {"contrast": "ExPEC(Salipante) vs EPEC(Hazen)", "n": len(ids), "k": K,
           "cov_threshold": COV_THRESH, "n_marker_families": len(families),
           "marker_families": {f: fam_group[f] for f in families},
           "eae_in_epec_of12": eae_in_epec, "eae_in_expec_of12": eae_in_expec,
           "auroc_eae_predicts_epec": auroc_eae,
           "auroc_expec_score_predicts_expec": auroc_expec,
           "mean_expec_score_expec": float(expec_score[y == 0].mean()),
           "mean_expec_score_epec": float(expec_score[y == 1].mean()),
           "verdict": verdict, "per_genome": rows,
           "coverage_matrix": {sid.split("|")[0]: cov[sid] for sid in ids}}
    out = REPO / "research_outputs/pathotype_v0_marker_screen_2026-06-02.json"
    out.write_text(json.dumps(res, indent=2), encoding="utf-8")
    print(f"[v0] wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
