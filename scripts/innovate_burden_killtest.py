"""Kill-test helper for the /innovate 'determinant-burden rescue' candidate. Exits 0 (DISPROOF FOUND =
KILLED) iff determinant burden does NOT separate R from S for <drug> in the gono NCBI-PD cohort (R-mean
minus S-mean burden < --min-gap); exits 1 (SURVIVED) iff it does separate."""
import argparse, csv, statistics as st, sys
ap = argparse.ArgumentParser()
ap.add_argument("--drug", required=True)
ap.add_argument("--min-gap", type=float, default=2.0)
ap.add_argument("--cohort", default="data/raw/gono_ncbipd_extval")
a = ap.parse_args()
lab = {r["biosample"]: r for r in csv.DictReader(open(f"{a.cohort}/cohort.tsv"), delimiter="\t")}
det = {r["biosample"]: [s for s in (r["determinants"] or "").split(";") if s]
       for r in csv.DictReader(open(f"{a.cohort}/determinants.tsv"), delimiter="\t")}
R = [len(det.get(b, [])) for b, r in lab.items() if r.get(a.drug) == "R"]
S = [len(det.get(b, [])) for b, r in lab.items() if r.get(a.drug) == "S"]
gap = (st.mean(R) - st.mean(S)) if R and S else 0.0
sys.exit(0 if gap < a.min_gap else 1)  # exit 0 = disproof (no separation) = KILLED
