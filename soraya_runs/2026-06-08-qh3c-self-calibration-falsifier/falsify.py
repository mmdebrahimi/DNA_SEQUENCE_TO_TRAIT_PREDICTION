"""Dual falsifier for the self-calibrating-tool strands (hypothesise #1 + #2).

#1 AUTO-THRESHOLD: can leave-one-out balanced-accuracy threshold selection recover the
   biologically-correct per-organism cipro QRDR-point threshold WITHOUT hand-tuning?
   Ground truth: Campylobacter -> 1 (single gyrA T86I), Klebsiella -> 2 (E.coli-style).
   KILL if LOO-selected threshold disagrees with ground truth on >=2 of the GT organisms.

#2 AUTO-INTRINSIC-FLAG: does "determinant present in >=90% of BOTH R and S strains" flag
   the known intrinsic (blaOXA-51-family in Acinetobacter) WITHOUT falsely flagging a
   discriminative acquired determinant (OXA-23)?
   KILL if OXA-51-family is not flagged, OR a strong discriminative determinant IS flagged.

Pure logic over cached AMRFinder main.tsv. No Docker, no downloads.
"""
import sys, json
from pathlib import Path
sys.path.insert(0, ".")
from dna_decode.eval.amr_rules import qrdr_point_determinants, cipro_determinants_from_main

RAW = Path("data/raw")

def load_labels(slug):
    out = {}
    for ln in (RAW / slug / "selected.tsv").read_text().splitlines():
        if "\t" in ln:
            a, rs = ln.split("\t")
            out[a] = rs.strip()
    return out

def confusion(counts_labels, thr):
    tp = fp = tn = fn = 0
    for n, lab in counts_labels:
        pred = "R" if n >= thr else "S"
        if (lab, pred) == ("R", "R"): tp += 1
        elif (lab, pred) == ("S", "S"): tn += 1
        elif (lab, pred) == ("S", "R"): fp += 1
        else: fn += 1
    N = len(counts_labels)
    sens = tp / (tp + fn) if (tp + fn) else 0.0
    spec = tn / (tn + fp) if (tn + fp) else 0.0
    return {"acc": (tp + tn) / N, "bal_acc": 0.5 * (sens + spec),
            "sens": sens, "spec": spec, "tp": tp, "tn": tn, "fp": fp, "fn": fn}

def loo_select_threshold(counts_labels, thr_grid=(1, 2, 3)):
    """For each held-out strain, pick the threshold maximizing balanced accuracy on the
    REST, then predict the held-out. Report the modal selected threshold + LOO accuracy."""
    sel = []
    correct = 0
    for i in range(len(counts_labels)):
        train = counts_labels[:i] + counts_labels[i+1:]
        best_thr, best_ba = None, -1.0
        for thr in thr_grid:
            ba = confusion(train, thr)["bal_acc"]
            if ba > best_ba:
                best_ba, best_thr = ba, thr
        sel.append(best_thr)
        n, lab = counts_labels[i]
        pred = "R" if n >= best_thr else "S"
        if pred == lab: correct += 1
    from collections import Counter
    modal = Counter(sel).most_common(1)[0][0]
    return {"modal_threshold": modal, "selected_thresholds": dict(Counter(sel)),
            "loo_acc": correct / len(counts_labels)}

# ---------- #1 AUTO-THRESHOLD ----------
print("=" * 70)
print("FALSIFIER #1 — AUTO-THRESHOLD (LOO over cipro QRDR-point counts)")
print("=" * 70)
GROUND_TRUTH = {"campylobacter_ciprofloxacin": 1, "klebsiella_cipro": 2}
cipro_cohorts = ["campylobacter_ciprofloxacin", "klebsiella_cipro",
                 "pseudomonas_aeruginosa_ciprofloxacin"]
res1 = {}
agree = disagree = 0
for slug in cipro_cohorts:
    runs = RAW / slug / "amrfinder_runs"
    if not runs.exists():
        continue
    labels = load_labels(slug)
    cl = []
    for acc, lab in labels.items():
        mt = runs / acc / "main.tsv"
        n = len(qrdr_point_determinants(mt)) if mt.exists() else 0
        cl.append((n, lab))
    loo = loo_select_threshold(cl)
    # also show fixed-threshold accuracy at the modal pick
    fixed = confusion(cl, loo["modal_threshold"])
    gt = GROUND_TRUTH.get(slug)
    verdict = ""
    if gt is not None:
        ok = (loo["modal_threshold"] == gt)
        verdict = f"GT={gt} -> {'AGREE' if ok else 'DISAGREE'}"
        agree += ok; disagree += (not ok)
    res1[slug] = {"n": len(cl), **loo, "fixed_acc_at_modal": fixed["acc"],
                  "ground_truth": gt, "verdict": verdict}
    print(f"\n{slug} (N={len(cl)})")
    print(f"  LOO modal threshold = {loo['modal_threshold']}  (picks: {loo['selected_thresholds']})")
    print(f"  LOO accuracy = {loo['loo_acc']:.3f} | fixed-acc@modal = {fixed['acc']:.3f}")
    if gt is not None:
        print(f"  ground truth = {gt}  ->  {verdict}")

kill1 = disagree >= 2
print(f"\n#1 RESULT: agree={agree} disagree={disagree} on ground-truth organisms "
      f"-> {'KILLED' if kill1 else 'SURVIVES'}")

# ---------- #2 AUTO-INTRINSIC-FLAG ----------
print("\n" + "=" * 70)
print("FALSIFIER #2 — AUTO-INTRINSIC-FLAG (R/S prevalence; flag intrinsic, spare acquired)")
print("=" * 70)
def determinant_prevalence(slug, drug):
    labels = load_labels(slug)
    runs = RAW / slug / "amrfinder_runs"
    R = {a for a, l in labels.items() if l == "R"}
    S = {a for a, l in labels.items() if l == "S"}
    from collections import defaultdict
    inR = defaultdict(int); inS = defaultdict(int)
    for acc, lab in labels.items():
        mt = runs / acc / "main.tsv"
        if not mt.exists():
            continue
        dets = cipro_determinants_from_main(mt, drug)  # broad drug-class determinants
        syms = {d["symbol"] for d in dets}
        for s in syms:
            if lab == "R": inR[s] += 1
            else: inS[s] += 1
    rows = []
    for s in set(inR) | set(inS):
        pr = inR[s] / len(R) if R else 0
        psv = inS[s] / len(S) if S else 0
        rows.append({"symbol": s, "prev_R": round(pr, 3), "prev_S": round(psv, 3),
                     "intrinsic_flag": (pr >= 0.90 and psv >= 0.90)})
    return sorted(rows, key=lambda r: (-r["intrinsic_flag"], -(r["prev_R"] + r["prev_S"])))

# Acinetobacter meropenem: OXA-51-family intrinsic; OXA-23 discriminative acquired.
ac = determinant_prevalence("acinetobacter_meropenem", "meropenem")
print("\nAcinetobacter x meropenem — determinant R/S prevalence (top 15):")
print(f"  {'symbol':28} {'prev_R':>7} {'prev_S':>7}  intrinsic_flag")
flagged, discrim_flagged = [], []
for r in ac[:15]:
    sym = r["symbol"]
    is_oxa51fam = any(t in sym for t in ("OXA-51","OXA-66","OXA-68","OXA-69","OXA-65",
                                          "OXA-90","OXA-100","OXA-113","OXA-312","OXA-64"))
    is_strong_acq = any(t in sym for t in ("OXA-23","OXA-24","OXA-40","OXA-72","NDM","IMP","VIM","KPC"))
    print(f"  {sym:28} {r['prev_R']:>7} {r['prev_S']:>7}  {r['intrinsic_flag']}"
          + ("   <- OXA-51-fam(intrinsic)" if is_oxa51fam else "")
          + ("   <- strong acquired" if is_strong_acq else ""))
    if r["intrinsic_flag"]:
        flagged.append(sym)
        if is_strong_acq: discrim_flagged.append(sym)

oxa51_flagged = any(any(t in s for t in ("OXA-51","OXA-66","OXA-68","OXA-69","OXA-65",
                    "OXA-90","OXA-100","OXA-113","OXA-312","OXA-64")) for s in flagged)
kill2 = (not oxa51_flagged) or bool(discrim_flagged)
print(f"\n  intrinsic-flagged determinants: {flagged}")
print(f"  OXA-51-family flagged (want True): {oxa51_flagged}")
print(f"  strong-acquired falsely flagged (want empty): {discrim_flagged}")
print(f"\n#2 RESULT: -> {'KILLED' if kill2 else 'SURVIVES'}")

# ---------- summary ----------
out = {"falsifier_1_auto_threshold": {"per_organism": res1, "agree": agree,
        "disagree": disagree, "verdict": "KILLED" if kill1 else "SURVIVES"},
       "falsifier_2_auto_intrinsic_flag": {"acinetobacter_prevalence": ac,
        "intrinsic_flagged": flagged, "oxa51_flagged": oxa51_flagged,
        "strong_acquired_falsely_flagged": discrim_flagged,
        "verdict": "KILLED" if kill2 else "SURVIVES"}}
Path("soraya_runs/2026-06-08-qh3c-self-calibration-falsifier/falsify_result.json").write_text(
    json.dumps(out, indent=2))
print("\n" + "=" * 70)
print(f"BOTH-SURVIVE: {(not kill1) and (not kill2)}")
print("wrote falsify_result.json")
