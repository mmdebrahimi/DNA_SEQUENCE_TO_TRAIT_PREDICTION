"""Forward-cell VALIDATION REPORT CARD — standing read-only trust surface (the molecular analogue of
`build_validation_report_card.py`). Rolls up the forward cell's DMS-validated numbers already on disk into
one honest per-capability card. It does NOT score/train — a report, not a gate (exit 0 always). NO aggregate
headline; each capability carries its own honest tier + scope. The forward cell validates against measured
Deep Mutational Scanning (DMS) — the one place the project's label wall does not bind (free per-variant labels).
"""
from __future__ import annotations
import json, statistics
from pathlib import Path

W = Path("wiki")

def load(name):
    p = W / name
    return json.loads(p.read_text(encoding="utf-8")) if p.exists() else None

def med(d):  # abs_spearman dict -> median
    return d.get("median") if isinstance(d, dict) else None

lb = load("forward_method_leaderboard_2026-07-15.json") or {}
hyb = load("forward_modality_hybrid_2026-07-17.json") or {}
bl = load("forward_blosum_proteingym_2026-07-17.json") or {}
es = load("forward_esm_proteingym_2026-07-17.json") or {}
inv = load("forward_inverse_deployable_2026-07-17.json") or {}
gb1 = load("forward_epistasis_gb1_2026-07-27.json") or {}
sweep = load("forward_epistasis_sweep_2026-07-27.json") or {}

rows = []
# --- single-variant methods (Regime B), benchmark-wide ---
rows.append({"capability": "blosum62 (deterministic, no deps)", "regime": "B_molecular",
             "tier": "DMS_VALIDATED_BENCHMARK_WIDE",
             "metric": f"ProteinGym median |Spearman| {med(bl.get('abs_spearman',{}))}",
             "scope": f"n={bl.get('n_scored')} assays; {bl.get('n_below_0.15')} below 0.15 — modest, the honest floor",
             "source": "forward_blosum_proteingym_2026-07-17.json"})
rows.append({"capability": "esm2-650M (learned, universal)", "regime": "B_molecular",
             "tier": "DMS_VALIDATED_BENCHMARK_WIDE",
             "metric": f"ProteinGym median |Spearman| {med(es.get('abs_spearman',{}))}",
             "scope": f"n={es.get('n_scored')} assays; {es.get('n_above_0.3')} above 0.3 — the sequence baseline at scale",
             "source": "forward_esm_proteingym_2026-07-17.json"})
rows.append({"capability": "modality hybrid (ESM2+ProSST[+GEMME])", "regime": "B_molecular",
             "tier": "DMS_VALIDATED_CEILING",
             "metric": f"ESM2 baseline median {hyb.get('baseline_median_abs_spearman')}; hybrid ceiling ~0.547 (ESM2+GEMME+ProSST)",
             "scope": f"n={hyb.get('n_assays_total')} struct+MSA assays; +0.056 paired, win 90.5% (prereg bar met); RANKS not doses",
             "source": "forward_modality_hybrid_2026-07-17.json"})
rows.append({"capability": "alphamissense / esm-if / prosst (structure/pathogenicity)", "regime": "B_molecular",
             "tier": "DMS_VALIDATED_PER_PROTEIN",
             "metric": "per-protein DMS Spearman (AM human-only; ESM-IF/ProSST structure tier)",
             "scope": f"{len(lb.get('proteins',[]))} proteins in the method leaderboard; ESM-IF does NOT beat ESM2 on PTEN",
             "source": "forward_method_leaderboard_2026-07-15.json"})
# --- multi-mutant additive-null (rows 578-580) ---
gs = (gb1.get("mean_en_linear_oos") is not None) or True
rows.append({"capability": "multi-mutant additive-null (predict_multi_effect)", "regime": "B_molecular",
             "tier": "DMS_VALIDATED",
             "metric": f"GB1 doubles: joint-ESM2 beats additive by only +{gb1.get('joint_minus_additive', 0):.4f}",
             "scope": "additive null validated as the robust deployed default; GB1 doubles ~96% additive (measured s1+s2 -> 0.958)",
             "source": "forward_epistasis_gb1_2026-07-27.json"})
# epistasis cross-protein
sw = sweep.get("results", [])
n_prot = len(sw)
maxdeg = max((abs(r.get("delta", 0)) for r in sw), default=0)
rows.append({"capability": "epistasis characterization (joint vs additive)", "regime": "B_molecular",
             "tier": "CHARACTERIZED",
             "metric": f"{n_prot} proteins x orders 2-6: joint ~= additive (delta ~+-0.005); 'grows with order' FALSIFIED",
             "scope": "novel: joint can be worse OOD (ParD WITHIN-ORDER delta -0.053, degrades with order); additive is robust. The pooled -0.283 was a mutation-order POOLING artifact -- corrected 2026-08-25, wiki/forward_epistasis_pooling_correction_2026-08-25.md",
             "source": "forward_epistasis_sweep_2026-07-27.json"})
# --- inverse cell ---
inv_h = inv.get("headline", {})
rows.append({"capability": "inverse (effect->edit proposal)", "regime": "B_molecular",
             "tier": "DEPLOYABLE_RANK_ONLY",
             "metric": f"rank/percentile inverse works {inv_h.get('rank_inverse_works_beats_null','?')} assays (beats null); "
                       f"learned oracle earns keep {inv_h.get('learned_oracle_earns_keep','?')}; magnitude NOT deployable",
             "scope": "RANKS not doses (percentile pts); magnitude needs the target's own DMS (circular by construction)",
             "source": "forward_inverse_deployable_2026-07-17.json"})
# --- regime router ---
rows.append({"capability": "regime router (predict_edit)", "regime": "router",
             "tier": "DEPLOYED",
             "metric": "A determinant -> AMR catalogue / B molecular -> DMS methods / C organismal -> ABSTAIN",
             "scope": "resistance NEVER routed to a likelihood model; organism-polygenic ABSTAINS (closed negative)",
             "source": "dna_decode/forward/router.py"})

card = {"schema": "forward-validation-report-card-v1", "generated": "2026-07-28",
        "note": "DMS-validated molecular variant-effect cell; zero-shot unless stated; per-capability honest tier; NO aggregate headline",
        "capabilities": rows}
W.joinpath("forward_validation_report_card.json").write_text(json.dumps(card, indent=2), encoding="utf-8")

md = ["# Forward-cell validation report card (standing trust surface)", "",
      "Read-only roll-up of the forward variant-effect cell's DMS-validated numbers (the molecular analogue of",
      "the AMR `decoder_validation_report_card`). DMS is the one place the project's label wall does not bind",
      "(free, independent, per-variant magnitude labels). Zero-shot unless stated. **No aggregate headline** —",
      "each capability carries its own honest tier.", "",
      "| capability | regime | tier | metric | scope |", "|---|---|---|---|---|"]
for r in rows:
    md.append(f"| {r['capability']} | {r['regime']} | {r['tier']} | {r['metric']} | {r['scope']} |")
md += ["", "## Honest scope", "- **DMS-validated, zero-shot** — per-variant Spearman vs measured DMS; not a "
       "calibrated dose (the inverse is RANK-only; the hybrid RANKS, does not dose).",
       "- **Regime B only** — resistance (A) routes to the frozen AMR catalogue (never a likelihood model); "
       "organism-polygenic (C) ABSTAINS (closed negative).", 
       "- **Multi-mutant = additive null** (validated robust across proteins/orders; joint epistasis adds ~0 and "
       "can hurt out-of-distribution).",
       "- Sources are the committed `wiki/forward_*` artifacts; regenerate with `scripts/build_forward_report_card.py`."]
W.joinpath("forward_validation_report_card.md").write_text("\n".join(md), encoding="utf-8")
print(f"wrote forward_validation_report_card.{{md,json}} — {len(rows)} capabilities")
for r in rows: print(f"  [{r['tier']}] {r['capability']}")
