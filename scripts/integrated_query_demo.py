"""Demo: the rung-2 + rung-3 fused query on real data (2026-07-13).

Runs the integrated query for HIV RT mutations (using the committed RT sequence + the cached ESM masked
marginals — instant) and an E. coli honest-fallback case, and writes a result packet. The headline the demo
makes visible: K103N is ESM-molecularly 'uncertain/benign' (damage_llr ~ 0) yet a MAJOR resistance mutation
(the catalog knows it) — the two rungs DISAGREE, and that disagreement is exactly why a curated catalog is
irreplaceable for antagonistically-selected resistance. Frozen surface untouched.
"""
from __future__ import annotations

import json
import sys
from datetime import date as _date
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO)); sys.path.insert(0, str(REPO / "scripts"))

from dna_decode.protein_effect import integration as IQ  # noqa: E402

RT_CACHE = REPO / "data" / "processed" / "hiv_rt_esm650m_masked_marginals.json"


def main():
    import hiv_esm_vs_catalog as ev
    seq = ev.rt_protein()
    logp = IQ.load_logp(RT_CACHE)
    queries = []
    for mut in ("K103N", "L100I", "Y181C", "M184V"):     # famous NNRTI/NRTI DRMs
        queries.append(IQ.integrated_query(mut, gene="RT", organism="HIV-1", sequence=seq, logp=logp))
    ecoli = {"query": {"organism": "Escherichia coli", "gene": "gyrA", "mutation": "S83L"},
             "known_phenotype_rung3": IQ.known_phenotype("S83L", gene="gyrA", organism="Escherichia coli"),
             "note": "molecular rung would need --gene gyrA (UniProt fetch + ESM); shown here for the rung-3 fallback"}
    res = {
        "artifact": "integrated_query_demo", "schema": "integrated-query-demo-v1", "date": str(_date.today()),
        "headline": ("A single query returns BOTH the molecular effect (universal ESM rank) and the known "
                     "AMR phenotype (catalog lookup). K103N is the load-bearing example: ESM rates it "
                     "molecularly benign/uncertain, yet the catalog knows it is a major NNRTI resistance "
                     "mutation — the rungs DISAGREE, which is precisely why the curated catalog is "
                     "irreplaceable for antagonistically-selected resistance (a chemically-innocent "
                     "substitution that confers resistance)."),
        "hiv_queries": queries, "ecoli_fallback": ecoli,
    }
    out = REPO / "wiki" / f"integrated_query_demo_{_date.today()}.json"
    out.write_text(json.dumps(res, indent=2), encoding="utf-8")
    print("=== rung-2 (molecular) + rung-3 (phenotype) fused query — HIV RT ===")
    print(f"{'mut':>6} {'damage_llr':>11} {'pctile':>7} {'mol_hint':>17} | {'known_DRM':>9} resistant_drugs")
    for q in queries:
        m, p = q["molecular_effect_rung2"], q["known_phenotype_rung3"]
        drugs = ",".join(d[:3] for d in (p["resistant_drugs"] or []))
        print(f"{q['query']['mutation']:>6} {m['damage_llr']:>11} {m['position_percentile']:>7} "
              f"{m['direction_hint']:>17} | {str(p['is_known_resistance_mutation']):>9} {drugs}")
    print(f"\nE. coli gyrA S83L rung-3: catalog={ecoli['known_phenotype_rung3']['catalog']} (honest fallback)")
    print(f"[wrote {out}]")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
