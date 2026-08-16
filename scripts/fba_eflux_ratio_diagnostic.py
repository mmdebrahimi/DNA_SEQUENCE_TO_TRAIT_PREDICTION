import sys; sys.path.insert(0,'.'); sys.path.insert(0,'scripts')
import numpy as np
from cobra.flux_analysis import single_gene_deletion
from dna_decode.fba.model import load_model, wildtype_growth
from dna_decode.fba.fitness_browser import open_db, carbon_conditions, load_records, apply_carbon_condition, ESSENTIAL_FITNESS
from dna_decode.fba.conditional_essentiality import conditionally_essential_genes
from fba_eflux_bridge import build_condition_expression, apply_eflux, FRAC

def main():
    model=load_model(); conn=open_db(); conds_all=carbon_conditions(conn,model)
    expr,_=build_condition_expression(conds_all)
    keys=sorted(expr); conds={k:conds_all[k] for k in keys}
    recs=load_records(conn,conds,gene_filter={g.id for g in model.genes},threshold=ESSENTIAL_FITNESS)
    genes=[r.gene_id for r in conditionally_essential_genes(recs)][:40]
    allex=tuple(conds_all.values())
    for cond in ("D-Galactose","Potassium acetate"):
        out={}
        for arm in ("baseline","eflux"):
            with model:
                apply_carbon_condition(model,conds[cond],all_carbon=allex)
                if arm=="eflux": apply_eflux(model,expr[cond])
                wt=wildtype_growth(model)
                res=single_gene_deletion(model,gene_list=[model.genes.get_by_id(g) for g in genes],processes=1)
                out[arm]={next(iter(r["ids"])):(0.0 if r["growth"]!=r["growth"] else r["growth"]/wt)
                          for _,r in res.iterrows()}
        b,e=out["baseline"],out["eflux"]
        d=np.array([abs(b[g]-e[g]) for g in genes])
        print(f"\n{cond} (n={len(genes)} genes)")
        print(f"  ratio |delta| : max={d.max():.4f}  mean={d.mean():.5f}  n>0.01: {(d>0.01).sum()}/{len(genes)}")
        print(f"  near {FRAC} thresh: {sum(1 for g in genes if abs(b[g]-FRAC)<0.05)}")
        print(f"  base->eflux   : {[(round(b[g],3),round(e[g],3)) for g in genes[:6]]}")
main()
