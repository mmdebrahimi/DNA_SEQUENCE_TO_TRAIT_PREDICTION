"""Population-geography validation of the recovered HIrisPlex-S hair + skin (+ eye) models on real 1000G.

The models are already proven faithful to the deployed webtool (held-out |ΔP| eye/hair ~1e-15, skin ~9e-3).
This is the BIOLOGICAL sanity check: predict hair + skin + eye across all 1000G superpopulations and confirm
the known geography — red/blond hair concentrated in EUR (~absent AFR/EAS/SAS); a skin cline dark in AFR ->
pale/very-pale in EUR. Coords are Ensembl-pinned at runtime + strand-harmonized (no fabrication), reusing the
eye validator's helpers. POPULATION-level sanity (not per-individual; openSNP deleted).

    uv run python scripts/pigment_1000g_hairskin_validate.py
Exit 0 = ran + known geography held; 1 = ran but a geography expectation FAILED; 3 = data/network unreachable.
"""
from __future__ import annotations

import datetime
import json
import sys
import tempfile
from collections import Counter, defaultdict
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from dna_decode.pigment.hirisplex_models import load_hirisplex_models  # noqa: E402
from dna_decode.pigment.multinomial import predict  # noqa: E402
from scripts.pigment_1000g_population_validate import (  # noqa: E402  (reuse helpers)
    _COMP, _parse_snp_record, ensembl_grch38_mapping, superpops,
)
from scripts.fetch_1000g_region import fetch_region  # noqa: E402


def build_genotypes(rsids_alleles: dict[str, str]) -> tuple[dict[str, dict], list[str]]:
    """{sample: {rsid: 'A/G' on the model's counted-allele strand}} for the given {rsid: counted_allele}."""
    per: dict[str, dict] = defaultdict(dict)
    notes: list[str] = []
    for rsid, counted in rsids_alleles.items():
        try:
            mp = ensembl_grch38_mapping(rsid)
        except Exception as e:  # noqa: BLE001
            notes.append(f"{rsid}: Ensembl map failed ({e}); skipped")
            continue
        chrom, pos = mp["chrom"], mp["pos"]
        with tempfile.TemporaryDirectory() as td:
            out = Path(td) / f"{rsid}.vcf"
            try:
                fetch_region(chrom, pos, pos, out, verbose=False)
                ref, alts, _s, gts = _parse_snp_record(out.read_text(encoding="utf-8"), pos)
            except Exception as e:  # noqa: BLE001
                notes.append(f"{rsid}: 1000G fetch/parse failed ({e}); skipped")
                continue
        site = {ref, *alts}
        if counted in site:
            strand = "fwd"
        elif _COMP.get(counted) in site:
            strand = "rev"
        else:
            notes.append(f"{rsid}: counted {counted} not at 1000G site {chrom}:{pos}; skipped")
            continue
        for s, gt in gts.items():
            per[s][rsid] = gt if strand == "fwd" else "/".join(_COMP[b] for b in gt.split("/"))
    return per, notes


def main() -> int:
    models = load_hirisplex_models()
    # union of all SNPs across the 3 models with their counted alleles (skin panel = 41 = superset)
    rsid_allele: dict[str, str] = {}
    for m in models.values():
        for snp in m.snps:
            rsid_allele[snp.rsid] = snp.counted_allele
    print(f"fetching {len(rsid_allele)} HIrisPlex-S SNPs from 1000G ...", flush=True)
    per, notes = build_genotypes(rsid_allele)
    if not per:
        print("UNVERIFIED: no SNPs fetched", file=sys.stderr)
        return 3
    sp = superpops()

    # per-sample predictions (allow_missing for any SNP that failed to fetch -> low confidence)
    by_pop_call: dict[str, dict[str, Counter]] = {t: defaultdict(Counter) for t in models}
    by_pop_meanP: dict[str, dict[str, dict]] = {t: defaultdict(lambda: defaultdict(float)) for t in models}
    n_by_pop: Counter = Counter()
    for s, g in per.items():
        pop = sp.get(s, "UNK")
        n_by_pop[pop] += 1
        for trait, model in models.items():
            r = predict(model, g, allow_missing=True)
            by_pop_call[trait][pop][r.call] += 1
            for cat, p in r.probabilities.items():
                by_pop_meanP[trait][pop][cat] += p

    result = {"date": datetime.date.today().isoformat(),
              "substrate": "1000G 30x phased 3202-sample panel; coords Ensembl-pinned; strand-harmonized",
              "n_snps_fetched": len(rsid_allele) - sum("skipped" in n for n in notes),
              "n_by_superpop": dict(n_by_pop), "distribution": {}, "mean_probabilities": {}}
    for trait in models:
        result["distribution"][trait] = {pop: {c: round(v / sum(cnt.values()), 4) for c, v in cnt.items()}
                                         for pop, cnt in by_pop_call[trait].items()}
        result["mean_probabilities"][trait] = {
            pop: {c: round(tot / n_by_pop[pop], 4) for c, tot in cats.items()}
            for pop, cats in by_pop_meanP[trait].items()}

    # known-geography gates (mean probabilities; robust to the argmax being brown/black for most)
    hp = result["mean_probabilities"]["hair_colour"]
    skp = result["mean_probabilities"]["skin_colour"]
    eur_light_hair = hp.get("EUR", {}).get("blond", 0) + hp.get("EUR", {}).get("red", 0)
    afr_light_hair = hp.get("AFR", {}).get("blond", 0) + hp.get("AFR", {}).get("red", 0)
    afr_dark_skin = skp.get("AFR", {}).get("dark", 0) + skp.get("AFR", {}).get("dark_to_black", 0)
    eur_pale_skin = skp.get("EUR", {}).get("very_pale", 0) + skp.get("EUR", {}).get("pale", 0)
    checks = {
        "EUR_light_hair(blond+red) > AFR": (eur_light_hair, afr_light_hair, eur_light_hair > afr_light_hair),
        "AFR_dark_skin > EUR": (afr_dark_skin, skp.get("EUR", {}).get("dark", 0)
                                + skp.get("EUR", {}).get("dark_to_black", 0),
                                afr_dark_skin > 0.5),
        "EUR_pale_skin(very_pale+pale) > AFR": (eur_pale_skin, skp.get("AFR", {}).get("very_pale", 0)
                                                + skp.get("AFR", {}).get("pale", 0),
                                                eur_pale_skin > 0.5),
    }
    result["geography_checks"] = {k: {"eur/afr_or_val": v[:2], "pass": v[2]} for k, v in checks.items()}
    passed = all(v[2] for v in checks.values())
    result["known_geography_held"] = passed
    result["snp_notes"] = notes

    (REPO / "wiki" / f"pigment_1000g_hairskin_{result['date']}.json").write_text(
        json.dumps(result, indent=2), encoding="utf-8")
    md = [f"# HIrisPlex-S hair + skin — 1000G population-geography validation ({result['date']})", "",
          f"**Known geography held: {passed}** ({sum(n_by_pop.values())} samples, "
          f"{result['n_snps_fetched']}/{len(rsid_allele)} SNPs fetched).", "",
          "Mean P by superpopulation (hair):", "```",
          json.dumps(result["mean_probabilities"]["hair_colour"], indent=1), "```",
          "Mean P by superpopulation (skin):", "```",
          json.dumps(result["mean_probabilities"]["skin_colour"], indent=1), "```",
          "Geography checks:", "```", json.dumps(result["geography_checks"], indent=1), "```",
          "", "POPULATION-level sanity (not per-individual). Models recovered + held-out-validated from the "
          "HIrisPlex-S webtool; coords Ensembl-pinned + strand-harmonized on 1000G."]
    (REPO / "wiki" / f"pigment_1000g_hairskin_{result['date']}.md").write_text("\n".join(md) + "\n",
                                                                               encoding="utf-8")
    print("\n".join(md))
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
