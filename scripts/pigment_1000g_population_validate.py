"""Population-frequency VALIDATION of the IrisPlex eye-colour cell on real 1000G genomes.

The IrisPlex cell (dna_decode/pigment) is "faithful-to-published-model, NOT scored" — its per-individual
label source (openSNP) was permanently deleted 2025-04-30. This runs the ONE measured check still free:
predict eye colour for all 1000G samples and confirm the PER-POPULATION distribution matches known
biology — blue eyes are common in EUROPEANS and near-absent in AFR/EAS/SAS. It is a POPULATION-level
sanity check (not per-individual concordance), and that scope is stated in the artifact.

ANTI-FABRICATION rails (load-bearing):
  - The 6 IrisPlex SNP GRCh38 coordinates are PINNED FROM ENSEMBL AT RUNTIME (reusing the project's
    verify_sentinel_coords Ensembl fetch), never hardcoded from memory — a wrong coord is invisible.
  - STRAND is auto-harmonized against each site's real 1000G alleles: the coefficient table's counted
    allele must appear (forward) at the fetched site, else its COMPLEMENT must — if neither, the SNP is a
    hard mismatch and the run ABORTS (never a silent wrong count). The genotype is normalized to the
    table's strand so `predict_eye_color`'s internal counted allele is correct unchanged.
  - Nothing is fabricated: if Ensembl or the 1000G fetch is unreachable the run reports UNVERIFIED and
    exits non-zero rather than emitting a made-up distribution.

Uses the Docker-free fetcher (scripts/fetch_1000g_region.py) + the committed 1000G_3202_samples.ped
(Superpopulation column). Reversible, no money.

    uv run python scripts/pigment_1000g_population_validate.py
Exit: 0 = ran + expected pattern held, 1 = ran but pattern FAILED (investigate), 3 = data/network unreachable.
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

from dna_decode.pigment.irisplex import IRISPLEX_SNPS, predict_eye_color  # noqa: E402
from scripts.fetch_1000g_region import fetch_region  # noqa: E402
from scripts.verify_sentinel_coords import _http_fetch as ensembl_fetch  # noqa: E402

_COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}
PED = REPO / "data" / "pgx_1000g" / "1000G_3202_samples.ped"


def ensembl_grch38_mapping(rsid: str) -> dict:
    """{chrom, pos, alleles:[...]} for rsid on GRCh38 (from Ensembl). Raises on no unique chr mapping."""
    data = ensembl_fetch(rsid)
    quals = [m for m in data.get("mappings", [])
             if m.get("assembly_name") == "GRCh38" and m.get("coord_system") == "chromosome"]
    if len(quals) != 1:
        raise RuntimeError(f"{rsid}: {len(quals)} GRCh38 chromosome mappings (need exactly 1)")
    m = quals[0]
    return {"chrom": "chr" + str(m["seq_region_name"]), "pos": int(m["start"]),
            "alleles": (m.get("allele_string") or "").split("/")}


def _parse_snp_record(vcf_text: str, pos: int) -> tuple[str, list[str], list[str], dict[str, str]]:
    """From a fetched single-region VCF, return (ref, alts, samples, {sample: 'A/G'}) for the record at `pos`."""
    samples: list[str] = []
    for line in vcf_text.splitlines():
        if line.startswith("##"):
            continue
        if line.startswith("#CHROM"):
            samples = line.rstrip("\n").split("\t")[9:]
            continue
        cols = line.rstrip("\n").split("\t")
        if len(cols) < 10 or not cols[1].isdigit() or int(cols[1]) != pos:
            continue
        ref, alts = cols[3], cols[4].split(",")
        if len(ref) != 1 or any(len(a) != 1 for a in alts):
            continue  # SNP only; skip a co-located indel
        code = {"0": ref, **{str(i + 1): a for i, a in enumerate(alts)}}
        fmt = cols[8].split(":")
        gi = fmt.index("GT") if "GT" in fmt else 0
        gts: dict[str, str] = {}
        for s, cell in zip(samples, cols[9:]):
            raw = cell.split(":")[gi].replace("|", "/").split("/")
            bases = [code.get(a) for a in raw if a in code]
            if len(bases) == 2:
                gts[s] = "/".join(bases)
        return ref, alts, samples, gts
    raise RuntimeError(f"no SNP record at pos {pos} in the fetched slice")


def build_genotypes() -> tuple[dict[str, dict], list[str]]:
    """Fetch all 6 IrisPlex SNPs; return ({sample: {rsid: 'A/G' on the table's strand}}, notes)."""
    per_sample: dict[str, dict] = defaultdict(dict)
    notes: list[str] = []
    for rsid, counted, _b_int, _b_brown in IRISPLEX_SNPS:
        mp = ensembl_grch38_mapping(rsid)
        chrom, pos = mp["chrom"], mp["pos"]
        with tempfile.TemporaryDirectory() as td:
            out = Path(td) / f"{rsid}.vcf"
            fetch_region(chrom, pos, pos, out, verbose=False)
            ref, alts, _samples, gts = _parse_snp_record(out.read_text(encoding="utf-8"), pos)
        site = {ref, *alts}
        if counted in site:
            strand = "fwd"
        elif _COMP[counted] in site:
            strand = "rev"
        else:
            raise RuntimeError(f"{rsid}: counted allele {counted!r} nor its complement is at the 1000G site "
                               f"{chrom}:{pos} (alleles {sorted(site)}) — coord/allele mismatch, ABORT")
        notes.append(f"{rsid} {chrom}:{pos} site={ref}/{','.join(alts)} counted={counted} strand={strand}")
        for s, gt in gts.items():
            per_sample[s][rsid] = gt if strand == "fwd" else "/".join(_COMP[b] for b in gt.split("/"))
    return per_sample, notes


def superpops() -> dict[str, str]:
    out = {}
    for line in PED.read_text(encoding="utf-8").splitlines()[1:]:
        f = line.split()
        if len(f) >= 7:
            out[f[1]] = f[6]
    return out


def main(argv=None) -> int:
    try:
        per_sample, notes = build_genotypes()
    except Exception as e:  # noqa: BLE001
        print(f"UNVERIFIED: {type(e).__name__}: {e}", file=sys.stderr)
        return 3
    sp = superpops()
    by_pop: dict[str, Counter] = defaultdict(Counter)
    n_all6 = 0
    for s, genos in per_sample.items():
        if len(genos) != len(IRISPLEX_SNPS):
            continue  # need all 6 (rs12913832 mandatory anyway)
        n_all6 += 1
        try:
            call = predict_eye_color(genos).call
        except Exception:  # noqa: BLE001
            continue
        by_pop[sp.get(s, "UNK")][call] += 1

    dist = {}
    for pop, c in sorted(by_pop.items()):
        tot = sum(c.values())
        dist[pop] = {"n": tot, "blue": round(c["blue"] / tot, 4), "intermediate": round(c["intermediate"] / tot, 4),
                     "brown": round(c["brown"] / tot, 4)}

    # KNOWN-biology expectation (the sanity gate): EUR carries substantial blue; AFR/EAS/SAS ~all brown.
    eur_blue = dist.get("EUR", {}).get("blue", 0.0)
    afr_brown = dist.get("AFR", {}).get("brown", 0.0)
    eas_brown = dist.get("EAS", {}).get("brown", 0.0)
    sas_brown = dist.get("SAS", {}).get("brown", 0.0)
    passed = (eur_blue >= 0.30 and afr_brown >= 0.90 and eas_brown >= 0.90 and sas_brown >= 0.85)

    res = {
        "cell": "typing:human:pigment (IrisPlex eye colour)", "date": datetime.date.today().isoformat(),
        "substrate": "1000G 30x phased 3202-sample panel (GRCh38); coords Ensembl-pinned at runtime",
        "n_scored_all6": n_all6, "distribution_by_superpopulation": dist,
        "expectation": {"EUR_blue>=0.30": eur_blue, "AFR_brown>=0.90": afr_brown,
                        "EAS_brown>=0.90": eas_brown, "SAS_brown>=0.85": sas_brown},
        "known_biology_pattern_held": passed,
        "snp_provenance": notes,
        "honesty": ("POPULATION-level sanity check, NOT per-individual concordance (openSNP, the free "
                    "per-individual label source, was deleted 2025-04-30). Confirms the published IrisPlex "
                    "model reproduces the known geography of eye colour on real genomes; it does NOT score "
                    "individual accuracy. Includes related samples (3202 panel)."),
    }
    out_json = REPO / "wiki" / f"pigment_1000g_population_{res['date']}.json"
    out_md = REPO / "wiki" / f"pigment_1000g_population_{res['date']}.md"
    out_json.write_text(json.dumps(res, indent=2), encoding="utf-8")
    lines = [f"# IrisPlex eye-colour — 1000G population-frequency validation ({res['date']})", "",
             f"**Known-biology pattern held: {passed}** (n={n_all6} scored on all 6 SNPs).", "",
             "| superpop | n | P(blue) | P(intermediate) | P(brown) |", "|---|---|---|---|---|"]
    for pop, d in dist.items():
        lines.append(f"| {pop} | {d['n']} | {d['blue']} | {d['intermediate']} | {d['brown']} |")
    lines += ["", res["honesty"], "", "SNP provenance (Ensembl-pinned GRCh38 + strand):"]
    lines += [f"- {n}" for n in notes]
    out_md.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print("\n".join(lines))
    print(f"\n[-> {out_json}]\n[-> {out_md}]")
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
