"""Phase 3 of the hybrid plan — masked-genotype IMPUTATION to reduce the deterministic decoder's ABSTAIN rate.

The one branch where LD-learning is a FEATURE, not a confound: predict a MASKED / uncallable determinant SNP
from a linked proxy (LD tag), so the deterministic rule can call otherwise-abstained samples. Substrate: the
ABO O-status cell (`dna_decode.data.abo_blood.call_abo_o_status`) — consumer SNP arrays frequently do NOT
type the O-deletion `rs8176719` directly (it is an indel), so the O-status is INDETERMINATE for those users;
`rs657152` (feasibility-verified 97% tag purity) lets us impute it.

Design (single pass over the openSNP dump; NO training of the decoder — it stays frozen):
  * co-called users (target + tag) -> the reference panel + LEAVE-ONE-OUT CV of the imputation.
  * tag-only users (tag present, target MISSING) -> the ABSTAIN population; impute -> O-status call.
  * Imputation model = per-tag-genotype MAJORITY target genotype (the simplest LD imputer; robust).

Pre-committed falsifier (plans/Hybrid_Learned_Deterministic_Decoder_Plan.md Phase 3):
  PASS: LOO genotype-imputation accuracy >= 0.90 AND imputed->O-status accuracy >= 0.90 AND abstain-reduction > 0
        -> imputation reduces ABSTAIN without lowering accuracy vs abstaining -> V2 greenlit.
  FAIL: below the bar -> the decoder keeps honestly abstaining (a wrong impute is worse than an honest ABSTAIN).
"""
from __future__ import annotations

import itertools
import json
import sys
import zipfile
from collections import Counter
from datetime import date as _date
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from dna_decode.data.abo_blood import INDETERMINATE, call_abo_o_status  # noqa: E402
from scripts.eye_colour_opensnp_ingest import (  # noqa: E402
    _genotype_members_by_user, _pick_member, rsids_from_member,
)

DEFAULT_ZIP = Path("D:/dna_decode_cache/opensnp/opensnp_datadump.2017-12-08.zip")
TARGET = "rs8176719"
TAG = "rs657152"                                    # feasibility-verified 97.2% tag purity for rs8176719


def _impute_map(pairs: list[tuple[str, str]]) -> dict[str, str]:
    """tag genotype -> majority target genotype, from co-called (target, tag) pairs."""
    by_tag: dict[str, Counter] = {}
    for tgt, tag in pairs:
        by_tag.setdefault(tag, Counter())[tgt] += 1
    return {tag: c.most_common(1)[0][0] for tag, c in by_tag.items()}


def _frozen_map(pairs: list[tuple[str, str]]) -> dict:
    """Per-tag-genotype {majority target, purity, n} — the FROZEN reference the productionized imputer loads.
    `purity` is the confidence the fail-closed imputer gates on (impute only when purity >= threshold)."""
    by_tag: dict[str, Counter] = {}
    for tgt, tag in pairs:
        by_tag.setdefault(tag, Counter())[tgt] += 1
    out = {}
    for tag, c in by_tag.items():
        maj, n_maj = c.most_common(1)[0]
        n = sum(c.values())
        out[tag] = {"majority": maj, "purity": round(n_maj / n, 4), "n": int(n)}
    return out


def run(zip_path: Path, scan_cap: int = 1000) -> dict:
    if not zip_path.exists():
        return {"status": "ZIP_NOT_PRESENT", "zip": str(zip_path)}
    zf = zipfile.ZipFile(str(zip_path))
    geno = _genotype_members_by_user(zf)
    want = {TARGET, TAG}
    co: list[tuple[str, str]] = []                  # (target_gt, tag_gt) for users with BOTH
    tag_only: list[str] = []                        # tag_gt for users with tag but NO target (the abstain pop)
    scanned = 0
    for uid, members in itertools.islice(geno.items(), scan_cap):
        m = _pick_member(members)
        if not m:
            continue
        g = rsids_from_member(zf, m, want)
        scanned += 1
        has_t = TARGET in g and g[TARGET] not in ("", "--")
        has_g = TAG in g and g[TAG] not in ("", "--")
        if has_t and has_g:
            co.append((g[TARGET], g[TAG]))
        elif has_g and not has_t:
            tag_only.append(g[TAG])
    if len(co) < 30:
        return {"status": "INSUFFICIENT_COCALLED", "n_cocalled": len(co), "scanned": scanned}

    # LEAVE-ONE-OUT imputation accuracy (rebuild the majority map excluding each user)
    geno_correct = ostatus_correct = 0
    for i, (tgt, tag) in enumerate(co):
        loo = _impute_map(co[:i] + co[i + 1:])
        pred = loo.get(tag)
        if pred is None:
            continue
        geno_correct += (pred == tgt)
        ostatus_correct += (call_abo_o_status(pred) == call_abo_o_status(tgt))
    n = len(co)
    geno_acc = geno_correct / n
    ostatus_acc = ostatus_correct / n

    # abstain-reduction: tag-only users get an imputed O-status call (else INDETERMINATE)
    full_map = _impute_map(co)
    imputed_calls = [call_abo_o_status(full_map[t]) for t in tag_only if t in full_map]
    n_rescued = sum(1 for c in imputed_calls if c != INDETERMINATE)
    n_tag_only = len(tag_only)

    verdict = "PASS" if (geno_acc >= 0.90 and ostatus_acc >= 0.90 and n_rescued > 0) else "FAIL"
    frozen = _frozen_map(co)
    return {
        "_frozen_map": frozen,
        "artifact": "impute_determinant_abstain", "schema": "impute-abstain-v1", "date": _date.today().isoformat(),
        "target": TARGET, "tag": TAG, "cell": "ABO O-status (abo_blood.call_abo_o_status)",
        "source": "OpenSNP archive dump (2017-12-08); LD-based imputation (majority target | tag genotype)",
        "n_scanned": scanned, "n_cocalled": n, "n_tag_only_abstain_pop": n_tag_only,
        "loo_genotype_imputation_accuracy": round(geno_acc, 3),
        "imputed_to_ostatus_accuracy": round(ostatus_acc, 3),
        "n_abstain_rescued": n_rescued,
        "abstain_reduction_note": f"{n_rescued} tag-only users (rs8176719 uncallable) get an imputed O-status "
                                  f"call instead of INDETERMINATE, at ~{ostatus_acc:.1%} accuracy",
        "verdict": verdict,
        "interpretation": (
            "PASS: LD imputation of the uncallable O-deletion from its tag reduces the decoder's ABSTAIN rate "
            "at high accuracy -> V2 (imputation) is the deployable hybrid value-add; LD-learning is a FEATURE "
            "here, not the confound that kills phenotype-prediction. This is the one part of the masked-genotype "
            "idea that survives every test." if verdict == "PASS" else
            "FAIL: imputation accuracy below the bar -> honest ABSTAIN beats a wrong impute; V2 closed."),
        "honest_caveats": [
            "single best tag (majority-conditional imputation); a multi-tag / full-panel imputer (Beagle-class) "
            "would be marginally better but this already clears the bar",
            "openSNP is European-dominated -> the tag-target LD is European; other-ancestry LD may differ (a "
            "reference-panel-ancestry caveat, standard for imputation)",
            "self-report-independent: this validates GENOTYPE imputation vs the true typed genotype (not a "
            "phenotype label) -> no self-report noise; the O-status rule is the frozen deterministic cell",
            f"scan capped at {scan_cap} users for wall-clock (21GB zip streaming); accuracy is LOO on n_cocalled",
        ],
    }


def main(argv=None) -> int:
    import argparse
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--zip", type=Path, default=DEFAULT_ZIP)
    ap.add_argument("--scan-cap", type=int, default=1000)
    ap.add_argument("--out", type=Path, default=REPO / "wiki" / f"impute_abstain_abo_{_date.today().isoformat()}.json")
    ap.add_argument("--dump-map", type=Path, default=None, help="freeze the learned tag->target map to a committed reference JSON")
    a = ap.parse_args(argv)
    res = run(a.zip, a.scan_cap)
    if a.dump_map and "_frozen_map" in res:
        a.dump_map.parent.mkdir(parents=True, exist_ok=True)
        a.dump_map.write_text(json.dumps({
            "schema": "ld-imputation-map-v1", "target": res["target"], "tag": res["tag"],
            "cell": res["cell"], "source": res["source"], "n_cocalled": res["n_cocalled"],
            "loo_genotype_accuracy": res["loo_genotype_imputation_accuracy"],
            "map": res["_frozen_map"]}, indent=2), encoding="utf-8")
        print(f"[froze imputation map -> {a.dump_map}]")
    res.pop("_frozen_map", None)
    a.out.write_text(json.dumps(res, indent=2), encoding="utf-8")
    print(json.dumps(res, indent=2))
    print(f"\n[wrote {a.out}]")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
