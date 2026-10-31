"""Emit `wiki/glm_generator_sampling_sweep.json` -- G-D's sampling-sweep unblock trigger.

WHY THIS IS A SCRIPT AND NOT A FILE COPY. That exact filename is the `file-exists` expiry test
for constraint `gd-finetune-gated-on-sampling-sweep-not-compute`, and both CLAUDE.md and the
2026-10-30 handoff say in as many words: do NOT create it BY HAND. Writing it on an unverified
sweep would fire a LoRA fine-tune's unblock trigger on a grid this project's own gate excludes --
the same class of error as discharging a bar without spending the lever.

So the discharge CONDITION is asserted here in code, against the real pulled artifacts, and the
file is refused unless it holds:

  every grid point is plan-eligible under `sweep_gates.verdict_eligibility`

which for the frozen Step-4 grid means `eligible:verified` at the two NON-DEFAULT points
(T=0.8, T=0.6) and `eligible:default-sampling` at the T=1.0 anchor. The anchor is NOT
"unverified" -- it requested no non-default knob, so it has nothing to verify, and excluding it
would leave the sharpened points with nothing to be sharper than.

WHAT EMITTING THIS DOES AND DOES NOT DO. It discharges the SWEEP half of G-D's unblock trigger.
Do NOT read the remaining "8 gb gpu" half as unmet hardware: the constraint record's own analysis
is that the compute half has been satisfied since 2026-07-29 (free Kaggle T4s are 16 GB, and this
repo runs its kernels on them), and that the SWEEP was the binding half -- it had to report first,
because sampling and weights are two generator levers and running them together means neither can
be attributed. So on the EVIDENCE the trigger is now discharged.

G-D nonetheless stays BLOCKED, and for a reason that is not evidential: the user's standing call
on 2026-10-31 was to HOLD that edge, and the frozen plan's Step 6 forbids changing it. Lifting it
is an AUTHORITY call, surfaced and not taken. Nothing here launches a fine-tune.

    uv run python scripts/emit_glm_sweep_unblock_artifact.py \
        --summary wiki/glm_sampling_sweep_2026-10-31_vk_glm_sampling_sweep.json \
        --cells wiki/glm_sampling_sweep_2026-10-31_vk_cell_t1.0_k0.json ... \
        --source-kernel glm-sampling-sweep-vk
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from dna_decode.glm.sweep_gates import ELIGIBLE, verdict_eligibility  # noqa: E402

CANONICAL = ROOT / "wiki" / "glm_generator_sampling_sweep.json"


class NotDischargeable(SystemExit):
    pass


def assess(cell_paths: list) -> list:
    """One row per grid point: its sampling stamp and its code-owned eligibility."""
    rows = []
    for p in cell_paths:
        d = json.loads(Path(p).read_text(encoding="utf-8"))
        knob = d.get("knob_contract", {})
        stamp = d.get("sampling", {})
        rows.append({
            "artifact": str(Path(p).as_posix()),
            "temperature": stamp.get("temperature"),
            "top_k": stamp.get("top_k"),
            "is_recorded_default": stamp.get("is_recorded_default"),
            "requested_non_default": bool(knob.get("requested_non_default")),
            "knob_wired": bool(knob.get("knob_wired")),
            "knob_verified_end_to_end": bool(knob.get("knob_verified_end_to_end")),
            "verification_verdict": (knob.get("verification") or {}).get("verdict"),
            "eligibility": verdict_eligibility(knob),
        })
    return rows


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--summary", required=True)
    ap.add_argument("--cells", required=True, nargs="+")
    ap.add_argument("--source-kernel", required=True)
    ap.add_argument("--force", action="store_true",
                    help="overwrite an existing canonical artifact")
    ap.add_argument("--out", default=None,
                    help="where to write; defaults to the canonical trigger path. Exists so the "
                         "refusals can be TESTED without writing the real trigger file -- a "
                         "guard whose only exercise is the live path is a guard nobody runs.")
    a = ap.parse_args(argv)
    out = Path(a.out) if a.out else CANONICAL

    rows = assess(a.cells)
    not_eligible = [r for r in rows if r["eligibility"] not in ELIGIBLE]
    if not_eligible:
        for r in not_eligible:
            print("  T=%s -> %s" % (r["temperature"], r["eligibility"]), file=sys.stderr)
        raise NotDischargeable(
            "REFUSING to emit %s: %d of %d grid point(s) are not plan-eligible. The trigger is "
            "not dischargeable on this grid." % (out.name, len(not_eligible), len(rows)))

    verified = [r for r in rows if r["eligibility"] == "eligible:verified"]
    if not verified:
        raise NotDischargeable(
            "REFUSING to emit %s: no grid point carries a VERIFIED non-default knob. A grid of "
            "default-sampling points attests nothing about the knob and cannot discharge a "
            "trigger that exists to prove the knob is live." % out.name)

    if out.is_file() and not a.force:
        raise NotDischargeable("%s already exists; pass --force to overwrite" % out.name)

    summary = json.loads(Path(a.summary).read_text(encoding="utf-8"))
    payload = {
        "record": "glm-generator-sampling-sweep-unblock-v1",
        "emitted_by": "scripts/emit_glm_sweep_unblock_artifact.py",
        "source_kernel": a.source_kernel,
        "source_summary": str(Path(a.summary).as_posix()),
        "discharges": {
            "half": "sampling-sweep",
            "constraint": "gd-finetune-gated-on-sampling-sweep-not-compute",
            "family": "glm-generator-finetune-2026-10-07",
            "trigger_state": "DISCHARGED ON EVIDENCE. The constraint record's own analysis is "
                             "that the compute half was already satisfied (Kaggle T4 = 16 GB, "
                             "reachable since 2026-07-09) and that the sweep was the binding "
                             "half. Do NOT read the edge's '8 gb gpu' clause as unmet hardware.",
            "still_blocked_because": "G-D's blocked_by edge is HELD by explicit user call "
                                     "2026-10-31 and by the frozen plan's Step 6, NOT by missing "
                                     "hardware. Lifting it is an AUTHORITY call, not an executor "
                                     "task; this artifact does not unblock G-D.",
        },
        "knob_attestation": {
            "n_points": len(rows),
            "n_eligible": len(rows),
            "n_verified_non_default": len(verified),
            "points": rows,
        },
        "grid_summary": summary,
        "honest_limits": [
            "The grid REPRODUCES the 2026-10-29 sweep exactly (all 6 falsifier rows per point, "
            "all three points) -- the attestation added a self-check without perturbing the "
            "measurement. It is NOT a byte-exact reproduction of the RECORDED 2026-10-07 "
            "baseline: the anchor reads kmer 0.7355 vs recorded 0.7237 (+0.0118) and both "
            "0.7787 vs 0.7706 (+0.0081), because that run's seed is unknown.",
            "NO grid point beats the Markov null in either feature mode. Every cell is "
            "WORSE_THAN or NO_BETTER_THAN_MARKOV_NULL. This artifact attests that the "
            "temperature knob is LIVE and that all three points are plan-eligible; it does NOT "
            "report a generator that works.",
            "Temperature is a real lever on the kmer gap over the null (0.1365 -> 0.0401 -> "
            "0.0733 across T=1.0 -> 0.8 -> 0.6, best at T=0.8) but no setting clears the null.",
            "top_k was NOT swept; it is 0 at every point.",
        ],
    }
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print("wrote %s" % out)
    for r in rows:
        print("  T=%-4s %-28s %s" % (r["temperature"], r["eligibility"],
                                     r["verification_verdict"] or "(nothing to verify)"))
    print("  verified non-default points: %d" % len(verified))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
