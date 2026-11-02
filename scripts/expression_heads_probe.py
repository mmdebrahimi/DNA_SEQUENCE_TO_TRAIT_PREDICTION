"""Feasibility, geometry and circularity probe for the shipped-expression-heads arm.

This probe has more to establish than its siblings, because the arm is NOT a drop-in. It answers
four questions, and the arm should not be run until all four have answers on disk:

1. **Are the heads what the arm assumes?** Every head's `fc1.weight` must be `(64, 11008)` and every
   head's embedded base weights must be byte-identical to the frozen base. A head that has drifted
   turns this from a FROZEN-encoder arm into a FINE-TUNED-encoder arm, which is a different claim
   and a different arm, so a drift is an exit-non-zero refusal rather than a footnote.

2. **Which windowing route should the bar freeze?** `11008 = 128 * 86` locks the heads to 86
   positions against a 150-position substrate, so a head must be run over 65 windows. The cheap
   route convolves once over the full latent and slices after; the true route slices first and
   recomputes per window. They are NOT the same function -- `padding='same'` means a window's edge
   columns see real neighbouring context in one and zero padding in the other. This probe measures
   the deviation AND the cost of both, per head, so Step 4 can choose from a number.

3. **Can the cheap route's error be distinguished from arithmetic noise?** Measured: in float32 the
   two routes differ at every column by order 1e-7 regardless of geometry, because `conv1d` blocks
   a 150-position input differently from an 86-position one. So this probe reports the deviation
   restricted to columns the closed form says are UNCONTAMINATED -- that figure is the noise floor,
   and a bar tolerance below it would be measuring arithmetic rather than geometry.

4. **What does a null here actually bound?** These heads were trained on OTHER libraries, so the arm
   measures TRANSFER to the tile rung and a null bounds transfer rather than bounding the heads.
   Recorded here so the framing is fixed before any score exists.

**Circularity, and it is narrower than it first looks.** The `urtecho2018` head IS circular on the
sigma70 rung (GSE108535) and must never be used as a comparator there. This substrate is GSE144621
peak tiles, and this project's own withdrawn lead established that no genome-wide or tiling Urtecho
head exists -- so the head is not circular HERE. Recorded as a scoped limit precisely so nobody
carries it to the wrong rung.

**The `tslpred-kosuri2013` endpoint does not match this rung, and it is included anyway.** It
predicts TRANSLATION; the other four predict transcription, which is what a promoter-activity rung
asks about. It is declared under the standing instruction to leave no stones unturned, with the
mismatch recorded as a stated limit rather than resolved by exclusion. Including it also costs
interval width: five heads mean an alpha/5 correction instead of alpha/4, which makes the arm
HARDER to pass, not easier.

Read-only. Writes one artifact. Exits non-zero on any refusal.
"""

from __future__ import annotations

import argparse
import json
import os
import platform
import time
from pathlib import Path

import numpy as np

from dna_decode.glm import atlas_heads as AH
from dna_decode.glm.atlas_embed import DEFAULT_PA_REPO, EXPECTED_PARAMS, audit_checkpoint
from scripts.positional_arm_probe import (
    RAM_HEADROOM,
    extrapolate_wall_seconds,
    fits_in_ram,
    free_ram_bytes,
    load_arm_cohort,
    sha256_of,
    stride_sample,
)

REPO = Path(__file__).resolve().parents[1]
OUT = REPO / "wiki/expression_heads_feasibility_2026-11-02.json"

ARM = "shipped-expression-heads-zeroshot-150bp"
CRITERION_HINT = "ceiling"

#: Pinned by the positional arm's own feasibility artifact. The clone lives in a volatile temp
#: directory, so a moved hash means a different checkpoint and every inherited number is void.
EXPECTED_SHA256 = "bb413bc8602f05b6903ee06da9cebecb461990aebfa129c62bb3cb56cd7bf624"
EXPECTED_BYTES = 22_313_814

#: Tiles in the timing sample. Small on purpose: the true route costs ~65 trunk evaluations per tile
#: per head, so a large sample would itself be the expensive run this probe exists to size.
DEFAULT_SAMPLE = 48

SUBSTRATE_ACCESSION = "GSE144621"
CIRCULAR_RUNG_ACCESSION = "GSE108535"


class ProbeRefused(RuntimeError):
    """A precondition the arm rests on is not satisfied."""


# ---------------------------------------------------------------------------------------------------
# pure helpers
# ---------------------------------------------------------------------------------------------------
def route_cost_projection(n_total: int, n_sample: int, elapsed_seconds: float,
                          n_heads: int) -> dict:
    """Project a sampled route timing to the full cohort for every head.

    Reuses the sibling arm's `extrapolate_wall_seconds` rather than restating the arithmetic, then
    multiplies by the head count, because each head runs its own trunk and classifier.
    """
    if int(n_heads) < 1:
        raise ValueError("n_heads must be >= 1")
    per_head = extrapolate_wall_seconds(n_total, n_sample, elapsed_seconds)
    total = per_head * int(n_heads)
    return {"seconds_per_head": round(per_head, 1), "n_heads": int(n_heads),
            "seconds_all_heads": round(total, 1), "hours_all_heads": round(total / 3600.0, 2)}


def uncontaminated_deviation(cheap_cols, true_cols, start: int, length: int,
                             width: int = AH.HEAD_SEQ_LENGTH) -> dict:
    """Deviation restricted to columns the closed form says CANNOT differ -- the noise floor.

    `cheap_cols` / `true_cols` are `(C, width)` pooled windows. Any disagreement here is arithmetic,
    not geometry, so this is the figure a bar's tolerance must sit above.
    """
    bad = AH.contaminated_columns(start, length, width)
    keep = [c for c in range(width) if c not in bad]
    if not keep:
        return {"n_columns": 0, "max_abs_dev": None,
                "why": "every column is contaminated at this offset, so there is no floor to read"}
    a = np.asarray(cheap_cols, dtype=np.float64)[..., keep]
    b = np.asarray(true_cols, dtype=np.float64)[..., keep]
    d = np.abs(a - b)
    return {"n_columns": len(keep), "max_abs_dev": float(d.max()),
            "median_abs_dev": float(np.median(d)),
            "why": ("columns the closed form says are identical; any deviation here is float32 "
                    "conv-blocking noise amplified by max-pool argmax ties, NOT the edge effect")}


def spearman_between(a, b) -> float:
    """Rank agreement of the two routes' reduced per-tile values.

    Uses the SAME `genomewide.spearman` the scorer uses, which maps NaN to 0.0 -- which is exactly
    why `atlas_heads` refuses non-finite input before this is ever reached. A 0.0 arriving here
    means uncorrelated, not broken.
    """
    from dna_decode.glm.genomewide import spearman

    return spearman(list(np.asarray(a, dtype=float)), list(np.asarray(b, dtype=float)))


# ---------------------------------------------------------------------------------------------------
# verification
# ---------------------------------------------------------------------------------------------------
def verify_checkpoint(pa_repo: Path | str | None = None) -> dict:
    """Re-measure the base checkpoint's identity. Refuses on any mismatch.

    The parameter count goes through `model.parameters()` and NOT through the state dictionary. The
    state dictionary sums to 1,847,492; `EXPECTED_PARAMS` is 1,847,428; the 64-element gap is eight
    8-element rotary `inv_freq` BUFFERS. Comparing the state-dict sum to `EXPECTED_PARAMS` reports a
    checkpoint mismatch that does not exist.
    """
    from dna_decode.glm.atlas_embed import load_atlas

    pa = Path(pa_repo if pa_repo is not None else DEFAULT_PA_REPO)
    ck = pa / "trained_weights/base_model/promoteratlas-base.pt"
    if not ck.is_file():
        raise ProbeRefused("base checkpoint absent: %s" % ck)

    size = ck.stat().st_size
    digest = sha256_of(ck)
    if digest != EXPECTED_SHA256 or size != EXPECTED_BYTES:
        raise ProbeRefused(
            "checkpoint MOVED: sha256 %s (expected %s), %d bytes (expected %d). The clone lives in "
            "a volatile temp directory, so every figure inherited from the positional arm is void "
            "until this matches." % (digest[:16], EXPECTED_SHA256[:16], size, EXPECTED_BYTES)
        )

    globals_named = audit_checkpoint(ck)
    model = load_atlas(pa)                      # itself refuses on a parameter-count mismatch
    n_params = sum(p.numel() for p in model.parameters())
    n_state = sum(v.numel() for v in model.state_dict().values())
    if n_params != EXPECTED_PARAMS:
        raise ProbeRefused("%d parameters, expected %d" % (n_params, EXPECTED_PARAMS))

    return {
        "path": str(ck), "bytes": size, "sha256": digest,
        "pickle_globals": globals_named,
        "n_parameters_via_model_parameters": n_params,
        "n_elements_via_state_dict": n_state,
        "buffer_elements_excluded_by_parameters": n_state - n_params,
        "count_trap": ("the state dictionary sums to %d while model.parameters() gives %d; the %d "
                       "element gap is 8 rotary inv_freq BUFFERS at 8 elements each. The repo "
                       "enforces the PARAMETER count, so a probe comparing the state-dict sum to it "
                       "reports a mismatch that does not exist."
                       % (n_state, n_params, n_state - n_params)),
    }, model


def verify_heads(pa_repo: Path | str | None = None) -> tuple[list, dict]:
    """Load every head, re-measuring geometry and base identity. Refuses on drift."""
    pa = Path(pa_repo if pa_repo is not None else DEFAULT_PA_REPO)
    base_sd = AH.load_base_state_dict(pa)
    expr_dir = pa / "trained_weights/expression_prediction"
    heads = AH.load_heads(expr_dir, base_state_dict=base_sd)   # refuses on shape or drift

    rows = {}
    for h in heads:
        ident = AH.base_identity(h.state_dict, base_sd)
        head_only = sum(v.numel() for k, v in h.state_dict.items()
                        if not k.startswith("base_model."))
        rows[h.name] = {
            "path": h.path, "endpoint": h.endpoint,
            "fc1_shape": list(AH.FC1_SHAPE),
            "implied_seq_length": AH.FC1_SHAPE[1] // AH.LATENT_CHANNELS,
            "base_tensors_checked": ident["n_checked"],
            "base_tensors_identical": ident["n_identical"],
            "base_max_abs_dev": ident["max_abs_dev"],
            "n_elements_total": sum(v.numel() for v in h.state_dict.values()),
            "n_elements_head_only": head_only,
            "normalization_params": {k: (float(v) if isinstance(v, (int, float, np.floating))
                                         else v)
                                     for k, v in h.normalization.items()},
            "log_transform": h.log_transform,
            "spearman_is_invariant_to_the_denormalisation": True,
            "denormalisation_note": ("mean/std (and the log transform where present) are MONOTONE, "
                                     "so they cannot move a rank-based verdict. An ABSOLUTE "
                                     "predicted expression value would be wrong without them."),
        }
    return heads, rows


# ---------------------------------------------------------------------------------------------------
# the measurement
# ---------------------------------------------------------------------------------------------------
def measure_routes(head, latents, *, length: int, probe_offset: int) -> dict:
    """Both routes on one head, with their costs and their disagreement."""
    import torch

    sd = head.state_dict
    z = torch.from_numpy(np.asarray(latents, dtype=np.float32))

    t0 = time.time()
    cheap = AH.head_forward_full(z, sd)
    t_cheap = time.time() - t0

    t0 = time.time()
    true_ = AH.head_forward_per_window(z, sd)
    t_true = time.time() - t0

    dev_windows = AH.route_deviation(cheap, true_)
    red = {}
    for how in ("mean", "max"):
        c = AH.reduce_windows(cheap, how).numpy()
        t = AH.reduce_windows(true_, how).numpy()
        red[how] = {"deviation": AH.route_deviation(AH.reduce_windows(cheap, how),
                                                    AH.reduce_windows(true_, how)),
                    "spearman_cheap_vs_true": spearman_between(c, t),
                    "cheap_is_constant": bool(np.std(c) == 0.0),
                    "true_is_constant": bool(np.std(t) == 0.0)}

    pooled_cheap = AH._trunk(z, sd)[:, :, probe_offset:probe_offset + AH.HEAD_SEQ_LENGTH]
    pooled_true = AH._trunk(z[:, :, probe_offset:probe_offset + AH.HEAD_SEQ_LENGTH], sd)
    floor = uncontaminated_deviation(pooled_cheap.detach().numpy(),
                                     pooled_true.detach().numpy(),
                                     probe_offset, length)

    return {
        "n_windows": int(cheap.shape[-1]),
        "seconds_cheap_route": round(t_cheap, 3),
        "seconds_true_route": round(t_true, 3),
        "true_route_cost_multiple": (round(t_true / t_cheap, 1) if t_cheap > 0 else None),
        "per_window_deviation": dev_windows,
        "per_reduction": red,
        "noise_floor_at_offset_%d" % probe_offset: floor,
        "contaminated_columns_at_that_offset": sorted(
            AH.contaminated_columns(probe_offset, length)),
    }


def probe(*, pa_repo: Path | str | None = None, n_sample: int = DEFAULT_SAMPLE,
          probe_offset: int = 30, cohort=None) -> dict:
    t_start = time.time()
    ck, model = verify_checkpoint(pa_repo)
    heads, head_rows = verify_heads(pa_repo)

    tiles, st = cohort if cohort is not None else load_arm_cohort()
    idx = stride_sample(len(tiles), min(n_sample, len(tiles)))
    seqs = [tiles[i].seq for i in idx]
    length = len(seqs[0])
    if length != 150:
        raise ProbeRefused("expected a 150 bp substrate, got %d bp" % length)

    from scripts.positional_arm_extract import make_atlas_encoder

    encode = make_atlas_encoder()
    t0 = time.time()
    latents = encode(seqs)
    t_encode = time.time() - t0
    if latents.shape[1:] != (AH.LATENT_CHANNELS, length):
        raise ProbeRefused("latent is %r, expected (B, %d, %d)"
                           % (latents.shape, AH.LATENT_CHANNELS, length))
    if not np.isfinite(latents).all():
        raise ProbeRefused("the encoder produced non-finite latents")

    per_head = {}
    for h in heads:
        per_head[h.name] = measure_routes(h, latents, length=length, probe_offset=probe_offset)

    n_tiles = len(tiles)
    agg_cheap = sum(v["seconds_cheap_route"] for v in per_head.values())
    agg_true = sum(v["seconds_true_route"] for v in per_head.values())
    projection = {
        "base_forward_pass": {"seconds_sampled": round(t_encode, 2),
                              "seconds_full_cohort": round(
                                  extrapolate_wall_seconds(n_tiles, len(idx), t_encode), 1),
                              "shared_across_heads": True,
                              "why_shared": ("all %d heads carry base weights byte-identical to the "
                                             "frozen base, so one forward pass feeds every head"
                                             % len(heads))},
        "cheap_route": route_cost_projection(n_tiles, len(idx), agg_cheap / len(heads), len(heads)),
        "true_route": route_cost_projection(n_tiles, len(idx), agg_true / len(heads), len(heads)),
    }

    n_windows = AH.n_windows(length)
    need = 8 * n_tiles * len(heads) * 2          # reduced scalars, both reductions, float64
    freeb = free_ram_bytes()
    resources = {
        "n_tiles": n_tiles, "n_windows": n_windows, "n_heads": len(heads),
        "reduced_values_bytes": need,
        "free_ram_bytes": freeb,
        "fits_in_ram": fits_in_ram(need, freeb, headroom=RAM_HEADROOM),
        "headroom_factor": RAM_HEADROOM,
        "ram_is_not_the_binding_constraint": True,
        "binding_constraint": ("WALL CLOCK on the true route, not memory: the arm reduces to one "
                               "scalar per tile per head, so the design is tiny, while the true "
                               "route costs %d trunk evaluations per tile per head." % n_windows),
        "fits_in_ram_semantics": "None means UNMEASURED and must never be read as True",
    }

    verdict, detail = classify_feasibility(per_head, projection, head_rows)

    return {
        "schema": "expression-heads-feasibility-v1",
        "date": "2026-11-02",
        "arm": ARM,
        "criterion_hint": CRITERION_HINT,
        "host": {"platform": platform.platform(), "python": platform.python_version(),
                 "cpu_count": os.cpu_count()},
        "checkpoint": ck,
        "heads": head_rows,
        "n_heads": len(heads),
        "multiplicity_consequence": (
            "The head count IS the arm's Bonferroni divisor. FIVE heads were found, not the four an "
            "earlier note recorded, so the bar must correct at alpha/5. That makes the arm HARDER "
            "to pass than alpha/4, which is the honest direction for an inherited miscount."),
        "substrate": {"n_tiles": n_tiles, "bp": length, "path": st.path,
                      "n_peaks": len({t.peak for t in tiles}),
                      "sample_indices_n": len(idx), "sample_is_strided": True},
        "geometry": {
            "fc1_in_features": AH.FC1_SHAPE[1],
            "arithmetic": "11008 = 128 x 86",
            "head_seq_length": AH.HEAD_SEQ_LENGTH,
            "substrate_length": length,
            "substrate_would_flatten_to": AH.LATENT_CHANNELS * length,
            "direct_feed_is_impossible": True,
            "n_windows": n_windows,
            "window_arithmetic": "%d - %d + 1 = %d" % (length, AH.HEAD_SEQ_LENGTH, n_windows),
            "edge_columns_left": AH.EDGE_COLS_LEFT,
            "edge_columns_right": AH.EDGE_COLS_RIGHT,
            "edge_asymmetry_note": ("the right edge is wider than the left because padding='same' "
                                    "on an even kernel splits (%d, %d), not (%d, %d)"
                                    % (AH.SAME_PAD_LEFT, AH.SAME_PAD_RIGHT,
                                       AH.SAME_PAD_RIGHT, AH.SAME_PAD_LEFT)),
        },
        "routes": per_head,
        "cost_projection": projection,
        "resources": resources,
        "zero_shot_framing": (
            "These heads were trained on OTHER libraries (hossain2020, lafleur2022, urtecho2018, "
            "yu2021, kosuri2013). This arm therefore measures TRANSFER to the GSE144621 tile rung. "
            "A null result bounds that transfer; it does NOT bound the heads."),
        "circularity_scope": {
            "circular_head": "promoteratlas-trspred-urtecho2018",
            "circular_on_rung": CIRCULAR_RUNG_ACCESSION,
            "this_substrate": SUBSTRATE_ACCESSION,
            "is_circular_here": False,
            "why": ("the urtecho2018 head IS circular on the sigma70 rung (%s) and must never be a "
                    "comparator there. This substrate is %s peak tiles, and this project's own "
                    "withdrawn lead established that no genome-wide or tiling Urtecho head exists, "
                    "so it is not circular HERE."
                    % (CIRCULAR_RUNG_ACCESSION, SUBSTRATE_ACCESSION)),
            "recorded_because": "so nobody reuses this head on the rung where it IS circular",
        },
        "endpoint_scope": {
            "transcription_heads": sorted(k for k, v in head_rows.items()
                                          if v["endpoint"] == AH.TRANSCRIPTION),
            "translation_heads": sorted(k for k, v in head_rows.items()
                                        if v["endpoint"] == AH.TRANSLATION),
            "limit": ("promoteratlas-tslpred-kosuri2013 predicts TRANSLATION, which is a different "
                      "quantity from the promoter activity this rung measures, and it uniquely "
                      "carries log_transform=True. It is declared anyway under the standing "
                      "instruction to leave no stones unturned, with the mismatch stated as a limit "
                      "rather than resolved by exclusion."),
        },
        "verdict": verdict,
        "verdict_detail": detail,
        "wall_seconds": round(time.time() - t_start, 1),
    }


def classify_feasibility(per_head: dict, projection: dict, head_rows: dict) -> tuple[str, str]:
    """Closed verdict set. A route choice is NOT made here -- that is the bar's job in Step 4."""
    if not per_head:
        return ("INDETERMINATE_NO_HEAD_MEASURED",
                "no head produced a route measurement, so nothing about the arm is established")
    constant = [n for n, v in per_head.items()
                if v["per_reduction"]["mean"]["cheap_is_constant"]
                or v["per_reduction"]["mean"]["true_is_constant"]]
    if constant:
        return ("INDETERMINATE_COLLAPSED_HEAD",
                "head(s) %r produced constant output over the sample; a collapsed head would score "
                "as a clean null rather than as a fault, so the arm must not run until this is "
                "understood" % sorted(constant))
    if any(v["base_tensors_identical"] != v["base_tensors_checked"] for v in head_rows.values()):
        return ("INDETERMINATE_BASE_DRIFT",
                "a head's embedded base weights are not byte-identical to the frozen base")
    return ("FEASIBLE_ROUTE_CHOICE_IS_THE_BARS",
            "the geometry, the base identity and both routes' cost and deviation are measured for "
            "every head; which route to use is an acceptance-bar decision and is deliberately NOT "
            "made here. Cheap route projects to %s h over all heads, true route to %s h."
            % (projection["cheap_route"]["hours_all_heads"],
               projection["true_route"]["hours_all_heads"]))


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--pa-repo", default=None)
    ap.add_argument("--n-sample", type=int, default=DEFAULT_SAMPLE)
    ap.add_argument("--probe-offset", type=int, default=30)
    ap.add_argument("--out", default=str(OUT))
    args = ap.parse_args(argv)

    try:
        rec = probe(pa_repo=args.pa_repo, n_sample=args.n_sample, probe_offset=args.probe_offset)
    except (ProbeRefused, AH.AtlasHeadsError) as e:
        print("REFUSED (%s): %s" % (type(e).__name__, e), flush=True)
        return 2
    except Exception as e:                                     # noqa: BLE001
        print("PROBE FAILED (%s): %s" % (type(e).__name__, e), flush=True)
        return 1

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(rec, indent=2, sort_keys=True), encoding="utf-8")

    print("HEADS    %d found, all geometry-locked to %d bp, all base-identical"
          % (rec["n_heads"], AH.HEAD_SEQ_LENGTH))
    for name, v in sorted(rec["routes"].items()):
        m = v["per_reduction"]["mean"]
        print("  %-38s cheap %6.2fs  true %7.2fs (%sx)  mean-dev %.3e  rho %.4f"
              % (name, v["seconds_cheap_route"], v["seconds_true_route"],
                 v["true_route_cost_multiple"], m["deviation"]["max_abs_dev"],
                 m["spearman_cheap_vs_true"]))
    print("PROJECT  cheap %s h | true %s h (all heads, %d tiles)"
          % (rec["cost_projection"]["cheap_route"]["hours_all_heads"],
             rec["cost_projection"]["true_route"]["hours_all_heads"],
             rec["substrate"]["n_tiles"]))
    print("VERDICT  %s" % rec["verdict"])
    print("         %s" % rec["verdict_detail"])
    print("wall     %.1f s" % rec["wall_seconds"])
    print("-> %s" % out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
