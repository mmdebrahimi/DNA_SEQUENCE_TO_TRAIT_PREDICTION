"""Standalone graphical browser for a genome-map JSON — the deferred v1 "graphical browser".

Pure + dependency-free: `build_genome_map_html(genome_map, ...)` turns the same map dict that
`scripts/genome_map.py` produces (and `render_genome_summary_md` renders as text) into a single
self-contained HTML file — one horizontal track per contig, every feature a block positioned by
`start`/`end` and coloured by `primary_tier`, with the AMR/virulence determinant features surfaced and
labelled. No external JS/CSS libraries, no network, no new dependency; the map JSON is the only input.

THE HONESTY WALL CARRIES INTO THE VISUAL (this is the whole point of the map):
  - Only `determinant-phenotype` (and, separately, `virulence-determinant`) features show a phenotype /
    pathotype claim. Every other tier is a molecular-function / homology / unknown annotation with NO
    phenotype — the browser never colours a phenotype onto them.
  - The `unknown_under_bakta_db_light` rate, `overlay_status`, and `all_joins_symbol_fallback` render as a
    prominent caveat banner, not buried.
  - The whole-contig `region` backbone feature defines each track's length and is NOT drawn as a block
    (it would otherwise cover the entire track — the same feature the coord-join excludes).
"""
from __future__ import annotations

import html
from datetime import date as _date

# tier -> (colour, track-height px, z-order, legend label). Determinant/virulence are tallest + on top.
_TIER_STYLE: dict[str, tuple[str, int, int, str]] = {
    "determinant-phenotype": ("#c0392b", 20, 6, "determinant-phenotype — AMR phenotype claim (behind the wall)"),
    "virulence-determinant": ("#8e44ad", 18, 5, "virulence-determinant — curated VF present (presence only)"),
    "curated-molecular-function": ("#2980b9", 9, 3, "curated-molecular-function"),
    "homology-only-hypothesis": ("#e67e22", 7, 2, "homology-only-hypothesis"),
    "unknown": ("#95a5a6", 6, 1, "unknown (db-light coverage caveat)"),
    "secondary-evidence": ("#bdc3c7", 6, 1, "secondary-evidence"),
}
_DEFAULT_STYLE = ("#7f8c8d", 6, 1, "other")
_BACKBONE_TYPES = {"region"}


def _e(x) -> str:
    return html.escape(str(x if x is not None else ""), quote=True)


def contig_lengths(features: list[dict]) -> dict[str, int]:
    """{seqid: length}. Prefer the whole-contig `region` backbone; else fall back to the max `end`."""
    region: dict[str, int] = {}
    max_end: dict[str, int] = {}
    for f in features:
        sid = f.get("seqid")
        if sid is None:
            continue
        end = int(f.get("end") or 0)
        max_end[sid] = max(max_end.get(sid, 0), end)
        if f.get("raw_feature_type") in _BACKBONE_TYPES:
            region[sid] = max(region.get(sid, 0), end)
    # prefer the region backbone length; fall back to the max feature end for contigs without one
    return {sid: region.get(sid, max_end[sid]) for sid in max_end}


def _feature_label(f: dict) -> str:
    if f.get("raw_gene_symbol"):
        return f["raw_gene_symbol"]
    # for a determinant-phenotype feature with no gene symbol, the curated determinant symbol
    # (e.g. `parC_E84V`) is far more informative than the strain-unique locus tag.
    dets = [p.get("determinant_symbol") for p in (f.get("phenotype") or []) if p.get("determinant_symbol")]
    if dets:
        return ", ".join(dict.fromkeys(dets))   # dedup, preserve order
    return f.get("raw_locus_tag") or (f.get("raw_product") or "")[:60] or f.get("raw_feature_type") or "feature"


def _feature_tooltip(f: dict) -> str:
    parts = [
        f"{_feature_label(f)}",
        f"tier: {f.get('primary_tier')}",
        f"{f.get('seqid')}:{f.get('start')}-{f.get('end')} ({f.get('strand')})",
        f"type: {f.get('raw_feature_type')}",
    ]
    if f.get("raw_product"):
        parts.append(f"product: {f['raw_product']}")
    pheno = f.get("phenotype") or []
    for p in pheno:
        det = p.get("determinant_symbol") or p.get("gene") or ""
        drug = p.get("drug") or p.get("phenotype") or ""
        pred = p.get("genome_prediction") or p.get("phenotype") or ""
        if det or drug:
            parts.append(f"  {det} → {drug} = {pred}")
    return " | ".join(str(x) for x in parts)


def _blocks_for_contig(features: list[dict], seqid: str, length: int) -> str:
    """Absolutely-positioned feature blocks for one contig, scaled to [0, length]."""
    if length <= 0:
        return ""
    rows: list[dict] = [
        f for f in features
        if f.get("seqid") == seqid and f.get("raw_feature_type") not in _BACKBONE_TYPES
        and f.get("start") is not None and f.get("end") is not None
    ]
    rows.sort(key=lambda f: _TIER_STYLE.get(f.get("primary_tier"), _DEFAULT_STYLE)[2])  # tallest last (on top)
    out: list[str] = []
    for f in rows:
        tier = f.get("primary_tier") or "unknown"
        color, h, _z, _lbl = _TIER_STYLE.get(tier, _DEFAULT_STYLE)
        s, e = int(f["start"]), int(f["end"])
        left = max(0.0, min(100.0, 100.0 * s / length))
        width = max(0.06, min(100.0 - left, 100.0 * (e - s) / length))  # min width so a point is visible
        tier_class = "t-" + tier.replace(" ", "_")
        labelled = tier in ("determinant-phenotype", "virulence-determinant")
        label_html = (f'<span class="flabel">{_e(_feature_label(f))}</span>' if labelled else "")
        out.append(
            f'<div class="feat {tier_class}" style="left:{left:.4f}%;width:{width:.4f}%;'
            f'height:{h}px;background:{color};" title="{_e(_feature_tooltip(f))}">{label_html}</div>'
        )
    return "\n".join(out)


def _metrics_header(gm: dict) -> str:
    m = gm.get("metrics", {})
    per_tier = m.get("per_tier_counts", {})
    chips = "".join(
        f'<span class="chip" style="background:{_TIER_STYLE.get(k, _DEFAULT_STYLE)[0]}">{_e(k)}: {_e(v)}</span>'
        for k, v in per_tier.items()
    )
    calls = m.get("genome_level_calls", {}) or {}
    call_str = ", ".join(f"{_e(d)}={_e(v.get('prediction'))}" for d, v in calls.items()) or "(none)"
    jq = m.get("join_quality", {}) or {}
    lines = [
        f'<div class="chips">{chips}</div>',
        f'<p><b>total features:</b> {_e(m.get("total_features"))} · '
        f'<b>determinant-phenotype:</b> {_e(m.get("determinant_phenotype_feature_count"))} · '
        f'<b>organism (-O):</b> <code>{_e(gm.get("amrfinder_organism"))}</code></p>',
        f'<p><b>genome-level R/S calls</b> (separate from features): {call_str}</p>',
        f'<p><b>determinant join quality:</b> n_main_rows={_e(jq.get("n_main_rows"))} '
        f'high_confidence={_e(jq.get("n_high_confidence_join"))} '
        f'symbol_fallback={_e(jq.get("n_symbol_fallback"))} unjoined={_e(jq.get("n_unjoined"))}</p>',
    ]
    vstatus = gm.get("virulence_status")
    if vstatus:
        pc = m.get("genome_pathotype_call") or {}
        dc = pc.get("derived_call") or {}
        pc_str = (f'{_e(dc.get("primary"))} [{_e(dc.get("confidence_tier"))}]'
                  if pc.get("status") == "ok" and dc else _e(pc.get("status")))
        lines.append(f'<p><b>virulence_status:</b> <code>{_e(vstatus)}</code> · '
                     f'<b>virulence-determinant features:</b> {_e(m.get("virulence_determinant_feature_count", 0))} · '
                     f'<b>pathotype call</b> (presence-only): {pc_str}</p>')
    return "\n".join(lines)


def _honesty_banner(gm: dict, gate_result: dict | None) -> str:
    m = gm.get("metrics", {})
    unk = m.get("unknown_under_bakta_db_light")
    overlay = gm.get("overlay_status")
    all_fb = m.get("all_joins_symbol_fallback")
    bits = [
        "<b>Honesty wall.</b> A phenotype/pathotype claim appears ONLY on <span class='k-dp'>"
        "determinant-phenotype</span> and <span class='k-vf'>virulence-determinant</span> features "
        "(behind a validated-determinant coordinate join). Every other tier is annotation with NO "
        "phenotype.",
    ]
    if unk is not None:
        bits.append(f"<code>unknown_under_bakta_db_light = {unk:.3f}</code> — the db-light coverage caveat "
                    "is IN the field name, NOT a biological-unknown claim.")
    if overlay and overlay != "FULL":
        bits.append(f"<b>overlay_status = {_e(overlay)}</b> — no determinant overlay; tiers from the GFF "
                    "only; phenotype claims require a full overlay.")
    if all_fb:
        bits.append("<b>all_joins_symbol_fallback = True</b> — determinant joins are ALL symbol-fallback "
                    "(NO-GO for a confident determinant claim).")
    if gate_result is not None:
        bits.append(f"G1 features: {_e(len(gate_result.get('g1_features', [])))} "
                    f"(demote={_e(gate_result.get('g1_demote_count'))}, "
                    f"surface={_e(gate_result.get('g1_surface_count'))}) · "
                    f"G2 phenotype-wall pass={_e((gate_result.get('g2_spotcheck') or {}).get('pass'))}")
    return '<div class="banner">' + "<br>".join(bits) + "</div>"


def _legend() -> str:
    items = []
    seen = set()
    for tier, (color, _h, _z, label) in sorted(_TIER_STYLE.items(), key=lambda kv: -kv[1][2]):
        if label in seen:
            continue
        seen.add(label)
        tier_class = "t-" + tier.replace(" ", "_")
        items.append(
            f'<label class="leg"><input type="checkbox" checked data-tier="{tier_class}"> '
            f'<span class="sw" style="background:{color}"></span>{_e(label)}</label>'
        )
    return '<div class="legend">' + "".join(items) + "</div>"


_CSS = """
body{font-family:-apple-system,Segoe UI,Roboto,Helvetica,Arial,sans-serif;margin:18px;color:#222;background:#fafafa}
h1{font-size:20px;margin:0 0 2px} .sub{color:#666;font-size:12px;margin:0 0 12px}
.banner{background:#fff7e6;border:1px solid #f0c674;border-radius:6px;padding:10px 12px;font-size:12.5px;margin:10px 0}
.k-dp{color:#c0392b;font-weight:600}.k-vf{color:#8e44ad;font-weight:600}
.chips{margin:6px 0}.chip{display:inline-block;color:#fff;border-radius:10px;padding:2px 9px;font-size:11.5px;margin:2px 4px 2px 0}
.legend{margin:10px 0 14px}.leg{display:inline-block;font-size:12px;margin-right:14px;cursor:pointer;user-select:none}
.sw{display:inline-block;width:11px;height:11px;border-radius:2px;vertical-align:middle;margin:0 4px}
.contig{margin:16px 0}.contig h3{font-size:13px;margin:0 0 4px;color:#444}
.track{position:relative;height:26px;background:#eee;border:1px solid #ddd;border-radius:3px}
.feat{position:absolute;bottom:0;border-radius:1px;box-sizing:border-box;opacity:.9;cursor:pointer}
.feat:hover{outline:1px solid #111;opacity:1;z-index:99}
.flabel{position:absolute;bottom:20px;left:0;white-space:nowrap;font-size:9.5px;color:#111;background:rgba(255,255,255,.75);padding:0 2px;border-radius:2px}
code{background:#eee;padding:0 3px;border-radius:2px;font-size:12px}
"""

_JS = """
document.querySelectorAll('.leg input').forEach(function(cb){
  cb.addEventListener('change',function(){
    document.querySelectorAll('.'+cb.dataset.tier).forEach(function(el){
      el.style.display = cb.checked ? '' : 'none';
    });
  });
});
"""


def build_genome_map_html(genome_map: dict, *, gate_result: dict | None = None,
                          generated: str | None = None) -> str:
    """Render a genome-map dict to a single standalone HTML string. Pure — no I/O, no network."""
    generated = generated or _date.today().isoformat()
    features = genome_map.get("features", []) or []
    acc = genome_map.get("genome_accession")
    lengths = contig_lengths(features)

    tracks: list[str] = []
    for sid in sorted(lengths, key=lambda s: -lengths[s]):  # largest contig first
        length = lengths[sid]
        n = sum(1 for f in features if f.get("seqid") == sid
                and f.get("raw_feature_type") not in _BACKBONE_TYPES)
        tracks.append(
            f'<div class="contig"><h3>{_e(sid)} — {length:,} bp · {n} features</h3>'
            f'<div class="track">{_blocks_for_contig(features, sid, length)}</div></div>'
        )

    return (
        "<!DOCTYPE html><html><head><meta charset='utf-8'>"
        f"<title>Genome map — {_e(acc)}</title><style>{_CSS}</style></head><body>"
        f"<h1>Genome map — {_e(acc)}</h1>"
        f"<p class='sub'>generated {_e(generated)} · rendered by dna_decode.genome_map.browser · "
        "NOT a learned phenotype predictor</p>"
        f"{_honesty_banner(genome_map, gate_result)}"
        f"{_metrics_header(genome_map)}"
        f"{_legend()}"
        f"{''.join(tracks)}"
        f"<script>{_JS}</script></body></html>"
    )
