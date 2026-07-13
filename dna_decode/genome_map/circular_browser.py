"""Standalone CIRCULAR genome-ring browser for a genome-map JSON (2026-07-13).

Pure + dependency-free: `build_circular_html(genome_map, ...)` renders the SAME map dict that
`scripts/genome_map.py` produces (and the LINEAR `browser.py` renders as horizontal tracks) as a
Circos/Proksee-style circular ring — one arc per contig laid end-to-end around the circle, features drawn
as coloured arcs in concentric strand bands, and every AMR/virulence determinant surfaced as a radial
tick + curated-symbol label pointing outward. The most direct "which parts of the DNA are doing what"
view. No external JS/CSS libraries, no CDN, no network; the map JSON is the only input.

Reuses `browser.py`'s `_TIER_STYLE` + `contig_lengths` (single source of truth for tier colours + the
whole-contig `region` backbone handling) so the linear + circular views never drift.

THE HONESTY WALL CARRIES INTO THE VISUAL (identical to the linear browser):
  - only `determinant-phenotype` + `virulence-determinant` features carry a phenotype/pathotype claim
    (behind the coordinate-join wall); every other tier is annotation with NO phenotype.
  - `unknown_under_bakta_db_light`, `overlay_status`, `all_joins_symbol_fallback` render as a banner.
  - the whole-contig `region` backbone defines each arc's span and is NOT drawn as a feature.
"""
from __future__ import annotations

import html
import math
from datetime import date as _date

from .browser import _BACKBONE_TYPES, _DEFAULT_STYLE, _TIER_STYLE, _feature_label, contig_lengths

_TWO_PI = 2 * math.pi


def _e(x) -> str:
    return html.escape(str(x if x is not None else ""), quote=True)


def _pt(cx: float, cy: float, r: float, frac: float) -> tuple[float, float]:
    """Point on the circle at fractional position `frac` in [0,1), 0 = top, clockwise."""
    ang = _TWO_PI * frac - math.pi / 2
    return (cx + r * math.cos(ang), cy + r * math.sin(ang))


def _arc(cx: float, cy: float, ri: float, ro: float, f0: float, f1: float) -> str:
    """SVG path for an annulus segment (ring band) between fracs f0..f1, radii ri..ro."""
    if f1 <= f0:
        f1 = f0 + 1e-6
    large = 1 if (f1 - f0) > 0.5 else 0
    x0o, y0o = _pt(cx, cy, ro, f0)
    x1o, y1o = _pt(cx, cy, ro, f1)
    x1i, y1i = _pt(cx, cy, ri, f1)
    x0i, y0i = _pt(cx, cy, ri, f0)
    return (f"M{x0o:.2f},{y0o:.2f} A{ro:.2f},{ro:.2f} 0 {large} 1 {x1o:.2f},{y1o:.2f} "
            f"L{x1i:.2f},{y1i:.2f} A{ri:.2f},{ri:.2f} 0 {large} 0 {x0i:.2f},{y0i:.2f} Z")


def _layout(features: list[dict], lengths: dict[str, int], gap_frac: float = 0.012):
    """Assign each contig an angular span; return (contig_spans, total_len, drawable_features).

    contig_spans[seqid] = (start_frac, span_frac, length). Positions within a contig map linearly onto
    its span. A small `gap_frac` separates contigs so boundaries are visible.
    """
    order = sorted(lengths, key=lambda s: -lengths[s])
    total = sum(lengths.values()) or 1
    n = len(order)
    usable = 1.0 - gap_frac * n
    spans: dict[str, tuple[float, float, int]] = {}
    cur = 0.0
    for sid in order:
        span = usable * (lengths[sid] / total)
        spans[sid] = (cur, span, lengths[sid])
        cur += span + gap_frac
    return spans, total, order


def _feature_frac(spans, seqid, pos) -> float | None:
    sp = spans.get(seqid)
    if sp is None or sp[2] <= 0:
        return None
    start_frac, span, length = sp
    return start_frac + span * max(0.0, min(1.0, pos / length))


def build_circular_html(genome_map: dict, *, gate_result: dict | None = None,
                        generated: str | None = None, size: int = 820) -> str:
    """Render a genome-map dict to one standalone circular-ring HTML string. Pure — no I/O, no network."""
    generated = generated or _date.today().isoformat()
    features = genome_map.get("features", []) or []
    acc = genome_map.get("genome_accession")
    lengths = contig_lengths(features)
    spans, total, order = _layout(features, lengths)

    cx = cy = size / 2
    R = size / 2 - 60
    # concentric bands: backbone, forward CDS, reverse CDS
    r_back_i, r_back_o = R, R + 12
    r_fwd_i, r_fwd_o = R - 46, R - 6
    r_rev_i, r_rev_o = R - 92, R - 52
    det_tick_o = R + 34

    svg: list[str] = []

    # contig backbone arcs (alternating shade) + labels
    for i, sid in enumerate(order):
        f0, span, length = spans[sid]
        shade = "#34495e" if i % 2 == 0 else "#5d6d7e"
        svg.append(f'<path d="{_arc(cx, cy, r_back_i, r_back_o, f0, f0 + span)}" fill="{shade}"/>')
        mid = f0 + span / 2
        lx, ly = _pt(cx, cy, r_back_o + 16, mid)
        svg.append(f'<text x="{lx:.1f}" y="{ly:.1f}" class="cl" text-anchor="middle">'
                   f'{_e(sid)} · {length:,}bp</text>')

    # feature arcs in strand bands
    dets: list[tuple[float, dict]] = []
    for f in features:
        ftype = f.get("raw_feature_type")
        if ftype in _BACKBONE_TYPES:
            continue
        sid = f.get("seqid")
        s, e = f.get("start"), f.get("end")
        if sid is None or s is None or e is None:
            continue
        f0 = _feature_frac(spans, sid, int(s))
        f1 = _feature_frac(spans, sid, int(e))
        if f0 is None or f1 is None:
            continue
        if f1 - f0 < 3e-4:                       # min angular width so a point stays visible
            f1 = f0 + 3e-4
        tier = f.get("primary_tier") or "unknown"
        color, _h, _z, _lbl = _TIER_STYLE.get(tier, _DEFAULT_STYLE)
        fwd = str(f.get("strand", "+")).startswith("+") or f.get("strand") in (1, "1")
        ri, ro = (r_fwd_i, r_fwd_o) if fwd else (r_rev_i, r_rev_o)
        tier_class = "t-" + tier.replace(" ", "_")
        tip = _e(f"{_feature_label(f)} | {tier} | {sid}:{s}-{e} ({f.get('strand')})")
        svg.append(f'<path class="feat {tier_class}" d="{_arc(cx, cy, ri, ro, f0, f1)}" '
                   f'fill="{color}" data-tip="{tip}"/>')
        if tier in ("determinant-phenotype", "virulence-determinant"):
            dets.append(((f0 + f1) / 2, f))

    # determinant callouts: radial tick + curated-symbol label pointing outward
    for mid, f in sorted(dets, key=lambda t: t[0]):
        color = _TIER_STYLE.get(f.get("primary_tier"), _DEFAULT_STYLE)[0]
        x0, y0 = _pt(cx, cy, r_back_o, mid)
        x1, y1 = _pt(cx, cy, det_tick_o, mid)
        lx, ly = _pt(cx, cy, det_tick_o + 3, mid)
        anchor = "start" if math.cos(_TWO_PI * mid - math.pi / 2) >= 0 else "end"
        svg.append(f'<line x1="{x0:.1f}" y1="{y0:.1f}" x2="{x1:.1f}" y2="{y1:.1f}" '
                   f'stroke="{color}" stroke-width="1.4"/>')
        svg.append(f'<text x="{lx:.1f}" y="{ly:.1f}" class="dl" text-anchor="{anchor}" '
                   f'fill="{color}">{_e(_feature_label(f))}</text>')

    m = genome_map.get("metrics", {}) or {}
    center = (f'<text x="{cx}" y="{cy-8}" class="ct" text-anchor="middle">{_e(acc)}</text>'
              f'<text x="{cx}" y="{cy+12}" class="cts" text-anchor="middle">'
              f'{total:,} bp · {len(order)} contig(s)</text>'
              f'<text x="{cx}" y="{cy+30}" class="cts" text-anchor="middle">'
              f'{_e(m.get("determinant_phenotype_feature_count", len(dets)))} determinant-phenotype</text>')

    from .browser import _honesty_banner, _legend, _metrics_header  # reuse the linear header/banner/legend
    return (
        "<!DOCTYPE html><html><head><meta charset='utf-8'>"
        f"<title>Circular genome map — {_e(acc)}</title><style>{_CSS}</style></head><body>"
        f"<h1>Circular genome map — {_e(acc)}</h1>"
        f"<p class='sub'>generated {_e(generated)} · rendered by dna_decode.genome_map.circular_browser · "
        "NOT a learned phenotype predictor</p>"
        f"{_honesty_banner(genome_map, gate_result)}{_metrics_header(genome_map)}{_legend()}"
        f"<div id='wrap'><svg id='ring' width='{size}' height='{size}'>{''.join(svg)}{center}</svg></div>"
        "<div id='tt' class='tip'></div>"
        f"<script>{_JS}</script></body></html>"
    )


_CSS = """
body{font-family:-apple-system,Segoe UI,Roboto,Helvetica,Arial,sans-serif;margin:16px;color:#222;background:#fafafa}
h1{font-size:20px;margin:0 0 2px}.sub{color:#666;font-size:12px;margin:0 0 10px}
.banner{background:#fff7e6;border:1px solid #f0c674;border-radius:6px;padding:10px 12px;font-size:12.5px;margin:8px 0}
.k-dp{color:#c0392b;font-weight:600}.k-vf{color:#8e44ad;font-weight:600}
.chips{margin:6px 0}.chip{display:inline-block;color:#fff;border-radius:10px;padding:2px 9px;font-size:11.5px;margin:2px 4px 2px 0}
.legend{margin:10px 0 14px}.leg{display:inline-block;font-size:12px;margin-right:14px;cursor:pointer;user-select:none}
.sw{display:inline-block;width:11px;height:11px;border-radius:2px;vertical-align:middle;margin:0 4px}
#wrap{border:1px solid #ddd;border-radius:4px;background:#fff;display:inline-block}
.feat{cursor:pointer}.feat:hover{opacity:1;stroke:#111;stroke-width:0.6}
.cl{font-size:10px;fill:#333}.dl{font-size:9px}.ct{font-size:15px;font-weight:600;fill:#222}.cts{font-size:11px;fill:#666}
.tip{position:fixed;pointer-events:none;background:#111;color:#fff;font-size:11.5px;padding:5px 8px;border-radius:4px;max-width:300px;display:none;z-index:9}
code{background:#eee;padding:0 3px;border-radius:2px;font-size:12px}
"""

_JS = """
const t=document.getElementById('tt');
document.querySelectorAll('.feat').forEach(function(el){
  el.addEventListener('mousemove',function(ev){t.innerHTML=el.dataset.tip;t.style.display='block';
    t.style.left=(ev.clientX+12)+'px';t.style.top=(ev.clientY+12)+'px';});
  el.addEventListener('mouseout',function(){t.style.display='none';});
});
document.querySelectorAll('.leg input').forEach(function(cb){
  cb.addEventListener('change',function(){document.querySelectorAll('.'+cb.dataset.tier).forEach(function(el){
    el.style.display=cb.checked?'':'none';});});
});
"""
