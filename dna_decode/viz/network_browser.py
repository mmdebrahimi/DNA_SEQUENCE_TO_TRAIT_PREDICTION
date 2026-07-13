"""Standalone force-directed browser for a co-occurrence graph (2026-07-13).

Pure + dependency-free: `build_network_html(graph, ...)` turns the {nodes, edges} dict that
`dna_decode.viz.network_adapter.build_graph` produces into ONE self-contained HTML file — a
force-directed network drawn on inline SVG with a tiny vanilla-JS spring layout. No external JS/CSS
libraries, no CDN, no network; the graph dict is the only input (opens offline).

THE DE-CONFOUND CARRIES INTO THE PIXELS (this is the whole point):
  - node BORDER encodes the leave-one-clade-out lineage status:
      solid  = generalizes-beyond-lineage (survives clade-grouped CV -> de-confounded association)
      dashed = lineage-mediated (clade-concentrated -> the link may be CLONAL population structure)
      dotted = untested (no cross-axis entry)
  - an EDGE is drawn DASHED iff either endpoint is lineage-mediated -> a possibly-clonal co-occurrence
    never renders as a solid line.
  - a persistent banner states: associational co-occurrence, NOT causal; dashed = clonal-structure
    caveat; organism-scoped; frozen decoder surface untouched. No edge is ever an arrow.
"""
from __future__ import annotations

import html
import json
from datetime import date as _date

from .network_adapter import GENERALIZES, LINEAGE_MEDIATED, UNTESTED

# lineage_status -> (border-css, legend label)
_BORDER = {
    GENERALIZES: ("solid", "generalizes beyond lineage (de-confounded — survives clade-grouped CV)"),
    LINEAGE_MEDIATED: ("dashed", "lineage-mediated (clade-concentrated — may be clonal structure)"),
    UNTESTED: ("dotted", "untested (no cross-axis leave-one-clade-out entry)"),
}
# node fill by gene family prefix (legible grouping; NOT a claim — colour is a grouping aid)
_FAMILY_COLORS = {
    "gyrA": "#c0392b", "gyrB": "#c0392b", "parC": "#e74c3c", "parE": "#e74c3c",  # QRDR (cipro)
    "bla": "#2980b9", "CTX": "#2980b9", "TEM": "#2980b9", "SHV": "#2980b9",       # beta-lactam
    "aac": "#27ae60", "aad": "#27ae60", "aph": "#27ae60", "ant": "#27ae60",       # aminoglycoside
    "sul": "#8e44ad", "dfr": "#8e44ad",                                            # sulfa/trimethoprim
    "tet": "#e67e22", "erm": "#e67e22", "mph": "#e67e22",                          # tet / macrolide
}
_DEFAULT_COLOR = "#7f8c8d"


def _e(x) -> str:
    return html.escape(str(x if x is not None else ""), quote=True)


def _family_color(det: str) -> str:
    for pref, col in _FAMILY_COLORS.items():
        if det.startswith(pref):
            return col
    return _DEFAULT_COLOR


def _prep_nodes(nodes: list[dict]) -> list[dict]:
    """Attach render attributes (colour, radius, border) — pure, deterministic."""
    prevs = [n.get("prevalence", 0) or 0 for n in nodes] or [0]
    pmax = max(prevs) or 1
    out = []
    for n in nodes:
        prev = n.get("prevalence", 0) or 0
        out.append({
            "id": n["id"],
            "color": _family_color(n["id"]),
            "r": round(6 + 14 * (prev / pmax) ** 0.5, 2),   # sqrt scale so big prevalences don't dominate
            "border": _BORDER.get(n.get("lineage_status", UNTESTED), _BORDER[UNTESTED])[0],
            "lineage_status": n.get("lineage_status", UNTESTED),
            "prevalence": prev,
            "auc": n.get("auc"),
            "linked": bool(n.get("linked", False)),
        })
    return out


def _banner(meta: dict) -> str:
    bits = [
        "<b>What this is.</b> An <b>associational</b> co-occurrence network of AMR determinants within "
        f"<code>{_e(meta.get('organism'))}</code> ({_e(meta.get('n_genomes'))} genomes). An edge means two "
        "determinants tend to appear in the same genome — <b>NOT</b> that one causes the other. No edge is a "
        "causal arrow.",
        "<b>The de-confound is in the drawing.</b> A <b>solid</b> node border = the association survives "
        "leave-one-clade-out CV (de-confounded). A <b>dashed</b> border / dashed edge = clade-concentrated: "
        "the co-occurrence may just be <b>clonal population structure</b>, not a real linkage. Dotted = untested.",
        f"<b>Scope + provenance.</b> verdict <code>{_e(meta.get('verdict'))}</code>; edges pruned to "
        f"cooc &ge; {_e(meta.get('min_cooc'))}; read-only over committed artifacts "
        f"(<code>{_e(meta.get('cooc_artifact'))}</code>"
        + (f", <code>{_e(meta.get('crossaxis_artifact'))}</code>" if meta.get("crossaxis_artifact") else "")
        + "). The frozen decoder surface is untouched.",
    ]
    for c in (meta.get("honest_caveats") or [])[:3]:
        bits.append(f"<span class='cav'>caveat:</span> {_e(c)}")
    return '<div class="banner">' + "<br>".join(bits) + "</div>"


def _legend() -> str:
    border_items = "".join(
        f'<span class="leg"><span class="swb" style="border-style:{css}"></span>{_e(lbl)}</span>'
        for css, lbl in _BORDER.values()
    )
    fams = {"QRDR (cipro)": "#c0392b", "beta-lactam": "#2980b9", "aminoglycoside": "#27ae60",
            "sulfa/trimethoprim": "#8e44ad", "tet/macrolide": "#e67e22", "other": _DEFAULT_COLOR}
    fam_items = "".join(
        f'<span class="leg"><span class="swf" style="background:{col}"></span>{_e(name)}</span>'
        for name, col in fams.items()
    )
    return (f'<div class="legend"><div class="legrow"><b>node border (de-confound):</b> {border_items}</div>'
            f'<div class="legrow"><b>fill (gene family):</b> {fam_items}</div>'
            '<div class="legrow"><b>node size</b> = prevalence (&radic; scale) · '
            '<b>edge width</b> = co-occurrence count</div></div>')


_CSS = """
body{font-family:-apple-system,Segoe UI,Roboto,Helvetica,Arial,sans-serif;margin:16px;color:#222;background:#fafafa}
h1{font-size:20px;margin:0 0 2px}.sub{color:#666;font-size:12px;margin:0 0 10px}
.banner{background:#fff7e6;border:1px solid #f0c674;border-radius:6px;padding:10px 12px;font-size:12.5px;margin:8px 0}
.cav{color:#b9770e;font-weight:600}
.legend{font-size:11.5px;margin:8px 0}.legrow{margin:3px 0}.leg{display:inline-block;margin-right:16px}
.swb{display:inline-block;width:13px;height:13px;border-radius:50%;border:2px solid #333;vertical-align:middle;margin:0 4px}
.swf{display:inline-block;width:12px;height:12px;border-radius:2px;vertical-align:middle;margin:0 4px}
#wrap{border:1px solid #ddd;border-radius:4px;background:#fff}
.tip{position:fixed;pointer-events:none;background:#111;color:#fff;font-size:11.5px;padding:5px 8px;border-radius:4px;max-width:280px;display:none;z-index:9}
text.nl{font-size:9px;fill:#111;pointer-events:none}
"""


def build_network_html(graph: dict, *, generated: str | None = None,
                       width: int = 1100, height: int = 720) -> str:
    """Render a {nodes, edges} graph dict to one standalone HTML string. Pure — no I/O, no network."""
    generated = generated or _date.today().isoformat()
    meta = graph.get("meta", {}) or {}
    rnodes = _prep_nodes(graph.get("nodes", []) or [])
    edges = graph.get("edges", []) or []
    coocs = [e.get("cooc", 0) or 0 for e in edges] or [1]
    cmax = max(coocs) or 1
    redges = [{
        "source": e["source"], "target": e["target"], "cooc": e.get("cooc", 0) or 0,
        "w": round(0.6 + 4.0 * ((e.get("cooc", 0) or 0) / cmax), 2),
        "dashed": bool(e.get("lineage_mediated")),
        "lift": e.get("lift"),
    } for e in edges]

    data = json.dumps({"nodes": rnodes, "edges": redges, "w": width, "h": height})
    # inline force layout (deterministic seed via index-based initial positions; no Math.random)
    js = """
const G = %s;
const svg = document.getElementById('net'), W = G.w, H = G.h;
const NS='http://www.w3.org/2000/svg';
const idx = {}; G.nodes.forEach((n,i)=>{idx[n.id]=i; n.x = W/2 + Math.cos(i*2.399)*Math.min(W,H)*0.34;
  n.y = H/2 + Math.sin(i*2.399)*Math.min(W,H)*0.34; n.vx=0; n.vy=0;});
const E = G.edges.map(e=>({...e, s:idx[e.source], t:idx[e.target]})).filter(e=>e.s!=null&&e.t!=null);
// spring-electrical, fixed iteration count (deterministic, no RNG)
for(let it=0; it<300; it++){
  const k = 0.045*(1 - it/400);
  for(let i=0;i<G.nodes.length;i++){for(let j=i+1;j<G.nodes.length;j++){
    const a=G.nodes[i],b=G.nodes[j]; let dx=a.x-b.x,dy=a.y-b.y; let d2=dx*dx+dy*dy+0.01;
    const rep=1200/d2; const d=Math.sqrt(d2); dx/=d;dy/=d;
    a.vx+=dx*rep;a.vy+=dy*rep;b.vx-=dx*rep;b.vy-=dy*rep;}}
  E.forEach(e=>{const a=G.nodes[e.s],b=G.nodes[e.t]; let dx=b.x-a.x,dy=b.y-a.y; const d=Math.sqrt(dx*dx+dy*dy)+0.01;
    const f=(d-90)*0.01*(1+e.cooc/40); dx/=d;dy/=d; a.vx+=dx*f;a.vy+=dy*f;b.vx-=dx*f;b.vy-=dy*f;});
  G.nodes.forEach(n=>{n.x+=Math.max(-8,Math.min(8,n.vx*k*10));n.y+=Math.max(-8,Math.min(8,n.vy*k*10));
    n.x=Math.max(24,Math.min(W-24,n.x));n.y=Math.max(24,Math.min(H-24,n.y));n.vx*=0.85;n.vy*=0.85;});
}
function mk(t,a){const el=document.createElementNS(NS,t);for(const k in a)el.setAttribute(k,a[k]);return el;}
E.forEach(e=>{const a=G.nodes[e.s],b=G.nodes[e.t];
  const ln=mk('line',{x1:a.x,y1:a.y,x2:b.x,y2:b.y,stroke:'#9aa',
    'stroke-width':e.w,'stroke-opacity':0.55,'stroke-dasharray':e.dashed?'5,4':''});
  ln.addEventListener('mousemove',ev=>tip(ev, e.source+' &harr; '+e.target+'<br>cooc='+e.cooc+' lift='+e.lift+
    (e.dashed?'<br><b>dashed:</b> lineage-mediated endpoint (possible clonal structure)':'')));
  ln.addEventListener('mouseout',hide); svg.appendChild(ln);});
const t=document.getElementById('tt');
function tip(ev,htmlv){t.innerHTML=htmlv;t.style.display='block';t.style.left=(ev.clientX+12)+'px';t.style.top=(ev.clientY+12)+'px';}
function hide(){t.style.display='none';}
G.nodes.forEach(n=>{const c=mk('circle',{cx:n.x,cy:n.y,r:n.r,fill:n.color,'fill-opacity':0.85,
    stroke:'#222','stroke-width':2,'stroke-dasharray':n.border==='dashed'?'4,3':(n.border==='dotted'?'1,3':'')});
  c.addEventListener('mousemove',ev=>tip(ev,'<b>'+n.id+'</b><br>prevalence='+n.prevalence+
    '<br>lineage: '+n.lineage_status+'<br>impute AUC='+(n.auc==null?'n/a':n.auc)+' linked='+n.linked));
  c.addEventListener('mouseout',hide); svg.appendChild(c);
  if(n.r>=11){const lbl=mk('text',{x:n.x+n.r+1,y:n.y+3,class:'nl'});lbl.textContent=n.id;svg.appendChild(lbl);}
});
""" % data

    return (
        "<!DOCTYPE html><html><head><meta charset='utf-8'>"
        f"<title>Co-occurrence network — {_e(meta.get('organism'))}</title><style>{_CSS}</style></head><body>"
        f"<h1>AMR determinant co-occurrence network — {_e(meta.get('organism'))}</h1>"
        f"<p class='sub'>generated {_e(generated)} · rendered by dna_decode.viz.network_browser · "
        f"{_e(meta.get('n_nodes'))} nodes / {_e(meta.get('n_edges'))} edges · "
        "associational, NOT causal · frozen decoder surface untouched</p>"
        f"{_banner(meta)}{_legend()}"
        f"<div id='wrap'><svg id='net' width='{width}' height='{height}'></svg></div>"
        "<div id='tt' class='tip'></div>"
        f"<script>{js}</script></body></html>"
    )
