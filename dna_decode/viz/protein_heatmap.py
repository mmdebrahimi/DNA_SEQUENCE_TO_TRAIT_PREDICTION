"""Standalone DMS-style protein damage heatmap (2026-07-13).

Pure + dependency-free: turns a protein sequence + its ESM masked-marginal log-probs (the SAME cache
`dna_decode.protein_effect.predictor.masked_marginals` writes) into ONE self-contained HTML heatmap —
position (x) x 20 amino acids (y), each cell coloured by the rung-2 `damage_llr` (= logP(wt) - logP(sub);
higher = more damaging), the dms-viz / MaveVis idiom. No external JS/CSS, no CDN, no network.

THE HONESTY RAIL CARRIES INTO THE VISUAL:
  - the colour is a ZERO-SHOT ESM molecular RANK (~0.5 Spearman on ProteinGym stability), NOT a
    probability and NOT a resistance/phenotype call. A banner says so; a supplied known-DRM marker is
    labelled "catalog phenotype" separately so the molecular vs. phenotype distinction stays explicit
    (the K103N lesson: molecularly benign, catalogued as major resistance).
  - the WT residue at each position is outlined (its own damage_llr is 0 by construction).
"""
from __future__ import annotations

import html
from datetime import date as _date

from ..protein_effect.predictor import AA, damage_llr

# diverging palette anchors (benign/tolerated -> damaging). Interpolated in-cell; no library.
_NEG = (43, 108, 176)    # blue  (sub is FAVOURED over wt -> negative llr)
_MID = (245, 245, 245)   # white (~0)
_POS = (176, 42, 55)     # red   (sub is DISFAVOURED vs wt -> positive llr, "damaging")


def _e(x) -> str:
    return html.escape(str(x if x is not None else ""), quote=True)


def _lerp(a, b, t):
    return tuple(round(a[i] + (b[i] - a[i]) * t) for i in range(3))


def _color(llr: float | None, scale: float) -> str:
    if llr is None:
        return "#dddddd"
    t = max(-1.0, min(1.0, llr / scale))
    rgb = _lerp(_MID, _POS, t) if t >= 0 else _lerp(_MID, _NEG, -t)
    return "#%02x%02x%02x" % rgb


def build_matrix(seq: str, logp: dict[int, dict[str, float]]) -> dict:
    """Pure position x AA damage_llr matrix + summary. `logp` keyed by 1-based position.

    cell[pos][aa] = damage_llr(logp, pos, wt, aa) = logP(wt) - logP(aa); None where logp is missing.
    """
    seq = seq.upper()
    cells: dict[int, dict[str, float | None]] = {}
    allv: list[float] = []
    for pos in range(1, len(seq) + 1):
        wt = seq[pos - 1]
        row: dict[str, float | None] = {}
        col = logp.get(pos)
        for aa in AA:
            if col is None or wt not in col or aa not in col:
                row[aa] = None
            elif aa == wt:
                row[aa] = 0.0
            else:
                row[aa] = round(damage_llr(logp, pos, wt, aa), 3)
                allv.append(row[aa])
        cells[pos] = row
    scale = max((abs(v) for v in allv), default=1.0) or 1.0
    return {"seq": seq, "length": len(seq), "cells": cells, "scale": round(scale, 3),
            "n_scored_positions": sum(1 for p in cells if logp.get(p) is not None)}


def _parse_markers(markers: list[str] | None, seq: str) -> list[dict]:
    """['K103N', ...] -> [{pos, wt, mut, label}] for overlay pins (catalog phenotype, shown separately)."""
    out = []
    for mk in markers or []:
        mk = mk.strip().upper()
        if len(mk) < 3 or not mk[1:-1].isdigit():
            continue
        pos = int(mk[1:-1])
        wt_ok = 1 <= pos <= len(seq) and seq[pos - 1].upper() == mk[0]
        out.append({"pos": pos, "wt": mk[0], "mut": mk[-1], "label": mk, "wt_matches_seq": wt_ok})
    return out


_CSS = """
body{font-family:-apple-system,Segoe UI,Roboto,Helvetica,Arial,sans-serif;margin:16px;color:#222;background:#fafafa}
h1{font-size:19px;margin:0 0 2px}.sub{color:#666;font-size:12px;margin:0 0 10px}
.banner{background:#fff7e6;border:1px solid #f0c674;border-radius:6px;padding:10px 12px;font-size:12.5px;margin:8px 0}
.grid{overflow-x:auto;border:1px solid #ddd;border-radius:4px;background:#fff;padding:6px}
table{border-collapse:collapse;font-size:9px}
td.c{width:12px;height:12px;padding:0;border:0.5px solid #f4f4f4}
td.c.wt{outline:1.4px solid #111;outline-offset:-1.4px}
td.c.mk{outline:2px solid #111;outline-offset:-1px}
th.aa{width:14px;text-align:center;color:#333;font-weight:600}
th.pos{font-size:7px;color:#888;height:34px;white-space:nowrap;writing-mode:vertical-rl;transform:rotate(180deg)}
.scalebar{display:inline-block;height:12px;width:150px;background:linear-gradient(90deg,#2b6cb0,#f5f5f5,#b02a37);vertical-align:middle;border:1px solid #ccc}
.tip{position:fixed;pointer-events:none;background:#111;color:#fff;font-size:11.5px;padding:5px 8px;border-radius:4px;display:none;z-index:9}
.mklist{font-size:11.5px;margin:6px 0}.mkpin{display:inline-block;background:#eef;border:1px solid #99c;border-radius:10px;padding:1px 8px;margin:2px 4px 2px 0}
code{background:#eee;padding:0 3px;border-radius:2px;font-size:12px}
"""

_JS = """
const t=document.getElementById('tt');
document.querySelectorAll('td.c').forEach(function(td){
  td.addEventListener('mousemove',function(ev){
    const v=td.dataset.v; t.innerHTML='pos '+td.dataset.p+' '+td.dataset.wt+'&rarr;'+td.dataset.aa+
      '<br>damage_llr = '+(v===''?'n/a (unscored)':v)+(td.classList.contains('wt')?' (wild-type)':'');
    t.style.display='block';t.style.left=(ev.clientX+12)+'px';t.style.top=(ev.clientY+12)+'px';});
  td.addEventListener('mouseout',function(){t.style.display='none';});
});
"""


def build_heatmap_html(seq: str, logp: dict[int, dict[str, float]], *, title: str = "protein",
                       markers: list[str] | None = None, generated: str | None = None,
                       max_positions: int = 400) -> str:
    """Render a position x AA damage-heatmap to one standalone HTML string. Pure — no I/O, no network."""
    generated = generated or _date.today().isoformat()
    mtx = build_matrix(seq, logp)
    seqU = mtx["seq"]
    scale = mtx["scale"]
    pins = _parse_markers(markers, seqU)
    pin_pos = {p["pos"]: p for p in pins}

    n = min(mtx["length"], max_positions)
    truncated = mtx["length"] > n
    # header row of positions
    pos_hdr = "".join(f'<th class="pos">{_e(seqU[i-1])}{i}</th>' for i in range(1, n + 1))
    rows = [f'<tr><th class="aa"></th>{pos_hdr}</tr>']
    for aa in AA:
        tds = [f'<th class="aa">{aa}</th>']
        for i in range(1, n + 1):
            v = mtx["cells"][i][aa]
            wt = seqU[i - 1]
            cls = "c" + (" wt" if aa == wt else "") + (" mk" if (i in pin_pos and pin_pos[i]["mut"] == aa) else "")
            tds.append(
                f'<td class="{cls}" style="background:{_color(v, scale)}" '
                f'data-p="{i}" data-wt="{_e(wt)}" data-aa="{aa}" data-v="{_e("" if v is None else v)}"></td>'
            )
        rows.append("<tr>" + "".join(tds) + "</tr>")

    pins_html = ""
    if pins:
        chips = "".join(
            f'<span class="mkpin">{_e(p["label"])}'
            + ("" if p["wt_matches_seq"] else " ⚠ wt≠seq")
            + f' · llr={_e(mtx["cells"].get(p["pos"], {}).get(p["mut"]))}</span>'
            for p in pins)
        pins_html = (f'<div class="mklist"><b>catalog markers</b> (phenotype from the AMR catalog — a '
                     f'SEPARATE axis from the molecular colour): {chips}</div>')

    banner = (
        '<div class="banner"><b>What the colour means.</b> Each cell is the rung-2 '
        '<code>damage_llr = logP(wt) − logP(substitution)</code> from ESM2-650M — a <b>zero-shot molecular '
        'RANK</b> (~0.5 Spearman on ProteinGym stability), <b>NOT a probability and NOT a resistance call</b>. '
        'Red = the substitution is disfavoured vs. wild-type (more likely disruptive); blue = favoured. '
        'The wild-type residue (outlined) is 0 by construction.<br>'
        '<b>Molecular ≠ phenotype.</b> A molecularly-benign residue can still be a catalogued resistance '
        'mutation (e.g. HIV K103N) — catalog markers are shown on a separate axis, never inferred from colour. '
        'Frozen decoder surface untouched.</div>'
    )
    trunc = (f'<p class="sub">⚠ showing first {n} of {mtx["length"]} positions (max_positions)</p>'
             if truncated else "")
    return (
        "<!DOCTYPE html><html><head><meta charset='utf-8'>"
        f"<title>Damage heatmap — {_e(title)}</title><style>{_CSS}</style></head><body>"
        f"<h1>Protein damage heatmap — {_e(title)}</h1>"
        f"<p class='sub'>generated {_e(generated)} · rendered by dna_decode.viz.protein_heatmap · "
        f"{_e(mtx['length'])} aa · {_e(mtx['n_scored_positions'])} scored positions · "
        "zero-shot ESM molecular rank, NOT a phenotype call</p>"
        f"{banner}{pins_html}"
        f'<p class="sub">tolerated <span class="scalebar"></span> damaging · '
        f'scale ±{_e(scale)} llr</p>{trunc}'
        f'<div class="grid"><table>{"".join(rows)}</table></div>'
        "<div id='tt' class='tip'></div>"
        f"<script>{_JS}</script></body></html>"
    )
