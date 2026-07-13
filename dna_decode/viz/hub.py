"""Unified visual-decoder hub — the semantic-zoom entry point (2026-07-13).

Pure + dependency-free: `build_hub_html(views)` ties the separate visual-decoder views into ONE landing
page organized by the three scales the decoder operates at ("genome -> network -> protein" — the
Google-Maps-for-the-genome semantic zoom, done honestly). No external JS/CSS, no CDN; the view
descriptors are the only input. Each card links a generated standalone HTML view and carries the SAME
honesty framing the view itself does — the hub never makes a claim the underlying view walls off.

A `View` descriptor = {scale, title, href, blurb, honesty}. `SCALES` is the fixed three-tier model.
"""
from __future__ import annotations

import html
from datetime import date as _date

# fixed scale model: id -> (ordinal, title, one-line "what question it answers")
SCALES: dict[str, tuple[int, str, str]] = {
    "genome": (1, "Whole genome — which parts do what",
               "point at one genome -> an evidence-tiered map of every feature (circular ring / linear browser)"),
    "network": (2, "The network — which parts relate to which",
                "associational co-occurrence of AMR determinants, with the leave-one-clade-out de-confound drawn in"),
    "protein": (3, "The protein — what a specific edit does",
                "a position x amino-acid ESM molecular-damage heatmap for one protein (zero-shot rank, not a phenotype)"),
}


def _e(x) -> str:
    return html.escape(str(x if x is not None else ""), quote=True)


# filename-prefix -> (scale, title-stem, blurb, honesty) for auto-discovery. Pure + testable.
_VIEW_PATTERNS: list[tuple[str, str, str, str, str]] = [
    ("network_cooccurrence_", "network", "Co-occurrence network",
     "force-directed AMR-determinant co-occurrence graph",
     "associational, NOT causal; dashed = lineage-mediated (possible clonal structure)"),
    ("circular_genome_", "genome", "Circular genome ring",
     "Circos-style ring; features by strand band, determinants called out",
     "phenotype claim only on determinant-phenotype features"),
    ("genome_map_", "genome", "Linear genome browser",
     "horizontal per-contig feature tracks",
     "phenotype claim only behind the determinant coordinate-join wall"),
    ("heatmap_", "protein", "Protein damage heatmap",
     "position x amino-acid ESM molecular-damage heatmap",
     "zero-shot molecular rank, NOT a resistance call; catalog markers on a separate axis"),
]


def classify_view_file(name: str) -> dict | None:
    """Map a generated view HTML filename -> a View descriptor, or None if unrecognized. Pure."""
    if not name.endswith(".html"):
        return None
    for prefix, scale, title, blurb, honesty in _VIEW_PATTERNS:
        if name.startswith(prefix):
            stem = name[len(prefix):-len(".html")]
            label = stem.replace("_", " ").strip()
            return {"scale": scale, "title": f"{title}" + (f" — {label}" if label else ""),
                    "href": name, "blurb": blurb, "honesty": honesty}
    return None


def _card(v: dict) -> str:
    scale = v.get("scale", "genome")
    ordv = SCALES.get(scale, (0, scale, ""))[0]
    return (
        f'<a class="card s-{_e(scale)}" href="{_e(v.get("href"))}">'
        f'<div class="badge">scale {ordv} · {_e(scale)}</div>'
        f'<div class="ctitle">{_e(v.get("title"))}</div>'
        f'<div class="cblurb">{_e(v.get("blurb"))}</div>'
        f'<div class="chonesty">{_e(v.get("honesty"))}</div></a>'
    )


def build_hub_html(views: list[dict], *, generated: str | None = None, title: str = "DNA decoder — visual") -> str:
    """Render the hub landing page. Pure — no I/O, no network. `views` = list of View descriptors."""
    generated = generated or _date.today().isoformat()
    by_scale: dict[str, list[dict]] = {}
    for v in views:
        by_scale.setdefault(v.get("scale", "genome"), []).append(v)

    sections = []
    for scale, (_ord, stitle, sq) in sorted(SCALES.items(), key=lambda kv: kv[1][0]):
        cards = by_scale.get(scale, [])
        cards_html = "".join(_card(v) for v in cards) if cards else \
            '<div class="empty">— no view generated for this scale yet —</div>'
        sections.append(
            f'<section><h2><span class="snum">{_ord}</span> {_e(stitle)}</h2>'
            f'<p class="sq">{_e(sq)}</p><div class="cards">{cards_html}</div></section>'
        )

    banner = (
        '<div class="banner"><b>What this is — and is not.</b> A legibility layer over the deterministic '
        'decoder + its findings. Every view is READ-ONLY over committed artifacts; the frozen decoder '
        'surface is untouched. The honesty rails carry into the pixels: the network is <b>associational, '
        'not causal</b> (the clade de-confound is drawn as solid vs. dashed); the protein heatmap is a '
        '<b>zero-shot molecular rank, not a resistance call</b> (catalog markers live on a separate axis); '
        'a phenotype claim appears only behind a validated-determinant wall. Nothing here is a learned '
        'phenotype predictor.</div>'
    )
    return (
        "<!DOCTYPE html><html><head><meta charset='utf-8'>"
        f"<title>{_e(title)}</title><style>{_CSS}</style></head><body>"
        f"<h1>{_e(title)}</h1>"
        f"<p class='sub'>generated {_e(generated)} · dna_decode.viz.hub · semantic zoom: genome &rarr; network &rarr; protein</p>"
        f"{banner}{''.join(sections)}"
        "</body></html>"
    )


_CSS = """
body{font-family:-apple-system,Segoe UI,Roboto,Helvetica,Arial,sans-serif;margin:20px auto;max-width:1000px;color:#222;background:#fafafa}
h1{font-size:23px;margin:0 0 2px}.sub{color:#666;font-size:12.5px;margin:0 0 12px}
.banner{background:#fff7e6;border:1px solid #f0c674;border-radius:6px;padding:11px 13px;font-size:12.5px;margin:10px 0 18px}
section{margin:20px 0}h2{font-size:16px;margin:0 0 2px;color:#2c3e50}
.snum{display:inline-block;width:22px;height:22px;line-height:22px;text-align:center;background:#2c3e50;color:#fff;border-radius:50%;font-size:12px;margin-right:6px}
.sq{color:#666;font-size:12.5px;margin:0 0 10px}
.cards{display:flex;flex-wrap:wrap;gap:12px}
.card{display:block;width:300px;text-decoration:none;color:inherit;background:#fff;border:1px solid #ddd;border-left:5px solid #7f8c8d;border-radius:6px;padding:11px 13px;transition:box-shadow .12s}
.card:hover{box-shadow:0 2px 10px rgba(0,0,0,.12)}
.card.s-genome{border-left-color:#c0392b}.card.s-network{border-left-color:#2980b9}.card.s-protein{border-left-color:#8e44ad}
.badge{font-size:10.5px;color:#888;text-transform:uppercase;letter-spacing:.4px}
.ctitle{font-size:14px;font-weight:600;margin:3px 0}.cblurb{font-size:12px;color:#444;margin:2px 0}
.chonesty{font-size:11px;color:#b9770e;margin-top:5px}
.empty{color:#aaa;font-size:12px;font-style:italic}
"""
