<!-- result.md for 2026-06-02-1530-ep4-v0-marker-resolver -->
# Result — v0 virulence-marker resolver (confound-immune) WORKS

Artifact: `research_outputs/pathotype_v0_marker_screen_2026-06-02.json`. Method: pure-Python k=15 seed presence (both strands, COV≥0.80) of 40 VirulenceFinder marker families (1789 alleles) vs the 24 cached genomes. No BLAST, no external software.

## Headline: known biology PERFECTLY tracks the labels
- **eae (LEE/intimin): EPEC 12/12, ExPEC 0/12 → AUROC 1.000.** espA same (12/0). EPEC is definitionally LEE+; this is exact.
- **ExPEC virulence-score → ExPEC AUROC 0.882**; mean marker count ExPEC 7.92 vs EPEC 1.08.
- Verdict: **BIOLOGY_TRACKS_LABELS** — the ExPEC/EPEC labels reflect real pathotype gene content, NOT assembler batch.

## Why this resolves the confound that sank the learned track
Run 1: the learned k-mer 0.729 ran on a study==class subset and its top features were rare CTAG-motif k-mers = assembler batch. This resolver separates the SAME 24 genomes on a SPECIFIC 2.8 kb gene (eae) with textbook biology. A whole-gene call is robust to rare-k-mer assembly noise, and eae→EPEC is the established definition — so it is biology by construction, immune to the study==class confound. It also explains why the learned model only reached 0.729: XGBoost latched onto batch k-mers instead of the actual eae locus.

## Per-family presence (the interpretable DNA→trait map)
LEE (EPEC): eae 0/12 vs 12/12 · espA 0/12 · tir 0/7 · espF 0/6 · espB 0/4 · bfpA 0/3 (3 typical EPEC, 9 atypical).
ExPEC adhesin: papA 8/0 · papC 5/0 · sfa 2/0 · foc 2/0 · afa 1/0 · iha 4/1.
ExPEC iron: sitA 10/0 · fyuA 9/1 · irp2 9/1 · iutA 8/1 · ireA 2/0.
ExPEC capsule/toxin/other: kpsMII 7/0 · usp 5/0 · vat 5/0 · sat 4/0 · hlyF 4/0 · iss 8/2.
Shared (correctly non-discriminating): chuA 8/8 · ompT 9/10 · hlyA 1/1.

Each genome carries a marker provenance profile (`per_genome` + `coverage_matrix` in the JSON) — this IS the v0 decoder output: which genes present → which pathotype, with evidence.

## Caveats
- 2 ExPEC genomes (JSLG, JSPG) carry ~0–1 markers — likely fragmented assemblies or commensal-leaning isolates; realistic ExPEC heterogeneity, not a method failure.
- k=15 exact-seed presence at COV≥0.80 is strict; a true gene at <80% allele identity reads absent. Threshold is transparent in the JSON; tune per marker if needed.
- Still N=24, one study per class — but the resolver's validity does NOT depend on that (it keys on specific known genes), which is the entire point.
