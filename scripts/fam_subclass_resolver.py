"""Resolve a gene/allele symbol -> AMRFinder (Class, Subclass) via the deployed fam.tsv.

The KEYSTONE for scoring cohorts whose genotype ships as a gene-presence table of allele NAMES
(e.g. Sci234 'Supplementary data 1': blaCTX-M-15, aac(3)-IIa, ...) rather than a real AMRFinder run.
The frozen decoder rule (`amr_rules.cipro_determinants_from_main`) keys off AMRFinder's Class/Subclass;
this resolver reproduces, per allele, exactly the Class/Subclass the DEPLOYED v4.2.7 AMRFinder DB
(`data/amrfinder_db/2026-03-24.1/fam.tsv`) curates for that gene family — so a synthesized main.tsv
row carries the same Subclass AMRFinder would have emitted, and the frozen rule is what's under test.

fam.tsv shape (TSV, header row): col1 #node_id, col2 parent_node_id, col3 gene_symbol, col16 class,
col17 subclass, col18 family_name. The tree: a node's class/subclass may be BLANK and inherited from
its parent_node_id chain. The reported FAMILY for a hit is col3 (gene_symbol); specific alleles are
distinct node_ids OR not present at all (then the family node_id/gene_symbol is the resolution target).

Resolution for a query symbol Q (case-insensitive):
  1. exact node_id match  -> that node's (class, subclass), inheriting up parents when blank
  2. exact gene_symbol match -> the FIRST node whose gene_symbol==Q (family node), inheriting
  3. longest-prefix match: the longest node_id OR gene_symbol S such that Q==S or Q starts with S at
     a token boundary (next char is '-', digit, or end) -> that node, inheriting
  4. unresolved -> (None, None); caller decides (conservatively: not counted)

This is faithful for the rule's PURPOSE (does Subclass carry GENTAMICIN / CEPHALOSPORIN/CARBAPENEM),
which is family-level curation. Allele-level subclass refinements that fam.tsv encodes only on rare
sub-nodes (e.g. aac(6')-Ib-G = GENTAMICIN vs the aac(6')-Ib family = generic AMINOGLYCOSIDE) resolve to
the FAMILY value — the conservative, defensible call, since the source caller is allele-level at best
and AMRFinder itself reports the family symbol. Documented as a scope-limit in the cohort artifact.
"""
from __future__ import annotations

import csv
from functools import lru_cache
from pathlib import Path

FAM_TSV_DEPLOYED = Path("data/amrfinder_db/2026-03-24.1/fam.tsv")


class FamResolver:
    """Loaded fam.tsv tree with node/gene_symbol indexes + parent-inheritance resolution."""

    def __init__(self, fam_tsv: Path = FAM_TSV_DEPLOYED):
        self.fam_tsv = Path(fam_tsv)
        # node_id -> row dict; gene_symbol(lower) -> first node_id; sorted key lists for prefix match
        self._by_node: dict[str, dict] = {}
        self._by_symbol: dict[str, str] = {}
        with self.fam_tsv.open(encoding="utf-8") as fh:
            rd = csv.DictReader(fh, delimiter="\t")
            for r in rd:
                node = (r.get("#node_id") or "").strip()
                if not node:
                    continue
                self._by_node[node] = {
                    "parent": (r.get("parent_node_id") or "").strip(),
                    "gene_symbol": (r.get("gene_symbol") or "").strip(),
                    "class": (r.get("class") or "").strip(),
                    "subclass": (r.get("subclass") or "").strip(),
                }
                sym = (r.get("gene_symbol") or "").strip()
                if sym and sym != "-" and sym.lower() not in self._by_symbol:
                    self._by_symbol[sym.lower()] = node
        # prefix-match candidates: every node_id + every gene_symbol, longest first
        self._lower_node = {n.lower(): n for n in self._by_node}
        keys = set(self._lower_node) | set(self._by_symbol)
        self._prefix_keys = sorted(keys, key=len, reverse=True)

    def _inherit(self, node: str) -> tuple[str, str]:
        """Walk node->parent until a non-blank (class, subclass) is found. Returns ('','') if none."""
        seen = set()
        cur = node
        cls = sub = ""
        while cur and cur in self._by_node and cur not in seen:
            seen.add(cur)
            row = self._by_node[cur]
            cls = cls or row["class"]
            sub = sub or row["subclass"]
            if cls and sub:
                break
            cur = row["parent"]
        return cls, sub

    @staticmethod
    def _boundary_ok(q: str, s: str) -> bool:
        if q == s:
            return True
        if not q.startswith(s):
            return False
        nxt = q[len(s)]
        return nxt in "-_" or nxt.isdigit() or nxt == "("

    def resolve(self, symbol: str) -> tuple[str | None, str | None, str]:
        """Resolve a gene/allele symbol -> (class, subclass, match_kind).

        match_kind in {node, gene_symbol, prefix, unresolved}. class/subclass are None if unresolved."""
        q = (symbol or "").strip()
        if not q:
            return None, None, "unresolved"
        ql = q.lower()
        if ql in self._lower_node:
            cls, sub = self._inherit(self._lower_node[ql])
            return cls or None, sub or None, "node"
        if ql in self._by_symbol:
            cls, sub = self._inherit(self._by_symbol[ql])
            return cls or None, sub or None, "gene_symbol"
        for k in self._prefix_keys:
            if self._boundary_ok(ql, k):
                node = self._lower_node.get(k) or self._by_symbol.get(k)
                cls, sub = self._inherit(node)
                return cls or None, sub or None, "prefix"
        return None, None, "unresolved"


@lru_cache(maxsize=1)
def deployed_resolver() -> FamResolver:
    return FamResolver(FAM_TSV_DEPLOYED)
