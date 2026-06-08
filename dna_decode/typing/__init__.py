"""Shared engine for curated-DB BLAST typing decoders (CGE-family: PlasmidFinder / SerotypeFinder /
ResFinder / VirulenceFinder-style). One generic `call_alleles` blastn-best-hit core; each decoder is a thin
config (DB + header parser + report shape) on top. Reuses the pathotype blastn resolvers. Offline-safe.
"""
