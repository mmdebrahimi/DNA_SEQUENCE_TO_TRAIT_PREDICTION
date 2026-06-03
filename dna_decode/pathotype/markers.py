"""v0 virulence marker-cluster catalog (EP-4 pathotype compatibility resolver).

The ledger-locked >=23 marker clusters (project_state v0 Output Contract). Each
cluster maps to a list of gene-name PREFIXES as they appear in the VirulenceFinder
E. coli DB headers (`>gene:accession:...`), matched case-insensitively by startswith.
This is the single source of truth for both detection (which alleles belong to a
cluster) and resolution (which clusters drive which pathotype call).

Reframing per ledger v5 (2026-05-30): v0 is an auditable marker-based pathotype
COMPATIBILITY resolver with abstention, NOT a predictor. Supported (H1-passing)
classes = ExPEC / EPEC / ETEC; EAEC / commensal / clean-EHEC are a documented
scope-limit (the resolver still reports their modules but flags low external validity).
"""
from __future__ import annotations

# cluster -> list of lowercased gene-name prefixes in the VirulenceFinder DB.
CLUSTER_MARKERS: dict[str, list[str]] = {
    # --- STEC / EHEC ---
    "STX1": ["stx1"],
    "STX2": ["stx2"],
    # --- attaching-effacing (EPEC/EHEC) ---
    "LEE": ["eae", "tir", "espa", "espb", "espf"],   # eae is the anchor; esp/tir support
    "BFP_EAF": ["bfpa", "bfpb", "pera"],             # typical EPEC bundle-forming pilus / EAF
    # --- ETEC enterotoxins ---
    "LT": ["elti"],                                  # eltIAB / eltIIAB / eltI-
    "ST": ["estap", "estah", "estb", "sta", "stb"],  # STh / STp / STb
    # --- EAEC modules ---
    "EAEC_REG": ["aggr"],
    "EAEC_TRANSPORT": ["aata"],
    "EAEC_T6SS": ["aaic"],
    "AAF_I": ["agga"],
    "AAF_II": ["aafa"],
    "AAF_III": ["agg3a"],
    "AAF_IV": ["agg4a"],
    "AAF_V": ["agg5a"],
    # --- ExPEC / UPEC compatible ---
    "P_FIMBRIAE": ["papc", "papg", "papa"],
    "S_FIMBRIAE": ["sfa", "foc"],
    "AFA_DRA": ["afa", "dra"],
    "HEMOLYSIN": ["hlya"],
    "CNF1": ["cnf"],
    "SIDEROPHORES": ["iuta", "iron", "fyua", "irp2", "irea", "chua", "sita"],
    "CAPSULE_SERUM": ["kpsmii", "kpsm", "iss", "trat"],
    # --- out-of-scope flags ---
    "EIEC_FLAG": ["ipah"],
    "DAEC_FLAG": ["daad", "daae"],
}

# Primary diarrheagenic E. coli (DEC) modules that define a pathotype call.
PRIMARY_DEC_MODULES = ["STX1", "STX2", "LEE", "BFP_EAF", "LT", "ST",
                       "EAEC_REG", "EAEC_TRANSPORT", "EAEC_T6SS",
                       "AAF_I", "AAF_II", "AAF_III", "AAF_IV", "AAF_V"]

# Strong ExPEC markers (>=2 => UPEC/ExPEC-compatible). SIDEROPHORES/CAPSULE are support-only.
EXPEC_STRONG = ["P_FIMBRIAE", "S_FIMBRIAE", "AFA_DRA", "HEMOLYSIN", "CNF1"]
EXPEC_SUPPORT = ["SIDEROPHORES", "CAPSULE_SERUM"]

# Supported (externally-valid) v0 surface vs documented scope-limit (abstain-leaning).
SUPPORTED_CLASSES = {"tEPEC_COMPATIBLE", "aEPEC_COMPATIBLE", "ETEC_COMPATIBLE",
                     "UPEC_COMPATIBLE"}  # ExPEC=UPEC-compatible; EPEC; ETEC
SCOPE_LIMITED_CLASSES = {"EHEC_COMPATIBLE", "STEC_NON_LEE", "EAEC_COMPATIBLE",
                         "COMMENSAL_LOW_MARKER_BURDEN"}

RULES_VERSION = "pathotype-rules-v0.1.0"
