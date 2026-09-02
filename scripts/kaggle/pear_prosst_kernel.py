"""Kaggle T4 kernel: ProSST + ESM2 variant tables for blaCTX-M-14, for the PEAR hybrid.

WHY KAGGLE. ProSST's structure-quantizer needs torch_geometric plus a repo clone, and the local host is
Windows where the quantizer's own docstring records `spawn` as unusable. Kaggle is Linux + GPU + internet,
which removes the toolchain wall rather than working around it. Free tier: no money gate.

WHAT IT EMITS. Two score tables over EVERY single-AA substitution of the mature protein:
  esm2_table   {pos: {aa: log-prob}}          -- masked marginals, ESM2-650M
  prosst_table {"<wt><pos><alt>": log-ratio}  -- ProSST-2048, structure-conditioned
Joining them to PEAR's measured fitness and computing the rank-average hybrid happens LOCALLY; this
kernel does not see the phenotype, so it cannot be tuned to it.

THE NUMBERING TRAP, HANDLED EXPLICITLY. PEAR's reference is the MATURE protein (the gene's 81-nt signal
peptide is trimmed) while AlphaFold models the FULL UniProt sequence. A silent off-by-27 would produce a
structure-conditioned score for the wrong residue at every position -- plausible numbers, entirely wrong.
So the offset is DERIVED by exact substring match of the mature sequence inside the UniProt sequence, the
match is asserted unique, and the run ABORTS if the mature sequence is not found. Nothing is assumed.

Everything is printed so the log itself is the audit trail.
"""
import json
import os
import subprocess
import sys
import traceback

os.environ.setdefault("PYTHONUTF8", "1")
OUT = "/kaggle/working"

# The PEAR reference, mature blaCTX-M-14 (795 nt -> 264 aa). Inlined so the kernel needs no upload.
MATURE_CDS = (
    "GCGCAGACGAGTGCGGTGCAGCAAAAGCTGGCGGCGCTGGAGAAAAGCAGCGGAGGGCGGCTGGGCGTCGCGCTCATCGATACCGCAGATAATACGCAGG"
    "TGCTTTATCGCGGTGATGAACGCTTTCCAATGTGCAGTACCAGTAAAGTTATGGCGGCCGCGGCGGTGCTTAAGCAGAGTGAAACGCAAAAGCAGCTGCT"
    "TAATCAGCCTGTCGAGATCAAGCCTGCCGATCTGGTTAACTACAATCCGATTGCCGAAAAACACGTCAACGGCACAATGACGCTGGCAGAACTGAGCGCG"
    "GCCGCGTTGCAGTACAGCGACAATACCGCCATGAACAAATTGATTGCCCAGCTCGGTGGCCCGGGAGGCGTGACGGCTTTTGCCCGCGCGATCGGCGATG"
    "AGACGTTTCGTCTGGATCGCACTGAACCTACGCTGAATACCGCCATTCCCGGCGACCCGAGAGACACCACCACGCCGCGGGCGATGGCGCAGACGTTGCG"
    "TCAGCTTACGCTGGGTCATGCGCTGGGCGAAACCCAGCGGGCGCAGTTGGTGACGTGGCTCAAAGGCAATACGACCGGCGCAGCCAGCATTCGGGCCGGC"
    "TTACCGACGTCGTGGACTGTGGGTGATAAGACCGGCAGCGGCGACTACGGCACCACCAATGATATTGCGGTGATCTGGCCGCAGGGTCGTGCGCCGCTGG"
    "TTCTGGTGACCTATTTTACCCAGCCGCAACAGAACGCAGAGAGCCGCCGCGATGTGCTGGCTTCAGCGGCGAGAATCATCGCCGAAGGGCTGTAA"
)
AA = "ACDEFGHIKLMNPQRSTVWY"
CODON = {}
for i, b1 in enumerate("TCAG"):
    for b2 in "TCAG":
        for b3 in "TCAG":
            CODON[b1 + b2 + b3] = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"[
                len(CODON)]


def translate(cds):
    return "".join(CODON[cds[i:i + 3]] for i in range(0, len(cds) - len(cds) % 3, 3))


def sh(cmd, check=True):
    print(f"$ {cmd}", flush=True)
    r = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if r.stdout:
        print(r.stdout[-3000:], flush=True)
    if r.returncode != 0:
        print("STDERR:", r.stderr[-3000:], flush=True)
        if check:
            raise SystemExit(f"command failed: {cmd}")
    return r


def resolve_uniprot(mature_seq):
    """Find the UniProt entry whose sequence CONTAINS the mature sequence. Returns (acc, seq, offset).

    offset = 0-based index of the mature sequence inside the full sequence, so
    full_pos = mature_pos + offset  (both 1-based -> full_pos = mature_pos + offset).
    """
    import urllib.parse
    import urllib.request
    q = urllib.parse.quote('protein_name:"beta-lactamase CTX-M-14" OR gene:blaCTX-M-14 OR '
                           '(protein_name:"CTX-M" AND organism_id:562)')
    url = (f"https://rest.uniprot.org/uniprotkb/search?query={q}"
           f"&fields=accession,protein_name,sequence,organism_name&format=json&size=200")
    print("querying UniProt...", flush=True)
    with urllib.request.urlopen(url, timeout=120) as fh:
        data = json.load(fh)
    hits = []
    for r in data.get("results", []):
        seq = (r.get("sequence") or {}).get("value", "")
        if mature_seq in seq:
            hits.append((r["primaryAccession"], seq, seq.index(mature_seq)))
    print(f"UniProt returned {len(data.get('results', []))} entries; "
          f"{len(hits)} CONTAIN the mature sequence exactly", flush=True)
    for acc, seq, off in hits[:10]:
        print(f"   {acc}  len={len(seq)}  offset={off}", flush=True)
    if not hits:
        raise SystemExit("ABORT: no UniProt entry contains the mature sequence. Refusing to guess a "
                         "structure/numbering correspondence.")
    offs = {o for _, _, o in hits}
    if len(offs) != 1:
        raise SystemExit(f"ABORT: candidate entries disagree on the offset ({offs}); ambiguous numbering.")
    # Every candidate agrees on the offset, so ANY of them gives the same mature<->full numbering.
    # That is what makes it safe to fall through to a different accession when AlphaFold has no model
    # for the first one (it 404s on many unreviewed entries).
    return hits


def fetch_alphafold(hits):
    """Try each candidate accession; AlphaFold does not model every UniProt entry."""
    import urllib.request
    last = None
    for acc, seq, off in hits:
        url = f"https://alphafold.ebi.ac.uk/files/AF-{acc}-F1-model_v4.pdb"
        dest = f"{OUT}/AF-{acc}.pdb"
        try:
            urllib.request.urlretrieve(url, dest)
            print(f"AlphaFold model found for {acc}: {os.path.getsize(dest)} bytes", flush=True)
            return acc, seq, off, dest
        except Exception as e:
            print(f"  no model for {acc} ({e})", flush=True)
            last = e
    raise SystemExit(f"ABORT: no AlphaFold model for any of the {len(hits)} candidates ({last}).")


def pdb_sequence(pdb_path):
    three = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
             "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
             "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}
    seen, seq = set(), []
    for line in open(pdb_path):
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            rid = int(line[22:26])
            if rid not in seen:
                seen.add(rid)
                seq.append(three.get(line[17:20].strip(), "X"))
    return "".join(seq)


def main():
    mature = translate(MATURE_CDS).rstrip("*")
    print(f"mature protein: {len(mature)} aa\n{mature}\n", flush=True)

    result = {"mature_len": len(mature)}

    # ---- 1. ESM2-650M masked marginals (GPU) --------------------------------------------------
    # DO NOT upgrade transformers/torch here. Kaggle ships a working GPU stack, and `pip install
    # --upgrade transformers torch` broke it on the first attempt with
    # "ModuleNotFoundError: Could not import module 'EsmConfig'" -- the registry mapping goes stale
    # against the already-loaded package set. Use what the image provides and PRINT the versions so a
    # later drift is visible in the log rather than silent.
    import torch
    import transformers
    print("preinstalled: torch", torch.__version__, "| transformers", transformers.__version__, flush=True)
    from transformers import AutoModelForMaskedLM, AutoTokenizer
    print("torch", torch.__version__, "cuda", torch.cuda.is_available(),
          torch.cuda.get_device_name(0) if torch.cuda.is_available() else "-", flush=True)
    dev = "cuda" if torch.cuda.is_available() else "cpu"
    name = "facebook/esm2_t33_650M_UR50D"
    tok = AutoTokenizer.from_pretrained(name)
    mdl = AutoModelForMaskedLM.from_pretrained(name).to(dev).eval()
    ids = tok(mature, return_tensors="pt")["input_ids"][0]
    aa_ids = tok.convert_tokens_to_ids(list(AA))
    esm = {}
    with torch.no_grad():
        for s in range(0, len(mature), 8):
            chunk = list(range(s + 1, min(s + 9, len(mature) + 1)))
            stack = ids.repeat(len(chunk), 1).clone()
            for r, p in enumerate(chunk):
                stack[r, p] = tok.mask_token_id
            lg = mdl(stack.to(dev)).logits
            for r, p in enumerate(chunk):
                lp = torch.log_softmax(lg[r, p].float(), dim=-1)
                esm[p] = {a: float(lp[i]) for a, i in zip(AA, aa_ids)}
    print(f"ESM2 table: {len(esm)} positions", flush=True)
    result["esm2_positions"] = len(esm)
    json.dump({str(k): v for k, v in esm.items()}, open(f"{OUT}/esm2_ctxm14_mature.json", "w"))

    # FREE ESM2 before loading the next model. Leaving it resident held 14.35 of the T4's 14.56 GiB and
    # made ESMFold OOM on a 16 MiB allocation. Three models run in this kernel; each releases the GPU.
    del mdl
    torch.cuda.empty_cache()
    print(f"freed ESM2; cuda reserved={torch.cuda.memory_reserved() / 2**30:.2f} GiB", flush=True)

    # ---- 2. resolve structure + numbering (ABORTS rather than guessing) ------------------------
    # AlphaFold has NO model for any CTX-M-14 UniProt entry -- all 9 that contain the mature sequence
    # 404, because AlphaFold DB covers reference proteomes and these are plasmid-borne entries from
    # clinical isolates. An experimental PDB is the other option, but crystal structures have missing
    # loops, and ProSST needs one structure token per residue of the sequence it scores -- a gapped
    # chain reintroduces exactly the alignment risk this kernel exists to avoid.
    #
    # So: FOLD THE MATURE SEQUENCE DIRECTLY with ESMFold. The structure is then of precisely the
    # sequence being scored -- offset 0, no gaps, no cross-database numbering to verify. The UniProt
    # lookup still runs because its unanimous offset=27 independently corroborates the signal-peptide
    # trim, but nothing downstream depends on it.
    try:
        hits = resolve_uniprot(mature)
        result["uniprot_candidates"] = [h[0] for h in hits]
        result["uniprot_offset_unanimous"] = hits[0][2]
        print(f"(UniProt corroborates a {hits[0][2]}-residue signal-peptide trim; not used downstream)",
              flush=True)
    except SystemExit as e:
        print(f"UniProt corroboration unavailable ({e}); continuing -- ESMFold needs no accession.",
              flush=True)

    from transformers import EsmForProteinFolding
    print("folding the mature sequence with ESMFold v1...", flush=True)
    fmodel = EsmForProteinFolding.from_pretrained("facebook/esmfold_v1", low_cpu_mem_usage=True)
    fmodel = fmodel.to(dev).eval()
    fmodel.trunk.set_chunk_size(64)                       # T4 memory headroom
    with torch.no_grad():
        pdb_text = fmodel.infer_pdb(mature)
    pdb = f"{OUT}/esmfold_ctxm14_mature.pdb"
    open(pdb, "w").write(pdb_text)
    print(f"  -> {pdb} ({os.path.getsize(pdb)} bytes)", flush=True)
    del fmodel
    torch.cuda.empty_cache()

    pseq = pdb_sequence(pdb)
    print(f"folded CA sequence length {len(pseq)}; matches the scored sequence: {pseq == mature}",
          flush=True)
    if pseq != mature:
        raise SystemExit("ABORT: folded model sequence != the mature sequence being scored.")
    full_seq, offset = mature, 0                          # by construction
    result.update({"structure_source": "ESMFold v1, folded from the mature sequence",
                   "full_len": len(full_seq), "offset": offset, "pdb_matches_scored_seq": True})

    # ---- 3. ProSST structure tokens + variant table --------------------------------------------
    # ProSST-only deps. Deliberately NOT touching torch/transformers (see the note above).
    # torch_scatter is imported by the ProSST GVP encoder but only its scatter ops are used; the wheel
    # rarely builds. Install what does build, and shim torch_scatter in pure python if it is absent.
    sh("pip -q install torch_geometric biotite pathos biopython joblib", check=False)
    sh("pip -q install torch_cluster -f https://data.pyg.org/whl/torch-2.6.0+cu124.html", check=False)
    try:
        import torch_scatter  # noqa: F401
        print("torch_scatter present", flush=True)
    except ImportError:
        import types
        m = types.ModuleType("torch_scatter")

        # The repo imports scatter_mean, scatter_sum AND scatter_max -- all three must exist, and
        # scatter_max returns a (values, argmax) TUPLE in the real package.
        def _expand(src, index, dim):
            idx = index
            for _ in range(src.dim() - index.dim()):
                idx = idx.unsqueeze(-1)
            return idx.expand_as(src)

        def _reduce(src, index, dim, dim_size, how, fill):
            import torch as _t
            n = dim_size if dim_size is not None else int(index.max()) + 1
            shape = list(src.shape)
            shape[dim] = n
            out = _t.full(shape, fill, dtype=src.dtype, device=src.device)
            return out.scatter_reduce(dim, _expand(src, index, dim), src, reduce=how,
                                      include_self=False)

        def scatter_sum(src, index, dim=0, out=None, dim_size=None):
            return _reduce(src, index, dim, dim_size, "sum", 0)

        def scatter_mean(src, index, dim=0, out=None, dim_size=None):
            return _reduce(src, index, dim, dim_size, "mean", 0)

        def scatter_max(src, index, dim=0, out=None, dim_size=None):
            import torch as _t
            vals = _reduce(src, index, dim, dim_size, "amax", float("-inf"))
            return vals, _t.zeros_like(vals, dtype=_t.long)   # argmax unused by the GVP encoder

        m.scatter = m.scatter_add = scatter_sum
        m.scatter_sum, m.scatter_mean, m.scatter_max = scatter_sum, scatter_mean, scatter_max
        sys.modules["torch_scatter"] = m
        print("torch_scatter ABSENT -> pure-python shim installed "
              "(scatter_sum/mean/max)", flush=True)
    sh("git clone -q https://github.com/ai4protein/ProSST.git /kaggle/working/ProSST", check=False)
    sys.path.insert(0, "/kaggle/working/ProSST")
    os.environ["PROSST_REPO"] = "/kaggle/working/ProSST"

    prosst = {}
    try:
        # Mirror dna_decode/forward/prosst_scorer.py::quantize_structure -- the path already validated
        # locally (self-quantized GRB2 == ProteinGym's pre-quantized tokens, 217/217). The entry point is
        # SSTPredictor from prosst.structure.get_sst_seq; `prosst.structure.quantizer.PdbQuantizer` does
        # NOT exist in this repo and was the previous run's ModuleNotFoundError.
        class _SerialPool:
            def __init__(self, *a, **k):
                pass

            def map(self, f, *it):
                return list(map(f, *it))
            imap = imap_unordered = map

            def close(self):
                pass
            join = terminate = restart = clear = close

        import pathos.multiprocessing as _mp
        import pathos.threading as _th
        _mp.Pool = _mp.ProcessPool = _SerialPool           # patch BEFORE the repo binds them at import
        _th.ThreadPool = _SerialPool

        # biotite 1.x renamed filter_backbone -> filter_peptide_backbone; the repo still imports the
        # old name. Alias it back rather than pinning an old biotite (which would fight Kaggle's image).
        import biotite.structure as _bs
        if not hasattr(_bs, "filter_backbone"):
            _bs.filter_backbone = _bs.filter_peptide_backbone
            print("shimmed biotite.structure.filter_backbone -> filter_peptide_backbone", flush=True)

        # PREFLIGHT: report EVERY missing dependency in one run. Discovering them one per push costs
        # ~10 minutes each; this lists the whole set at once. Derived from the repo's own imports.
        import importlib
        missing = []
        for mod in ("Bio", "biotite", "joblib", "numpy", "pandas", "pathos", "scipy",
                    "torch_geometric", "torch_scatter", "tqdm", "torch_cluster"):
            try:
                importlib.import_module(mod)
            except Exception as e:
                missing.append(f"{mod} ({type(e).__name__})")
        print(f"ProSST dependency preflight -- missing: {missing or 'none'}", flush=True)
        result["prosst_missing_deps"] = missing

        from prosst.structure.get_sst_seq import SSTPredictor

        predictor = SSTPredictor(structure_vocab_size=2048, num_processes=0, num_threads=1)
        res = predictor.predict_from_pdb(pdb)
        rec = res[0] if isinstance(res, list) else res
        if isinstance(rec, dict):
            key = next((k for k in rec if str(k).endswith("2048_sst_seq")), None)
            if key is None:
                raise KeyError(f"no '*2048_sst_seq' key in quantizer result (keys: {list(rec)})")
            toks = list(rec[key])
        else:
            toks = list(rec)
        toks = [int(t) for t in toks]
        print(f"structure tokens: {len(toks)} (full protein, expected {len(full_seq)})", flush=True)
        if len(toks) != len(full_seq):
            raise SystemExit(f"ABORT: {len(toks)} structure tokens != {len(full_seq)} residues.")

        from transformers import AutoModelForMaskedLM as AM, AutoTokenizer as AT
        pt = AT.from_pretrained("AI4Protein/ProSST-2048", trust_remote_code=True)
        pm = AM.from_pretrained("AI4Protein/ProSST-2048", trust_remote_code=True).to(dev).eval()
        enc = pt(full_seq, return_tensors="pt")
        st = torch.tensor([[1] + [t + 3 for t in toks] + [2]])
        with torch.no_grad():
            logits = pm(input_ids=enc["input_ids"].to(dev),
                        attention_mask=enc["attention_mask"].to(dev),
                        ss_input_ids=st.to(dev)).logits
        lp = torch.log_softmax(logits[0].float(), dim=-1).cpu()
        vocab = {a: pt.convert_tokens_to_ids(a) for a in AA}
        for i, wt in enumerate(mature, start=1):          # MATURE numbering out
            full_i = i + offset                            # -> FULL numbering for the table lookup
            row = lp[full_i]                               # +1 for <cls> is cancelled by 1-based full_i
            for alt in AA:
                if alt != wt:
                    prosst[f"{wt}{i}{alt}"] = float(row[vocab[alt]] - row[vocab[wt]])
        print(f"ProSST table: {len(prosst)} mutants", flush=True)
        result["prosst_mutants"] = len(prosst)
    except Exception:
        traceback.print_exc()
        result["prosst_error"] = traceback.format_exc()[-2000:]

    json.dump(prosst, open(f"{OUT}/prosst_ctxm14_mature.json", "w"))
    json.dump(result, open(f"{OUT}/pear_prosst_result.json", "w"), indent=2)
    print("\nRESULT:", json.dumps(result, indent=2)[:2000], flush=True)


if __name__ == "__main__":
    main()
