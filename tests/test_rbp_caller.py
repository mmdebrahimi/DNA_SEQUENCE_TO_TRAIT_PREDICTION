"""Offline tests for the RBP-level receptor caller (dna_decode/phage/rbp_caller.py)."""
from __future__ import annotations

from dna_decode.phage import rbp_caller as rc


def test_protein_kmers():
    assert rc.protein_kmers("MKVL", k=4) == frozenset({"MKVL"})
    assert rc.protein_kmers("MKV", k=4) == frozenset()          # too short
    assert rc.protein_kmers("MK-V L", k=2) == rc.protein_kmers("MKVL", k=2)  # non-AA stripped


def test_kmer_similarity_jaccard():
    a = frozenset({"AAAA", "BBBB"}); b = frozenset({"AAAA", "CCCC"})
    assert rc.kmer_similarity(a, b) == 1 / 3       # 1 shared / 3 union
    assert rc.kmer_similarity(a, frozenset()) == 0.0
    assert rc.kmer_similarity(a, a) == 1.0


def test_nearest_rbp_transfers_receptor():
    ref = {"p1": rc.protein_kmers("MKVLAAAAWWWW", 4), "p2": rc.protein_kmers("QQQQYYYYZZZZ", 4)}
    recep = {"p1": "Tsx", "p2": "OmpC"}
    q = rc.protein_kmers("MKVLAAAAWWWW", 4)         # identical to p1
    call = rc.nearest_rbp_receptor(q, ref, recep, exclude=None)
    assert call.status == "CALLED" and call.predicted_receptor == "Tsx" and call.nearest_phage == "p1"


def test_nearest_abstains_below_threshold():
    ref = {"p1": rc.protein_kmers("AAAAAAAA", 4)}
    recep = {"p1": "Tsx"}
    q = rc.protein_kmers("WWWWWWWW", 4)             # no shared k-mers
    call = rc.nearest_rbp_receptor(q, ref, recep, min_similarity=0.05)
    assert call.status == "INDETERMINATE" and call.predicted_receptor is None


def test_leave_one_out_rbp_accounting():
    # 4 phages: two Tsx with near-identical RBPs (recover each other), two OmpC likewise
    km = {
        "t1": rc.protein_kmers("MKVLAAAAWWWWTTTT", 4), "t2": rc.protein_kmers("MKVLAAAAWWWWTTTS", 4),
        "o1": rc.protein_kmers("QQQQYYYYZZZZGGGG", 4), "o2": rc.protein_kmers("QQQQYYYYZZZZGGGH", 4),
    }
    recep = {"t1": "Tsx", "t2": "Tsx", "o1": "OmpC", "o2": "OmpC"}
    res = rc.leave_one_out_rbp(km, recep, min_similarity=0.05)
    assert res.n_total == 4 and res.n_called == 4 and res.n_correct == 4
    assert res.accuracy == 1.0
    assert res.per_receptor["Tsx"] == [2, 2] and res.per_receptor["OmpC"] == [2, 2]


def test_load_committed_rbp_reference():
    # the committed MIT-attributed reference is loadable + self-consistent
    km, rec = rc.load_rbp_reference()
    assert len(km) >= 100 and len(set(rec.values())) >= 10
    from dna_decode.data.phage_receptor import is_receptor_class
    assert all(is_receptor_class(r) for r in rec.values())


def test_call_rbp_from_protein_roundtrip(tmp_path):
    # take one reference RBP protein back out of the .faa and confirm it calls its own receptor
    faa = rc.DEFAULT_RBP_REFERENCE
    label = receptor = None
    seq = []
    with open(faa, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith(";"):
                continue
            if line.startswith(">"):
                if label:
                    break
                label, _, receptor = line[1:].strip().partition("|")
            elif label:
                seq.append(line.strip())
    call = rc.call_rbp_from_protein("".join(seq))
    assert call.status == "CALLED" and call.predicted_receptor == receptor
