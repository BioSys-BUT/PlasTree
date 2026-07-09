import pytest
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from plastree.kmertree import kmer_tree


def make_seq(seq_str, seq_id="seq"):
    return SeqRecord(Seq(seq_str), id=seq_id, description="")


def test_identical_sequences_distance_zero():
    seq = "ATCGATCGATCGATCGATCGATCGATCG" * 5
    seqs = [make_seq(seq, "A"), make_seq(seq, "B")]
    tree, dist_mat = kmer_tree(seqs)
    assert dist_mat["A", "B"] == pytest.approx(0.0, abs=1e-6)


def test_distance_symmetry():
    seq1 = "ATCGATCGATCGATCGATCG" * 5
    seq2 = "GCTAGCTAGCTAGCTAGCTA" * 5
    seqs = [make_seq(seq1, "A"), make_seq(seq2, "B")]
    tree, dist_mat = kmer_tree(seqs)
    assert dist_mat["A", "B"] == pytest.approx(dist_mat["B", "A"], abs=1e-10)


def test_distance_range():
    seq1 = "ATCGATCGATCG" * 10
    seq2 = "GCTAGCTAGCTA" * 10
    seqs = [make_seq(seq1, "A"), make_seq(seq2, "B")]
    tree, dist_mat = kmer_tree(seqs)
    assert 0.0 <= dist_mat["A", "B"] <= 10.0  # normalizováno, ale závisí na k


def test_correct_leaf_names():
    seq = "ATCGATCG" * 10
    names = ["seq1", "seq2", "seq3"]
    seqs = [make_seq(seq, n) for n in names]
    tree, _ = kmer_tree(seqs)
    leaf_names = {c.name for c in tree.get_terminals()}
    assert leaf_names == set(names)


def test_similar_sequences_closer_than_different():
    # K-mer metoda rozlišuje sekvence podle frekvenčních profilů k-merů.
    # Použijeme deterministické pseudonáhodné sekvence s různou GC kompozicí.
    import random
    rng = random.Random(42)

    def rand_seq(length, weights):
        return "".join(rng.choices("ACGT", weights=weights, k=length))

    # base: rovnoměrná kompozice (25 % každý nukleotid)
    base      = rand_seq(500, [1, 1, 1, 1])
    # similar: téměř stejná kompozice (mírně AT-bohatší)
    similar   = rand_seq(500, [1.2, 0.8, 0.8, 1.2])
    # different: extrémně GC-bohatá (80 % G+C)
    different = rand_seq(500, [0.1, 0.9, 0.9, 0.1])

    seqs = [make_seq(base, "base"), make_seq(similar, "similar"), make_seq(different, "diff")]
    tree, dist_mat = kmer_tree(seqs)
    assert dist_mat["base", "similar"] < dist_mat["base", "diff"]


def test_three_sequences_tree_topology():
    # A a B jsou identické, C je odlišné → A a B by se měly shlukovat dohromady
    seq_ab = "ATCGATCGATCGATCGATCGATCG" * 5
    seq_c  = "GCTAGCTAGCTAGCTAGCTAGCTA" * 5
    seqs = [make_seq(seq_ab, "A"), make_seq(seq_ab, "B"), make_seq(seq_c, "C")]
    tree, dist_mat = kmer_tree(seqs)
    assert dist_mat["A", "B"] < dist_mat["A", "C"]
    assert dist_mat["A", "B"] < dist_mat["B", "C"]