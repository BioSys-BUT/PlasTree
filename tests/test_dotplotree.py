import pytest
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from plastree.dotplotree import _find_long_segments, _dot_mat_distance_core, dotplot_tree


def make_seq(seq_str, seq_id="seq"):
    return SeqRecord(Seq(seq_str), id=seq_id, description="")


# --- _find_long_segments ---

def test_find_long_segments_empty():
    assert _find_long_segments([], window=5, gap=10) == []


def test_find_long_segments_perfect_diagonal():
    # Bodová diagonála: (0,0),(1,1),(2,2),...,(29,29)
    points = [(i, i) for i in range(30)]
    segs = _find_long_segments(points, window=5, gap=10)
    assert len(segs) == 1
    assert len(segs[0]) == 30


def test_find_long_segments_too_short():
    # Segment kratší než window → nemá se vrátit
    points = [(i, i) for i in range(4)]
    segs = _find_long_segments(points, window=5, gap=10)
    assert segs == []


def test_find_long_segments_two_diagonals_merged():
    # Dvě části stejné diagonály s mezerou menší než gap → merge
    points = [(i, i) for i in range(20)] + [(i, i) for i in range(25, 45)]
    segs = _find_long_segments(points, window=5, gap=10)
    assert len(segs) == 1


def test_find_long_segments_two_diagonals_not_merged():
    # Dvě části stejné diagonály s mezerou větší než gap → zůstanou oddělené
    points = [(i, i) for i in range(20)] + [(i, i) for i in range(100, 130)]
    segs = _find_long_segments(points, window=5, gap=10)
    assert len(segs) == 2


def test_find_long_segments_sort_by_diagonal_offset():
    # Ověří, že nové třídění (offset y-x, pak x) správně detekuje diagonálu
    # Diagonála na offsetu 100: body (0,100),(1,101),(2,102),...
    # Přimíchané body z jiné diagonály (offset 0): (0,0),(1,1),...
    diagonal_main = [(i, i + 100) for i in range(30)]
    diagonal_noise = [(i, i) for i in range(30)]
    all_points = diagonal_main + diagonal_noise
    # Seřazení podle (y-x, x) — jako v opravené verzi
    sorted_points = sorted(all_points, key=lambda p: (p[1] - p[0], p[0]))
    segs = _find_long_segments(sorted_points, window=5, gap=10)
    # Musí najít obě diagonály (každá >= window=5 bodů)
    assert len(segs) == 2


# --- _dot_mat_distance_core ---

def test_distance_identical_sequences():
    seq = "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG" * 10
    dist = _dot_mat_distance_core(seq, seq, window=22, gap=44)
    assert dist == pytest.approx(0.0, abs=1e-6)


def test_distance_completely_different():
    seq1 = "A" * 500
    seq2 = "C" * 500
    dist = _dot_mat_distance_core(seq1, seq2, window=22, gap=44)
    assert dist > 0.5


def test_distance_symmetry():
    seq1 = "ATCGATCGATCG" * 20
    seq2 = "TAGCTAGCTAGC" * 20
    d12 = _dot_mat_distance_core(seq1, seq2, window=22, gap=44)
    d21 = _dot_mat_distance_core(seq2, seq1, window=22, gap=44)
    assert d12 == pytest.approx(d21, abs=1e-10)


def test_distance_reverse_complement():
    # Reverzní komplement by měl dát vzdálenost 0 (nebo velmi malou)
    seq = "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG" * 10
    seq_rc = str(Seq(seq).reverse_complement())
    dist = _dot_mat_distance_core(seq, seq_rc, window=22, gap=44)
    assert dist == pytest.approx(0.0, abs=1e-6)


def test_distance_circular_permutation():
    # Cirkulárně permutovaná sekvence — vzdálenost by měla být blízká 0
    base = "ATCGATCGATCGATCGATCGATCGATCGATCG" * 30   # 960 bp
    offset = 200
    permuted = base[offset:] + base[:offset]
    dist = _dot_mat_distance_core(base, permuted, window=22, gap=44)
    assert dist < 0.05


def test_distance_range():
    seq1 = "ATCGATCG" * 50
    seq2 = "GCTAGCTA" * 50
    dist = _dot_mat_distance_core(seq1, seq2, window=22, gap=44)
    assert 0.0 <= dist <= 1.0


# --- dotplot_tree ---

def test_tree_identical_sequences():
    seq = "ATCGATCGATCGATCGATCGATCGATCGATCG" * 30
    seqs = [
        make_seq(seq, "A"),
        make_seq(seq, "B"),
        make_seq(seq, "C"),
    ]
    tree, dist_mat = dotplot_tree(seqs, window=22, threads=2)
    assert dist_mat["A", "B"] == pytest.approx(0.0, abs=1e-6)
    assert dist_mat["A", "C"] == pytest.approx(0.0, abs=1e-6)
    assert dist_mat["B", "C"] == pytest.approx(0.0, abs=1e-6)


def test_tree_has_correct_leaves():
    seq = "ATCGATCGATCGATCGATCGATCG" * 20
    seqs = [make_seq(seq, name) for name in ["X", "Y", "Z"]]
    tree, _ = dotplot_tree(seqs, window=22, threads=1)
    leaf_names = {c.name for c in tree.get_terminals()}
    assert leaf_names == {"X", "Y", "Z"}


def test_tree_circular_permutation_clusters_together():
    # Tři cirkulárně permutované varianty by se měly shlukovat blízko k sobě
    base = "ATCGATCGATCGATCGATCGATCGATCGATCG" * 30
    seqs = [
        make_seq(base, "orig"),
        make_seq(base[100:] + base[:100], "perm1"),
        make_seq(base[200:] + base[:200], "perm2"),
    ]
    tree, dist_mat = dotplot_tree(seqs, window=22, threads=2)
    assert dist_mat["orig", "perm1"] < 0.05
    assert dist_mat["orig", "perm2"] < 0.05
    assert dist_mat["perm1", "perm2"] < 0.05

def test_window_must_match_the_seed_mask():
    """Only the 22 bp window is supported; the mask fixes the window length."""
    from plastree.dotplotree import DOTPLOT_WINDOW, SEED_MASK, _dot_mat_distance_core

    assert DOTPLOT_WINDOW == len(SEED_MASK) == 22
    seq = "ATCGATCGATCGATCGATCGATCGATCGATCG" * 4
    with pytest.raises(ValueError, match="window must be 22"):
        _dot_mat_distance_core(seq, seq, window=19, gap=38)
