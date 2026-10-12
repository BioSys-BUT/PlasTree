import io

import pytest
from Bio import Phylo
from Bio.Phylo.Consensus import majority_consensus

from plastree.consensus import weighted_greedy_consensus


def parse(newick):
    return Phylo.read(io.StringIO(newick), "newick")


def clade_sets(tree):
    return {frozenset(t.name for t in n.get_terminals())
            for n in tree.get_nonterminals()}


def test_majority_clade_is_kept_and_minority_dropped():
    trees = [parse("((A:1,B:1):1,(C:1,D:1):1);"),
             parse("((A:1,B:1):1,(C:1,D:1):1);"),
             parse("((A:1,C:1):1,(B:1,D:1):1);")]
    sets = clade_sets(weighted_greedy_consensus(trees))
    assert frozenset({"A", "B"}) in sets
    assert frozenset({"A", "C"}) not in sets


@pytest.mark.parametrize("newick", [
    "((A:1,B:1):1,(C:1,D:1):1);",
    "(((A:1,B:1):1,(C:1,D:1):1):1,E:1);",
    "(((A:1,B:1):1,(C:1,D:1):1):1,(E:1,F:1):1);",
])
def test_matches_biopython_when_leaf_orders_agree(newick):
    """Same result as the reference when the counting difference cannot bite.

    With identical input trees the terminal orderings coincide, so Biopython's
    bitstring counting and counting by taxon identity agree. This pins the
    selection, tie-breaking and stopping conventions to the reference,
    including its early stop, which can leave out a unanimously supported
    clade on very small trees.
    """
    mine = clade_sets(weighted_greedy_consensus([parse(newick), parse(newick)]))
    theirs = clade_sets(majority_consensus([parse(newick), parse(newick)]))
    assert mine == theirs


def test_counts_clades_by_taxon_not_by_leaf_order():
    """A clade shared by two differently shaped trees is counted as shared."""
    # {A,B} is present in both, but the trees list their terminals differently.
    t1 = parse("(((A:1,B:1):1,C:1):1,D:1);")
    t2 = parse("(D:1,(C:1,(A:1,B:1):1):1);")
    cons = weighted_greedy_consensus([t1, t2])
    ab = [n for n in cons.get_nonterminals()
          if {t.name for t in n.get_terminals()} == {"A", "B"}]
    assert len(ab) == 1
    assert ab[0].confidence == pytest.approx(100.0)


def test_weighting_by_repeated_trees():
    """Three copies of one tree outvote two copies of another."""
    major = parse("((A:1,B:1):1,(C:1,D:1):1);")
    minor = parse("((A:1,C:1):1,(B:1,D:1):1);")
    trees = [major, major, major, minor, minor]
    sets = clade_sets(weighted_greedy_consensus(trees))
    assert frozenset({"A", "B"}) in sets


def test_all_taxa_present_exactly_once():
    trees = [parse("(((A:1,B:1):1,C:1):1,D:1);"), parse("((A:1,B:1):1,(C:1,D:1):1);")]
    names = sorted(t.name for t in weighted_greedy_consensus(trees).get_terminals())
    assert names == ["A", "B", "C", "D"]


def test_deterministic_under_permutation_of_input_trees():
    t1 = parse("(((A:1,B:1):1,C:1):1,D:1);")
    t2 = parse("((A:1,B:1):1,(C:1,D:1):1);")
    order = ["A", "B", "C", "D"]
    a = clade_sets(weighted_greedy_consensus([t1, t2, t1], taxon_order=order))
    b = clade_sets(weighted_greedy_consensus([t1, t1, t2], taxon_order=order))
    c = clade_sets(weighted_greedy_consensus([t2, t1, t1], taxon_order=order))
    assert a == b == c


def test_internal_branch_length_is_mean_over_occurrences():
    t1 = parse("((A:1,B:1):0.4,(C:1,D:1):1);")
    t2 = parse("((A:1,B:1):0.6,(C:1,D:1):1);")
    cons = weighted_greedy_consensus([t1, t2])
    ab = next(n for n in cons.get_nonterminals()
              if {t.name for t in n.get_terminals()} == {"A", "B"})
    assert ab.branch_length == pytest.approx(0.5)


@pytest.mark.parametrize("trees, message", [
    ([], "At least one input tree"),
    ([parse("((A:1,B:1):1,(C:1,D:1):1);"), parse("((A:1,B:1):1,(C:1,E:1):1);")],
     "missing, duplicate or mismatched taxa"),
])
def test_invalid_input_is_rejected(trees, message):
    with pytest.raises(ValueError, match=message):
        weighted_greedy_consensus(trees)


def test_taxon_order_must_cover_every_taxon():
    trees = [parse("((A:1,B:1):1,(C:1,D:1):1);")]
    with pytest.raises(ValueError, match="taxon_order"):
        weighted_greedy_consensus(trees, taxon_order=["A", "B", "C"])
