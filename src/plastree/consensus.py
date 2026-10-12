"""Weighted greedy consensus with a shared mapping of taxon names.

The selection, tie-breaking, branch-length and stopping conventions follow
Bio.Phylo.Consensus.majority_consensus in Biopython 1.87 with cutoff=0.
Clades are counted by taxon identity, not by each tree's local leaf order.
This implementation uses public Bio.Phylo objects and does not patch Biopython.
"""

from copy import deepcopy

from Bio.Phylo.BaseTree import Clade, Tree


def weighted_greedy_consensus(trees, *, taxon_order=None):
    """Return a greedy consensus of identically labelled, rooted input trees.

    Integer weights are represented by repeated input trees. For example,
    ``[plotree, pangtree, plotree, kmertree, plotree, pangtree]`` implements
    3:2:1. A zero-weight method is omitted. Input trees are never modified.

    Clades are ranked by decreasing occurrence count, number of taxa and
    membership bit pattern. One shared taxon order is used throughout:
    the first tree's terminal order, or the explicit ``taxon_order``. Fix that
    order when testing invariance to rotations or permutations of input trees.
    Ties can depend on this order. Compatible clades are admitted greedily,
    without a minimum support threshold; this is not a strict majority rule.

    Internal confidence is 100 times weighted clade frequency. Internal branch
    lengths are averaged over occurrences containing that clade (missing
    lengths count as zero); terminal lengths/metadata come from the first
    tree. No branch-length normalization or rerooting is performed. As in the
    original builder, selection stops at a binary rooted tree or a fully
    resolved tree with three children at its root. The returned Tree uses
    Bio.Phylo's default rooted flag (False), matching the original builder,
    although the selection itself compares rooted clades.

    Raises ValueError for empty input, unnamed/duplicate/mismatched taxa,
    invalid taxon_order, or unary internal nodes (ambiguous repeated clades).
    """
    trees = list(trees)
    if not trees:
        raise ValueError("At least one input tree is required.")

    first_tips = trees[0].get_terminals()
    first_names = [tip.name for tip in first_tips]
    if any(not isinstance(n, str) or not n for n in first_names):
        raise ValueError("Every terminal must have a non-empty string name.")
    if len(set(first_names)) != len(first_names):
        raise ValueError("Duplicate taxon names in the first tree.")
    order = first_names if taxon_order is None else list(taxon_order)
    if len(order) != len(first_names) or set(order) != set(first_names):
        raise ValueError("taxon_order must contain every taxon exactly once.")
    universe = frozenset(order)
    counts = {}
    lengths = {}
    for index, tree in enumerate(trees):
        names = [tip.name for tip in tree.get_terminals()]
        if len(names) != len(universe) or set(names) != universe:
            raise ValueError(f"Tree {index + 1} has missing, duplicate or mismatched taxa.")
        for node in tree.get_nonterminals():
            if len(node.clades) < 2:
                raise ValueError(f"Tree {index + 1} contains a unary internal node.")
            members = frozenset(tip.name for tip in node.get_terminals())
            counts[members] = counts.get(members, 0) + 1
            lengths[members] = lengths.get(members, 0.0) + (node.branch_length or 0.0)

    tips_by_name = {tip.name: deepcopy(tip) for tip in first_tips}
    if len(universe) == 1:
        return Tree(root=tips_by_name[order[0]])
    tips = [tips_by_name[name] for name in order]
    root = Clade(clades=list(tips))
    accepted = {universe: root}
    ranked = sorted(
        (s for s in counts if s != universe),
        key=lambda s: (counts[s], len(s), tuple(name in s for name in order)),
        reverse=True,
    )
    for members in ranked:
        if any(not (members <= s or s <= members or members.isdisjoint(s))
               for s in accepted):
            continue
        parent_members = min((s for s in accepted if members < s), key=len)
        parent = accepted[parent_members]
        # Move only maximal accepted descendants; smaller ones are already
        # nested inside them. Size ordering preserves the original convention.
        children = []
        for s in sorted(accepted, key=len, reverse=True):
            if s < members and all(s.isdisjoint(c) for c in children):
                children.append(s)
        direct_tips = [tips_by_name[name] for name in order if name in members]
        node = Clade(
            branch_length=lengths[members] / counts[members],
            confidence=100.0 * counts[members] / len(trees),
            clades=direct_tips,
        )
        parent.clades = [c for c in parent.clades if c not in direct_tips]
        parent.clades.append(node)
        covered = frozenset().union(*children)
        for child_members in children:
            child = accepted[child_members]
            parent.clades.remove(child)
            node.clades.append(child)
        node.clades = [c for c in node.clades
                       if not (c.is_terminal() and c.name in covered)]
        # Match the builder's insertion-order convention for subsequent ties.
        del accepted[parent_members]
        accepted[parent_members] = parent
        accepted[members] = node
        if len(accepted) == len(tips) - 1 or (
            len(accepted) == len(tips) - 2 and len(root.clades) == 3
        ):
            break
    return Tree(root=root)
