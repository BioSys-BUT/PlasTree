# plasmid_tree_tool.py
# Unified script to compute plasmid phylogenetic trees using dotplot, k-mer, and pangenome presence/absence methods
# Generates individual and consensus trees and identifies clusters

import os
import argparse
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend for headless environments
import matplotlib.pyplot as plt
import pandas as pd
from Bio import SeqIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceMatrix
from Bio.Phylo.Consensus import majority_consensus
from dotplotree import dotplot_tree
from pangtree import pang_tree
from kmertree import kmer_tree

def get_unique_filename(base_name, ext, max_digits=2):
    """
    Generate a unique filename by appending a numeric suffix if file exists.

    Args:
        base_name (str): Base filename without extension.
        ext (str): File extension (e.g. 'tsv', 'newick').
        max_digits (int): Number of digits for numbering.

    Returns:
        str: Unique filename with extension.
    """
    if not os.path.exists(f"{base_name}.{ext}"):
        return f"{base_name}.{ext}"
    i = 1
    while True:
        candidate = f"{base_name}{str(i).zfill(max_digits)}.{ext}"
        if not os.path.exists(candidate):
            return candidate
        i += 1

def get_sequences(fasta_path):
    """
    Load sequences from FASTA file.

    Args:
        fasta_path (str): Path to FASTA file.

    Returns:
        list: List of SeqRecord objects.
    """
    return list(SeqIO.parse(fasta_path, "fasta"))

def save_tree(tree, filename):
    """
    Save phylogenetic tree to Newick file.

    Args:
        tree (Tree): Bio.Phylo tree object.
        filename (str): Path to output file.
    """
    Phylo.write(tree, filename, "newick")

def save_distance_matrix(dist_mat, filename):
    """
    Save distance matrix to TSV file.

    Args:
        dist_mat (DistanceMatrix): Bio.Phylo DistanceMatrix object.
        filename (str): Path to output file.
    """
    with open(filename, "w") as f:
        f.write("\t" + "\t".join(dist_mat.names) + "\n")
        for i, row_name in enumerate(dist_mat.names):
            row = [str(dist_mat[i, j]) if j < i else "" for j in range(len(dist_mat.names))]
            f.write(row_name + "\t" + "\t".join(row) + "\n")

def largest_distant_nodes_with_table(tree, max_branches_in_cluster=None):
    """
    Identify and group most distinct clusters from a tree.

    Args:
        tree (Tree): Bio.Phylo tree object.
        max_branches_in_cluster (int): Max OTUs per cluster.

    Returns:
        tuple: (DataFrame with cluster IDs, list of cluster OTUs)
    """
    def get_sorted_branches(tree):
        branches = []
        for node in tree.get_nonterminals():
            branch = [leaf.name for leaf in node.get_terminals()]
            branch.sort()
            branches.append(branch)
        return set(tuple(branch) for branch in branches)

    branches_t1 = get_sorted_branches(tree)
    total_branches = len(branches_t1) + 1

    if max_branches_in_cluster is None:
        max_branches_in_cluster = total_branches // 2

    sorted_clusters = sorted(branches_t1, key=len, reverse=True)
    selected_clusters = []
    selected_branches = set()
    cluster_to_branch_map = {}
    cluster_id = 1

    for cluster in sorted_clusters:
        if len(cluster) <= max_branches_in_cluster and not selected_branches.intersection(cluster):
            selected_clusters.append(cluster)
            selected_branches.update(cluster)
            for branch in cluster:
                cluster_to_branch_map[branch] = cluster_id
            cluster_id += 1

    all_branches = {branch for cluster in branches_t1 for branch in cluster}
    unassigned_branches = all_branches - selected_branches
    for branch in unassigned_branches:
        cluster_to_branch_map[branch] = 0

    table = [{"OTU": branch, "Cluster ID": cluster_to_branch_map[branch]} for branch in all_branches]
    df = pd.DataFrame(table).sort_values(by="Cluster ID")
    return df, list(selected_clusters)

def main():
    parser = argparse.ArgumentParser(
        description="Plasmid Tree Tool: Build phylogenetic trees from plasmid sequences using dotplot, pangenome, and k-mer methods."
    )
    parser.add_argument("--fasta", required=True, help="Path to input FASTA file with plasmid sequences.")
    parser.add_argument("--pang_csv", default=None, help="Roary output CSV file for pangenome presence/absence matrix.")
    parser.add_argument("--methods", nargs="*", choices=["dotplot", "pang", "kmer"], default=["dotplot"], help="Tree construction methods to apply. Default: dotplot only.")
    parser.add_argument("--outdir", default="output", help="Directory to store output files. Default: /output.")
    parser.add_argument("--max_cluster_size", type=int, default=None, help="Maximum OTUs per cluster in consensus tree. Default: 1/2 number of sequences.")
    parser.add_argument("--dotplot_window", type=int, default=21, help="Size of the sliding window for dotplot comparison (e.g., 19 or 21).")
    parser.add_argument("--dotplot_png", action='store_true', help="Save dotplots for all pairwise comparisons as PNGs. Only for ≤20 sequences. Default: off.")
    parser.add_argument("--result_tree_png", action='store_true', help="Save consensus/result tree as PNG with node support and cluster coloring. Default: off.")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    sequences = get_sequences(args.fasta)
    trees = []

    def color_cluster_tree(tree, common_clusters):
        """
        Render and save tree with colored clusters.

        Args:
            tree (Tree): Bio.Phylo tree object.
            common_clusters (list): List of clusters (lists of OTU names).
        """
        # colors = ["#FF0000", "#0000FF", "#00FF00", "#800080", "#008080", "#FFA500", "#00FF85", "#FF00A5"]
        colors = ["#FF0000", "#0000FF", "#00FF00", "#800080", "#008080", "#FFA500", "#00FF85", "#FF00A5", "#A52A2A", "#FFFF00", "#00CED1", "#FF1493", "#7FFF00", "#DC143C", "#000080", "#FFD700", "#20B2AA", "#C71585", "#ADFF2F", "#FF6347", "#1E90FF", "#B22222", "#DA70D6", "#5F9EA0", "#40E0D0", "#FF4500", "#9ACD32", "#BA55D3", "#3CB371", "#8B4513"]
        for i, cluster in enumerate(common_clusters):
            if len(cluster) > 1:
                nodes = [{"name": name} for name in cluster]
                mrca = tree.common_ancestor(*nodes)
                mrca.color = colors[i % len(colors)]
        import warnings
        warnings.filterwarnings("ignore", category=UserWarning, module="Bio.Phylo._utils")
        fig = plt.figure(figsize=(10, 15))
        axes = fig.add_subplot(1, 1, 1)
        Phylo.draw(tree, axes=axes, branch_labels=lambda c: "{:.2f}".format(c.confidence) if c.confidence else '')
        # Save to PNG instead of showing
        plt.tight_layout()
        plt.savefig(get_unique_filename(os.path.join(args.outdir, "cluster_tree"), "png"), dpi=300)
        plt.close()

    # DOTPLOT METHOD
    if "dotplot" in args.methods:
        print("Running dotplot method...")
        dotmat_tree, dotmat_dist = dotplot_tree(sequences, args.dotplot_window, args.dotplot_png, args.outdir)
        trees.append(dotmat_tree)
        save_tree(dotmat_tree, get_unique_filename(os.path.join(args.outdir, "dotplot_tree"), "newick"))
        save_distance_matrix(dotmat_dist, get_unique_filename(os.path.join(args.outdir, "dotplot_distmat"), "tsv"))

    # PANGTREE METHOD
    if "pang" in args.methods:
        if args.pang_csv:
            print("Running pangenome binary method...")
            pang_tree_result, pang_dist = pang_tree(args.pang_csv)
            trees.append(pang_tree_result)
            if 'dotmat_tree' in locals() and dotmat_tree is not None:
                trees.append(dotmat_tree)
            save_tree(pang_tree_result, get_unique_filename(os.path.join(args.outdir, "pang_tree"), "newick"))
            save_distance_matrix(pang_dist, get_unique_filename(os.path.join(args.outdir, "pang_distmat"), "tsv"))
        else:
            print("Warning: Method 'pang' was selected but no --pang_csv provided. Skipping pangtree method.")

    # K-MER METHOD
    if "kmer" in args.methods:
        print("Running k-mer method...")
        kmer_result_tree, kmer_dist = kmer_tree(sequences)
        trees.append(kmer_result_tree)
        if 'dotmat_tree' in locals() and dotmat_tree is not None:
            trees.append(dotmat_tree)
        if 'pang_tree_result' in locals() and pang_tree_result is not None:
            trees.append(pang_tree_result)
        save_tree(kmer_result_tree, get_unique_filename(os.path.join(args.outdir, "kmer_tree"), "newick"))
        save_distance_matrix(kmer_dist, get_unique_filename(os.path.join(args.outdir, "kmer_distmat"), "tsv"))

    # CONSENSUS TREE
    if len(trees) > 1:
        print("Computing consensus tree...")
        result_tree = majority_consensus(trees)
        save_tree(result_tree, get_unique_filename(os.path.join(args.outdir, "consensus_tree"), "newick"))
        print("Consensus tree:")
    else:
        print("Only one method used — consensus tree not computed.")
        print(", ".join(args.methods) + " tree:")
        result_tree = trees[0]

    Phylo.draw_ascii(result_tree)

    # CLUSTERING + SAVE
    cluster_table, clusters = largest_distant_nodes_with_table(result_tree, args.max_cluster_size)
    cluster_table.to_csv(get_unique_filename(os.path.join(args.outdir, "plasmid_clusters"), "tsv"), sep="\t", index=False)

    if args.result_tree_png:
        color_cluster_tree(result_tree, clusters)

if __name__ == "__main__":
    main()
