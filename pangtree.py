# pangtree.py
# Tree construction based on binary pangenome presence/absence matrix
# Input is a CSV file in Roary or Panaroo output format

import pandas as pd
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, DistanceMatrix
from scipy.spatial.distance import hamming
import math

def pang_tree(input_file, expected_samples=None):
    """
    Construct a UPGMA tree from pangenome presence/absence matrix.

    This function reads a pangenome CSV file (Roary or Panaroo format),
    extracts the binary matrix of gene presence across samples, calculates pairwise
    Hamming distances between samples, and builds a UPGMA phylogenetic tree.

    Args:
        input_file (str): Path to pangenome CSV (Roary 'gene_presence_absence_roary.csv' 
                          or Panaroo 'gene_presence_absence.csv').
        expected_samples (list[str] | None): Optional ordered list of sample names
            expected in output (typically FASTA IDs). Missing samples are added
            as all-zero columns.

    Returns:
        tuple: (Bio.Phylo.BaseTree.Tree, DistanceMatrix)
            - The UPGMA phylogenetic tree.
            - The corresponding distance matrix.
    """
    # Load the CSV file
    csv_file = pd.read_csv(input_file, low_memory=False)

    # Identify format (Roary vs Panaroo) and find presence/absence start
    if 'Avg group size nuc' in csv_file.columns:
        # Roary format: find 'Avg group size nuc' and start after it
        presence_start_idx = csv_file.columns.get_loc('Avg group size nuc') + 1
    elif 'Annotation' in csv_file.columns:
        # Panaroo format: find 'Annotation' and start after it
        presence_start_idx = csv_file.columns.get_loc('Annotation') + 1
    else:
        raise ValueError("Could not determine pangenome format. Neither 'Avg group size nuc' (Roary) nor 'Annotation' (Panaroo) column found.")

    # Extract presence/absence matrix: columns after the metadata section
    df = csv_file.iloc[:, presence_start_idx:]

    # Convert all non-empty cells (presence) to 1, and missing values (absence) to 0
    df = df.notna().astype(int)

    if expected_samples is not None:
        # Map Roary columns like "sample_bakta" back to FASTA IDs when possible.
        renamed_columns = {}
        for col in df.columns:
            if col.endswith("_bakta"):
                base = col[:-6]
                if base in expected_samples:
                    renamed_columns[col] = base
        if renamed_columns:
            df = df.rename(columns=renamed_columns)

        # Add missing expected samples as all-zero presence/absence columns.
        missing_samples = [name for name in expected_samples if name not in df.columns]
        for name in missing_samples:
            df[name] = 0

        # Keep only expected samples and in FASTA order for cross-method consistency.
        df = df[[name for name in expected_samples]]

    # Extract sample names
    names = df.columns.tolist()
    total_pairs = len(names) * (len(names) - 1) // 2
    done_pairs = 0
    progress_step = max(1, math.ceil(total_pairs / 100)) if total_pairs else 1

    if total_pairs:
        print(f"Pang progress: 0/{total_pairs} pairs", flush=True)

    # Compute lower-triangular Hamming distance matrix
    matrix = []
    for i in range(len(names)):
        row = []
        for j in range(i):
            vec1 = df[names[i]].values
            vec2 = df[names[j]].values
            dist = hamming(vec1, vec2)
            row.append(dist)
            done_pairs += 1
            if done_pairs % progress_step == 0:
                print(f"Pang progress: {done_pairs}/{total_pairs} pairs", flush=True)
        matrix.append(row)

    # Add diagonal (0s)
    lower_triangle = [
        [matrix[i][j] for j in range(i)] + [0]  # Diagonal value = 0
        for i in range(len(names))
    ]

    # Build and return the tree
    dist_mat = DistanceMatrix(names, lower_triangle)
    constructor = DistanceTreeConstructor()
    tree = constructor.upgma(dist_mat)
    return tree, dist_mat
