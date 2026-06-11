# kmertree.py
# Tree construction based on ranked k-mer frequencies

from itertools import product, combinations
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor
import numpy as np
import math

def kmer_tree(sequences):
    """
    Constructs a phylogenetic tree based on ranked k-mer frequency profiles
    averaged across multiple k-mer lengths (3 to 8).

    Args:
        sequences (list): List of Bio.SeqRecord objects representing the plasmid sequences.

    Returns:
        tuple: A tuple (tree, distance_matrix), where tree is a UPGMA phylogenetic tree,
               and distance_matrix is a Bio.Phylo DistanceMatrix object.
    """

    def generate_kmers(k):
        """
        Generate all possible k-mers of length k composed of A, C, G, T.

        Args:
            k (int): Length of the k-mers.

        Returns:
            list: All possible k-mers of given length.
        """
        nucleotides = ['A', 'C', 'G', 'T']
        kmers = [''.join(p) for p in product(nucleotides, repeat=k)]
        return kmers

    def count_kmers(k, seq, kmers_set):
        """
        Count k-mer occurrences in a sequence and rank them by frequency.

        Args:
            k (int): Length of the k-mers.
            seq (str): Nucleotide sequence.
            kmers_set (set): Set of all valid k-mers.

        Returns:
            list: Normalized rank values for each k-mer in original order.
        """
        kmer_counts = {kmer: 0 for kmer in kmers_set}
        seq_len = len(seq)
        for i in range(seq_len - k + 1):
            kmer = seq[i: i + k]
            if kmer in kmer_counts:
                kmer_counts[kmer] += 1

        sorted_kmer_counts = sorted(kmer_counts.items(), key=lambda item: item[1], reverse=True)
        ranks = {}
        current_rank = 1
        previous_count = None

        for kmer, count in sorted_kmer_counts:
            if count != previous_count:
                ranks[kmer] = current_rank
                current_rank += 1
            else:
                ranks[kmer] = current_rank - 1
            previous_count = count

        num_kmers = len(kmers_set)
        normalized_ranks = {kmer: rank / num_kmers for kmer, rank in ranks.items()}
        output_ranks = [normalized_ranks[kmer] for kmer in kmers_set]
        return output_ranks

    def compute_all_kmer_ranks(sequences, k, rev):
        """
        Compute normalized k-mer frequency ranks for all sequences.

        Args:
            sequences (list): List of Bio.SeqRecord objects.
            k (int): Length of k-mers.
            rev (int): Unused parameter (placeholder).

        Returns:
            dict: Mapping from sequence ID to rank vector.
        """
        kmers = generate_kmers(k)
        kmers_set = set(kmers)
        kmer_ranks = {}
        total = len(sequences)
        step = max(1, math.ceil(total / 20)) if total else 1
        for idx, seq in enumerate(sequences, start=1):
            kmer_ranks[seq.id] = count_kmers(k, str(seq.seq) + str(seq.seq.reverse_complement()), kmers_set)
            if idx % step == 0 or idx == total:
                print(f"K-mer ranks (k={k}): {idx}/{total} sequences", flush=True)
        return kmer_ranks

    def pairwise_distance(ranks1, ranks2):
        """
        Compute Euclidean distance between two k-mer rank vectors.

        Args:
            ranks1 (list): Rank vector of first sequence.
            ranks2 (list): Rank vector of second sequence.

        Returns:
            float: Euclidean distance.
        """
        return sum((r1 - r2) ** 2 for r1, r2 in zip(ranks1, ranks2)) ** 0.5

    def distance_matrix(sequences, k):
        """
        Construct a distance matrix using k-mer rank vectors.

        Args:
            sequences (list): List of Bio.SeqRecord objects.
            k (int): Length of k-mers.

        Returns:
            DistanceMatrix: Lower triangular distance matrix.
        """
        names = [s.id for s in sequences]
        matrix = DistanceMatrix(names)
        kmer_ranks = compute_all_kmer_ranks(sequences, k, 0)
        total_pairs = len(names) * (len(names) - 1) // 2
        done_pairs = 0
        step = max(1, math.ceil(total_pairs / 100)) if total_pairs else 1
        if total_pairs:
            print(f"K-mer distance (k={k}): 0/{total_pairs} pairs", flush=True)
        for (tax1, tax2) in combinations(names, 2):
            matrix[tax1, tax2] = pairwise_distance(kmer_ranks[tax1], kmer_ranks[tax2])
            done_pairs += 1
            if done_pairs % step == 0 or done_pairs == total_pairs:
                print(f"K-mer distance (k={k}): {done_pairs}/{total_pairs} pairs", flush=True)
        return matrix

    def distance_matrix_to_full(dm):
        """
        Convert a Bio.Phylo DistanceMatrix to a symmetric NumPy matrix.

        Args:
            dm (DistanceMatrix): Input distance matrix.

        Returns:
            np.ndarray: Symmetric matrix of distances.
        """
        full = np.zeros((n, n))
        for i in range(n):
            for j in range(i):
                dist = dm[i, j]
                full[i][j] = dist
                full[j][i] = dist
        return full

    normalized_matrices = []
    names = [s.id for s in sequences]
    n = len(names)

    # Compute distance matrices for k from 3 to 8 and normalize them
    for k in range(3, 9):
        print(f"K-mer progress: starting k={k}", flush=True)
        dm = distance_matrix(sequences, k)
        full_matrix = distance_matrix_to_full(dm)
        max_val = np.max(full_matrix)
        normalized = full_matrix / max_val if max_val != 0 else full_matrix
        normalized_matrices.append(normalized)
        print(f"K-mer progress: finished k={k}", flush=True)

    # Average the normalized matrices
    average_matrix = np.mean(normalized_matrices, axis=0)
    lower_triangle = [
        [average_matrix[i][j] for j in range(i)] + [0]
        for i in range(n)
    ]

    dist_mat = DistanceMatrix(names, lower_triangle)
    constructor = DistanceTreeConstructor()
    tree = constructor.upgma(dist_mat)
    return tree, dist_mat
