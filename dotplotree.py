# -*- coding: utf-8 -*-
# dotplotree.py
# Tree construction using dotplot-based similarity between DNA sequences

from Bio.Seq import Seq
from Bio.Phylo.TreeConstruction import DistanceMatrix
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
import numpy as np
import math
import os
import matplotlib.pyplot as plt
from concurrent.futures import ProcessPoolExecutor, as_completed


def _replace_with_mask(sequence, mask, replacement_char='N'):
    return Seq(''.join(replacement_char if mask_bit == '0' else char for char, mask_bit in zip(sequence, mask)))


def _find_long_segments(points, window, gap):
    segments = []
    if not points:
        return segments

    current_segment = [points[0]]
    half_window = (window - 1) / 2

    for i in range(1, len(points)):
        distance_x = abs(points[i][0] - current_segment[-1][0])
        distance_y = abs(points[i][1] - current_segment[-1][1])

        if distance_x <= half_window and distance_y <= half_window and distance_x == distance_y:
            current_segment.append(points[i])
        else:
            if len(current_segment) >= window:
                segments.append(current_segment)
            current_segment = [points[i]]

    if len(current_segment) >= window:
        segments.append(current_segment)

    final_segments = []
    current_final_segment = []

    for segment in segments:
        if not current_final_segment:
            current_final_segment = segment
        else:
            distance_x = abs(current_final_segment[-1][0] - segment[0][0])
            distance_y = abs(current_final_segment[-1][1] - segment[0][1])

            if distance_x <= gap and distance_y <= gap and abs(distance_x - distance_y) <= round(window / 2):
                current_final_segment.extend(segment[1:])
            else:
                if len(current_final_segment) >= 2 * window:
                    final_segments.append(current_final_segment)
                current_final_segment = segment

    if current_final_segment and len(current_final_segment) >= 2 * window:
        final_segments.append(current_final_segment)

    return final_segments


def _dot_mat_distance_core(seq1_str, seq2_str, window, gap):
    if window == 21:
        mask = '1110101100110010101111'
    else:
        mask = '1110100110010101111'

    seq1 = seq1_str.upper()
    seq2 = seq2_str.upper()
    seq2rc = str(Seq(seq2).reverse_complement())
    dict_one, dict_two, dict_two2 = {}, {}, {}

    for seq, section_dict in [(seq1, dict_one), (seq2, dict_two), (seq2rc, dict_two2)]:
        for i in range(len(seq) - window + 1):
            section = _replace_with_mask(seq[i: i + window], mask)
            section_dict.setdefault(section, []).append(i)

    def get_match_coordinates(matches, one, two):
        x, y = [], []
        for section in matches:
            for i in one[section]:
                for j in two[section]:
                    x.append(i)
                    y.append(j)
        return x, y

    matches = set(dict_one).intersection(dict_two)
    matches2 = set(dict_one).intersection(dict_two2)
    x, y = get_match_coordinates(matches, dict_one, dict_two)
    x2, y2 = get_match_coordinates(matches2, dict_one, dict_two2)

    segments_direct = _find_long_segments(sorted(zip(x, y)), window, gap)
    segments_reverse = _find_long_segments(sorted(zip(x2, y2)), window, gap)

    bin_dir1 = np.zeros(len(seq1), dtype=int)
    bin_dir2 = np.zeros(len(seq2), dtype=int)
    bin_rev1 = np.zeros(len(seq1), dtype=int)
    bin_rev2 = np.zeros(len(seq2), dtype=int)

    for segment in segments_direct:
        bin_dir1[segment[0][0]:segment[-1][0] + window] = 1
        bin_dir2[segment[0][1]:segment[-1][1] + window] = 1

    for segment in segments_reverse:
        bin_rev1[segment[0][0]:segment[-1][0] + window] = 1
        bin_rev2[segment[0][1]:segment[-1][1] + window] = 1

    bin_rrev2 = bin_rev2[::-1]
    merged_bin1 = np.logical_or(bin_dir1, bin_rev1).astype(int)
    merged_bin2 = np.logical_or(bin_dir2, bin_rrev2).astype(int)

    bin_match = np.sum(merged_bin1) + np.sum(merged_bin2)
    return 1 - (bin_match / (len(seq1) + len(seq2)))


def dotplot_tree(sequences, window=21, dotplot_view=False, output_dir=None, workers=1):
    """
    Construct a UPGMA tree based on dotplot-derived similarity between sequences.

    Args:
        sequences (list of SeqRecord): Input sequences.
        window (int): Size of the sliding window for dotplot comparison.
        dotplot_view (bool): Whether to render dotplots for each pair.
        workers (int): Number of worker processes for pairwise distance calculation.

    Returns:
        tuple: (Phylogenetic tree, Distance matrix)
    """
    def get_unique_filename(base_name, ext, max_digits=3):
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
            
    def dotplotxy_segments(x, y, x2, y2, segments, segments2, window, mask):
        """
        Visualize the dotplot and matching segments.

        Args:
            x, y, x2, y2 (list): Match coordinates.
            segments, segments2 (list): Direct and reverse matching segments.
            window (int): Window size.
            mask (str): Binary mask used for filtering.
        """
        plt.scatter(x, y, color='black', marker='.', linewidths=0.05)
        plt.scatter(x2, len(rec_two)-np.array(y2), color='gray', marker='.', linewidths=0.05)

        for segment in segments:
            x_segment = [point[0] for point in segment]
            y_segment = [point[1] for point in segment]
            plt.plot(x_segment, y_segment, color='red')

        for segment in segments2:
            x_segment = [point[0] for point in segment]
            y_segment = [point[1] for point in segment]
            plt.plot(x_segment, len(rec_two)-np.array(y_segment), color='#00FF00')

        plt.xlabel("%s (length %i bp)" % (rec_one.id, len(rec_one)))
        plt.ylabel("%s (length %i bp)" % (rec_two.id, len(rec_two)))
        plt.grid(True)
        plt.xlim(0, len(rec_one) - window)
        plt.ylim(0, len(rec_two) - window)
        plt.gca().set_aspect('equal', adjustable='box')
        plt.tight_layout()
        plt.savefig(get_unique_filename(os.path.join(output_dir, "dotplot"), "png"), dpi=150)
        plt.close()

    def dot_mat_distance(seq1, seq2, window, gap, dotplot):
        """
        Calculate dotplot-based distance between two sequences.

        Args:
            seq1, seq2 (SeqRecord): DNA sequences.
            window (int): Size of sliding window.
            gap (int): Gap threshold for segment merging.
            dotplot (bool): Whether to visualize.

        Returns:
            float: Dotplot-based distance.
        """
        if window == 21:
            mask = '1110101100110010101111'
        else:
            mask = '1110100110010101111'

        seq2rc = seq2.seq.reverse_complement()
        dict_one, dict_two, dict_two2 = {}, {}, {}

        for seq, section_dict in [(seq1.seq.upper(), dict_one),
                                  (seq2.seq.upper(), dict_two),
                                  (seq2rc, dict_two2)]:
            for i in range(len(seq) - window + 1):
                section = _replace_with_mask(seq[i: i + window], mask)
                section_dict.setdefault(section, []).append(i)

        def get_match_coordinates(matches, dict_one, dict_two):
            x, y = [], []
            for section in matches:
                for i in dict_one[section]:
                    for j in dict_two[section]:
                        x.append(i)
                        y.append(j)
            return x, y

        matches = set(dict_one).intersection(dict_two)
        matches2 = set(dict_one).intersection(dict_two2)
        x, y = get_match_coordinates(matches, dict_one, dict_two)
        x2, y2 = get_match_coordinates(matches2, dict_one, dict_two2)

        segments_direct = _find_long_segments(sorted(zip(x, y)), window, gap)
        segments_reverse = _find_long_segments(sorted(zip(x2, y2)), window, gap)

        if dotplot:
            dotplotxy_segments(x, y, x2, y2, segments_direct, segments_reverse, window, mask)

        bin_dir1 = np.zeros(len(seq1), dtype=int)
        bin_dir2 = np.zeros(len(seq2), dtype=int)
        bin_rev1 = np.zeros(len(seq1), dtype=int)
        bin_rev2 = np.zeros(len(seq2), dtype=int)

        for segment in segments_direct:
            bin_dir1[segment[0][0]:segment[-1][0]+window] = 1
            bin_dir2[segment[0][1]:segment[-1][1]+window] = 1

        for segment in segments_reverse:
            bin_rev1[segment[0][0]:segment[-1][0]+window] = 1
            bin_rev2[segment[0][1]:segment[-1][1]+window] = 1

        bin_rrev2 = bin_rev2[::-1]
        merged_bin1 = np.logical_or(bin_dir1, bin_rev1).astype(int)
        merged_bin2 = np.logical_or(bin_dir2, bin_rrev2).astype(int)

        bin_match = np.sum(merged_bin1) + np.sum(merged_bin2)
        bin_dist = 1 - (bin_match / (len(seq1) + len(seq2)))
        return bin_dist

    gap = 2 * window
    names = [s.id for s in sequences]
    dist_mat = DistanceMatrix(names)
    n = len(sequences)
    total_pairs = n * (n - 1) // 2
    done_pairs = 0
    progress_step = max(1, math.ceil(total_pairs / 100)) if total_pairs else 1

    if workers < 1:
        raise ValueError("workers must be >= 1")

    if dotplot_view and len(sequences) > 20:
        dotplot_view = False
        print("Too many sequences. Dot plot rendering has been disabled. Maximum 20 sequences.")

    if workers > 1 and dotplot_view:
        dotplot_view = False
        print("Dot plot rendering has been disabled for parallel execution (workers > 1).")

    if total_pairs:
        print(f"Dotplot progress: 0/{total_pairs} pairs", flush=True)

    if workers == 1:
        for i, rec_one in enumerate(sequences):
            for j in range(i + 1, len(sequences)):
                rec_two = sequences[j]
                dist_mat[i, j] = dot_mat_distance(rec_one, rec_two, window, gap, dotplot_view)
                done_pairs += 1
                if done_pairs % progress_step == 0 or done_pairs == total_pairs:
                    print(f"Dotplot progress: {done_pairs}/{total_pairs} pairs", flush=True)
    else:
        future_to_pair = {}
        with ProcessPoolExecutor(max_workers=workers) as executor:
            for i, rec_one in enumerate(sequences):
                for j in range(i + 1, len(sequences)):
                    rec_two = sequences[j]
                    future = executor.submit(_dot_mat_distance_core, str(rec_one.seq), str(rec_two.seq), window, gap)
                    future_to_pair[future] = (i, j)

            for future in as_completed(future_to_pair):
                i, j = future_to_pair[future]
                dist_mat[i, j] = future.result()
                done_pairs += 1
                if done_pairs % progress_step == 0 or done_pairs == total_pairs:
                    print(f"Dotplot progress: {done_pairs}/{total_pairs} pairs", flush=True)

    constructor = DistanceTreeConstructor()
    tree = constructor.upgma(dist_mat)
    return tree, dist_mat