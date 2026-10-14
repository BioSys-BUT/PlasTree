# PlasTree

A tree-based pipeline for clustering complete bacterial plasmids from long-read assemblies, combining structural, k-mer, and gene-content similarity into a single weighted consensus tree.

![PlasTree overview](https://raw.githubusercontent.com/BioSys-BUT/PlasTree/main/plastree_overview.png "PlasTree pipeline and an example consensus tree with color-coded clusters")

## Overview

PlasTree integrates three complementary similarity methods, DotPloTree (pairwise dotplot-based structural similarity), KmerTree (ranked k-mer frequency similarity), and PanGTree (pangenome presence/absence similarity), and combines them into a single weighted consensus tree. Plasmids are then automatically grouped into clusters based on branch length and cluster size thresholds. PlasTree is intended for complete plasmid assemblies from long-read sequencing (e.g. Oxford Nanopore Technologies), where rearrangements, recombination, and accessory gene turnover can obscure relationships found by a single similarity measure.

## Content

- [Key Features](#key-features)
- [Requirements](#requirements)
- [Installation](#installation)
- [Usage](#usage)
  - [Parameters](#parameters)
  - [Example](#example)
  - [Output](#output)
- [Changes in 0.2.1](#changes-in-021)
- [Changes in 0.2.0](#changes-in-020)
- [Generating pangenome input (Bakta + Roary/Panaroo)](#generating-pangenome-input-bakta--roarypanaroo)
- [License](#license)
- [Citation](#citation)
- [Contact](#contact)

## Key Features

- **Three independent tree-building methods**: DotPloTree (structural), KmerTree (nucleotide composition), and PanGTree (gene content, Roary or Panaroo format)
- **Weighted consensus tree** (2:1:1 DotPloTree:PanGTree:KmerTree) when two or more methods are combined; clades are counted by plasmid identity, so the result does not depend on the order in which the methods are run
- **Automatic clustering** of plasmids by branch length, with support for multiple thresholds computed in a single run
- **Installable as a standard Python package**, exposing a single `plastree` command
- Optional PNG rendering of pairwise dotplots and the final cluster-colored consensus tree

## Requirements

- Python 3.10+
- All other dependencies (Biopython, NumPy, Pandas, SciPy, Matplotlib) are installed automatically with the package.

## Installation

First create an isolated environment named `plastree` (pick one option below), then install PlasTree into it.

**Option A: venv**
```bash
python3 -m venv plastree
source plastree/bin/activate
pip install --upgrade pip
```

**Option B: conda / miniforge**
```bash
conda create -n plastree python=3.12
conda activate plastree
```

With the environment active, install PlasTree from PyPI:
```bash
pip install plastree
```

Or, install directly from GitHub (e.g. to get unreleased changes):
```bash
pip install git+https://github.com/BioSys-BUT/PlasTree.git
```

Or, if you want a local copy of the source:
```bash
git clone https://github.com/BioSys-BUT/PlasTree.git
cd PlasTree
pip install .
```

Verify the install:
```bash
plastree --help
```

## Usage

```bash
plastree --fasta <fasta_path> [--pang_format {roary,panaroo}] [--pang_csv <csv_path>] [--methods dotplot kmer pang] [--outdir <out_dir>] [--max_cluster_size <int>] [--max_branch_length <float> [<float> ...]] [--dotplot_threads <int>] [--dotplot_png] [--result_tree_png]
```

### Parameters

- `--fasta`: Path to a single FASTA file containing all plasmid sequences (one record per plasmid). **Required.**
- `--pang_format {roary,panaroo}`: Format of the pangenome CSV. Required when using `--methods pang`.
- `--pang_csv`: Path to the pangenome presence/absence CSV file. Required when using `--methods pang`.
- `--methods [dotplot pang kmer ...]`: Tree construction methods to apply. Default: `dotplot` only.
- `--outdir`: Directory to store output files. Default: `output`.
- `--max_cluster_size`: Maximum number of plasmid sequences per cluster. Default: half the number of input sequences.
- `--max_branch_length`: Maximum allowed branch length within a cluster (0-1). Accepts one or more values (e.g. `0.1 0.2 0.05`); one cluster-assignment table is produced per value. Default: `0.02 0.1`.
- `--dotplot_threads`: Number of parallel worker processes for dotplot pairwise distance computation. Default: `1`.
- `--dotplot_png`: Save pairwise dotplots as PNGs. Only for <=20 sequences. Default: off.
- `--result_tree_png`: Save the consensus tree as a PNG with node support and cluster coloring. Default: off.

### Example

```bash
plastree \
  --fasta sequences.fasta \
  --methods dotplot kmer pang \
  --pang_format panaroo \
  --pang_csv gene_presence_absence.csv \
  --dotplot_threads 8 \
  --max_cluster_size 30 \
  --max_branch_length 0.05 0.1 0.2 0.4 \
  --result_tree_png \
  --outdir output
```

with example data (8 small real plasmids from *Enterococcus faecium*, included in [`examples/`](examples)):
```bash
plastree \
  --fasta examples/plasmids.fasta \
  --methods dotplot kmer pang \
  --pang_format panaroo \
  --pang_csv examples/gene_presence_absence.csv \
  --result_tree_png \
  --outdir examples/output
```

### Output

All output files are written to `--outdir`.

| File | Description |
|---|---|
| `dotplot_tree.newick` | DotPloTree structural similarity tree (`--methods dotplot`) |
| `dotplot_distmat.tsv` | DotPloTree distance matrix |
| `kmer_tree.newick` | KmerTree nucleotide composition tree (`--methods kmer`) |
| `kmer_distmat.tsv` | KmerTree distance matrix |
| `pang_tree.newick` | PanGTree gene content tree (`--methods pang`) |
| `pang_distmat.tsv` | PanGTree distance matrix |
| `consensus_tree.newick` | weighted consensus tree (only when 2+ methods are used) |
| `plasmid_clusters.tsv` | cluster assignment per plasmid (cluster 0 = outlier) |
| `cluster_tree.png` | consensus tree rendered with clusters color-coded (only with `--result_tree_png`) |

If multiple `--max_branch_length` values are given, `plasmid_clusters.tsv` and `cluster_tree.png` are produced once per threshold, suffixed with the value (dots replaced by underscores), e.g. `plasmid_clusters_bl0_2.tsv` and `cluster_tree_0_2.png` for threshold `0.2`.

## Changes in 0.2.1

The sequences are sorted by identifier after the FASTA file is read, and clades of equal size are ordered by their members, so the clusters and their numbering no longer depend on the order of the sequences in the input file (UPGMA resolves tied distances by the order of the matrix). The results for the example data and for the datasets of the manuscript are unchanged.

## Changes in 0.2.0

Version 0.2.0 changes the results, so trees and clusters will differ from those produced by 0.1.x for the same input.

- **PanGTree** uses Jaccard instead of Hamming distance. Gene groups absent from both plasmids no longer count as agreement, which matters because the presence/absence matrix is sparse.
- **KmerTree** divides each distance matrix by the largest value attainable for a given k (`2**k`) instead of the largest value observed in the input. The same pair of plasmids now receives the same distance regardless of which other plasmids are in the file, so a `--max_branch_length` threshold has a fixed meaning across analyses.
- **Consensus tree** is computed by a weighted greedy consensus that counts clades by plasmid identity (DotPloTree : PanGTree : KmerTree = 2:1:1). Previously Biopython's `majority_consensus` was used, which compared clades through each tree's own leaf order.
- **`--dotplot_window` was removed.** The window length is now fixed at 22 bp by the spaced seed. Remove this option from existing commands, otherwise PlasTree exits with an error.
- **Default `--max_branch_length`** changed from `0.2` to `0.02 0.1`. Pass the values explicitly (e.g. `--max_branch_length 0.05 0.1 0.2 0.4`) to reproduce earlier thresholds.

## Generating pangenome input (Bakta + Roary/Panaroo)

The `--pang_csv` file used by `--methods pang` is not produced by PlasTree itself. It comes from annotating the same plasmid FASTA files and clustering the predicted genes into orthologous groups, using external tools:

1. Annotate each plasmid with [Bakta](https://github.com/oschwengers/bakta) (or another tool that outputs GFF3 annotations), producing one GFF3 file per plasmid.
2. Collect all resulting GFF3 files into one directory.
3. Run [Roary](https://sanger-pathogens.github.io/Roary/) and/or [Panaroo](https://github.com/gtonkinhill/panaroo) on the collected GFF3 files:
```bash
roary *.gff3 -f roary_out
```
```bash
panaroo -i *.gff3 --clean-mode sensitive -o panaroo_out
```

Both produce a `gene_presence_absence.csv`. Pass this file as `--pang_csv`, together with the matching `--pang_format roary` or `--pang_format panaroo`. See each tool's own documentation for annotation database setup, threading, and other options. Bakta, Roary, and Panaroo have their own, sometimes conflicting, dependencies, so it is common to install them into separate environments (e.g. separate conda environments).

## License

See [LICENSE](LICENSE) for details (MIT).

## Citation

The manuscript describing PlasTree is currently in submission. Until it is published, if you use PlasTree in your research, please cite this repository:

```bibtex
@software{plastree2026,
  title={PlasTree},
  author={Vitkova, Helena and Jakubickova, Marketa and Bezdicek, Matej and Holubova, Ema and Lengerova, Martina},
  year={2026},
  url={https://github.com/BioSys-BUT/PlasTree}
}
```

## Contact

For questions and feedback, please open an issue on GitHub or contact us via email at [jakubickova@vut.cz](mailto:jakubickova@vut.cz).