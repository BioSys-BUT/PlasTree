# PlasTree

Build phylogenetic trees for plasmid sequences using dotplot, k-mer, and
pangenome presence/absence (Roary/Panaroo) methods, then combine them into a
consensus tree and split plasmids into clusters.

## Features

- Three independent tree-building methods: **dotplot**, **k-mer**, and
  **pangenome** presence/absence (Roary or Panaroo format)
- Weighted **consensus tree** when two or more methods are combined
- Automatic **clustering** of plasmids based on branch length, with support
  for multiple thresholds in a single run
- Optional PNG rendering of dotplots and the final cluster-colored tree

## Requirements

- Python 3.10+
- All other dependencies (biopython, numpy, pandas, scipy, matplotlib) are
  installed automatically with the package.

## Installation

Extract the archive and move into the `tool/` directory (it contains
`pyproject.toml`):

```bash
tar -xzf plastree.tar.gz
cd tool
```

Pick one of the environment options below. Both create an isolated
environment named `plastree`.

### Option A: venv

```bash
python3 -m venv plastree
source plastree/bin/activate
pip install --upgrade pip
pip install .
```

### Option B: conda / miniforge

```bash
conda create -n plastree python=3.12
conda activate plastree
pip install .
```

### Verify

```bash
plastree --help
```

If the help text prints with all the parameters (`--fasta`, `--methods`,
...), the install worked. Remember to (re)activate the environment
(`source plastree/bin/activate` or `conda activate plastree`) in every new
terminal session before running `plastree`.

## Usage

Basic run (dotplot method only):

```bash
plastree --fasta sequences.fasta --outdir output
```

All three methods, plus pangenome analysis via Panaroo, with multiple
clustering thresholds computed in one go:

```bash
plastree \
  --fasta sequences.fasta \
  --methods dotplot kmer pang \
  --pang_format panaroo \
  --pang_csv gene_presence_absence.csv \
  --dotplot_window 22 \
  --dotplot_threads 8 \
  --max_cluster_size 30 \
  --max_branch_length 0.05 0.1 0.2 0.4 \
  --outdir output
```

See `plastree --help` for the full list of options.

### Output (in `--outdir`)

| File | Description |
|---|---|
| `dotplot_tree.newick`, `kmer_tree.newick`, `pang_tree.newick` | tree from each selected method |
| `*_distmat.tsv` | corresponding distance matrix |
| `consensus_tree.newick` | consensus tree (only when 2+ methods are used) |
| `plasmid_clusters.tsv` (or `plasmid_clusters_bl<value>.tsv` when multiple `--max_branch_length` values are given) | cluster assignment per plasmid (cluster 0 = outlier) |
| `cluster_tree.png` (only with `--result_tree_png`) | tree rendered with clusters color-coded |

## Testing (optional, for development)

```bash
pip install pytest
pytest tests/ -q
```
