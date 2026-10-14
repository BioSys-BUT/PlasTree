"""The result must not depend on the order of the sequences in the input file."""

import os
import random
import sys

from Bio import SeqIO

from plastree import plasmid_tree_tool

EXAMPLES = os.path.join(os.path.dirname(__file__), "..", "examples", "plasmids.fasta")


def run(fasta, outdir, monkeypatch):
    monkeypatch.setattr(sys, "argv", [
        "plastree", "--fasta", str(fasta), "--methods", "dotplot", "kmer",
        "--max_branch_length", "0.1", "0.02", "--outdir", str(outdir)])
    plasmid_tree_tool.main()
    tables = {}
    for name in ("plasmid_clusters_bl0_1.tsv", "plasmid_clusters_bl0_02.tsv"):
        with open(os.path.join(outdir, name)) as fh:
            tables[name] = fh.read()
    return tables


def test_clusters_and_ids_do_not_depend_on_input_order(tmp_path, monkeypatch):
    records = list(SeqIO.parse(EXAMPLES, "fasta"))
    reference = None
    for seed in (None, 1, 2, 3):
        recs = records[:]
        if seed is not None:
            random.Random(seed).shuffle(recs)
        fasta = tmp_path / f"in_{seed}.fasta"
        SeqIO.write(recs, fasta, "fasta")
        tables = run(fasta, tmp_path / f"out_{seed}", monkeypatch)
        if reference is None:
            reference = tables
        assert tables == reference
