import os
import pytest
import tempfile
import pandas as pd
from plastree.pangtree import pang_tree


def write_panaroo_csv(tmp_path, data: dict) -> str:
    """Zapíše CSV ve formátu Panaroo a vrátí cestu."""
    # Sloupce: Gene, Non-unique Gene name, Annotation, sample1, sample2, ...
    samples = list(data.keys())
    genes = list(data[samples[0]].keys())
    rows = []
    for gene in genes:
        row = {"Gene": gene, "Non-unique Gene name": "", "Annotation": ""}
        for sample in samples:
            row[sample] = gene if data[sample][gene] else ""
        rows.append(row)
    df = pd.DataFrame(rows)
    path = os.path.join(str(tmp_path), "gene_presence_absence.csv")
    df.to_csv(path, index=False)
    return path


def write_roary_csv(tmp_path, data: dict) -> str:
    """Zapíše CSV ve formátu Roary a vrátí cestu."""
    samples = list(data.keys())
    genes = list(data[samples[0]].keys())
    rows = []
    for gene in genes:
        row = {
            "Gene": gene,
            "Non-unique Gene name": "",
            "Annotation": "",
            "No. isolates": "",
            "No. sequences": "",
            "Avg sequences per isolate": "",
            "Genome Fragment": "",
            "Order within Fragment": "",
            "Accessory Fragment": "",
            "Accessory Order with Fragment": "",
            "QC": "",
            "Min group size nuc": "",
            "Max group size nuc": "",
            "Avg group size nuc": "",
        }
        for sample in samples:
            row[sample] = gene if data[sample][gene] else ""
        rows.append(row)
    df = pd.DataFrame(rows)
    path = os.path.join(str(tmp_path), "gene_presence_absence_roary.csv")
    df.to_csv(path, index=False)
    return path


@pytest.fixture
def tmp_path():
    with tempfile.TemporaryDirectory() as d:
        yield d


def test_panaroo_identical_samples(tmp_path):
    # Dvě identická vzorky → vzdálenost 0
    data = {
        "A": {"gene1": True, "gene2": True, "gene3": False},
        "B": {"gene1": True, "gene2": True, "gene3": False},
    }
    csv = write_panaroo_csv(tmp_path, data)
    tree, dist_mat = pang_tree(csv, expected_format="panaroo")
    assert dist_mat["A", "B"] == pytest.approx(0.0, abs=1e-6)


def test_panaroo_completely_different_samples(tmp_path):
    # Vzorky bez sdílených genů → vzdálenost 1
    data = {
        "A": {"gene1": True,  "gene2": True,  "gene3": False, "gene4": False},
        "B": {"gene1": False, "gene2": False, "gene3": True,  "gene4": True},
    }
    csv = write_panaroo_csv(tmp_path, data)
    tree, dist_mat = pang_tree(csv, expected_format="panaroo")
    assert dist_mat["A", "B"] == pytest.approx(1.0, abs=1e-6)


def test_roary_format(tmp_path):
    data = {
        "A": {"gene1": True, "gene2": False},
        "B": {"gene1": True, "gene2": True},
        "C": {"gene1": False, "gene2": True},
    }
    csv = write_roary_csv(tmp_path, data)
    tree, dist_mat = pang_tree(csv, expected_format="roary")
    leaf_names = {c.name for c in tree.get_terminals()}
    assert leaf_names == {"A", "B", "C"}


def test_format_mismatch_raises(tmp_path):
    data = {
        "A": {"gene1": True},
        "B": {"gene1": False},
    }
    csv = write_panaroo_csv(tmp_path, data)
    with pytest.raises(ValueError, match="format mismatch"):
        pang_tree(csv, expected_format="roary")


def test_missing_sample_added_as_zeros(tmp_path):
    # expected_samples obsahuje vzorek, který v CSV není → přidá se jako all-zero
    data = {
        "A": {"gene1": True, "gene2": True},
        "B": {"gene1": True, "gene2": False},
    }
    csv = write_panaroo_csv(tmp_path, data)
    tree, dist_mat = pang_tree(csv, expected_samples=["A", "B", "C"], expected_format="panaroo")
    leaf_names = {c.name for c in tree.get_terminals()}
    assert "C" in leaf_names
    # C je all-zero, A má všechny geny → vzdálenost = 1
    assert dist_mat["A", "C"] == pytest.approx(1.0, abs=1e-6)


def test_expected_samples_order(tmp_path):
    data = {
        "A": {"gene1": True},
        "B": {"gene1": False},
        "C": {"gene1": True},
    }
    csv = write_panaroo_csv(tmp_path, data)
    tree, dist_mat = pang_tree(csv, expected_samples=["C", "A", "B"], expected_format="panaroo")
    leaf_names = {c.name for c in tree.get_terminals()}
    assert leaf_names == {"A", "B", "C"}


def test_distance_symmetry(tmp_path):
    data = {
        "A": {"gene1": True, "gene2": False, "gene3": True},
        "B": {"gene1": False, "gene2": True, "gene3": True},
    }
    csv = write_panaroo_csv(tmp_path, data)
    tree, dist_mat = pang_tree(csv, expected_format="panaroo")
    assert dist_mat["A", "B"] == pytest.approx(dist_mat["B", "A"], abs=1e-10)

def test_jaccard_ignores_shared_absences(tmp_path):
    """Sharing most of a repertoire must beat sharing nothing.

    Under Hamming the shared absences dominate a sparse matrix, so C and D,
    which have no gene in common, came out closer than A and B, which share
    eight of twelve. Jaccard ignores positions where both are zero.
    """
    genes = [f"g{i}" for i in range(12)]
    data = {
        "A": {g: i < 10 for i, g in enumerate(genes)},
        "B": {g: (i < 8 or i >= 10) for i, g in enumerate(genes)},
        "C": {g: i == 8 for i, g in enumerate(genes)},
        "D": {g: i == 9 for i, g in enumerate(genes)},
    }
    _, dm = pang_tree(write_panaroo_csv(tmp_path, data))

    assert dm["C", "D"] == pytest.approx(1.0)          # disjoint repertoires
    assert dm["A", "B"] == pytest.approx(1 - 8 / 12)   # 8 shared of 12 in union
    assert dm["A", "B"] < dm["C", "D"]


def test_identical_gene_content_is_zero(tmp_path):
    genes = [f"g{i}" for i in range(5)]
    data = {"A": {g: True for g in genes}, "B": {g: True for g in genes}}
    _, dm = pang_tree(write_panaroo_csv(tmp_path, data))
    assert dm["A", "B"] == pytest.approx(0.0)


def test_plasmids_without_genes_are_not_distant(tmp_path):
    """Samples absent from the pangenome table become all-zero columns."""
    genes = [f"g{i}" for i in range(4)]
    data = {"A": {g: False for g in genes}, "B": {g: False for g in genes}}
    _, dm = pang_tree(write_panaroo_csv(tmp_path, data))
    assert dm["A", "B"] == pytest.approx(0.0)
