from hidden_genes_phd.OrthologyGroup import OrthologyGroup
from hidden_genes_phd.OrthologyTable import OrthologyTable
from hidden_genes_phd.Genome import Genome
from pathlib import Path
import pytest
import polars as pl

orthology_table_csv = Path("orthology_table-68species_ensembl_ids.csv")
orthology_table_tsv = Path("orthology_table-68species_ensembl_ids.tsv")
orthology_table_excel = Path("orthology_table-68species_ensembl_ids.xlsx")
orthology_table_invalid_suffix = Path("orthology_table-68species_ensembl_ids.foo")
nonexistent_orthology_table = ""
tiny_orthology_table = Path("test_orthology_table.csv")


def test_class_attributes():
    ot = OrthologyTable(orthology_table_csv)
    assert ot.table_file == orthology_table_csv
    assert ot.orthology_df.is_empty() is False


def test_parse_orthology_table():
    with pytest.raises(FileNotFoundError):
        OrthologyTable(nonexistent_orthology_table)

    ot_csv = OrthologyTable(orthology_table_csv)
    ot_tsv = OrthologyTable(orthology_table_tsv)
    ot_excel = OrthologyTable(orthology_table_excel)
    assert ot_csv.orthology_df.is_empty() is False
    assert ot_tsv.orthology_df.is_empty() is False
    assert ot_excel.orthology_df.is_empty() is False

    correct_dims = (54120, 68)
    assert ot_csv.orthology_df.shape == correct_dims
    assert ot_tsv.orthology_df.shape == correct_dims
    assert ot_excel.orthology_df.shape == correct_dims

    with pytest.raises(ValueError):
        OrthologyTable(orthology_table_invalid_suffix)


def test_species_list():
    ot = OrthologyTable(orthology_table_csv)
    assert type(ot.get_species_list()) == list
    assert len(ot.get_species_list()) == len(ot.orthology_df.columns)


def test_annotation_names():
    ot = OrthologyTable(orthology_table_csv)
    assert type(ot.get_annotation_names()) == list
    assert len(ot.get_annotation_names()) == len(ot.orthology_df.columns)


def test_add_taxonomy_class_to_df():
    ot = OrthologyTable(tiny_orthology_table)
    gs = ot.create_genomes_for_species_in_table()
    ot = ot.add_taxonomy_class_to_df(gs)
    assert ot.orthology_taxonomy_df.null_count().select("class").item() == 0
    new_df_shape = ot.orthology_df.shape[::-1]
    new_df_shape = list(new_df_shape)
    new_df_shape[1] += 2
    new_df_shape = tuple(new_df_shape)
    assert ot.orthology_taxonomy_df.shape == new_df_shape


def test_create_genomes_for_species_in_table():
    ot = OrthologyTable(tiny_orthology_table)
    gs = ot.create_genomes_for_species_in_table()
    assert len(gs) == ot.orthology_df.select(pl.count()).item()
    for name, obj in gs.items():
        assert isinstance(obj, Genome)


def test_create_orthology_groups():
    ot = OrthologyTable(tiny_orthology_table)
    ogs = ot.create_orthology_groups()
    assert len(ogs) == len(ot.orthology_df.columns)
    for og in ogs:
        assert isinstance(og, OrthologyGroup)
