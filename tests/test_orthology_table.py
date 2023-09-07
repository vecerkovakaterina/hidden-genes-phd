from hidden_genes_phd.OrthologyTable import OrthologyTable
from pathlib import Path
import pytest

orthology_table_csv = Path("orthology_table-68species_ensembl_ids.csv")
orthology_table_tsv = Path("orthology_table-68species_ensembl_ids.tsv")
orthology_table_excel = Path("orthology_table-68species_ensembl_ids.xlsx")
orthology_table_invalid_suffix = Path("orthology_table-68species_ensembl_ids.foo")
nonexistent_orthology_table = ""
tiny_orthology_table = Path("test_orthology_table.csv")


def test_class_attributes():
    ot = OrthologyTable(orthology_table_csv)
    assert ot.table_file == orthology_table_csv
    assert ot.orthology_df.empty is False


def test_parse_orthology_table():
    with pytest.raises(FileNotFoundError):
        OrthologyTable(nonexistent_orthology_table)

    ot_csv = OrthologyTable(orthology_table_csv)
    ot_tsv = OrthologyTable(orthology_table_tsv)
    ot_excel = OrthologyTable(orthology_table_excel)
    assert ot_csv.orthology_df.empty is False
    assert ot_tsv.orthology_df.empty is False
    assert ot_excel.orthology_df.empty is False

    correct_dims = (54120, 68)
    assert ot_csv.orthology_df.shape == correct_dims
    assert ot_tsv.orthology_df.shape == correct_dims
    assert ot_excel.orthology_df.shape == correct_dims

    with pytest.raises(ValueError):
        OrthologyTable(orthology_table_invalid_suffix)


def test_species_list():
    ot = OrthologyTable(orthology_table_csv)
    assert type(ot.get_species_list()) == list
    assert len(ot.get_species_list()) == 68


def test_annotation_names():
    ot = OrthologyTable(orthology_table_csv)
    assert type(ot.get_annotation_names()) == list
    assert len(ot.get_annotation_names()) == 68


def test_add_taxonomy_class_to_df():
    ot = OrthologyTable(tiny_orthology_table)
    gs = ot.create_genomes_for_species_in_table()
    ot = ot.add_taxonomy_class_to_df(gs)
    assert ot.orthology_taxonomy_df["class"].isna().sum() == 0
    new_df_shape = ot.orthology_df.shape[::-1]
    new_df_shape = list(new_df_shape)
    new_df_shape[1] += 1
    new_df_shape = tuple(new_df_shape)
    assert ot.orthology_taxonomy_df.shape == new_df_shape
