from hidden_genes_phd.OrthologyTable import OrthologyTable
from pathlib import Path
import pytest

orthology_table_csv = Path("orthology_table-68species_ensembl_ids.csv")
orthology_table_tsv = Path("orthology_table-68species_ensembl_ids.tsv")
orthology_table_excel = Path("orthology_table-68species_ensembl_ids.xlsx")
orthology_table_invalid_suffix = Path("orthology_table-68species_ensembl_ids.foo")
nonexistent_orthology_table = ""


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
