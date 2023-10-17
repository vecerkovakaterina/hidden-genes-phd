from pathlib import Path
import polars as pl
from hidden_genes_phd.OrthologyGroup import OrthologyGroup
from hidden_genes_phd.OrthologyTable import OrthologyTable

orthology_table_file = Path("test_orthology_table.csv")
orthology_table = OrthologyTable(orthology_table_file)

orthology_taxonomy_df = pl.DataFrame([["A1", "B1"], ["A2", None], ["Mammalia", "Aves"]])

orthology_table.orthology_taxonomy_df = orthology_taxonomy_df
orthology_table.orthology_taxonomy_df.columns = ["0", "1", "class"]

orthology_group_1 = OrthologyGroup(["A1", "B1"], ["homo_sapiens", "gallus_gallus"])
orthology_group_2 = OrthologyGroup(["nan", "B2"], ["gallus_gallus"])

target_taxon = "Mammalia"


def test_instance_attributes():
    assert orthology_group_1.orthologs == ["A1", "B1"]
    assert len(OrthologyGroup.orthology_groups_list) == 2


def test_calculate_score_without_taxon():
    orthology_group_1.calculate_score_all_groups(orthology_table)
    assert orthology_group_1.score == 1.0
    assert orthology_group_2.score == 0.5


def test_calculate_score_with_taxon():
    orthology_group_1.calculate_score_all_groups(orthology_table, target_taxon)
    assert orthology_group_1.score == 2.0
    assert orthology_group_2.score == 1.5


def test_get_max_possible_score_without_taxon():
    assert OrthologyGroup.get_max_possible_score(orthology_table) == 1.0


def test_get_max_possible_score_with_taxon():
    assert OrthologyGroup.get_max_possible_score(orthology_table, target_taxon) == 2.0


def test_rank_groups_by_score():
    OrthologyGroup.rank_groups_by_score()
    expected_scores = [2.0, 1.5]
    actual_scores = [x.score for x in OrthologyGroup.orthology_groups_list]
    assert actual_scores == expected_scores


def test_delete_full_groups_without_taxon():
    orthology_group_1.calculate_score_all_groups(orthology_table)
    orthology_group_1.delete_full_groups(orthology_table)
    assert len(OrthologyGroup.orthology_groups_list) == 1


def test_delete_full_groups_with_taxon():
    orthology_group_1.calculate_score_all_groups(orthology_table)
    orthology_group_1.delete_full_groups(orthology_table)
    assert len(OrthologyGroup.orthology_groups_list) == 1


# test static methods
def test_drop_nans_from_list():
    list_wo_nans = OrthologyGroup.drop_nans_from_list(
        ["nan", "nun", "nen", 0, "nan", "nan"]
    )
    assert list_wo_nans == ["nun", "nen", 0]


def test_taxon_score():
    assert orthology_group_1.taxon_score(1, orthology_table, "Mammalia") == 1
    assert orthology_group_1.taxon_score(1, orthology_table, "Aves") == 0
    assert orthology_group_2.taxon_score(2, orthology_table, "Mammalia") == 1
    assert orthology_group_2.taxon_score(2, orthology_table, "Aves") == 1
