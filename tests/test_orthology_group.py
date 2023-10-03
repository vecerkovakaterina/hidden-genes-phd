from pathlib import Path
import pandas as pd
from hidden_genes_phd.OrthologyGroup import OrthologyGroup, get_max_possible_score
from hidden_genes_phd.OrthologyTable import OrthologyTable

orthology_table_file = Path("test_orthology_table.csv")
orthology_table = OrthologyTable(orthology_table_file)

orthology_taxonomy_df = pd.DataFrame(
    [["A1", "B1"], ["A2", None], ["Mammalia", "Aves"]]
).T

orthology_table.orthology_taxonomy_df = orthology_taxonomy_df
orthology_table.orthology_taxonomy_df.columns = [0, 1, "class"]

orthology_group_1 = OrthologyGroup(["A1", "B1"])
orthology_group_2 = OrthologyGroup(["B2"])

target_taxon = "Mammalia"


def test_instance_attributes():
    assert orthology_group_1.orthologs == ["A1", "B1"]
    assert len(OrthologyGroup.orthology_groups_list) == 2


def test_score_without_taxon():
    orthology_group_1.calculate_score(orthology_table)
    assert orthology_group_1.score == 1.0
    assert orthology_group_2.score == 0.5


def test_score_with_taxon():
    orthology_group_1.calculate_score(orthology_table, target_taxon)
    assert orthology_group_1.score == 2.0
    assert orthology_group_2.score == 1.5


def test_get_max_possible_score_without_taxon():
    assert get_max_possible_score(orthology_table) == 1.0


def test_get_max_possible_score_with_taxon():
    assert get_max_possible_score(orthology_table, target_taxon) == 2.0


def test_rank_groups_by_score():
    OrthologyGroup.rank_groups_by_score()
    expected_scores = [2.0, 1.5]
    actual_scores = [x.score for x in OrthologyGroup.orthology_groups_list]
    assert actual_scores == expected_scores


def test_delete_full_groups_without_taxon():
    orthology_group_1.calculate_score(orthology_table)
    orthology_group_1.delete_full_groups(orthology_table)
    assert len(OrthologyGroup.orthology_groups_list) == 1


def test_delete_full_groups_with_taxon():
    orthology_group_1.calculate_score(orthology_table)
    orthology_group_1.delete_full_groups(orthology_table)
    assert len(OrthologyGroup.orthology_groups_list) == 1
