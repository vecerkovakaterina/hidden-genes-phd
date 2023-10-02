from hidden_genes_phd.OrthologyGroup import OrthologyGroup

mock_orthologs_list = ["A", "B", "C"]


def test_instance_attributes():
    og = OrthologyGroup(mock_orthologs_list)
    assert og.orthologs == mock_orthologs_list
    assert len(OrthologyGroup.orthology_groups_list) == 1
