import pytest

from hidden_genes_phd.HiddenGene import HiddenGene


def test_validate_side():
    assert HiddenGene.validate_side("left") == None
    assert HiddenGene.validate_side("right") == None
    with pytest.raises(ValueError) as ex:
        HiddenGene.validate_side("middle")
        assert str(ex.value) == "Side must be of ['left', 'right']"


def test_get_n_neighbors_side():
    # TODO
    pass


def test_most_common_neighbor_side():
    # TODO
    pass
