from pathlib import Path
from hidden_genes_phd.Annotation import Annotation

test_gtf = Path("genomes", "homo_sapiens.gtf")
non_existent_gtf = Path("genomes", "gnomus_hortensis.gtf")  # garden gnome in Latin


def test_instance_attributes():
    assert Annotation(test_gtf).gtf_file == test_gtf


def test_gtf_file_exists():
    assert Annotation(test_gtf).gtf_file_exists() is True
    assert Annotation(non_existent_gtf).gtf_file_exists() is False


def test_gtf_count_header_lines():
    assert Annotation("cavia_porcellus.gtf").gtf_count_header_lines() == 5
    assert Annotation("maylandia_zebra.gtf").gtf_count_header_lines() == 5
    assert Annotation("octodon_degus.gtf").gtf_count_header_lines() == 5


def test_gtf_to_df():
    an = Annotation(test_gtf)
    assert isinstance(an.gtf_to_df(), Annotation)
    assert an.df.is_empty() is False


def test_sort_df_by_coordinates():
    # TODO
    pass


def test_filter_only_genes():
    # TODO
    pass
