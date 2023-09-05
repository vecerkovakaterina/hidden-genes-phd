from hidden_genes_phd.Genome import Genome
from pathlib import Path

valid_annotation_name = "heterocephalus_glaber_female"
invalid_annotation_name = "heterocephalus_glaber"


def test_class_attributes():
    g = Genome(valid_annotation_name)
    assert g.kind == "vertebrate"
    assert g.ensembl_release == 109


def test_class_attribute_genomes_dict():
    assert len(Genome.genomes_dict) == 1
    assert "heterocephalus_glaber_female" in Genome.genomes_dict


def test_instance_attribute_species():
    g = Genome(valid_annotation_name)
    assert g.species == valid_annotation_name


def test_instance_attribute_ftp_link():
    g = Genome(valid_annotation_name)
    assert (
        g.ftp_link
        == "http://ftp.ensembl.org/pub/release-109/gtf/heterocephalus_glaber_female/Heterocephalus_glaber_female.HetGla_female_1.0.109.gtf.gz"
    )


def test_instance_attribute_genome_file():
    g = Genome(valid_annotation_name)
    assert Path("genomes", g.species + ".gtf.gz").exists() is True


def test_genome_annotation_unpacked():
    g = Genome(valid_annotation_name)
    assert Path("genomes", g.species + ".gtf").exists() is True


def test_class_method_annotation_name_exists():
    assert Genome.annotation_name_exists(valid_annotation_name) is True
    assert Genome.annotation_name_exists(invalid_annotation_name) is False
