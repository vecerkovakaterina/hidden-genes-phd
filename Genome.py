import gget
import pathlib
import urllib.request
import gzip
import shutil


class Genome:
    kind = "vertebrate"
    ensembl_release = 109

    genomes_dict = {}

    def __new__(cls, species):
        if not cls.assembly_name_exists(species):
            return None

        instance = super().__new__(cls)
        instance.species = species
        instance.ftp_link = ""
        instance.genome_file = pathlib.Path("genomes/", instance.species + ".gtf")
        instance.download_annotation()
        Genome.genomes_dict[instance.species] = instance

        return instance

    # def __init__(self, species):
    #     self.species = species
    #     self.ftp_link = ""
    #     self.genome_file = pathlib.Path("genomes/", self.species + ".gtf")
    #
    #     self.download_annotation()
    #     Genome.genomes_dict[self.species] = self

    def unpack_gz(self, gzip_file):
        with gzip.open(gzip_file, "rb") as packed:
            with open(self.genome_file, "wb") as unpacked:
                shutil.copyfileobj(packed, unpacked)

    @classmethod
    def assembly_name_exists(cls, species):
        try:
            gget.ref(species, which=["gtf"], release=cls.ensembl_release)[species][
                "annotation_gtf"
            ]["ftp"]
        except ValueError:
            print(f"Invalid assembly name: {species}")
            return False
        else:
            return True

    def download_annotation(self):
        self.ftp_link = gget.ref(
            self.species, which=["gtf"], release=self.ensembl_release
        )[self.species]["annotation_gtf"]["ftp"]

        genome_file_gz = pathlib.Path("genomes", self.species + ".gtf.gz")
        urllib.request.urlretrieve(self.ftp_link, genome_file_gz)
        self.unpack_gz(genome_file_gz)
