import gget
import pathlib
import urllib.request
import gzip
import shutil
import requests
import json


def access_ensembl_api(ext):
    server = "https://rest.ensembl.org"
    response = requests.get(server + ext, headers={"Content-Type": "application/json"})
    response = json.dumps(response.json())
    return response


class Genome:
    kind = "vertebrate"
    ensembl_release = 109

    genomes_dict = {}

    def __init__(self, species_name, annotation_name):
        self.taxonomy_class = None
        self.species_name = species_name
        self.annotation_name = annotation_name
        self.ftp_link = None
        self.genome_file = pathlib.Path("genomes/", self.species_name + ".gtf")
        self.download_annotation()
        Genome.genomes_dict[self.species_name] = self

    @classmethod
    def unpack_gz(cls, gzip_file, output_file):
        with gzip.open(gzip_file, "rb") as packed:
            with open(output_file, "wb") as unpacked:
                shutil.copyfileobj(packed, unpacked)

    def annotation_name_exists(self):
        try:
            gget.ref(self.annotation_name, which=["gtf"], release=self.ensembl_release)[
                self.annotation_name
            ]["annotation_gtf"]["ftp"]
        except ValueError:
            raise ValueError(f"Invalid annotation name: {self.annotation_name}")
        else:
            return True

    def annotation_downloaded(self):
        return pathlib.Path("genomes", self.annotation_name + ".gtf").exists()

    def download_annotation(self):
        if self.annotation_name_exists():
            self.ftp_link = gget.ref(
                self.annotation_name, which=["gtf"], release=self.ensembl_release
            )[self.annotation_name]["annotation_gtf"]["ftp"]

        if not self.annotation_downloaded():
            gtf_gz = pathlib.Path("genomes", self.annotation_name + ".gtf.gz")
            urllib.request.urlretrieve(self.ftp_link, gtf_gz)
            self.unpack_gz(gtf_gz, self.genome_file)

    def get_species_class(self):
        ens_api_response = access_ensembl_api(
            f"/taxonomy/classification/{self.species_name}?"
        )
        classes = ["Mammalia", "Aves", "Reptilia", "Actinopteri", "Amphibia"]

        for taxon in classes:
            if taxon in ens_api_response:
                self.taxonomy_class = taxon
                return taxon
