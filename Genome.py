class Genome:

    kind = "vertebrate"
    ensembl_release = 109

    def __init__(self, species):
        self.species = species

    def download_anotation(self, ftp_link):
        pass