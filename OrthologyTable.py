from pathlib import Path

import polars as pl

from hidden_genes_phd.Genome import Genome
from hidden_genes_phd.OrthologyGroup import OrthologyGroup

known_problematic_annotation_names = {
    "heterocephalus_glaber": "heterocephalus_glaber_female",
    "gorilla_gorilla_gorilla": "gorilla_gorilla",
    "cricetulus_griseus": "cricetulus_griseus_chok1gshd",
    "ovis_aries": "ovis_aries_rambouillet",
}


class OrthologyTable:
    def __init__(self, path_to_orthology_table):
        self.table_file = path_to_orthology_table
        self.orthology_df = self.parse_orthology_table()
        self.orthology_taxonomy_df = None

    @classmethod
    def file_exists(cls, file):
        return Path(file).is_file()

    def parse_orthology_table(self):
        if self.file_exists(self.table_file):
            suffix = Path(self.table_file).suffix
            if suffix == ".csv":
                return pl.read_csv(source=self.table_file, has_header=True).drop(
                    columns=""
                )
            elif suffix == ".tsv":
                return pl.read_csv(
                    source=self.table_file, separator="\t", has_header=True
                ).drop(columns="")
            elif suffix == ".xlsx":
                return pl.read_excel(
                    source=self.table_file, read_csv_options={"has_header": True}
                ).drop(columns="")
            else:
                raise ValueError(
                    "Invalid orthology table format. Accepted formats .csv, .tsv, .xlsx."
                )
        else:
            raise FileNotFoundError(f"File {self.table_file} does not exist.")

    def get_species_list(self):
        return self.orthology_df.columns

    def get_annotation_names(self):
        annotation_names = [
            "_".join(x.split(" ")).lower() for x in self.get_species_list()
        ]
        return annotation_names

    @classmethod
    def transpose_df(cls, df):
        return df.transpose(include_header=True)

    def add_taxonomy_class_to_df(self, genomes_dict):
        annotation_names = self.get_annotation_names()
        transposed_orthology_table = OrthologyTable.transpose_df(self.orthology_df)

        classes = [genomes_dict[name].get_species_class() for name in annotation_names]

        transposed_orthology_table = transposed_orthology_table.with_columns(
            pl.Series(name="class", values=classes)
        )

        self.orthology_taxonomy_df = transposed_orthology_table
        return self

    def create_genomes_for_species_in_table(self):
        annotation_names = self.get_annotation_names()

        for annotation in annotation_names:
            if annotation in known_problematic_annotation_names:
                Genome(annotation, known_problematic_annotation_names[annotation])
            else:
                Genome(annotation, annotation)

        return Genome.genomes_dict

    def create_orthology_groups(self):
        OrthologyGroup.max_orthologs = self.orthology_df.shape[1]
        for row in self.orthology_df.rows():
            OrthologyGroup(list(row))

        return OrthologyGroup.orthology_groups_list
