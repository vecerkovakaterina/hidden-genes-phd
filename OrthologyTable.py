from pathlib import Path

import pandas as pd

from Genome import Genome
from OrthologyGroup import OrthologyGroup

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
                return pd.read_csv(
                    self.table_file, sep=",", header=0, index_col=0, low_memory=False
                )
            elif suffix == ".tsv":
                return pd.read_csv(
                    self.table_file, sep="\t", header=0, index_col=0, low_memory=False
                )
            elif suffix == ".xlsx":
                return pd.read_excel(self.table_file, header=0, index_col=0)
            else:
                raise ValueError(
                    "Invalid orthology table format. Accepted formats .csv, .tsv, .xlsx."
                )
        else:
            raise FileNotFoundError(f"File {self.table_file} does not exist.")

    def get_species_list(self):
        return self.orthology_df.columns.tolist()

    def get_annotation_names(self):
        annotation_names = [
            "_".join(x.split(" ")).lower() for x in self.get_species_list()
        ]
        return annotation_names

    @classmethod
    def transpose_df(cls, df):
        return df.transpose()

    def add_taxonomy_class_to_df(self, genomes_dict):
        annotation_names = self.get_annotation_names()
        transposed_orthology_table = OrthologyTable.transpose_df(self.orthology_df)
        transposed_orthology_table.index = self.get_annotation_names()

        for annotation in annotation_names:
            transposed_orthology_table.loc[annotation, "class"] = genomes_dict[
                annotation
            ].get_species_class()

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
        for index, row in self.orthology_df.iterrows():
            OrthologyGroup(row.dropna().values.tolist())

        return OrthologyGroup.orthology_groups_list
