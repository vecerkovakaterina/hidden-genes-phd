import polars as pl
import gget
from pathlib import Path
import subprocess


def drop_nans_form_list(lst):
    return [x for x in lst if x != "nan"]


def taxon_score(col, orthology_table, taxon):
    nans_in_taxon = (
        orthology_table.orthology_taxonomy_df.filter(pl.col("class") == taxon)
        .select(str(orthology_table.orthology_taxonomy_df.columns[col]))
        .null_count()
        .item()
    )
    species_in_taxon = (
        orthology_table.orthology_taxonomy_df.filter(pl.col("class") == taxon)
        .select(str(orthology_table.orthology_taxonomy_df.columns[col]))
        .select(pl.count())
        .item()
    )
    orthologs_in_taxon = species_in_taxon - nans_in_taxon
    return orthologs_in_taxon


def get_max_possible_score(orthology_table, taxon=None):
    full_orthology_group = 1.0
    if taxon:
        full_orthology_group += (
            orthology_table.orthology_taxonomy_df.filter(pl.col("class") == taxon)
            .select(pl.count())
            .item()
        )
    return full_orthology_group


def get_max_size_orthology_group(orthology_table):
    return len(orthology_table.orthology_taxonomy_df.columns) - 1


class OrthologyGroup:
    orthology_groups_list = []

    def __init__(self, orthologs):
        self.orthologs = orthologs
        self.score = None
        self.sequences_fasta = None
        OrthologyGroup.orthology_groups_list.append(self)

    def calculate_score(self, orthology_table, taxon=None):
        max_orthologs = get_max_size_orthology_group(orthology_table)

        for i, og in enumerate(self.orthology_groups_list):
            number_orthologs = len([x for x in og.orthologs if str(x) != "nan"])
            og.score = number_orthologs / max_orthologs
            if taxon:
                og.score += taxon_score(i, orthology_table, taxon)

    def delete_full_groups(self, orthology_table):
        max_orthologs = get_max_size_orthology_group(orthology_table)
        for og in self.orthology_groups_list:
            if len(og.orthologs) == max_orthologs:
                self.orthology_groups_list.remove(og)

    @classmethod
    def rank_groups_by_score(cls):
        cls.orthology_groups_list.sort(key=lambda x: x.score, reverse=True)

    def create_ortholog_sequences_list(self):
        orthologs = drop_nans_form_list(self.orthologs)
        ortholog_aa_sequences = []
        for ortholog in orthologs:
            try:
                sequence = gget.seq(ortholog, translate=True, verbose=False)
                if len(sequence) > 0:
                    ortholog_aa_sequences.append(sequence)
            except TypeError:
                continue

        return ortholog_aa_sequences

    def write_ortholog_sequences_to_fasta(self):
        aa_sequences = self.create_ortholog_sequences_list()
        aa_sequences = "\n".join(["\n".join(sequence) for sequence in aa_sequences])
        aa_sequences_file = Path(
            "ortholog_sequences",
            f"{'_'.join(drop_nans_form_list(self.orthologs)) + '.fa'}",
        )

        with open(aa_sequences_file, "w") as fasta_file:
            fasta_file.write(aa_sequences)

        self.sequences_fasta = aa_sequences_file

    def fasta_to_custom_blast_db(self):
        command = f"makeblastdb -in {self.sequences_fasta} -parse_seqids -dbtype prot"
        result = subprocess.run(
            command,
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=True,
            text=True,
        )
        if result.returncode != 0:
            raise Exception(
                f"Error when creating custom BLAST db: {self.sequences_fasta}!"
            )
